#!/usr/bin/env python3

# JM: 05 Sep 2018
# plot the eke variable in log scale
# (designed for the GEOMETRIC depth-integrated eddy energy but ok for NEMO 
#  generated too probably)
# styling is default and this script is intended to be used for quick and dirty 
# visualisations

import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt

from numpy import maximum, sum, nan, linspace, log10, arange, newaxis
from orca_plotting_commands import *
from midpointnorm import *

import netCDF4, sys
# need iris if wanting to use cartopy_command
import iris
import iris.analysis.cartography

# for editing
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

# style settings

plt.rcParams["font.family"] = "DejaVu Serif"
plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams["mathtext.rm"] = "serif"
plt.rcParams["image.cmap"] = "RdBu_r" # "*_r" is reverse of standard colour

#--------------------------------------------------------
# define the argument parser
import argparse

parser = argparse.ArgumentParser(description = """Plot the eke variable in log scale with Plate Carree projection
                                                  (needs Iris and Cartopy package)""")

# fixed arguments
parser.add_argument("data_dir", type = str, 
                    help = "specify data directory")
parser.add_argument("fileT",    type = str, 
                    help = "specify data filename")

# optional arguments
parser.add_argument("--t_var", type = str,
                    help = "specify T variable name (default = eke)", default = "eke")
parser.add_argument("--lprint", 
                    help = "print out the variables available in fileT", action = "store_true")
parser.add_argument("--kt", type = int,
                    help = "plot a specified time slice (default = last time entry in variable)")
parser.add_argument("--level", nargs = 2, type = float,
                    help = "specify the limits of levels to plot if any (default: 1000 to 100000)")
parser.add_argument("--cshift", type = float,
                    help = "specify a shift of the data to emphasis high/low values (set to low value to emphasise high colours)")
                    
parser.add_argument("--file_out",  type = str, 
                    help = "specify output name (default = fileT + _eke.png)") 
                    
# collect arguments
args = parser.parse_args()

if args.level is None:
  args.level = []
  if args.lzavg:
    args.level.append(1e-4)
    args.level.append(1e-1)
  else:
    args.level.append(1e3)
    args.level.append(1e6)

#--------------------------------------------------------
# load files if necessary
data_netcdf4 = netCDF4.Dataset(args.data_dir + args.fileT)
if args.lprint:
  print(data_netcdf4)
  data_netcdf4.close()
  sys.exit("finished displaying data in file, exiting...")
if args.kt is None:
  kt = data_netcdf4.dimensions["time_counter"].size - 1 # default load the last time level
eE_raw = data_netcdf4.variables[args.t_var][kt, :, :]
latT = data_netcdf4.variables["nav_lat"][:, :]
lonT = data_netcdf4.variables["nav_lon"][:, :]
depthT = data_netcdf4.variables["deptht"][:]
e3T = data_netcdf4.variables["e3t"][kt, :, :, :]
data_netcdf4.close()

data_netcdf4 = netCDF4.Dataset(args.data_dir + "mesh_mask.nc")
tmask = data_netcdf4.variables["tmask"][0, :, : ,:]
data_netcdf4.close()

#--------------------------------------------------------
# Main plotting commands

# process and projection step

depth = sum(depthT[:, newaxis, newaxis] * tmask, axis = 0)

#eE_raw /= maximum(depth, 1e-16)
rho0 = 1024.0
eE_raw *= rho0

iris.FUTURE.netcdf_promote = True

pcarree = ccrs.PlateCarree()
target_proj = pcarree

lat = iris.coords.AuxCoord(latT, standard_name = "latitude", units = "degrees")
lon = iris.coords.AuxCoord(lonT, standard_name = "longitude", units = "degrees")
data_cube = iris.cube.Cube(eE_raw, 
                           long_name = "geom_eke_depth_avg", 
                           units = "m2 s-2",
                           aux_coords_and_dims = [(lat, (0, 1)), (lon, (0,1))])
data_proj, extent = iris.analysis.cartography.project(data_cube[:, :], pcarree, nx = 600, ny = 300)
x = data_proj.coord('projection_x_coordinate').points
y = data_proj.coord('projection_y_coordinate').points
plot_data = data_proj.data

# plot

plot_data[(plot_data == 0)] = nan

fig = plt.figure(figsize=(12, 7))

misc_format = {"levels" : linspace(log10(args.level[0]), log10(args.level[1]), 21),
               "extend" : "both",
               "cmap"   : "Spectral_r"}
               
if args.cshift is not None:
  misc_format["norm"] = MidPointNorm(midpoint = log10(args.cshift))

ax, mesh = cartopy_contourf(x, y, log10(plot_data), proj = target_proj, **misc_format)

ax.set_extent([-180, 180, -75, 80], crs = ccrs.PlateCarree())

# set axes, title and add gridlines

gl = ax.gridlines(crs=ccrs.PlateCarree(),
              draw_labels = True, linewidth = 1, linestyle = '--')
gl.xlabels_top = False
gl.ylabels_right = False
gl.ylocator = mpl.ticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

ax.text(-0.05, 0.5, r'Lat $\left( {}^\circ \right)$', 
        va='bottom', ha='center',
        rotation=90, rotation_mode='anchor',
        transform=ax.transAxes)
ax.text(0.5, -0.1, r'Lon $\left( {}^\circ \right)$', 
        va='bottom', ha='center',
        rotation='horizontal', rotation_mode='anchor',
        transform=ax.transAxes)

# add colorbar
divider = make_axes_locatable(ax)
ax_cb = divider.new_horizontal(size="1%", pad=0.2, axes_class=plt.Axes)
fig.add_axes(ax_cb)
cb = plt.colorbar(mesh, cax=ax_cb, orientation = "vertical")
cb.ax.set_title(r"$\mathrm{J}\ \mathrm{m}^{-2}$")
  
cb.set_ticks(arange(int(log10(args.level[0])), int(log10(args.level[1])) + 1, 1))
cb.set_ticklabels([r"$10^{%s}$" % i for i in range(int(log10(args.level[0])), int(log10(args.level[1])) + 1, 1)])
  
#--------------------------------------------------------
# saving commands

if args.file_out is None:
  args.file_out = args.fileT.replace(".nc", "") + "_eke.png"

fig.savefig(args.file_out, dpi = 300, bbox_inches = "tight")
plt.close(fig)

print("generated %s , exiting..." % args.file_out)
