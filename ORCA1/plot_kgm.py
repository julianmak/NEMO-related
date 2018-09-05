#!/usr/bin/env python3

# JM: 05 Sep 2018
# plots the kapgm variable (aeiw in NEMO)
# assumes it's on a W grid
# styling is default and this script is intended to be used for quick and dirty 
# visualisations

import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt

from numpy import sum, newaxis, maximum, nan, linspace, log10, arange
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

parser = argparse.ArgumentParser(description = """Plot the z averaged kapgm variable with Plate Carree projection
                                                  (needs Iris and Cartopy package)""")

# fixed arguments
parser.add_argument("data_dir", type = str, 
                    help = "specify data directory")
parser.add_argument("fileW",    type = str, 
                    help = "specify data filename")

# optional arguments
parser.add_argument("--w_var", type = str,
                    help = "specify w variable name (default = aeiw)", default = "aeiw")
parser.add_argument("--lprint", 
                    help = "print out the variables available in fileU", action = "store_true")
parser.add_argument("--kt", type = int,
                    help = "plot a specified time slice (default = last time entry in variable)")
parser.add_argument("--level", nargs = 2, type = float,
                    help = "specify the limits of levels to plot if any (default: 10 to 3000)")
parser.add_argument("--cshift", type = float,
                    help = "specify a shift of the data to emphasis high/low values (set to low value to emphasise high colours)")
                    
parser.add_argument("--file_out",  type = str, 
                    help = "specify output name (default = fileW + _kgm.png)") 
                    
# collect arguments
args = parser.parse_args()

if args.level is None:
  args.level = []
  args.level.append(1e2)
  args.level.append(3e3)

#--------------------------------------------------------
# load files if necessary
data_netcdf4 = netCDF4.Dataset(args.data_dir + args.fileW)
if args.lprint:
  print(data_netcdf4)
  data_netcdf4.close()
  sys.exit("finished displaying data in file, exiting...")
if args.kt is None:
  kt = data_netcdf4.dimensions["time_counter"].size - 1 # default load the last time level
kgm_raw = data_netcdf4.variables[args.w_var][kt, :, :, :]
depthW = data_netcdf4.variables["depthw"][:]
latW = data_netcdf4.variables["nav_lat"][:, :]
lonW = data_netcdf4.variables["nav_lon"][:, :]
e3W = data_netcdf4.variables["e3w"][kt, :, :, :]
data_netcdf4.close()

data_netcdf4 = netCDF4.Dataset(args.data_dir + "mesh_mask.nc")
tmask = data_netcdf4.variables["tmask"][0, :, : ,:]
data_netcdf4.close()

#--------------------------------------------------------
# Main plotting commands

# process and projection step

depth = sum(depthW[:, newaxis, newaxis] * tmask, axis = 0) # broadcast the 1d array to 3d

kgm_zavg = kgm_raw[0, :, :]
#kgm_zavg = sum(kgm_raw * e3W * tmask, axis = 0) / maximum(depth, 1e-16)

iris.FUTURE.netcdf_promote = True

pcarree = ccrs.PlateCarree()
target_proj = pcarree

lat = iris.coords.AuxCoord(latW, standard_name = "latitude", units = "degrees")
lon = iris.coords.AuxCoord(lonW, standard_name = "longitude", units = "degrees")
data_cube = iris.cube.Cube(kgm_zavg, 
                           long_name = "kgm_depth_avg", 
                           units = "m s-2",
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
cb.ax.set_title(r"$\mathrm{m}\ \mathrm{s}^{-2}$") # title on colourbar

clevel = [10, 100, 500]
for i in range(1000, int(args.level[1] + 1), 1000):
    clevel.append(i)

cb.set_ticks(log10(clevel))
cb.set_ticklabels(clevel)
  
#--------------------------------------------------------
# saving commands

if args.file_out is None:
  args.file_out = args.fileW.replace(".nc", "") + "_kgm.png"

fig.savefig(args.file_out, dpi = 300, bbox_inches = "tight")
plt.close(fig)

print("generated %s , exiting..." % args.file_out)
