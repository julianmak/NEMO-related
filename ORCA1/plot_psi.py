#!/usr/bin/env python3

# JM: 03 Sep 2018
# plots barotropic streamfunction using iris projection
# needs the Iris and/or Cartopy

import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt

import numpy as np
from orca_plotting_commands import *

import netCDF4
# need iris if wanting to use cartopy_command
import iris
import iris.analysis.cartography

# for editing
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from pyCDFTOOLS.cdfpsi import *

# style settings

plt.rcParams["font.family"] = "DejaVu Serif"
plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams["mathtext.rm"] = "serif"
plt.rcParams["image.cmap"] = "RdBu_r" # "*_r" is reverse of standard colour

#--------------------------------------------------------
# define the argument parser
import argparse

parser = argparse.ArgumentParser(description = """Plot the barotropic streamfunction with Plate Carree 
                                                  (needs Iris and Cartopy package)""")

# fixed arguments
parser.add_argument("data_dir", type = str, 
                    help = "specify data directory")
parser.add_argument("fileU",    type = str, 
                    help = "specify data filename")

# optional arguments
parser.add_argument("--u_var", type = str,
                    help = "specify u variable name (default = uo)", default = "uo")
parser.add_argument("--lprint", 
                    help = "print out the variables available in fileU", action = "store_true")
parser.add_argument("--ll_v",   
                    help = "use the v variable instead", action = "store_true")
parser.add_argument("--v_var",
                    help = "specify v variable name (default = vo)", default = "vo")
parser.add_argument("--lg_vvl",   
                    help = "use time varying metric e3u/v", action = "store_true")
parser.add_argument("--level", nargs = 2, type = float,
                    help = "specify the limits of levels to plot")
                    
# collect arguments
args = parser.parse_args()

args.fileV = args.fileU.replace("_U", "_V")

#--------------------------------------------------------
# Main plotting commands

kwargs = {"lprint" : args.lprint,
          "ll_v"   : args.ll_v,
          "lg_vvl" : args.lg_vvl}

if args.ll_v:
  print(" ")
  print("using v variable %s" % args.v_var)
else:
  print(" ")
  print("using u variable %s" % args.u_var)
  
if args.lg_vvl:
  print(" ")
  print("using time varying metric")
  
if args.level is None:
  args.level = []
  args.level.append(-60)
  args.level.append(190)
  
lonT, latT, psi, opt_dic = cdfpsi(args.data_dir, args.fileU, args.u_var, args.fileV, args.v_var, **kwargs)

# projection step
iris.FUTURE.netcdf_promote = True

pcarree = ccrs.PlateCarree()
target_proj = pcarree

lat = iris.coords.AuxCoord(latT, standard_name = "latitude", units = "degrees")
lon = iris.coords.AuxCoord(lonT, standard_name = "longitude", units = "degrees")
data_cube = iris.cube.Cube(psi, 
                           long_name = "Psi", 
                           units = "Sv",
                           aux_coords_and_dims = [(lat, (0, 1)), (lon, (0,1))])
data_proj, extent = iris.analysis.cartography.project(data_cube[:, :], pcarree, nx = 400, ny = 200)
x = data_proj.coord('projection_x_coordinate').points
y = data_proj.coord('projection_y_coordinate').points
plot_data = data_proj.data / 1e6 # in Sv

# plot

# touch up the data and set levels
# "mask" the data by setting the land point data to nan

plot_data[(plot_data == 0)] = np.nan

fig = plt.figure(figsize=(10, 7))

# plot denser colour contours
misc_format = {"levels" : np.arange(args.level[0], args.level[1], 10),
               "extend" : "both"}
ax, mesh = cartopy_contourf(x, y, plot_data, proj = target_proj, **misc_format)

# plot more sparse line contours
misc_format = {"levels" : np.arange(args.level[0], args.level[1], 30),
               "extend" : "both"}
line = plt.contour(x, y, plot_data, transform = pcarree, colors = "k", **misc_format)
plt.clabel(line, fmt = r"%2.0f", colors = 'k')

# set axes, title and add gridlines

ax.gridlines(linewidth = 1, linestyle = '--')
gl = ax.gridlines(crs=ccrs.PlateCarree(),
              draw_labels = True, linewidth = 1, linestyle = '--')
gl.xlabels_top = False
gl.ylabels_right = False
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
cb.ax.set_title(r'$\mathrm{Sv}$') # title on colourbar
  
#--------------------------------------------------------
# saving commands

save_filename = args.fileU.replace(".nc", "") + "_psi.png"

fig.savefig(save_filename, dpi = 300, bbox_inches = "tight")
plt.close(fig)

print("generated %s , exiting..." % save_filename)
