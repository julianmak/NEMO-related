#!/usr/bin/env python3

# JM: 03 Sep 2018
# plots some meridional sections
# note that "meridional" here means plot along a fixed ORCA index, which is 
# slightly curved towards the Arctic

import matplotlib
matplotlib.use('agg')

from orca_plotting_commands import *
from midpointnorm import *

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

from pyCDFTOOLS import eos

# style settings

plt.rcParams["font.family"] = "DejaVu Serif"
plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams["mathtext.rm"] = "serif"
plt.rcParams["image.cmap"] = "RdBu_r" # "*_r" is reverse of standard colour

#--------------------------------------------------------
# define the argument parser
import argparse

parser = argparse.ArgumentParser(description = "Plot the MOC averaged at fixed height")

# fixed arguments
parser.add_argument("data_dir", type = str, 
                    help = "specify data directory")
parser.add_argument("fileT",    type = str, 
                    help = "specify data filename")

# optional arguments
parser.add_argument("--t_var",  type = str, 
                    help = "specify temperature variable name (default = thetao)", default = "thetao")
parser.add_argument("--s_var",  type = str, 
                    help = "specify salinity variable name (default = so)", default = "so")
parser.add_argument("--lprint", 
                    help = "print out the variables available in fileT", action = "store_true")
parser.add_argument("--lon",  type = float, 
                    help = "specify the target longitude to plot (default = 0.0)", default = 0.0)

# collect arguments
args = parser.parse_args()

#--------------------------------------------------------
# load the data

data   = Dataset(args.data_dir + args.fileT)
depth  =-data.variables["deptht"][:]
lonT   = data.variables["nav_lon"][:, :]
latT   = data.variables["nav_lat"][:, :]
thetaT = data.variables["thetao"][0, :, :, :]
salinT = data.variables["so"][0, :, :, :]

if args.lprint:
  print(data)
  data.close()
  sys.exit("printing variables available in file, exiting...")
  
data.close()

#mask out land points
tmask = (salinT == 0)
salinT[tmask] = np.nan
thetaT[tmask] = np.nan

#--------------------------------------------------------
# plotting

# find the equator index and make that the reference
eq_index = np.where(latT[:, 0] == 0)[0][0]

if args.lon > 0:
    lon_str = "%g E" % abs(args.lon)
else:
    lon_str = "%g W" % abs(args.lon)
dummy = lonT[eq_index, :] - args.lon
lon_index = np.where(dummy == np.amin(abs(dummy)))[0][0]

fig = plt.figure(figsize = (8, 16))

plt.subplot(3, 1, 1)
mesh = plt.contourf(latT[:, lon_index], depth, thetaT[:, :, lon_index],
                    np.arange(-2, 30, 1), cmap = "RdBu_r", extend = "both")
plt.gca().set_facecolor('gray')
line = plt.contour(latT[:, lon_index], depth, thetaT[:, :, lon_index], levels = [0, 10, 20, 30],
                   colors = "k")
plt.clabel(line, fmt = "%2.0f", colors = 'k')
line = plt.contour(latT[:, lon_index], depth, thetaT[:, :, lon_index], levels = [4],
                   colors = "k", linestyles = "--")
plt.clabel(line, fmt = "%2.1f", colors = 'k')
plt.xlim(-75, 60)
plt.ylabel(r"$z$ ($\mathrm{m}$)")
plt.title(r"temeprature at Lon = %s" % lon_str)
cb = plt.colorbar(mesh)
cb.ax.set_title(r"${}^\circ\ \mathrm{C}$")

plt.subplot(3, 1, 2)
mesh = plt.contourf(latT[:, lon_index], depth, salinT[:, :, lon_index],
                    np.arange(33.5, 37.5, 0.1), cmap = "gist_rainbow_r", extend = "both")
plt.gca().set_facecolor('gray')
line = plt.contour(latT[:, lon_index], depth, salinT[:, :, lon_index], levels = [34, 35, 36],
                   colors = "k")
plt.clabel(line, fmt = "%2.0f", colors = 'k')
line = plt.contour(latT[:, lon_index], depth, salinT[:, :, lon_index], levels = [34.5],
                   colors = "k", linestyles = "--")
plt.clabel(line, fmt = "%2.2f", colors = 'k')
plt.ylabel(r"$z$ ($\mathrm{m}$)")
plt.xlim(-75, 60)
plt.title(r"salinity at Lon = %s" % lon_str)
cb = plt.colorbar(mesh)
cb.ax.set_title(r"$\mathrm{g}\ \mathrm{kg}^{-1}$")

# compute the relevant sigma2 on the fly
# potential density referenced to 2000m depth
sigma2 = eos.sigmai_dep(thetaT[:, :, lon_index], salinT[:, :, lon_index], 2000)

plt.subplot(3, 1, 3)
mesh = plt.contourf(latT[:, lon_index], depth, sigma2,
                    np.arange(30, 37.5, 0.25), cmap = "gist_rainbow_r", extend = "both", norm = MidPointNorm(midpoint = 35))
plt.gca().set_facecolor('gray')
line = plt.contour(latT[:, lon_index], depth, sigma2, levels = [36.9],
                   colors = "k")
plt.clabel(line, fmt = r"$\sigma_2$ = %2.1f", colors = 'k')
plt.xlim(-75, 60)
plt.xlabel(r"Lat (${}^\circ$)")
plt.ylabel(r"$z$ ($\mathrm{m}$)")
plt.title(r"$\sigma_2 - 1000$ at Lon = %s" % lon_str)
cb = plt.colorbar(mesh)
cb.ax.invert_yaxis()
cb.ax.set_title(r"$\mathrm{kg}\ \mathrm{m}^{-3}$")
  
#--------------------------------------------------------
# saving commands

save_filename = args.fileT.replace(".nc", "") + "_lat_slice.png"

fig.savefig(save_filename, dpi = 300, bbox_inches = "tight")
plt.close(fig)

print("generated %s , exiting..." % save_filename)
