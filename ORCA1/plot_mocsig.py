#!/usr/bin/env python3

# JM: 04 Sep 2018
# plots the processed MOC averaged at fixed density in density space
# (use plot_mocsig_isodep for plotting the MOC mapped back to depth space)
# styling is default and this script is intended to be used for quick and dirty 
# visualisations

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from numpy import arange, amax, amin, where
import sys
import netCDF4

# style settings

plt.rcParams["font.family"] = "DejaVu Serif"
plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams["mathtext.rm"] = "serif"
plt.rcParams["image.cmap"] = "RdBu_r" # "*_r" is reverse of standard colour

#--------------------------------------------------------
# define the argument parser
import argparse

parser = argparse.ArgumentParser(description = "Plot the MOC averaged at fixed density in density space")

# fixed arguments                    
parser.add_argument("data_dir", type = str, 
                    help = "specify data directory")
parser.add_argument("fileMOC", type = str, 
                    help = "specify data filename")

# optional arguments
parser.add_argument("--MOC_var", type = str,
                    help = "specify MOC variable name (default = rmoc_glob_tave)", default = "rmoc_glob_tave")
parser.add_argument("--region", type = str,
                    help = """specify the region of interest: global (default),
                                                              so,
                                                              atlantic """, default = "global")
parser.add_argument("--lprint", 
                    help = "print out the variables available in fileMOC", action = "store_true")
parser.add_argument("--clim", nargs = 2, type = float,
                    help = "specify the MOC range to plot over if any (default: -30 to 30)")
                    
parser.add_argument("--file_out",  type = str, 
                    help = "specify output name (default = fileMOC + _sigmaMOC.png)")

# collect arguments
args = parser.parse_args()

if args.clim is None:
  args.clim = []
  args.clim.append(-30)
  args.clim.append(30)
# plot roughly 20 contours
args.clim.append(int((args.clim[1] - args.clim[0]) / 20))

#--------------------------------------------------------
# load files if necessary
data_netcdf4 = netCDF4.Dataset(args.data_dir + args.fileMOC)
if args.lprint:
  print(data_netcdf4)
  data_netcdf4.close()
  sys.exit("finished displaying data in file, exiting...")
sigma = data_netcdf4.variables["sigma"][:]
latV  = data_netcdf4.variables["latV"][:]
rmoc  = data_netcdf4.variables[args.MOC_var][0, :, :]
data_netcdf4.close()

#--------------------------------------------------------
# define some options and markers

if args.region.lower() == "so":
  print(" ")
  print("plotting with SO focus...")
  
  ylim = []
  ylim.append(-75)
  ylim.append(-30)
  slim = []
  slim.append(34)
  slim.append(37.2)

  # define some masks to pull out some values
  # 1) NADW as max of selected domain
  latV_mask  = (latV > ylim[0]) & (latV < ylim[1])
  sigma_mask = (sigma > slim[0]) & (sigma < slim[1])
  rmoc_select = rmoc[sigma_mask, :][:, latV_mask]
  NADW_strength = amax(rmoc_select)
  NADW_index = where(rmoc_select == NADW_strength)

  NADW_lat = latV[latV_mask][NADW_index[1]][0]
  NADW_sig = sigma[sigma_mask][NADW_index[0]][0]

  # 2) SAMW as min of lighter water above sigma = 35.5
#  sigma_mask = (sigma < 35.5) & (sigma > slim[0])
#  rmoc_select = rmoc[sigma_mask, :][:, latV_mask]
#  SAMW_strength = amin(rmoc_select)
#  SAMW_index = where(rmoc_select == SAMW_strength)
#  SAMW_lat = latV[latV_mask][SAMW_index[1]][0]
#  SAMW_sig = sigma[sigma_mask][SAMW_index[0]][0]

  # 3) AABW as min of denser water below sigma = 35.5
  sigma_mask = (sigma > 35.5)
  rmoc_select = rmoc[sigma_mask, :][:, latV_mask]
  AABW_strength = amin(rmoc_select)
  AABW_index = where(rmoc_select == AABW_strength)
  AABW_lat = latV[latV_mask][AABW_index[1]][0]
  AABW_sig = sigma[sigma_mask][AABW_index[0]][0]
  
elif args.region.lower() == "atlantic":
  print(" ")
  print("plotting with Atlantic focus...")
  
  ylim = []
  ylim.append(-30)
  ylim.append(80)
  slim = []
  slim.append(35)
  slim.append(37.5)
  
  latV_mask  = (latV > ylim[0]) & (latV < ylim[1])
  sigma_mask = (sigma > slim[0]) & (sigma < slim[1])
  rmoc_select = rmoc[sigma_mask, :][:, latV_mask]
  NADW_strength = amax(rmoc_select)
  NADW_index = where(rmoc_select == NADW_strength)

  NADW_lat = latV[latV_mask][NADW_index[1]][0]
  NADW_sig = sigma[sigma_mask][NADW_index[0]][0]
  
else:
  print(" ")
  print("plotting the global MOC...")
  
  ylim = []
  ylim.append(amin(latV))
  ylim.append(amax(latV))
  slim = []
  slim.append(amin(sigma))
  slim.append(amax(sigma))

#--------------------------------------------------------
# Main plotting commands

fig = plt.figure(figsize=(6, 6))
ax = plt.axes()
mesh = ax.contourf(latV, sigma, rmoc, arange(args.clim[0], args.clim[1] + 1, args.clim[2]), 
                   extend = "both", cmap = "Spectral_r")
ax.set_xlim(ylim[0], ylim[1])
ax.set_ylim(slim[0], slim[1])
ax.set_xlabel(r"Lat (${}^\circ$)")
ax.set_ylabel(r"$\sigma_{2000}$ ($\mathrm{kg}\ \mathrm{m}^{-3}$)")
ax.invert_yaxis() # flip arrays upside down (lighter density class = higher up)
ax.set_title(args.MOC_var)

if args.region.lower() in ["so", "atlantic"]:
  ax.plot(NADW_lat, NADW_sig, "k+")
  ax.text(NADW_lat - 20, NADW_sig, "NADW = %.1f Sv" % NADW_strength)
if args.region.lower() == "so":
#  ax.plot(SAMW_lat, SAMW_sig, "k.")
#  ax.text(SAMW_lat - 20, SAMW_sig, "SAMW = %.1f Sv" % SAMW_strength)
  ax.plot(AABW_lat, AABW_sig, "kv")
  ax.text(AABW_lat - 7, AABW_sig - 0.2, "AABW = %.1f Sv" % AABW_strength)

cb = plt.colorbar(mesh)
cb.ax.set_title(r"Sv")
  
#--------------------------------------------------------
# saving commands

if args.file_out is None:
  args.file_out = args.fileMOC.replace(".nc", "") + "_sigmaMOC.png"

fig.savefig(args.file_out, dpi = 300, bbox_inches = "tight")
plt.close(fig)

print("generated %s , exiting..." % args.file_out)
