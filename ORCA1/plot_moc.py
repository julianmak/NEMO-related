#!/usr/bin/env python3

# JM: 30 Aug 2018
# plots the MOC averaged at fixed height, styling is default and this script
# is intended to be used for quick and dirty visualisations

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from numpy import linspace, amax, where
from pyCDFTOOLS.cdfmoc import *

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
parser.add_argument("fileV",    type = str, 
                    help = "specify data filename")

# optional arguments
parser.add_argument("--v_var",  type = str, 
                    help = "specify v variable name (default = vo)", default = "vo")
parser.add_argument("--lprint", 
                    help = "print out the variables available in fileV", action = "store_true")
parser.add_argument("--lg_vvl",   
                    help = "using time-varying metric e3v (assume it's in fileV)", action = "store_true")
parser.add_argument("--lbas",   
                    help = "use basin decomposition (assumes there is a new_maskglo.nc in data_dir)", action = "store_true")
parser.add_argument("--eivv_var",   type = str,
                    help = "add the eddy induced velocity (give the variable name here)")

# collect arguments
args = parser.parse_args()

#--------------------------------------------------------
# Main plotting commands

kwargs = {"lprint" : args.lprint,
          "lg_vvl" : args.lg_vvl,
          "lbas"   : args.lbas}

if args.lg_vvl:
  print(" ")
  print("using time-varying metric e3v")
  
if args.lbas:
  print(" ")
  print("using basin decomposition")

if args.eivv_var is not None:
  kwargs["leiv"] = True
  kwargs["eivv_var"] = args.eivv_var
  print("including eddy induced velocity...")

zW, latV, dmoc, opt_dic = cdfmoc(args.data_dir, args.fileV, args.v_var, **kwargs)

if args.lbas:

  fig = plt.figure(figsize=(10, 10))
  
  # find the max north of 30N below 500m
  latV_mask  = (latV > 30)
  zW_mask = (zW < -500)
  rmoc_select = dmoc[0, zW_mask, :][:, latV_mask]
  NADW_strength = amax(rmoc_select)
  NADW_index = where(rmoc_select == NADW_strength)
  NADW_lat = latV[latV_mask][NADW_index[1]][0]
  NADW_dep = zW[zW_mask][NADW_index[0]][0]

  ax1 = plt.subplot(3, 1, 1)
  mesh1 = ax1.contourf(latV, zW, dmoc[0, :, :], linspace(-20, 20, 21), cmap = "RdBu_r", extend = "both")
  ax1.plot(NADW_lat, NADW_dep, "k+")
  ax1.text(NADW_lat - 10, NADW_dep - 500, "NADW = %.1f Sv" % NADW_strength)
  ax1.set_ylabel(r"z ($\mathrm{m}$)")
  ax1.set_title("Global")
  cb = plt.colorbar(mesh1)
  cb.ax.set_title(r"$\mathrm{Sv}$")
  
  rmoc_select = dmoc[1, zW_mask, :][:, latV_mask]
  NADW_strength = amax(rmoc_select)
  NADW_index = where(rmoc_select == NADW_strength)
  NADW_lat = latV[latV_mask][NADW_index[1]][0]
  NADW_dep = zW[zW_mask][NADW_index[0]][0]
  NADW_lat = latV[latV_mask][NADW_index[1]][0]
  NADW_dep = zW[zW_mask][NADW_index[0]][0]

  ax2 = plt.subplot(3, 1, 2)
  mesh2 = ax2.contourf(latV, zW, dmoc[1, :, :], linspace(-20, 20, 21), cmap = "RdBu_r", extend = "both")
  ax2.plot(NADW_lat, NADW_dep, "k+")
  ax2.text(NADW_lat - 10, NADW_dep - 500, "NADW = %.1f Sv" % NADW_strength)
  ax2.set_ylabel(r"z ($\mathrm{m}$)")
  ax2.set_title("Atlantic")
  cb = plt.colorbar(mesh2)

  ax3 = plt.subplot(3, 1, 3)
  mesh3 = ax3.contourf(latV, zW, dmoc[0, :, :] - dmoc[1, :, :], linspace(-20, 20, 21), cmap = "RdBu_r", extend = "both")
  ax3.set_xlabel(r"Lat (${}^\circ$)")
  ax3.set_ylabel(r"z ($\mathrm{m}$)")
  ax3.set_title("Not Atlantic")
  cb = plt.colorbar(mesh3)
  
else:

  # find the max north of 30N below 500m
  latV_mask  = (latV > 30)
  zW_mask = (zW < -500)
  rmoc_select = dmoc[0, zW_mask, :][:, latV_mask]
  NADW_strength = amax(rmoc_select)
  NADW_index = where(rmoc_select == NADW_strength)
  NADW_lat = latV[latV_mask][NADW_index[1]][0]
  NADW_dep = zW[zW_mask][NADW_index[0]][0]


  # plot
  fig = plt.figure(figsize=(10, 3))
  ax = plt.axes()
  mesh = ax.contourf(latV, zW, dmoc[0, :, :], linspace(-20, 20, 21), cmap = "RdBu_r", extend = "both")
  ax.plot(NADW_lat, NADW_dep, "k+")
  ax.text(NADW_lat - 10, NADW_dep - 500, "NADW = %.1f Sv" % NADW_strength)
  ax.set_xlabel(r"Lat (${}^\circ$)")
  ax.set_ylabel(r"z ($\mathrm{m}$)")
  ax.set_title("Global")
  cb = plt.colorbar(mesh)
  cb.ax.set_title(r"$\mathrm{Sv}$")
  
#--------------------------------------------------------
# saving commands

save_filename = args.fileV.replace(".nc", "") + "_zMOC.png"

fig.savefig(save_filename, dpi = 300, bbox_inches = "tight")
plt.close(fig)

print("generated %s , exiting..." % save_filename)
