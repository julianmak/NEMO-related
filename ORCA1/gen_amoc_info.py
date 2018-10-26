#!/usr/bin/env python3

# JM: 13 Oct 2018
# process the AMOC from the *grid_V.nc files and write them to a text file

import glob, sys

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from numpy import zeros, argmax, unravel_index, amax, where, arange
from netCDF4 import Dataset

from pyCDFTOOLS.cdfmoc_atl import *

# style settings

plt.rcParams["font.family"] = "DejaVu Serif"
plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams["mathtext.rm"] = "serif"
plt.rcParams["image.cmap"] = "RdBu_r" # "*_r" is reverse of standard colour

#--------------------------------------------------------
# define the argument parser
import argparse

parser = argparse.ArgumentParser(description = "Process the scalar netcdf files and output as text")

# fixed arguments
parser.add_argument("data_dir", type = str, 
                    help = "specify data directory")
                    
# optional arguments
parser.add_argument("--lquery", 
                    help = "print out the variables available", action = "store_true")
parser.add_argument("--keys", type = str,
                    help = "grab a specified set of matching strs (enter string in quotes, default grab everything)",
                    default = "*")
parser.add_argument("--v_var", type = str,
                    help = "variable name of v velocity",
                    default = "vo")
parser.add_argument("--lpng",
                    help = "output png files", action = "store_true")

# collect arguments
args = parser.parse_args()

#--------------------------------------------------------

# plotting subroutine

def plot_amoc(NADW_info, zW, latV, dmoc, filename):
  fig = plt.figure(figsize = (8, 3))
  ax = plt.axes()
  mesh = ax.contourf(latV, zW, dmoc, arange(-25, 26, 2), cmap = "RdBu_r", extend = "both")
  ax.set_xlim(-30, 90)
  ax.plot(NADW_info[1], NADW_info[2], "k+")
  ax.text(NADW_info[1] - 10, NADW_info[2] - 500, "NADW = %.1f Sv" % NADW_info[0])
  ax.set_xlabel(r"Lat (${}^\circ$)")
  ax.set_ylabel(r"z ($\mathrm{m}$)")
  ax.set_title("Atlantic (no Med sea)")
  cb = plt.colorbar(mesh)
  cb.ax.set_title(r"$\mathrm{Sv}$")
  
  fig.savefig(filename, dpi = 150, bbox_inches = "tight")
  plt.close(fig)

  print("generated %s , exiting..." % filename)

#--------------------------------------------------------

# Main commands

# grab the relevant filenames
file_list = []
for file in glob.glob(args.data_dir + args.keys + "grid_V.nc"):
  file_list.append(file) 
  
if not file_list:
  print("no files grabbed, are you in the right directory?")
  print("no files grabbed, are you in the right directory?")
  print("no files grabbed, are you in the right directory?")

# sort it according to the timestamps
file_list.sort()

# cycle through the files
for i in range(len(file_list)):
  # grab output time in years
  # assumes file format is $EXP_$PERIOD_$START_$END_scalar.nc
  # so pulls out the $START and $END and keeps only the first four entries
  # string here for use in output
  start_time = file_list[i].replace(args.data_dir, "").split("_")[2][0:4]
  end_time   = file_list[i].replace(args.data_dir, "").split("_")[3][0:4]

  data = Dataset(file_list[i])
  t = data.variables["time_centered"][:]
  
  if args.lquery:
    for name, variable in data.variables.items():
      for attrname in variable.ncattrs():
        if attrname == "standard_name":
          print("{} -- {}".format(name, getattr(variable, attrname)))
    
    data.close()
    
    print(" ")
    
    sys.exit("finished query, exiting gen_amoc_info...")
  
  else:
  
    # ?? could do something like the query loop above to pull out all keys
    #    and dump accordingly; leaving it as manual for now
    
    # pull out the data written in 2d field for some reason and write it out

    txt_filename = args.data_dir + file_list[i].replace(args.data_dir, "").replace("grid_V.nc", "amoc.txt")
    png_filename = args.data_dir + file_list[i].replace(args.data_dir, "").replace("grid_V.nc", "amoc.png")
    
    txt_file = open(txt_filename, "w")
    txt_file.write( "%s %s %s %s\n" % ("time", "NADW_str", "NADW_lat", "NADW_dep") )
    for kt in range(len(t)):
      print("processing %s at index %i / %i..." 
        % (file_list[i].replace(args.data_dir, ""), kt, len(t))
           )
      # global mean/totals
      time     = (t[kt] - t[0]) / (3600 * 24 * 365) + int(start_time)
      kwargs = {"lprint" : False,
                "lg_vvl" : True,
                "leiv"   : True, "eivv_var" : "vo_eiv"}

      NADW_info, zW, latV, dmoc, _ = cdfmoc_atl(args.data_dir, file_list[i].replace(args.data_dir, ""), args.v_var, **kwargs)
      txt_file.write( "%.2f %.8f %.8f %.8f\n" % (time, NADW_info[0], NADW_info[1], NADW_info[2]) )
      
      if args.lpng:
        plot_amoc(NADW_info, zW, latV, dmoc[1, :, :], png_filename)
          
    txt_file.close()
    data.close()

  print("finished processing, exiting gen_amoc_info...")
