#!/usr/bin/env python3

# JM: 26 Oct 2018
# process the transport from *grid_U.nc files and write them out

import glob, sys

from numpy import zeros, int8, where, sum, mean
from netCDF4 import Dataset
from pyCDFTOOLS.cdfpsi import *

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
parser.add_argument("--trans_var", type = str,
                    help = "variable name of transport diagnostic",
                    default = "uocetr_eff")
parser.add_argument("--uo_var", type = str,
                    help = "variable name of velocity",
                    default = "uo")                    
parser.add_argument("--lpng",
                    help = "output png files", action = "store_true")

# collect arguments
args = parser.parse_args()

#--------------------------------------------------------

# Main commands

# grab the relevant filenames
file_list = []
for file in glob.glob(args.data_dir + args.keys + "grid_U.nc"):
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
    
    sys.exit("finished query, exiting gen_ubt_info...")
  
  else:
  
    # ?? could do something like the query loop above to pull out all keys
    #    and dump accordingly; leaving it as manual for now
    
    # pull out the data written in 2d field for some reason and write it out

    txt_filename = args.data_dir + file_list[i].replace(args.data_dir, "").replace("grid_U.nc", "utrans.txt")
    
    txt_file = open(txt_filename, "w")
    txt_file.write( "%s %s %s\n" % ("time", "utot", "ubot") )
    
    # load the mask files outside
    data = Dataset(args.data_dir + "mesh_mask.nc")
    jpiglo = data.dimensions["x"].size
    jpjglo = data.dimensions["y"].size
    latU = data.variables["gphiu"][0, :, :]
    lonU = data.variables["glamu"][0, :, :]
    umask = data.variables["umask"][0, :, :, :]
    tmask = data.variables["tmask"][0, :, :, :]
    data.close()
    
    # find the bottom index (don't use mbathy here because that's for T-grid)
    k_bot = zeros((jpjglo, jpiglo))
    for jj in range(jpjglo):
        for ji in range(jpiglo):
            k_bot[jj, ji] = where((umask[:, jj, ji] == 0))[0][0]
    k_bot = int8(k_bot)
    
    for kt in range(len(t)):
      print("processing %s at index %i / %i..." 
        % (file_list[i].replace(args.data_dir, ""), kt, len(t))
           )
      # global mean/totals
      time     = (t[kt] - t[0]) / (3600 * 24 * 365) + int(start_time)
      
      data = Dataset(file_list[i])
      u_trans = data.variables[args.trans_var][0, :, :, :] # this should have all the length factors in already
      e3u = data.variables["e3u"][0, :, :, :]
      data.close()

      real_depth = sum(e3u * umask, axis = 0)

      u_tot = sum(u_trans, axis = 0) # vertically integrate
      u_bot = zeros((jpjglo, jpiglo))
      for jj in range(jpjglo):
          for ji in range(jpiglo):
              if k_bot[jj, ji] > 0:
                  u_bot[jj, ji] = u_trans[k_bot[jj, ji] - 1, jj, ji] / e3u[k_bot[jj, ji] - 1, jj, ji] * real_depth[jj, ji]
                
      # calculate the barotropic streamfunction as a way to define the ACC
      file_U = file_list[i].replace(args.data_dir, "")
      file_V = file_U.replace("_U", "_V")
      kwargs = {"lprint" : False,
                "ll_v"   : False,
                "lg_vvl" : True}

      lonT, latT, psi, opt_dic = cdfpsi(args.data_dir, file_U, args.uo_var, file_V, "vo", **kwargs)
      psi *= tmask[0, :, :] / 1e6 # mask and normalise
      ACC_mask = (latT < - 40) & (psi > 1.0) # define the ACC grid points
                
      # calculate the transports over the ACC defined by the barotropic streamfunction
      T_tot = mean(sum(u_tot * ACC_mask, axis = 0)) / 1e6
      T_bot = mean(sum(u_bot * ACC_mask, axis = 0)) / 1e6
      
      txt_file.write( "%.2f %.8f %.8f\n" % (time, T_tot, T_bot) )
          
    txt_file.close()

  print("finished processing, exiting gen_ubt_info...")
