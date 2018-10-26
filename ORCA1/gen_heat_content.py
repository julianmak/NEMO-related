#!/usr/bin/env python3

# JM: 14 Oct 2018
# calculate the domain ocean heat content from the *grid_T.nc files and write 
# them to a text file

import glob, sys

from numpy import zeros, sum
from netCDF4 import Dataset

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
parser.add_argument("--t_var", type = str,
                    help = "variable name of t temperature",
                    default = "heatc")
parser.add_argument("--lpng",
                    help = "output png files", action = "store_true")

# collect arguments
args = parser.parse_args()

#--------------------------------------------------------

# subroutine in question

def calc_heatc(data_dir, t_file, t_var, **kwargs):
  """
  Compute the domain integrated heat content in different basins

  Needs associated mesh_mask.nc and new_maskglo.nc file in the same data folder
  
  Inputs:
    data_dir = string for data directory
    t_file   = string for file with t in
    t_var    = string for t variable name
  
  Optional arguments (as of 16 Apr 2018):
    lprint   = True   for printing out variable names in netcdf file
    lg_vvl   = True   for using s-coord (time-varying metric)
    lbas     = True   for computing heat content in different basins
    
  Returns:
    heatc (with basin entries), opt_dic for record
  """
  # some defaults for optional arguments
  opt_dic = {"kt"     : 0,
             "lprint" : False,
             "lg_vvl" : False,
             "lbas"   : False}

  # overwrite the options by cycling through the input dictionary
  for key in kwargs:
    opt_dic[key] = kwargs[key]

  # open some files and pull variables out
  cf_tfil = Dataset(data_dir + t_file)
  if opt_dic["lprint"]:
    print(cf_tfil)
  zheatc    = cf_tfil.variables[t_var][opt_dic["kt"], :, :]
  npiglo  = cf_tfil.dimensions["x"].size
  npjglo  = cf_tfil.dimensions["y"].size
  cf_tfil.close()
  
  cn_mask = Dataset(data_dir + "mesh_mask.nc")
  e1t     = cn_mask.variables["e1t"][0, :, :]
  e2t     = cn_mask.variables["e2t"][0, :, :]
  gphit   = cn_mask.variables["gphit"][0, :, :]
  tmask   = cn_mask.variables["tmask"][0, 0, :, :]
  cn_mask.close()
  
  # 0 : global ; 1/2 : Atlantic ; 3 : Indo-Pacif ; 4 : Indian ; 5 : Pacif
  cn_basin = Dataset(data_dir + "new_maskglo.nc")
  nbasins = 6
  ibmask = zeros((nbasins, npjglo, npiglo))
  ibmask[0, :, :] = tmask
  ibmask[1, :, :] = cn_basin.variables["atlmsk"][:, :]
  ibmask[2, :, :] = cn_basin.variables["atlmsk_nomed"][:, :]
  ibmask[3, :, :] = cn_basin.variables["indpacmsk"][:, :]
  ibmask[4, :, :] = cn_basin.variables["indmsk"][:, :]
  ibmask[5, :, :] = cn_basin.variables["pacmsk"][:, :]
  cn_basin.close()
  
  # initialise
  e1e2t = e1t * e2t * tmask
  
  heatc = zeros(nbasins)
  for ib in range(nbasins):
    heatc[ib] = sum(sum(zheatc * e1e2t * ibmask[ib, :, :]))

  return (heatc, opt_dic)

#--------------------------------------------------------

# Main commands

# grab the relevant filenames
file_list = []
for file in glob.glob(args.data_dir + args.keys + "grid_T.nc"):
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
    
    sys.exit("finished query, exiting gen_heat_content...")
  
  else:
  
    # ?? could do something like the query loop above to pull out all keys
    #    and dump accordingly; leaving it as manual for now
    
    # pull out the data written in 2d field for some reason and write it out

    txt_filename = args.data_dir + file_list[i].replace(args.data_dir, "").replace("grid_T.nc", "heatc.txt")
    
    txt_file = open(txt_filename, "w")
    txt_file.write( "%s %s %s %s %s %s %s\n" 
                   % ("time", "glob", "atl", "atl_nomed", "indpac", "ind", "pac") 
                  )
    for kt in range(len(t)):
      print("processing %s at index %i / %i..." 
        % (file_list[i].replace(args.data_dir, ""), kt, len(t))
           )
      # global mean/totals
      time     = (t[kt] - t[0]) / (3600 * 24 * 365) + int(start_time)
      kwargs = {"lprint" : False,
                "lg_vvl" : False,
                "lbas"   : True}

      heatc, _ = calc_heatc(args.data_dir, file_list[i].replace(args.data_dir, ""), args.t_var, **kwargs)
      txt_file.write( "%.2f %.6e %.6e %.6e %.6e %.6e %.6e\n" 
                     % (time, heatc[0], heatc[1], heatc[2], heatc[3], heatc[4], heatc[5]) )
          
    txt_file.close()
    data.close()

  print("finished processing, exiting gen_heat_content...")
