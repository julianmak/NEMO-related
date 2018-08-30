#!/usr/bin/env python3

# JM: 30 Aug 2018
# process the *scalar.nc files (because I couldn't get XIOS to spit out just
# a number for whatever reason...)

import netCDF4
import glob, sys

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

# collect arguments
args = parser.parse_args()

#--------------------------------------------------------

# Main commands

# grab the relevant filenames
file_list = []
for file in glob.glob(args.data_dir + "*scalar*.nc"):
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

  data = netCDF4.Dataset(file_list[i])
  t = data.variables["time_centered"][:]
  
  if args.lquery:
    for name, variable in data.variables.items():
      for attrname in variable.ncattrs():
        if attrname == "standard_name":
          print("{} -- {}".format(name, getattr(variable, attrname)))
    
    data.close()
    
    print(" ")
    
    sys.exit("finished query, exiting scalar_to_txt...")
  
  else:
  
    # ?? could do something like the query loop above to pull out all keys
    #    and dump accordingly; leaving it as manual for now
    
    # pull out the data written in 2d field for some reason and write it out
    
    txt_file = open(args.data_dir + file_list[i].replace(args.data_dir, "").replace(".nc", ".txt"), "w")
    txt_file.write("%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n"
                  % ("time",
                     "scvoltot", "scsshtot", "scsshste", "scsshtst", "sctemtot", 
                     "scsaltot",
                     "bgtemper", "bgsaline", "bgheatco", "bgsaltco", "bgvolssh", 
                     "bgvole3t", "bgfrcvol", "bgfrctem", "bgfrcsal"
                    )
                  )
    for kt in range(len(t)):
      print("processing %s at index %i / %i..." 
        % (file_list[i].replace(args.data_dir, ""), kt, len(t))
           )
      # global mean/totals
      time     = (t[kt] - t[0]) / (3600 * 24 * 365) + int(start_time)
      scvoltot = data.variables["scvoltot"][kt, 0, 0] # sea water volume,                     m^3
      scsshtot = data.variables["scsshtot"][kt, 0, 0] # mean ssh,                             m
      scsshste = data.variables["scsshste"][kt, 0, 0] # mean ssh steric,                      m
      scsshtst = data.variables["scsshtst"][kt, 0, 0] # mean ssh thermosteric,                m
      sctemtot = data.variables["sctemtot"][kt, 0, 0] # mean temperature,                     C
      scsaltot = data.variables["scsaltot"][kt, 0, 0] # mean salinity,                        psu
      #global drifts wrt time step 1
      bgtemper = data.variables["bgtemper"][kt, 0, 0] # mean temperature,                     C
      bgsaline = data.variables["bgsaline"][kt, 0, 0] # mean saline,                          psu    
      bgheatco = data.variables["bgheatco"][kt, 0, 0] # mean heat content,              1e+20 J
      bgsaltco = data.variables["bgsaltco"][kt, 0, 0] # mean salt content,              1e-3  km^3
      bgvolssh = data.variables["bgvolssh"][kt, 0, 0] # mean ssh volume,                      km^3
      bgvole3t = data.variables["bgvole3t"][kt, 0, 0] # mean volume variation,                km^3
      bgfrcvol = data.variables["bgfrcvol"][kt, 0, 0] # mean volume from forcing,             km^3
      bgfrctem = data.variables["bgfrctem"][kt, 0, 0] # mean heat content from forcing, 1e+20 J
      bgfrcsal = data.variables["bgfrcsal"][kt, 0, 0] # mean salt content from forcing, 1e-3  km^3
      txt_file.write("%.2f %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e\n"
                    % (time,
                       scvoltot, scsshtot, scsshste, scsshtst, sctemtot, 
                       scsaltot,
                       bgtemper, bgsaline, bgheatco, bgsaltco, bgvolssh, 
                       bgvole3t, bgfrcvol, bgfrctem, bgfrcsal
                      )
                    )
          
    txt_file.close()
    data.close()

  print("finished processing, exiting scalar_to_txt...")
