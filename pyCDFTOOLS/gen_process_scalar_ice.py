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
for file in glob.glob(args.data_dir + "*[1][2][3][1]_SBC_scalar.nc"): # distinguish this with *SBC_scalar.nc designed for ORCA
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
    
    sys.exit("finished query, exiting scalar_to_txt_ice...")
  
  else:
  
    # ?? could do something like the query loop above to pull out all keys
    #    and dump accordingly; leaving it as manual for now
    
    # pull out the data written in 2d field for some reason and write it out
    
    txt_file = open(args.data_dir + file_list[i].replace(args.data_dir, "").replace(".nc", ".txt"), "w")
    txt_file.write("%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n"
                  % ("time",
                     "ibgvol_tot", "sbgvol_tot", "ibgarea_tot", "ibgsalt_tot", "ibgheat_tot", "sbgheat_tot",
                     "ibgvolume", "ibgsaltco", "ibgheatco", "ibgheatfx", "ibgfrcvoltop", 
                     "ibgfrcvolbot", "ibgfrctemtop", "ibgfrctembot", "ibgfrcsal", "ibgfrchfxtop", "ibgfrchfxbot"
                    )
                  )
    for kt in range(len(t)):
      print("processing %s at index %i / %i..." 
        % (file_list[i].replace(args.data_dir, ""), kt, len(t))
           )
      # global mean/totals
      time        = (t[kt] - t[0]) / (3600 * 24 * 365) + int(start_time)
      ibgvol_tot  = data.variables["ibgvol_tot"][kt, 0, 0]  # global mean ice volume,            km^3
      sbgvol_tot  = data.variables["sbgvol_tot"][kt, 0, 0]  # global mean snow volume,           km^3
      ibgarea_tot = data.variables["ibgarea_tot"][kt, 0, 0] # global mean ice area,              km^2
      ibgsalt_tot = data.variables["ibgsalt_tot"][kt, 0, 0] # global mean ice salt content  1e-3*km^3
      ibgheat_tot = data.variables["ibgheat_tot"][kt, 0, 0] # global mean ice heat content  1e20* J
      sbgheat_tot = data.variables["sbgheat_tot"][kt, 0, 0] # global mean snow heat content 1e20* J
      #global drifts wrt time step 1
      ibgvolume = data.variables["ibgvolume"][kt, 0, 0] # drift in ice/snow volume               km3
      ibgsaltco = data.variables["ibgsaltco"][kt, 0, 0] # drift in ice salt content          pss*km3
      ibgheatco = data.variables["ibgheatco"][kt, 0, 0] # drift in ice/snow heat content    1e20* J
      ibgheatfx = data.variables["ibgheatfx"][kt, 0, 0] # drift in ice/snow heat flux          W/ m^2
      
      ibgfrcvoltop = data.variables["ibgfrcvoltop"][kt, 0, 0] # global mean ice/snow forcing at interface ice/snow-atm          km^3
      ibgfrcvolbot = data.variables["ibgfrcvolbot"][kt, 0, 0] # global mean ice/snow forcing at interface ice/snow-ocean        km^3
      ibgfrctemtop = data.variables["ibgfrctemtop"][kt, 0, 0] # global mean heat on top of ice/snw/ocean-atm               1e20* J
      ibgfrctembot = data.variables["ibgfrctembot"][kt, 0, 0] # global mean heat below ice (on top of ocean)               1e20* J
      ibgfrcsal    = data.variables["ibgfrcsal"   ][kt, 0, 0] # global mean ice/snow forcing (salt equivalent ocean volume) pss*km^3
      ibgfrchfxtop = data.variables["ibgfrchfxtop"][kt, 0, 0] # global mean heat flux on top of ice/snw/ocean-atm             W/ m^2
      ibgfrchfxbot = data.variables["ibgfrchfxbot"][kt, 0, 0] # global mean heat flux below ice (on top of ocean)             W/ m^2
      txt_file.write("%.2f %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e\n"
                    % (time,
                       ibgvol_tot, sbgvol_tot, ibgarea_tot, ibgsalt_tot, ibgheat_tot, sbgheat_tot,
                       ibgvolume, ibgsaltco, ibgheatco, ibgheatfx, ibgfrcvoltop, 
                       ibgfrcvolbot, ibgfrctemtop, ibgfrctembot, ibgfrcsal,
                       ibgfrchfxtop, ibgfrchfxbot
                      )
                    )
          
    txt_file.close()
    data.close()

  print("finished processing, exiting scalar_to_txt_ice...")
