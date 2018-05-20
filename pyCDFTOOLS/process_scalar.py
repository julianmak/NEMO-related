#!/usr/bin/env python3
#
# subfunction for processing scalar data
#  !!======================================================================
#  !!                     ***  PROGRAM  process_scalar  ***
#  !!=====================================================================
#  !!  ** Purpose : (i)  read "scalar" data and dump it out into a text file
#  !!               (ii) read the text file accordingly
#  !!
#  !!  ** Method  : python native, brute force coded at the moment
#  !!
#  !! History : new  : 05/2018  : J. Mak : Original code
#  !!----------------------------------------------------------------------
#

import netCDF4
import os, glob

def scalar_to_txt(data_dir):
  """
  Given a data_dir, process all files with the "*scalar*.nc*" match
  and spit out some of its contents 
  
  (because for whatever reason I cannot get XIOS to just spit out a number, and
  currently it spits out a 2d field...)

  variable names and order hardcoded in here at the moment
  
  Inputs:
    data_dir = string for data directory
    
  Returns:
    nothing except the text file
  """

  # grab the relevant filenames
  file_list = []
  for file in glob.glob(data_dir + "*scalar*.nc"):
    file_list.append(file) 

  # sort it according to the timestamps
  file_list.sort()

  # cycle through the files
  for i in range(len(file_list)):
    # grab output time in years
    # assumes file format is $EXP_$PERIOD_$START_$END_scalar.nc
    # so pulls out the $START and $END and keeps only the first four entries
    # string here for use in output
    start_time = file_list[i].replace(data_dir, "").split("_")[2][0:4]
    end_time   = file_list[i].replace(data_dir, "").split("_")[3][0:4]

    data = netCDF4.Dataset(file_list[i])
    t = data.variables["time_centered"][:]
      
    # pull out the data written in 2d field for some reason and write it out
    txt_file = open(data_dir + file_list[i].replace(data_dir, "").replace(".nc", ".txt"), "w")
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
        % (file_list[i].replace(data_dir, ""), kt, len(t))
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
