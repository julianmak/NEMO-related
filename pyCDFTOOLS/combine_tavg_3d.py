#!/usr/bin/env python3

# JM: 16 Nov 2018
# generic script to combine multiple netCDF4 files into one average
# do this by opening a grabbed list of files and then averaging it by the length
# of the file list
# only does it per variable, option to output

import numpy as np
import glob, netCDF4, copy, sys

#--------------------------------------------------------
# define the argument parser
import argparse

parser = argparse.ArgumentParser(description = 
    """
    Generic script to combine multiple netCDF4 files into one file as an average.
    See flags for options.
    """)

# fixed arguments
parser.add_argument("data_dir", type = str, 
                    help = "specify data directory")
parser.add_argument("file",     type = str, 
                    help = "specify string in file that matches the data (e.g. \"*MOC_V*\", note the quote marks)")
parser.add_argument("var",      type = str, 
                    help = "specify variable name in file (e.g. votemper)")
parser.add_argument("depth_var",    type = str, 
                    help = "specify depth variable name in file (e.g. usually deptht/w/u/v)")
parser.add_argument("file_out", type = str, 
                    help = "specify output name (e.g. kgm_tave.nc)")

# collect arguments
args = parser.parse_args()

#--------------------------------------------------------

# initialisation

# grab file names  
filenames = args.data_dir + args.file + ".nc"
file_list = []
for file in glob.glob(filenames):
    file_list.append(file)
    
if len(file_list) == 0:
  sys.exit("no files grabbed in %s, are you in the right folder?" % filenames)

# sort it according to the timestamps
file_list.sort()

#--------------------------------------------------------
# cycle through the file lists and compile the moc
print("%g files found, cycling through them..." % len(file_list))

for i in range(len(file_list)):
  file = file_list[i].replace(args.data_dir, "") # strip out the data_dir
  print(" ")
  print("working in file = %g / %g" % (i + 1, len(file_list)))
  data = netCDF4.Dataset(file_list[i])
  try:
    var = data.variables[args.var]
  except KeyError:
    print(data)
    data.close()
    sys.exit("key error found, look up the list printed above")
    
  # work on the assumption that there is a time-dimension
  if i == 0:
    long_name = copy.deepcopy(var.long_name)
    units     = copy.deepcopy(var.units)
    nav_lat   = data.variables["nav_lat"][:, :]
    nav_lon   = data.variables["nav_lon"][:, :]
    depth     = data.variables[args.depth_var][:]
    array  = var[0, :, :, :] / len(file_list)
  else:
    array += var[0, :, :, :] / len(file_list)

  data.close()
        
print(" ")
print("returning final time-averaged field")

#--------------------------------------------------------
# write the file

ncfile = netCDF4.Dataset(args.data_dir + args.file_out, "w", format = "NETCDF4") 
ncfile.title = "combined variable %s" % args.var

# create the dimensions.

ncfile.createDimension("x", nav_lat.shape[1])
ncfile.createDimension("y", nav_lat.shape[0])
ncfile.createDimension("z", depth.shape[0])
ncfile.createDimension("time", len(np.asarray([0.0])))

# first argument is name of variable, 
# second is datatype,
# third is a tuple with the names of dimensions.

lon_netcdf = ncfile.createVariable("nav_lon", np.dtype("float32").char, ("y", "x"))
lon_netcdf[:] = nav_lon
lon_netcdf.units = "deg"
lon_netcdf.long_name = "x"

lat_netcdf = ncfile.createVariable("nav_lat", np.dtype("float32").char, ("y", "x"))
lat_netcdf[:] = nav_lat
lat_netcdf.units = "deg"
lat_netcdf.long_name = "y"

dep_netcdf = ncfile.createVariable("depth", np.dtype("float32").char, ("z"))
dep_netcdf[:] = depth
dep_netcdf.units = "m"
dep_netcdf.long_name = "z"

var_netcdf = ncfile.createVariable(args.var, np.dtype("float32").char, ("time", "z", "y", "x"))
var_netcdf[:] = array
var_netcdf.units = units
var_netcdf.long_name = long_name

print("*** SUCCESS writing to %s! ***" % args.file_out)

