#!/usr/bin/env python3

# JM: 01 Sep 2018
# generates the MOC in sigma co-ordinates; see flags for options below
#
# input a data_dir and it tries to grab *all* the files matching a string within
# the folder, so may want to do some pre-processing first
#
# by default generates a "mocsig.nc", there is a specify output name

from pyCDFTOOLS.cdfmocsig import *
import numpy as np
import glob, netCDF4, copy

#--------------------------------------------------------
# define the argument parser
import argparse

parser = argparse.ArgumentParser(description = 
    """
    Generate the time-averaged MOC in sigma-coordinates, see flags for options below.
    """)

# fixed arguments
parser.add_argument("data_dir", type = str, 
                    help = "specify data directory")
parser.add_argument("fileV",    type = str, 
                    help = "specify string in V-grid file that matches the data (e.g. \"*MOC_V*\", note the quote marks)")
                    
# optional arguments
parser.add_argument("--v_var",  type = str, 
                    help = "specify v variable name (default = voce)", default = "voce")
parser.add_argument("--t_var",  type = str, 
                    help = "specify temperature variable name (default = toce)", default = "toce")
parser.add_argument("--s_var",  type = str, 
                    help = "specify salinity variable name (default = soce)", default = "soce")

parser.add_argument("--lprint", 
                    help = "print out the variables available in fileV", action = "store_true")
parser.add_argument("--lg_vvl",   
                    help = "using time-varying metric e3v (assume it's in fileV)", action = "store_true")
parser.add_argument("--lbas",
                    help = "use basin decomposition (assumes there is a new_maskglo.nc in data_dir)", action = "store_true")
parser.add_argument("--eivv_var",   type = str,
                    help = "add the eddy induced velocity (give the variable name here)")
parser.add_argument("--lisodep",   
                    help = "diagonse mean depth of isopycnal", action = "store_true")
parser.add_argument("--lntr",   
                    help = "use netural density for binning", action = "store_true")
                    
parser.add_argument("--file_out",  type = str, 
                    help = "specify output name (default = mocsig_tave.nc)", default = "mocsig_tave.nc")                    

# collect arguments
args = parser.parse_args()

#--------------------------------------------------------

# initialisation

kwargs = {"lprint" : args.lprint,
          "lg_vvl" : args.lg_vvl,
          "lbas"   : args.lbas,
          "lisodep": args.lisodep,
          "lntr"   : args.lntr}

if args.lg_vvl:
  print(" ")
  print("using time-varying metric e3v")
  print(" ")
  
if args.lbas:
  print(" ")
  print("using basin decomposition")
  print(" ")
  
if args.lisodep:
  print(" ")
  print("diagnosing mean isopycnal depth")
  print(" ")

if args.eivv_var is not None:
  kwargs["leiv"] = True
  kwargs["eivv_var"] = args.eivv_var
  print(" ")
  print("including eddy induced velocity...")
  print(" ")

# grab file names  
file_list = []
for file in glob.glob(args.data_dir + args.fileV + ".nc"):
    file_list.append(file)

# sort it according to the timestamps
file_list.sort()

if args.lprint:
  for file in file_list:
    print(file)

# define bins
if args.lntr:
  sys.exit("neutral density option not implemented yet, stopping...")
else:
  # generate bins from presets (available for 0, 1000, 2000m depth), or define it in a dictionary as
  bins = {"pref"   : 2000,
          "nbins"  : 316,
          "sigmin" : 30.0,
          "sigstp" : 0.025}

  #bins = sigma_bins(2000)
  
  print(" ")
  print("using potential density referenced to %g m" % bins["pref"])
  print("%g bins from sigmin = %g to sigmax = %g" 
    % (bins["nbins"], bins["sigmin"], bins["sigmin"] + (bins["nbins"] - 1) * bins["sigstp"])
       )
  print(" ")

#--------------------------------------------------------
# cycle through the file lists and compile the moc
print("%g files found, cycling through them..." % len(file_list))

for i in range(len(file_list)):
    fileV = file_list[i].replace(args.data_dir, "") # strip out the data_dir
    fileT = fileV.replace("_V", "_T") # replace V with T
    print(" ")
    print("working in file = %g / %g" % (i + 1, len(file_list)))
    if i == 0:
        sigma, depi_temp, latV, dmoc_temp, opt_dic = \
          cdfmocsig_tave(args.data_dir, fileV, args.v_var, fileT, args.t_var, args.s_var, bins, **kwargs)
        dmoc = dmoc_temp / len(file_list)
        depi = depi_temp / len(file_list)
    else:
        _, depi_temp, _, dmoc_temp, _ = \
          cdfmocsig_tave(args.data_dir, fileV, args.v_var, fileT, args.t_var, args.s_var, bins, **kwargs)
        dmoc += dmoc_temp / len(file_list)
        depi += depi_temp / len(file_list)
        
print(" ")
print("returning final time-averaged field")

#--------------------------------------------------------
# write the file

latV_mesh = np.zeros(depi[0, :, :].shape)
for j in range(latV_mesh.shape[0]):
    latV_mesh[j, :] = latV

# open a new netCDF file for writing.
ncfile = netCDF4.Dataset(args.data_dir + args.file_out, "w", format = "NETCDF4") 
ncfile.title = "diagnosed MOC in density co-ordinates"

# create the dimensions.
ncfile.createDimension("sigma", sigma.shape[0])
ncfile.createDimension("y", latV.shape[0])
ncfile.createDimension("time", len(np.asarray([0.0])))

# first argument is name of variable, 
# second is datatype,
# third is a tuple with the names of dimensions.

lat_netcdf = ncfile.createVariable("latV", np.dtype("float32").char, "y")
lat_netcdf[:] = latV
lat_netcdf.units = "deg"
lat_netcdf.long_name = "y"

sigma_netcdf = ncfile.createVariable("sigma", np.dtype("float32").char, "sigma")
sigma_netcdf[:] = sigma
sigma_netcdf.units = "kg m-3"
sigma_netcdf.long_name = "sigma"

# global moc
moc_glob_tave_netcdf = ncfile.createVariable("rmoc_glob_tave", np.dtype("float32").char, 
                                            ("time", "sigma", "y"), fill_value = 1e20)
moc_glob_tave_netcdf[:] = dmoc[0, :, :]
moc_glob_tave_netcdf.units = "Sv"
moc_glob_tave_netcdf.long_name = "global RMOC in density co-ordinates"

# add the other MOCs and variables here after figuring out how to split them

if opt_dic["lbas"]:
    moc_atl_tave_netcdf = ncfile.createVariable("rmoc_atl_tave", np.dtype("float32").char, 
                                               ("time", "sigma", "y"), fill_value = 1e20)
    moc_atl_tave_netcdf[:] = dmoc[1, :, :]
    moc_atl_tave_netcdf.units = "Sv"
    moc_atl_tave_netcdf.long_name = "Atlantic RMOC in density co-ordinates"
    
    moc_indpac_tave_netcdf = ncfile.createVariable("rmoc_indpac_tave", np.dtype("float32").char, 
                                                  ("time", "sigma", "y"), fill_value = 1e20)
    moc_indpac_tave_netcdf[:] = dmoc[2, :, :]
    moc_indpac_tave_netcdf.units = "Sv"
    moc_indpac_tave_netcdf.long_name = "Indo-Pacific RMOC in density co-ordinates"
    
    moc_ind_tave_netcdf = ncfile.createVariable("rmoc_ind_tave", np.dtype("float32").char, 
                                               ("time", "sigma", "y"), fill_value = 1e20)
    moc_ind_tave_netcdf[:] = dmoc[3, :, :]
    moc_ind_tave_netcdf.units = "Sv"
    moc_ind_tave_netcdf.long_name = "Indian RMOC in density co-ordinates"
    
    moc_pac_tave_netcdf = ncfile.createVariable("rmoc_pac_tave", np.dtype("float32").char, 
                                               ("time", "sigma", "y"), fill_value = 1e20)
    moc_pac_tave_netcdf[:] = dmoc[4, :, :]
    moc_pac_tave_netcdf.units = "Sv"
    moc_pac_tave_netcdf.long_name = "Pacific RMOC in density co-ordinates"

if opt_dic["lisodep"]:
    latV_mesh_netcdf = ncfile.createVariable("latV_mesh", np.dtype("float32").char, ("sigma", "y"))
    latV_mesh_netcdf[:] = latV_mesh
    latV_mesh_netcdf.units = "deg"
    latV_mesh_netcdf.long_name = "latV in a mesh"

    isodep_glob_tave_netcdf = ncfile.createVariable("isodep_glob_tave", np.dtype("float32").char, 
                                                   ("time", "sigma", "y"), fill_value = 1e20)
    isodep_glob_tave_netcdf[:] = depi[0, :, :]
    isodep_glob_tave_netcdf.units = "m"
    isodep_glob_tave_netcdf.long_name = "global mean isopycnal depth"
    
    if opt_dic["lbas"]:
        isodep_atl_tave_netcdf = ncfile.createVariable("isodep_atl_tave", np.dtype("float32").char, 
                                                      ("time", "sigma", "y"), fill_value = 1e20)
        isodep_atl_tave_netcdf[:] = depi[1, :, :]
        isodep_atl_tave_netcdf.units = "m"
        isodep_atl_tave_netcdf.long_name = "Atlantic mean isopycnal depth"
        
        isodep_indpac_tave_netcdf = ncfile.createVariable("isodep_indpac_tave", np.dtype("float32").char, 
                                                         ("time", "sigma", "y"), fill_value = 1e20)
        isodep_indpac_tave_netcdf[:] = depi[2, :, :]
        isodep_indpac_tave_netcdf.units = "m"
        isodep_indpac_tave_netcdf.long_name = "Indo-Pacific mean isopycnal depth"
        
        isodep_ind_tave_netcdf = ncfile.createVariable("isodep_ind_tave", np.dtype("float32").char, 
                                                    ("time", "sigma", "y"), fill_value = 1e20)
        isodep_ind_tave_netcdf[:] = depi[3, :, :]
        isodep_ind_tave_netcdf.units = "m"
        isodep_ind_tave_netcdf.long_name = "Indian mean isopycnal depth"
        
        isodep_pac_tave_netcdf = ncfile.createVariable("isodep_pac_tave", np.dtype("float32").char, 
                                                      ("time", "sigma", "y"), fill_value = 1e20)
        isodep_pac_tave_netcdf[:] = depi[4, :, :]
        isodep_pac_tave_netcdf.units = "m"
        isodep_pac_tave_netcdf.long_name = "global mean isopycnal depth"

# close the file.
ncfile.close()

print("*** SUCCESS writing to %s! ***" % args.file_out)

