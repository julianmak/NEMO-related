#!/usr/bin/env python3
#
# subfunction for generating the modulus of the velocity
# partly adapted from CDFTOOLS/cdfcurl and CDFTOOLS/cdfspeed:
#  !!======================================================================
#  !!                     ***  PROGRAM  cdfspeed  ***
#  !!=====================================================================
#  !!  ** Purpose : Compute current speed
#  !!
#  !!  ** Method  : Compute the 2D field absu from U and V velocity (on T point)
#  !!
#  !! History : new  : 04/2018  : J. Mak : "new" code though it turns out
#  !!                                      cdfspeed does something similar...
#  !!=====================================================================
#

import numpy as np
from netCDF4 import Dataset
import copy
from sys import exit

def cdfspeed(data_dir, u_file, u_var, v_file, v_var, **kwargs):
  """
  Compute (Absolute) Current Speed

  Needs associated mesh_mask.nc file in the same data folder
  
  Inputs:
    data_dir = string for data directory
    u_file   = string for file with u in
    u_var    = string for u variable name
    v_file   = string for file with v in
    v_var    = string for v variable name
  
  Optional arguments (as of 08 Apr 2018):
    kt       = number for using a specified time entry
    kz       = vertical level to use
    lprint   = True   for printing out variable names in netcdf file
    lperiod  = True   for imposing E-W periodicity
    
  Returns:
    lonT (glamt), latT (gphit), absu for plotting, opt_dic for record
  """
  # some defaults for optional arguments
  opt_dic = {"kt"     : 0, 
             "kz"     : 0,
             "lperio" : False, 
             "lprint" : False}

  # overwrite the options by cycling through the input dictionary
  for key in kwargs:
    opt_dic[key] = kwargs[key]

  # open some files and pull variables out
  cn_mask = Dataset(data_dir + "mesh_mask.nc")
  npiglo  = cn_mask.dimensions["x"].size
  npjglo  = cn_mask.dimensions["y"].size
  glamt   = cn_mask.variables["glamt"][0, :, :]
  gphit   = cn_mask.variables["gphit"][0, :, :]
  e1u     = cn_mask.variables["e1u"][0, :, :]  
  e2v     = cn_mask.variables["e2v"][0, :, :]
  cn_mask.close()

  cf_ufil = Dataset(data_dir + u_file)
  if opt_dic["lprint"]:
    print(cf_ufil)
  zu      = cf_ufil.variables[u_var][opt_dic["kt"], opt_dic["kz"], :, :]
  cf_ufil.close()  
  
  cf_vfil = Dataset(data_dir + v_file)
  if opt_dic["lprint"]:
    print(cf_vfil)
  zv      = cf_vfil.variables[v_var][opt_dic["kt"], opt_dic["kz"], :, :]
  cf_vfil.close()
  
  # average the velocities onto the T points
  absu = np.zeros((npjglo, npiglo))
  
  for ji in range(npiglo-1):
    for jj in range(npjglo-1):
        absu[jj, ji] = np.sqrt(  ( (e1u[jj, ji+1] * zu[jj, ji  ] 
                                  + e1u[jj, ji  ] * zu[jj, ji+1]) 
                                 / (e1u[jj, ji] + e1u[jj, ji+1]) ) ** 2
                              +  ( (e2v[jj+1, ji] * zv[jj  , ji] 
                                  + e2v[jj, ji  ] * zv[jj+1, ji]) 
                                 / (e2v[jj, ji] + e2v[jj+1, ji]) ) ** 2
                              )

  # check periodicity
  if (np.count_nonzero(glamt[0, 0] - glamt[0, npiglo-2]) == 0): # python indexing
    print("E-W periodicity detected!")
    if not opt_dic["lperio"]:
      print("--> forcing lperio = True")
      opt_dic["lperio"] = True
  
  if opt_dic["lperio"]:
    absu[:, npiglo-1] = absu[:, 1] # python indexing

  return (glamt, gphit, absu, opt_dic)
