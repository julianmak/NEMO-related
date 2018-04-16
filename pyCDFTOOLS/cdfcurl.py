#!/usr/bin/env python3
#
# subfunction for generating the curl of a variable
# adapted from CDFTOOLS/cdfcurl:
#  !!======================================================================
#  !!                     ***  PROGRAM  cdfcurl  ***
#  !!=====================================================================
#  !!  ** Purpose : Compute the curl on F-points for given gridU gridV 
#  !!               files and variables
#  !!
#  !!  ** Method  : Use the same algorithm than NEMO
#  !!
#  !! History : 2.1  : 05/2005  : J.M. Molines : Original code
#  !!         : 2.1  : 06/2007  : P. Mathiot   : for use with forcing fields
#  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
#  !!         : 4.0  : 03/2017  : J.M. Molines  
#  !!----------------------------------------------------------------------
#

from numpy import zeros, pi, sin
from netCDF4 import Dataset
import copy
from sys import exit

def cdfcurl(data_dir, u_file, u_var, v_file, v_var, **kwargs):
  """
  Compute the curl on something F-points from U and V variables
  
  EDIT: unlike cdfcurl.f90 I did not put an option to work out 
        things like wind-stress curl for A-grid data
        (e.g. CORE-II forcing, where both components are on the same point)

  Needs associated mesh_mask.nc file in the same data folder
  
  Inputs:
    data_dir = string for data directory
    u_file   = string for file with u in
    u_var    = string for u variable name
    v_file   = string for file with v in
    v_var    = string for v variable name
  
  Optional arguments (as of 10 Apr 2018):
    kt       = number for using a specified time entry
    kz       = vertical level to use
    lprint   = True   for printing out variable names in netcdf file
    lperio   = True   for imposing E-W periodicity
    loverf   = True   for having curl / f 
    
  Returns:
    lonF (glamf), latF (gphif), curl for plotting, opt_dic for record
  """
  # some defaults for optional arguments
  opt_dic = {"kt"     : 0, 
             "kz"     : 0,
             "lprint" : False,
             "lperio" : False,
             "loverf" : False}

  # overwrite the options by cycling through the input dictionary
  for key in kwargs:
    opt_dic[key] = kwargs[key]

  # open some files and pull variables out
  cn_mask = Dataset(data_dir + "mesh_mask.nc")
  npiglo  = cn_mask.dimensions["x"].size
  npjglo  = cn_mask.dimensions["y"].size
  glamf   = cn_mask.variables["glamf"][0, :, :]
  gphif   = cn_mask.variables["gphif"][0, :, :]
  e1u     = cn_mask.variables["e1u"][0, :, :]
  e1f     = cn_mask.variables["e1f"][0, :, :]
  e2v     = cn_mask.variables["e2v"][0, :, :]
  e2f     = cn_mask.variables["e2f"][0, :, :]
  fmask   = cn_mask.variables["fmask"][0, opt_dic["kz"], :, :]
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
  
  # check periodicity
  if (glamf[0, 0] - glamf[0, npiglo-2] == 0): # python indexing
    print("E-W periodicity detected!")
    if not opt_dic["lperio"]:
      print("--> forcing lperio = True")
      opt_dic["lperio"] = True
  
  curl = zeros((npjglo, npiglo))
  
  for ji in range(npiglo-1):
    for jj in range(npjglo-1):
        curl[jj, ji] = (  e2v[jj  , ji+1] * zv[jj  , ji+1] 
                        - e2v[jj  , ji  ] * zv[jj  , ji  ]
                        - e1u[jj+1, ji  ] * zu[jj+1, ji  ]
                        + e1u[jj  , ji  ] * zu[jj  , ji  ]
                       ) * fmask[jj, ji] / (e1f[jj, ji] * e2f[jj, ji])
  
  if opt_dic["lperio"]:
    curl[:, npiglo-1] = curl[:, 1] # python indexing
    
  if opt_dic["loverf"]:
    dl_omega = (2 * pi / 86400.0)
    dl_ff = 2 * dl_omega * sin(gphif * pi / 180.0)
    # avoid divide by zero
    curl[(dl_ff == 0)] = 0.0
    dl_ff[(dl_ff == 0)] = 1.0
    curl = curl / dl_ff

  return (glamf, gphif, curl, opt_dic)
