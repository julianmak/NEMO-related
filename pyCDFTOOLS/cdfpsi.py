#!/usr/bin/env python3
#
# subfunction for generating barotropic streamfunction from NEMO data
# adapted from CDFTOOLS/cdfpsi:
#  !!======================================================================
#  !!                     ***  PROGRAM  cdfpsi  ***
#  !!=====================================================================
#  !!  ** Purpose : Compute Barotropic Stream Function
#  !!
#  !!  ** Method  : Compute the 2D fields dtrpu, dtrpv as the integral on 
#  !!               the vertical of u, v on their respective points.
#  !!               Then integrate from south to north : ==> dpsiu
#  !!               Then integrate from West to East   : ==> dpsiv
#  !!                  (should be almost the same (if no error ))
#  !!               Default (appropriate for global model): output dpsiu;
#  !!               normalizes the values setting psi (jpi,jpj) = 0
#  !!               If option "V" is given as last argument, output dpsiv,
#  !!               normalizes values setting psi(jpi,1) = 0.
#  !!               This is appropriate for North Atlantic
#  !!
#  !! History : 2.1  : 05/2005  : J.M. Molines : Original code
#  !!           3.0  : 05/2011  : J.M. Molines : Doctor norm + Lic.
#  !!         : 4.0  : 03/2017  : J.M. Molines  
#  !!=====================================================================
#
# warning: this is stripped down and does not yet include full functionality

import numpy as np
from netCDF4 import Dataset
import copy
from sys import exit

def cdfpsi(data_dir, u_file, u_var, v_file, v_var, **kwargs):
  """
  Compute Barotropic Stream Function

  Needs associated mesh_mask.nc file in the same data folder
  
  Inputs:
    data_dir = string for data directory
    u_file   = string for file with u in
    u_var    = string for u variable name
    v_file   = string for file with v in
    v_var    = string for v variable name
  
  Optional arguments (as of 21 Mar 2018):
    kt       = number for using a specified time entry
    lprint   = True   for printing out variable names in netcdf file
    ll_v     = True   for using V instead of the default U
    lg_vvl   = True   for using s-coord (time-varying metric)
    i[ij]ref = number for using a different normalising reference (need both to be specified)
    
  Returns:
    lonF (glamf), latF (gphif), dpsi for plotting, opt_dic for record
  """
  # some defaults for optional arguments
  opt_dic = {"kt"     : 0, 
             "ll_v"   : False, 
             "lprint" : False, 
             "lg_vvl" : False,
             "iiref"  : -1,
             "ijref"  : -1}

  # overwrite the options by cycling through the input dictionary
  for key in kwargs:
    opt_dic[key] = kwargs[key]

  # open some files and pull variables out
  if opt_dic["ll_v"]:
    cf_vfil = Dataset(data_dir + v_file)
    if opt_dic["lprint"]:
      print(cf_vfil)
    npiglo  = cf_vfil.dimensions["x"].size
    npjglo  = cf_vfil.dimensions["y"].size
    npk     = cf_vfil.dimensions["depthv"].size
    npt     = cf_vfil.dimensions["time_counter"].size
    zv      = cf_vfil.variables[v_var][opt_dic["kt"], :, :, :]
    if opt_dic["lg_vvl"]:
      e3v     = cf_vfil.varaibles["e3v"][opt_dic["kt"], :, :, :]
    cf_vfil.close()

    cn_mask = Dataset(data_dir + "mesh_mask.nc")
    glamf   = cn_mask.variables["glamf"][0, :, :]
    gphif   = cn_mask.variables["gphif"][0, :, :]
    e1v     = cn_mask.variables["e1v"][0, :, :]
    if not opt_dic["lg_vvl"]:
      e3v     = cn_mask.variables["e3v_0"][opt_dic["kt"], :, :, :]
    cn_mask.close()
  else:
    cf_ufil = Dataset(data_dir + u_file)
    if opt_dic["lprint"]:
      print(cf_ufil)
    npiglo  = cf_ufil.dimensions["x"].size
    npjglo  = cf_ufil.dimensions["y"].size
    npk     = cf_ufil.dimensions["depthu"].size
    npt     = cf_ufil.dimensions["time_counter"].size
    zu      = cf_ufil.variables[u_var][opt_dic["kt"], :, :, :]
    if opt_dic["lg_vvl"]:
      e3u     = cf_ufil.variables["e3u"][opt_dic["kt"], :, :, :]
    cf_ufil.close()

    cn_mask = Dataset(data_dir + "/mesh_mask.nc")
    glamf   = cn_mask.variables["glamf"][0, :, :]
    gphif   = cn_mask.variables["gphif"][0, :, :]
    e2u     = cn_mask.variables["e2u"][0, :, :]
    if not opt_dic["lg_vvl"]:
      e3u     = cn_mask.variables["e3u_0"][opt_dic["kt"], :, :, :]
    cn_mask.close()

  if (opt_dic["iiref"] == -1 or opt_dic["ijref"] == -1):
    opt_dic["iiref"] = npiglo
    opt_dic["ijref"] = npjglo

  # cumulative sum over depth
  if opt_dic["ll_v"]:
    dtrpv = np.zeros(e1v.shape)
    for jk in range(npk):
        dtrpv += zv[jk, :, :] * e1v[:, :] * e3v[jk, :, :]

    # do zonal integration
    dpsiv = np.zeros(e1v.shape)
    for ji in range(npiglo-2, 0, -1):
        dpsiv[:, ji] = dpsiv[:, ji+1] - dtrpv[:, ji]

    dpsi = copy.deepcopy(dpsiv)
  else:
    dtrpu = np.zeros(e2u.shape)
    for jk in range(npk):
        dtrpu += zu[jk, :, :] * e2u[:, :] * e3u[jk, :, :]

    # do meridional integration
    dpsiu = np.zeros(e2u.shape)
    for jj in range(1, npjglo):
        dpsiu[jj, :] = dpsiu[jj-1, :] - dtrpu[jj, :]

    dpsi = copy.deepcopy(dpsiu)

  # normalise
  dpsi = dpsi - dpsi[opt_dic["ijref"]-1, opt_dic["iiref"]-1] # python indexing

  return (glamf, gphif, dpsi, opt_dic)
