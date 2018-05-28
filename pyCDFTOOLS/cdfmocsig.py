#!/usr/bin/env python3
#
# subfunction for generating the MOC in density co-ordinates
# adapted from CDFTOOLS/cdfmocsig:
#  !!======================================================================
#  !!                     ***  PROGRAM  cdfmocsig  ***
#  !!=====================================================================
#  !!  ** Purpose : Compute the Meridional Overturning Cell (MOC)
#  !!               using density bins. 
#  !!
#  !!  ** Method  : The MOC is computed from the V velocity field, collected in density bins,
#  !!               (reference depth is given as the 3rd argument) and integrated
#  !!               throughout the density bins, then zonally averaged with
#  !!               eventual masking for oceanic basins.
#  !!               In the present version the masking corresponds to the global
#  !!               configuration. MOC for Global, Atlantic, Indo-Pacific, Indian,Pacific ocean
#  !!               Results are saved on mocsig.nc file with variables name respectively
#  !!               zomsfglo, zomsfatl, zomsfinp, zomsfind, zomsfpac.
#  !!               If no new_maskglo.nc file found, then the mask.nc file is used and
#  !!               only zomsfglo is computed.
#  !!
#  !! History : 2.1  : 11/2005  : A.M. Treguier : Original code from cdfmoc
#  !!                : 03/2010  : C. Dufour     : Choice of depth reference
#  !!                                             improvements 
#  !!           3.0  : 04/2011  : J.M. Molines  : Doctor norm + Lic.
#  !!           3.0  : 06/2013  : J.M. Molines  : add neutral density
#  !!           3.0  : 06/2013  : J.M. Molines  : add bin mean depth calculation and OpenMp directives
#  !!         : 4.0  : 03/2017  : J.M. Molines  
#  !!----------------------------------------------------------------------

from pyCDFTOOLS import eos

from numba import jit, int32 # need explicit type for use within jit
from numpy import zeros, argmax, unravel_index, linspace, maximum, minimum
from netCDF4 import Dataset

import numpy as np

def cdfmocsig(data_dir, v_file, v_var, t_file, t_var, s_var, bins, **kwargs):
  """
  Compute the MOC in density co-ordinates

  Needs associated mesh_mask.nc file in the same data folder
  
  Inputs:
    data_dir  = string for data directory
    v_file    = string for file with v in
    v_var     = string for v variable name
    t_file    = string for file with t-variables in
    t_var     = string for temperature variable name
    s_var     = string for s variable name
    bins      = specify density bins (as a dictionary)
  
  Optional arguments (as of 22 Apr 2018):
    lprint   = True   for printing out variable names in netcdf file
    lverb    = True   for printing out more information
    lg_vvl   = True   for using s-coord (time-varying metric)
    ldec     = True   decompose the MOC into some components
    leiv     = True   for adding the eddy induced velocity component
      eivv_var = string for EIV-v variable name
    lisodep  = True   (not yet implemented) output zonal averaged isopycnal depth
    lntr     = True   (not yet implemented) do binning with neutral density
    
  Returns:
    sigma (density bins), latV (rdumlat), dmoc for plotting, opt_dic for record
  """
  # some defaults for optional arguments
  opt_dic = {"kt"     : 0,
             "lprint" : False,
             "lverb"  : False,
             "lg_vvl" : False,
             "ldec"   : False,
             "leiv"   : False,
             "lisodep": False,
             "lntr"   : False}

  # overwrite the options by cycling through the input dictionary
  for key in kwargs:
    opt_dic[key] = kwargs[key]

  # open some files and pull variables out
  cf_vfil = Dataset(data_dir + v_file)
  if opt_dic["lprint"]:
    print(cf_vfil)
  npiglo  = cf_vfil.dimensions["x"].size
  npjglo  = cf_vfil.dimensions["y"].size
  npk     = cf_vfil.dimensions["depthv"].size
  zv      = cf_vfil.variables[v_var][opt_dic["kt"], :, :, :]
  if opt_dic["lg_vvl"]:
    e3v     = cf_vfil.variables["e3v"][opt_dic["kt"], :, :, :]
  if opt_dic["leiv"]:
    if opt_dic["lverb"]:
      print("using eddy induced velocity too")
    # load and add the contribution to zv
    zeiv    = cf_vfil.variables[opt_dic["eivv_var"]][opt_dic["kt"], :, :, :]
    zv += zeiv
  cf_vfil.close()
  
  cf_tfil = Dataset(data_dir + t_file)
  if opt_dic["lprint"]:
    print(cf_tfil)
  zt      = cf_tfil.variables[t_var][opt_dic["kt"], :, :, :]
  zs      = cf_tfil.variables[s_var][opt_dic["kt"], :, :, :]
  cf_tfil.close()
  
  cn_mask = Dataset(data_dir + "mesh_mask.nc")
  e1v     = cn_mask.variables["e1v"][0, :, :]
  gphiv   = cn_mask.variables["gphiv"][0, :, :]
  vmask   = cn_mask.variables["vmask"][0, :, :, :]
  tmask   = cn_mask.variables["tmask"][0, :, :, :]
  if not opt_dic["lg_vvl"]:
    e3v     = cn_mask.variables["e3v_0"][0, :, :, :]
  cn_mask.close()
  
  if opt_dic["ldec"]:
    print("NOT DONE THIS YET! -- 16 APR 2018")

#  --------------------------
#  0) Initialisations
#  --------------------------
  pref, nbins, sigmin, sigstp = bins["pref"], bins["nbins"], bins["sigmin"], bins["sigstp"]
  
  # make the density vector
  sigma = sigmin + (linspace(1, nbins, nbins) - 0.5) * sigstp
  if opt_dic["lverb"]:
    print("min density for binning = %g" % sigma[0])
    print("max density for binning = %g" % sigma[nbins - 1])
  
  # grab a latitude vector
  iloc = unravel_index(argmax(gphiv, axis = None), gphiv.shape)
  rdumlat = gphiv[:, iloc[1]]
  
  # work out areas
  zarea = e1v * e3v
  
#  --------------------------
#  1) Bin and integrate across density classes to form dmoc
#     [horizontal loop done using jit to speed it up (grindingly slow otherwise)]
#  --------------------------
  dmoc = zeros((nbins, npjglo))
  
  for jk in range(0, npk - 1):
    if opt_dic["lverb"]:
      if jk == 0:
        print("slow k loop, progress = ", end = "")
      if (jk % 10) == 0:
        print("%.2f %%..." % (jk / npk), end = "")
      
    #TODO (22 Apr): average T and S onto V point? (skip for now)
    if opt_dic["lntr"]:
      print("neutral density option not done yet! (22 Apr)")
      return (0.0, 0.0, 0.0, opt_dic)
    else:
      dens = eos.sigmai_dep(zt[jk, :, :], zs[jk, :, :], pref)
    
    zttmp = dens * tmask[jk, :, :]
    
    # find bin numbers
    ibin = (zttmp - sigmin) / sigstp
    ibin = minimum( maximum(ibin.astype(int), 1), nbins)
    
    # used in case broken line transects are required (where npjglo = 1)
    # NOT IMPLEMENTED YET
    ij1, ij2 = 1, npjglo - 1 # python indexing
    
    zv_weight = zv[jk, :, :] * zarea[jk, :, :]
    
    # sum up over the density classes
    # routine takes dmoc as an input and adds to it, so no += required
    dmoc = dmoc_loop(dmoc, zv_weight, ibin, nbins, ij1, ij2, npjglo, npiglo)

  # Integrate across the bins from high to low density
  dmoc[nbins - 1, :] /= 1.0e6
  for jbin in range(nbins - 1, 0, -1):
    dmoc[jbin - 1, :] = dmoc[jbin, :] + dmoc[jbin - 1, :] / 1.0e6

  return (sigma, rdumlat, dmoc, opt_dic)
  
#-------------------------------------------------------------------------------
@jit(nopython = True)
def dmoc_loop(dmoc, zv_weight, ibin, nbins, ij1, ij2, npjglo, npiglo):
  
  # Do the binning and convert k into sigma
  
  for jj in range(ij1, ij2):
    dmoc_tmp = zeros((nbins, npiglo))
    for ji in range(1, npiglo - 1):
      # cycle through the indices and bin according to density classes
      ib = ibin[jj, ji] - 1 # python indexing
      
      dmoc_tmp[ib, ji] -= zv_weight[jj, ji]
    
#      if opt_dic["lisodep"]:
#        print("diagnosing isopycnal depth not implemented yet! (22 Apr)")
#        return (0.0, 0.0, 0.0, opt_dic)
      
    # integrate 'zonally' (along i-coordingate)
    # add to dmoc
    for jbin in range(0, nbins):
      for ji in range(1, npiglo - 1):
        dmoc[jbin, jj] += dmoc_tmp[jbin, ji]
  
  return dmoc
  
#-------------------------------------------------------------------------------

def cdfmocsig_tave(data_dir, v_file, v_var, t_file, t_var, s_var, bins, **kwargs):
  """
  Compute the MOC in density co-ordinates as cdfmocsig, but grabs every time 
  entry (assumed to be equally spaced) and time-average it
    
  Returns:
    sigma (density bins), latV (rdumlat), dmoc for plotting, opt_dic for record
  """
  
  # load the relevant V file and see how many entries there are
  cf_vfil = Dataset(data_dir + v_file)
  nt = cf_vfil.dimensions["time_counter"].size
  cf_vfil.close()
  
  print("%g frames found, cycling through them..." % nt)
  print(" ")
  
  # cycle through every single time entry
  for kt in range(nt):
    kwargs["kt"] = kt
    if kt == 0:
      print("working at frame %g..." % (kt + 1), end = "")
      sigma, rdumlat, dmoc_temp, opt_dic = cdfmocsig(data_dir, v_file, v_var, t_file, t_var, s_var, bins, **kwargs)
      dmoc = dmoc_temp / nt
    else:
      print("%g..." % (kt + 1), end = "")
      _, _, dmoc_temp, _ = cdfmocsig(data_dir, v_file, v_var, t_file, t_var, s_var, bins, **kwargs)
      dmoc += dmoc_temp / nt
  
  print(" ")
  print("returning time-averaged field")
  
  return (sigma, rdumlat, dmoc, opt_dic)
  
#-------------------------------------------------------------------------------
  
def sigma_bins(pref):
  """
  Some pre-defined density bin options. Add to here as appropriate
  
  Inputs:
    pref = a reference depth
    
  Outputs:
    bins = a dictionary with some entries
  """
  if int(pref) == 0:
    bins = {"pref"   : int(pref),
            "nbins"  : 52,
            "sigmin" : 23.0,
            "sigstp" : 0.1}
  elif int(pref) == 1000:
    bins = {"pref"   : int(pref),
            "nbins"  : 88,
            "sigmin" : 24.0,
            "sigstp" : 0.1}
  elif int(pref) == 2000:
    bins = {"pref"   : int(pref),
            "nbins"  : 158,
            "sigmin" : 30.0,
            "sigstp" : 0.05}
  else:
    print("preset not defined, stopping...")
    return 0.0
    
  return bins
