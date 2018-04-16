#!/usr/bin/env python3
#
# subfunction for generating the MOC
# adapted from CDFTOOLS/cdfmoc (see also cdfmocsig):
#  !!======================================================================
#  !!                     ***  PROGRAM  cdfmoc  ***
#  !!=====================================================================
#  !!  ** Purpose : Compute the Meridional Overturning Cell (MOC)
#  !!
#  !!  ** Method  : The MOC is computed from the V velocity field, integrated
#  !!               from the bottom to the surface, then zonally averaged with
#  !!               eventual masking for oceanic basins.
#  !!               The program looks for the file "new_maskglo.nc". If it 
#  !!               does not exist, only the calculation over all the domain
#  !!               is performed (this is adequate for a basin configuration).
#  !!               In new_maskglo.nc the masking corresponds to the global
#  !!               configuration. MOC for Global, Atlantic, Indo-Pacific, 
#  !!               Indian, Pacific ocean, inp0=Global-Atlantic
#  !!               Results are saved on moc.nc file with variables name 
#  !!               respectively zomsfglo, zomsfatl, zomsfinp, zomsfind, zomsfpac, zomsinp0
#  !!
#  !! History : 2.1  : 07/2005  : J.M. Molines  : Original code
#  !!                : 04/2006  : A.M. Treguier : Adaptation to NATL4 case
#  !!                : 09/2007  : G. Smith      : MOC decomposition
#  !!                : 01/2008  : A. Lecointre  : MOC decomposition adaptation 
#  !!           3.0  : 03/2011  : J.M. Molines  : Merge all MOC prog, Doctor norm + Lic.
#  !!                : 10/2012  : M.A. Balmaseda: it adds basin INP0=GLOBAL-ATL, different from INP.
#  !!         : 4.0  : 03/2017  : J.M. Molines  
#  !!         
#  !!
#  !! References :  For MOC decomposition : Lee & Marotzke (1998), 
#  !!               Baehr, Hirschi, Beismann &  Marotzke (2004),
#  !!               Cabanes, Lee, & Fu (2007),  Koehl & Stammer (2007).
#  !!               See also the powerpoint presentation by Tony Lee at the third 
#  !!               CLIVAR-GSOP intercomparison  available at : 
#  !!    http://www.clivar.org/organization/gsop/synthesis/mit/talks/lee_MOC_comparison.ppt
#  !!               See :   AMOC Metrics guidelines availablke on:
#  !!    https://www.godae-oceanview.org/documents/q/action-edit/ref-264/parent-261/
#  !!----------------------------------------------------------------------

from numpy import zeros, argmax, unravel_index
from netCDF4 import Dataset

def cdfmoc(data_dir, v_file, v_var, **kwargs):
  """
  Compute (Absolute) Current Speed

  Needs associated mesh_mask.nc file in the same data folder
  
  Inputs:
    data_dir = string for data directory
    v_file   = string for file with v in
    v_var    = string for v variable name
  
  Optional arguments (as of 16 Apr 2018):
    lprint   = True   for printing out variable names in netcdf file
    lg_vvl   = True   for using s-coord (time-varying metric)
    ldec     = True   decompose the MOC into some components
    
  Returns:
    zW (gdepw_1d), latV (rdumlat), dmoc for plotting, opt_dic for record
  """
  # some defaults for optional arguments
  opt_dic = {"kt"     : 0,
             "lprint" : False,
             "lg_vvl" : False,
             "ldec"   : False}

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
    e3v     = cf_vfil.varaibles["e3v"][opt_dic["kt"], :, :, :]
  cf_vfil.close()
  
  cn_mask = Dataset(data_dir + "mesh_mask.nc")
  e1v     = cn_mask.variables["e1v"][0, :, :]
  gphiv   = cn_mask.variables["gphiv"][0, :, :]
  gdepw   =-cn_mask.variables["gdepw_1d"][0, :]
  vmask   = cn_mask.variables["vmask"][0, :, :, :]
  if not opt_dic["lg_vvl"]:
    e3v     = cn_mask.variables["e3v_0"][opt_dic["kt"], :, :, :]
  cn_mask.close()
  
  if opt_dic["ldec"]:
    print("NOT DONE THIS YET! -- 16 APR 2018")

#  --------------------------
#  0) Create a dummy latitude to output
#  --------------------------
  iloc = unravel_index(argmax(gphiv, axis = None), gphiv.shape)
  rdumlat = gphiv[:, iloc[1]]
  
#  --------------------------
#  1) Compute total MOC: dmoc
#  --------------------------
  dmoc = zeros((npk, npjglo)) # change this if having the multiple basin decomposition
  
  # integrate 'zonally' (along i-coordinate)
  for jk in range(npk):
    if (jk % 10) == 0:
      print("processing in the slow loop, jk = %i / %g" % (jk, npk))
    for jj in range(npjglo):
      for ji in range(npiglo):
        dmoc[jk, jj] -= e1v[jj, ji] * e3v[jk, jj, ji] * vmask[jk, jj, ji] * zv[jk, jj, ji]
        
  # integrate vertically from bottom to surface
  for jj in range(npjglo):
    for jk in range(npk-2, 0, -1): # python indexing
      dmoc[jk, jj] = dmoc[jk + 1, jj] + dmoc[jk, jj] / 1.0e6

  return (gdepw, rdumlat, dmoc, opt_dic)
