#!/usr/bin/env python3
#
# module for doing density relating things
# adapted from CDFTOOLS/eos:
#  !!======================================================================
#  !!                     ***  MODULE  eos  ***
#  !! All routines dealing with the Equation Of State of sea water
#  !!=====================================================================
#  !! History : 2.1  !  2004   : J.M. Molines : Original code ported
#  !!                                           from NEMO
#  !!           3.0    12/2010 : J.M. Molines : Doctor norm + Lic.
#  !!           4.0  : 03/2017 : J.M. Molines
#  !!           4.0  : 06/2018 : J.   Mak     : split out another TEOS-10 
#  !!                                           file + adapt from GSW toolbox
#  !!----------------------------------------------------------------------
#  !! Notes: splitting out a eos_pre_TEOS10 and eos_TEOS10
#  !!        tried to have same routine names to keep consistency, 
#  !!        ***should*** be able to just overwrite the existing "eos"
#  !!
#  !!        sigmantr requires a conversion, maybe one that doesn't will
#  !!        exist at some point
#  !!----------------------------------------------------------------------
#  !!   routines      : description
#  !!   sigmai_dep    : compute potential density referenced to some depth
#  !!   sigmantr      : compute neutral density
#  !!                   NOTE! this needs inputs as potential temperature
#  !!                                              practical salinity
#  !!----------------------------------------------------------------------

from numpy import zeros, ones, sqrt, abs, floor, int32, reshape, nan, asarray
from numpy import mean, nanmax, sum, isfinite, isnan, sin, pi
from scipy.io import loadmat
from copy import deepcopy
from scipy.interpolate import interp1d

################################################################################
# sigmai_dep (was gsw_rho)                    in-situ density (75-term equation)
#==========================================================================
# 
# USAGE:  
#  rho = gsw_rho(SA,CT,p)
#
# DESCRIPTION:
#  Calculates in-situ density from Absolute Salinity and Conservative 
#  Temperature, using the computationally-efficient expression for
#  specific volume in terms of SA, CT and p  (Roquet et al., 2015).
#
#  Note that potential density with respect to reference pressure, pr, is
#  obtained by calling this function with the pressure argument being pr
#  (i.e. "gsw_rho(SA,CT,pr)").
#
#  Note that this 75-term equation has been fitted in a restricted range of 
#  parameter space, and is most accurate inside the "oceanographic funnel" 
#  described in McDougall et al. (2003).  The GSW library function 
#  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if 
#  some of one's data lies outside this "funnel".  
#
# INPUT:
#  SA  =  Absolute Salinity                                        [ g/kg ]
#  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
#  p   =  sea pressure                                             [ dbar ]
#         ( i.e. absolute pressure - 10.1325 dbar )
#
#  SA & CT need to have the same dimensions.
#  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
#
# OUTPUT:
#  rho  =  in-situ density                                         [ kg/m^3 ]
#
# AUTHOR: 
#  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
#
# VERSION NUMBER: 3.05 (27th November, 2015)
#
# REFERENCES:
#  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
#   seawater - 2010: Calculation and use of thermodynamic properties.  
#   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
#   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
#    See appendix A.20 and appendix K of this TEOS-10 Manual. 
#
#  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003: 
#   Accurate and computationally efficient algorithms for potential 
#   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
#   pp. 730-741.
#
#  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
#   polynomial expressions for the density and specifc volume of seawater
#   using the TEOS-10 standard. Ocean Modelling, 90, pp. 29-43.
#
# The software is available from http://www.TEOS-10.org

#==========================================================================

def sigmai_dep(CT, SA, p):
  """
  Compute the in-situ (or potential) density from CONSERVATIVE Temperature, 
  ABSOLUTE Salinity and pressure fields using the TEOS-10 EOS
  
  p can be a field (computed with subroutine p_from_z say) to get in-situ 
  density or can be a number p = p_ref, then it is a potential density 
  referenced to p_ref
  
  Adapted from the MATLAB GSW toolbox (http://www.TEOS-10.org)
  
  Inputs:
    CT             = Conservative Temperature      t        deg celsius
    SA             = Absolute Salinity             s        g / kg
    p              = (reference) pressure          p        dbar
    
  Returns:
    sigmai_dep_out = (potential density            rho      kg / m^3
    
  """
  # ensures that SA is non-negative.
  SA = abs(SA)

  # deltaS = 24
  sfac = 0.0248826675584615                 # sfac   = 1/(40*(35.16504/35)).
  offset = 5.971840214030754e-1             # offset = deltaS*sfac.

  x2 = sfac * SA
  xs = sqrt(x2 + offset)
  ys = CT * 0.025
  z  = p * 1e-4

  v000 =  1.0769995862e-3
  v001 = -6.0799143809e-5
  v002 =  9.9856169219e-6
  v003 = -1.1309361437e-6
  v004 =  1.0531153080e-7
  v005 = -1.2647261286e-8
  v006 =  1.9613503930e-9
  v010 = -1.5649734675e-5
  v011 =  1.8505765429e-5
  v012 = -1.1736386731e-6
  v013 = -3.6527006553e-7
  v014 =  3.1454099902e-7
  v020 =  2.7762106484e-5
  v021 = -1.1716606853e-5
  v022 =  2.1305028740e-6
  v023 =  2.8695905159e-7
  v030 = -1.6521159259e-5
  v031 =  7.9279656173e-6
  v032 = -4.6132540037e-7
  v040 =  6.9111322702e-6
  v041 = -3.4102187482e-6
  v042 = -6.3352916514e-8
  v050 = -8.0539615540e-7
  v051 =  5.0736766814e-7
  v060 =  2.0543094268e-7
  v100 = -3.1038981976e-4
  v101 =  2.4262468747e-5
  v102 = -5.8484432984e-7
  v103 =  3.6310188515e-7
  v104 = -1.1147125423e-7
  v110 =  3.5009599764e-5
  v111 = -9.5677088156e-6
  v112 = -5.5699154557e-6
  v113 = -2.7295696237e-7
  v120 = -3.7435842344e-5
  v121 = -2.3678308361e-7
  v122 =  3.9137387080e-7
  v130 =  2.4141479483e-5
  v131 = -3.4558773655e-6
  v132 =  7.7618888092e-9
  v140 = -8.7595873154e-6
  v141 =  1.2956717783e-6
  v150 = -3.3052758900e-7
  v200 =  6.6928067038e-4
  v201 = -3.4792460974e-5
  v202 = -4.8122251597e-6
  v203 =  1.6746303780e-8
  v210 = -4.3592678561e-5
  v211 =  1.1100834765e-5
  v212 =  5.4620748834e-6
  v220 =  3.5907822760e-5
  v221 =  2.9283346295e-6
  v222 = -6.5731104067e-7
  v230 = -1.4353633048e-5
  v231 =  3.1655306078e-7
  v240 =  4.3703680598e-6
  v300 = -8.5047933937e-4
  v301 =  3.7470777305e-5
  v302 =  4.9263106998e-6
  v310 =  3.4532461828e-5
  v311 = -9.8447117844e-6
  v312 = -1.3544185627e-6
  v320 = -1.8698584187e-5
  v321 = -4.8826139200e-7
  v330 =  2.2863324556e-6
  v400 =  5.8086069943e-4
  v401 = -1.7322218612e-5
  v402 = -1.7811974727e-6
  v410 = -1.1959409788e-5
  v411 =  2.5909225260e-6
  v420 =  3.8595339244e-6
  v500 = -2.1092370507e-4
  v501 =  3.0927427253e-6
  v510 =  1.3864594581e-6
  v600 =  3.1932457305e-5

  v = v000 + ( 
        xs * (v100 + xs * (v200 + xs * (v300 + xs * (v400 + xs * (v500 
      + v600 * xs))))) + ys * (v010 + xs * (v110 + xs * (v210 + xs * (v310 
      + xs * (v410 + v510 * xs)))) + ys * (v020 + xs * (v120 + xs * (v220 
      + xs * (v320 + v420 * xs))) + ys * (v030 + xs * (v130 + xs * (v230 
      + v330 * xs)) + ys * (v040 + xs * (v140 + v240*xs) + ys * (v050 
      + v150 * xs + v060 * ys))))) + z * (v001 + xs * (v101 + xs * (v201 
      + xs * (v301 + xs * (v401 + v501 * xs)))) + ys * (v011 + xs * (v111
      + xs * (v211 + xs * (v311 + v411 * xs))) + ys * (v021 + xs * (v121 
      + xs * (v221 + v321 * xs)) + ys * (v031 + xs * (v131 + v231 * xs) 
      + ys * (v041 + v141 * xs + v051 * ys)))) + z * (v002 + xs * (v102 
      + xs * (v202 + xs * (v302 + v402 * xs))) + ys * (v012 + xs * (v112 
      + xs * (v212 + v312 * xs)) + ys * (v022 + xs * (v122 + v222 * xs) 
      + ys * (v032 + v132 * xs + v042 * ys))) + z * (v003 + xs * (v103 
      + v203 * xs) + ys * (v013 + v113 * xs + v023 * ys) + z * (v004 
      + v104 * xs + v014 * ys + z * (v005 + v006 * z)))))
              )

  sigmai_dep_out = (1 / v) - 1000
  
  return sigmai_dep_out

################################################################################
#    !! --------------------------------------------------------------------
#    !! ** Purpose :   Compute the  neutral volumic mass (kg/m3) from known
#    !!      potential temperature and salinity fields using an equation of
#    !!      state. 
#    !!
#    !! ** Method  :
#    !!       McDougall and Jackett (2005) equation of state.
#    !!              potential temperature         t        deg celsius
#    !!              salinity                      s        psu
#    !!              neutral density               rho      kg/m**3
#    !!       result is not masked at this stage.
#    !!       Check value: rho(20,35) = 1024.59416751197 kg/m**3  -1000.
#    !!       t = 20 deg celcius, s=35 psu
#    !!
#    !! ** References : McDougall and Jackett, The material derivative of neutral density
#    !!        Journal of Marine Research, 63, 159-185, 2005
#    !! --------------------------------------------------------------------
def sigmantr(CT, SA, z, lon, lat):
  """
  Compute the neutral volumic mass (kg/m3) from known CONSERVATIVE Temperature 
  and PRACTICAL salinity fields from the Jackett and McDougall (2005) EOS
  
  Done via a conversion of CT to pt and SA to ps routines. The CT and SA fields
  should be massaged into a 3d field in (depth, lat, lon) arrangement even if
  one of the dimensions is only of length 1 so that the cycling through the
  depth index is ok
  
  Inputs:
    CT           = Conservative Temperature      t       deg celsius
    SA           = Absolute Salinity             s         g / kg
    z            = z-location of data                      m
                   (negative z is ocean depth)
    lon          = longitude of data                     deg
    lat          = latitude  of data                     deg
    
  Returns:
    sigmantr_out = neutral density               rho      kg / m^3
    
  """
  zrau0 = 1000.0
  
  sigmantr_out = zeros(SA.shape)
  
  # cycle through the depth index
  for k in range(len(z)):
  
    pt = pt_from_CT(CT[k, :, :], SA[k, :, :])   # convert CT to pt
    p  = p_from_z(z[k], lat)                    # convert z  to p
    ps = sp_from_SA(SA[k, :, :], p, lon, lat)   # convert SA to ps
    
    zws = sqrt( abs(ps) )
    
    # Numerator
    # T-Polynome
    zr1 = ( ( (-4.3159255086706703e-4 * pt + 8.1157118782170051e-2) * pt 
             + 2.2280832068441331e-1 ) * pt + 1002.3063688892480e0
           )
    # S-T Polynome
    zr2 = ( (-1.7052298331414675e-7 * ps - 3.1710675488863952e-3 * pt 
          - 1.0304537539692924e-4) * ps
          )
    # Denominator
    # T-Polynome
    zr3 = ( ( ( (-2.3850178558212048e-9 * pt -1.6212552470310961e-7) * pt
          + 7.8717799560577725e-5 ) * pt + 4.3907692647825900e-5  ) * pt + 1.0e0
          )
    # S-T Polynome
    zr4 = ( ( ( -2.2744455733317707e-9 * pt * pt + 6.0399864718597388e-6) * pt 
          - 5.1268124398160734e-4 ) * ps
          )
    # S-T Polynome
    zr5 = ( -1.3409379420216683e-9 * pt * pt - 3.6138532339703262e-5) * ps * zws

    # Neutral density
    sigmantr_out[k, :, :] = ( zr1 + zr2 ) / ( zr3 + zr4 + zr5 ) - zrau0
  
  return sigmantr_out

################################################################################
# pt_from_CT (was gsw_pt_from_CT)                 potential temperature --> CT
#==========================================================================
#
# USAGE:  
#  pt = gsw_pt_from_CT(SA,CT)
#
# DESCRIPTION:
#  Calculates potential temperature (with a reference sea pressure of
#  zero dbar) from Conservative Temperature.  This function uses 1.5 
#  iterations through a modified Newton-Raphson (N-R) iterative solution 
#  proceedure, starting from a rational-function-based initial condition 
#  for both pt and dCT_dpt. 
#
# INPUT:
#  SA  =  Absolute Salinity                                        [ g/kg ]
#  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
#
#  SA & CT need to have the same dimensions.
#
# OUTPUT:
#  pt  =  potential temperature referenced to a sea pressure 
#         of zero dbar (ITS-90)                                   [ deg C ]
#
# AUTHOR: 
#  Trevor McDougall, David Jackett, Claire Roberts-Thomson and Paul Barker. 
#                                                      [ help@teos-10.org ]
#
# VERSION NUMBER: 3.05 (27th January 2015)
#
# REFERENCES:
#  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
#   seawater - 2010: Calculation and use of thermodynamic properties.  
#   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
#   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
#    See sections 3.1 and 3.3 of this TEOS-10 Manual.
#
#  McDougall, T.J., and S.J. Wotherspoon, 2014: A simple modification of 
#   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied 
#   Mathematics Letters, 29, 20-25.  
#
#  The software is available from http://www.TEOS-10.org

#==========================================================================

def pt_from_CT(CT, SA):
  """
  Compute POTENTIAL temperature from Conservative Temperature and 
  Absolute Salinity
  
  Adapted from the MATLAB GSW toolbox (http://www.TEOS-10.org)
  
  Inputs:
    CT           = Conservative Temperature      t        deg celsius
    SA           = Absolute Salinity             s        g / kg
    
  Returns:
    pt_out     = potential temperature         t        deg celsius
    
  """

  # ensures that SA is non-negative.
  SA = abs(SA)

  s1 = SA*0.995306702338459   # Note that 0.995306702338459 = (35./35.16504) 

  a0 = -1.446013646344788e-2
  a1 = -3.305308995852924e-3
  a2 =  1.062415929128982e-4
  a3 =  9.477566673794488e-1
  a4 =  2.166591947736613e-3
  a5 =  3.828842955039902e-3

  b0 =  1
  b1 =  6.506097115635800e-4
  b2 =  3.830289486850898e-3
  b3 =  1.247811760368034e-6

  a5CT = a5 * CT
  b3CT = b3 * CT
  CT_factor   = (a3 + a4 * s1 + a5CT)
  pt_num    = a0 + s1 * (a1 + a2 * s1) + CT * CT_factor
  pt_recden = 1.0 / ( b0 + b1 * s1 + CT * (b2 + b3CT) )
  pt        = pt_num * pt_recden      
                               # At this point the abs max error is 1.5e-2 deg C

  dpt_dCT = ( CT_factor + a5CT - (b2 + b3CT + b3CT) * pt ) * pt_recden
                            
  # start the 1.5 iterations through the modified Newton-Rapshon iterative 
  # method (McDougall and Wotherspoon, 2014). 

  CT_diff = CT_from_pt(pt, SA) - CT
  pt_old = pt
  pt = pt_old - CT_diff * dpt_dCT    # 1/2-way through the 1st modified N-R loop
                               # At this point the abs max error is 6.6e-5 deg C

  ptm = 0.5 * (pt + pt_old);

  # This routine calls gibbs_pt0_pt0(SA,pt0) to get the second derivative 
  # of the Gibbs function with respect to temperature at zero sea pressure.  
  
  dpt_dCT = -3991.86795711963 / ( (ptm + 273.15) * gibbs_pt0_pt0(ptm, SA) )
  pt = pt_old - CT_diff * dpt_dCT       # end of 1st full modified N-R iteration
                               # At this point the abs max error is 1.e-10 deg C

  CT_diff = CT_from_pt(pt, SA) - CT
  pt_old = pt
  pt_out = pt_old - CT_diff * dpt_dCT
                                     # 1.5 iterations of the modified N-R method

  # The abs max error of the result is 1.42e-14 deg C
  
  return pt_out

################################################################################
# CT_from_pt (was gsw_CT_from_pt)                   CT --> potential temperature
#==========================================================================
#
# USAGE:
#  CT = gsw_CT_from_pt(SA,pt)
#
# DESCRIPTION:
#  Calculates Conservative Temperature of seawater from potential 
#  temperature (whose reference sea pressure is zero dbar).
#
# INPUT:
#  SA  =  Absolute Salinity                                        [ g/kg ]
#  pt  =  potential temperature (ITS-90)                          [ deg C ]
#
#  SA & pt need to have the same dimensions.
#
# OUTPUT:
#  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
#
# AUTHOR: 
#  David Jackett, Trevor McDougall and Paul Barker     [ help@teos-10.org ]
#  
# VERSION NUMBER: 3.05 (27th January 2015)
#
# REFERENCES:
#  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
#   seawater - 2010: Calculation and use of thermodynamic properties.  
#   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
#   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
#    See section 3.3 of this TEOS-10 Manual. 
#
#  The software is available from http://www.TEOS-10.org
#
#==========================================================================

def CT_from_pt(pt, SA):
  """
  Compute CONSERVATIVE temperature from potential emperature and 
  Absolute Salinity
  
  Adapted from the MATLAB GSW toolbox (http://www.TEOS-10.org)
  
  Inputs:
    pt           = potential temperature         t        deg celsius
    SA           = Absolute Salinity             s          g / kg
    
  Returns:
    CT_out       = Conservative Temperature      t        deg celsius
    
  """

  # ensures that SA is non-negative.
  SA = abs(SA)

  sfac = 0.0248826675584615;                  # sfac = 1/(40.*(35.16504/35))

  x2 = sfac * SA
  x = sqrt(x2)
  y = pt * 0.025                           # normalize for F03 and F08.
  
  #--------------------------------------------------------------------------
  # The below polynomial for pot_enthalpy is the full expression for 
  # potential entahlpy in terms of SA and pt, obtained from the Gibbs 
  # function as below.  The above polynomial has simply collected like powers
  # of x and y so that it is computationally faster than calling the Gibbs 
  # function twice as is done in the commented code below.  When this code 
  # below is run, the results are identical to calculating pot_enthalpy as 
  # above, to machine precision.  
  #
  #  pr0 = zeros(size(SA));
  #  pot_enthalpy = gsw_gibbs(0,0,0,SA,pt,pr0) - ...
  #                       (273.15 + pt).*gsw_gibbs(0,1,0,SA,pt,pr0);
  #
  #-----------------This is the end of the alternative code------------------

  pot_enthalpy =  61.01362420681071 + (
      y * (168776.46138048015 +
      y * (-2735.2785605119625 + y * (2574.2164453821433 +
       y * (-1536.6644434977543 + y * (545.7340497931629 +
      (-50.91091728474331 - 18.30489878927802 * y) * y))))) +
      x2 * (268.5520265845071 + y * (-12019.028203559312 +
      y * (3734.858026725145 + y * (-2046.7671145057618 +
      y * (465.28655623826234 + (-0.6370820302376359 -
      10.650848542359153 * y) * y)))) +
      x * (937.2099110620707 + y * (588.1802812170108 +
      y * (248.39476522971285 + (-3.871557904936333 -
      2.6268019854268356 * y) * y)) + x * (-1687.914374187449 + 
      x * (246.9598888781377 +
      x * (123.59576582457964 - 48.5891069025409 * x)) +
      y * (936.3206544460336 +
      y * (-942.7827304544439 + y * (369.4389437509002 +
      (-33.83664947895248 - 9.987880382780322 * y) * y))))))
                                       )
      
  CT_out = pot_enthalpy / 3991.86795711963 # gsw_cp0 = 3991.86795711963

  return CT_out

################################################################################
# gibbs_pt0_pt0 (was gsw_gibbs_pt0_pt0)                    gibbs_tt at (SA,pt,0)
#==========================================================================
# This function calculates the second derivative of the specific Gibbs 
# function with respect to temperature at zero sea pressure.  The inputs 
# are Absolute Salinity and potential temperature with reference sea 
# pressure of zero dbar.  This library function is called by both 
# "gsw_pt_from_CT(SA,CT)" ,"gsw_pt0_from_t(SA,t,p)" and
# "gsw_pt_from_entropy(SA,entropy)".
#
# VERSION NUMBER: 3.05 (27th January 2015)
#
#==========================================================================

def gibbs_pt0_pt0(pt0, SA):
  """
  Compute second derivative of the specific Gibbs function with respect to 
  temperature at zero sea pressure, from Absolute Salinity and potential 
  temperature referenced to zero dbar
  
  (Used in the pt_from_CT function)
  
  Adapted from the MATLAB GSW toolbox (http://www.TEOS-10.org)
  
  Inputs:
    pt0           = potential temperature ref to z = 0       t       deg celsius
    SA            = Absolute Salinity                        s       g / kg
    
  Returns:
    gibbs_pt0_pt0 = second derivative of Gibbs function         (J/kg) / K^2
    
  """

  # ensures that SA is non-negative.
  SA = abs(SA)

  sfac = 0.0248826675584615                # sfac = 1/(40*(35.16504/35))

  x2 = sfac * SA
  x  = sqrt(x2)
  y  = pt0 * 0.025

  g03 = -24715.571866078 + (
      y * (4420.4472249096725 +
      y * (-1778.231237203896 +
      y * (1160.5182516851419 +
      y * (-569.531539542516 + y * 128.13429152494615))))
                            )

  g08 = x2 * (1760.062705994408 + x * (-86.1329351956084 +
      x * (-137.1145018408982 + y * (296.20061691375236 +
      y * (-205.67709290374563 + y * 49.9394019139016))) +
      y * (-60.136422517125 + y * 10.50720794170734)) +
      y * (-1351.605895580406 + y * (1097.1125373015109 +
      y * (-433.20648175062206 + y * 63.905091254154904 * y))))

  gibbs_pt0_pt0_out = (g03 + g08) * 0.000625

  return gibbs_pt0_pt0_out

################################################################################  
# sp_from_sa (was gsw_SP_from_SA)                                      SA --> ps
#==========================================================================
#
# USAGE:
#  [SP, in_ocean] = gsw_SP_from_SA(SA,p,long,lat)
#
# DESCRIPTION:
#  Calculates Practical Salinity from Absolute Salinity. 
#
# INPUT:
#  SA    =  Absolute Salinity                                      [ g/kg ]
#  p     =  sea pressure                                           [ dbar ]
#           ( i.e. absolute pressure - 10.1325 dbar )
#  long  =  longitude in decimal degrees                     [ 0 ... +360 ]
#                                                     or  [ -180 ... +180 ]
#  lat   =  latitude in decimal degrees north               [ -90 ... +90 ]
#
#  p, lat and long may have dimensions 1x1 or Mx1 or 1xN or MxN,
#  where SA is MxN.
#
# OUTPUT:
#  SP        =  Practical Salinity  (PSS-78)                   [ unitless ]
#  in_ocean  =  0, if long and lat are a long way from the ocean 
#            =  1, if long and lat are in the ocean
#  Note. This flag is only set when the observation is well and truly on
#    dry land; often the warning flag is not set until one is several 
#    hundred kilometres inland from the coast. 
#
# AUTHOR: 
#  David Jackett, Trevor McDougall and Paul Barker    [ help_gsw@csiro.au ]
#
# VERSION NUMBER: 3.05 (27th January 2015)
#
# REFERENCES:
#  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
#   seawater - 2010: Calculation and use of thermodynamic properties.  
#   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
#   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
#
#  McDougall, T.J., D.R. Jackett, F.J. Millero, R. Pawlowicz and 
#   P.M. Barker, 2012: A global algorithm for estimating Absolute Salinity.
#   Ocean Science, 8, 1123-1134.  
#   http://www.ocean-sci.net/8/1123/2012/os-8-1123-2012.pdf 
#
#  The software is available from http://www.TEOS-10.org
#
#==========================================================================

def sp_from_SA(SA, p, lon, lat):
  """
  Compute PRACTICAL salinity from Absolute Salinity, pressure, longitude and 
  latitude. Uses a look up table.
  
  Adapted from the MATLAB GSW toolbox (http://www.TEOS-10.org)
  
  Inputs:
    SA           = Absolute Salinity             s         g / kg
    p            = pressure                      p      dbar
    lon          = longitude of data                     deg
    lat          = latitude  of data                     deg
    
  Returns:
    sp_out       = practical salinity                    psu
    
  """                   
                        
  # turn things into 1d vectors here for SAAR_mod, then reshape it back
  orig_shape = SA.shape
  lat = lat.flatten()
  lon = lon.flatten()
  SA  = SA.flatten()
  
  # move the longitude
  lon[lon < 0] += 360 
  
  # lat and lon should have the same size as SA already
  p   *= ones(lat.shape) 
  SP   = zeros(lat.shape)
  SAAR = zeros(lat.shape)

  SAAR = SAAR_mod(p, lon, lat)

  sp_out = (35 / 35.16504) * SA / (1 + SAAR)

  SP_baltic = SP_from_SA_Baltic(SA, lon, lat)

  dum_cond = (~isnan(SP_baltic))
  if dum_cond.any():
    Ibaltic = asarray([i for i in range(len(dum_cond)) if dum_cond[i]])
    sp_out[Ibaltic] = SP_baltic[Ibaltic]
  
  sp_out = reshape(sp_out, orig_shape)
  
  return sp_out
  
################################################################################
# SP_from_SA_Baltic (was gsw_SP_from_SA_Baltic)
#==========================================================================
#
# USAGE:  
#  SP_baltic = gsw_SP_from_SA_Baltic(SA,long,lat)
#
# DESCRIPTION:
#  Calculates Practical Salinity for the Baltic Sea, from a value computed
#  analytically from Absolute Salinity.
#  Note. This programme will only produce Practical Salinty values for the
#    Baltic Sea.
#
# INPUT:
#  SA    =  Absolute Salinity in the Baltic Sea                 [ g kg^-1 ]
#  long  =  Longitude in decimal degress east                [ 0 ... +360 ]    
#  lat   =  Latitude in decimal degress north               [ -90 ... +90 ]  
#
# OUTPUT:
#  SP_baltic  =  Practical Salinity                            [ unitless ]
#
# AUTHOR: 
#  David Jackett, Trevor McDougall & Paul Barker       [ help@teos-10.org ]
#
# VERSION NUMBER: 3.05 (27th January 2015)
#
# REFERENCES:
#  Feistel, R., S. Weinreben, H. Wolf, S. Seitz, P. Spitzer, B. Adel, 
#   G. Nausch, B. Schneider and D. G. Wright, 2010c: Density and Absolute 
#   Salinity of the Baltic Sea 2006-2009.  Ocean Science, 6, 3-24.
#   http://www.ocean-sci.net/6/3/2010/os-6-3-2010.pdf 
#
#  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
#   seawater - 2010: Calculation and use of thermodynamic properties.  
#   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
#   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
#
#  McDougall, T.J., D.R. Jackett, F.J. Millero, R. Pawlowicz and 
#   P.M. Barker, 2012: A global algorithm for estimating Absolute Salinity.
#   Ocean Science, 8, 1123-1134.  
#   http://www.ocean-sci.net/8/1123/2012/os-8-1123-2012.pdf 
#
#  The software is available from http://www.TEOS-10.org
#
#==========================================================================

def SP_from_SA_Baltic(SA, lon, lat):
  """
  Calculates Practical Salinity for the Baltic Sea, from a value computed
  analytically from Absolute Salinity. Note this program will only produce 
  Practical Salinty values for the Baltic Sea.
  
  Adapted from the MATLAB GSW toolbox (http://www.TEOS-10.org)
  
  Inputs (as vectors!):
    SA           = absolute salinity                  g / kg
    lon          = longitude of data                     deg
    lat          = latitude  of data                     deg
    
  Returns (as vector):
    SP           = practical salinity                    psu
    
  """
  
  # Baltic sea co-ordinates
  xb1 = 12.6
  xb2 = 7
  xb3 = 26
  xb1a = 45
  xb3a = 26

  yb1 = 50 
  yb2 = 59 
  yb3 = 69

  SP_baltic = nan * ones(SA.shape)

  dum_cond1 = ((xb2 < lon) & (lon < xb1a) & (yb1 < lat) & (lat < yb3))
  if dum_cond1.any():
    inds_baltic = asarray([i for i in range(len(dum_cond1)) if dum_cond1[i]])
    f = interp1d(asarray([yb1, yb2, yb3]), asarray([xb1, xb2, xb3]))
    xx_left = f(lat[inds_baltic])
    f = interp1d(asarray([yb1, yb3]), asarray([xb1a, xb3a]))
    xx_right = f(lat[inds_baltic])
    dum_cond2 = ((xx_left <= lon[inds_baltic]) & (lon[inds_baltic] <= xx_right))
    if dum_cond2.any():
        inds_baltic1 = asarray([i for i in range(len(dum_cond2)) if dum_cond2[i]])
        SP_baltic[inds_baltic[inds_baltic1]] = (
            35 / (35.16504 - 0.087) * (SA[inds_baltic[inds_baltic1]] - 0.087)
                                                )
  
  return SP_baltic
  
################################################################################
# SAAR_mod (was gsw_SAAR)     ASalinity Anomaly Ratio (excluding the Baltic Sea)
#==========================================================================
#
# USAGE:  
#  [SAAR, in_ocean] = gsw_SAAR(p,long,lat)
#
# DESCRIPTION:
#  Calculates the Absolute Salinity Anomaly Ratio, SAAR, in the open ocean
#  by spatially interpolating the global reference data set of SAAR to the
#  location of the seawater sample.  
# 
#  This function uses version 3.0 of the SAAR look up table (15th May 2011). 
#
#  The Absolute Salinity Anomaly Ratio in the Baltic Sea is evaluated 
#  separately, since it is a function of Practical Salinity, not of space. 
#  The present function returns a SAAR of zero for data in the Baltic Sea. 
#  The correct way of calculating Absolute Salinity in the Baltic Sea is by 
#  calling gsw_SA_from_SP.  
#
# INPUT:
#  p     =  sea pressure                                           [ dbar ] 
#          ( i.e. absolute pressure - 10.1325 dbar )
#  long  =  Longitude in decimal degrees                     [ 0 ... +360 ]
#                                                      or [ -180 ... +180 ]
#  lat   =  Latitude in decimal degrees north               [ -90 ... +90 ]
#
#  p, long & lat need to be vectors and have the same dimensions.
#
# OUTPUT:
#  SAAR      =  Absolute Salinity Anomaly Ratio                [ unitless ]
#  in_ocean  =  0, if long and lat are a long way from the ocean 
#            =  1, if long and lat are in or near the ocean
#  Note. This flag is only set when the observation is well and truly on
#    dry land; often the warning flag is not set until one is several 
#    hundred kilometres inland from the coast. 
#
# AUTHOR: 
#  David Jackett                                       [ help@teos-10.org ]
#
# MODIFIED:
#  Paul Barker and Trevor McDougall 
#  Acknowledgment. Matlab programming assisance from Sunke Schmidtko.
#
# VERSION NUMBER: 3.06.10 (22nd March 2018)
#
# REFERENCES:
#  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
#   seawater - 2010: Calculation and use of thermodynamic properties.  
#   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
#   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
#
#  McDougall, T.J., D.R. Jackett, F.J. Millero, R. Pawlowicz and 
#   P.M. Barker, 2012: A global algorithm for estimating Absolute Salinity.
#   Ocean Science, 8, 1123-1134.  
#   http://www.ocean-sci.net/8/1123/2012/os-8-1123-2012.pdf 
#
#  See also gsw_SA_from_SP, gsw_deltaSA_atlas
#
#  Reference page in Help browser
#       <a href="matlab:doc gsw_SAAR">doc gsw_SAAR</a>
#  Note that this reference page includes the code contained in gsw_SAAR.
#  We have opted to encode this programme as it is a global standard and 
#  such we cannot allow anyone to change it.
#
#==========================================================================

def SAAR_mod(p, lon, lat):
  """
  Calculates the Absolute Salinity Anomaly Ratio, SAAR, in the open ocean by
  spatially interpolating the global reference data set of SAAR to the location
  of the seawater sample. The Absolute Salinity Anomaly Ratio in the Baltic Sea
  is evaluated separately, since it is a function of Practical Salinity, not of
  space. The present function returns a SAAR of zero for data in the Baltic Sea. 
 
  This function uses version 3.0 of the SAAR look up table (15th May 2011). The
  "gsw_data_v3_0.mat" file should sit in the working directory. Find it in the
  GSW toolbox download if need be. 
  
  *** All entries should be 1d vectors, flatten or reshape it before input (be
  cafeful with ordering!) ***
  
  Adapted from the MATLAB GSW toolbox (http://www.TEOS-10.org)
  
  Inputs (as vectors!):
    p            = pressure                      p      dbar
    lon          = longitude of data                     deg
    lat          = latitude  of data                     deg
    
  Returns (as vector):
    SAAR         = Absolute Salinity Anomaly Ratio       ---
    
  """

  # set any pressures between 0 and -1.5 to be equal to 0 (i.e. the surface)
  p[(p < 0)] = 0

  #--------------------------------------------------------------------------
  # Start of the calculation (extracting from a look up table)
  #--------------------------------------------------------------------------
  data = loadmat('gsw_data_v3_0.mat')
  SAAR_ref   = data["SAAR_ref"][:, :, :] # (p, lat, lon) indexing
  lats_ref   = data["lats_ref"][:, 0]
  long_ref   = data["longs_ref"][:, 0]
  p_ref      = data["p_ref"][:, 0]
  ndepth_ref = data["ndepth_ref"][:, :]  # (lat, lon)    indexing

  # precalculate constants 
  nx  = len(long_ref)
  ny  = len(lats_ref)
  nz  = len(p_ref)
  nyz = ny * nz

  n0  = len(p)

  dlong_ref = long_ref[1] - long_ref[0] 
  dlats_ref  = lats_ref[1] - lats_ref[0]

  indsx0 = floor(1 + (nx - 1) * (lon - long_ref[0]) 
                              / (long_ref[nx - 1] - long_ref[0])
                 )
  indsx0[(indsx0 == nx)] = nx - 1
  indsx0 = int32(indsx0 - 1)

  indsy0 = floor(1 + (ny - 1) * (lat - lats_ref[0])
                              / (lats_ref[ny - 1] - lats_ref[0])
                 )
  indsy0[(indsy0 == ny)] = ny - 1
  indsy0 = int32(indsy0 - 1)

  indsz0 = ones(n0)
  for i in range(1, nz - 2):
    indsz0[(p >= p_ref[i - 1]) & (p < p_ref[i])] = i - 1
  indsz0[(p >= p_ref[nz - 1])] = nz - 1
  indsz0 = int32(indsz0)
       
  indsy0_indsx0_ny = int32(indsy0 + indsx0 * ny + ny)       # python indexing reponsible for extra +ny
  indsn1           = int32(indsy0_indsx0_ny - ny)           # 4 xy grid points surrounding the data
  indsn2           = int32(indsy0_indsx0_ny)
  indsn3           = int32(indsy0_indsx0_ny + 1)
  indsn4           = int32(indsy0_indsx0_ny + (1 - ny))

  # flatten but do it like MATLAB
  nmax_dum       = zeros((4, n0))
  nmax_dum[0, :] = reshape(ndepth_ref, ny * nx, order = "F")[indsn1]
  nmax_dum[1, :] = reshape(ndepth_ref, ny * nx, order = "F")[indsn2]
  nmax_dum[2, :] = reshape(ndepth_ref, ny * nx, order = "F")[indsn3]
  nmax_dum[3, :] = reshape(ndepth_ref, ny * nx, order = "F")[indsn4]

  nmax = nanmax(nmax_dum, axis = 0) - 1    #?? -1 for python index
  
  del nmax_dum

  if any(indsz0 > nmax):                                          # python will complain aboout comparing NaNs
    p[(indsz0 > nmax)] = p_ref[int32(nmax[(indsz0 > nmax)])]      # casts deeper than GK maximum
    indsz0[(indsz0 > nmax)] = nmax[(indsz0 > nmax)] - 1           # have reset p here so have to reset indsz0

  indsyx_tmp = indsy0_indsx0_ny * nz + nz - 1         # python indexing reponsible for extra +nz - 1
  inds0 =  indsz0 + indsyx_tmp - (nyz + nz) + 1       # python indexing took too many off so +1 back in

  data_indices       = zeros((n0, 4))
  data_indices[:, 0] = indsx0
  data_indices[:, 1] = indsy0
  data_indices[:, 2] = indsz0
  data_indices[:, 3] = inds0
  data_inds          = data_indices[:, 2]

  r1 = (lon - long_ref[indsx0]) / (long_ref[indsx0 + 1] - long_ref[indsx0])
  s1 = (lat - lats_ref[indsy0]) / (lats_ref[indsy0 + 1] - lats_ref[indsy0])
  t1 = (p   - p_ref[indsz0])    / (p_ref[indsz0 + 1]    - p_ref[indsz0])
      
  sa_upper = nan * ones(data_inds.shape)
  sa_lower = nan * ones(data_inds.shape)
  SAAR     = nan * ones(data_inds.shape)

  saar_nan = nan * ones((4, n0))

  SAAR_ref = reshape(SAAR_ref, (nx * ny * nz), order = "F")   # flatten but do it like MATLAB
  
  #-----------------------------------------------------------------------------
  # start loop

  for k in range(nz):
    
    # only do something and find indices if there is something that needs doing
    if (indsz0 == k).any():
      inds_k = asarray([i for i in range(len(indsz0)) if indsz0[i] == k])
      indsXYZ = k + indsyx_tmp[inds_k] + 1   # python indexing took too many off so +1 back in
      
      # level k interpolation
      saar = deepcopy(saar_nan)
      saar[0, inds_k] = SAAR_ref[indsXYZ - (nz + nyz)]
      saar[1, inds_k] = SAAR_ref[indsXYZ - (nz      )]
      saar[2, inds_k] = SAAR_ref[indsXYZ             ]
      saar[3, inds_k] = SAAR_ref[indsXYZ - (     nyz)]
      
      inds_pan = asarray([i for i in inds_k 
                      if (abs(lon[i] - 277.6085) <= 17.6085) 
                       & (abs(lat[i] - 9.775) <= 9.775)]
                         )
      
      if inds_pan != []:
        inds = inds_k[inds_pan]
        saar[:, inds] = saar_add_barrier(saar[:, inds],
                                         lon[inds], lat[inds],
                                         long_ref[indsx0[inds]],
                                         lats_ref[indsy0[inds]],
                                         dlong_ref, dlats_ref
                                         )

      if isnan(saar[:, inds_k]).any():
        dummy = isnan(sum(saar[:, inds_k], axis = 0))
        dummy = [i for i in range(len(dummy)) if dummy[i]]
        inds = asarray(inds_k)[dummy]
        saar[:, inds] = saar_add_mean(saar[:, inds])
      
      # level k+1 interpolation
      sa_upper[inds_k] = (
          (1 - s1[inds_k]) * (saar[0, inds_k] + r1[inds_k] * (saar[1, inds_k] - saar[0, inds_k])) +
               s1[inds_k]  * (saar[3, inds_k] + r1[inds_k] * (saar[2, inds_k] - saar[3, inds_k]))
                          )
      
      saar = saar_nan
      saar[0, inds_k] = SAAR_ref[indsXYZ + (1 - nz - nyz)]
      saar[1, inds_k] = SAAR_ref[indsXYZ + (1 - nz      )]
      saar[2, inds_k] = SAAR_ref[indsXYZ +  1            ]
      saar[3, inds_k] = SAAR_ref[indsXYZ + (1 -      nyz)]

      if inds_pan != []:
        inds = inds_k[inds_pan]
        saar[:, inds] = saar_add_barrier(saar[:, inds],
                                         lon[inds], lat[inds],
                                         long_ref[indsx0[inds]],
                                         lats_ref[indsy0[inds]],
                                         dlong_ref, dlats_ref
                                         )
      
      if isnan(saar[:, inds_k]).any():
        dummy = isnan(sum(saar[:, inds_k], axis = 0))
        dummy = [i for i in range(len(dummy)) if dummy[i]]
        inds = asarray(inds_k)[dummy]
        saar[:, inds] = saar_add_mean(saar[:, inds])
          
      sa_lower[inds_k] = (
          (1 - s1[inds_k]) * (saar[0, inds_k] 
                           + r1[inds_k] * (saar[1, inds_k] - saar[0, inds_k])) 
             + s1[inds_k]  * (saar[3, inds_k] 
                           + r1[inds_k] * (saar[2, inds_k] - saar[3, inds_k]))
                          )
      
      if (isfinite(sa_upper[inds_k]) & isnan(sa_lower[inds_k])).any():
        dummy = (isfinite(sa_upper[inds_k]) & isnan(sa_lower[inds_k]))
        inds_diff = asarray([i for i in range(len(dummy)) if dummy[i]])
        sa_lower[inds_k[inds_diff]] = sa_upper[inds_k[inds_diff]]

      SAAR[inds_k] = sa_upper[inds_k] + t1[inds_k] * (sa_lower[inds_k] - sa_upper[inds_k])
      
  # end loop
  #-----------------------------------------------------------------------------
    
  # set the non-finite things to zero
  SAAR[~isfinite(SAAR)] = 0
  
  return SAAR
  
###########################################################################
# saar_add_mean (was gsw_saar_add_mean)
#==========================================================================
#
# USAGE:
#  SAAR = gsw_saar_add_mean(saar)
#
# DESCRIPTION:
#  Replaces NaN's with nanmean of the 4 adjacent neighbours
#
# INPUT:
#  saar  =  Absolute Salinity Anomaly Ratio of the 4 adjacent neighbours  
#                                                              [ unitless ]
#
# OUTPUT:
#  SAAR  =  nanmean of the 4 adjacent neighbours               [ unitless ]
#
# AUTHOR: 
#  David Jackett
#
# MODIFIED:
#  Paul Barker and Trevor McDougall
#  Aknowlegments. Matlab programming assisance from Sjoerd Groeskamp.
#
# VERSION NUMBER: 3.05 (27th January 2015)
#
# REFERENCES:
#  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
#   seawater - 2010: Calculation and use of thermodynamic properties.  
#   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
#   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
#
#  McDougall, T.J., D.R. Jackett, F.J. Millero, R. Pawlowicz and 
#   P.M. Barker, 2012: A global algorithm for estimating Absolute Salinity.
#   Ocean Science, 8, 1123-1134.  
#   http://www.ocean-sci.net/8/1123/2012/os-8-1123-2012.pdf 
#
#  The software is available from http://www.TEOS-10.org
#
#==========================================================================

def saar_add_mean(saar):
  """
  Replaces NaN's with nanmean of the 4 adjacent neighbours
  
  Inputs (as vectors!):
    SAAR         = Absolute Salinity Anomaly Ratio       ---
    
  Returns (as vector):
    SAAR         = Absolute Salinity Anomaly Ratio       ---
    
  """
  saar_mean = mean(saar, axis = 0)
  inds_nan = asarray([i for i in range(len(saar_mean)) if isnan(saar_mean[i])])
  no_nan = len(inds_nan)
  for kk in range(no_nan):
      col = inds_nan[kk]
      Inn = [i for i in range(4) if ~isnan(saar[i, col])]
      if Inn != []:
          saar[isnan(saar[:, col]), col] = sum(saar[Inn,col]) / len(Inn)
          
  return saar

###########################################################################
# saar_add_barrier
#==========================================================================
#
# USAGE:
#  SAAR = gsw_saar_add_barrier(saar,long,lat,longs_ref,lats_ref,dlong_ref,dlats_ref)
#
# DESCRIPTION:
#  Adds a barrier through Central America (Panama) and then averages
#  over the appropriate side of the barrier
#
# INPUT:
#  saar        =  Absolute Salinity Anomaly Ratio                          [ unitless ]
#  long        =  Longitudes of data in decimal degrees east               [ 0 ... +360 ]
#  lat         =  Latitudes of data in decimal degrees north               [ -90 ... +90 ]
#  longs_ref   =  Longitudes of regular grid in decimal degrees east       [ 0 ... +360 ]
#  lats_ref    =  Latitudes of regular grid in decimal degrees north       [ -90 ... +90 ]
#  dlong_ref  =  Longitude difference of regular grid in decimal degrees   [ deg longitude ]
#  dlats_ref   =  Latitude difference of regular grid in decimal degrees   [ deg latitude ]
#
# OUTPUT:
#  SAAR        =  Absolute Salinity Anomaly Ratio                          [ unitless ]
#
# AUTHOR: 
#  David Jackett
#
# MODIFIED:
#  Paul Barker and Trevor McDougall
#
# VERSION NUMBER: 3.05 (27th January 2015)
#
# REFERENCES:
#  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
#   seawater - 2010: Calculation and use of thermodynamic properties.  
#   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
#   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
#
#  McDougall, T.J., D.R. Jackett, F.J. Millero, R. Pawlowicz and 
#   P.M. Barker, 2012: A global algorithm for estimating Absolute Salinity.
#   Ocean Science, 8, 1123-1134.  
#   http://www.ocean-sci.net/8/1123/2012/os-8-1123-2012.pdf 
#
#  The software is available from http://www.TEOS-10.org
#
#==========================================================================

def saar_add_barrier(saar, lon, lat, long_ref, lats_ref, dlong_ref, dlats_ref):
  """
  Adds a barrier through Central America (Panama) and then averages over the
  appropriate side of the barrier
  
  Inputs (as vectors!):
    SAAR         = Absolute Salinity Anomaly Ratio       ---
    lon          = longitude                             deg
    lat          = latitude                              deg
    long_ref     = longitude reference from lookup       deg
    lats_ref     = latitutde reference from lookup       deg
    dlong_ref    = spacing of long_ref                   deg
    dlats_ref    = spacing of lats_ref                   deg
    
  Returns (as vector):
    SAAR         = Absolute Salinity Anomaly Ratio       ---
    
  """

  long_pan = asarray([260.0000, 272.5900, 276.5000, 278.6500, 280.7300, 295.2170])
  lats_pan = asarray([ 19.5500,  13.9700,   9.6000,   8.1000,   9.3300,   0])

  f = interp1d(long_pan, lats_pan, bounds_error = False, fill_value = nan)
  lats_lines0 = f(lon)
  lats_lines1 = f(long_ref)
  lats_lines2 = f(long_ref + dlong_ref)
  
  above_line = zeros(4)

  for k0 in range(len(lon)):
    if (lats_lines0[k0] <= lat[k0]):
        above_line0 = 1
    else:
        above_line0 = 0

    if (lats_lines1[k0] <= lats_ref[k0]):
        above_line[0] = 1

    if (lats_lines1[k0] <= (lats_ref[k0] + dlats_ref)):
        above_line[3] = 1

    if (lats_lines2[k0] <= lats_ref[k0]):
        above_line[1] = 1

    if (lats_lines2[k0] <= (lats_ref[k0] + dlats_ref)):
        above_line[2] = 1

    # indices of different sides of CA line
    saar[above_line != above_line0, k0] = nan

  if isnan(saar).any():
    saar = saar_add_mean(saar)

  return saar
  
################################################################################
# p_from_z (was gsw_p_from_z)            pressure from height (75-term equation)
#==========================================================================
#
# USAGE:
#  p = gsw_p_from_z(z,lat,{geo_strf_dyn_height},{sea_surface_geopotental})
#
# DESCRIPTION:
#  Calculates sea pressure from height using computationally-efficient 
#  75-term expression for density, in terms of SA, CT and p (Roquet et al.,
#  2015).  Dynamic height anomaly, geo_strf_dyn_height, if provided,
#  must be computed with its p_ref = 0 (the surface). Also if provided,
#  sea_surface_geopotental is the geopotential at zero sea pressure. This 
#  function solves Eqn.(3.32.3) of IOC et al. (2010) iteratively for p.  
#
#  Note. Height (z) is NEGATIVE in the ocean.  Depth is -z.  
#    Depth is not used in the GSW computer software library. 
#
#  Note that this 75-term equation has been fitted in a restricted range of 
#  parameter space, and is most accurate inside the "oceanographic funnel" 
#  described in McDougall et al. (2003).  The GSW library function 
#  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if 
#  some of one's data lies outside this "funnel".  
#
# INPUT:
#  z  =  height                                                       [ m ]
#   Note. At sea level z = 0, and since z (HEIGHT) is defined 
#     to be positive upwards, it follows that while z is 
#     positive in the atmosphere, it is NEGATIVE in the ocean.
#  lat  =  latitude in decimal degrees north                [ -90 ... +90 ]
#   
# OPTIONAL:
#  geo_strf_dyn_height = dynamic height anomaly                 [ m^2/s^2 ]
#    Note that the reference pressure, p_ref, of geo_strf_dyn_height must
#     be zero (0) dbar.
#  sea_surface_geopotental = geopotential at zero sea pressure  [ m^2/s^2 ]
#
#  lat may have dimensions 1x1 or Mx1 or 1xN or MxN, where z is MxN.
#  geo_strf_dyn_height and geo_strf_dyn_height, if provided, must have 
#  dimensions MxN, which are the same as z.
#
# OUTPUT:
#   p  =  sea pressure                                             [ dbar ]
#         ( i.e. absolute pressure - 10.1325 dbar )
#
# AUTHOR: 
#  Trevor McDougall, Claire Roberts-Thomson and Paul Barker. 
#                                                      [ help@teos-10.org ]
#
# VERSION NUMBER: 3.05 (27th January 2015)
#
# REFERENCES:
#  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
#   seawater - 2010: Calculation and use of thermodynamic properties.
#   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
#   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
#
#  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003: 
#   Accurate and computationally efficient algorithms for potential 
#   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
#   pp. 730-741.
#
#  McDougall, T.J., and S.J. Wotherspoon, 2013: A simple modification of 
#   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied 
#   Mathematics Letters, 29, pp. 20-25.  
#
#  Moritz, H., 2000: Geodetic reference system 1980. J. Geodesy, 74, 
#   pp. 128-133.
#
#  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
#   polynomial expressions for the density and specifc volume of seawater
#   using the TEOS-10 standard. Ocean Modelling, 90, pp. 29-43.
#
#  Saunders, P.M., 1981: Practical conversion of pressure to depth. 
#   Journal of Physical Oceanography, 11, pp. 573-574.
#
#  This software is available from http://www.TEOS-10.org
#
#==========================================================================

def p_from_z(z, lat):

  #--------------------------------------------------------------------------
  # Start of the calculation
  #--------------------------------------------------------------------------

  geo_strf_dyn_height = 0
  sea_surface_geopotential = 0
  
  db2Pa = 1e4
#  gamma = 2.26e-7 # If the graviational acceleration were to be regarded as 
                   # being depth-independent, which is often the case in 
                   # ocean models, then gamma would be set to be zero here,
                   # and the code below works perfectly well.  
  gamma = 0
  deg2rad = pi / 180
  sinlat = sin(lat * deg2rad)
  sin2 = sinlat* sinlat
  gs = 9.780327 * (  1.0 + ( 5.2792e-3 + (2.32e-5 * sin2) ) * sin2  )

  # get the first estimate of p from Saunders (1981)
  c1 =  5.25e-3 * sin2 + 5.92e-3
  p  = -2 * z / (  (1 - c1) + sqrt( (1 - c1) *(1 - c1) + 8.84e-6 * z)  )
  # end of the first estimate from Saunders (1981)

  df_dp = db2Pa * specvol_SSO_0(p) # initial value of the derivative of f

  f = enthalpy_SSO_0(p) + gs * ( (z - 0.5 * gamma * z * z)
                               - geo_strf_dyn_height
                               - sea_surface_geopotential
                               )
  p_old = p
  p     = p_old - f / df_dp
  p_mid = 0.5 * (p + p_old)
  df_dp = db2Pa * specvol_SSO_0(p_mid)
  p     = p_old - f / df_dp

  # After this one iteration through this modified Newton-Raphson iterative
  # procedure (McDougall and Wotherspoon, 2013), the remaining error in p is 
  # at computer machine precision, being no more than 1.6e-10 dbar.
  
  return p
  
################################################################################
# enthalpy_SSO_0 (was gsw_enthalpy_SSO_0)               enthalpy at (SSO,CT=0,p)
#                                                            (75-term eqn.)
#==========================================================================
#  This function calculates enthalpy at the Standard Ocean Salinity, SSO, 
#  and at a Conservative Temperature of zero degrees C, as a function of
#  pressure, p, in dbar, using a streamlined version of the 75-term 
#  computationally-efficient expression for specific volume, that is, a 
#  streamlined version of the code "gsw_enthalpy(SA,CT,p)".
#
# VERSION NUMBER: 3.05 (27th January 2015)
#
# REFERENCES:
#  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
#   polynomial expressions for the density and specifc volume of seawater
#   using the TEOS-10 standard. Ocean Modelling.
#
#==========================================================================

def enthalpy_SSO_0(p):

  z = p * 1e-4

  h006 = -2.1078768810e-9
  h007 =  2.8019291329e-10

  dynamic_enthalpy_SSO_0_p = (
        z * (9.726613854843870e-4  + z * (-2.252956605630465e-5
      + z * (2.376909655387404e-6  + z * (-1.664294869986011e-7
      + z * (-5.988108894465758e-9 + z * (h006 + h007 * z))))))
                              )

  enthalpy_SSO_0_out = dynamic_enthalpy_SSO_0_p * 1e8     #Note. 1e8 = db2Pa*1e4
  
  return enthalpy_SSO_0_out
  
################################################################################  
# specvol_SSO_0 (was gsw_specvol_SSO_0)          specific volume at (SSO,CT=0,p)
#                                                        (75-term equation)
#==========================================================================
#  This function calculates specifc volume at the Standard Ocean Salinity,
#  SSO, and at a Conservative Temperature of zero degrees C, as a function 
#  of pressure, p, in dbar, using a streamlined version of the 75-term CT
#  version of specific volume, that is, a streamlined version of the code
#  "gsw_specvol(SA,CT,p)".
#
# VERSION NUMBER: 3.05 (27th January 2015)
#
# REFERENCES:
#  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
#   polynomial expressions for the density and specifc volume of seawater
#   using the TEOS-10 standard. Ocean Modelling.
#
#==========================================================================
 
def specvol_SSO_0(p):
  z = p * 1e-4

  v005 = -1.2647261286e-8
  v006 =  1.9613503930e-9
          
  specvol_SSO_0_out = (
              9.726613854843870e-04 + z * (-4.505913211160929e-05
      + z * ( 7.130728965927127e-06 + z * (-6.657179479768312e-07
      + z * (-2.994054447232880e-08 + z * (v005 + v006 * z)))))
                       )

  return specvol_SSO_0_out
  
  
  
  
  
  
  
  
  
