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
#  !!----------------------------------------------------------------------
#  !!----------------------------------------------------------------------
#  !!   routines      : description
#  !!   sigma0     : compute sigma-0 
#  !!   eosbn2     : compute Brunt Vaissala Frequency
#  !!   sigmai     : compute sigma-i ( refered to a depth given in argument
#  !!   albet      : Compute the ratio alpha/beta ( Thermal/haline exapnsion)
#  !!   beta       : compute beta (haline expension)
#  !!----------------------------------------------------------------------

from numba import jit
from numpy import zeros, sqrt, abs

################################################################################
#    !!---------------------------------------------------------------------
#    !!                  ***  FUNCTION sigma0  ***
#    !!
#    !! ** Purpose : Compute the in situ density (ratio rho/rau0) and the
#    !!              potential volumic mass (Kg/m3) from potential temperature 
#    !!              and salinity fields using an equation of state defined 
#    !!              through the namelist parameter neos.
#    !!
#    !! ** Method  : Jackett and McDougall (1994) equation of state.
#    !!              The in situ density is computed directly as a function of
#    !!              potential temperature relative to the surface (the opa t
#    !!              variable), salt and pressure (assuming no pressure variation
#    !!              along geopotential surfaces, i.e. the pressure p in decibars
#    !!              is approximated by the depth in meters.
#    !!              prd(t,s,p) = ( rho(t,s,p) - rau0 ) / rau0
#    !!              rhop(t,s)  = rho(t,s,0)
#    !!              with pressure                      p        decibars
#    !!                   potential temperature         t        deg celsius
#    !!                   salinity                      s        psu
#    !!                   reference volumic mass        rau0     kg/m**3
#    !!                   in situ volumic mass          rho      kg/m**3
#    !!                   in situ density anomalie      prd      no units
#    !!
#    !!----------------------------------------------------------------------
def sigma0(ptem, psal):
  """
  Compute the in situ density (ratio rho/rau0) and the potential volumic mass
  (Kg/m3) from potential temperature and salinity fields using an equation of
  state defined through the namelist parameter neos.
  
  JM: changed input so that array sizes by default match the input array sizes
  
  Inputs:
    ptem       = potential temperature         t        deg celsius
    psal       = salinity                      s        psu
    
  Returns:
    sigma0_out = in situ density using Jackett and McDougall (1994) equation of state
    
  """
  zrau0 = 1000.0

  # compute volumic mass pure water at atm pressure
  zr1 = ( ( ( ( ( 6.536332e-9 * ptem - 1.120083e-6 ) * ptem + 1.001685e-4) * ptem   
       -9.095290e-3 ) * ptem + 6.793952e-2 ) * ptem + 999.842594
        )
  # seawater volumic mass atm pressure
  zr2 = ( ( ( ( 5.3875e-9 * ptem - 8.2467e-7 ) * ptem + 7.6438e-5 ) * ptem
       -4.0899e-3 ) * ptem + 0.824493
        )
  zr3 = ( -1.6546e-6 * ptem + 1.0227e-4 ) * ptem - 5.72466e-3
  zr4 = 4.8314e-4
  
  # potential volumic mass (reference to the surface)
  sigma0_out = ( zr4 * psal + zr3 * zws + zr2 ) * psal + zr1 - zrau0

  return sigma0_out
  
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
def sigmantr(ptem, psal):
  """
  Compute the neutral volumic mass (kg/m3) from known potential temperature and
  salinity fields using an equation of state. 
  
  JM: Uchanged input so that array sizes by default match the input array sizes
  
  Inputs:
    ptem         = potential temperature         t        deg celsius
    psal         = salinity                      s        psu
    
  Returns:
    sigmantr_out = neutral density using Jackett and McDougall (2005) equation of state
    
  """
  zrau0 = 1000.0
  
  zws = sqrt( abs(psal) )
  
  # Numerator
  # T-Polynome
  zr1 = ( ( (-4.3159255086706703e-4 * ptem + 8.1157118782170051e-2) * ptem 
           + 2.2280832068441331e-1 ) * ptem + 1002.3063688892480e0
         )
  # S-T Polynome
  zr2 = ( (-1.7052298331414675e-7 * psal - 3.1710675488863952e-3 * ptem 
        - 1.0304537539692924e-4) * psal
        )
  # Denominator
  # T-Polynome
  zr3 = ( ( ( (-2.3850178558212048e-9 * ptem -1.6212552470310961e-7) * ptem
        + 7.8717799560577725e-5 ) * ptem + 4.3907692647825900e-5  ) * ptem + 1.0e0
        )
  # S-T Polynome
  zr4 = ( ( ( -2.2744455733317707e-9 * ptem * ptem + 6.0399864718597388e-6) * ptem 
        - 5.1268124398160734e-4 ) * psal
        )
  # S-T Polynome
  zr5 = ( -1.3409379420216683e-9 * ptem * ptem - 3.6138532339703262e-5) * psal * zws

  # Neutral density
  sigmantr_out = ( zr1 + zr2 ) / ( zr3 + zr4 + zr5 ) - zrau0
  
  return sigmantr_out
  
################################################################################
#    !! --------------------------------------------------------------------
#    !! ** Purpose :   Compute the  density referenced to pref (ratio rho/rau0) 
#    !!       from potential temperature and
#    !!      salinity fields using an equation of state defined through the
#    !!     namelist parameter neos.
#    !!
#    !! ** Method  :
#    !!       Jackett and McDougall (1994) equation of state.
#    !!         the in situ density is computed directly as a function of
#    !!         potential temperature relative to the surface (the opa t
#    !!         variable), salt and pressure (assuming no pressure variation
#    !!         along geopotential surfaces, i.e. the pressure p in decibars
#    !!         is approximated by the depth in meters.
#    !!              prd(t,s,p) = ( rho(t,s,p) - rau0 ) / rau0
#    !!              rhop(t,s)  = rho(t,s,0)
#    !!         with pressure                      p        decibars
#    !!              potential temperature         t        deg celsius
#    !!              salinity                      s        psu
#    !!              reference volumic mass        rau0     kg/m**3
#    !!              in situ volumic mass          rho      kg/m**3
#    !!              in situ density anomalie      prd      no units
#    !! --------------------------------------------------------------------
@jit(nopython = True)
def sigmai_dep(ptem, psal, pref):
  """
  Compute the density referenced to pref (ratio rho/rau0) from potential
  temperature and salinity fields using an equation of state defined through the
  namelist parameter neos.
  
  JM: changed input so that array sizes by default match the input array sizes
      combined sigmai_dep2d and sigmai_dep (here it doesn't care, and reference
      can be 2d array matching input size or just a number)
  
  Inputs:
    ptem           = potential temperature         t        deg celsius
    psal           = salinity                      s        psu
    
  Returns:
    sigmai_dep_out = density referenced to pref using Jackett and McDougall
                     (2005) equation of state
    
  """
  zr4   = 4.8313e-4
  zd    =-2.042967e-2
  zrau0 = 1000.e0
  
  sigmai_dep_out = zeros(psal.shape)
  
  # ?? for whatever reason sqrt(abs(psal)) seems to kick up a fuss when arrays
  #    exceed a certain size...??? otherwise this could be vectorised
  # TODO: if pref is a number, broadcast it into a 2d field
  
  for jj in range(psal.shape[0]): # python indexing
    for ji in range(psal.shape[1]):
    
      ztem = ptem[jj, ji]
      zsal = psal[jj, ji]
      zws = sqrt( abs(psal[jj, ji]) )
    
      # Compute the volumic mass of pure water at atmospheric pressure.
      zr1 = ( ( ( ( (6.536332e-9 * ztem - 1.120083e-6) * ztem + 1.001685e-4 )
                * ztem - 9.095290e-3 ) * ztem + 6.793952e-2 ) * ztem + 999.842594e0
            )

      # Compute the seawater volumic mass at atmospheric pressure.
      zr2 = ( ( ( ( 5.3875e-9 * ztem - 8.2467e-7) * ztem + 7.6438e-5)
           * ztem - 4.0899e-3) * ztem + 0.824493e0
            )

      zr3 = (-1.6546e-6 * ztem + 1.0227e-4) * ztem - 5.72466e-3

      # Compute the potential volumic mass (referenced to the surface).
      zrhop = (zr4 * zsal + zr3 * zws + zr2) * zsal + zr1

      # Compute the compression terms.
      ze  = (-3.508914e-8 * ztem - 1.248266e-8) * ztem - 2.595994e-6

      zbw = (1.296821e-6 * ztem - 5.782165e-9) * ztem + 1.045941e-4

      zb  = zbw + ze * zsal

      zc  = (-7.267926e-5 * ztem + 2.598241e-3) * ztem + 0.1571896e0

      zaw = ( ( (5.939910e-6 * ztem + 2.512549e-3) * ztem - 0.1028859e0 ) 
             * ztem - 4.721788e0
             )

      za  = (zd * zws + zc) * zsal + zaw

      zb1 = (-0.1909078e0 * ztem + 7.390729e0) * ztem - 55.87545e0

      za1 = ( ( (2.326469e-3 * ztem + 1.553190e0) * ztem - 65.00517e0)
             * ztem + 1044.077e0
             )

      zkw = ( ( ( (-1.361629e-4 * ztem - 1.852732e-2) * ztem - 30.41638e0)
             * ztem + 2098.925e0) * ztem + 190925.60
             )

      zk0 = (zb1 * zws + za1) * zsal + zkw

      # Compute the potential density anomaly.
      sigmai_dep_out[jj, ji] = ( zrhop / (1.0e0 - pref / 
                                         ( zk0 - pref * (za - pref * zb) ) )
                               - zrau0
                               )
  
  return sigmai_dep_out
