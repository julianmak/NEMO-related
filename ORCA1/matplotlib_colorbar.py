#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""
Convenience routines for maptlotlib colour bars. Also includes some ParaView
colour maps.
"""

import copy
import matplotlib.ticker
import unittest
import numpy as np

__all__ = \
  [
    "BoundedLocator",
    "SciFormatter",
    "cold_and_hot",
    "cold_and_hot_white",
    "cool_to_warm",
    "grayscale",
    "rainbow_desaturated",
    "rainbow_desaturated_white",
    "BluWhiOraRed"
  ]

def _cdict(rgb):
  return {"red"  :tuple([(v[0], v[1][0] / 255, v[1][0] / 255) for i, v in enumerate(rgb)]),
          "green":tuple([(v[0], v[1][1] / 255, v[1][1] / 255) for i, v in enumerate(rgb)]),
          "blue" :tuple([(v[0], v[1][2] / 255, v[1][2] / 255) for i, v in enumerate(rgb)])}

################################################################################

# Numbers extracted from "Rainbow Desaturated" ParaView 4.0.1 Preset
rainbow_desaturated = [(0, (71, 71, 219)),
                       (0.143, (0, 0, 92)),
                       (0.285, (0, 255, 255)),
                       (0.429, (0, 128, 0)),
                       (0.571, (255, 255, 0)),
                       (0.714, (255, 97, 0)),
                       (0.857, (107, 0, 0)),
                       (1, (224, 77, 77))]
rainbow_desaturated = matplotlib.colors.LinearSegmentedColormap("rainbow_desaturated", _cdict(rainbow_desaturated), 2 ** 12)

# Numbers extracted from "Rainbow Desaturated" ParaView 4.0.1 Preset
rainbow_desaturated_white = [(0, (71, 71, 219)),
                             (0.143, (0, 0, 92)),
                             (0.285, (0, 255, 255)),
                             (0.429, (0, 128, 0)),
                             (0.5, (255, 255, 255)),
                             (0.571, (255, 255, 0)),
                             (0.714, (255, 97, 0)),
                             (0.857, (107, 0, 0)),
                             (1, (224, 77, 77))]
rainbow_desaturated_white = matplotlib.colors.LinearSegmentedColormap("rainbow_desaturated_white", _cdict(rainbow_desaturated_white), 2 ** 12)

# Numbers extracted from "Cold and Hot" ParaView 4.0.1 Preset
cold_and_hot = [(0, (0, 255, 255)),
                (0.45, (0, 0, 255)),
                (0.5, (0, 0, 128)),
                (0.55, (255, 0, 0)),
                (1, (255, 255, 0))]
cold_and_hot = matplotlib.colors.LinearSegmentedColormap("cold_and_hot", _cdict(cold_and_hot), 2 ** 12)

# Modified version of "Cold and Hot" ParaView 4.0.1 Preset
cold_and_hot_white = [(0, (0, 255, 0)),
                      (0.3, (0, 255, 255)),
                      (0.47, (0, 0, 255)),
                      (0.5, (255, 255, 255)),
                      (0.53, (255, 0, 0)),
                      (0.7, (255, 0, 255)),
                      (1, (255, 255, 0))]
cold_and_hot_white = matplotlib.colors.LinearSegmentedColormap("cold_and_hot_white", _cdict(cold_and_hot_white), 2 ** 12)

# Numbers extracted from "Cool to Warm" 4.1.0 Paraview Preset
cool_to_warm = [(0, (59, 76, 192)),
                (0.5, (221, 221, 221)),
                (1.0, (180, 4, 38))]
cool_to_warm = matplotlib.colors.LinearSegmentedColormap("cool_to_warm", _cdict(cool_to_warm), 2 ** 12)

# Numbers extracted from "Grayscale" ParaView 4.0.1 Preset
grayscale = [(0, (0, 0, 0)),
             (1, (255, 255, 255))]
grayscale = matplotlib.colors.LinearSegmentedColormap("grayscale", _cdict(grayscale), 2 ** 12)

# Numbers extracted from "BluWhiOraRed" NCL UCAR
BluWhiOraRed = [(0, (27, 44, 98)),
             (0.09, (36, 89, 167)),
             (0.18, (87, 164, 216)),
             (0.27, (157, 218, 247)),
             (0.36, (215, 239, 249)),
             (0.49, (255, 255, 255)),
             (0.51, (255, 255, 255)),
             (0.63, (252, 235, 155)),
             (0.72, (254, 176, 50)),
             (0.81, (245, 106, 41)),
             (0.90, (211, 31, 40)),
             (1, (141, 21, 25))]

BluWhiOraRed = matplotlib.colors.LinearSegmentedColormap("BluWhiOraRed", _cdict(BluWhiOraRed), 2 ** 12)

################################################################################

class SciFormatter(matplotlib.ticker.Formatter):
  """
  Scientific notation Ticker Formatter.

  Constructor arguments:

    sf: (Optional) Number of significant figures
    strip_zeros: (Optional) Whether to strip trailing zeros
  """
  
  def __init__(self, sf = 7, strip_zeros = True):
    self.__sf = sf
    self.__sz = strip_zeros

    return

  def __call__(self, x, pos = None):
    if x == 0:
      return "$0$"
    elif abs(x) < 1:
      exp = int(np.log10(abs(x))) - 1
    else:
      exp = int(np.log10(abs(x)))

    v = ("%%.%if" % (self.__sf - 1)) % (x * (10 ** -exp))
    if v.startswith("10.") or v.startswith("-10."):
      exp += 1
      v = ("%%.%if" % (self.__sf - 1)) % (x * (10 ** -exp))
    assert(len(v.split(".")[0]) == 2 if v.startswith("-") else 1)
    
    if self.__sz:
      i = len(v) - 1
      while v[i] == "0":
        i -= 1
      v = v[:i + 1]
    if v.endswith("."):
      v = v[:-1]
      
    if v == "0":
      return "$0$"
    else:
      return "$%s \\times 10^{%i}$" % (v, exp)

class BoundedLocator(matplotlib.ticker.Locator):
  """
  Ticker Locator which includes field bounds.

  Constructor arguments:
    pad: Fraction of colour bar over which interior ticks should be removed
  """
  
  def __init__(self, pad = 0.1):
    matplotlib.ticker.Locator.__init__(self)
    self.__pad = pad

    return

  def __call__(self):
    return self.tick_values(*self.axis.get_data_interval())

  def tick_values(self, vmin, vmax):    
    # Round the lower bound up and the upper bound down to 1 sf. This determines
    # our initial internal tick value bounds.
    vrange = vmax - vmin
    if abs(vrange) < 1:
      exp = int(np.log10(abs(vrange))) - 1
    else:
      exp = int(np.log10(abs(vrange)))
    v1 = np.ceil(vmin * (10 ** -exp))
    v2 = np.floor(vmax * (10 ** -exp))
    vn = int(v2 - v1 + 1.5)

    # If we have too few tick values, round the bounds to 2 sf
    def vp():
      vp = vn
      if abs(v1 - vmin) < (vmax - vmin) * self.__pad:
        vp -= 1
      if abs(v2 - vmax) < (vmax - vmin) * self.__pad:
        vp -= 1
      return vp
    if vp() < 4:
      exp -= 1
      v1 = np.ceil(vmin * (10 ** -exp))
      v2 = np.floor(vmax * (10 ** -exp))
    vn = int(v2 - v1 + 1.5)

    # If we have too many tick values, round more aggressively
    d = 2
    while vp() > 8:
      v1 = d * np.ceil(v1 / d)
      v2 = d * np.floor(v2 / d)
      vn = int(((v2 - v1) / d) + 1.5)
      d *= 2

    # The internal tick values
    vi = (np.linspace(v1, v2, vn) * (10 ** exp)).tolist()

    # Remove tick values that are very close to the value minima and maxima,
    # as these may clash with the value bounds
    if abs(vi[0] - vmin) < (vmax - vmin) * self.__pad:
      del(vi[0])
    if abs(vi[-1] - vmax) < (vmax - vmin) * self.__pad:
      del(vi[-1])

    # Add ticks for the value bounds
    return [vmin] + vi + [vmax]

class matplotlib_colorbar_UnitTests(unittest.TestCase):
  def testSciFormatter(self):
    s = SciFormatter()
    self.assertEquals(s(3.101e-6), "$3.101 \\times 10^{-6}$")
    self.assertEquals(s(-3e-6), "$-3 \\times 10^{-6}$")
    self.assertEquals(s(3e6), "$3 \\times 10^{6}$")
    self.assertEquals(s(1e6), "$1 \\times 10^{6}$")
    self.assertEquals(s(-1e6), "$-1 \\times 10^{6}$")
    self.assertEquals(s(0), "$0$")

    s = SciFormatter(3)
    self.assertEquals(s(3.101e-6), "$3.1 \\times 10^{-6}$")
    self.assertEquals(s(-3e-6), "$-3 \\times 10^{-6}$")
    self.assertEquals(s(3e6), "$3 \\times 10^{6}$")
    self.assertEquals(s(1e6), "$1 \\times 10^{6}$")
    self.assertEquals(s(-1e6), "$-1 \\times 10^{6}$")
    self.assertEquals(s(0), "$0$")

    return

if __name__ == "__main__":
  unittest.main()
