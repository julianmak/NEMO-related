#!/usr/bin/env python3

# JM: 27 Nov 2018
# draw a clock

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

# units in radians

def draw_clock(yyyymmdd, clock_color = "xkcd:grey", progress_color = "Spectral", 
               fontsize = 14, ax = None, uniform_months = False):
  """
  fairly dumb way of drawing a clock, takes input as yyyymmdd and draws a clock
  through drawing some filled circles with contourf in polar plot
  
  the plot will not automatically scale but the "fontsize" number can be modified
  accordingly to make it scale, so test this with some sample images first
  
  input:
     yyyymmdd         string of yyyy/mm/dd, how you grab from data is up to you
     clock_color      default is "xkcd:grey", modify accordingly as RGB, hex, python words etc.
     progress_color   default is orange progress lime background from "Spectral", change it
                        by inputing a colormap if you like
     fontsize         default 14, modify this depending on clock size
     ax               subplot axes to plot it, suggestion is in the parent axes do say
     
                        a = plt.axes([0.95, .6, .2, .2], polar = True)
                        draw_clock("19510630", ax = a, fontsize = 10)
            
     uniform_months   if False then assumes no leap years, otherwise 30 days a month
  """
  
  if ax is None:
    ax = plt.axes(projection = 'polar')

  # set up the clock as a circle
  ax.set_theta_offset(np.pi / 2.0)    # start counting at 12 o'clock
  ax.set_theta_direction("clockwise") # go clockwise
  ax.set_xticks([])                   # kill all the ticks
  ax.set_rticks([])
  ax.set_rlim(0, 1)                   # set the clockface

  ax.set_facecolor(clock_color)

  # set up an array to plot the invariant parts
  outer_line = 0.70
  inner_line = 0.45
  theta_vec = np.linspace(0, 2 * np.pi, 71)
  r_vec = np.linspace(inner_line, outer_line, 31)
  theta, r = np.meshgrid(theta_vec, r_vec)

  # set up some settings (hand tuned for now...)
  months = {}
  months["label"] = ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"]
  months["theta"] = np.arange(0, 12, 1) * np.pi / 6.0 + np.pi / 12.0
  if uniform_months:
    months["days"] = [30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30] # assume 30 days a month
  else:
    months["days"] = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31] # assume no leap years
  # work out the date in radians
  year  = int(yyyymmdd[0:4])
  month = int(yyyymmdd[4:6])
  day   = int(yyyymmdd[6::])
  if (year > 9999) or (year < 0):
    print("year grabbed is %.4d and is out of bounds?" % year)
    return
  if (month > 12) or (month < 0):
    print("month grabbed is %.2d and is out of bounds?" % month)
    return
  if (day > months["days"][month - 1]) or (day < 0):
    print("month grabbed is %.2d" % month)
    print("but date grabbed is %.2d so is out of bounds?" % day)
    return
  date_in_rad = (month - 1) * np.pi / 6.0 + (day / months["days"][month - 1]) * np.pi / 6.0

  ax.plot(theta_vec, inner_line * np.ones(theta_vec.shape), 'k')
  ax.plot(theta_vec, outer_line * np.ones(theta_vec.shape), 'k')
  ax.plot(theta_vec, 1.0 * np.ones(theta_vec.shape), 'k', linewidth = 2)
  ax.text(3 * np.pi / 2.0, 0.0, "%.4d" % year, 
          ha = 'center', va = 'center', fontsize = fontsize)
  for month in range(12):
    ax.plot([month * np.pi / 6.0, month * np.pi / 6.0], [outer_line, 1.0], 'k-')
    ax.text(months["theta"][month], 0.85, months["label"][month], 
            ha = 'center', va = 'center', fontsize = fontsize)

  filled_region = np.where(theta < date_in_rad + 0.01, -1, 1) # ad a little increment to push the contour over
  ax.contourf(theta, r, filled_region, levels = np.linspace(-2, 2, 3), cmap = progress_color)

def hex_to_rgb(color_hex):
    
  color_rgb = tuple(int(color_hex[i:i+2], 16) / 255 for i in (0, 2 ,4))
  
  return color_rgb

def hex_duple_colormap(color_hex1, color_hex2, sample = False):
  
  colors = [hex_to_rgb(color_hex1), hex_to_rgb(color_hex2)]
  cmap = LinearSegmentedColormap.from_list("murp", colors, N = 2)
  
  if sample:
    x = np.arange(0, np.pi, 0.1)
    y = np.arange(0, 2*np.pi, 0.1)
    X, Y = np.meshgrid(x, y)
    Z = np.cos(X) * np.sin(Y) * 10

    ax = plt.axes()
    im = ax.imshow(Z, interpolation='nearest', origin='lower', cmap = cmap)
    fig.colorbar(im, ax = ax)
    plt.show()
    
  return cmap
  
def convert_nemo_times(time_in_sec, uniform_months = False):
  """
  Converts NEMO outputs of seconds since 1900 Jan 01 to a yyyymmdd date
  
  Can take a vector or number and returns a list
  
  NOTE: NOT TESTED WITH NON-UNIFORM MONTH DATA [Apr 19 2020]
  """
  if uniform_months: # assume 30 days per month
    days_in_month = np.array([30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30])
    days_in_year = 360
  else: # assume no leap years
    days_in_month = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
    days_in_year = 365
    
  cumsum_days = np.cumsum(days_in_month)
  sec_per_yr = 3600 * 24 * days_in_year
  time_in_sec = (time_in_sec + sec_per_yr * 1900) / sec_per_yr
  
  yyyymmdd = []
  for kt in range(len(time_in_sec)):
    year  = int(time_in_sec[kt])
    month = np.argmax((time_in_sec[kt] - year) * days_in_year < cumsum_days) + 1
    day = int((time_in_sec[kt] - year) * days_in_year - cumsum_days[month - 1] + cumsum_days[0])
    yyyymmdd.append(f"{year:04d}{month:02d}{day:02d}")
    
  return yyyymmdd
  
