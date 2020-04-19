#!/usr/bin/env python3

# JM: 19 Apr 2020
# some generic plotting commands

import matplotlib.pyplot as plt
import numpy as np

# define some defaults\n",
plt.rcParams["font.family"] = "DejaVu Serif"
plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams["mathtext.rm"] = "serif"
plt.rcParams["axes.formatter.limits"] = [-4, 4]
plt.rcParams["font.size"] = 12.0

# define a colorset (just the default but written out in full)
color_format = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']

def plot_2d_contour(x, y, data, 
                    ax = None, ax_labels = None,
                    **plot_opts):
  
  """
  Generic command to plot a 2-dimensional field
  """
  if ax is None:
    ax = plt.axes()
  
  if ax_labels is None:
    ax_labels = {"title" : "", "xlabel" : "", "ylabel" : ""}
    
  cs = ax.contour(x, y, data, **plot_opts)
  ax.set_title(ax_labels["title"])
  ax.set_xlabel(ax_labels["xlabel"])
  ax.set_ylabel(ax_labels["ylabel"])
  ax.grid()
  
  return (ax, cs)
    
def plot_2d_contourf(x, y, data, 
                     ax = None, ax_labels = None,
                     **plot_opts):
  
  """
  Generic command to plot a 2-dimensional field
  """
  if ax is None:
    ax = plt.axes()
  
  if ax_labels is None:
    ax_labels = {"title" : "", "xlabel" : "", "ylabel" : ""}
    
  cs = ax.contourf(x, y, data, **plot_opts)
  ax.set_title(ax_labels["title"])
  ax.set_xlabel(ax_labels["xlabel"])
  ax.set_ylabel(ax_labels["ylabel"])
  ax.grid()
  
  return (ax, cs)
