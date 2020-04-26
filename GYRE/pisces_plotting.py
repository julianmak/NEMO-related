#!/usr/bin/env python3

# JM: 25 Apr 2020
# some generic plotting commands for PISCES

from pyCDFTOOLS.generic_plotting import *
from netCDF4 import Dataset

def plot_bgc_all_latlon(data_dir, filename, var_list = None,
                        kt = -1, kz = 0, 
                        axes = None, axes_cb = None, plot_opts = None):
  """
  Default command to plot all BGC fields
  """
  if plot_opts is None:
    plot_opts = [{}] * len(axes)
    
  if var_list is None:
    var_list = ["DET", "ZOO", "PHY", "NO3", "NH4", "DOM"]
  
  # keep file open to load variables
  data    = Dataset(data_dir + "/" + filename)
  xT      = data.variables["nav_lon"][:, :]
  yT      = data.variables["nav_lat"][:, :]
  
  lines = []
  for i in range(len(axes)):
    
    if plot_opts[i].get("extend") is None:
      plot_opts[i]["extend"] = "both"
    if plot_opts[i].get("cmap") is None:
      plot_opts[i]["cmap"] = cm.get_cmap("gist_ncar")
    
    var = data.variables[var_list[i]][kt, kz, :, :]
    # don't plot the boundary values to appease the contour algorithms
    var_field = var[1:-1, 1:-1]
    ax, cs = plot_2d_contourf(xT[1:-1, 1:-1], yT[1:-1, 1:-1], var_field, 
                              ax = axes[i], ax_labels = None,
                              **plot_opts[i])
    ax.set_aspect("equal")
    
    if plot_opts[i].get("vmin") is None:
      plot_opts[i]["vmin"] = np.min(var_field)
    if plot_opts[i].get("vmax") is None:
      plot_opts[i]["vmax"] = np.max(var_field)
    norm = Normalize(vmin = plot_opts[i]["vmin"], vmax = plot_opts[i]["vmax"])
    colors = plot_opts[i]["cmap"](np.linspace(0, 1, plot_opts[i]["cmap"].N))
    cmap2 = LinearSegmentedColormap.from_list('dummy', colors)
    cb = ColorbarBase(axes_cb[i], cmap = cmap2, norm = norm)
    
    lines.append(cs)
  
  data.close()
  
  return (axes, lines, axes_cb)
