#!/usr/bin/env python3
#
# JM, 28 Dec 2017
# some default commands for plotting ORCA data
# feel free to adapt as appropriate

#### load modules

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import netCDF4

from matplotlib_colorbar import *

# needed for adding colorbar in certain transforms
from mpl_toolkits.axes_grid1 import make_axes_locatable

# cartopy puts things on maps by transforming data co-ordinates etc
# though it seems it only works for certain projections with raw ORCA2 data
import cartopy
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

import iris
import iris.analysis.cartography

# style settings

plt.rcParams["font.family"] = "DejaVu Serif"
plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams["mathtext.rm"] = "serif"
plt.rcParams["image.cmap"] = "RdBu_r" # "*_r" is reverse of standard colour

################################################################################
#### some specific plotting commands adapted for ORCA

def plot_eke(data_dir, filename, var_name, kt = 0, lprint = False, ax = None,
             target_proj = None, **misc_format):
  """
  plot the depth-integrated eddy energy that comes out from the GEOMETRIC scheme
  
  Inputs:
    data_dir     = data directory
    filename     = file name
    var_name     = variable name in filename
    kt           = 0 by default
    lprint       = false by default, print out filename if need be
    misc_format  = add in formats for plotting commands
  
  Returns:
    ax           = diagram axes
    mesh         = data mesh
    target_proj  = projection used
  """
  # load files
  eE_raw, lonT, latT, _ = load_var_2d(data_dir, filename, var_name)
  #e3T, _, _, _ = load_var_3d(data_dir, filename, "e3t")

  #tmask = load_mask(data_dir, "../mesh_mask.nc", "tmask")
  
  # depth averaging
  #eE_zavg = eE_raw / np.maximum(np.sum(e3T * tmask, axis = 0), 1e-16)
  eE_zavg = eE_raw
  
  # load in iris and defined some generic things
  iris.FUTURE.netcdf_promote = True

  pcarree = ccrs.PlateCarree()
  
  if target_proj is None:
    print("no projection specified, using pcarree")
    target_proj = pcarree

  lat = iris.coords.AuxCoord(latT, standard_name = "latitude", units = "degrees")
  lon = iris.coords.AuxCoord(lonT, standard_name = "longitude", units = "degrees")
  data_cube = iris.cube.Cube(eE_zavg, 
                             long_name = "geom_eke_depth_avg", 
                             units = "m2 s-2",
                             aux_coords_and_dims = [(lat, (0, 1)), (lon, (0,1))])
  data_proj, extent = iris.analysis.cartography.project(data_cube[:, :], pcarree, nx = 600, ny = 300)
  x = data_proj.coord('projection_x_coordinate').points
  y = data_proj.coord('projection_y_coordinate').points
  plot_data = data_proj.data
  
  # touch up the data and set levels
  # "mask" the data by setting the land point data to be the minimum of the cmin
  cmin = plot_data.min()
  cmax = np.ceil(plot_data.max())
  plot_data[(plot_data == 0)] = cmin
  cmin = np.floor(cmin)

  # multiply by rho_0 = 1024 here (or change it if you want)
  ax, mesh = cartopy_contourf(x, y, np.log10(1024 * plot_data), ax, target_proj, **misc_format)

  if type(target_proj).__name__ in ["PlateCarree", "Mercator"]:
    gl = ax.gridlines(crs=ccrs.PlateCarree(),
                  draw_labels = True, linewidth = 1, linestyle = '--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
  else:
    ax.gridlines(linewidth = 1, linestyle = '--')
  
  return (ax, mesh, target_proj)
  
#-------------------------------------------------------------------------------

def plot_kgm(data_dir, filename, var_name, kt = 0, lprint = False, ax = None,
             target_proj = None, **misc_format):
  """
  plot the depth-averaged kap_gm that comes out from the GEOMETRIC scheme
  
  Inputs:
    data_dir     = data directory
    filename     = file name
    var_name     = variable name in filename
    kt           = 0 by default
    lprint       = false by default, print out filename if need be
    misc_format  = add in formats for plotting commands
  
  Returns:
    ax           = diagram axes
    mesh         = data mesh
    target_proj  = projection used
  """
  # load files
  kgm_raw, lonW, latW, _ = load_var_3d(data_dir, filename, var_name)
  #e3W, _, _, _ = load_var_3d(data_dir, filename, "e3w")
  #e3T, _, _, _ = load_var_3d(data_dir, filename.replace("_W", "_T"), "e3t")

  #tmask = load_mask(data_dir, "../mesh_mask.nc", "tmask")
  
  # depth averaging
  #kgm_zavg = np.sum(kgm_raw * e3W * tmask, axis = 0) / np.maximum(np.sum(e3T * tmask, axis = 0), 1e-16)
  kgm_zavg = kgm_raw[0, :, :]
  
  # load in iris and defined some generic things
  iris.FUTURE.netcdf_promote = True

  pcarree = ccrs.PlateCarree()
  
  if target_proj is None:
    print("no projection specified, using pcarree")
    target_proj = pcarree

  lat = iris.coords.AuxCoord(latW, standard_name = "latitude", units = "degrees")
  lon = iris.coords.AuxCoord(lonW, standard_name = "longitude", units = "degrees")
  data_cube = iris.cube.Cube(kgm_zavg, 
                             long_name = "kgm_depth_avg", 
                             units = "m s-2",
                             aux_coords_and_dims = [(lat, (0, 1)), (lon, (0,1))])
  data_proj, extent = iris.analysis.cartography.project(data_cube[:, :], pcarree, nx = 600, ny = 300)
  x = data_proj.coord('projection_x_coordinate').points
  y = data_proj.coord('projection_y_coordinate').points
  plot_data = data_proj.data
  
  # touch up the data and set levels
  # "mask" the data by setting the land point data to be the minimum of the cmin
  cmin = plot_data.min()
  cmax = np.ceil(plot_data.max())
  plot_data[(plot_data == 0)] = cmin
  cmin = np.floor(cmin)

  ax, mesh = cartopy_contourf(x, y, np.log10(plot_data), ax, target_proj, **misc_format)

  if type(target_proj).__name__ in ["PlateCarree", "Mercator"]:
    gl = ax.gridlines(crs=ccrs.PlateCarree(),
                  draw_labels = True, linewidth = 1, linestyle = '--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
  else:
    ax.gridlines(linewidth = 1, linestyle = '--')
  
  return (ax, mesh, target_proj)

################################################################################
################################################################################
#### some general plotting / loading commands

def raw_plot_tricontourf(x_coord, y_coord, plot_data, ax = None, levels = 20):
  """
  plots the raw data with tricontourf
  
  Returns:
    ax, mesh for editing/formatting purposes
  """
  if ax is None:
    ax = plt.axes(projection = proj)
    
  mesh = plt.tricontourf(x_coord.flatten(), y_coord.flatten(), plot_data.flatten(), levels)
    
  # magic command to get rid of white lines in contourf outputs
  for c in mesh.collections:
      c.set_edgecolor("face")

  return (ax, mesh)
  
#-------------------------------------------------------------------------------

def carree_tricontourf(x_coord, y_coord, plot_data, ax = None, levels = 20):
  """
  plots the data with tricontourf on the PlateCarree projection
  
  Returns:
    ax, mesh for editing/formatting purposes
  """
  if ax is None:
    ax = plt.axes(projection = ccrs.PlateCarree())
    
  mesh = ax.tricontourf(x_coord.flatten(), y_coord.flatten(), plot_data.flatten(), 
                        levels, transform = ccrs.PlateCarree())
  ax.set_global()
  ax.add_feature(cartopy.feature.LAND, zorder = 10, edgecolor = 'k')

  # magic command to get rid of white lines in contourf outputs
  for c in mesh.collections:
    c.set_edgecolor("face")
    
  return (ax, mesh)
  
#-------------------------------------------------------------------------------
  
def cartopy_pcolormesh(x_coord, y_coord, plot_data, ax = None,
                       proj = ccrs.PlateCarree(),
                       **kwargs):
  """
  plots the data with pcolormesh on a specified projection
  
  Returns:
    ax, mesh for editing/formatting purposes
  """
  if ax is None:
    ax = plt.axes(projection = proj)
    
  mesh = ax.pcolormesh(x_coord, y_coord, plot_data,
                       transform = ccrs.PlateCarree(),
                       **kwargs)
  ax.set_global()
  ax.add_feature(cartopy.feature.LAND, zorder = 10, edgecolor = 'k')

  # magic command to get rid of some box lines in output in pcolormesh
  mesh.set_edgecolor('face')
  
  return (ax, mesh)
  
#-------------------------------------------------------------------------------
  
def cartopy_contourf(x_coord, y_coord, plot_data, ax = None,
                     proj = ccrs.PlateCarree(),
                     **kwargs):
  """
  plots the data with contourf on the specified projection
  NOTE: needs an interpolation with iris first to get it on a regular grid
  
  Returns:
    ax, mesh for editing/formatting purposes
  """
  if ax is None:
    ax = plt.axes(projection = proj)
    
  mesh = ax.contourf(x_coord, y_coord, plot_data,
                     transform = ccrs.PlateCarree(),
                     **kwargs)
  ax.set_global()
  ax.add_feature(cartopy.feature.LAND, zorder = 10, edgecolor = 'k')
  
  # magic command to get rid of white lines in contourf outputs
  for c in mesh.collections:
    c.set_edgecolor("face")
    
  return (ax, mesh)

#-------------------------------------------------------------------------------

def load_var_2d(data_dir, filename, var_name, kt = 0, lprint = False):
  """
  load a 2d variable, defaults to loading the first time slice
  
  Inputs:
    data_dir     = data directory
    filename     = file name
    var_name     = variable name in filename
    kt           = 0 by default
    lprint       = false by default, print out filename if need be
  
  Returns:
    var          = variable data array
    lon          = lon grid
    lat          = lat grid
    time_vec     = time stamp according to file
  """
  files_netcdf4 = netCDF4.Dataset(data_dir + filename)
  if lprint:
    print(file_netcdf4)
    
  time_vec = files_netcdf4.variables["time_instant"][:]
  time_vec = (time_vec - time_vec[0]) / (3600 * 24) # time vector in days
  var = files_netcdf4.variables[var_name][kt, :, :]
  lat = files_netcdf4.variables["nav_lat"][:, :]
  lon = files_netcdf4.variables["nav_lon"][:, :]
  files_netcdf4.close()
  
  return (var, lon, lat, time_vec)

#-------------------------------------------------------------------------------

def load_var_3d(data_dir, filename, var_name, kt = 0, lprint = False):
  """
  load a 3d variable, defaults to loading the first time slice
  
  Inputs:
    data_dir     = data directory
    filename     = file name
    var_name     = variable name in filename
    kt           = 0 by default
    lprint       = false by default, print out filename if need be
  
  Returns:
    var          = variable data array
    lon          = lon grid
    lat          = lat grid
    time_vec     = time stamp according to file
  """
  files_netcdf4 = netCDF4.Dataset(data_dir + filename)
  if lprint:
    print(file_netcdf4)
    
  time_vec = files_netcdf4.variables["time_instant"][:]
  time_vec = (time_vec - time_vec[0]) / (3600 * 24) # time vector in days
  var = files_netcdf4.variables[var_name][kt, :, :, :]
  lat = files_netcdf4.variables["nav_lat"][:, :]
  lon = files_netcdf4.variables["nav_lon"][:, :]
  files_netcdf4.close()
  
  return (var, lon, lat, time_vec)
  
#-------------------------------------------------------------------------------

def load_mask(data_dir, filename, var_name, lprint = False):
  """
  load the 3d mask from the mesh_mask.nc file
  
  Inputs:
    data_dir     = data directory
    filename     = file name
    var_name     = variable name in filename
    kt           = 0 by default
    lprint       = false by default, print out filename if need be
  
  Returns:
    var          = variable data array
  """
  files_netcdf4 = netCDF4.Dataset(data_dir + filename)
  if lprint:
    print(file_netcdf4)
    
  var = files_netcdf4.variables[var_name][0, :, : ,:]
  
  files_netcdf4.close()
  
  return var
