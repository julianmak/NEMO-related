#!/usr/bin/env python3

# JM: 03 Sep 2018
# plots the time series of the transport of *one* experiment

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from pyCDFTOOLS.process_scalar import txt_to_array, read_vol_transport

# style settings

plt.rcParams["font.family"] = "DejaVu Serif"
plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams["mathtext.rm"] = "serif"
plt.rcParams["image.cmap"] = "RdBu_r" # "*_r" is reverse of standard colour

#--------------------------------------------------------
# define the argument parser
import argparse

parser = argparse.ArgumentParser(description = "Plots the time series of the transport of one experiment")

# fixed arguments
parser.add_argument("data_dir", type = str, 
                    help = "specify data directory")

# optional arguments
parser.add_argument("--str_key", type = str, 
                    help = "specify transport key to search for (has to match string exactly! default = \"ACC_Drake_Passage\")", default = "ACC_Drake_Passage")
parser.add_argument("--xlim", nargs = 2, type = float,
                    help = "specify the time window to plot over if any")
parser.add_argument("--ylim", nargs = 2, type = float,
                    help = "specify the transport range to plot over if any")
parser.add_argument("--l_temp",   
                    help = "plot the mean temperature as well on the same plot", action = "store_true")
parser.add_argument("--tylim", nargs = 2, type = float,
                    help = "specify the temperature range to plot over if any")

# collect arguments
args = parser.parse_args()

#--------------------------------------------------------
# Main plotting commands

time, transport = read_vol_transport(args.data_dir, args.str_key)

if args.xlim is None:
  args.xlim = []
  args.xlim.append(time[0])
  args.xlim.append(time[-1])
  
if args.ylim is None:
  args.ylim = []
  args.ylim.append(min(transport) - 10)
  args.ylim.append(max(transport) + 10)

fig = plt.figure(figsize=(10, 4)) # dpi shouldn't matter for pdf outputs

ax1 = plt.axes()
line1 = plt.plot(time, transport, 'b', label = r"Transport")
plt.xlabel(r"$t$ ($\mathrm{yrs}$)")
plt.ylabel(r"Transport ($\mathrm{Sv}$)")
plt.xlim(args.xlim[0], args.xlim[1])
plt.ylim(args.ylim[0], args.ylim[1])
plt.title(args.str_key)
plt.grid()

lines = line1

if args.l_temp:

  time, mean_temp = txt_to_array(args.data_dir, "sctemtot", query = False)
  
  if args.tylim is None:
    args.tylim = []
    args.tylim.append(min(mean_temp) - 0.2)
    args.tylim.append(max(mean_temp) + 0.2)

  ax2 = ax1.twinx()
  line2 = plt.plot(time, mean_temp, 'r', label = r"$\left\langle\theta\right\rangle$")
  plt.ylim(args.tylim[0], args.tylim[1])
  plt.ylabel(r"$\left\langle\theta\right\rangle$ (${}^\circ\mathrm{C}$)")

  lines += line2

plt.legend(handles = lines, loc=4)

#--------------------------------------------------------
# Saving commands

save_filename = args.str_key + "_transport.pdf"

fig.savefig(save_filename, bbox_inches = "tight")

print("generated %s , exiting..." % save_filename)


