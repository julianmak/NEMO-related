#!/usr/bin/env python3

# JM: 03 Sep 2018
# plots the time series of the transport of *a batch of* experiments

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
                    help = "specify control data directory")

# optional arguments
parser.add_argument("--str_key", type = str, 
                    help = "specify transport key to search for (has to match string exactly! default = \"ACC_Drake_Passage\")", default = "ACC_Drake_Passage")
parser.add_argument("--exp_key", nargs = "+", type = str,
                    help = "specify the folder names of other experiments to plot (e.g. GEOM_x200 to replace GEOM_x100)")
parser.add_argument("--xlim", nargs = 2, type = float,
                    help = "specify the time window to plot over if any")
parser.add_argument("--ylim", nargs = 2, type = float,
                    help = "specify the transport range to plot over if any")

# collect arguments
args = parser.parse_args()

#--------------------------------------------------------
# Main plotting commands

time, transport = read_vol_transport(args.data_dir, args.str_key)

match = "GEOM"
for key in args.data_dir.split("/"):
  if match in key:
    exp_name = key.replace("EXP_", "")
    print(" ")
    print("working in %s" % args.data_dir)
    print(" ")
    print("control experiment name grabbed = %s ..." % exp_name)
    print(" ")

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
lines = ax1.plot(time, transport, label = exp_name)

if args.exp_key is not None:
  for i in range(len(args.exp_key)):
    print("exchanging %s for %s ..." % (exp_name, args.exp_key[i]))
    print(" ")
    time, transport = read_vol_transport(args.data_dir.replace(exp_name, args.exp_key[i]), args.str_key)
    line = plt.plot(time, transport, label = args.exp_key[i])
    lines += line
    
ax1.set_xlabel(r"$t$ ($\mathrm{yrs}$)")
ax1.set_ylabel(r"Transport ($\mathrm{Sv}$)")
ax1.set_xlim(args.xlim[0], args.xlim[1])
ax1.set_ylim(args.ylim[0], args.ylim[1])
ax1.set_title(args.str_key)
ax1.grid()
ax1.legend(handles = lines, loc = 4, ncol = 2, fontsize = "x-small")

#--------------------------------------------------------
# Saving commands

save_filename = args.str_key + "_transport_batch.pdf"

fig.savefig(save_filename, bbox_inches = "tight")

print("generated %s , exiting..." % save_filename)


