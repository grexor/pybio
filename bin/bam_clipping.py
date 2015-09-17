#!/usr/bin/python

import os
import apa
import pybio
import pysam
import argparse
import gzip
import numpy
import glob
import sys

import matplotlib
matplotlib.use("Agg", warn=False)
import matplotlib.pyplot as plt
import math
from matplotlib import cm as CM
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as mcolors
c = mcolors.ColorConverter().to_rgb

# styling
matplotlib.rcParams['axes.labelsize'] = 10
matplotlib.rcParams['axes.titlesize'] = 10
matplotlib.rcParams['xtick.labelsize'] = 10
matplotlib.rcParams['ytick.labelsize'] = 10
matplotlib.rcParams['legend.fontsize'] = 10
matplotlib.rc('axes',edgecolor='gray')
matplotlib.rcParams['axes.linewidth'] = 0.3
matplotlib.rcParams['legend.frameon'] = 'False'
from matplotlib import gridspec

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-bam', action="store", dest="bam", default=None)
parser.add_argument('-image', action="store", dest="image", default=None)
args = parser.parse_args()

seq_len = 0
seqs = []
all_files = glob.glob(args.bam)

for bam_filename in all_files:
    bam_file = pysam.Samfile(bam_filename)
    count = 0
    for a in bam_file.fetch():
        seq_len = max(seq_len, len(a.seq))
        sc_3 = 0
        sc_5 = 0
        if a.is_reverse:
            last_cigar = a.cigar[0]
            first_cigar = a.cigar[-1]
        else:
            last_cigar = a.cigar[-1]
            first_cigar = a.cigar[0]
        removed_3 = 0
        removed_5 = 0
        if last_cigar[0]==4:
            removed_3 = last_cigar[1]
        if first_cigar[0]==4:
            removed_5 = first_cigar[1]
        if a.is_reverse:
            read_sequence = pybio.sequence.reverse_complement(a.seq)
            mapped_sequence = pybio.sequence.reverse_complement(a.query)
        else:
            read_sequence = a.seq
            mapped_sequence = a.query

        flank_5 = ""
        flank_3 = ""
        if removed_5>0:
            flank_5 = read_sequence[:removed_5]
        if removed_3>0:
            flank_3 = read_sequence[-removed_3:]

        read_reconstructed =  "%s%s%s" % (flank_5,mapped_sequence,flank_3)
        assert(read_sequence==read_reconstructed)
        if removed_3==0:
            assert(read_sequence[removed_5:]==mapped_sequence)
        else:
            assert(read_sequence[removed_5:-removed_3]==mapped_sequence)

        row = [flank_5, mapped_sequence,flank_3]
        seqs.append(row)
        count += 1
        if count%10000==0:
            print "reading %s, %sM sequences read [read_len=%s]" % (bam_filename, count/1e6, seq_len)
        #if count>10000:
        #    break

print len(seqs), count

map_5 = [0]*seq_len
map_3 = [0]*seq_len
map_read = [0]*seq_len

count = 0
for r in seqs:
    count += 1
    read_5 = r[0]
    read_mapped = r[1]
    read_3 = r[2]
    for i in range(0, len(read_5)):
        map_5[i] = map_5[i]+1
    for i in range(len(read_5), len(read_5)+len(read_mapped)):
        map_read[i] = map_read[i]+1
    for i in range(len(read_5)+len(read_mapped), len(read_5)+len(read_mapped)+len(read_3)):
        map_3[i] = map_3[i]+1
    if count%10000==0:
        print "%.3fM reads processed" % (count/1e6)

map_5 = [x/float(count) for x in map_5]
map_read = [x/float(count) for x in map_read]
map_3 = [x/float(count) for x in map_3]

fig = plt.figure(figsize=(20, 4))
gs = gridspec.GridSpec(2, 1, height_ratios=[10, 1])
ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1], sharex=ax1)
#fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(20, 4))
ax1.set_xlim(0, seq_len)
ax1.set_ylim(0, 1)

ax1.set_title("Soft clipping distribution [%.3fM sites]" % (count/1e6))
ax1.set_ylabel("fraction")
plt.xlabel("position in read")

legend = []
for name, data in [("clipped_5", map_5), ("mapped to genome", map_read), ("clipped_3", map_3)]:
    legend.append(name)
    ax1.plot(range(0, seq_len), data, alpha=0.6, linewidth=2)

ax1.legend(legend, loc='upper left')

# heatmaps
import matplotlib.cm as cm
import matplotlib.colors as mcolors

heatmap = ax2.pcolor(numpy.array([map_read]), cmap="BuGn")
ax2.set_yticks([]) # remove y scale ticks for heatmap

plt.savefig(args.image+".png", dpi=150)
plt.close()

f = open(args.image+".tab", "wt")
f.write("count=%s\n" % count)
f.write("map_5=%s\n" % map_5)
f.write("map_read=%s\n" % map_read)
f.write("map_3=%s\n" % map_3)
f.close()
