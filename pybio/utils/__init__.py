import os
import sys
import subprocess
import gzip
import pybio
import math
import numpy as np
import pickle

endings = [".gzip", ".gz", ".bz2"]

def gauss(hw=2, sigma=1):
    """
    Gauss
    """
    hw = 2*hw+1
    r = range(-int(hw/2),int(hw/2)+1)
    return [1 / (sigma * math.sqrt(2*math.pi)) * math.exp(-float(x)**2/(2*sigma**2)) for x in r]

def smooth(data, hw=2):
    """
    Smooth (convolve) with gauss.
    """
    g = gauss(hw)
    return np.convolve(data, g, "same")

def density(x, x_grid, bandwidth=0.2, **kwargs):
    """
    Kernel Density Estimation with Scikit-learn
    """
    kde_skl = KernelDensity(bandwidth=bandwidth, **kwargs)
    kde_skl.fit(x[:, np.newaxis])
    # score_samples() returns the log-likelihood of the samples
    log_pdf = kde_skl.score_samples(x_grid[:, np.newaxis])
    return np.exp(log_pdf)

def merge_intervals(ints):
    """
    | Each interval defined by: (from, to), e.g. (10,100); the interval region is closed, meaning (1,10) includes positions [1,10]
    | Input: list of intervals, e.g. (10,20), (15,25), (100,120)
    | Output: [(10,25), (100,120)]
    """

    def take_next(ints):
        to_return = ints[0]
        del ints[0]
        return to_return

    def add_new(ints, start, stop):
        if len(ints)==0:
            ints.append((start, stop))
            return
        (start_0, stop_0) = ints[-1]
        if start_0 <= start <= stop_0:
            del ints[-1]
            ints.append((start_0, max(stop_0, stop)))
        else:
            ints.append((start, stop))

    ints.sort()
    ints_merged = []
    while len(ints)>0:
        (start, stop) = take_next(ints)
        add_new(ints_merged, start, stop)
    ints_merged.sort()
    return ints_merged

def join_intervals(intervals):
    """
    | Each interval defined by: (from, to, value), e.g. (10,100,1); the interval region is right open, meaning (1,10) includes positions [1,9] = [1, 10)
    | Input: list of intervals, e.g. (100,120,1), (200,300,1), (200,220,2), (215,250,10)
    | Output: merged and spliced intervals, according to value
    | E.g.:[(100, 120, 1), (200, 215, 3), (215, 220, 13), (220, 250, 11), (250, 300, 1)]
    """
    def merge_overlapping(data):
        data.sort()
        starts = {}
        stops = {}
        ranges = set()
        for start, stop, value in data:
            starts[start] = starts.setdefault(start, 0) + value
            stops[stop] = stops.setdefault(stop, 0) + value
            ranges.add(start)
            ranges.add(stop)
        ranges = list(ranges)
        ranges.sort()
        value = 0
        result = []
        for r1, r2 in zip(ranges, ranges[1:]):
            value += starts.get(r1, 0)
            if starts.get(r1, None)!=None:
                if starts.get(r2, None)!=None:
                    result.append((r1, r2, value))
                elif stops.get(r2, None)!=None:
                    result.append((r1, r2, value))
            elif stops.get(r1, None)!=None:
                if stops.get(r2, None)!=None:
                    result.append((r1, r2, value))
                elif starts.get(r2, None)!=None:
                    result.append((r1, r2, value))
            value -= stops.get(r2, 0)
        return result

    intervals.sort()
    intervals.append((None, None, None)) # a trick to process all items in the list

    result = []
    overlaps = []
    overlap_flag = False
    overlap_start, overlap_stop, _ = intervals[0]
    for (start1, stop1, value1), (start2, stop2, value2) in zip(intervals, intervals[1:]):
            overlap = pybio.utils.interval_overlap(overlap_start, overlap_stop, start2, stop2)
            if overlap==0:
                if overlap_flag:
                    overlaps.append((start1, stop1, value1))
                    result = result + merge_overlapping(overlaps)
                    overlap_flag = False
                    overlaps = []
                else:
                    result.append((start1, stop1, value1))
                overlap_start = start2
                overlap_stop = stop2
            else:
                overlaps.append((start1, stop1, value1))
                overlap_start = min(overlap_start, start1)
                overlap_stop = max(overlap_stop, stop1, stop2)
                overlap_flag = True
    result.sort()
    return result

def coverage_to_intervals(coverage):
    """
    | input: coverage={10:"5", 12:"5", 13:"o", 14:"o", 25:'3', 26:'3',27:'5'}
    | output: [(10, 10, '5'), (12, 12, '5'), (13, 14, 'o'), (25, 26, '3'), (27, 27, '5')]
    | intervals are inclusive left and right
    """
    if len(coverage)==0:
        return [] # no coverage? no intervals
    L = []
    i = None
    val_index = min(coverage.keys())
    val = coverage[val_index]
    for i in range(min(coverage.keys())+1, max(coverage.keys())+1):
        new_val = coverage.get(i, None)
        if val!=None and new_val!=None and val!=new_val:
            L.append((val_index, i-1, val))
            val = new_val
            val_index = i
        if val!=None and new_val==None:
            L.append((val_index, i-1, val))
            val = new_val
            val_index = i
        if val==None and new_val!=None:
            val = new_val
            val_index = i
    if i!=None:
        L.append((val_index, i, val))
    return L

def interval_overlap(start1, stop1, start2, stop2):
    """
    | Return how much two intervals overlap:
    | (start1, stop1) vs (start2, stop2)
    | Return 0 if they don't overlap
    """
    if stop1 < start2 or stop2 < start1:
        return 0
    else:
        return max(0, min(stop1, stop2) - max(start1, start2)) + 1

class Cmd():
    """
    Interface for running commands from python.
    """

    def __init__(self, command):
        self.command = command
        self.returncode = None
        self.process = subprocess.Popen(['/bin/sh', '-cl', command], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.pid = self.process.pid

    def run(self):
        output, error = self.process.communicate()
        self.returncode = self.process.returncode
        return output.decode("utf-8"), error.decode("utf-8")

def process_isrunning(pid, os="linux"):
    """
    Check if the process with pid is running.
    """
    if os=="linux":
        output, error = Cmd("ps -p %s" % pid).run()
        output = output.split("\n")
        return len(output)>=3

    if os=="windows":
        from win32com.client import GetObject
        WMI = GetObject('winmgmts:')
        processes = WMI.InstancesOf('Win32_Process')
        plist = [process.Properties_('ProcessID').Value for process in processes]
        return pid in plist

def decompress(source, dest=None):
    """
    Decompress gzip or bzip2 files. The original compressed file is left untouched.

    :param dest: the uncompressed output is written to this file. If omitted, the uncompressed content is written to source name but without ".gz" or ".bz2". The original compressed file is left untouched.
    """
    if dest==None:
        for end in endings:
            if source.endswith(end):
                dest = source[:-len(end)]
                break

    if source.endswith(".gzip") or source.endswith(".gz"):
        command = "gunzip -c %s > %s" % (source, dest)
        out, err = cmd(command)
        return dest

    if source.endswith(".bz2"):
        command = "bunzip2 -c %s > %s" % (source, dest)
        out, err = cmd(command)
        return dest

    return source # no decompression

def gzip(source):
    """
    Compress (gzip) input.
    """

    command = "gzip -f %s" % (source)
    out, err = cmd(command)
    return source+".gz"

def split_string(seq, length):
    """
    Splitting strings (used by FASTA module).
    """
    return [seq[i:i+length] for i in range(0, len(seq), length)]

def is_sorted(l):
    return all(l[i] <= l[i+1] for i in range(len(l)-1))

# from Orange bioinformatics stats library
# http://orange.biolab.si
def FDR(p_values, dependent=False, m=None, ordered=False):

    if not ordered:
        ordered = is_sorted(p_values)

    if not ordered:
        joined = [ (v,i) for i,v in enumerate(p_values) ]
        joined.sort()
        p_values = [ p[0] for p in joined ]
        indices = [ p[1] for p in joined ]

    if not m:
        m = len(p_values)
    if m <= 0 or not p_values:
        return []

    if dependent: # correct q for dependent tests
        k = c[m-1] if m <= len(c) else math.log(m) + 0.57721566490153286060651209008240243104215933593992
        m = m * k

    tmp_fdrs = [p*m/(i+1.0) for (i, p) in enumerate(p_values)]
    fdrs = []
    cmin = tmp_fdrs[-1]
    for f in reversed(tmp_fdrs):
        cmin = min(f, cmin)
        fdrs.append( cmin)
    fdrs.reverse()

    if not ordered:
        new = [ None ] * len(fdrs)
        for v,i in zip(fdrs, indices):
            new[i] = v
        fdrs = new

    return fdrs

# replaces column cname in file fname with corrected FDR p-values
def FDR_tab(fname, cname):
    L = []
    f = open(fname)
    header = f.readline().replace("\n", "").replace("\r", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\n", "").replace("\r", "").split("\t")
        data = dict(zip(header, r))
        L.append(float(data[cname]))
        r = f.readline()
    f.close()

    os.rename(fname, fname+"fdr_temp")

    L2 = pybio.utils.FDR(L)
    f1 = open(fname+"fdr_temp")
    f2 = open(fname, "wt")
    header = f1.readline().replace("\n", "").replace("\r", "").split("\t")
    f2.write("\t".join(header)+"\n")
    r1 = f1.readline()
    while r1:
        r1 = r1.replace("\n", "").replace("\r", "").split("\t")
        data = dict(zip(header, r))
        r1[header.index(cname)] = L2.pop(0)
        f2.write("\t".join(str(x) for x in r1) + "\n")
        r1 = f1.readline()
    f1.close()
    f2.close()
    assert(len(L2)==0) # we used all corrected value

    #os.remove(fname+"fdr_temp")

def save_pickle(obj, fname):
    pickle.dump(obj, open(fname, "wb"), -1)

def load_pickle(fname):
    return pickle.load(open(fname, "rb"))

# is a specific software installed on the system? # uses which to check
def is_tool(name):
    from shutil import which
    return which(name) is not None