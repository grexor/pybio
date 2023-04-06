import pybio
import numpy

reverse_code = { \
                'A': 'T', 'T': 'A', 'U': 'A', 'G': 'C', 'C': 'G', 'R': 'Y', 'Y': 'R', 'K': 'M', 'M': 'K', 'S': 'S', 'W': 'W', 'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B', 'N': 'N', \
                'a': 't', 't': 'a', 'u': 'a', 'g': 'c', 'c': 'g', 'r': 'y', 'y': 'r', 'k': 'm', 'm': 'k', 's': 's', 'w': 'w', 'b': 'v', 'd': 'h', 'h': 'd', 'v': 'b', 'n': 'n' \
                }
reverse_code_rysw = {'R' : 'Y', 'Y' : 'R', 'S' : 'S', 'W' : 'W'}
expand_code = { "R": ["A", "G"], "Y" : ["C", "T"], "S" : ["G", "C"], "W" : ["A", "T"], "A":"A", "T":"T", "C":"C", "G":"G" }

def reverse_complement(s):
    """
    :param s: sequence to reverse complement

    Reverse complement sequence
    """
    reverse_str = lambda s: ''.join([s[i] for i in range(len(s)-1, -1, -1)])
    rs = reverse_str(s)
    new_sequence = ""
    for c in rs:
        rc = reverse_code.get(c.upper(), c.upper())
        if c.islower():
            rc = rc.lower()
        new_sequence += rc
    return new_sequence

def turn_char(c1):
    if c1=="C":
        return "T"
    if c1=="T":
        return "C"
    if c1=="A":
        return "G"
    if c1=="G":
        return "A"

def all_indices(string, sub,offset=0):
    listindex = []
    i = string.find(sub, offset)
    while i >= 0:
        listindex.append(i)
        i = string.find(sub, i + 1)
    return listindex

def overlap(item1, item2):
    start1 = item1[0]
    start2 = item2[0]
    stop1 = item1[1]
    stop2 = item2[1]
    if start2>=start1+1 and start2<=stop1+1:
        return ((start1, stop2), True)
    else:
        return (item1, False)

def expand(motif):
    """
    Expands motif with beginning and ending nucleotide R,Y,S,W to a list of 4 motifs with alphabet A,T,C,G

    Example: RAAY -> ['AAAC', 'AAAT', 'GAAC', 'GAAT']
    """
    start_n = motif[0]
    end_n = motif[-1]
    if expand_code.get(start_n, None)==None or expand_code.get(end_n, None)==None:
        return [motif]
    motifs = set()
    for x in expand_code[start_n]:
        for y in expand_code[end_n]:
            motifs.add(x+motif[1:-1]+y)
    return list(motifs)

def search(input_string, motif_list, strict=False):
    vector = [0] * len(input_string)
    if type(motif_list) is not list:
        motif_list = motif_list.split("_")
    for m in motif_list:
        expanded = expand(m)
        pool = []
        for em in expanded:
            z = all_indices(input_string, em)
            if len(z)==0 and strict:
                return [], [0] * len(input_string)
            for i in z:
                pool.append((i, i+len(em)))
            for (start, stop) in pool:
                vector[start:stop] = [1] * (abs(start-stop))
    return pool, vector

def filter(vector, hw=30, hwt=1):
    vector = numpy.convolve(vector, [1]*(hw*2+1), "same")
    vector = [1 if x>=hwt else 0 for x in vector]
    return vector

def convolve(vector, hw):
    vector = numpy.convolve(vector, [1]*(hw*2+1), "same")
    return vector

def extend(vector, window_size):
    new_vector = [sum(vector[max(0, i-window_size/2):i+window_size/2+1]) for i in range(0, len(vector))]
    return new_vector

def draw(sequences, fname):

    import matplotlib
    matplotlib.use("Agg", warn=False)
    import matplotlib.pyplot as plt
    import math
    import gzip
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

    def compute_f(sequences, max_l):
        dist = {}
        for i in range(0, max_l):
            dist[i] = {"A":0, "C":0, "U":0, "G":0, "N":0}
        for seq in sequences:
            for i in range(0, len(seq)):
                c = {"A":"A", "C":"C", "T":"U", "G":"G", "N":"N"}[seq[i]]
                dist[i][c] = dist[i].get(c, 0) + 1
        v_a = []
        v_u = []
        v_c = []
        v_g = []
        for i in range(max_l):
            v_sum = dist[i]["A"]+dist[i]["U"]+dist[i]["C"]+dist[i]["G"]
            v_a.append(float(dist[i]["A"])/v_sum*100)
            v_u.append(float(dist[i]["U"])/v_sum*100)
            v_c.append(float(dist[i]["C"])/v_sum*100)
            v_g.append(float(dist[i]["G"])/v_sum*100)
        return v_a, v_u, v_c, v_g, len(sequences)

    max_l = 0
    for s in sequences:
        max_l = max(max_l, len(s))
    fig, ax = plt.subplots(1, 1, figsize=(12, 3))
    ax.set_ylim(0, 100)
    ax.set_xlim(0, max_l)

    legend = []
    v_a, v_u, v_c, v_g, count = compute_f(sequences, max_l)
    ax.set_title("nucleotide composition [%s sites]" % format(count, ','))
    ax.set_ylabel("nucleotide percentage")
    plt.xlabel("nucleotide position")
    for name, data, color in [("A", v_a, "red"), ("U", v_u, "blue"), ("C", v_c, "green"), ("G", v_g, "orange")]:
        legend.append(name)
        ax.plot(range(0, max_l), data, alpha=0.7, linewidth=3, color=color)
    ax.set_xticks(range(0, max_l, 10))
    ax.set_xticklabels(range(0, max_l, 10))
    ax.legend(legend, loc='upper left')
    ax.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig(fname)
    plt.close()
