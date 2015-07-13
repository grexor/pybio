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
    reverse_str = lambda s: ''.join([s[i] for i in xrange(len(s)-1, -1, -1)])
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
        return motif
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