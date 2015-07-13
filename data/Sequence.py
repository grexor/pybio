bc = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N', 'a': 't', 'c': 'g', 't': 'a', 'g': 'c', 'n': 'n'}

def rev_comp(seq):
    seq = [bc[s] for s in seq]
    return "".join(seq[::-1])