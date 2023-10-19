import pybio
from os.path import join as pjoin


# modules
from pybio.data.Fastq import *
from pybio.data.Fasta import *
from pybio.data.Bedgraph import *

def bedgraph_bigwig(filename_bed, filename_bw, genome):
    """
    BigWig can't display + and - strands together, always filter bedGraph before
    """
    filename_genome_chrs = pjoin(pybio.path.folder_genomes, genome, "sequence", "%s.chrs" % genome)
    command = "sort -k1,1 -k2,2n {filename_bed} > {filename_bed}.temp".format(filename_bed=filename_bed) # sort bedGraph
    output, error = pybio.utils.cmd(command)
    command = "bedGraphToBigWig {filename_bed}.temp {filename_genome_chrs} {filename_bw}".format(filename_bed=filename_bed, filename_genome_chrs=filename_genome_chrs, filename_bw=filename_bw)
    output, error = pybio.utils.cmd(command)
    command = "rm {filename_bed}.temp".format(filename_bed=filename_bed) # delete sorted temo bedGraph file
    output, error = pybio.utils.cmd(command)

def fastq_qminmax(filename):
    """
    Returns min and max quality ord value from fastq file.
    """
    qmin = ()
    qmax = 0
    f = pybio.data.Fastq(filename)
    while f.read():
        qmin = min(qmin, ord(f.quality[-1]))
        qmax = max(qmax, ord(f.quality[0]))
    return qmin, qmax

def fasta_check(filename, allowed_chars=["A", "C", "T", "G", "N"]):
    """
    Checks if FASTA file format is valid.
    """
    f = pybio.data.Fasta(filename)
    valid_fasta = True
    while f.read():
        seq_len = 0
        if len(f.sequence)==0:
            return False, "Sequence with ID %s has length 0" % (f.id)
        for allowed in allowed_chars:
            seq_len += f.sequence.upper().count(allowed)
        if seq_len!=len(f.sequence):
            return False, "Sequence with characters other than %s" % str(allowed_chars)
    return True, "File in FASTA format"
