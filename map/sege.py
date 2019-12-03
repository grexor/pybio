import pybio
from os.path import join as pjoin
import os
import sys
import locale

def sege(genome, fastq_file, output_folder, name, cpu=1, verbose=True):
    """
    Use the `segemehl.x aligner <http://www.bioinf.uni-leipzig.de/Software/segemehl/>`_ to map single-end reads to the reference genome.

    :param genome: the genome ID to which to map the reads to
    :param fastq_file: path to fastq file (can be gzipped)
    :param output_folder: results will be stored in this folder
    :param name: name of results
    :param cpu: number of cores to use
    :param verbose: print executed STAR commands
    """

    genome_index = os.path.join(pybio.path.genomes_folder, "%s_indices/%s_sege.idx" % (genome, genome))
    genome_fasta = os.path.join(pybio.path.genomes_folder, "%s.assembly/%s.fasta" % (genome, genome))
    command = "pybio.sege %s %s %s %s %s %s" % (output_folder, genome_index, genome_fasta, os.path.abspath(fastq_file), name, cpu)
    if verbose:
        print(command)
    output, error = pybio.utils.Cmd(command).run()
    return output
