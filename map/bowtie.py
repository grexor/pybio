import pybio
from os.path import join as pjoin
import os
import sys
import locale

def bowtie(genome, fastq_file, output_folder, name, cpu=1, verbose=True):
    """
    Use the `bowtie aligner <http://bowtie-bio.sourceforge.net/manual.shtml>`_ to map single-end reads to the reference genome.

    :param genome: the genome ID to which to map the reads to
    :param fastq_file: path to fastq file (can be gzipped)
    :param output_folder: results will be stored in this folder
    :param name: name of results
    :param cpu: number of cores to use
    :param verbose: print executed STAR commands
    """

    f = pybio.data.Fastq(fastq_file)
    f.read()
    seq_len = len(f.sequence)

    genome_index = os.path.join(pybio.path.genomes_folder, "%s_indices/%s.bowtie" % (genome, genome))
    command = "pybio.bowtie %s %s %s %s %s %s" % (output_folder, genome_index, os.path.abspath(fastq_file), name, cpu, seq_len)
    if verbose:
        print(command)
    output, error = pybio.utils.Cmd(command).run()
    return output

def bowtie2(genome, fastq_file, output_folder, name, cpu=1, verbose=True):
    genome_index = os.path.join(pybio.path.genomes_folder, "%s.assembly.ensembl90/%s" % (genome, genome))
    command = "pybio.bowtie2 %s %s %s %s %s %s" % (output_folder, genome_index, os.path.abspath(fastq_file), name, cpu)
    if verbose:
        print(command)
    output, error = pybio.utils.Cmd(command).run()
    return output
