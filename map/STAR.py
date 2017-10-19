import pybio
from os.path import join as pjoin
import os
import sys
import locale

def star(genome, fastq_file, output_folder, name, cpu=1, verbose=True, minlen=0.66):
    """
    Use the `STAR aligner <https://code.google.com/p/rna-star>`_ to map single-end reads to the reference genome.

    :param genome: the genome ID to which to map the reads to
    :param fastq_file: path to fastq file (can be gzipped)
    :param output_folder: results will be stored in this folder
    :param name: name of results
    :param cpu: number of cores to use
    :param verbose: print executed STAR commands
    :param minlen: outFilterMatchNminOverLread, if -1, use default
    """

    version = pybio.genomes.get_latest_version(genome)
    genome_folder = os.path.join(pybio.path.genomes_folder, "%s.assembly.%s.star" % (genome, version))
    command = "pybio.star %s %s %s %s %s %s" % (output_folder, genome_folder, os.path.abspath(fastq_file), name, cpu, minlen)
    if verbose:
        print command
    output, error = pybio.utils.Cmd(command).run()
    return output

def star_pair(genome, file1, file2, output_folder, name, cpu=1, verbose=True):
    """
    Use the `STAR aligner <https://code.google.com/p/rna-star>`_ to map paired-end reads to the reference genome.

    :param genome: the genome ID to which to map the reads to
    :param file1: path to fastq file 1 (can be gzipped)
    :param file2: path to fastq file 2 (can be gzipped)
    :param output_folder: results will be stored in this folder
    :param name: name of results
    :param cpu: number of cores to use
    :param verbose: print executed STAR commands
    """

    version = pybio.genomes.get_latest_version(genome)
    genome_folder = os.path.join(pybio.path.genomes_folder, "%s.assembly.%s.star" % (genome, version))
    command = "pybio.star.pair %s %s %s %s %s %s" % (output_folder, genome_folder, os.path.abspath(file1), os.path.abspath(file2), name, cpu)
    if verbose:
        print command
    output, error = pybio.utils.Cmd(command).run()
    return output
