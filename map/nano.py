import pybio
from os.path import join as pjoin
import os
import sys
import locale

def nano(genome, fastq_file, output_folder, name, cpu=1, verbose=True, minlen=0.66):
    version = pybio.genomes.get_latest_version(genome)
    genome_fasta_file = os.path.join(pybio.path.genomes_folder, "%s.assembly.%s/%s.fasta" % (genome, version, genome))
    command = "pybio.nano %s %s %s %s %s %s" % (output_folder, genome_fasta_file, os.path.abspath(fastq_file), name, cpu, minlen)
    if verbose:
        print(command)
    output, error = pybio.utils.Cmd(command).run()
    return output
