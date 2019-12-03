import pybio
import pysam
from os.path import join as pjoin

def bam_reads(filename):
    print("bam_reads {filename}".format(filename=filename))
    return reduce(lambda x, y: x + y, [ eval('+'.join(l.rstrip('\n').split('\t')[2:]) ) for l in pysam.idxstats(filename) ])

class Bam():

    def __init__(self, filename):
        self.filename = filename

    def get_coverage(self, chr, strand, start, stop):
        raw_count = 0
        command = "samtools view -F 4 -q 30 {filename} {chr}:{start}-{stop}".format(filename = self.filename, chr = chr, start = start, stop = stop)
        output, err = pybio.utils.cmd(command)
        output = output.split("\n")
        for line in output:
            line = line.split("\t")
            if len(line)>3:
                if int(line[3])<=start and int(line[7])<=stop:
                    raw_count += 1
        return raw_count
