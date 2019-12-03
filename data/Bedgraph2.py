import gzip
import bisect
import pybio
import os
import json

def mix_sort(x):
    if x.isdigit():
        return int(x)
    else:
        return str(x)

class Bedgraph2():
    """
    Class for reading bedgraph files with region-type storage and search.
    """

    def is_int(self, val):
        if type(val) not in [float, int]:
            return False
        return val == int(val)

    def __init__(self, filename=None, genome="ensembl", fixed_cDNA=False, fast=False):
        self.clear()
        if filename!=None:
            self.load(filename, genome=genome, fast=fast, fixed_cDNA=fixed_cDNA)

    def clear(self):
        self.raw = {} # raw counts
        self.rawx = {} # raw counts
        self.filename = ""

    def load(self, filename, track_id=None, meta=None, min_cpm=0, min_raw=0, fixed_cDNA=False, compute_cpm=True, genome="ensembl", fast=False, force_strand=None):
        """
        Load Bedgraph file (can also be gzipped). A Bedgraph object load method can be called multiple times on various files, the content is added up.
        """
        print("loading : {filename}".format(filename=filename))
        self.filename = filename

        if track_id==None:
            track_id = filename

        temp_raw = {}
        temp_rawx = {}
        raw_sum = 0
        if filename.endswith(".gz"):
            f = gzip.open(filename, "rb")
        else:
            f = open(filename)
        r = f.readline()
        not_converted = 0
        while r:
            if r.startswith("track"):
                r = f.readline()
                continue
            if r.startswith("#"):
                r = f.readline()
                continue
            r = r.rstrip("\r").rstrip("\n").split("\t")
            if r==[""]:
                r = f.readline()
                continue
            chr = r[0]
            if genome!="ensembl":
                chr = pybio.genomes.chr_uscs_ensembl.get(genome, {}).get(chr, None) # always convert ucsc to ensembl
                if chr==None: # the mapping of chromosome from ucsc to ensembl didnt success
                    not_converted += 1
                    r = f.readline()
                    continue
            raw = float(r[3])
            # if reading integer numbers, make them integer
            if self.is_int(raw):
                raw = int(raw)
            strand = "+" if raw>=0 else "-"
            if fixed_cDNA!=False:
                raw = fixed_cDNA
            if force_strand!=None:
                assert(force_strand in ["+", "-"])
                strand = force_strand

            temp_raw.setdefault(chr, {}).setdefault(strand, [])
            temp_raw[chr][strand].append((int(r[1]), int(r[2])-1, abs(raw)))
            temp_rawx.setdefault(chr, {}).setdefault(strand, [])
            r = f.readline()
        f.close()
        self.genome = "ensembl" # we converted to ensembl

        if not_converted>0:
            print("{not_converted} positions were skipped; unable to convert chromosome name to ensembl format".format(not_converted=not_converted))

        self.raw = temp_raw
        self.rawx = temp_rawx
        # sort things
        for chr, strand_data in self.raw.items():
            for strand, L in strand_data.items():
                self.raw[chr][strand].sort()
                self.rawx[chr][strand] = [start for (start, stop, val) in self.raw[chr][strand]]
        return

    def get_value(self, chr, strand, pos, db="raw"):
        """
        Get value at chromosome (chr), strand, position.
        """
        return self.get_vector(chr, strand, pos, 0, 1)[0]

    def get_vector(self, chr, strand, pos, start, stop, db="raw"):
        # returns vector of values around pos
        # considers strand but does not reverse the vector
        # example: strand=-, start=-10, stop=20, would return a vector of pos-20 .. pos+10
        # example: strand=+, start=-10, stop=20, would return a vector of pos-10 .. pos+20
        if strand=="-":
            start, stop = -start, -stop
        start, stop = min(start, stop), max(start, stop)
        start = pos+start
        stop = pos+stop
        R = [0] * (stop-start+1)
        i = bisect.bisect(self.rawx.get(chr, {}).get(strand, []), start)
        i -= 1
        x = stop
        len_L = len(self.rawx.get(chr, {}).get(strand, []))
        if len_L==0:
            return R
        while x<=stop and i<len_L:
            x, y, v = self.raw.get(chr, {}).get(strand, [])[max(0, i)]
            r_start = max(start, x)
            r_stop = min(stop, y)
            insert = [v] * (r_stop-r_start+1)
            R[r_start-start:r_start-start+len(insert)] = insert
            i+=1
        return R
