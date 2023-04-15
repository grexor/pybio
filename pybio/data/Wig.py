import gzip
import pybio

class Wig():
    def __init__(self, filename=None, min_value=None, only_positions=None):
        self.trackname = {}
        self.data = {}
        if not filename==None:
            self.load(filename, min_value=min_value, only_positions = only_positions)

    def load(self, filename, min_value=0, only_positions=None):
        if filename.endswith(".gz"):
            f = gzip.open(filename, "rb")
        else:
            f = open(filename, "rb")
        r = f.readline()
        r = f.readline()
        chr = None
        while r:
            r = r.rstrip("\r").rstrip("\n").split("\t")
            if r[0].startswith("variableStep"):
                r = r[0].split(" ")
                if chr!=None:
                    self.data[chr] = L
                chr = r[1].split("=")[1]
                span = int(r[2].split("=")[1])
                self.data[chr] = {'+':[], '-':[]}
                L = self.data[chr]
                r = f.readline()
                continue
            chr = r[0]
            start = int(r[0])
            value = int(r[1])
            strand = "+" if value>=0 else "-"
            L[strand].append((span, start, value))
            r = f.readline()
        f.close()
        for chr in self.data.keys():
            self.data[chr]["+"].sort()
            self.data[chr]["-"].sort()

    def chromosomes(self):
        return self.data.keys()

    def value(self, chr, strand, pos):
        L = self.data.get(chr, {}).get(strand, [])
        for (span, start, value) in L:
            if start<=pos<=start+span-1:
                return value
        return 0

    def region(self, chr, strand, pos_from, pos_to):
        L = self.data.get(chr, {}).get(strand, [])
        count = 0
        for (span, start, value) in L:
            if start>pos_to:
                break
            if start+span-1 < pos_from:
                continue
            count += pybio.utils.interval_overlap(start, start+span-1, pos_from, pos_to) * value
        return count
