import gzip
import pybio

class Sissrs():
    def __init__(self, filename):
        self.data = []
        f = open(filename, "rt")
        r = f.readline()
        while not r.startswith("Chr	cStart	cEnd	NumTags"):
            r = f.readline()
        r = f.readline()
        r = f.readline()
        while r:
            if r.startswith("==========="):
                r = f.readline()
                continue
            r = r.replace("\r", "").replace("\n", "").split("\t")
            self.data.append(r)
            r = f.readline()
        f.close()

    def bedgraph(self, filename, p_value = None):
        f = open(filename, "wt")
        for r in self.data:
            if len(r)>4 and p_value!=None:
                if float(r[-1])>p_value:
                    r = f.readline()
                    continue
            f.write("%s\t%s\t%s\t%s\n" % (r[0], r[1], r[2], r[3]))
        f.close()
