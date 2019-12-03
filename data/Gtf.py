import pybio
import os
import sys
import pickle

cache_data = {}
def cache_string(string):
    if cache_data.get(string, None)==None:
        cache_data[string] = string
        return string
    else:
        return cache_data[string]

class Gtf():

    def __init__(self, filename):
        self.genes = {}
        self.filename = filename
        f = pybio.data.TabReader(filename)
        while f.readline():
            chr = f.r[0]
            gene_type = f.r[2]
            start = int(f.r[3])
            stop = int(f.r[4])
            strand = f.r[6]
            attrs = {}
            temp = f.r[-1].split(";")
            for att in temp:
                att = att.replace("\"", "")
                att = att.lstrip(" ")
                att = att.split(" ")
                attrs[att[0]] = " ".join(att[1:])
            if attrs.get("gene_id", None)==None:
                continue
            gene = self.genes.get(attrs["gene_id"], pybio.data.Gene(attrs["gene_id"], chr, strand, attrs=attrs))
            feature = pybio.data.GeneFeature(start, stop, gene_type, gene)
            gene.add_feature(feature)
            self.genes[gene.id] = gene

    def get_genes(self, chr, pos):
        bin = pos/self.bin_size
        candidate_genes = self.pindex.get(chr, {}).get(bin, [])
        position_genes = set()
        for gene_id in candidate_genes:
            for feature in self.genes[gene_id].features:
                if feature.type!="exon":
                    continue
                if feature.start<=pos<=feature.stop:
                    position_genes.add(gene_id)
        return position_genes

    def write_gff3(self, filename):
        f = open(filename, "wt")
        for gene_id, gene in self.genes.iteritems():
            row = [gene.chr, "ap", "gene", gene.start, gene.stop, "", gene.strand, ".", "ID=%s;Name=%s" % (gene_id, gene_id)] # gene
            f.write("\t".join(str(x) for x in row) + "\n")
            row = [gene.chr, "ap", "mRNA", gene.start, gene.stop, "", gene.strand, ".", "ID=%s.t1;Parent=%s" % (gene_id, gene_id)] # mRNA
            f.write("\t".join(str(x) for x in row) + "\n")
            for exon_index, feature in enumerate(gene.features):
                if feature.type not in ["exon", "CDS"]:
                    continue
                row = [gene.chr, "ap", "CDS", feature.start, feature.stop, "", gene.strand, ".", "ID=%s.t1.cds;Parent=%s.t1" % (gene_id, gene_id)] # mRNA
                f.write("\t".join(str(x) for x in row) + "\n")
                row = [gene.chr, "ap", "exon", feature.start, feature.stop, "", gene.strand, ".", "ID=%s.t1.exon%s;Parent=%s.t1" % (gene_id, exon_index+1, gene_id)] # mRNA
                f.write("\t".join(str(x) for x in row) + "\n")
            f.write("\n")
        f.close()
