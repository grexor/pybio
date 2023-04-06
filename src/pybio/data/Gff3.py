import pybio
import os
import sys

class Gff3():

    def __init__(self, filename):
        fasta_part = 0
        mRNA_part = 0
        f = open(filename, "rt")
        r = f.readline()
        self.genes = {}
        self.mRNA_genes = {}
        l = 1
        while r:
            if r.startswith("##FASTA") or fasta_part==1:
                r = f.readline()
                fasta_part = 1
                l+=1
                continue
            if r.startswith("#"):
                r = f.readline()
                l+=1
                continue
            r = r.rstrip("\r\n").split("\t")
            if len(r)==1:
                r = f.readline()
                continue
            if r[0]=="##gff-version   3":
                fasta_part = 0
                r = f.readline()
                l+=1
                continue
            seqid = r[0]
            source = r[1]
            type = r[2]
            start = int(r[3])
            stop = int(r[4])
            strand = r[6]
            attributes = {}
            chromosome = r[0]
            for att in r[-1].split(";"):
                att = att.split("=")
                attributes[att[0]] = att[1]
            if type=="gene":
                mRNA_part = 0
                self.genes[attributes["ID"]] = {'chromosome':chromosome, 'strand':strand, 'data':{}, 'attributes':attributes}
            if type=="mRNA":
                mRNA_part = 1
                gene_id = attributes["Parent"]
                mRNA_id = attributes["ID"]
                gene_data = self.genes.get(gene_id)["data"]
                gene_data[mRNA_id] = {'exons':[], 'CDS':[], 'attributes':attributes}
                self.genes[gene_id]["data"] = gene_data
                self.mRNA_genes[attributes["ID"]] = attributes["Parent"]
            if type=="pseudogene" or type=="tRNA":
                mRNA_part = 0
            if mRNA_part==0:
                r = f.readline()
                continue
            if type=="CDS":
                gene_id = self.mRNA_genes[attributes["Parent"]]
                mRNA_id = attributes["Parent"]
                self.genes[gene_id]["data"][mRNA_id]["exons"].append((start, stop))
            r = f.readline()
            l+=1

    def write_gtf(self, filename):
        f = open(filename, "wt")
        gene_names = self.genes.keys()
        for gene_id, gene_data in self.genes.items():
            gene_strand = gene_data["strand"]
            gene_chromosome = gene_data["chromosome"]
            transcripts = gene_data["data"]
            for mRNA_id, mRNA_data in transcripts.items():
                for (exon_start, exon_stop) in mRNA_data["exons"]:
                    f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (gene_chromosome, "", "exon", exon_start, exon_stop, ".", gene_strand, ".", "gene_id \"%s\"; transcript_id \"%s\";" % (gene_id, mRNA_id)))
                for mRNA_id, mRNA_data in transcripts.items():
                    for (CDS_start, CDS_stop) in mRNA_data["CDS"]:
                        f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (gene_chromosome, "", "CDS", CDS_start, CDS_stop, ".", gene_strand, ".", "gene_id \"%s\"; transcript_id \"%s\";" % (gene_id, mRNA_id)))

    def return_genes(self):
        pass
