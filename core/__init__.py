import os
import sys
import gzip
import psutil
import math
import pickle

genes_db = {}
transcripts_db = {}
exons_db = {}
gene_bins_db = {}
gene_bin_size = 100000

process = psutil.Process(os.getpid())

class Gene:

  def __init__(self, gene_id, gene_name, chr, strand, start, stop, atts):
    self.gene_id = gene_id
    self.gene_name = gene_name
    self.chr = chr
    self.strand = strand
    self.start = int(start)-1
    self.stop = int(stop)-1
    self.transcripts = set()
    self.atts = atts

class Transcript:

  def __init__(self, transcript_id, gene_id, start, stop):
    self.transcript_id = transcript_id
    self.gene = genes_db[gene_id]
    self.start = int(start)-1
    self.stop = int(stop)-1
    self.utr3 = None
    self.utr5 = None
    self.exons = set()
    self.gene.transcripts.add(self)

class Exon:

  def __init__(self, exon_id, transcript_id, start, stop):
    self.exon_id = exon_id
    self.transcript = transcripts_db[transcript_id]
    self.start = int(start)-1
    self.stop = int(stop)-1
    self.transcript.exons.add(self)

class Utr3:

  def __init__(self, transcript_id, start, stop):
    self.transcript = transcripts_db[transcript_id]
    self.start = int(start)-1
    self.stop = int(stop)-1
    self.transcript.utr3 = self

class Utr5:

  def __init__(self, transcript_id, start, stop):
    self.transcript = transcripts_db[transcript_id]
    self.start = int(start)-1
    self.stop = int(stop)-1
    self.transcript.utr5 = self

def load():
    global gene_bins_db
    global genes_db
    global transcripts_db
    global exons_db

    print("[pybio] reading gene_bins_db.pickle")
    gene_bins_db = pickle.load(open("gene_bins_db.pickle", "rb"))
    print("[pybio] reading genes_db.pickle")
    genes_db = pickle.load(open("genes_db.pickle", "rb"))

    """
    print("[pybio] reading transcripts_db.pickle")
    #transcripts_db = pickle.load(open("transcripts_db.pickle", "rb"))
    print("[pybio] reading exons_db.pickle")
    #exons_db = pickle.load(open("exons_db.pickle", "rb"))
    """

def annotate(chr, strand, pos):
    genes = set()
    transcripts = set()
    exons = set()
    utr3 = set()
    utr5 = set()
    gene_bin = math.floor(pos/gene_bin_size)
    gene_list = gene_bins_db[chr].get(gene_bin, set())
    for gene_id in gene_list:
        gene = genes_db[gene_id]
        if gene.start<=pos<=gene.stop and gene.chr==chr and gene.strand==strand:
            genes.add(gene)
            for transcript in gene.transcripts:
                if transcript.start<=pos<=transcript.stop:
                    transcripts.add(transcript)
                    if transcript.utr3!=None:
                        if transcript.utr3.start<=pos<=transcript.utr3.stop:
                            utr3.add(transcript.utr3)
                    if transcript.utr5!=None:
                        if transcript.utr5.start<=pos<=transcript.utr5.stop:
                            utr5.add(transcript.utr5)
                    for exon in transcript.exons:
                        if exon.start<=pos<=exon.stop:
                            exons.add(exon)       
    return list(genes), list(transcripts), list(exons), list(utr3), list(utr5)

def prepare():
    chr_list = set()
    gtf_fields = ["chr", "source", "feature", "start", "stop", "a1", "strand", "a2", "atts"]

    f = gzip.open("/projects/site/pred/spliceosome/pybio/genomes/hg38.annotation.ensembl109/Homo_sapiens.GRCh38.109.chr.gtf.gz", "rt")
    clines = 0
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        if r[0].startswith("#"):
            r = f.readline()
            continue
        data = dict(zip(gtf_fields, r))
        atts = {}
        temp = data["atts"].split("; ")
        for el in temp:
            el1, el2 = el.split(" ")[0], " ".join(el.split(" ")[1:])
            atts[el1] = el2[1:-1]
        if data["feature"] == "gene":
            gene = Gene(atts["gene_id"], atts.get("gene_name", ""), data["chr"], data["strand"], data["start"], data["stop"], atts)
            genes_db[atts["gene_id"]] = gene
            chr_list.add(data["chr"])
            gene_bins = set()
            gene_start = int(data["start"])-1
            gene_stop = int(data["stop"])-1
            bin_min = math.floor(gene_start/gene_bin_size)
            bin_max = math.floor(gene_stop/gene_bin_size)
            for gene_bin_number in range(bin_min, bin_max+1):
                gene_bins.add(gene_bin_number)
            for gene_bin in gene_bins:
                gene_bins_db.setdefault(data["chr"], {}).setdefault(gene_bin, set()).add(atts["gene_id"])
        if data["feature"] == "transcript":
            transcript = Transcript(atts["transcript_id"], atts["gene_id"], data["start"], data["stop"])
            transcripts_db[atts["transcript_id"]] = transcript
        if data["feature"] == "exon":
            exon = Exon(atts["exon_id"], atts["transcript_id"], data["start"], data["stop"])
            exons_db[atts["exon_id"]] = exon
        if data["feature"] == "three_prime_utr":
            utr3 = Utr3(atts["transcript_id"], data["start"], data["stop"])
        if data["feature"] == "five_prime_utr":
            utr5 = Utr5(atts["transcript_id"], data["start"], data["stop"])
        r = f.readline()
        clines += 1
        if clines%100000==0:
            print("%.2f M lines parsed, RAM GB used: %.3f" % (clines/1000000.0, process.memory_info().rss/1000000000))
    f.close()
    pickle.dump(gene_bins_db, open("gene_bins_db.pickle", "wb"))

    print("#exons = ", len(exons_db.keys()))
    print("#transcripts = ", len(transcripts_db.keys()))
    print("#genes = ", len(genes_db.keys()))
    print("chromosome list = ", "; ".join(list(chr_list)))

    pickle.dump(exons_db, open("exons_db.pickle", "wb"))
    pickle.dump(transcripts_db, open("transcripts_db.pickle", "wb"))
    pickle.dump(genes_db, open("genes_db.pickle", "wb"))

def test():
    load()
    chr_list = set()
    gtf_fields = ["chr", "source", "feature", "start", "stop", "a1", "strand", "a2", "atts"]

    f = gzip.open("/projects/site/pred/spliceosome/pybio/genomes/hg38.annotation.ensembl109/Homo_sapiens.GRCh38.109.chr.gtf.gz", "rt")
    clines = 0
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        if r[0].startswith("#"):
            r = f.readline()
            continue
        data = dict(zip(gtf_fields, r))
        atts = {}
        temp = data["atts"].split("; ")
        for el in temp:
            el1, el2 = el.split(" ")[0], " ".join(el.split(" ")[1:])
            atts[el1] = el2[1:-1]
        if data["feature"] == "gene":
            gene_len = round((int(data["stop"])-int(data["start"]))/1000)
            if gene_len<2: # test short genes only
                print("test gene", data["chr"], data["strand"], atts["gene_id"], str(gene_len) + "K")
                # positive test
                for pos in range(int(data["start"])-1, int(data["stop"])):
                    genes, _, _, _, _ = annotate(data["chr"], data["strand"], pos)
                    genes = [gene.gene_id for gene in genes]
                    assert(atts["gene_id"] in genes)
                # negative test
                for pos in range(int(data["start"])-100, int(data["start"])-1):
                    genes, _, _, _, _ = annotate(data["chr"], data["strand"], pos)
                    genes = [gene.gene_id for gene in genes]
                    assert(atts["gene_id"] not in genes)
                # negative test
                for pos in range(int(data["stop"]), int(data["stop"])+100):
                    genes, _, _, _, _ = annotate(data["chr"], data["strand"], pos)
                    genes = [gene.gene_id for gene in genes]
                    assert(atts["gene_id"] not in genes)
        if data["feature"] == "transcript":
            transcript_len = round((int(data["stop"])-int(data["start"]))/1000)
            if transcript_len<=2:
                print("test transcript", data["chr"], data["strand"], atts["transcript_id"], str(transcript_len) + "K")
                # positive test
                for pos in range(int(data["start"])-1, int(data["stop"])):
                    _, transcripts, _, _, _ = annotate(data["chr"], data["strand"], pos)
                    transcripts = [transcript.transcript_id for transcript in transcripts]
                    assert(atts["transcript_id"] in transcripts)
                # negative test
                for pos in range(int(data["start"])-100, int(data["start"])-1):
                    _, transcripts, _, _, _ = annotate(data["chr"], data["strand"], pos)
                    transcripts = [transcript.transcript_id for transcript in transcripts]
                    assert(atts["transcript_id"] not in transcripts)
                # negative test
                for pos in range(int(data["stop"]), int(data["stop"])+100):
                    _, transcripts, _, _, _ = annotate(data["chr"], data["strand"], pos)
                    transcripts = [transcript.transcript_id for transcript in transcripts]
                    assert(atts["transcript_id"] not in transcripts)
        if data["feature"] == "exon":
            exon_len = round((int(data["stop"])-int(data["start"]))/1000)
            if exon_len<=2:
                print("test exon", data["chr"], data["strand"], atts["exon_id"], str(exon_len) + "K")
                # positive test
                for pos in range(int(data["start"])-1, int(data["stop"])):
                    _, _, exons, _, _ = annotate(data["chr"], data["strand"], pos)
                    exons = [exon.exon_id for exon in exons]
                    assert(atts["exon_id"] in exons)
                # negative test
                for pos in range(int(data["start"])-100, int(data["start"])-1):
                    _, _, exons, _, _ = annotate(data["chr"], data["strand"], pos)
                    exons = [exon.exon_id for exon in exons]
                    assert(atts["exon_id"] not in exons)
                # negative test
                for pos in range(int(data["stop"]), int(data["stop"])+100):
                    _, _, exons, _, _ = annotate(data["chr"], data["strand"], pos)
                    exons = [exon.exon_id for exon in exons]
                    assert(atts["exon_id"] not in exons)
        r = f.readline()
    f.close()

#prepare()
test()