# Ensembl REST API, https://rest.ensembl.org

import os
import bisect
import pybio
from os.path import join as pjoin
import glob
import json
import struct
import itertools
import copy
import gzip
import psutil
import pickle
import math
import requests

process = psutil.Process(os.getpid())

genes_db = {}
transcripts_db = {}
exons_db = {}
gene_bins_db = {}
gene_bin_size = 100000
species_db = {}
genome_loaded = None

code = {"R": ["A", "G"], "Y": ["C", "T"], "S": ["G", "C"], "W": ["A", "T"]}
revCode = {'A': 'T', 'T': 'A', 'U': 'A', 'G': 'C', 'C': 'G', 'R': 'Y', 'Y': 'R', 'K': 'M', 'M': 'K', 'S': 'S', 'W': 'W', 'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B', 'N': 'N'}
revCodeRYSW = {'R' : 'Y', 'Y' : 'R', 'S' : 'S', 'W' : 'W'}
chr_uscs_ensembl = {}
chr_ensembl_ucsc = {}

# ensembl names
ensembl_longA = {"hg38":"homo_sapiens"}
ensembl_longB = {"hg38":"Homo_sapiens"}
ensembl_longC = {"hg38":"GRCh38"}

# load existing annotations
genes = {}
intervals = {}
linear = {}

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

def init():
    #load_chr_ucsc_ensembl()
    ensembl_species_fname = os.path.join(pybio.config.genomes_folder, "ensembl_species.tab")
    if not os.path.exists(ensembl_species_fname):
        pybio.core.genomes.list_species()
    f = open(os.path.join(pybio.config.genomes_folder, "ensembl_species.tab"), "rt")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        species_db[data["species"]] = {"assembly": data["assembly"]}
        r = f.readline()
    f.close()

def prepare(species="hg38", version=None):
    if version==None:
        version = get_latest_version(species)
    print("[pybio] {species}: processing annotation".format(species=species))

    annotation_folder = os.path.join(pybio.config.genomes_folder, "%s.annotation.ensembl%s" % (species, version))
    chr_list = set()
    gtf_fields = ["chr", "source", "feature", "start", "stop", "a1", "strand", "a2", "atts"]

    gtf_files = glob.glob(os.path.join(annotation_folder, "*.gtf.gz"))
    if len(gtf_files)==0:
        gtf_files = glob.glob(os.path.join(annotation_folder, "*.gtf"))
    gtf_fname = gtf_files[0]
    if gtf_fname.endswith(".gz"):
        f = gzip.open(gtf_fname, "rt")
    else:
        f = open(gtf_fname, "rt")

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
    pickle.dump(gene_bins_db, open(os.path.join(annotation_folder, "gene_bins_db.pickle"), "wb"))

    print("#exons = ", len(exons_db.keys()))
    print("#transcripts = ", len(transcripts_db.keys()))
    print("#genes = ", len(genes_db.keys()))
    print("chromosome list = ", "; ".join(list(chr_list)))

    pickle.dump(exons_db, open(os.path.join(annotation_folder, "exons_db.pickle"), "wb"))
    pickle.dump(transcripts_db, open(os.path.join(annotation_folder, "transcripts_db.pickle"), "wb"))
    pickle.dump(genes_db, open(os.path.join(annotation_folder, "genes_db.pickle"), "wb"))

def seq_direct(genome, chr, strand, start, stop, flank="N", version=None):
    """
    Returns chromosome sequence from [start..stop]
    Note: negative strand returns reverse complement; coordinates are 0-based, left and right closed; start must be < stop;
    """

    if version==None:
        version = get_latest_version(genome)

    if start>stop:
        start, stop = stop, start
    assert(start<=stop)
    seq = ""
    if start<0:
        seq = flank * abs(start)
    start = max(0, start) # start can only be non-negative
    fname = os.path.join(pybio.path.genomes_folder, "%s.assembly.%s" % (genome, version), "%s.string" % chr)
    if not os.path.exists(fname):
        return ""
    f = open(fname, "rt")
    f.seek(start, 0)
    seq += f.read(stop-start+1)
    diff = (stop-start+1) - len(seq)
    seq += flank * diff
    if strand=="-":
        return pybio.sequence.reverse_complement(seq)
    else:
        return seq

# get genomic sequence
def seq(genome, chr, strand, pos, start=0, stop=0, flank="N"):
    """
    Returns chromosome sequence from [pos-start..pos+stop]
    Note: negative strand returns reverse complement; coordinates are 0-based, left and right closed; start must be < stop;
    """
    assert(start<=stop)
    if strand=="+":
        gstart = pos+start
        gstop = pos+stop
    else:
        gstart = pos-stop
        gstop = pos-start
    return seq_direct(genome, chr, strand, gstart, gstop)

def make_motifs(motif_size):
    motifs = []
    z = itertools.product("ATCG", repeat=motif_size)
    for el in z:
        motifs.append("".join(el))
    temp = []
    z = itertools.product("ATCG", repeat=motif_size-2)
    for el in z:
        temp.append("".join(el))
    temp2 = []
    for B in ["R", "Y", "S", "W"]:
        for E in ["R", "Y", "S", "W"]:
            for m in temp:
                new_m = B+m+E
                temp2.append(new_m)
    motifs.extend(temp2)
    return motifs

def make_motifs_nr(motif_size):
    motifs = []
    z = itertools.product("ATCG", repeat=motif_size)
    for el in z:
        motifs.append("".join(el))
    return motifs

def url_exists(path):
    r = requests.head(path)
    return r.status_code == requests.codes.ok

def download_assembly(species="homo_sapiens", ensembl_version=None):
    if ensembl_version==None:
        return False
    species_capital = species.capitalize()
    assembly = species_db[species]["assembly"]

    # first, try primary assembly
    fasta_url = f"https://ftp.ensembl.org/pub/release-{ensembl_version}/fasta/{species}/dna/{species_capital}.{assembly}.dna.primary_assembly.fa.gz"
    # no? download the toplevel
    if not url_exists(fasta_url):
        fasta_url = f"https://ftp.ensembl.org/pub/release-{ensembl_version}/fasta/{species}/dna/{species_capital}.{assembly}.dna.toplevel.fa.gz"

    script = """
cd {gdir}
rm {species}.assembly.ensembl{ensembl_version}/*
mkdir {species}.assembly.ensembl{ensembl_version}
cd {species}.assembly.ensembl{ensembl_version}
wget {fasta_url} -O {species}.fasta.gz
gunzip -f {species}.fasta.gz
python3 -c "import pybio; pybio.data.Fasta('{species}.fasta').split()"
"""

    # download FASTA and split
    print(f"[pybio.core.genomes] download FASTA for {species}.ensembl{ensembl_version}")
    command = script.format(gdir=pybio.config.genomes_folder, fasta_url=fasta_url, species=species, species_capital=species_capital, assembly=assembly, ensembl_version=ensembl_version)
    os.system(command)

def download_annotation(species="homo_sapiens", ensembl_version=None):
    if ensembl_version==None:
        return False
    species_capital = species.capitalize()
    assembly = species_db[species]["assembly"]

    script = """
cd {gdir}
rm {species}.annotation.ensembl{ensembl_version}/*
mkdir {species}.annotation.ensembl{ensembl_version}
cd {species}.annotation.ensembl{ensembl_version}
# https://www.biostars.org/p/279235/#279238
wget ftp://ftp.ensembl.org/pub/release-{ensembl_version}/gtf/{species}/{species_capital}.{assembly}.{ensembl_version}.chr.gtf.gz -O {species_capital}.{assembly}.{ensembl_version}.chr.gtf.gz
"""
    # download GTF
    print(f"[pybio.core.genomes] download annotation GTF for {species}.ensembl{ensembl_version}")
    command = script.format(gdir=pybio.config.genomes_folder, species=species, species_capital=species_capital, assembly=assembly, ensembl_version=ensembl_version)
    os.system(command)

def star_index(species="hg38", ensembl_version=None):
    species_capital = species.capitalize()
    assembly = species_db[species]["assembly"]
    script = """
cd {gdir}
rm {species}.assembly.ensembl{ensembl_version}.star/*
mkdir {species}.assembly.ensembl{ensembl_version}.star
gunzip {species}.annotation.ensembl{ensembl_version}/{species_capital}.{assembly}.{ensembl_version}.chr.gtf.gz
STAR --runMode genomeGenerate --genomeDir {species}.assembly.ensembl{ensembl_version}.star --genomeFastaFiles {species}.assembly.ensembl{ensembl_version}/{species}.fasta --runThreadN 4 --sjdbGTFfile {species}.annotation.ensembl{ensembl_version}/{species_capital}.{assembly}.{ensembl_version}.chr.gtf
gzip -f {species}.annotation.ensembl{ensembl_version}/{species_capital}.{assembly}.{ensembl_version}.chr.gtf
"""
    command = script.format(gdir=pybio.config.genomes_folder, species=species, species_capital=species_capital, assembly=assembly, ensembl_version=ensembl_version)
    os.system(command)

def salmon_index(species="hg38", ensembl_version=None):
    species_capital = species.capitalize()
    assembly = species_db[species]["assembly"]
    script = """
cd {gdir}
rm {species}.transcripts.ensembl{ensembl_version}/*
mkdir {species}.transcripts.ensembl{ensembl_version}
cd {species}.transcripts.ensembl{ensembl_version}
wget ftp://ftp.ensembl.org/pub/release-{ensembl_version}/fasta/{species}/cdna/{species_capital}.{assembly}.cdna.all.fa.gz

cd {gdir}
salmon index -t {species}.transcripts.ensembl{ensembl_version}/{species_capital}.{assembly}.cdna.all.fa.gz -i {species}.transcripts.ensembl{ensembl_version}.salmon
"""
    command = script.format(gdir=pybio.config.genomes_folder, species=species, species_capital=species_capital, assembly=assembly, ensembl_version=ensembl_version)
    os.system(command)

def list_species(ensembl_version=109):
    print("[pybio.core.genomes] getting information about all available species from Ensembl, this is done once and takes 1 minute")
    from bs4 import BeautifulSoup
    import requests
    species_db = {}
    def listFD(url, ext=''):
        page = requests.get(url).text
        soup = BeautifulSoup(page, 'html.parser')
        return [node.get('href')[:-1] for node in soup.find_all('a') if node.get('href').endswith(ext)]
    f = open(os.path.join(pybio.config.genomes_folder, "ensembl_species.tab"), "wt")
    f.write("species\tassembly\tensembl_version\n")
    for species in listFD(f"https://ftp.ensembl.org/pub/release-{ensembl_version}/fasta/", "/")[1:]:
        dna_folder_url = f"https://ftp.ensembl.org/pub/release-{ensembl_version}/fasta/{species}/dna/"
        #fasta_file = listFD(dna_folder_url, ext='.dna.primary_assembly.fa.gz')
        fasta_file = listFD(dna_folder_url, ext='.dna.toplevel.fa.gz')
        if len(fasta_file)==0:
            print(f"[pybio.core.genomes] skipping {species}, no fasta file found")
            continue
        fasta_file = fasta_file[0]
        species_assembly = fasta_file.replace(".dna.toplevel.fa.g", "").split(".")
        species_long = species_assembly[0]
        species_assembly = ".".join(species_assembly[1:])
        species_db[species] = (species, species_assembly)
        f.write(f"{species}\t{species_assembly}\t{ensembl_version}\n")
        assert(species.capitalize()==species_long)
    f.close()
    print("[pybio.core.genomes] species list downloaded to:", os.path.join(pybio.config.genomes_folder, "ensembl_species.tab"))
    print("[pybio.core.genomes] example download of fasta/gtf/annotation for homo_sapiens:")
    print("pybio ensembl homo_sapiens [ensembl_version]")

def load(species, ensembl_version=None):
    global gene_bins_db
    global genes_db
    global transcripts_db
    global exons_db
    global genome_loaded
    if ensembl_version==None:
        ensembl_version = pybio.config.ensembl_version_latest
    print(f"[pybio] loading genome annotation for {species} with Ensembl version {ensembl_version}")
    annotation_folder = os.path.join(pybio.config.genomes_folder, f"{species}.annotation.ensembl{ensembl_version}")
    gene_bins_db = pickle.load(open(os.path.join(annotation_folder, "gene_bins_db.pickle"), "rb"))
    genes_db = pickle.load(open(os.path.join(annotation_folder, "genes_db.pickle"), "rb"))
    genome_loaded = (species, ensembl_version)

def annotate(species, chr, strand, pos, ensembl_version=None):
    global genome_loaded
    if ensembl_version==None:
        ensembl_version = pybio.config.ensembl_version_latest
    if genome_loaded!=(species, ensembl_version):
        load(species, ensembl_version)
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

def test():
    load("homo_sapiens")
    chr_list = set()
    gtf_fields = ["chr", "source", "feature", "start", "stop", "a1", "strand", "a2", "atts"]
    f = gzip.open("genomes/homo_sapiens.annotation.ensembl109/Homo_sapiens.GRCh38.109.chr.gtf.gz", "rt")
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
                    genes, _, _, _, _ = annotate("homo_sapiens", data["chr"], data["strand"], pos)
                    genes = [gene.gene_id for gene in genes]
                    assert(atts["gene_id"] in genes)
                # negative test
                for pos in range(int(data["start"])-100, int(data["start"])-1):
                    genes, _, _, _, _ = annotate("homo_sapiens", data["chr"], data["strand"], pos)
                    genes = [gene.gene_id for gene in genes]
                    assert(atts["gene_id"] not in genes)
                # negative test
                for pos in range(int(data["stop"]), int(data["stop"])+100):
                    genes, _, _, _, _ = annotate("homo_sapiens", data["chr"], data["strand"], pos)
                    genes = [gene.gene_id for gene in genes]
                    assert(atts["gene_id"] not in genes)
        if data["feature"] == "transcript":
            transcript_len = round((int(data["stop"])-int(data["start"]))/1000)
            if transcript_len<=2:
                print("test transcript", data["chr"], data["strand"], atts["transcript_id"], str(transcript_len) + "K")
                # positive test
                for pos in range(int(data["start"])-1, int(data["stop"])):
                    _, transcripts, _, _, _ = annotate("homo_sapiens", data["chr"], data["strand"], pos)
                    transcripts = [transcript.transcript_id for transcript in transcripts]
                    assert(atts["transcript_id"] in transcripts)
                # negative test
                for pos in range(int(data["start"])-100, int(data["start"])-1):
                    _, transcripts, _, _, _ = annotate("homo_sapiens", data["chr"], data["strand"], pos)
                    transcripts = [transcript.transcript_id for transcript in transcripts]
                    assert(atts["transcript_id"] not in transcripts)
                # negative test
                for pos in range(int(data["stop"]), int(data["stop"])+100):
                    _, transcripts, _, _, _ = annotate("homo_sapiens", data["chr"], data["strand"], pos)
                    transcripts = [transcript.transcript_id for transcript in transcripts]
                    assert(atts["transcript_id"] not in transcripts)
        if data["feature"] == "exon":
            exon_len = round((int(data["stop"])-int(data["start"]))/1000)
            if exon_len<=2:
                print("test exon", data["chr"], data["strand"], atts["exon_id"], str(exon_len) + "K")
                # positive test
                for pos in range(int(data["start"])-1, int(data["stop"])):
                    _, _, exons, _, _ = annotate("homo_sapiens", data["chr"], data["strand"], pos)
                    exons = [exon.exon_id for exon in exons]
                    assert(atts["exon_id"] in exons)
                # negative test
                for pos in range(int(data["start"])-100, int(data["start"])-1):
                    _, _, exons, _, _ = annotate("homo_sapiens", data["chr"], data["strand"], pos)
                    exons = [exon.exon_id for exon in exons]
                    assert(atts["exon_id"] not in exons)
                # negative test
                for pos in range(int(data["stop"]), int(data["stop"])+100):
                    _, _, exons, _, _ = annotate("homo_sapiens", data["chr"], data["strand"], pos)
                    exons = [exon.exon_id for exon in exons]
                    assert(atts["exon_id"] not in exons)
        r = f.readline()
    f.close()
