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
import sys
requests.packages.urllib3.disable_warnings()

process = psutil.Process(os.getpid())

providers_ftp = {}
providers_ftp["ensembl"] = "https://ftp.ensembl.org/pub"
providers_ftp["ensemblgenomes"] = "https://ftp.ensemblgenomes.ebi.ac.uk/pub"

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
    self.exons_list = [] # list of exons, (start, stop), sorted in 5->3 reading direction, start<=stop
    self.gene.transcripts.add(self)

class Exon:

  def __init__(self, exon_id, transcript_id, start, stop):
    self.exon_id = exon_id
    self.transcript = transcripts_db.get(transcript_id, None)
    self.start = int(start)-1
    self.stop = int(stop)-1
    if self.transcript!=None:
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
    genome_species_fname_finished = os.path.join(pybio.config.genomes_folder, "genome_species.tab.finished")
    if not os.path.exists(genome_species_fname_finished):
        pybio.core.genomes.list_species_ensembl()
    f = open(os.path.join(pybio.config.genomes_folder, "genome_species.tab"), "rt")
    header = f.readline().replace("\r", "").replace("\n", "").split("\t")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        data = dict(zip(header, r))
        species_db[data["species"]] = {"display_name":data["display_name"], "assembly": data["assembly"], "genome_version": data["genome_version"], "provider": data["provider"], "provider_subfolder": data["provider_subfolder"]}
        r = f.readline()
    f.close()

def prepare(species="homo_sapiens", genome_version=None):
    print(f"[pybio.core.genomes] {species}: processing annotation".format(species=species))
    provider = species_db[species]["provider"]
    if genome_version==None:
        genome_version = species_db[species]["ensembl_version"]
    annotation_folder = os.path.join(pybio.config.genomes_folder, "%s.annotation.%s" % (species, genome_version))
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
    print(f"[pybio] reading {gtf_fname}")
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
            if atts.get("exon_id", None)==None:
                atts["exon_id"] = f"exon_{data['start']}_{data['stop']}"
            # note: some annotations don't have "transcript" records, just "transcript_id" in their "exon" records
            if transcripts_db.get(atts["transcript_id"], None)==None:
                transcript = Transcript(atts["transcript_id"], atts["gene_id"], data["start"], data["stop"])
                transcripts_db[atts["transcript_id"]] = transcript
            # we also need to update the transcript start/stop if its already there
            else:
                transcript = transcripts_db[atts["transcript_id"]]
                transcript.start = min(int(data["start"])-1, transcript.start)
                transcript.stop = max(int(data["stop"])-1, transcript.stop)
                transcripts_db[atts["transcript_id"]] = transcript
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

    # sort exons
    print("pybio: sorting exons of transcripts to be given in 5'->3' order")
    for gene_id, gene in genes_db.items():
        for transcript in gene.transcripts:
            exons_list = [(exon.start, exon.stop) for exon in transcript.exons]
            if gene.strand=="+":
                exons_list.sort(key = lambda x: x[0])
            else:
                exons_list.sort(key = lambda x: x[1], reverse=True)
            transcript.exons_list = exons_list

    pickle.dump(exons_db, open(os.path.join(annotation_folder, "exons_db.pickle"), "wb"))
    pickle.dump(transcripts_db, open(os.path.join(annotation_folder, "transcripts_db.pickle"), "wb"))
    pickle.dump(genes_db, open(os.path.join(annotation_folder, "genes_db.pickle"), "wb"))
    return 0

def seq_direct(species, chr, strand, start, stop, flank="N", genome_version=None):
    """
    Returns chromosome sequence from [start..stop]
    Note: negative strand returns reverse complement; coordinates are 0-based, left+right inclusive; start must be < stop;
    """
    if genome_version==None:
        genome_version = pybio.core.genomes.species_db.get(species, {}).get("genome_version", None)
    if genome_version==None:
        return ""
    if start>stop:
        start, stop = stop, start
    assert(start<=stop)
    seq = ""
    if start<0:
        seq = flank * abs(start)
    start = max(0, start) # start can only be non-negative
    fasta_fname = os.path.join(pybio.config.genomes_folder, "%s.assembly.%s" % (species, genome_version), "%s.string" % chr)
    if not os.path.exists(fasta_fname):
        return ""
    f = open(fasta_fname, "rt")
    f.seek(start, 0)
    seq += f.read(stop-start+1)
    diff = (stop-start+1) - len(seq)
    seq += flank * diff
    if strand=="-":
        return pybio.sequence.reverse_complement(seq)
    else:
        return seq

# get genomic sequence
def seq(species, chr, strand, pos, start=0, stop=0, flank="N", genome_version=None):
    """
    Returns chromosome sequence from [pos-start..pos+stop]
    Note: negative strand returns reverse complement; coordinates are 0-based, left+right inclusive; start must be < stop;
    """
    assert(start<=stop)
    if strand=="+":
        gstart = pos+start
        gstop = pos+stop
    else:
        gstart = pos-stop
        gstop = pos-start
    return seq_direct(species, chr, strand, gstart, gstop, genome_version=genome_version)

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
    r = requests.head(path, verify=False)
    return r.status_code == requests.codes.ok

def download_assembly(species, genome_version):
    if genome_version==None:
        return False
    species_capital = species.capitalize()
    assembly = species_db[species]["assembly"]
    ensembl_version = genome_version.replace("ensembl", "")
    ensembl_version = ensembl_version.replace("genomes", "")

    ftp_address = providers_ftp[species_db[species]["provider"]]
    subfolder = species_db[species]["provider_subfolder"]
    if subfolder!="":
        ftp_address = f"{ftp_address}/{subfolder}"
    # first, try primary assembly
    fasta_url = f"{ftp_address}/release-{ensembl_version}/fasta/{species}/dna/{species_capital}.{assembly}.dna.primary_assembly.fa.gz"
    # no? download the toplevel
    if not url_exists(fasta_url):
        fasta_url = f"{ftp_address}/release-{ensembl_version}/fasta/{species}/dna/{species_capital}.{assembly}.dna.toplevel.fa.gz"
    # no? download nonchromosomal
    if not url_exists(fasta_url):
        fasta_url = f"{ftp_address}/release-{ensembl_version}/fasta/{species}/dna/{species_capital}.{assembly}.dna.nonchromosomal.fa.gz"
    print(fasta_url)
    script = """
cd {gdir}
rm {species}.assembly.{genome_version}/* >/dev/null 2>&1
mkdir {species}.assembly.{genome_version} >/dev/null 2>&1
cd {species}.assembly.{genome_version}
wget {fasta_url} -O {species}.fasta.gz --no-check-certificate
gunzip -f {species}.fasta.gz
python3 -c "import pybio; pybio.data.Fasta('{species}.fasta').split()"
"""

    # download FASTA and split
    print(f"[pybio.core.genomes] download FASTA for {species}.{genome_version}")
    command = script.format(gdir=pybio.config.genomes_folder, fasta_url=fasta_url, species=species, species_capital=species_capital, assembly=assembly, genome_version=genome_version, ensembl_version=ensembl_version)
    return os.system(command)

def download_annotation(species, genome_version):
    if genome_version==None:
        return False
    species_capital = species.capitalize()
    assembly = species_db[species]["assembly"]
    ensembl_version = genome_version.replace("ensembl", "")
    ensembl_version = ensembl_version.replace("genomes", "")

    ftp_address = providers_ftp[species_db[species]["provider"]]
    subfolder = species_db[species]["provider_subfolder"]
    if subfolder!="":
        ftp_address = f"{ftp_address}/{subfolder}"

    # prepare GFF3 download
    gff3_url = f"{ftp_address}/release-{ensembl_version}/gff3/{species}/{species_capital}.{assembly}.{ensembl_version}.chr.gff3.gz"
    # no? download the toplevel
    if not url_exists(gff3_url):
        gff3_url = f"{ftp_address}/release-{ensembl_version}/gff3/{species}/{species_capital}.{assembly}.{ensembl_version}.gff3.gz"
    script = """
cd {gdir}
rm {species}.annotation.{genome_version}/* >/dev/null 2>&1
mkdir {species}.annotation.{genome_version} >/dev/null 2>&1
cd {species}.annotation.{genome_version}
# https://www.biostars.org/p/279235/#279238
wget {gff3_url} -O {species}.gff3.gz --no-check-certificate
"""
    # download GFF3
    print(f"[pybio.core.genomes] download annotation GFF3 for {species}.{genome_version}")
    command = script.format(gff3_url=gff3_url, gdir=pybio.config.genomes_folder, species=species, genome_version=genome_version)
    print(command)
    os.system(command)

    # prepare GTF download
    gtf_url = f"{ftp_address}/release-{ensembl_version}/gtf/{species}/{species_capital}.{assembly}.{ensembl_version}.chr.gtf.gz"
    # no? download the toplevel
    if not url_exists(gtf_url):
        gtf_url = f"{ftp_address}/release-{ensembl_version}/gtf/{species}/{species_capital}.{assembly}.{ensembl_version}.gtf.gz"

    script = """
cd {gdir}
mkdir {species}.annotation.{genome_version} >/dev/null 2>&1
cd {species}.annotation.{genome_version}
# https://www.biostars.org/p/279235/#279238
wget {gtf_url} -O {species}.gtf.gz --no-check-certificate
"""
    # download GTF
    print(f"[pybio.core.genomes] download annotation GTF for {species}.{genome_version}")
    command = script.format(gtf_url=gtf_url, gdir=pybio.config.genomes_folder, species=species, genome_version=genome_version)
    return os.system(command)

def star_index(species, genome_version, threads=1):
    species_capital = species.capitalize()
    assembly = species_db[species]["assembly"]
    ensembl_version = genome_version.replace("ensembl", "")
    ensembl_version = ensembl_version.replace("genomes", "")
    script = """
cd {gdir}
rm {species}.assembly.{genome_version}.star/* >/dev/null 2>&1
mkdir {species}.assembly.{genome_version}.star >/dev/null 2>&1
cd {species}.assembly.{genome_version}.star
gunzip -k ../{species}.annotation.{genome_version}/{species}.gtf.gz # -k to keep both .gz and uncompressed GTF, some tools require uncompressed GTF
STAR --runMode genomeGenerate --genomeSAindexNbases {genomeSAindexNbases} --genomeDir ../{species}.assembly.{genome_version}.star --genomeFastaFiles ../{species}.assembly.{genome_version}/{species}.fasta --runThreadN {threads} --sjdbGTFfile ../{species}.annotation.{genome_version}/{species}.gtf
"""
    fasta_file = f"{pybio.config.genomes_folder}/{species}.assembly.{genome_version}/{species}.fasta"
    fasta_size = os.path.getsize(fasta_file)
    genomeSAindexNbases = int(min(14, math.log(fasta_size, 2)/2 - 1))
    command = script.format(threads=threads, genomeSAindexNbases=genomeSAindexNbases, gdir=pybio.config.genomes_folder, species=species, species_capital=species_capital, assembly=assembly, ensembl_version=ensembl_version, genome_version=genome_version)
    return os.system(command)

def salmon_index(species, genome_version):
    species_capital = species.capitalize()
    assembly = species_db[species]["assembly"]
    ensembl_version = genome_version.replace("ensembl", "")
    ensembl_version = ensembl_version.replace("genomes", "")
    ftp_address = providers_ftp[species_db[species]["provider"]]
    subfolder = species_db[species]["provider_subfolder"]
    if subfolder!="":
        ftp_address = f"{ftp_address}/{subfolder}"
    
    script = """
cd {gdir}
rm {species}.transcripts.{genome_version}/* >/dev/null 2>&1
mkdir {species}.transcripts.{genome_version} >/dev/null 2>&1
cd {species}.transcripts.{genome_version}
wget {ftp_address}/release-{ensembl_version}/fasta/{species}/cdna/{species_capital}.{assembly}.cdna.all.fa.gz -O {species}.transcripts.fasta.gz --no-check-certificate

cd {gdir}
salmon index -t {species}.transcripts.{genome_version}/{species}.transcripts.fasta.gz -i {species}.transcripts.{genome_version}.salmon
"""
    command = script.format(ftp_address=ftp_address, gdir=pybio.config.genomes_folder, species=species, species_capital=species_capital, assembly=assembly, ensembl_version=ensembl_version, genome_version=genome_version)
    return os.system(command)

def get_latest_ensembl():
    server = "https://rest.ensembl.org"
    ext = "/info/data/?"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"}, verify=False)
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    decoded = r.json()
    return str(decoded["releases"][0])

def get_latest_ensemblgenomes():
    server = "https://rest.ensembl.org"
    ext = "/info/eg_version?"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"}, verify=False)
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    decoded = r.json()
    return str(decoded["version"])

def get_genome_info(genome_id): 
  server = "https://rest.ensembl.org"
  ext = f"/info/genomes/{genome_id}?"
  r = requests.get(server+ext, headers={ "Content-Type" : "application/json"}, verify=False)
  if not r.ok:
      return {}
  decoded = r.json()
  return decoded

def list_species_ensembl(prepared=True):
    from bs4 import BeautifulSoup
    genome_species_fname_finished = os.path.join(pybio.config.genomes_folder, "genome_species.tab.finished")
    os.system(f"rm {genome_species_fname_finished} >/dev/null 2>&1")
    # download prepared file? (do not query Ensembl)
    if prepared:
        if not os.path.exists(pybio.config.genomes_folder):
            os.makedirs(pybio.config.genomes_folder)
        genome_species_online = "https://raw.githubusercontent.com/grexor/pybio/master/ensembl/genome_species.tab"
        r = requests.get(genome_species_online, allow_redirects=True, verify=False)
        open(os.path.join(pybio.config.genomes_folder, "genome_species.tab"), "wb").write(r.content)
    else:
        print("[pybio.core.genomes] Species list from Ensembl; done once and takes ~ 1 minute")
        ensembl_version = get_latest_ensembl()
        ensemblgenomes_version = get_latest_ensemblgenomes()
        species_db = {}
        def listFD(url, ext=''):
            page = requests.get(url, verify=False).text
            soup = BeautifulSoup(page, 'html.parser')
            return [node.get('href')[:-1] for node in soup.find_all('a') if node.get('href').endswith(ext)]
        if not os.path.exists(pybio.config.genomes_folder):
            os.makedirs(pybio.config.genomes_folder)
        f = open(os.path.join(pybio.config.genomes_folder, "genome_species.tab"), "wt")
        f.write("species\tassembly\tprovider\tprovider_subfolder\tgenome_version\tdisplay_name\n")
        for species in listFD(f"https://ftp.ensembl.org/pub/release-{ensembl_version}/fasta/", "/")[1:]:
            print(f"[pybio.core.genomes] Checking {species}" + " "*30, end="\r")
            dna_folder_url = f"https://ftp.ensembl.org/pub/release-{ensembl_version}/fasta/{species}/dna/"
            fasta_file = listFD(dna_folder_url, ext='.dna.primary_assembly.fa.gz')
            if len(fasta_file)==0:
                fasta_file = listFD(dna_folder_url, ext='.dna.toplevel.fa.gz')
            if len(fasta_file)==0:
                fasta_file = listFD(dna_folder_url, ext='.dna.nonchromosomal.fa.gz')
            if len(fasta_file)==0:
                print(f"[pybio.core.genomes] skipping {species}, no fasta file found")
                continue
            fasta_file = fasta_file[0]
            species_assembly = fasta_file.replace(".dna.toplevel.fa.g", "").replace(".dna.primary_assembly.fa.g", "").replace(".dna.nonchromosomal.fa.g", "").split(".")
            species_long = species_assembly[0]
            species_assembly = ".".join(species_assembly[1:])
            species_db[species] = (species, species_assembly)
            genome_data = get_genome_info(species)
            f.write(f"{species}\t{species_assembly}\tensembl\t\tensembl{ensembl_version}\t{genome_data.get('display_name', '')}\n")
            assert(species.capitalize()==species_long)

        for subfolder in ["plants", "fungi", "protists", "metazoa"]:
            print(f"[pybio.core.genomes] Species list from Ensembl Genomes: {subfolder}; done once and takes ~ 1 minute")
            for species in listFD(f"https://ftp.ensemblgenomes.ebi.ac.uk/pub/{subfolder}/release-{ensemblgenomes_version}/fasta/", "/")[1:]:
                print(f"[pybio.core.genomes] Checking {species}" + " "*30, end="\r")
                dna_folder_url = f"https://ftp.ensemblgenomes.ebi.ac.uk/pub/{subfolder}/release-{ensemblgenomes_version}/fasta/{species}/dna/"
                fasta_file = listFD(dna_folder_url, ext='.dna.toplevel.fa.gz')
                if len(fasta_file)==0:
                    print(f"[pybio.core.genomes] Skipping {species}, no fasta file found")
                    continue
                fasta_file = fasta_file[0]
                species_assembly = fasta_file.replace(".dna.toplevel.fa.g", "").split(".")
                species_long = species_assembly[0]
                species_assembly = ".".join(species_assembly[1:])
                species_db[species] = (species, species_assembly)
                genome_data = get_genome_info(species)
                f.write(f"{species}\t{species_assembly}\tensemblgenomes\t{subfolder}\tensemblgenomes{ensemblgenomes_version}\t{genome_data.get('display_name', '')}\n")
                assert(species.capitalize()==species_long)
        f.close()

    os.system(f"touch {genome_species_fname_finished} >/dev/null 2>&1")
    print()
    print("[pybio.core.genomes] Complete species list downloaded to:", os.path.join(pybio.config.genomes_folder, "genome_species.tab"))
    print()
    print("[pybio.core.genomes] Example command to download and process homo_sapiens genome:")
    print("$ pybio genome homo_sapiens 109")

def load(species, genome_version=None):
    global gene_bins_db
    global genes_db
    global transcripts_db
    global exons_db
    global genome_loaded
    if genome_version==None:
        genome_version = pybio.config.ensembl_version_latest
    print(f"[pybio] loading genome annotation for {species} with Ensembl version {genome_version}")
    annotation_folder = os.path.join(pybio.config.genomes_folder, f"{species}.annotation.{genome_version}")
    gene_bins_db = pickle.load(open(os.path.join(annotation_folder, "gene_bins_db.pickle"), "rb"))
    genes_db = pickle.load(open(os.path.join(annotation_folder, "genes_db.pickle"), "rb"))
    genome_loaded = (species, genome_version)

def annotate(species, chr, strand, pos, genome_version=None):
    global genome_loaded
    if genome_version==None:
        genome_version = pybio.config.ensembl_version_latest
    if genome_loaded!=(species, genome_version):
        load(species, genome_version)
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

# find genes (gene_id, gene_name) by exact match of a search term
def find_genes(search):
    results = set()
    for gene_id, gene in genes_db.items():
        if gene.gene_name==search:
            results.add((0, gene_id, gene.gene_name))
        elif gene.gene_name.find(search)!=-1:
            results.add((1, gene_id, gene.gene_name))
    results = list(results)
    results.sort()
    results = [(e1, e2) for (_, e1, e2) in results]
    return results
