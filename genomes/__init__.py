import os
import bisect
import pybio
from os.path import join as pjoin
import glob
import json
import struct
import itertools
import copy

code = {"R": ["A", "G"], "Y": ["C", "T"], "S": ["G", "C"], "W": ["A", "T"]}
revCode = {'A': 'T', 'T': 'A', 'U': 'A', 'G': 'C', 'C': 'G', 'R': 'Y', 'Y': 'R', 'K': 'M', 'M': 'K', 'S': 'S', 'W': 'W', 'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B', 'N': 'N'}
revCodeRYSW = {'R' : 'Y', 'Y' : 'R', 'S' : 'S', 'W' : 'W'}
chr_uscs_ensembl = {}
chr_ensembl_ucsc = {}

# load existing annotations
genes = {}
intervals = {}
linear = {}

def init():
    load_chr_ucsc_ensembl()

def get_latest_version(species):
    if species=="dm6":
        return "ensembl90"
    if species=="at":
        return "ensembl39"
    if species=="mm9":
        return "ensembl67"
    if species=="mm10":
        return "ensembl90"
    if species=="hg19":
        return "ensembl75"
    if species=="hg38":
        return "ensembl97"
    if species=="rn6":
        return "ensembl91"
    if species=="PR8":
        return ""
    if species=="SC35M":
        return ""
    return "ensembl90" # for all other species

def genomes_list(version="ensembl90"):
    genomes = glob.glob(os.path.join(pybio.path.genomes_folder, "*.assembly.%s" % version))
    return [os.path.basename(x.rstrip(".assembly.%s" % version)) for x in genomes]

def load_chr_ucsc_ensembl():
    for species in ["hg19", "mm10"]:
        version = get_latest_version(species)
        map_filename = os.path.join(pybio.path.genomes_folder, "%s.assembly.%s" % (species, version), "%s.chr.ucsc.ensembl" % species)
        if os.path.exists(map_filename):
            f = open(map_filename, "rt")
            r = f.readline()
            while r:
                r = r.replace("\r", "").replace("\n", "").split("\t")
                ucsc = r[0]
                ensembl = r[1]
                if species=="hg19" and ensembl.find("GL")!=-1:
                    ensembl = ensembl+".1"
                if chr_uscs_ensembl.get(ucsc, None)!=None:
                    assert(chr_uscs_ensembl[ucsc]==ensembl)
                chr_uscs_ensembl[ucsc] = ensembl
                if chr_ensembl_ucsc.get(ensembl, None)!=None:
                    assert(chr_ensembl_ucsc[ensembl]==ucsc)
                chr_ensembl_ucsc[ensembl] = ucsc
                r = f.readline()

def load(species, version=None, force=False):
    """
    Loads the genome annotation into memory

    :param species: the genome id
    :param version: the Ensembl version of the annotation, if applicable
    """

    if genes.get(species, None)!=None and not force: # already loaded?
        return
    print "%s: loading genome" % species
    if version==None:
        version = get_latest_version(species)

    genes_filename = os.path.join(pybio.path.genomes_folder, "%s.annotation.%s" % (species, version), "genes.json")
    genes[species] = json.loads(open(genes_filename).read())
    intervals_filename = os.path.join(pybio.path.genomes_folder, "%s.annotation.%s" % (species, version), "intervals.json")
    intervals[species] = json.loads(open(intervals_filename).read())
    linear_filename = os.path.join(pybio.path.genomes_folder, "%s.annotation.%s" % (species, version), "linear.json")
    linear[species] = json.loads(open(linear_filename).read())

    # json stores [] instead of (), bisect doesnt work with []
    for chrstrand, gi in intervals[species].items():
        new_gi = []
        for e in gi:
            new_gi.append(tuple(e))
        intervals[species][chrstrand] = new_gi

    for chrstrand, gi in linear[species].items():
        new_gi = []
        for e in gi:
            new_gi.append(tuple(e))
        linear[species][chrstrand] = new_gi

def adjust_gene(gid, limit_intervals, dbgenes):
    # limit_intervals must be defined inside gene
    if len(limit_intervals)>1:
        dbgenes[gid]["gene_start"] = [start for (start, stop) in limit_intervals]
        dbgenes[gid]["gene_stop"] = [stop for (start, stop) in limit_intervals]
    else:
        dbgenes[gid]["gene_start"], dbgenes[gid]["gene_stop"] = limit_intervals[0][0], limit_intervals[0][1]
    # now adjust gene intervals
    new_gene_intervals = []
    for (limit_start, limit_stop) in limit_intervals:
        for (start, stop, t) in dbgenes[gid]["gene_intervals"]:
            overlap = pybio.utils.interval_overlap(start, stop, limit_start, limit_stop)
            if overlap>0:
                new_gene_intervals.append([max(start, limit_start), min(stop, limit_stop), t])
    new_gene_intervals.sort()
    dbgenes[gid]["gene_intervals"] = new_gene_intervals

def add_cluster(gid, dbgenes, dbclusters):
    overlaps = []
    for (start, stop), _ in dbclusters.items():
        val = pybio.utils.interval_overlap(start, stop, dbgenes[gid]["gene_start"], dbgenes[gid]["gene_stop"])
        if val>0:
            overlaps.append((start, stop))
    if len(overlaps)>0:
        # join all clusters which the new gene overlaps
        cluster_genes = [dbgenes[gid]["gene_id"]]
        cluster_start, cluster_stop = overlaps[0]
        for (start, stop) in overlaps:
            gene_list = dbclusters[(start, stop)]
            cluster_start = min(start, cluster_start, dbgenes[gid]["gene_start"])
            cluster_stop = max(stop, cluster_stop, dbgenes[gid]["gene_stop"])
            cluster_genes = cluster_genes + dbclusters[(start, stop)]
            del dbclusters[(start, stop)] # delete old cluster
        dbclusters[(cluster_start, cluster_stop)] = cluster_genes # newly created cluster
    else:
        cluster_genes = [dbgenes[gid]["gene_id"]]
        cluster_start, cluster_stop = dbgenes[gid]["gene_start"], dbgenes[gid]["gene_stop"]
        dbclusters[(cluster_start, cluster_stop)] = cluster_genes # newly created cluster

def prepare(species="hg19", version=None):

    if version==None:
        version = get_latest_version(species)

    print "%s: processing annotation" % species
    annotation_folder = os.path.join(pybio.path.genomes_folder, "%s.annotation.%s" % (species, version))
    f_log = open(os.path.join(annotation_folder, "log.txt"), "wt")
    chrs_list = pybio.genomes.chr_list(species, version)
    assert(len(chrs_list)>0)
    chrs_names = [name for (name, size) in chrs_list]
    temp_genes = {}
    temp_intervals = {}

    f_log.write("reading annotation\n")
    # A: read all data; merge transcript records and assign to genes
    # make all coordinates 0-based (Ensembl is 1-based)
    # -1 on all coordinates
    # convert strand: 1 = +, -1 = "-"
    f = pybio.data.TabReader(os.path.join(annotation_folder, "%s.annotation.%s.tab" % (species, version)))
    cline = 0
    while f.readline():
        cline += 1
        if cline%100000==0:
            print "%s: processed %sM annotation rows" % (species, cline/100000)
        for k, item in f.data.items():
            f.data[k.lower()] = item
        utr5_start = utr5_stop = utr3_start = utr3_stop = ""
        if f.data.get("Ensembl Gene ID", None)!=None:
            gene_id = f.data["Ensembl Gene ID"]
        elif f.data.get("Gene stable ID", None)!=None:
            gene_id = f.data["Gene stable ID"]
        if "Gene start (bp)" in f.data.keys():
            gene_start = int(f.data["Gene start (bp)"])-1
        if "Gene Start (bp)" in f.data.keys():
            gene_start = int(f.data["Gene Start (bp)"])-1
        if "Gene End (bp)" in f.data.keys():
            gene_stop = int(f.data["Gene End (bp)"])-1
        if "Gene end (bp)" in f.data.keys():
            gene_stop = int(f.data["Gene end (bp)"])-1
        if f.data.get("Chromosome Name")!=None:
            gene_chr = f.data["Chromosome Name"]
        elif f.data.get("Chromosome/scaffold name")!=None:
            gene_chr = f.data["Chromosome/scaffold name"]
        gene_strand = int(f.data["Strand"])
        gene_strand = "+" if gene_strand==1 else "-"
        if gene_chr not in chrs_names:
            continue
        if "5' UTR start" in f.data:
            if f.data["5' UTR start"]!="":
                utr5_start = int(f.data["5' UTR start"])-1
        if "5' UTR end" in f.data:
            if f.data["5' UTR end"]!="":
                utr5_stop = int(f.data["5' UTR end"])-1
        if "3' UTR start" in f.data:
            if f.data["3' UTR start"]!="":
                utr3_start = int(f.data["3' UTR start"])-1
        if "3' UTR end" in f.data:
            if f.data["3' UTR end"]!="":
                utr3_stop = int(f.data["3' UTR end"])-1

        if "5' UTR Start" in f.data:
            if f.data["5' UTR Start"]!="":
                utr5_start = int(f.data["5' UTR Start"])-1
        if "5' UTR End" in f.data:
            if f.data["5' UTR End"]!="":
                utr5_stop = int(f.data["5' UTR End"])-1
        if "3' UTR Start" in f.data:
            if f.data["3' UTR Start"]!="":
                utr3_start = int(f.data["3' UTR Start"])-1
        if "3' UTR End" in f.data:
            if f.data["3' UTR End"]!="":
                utr3_stop = int(f.data["3' UTR End"])-1

        biotype = f.data.get("Gene Biotype", None)
        if biotype==None:
            biotype = f.data.get("Gene type", "")

        if f.data.get("Exon Chr Start (bp)", None)!=None:
            exon_start = int(f.data["Exon Chr Start (bp)"])-1
        elif f.data.get("Exon region start (bp)", None)!=None:
            exon_start = int(f.data["Exon region start (bp)"])-1

        if f.data.get("Exon Chr End (bp)", None)!=None:
            exon_stop = int(f.data["Exon Chr End (bp)"])-1
        elif f.data.get("Exon region end (bp)", None)!=None:
            exon_stop = int(f.data["Exon region end (bp)"])-1

        if "Constitutive exon" in f.data.keys():
            constitutive = f.data["Constitutive exon"]
        if "Constitutive Exon" in f.data.keys():
            constitutive = f.data["Constitutive Exon"]
        # define gene record dictionary
        geneD = temp_genes.get(gene_id, {})
        geneD["record"] = "gene"
        geneD["gene_id"] = gene_id
        geneD["gene_chr"] = gene_chr
        geneD["gene_strand"] = gene_strand
        geneD["gene_start"] = gene_start
        geneD["gene_stop"] = gene_stop
        geneD.setdefault("utr5", [])
        geneD.setdefault("utr3", [])
        geneD.setdefault("exons", [])
        if f.data.get("Associated Gene Name", None)!=None:
            geneD["gene_name"] = f.data["Associated Gene Name"].replace(" ", "")
        elif f.data.get("Gene name", None)!=None:
            geneD["gene_name"] = f.data["Gene name"].replace(" ", "")
        geneD["gene_biotype"] = biotype
        geneD["gene_length"] = gene_stop - gene_start + 1
        if utr3_start!="" and utr3_stop!="":
            utr3 = geneD.get("utr3")
            utr3.append((utr3_start, utr3_stop))
            geneD["utr3"] = utr3
        if utr5_start!="" and utr5_stop!="":
            utr5 = geneD.get("utr5")
            utr5.append((utr5_start, utr5_stop))
            geneD["utr5"] = utr5
        exons = geneD.get("exons")
        exons.append((exon_start, exon_stop))
        geneD["exons"] = exons
        temp_genes[gene_id] = geneD

    max_intervals = 0
    all_genes = len(temp_genes.keys())

    current_gene = 0
    for gene_id, geneD in temp_genes.items():
        current_gene += 1
        if current_gene%100==0:
            f_log.write("%s: creating intervals: %.2f done\n" % (species, current_gene / float(all_genes)))
        coverage = {}
        coverage_utrs = {} # only 5' and 3'
        geneD["exons"] = pybio.utils.merge_intervals(geneD["exons"])
        geneD["utr5"] = pybio.utils.merge_intervals(geneD["utr5"])
        geneD["utr3"] = pybio.utils.merge_intervals(geneD["utr3"])

        # precedence: 3utr -> 5utr -> exon
        if geneD.get("exons", None)!=None:
            for (exon_start, exon_stop) in geneD["exons"]:
                for i in xrange(exon_start, exon_stop+1):
                    coverage[i] = 'o'
        if geneD.get("utr5", None)!=None:
            for (utr5_start, utr5_stop) in geneD["utr5"]:
                for i in xrange(utr5_start, utr5_stop+1):
                    coverage_utrs[i] = '5'
                    coverage[i] = 'o'
        if geneD.get("utr3", None)!=None:
            for (utr3_start, utr3_stop) in geneD["utr3"]:
                for i in xrange(utr3_start, utr3_stop+1):
                    coverage_utrs[i] = '3'
                    coverage[i] = 'o'

        ints = pybio.utils.coverage_to_intervals(coverage)
        utr_intervals = pybio.utils.coverage_to_intervals(coverage_utrs)
        # add intronic intervals
        all_intervals = [ints[0]]
        for (i1_start, i1_stop, i1_value), (i2_start, i2_stop, i2_value) in zip(ints, ints[1:]):
            if i1_stop+1<i2_start:
                all_intervals.append((i1_stop+1, i2_start-1, 'i'))
            all_intervals.append((i2_start, i2_stop, i2_value))
        assert(all_intervals[0][0]==geneD["gene_start"])
        assert(all_intervals[-1][1]==geneD["gene_stop"])
        geneD["gene_intervals"] = all_intervals
        geneD["gene_utrs"] = utr_intervals

        # delete original annotation intervals
        for feature in ["exons", "utr5", "utr3"]:
            if feature in geneD:
                del geneD[feature]

    # B: all genes are read, now we need to resolve overlapping clusters of genes
    genes_by_chrstrand = {}
    for gene_id, geneD in temp_genes.items():
        chrstrand = "%s:%s" % (geneD["gene_chr"], geneD["gene_strand"])
        L = genes_by_chrstrand.get(chrstrand, [])
        L.append((geneD["gene_start"], geneD["gene_id"]))
        genes_by_chrstrand[chrstrand] = L

    for chrstrand, gene_list in genes_by_chrstrand.items():
        f_log.write("making clusters %s\n" % chrstrand)
        gene_list.sort()
        clusters = {}
        for gene_start, gid in gene_list:
            add_cluster(gid, temp_genes, clusters)

        f_log.write("overlaying genes in clusters %s\n" % chrstrand)
        # superimpose genes in clusters (shorter genes have preference)
        for (start, stop), cluster_genes in clusters.items():
            # change genes that are overlaping (they are inside |clusters|>1)
            if len(cluster_genes)>1: # sort gene in cluster by length
                gene_list_bylen = [(temp_genes[gid]["gene_length"], gid) for gid in cluster_genes]
                gene_list_bylen.sort(reverse=True)
                # give preference to genes that are protein coding
                f_log.write("cluster len = %s, chrstrand = %s\n" % (len(cluster_genes), chrstrand))
                coverage = {}
                for (_, gid) in gene_list_bylen: # first place non protein coding genes
                    if temp_genes[gid]["gene_biotype"]!="protein_coding":
                        f_log.write("%s %s %s\n" % (gid, temp_genes[gid]["gene_start"], temp_genes[gid]["gene_stop"]))
                        for i in xrange(temp_genes[gid]["gene_start"], temp_genes[gid]["gene_stop"]+1):
                            coverage[i] = gid
                for (_, gid) in gene_list_bylen: # second, place protein_coding genes (they have preference)
                    if temp_genes[gid]["gene_biotype"]=="protein_coding":
                        f_log.write("%s %s %s\n" % (gid, temp_genes[gid]["gene_start"], temp_genes[gid]["gene_stop"]))
                        for i in xrange(temp_genes[gid]["gene_start"], temp_genes[gid]["gene_stop"]+1):
                            coverage[i] = gid
                f_log.write("\n")
                coverage = pybio.utils.coverage_to_intervals(coverage)
                # check coverage intervals
                for (i1_start, i1_stop, gid1), (i2_start, i2_stop, gid2) in zip(coverage, coverage[1:]):
                    assert(i1_stop+1==i2_start)

                # it may happen that one gene completelly (exactly) overlaps another
                # get new gene list from coverage
                new_gene_list = [gid for (i_start, i_stop, gid) in coverage]
                delete_genes = set(cluster_genes).difference(set(new_gene_list))
                if len(delete_genes)>0:
                    f_log.write("\n")
                    f_log.write("deleted (lost) genes = %s\n" % (str(delete_genes)))
                    f_log.write("these genes completelly (start to stop) overlap with some other gene from this cluster\n")
                    f_log.write("\n")
                # adjust each gene (recover intervals where the gene is defined
                for gid in new_gene_list:
                    limit_intervals = [(start, stop) for (start, stop, interval_gid) in coverage if interval_gid==gid]
                    adjust_gene(gid, limit_intervals, temp_genes)

        f_log.write("making genome intervals from resolved clusters %s\n" % chrstrand)
        gi_list = temp_intervals.get(chrstrand, [])
        for gid, geneD in temp_genes.items():
            k = "%s:%s" % (geneD["gene_chr"], geneD["gene_strand"])
            if k!=chrstrand:
                continue
            if type(geneD["gene_start"])==list:
                f_log.write("gene %s is split (defined) on intervals:\n" % gid)
                for (start, stop) in zip(geneD["gene_start"], geneD["gene_stop"]):
                    f_log.write("%s %s\n" % (start, stop))
                f_log.write("\n")
            if type(geneD["gene_start"])==list:
                for (gene_start, gene_stop) in zip(geneD["gene_start"], geneD["gene_stop"]):
                    gi_list.append((gene_start, gene_stop, gid))
            else:
                gi_list.append((geneD["gene_start"], geneD["gene_stop"], gid))
        gi_list.sort()
        temp_intervals[chrstrand] = gi_list

    f_log.write("saving genome intervals\n")
    f_log.close()

    f = open(os.path.join(annotation_folder, "intervals.json"), "wt")
    f.write(json.dumps(temp_intervals))
    f.close()

    f = open(os.path.join(annotation_folder, "genes.json"), "wt")
    f.write(json.dumps(temp_genes))
    f.close()

    linear = {}
    f = open(os.path.join(annotation_folder, "linear.json"), "wt")
    for gene_id, gene in temp_genes.items():
        key = "%s:%s" % (gene["gene_chr"], gene["gene_strand"])
        chrstrand = linear.get(key, [])
        for (i_start, i_stop, i_type) in gene["gene_intervals"]:
            chrstrand.append((i_start, i_stop, i_type, gene_id))
        linear[key] = chrstrand
    for key, chrstrand in linear.items():
        chrstrand.sort()
        linear[key] = chrstrand
    f.write(json.dumps(linear))

def annotate(species, chr, strand, pos, extension = 0):
    """
    | Annotates given genomic position.
    | Returns triple (upstream_gene, position_gene, downstream_gene).
    """
    load(species)
    chrstrand = "%s:%s" % (chr, strand)
    genome_linear = linear[species].get(chrstrand, None)
    if genome_linear==None: # no genes on this chrstrand
        return (None, None, None, None, None)
    closest_index, gene_id = find_linear(genome_linear, pos) # if gene_id==None, pos is intergenic, closest_index is upstream linear
    if gene_id!=None:
        _, gid_up = find_left_linear(genome_linear, closest_index)
        _, gid_down = find_right_linear(genome_linear, closest_index)
        interval = (genome_linear[closest_index][0], genome_linear[closest_index][1], genome_linear[closest_index][2]) # convert linear (start, stop, i/o, gene_id) to (start, stop, i/o)
        if strand=="+":
            return (gid_up, gene_id, gid_down, interval, find_feature_type(species, gene_id, interval[2], pos))
        else:
            return (gid_down, gene_id, gid_up, interval, find_feature_type(species, gene_id, interval[2], pos))
    else: # position is intergenic
        gene_up = genome_linear[closest_index]
        index_down, _ = find_right_linear(genome_linear, closest_index)
        gene_down = genome_linear[index_down]

        if strand=="+":
            distance_to_up = abs(pos-gene_up[1])
            distance_to_down = abs(pos-gene_down[0])
        if strand=="-":
            distance_to_up = abs(pos-gene_down[0])
            distance_to_down = abs(pos-gene_up[1])

        # assign intergenic position to upstream gene or return current gene as none and return upstream and downstream genes
        if distance_to_up<=extension and distance_to_up<=distance_to_down:
            _, gid_up = find_left_linear(genome_linear, closest_index)
            gene_id = genome_linear[closest_index][3]
            gene_linear = genome_linear[closest_index]
            interval = (gene_linear[0], gene_linear[1], gene_linear[2]) # convert linear (start, stop, i/o, gene_id) to (start, stop, i/o)
            _, gid_down = find_right_linear(genome_linear, closest_index)
        else:
            gid_up = genome_linear[closest_index][3]
            gene_linear = None
            interval = (None, None, None)
            _, gid_down = find_right_linear(genome_linear, closest_index)

        if strand=="+":
            return (gid_up, gene_id, gid_down, interval, find_feature_type(species, gene_id, interval[2], pos))
        else:
            return (gid_down, gene_id, gid_up, interval, find_feature_type(species, gene_id, interval[2], pos))

def find_feature_type(species, gene_id, interval, pos):
    if interval==None:
        return None
    feature_type = {"o":"exon", "i":"intron"}[interval]
    if interval=="o":
        # detect if it's 3'utr or 5'utr
        for (start, stop, ftype) in genes[species][gene_id]["gene_utrs"]:
            if start<=pos<=stop:
                feature_type = ftype
    return {"exon":"exon", "3":"utr3", "5":"utr5", "intron":"intron"}[feature_type]

def find_linear(genome_linear, pos):
    index = bisect.bisect_left(genome_linear, (pos, "", ""))
    in_gene = False
    gene_id = None
    if (genome_linear[index-1][0]<=pos<=genome_linear[index-1][1]):
        gene_id = genome_linear[index-1][3]
    return index-1, gene_id

def find_left_linear(genome_linear, index):
    gene_id = genome_linear[index][3]
    left_gene_id = gene_id
    left_index = index
    while gene_id==left_gene_id and left_index>0:
        left_index = left_index - 1
        left_gene_id = genome_linear[left_index][3]
    return left_index, left_gene_id

def find_right_linear(genome_linear, index):
    gene_id = genome_linear[index][3]
    right_gene_id = gene_id
    right_index = index
    while gene_id==right_gene_id and right_index<(len(genome_linear)-1):
        right_index = right_index + 1
        right_gene_id = genome_linear[right_index][3]
    return right_index, right_gene_id

# get chromosome list
def chr_list(genome, version=None):
    """
    :param genome: mm9, hg19, ...

    Returns chromosome list (names) of the given genome.
    """

    if version==None:
        version = get_latest_version(genome)

    files = glob.glob(pybio.path.root_folder+"/genomes/%s.assembly.%s/*.string" % (genome, version))
    R = []
    for f in files:
        chr_name = f.rstrip(".string").split("/")[-1]
        chr_size = os.path.getsize(f)
        R.append((chr_name, chr_size))
    files = glob.glob(pybio.path.root_folder+"/genomes/%s.assembly.%s/*.raw" % (genome, version))
    for f in files:
        chr_name = f.rstrip(".raw").split("/")[-1]
        chr_size = os.path.getsize(f)
        R.append((chr_name, chr_size))
    return R

# get genomic sequence
def seq_direct(genome, chr, strand, start, stop, flank="N", version=None):
    """
    :param genome: mm9, hg19, ...
    :param chr: chr1, chrX, ...
    :param start: genomic position, 0-based, right closed, positive integer
    :param stop: genomic position, 0-based, right closed, positive integer

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
    fname = os.path.join(pybio.path.root_folder, "genomes", "%s.assembly.%s" % (genome, version), "%s.string" % chr)
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
    :param genome: mm9, hg19, ...
    :param chr: chr1, chrX, ...
    :param pos: genomic position
    :param start: offset relative to pos
    :param stop: offset relative to pos

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
    #for B in ["R", "Y", "S", "W"]:
    #    for E in ["R", "Y", "S", "W"]:
    for B in ["R", "Y", "S", "W"]:
        for E in ["R", "Y", "S", "W"]:
            for m in temp:
                new_m = B+m+E
                temp2.append(new_m)
    motifs.extend(temp2)
    #return motifs[:4]
    return motifs

def make_motifs_nr(motif_size):
    motifs = []
    z = itertools.product("ATCG", repeat=motif_size)
    for el in z:
        motifs.append("".join(el))
    return motifs

def ensembl_species():
    #https://rest.ensembl.org/info/species?content-type=application/json
    return True
