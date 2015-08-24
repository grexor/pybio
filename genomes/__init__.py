import os
import bisect
import pybio
from os.path import join as pjoin
import glob
import json
import struct
import itertools

code = {"R": ["A", "G"], "Y": ["C", "T"], "S": ["G", "C"], "W": ["A", "T"]}
revCode = {'A': 'T', 'T': 'A', 'U': 'A', 'G': 'C', 'C': 'G', 'R': 'Y', 'Y': 'R', 'K': 'M', 'M': 'K', 'S': 'S', 'W': 'W', 'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B', 'N': 'N'}
revCodeRYSW = {'R' : 'Y', 'Y' : 'R', 'S' : 'S', 'W' : 'W'}
chr_uscs_ensembl = {}
chr_ensembl_ucsc = {}

# load existing annotations
genes = {}
intervals = {}

def init():
    load_chr_ucsc_ensembl("hg19")
    load_chr_ucsc_ensembl("mm10")

def genomes_list():
    genomes = glob.glob(os.path.join(pybio.path.genomes_folder, "*.assembly"))
    return [os.path.basename(x.rstrip(".assembly")) for x in genomes]

def load_chr_ucsc_ensembl(species):
    if chr_uscs_ensembl.get(species, None)==None:
        chr_uscs_ensembl[species] = {}
        chr_ensembl_ucsc[species] = {}
        map_filename = os.path.join(pybio.path.genomes_folder, "%s.assembly" % species, "%s.chr.ucsc.ensembl" % species)
        if os.path.exists(map_filename):
            f = open(map_filename, "rt")
            r = f.readline()
            while r:
                r = r.replace("\r", "").replace("\n", "").split("\t")
                ucsc = r[0]
                ensembl = r[1]
                if species=="hg19" and ensembl.find("GL")!=-1:
                    ensembl = ensembl+".1"
                chr_uscs_ensembl[species][ucsc] = ensembl
                chr_ensembl_ucsc[species][ensembl] = ucsc
                r = f.readline()

def load(species, version="ensembl74"):
    """
    Loads the genome annotation into memory

    :param species: the genome id
    :param version: the Ensembl version of the annotation, if applicable
    """

    if genes.get(species, None)!=None: # already loaded?
        return
    genes_filename = os.path.join(pybio.path.genomes_folder, "%s.annotation.%s" % (species, version), "genes.json")
    genes[species] = json.loads(open(genes_filename).read())
    intervals_filename = os.path.join(pybio.path.genomes_folder, "%s.annotation.%s" % (species, version), "intervals.json")
    intervals[species] = json.loads(open(intervals_filename).read())
    # json stores [] instead of (), bisect doesnt work with []
    for chrstrand, gi in intervals[species].items():
        new_gi = []
        for e in gi:
            new_gi.append(tuple(e))
        intervals[species][chrstrand] = new_gi

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

def prepare(species="hg19", version="ensembl74"):
    if species=="mm9":
        version = "ensembl67" # latest ensembl annotation for mm9
    print "preparing annotation for %s" % species
    annotation_folder = os.path.join(pybio.path.genomes_folder, "%s.annotation.%s" % (species, version))
    f_log = open(os.path.join(annotation_folder, "log.txt"), "wt")
    chrs_list = pybio.genomes.chr_list(species)
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
    while f.readline():
        utr5_start = utr5_stop = utr3_start = utr3_stop = ""
        gene_id = f.data["Ensembl Gene ID"]
        gene_start = int(f.data["Gene Start (bp)"])-1
        gene_stop = int(f.data["Gene End (bp)"])-1
        gene_chr = f.data["Chromosome Name"]
        gene_strand = int(f.data["Strand"])
        gene_strand = "+" if gene_strand==1 else "-"
        if gene_chr not in chrs_names:
            continue
        if f.data["5' UTR Start"]!="":
            utr5_start = int(f.data["5' UTR Start"])-1
        if f.data["5' UTR End"]!="":
            utr5_stop = int(f.data["5' UTR End"])-1
        if f.data["3' UTR Start"]!="":
            utr3_start = int(f.data["3' UTR Start"])-1
        if f.data["3' UTR End"]!="":
            utr3_stop = int(f.data["3' UTR End"])-1
        biotype = f.data["Gene Biotype"]
        exon_start = int(f.data["Exon Chr Start (bp)"])-1
        exon_stop = int(f.data["Exon Chr End (bp)"])-1
        constitutive = f.data["Constitutive Exon"]
        # define gene record dictionary
        geneD = temp_genes.get(gene_id, {})
        geneD["record"] = "gene"
        geneD["gene_id"] = gene_id
        geneD["gene_chr"] = gene_chr
        geneD["gene_strand"] = gene_strand
        geneD["gene_start"] = gene_start
        geneD["gene_stop"] = gene_stop
        geneD["gene_name"] = f.data.get("Associated Gene Name", "")
        geneD["gene_biotype"] = f.data["Gene Biotype"]
        geneD["gene_length"] = gene_stop - gene_start + 1
        if utr3_start!="" and utr3_stop!="":
            utr3 = geneD.get("utr3", [])
            utr3.append((utr3_start, utr3_stop))
            geneD["utr3"] = utr3
        if utr5_start!="" and utr5_stop!="":
            utr5 = geneD.get("utr5", [])
            utr5.append((utr5_start, utr5_stop))
            geneD["utr5"] = utr5
        exons = geneD.get("exons", [])
        exons.append((exon_start, exon_stop))
        geneD["exons"] = exons
        temp_genes[gene_id] = geneD

    max_intervals = 0
    all_genes = len(temp_genes.keys())

    current_gene = 0
    for gene_id, geneD in temp_genes.items():
        current_gene += 1
        if current_gene%100==0:
            f_log.write("creating intervals: %.2f done\n" % (current_gene / float(all_genes)))
        coverage = {}
        geneD["exons"] = pybio.utils.merge_intervals(geneD["exons"])
        if geneD.get("utr3", None)!=None:
            for (utr3_start, utr3_stop) in geneD["utr3"]:
                for i in xrange(utr3_start, utr3_stop+1):
                    #coverage[i] = '3'
                    coverage[i] = 'o'
        if geneD.get("utr5", None)!=None:
            for (utr5_start, utr5_stop) in geneD["utr5"]:
                for i in xrange(utr5_start, utr5_stop+1):
                    if coverage.get(i, None)!=None:
                        continue
                    #coverage[i] = '5'
                    coverage[i] = 'o'
        if geneD.get("exons", None)!=None:
            for (exon_start, exon_stop) in geneD["exons"]:
                for i in xrange(exon_start, exon_stop+1):
                    if coverage.get(i, None)!=None:
                        continue
                    coverage[i] = 'o'

        ints = pybio.utils.coverage_to_intervals(coverage)
        # add intronic intervals
        all_intervals = [ints[0]]
        for (i1_start, i1_stop, i1_value), (i2_start, i2_stop, i2_value) in zip(ints, ints[1:]):
            if i1_stop+1<i2_start:
                all_intervals.append((i1_stop+1, i2_start-1, 'i'))
            all_intervals.append((i2_start, i2_stop, i2_value))
        assert(all_intervals[0][0]==geneD["gene_start"])
        assert(all_intervals[-1][1]==geneD["gene_stop"])
        geneD["gene_intervals"] = all_intervals
        # delete old intervals
        del geneD["exons"]
        if geneD.get("utr5", None)!=None:
            del geneD["utr5"]
        if geneD.get("utr3", None)!=None:
            del geneD["utr3"]

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

def annotate_position(species, chr, strand, pos):
    """
    | Annotages given genomic position.
    | Returns triple (upstream_gene, position_gene, downstream_gene).
    """
    chrstrand = "%s:%s" % (chr, strand)
    gi = intervals[species].get(chrstrand, None)
    if gi==None: # no genes on this chrstrand
        return (None, None, None)
    i, gid = find_gene(gi, pos)
    # bisect returns the point of insertion, we still need to check if pos is inside gene
    # if intergenic, bisect returns the upstream gene
    pos_in_gene = False
    if type(genes[species][gid]["gene_start"])==list:
        for (start, stop) in zip(genes[species][gid]["gene_start"], genes[species][gid]["gene_stop"]):
            if start<=pos<=stop:
                pos_in_gene = True
    if genes[species][gid]["gene_start"]<=pos<=genes[species][gid]["gene_stop"]:
        pos_in_gene = True
    if pos_in_gene:
        if i>0:
            gid_up = gi[i-1][2]
        else:
            gid_up = None
        if i<len(gi)-1:
            gid_down = gi[i+1][2]
        else:
            gid_down = None
        if strand=="+":
            return (gid_up, gid, gid_down)
        else:
            return (gid_down, gid, gid_up)
    else: # position not in gene
        gid = None
        if i>=0:
            gid_up = gi[i][2] # bisect point of insertion
        else:
            gid_up = None
        if i<len(gi)-1:
            gid_down = gi[i+1][2]
        else:
            gid_down = None
        if strand=="+":
            return (gid_up, gid, gid_down)
        else:
            return (gid_down, gid, gid_up)

def find_gene(gi, pos):
    # print gi[:12]
    i = bisect.bisect_left(gi, (pos, "", ""))
    return i-1, gi[i-1][2]

def find_feature(species, gid, pos):
    intervals = genes[species][gid]["gene_intervals"]
    for (start, stop, t) in intervals:
        if start<=pos<=stop:
            return (start, stop, t)
    return None

def annotate(species, chr, strand, pos):
    load(species)
    pos = int(pos)
    up_gene, gid, down_gene = annotate_position(species, chr, strand, pos)
    interval = None
    if gid:
        interval = find_feature(species, gid, pos)
    return (up_gene, gid, down_gene, interval)

# get chromosome list
def chr_list(genome):
    """
    :param genome: mm9, hg19, ...

    Returns chromosome list (names) of the given genome.
    """
    files = glob.glob(pybio.path.root_folder+"/genomes/%s.assembly/*.string" % (genome))
    R = []
    for f in files:
        chr_name = f.rstrip(".string").split("/")[-1]
        chr_size = os.path.getsize(f)
        R.append((chr_name, chr_size))
    return R

# get genomic sequence
def seq(genome, chr, strand, start, stop, flank="N"):
    """
    :param genome: mm9, hg19, ...
    :param chr: chr1, chrX, ...
    :param start: 0-based
    :param stop: 0-based, right closed

    Returns chromosome sequence from start_pos to end_pos. If strand is negative, return reverse complement
    Coordinates are 0-based, inclusive
    """
    if start>stop: # always positive orientation
        start, stop = stop, start
    chrs = chr_list(genome)
    seq = ""
    if start<0:
        seq = flank * abs(start)
    start = max(0, start) # start can only be non-negative
    fname = os.path.join(pybio.path.root_folder, "genomes", "%s.assembly" % genome, "%s.string" % chr)
    if not os.path.exists(fname):
        return ""
    f = open(fname, "rt")
    f.seek(start, 0)
    seq += f.read(stop-start+1)
    diff = (stop-start+1) - len(seq)
    seq += flank * diff
    if strand=="-":
        return pybio.sequence.reverse_complement(seq)
    return seq

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
