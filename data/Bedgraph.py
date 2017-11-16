import gzip
import bisect
import pybio
import os
import json
import pysam
import numpy as np

def mix_sort(x):
    if x.isdigit():
        return int(x)
    else:
        return str(x)

class Bedgraph():
    """
    Class for reading bedgraph files.
    Scaling is done by / total_raw * 1e6.
    Each position can carry metadata, weighted by cpm value.
    """

    def is_int(self, val):
        if type(val) not in [float, int]:
            return False
        return val == int(val)

    def __init__(self, filename=None, genome="ensembl", fixed_cDNA=False, fast=False, min_cDNA=0):
        self.clear()
        if filename!=None:
            self.load(filename, genome=genome, fast=fast, fixed_cDNA=fixed_cDNA, min_cDNA=min_cDNA)

    def clear(self):
        self.raw = {} # raw counts
        self.cpm = {} # cpm (count per million) scaled raw counts: raw / self.total_raw * 1e6
        self.support = {} # set of ids: if cpm >= self.min_cpm set value to id of file
        self.meta = {} # metadata: weighted by cpm data, fraction of various features (tissue, experiment, etc)
        self.total_raw = 0
        self.total_pos = 0
        self.filename = ""

    def overlay(self, template_filename, source_filename, start=-100, stop=25, db_template="raw", db_source="raw", fast=False, genome_template="ensembl", genome_source="ensembl"):
        self.clear()
        bg_template = pybio.data.Bedgraph(template_filename, fast=fast, genome=genome_template)
        bg_source = pybio.data.Bedgraph(source_filename, fast=fast, genome=genome_source)
        data_template = getattr(bg_template, db_template)
        for chr, strand_data in data_template.items():
            for strand, pos_data in strand_data.items():
                for pos in pos_data.keys():
                    if (start<0) and (stop>0): # all this to consider nearby sites and only add /2 distance if closer than start, stop
                        check_start = 2*start
                        check_stop = 2*stop
                        min_start = start
                        min_stop = stop
                        vector = bg_template.get_vector(chr, strand, pos, start=check_start, stop=-1)
                        vector.reverse()
                        for i in range(len(vector)):
                            if vector[i]>0:
                                min_start = max(min_start, -i/2)
                        vector = bg_template.get_vector(chr, strand, pos, start=1, stop=check_stop)
                        for i in range(len(vector)):
                            if vector[i]>0:
                                min_stop = min(min_stop, i/2)
                        cDNA = bg_source.get_region(chr, strand, pos, start=min_start, stop=min_stop, db = db_source)
                    else:
                        cDNA = bg_source.get_region(chr, strand, pos, start=start, stop=stop, db = db_source)
                    self.set_value(chr, strand, pos, cDNA)
                    #print chr, strand, pos, cDNA
                    self.total_raw += cDNA

    def overlay_old(self, template_filename, source_filename, start=-100, stop=25, db_template="raw", db_source="raw", fast=False, genome_template="ensembl", genome_source="ensembl"):
        self.clear()
        bg_template = pybio.data.Bedgraph(template_filename, fast=fast, genome=genome_template)
        bg_source = pybio.data.Bedgraph(source_filename, fast=fast, genome=genome_source)
        data_template = getattr(bg_template, db_template)
        for chr, strand_data in data_template.items():
            for strand, pos_data in strand_data.items():
                for pos in pos_data.keys():
                    cDNA = bg_source.get_region(chr, strand, pos, start=start, stop=stop, db = db_source)
                    self.set_value(chr, strand, pos, cDNA)
                    self.total_raw += cDNA

    def overlay2(self, template_filename, source_filename, start=-100, stop=25, db_template="raw", db_source="raw", fast=False, genome_template="ensembl", genome_source="ensembl", reverse_strand=False):
        self.clear()
        bg_template = pybio.data.Bedgraph(template_filename, fast=fast, genome=genome_template)
        data_template = getattr(bg_template, db_template)
        bam = pysam.AlignmentFile(source_filename, "rb")

        proc = 0
        for chr, strand_data in data_template.items():
            try:
                bam.pileup(chr, 0, 10, truncate=True)
            except:
                print "skipping", chr
                continue
            for strand, pos_data in strand_data.items():
                if strand=="-":
                    start, stop = -start, -stop
                start, stop = min(start, stop), max(start, stop)
                match_strand = strand
                if reverse_strand:
                    match_strand = {"+":"-", "-":"+"}[match_strand]
                for pos in pos_data.keys():
                    max_coverage = 0
                    if proc%1000==0:
                        print template_filename, source_filename, "%sK processed" % (proc/1000)
                    proc += 1
                    for pup in bam.pileup(chr, pos+start, pos+stop+1, truncate=True):
                        present = set()
                        for read in pup.pileups:
                            cigar = read.alignment.cigar
                            cigar_types = [t for (t, v) in cigar]
                            if 3 in cigar_types:
                                continue
                            if (match_strand=="+" and not read.alignment.is_reverse) or (match_strand=="-" and read.alignment.is_reverse):
                                present.add(read.alignment.qname)
                        max_coverage = max(len(present), max_coverage)
                    #if max_coverage>0:
                    #    print chr, strand, pos, max_coverage
                    self.set_value(chr, strand, pos, max_coverage)
                    self.total_raw += max_coverage

    def get_bam_coverage(self, bam_filename, chr, strand, pos, start, stop, already_processed=set()):
        max_coverage = 0
        bam.close()
        return max_coverage

    def load(self, filename, track_id=None, meta=None, fixed_cDNA=False, compute_cpm=True, genome="ensembl", fast=False, force_strand=None, min_cDNA=0):
        """
        Load Bedgraph file (can also be gzipped). A Bedgraph object load method can be called multiple times on various files, the content is added up.
        """
        # if genome != "ensembl", data is loaded from UCSC genome (hg19, mm10) and converted to ensembl
        # load can be called multiple times on different bed files
        # track_id : unique identifier for this file
        # genome: ensembl/ucsc
        print "loading : %s" % filename
        self.filename = filename

        if track_id==None:
            track_id = filename

        temp_raw = {}
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
                chr = pybio.genomes.chr_uscs_ensembl.get(chr, None) # convert ucsc to ensembl
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
            for pos in range(int(r[1]), int(r[2])):
                if abs(raw)>=min_cDNA:
                    temp_raw.setdefault(chr, {}).setdefault(strand, {}).setdefault(pos, 0)
                    temp_raw[chr][strand][pos] += abs(raw)
                    self.total_raw += abs(raw)
                    self.total_pos += 1
            r = f.readline()
        f.close()
        self.genome = "ensembl" # we converted to ensembl

        if not_converted>0:
            print "%s positions were skipped; unable to convert chromosome name to ensembl format" % not_converted

        for chr, strand_data in temp_raw.items():
            for strand, pos_data in strand_data.items():
                for pos, raw in pos_data.items():
                    self.add_value(chr, strand, pos, raw, db="raw")
                    if meta!=None:
                        self.add_meta(chr, strand, pos, meta, raw)

        del temp_raw
        #self.add_value(chr, strand, pos, track_id, db="support")
        return

    # normalize (cpm) and consider thresholds (min_cpm, min_raw)
    def norm(self):
        for chr, strand_data in self.raw.items():
            for strand, pos_data in strand_data.items():
                for pos, raw in pos_data.items():
                    cpm = raw / float(self.total_raw) * 1e6
                    self.add_value(chr, strand, pos, cpm, db="cpm")

    def save(self, filename, db_save="raw", min_raw=0, min_support=0, min_cpm=0, filetype="bed", genome="ensembl", track_id=None, without_track=False):
        # IF genome != "ensembl", convert chromosome names to UCSC
        # Get positions from db_source, check min_raw agains db_source
        # Save raw value from db_save
        # Usually db_source and db_save are the same
        """
        :param without_track: do not include track names
        Save bedgraph data to file.
        """
        self.filename = filename # set the new filename
        print "save : %s, format : %s" % (filename, filetype)
        if track_id!=None:
            self.track_id = track_id
        else:
            self.track_id = os.path.basename(filename)
        not_converted = 0
        if filename.endswith(".gz"):
            f_out = gzip.open(filename, "wb")
        else:
            f_out = open(filename, "wt")
        if filetype=="complete":
            header = ["chr", "strand", "pos", "raw", "cpm", "support", "meta"]
            f_out.write("\t".join(header)+"\n")
        else:
            if not without_track:
                f_out.write('track type=bedGraph name="%s" description="%s" altColor="200,120,59" color="120,101,172" maxHeightPixels="100:50:0" visibility="full" priority="20"\n' % (self.track_id, self.track_id))
        data_save = getattr(self, db_save)
        chr_names = sorted(data_save.keys(), key=mix_sort)
        for chr in chr_names:
            chr_data = data_save[chr]
            for strand, strand_data in chr_data.items():
                positions = [(pos, raw) for pos, raw in strand_data.items()]
                positions.sort()
                cluster_start = None
                for (pos, val), (pos2, val2) in zip(positions, positions[1:]+[(None, None)]):
                    strand_str = "" if strand=="+" else "-"
                    if db_save=="support":
                        val = len(val)
                    raw = self.get_value(chr, strand, pos, db="raw")
                    cpm = self.get_value(chr, strand, pos, db="cpm")
                    support = len(self.get_value(chr, strand, pos, db="support"))
                    meta = self.get_value(chr, strand, pos, db="meta")
                    if raw>=min_raw and cpm>=min_cpm and support>=min_support:
                        chr_str = chr
                        if genome=="ucsc":
                            chr_str = pybio.genomes.chr_ensembl_ucsc.get(chr_str, None) # convert ensembl to ucsc
                            if chr_str==None:
                                not_converted += 1
                                continue
                        if filetype=="bed":
                            if pos2==pos+1 and val==val2:
                                if cluster_start==None:
                                    cluster_start = pos
                                continue
                            if cluster_start!=None:
                                row = [chr_str, str(cluster_start), str(pos+1)]
                                cluster_start = None
                            else:
                                row = [chr_str, str(pos), str(pos+1)]
                            if type(val)==float:
                                row.append("%s%.2f" % (strand_str, val))
                            else:
                                row.append("%s%s" % (strand_str, val))
                        if filetype=="complete":
                            for mtype, mvalue in meta.items():
                                meta[mtype] = "%.2f" % meta[mtype]
                            row = [chr_str, strand, str(pos), str(raw), "%.2f" % cpm, str(support), json.dumps(meta)]
                        f_out.write("\t".join(row)+"\n")
        if not_converted>0:
            print "%s not converted during save" % not_converted
            f_out.write("#not converted from ensembl chr names = %s" % not_converted)
        f_out.close()
        if genome!="ucsc":
            if filename.endswith(".bed"):
                self.save(filename.replace(".bed", "_ucsc.bed"), db_save=db_save, min_raw=min_raw, min_support=min_support, min_cpm=min_cpm, filetype=filetype, genome="ucsc")

    def set_value(self, chr, strand, pos, val, db="raw"):
        if val==0:
            return
        data = getattr(self, db)
        if db in ["raw", "cpm", "meta"]:
            data.setdefault(chr, {}).setdefault(strand, {}).setdefault(pos, abs(val))
            data[chr][strand][pos] = abs(val)
        if db in ["support"]:
            data.setdefault(chr, {}).setdefault(strand, {}).setdefault(pos, set(val))
            data[chr][strand][pos] = set(val)

    def add_value(self, chr, strand, pos, val, db="raw"):
        if val==0:
            return
        data = getattr(self, db)
        if db in ["raw", "cpm"]:
            data.setdefault(chr, {}).setdefault(strand, {}).setdefault(pos, 0)
            data[chr][strand][pos] += abs(val)
        if db in ["support"]:
            data.setdefault(chr, {}).setdefault(strand, {}).setdefault(pos, set())
            data[chr][strand][pos].add(val)

    def add_meta(self, chr, strand, pos, meta, val):
        # each meta type has its own value, additional parameter, so thats why a separate function compared to add_value
        if val==0 or meta=={} or meta==None:
            return
        data = getattr(self, "meta")
        data.setdefault(chr, {}).setdefault(strand, {}).setdefault(pos, {})
        for mt in meta:
            data[chr][strand][pos] = meta

    def get_value(self, chr, strand, pos, db="raw"):
        """
        Get value at chromosome (chr), strand, position.
        """
        data = getattr(self, db)
        if db in ["raw", "cpm"]:
            return data.get(chr, {}).get(strand, {}).get(pos, 0)
        if db in ["support"]:
            return data.get(chr, {}).get(strand, {}).get(pos, set())
        if db in ["meta"]:
            return data.get(chr, {}).get(strand, {}).get(pos, "")

    def fetch(self, db="raw"):
        data = getattr(self, db)
        for chr, chr_data in data.items():
            for strand, pos_data in chr_data.items():
                for pos, cDNA in pos_data.items():
                    yield (chr, strand, pos, cDNA)

    def set_region(self, chr, strand, start, stop, val, db="raw"):
        # only for raw and cpm
        data = getattr(self, db)
        data.setdefault(chr, {}).setdefault(strand, {})
        for i in range(start, stop+1):
            data[chr][strand][i] = val

    def get_region(self, chr, strand, pos, start=0, stop=0, db="raw"):
        # only for raw and cpm
        data = getattr(self, db)
        if strand=="-":
            start, stop = -start, -stop
        start, stop = min(start, stop), max(start, stop)
        region_sum = 0
        for i in range(pos+start, pos+stop+1):
            region_sum += data.get(chr, {}).get(strand, {}).get(i, 0)
        return region_sum

    def get_vector(self, chr, strand, pos, start, stop, db="raw"):
        # returns vector of values around pos
        # considers strand but does not reverse the vector
        # example: strand=-, start=-10, stop=20, would return a vector of pos-20 .. pos+10
        # example: strand=+, start=-10, stop=20, would return a vector of pos-10 .. pos+20
        vector = []
        # only for raw and cpm
        data = getattr(self, db)
        if strand=="-":
            start, stop = -start, -stop
        start, stop = min(start, stop), max(start, stop)
        for i in range(pos+start, pos+stop+1):
            vector.append(data.get(chr, {}).get(strand, {}).get(i, 0))
        return vector

    def cluster(self, region_up=150, region_down=150):
        # go down position list: highly expressed -> lowly expressed
        # sum [region_up, region_down] and assign to main position
        """
        Cluster positions in bedgraph. For each chromosome & strand, sort positions by value (descending), and assign sum(region_up, region_down) around each position.
        """
        for chr, strand_data in self.raw.items():
            for strand, pos_data in strand_data.items():
                print "clustering : %s %s" % (chr, strand)
                offset_up, offset_down = region_up, region_down
                if strand=="-":
                    offset_up, offset_down = offset_down, offset_up
                mp = 1 if strand=="+" else -1 # if same value, take downstream position first
                positions = [(raw, mp*pos, pos) for pos, raw in pos_data.items()]
                positions = sorted(positions, reverse=True)
                temp_raw = {}
                for _, _, main_pos in positions:
                    if self.get_value(chr, strand, main_pos, db="raw")==0: # has this position been deleted (added or ignored)? move on
                        continue
                    new_raw = 0
                    for i in range(main_pos-offset_up, main_pos+offset_down+1):
                        raw = self.get_value(chr, strand, i, db="raw")
                        if raw!=0:
                            new_raw += raw
                            del self.raw[chr][strand][i]
                    temp_raw[main_pos] = new_raw

                self.raw[chr][strand] = temp_raw

    def filter(self, min_distance=125):
        # filter on raw data
        # sort positions and go down the list
        # only keep positions that are min_distance apart (don't sum anything)
        """
        Sort chromosome + strand positions by value (descending). Go down the list and only keep positions that are min_distance apart.
        """
        for chr, strand_data in self.raw.items():
            for strand, pos_data in strand_data.items():
                print "filtering: %s%s (min distance = %snt)" % (strand, chr, min_distance)
                mp = 1 if strand=="+" else -1 # if same value, take downstream position first
                positions = [(raw, mp*pos, pos) for pos, raw in pos_data.items()]
                positions = sorted(positions, reverse=True)
                temp_raw = {}
                for _, _, main_pos in positions:
                    if self.get_value(chr, strand, main_pos)==0: # has this position been deleted (added or ignored)? move on
                        continue
                    temp_raw[main_pos] = self.get_value(chr, strand, main_pos)
                    for i in range(main_pos-min_distance, main_pos+min_distance+1): # delete main_pos and surroundings [-min_distance, min_distance]
                        if self.get_value(chr, strand, i)!=0:
                            del self.raw[chr][strand][i]
                self.raw[chr][strand] = temp_raw
