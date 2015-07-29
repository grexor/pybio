import gzip
import bisect
import pybio
import os
import json

def mix_sort(x):
    if x.isdigit():
        return int(x)
    else:
        return str(x)

class Bedgraph():
    """
    Class for reading bedgraph files.
    Scaling is done by / total_raw * 1e6
    Each position can carry metadata, weighted by cpm value.
    """

    def is_int(self, val):
        if type(val) not in [float, int]:
            return False
        return val == int(val)

    def __init__(self, filename=None, genome="ensembl", fast=False):
        self.clear()
        if filename!=None:
            self.load(filename, genome=genome, fast=fast)

    def clear(self):
        self.raw = {} # raw counts
        self.cpm = {} # cpm (count per million) scaled raw counts: raw / self.total_raw * 1e6
        self.support = {} # set of ids: if cpm >= self.min_cpm set value to id of file
        self.meta = {} # metadata: weighted by cpm data, fraction of various features (tissue, experiment, etc)
        self.total_raw = 0

    def overlay(self, template_filename, source_filename, region_up=100, region_down=25, db_template="raw", db_source="raw", fast=False, genome_template="ensembl", genome_source="ensembl"):
        self.clear()
        bg_template = pybio.data.Bedgraph(template_filename, fast=fast, genome=genome_template)
        bg_source = pybio.data.Bedgraph(source_filename, fast=fast, genome=genome_source)
        data_template = getattr(bg_template, db_template)
        for chr, strand_data in data_template.items():
            for strand, pos_data in strand_data.items():
                offset_up, offset_down = region_up, region_down
                if strand=="-":
                    offset_up, offset_down = offset_down, offset_up
                for pos, val in pos_data.items():
                    raw = 0
                    for i in range(pos-offset_up, pos+offset_down+1):
                        raw += bg_source.get_value(chr, strand, i, db = db_source)
                    self.set_value(chr, strand, pos, raw)

    def load(self, filename, track_id=None, meta=None, min_cpm=0, min_raw=0, compute_cpm=True, genome="ensembl", fast=False, force_strand=None):
        # if genome != "ensembl", data is loaded from UCSC genome (hg19, mm10) and converted to ensembl
        # load can be called multiple times on different bed files
        # track_id : unique identifier for this file
        # genome: ensembl/ucsc
        if track_id==None:
            track_id = filename

        temp_raw = {}
        raw_sum = 0
        print "loading : %s" % filename
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
                chr = pybio.genomes.chr_uscs_ensembl.get(genome, {}).get(chr, None) # always convert ucsc to ensembl
                if chr==None: # the mapping of chromosome from ucsc to ensembl didnt success
                    not_converted += 1
                    r = f.readline()
                    continue
            raw = float(r[3])
            # if reading integer numbers, make them integer
            if self.is_int(raw):
                raw = int(raw)
            strand = "+" if raw>=0 else "-"
            if force_strand!=None:
                assert(force_strand in ["+", "-"])
                strand = force_strand
            for pos in range(int(r[1]), int(r[2])):
                temp_raw.setdefault(chr, {}).setdefault(strand, {}).setdefault(pos, 0)
                temp_raw[chr][strand][pos] += abs(raw)
                raw_sum += abs(raw)
            self.total_raw += abs(raw)
            r = f.readline()
        f.close()
        self.genome = "ensembl" # we converted to ensembl

        if not_converted>0:
            print "%s positions were skipped; unable to convert chromosome name to ensembl format" % not_converted

        if fast:
            self.raw = temp_raw
            return

        # scaling (cpm) and consider thresholds (min_cpm, min_raw)
        for chr, strand_data in temp_raw.items():
            for strand, pos_data in strand_data.items():
                for pos, raw in pos_data.items():
                    cpm = raw / float(raw_sum) * 1e6
                    if raw>=min_raw and cpm>=min_cpm:
                        self.add_value(chr, strand, pos, raw, db="raw")
                        if compute_cpm: # add cpm, support, meta?
                            self.add_value(chr, strand, pos, cpm, db="cpm")
                            self.add_meta(chr, strand, pos, meta, cpm)
                            self.add_value(chr, strand, pos, track_id, db="support")
        del temp_raw

    def save(self, filename, db_save="raw", min_raw=0, min_support=0, min_cpm=0, filetype="bed", genome="ensembl", track_id=None):
        # IF genome != "ensembl", convert chromosome names to UCSC
        # Get positions from db_source, check min_raw agains db_source
        # Save raw value from db_save
        # Usually db_source and db_save are the same
        print "save : %s, format : %s" % (filename, filetype)
        if track_id!=None:
            self.track_id = track_id
        else:
            self.track_id = os.path.basename(filename)
        not_converted = 0
        f_out = open(filename, "wt")
        if filetype=="complete":
            header = ["chr", "strand", "pos", "raw", "cpm", "support", "meta"]
            f_out.write("\t".join(header)+"\n")
        else:
            f_out.write('track type=bedGraph name="%s" description="%s" altColor="200,120,59" color="120,101,172" maxHeightPixels="100:50:0" visibility="full" priority="20"\n' % (self.track_id, self.track_id))
        data_save = getattr(self, db_save)
        chr_names = sorted(data_save.keys(), key=mix_sort)
        for chr in chr_names:
            chr_data = data_save[chr]
            for strand, strand_data in chr_data.items():
                positions = [(pos, raw) for pos, raw in strand_data.items()]
                positions.sort()
                for (pos, val) in positions:
                    strand_str = "" if strand=="+" else "-"
                    if db_save=="support":
                        val = len(val)
                    raw = self.get_value(chr, strand, pos, db="raw")
                    cpm = self.get_value(chr, strand, pos, db="cpm")
                    support = len(self.get_value(chr, strand, pos, db="support"))
                    meta = self.get_value(chr, strand, pos, db="meta")
                    if raw>=min_raw and cpm>=min_cpm and support>=min_support:
                        chr_str = chr
                        if genome!="ensembl":
                            chr_str = pybio.genomes.chr_ensembl_ucsc.get(genome, {}).get(chr_str, None) # always convert ucsc to ensembl
                            if chr_str==None:
                                not_converted += 1
                                continue
                        if filetype=="bed":
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

    def set_value(self, chr, strand, pos, val, db="raw"):
        if val==0:
            return
        data = getattr(self, db)
        if db in ["raw", "cpm", "meta"]:
            data.setdefault(chr, {}).setdefault(strand, {}).setdefault(pos, val)
        if db in ["support"]:
            data.setdefault(chr, {}).setdefault(strand, {}).setdefault(pos, set(val))

    def add_value(self, chr, strand, pos, val, db="raw"):
        if val==0:
            return
        data = getattr(self, db)
        if db in ["raw", "cpm"]:
            data.setdefault(chr, {}).setdefault(strand, {}).setdefault(pos, 0)
            data[chr][strand][pos] += val
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
            data[chr][strand][pos][mt] = data[chr][strand][pos].get(mt, 0) + val

    def get_value(self, chr, strand, pos, db="raw"):
        data = getattr(self, db)
        if db in ["raw", "cpm"]:
            return data.get(chr, {}).get(strand, {}).get(pos, 0)
        if db in ["support"]:
            return data.get(chr, {}).get(strand, {}).get(pos, set())
        if db in ["meta"]:
            return data.get(chr, {}).get(strand, {}).get(pos, {})

    def set_region(self, chr, strand, start, stop, val, db="raw"):
        # only for raw and cpm
        data = getattr(self, db)
        for i in range(start, stop+1):
            data.setdefault(chr, {}).setdefault(strand, {})
            data[chr][strand][i] = val

    def get_region(self, chr, strand, start, stop, db="raw"):
        # only for raw and cpm
        data = getattr(self, db)
        region_sum = 0
        for i in range(start, stop+1):
            region_sum += data.get(chr, {}).get(strand, {}).get(i, 0)
        return region_sum

    def get_vector(self, chr, strand, start, stop, db="raw"):
        # only for raw and cpm
        data = getattr(self, db)
        vector = []
        for i in range(start, stop+1):
            vector.append(data.get(chr, {}).get(strand, {}).get(i, 0))
        return vector

    def cluster(self, region_up=150, region_down=150):
        # go down position list: highly expressed -> lowly expressed
        # sum [region_up, region_down] and assign to main position
        for chr, strand_data in self.raw.items():
            for strand, pos_data in strand_data.items():
                print "clustering : %s %s" % (chr, strand)
                offset_up, offset_down = region_up, region_down
                if strand=="-":
                    offset_up, offset_down = offset_down, offset_up
                mp = 1 if strand=="+" else -1
                positions = [(raw, mp*pos, pos) for pos, raw in pos_data.items()]
                positions = sorted(positions, reverse=True)
                self.cpm.setdefault(chr, {}).setdefault(strand, {})
                self.support.setdefault(chr, {}).setdefault(strand, {})
                self.meta.setdefault(chr, {}).setdefault(strand, {})
                temp_raw = {}
                temp_cpm = {}
                temp_support = {}
                temp_meta = {}
                for _, _, main_pos in positions:
                    if self.get_value(chr, strand, main_pos, db="raw")==0: # has this position been deleted (added or ignored)? move on
                        continue
                    new_raw = 0
                    new_cpm = 0
                    new_support = set()
                    new_meta = {}
                    for i in range(main_pos-offset_up, main_pos+offset_down+1):
                        raw = self.get_value(chr, strand, i, db="raw")
                        if raw!=0:
                            new_raw += raw
                            del self.raw[chr][strand][i]
                        cpm = self.get_value(chr, strand, i, db="cpm")
                        if cpm!=0:
                            new_cpm += cpm
                            del self.cpm[chr][strand][i]
                        support = self.get_value(chr, strand, i, db="support")
                        if support!=set():
                            new_support = new_support.union(support)
                            del self.support[chr][strand][i]
                        meta = self.get_value(chr, strand, i, db="meta")
                        if meta!={}:
                            for mt, mv in meta.items():
                                new_meta[mt] = new_meta.get(mt, 0) + mv
                            del self.meta[chr][strand][i]
                    temp_raw[main_pos] = new_raw
                    temp_cpm[main_pos] = new_cpm
                    temp_support[main_pos] = new_support
                    temp_meta[main_pos] = new_meta

                self.raw[chr][strand] = temp_raw
                self.cpm[chr][strand] = temp_cpm
                self.support[chr][strand] = temp_support
                self.meta[chr][strand] = temp_meta

    def filter(self, min_distance=125):
        # filter on raw data
        # sort positions and go down the list
        # only keep positions that are min_distance apart (don't sum anything)
        for chr, strand_data in self.raw.items():
            for strand, pos_data in strand_data.items():
                print "filtering: %s%s (MD = %snt)" % (strand, chr, min_distance)
                # if same value, take downstream position first
                mp = 1 if strand=="+" else -1
                positions = [(raw, mp*pos, pos) for pos, raw in pos_data.items()]
                positions = sorted(positions, reverse=True)
                temp_raw = {}
                temp_cpm = {}
                temp_support = {}
                temp_meta = {}
                for _, _, main_pos in positions:
                    if self.get_value(chr, strand, main_pos, db="raw")==0: # has this position been deleted (added or ignored)? move on
                        continue
                    temp_raw[main_pos] = self.get_value(chr, strand, main_pos, db="raw")
                    temp_cpm[main_pos] = self.get_value(chr, strand, main_pos, db="cpm")
                    temp_support[main_pos] = self.get_value(chr, strand, main_pos, db="support")
                    temp_meta[main_pos] = self.get_value(chr, strand, main_pos, db="meta")
                    # delete main_pos and surroundings [-min_distance, min_distance]
                    for i in range(main_pos-min_distance, main_pos+min_distance+1):
                        if self.get_value(chr, strand, i, db="raw")!=0:
                            del self.raw[chr][strand][i]
                        if self.get_value(chr, strand, i, db="cpm")!=0:
                            del self.cpm[chr][strand][i]
                        if self.get_value(chr, strand, i, db="support")!=set():
                            del self.support[chr][strand][i]
                        if self.get_value(chr, strand, i, db="meta")!={}:
                            del self.meta[chr][strand][i]
                self.raw[chr][strand] = temp_raw
                self.cpm[chr][strand] = temp_cpm
                self.support[chr][strand] = temp_support
                self.meta[chr][strand] = temp_meta
