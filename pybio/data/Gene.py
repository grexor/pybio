import pybio
from operator import attrgetter

"""
GTF format reader. Very simple for now.
"""

class Gene():

    def __init__(self, id = None, chr = None, strand = None, attrs = None, name = None, alias = None):
        self.id = id
        self.name = name
        self.chr = chr
        self.strand = strand
        self.attrs = attrs
        self.alias = alias
        self.features = []
        self.start = ()
        self.stop = 0

    def add_feature(self, feature_new, merge = True):
        """
        Add feature to Gene

        :param merge_overlapping: if the feature (exon, CDS) overlaps with an existing feature of the same type, simply adjust the start,stop positions
        """
        added = False
        # if feature_new is overlapping with an existing feature of the same type, adjust start,stop positions
        if merge:
            for feature in self.features:
                if feature.type!=feature_new.type:
                    continue
                if pybio.utils.interval_overlap(feature.start, feature.stop, feature_new.start, feature_new.stop)>0:
                    feature.start = min(feature.start, feature_new.start)
                    feature.stop = max(feature.stop, feature_new.stop)
                    added = True
                    break
        if not added:
            self.features.append(feature_new)

        # keep features sorted by start position
        self.features = sorted(self.features, key=attrgetter("start"))

        self.start = min(self.start, feature_new.start)
        self.stop = max(self.stop, feature_new.stop)
