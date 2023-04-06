import gzip
import bz2
import os

class Fasta:

    """
    FASTA file class
    """

    def __init__(self, fname):
        self.fname = fname
        self.open_file()

    def open_file(self):
        self.id = None
        self.next_id = None
        if self.fname.endswith(".gz"):
            self.f = gzip.open(self.fname, "rt")
        elif self.fname.endswith(".bz2"):
            self.f = bz2.BZ2File(file_name, "r")
        else:
            self.f = open(self.fname, "rt")

    def read(self):
        """
        Reads next ID and SEQUENCE from FASTA file
        """
        self.sequence = []
        if self.next_id!=None:
            self.id = self.next_id
            self.next_id = None
        r = self.f.readline()
        if not r:
            return False
        while r:
            r = r.replace("\r", "").replace("\n", "")
            if r=="":
                r = self.f.readline()
                continue
            if r[0]==">" and self.id==None:
                self.id = r[1:]
                r = self.f.readline()
                continue
            elif r[0]==">":
                self.next_id = r[1:]
                self.sequence = "".join(self.sequence)
                return True
            self.sequence.append(r)
            r = self.f.readline()
        self.sequence = "".join(self.sequence)
        return True

    def split(self):
        """
        Splits FASTA file
        """
        while self.read():
            folder_name = os.path.dirname(self.fname)
            seq_id = self.id.split(" ")[0] # take everything until first space
            print(f"[pybio] splitting fasta file {self.fname}, sequence {seq_id}")
            if folder_name!="":
                fout = open(f"{folder_name}/%s.string" % seq_id, "wt")
            else:
                fout = open("%s.string" % seq_id, "wt")
            fout.write(self.sequence)
            fout.close()
