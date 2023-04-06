class TabReader():

    """
    This class reads text data (.tab separated).
    """

    def __init__(self, filename="", header_at = 0):
        self.f = None
        self.header = None
        if filename!="":
            self.initialize(filename, header_at = header_at)
            
    def initialize(self, filename, header_at = 0):
        if filename.endswith(".gz"):
            self.f = gzip.open(filename)
        else:
            self.f = open(filename)
        while header_at>0:
            self.f.readline()
            header_at -= 1
        self.header = self.f.readline()
        self.header = self.header.replace("\"", "").replace("\r", "").replace("\n", "").split("\t")
    
    def readline(self):
        line = self.f.readline()
        if not line:
            return False
        line = line.replace("\r", "").replace("\n", "").split("\t")
        self.r = line
        self.data = dict(zip(self.header, line))
        return True
