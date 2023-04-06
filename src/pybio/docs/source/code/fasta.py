f = pybio.data.Fasta("data.fasta")
while f.read():
    print f.id
    print f.sequence
