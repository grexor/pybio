f = pybio.data.TabReader("data.tab")
while f.read():
    print f.Gene_name # string : Gene_name, note that "Gene name" from the tab file was renamed to "Gene_name" (space to _)
    print f.Alias     # string : Alias
    print f.Exp       # float : Exp
    print f.r         # list of all columns of current line
    print f.r[0]      # = f.Gene_name (first column)
