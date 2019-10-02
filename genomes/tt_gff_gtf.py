header = []
header.append("Gene stable ID")
header.append("Transcript stable ID")
header.append("Chromosome/scaffold name")
header.append("Gene start (bp)")
header.append("Gene end (bp)")
header.append("Transcript start (bp)")
header.append("Transcript end (bp)")
header.append("Strand")
header.append("5' UTR start")
header.append("5' UTR end")
header.append("3' UTR start")
header.append("3' UTR end")
header.append("Gene type")
header.append("Exon region start (bp)")
header.append("Exon region end (bp)")
header.append("Constitutive exon")
header.append("Exon rank in transcript")
header.append("Exon stable ID")
header.append("Gene name")

geneinfo = {}

f = open("tt.annotation.custom/tt.gff")
r = f.readline()
while r:
    r = r.replace("\n", "").replace("\r", "").split("\t")
    rtype = r[2]
    att = {}
    temp = r[-1].split(";")
    for rec in temp:
        aname, aval = rec.split("=")
        att[aname] = aval
    if rtype=="gene":
        gene_id = att["ID"]
    if rtype=="exon":
        exon_start = int(r[3])
        exon_stop = int(r[4])
        gene_start, gene_stop = geneinfo.get(gene_id, [(), 0])
        gene_start = min(gene_start, exon_start)
        gene_stop = max(gene_stop, exon_stop)
        geneinfo[gene_id] = [gene_start, gene_stop]
    r = f.readline()
f.close()

fout = open("tt.annotation.custom/tt.annotation.custom.tab", "wt")
fout.write("\t".join(header) + "\n")
fout2 = open("tt.annotation.custom/tt.gtf", "wt")
gene_id = None
f = open("tt.annotation.custom/tt.gff")
r = f.readline()
while r:
    r = r.replace("\n", "").replace("\r", "").split("\t")
    rtype = r[2]
    att = {}
    temp = r[-1].split(";")
    for rec in temp:
        aname, aval = rec.split("=")
        att[aname] = aval
    if rtype=="gene":
        gene_id = att["ID"]
        gene_name = att["Name"]
        exon_start = None
        exon_stop = None
        #print gene_id, gene_name
    if rtype=="exon":
        chr = r[0]
        strand = r[6]
        exon_start = r[3]
        exon_stop = r[4]
        gene_start, gene_stop = geneinfo[gene_id]
        row_ensembl = [gene_id, "", chr, gene_start, gene_stop, "", "", {"+":1, "-":-1}[strand], "", "", "", "", "protein_coding", exon_start, exon_stop, 1, 1, "", gene_name]
        fout.write("\t".join(str(e) for e in row_ensembl) + "\n")
        row_gtf = [r[0], "", "exon", exon_start, exon_stop, ".", r[6], ".", "gene_id \"%s\"; gene_name \"%s\"" % (gene_id, gene_name)]
        fout2.write("\t".join(str(e) for e in row_gtf) + "\n")
    r = f.readline()
f.close()
fout.close()
