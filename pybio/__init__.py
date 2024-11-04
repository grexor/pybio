"""
pybio

A compact python module developed to handle common bioinformatics file formats, especially
in the next-generation sequencing (NGS) field.

Github: https://github.com/grexor/pybio
"""

import os
import json
import gzip
import pybio.path
import pybio.data
import pybio.utils
import pybio.maths
import pybio.config
import pybio.sequence
import pybio.core

pybio_path = os.path.abspath(__file__)
pybio_folder = os.path.dirname(pybio_path)
version = open(os.path.join(pybio_folder, "version"), "rt").readlines()[0].replace("\n", "").replace("\r", "")

# initialize path module
pybio.config.init()
pybio.path.init()
pybio.core.genomes.init()

def genome_download(species, genome_version, args):
    print(f"pybio | genome | download species {species} and version {genome_version}")
    pybio.core.genomes.genomes_present.setdefault(species, {}).setdefault(genome_version, {}).setdefault("assembly", False)
    pybio.core.genomes.genomes_present.setdefault(species, {}).setdefault(genome_version, {}).setdefault("annotation", False)
    assembly_folder = os.path.join(pybio.config.genomes_folder, f"{species}.assembly.{genome_version}")
    if not pybio.core.genomes.genomes_present.get(species, {}).get(genome_version, {}).get("assembly", False) or not os.path.exists(assembly_folder):
        return_code = pybio.core.genomes.download_assembly(species, genome_version)
        if return_code==0:
            pybio.core.genomes.genomes_present[species][genome_version]["assembly"] = True
            # write gingo file
            ginfo_fname = "{gdir}/{species}.assembly.{genome_version}/genome.info".format(gdir=pybio.config.genomes_folder, species=species, genome_version=genome_version)
            pybio.core.genomes.write_ginfo_file(ginfo_fname, {species: {genome_version: {"assembly":True}}})
        else:
            pybio.core.genomes.genomes_present[species][genome_version]["assembly"] = False
    else:
        print(f"pybio | genome | FASTA ready at {pybio.config.genomes_folder}/{species}.assembly.{genome_version}")

    annotation_folder = os.path.join(pybio.config.genomes_folder, f"{species}.annotation.{genome_version}")
    if not pybio.core.genomes.genomes_present.get(species, {}).get(genome_version, {}).get("annotation", False) or not os.path.exists(annotation_folder):
        return_code = pybio.core.genomes.download_annotation(species, genome_version)
        if return_code==0:
            return_code = pybio.core.genomes.prepare(species, genome_version)
            # write gingo file
            ginfo_fname = "{gdir}/{species}.annotation.{genome_version}/genome.info".format(gdir=pybio.config.genomes_folder, species=species, genome_version=genome_version)
            pybio.core.genomes.write_ginfo_file(ginfo_fname, {species: {genome_version: {"annotation":True}}})
        if return_code==0:
            pybio.core.genomes.genomes_present[species][genome_version]["annotation"] = True
        else:
            pybio.core.genomes.genomes_present[species][genome_version]["annotation"] = False
    else:
        print(f"[pybio genome] genome annotation at {pybio.config.genomes_folder}/{species}.annotation.{genome_version}")

    json.dump(pybio.core.genomes.genomes_present, open(os.path.join(pybio.config.genomes_folder, "genomes_ready.json"), "wt"), indent=4)

def genome_import(species, genome_version, args):
    print(f"pybio | genome | species {species} and version {genome_version}")
    pybio.core.genomes.genomes_present.setdefault(species, {}).setdefault(genome_version, {}).setdefault("assembly", False)
    pybio.core.genomes.genomes_present.setdefault(species, {}).setdefault(genome_version, {}).setdefault("annotation", False)
    assembly_folder = os.path.join(pybio.config.genomes_folder, f"{species}.assembly.{genome_version}")
    if not pybio.core.genomes.genomes_present.get(species, {}).get(genome_version, {}).get("assembly", False) or not os.path.exists(assembly_folder):
        fasta_fname = args.fasta
        print(f"pybio | genome | importing {fasta_fname} to {assembly_folder}/{species}.fasta")
        if fasta_fname.endswith(".gz"):
            os.system(f"mkdir {assembly_folder}; cp {fasta_fname} {assembly_folder}/{species}.fasta.gz; gunzip -f {assembly_folder}/{species}.fasta.gz")
        else:
            os.system(f"mkdir {assembly_folder}; cp {fasta_fname} {assembly_folder}/{species}.fasta")
        return_code = os.system("python3 -c \"import pybio; pybio.data.Fasta('{assembly_folder}/{species}.fasta').split()\"".format(assembly_folder=assembly_folder, species=species))
        if return_code==0:
            pybio.core.genomes.genomes_present[species][genome_version]["assembly"] = True
            # write gingo file
            ginfo_fname = "{gdir}/{species}.assembly.{genome_version}/genome.info".format(gdir=pybio.config.genomes_folder, species=species, genome_version=genome_version)
            pybio.core.genomes.write_ginfo_file(ginfo_fname, {species: {genome_version: {"assembly":True}}})
        else:
            pybio.core.genomes.genomes_present[species][genome_version]["assembly"] = False
    else:
        print(f"pybio | genome | FASTA ready at {pybio.config.genomes_folder}/{species}.assembly.{genome_version}")
    json.dump(pybio.core.genomes.genomes_present, open(os.path.join(pybio.config.genomes_folder, "genomes_ready.json"), "wt"), indent=4)

    annotation_folder = os.path.join(pybio.config.genomes_folder, f"{species}.annotation.{genome_version}")
    if not pybio.core.genomes.genomes_present.get(species, {}).get(genome_version, {}).get("annotation", False) or not os.path.exists(annotation_folder):
        gtf_fname = args.gtf
        print(f"[pybio genome] importing {gtf_fname} to {annotation_folder}/{species}.gtf.gz")
        if gtf_fname.endswith(".gz"):
            os.system(f"mkdir {annotation_folder}; cp {gtf_fname} {annotation_folder}/{species}.gtf.gz; gunzip -f {annotation_folder}/{species}.gtf.gz")
        else:
            os.system(f"mkdir {annotation_folder}; cp {gtf_fname} {annotation_folder}/{species}.gtf; gzip -f {annotation_folder}/{species}.gtf")
        return_code = pybio.core.genomes.prepare(species, genome_version)
        if return_code==0:
            pybio.core.genomes.genomes_present[species][genome_version]["annotation"] = True
            # write gingo file
            ginfo_fname = "{gdir}/{species}.annotation.{genome_version}/genome.info".format(gdir=pybio.config.genomes_folder, species=species, genome_version=genome_version)
            pybio.core.genomes.write_ginfo_file(ginfo_fname, {species: {genome_version: {"annotation":True}}})
        else:
            pybio.core.genomes.genomes_present[species][genome_version]["annotation"] = False
            pybio.core.genomes.genomes_present[species][genome_version]["STAR"] = False
            pybio.core.genomes.genomes_present[species][genome_version]["salmon"] = False
    else:
        print(f"pybio | genome | annotation at {pybio.config.genomes_folder}/{species}.annotation.{genome_version}")
    json.dump(genomes_ready, open(os.path.join(pybio.config.genomes_folder, "genomes_ready.json"), "wt"), indent=4)

def genome_prepare(species, genome_version, args):
    pybio.core.genomes.genomes_present.setdefault(species, {}).setdefault(genome_version, {}).setdefault("assembly", False)
    pybio.core.genomes.genomes_present.setdefault(species, {}).setdefault(genome_version, {}).setdefault("annotation", False)
    if pybio.utils.is_tool("STAR") and not args.nostar:
        star_folder = os.path.join(pybio.config.genomes_folder, f"{species}.assembly.{genome_version}.star")
        if not pybio.core.genomes.genomes_present.get(species, {}).get(genome_version, {}).get("STAR", False) or not os.path.exists(star_folder):
            return_code = pybio.core.genomes.star_index(species, genome_version, threads=args.threads)
            if return_code==0:
                pybio.core.genomes.genomes_present[species][genome_version]["STAR"] = True
                # write gingo file
                ginfo_fname = "{gdir}/{species}.assembly.{genome_version}.star/genome.info".format(gdir=pybio.config.genomes_folder, species=species, genome_version=genome_version)
                pybio.core.genomes.write_ginfo_file(ginfo_fname, {species: {genome_version: {"STAR":True}}})
            else:
                pybio.core.genomes.genomes_present[species][genome_version]["STAR"] = False
        else:
            print(f"pybio | genome | STAR index ready at {pybio.config.genomes_folder}/{species}.assembly.{genome_version}.star")

    if pybio.utils.is_tool("salmon") and not args.nosalmon:
        salmon_folder = os.path.join(pybio.config.genomes_folder, f"{species}.transcripts.{genome_version}.salmon")
        if not pybio.core.genomes.genomes_present.get(species, {}).get(genome_version, {}).get("salmon", False) or not os.path.exists(salmon_folder):
            return_code = pybio.core.genomes.salmon_index(species, genome_version)
            if return_code==0:
                pybio.core.genomes.genomes_present[species][genome_version]["salmon"] = True
                # write gingo file
                ginfo_fname = "{gdir}/{species}.transcripts.{genome_version}/genome.info".format(gdir=pybio.config.genomes_folder, species=species, genome_version=genome_version)
                pybio.core.genomes.write_ginfo_file(ginfo_fname, {species: {genome_version: {"transcripts":True}}})
                ginfo_fname = "{gdir}/{species}.transcripts.{genome_version}.salmon/genome.info".format(gdir=pybio.config.genomes_folder, species=species, genome_version=genome_version)
                pybio.core.genomes.write_ginfo_file(ginfo_fname, {species: {genome_version: {"salmon":True}}})
            else:
                pybio.core.genomes.genomes_present[species][genome_version]["salmon"] = False
        else:
            print(f"pybio | genome | salmon index ready at {pybio.config.genomes_folder}/{species}.transcripts.{genome_version}.salmon")

    json.dump(pybio.core.genomes.genomes_present, open(os.path.join(pybio.config.genomes_folder, "genomes_ready.json"), "wt"), indent=4)

def gff4jbrowse(fname_input, fname_output):
    print(f"pybio | gff3 for JBrowse2 | {fname_input} {fname_output}")
    fin = gzip.open(fname_input) if fname_input.endswith(".gz") else open(fname_input)
    fout = gzip.open(fname_output, "wt") if fname_input.endswith(".gz") else open(fname_output, "wt")
    r = fin.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "").split("\t")
        if len(r)>2:
            if r[2]!="gene":
                atts = r[-1]
                atts = atts.split(";")
                new_atts = []
                gene_id = None
                for att in atts:
                    if att.find("Parent=gene:")!=-1:
                        gene_id = att.split("Parent=gene:")[1]
                    elif att.startswith("Name="):
                        if gene_id!=None:
                            new_atts.append(f"{att},{gene_id}")
                        else:
                            new_atts.append(att)
                    else:
                        new_atts.append(att)
                new_atts = ";".join(new_atts)
                r[-1] = new_atts
                fout.write("\t".join(r)+"\n")
        r = fin.readline()
    fin.close()
    fout.close()