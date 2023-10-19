"""
pybio

A compact python module developed to handle common bioinformatics file formats, especially
in the next-generation sequencing (NGS) field.

Github: https://github.com/grexor/pybio
"""

import os
import json
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
    print(f"[pybio genome_download] species {species} and version {genome_version}")
    genomes_ready_fname = os.path.join(pybio.config.genomes_folder, "genomes_ready.json")
    if os.path.exists(genomes_ready_fname):
        genomes_ready = json.load(open(genomes_ready_fname, "rt"))
    else:
        genomes_ready = {}
    genomes_ready.setdefault(species, {}).setdefault(genome_version, {}).setdefault("assembly", False)
    genomes_ready.setdefault(species, {}).setdefault(genome_version, {}).setdefault("annotation", False)
    assembly_folder = os.path.join(pybio.config.genomes_folder, f"{species}.assembly.{genome_version}")
    if not genomes_ready.get(species, {}).get(genome_version, {}).get("assembly", False) or not os.path.exists(assembly_folder):
        return_code = pybio.core.genomes.download_assembly(species, genome_version)
        if return_code==0:
            genomes_ready[species][genome_version]["assembly"] = True
        else:
            genomes_ready[species][genome_version]["assembly"] = False
    else:
        print(f"[pybio genome] FASTA ready at {pybio.config.genomes_folder}/{species}.assembly.{genome_version}")

    annotation_folder = os.path.join(pybio.config.genomes_folder, f"{species}.annotation.{genome_version}")
    if not genomes_ready.get(species, {}).get(genome_version, {}).get("annotation", False) or not os.path.exists(annotation_folder):
        return_code = pybio.core.genomes.download_annotation(species, genome_version)
        if return_code==0:
            return_code = pybio.core.genomes.prepare(species, genome_version)
        if return_code==0:
            genomes_ready[species][genome_version]["annotation"] = True
        else:
            genomes_ready[species][genome_version]["annotation"] = False
    else:
        print(f"[pybio genome] genome annotation at {pybio.config.genomes_folder}/{species}.annotation.{genome_version}")

    json.dump(genomes_ready, open(genomes_ready_fname, "wt"))

def genome_import(species, genome_version, args):
    print(f"[pybio genome] species {species} and version {genome_version}")
    genomes_ready_fname = os.path.join(pybio.config.genomes_folder, "genomes_ready.json")
    if os.path.exists(genomes_ready_fname):
        genomes_ready = json.load(open(genomes_ready_fname, "rt"))
    else:
        genomes_ready = {}
    genomes_ready.setdefault(species, {}).setdefault(genome_version, {}).setdefault("assembly", False)
    genomes_ready.setdefault(species, {}).setdefault(genome_version, {}).setdefault("annotation", False)
    assembly_folder = os.path.join(pybio.config.genomes_folder, f"{species}.assembly.{genome_version}")
    if not genomes_ready.get(species, {}).get(genome_version, {}).get("assembly", False) or not os.path.exists(assembly_folder):
        fasta_fname = args.fasta
        print(f"[pybio genome] importing {fasta_fname} to {assembly_folder}/{species}.fasta")
        if fasta_fname.endswith(".gz"):
            os.system(f"mkdir {assembly_folder}; cp {fasta_fname} {assembly_folder}/{species}.fasta.gz; gunzip -f {assembly_folder}/{species}.fasta.gz")
        else:
            os.system(f"mkdir {assembly_folder}; cp {fasta_fname} {assembly_folder}/{species}.fasta")
        return_code = os.system("python3 -c \"import pybio; pybio.data.Fasta('{assembly_folder}/{species}.fasta').split()\"".format(assembly_folder=assembly_folder, species=species))
        if return_code==0:
            genomes_ready[species][genome_version]["assembly"] = True
        else:
            genomes_ready[species][genome_version]["assembly"] = False
    else:
        print(f"[pybio genome] FASTA ready at {pybio.config.genomes_folder}/{species}.assembly.{genome_version}")
    json.dump(genomes_ready, open(genomes_ready_fname, "wt"))

    annotation_folder = os.path.join(pybio.config.genomes_folder, f"{species}.annotation.{genome_version}")
    if not genomes_ready.get(species, {}).get(genome_version, {}).get("annotation", False) or not os.path.exists(annotation_folder):
        gtf_fname = args.gtf
        print(f"[pybio genome] importing {gtf_fname} to {annotation_folder}/{species}.gtf.gz")
        if gtf_fname.endswith(".gz"):
            os.system(f"mkdir {annotation_folder}; cp {gtf_fname} {annotation_folder}/{species}.gtf.gz; gunzip -f {annotation_folder}/{species}.gtf.gz")
        else:
            os.system(f"mkdir {annotation_folder}; cp {gtf_fname} {annotation_folder}/{species}.gtf; gzip -f {annotation_folder}/{species}.gtf")
        return_code = pybio.core.genomes.prepare(species, genome_version)
        if return_code==0:
            genomes_ready[species][genome_version]["annotation"] = True
        else:
            genomes_ready[species][genome_version]["annotation"] = False
            genomes_ready[species][genome_version]["STAR"] = False
            genomes_ready[species][genome_version]["salmon"] = False
    else:
        print(f"[pybio genome] genome annotation at {pybio.config.genomes_folder}/{species}.annotation.{genome_version}")
    json.dump(genomes_ready, open(genomes_ready_fname, "wt"))

def genome_prepare(species, genome_version, args):
    genomes_ready_fname = os.path.join(pybio.config.genomes_folder, "genomes_ready.json")
    if os.path.exists(genomes_ready_fname):
        genomes_ready = json.load(open(genomes_ready_fname, "rt"))
    else:
        genomes_ready = {}
    genomes_ready.setdefault(species, {}).setdefault(genome_version, {}).setdefault("assembly", False)
    genomes_ready.setdefault(species, {}).setdefault(genome_version, {}).setdefault("annotation", False)
    if pybio.utils.is_tool("STAR") and not args.nostar:
        star_folder = os.path.join(pybio.config.genomes_folder, f"{species}.assembly.{genome_version}.star")
        if not genomes_ready.get(species, {}).get(genome_version, {}).get("STAR", False) or not os.path.exists(star_folder):
            return_code = pybio.core.genomes.star_index(species, genome_version, threads=args.threads)
            if return_code==0:
                genomes_ready[species][genome_version]["STAR"] = True
            else:
                genomes_ready[species][genome_version]["STAR"] = False
        else:
            print(f"[pybio genome] STAR index ready at {pybio.config.genomes_folder}/{species}.assembly.{genome_version}.star")

    if pybio.utils.is_tool("salmon") and not args.nosalmon:
        salmon_folder = os.path.join(pybio.config.genomes_folder, f"{species}.transcripts.{genome_version}.salmon")
        if not genomes_ready.get(species, {}).get(genome_version, {}).get("salmon", False) or not os.path.exists(salmon_folder):
            return_code = pybio.core.genomes.salmon_index(species, genome_version)
            if return_code==0:
                genomes_ready[species][genome_version]["salmon"] = True
            else:
                genomes_ready[species][genome_version]["salmon"] = False
        else:
            print(f"[pybio genome] salmon index ready at {pybio.config.genomes_folder}/{species}.transcripts.{genome_version}.salmon")

    json.dump(genomes_ready, open(genomes_ready_fname, "wt"))
