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
import pybio.map
import pybio.utils
import pybio.expression
import pybio.genomes
import pybio.maths
import pybio.config
import pybio.sequence
import pybio.core

# initialize path module
pybio.config.init()
pybio.path.init()
pybio.core.genomes.init()

def genome_download_prepare(species, ensembl_version, args):
    print(f"[pybio ensembl] species {species} and Ensembl version {ensembl_version}")

    genomes_ready_fname = os.path.join(pybio.config.genomes_folder, "genomes_ready.json")
    if os.path.exists(genomes_ready_fname):
        genomes_ready = json.load(open(genomes_ready_fname, "rt"))
    else:
        genomes_ready = {}

    assembly_folder = os.path.join(pybio.config.genomes_folder, f"{species}.assembly.ensembl{ensembl_version}")
    if not genomes_ready.get(species, {}).get(ensembl_version, {}).get("assembly", False) or not os.path.exists(assembly_folder):
        pybio.core.genomes.download_assembly(species, ensembl_version)
        genomes_ready.setdefault(species, {}).setdefault(ensembl_version, {}).setdefault("assembly", True)
    else:
        print(f"[pybio ensembl] FASTA ready at {pybio.config.genomes_folder}/{species}.assembly.ensembl{ensembl_version}")

    annotation_folder = os.path.join(pybio.config.genomes_folder, f"{species}.annotation.ensembl{ensembl_version}")
    if not genomes_ready.get(species, {}).get(ensembl_version, {}).get("annotation", False) or not os.path.exists(annotation_folder):
        pybio.core.genomes.download_annotation(species, ensembl_version)
        pybio.core.genomes.prepare(species, ensembl_version)
        genomes_ready.setdefault(species, {}).setdefault(ensembl_version, {}).setdefault("annotation", True)
    else:
        print(f"[pybio ensembl] genome annotation at {pybio.config.genomes_folder}/{species}.annotation.ensembl{ensembl_version}")

    if pybio.utils.is_tool("STAR") and not args.nostar:
        star_folder = os.path.join(pybio.config.genomes_folder, f"{species}.assembly.ensembl{ensembl_version}.star")
        if not genomes_ready.get(species, {}).get(ensembl_version, {}).get("STAR", False) or not os.path.exists(star_folder):
            pybio.core.genomes.star_index(species, ensembl_version)
            genomes_ready.setdefault(species, {}).setdefault(ensembl_version, {}).setdefault("STAR", True)
        else:
            print(f"[pybio ensembl] STAR index ready at {pybio.config.genomes_folder}/{species}.annotation.ensembl{ensembl_version}.star")

    if pybio.utils.is_tool("salmon") and not args.nosalmon:
        salmon_folder = os.path.join(pybio.config.genomes_folder, f"{species}.transcripts.ensembl{ensembl_version}.salmon")
        if not genomes_ready.get(species, {}).get(ensembl_version, {}).get("salmon", False) or not os.path.exists(salmon_folder):
            pybio.core.genomes.salmon_index(species, ensembl_version)
            genomes_ready.setdefault(species, {}).setdefault(ensembl_version, {}).setdefault("salmon", True)
        else:
            print(f"[pybio ensembl] salmon index ready at {pybio.config.genomes_folder}/{species}.annotation.ensembl{ensembl_version}.salmon")

    json.dump(genomes_ready, open(genomes_ready_fname, "wt"))