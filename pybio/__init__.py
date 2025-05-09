"""
pybio

A compact python module developed to handle common bioinformatics file formats, especially
in the next-generation sequencing (NGS) field.

Github: https://github.com/grexor/pybio
"""

import os
import sys
import json
import gzip
import pybio.path
import pybio.data
import pybio.utils
import pybio.maths
import pybio.config
import pybio.sequence
import pybio.core
import argparse
import pybio.aimux

pybio_path = os.path.abspath(__file__)
pybio_folder = os.path.dirname(pybio_path)
version = open(os.path.join(pybio_folder, "version"), "rt").readlines()[0].replace("\n", "").replace("\r", "")

# initialize path module
pybio.config.init()
pybio.path.init()
pybio.core.genomes.init()

def genome_download(species, genome_version, args):
    print(f"pybio | genome | download/check | {species} | {genome_version}")
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
        print(f"pybio | genome | assembly | {pybio.config.genomes_folder}/{species}.assembly.{genome_version}")

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
        print(f"pybio | genome | annotation | {pybio.config.genomes_folder}/{species}.annotation.{genome_version}")

    #json.dump(pybio.core.genomes.genomes_present, open(os.path.join(pybio.config.genomes_folder, "genomes_ready.json"), "wt"), indent=4)

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

    annotation_folder = os.path.join(pybio.config.genomes_folder, f"{species}.annotation.{genome_version}")
    if not pybio.core.genomes.genomes_present.get(species, {}).get(genome_version, {}).get("annotation", False) or not os.path.exists(annotation_folder):
        gtf_fname = args.gtf
        print(f"pybio | genome | importing {gtf_fname} to {annotation_folder}/{species}.gtf.gz")
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
    #json.dump(genomes_ready, open(os.path.join(pybio.config.genomes_folder, "genomes_ready.json"), "wt"), indent=4)

def genome_prepare(species, genome_version, args, unknown_args=""):
    pybio.core.genomes.genomes_present.setdefault(species, {}).setdefault(genome_version, {}).setdefault("assembly", False)
    pybio.core.genomes.genomes_present.setdefault(species, {}).setdefault(genome_version, {}).setdefault("annotation", False)
    if pybio.utils.is_tool("STAR") and not args.nostar:
        star_folder = os.path.join(pybio.config.genomes_folder, f"{species}.assembly.{genome_version}.star")
        if not pybio.core.genomes.genomes_present.get(species, {}).get(genome_version, {}).get("STAR", False) or not os.path.exists(star_folder):
            return_code = pybio.core.genomes.star_index(species, genome_version, args, unknown_args=unknown_args)
            if return_code==0:
                pybio.core.genomes.genomes_present[species][genome_version]["STAR"] = True
                # write gingo file
                ginfo_fname = "{gdir}/{species}.assembly.{genome_version}.star/genome.info".format(gdir=pybio.config.genomes_folder, species=species, genome_version=genome_version)
                pybio.core.genomes.write_ginfo_file(ginfo_fname, {species: {genome_version: {"STAR":True}}})
            else:
                pybio.core.genomes.genomes_present[species][genome_version]["STAR"] = False
        else:
            print(f"pybio | genome | STAR index | {pybio.config.genomes_folder}/{species}.assembly.{genome_version}.star")

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
            print(f"pybio | genome | salmon index | {pybio.config.genomes_folder}/{species}.transcripts.{genome_version}.salmon")

    #json.dump(pybio.core.genomes.genomes_present, open(os.path.join(pybio.config.genomes_folder, "genomes_ready.json"), "wt"), indent=4)

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

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, add_help=False)
    parser.add_argument('commands', help="which command to run?", nargs='*')
    parser.add_argument("-version", "--version", help="Print version", action="store_true")
    parser.add_argument("-nostar", "--nostar", help="Do not process genome with STAR (index)", action="store_true")
    parser.add_argument("-nosalmon", "--nosalmon", help="Do not process transcripts with salmon (index)", action="store_true")
    parser.add_argument("-species", "--species", help="Species name, e.g. homo_sapiens")
    parser.add_argument("-genome_version", "--genome_version", help="Version of the genome, e.g. ensembl109; if not provided, the latest Ensembl version is used")
    parser.add_argument("-fasta", "--fasta", help="Process sequences from fasta file")
    parser.add_argument("-gtf", "--gtf", help="GTF file to import")
    parser.add_argument("-threads", "--threads", '-t', '--t', default="1", help="Number of threads to use (default: 1)")
    parser.add_argument("-help", "-h", "--help", action="store_true")
    parser.add_argument("-xs", "--xs", help="Add '--outSAMstrandField intronMotif'", action="store_false")
    parser.add_argument("-alignIntronMax", "--alignIntronMax", help="STAR alignIntronMax", type=int, default=None)
    parser.add_argument("-genomeSAindexNbases", "--genomeSAindexNbases", help="STAR genomeSAindexNbases")
    parser.add_argument("-genomeChrBinNbits", "--genomeChrBinNbits", help="STAR genomeChrBinNbits")

    args, unknown_args = parser.parse_known_args()
    unknown_args = " ".join(unknown_args)

    help_0 = """
    Usage: pybio <species|command> [options]

    Commands:
    -- <species>             Directly provide the genome species (e.g. homo_sapiens) as the first parameter
                            Example: "pybio homo_sapiens"

    -- Genome
        genome species        Same as above, however use the command "genome"
                            Example: "pybio genome homo_sapiens"

    -- Search
        search search_text    Displays all available genomes matching search_text
                            Example: "pybio search homo_sapiens"

    -- Configuration
        config [folder]       If no folder provided, prompts for location of genomes data folder
        path species          Prints FASTA and GTF file path for provided species (+genome version)

    Options:

        -genome_version       Specify genome version (e.g. ensembl113)
        -nostar, -nosalmon    Avoid building STAR and samlmon index for the genome
        -threads n            Number of threads to use (STAR mapping, etc), default: 1
        -fasta fasta_file     Fasta file for custom genome
        -gtf gtf_file         GTF file for custom genome
        -version              Display pybio version

    """

    help_genome = """
    To download genomes, please specify species.

    Example:

    $ pybio genome homo_sapiens

    The above will download the latest Ensembl human genome.
    """

    help_star = """
    Aligning reads to a reference genome with STAR:

    Example:

    # single-end sequencing
    $ pybio star homo_sapiens r1.fastq.gz output.bam

    # paired-end sequencing
    $ pybio star homo_sapiens r1.fastq.gz r2.fastq.gz output.bam
    """

    help_gff4jbrowse = """
    Prepare GFF3 for JBrowse2, without gene records and moving the Parent=gene: to Name property of transcripts

    Example:

    $ pybio gff3jbrowse2 input.gff output.gff
    """

    help_sam2bam = """
    Convert .sam to .bam file:

    Example:

    $ pybio sam2bam file1.sam file1.bam
    """

    print(f"pybio | v{pybio.version}, https://github.com/grexor/pybio")
    print(f"pybio | config file: {pybio.config.config_fname()}")
    print(f"pybio | genomes folder: {pybio.config.genomes_folder}")
    print()

    if args.version:
        sys.exit(0)

    if len(args.commands)==0 and args.help:
        print(help_0)
        sys.exit(0)

    def is_known_species(species):
        species = species.lower()
        return species in pybio.core.genomes.species_db.keys()

    def determine_species(species):
        provided_species = species
        id_species = None
        species = species.lower()
        potential_hits = []
        for species_id, species_data in pybio.core.genomes.species_db.items():
            display_name = species_data["display_name"].lower()
            if display_name.find(species)!=-1 or species_id.find(species)!=-1:
                potential_hits.append((len(display_name), species_id, display_name))
        potential_hits.sort()
        if len(potential_hits)>0:
            return potential_hits
        return []

    def display_potential_hits(potential_hits, search_text):
        print(f"pybio found {len(potential_hits)} genome hits for your provided genome species `{search_text}`.\n")
        for hit in potential_hits:
            print(f"Species = '\033[32;1m{hit[1]}\033[0m', display name = '{hit[2]}'")
        if len(potential_hits)>0:
            print(f"\nFor example, to download the first hit from the list above, you could write:\n")
            print(f"$ \033[32;1mpybio genome {potential_hits[0][1]}\033[0m")
            print()
        sys.exit(1)

    def resolve_species_version(species, args, downloaded_only=False):
        known_species = is_known_species(species)
        potential_hits = []
        if not known_species:
            potential_hits = determine_species(species)
            if len(potential_hits)==1:
                species = potential_hits[0][1]
        if args.genome_version!=None:
            genome_version = args.genome_version
            return species.lower(), genome_version, potential_hits
        if len(args.commands)>=3:
            genome_version = args.commands[2]
            return species.lower(), genome_version, potential_hits
        if len(args.commands)>=2:
            if args.commands[0] not in ["search", "genome", "genomes", "path", "species"]:
                genome_version = args.commands[1]
                return species.lower(), genome_version, potential_hits
        genome_version = determine_genome_version(species, downloaded_only=downloaded_only)
        return species.lower(), genome_version, potential_hits

    def determine_genome_version(species, downloaded_only):
        species = species.lower()
        potential_versions = []
        if downloaded_only:
            L = pybio.core.genomes.genomes_present.get(species, {})
            for k,v in L.items():
                if isinstance(v, dict):
                    ke = k.replace("ensemblgenomes", "").replace("ensembl","")
                    try:
                        ke = int(ke)
                        potential_versions = [k] + potential_versions
                    except:
                        potential_versions.append(k)
            if len(potential_versions)>0:
                return potential_versions[0]
            else:
                return None
        return pybio.core.genomes.species_db.get(species, {}).get("genome_version", None)

    def handle_genome(args, unknown_args=""):
        search_string = args.commands[1] if args.commands[0] in ["genome", "genomes", "search"] else args.commands[0]
        species, genome_version, potential_hits = resolve_species_version(search_string, args)
        if len(potential_hits)>1:
            display_potential_hits(potential_hits, species)
        if genome_version==None:
            print("Could not identify species name.\n")
            sys.exit(1)
        if genome_version.find("ensembl")!=-1:
            pybio.genome_download(species, genome_version, args)
            pybio.genome_prepare(species, genome_version, args, unknown_args=unknown_args)
        elif args.fasta!=None and args.gtf!=None: # user provided genome, species, fasta, gtf, genome_version
            species = args.commands[1].lower()
            args.nosalmon = True
            pybio.genome_import(species, genome_version, args)
            pybio.genome_prepare(species, genome_version, args)
        elif (args.fasta==None or args.gtf==None) and (genome_version.find("ensembl")==-1):
            print(f"Could not find a match for the species '{species}'.")
            print("If you would like to import a custom genome from your own FASTA and GTF files,\nan example call would be:\n")
            print("$ pybio genome custom_species -fasta /path/to/fasta -gtf /path/to/gtf -genome_version custom_genome_v1\n")
            print("The above imported genome would be reachable under species 'custom_species' and genome_version 'custom_genome_v1'.\n")
            print("In case you would like to import an Ensembl ready genome:\n")
            print("$ pybio genome homo_sapiens ensembl113\n")
            print("Ommiting the genome version will download the latest Ensembl release of the genome.\n")

    if len(args.commands)>0:

        if args.commands[0] == "aimux":
            sub_args = sys.argv[2:]
            parser = argparse.ArgumentParser(prog="pybio aimux")
            parser.add_argument('-r1', required=True, help='R1 FASTQ file')
            parser.add_argument('-r2', required=True, help='R2 FASTQ file')
            parser.add_argument('-i1', required=True, help='I1 FASTQ file')
            parser.add_argument('-i2', required=True, help='I2 FASTQ file')
            parser.add_argument('-stats', required=True, help='stats file name')
            parser.add_argument('-annotation', required=True, help='samples annotation file with barcodes, TAB delimited, 3 columns: sample_name (1), barcode_forward (2), barcode_reverse (3)')
            parser.add_argument('-barcodes', required=True, help='barcodes definition')
            parser.add_argument('-output', required=False, default="demulti", help='output folder')
            args = parser.parse_args(sub_args)
            pybio.aimux.process(args)

        if args.commands[0]=="path":
            if len(args.commands)<2:
                print("Please provide species, e.g.:\n\n$ pybio path homo_sapiens\n")
                sys.exit(0)
            species = args.commands[1]
            # genome_version can also be the second command line parameter
            # $ pybio path homo_sapiens genome_version
            # or can be specified as named parameter
            # $ pybio path homo_sapiens -genome_version refseq
            species, genome_version, potential_hits = resolve_species_version(species, args, downloaded_only=True) 
            annotation_folder = os.path.join(pybio.config.genomes_folder, "%s.annotation.%s" % (species, genome_version))
            assembly_folder = os.path.join(pybio.config.genomes_folder, "%s.assembly.%s" % (species, genome_version))
            fasta_file = os.path.join(assembly_folder, f"{species}.fasta")
            gtf_file = os.path.join(annotation_folder, f"{species}.gtf.gz")
            gff3_file = os.path.join(annotation_folder, f"{species}.gff3.gz")
            if not os.path.exists(fasta_file) or genome_version==None:
                print(f"pybio couldn't find any genomes for species '{species}' with genome version {genome_version}.\nPlease provide a valid genome species to retrieve paths.")
                sys.exit(0)
            print(fasta_file)
            print(gtf_file)
            print(gff3_file)
            sys.exit()

        if args.commands[0] in ["species", "search"]:
            if len(args.commands)>1:
                potential_hits = determine_species(args.commands[1])
                display_potential_hits(potential_hits, args.commands[1])
            else:
                print(f"Please provide a search term for the species, otherwise all available Ensembl genome species are stored in the local database at:\n\n{pybio.config.genomes_folder}/ensembl.json")
            sys.exit()

        if args.commands[0] in ["genome", "genomes"]:
            if len(args.commands)<2:
                print("Please provide a search term for the species")
                sys.exit()
            handle_genome(args, unknown_args=unknown_args)
            sys.exit()

        if args.commands[0]=="config":
            if len(args.commands)>1:
                pybio.config.init(args.commands[1]) # pybio config genome_folder
            sys.exit()

        if args.commands[0]=="sam2bam":
            if len(args.commands) not in [3]:
                print(help_sam2bam)
                sys.exit()
            input_fname = args.commands[1]
            output_fname = args.commands[2]
            os.system(f"{pybio.config.shell} -c 'samtools view -@ {args.threads} -F 4 -bS {input_fname} -o {output_fname} # -F 4 means only mapped reads'")
            os.system(f"{pybio.config.shell} -c 'samtools sort {output_fname} -o {output_fname} -@ {args.threads}'")
            os.system(f"{pybio.config.shell} -c 'samtools index {output_fname} -@ {args.threads}'")
            sys.exit()

        if args.commands[0]=="gff4jbrowse":
            if len(args.commands) not in [3]:
                print(help_gff4jbrowse)
                sys.exit()
            file1=args.commands[1]
            file2=args.commands[2]
            pybio.gff4jbrowse(file1, file2)
            sys.exit()

        if args.commands[0]=="aimux":
            if len(args.commands) not in [4,5]:
                print(help_star)
                sys.exit()       

        if args.commands[0]=="star":
            if len(args.commands) not in [4,5]:
                print(help_star)
                sys.exit()       

            species = args.commands[1].lower()

            if args.genome_version!=None:
                genome_version = args.genome_version
            else:
                genome_version = pybio.core.genomes.species_db.get(species, {}).get("genome_version", None)
            if genome_version==None:
                print("Could not determine genome version")
                sys.exit()
            star_folder = os.path.join(pybio.config.genomes_folder, f"{species}.assembly.{genome_version}.star")
            if len(args.commands)==4:
                file1=args.commands[2]
                file2=None
                output = args.commands[3]
            else:
                file1=args.commands[2]
                file2=args.commands[3]
                output = args.commands[4]

            if output.endswith(".bam"): # remove .bam if provided
                output = output[:-4]

            add_params = []
            if args.xs: # --outFilterMultimapNmax 1 --outSAMstrandField intronMotif
                add_params.append("--outSAMstrandField intronMotif")
            if args.alignIntronMax is not None:
                add_params.append(f"--alignIntronMax {args.alignIntronMax}")
            add_params = " ".join(add_params)
            if file2!=None: # paired-end
                os.system(f"{pybio.config.shell} -c 'STAR --runThreadN {args.threads} --outFileNamePrefix {output}_  --genomeDir {star_folder} --readFilesIn {file1} {file2} --readFilesCommand zcat {add_params} {unknown_args}'")
            if file2==None: # single-end
                os.system(f"{pybio.config.shell} -c 'STAR --runThreadN {args.threads} --outFileNamePrefix {output}_  --genomeDir {star_folder} --readFilesIn {file1} --readFilesCommand zcat {add_params} {unknown_args}'")

            os.system(f"{pybio.config.shell} -c 'pybio sam2bam {output}_Aligned.out.sam {output}.bam -threads {args.threads}'")
            os.system(f"{pybio.config.shell} -c 'mv {output}_Log.final.out {output}.stats.txt'")
            os.system(f"{pybio.config.shell} -c 'mv {output}_Log.out {output}.log.txt'")
            os.system(f"{pybio.config.shell} -c 'mv {output}_Log.progress.out {output}.progress.txt'")
            os.system(f"{pybio.config.shell} -c 'rm {output}_Aligned.out.sam'")
            os.system(f"{pybio.config.shell} -c 'mv {output}_SJ.out.tab {output}.splice.tab'")
            sys.exit()

        handle_genome(args) # lastly: is the first parameter maybe a genome / species name?    