#!/usr/bin/env python3

"""
Illumina reads demultiplexer
Author: Gregor Rot, https://grexor.github.io
"""

import os
import sys
import argparse
import pybio
import gzip
import itertools

def process(args):

    # exmaple for "123" and max_mutations = 1
    # {'1N3', '1C3', 'G23', '12T', '12C', '12G', '1A3', 'C23', '1T3', 'A23', '123', '1G3', '12N', 'N23', 'T23', '12A'}
    # note the above also includes "123" (not mutated)
    def mutate_barcodes(s, max_mutations, charset):
        s = s.upper()
        mutated_strings = set()
        charset = set(charset)
        for num_mutations in range(1, max_mutations + 1):
            for positions in itertools.combinations(range(len(s)), num_mutations):
                original_chars = [s[i] for i in positions]
                replacement_options = [charset - {s[i]} for i in positions]
                for replacements in itertools.product(*replacement_options):
                    mutated = list(s)
                    for pos, new_char in zip(positions, replacements):
                        mutated[pos] = new_char
                    mutated_strings.add(''.join(mutated))
        mutated_strings.add(s)
        return mutated_strings

    try:
        os.makedirs(args.output)
    except:
        pass

    barcodes = []
    for barcode in args.barcodes.split(","):
        barcode_name, barcode_file, barcode_code = barcode.split(":")
        barcodes.append((barcode_name, barcode_file, barcode_code))

    annotation = {}
    fsamples = {}
    stats = {}
    stop_at = float('inf')
    match_db = {}

    if not os.path.exists(args.annotation):
        print(f"Provided annotation file {args.annotation} not present, please check.")
        sys.exit(1)

    num_samples = 0
    with open(args.annotation, 'r') as f:
        num_samples = sum(1 for _ in f)

    f = open(args.annotation, 'rt')
    header = f.readline()
    header = header.replace('\r', '').replace('\n', '').split('\t')
    r = f.readline()
    samples_encountered = set()
    csample = 0
    while r:
        csample += 1
        r = r.replace('\r', '').replace('\n', '').split('\t')
        data = dict(zip(header, r))
        sample_name = data["sample_name"]
        print(f"adding barcodes for sample: {sample_name} ({csample}/{num_samples})", flush=True)

        # sample name must be unique
        assert(sample_name not in samples_encountered)
        samples_encountered.add(sample_name)
        
        b = []
        for (barcode_name, barcode_file, barcode_code) in barcodes:
            bcode = data[barcode_name]
            if barcode_code.find("R")!=-1: # is it reverse? (RRRR...)
                bcode = pybio.sequence.reverse_complement(bcode)
            b.append(bcode) 

        all_barcodes = []
        for index, barcode_value in enumerate(b):
            allow_mismatches = barcodes[index][2].split("_m") # RRRRRRRRRRRR_0_m1 for example
            if len(allow_mismatches)==2:
                allow_mismatches = int(allow_mismatches[-1])
            else:
                allow_mismatches = 0
            bm = mutate_barcodes(barcode_value, allow_mismatches, ["A", "T", "C", "G", "N"])
            all_barcodes.append(bm)

        combs = list(itertools.product(*all_barcodes))
        for el in combs:
            match_db.setdefault(el, set()).add(sample_name)

        fsamples[sample_name+"_R1"] = gzip.open(f"{args.output}/{sample_name}_R1.fastq.gz", "wt")
        fsamples[sample_name+"_R2"] = gzip.open(f"{args.output}/{sample_name}_R2.fastq.gz", "wt")      
        r = f.readline()
    f.close()

    fsamples["unknown_R1"] = gzip.open(f"{args.output}/unknown_R1.fastq.gz", "wt")
    fsamples["unknown_R2"] = gzip.open(f"{args.output}/unknown_R2.fastq.gz", "wt")

    fsamples["ambiguous_R1"] = gzip.open(f"{args.output}/ambiguous_R1.fastq.gz", "wt")
    fsamples["ambiguous_R2"] = gzip.open(f"{args.output}/ambiguous_R2.fastq.gz", "wt")

    def get_matches(i1, i2, r1, r2):
        match = True
        files = {"i1": i1, "i2": i2, "r1": r1, "r2": r2}
        k = []
        for (barcode_name, barcode_file, barcode_code) in barcodes:
            barcode_pos = int(barcode_code.split("_")[1])
            barcode_code = barcode_code.split("_")[0]
            code = files[barcode_file][barcode_pos:barcode_pos+len(barcode_code)]
            k.append(code)
        return list(match_db.get(tuple(k), set()))

    r1f = pybio.data.Fastq(args.r1)
    r2f = pybio.data.Fastq(args.r2)
    i1f = pybio.data.Fastq(args.i1)
    i2f = pybio.data.Fastq(args.i2)

    current = 0
    while i1f.read() and i2f.read() and r1f.read() and r2f.read():
        current += 1
        i1 = i1f.sequence
        i2 = i2f.sequence
        r1 = r1f.sequence
        r2 = r2f.sequence
        matches = get_matches(i1, i2, r1, r2)
        if len(matches)==1:
            sample_name = matches[0]
            fsamples[sample_name+"_R1"].write(f"{r1f.id}\n{r1f.sequence}\n+\n{r1f.quality}\n")
            fsamples[sample_name+"_R2"].write(f"{r2f.id}\n{r2f.sequence}\n+\n{r2f.quality}\n")
            stats[sample_name] = stats.get(sample_name, 0) + 1
        elif len(matches)==0:
            sample_name = "unknown"
            stats[sample_name] = stats.get(sample_name, 0) + 1
            fsamples[sample_name+"_R1"].write(f"{r1f.id}:{i1}\n{r1f.sequence}\n+\n{r1f.quality}\n")
            fsamples[sample_name+"_R2"].write(f"{r2f.id}:{i2}\n{r2f.sequence}\n+\n{r2f.quality}\n")
        elif len(matches)>1:
            sample_name = "ambiguous"
            stats[sample_name] = stats.get(sample_name, 0) + 1
            fsamples[sample_name+"_R1"].write(f"{r1f.id}:{i1}\n{r1f.sequence}\n+\n{r1f.quality}\n")
            fsamples[sample_name+"_R2"].write(f"{r2f.id}:{i2}\n{r2f.sequence}\n+\n{r2f.quality}\n")

        if current>stop_at:
            break
        if current%10000==0:
            print(f"processed {current/1000000}M reads", flush=True)

    for f in fsamples.values():
        f.close()

    stats = [(count, sample) for sample, count in stats.items()]
    stats.sort(reverse=True)

    f = open(args.stats, "wt")
    f.write("sample_name\tread_count\n")
    for (count, sample) in stats:
        f.write(f"{sample}\t{count}\n")
    f.close()
