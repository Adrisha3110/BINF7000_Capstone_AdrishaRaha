#!/usr/bin/env python3

import os

species_file = "Data/Tree/rodent_species_list.txt"
input_dir = "Results/chrY_chunks_padded"
output_dir = "Results/chrY_chunks_padded_rodents"

with open(species_file) as f:
    keep_species = set(line.strip() for line in f if line.strip())

os.makedirs(output_dir, exist_ok=True)

for fname in os.listdir(input_dir):
    if not fname.endswith(".fa"):
        continue

    in_path = os.path.join(input_dir, fname)
    out_path = os.path.join(output_dir, fname)

    with open(in_path) as fin, open(out_path, "w") as fout:
        write_block = False
        for line in fin:
            if line.startswith(">"):
                header = line[1:].strip()
                species = header.split("-")[0]
                write_block = species in keep_species
                if write_block:
                    fout.write(line)
            else:
                if write_block:
                    fout.write(line)
