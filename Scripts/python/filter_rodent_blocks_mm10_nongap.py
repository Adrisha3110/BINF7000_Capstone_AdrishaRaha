#!/usr/bin/env python3

import os

input_dir = "Results/chrY_chunks_padded_rodents"
output_dir = "Results/chrY_chunks_padded_rodents_mm10ref"
report_file = "Results/mm10_block_filter_summary.txt"
ref_species = "mm10"

os.makedirs(output_dir, exist_ok=True)

kept = 0
removed = 0

for fname in sorted(os.listdir(input_dir)):
    if not fname.endswith(".fa"):
        continue

    in_path = os.path.join(input_dir, fname)
    out_path = os.path.join(output_dir, fname)

    found_mm10 = False
    mm10_seq = None

    with open(in_path) as fin:
        lines = fin.readlines()

    i = 0
    while i < len(lines):
        if lines[i].startswith(">"):
            header = lines[i][1:].strip()
            species = header.split("-")[0]
            seq = ""
            i += 1
            while i < len(lines) and not lines[i].startswith(">"):
                seq += lines[i].strip()
                i += 1
            if species == ref_species:
                found_mm10 = True
                mm10_seq = seq
        else:
            i += 1

    if found_mm10 and mm10_seq and set(mm10_seq) != {"-"}:
        with open(out_path, "w") as fout:
            fout.writelines(lines)
        kept += 1
    else:
        removed += 1

with open(report_file, "w") as rep:
    rep.write(f"Reference species: {ref_species}\n")
    rep.write(f"Kept blocks (mm10 has at least one non-gap base): {kept}\n")
    rep.write(f"Removed blocks (mm10 missing or all gaps): {removed}\n")

print(f"Kept: {kept}")
print(f"Removed: {removed}")
print(f"Report written to: {report_file}")
