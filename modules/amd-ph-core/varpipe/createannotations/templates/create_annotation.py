#!/usr/bin/env python3

"""Parses vcf files and creates an annotation file.

Inputs:
  - snpEff annotated vcf
  - sample_id

"""

import csv
import platform
from itertools import zip_longest

import yaml

FORMAT_KEYS = ["GT", "AD", "AF", "DP", "F1R2", "F2R1", "FAD"]
OPTIONAL_KEYS = ["PGT", "PID", "PS"]
LAST_KEY = ["SB"]
FORMAT_SHORT_KEYS = list(FORMAT_KEYS + LAST_KEY)
FORMAT_EXTENDED_KEYS = list(FORMAT_KEYS + OPTIONAL_KEYS + LAST_KEY)


def parse_format_field(format_string):
    values = format_string.split(":")
    keys = FORMAT_SHORT_KEYS if len(values) == 8 else FORMAT_EXTENDED_KEYS
    format_fields = dict(zip_longest(keys, values, fillvalue=""))

    return format_fields


def parse_info_field(info_string):
    info_arr = info_string.split(";")
    info_fields = {}

    for item in info_arr:
        # Skip items without key=value format
        if "=" not in item:
            continue

        # Split into key and value
        key, value = item.split("=", 1)
        info_fields[key] = value

    return info_fields


def process_annotation(annotation, sample_id, output_file):
    with open(annotation) as infile:
        lines = [line.rstrip("\\r\\n").split("\\t") for line in infile if not line.startswith("#")]
    # Predefined fieldnames
    fieldnames = [
        "Sample ID",
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "FILTER",
        "GT",
        "AD",
        "AF",
        "DP",
        "F1R2",
        "F2R1",
        "FAD",
        "PGT",
        "PID",
        "PS",
        "SB",
        "AS_FilterStatus",
        "AS_SB_TABLE",
        "rawDP",  # Changed from "DP"
        "ECNT",
        "MBQ",
        "MFRL",
        "MMQ",
        "MPOS",
        "POPAF",
        "RPA",
        "RU",
        "TLOD",
        "ANN",
    ]

    # Parse the VCF data
    clean_lines = []
    for line in lines:
        if len(line) < 10:  # Skip malformed lines
            continue

        # Initialize fields with basic VCF information
        fields = {
            "Sample ID": sample_id,
            "CHROM": line[0],
            "POS": line[1],
            "REF": line[3],
            "ALT": line[4],
            "FILTER": line[6],
        }

        # Parse FORMAT and INFO fields
        format_fields = parse_format_field(line[9])
        info_fields = parse_info_field(line[7])

        # Add FORMAT fields to the result
        fields.update(format_fields)

        # Add INFO fields to the result
        fields.update(info_fields)

        # Get DP values and select the lower one
        format_dp = format_fields["DP"]
        info_dp = info_fields.get("DP", format_dp)
        fields["DP"] = str(min(int(format_dp), int(info_dp)))

        # Add rawDP from INFO
        fields["rawDP"] = info_fields.get("DP", "")

        # Ensure ANN is prefixed with 'ANN='
        if "ANN" in fields:
            fields["ANN"] = f"ANN={fields['ANN']}"

        # Filter out any fields not in predefined fieldnames
        filtered_fields = {k: fields.get(k, "") for k in fieldnames}

        # Add to clean lines
        clean_lines.append(filtered_fields)

    # Write output
    with open(output_file, "w") as tsvfile:
        writer = csv.DictWriter(tsvfile, fieldnames=fieldnames, delimiter="\\t")
        writer.writeheader()

        for line in clean_lines:
            writer.writerow(line)


if __name__ == "__main__":
    prefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else "$meta.id"
    annotation = "${annotation}"

    input_type = "DR_loci" if "DR_loci" in "${annotation}" else "full"
    output = f"{prefix}_{input_type}_annotation.txt"

    task_process = "${task.process}"

    # Process the annotation file
    process_annotation(annotation, prefix, output)

    versions = {
        task_process: {
            "python": platform.python_version(),
        }
    }

    with open("versions.yml", "w", encoding="utf-8") as f:
        f.write(yaml.dump(versions))
