#!/usr/bin/env python3

"""
Cross-references coverage depth information with a set of annotated
structural regions and analyzes for deletions within these sites
"""

import csv
import platform
from collections import defaultdict

import yaml


def parse_genome_stats(stats_file):
    """Parse genome region coverage stats into dictionary"""
    stats = {}
    with open(stats_file, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\\t")
        for row in reader:
            stats[row["Reference"], row["Gene Name"]] = row
    return stats


def parse_bed_file(bed_file):
    """Parse BED file of annotated regions into a dictionary"""
    regions = {}
    with open(bed_file, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue  # Skip lines starting with '#'

            fields = line.strip().split("\\t")
            if fields[-1].isdigit():
                ref, start, stop, gene_id, gene_name, _ = fields
                regions[(ref, gene_name, gene_id)] = (int(start), int(stop))
            else:
                ref, start, stop, gene_name, gene_id, strand = fields
                regions[(ref, gene_name)] = (int(start), int(stop), strand)
    return regions


def parse_coverage_file(coverage_file, regions):
    """Parse coverage file into a dictionary and return with metrics"""
    coverage_data = defaultdict(dict)
    positions = 0
    total_count = 0
    total_coverage = 0
    region_positions = {p for v in regions.values() for p in range(v[0], v[1] + 1)}
    with open(coverage_file, encoding="utf-8") as f:
        for line in f:
            reference, position, coverage = line.strip().split("\\t")
            positions += 1
            coverage = int(coverage)
            if coverage > 0:
                total_count += 1
                total_coverage += coverage
            if int(position) in region_positions:
                coverage_data[reference][int(position)] = coverage
    return coverage_data, (positions, total_count, total_coverage)


def calculate_coordinates(strand, gene_start, gene_stop, del_start, del_stop):
    """Calculate CDS and amino acid coordinates based on strand orientation"""
    cds_start = 0
    cds_stop = 0
    aa_start = "NA"
    aa_stop = "NA"
    if strand == "forward":
        cds_start = (del_start + 1) - gene_start
        cds_stop = (del_stop + 1) - gene_start
        aa_start = int(cds_start / 3)
        aa_stop = int(cds_stop / 3)
    elif strand == "reverse":
        cds_start = (gene_stop + 1) - del_stop
        cds_stop = (gene_stop + 1) - del_start
        aa_start = int(cds_start / 3)
        aa_stop = int(cds_stop / 3)
    elif strand == "forward_promoter":
        cds_start = del_stop - gene_start
        cds_stop = del_start - gene_start
    elif strand == "reverse_promoter":
        cds_start = gene_stop - del_start
        cds_stop = gene_stop - del_stop

    return cds_start, cds_stop, aa_start, aa_stop


def analyze_deletions(coverage_data, gene_name, coordinates):
    """Identify potential deletions within a given region"""
    result = {}
    gene_start, gene_stop, gene_strand = coordinates
    deletions = [pos for pos in range(gene_start, gene_stop + 1) if pos not in coverage_data or coverage_data[pos] < 1]
    if deletions:
        deletion_start = deletions[0]
        deletion_stop = deletions[-1]
        cds_start, cds_stop, aa_start, aa_stop = calculate_coordinates(
            gene_strand, gene_start, gene_stop, deletion_start, deletion_stop
        )
        result = {
            "Gene": gene_name,
            "SV Length": deletion_stop - deletion_start,
            "Ref Start": deletion_start,
            "Ref Stop": deletion_stop,
            "CDS Start": cds_start,
            "CDS Stop": cds_stop,
            "Amino Acid Start": aa_start,
            "Amino Acid Stop": aa_stop,
        }
    return result


def output_results(sample_id, data):
    """Write results to structural variants txt file"""
    columns = [
        "Sample ID",
        "Gene",
        "SV Length",
        "Ref Start",
        "Ref Stop",
        "CDS Start",
        "CDS Stop",
        "Amino Acid Start",
        "Amino Acid Stop",
    ]
    outfile = f"{sample_id}_structural_variants.txt"
    with open(outfile, "w", encoding="utf-8") as outfile:
        writer = csv.DictWriter(
            outfile,
            fieldnames=columns,
            restval="NA",
            delimiter="\\t",
            lineterminator="\\n",
        )
        writer.writeheader()
        writer.writerows([{**d, "Sample ID": sample_id} for d in data])


def analyze_coverage_data(sample_id, coverage_file, regions_file, genome_stats_file):
    """Analyze coverage data across annotated regions"""
    results = []

    gene_regions = parse_bed_file(regions_file)
    coverage_stats = parse_genome_stats(genome_stats_file)
    coverage_data, _metrics = parse_coverage_file(coverage_file, gene_regions)

    for ref_gene, coords in gene_regions.items():
        if ref_gene not in coverage_stats:
            continue

        ref, gene_name = ref_gene
        avg_depth = float(coverage_stats[ref_gene]["Average Depth"])
        pct_coverage = float(coverage_stats[ref_gene]["Percentage Coverage"])
        if avg_depth > 2 and pct_coverage > 99:
            status = "No large deletion"
            results.append({"Gene": gene_name, "Ref Start": status, "Ref Stop": status})
        elif avg_depth < 2 or pct_coverage < 1:
            status = "Complete deletion"
            results.append({"Gene": gene_name, "Ref Start": status, "Ref Stop": status})
        else:
            deletion_data = analyze_deletions(coverage_data[ref], gene_name, coords)
            results.append(deletion_data)
    output_results(sample_id, results)


if __name__ == "__main__":

    def main():
        prefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else "$meta.id"
        coverage = "${coverage}"
        structures = "BED.txt"
        genome_stats = "${genome_region_coverage}"

        analyze_coverage_data(prefix, coverage, structures, genome_stats)

        versions = {
            "${task.process}": {
                "python": platform.python_version(),
            }
        }

        with open("versions.yml", "w", encoding="utf-8") as f:
            f.write(yaml.dump(versions))

    main()
