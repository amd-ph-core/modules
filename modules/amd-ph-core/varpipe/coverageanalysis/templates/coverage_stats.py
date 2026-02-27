#!/usr/bin/env python3

"""
Reads in coverage depth output from samtools depth along with mapped/unmapped
read counts, and cross-references against gene and target feature region
annotations to produce coverage statistics and QC logs for failed samples
"""

import csv
import datetime
import platform
from collections import defaultdict

import yaml


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


def check_coverage(coverage_data, reference, start, stop, threshold):
    """Check for low coverage within a given region"""
    for position in range(start, stop + 1):
        if position not in coverage_data[reference] or coverage_data[reference][position] < threshold:
            return True  # Low coverage detected
    return False


def get_region_coverage(coverage_data, reference, start, stop):
    """Obtain coverage data for a given region"""
    total_coverage = 0
    total_count = 0
    total_coverage = 0
    for position in range(start, stop + 1):
        if position in coverage_data[reference]:
            if coverage_data[reference][position] > 0:
                total_count += 1
                total_coverage += coverage_data[reference][position]
    return ((stop - start + 1), total_count, total_coverage)


def calculate_stats(metrics):
    """Calculate average depth and percentage coverage"""
    total_positions, total_count, total_coverage = metrics
    average_coverage = 0.0
    percentage_coverage = 0.0
    if total_count > 0:
        average_coverage = round(total_coverage / total_count, 2)
        percentage_coverage = round((total_count / total_positions) * 100, 2)
    return average_coverage, percentage_coverage


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


def check_qc_failed(stats_dict):
    """Determine if sample has failed to meet QC thresholds"""
    return (
        stats_dict["Percent Reads Mapped"] < 90
        or stats_dict["Average Genome Coverage Depth"] < 30
        or stats_dict["Percent Reference Genome Covered"] < 90
    )


def analyze_coverage_stats(sample_id, genes, features, coverage, unmapped, mapped):
    """Main script logic for coverage data analysis"""
    # Ensure single file paths
    if isinstance(coverage, list):
        coverage = coverage[0]
    if isinstance(unmapped, list):
        unmapped = unmapped[0]
    if isinstance(mapped, list):
        mapped = mapped[0]

    gene_regions = parse_bed_file(genes)

    feature_regions = {}
    for feature_file in features:
        feature_regions.update(parse_bed_file(feature_file))

    coverage_data, genome_metrics = parse_coverage_file(coverage, gene_regions | feature_regions)

    gene_flags = []
    for ref_gene, coords in gene_regions.items():
        ref, gene_name, gene_id = ref_gene
        start, stop = coords
        low_coverage_status = check_coverage(coverage_data, ref, start, stop, 10)
        flag = "Review" if low_coverage_status else "No deletion"
        gene_flags.append((sample_id, ref, start, stop, gene_name, gene_id, flag))

    feature_stats = []
    for ref_gene, coords in feature_regions.items():
        ref, gene_name, gene_id = ref_gene
        start, stop = coords
        avg, pct = calculate_stats(get_region_coverage(coverage_data, ref, start, stop))
        feature_stats.append((sample_id, ref, start, stop, gene_name, gene_id, avg, pct))

    # Calculate mapped/unmapped percentage
    with open(unmapped, encoding="utf-8") as f:
        unmapped_reads = int(f.readline().strip())
    with open(mapped, encoding="utf-8") as f:
        mapped_reads = int(f.readline().strip())
    pct_reads_mapped = round((mapped_reads / unmapped_reads) * 100, 2)

    # Calculate genome average coverage and percentage coverage
    genome_avg_depth, genome_percentage = calculate_stats(genome_metrics)

    # Calculate number of flagged genes
    total_flagged = len([flag for *_, flag in gene_flags if flag == "Review"])

    # Output coverage statistics
    coverage_stats_dict = {
        "Sample ID": sample_id,
        "Sample Name": sample_id,
        "Percent Reads Mapped": pct_reads_mapped,
        "Average Genome Coverage Depth": genome_avg_depth,
        "Percent Reference Genome Covered": genome_percentage,
        "Coverage Drop": total_flagged,
        "Pipeline Version": pipeline_version,
        "Date": datetime.datetime.now().astimezone().strftime("%Y/%m/%d %H:%M:%S %Z"),
    }
    outfile = f"{sample_id}_stats.txt"
    with open(outfile, "w", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\\t", lineterminator="\\n")
        writer.writerow(list(coverage_stats_dict.keys()))
        writer.writerow(list(coverage_stats_dict.values()))

    # Write to log file if the sample failed QC
    if check_qc_failed(coverage_stats_dict):
        outfile = f"{sample_id}_qc_log.txt"
        with open(outfile, "w", encoding="utf-8") as f:
            print(f"Sample {sample_id} failed QC checks", file=f)
    else:
        with open(".qc_passed.txt", "w", encoding="utf-8") as f:
            pass  # Additional log used for filtering within Nextflow

        # Output annotated coverage files if sample passed QC
        header = ["Sample ID", "Reference", "Start", "Stop", "Gene Name", "Gene ID"]

        # Output gene regions with flagged status
        outfile = f"{sample_id}_target_region_coverage.txt"
        with open(outfile, "w", encoding="utf-8") as f:
            writer = csv.writer(f, delimiter="\\t", lineterminator="\\n")
            writer.writerow(header + ["Flag"])
            # Sort target region output by start position
            gene_flags.sort(key=lambda x: x[2])
            writer.writerows(list(gene_flags))

        # Output feature regions with average depth and percentage covered
        outfile = f"{sample_id}_genome_region_coverage.txt"
        with open(outfile, "w", encoding="utf-8") as f:
            writer = csv.writer(f, delimiter="\\t", lineterminator="\\n")
            writer.writerow(header + ["Average Depth", "Percentage Coverage"])
            # Sort genome region output by start position
            feature_stats.sort(key=lambda x: x[2])
            writer.writerows(list(feature_stats))


if __name__ == "__main__":
    prefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else "$meta.id"
    genes = "${amp}"
    features = ["bed_1.txt", "bed_2.txt"]

    # Find files
    coverage = "${coverage}"
    unmapped = "${unmapped}"
    mapped = "${mapped}"
    pipeline_version = "Varpipeline: Varpipe_wgs_${pipeline_version}"

    analyze_coverage_stats(prefix, genes, features, coverage, unmapped, mapped)

    versions = {
        "${task.process}": {
            "python": platform.python_version(),
        }
    }

    with open("versions.yml", "w", encoding="utf-8") as f:
        f.write(yaml.dump(versions))
