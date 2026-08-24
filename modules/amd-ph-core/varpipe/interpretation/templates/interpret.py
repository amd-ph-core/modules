#!/usr/bin/env python3

import platform

import yaml


def load_variant_interpretations(input_file):
    # Parallel lists (rather than a dict) because dual-drug genes have two
    # reportable.txt entries for the same variant (one per drug)
    variants = []
    interpretations = []

    with open(input_file) as fh:
        for line in fh:
            lined = line.rstrip("\\r\\n").split("\\t")
            variants.append(lined[0])
            interpretations.append(lined[1])

    return variants, interpretations


def process_variants(input_file, variants, interpretations, primary_drug_dict, secondary_drug_dict):
    array_list = []
    mutations = set()
    position = set()
    precedence = False
    priority_interpretations = {"BDQ-R", "CFZ-R"}

    with open(input_file) as fh:
        for line in fh:
            lined = line.rstrip("\\r\\n").split("\\t")
            if not lined[0].isdigit():
                continue

            variant = f"{lined[1]}_{lined[3] if 'upstream' not in lined[1] else lined[2]}"

            mutations.add(variant)
            position.add(lined[0])

            if variant not in variants:
                continue

            index = variants.index(variant)
            primary_drug = primary_drug_dict[lined[1]]

            if interpretations[index] == f"{primary_drug}-R":
                array_list.append(f"{primary_drug}\\t{variant}\\t{interpretations[index]}")
                if interpretations[index] in priority_interpretations:
                    precedence = True

            # Dual-drug genes: the secondary drug's call immediately follows
            # the primary drug's call in reportable.txt
            if lined[1] in secondary_drug_dict and index + 1 < len(interpretations):
                secondary_drug = secondary_drug_dict[lined[1]]
                if interpretations[index + 1] == f"{secondary_drug}-R":
                    array_list.append(f"{secondary_drug}\\t{variant}\\t{interpretations[index + 1]}")
                    if interpretations[index + 1] in priority_interpretations:
                        precedence = True

            if interpretations[index] == "S":
                array_list.append(f"{primary_drug}\\t{variant}\\t{primary_drug}-S")
                if lined[1] in secondary_drug_dict:
                    secondary_drug = secondary_drug_dict[lined[1]]
                    array_list.append(f"{secondary_drug}\\t{variant}\\t{secondary_drug}-S")

    return array_list, mutations, position, precedence


def append_interpretation_summary(input_file):
    with open(input_file, "a") as fh:
        print("\\n", file=fh)
        print("Interpretations Summary:", file=fh)
        print("Drug\\tVariant\\tInterpretation", file=fh)


def process_deletions(input_file, array_list):
    # Map gene patterns to one or more (drug, call) pairs and index positions.
    # Genes reported against their own CDS coordinates use the Amino Acid
    # Start/Stop columns (7, 8); intragenic/promoter regions use the CDS
    # Start/Stop columns (5, 6) since no amino acid position is available.
    deletion_mapping = {
        "katG": {"calls": [("INH", "INH-R")], "start_idx": 7, "end_idx": 8},
        "Rv1909c": {"calls": [("INH", "INH-R")], "start_idx": 5, "end_idx": 6},
        "furA": {"calls": [("INH", "INH-R")], "start_idx": 7, "end_idx": 8},
        "pncA": {"calls": [("PZA", "PZA-R")], "start_idx": 7, "end_idx": 8},
        "Rv2043c": {"calls": [("PZA", "PZA-R")], "start_idx": 5, "end_idx": 6},
        "fgd1": {"calls": [("DLM", "DLM-R"), ("PMD", "PMD-R")], "start_idx": 7, "end_idx": 8},
        "Rv0406c": {"calls": [("DLM", "DLM-U"), ("PMD", "PMD-U")], "start_idx": 5, "end_idx": 6},
        "mmpR": {"calls": [("BDQ", "BDQ-R"), ("CFZ", "CFZ-R")], "start_idx": 7, "end_idx": 8},
        "mmpL5": {"calls": [("BDQ", "BDQ-S"), ("CFZ", "CFZ-S")], "start_idx": 7, "end_idx": 8},
        "mmpS5": {"calls": [("BDQ", "BDQ-S"), ("CFZ", "CFZ-S")], "start_idx": 7, "end_idx": 8},
        "Rv0677c": {"calls": [("BDQ", "BDQ-U"), ("CFZ", "CFZ-U")], "start_idx": 5, "end_idx": 6},
        "fbiC": {"calls": [("DLM", "DLM-R"), ("PMD", "PMD-R")], "start_idx": 7, "end_idx": 8},
        "Rv1172c": {"calls": [("DLM", "DLM-U"), ("PMD", "PMD-U")], "start_idx": 5, "end_idx": 6},
        "fbiA": {"calls": [("DLM", "DLM-R"), ("PMD", "PMD-R")], "start_idx": 7, "end_idx": 8},
        "fbiB": {"calls": [("DLM", "DLM-R"), ("PMD", "PMD-R")], "start_idx": 7, "end_idx": 8},
        "Rv2982c": {"calls": [("DLM", "DLM-U"), ("PMD", "PMD-U")], "start_idx": 5, "end_idx": 6},
        "Rv2983": {"calls": [("DLM", "DLM-R"), ("PMD", "PMD-R")], "start_idx": 7, "end_idx": 8},
        "Rv3260c": {"calls": [("DLM", "DLM-U"), ("PMD", "PMD-U")], "start_idx": 5, "end_idx": 6},
        "ddn": {"calls": [("DLM", "DLM-R"), ("PMD", "PMD-R")], "start_idx": 7, "end_idx": 8},
        "Rv3546": {"calls": [("DLM", "DLM-U"), ("PMD", "PMD-U")], "start_idx": 5, "end_idx": 6},
        "ethA": {"calls": [("ETO", "ETO-R")], "start_idx": 7, "end_idx": 8},
        "Rv3855": {"calls": [("ETO", "ETO-U")], "start_idx": 5, "end_idx": 6},
    }

    mmp_deletion_seen = False

    with open(input_file) as fh:
        for line in fh:
            lined = line.rstrip("\\r\\n").split("\\t")
            if line.startswith("Sample ID"):
                continue

            if not (lined[2].isdigit() or lined[3] == "Complete deletion"):
                continue

            # Find which gene pattern matches
            for gene, config in deletion_mapping.items():
                if gene in lined[1]:
                    if lined[3] == "Complete deletion":
                        variant = f"{lined[1]}_complete deletion"
                    else:
                        start_idx = config["start_idx"]
                        end_idx = config["end_idx"]
                        variant = f"{lined[1]}_deletion_{lined[start_idx]}_{lined[end_idx]}"

                    for drug, call in config["calls"]:
                        array_list.append(f"{drug}\\t{variant}\\t{call}")

                    if gene in ("mmpL5", "mmpS5"):
                        mmp_deletion_seen = True

    return array_list, mmp_deletion_seen


def process_non_cataloged_variants(
    input_file, position, mutations, variants, array_list, primary_drug_dict, secondary_drug_dict
):
    # Define LOF dictionaries by gene category
    lof_dicts = {
        "group1": {
            "frameshift_variant": "R",
            "frameshift_variant&stop_gained": "R",
            "frameshift_variant&stop_lost&splice_region_variant": "R",
            "stop_gained": "R",
            "start_lost": "R",
            "synonymous_variant": "S",
            "missense_variant": "U",
            "upstream_gene_variant": "U",
            "downstream_gene_variant": "U",
            "disruptive_inframe_insertion": "R",
            "disruptive_inframe_deletion": "R",
            "conservative_inframe_insertion": "R",
            "conservative_inframe_deletion": "R",
            "stop_lost&splice_region_variant": "R",
            "start_lost&conservative_inframe_deletion": "R",
        },
        "group2": {
            "frameshift_variant": "R",
            "frameshift_variant&stop_gained": "R",
            "frameshift_variant&stop_lost&splice_region_variant": "R",
            "stop_gained": "R",
            "start_lost": "R",
            "synonymous_variant": "S",
            "missense_variant": "U",
            "upstream_gene_variant": "U",
            "downstream_gene_variant": "U",
            "disruptive_inframe_insertion": "U",
            "disruptive_inframe_deletion": "U",
            "conservative_inframe_insertion": "U",
            "conservative_inframe_deletion": "U",
            "stop_lost&splice_region_variant": "R",
            "start_lost&conservative_inframe_deletion": "R",
        },
        "others": {
            "frameshift_variant": "U",
            "frameshift_variant&stop_gained": "U",
            "frameshift_variant&stop_lost&splice_region_variant": "U",
            "stop_gained": "U",
            "start_lost": "U",
            "synonymous_variant": "U",
            "missense_variant": "U",
            "upstream_gene_variant": "U",
            "downstream_gene_variant": "U",
            "disruptive_inframe_insertion": "U",
            "disruptive_inframe_deletion": "U",
            "conservative_inframe_insertion": "U",
            "conservative_inframe_deletion": "U",
            "stop_lost&splice_region_variant": "U",
            "start_lost&conservative_inframe_deletion": "U",
        },
    }

    # Genes with dedicated frameshift/inframe-indel interpretation rules
    group1_genes = ("rpoB", "atpE")
    group2_genes = (
        "katG",
        "pncA",
        "fbiA",
        "fbiB",
        "fbiC",
        "ddn",
        "fgd1",
        "Rv2983",
        "pepQ",
        "ethA",
        "mmpR",
        "mmpL5",
        "mmpS5",
    )

    mmp_lof_seen = False

    with open(input_file) as fh:
        for line in fh:
            lined = line.rstrip("\\r\\n").split("\\t")
            if line.startswith("Sample ID") or lined[2] not in position:
                continue

            annot = lined[29].split(",")
            for x in annot:
                subannot = x.split("|")
                gene = subannot[3]
                variant_type = subannot[1]

                # Skip if gene not in drug dict or variant type not recognized
                if gene not in primary_drug_dict:
                    continue

                # Determine which LOF dictionary to use
                if gene in group1_genes:
                    lof_dict = lof_dicts["group1"]
                elif gene in group2_genes:
                    lof_dict = lof_dicts["group2"]
                else:
                    lof_dict = lof_dicts["others"]

                # Skip if variant type not in LOF dict
                if variant_type not in lof_dict:
                    continue

                interpretation = lof_dict[variant_type]

                # Check both regular and upstream variants
                for interpret_suffix, pos_field in [("", "10"), ("upstream", "9")]:
                    interpret_key = (
                        f"{gene}{' ' if interpret_suffix else ''}{interpret_suffix}_{subannot[int(pos_field)]}"
                    )

                    if interpret_key not in mutations or interpret_key in variants:
                        continue

                    primary_drug = primary_drug_dict[gene]
                    array_list.append(f"{primary_drug}\\t{interpret_key}\\t{primary_drug}-{interpretation}")

                    # Genes with a secondary drug get calls recorded for both
                    if gene in group2_genes and gene in secondary_drug_dict:
                        secondary_drug = secondary_drug_dict[gene]
                        array_list.append(
                            f"{secondary_drug}\\t{interpret_key}\\t{secondary_drug}-{interpretation}"
                        )

                    if gene in ("mmpS5", "mmpL5") and interpret_suffix == "" and interpretation == "R":
                        mmp_lof_seen = True

    return array_list, mmp_lof_seen


def process_review_coverage(input_file, primary_drug_dict, secondary_drug_dict):
    primary_targets = {}
    secondary_targets = {}

    # Amplicon coverage rows for paired-gene targets use a combined gene label
    special_secondary_targets = {
        "mmpL5_S5": "CFZ",
        "fbiA_B": "PMD",
    }

    with open(input_file) as fh:
        for line in fh:
            lined = line.rstrip("\\r\\n").split("\\t")
            if line.startswith("Sample ID"):
                continue

            if "Review" not in lined:
                continue

            gene = lined[4]
            if gene in primary_drug_dict:
                primary_targets[gene] = primary_drug_dict[gene]
            if gene in secondary_drug_dict:
                secondary_targets[gene] = secondary_drug_dict[gene]
            if gene in special_secondary_targets:
                secondary_targets[gene] = special_secondary_targets[gene]

    return primary_targets, secondary_targets


def build_final_output(
    array_list,
    primary_targets,
    secondary_targets,
    sample_id,
    input_summary_file,
    mmp_lof_seen,
    mmp_deletion_seen,
    precedence,
):
    drug_order = ["INH", "RIF", "PZA", "FQ", "EMB", "CFZ", "DLM", "LZD", "BDQ", "PMD", "ETO"]
    second_line_drugs = {"CFZ", "DLM", "LZD", "BDQ", "PMD", "ETO"}

    # Dictionary to store drug information
    drug_data = {drug: {"variants": "", "interpretation": ""} for drug in drug_order}

    # Build interpretation priority map
    interp_priority = {"R": 3, "U": 2, "S": 1}

    # Process array list
    for entry in array_list:
        drug, variant, interpretation = entry.split("\\t")

        # Add variant to drug's variant list
        if not drug_data[drug]["variants"]:
            drug_data[drug]["variants"] = variant
        else:
            drug_data[drug]["variants"] += f",{variant}"

        # Update interpretation based on priority
        if not drug_data[drug]["interpretation"]:
            drug_data[drug]["interpretation"] = interpretation
        else:
            # Extract current and new status (R, U, S)
            current_status = drug_data[drug]["interpretation"].split("-")[-1]
            new_status = interpretation.split("-")[-1]

            # Keep the interpretation with higher priority
            if interp_priority.get(new_status, 0) > interp_priority.get(current_status, 0):
                drug_data[drug]["interpretation"] = interpretation

    # An mmpL5/mmpS5 loss-of-function call or deletion downgrades BDQ/CFZ to
    # uncertain unless a higher-precedence BDQ-R/CFZ-R call was also recorded
    if (mmp_lof_seen or mmp_deletion_seen) and not precedence:
        drug_data["BDQ"]["interpretation"] = "BDQ-U"
        drug_data["CFZ"]["interpretation"] = "CFZ-U"

    # Create interpretation output file
    interpretation_filename = f"{sample_id}_interpretation.txt"

    # Create summary output file
    summary_filename = f"{sample_id}_summary.txt"

    review_targets = set(primary_targets.values()) | set(secondary_targets.values())

    # Write CLI output to interpretation file
    with open(interpretation_filename, "w") as fh:
        # Write header
        fh.write("Sample ID\\tDrug\\tVariant\\tInterpretation\\n")

        # Write drug information
        for drug in drug_order:
            no_variant_msg = "No reportable variant detected"
            default_interpretation = f"{drug}-U" if drug in second_line_drugs else f"{drug}-S"

            if drug_data[drug]["variants"]:
                # Drug has variants
                line = f"{sample_id}\\t{drug}\\t{drug_data[drug]['variants']}\\t{drug_data[drug]['interpretation']}\\n"
                fh.write(line)
                print(line.rstrip())
            elif drug in review_targets:
                # Drug needs review
                line = f"{sample_id}\\t{drug}\\t{no_variant_msg}\\tReview coverage\\n"
                fh.write(line)
                print(line.rstrip())
            else:
                # Drug defaults to susceptible (first-line) or uncertain (second-line)
                line = f"{sample_id}\\t{drug}\\t{no_variant_msg}\\t{default_interpretation}\\n"
                fh.write(line)
                print(line.rstrip())

    # Append interpretations directly to the original summary file
    with open(input_summary_file, "a") as fh:
        for drug in drug_order:
            no_variant_msg = "No reportable variant detected"
            default_interpretation = f"{drug}-U" if drug in second_line_drugs else f"{drug}-S"

            if drug_data[drug]["variants"]:
                print(f"{drug}\\t{drug_data[drug]['variants']}\\t{drug_data[drug]['interpretation']}", file=fh)
            elif drug in review_targets:
                print(f"{drug}\\t{no_variant_msg}\\tReview coverage", file=fh)
            else:
                print(f"{drug}\\t{no_variant_msg}\\t{default_interpretation}", file=fh)

    # Copy the original summary file to the sample_id_summary.txt file
    with open(input_summary_file) as src, open(summary_filename, "w") as dst:
        dst.write(src.read())


def main():
    prefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else "$meta.id"
    reported = "${reported}"
    summary_in = "${summary}"
    structural_variants = "${structuralVariants}"
    dr_loci_annotation = "${drLociAnnotation}"
    target_region_coverage = "${targetRegionCoverage}"
    task_process = "${task.process}"

    # Primary drug assignment per gene (first-choice interpretation)
    primary_drug_dict = {
        "katG": "INH",
        "fabG1": "ETO",
        "fabG1 upstream": "ETO",
        "rpoB": "RIF",
        "pncA": "PZA",
        "pncA upstream": "PZA",
        "gyrA": "FQ",
        "gyrB": "FQ",
        "embB": "EMB",
        "fbiA": "DLM",
        "fbiA upstream": "DLM",
        "fbiB": "DLM",
        "fbiC": "DLM",
        "fbiC upstream": "DLM",
        "ddn": "DLM",
        "ddn upstream": "DLM",
        "fgd1": "DLM",
        "fgd1 upstream": "DLM",
        "Rv2983": "DLM",
        "Rv2983 upstream": "DLM",
        "rplC": "LZD",
        "rplC upstream": "LZD",
        "rrl": "LZD",
        "mmpR": "BDQ",
        "mmpR upstream": "BDQ",
        "mmpL5": "BDQ",
        "mmpS5": "BDQ",
        "atpE": "BDQ",
        "atpE upstream": "BDQ",
        "pepQ": "BDQ",
        "pepQ upstream": "BDQ",
        "ethA": "ETO",
        "ethA upstream": "ETO",
        "inhA": "ETO",
    }

    # Secondary drug assignment per gene, where applicable
    secondary_drug_dict = {
        "fabG1": "INH",
        "fabG1 upstream": "INH",
        "fbiA": "PMD",
        "fbiA upstream": "PMD",
        "fbiB": "PMD",
        "fbiC": "PMD",
        "fbiC upstream": "PMD",
        "ddn": "PMD",
        "ddn upstream": "PMD",
        "fgd1": "PMD",
        "fgd1 upstream": "PMD",
        "Rv2983": "PMD",
        "Rv2983 upstream": "PMD",
        "mmpR": "CFZ",
        "mmpR upstream": "CFZ",
        "mmpL5": "CFZ",
        "mmpS5": "CFZ",
        "pepQ": "CFZ",
        "pepQ upstream": "CFZ",
        "inhA": "INH",
    }

    # Load variants and interpretations
    variants, interpretations = load_variant_interpretations(reported)

    # Process variant file
    array_list, mutations, position, precedence = process_variants(
        summary_in, variants, interpretations, primary_drug_dict, secondary_drug_dict
    )

    # Append interpretation summary header to file
    append_interpretation_summary(summary_in)

    # Process deletion variants
    array_list, mmp_deletion_seen = process_deletions(structural_variants, array_list)

    # Process non-cataloged variants
    array_list, mmp_lof_seen = process_non_cataloged_variants(
        dr_loci_annotation, position, mutations, variants, array_list, primary_drug_dict, secondary_drug_dict
    )

    # Process review coverage targets
    primary_targets, secondary_targets = process_review_coverage(
        target_region_coverage, primary_drug_dict, secondary_drug_dict
    )

    # Build and write final output
    build_final_output(
        array_list,
        primary_targets,
        secondary_targets,
        prefix,
        summary_in,
        mmp_lof_seen,
        mmp_deletion_seen,
        precedence,
    )

    versions = {
        task_process: {
            "python": platform.python_version(),
        }
    }

    with open("versions.yml", "w", encoding="utf-8") as f:
        f.write(yaml.dump(versions))


if __name__ == "__main__":
    main()
