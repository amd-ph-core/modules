#!/usr/bin/env python3

import platform

import yaml


def load_variant_interpretations(input_file):
    variant_interp_map = {}

    with open(input_file) as fh:
        for line in fh:
            lined = line.rstrip("\\r\\n").split("\\t")
            variant_interp_map[lined[0]] = lined[1]

    return variant_interp_map


def process_variants(input_file, variant_interp_map, drug_dict):
    array_list = []
    mutations = set()
    position = set()

    with open(input_file) as fh:
        for line in fh:
            lined = line.rstrip("\\r\\n").split("\\t")
            if not lined[0].isdigit():
                continue

            variant = f"{lined[1]}_{lined[3] if 'upstream' not in lined[1] else lined[2]}"

            mutations.add(variant)
            position.add(lined[0])

            if variant in variant_interp_map:
                interpretation = variant_interp_map[variant]
                drug = drug_dict[lined[1]]

                if interpretation == "S":
                    array_list.append(f"{drug}\\t{variant}\\t{drug}-{interpretation}")
                else:
                    array_list.append(f"{drug}\\t{variant}\\t{interpretation}")

    return array_list, mutations, position


def append_interpretation_summary(input_file):
    with open(input_file, "a") as fh:
        print("\\n", file=fh)
        print("Interpretations Summary:", file=fh)
        print("Drug\\tVariant\\tInterpretation", file=fh)


def process_deletions(input_file, array_list):
    # Map gene patterns to drug and index positions
    deletion_mapping = {
        "katG": {"drug": "INH", "complete_idx": -1, "start_idx": 7, "end_idx": 8},
        "Rv1909c": {"drug": "INH", "complete_idx": -1, "start_idx": 5, "end_idx": 6},
        "furA": {"drug": "INH", "complete_idx": -1, "start_idx": 7, "end_idx": 8},
        "pncA": {"drug": "PZA", "complete_idx": -1, "start_idx": 7, "end_idx": 8},
        "Rv2043c": {"drug": "PZA", "complete_idx": -1, "start_idx": 5, "end_idx": 6},
    }

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
                    drug = config["drug"]
                    if lined[3] == "Complete deletion":
                        array_list.append(f"{drug}\\t{lined[1]}_complete deletion\\t{drug}-R")
                    elif lined[2].isdigit():
                        start_idx = config["start_idx"]
                        end_idx = config["end_idx"]
                        array_list.append(
                            f"{drug}\\t{lined[1]}_deletion_{lined[start_idx]}_{lined[end_idx]}\\t{drug}-R"
                        )

    return array_list


def process_non_cataloged_variants(input_file, position, mutations, variant_interp_map, array_list, drug_dict):
    # Define LOF dictionaries by gene category
    lof_dicts = {
        "rpoB": {
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
        "katG_pncA": {
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
                if gene not in drug_dict:
                    continue

                # Determine which LOF dictionary to use
                if gene == "rpoB":
                    lof_dict = lof_dicts["rpoB"]
                elif gene in ["katG", "pncA"]:
                    lof_dict = lof_dicts["katG_pncA"]
                else:
                    lof_dict = lof_dicts["others"]

                # Skip if variant type not in LOF dict
                if variant_type not in lof_dict:
                    continue

                # Check both regular and upstream variants
                for interpret_suffix, pos_field in [("", "10"), ("upstream", "9")]:
                    interpret_key = (
                        f"{gene}{' ' if interpret_suffix else ''}{interpret_suffix}_{subannot[int(pos_field)]}"
                    )

                    if interpret_key in mutations and interpret_key not in variant_interp_map:
                        drug = drug_dict[gene]
                        interpretation = lof_dict[variant_type]
                        array_list.append(f"{drug}\\t{interpret_key}\\t{drug}-{interpretation}")

    return array_list


def process_review_coverage(input_file, drug_dict):
    target_drugs = {}

    with open(input_file) as fh:
        for line in fh:
            lined = line.rstrip("\\r\\n").split("\\t")
            if line.startswith("Sample ID"):
                continue

            if "Review" in lined and lined[4] in drug_dict:
                target_drugs[lined[4]] = drug_dict[lined[4]]

    return target_drugs


def build_final_output(array_list, target_drugs, sample_id, input_summary_file):
    # Dictionary to store drug information
    drug_data = {drug: {"variants": "", "interpretation": ""} for drug in ["INH", "RIF", "PZA", "FQ", "EMB"]}

    # Build interpretation priority map
    interp_priority = {"R": 3, "U": 2, "S": 1}

    # Process array list
    for entry in array_list:
        drug, variant, interpretation = entry.split("\t")

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

    # Create interpretation output file
    interpretation_filename = f"{sample_id}_interpretation.txt"

    # Create summary output file
    summary_filename = f"{sample_id}_summary.txt"

    # Write CLI output to interpretation file
    with open(interpretation_filename, "w") as fh:
        # Write header
        fh.write("Sample ID\\tDrug\\tVariant\\tInterpretation\\n")

        # Write drug information
        for drug in drug_data:
            no_variant_msg = "No reportable variant detected"

            if drug_data[drug]["variants"]:
                # Drug has variants
                line = f"{sample_id}\\t{drug}\\t{drug_data[drug]['variants']}\\t{drug_data[drug]['interpretation']}\\n"
                fh.write(line)
                print(line.rstrip())
            elif drug in target_drugs.values():
                # Drug needs review
                line = f"{sample_id}\\t{drug}\\t{no_variant_msg}\\tReview coverage\\n"
                fh.write(line)
                print(line.rstrip())
            else:
                # Drug is susceptible
                line = f"{sample_id}\\t{drug}\\t{no_variant_msg}\\t{drug}-S\\n"
                fh.write(line)
                print(line.rstrip())

    # Append interpretations directly to the original summary file
    with open(input_summary_file, "a") as fh:
        for drug in drug_data:
            no_variant_msg = "No reported variant detected"

            if drug_data[drug]["variants"]:
                print(f"{drug}\\t{drug_data[drug]['variants']}\\t{drug_data[drug]['interpretation']}", file=fh)
            elif drug in target_drugs.values():
                print(f"{drug}\\t{no_variant_msg}\\tReview coverage", file=fh)
            else:
                print(f"{drug}\\t{no_variant_msg}\\t{drug}-S", file=fh)

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

    # Define drug dictionary
    drug_dict = {
        "katG": "INH",
        "fabG1": "INH",
        "fabG1 upstream": "INH",
        "rpoB": "RIF",
        "pncA": "PZA",
        "pncA upstream": "PZA",
        "gyrA": "FQ",
        "gyrB": "FQ",
        "embB": "EMB",
        "inhA": "INH",
    }

    # Load variants and interpretations
    variant_interp_map = load_variant_interpretations(reported)

    # Process variant file
    array_list, mutations, position = process_variants(summary_in, variant_interp_map, drug_dict)

    # Append interpretation summary header to file
    append_interpretation_summary(summary_in)

    # Process deletion variants
    array_list = process_deletions(structural_variants, array_list)

    # Process non-cataloged variants
    array_list = process_non_cataloged_variants(
        dr_loci_annotation, position, mutations, variant_interp_map, array_list, drug_dict
    )

    # Process review coverage targets
    target_drugs = process_review_coverage(target_region_coverage, drug_dict)

    # Build and write final output
    build_final_output(array_list, target_drugs, prefix, summary_in)

    versions = {
        task_process: {
            "python": platform.python_version(),
        }
    }

    with open("versions.yml", "w", encoding="utf-8") as f:
        f.write(yaml.dump(versions))


if __name__ == "__main__":
    main()
