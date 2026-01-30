#!/usr/bin/env python3

import csv
import platform

import yaml

tribe_map = {
    str(k): v
    for k, v in enumerate(
        [
            "Indo-Oceanic",
            "East-Asian",
            "East-African-Indian",
            "Euro-American",
            "West-Africa 1",
            "West-Africa 2",
            "Ethiopian",
        ],
        start=1,
    )
}

animal_map = {
    "BCG": ("Bovis-BCG", "M. bovis-BCG"),
    "BOV": ("Bovis", "M. bovis"),
    "ORYX": ("Oryx", "M. orygis"),
    "CAP": ("Caprae", "M. caprae"),
}
snp_tags = {
    "931123": "linfour",
    "1759252": "hrv37",
    "2831482": "CAP",
    "2289073": "BOV",
    "8624": "BCG",
    "2726378": "ORYX",
}
tag_snps = {v: k for k, v in snp_tags.items()}


def parse_lineage(lineage_marker_file, annotation_file, sample_id):
    def _get_lineage_detail():
        sublin, prevlin, prevsub = "", [], []
        for snp in snp_data:
            if snp["POS"] in position_map:
                lineage = position_map[snp["POS"]]["##Lineage"]  # FIXME the data should change, no ## needed
                alt = position_map[snp["POS"]]["Alt"]
                if alt == snp["ALT"]:
                    if len(lineage) > 1:
                        prevsub.append(lineage)
                        sublin = max(prevsub, key=len)
                        text_report.append(f"SNP {snp['POS']} suggests sub-lineage: {lineage}")
                    else:
                        prevlin.append(lineage)
                        text_report.append(f"SNP {snp['POS']} suggests lineage: {lineage}")
        return sublin, prevlin, prevsub

    def _check_discordance():
        discordance = False
        split_first = sublin.split(".")
        for sub in prevsub:
            split_lin = sub.split(".")
            if split_lin[0] != split_first[0] or split_lin[1] != split_first[1]:
                discordance = True
                break
        return discordance, split_first

    def _check_sublin_discordance():
        discordance = False
        for sub in prevsub:
            split_lin = sub.split(".")
            if split_lin[0] != prevlin[0] or split_lin[0] != split_first[0]:
                discordance = True
                break
        return discordance

    def _check_multiple_discordance():
        discordance = False
        if any(prevlin[0] != lin for lin in prevlin):
            discordance = True
        return discordance

    def _report_mixed_lineage():
        text_report.append("no precise lineage inferred")
        tsv_report.append([sample_id, "mixed lineage(s)", "mixed lineage(s)"])

    def _report_tribe(tribe):
        if discordance or animals:
            _report_mixed_lineage()
        else:
            text_report.append(f"Lineage: {tribe} {tribe_map[tribe]}")
            tsv_report.append([sample_id, split_first[0], tribe_map[tribe]])

    def _report_linfour():
        if animals:
            _report_mixed_lineage()
        else:
            text_report.append("Absence of SNP 931123 suggests lineage 4")
            if "hrv37" not in tags:
                text_report.append("Absence of SNP 1759252 suggests sublineage 4.9")
            text_report.append("Lineage: 4 Euro-American")
            tsv_report.append([sample_id, "4", "Euro-American"])

    def _report_animals():
        for tag, (short_species, species) in animal_map.items():
            if tag in tags:
                text_report.append(f"Lineage: {short_species}")
                tsv_report.append([sample_id, short_species, species])
                break
        else:
            text_report.append("No Informative SNPs detected")
            tsv_report.append(
                [
                    sample_id,
                    "No Informative SNPs detected",
                    "No Informative SNPs detected",
                ]
            )

    def _report_prevlin():
        text_report.append(f"Lineage: {prevlin[0]} {tribe_map[prevlin[0]]}")
        tsv_report.append([sample_id, prevlin[0], tribe_map[prevlin[0]]])

    def _report_multiple_prevlin():
        if discordance:
            text_report.append("no concordance between predicted lineage and sublineage(s)")
            tsv_report.append([sample_id, "mixed lineage(s)", "mixed lineage(s)"])
        else:
            _report_prevlin()

    with open(lineage_marker_file) as handle:
        marker_parser = csv.DictReader(handle, dialect="excel-tab")
        position_map = {marker["Position"]: marker for marker in marker_parser}

    with open(annotation_file) as handle:
        snp_data = list(csv.DictReader(handle, dialect="excel-tab"))

    snps = {snp["POS"] for snp in snp_data}
    tags = {snp_tags[pos] for pos in snps if pos in snp_tags}
    animals = any(tag in tags for tag in animal_map)

    tsv_report = [["Sample ID", "Lineage", "Lineage Name"]]
    text_report = []

    for tag, (short_species, species) in animal_map.items():
        if tag in tags:
            text_report.append(f"SNP {tag_snps[tag]} suggests {species}")

    sublin, prevlin, prevsub = _get_lineage_detail()
    split_first = sublin.split(".") if prevsub and sublin else ["NA"]
    if not prevlin:
        discordance, split_first = _check_discordance()
        if len(split_first) > 1:
            _report_tribe(split_first[0])
        elif "linfour" not in tags:
            _report_linfour()
        else:
            _report_animals()
    elif len(prevlin) == 1:
        if sublin:
            discordance = _check_sublin_discordance()
            _report_tribe(prevlin[0])
        else:
            _report_prevlin()
    else:
        discordance = _check_multiple_discordance()
        _report_multiple_prevlin()

    # FIXME if output filenames can be derived, we should not include lineage_report_file as a script param
    with open(f"{sample_id}.lineage_report.txt", "w") as handle:
        writer = csv.writer(handle, dialect="excel-tab", lineterminator="\\n")
        writer.writerows(tsv_report)

    with open(f"{sample_id}_Lineage.txt", "w") as handle:
        print("\\n".join(text_report), file=handle)


if __name__ == "__main__":

    def main():
        prefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else "$meta.id"
        lineage_marker_file = "${lineages}"
        annotation_file = "${fullFinalAnnotation}"

        parse_lineage(lineage_marker_file, annotation_file, prefix)

        versions = {
            "${task.process}": {
                "python": platform.python_version(),
            }
        }

        with open("versions.yml", "w", encoding="utf-8") as f:
            f.write(yaml.dump(versions))

    main()
