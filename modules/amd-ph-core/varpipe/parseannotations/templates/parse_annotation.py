#!/usr/bin/env python3

"""The script creates a final annotation file.

Inputs:
  - snpEff annotated file (tsv derived from vcf)
  - mutation loci CSV file
  - sample_id
  - output filename

"""

import csv
import platform
import re

import yaml


class SnpEffVcfParser:
    """Callable class to parse TB snpEff annotated VCF files"""

    @staticmethod
    def _get_loci(loci_path):
        """Parse the loci file into dictionary objects"""
        with open(loci_path) as handle:
            reader = csv.DictReader(handle, delimiter="\\t")
            loci = {row["H37Rv gene"]: row for row in reader}
            for locus in loci.values():
                locus["start"] = int(locus["start"])
                locus["stop"] = int(locus["stop"])
        return loci

    def __init__(self, loci_path):
        self.hgvs_re = re.compile(r"^p\\.(?P<wt>.*?)(?P<pos>\\d+)(?P<mut>.*?)\$")
        self.indel_re = re.compile(r"^(?P<indel_start>.*?)(?P<indel_pos>\\d+)(?P<indel_stop>.*?)\$")
        self.position_re = re.compile(r"-?\\d+")
        self.gene_map = {
            "erm_37_": "erm(37)",
        }
        self.variant_type_map = {
            "del": "Deletion",
            "ins": "Insertion",
            "dup": "Insertion",
        }
        self.ann_fields = [
            "Allele",
            "Annotation",
            "Annotation_Impact",
            "Gene_Name",
            "Gene_ID",
            "Feature_Type",
            "Feature_ID",
            "Transcript_BioType",
            "Rank",
            "HGVS_c",
            "HGVS_p",
            "cDNA_position",
            "CDS_position",
            "Protein_position",
            "Distance",
            "Errors_Warnings_Info",
        ]
        self.report_fields = [
            "Sample ID",
            "CHROM",
            "POS",
            "REF",
            "ALT",
            "Read Depth",
            "Percent Alt Allele",
            "Annotation",
            "Variant Type",
            "Nucleotide Change",
            "Position within CDS",
            "Amino acid Change",
            "REF Amino acid",
            "ALT Amino acid",
            "Codon Position",
            "Gene Name",
            "Gene ID",
        ]

        # FIXME should this be hard coded??
        self.loci = self._get_loci(loci_path)

    def _filter_vcf(self, sample_id, vcf_path):
        """Parse and filter VCF file"""
        records = []
        with open(vcf_path) as handle:
            reader = csv.DictReader(handle, dialect="excel-tab")
            records = [row for row in reader if row["FILTER"] == "PASS"]

        for record in records:
            record["Sample ID"] = sample_id
            record["Read Depth"] = record["DP"]
            support = float(record["AF"]) * 100
            record["Percent Alt Allele"] = f"{support:.2f}"
            record["support"] = support
            record["annotations"] = [
                dict(zip(self.ann_fields, a.split("|"))) for a in record["ANN"].split("=", 1)[1].split(",")
            ]

        return [record for record in records if record["support"] >= 5.0]

    def _get_variant_type(self, record, annotation):
        variant = "SNP"
        if len(record["REF"]) > 1 and (len(record["REF"]) == len(record["ALT"])):
            variant = "MNP"
        else:
            for short, long in self.variant_type_map.items():
                if short in annotation["HGVS_c"] or short in annotation["HGVS_p"]:
                    variant = long
                    break
        return variant

    def _update_cds_position(self, record, annotation):
        cds_positions = self.position_re.findall(annotation["HGVS_c"])
        if record["Variant Type"] == "MNP":
            record["Position within CDS"] = cds_positions[0]
        else:
            record["Position within CDS"] = (
                cds_positions[0] if "_" not in annotation["HGVS_c"] else f"{cds_positions[0]}-{cds_positions[1]}"
            )

    def _update_coding_fields(self, record, annotation):
        if not annotation["HGVS_p"]:
            return

        codon_positions = self.position_re.findall(annotation["HGVS_p"])
        record["Codon Position"] = (
            codon_positions[0] if "_" not in annotation["HGVS_p"] else f"{codon_positions[0]}-{codon_positions[1]}"
        )

        match = self.hgvs_re.match(annotation["HGVS_p"])
        if match:
            wt, mut = match.group("wt", "mut")
            record["Annotation"] = "Synonymous" if wt == mut else "Non-synonymous"
            record["REF Amino acid"] = wt or codon_positions[0]

            if record["Variant Type"] not in ("Insertion", "Deletion"):
                if any([len(aa) > 3 for aa in (wt, mut)]) and "-" not in record["Codon Position"]:
                    record["Codon Position"] = f"{record["Codon Position"]}-{int(record["Codon Position"])+1}"
                for non_aa in ("del", "*", "?", "fs"):
                    if non_aa in annotation["HGVS_p"]:
                        break
                else:
                    record["ALT Amino acid"] = mut or codon_positions[0]

    def _get_report_annotation(self, record):
        annotation = record["annotations"][0]
        variant_type = self._get_variant_type(record, annotation)
        record.update(
            {
                "Amino acid Change": annotation["HGVS_p"] or "NA",
                "Nucleotide Change": annotation["HGVS_c"],
                "Position within CDS": "NA",
                "Gene Name": self.gene_map.get(annotation["Gene_Name"], annotation["Gene_Name"]),
                "Gene ID": annotation["Gene_ID"],
                "Codon Position": "NA",
                "Annotation": "NA",
                "Variant Type": variant_type,
                "REF Amino acid": "NA",
                "ALT Amino acid": "NA",
            }
        )
        self._update_cds_position(record, annotation)
        self._update_coding_fields(record, annotation)

    @staticmethod
    def _get_modifier_variant_type(record):
        ref, alt = record["REF"], record["ALT"]
        if len(ref) == len(alt):
            variant_type = "MNP" if len(ref) > 1 else "SNP"
        else:
            variant_type = "Deletion" if len(ref) > len(alt) else "Insertion"
        return variant_type

    @staticmethod
    def _filter_modifier_annotations(record):
        filtered_annotations = {
            int(a["Distance"]): a
            for a in record["annotations"]
            if "downstream" not in a["Annotation"] and a["Distance"]
        }

        if filtered_annotations:
            annotation = filtered_annotations[sorted(filtered_annotations.keys())[0]]
            modifier = "upstream"
        else:
            annotation = record["annotations"][0]
            modifier = "downstream"

        record.update(
            {
                "Annotation": "Non-Coding",
                "Gene Name": f"{annotation["Gene_Name"]} {modifier}",
                "Gene ID": f"{annotation["Gene_ID"]} {modifier}",
                "Nucleotide Change": annotation["HGVS_c"],
                "Position within CDS": "NA",
                "Amino acid Change": "NA",
                "REF Amino acid": "NA",
                "ALT Amino acid": "NA",
                "Codon Position": "NA",
            }
        )

    @staticmethod
    def _annotate_known_loci(locus, record):
        mutation_position = int(record["POS"]) - locus["start"] + 1
        cds_position = "NA"
        if "SNP" == record["Variant Type"]:
            hgvs_c = f"c.{mutation_position}{record["REF"]}>{record["ALT"]}"
            cds_position = mutation_position
        elif "MNP" == record["Variant Type"]:
            hgvs_c = f"c.{mutation_position}_{mutation_position+1}del{record["REF"]}ins{record["ALT"]}"
        elif "Insertion" == record["Variant Type"]:
            hgvs_c = f"c.{mutation_position}_{mutation_position+1}ins{record["ALT"][len(record["REF"]):]}"
            cds_position = f"{mutation_position}-{mutation_position+1}"
        elif "Deletion" == record["Variant Type"]:
            if len(record["REF"]) - len(record["ALT"]) == 1:
                hgvs_c = f"c.{mutation_position+len(record["ALT"])}del{record['REF'][len(record["ALT"]):]}"
                cds_position = str(mutation_position + len(record["ALT"]))
            else:
                hgvs_c = f"c.{mutation_position+len(record["ALT"])}_{mutation_position+len(record["REF"])-1}del{record["REF"][len(record["ALT"]):]}"
                cds_position = f"{mutation_position+len(record["ALT"])}-{mutation_position+len(record["REF"])-1}"

        record.update(
            {
                "Annotation": locus["type"],
                "Nucleotide Change": hgvs_c,
                "Position within CDS": cds_position if locus["type"] == "CDS" else "NA",
                "Amino acid Change": "NA",
                "REF Amino acid": "NA",
                "ALT Amino acid": "NA",
                "Codon Position": "NA",
                "Gene Name": locus["H37Rv gene"],
                "Gene ID": locus["H37Rv gene id"],
            }
        )

    def _get_modifier_report_annotation(self, record):
        record["Variant Type"] = self._get_modifier_variant_type(record)
        for locus in self.loci.values():
            if locus["start"] <= int(record["POS"]) <= locus["stop"]:
                self._annotate_known_loci(locus, record)
                break
        else:
            # our location does not match a known loci
            # construct the information from other annotations
            self._filter_modifier_annotations(record)

    def _write_annotation_report(self, output_path, vcf_records):
        with open(output_path, "w") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=self.report_fields,
                extrasaction="ignore",
                delimiter="\\t",
                lineterminator="\\n",
            )
            writer.writeheader()
            for record in vcf_records:
                writer.writerow(record)

    def __call__(self, sample_id, vcf_path, output_path):
        vcf_records = self._filter_vcf(sample_id, vcf_path)

        for record in vcf_records:
            if not record["annotations"][0]["Annotation_Impact"] == "MODIFIER":
                self._get_report_annotation(record)
            else:
                self._get_modifier_report_annotation(record)

        self._write_annotation_report(output_path, vcf_records)


if __name__ == "__main__":

    def main():
        prefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else "$meta.id"
        loci_path = "${mutationloci}"
        vcf_path = "${annotation}"

        input_type = "DR_loci" if "DR_loci" in "${annotation}" else "full"
        output_path = f"{prefix}_{input_type}_Final_annotation.txt"

        task_process = "${task.process}"

        parser = SnpEffVcfParser(loci_path)
        parser(prefix, vcf_path, output_path)

        versions = {
            task_process: {
                "python": platform.python_version(),
            }
        }

        with open("versions.yml", "w", encoding="utf-8") as f:
            f.write(yaml.dump(versions))

    main()
