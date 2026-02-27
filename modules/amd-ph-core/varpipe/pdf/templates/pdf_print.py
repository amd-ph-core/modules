#!/usr/bin/env python3

import platform

import fpdf
import yaml


class SampleSummaryPdf(fpdf.FPDF):
    def __init__(self, summary_data, *args, **kwargs):
        self.summary_data = summary_data
        super().__init__(*args, **kwargs)

    def header(self):
        self.set_font("Arial", "", 9)
        self.ln(0.5)
        for k, v in self.summary_data.items():
            self.cell(30, 5, k)
            self.cell(30, 5, v)
            self.ln(5)
        self.ln(10)

    def footer(self):
        # Position at 1.5 cm from bottom
        self.set_y(-15)
        self.set_font("Arial", "I", 8)
        self.cell(0, 10, "Page " + str(self.page_no()) + " of {nb}", 0, 0, "C")


def parse_sample_summary(input_file):
    """Parse a TB sample summary text file"""

    def _parse_sample_summary():
        """Sample Summary section"""
        summary_data = {}
        in_sample_summary = False

        for line in infile:
            line = line.strip()
            if line.startswith("Sample Summary:"):
                in_sample_summary = True
                continue
            elif in_sample_summary and not line:
                break
            elif in_sample_summary:
                if ":" in line:
                    key, value = line.split(":", 1)
                    summary_data[key.strip()] = value.strip()
                else:
                    print(f"Warning: Line without colon in Sample Summary section: {line}")

        return summary_data

    def _parse_section(section_name):
        """named tsv sections"""
        rows = []
        in_section = False

        for line in infile:
            line = line.strip()
            if line and section_name in line:
                # locate section by title
                in_section = True
                continue
            elif in_section and not line:
                break
            elif in_section and line:
                rows.append(line.split("\\t"))

        return rows

    section_titles = (
        "Target Coverage Summary:",
        "Variant Summary:",
        "Interpretations Summary:",
    )
    with open(input_file) as infile:
        sample_summary = _parse_sample_summary()
        sections = {title: _parse_section(title) for title in section_titles}

    return sample_summary, sections


def create_pdf(input_file, output_file):
    """Format pdf with common header and one table per page"""

    def _init_pdf():
        pdf = SampleSummaryPdf(sample_summary, "L", "mm", "A4")
        pdf.alias_nb_pages()
        pdf.add_page()
        return pdf

    def _format_section_header():
        pdf.set_font("Arial", "B", 9)
        pdf.cell(0, 5, header, 0, 1)
        pdf.set_font("Arial", "", 9)
        pdf.ln(0.5)

    def _format_coverage():
        pdf.set_font("Arial", "", 9)
        for row in rows:
            for value in row:
                pdf.cell(40, 5, value, border=0)
            pdf.ln(5)
        pdf.ln(160)

    def _format_variants():
        for row in rows:
            for value in row:
                pdf.set_font("Arial", "", 8 if len(value) > 18 else 9)
                pdf.cell(
                    38,
                    5,
                    f"{value[:21]}*" if len(value) > 20 and "c." in value else value,
                    border=0,
                )
            pdf.ln(5)
        pdf.ln(160)

    def _format_interpretations():
        pdf.set_font("Arial", "", 9)
        for row in rows:
            for value in row:
                pdf.cell(120, 5, f"{value[:69]}*" if len(value) > 70 else value, border=0)
            pdf.ln(5)

    section_formatters = {
        "Target Coverage Summary:": _format_coverage,
        "Variant Summary:": _format_variants,
        "Interpretations Summary:": _format_interpretations,
    }
    sample_summary, sections = parse_sample_summary(input_file)

    pdf = _init_pdf()
    for header, rows in sections.items():
        _format_section_header()
        section_formatters[header]()
    pdf.output(output_file)


if __name__ == "__main__":

    def main():
        prefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else "$meta.id"
        input_file = "${summary}"
        output_file = f"{prefix}_report.pdf"

        create_pdf(input_file, output_file)

        versions = {
            "${task.process}": {
                "python": platform.python_version(),
            }
        }

        with open("versions.yml", "w", encoding="utf-8") as f:
            f.write(yaml.dump(versions))

    main()
