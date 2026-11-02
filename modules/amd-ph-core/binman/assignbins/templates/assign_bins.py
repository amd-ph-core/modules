#!/usr/bin/env python3
"""Fold tier-2/tier-3 enrichment hits into the tier-1 kraken2 family bins.

The policy here is measured, not assumed. On a 200,000-pair subsample of an rRNA-depleted
wastewater library, with all four candidate tiers reduced to "set of families per read":

  * no tier is a superset of any other. kraken2 -- the least sensitive -- still held 1,706
    reads nothing else recovered, the dense panel 841, mmseqs 1,442. So enrichment ADDS to the
    kraken2 bins and never replaces them, and a read already binned at tier 1 keeps that bin.
  * within-tier cross-recruitment is rare: 0 reads for kraken2 and the dense panel, 1,155
    (0.8%) for mmseqs. Across the union, 1,155 of 181,772 recruited reads (0.6%) carry more
    than one family.
  * between-tier family disagreement was ZERO across all six tier pairs, including 150,222
    reads shared by mmseqs and the assembly tier. When two methods both place a read, they
    agree. That is why a simple union of family labels is safe: it is not averaging away a
    conflict, because there is no conflict to average.

Hence: union the family labels a read earned from any tier; place it if exactly one family
survives; hold it in the residual if more than one does. Holding back cross-mappers costs 0.6%
of recruited reads and keeps every bin single-family by construction, which is what downstream
consumers need in order to treat a bin as one organism group.

Reads no tier recruited are written to the residual too, never dropped -- 9.1% of the pool as
measured, which is reasonable overhead to pass along for downstream tools to mine.
"""

import argparse
import collections
import csv
import gzip
import pathlib
import sys


def norm(read_id):
    """Strip mate suffixes to get the fragment (pair) name.

    kraken2 reports one row per PAIR while aligners report per MATE, and rammap appends its own
    mate tag on top of whatever the FASTQ already carried -- names arrive as `...:1364/1/1`.
    Stripping once leaves half the pairs double-counted, so strip repeatedly.
    """
    while read_id.endswith(("/1", "/2", ".1", ".2")):
        read_id = read_id[:-2]
    return read_id


def read_bin_manifest(path):
    """bin taxid -> scientific name, from the tool's own bins.tsv.

    Joining on the taxid rather than parsing bin filenames is deliberate: the filename carries a
    sanitised scientific name that is not round-trippable, while the taxid is authoritative.
    """
    out = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\\t"):
            taxid = row.get("bin_taxid", "").strip()
            if taxid and taxid != "0":
                out[taxid] = row.get("scientific_name", "").strip()
    return out


def load_mmseqs(path, families):
    """read -> {family}, from an mmseqs m8 whose target names are `<family>|<record>|<orf>`.

    rRNA control targets are counted separately rather than discarded. They are the
    false-positive channel a six-frame translated search opens up, and the measured rate (0-2
    best hits per 2,000 pairs, ~0.1%) is only knowable if they stay visible.
    """
    hits = collections.defaultdict(set)
    controls = set()
    if not path:
        return hits, controls
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\\n").split("\\t")
            if len(f) < 2:
                continue
            fam = f[1].split("|")[0].lower()
            read = norm(f[0])
            if fam == "rrna_control":
                controls.add(read)
            elif fam in families:
                hits[read].add(fam)
    return hits, controls


def load_rammap(path, blocks, min_len, min_ident):
    """read -> {family}, from a PAF filtered to a dense panel's blocks.tsv.

    No MAPQ filter, deliberately: on a dense panel the competing alignments are other variants
    of the same organism, so MAPQ is near zero for exactly the reads the panel exists to catch.
    Identity comes from the PAF match/block columns, which measures divergence, not ambiguity.
    """
    hits = collections.defaultdict(set)
    if not path or not blocks:
        return hits
    block_family = {}
    with open(blocks, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\\t"):
            bid = row.get("block_id")
            fam = (row.get("family") or row.get("bin") or "").strip().lower()
            if bid and fam:
                block_family[bid] = fam
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\\n").split("\\t")
            if len(f) < 11:
                continue
            try:
                nmatch, alen = int(f[9]), int(f[10])
            except ValueError:
                continue
            if alen < min_len or not alen or 100.0 * nmatch / alen < min_ident:
                continue
            fam = block_family.get(f[5])
            if fam:
                hits[norm(f[0])].add(fam)
    return hits


def write_fastq(records, path):
    with gzip.open(path, "wt") as out:
        for name, seq, qual in records:
            out.write(f"@{name}\\n{seq}\\n+\\n{qual}\\n")


def iter_fastq(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            seq = fh.readline().rstrip("\\n")
            fh.readline()
            qual = fh.readline().rstrip("\\n")
            yield h[1:].rstrip("\\n"), seq, qual


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin-manifest", required=True, help="bins.tsv emitted by kraken2-hunter")
    ap.add_argument("--residual-1", required=True, help="residual mate 1 (unresolved + unclassified)")
    ap.add_argument("--residual-2", help="residual mate 2, for paired input")
    ap.add_argument("--mmseqs-hits", help="tier 2 m8, optionally gzipped")
    ap.add_argument("--rammap-paf", help="tier 3 PAF, optionally gzipped")
    ap.add_argument("--rammap-blocks", help="blocks.tsv mapping panel records to families")
    ap.add_argument("--rammap-min-aln-len", type=int, default=60)
    ap.add_argument("--rammap-min-identity", type=float, default=80.0)
    ap.add_argument(
        "--hold-cross-mapped",
        choices=["true", "false"],
        default="true",
        help="hold reads recruited to >1 family in the residual (measured cost 0.6%%)",
    )
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--outdir", default="enriched")
    a = ap.parse_args()

    families = {v.lower(): k for k, v in read_bin_manifest(a.bin_manifest).items()}
    if not families:
        # An empty manifest means tier 1 opened no bins at all. That is a real possibility on a
        # library with no viral content, and it is not an error -- but enrichment has nowhere to
        # put anything, so say so rather than emit silence.
        print("no tier-1 bins in manifest; all residual reads stay residual", file=sys.stderr)

    mm, controls = load_mmseqs(a.mmseqs_hits, set(families))
    rm = load_rammap(a.rammap_paf, a.rammap_blocks, a.rammap_min_aln_len, a.rammap_min_identity)

    union = collections.defaultdict(set)
    for src in (mm, rm):
        for read, fams in src.items():
            union[read] |= fams

    hold = a.hold_cross_mapped == "true"
    assigned = {}
    cross = set()
    for read, fams in union.items():
        if len(fams) == 1:
            assigned[read] = next(iter(fams))
        elif hold:
            cross.add(read)
        else:
            # Without the hold rule a multi-family read has to go somewhere; lexicographic
            # choice at least makes it deterministic rather than dependent on set ordering.
            assigned[read] = sorted(fams)[0]

    outdir = pathlib.Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    mates = [("1", a.residual_1)] + ([("2", a.residual_2)] if a.residual_2 else [])
    counts = collections.Counter()
    control_hits = 0
    for mate, path in mates:
        per_family = collections.defaultdict(list)
        leftover = []
        for name, seq, qual in iter_fastq(path):
            read = norm(name.split()[0])
            fam = assigned.get(read)
            if fam:
                per_family[fam].append((name, seq, qual))
                counts[fam] += 1
            else:
                leftover.append((name, seq, qual))
                counts["_residual"] += 1
                if read in cross:
                    counts["_cross_mapped"] += 1
                elif read in controls:
                    control_hits += 1
        for fam, recs in per_family.items():
            taxid = families[fam]
            write_fastq(recs, outdir / f"taxid_{taxid}_{fam}_enriched_R{mate}.fastq.gz")
        write_fastq(leftover, f"{a.prefix}.residual_{mate}.fastq.gz")

    # `_cross_mapped` is a breakdown OF `_residual`, not a sibling of it, so it must not be
    # summed into the total -- doing so reported 10 mates for an 8-mate input.
    enriched_total = sum(v for k, v in counts.items() if not k.startswith("_"))
    with open(f"{a.prefix}.assignment_report.tsv", "w") as fh:
        fh.write("metric\\tvalue\\n")
        # Every count here is MATES, not fragments. The loop walks each mate file, so a paired
        # sample reports twice its fragment count; naming the unit in the metric avoids the
        # reader having to know that. Fragment-level counts live in kraken2's own bins.tsv.
        fh.write(f"residual_mates_in\\t{enriched_total + counts['_residual']}\\n")
        fh.write(f"enriched_mates_into_bins\\t{enriched_total}\\n")
        fh.write(f"held_cross_mapped_mates\\t{counts['_cross_mapped']}\\n")
        fh.write(f"rrna_control_only_mates\\t{control_hits}\\n")
        fh.write(f"residual_mates_out\\t{counts['_residual']}\\n")
        for fam in sorted(k for k in counts if not k.startswith("_")):
            fh.write(f"enriched_mates_{fam}\\t{counts[fam]}\\n")

    print(
        f"enriched {sum(v for k, v in counts.items() if not k.startswith('_')):,} reads into "
        f"{sum(1 for k in counts if not k.startswith('_')):,} bins; "
        f"{counts['_residual']:,} stay residual "
        f"({counts['_cross_mapped']:,} held for cross-mapping)",
        file=sys.stderr,
    )
    return 0


# ---------------------------------------------------------------------------------------------
# Nextflow `template` entry point.
#
# This file is used two ways on purpose. Nextflow interpolates the placeholder references below
# and runs it with no arguments; a developer or a test can run the very same file with real
# CLI flags. The injection fires only when argv is empty, so the standalone path is unaffected
# and the policy code above has exactly one copy.
#
# Only PROCESS-level references may appear here. Variables declared with `def` in the process
# script block are not part of a template's binding, so everything the script needs is either
# a process input, something off `task`, or derived here in Python.
if __name__ == "__main__":
    import shlex

    if len(sys.argv) == 1:
        prefix = "${task.ext.prefix}" if "${task.ext.prefix}" != "null" else "${meta.id}"
        single_end = "${meta.single_end}" == "true"

        # The residual mates are globbed rather than interpolated from the `residual` input,
        # because interpolating a Groovy List yields its toString -- `[a.gz, b.gz]` -- and
        # parsing that back is worse than globbing. Globbing is safe here in a way it was not
        # in BINMAN_RESIDUAL: the producer is BINMAN_RESIDUAL itself and its names are fixed,
        # and the count is asserted against meta.single_end below rather than assumed.
        res = sorted(pathlib.Path(".").glob("*.residual_pool_*.fastq.gz"))
        want = 1 if single_end else 2
        if len(res) != want:
            sys.exit(
                f"BINMAN_ASSIGNBINS: staged {len(res)} residual mate(s) "
                f"{[p.name for p in res]}, expected {want} for single_end={single_end}"
            )

        argv = ["--bin-manifest", "${bin_manifest}", "--residual-1", str(res[0])]
        if not single_end:
            argv += ["--residual-2", str(res[1])]

        # An absent optional input interpolates to an empty string. Omitting the flag entirely
        # is not the same as passing an empty path: a disabled tier must not look like a tier
        # that ran and found nothing.
        for flag, value in (
            ("--mmseqs-hits", "${mmseqs_hits}"),
            ("--rammap-paf", "${rammap_paf}"),
            ("--rammap-blocks", "${rammap_blocks}"),
        ):
            if value and value != "null":
                argv += [flag, value]

        argv += ["--prefix", prefix, "--outdir", "enriched"]

        # Thresholds and the cross-mapping policy arrive through ext.args rather than `params`.
        # A module in a shared repo cannot read pipeline params -- a consumer that has never
        # heard of `params.rammap_min_identity` would fail at script evaluation time.
        extra = "${task.ext.args}"
        if extra and extra != "null":
            argv += shlex.split(extra)

        sys.argv += argv

    # `eval` outputs are rejected on a non-Bash interpreter -- "Process output of type 'eval'
    # is only allowed with Bash process scripts" -- so a template module reports its versions
    # by writing the file and declaring it on the topic line. It is written BEFORE main() runs:
    # a declared output that is missing because the tool exited non-zero turns one failure into
    # two, and the second one hides the first.
    import platform

    with open("versions.yml", "w", encoding="utf-8") as fh:
        fh.write('"${task.process}":\\n    python: ' + platform.python_version() + "\\n")

    sys.exit(main())
