#!/usr/bin/env python3
"""Extract FASTA/FASTQ records selected by Kraken assignments.

The ordered fast path streams Kraken output and sequence records together.  If
the assignment order does not match the sequence order, the command restarts
with a disk-backed ID index so compatibility does not require an unbounded
in-memory dictionary.
"""

from __future__ import annotations

import argparse
import contextlib
import gzip
import io
import json
import os
from pathlib import Path
import re
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import textwrap
import time
from dataclasses import dataclass
from typing import BinaryIO, Iterator, TextIO


VERSION = "1.0.0"
ARTIFICIAL_SEQUENCE_TAXID = 81077
TAXID_IN_NAME = re.compile(r"taxid\s+(\d+)\)?")


class ExtractionError(RuntimeError):
    """Base error for invalid or inconsistent extraction inputs."""


class OrderMismatch(ExtractionError):
    """Raised when assignments and sequence records are not in the same order."""


@dataclass(frozen=True)
class FastxRecord:
    identifier: str
    description: str
    sequence: str
    quality: str | None


@dataclass(frozen=True)
class Assignment:
    identifier: str
    taxid: int


@dataclass
class ExtractionStats:
    mode: str
    assignments_processed: int = 0
    records_processed: int = 0
    records_written: int = 0
    elapsed_seconds: float = 0.0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extract reads using Kraken/Kraken2 assignment output."
    )
    parser.add_argument("-k", "--kraken", "--kraken-file", dest="kraken_file", required=True)
    parser.add_argument(
        "-s", "-s1", "-1", "-U", dest="seq_file1", required=True,
        help="First FASTA/FASTQ input; may be gzip-compressed.",
    )
    parser.add_argument("-s2", "-2", dest="seq_file2", default="")
    parser.add_argument("-t", "--taxid", "--taxids", nargs="+", required=True)
    parser.add_argument("-o", "--output", dest="output_file", required=True)
    parser.add_argument("-o2", "--output2", dest="output_file2", default="")
    append_group = parser.add_mutually_exclusive_group()
    append_group.add_argument("--append", dest="append", action="store_true")
    append_group.add_argument("--noappend", dest="append", action="store_false")
    parser.set_defaults(append=False)
    parser.add_argument("--max", dest="max_reads", type=int, default=100_000_000)
    parser.add_argument("-r", "--report", "--report-file", dest="report_file", default="")
    parser.add_argument("--include-parents", dest="parents", action="store_true")
    parser.add_argument("--include-children", dest="children", action="store_true")
    parser.add_argument("--exclude", action="store_true")
    parser.add_argument("--fastq-output", dest="fastq_output", action="store_true")
    parser.add_argument(
        "--threads", type=int, default=1,
        help="Total CPU budget used for streaming and compression.",
    )
    parser.add_argument(
        "--gzip-output", action="store_true",
        help="Compress output with deterministic gzip metadata.",
    )
    parser.add_argument(
        "--status-file", default="",
        help="Optional JSON status path for workflow integration.",
    )
    parser.add_argument("--version", action="version", version=VERSION)
    return parser


def normalized_id(identifier: str) -> str:
    if identifier.endswith("/1") or identifier.endswith("/2"):
        return identifier[:-2]
    return identifier


def parse_taxid(value: str) -> int:
    value = value.strip()
    if value == "A":
        return ARTIFICIAL_SEQUENCE_TAXID
    if "taxid" in value:
        match = TAXID_IN_NAME.search(value)
        if not match:
            raise ExtractionError(f"Unable to parse taxonomy ID from {value!r}")
        value = match.group(1)
    try:
        return int(value)
    except ValueError as exc:
        raise ExtractionError(f"Invalid taxonomy ID {value!r}") from exc


def iter_assignments(path: str) -> Iterator[Assignment]:
    with open(path, "rt", encoding="utf-8", newline=None) as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 5:
                raise ExtractionError(
                    f"Malformed Kraken assignment at {path}:{line_number}: expected 5 fields"
                )
            yield Assignment(fields[1], parse_taxid(fields[2]))


def parse_report_entry(line: str) -> tuple[int, int] | None:
    fields = line.rstrip("\r\n").split("\t")
    if len(fields) < 6:
        return None
    try:
        int(fields[1])
    except ValueError:
        return None

    taxid: int | None = None
    name = fields[-1]
    for candidate in (fields[-2], fields[-3] if len(fields) >= 7 else ""):
        try:
            taxid = int(candidate)
            break
        except ValueError:
            continue
    if taxid is None:
        return None
    indentation = len(name) - len(name.lstrip(" "))
    return taxid, indentation // 2


def expanded_taxids(
    requested: set[int], report_file: str, include_parents: bool, include_children: bool
) -> set[int]:
    if not include_parents and not include_children:
        return set(requested)
    if not report_file:
        raise ExtractionError(
            "--report is required with --include-parents or --include-children"
        )

    parent: dict[int, int | None] = {}
    children: dict[int, list[int]] = {}
    stack: list[tuple[int, int]] = []
    with open(report_file, "rt", encoding="utf-8", newline=None) as handle:
        for line in handle:
            entry = parse_report_entry(line)
            if entry is None:
                continue
            taxid, level = entry
            if taxid == 0:
                continue
            while stack and stack[-1][1] >= level:
                stack.pop()
            parent_taxid = stack[-1][0] if stack else None
            parent[taxid] = parent_taxid
            if parent_taxid is not None:
                children.setdefault(parent_taxid, []).append(taxid)
            stack.append((taxid, level))

    selected = set(requested)
    if include_parents:
        for taxid in requested:
            current = parent.get(taxid)
            while current is not None:
                selected.add(current)
                current = parent.get(current)
    if include_children:
        pending = [taxid for taxid in requested if taxid in parent]
        while pending:
            current = pending.pop()
            for child in children.get(current, []):
                if child not in selected:
                    selected.add(child)
                    pending.append(child)
    return selected


class InputHandle:
    def __init__(self, path: str, use_pigz: bool):
        self.path = path
        self.use_pigz = use_pigz and path.endswith(".gz") and shutil.which("pigz") is not None
        self.process: subprocess.Popen[bytes] | None = None
        self.handle: TextIO | None = None

    def __enter__(self) -> TextIO:
        if self.use_pigz:
            self.process = subprocess.Popen(
                ["pigz", "-dc", "-p", "1", "--", self.path],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            assert self.process.stdout is not None
            self.handle = io.TextIOWrapper(self.process.stdout, encoding="utf-8", newline=None)
        elif self.path.endswith(".gz"):
            self.handle = gzip.open(self.path, "rt", encoding="utf-8", newline=None)
        else:
            self.handle = open(self.path, "rt", encoding="utf-8", newline=None)
        return self.handle

    def __exit__(self, exc_type, exc, traceback) -> None:
        if self.handle is not None:
            self.handle.close()
        if self.process is not None:
            stderr = self.process.stderr.read().decode("utf-8", errors="replace") if self.process.stderr else ""
            return_code = self.process.wait()
            if return_code not in (0, -13, 141) and exc_type is None:
                raise ExtractionError(f"pigz failed while reading {self.path}: {stderr.strip()}")


def detect_format(path: str, use_pigz: bool) -> str:
    with InputHandle(path, use_pigz) as handle:
        first = handle.read(1)
    if first == ">":
        return "fasta"
    if first == "@":
        return "fastq"
    if not first:
        raise ExtractionError(f"Sequence file is empty: {path}")
    raise ExtractionError(f"Sequence file must be FASTA or FASTQ: {path}")


def iter_fasta(handle: TextIO, path: str) -> Iterator[FastxRecord]:
    description: str | None = None
    sequence_lines: list[str] = []
    for line_number, line in enumerate(handle, start=1):
        line = line.rstrip("\r\n")
        if line.startswith(">"):
            if description is not None:
                sequence = "".join(sequence_lines)
                identifier = description.split(None, 1)[0]
                yield FastxRecord(identifier, description, sequence, None)
            description = line[1:]
            if not description:
                raise ExtractionError(f"Blank FASTA header at {path}:{line_number}")
            sequence_lines = []
        elif description is None:
            raise ExtractionError(f"FASTA sequence precedes first header at {path}:{line_number}")
        else:
            sequence_lines.append("".join(line.split()))
    if description is not None:
        sequence = "".join(sequence_lines)
        identifier = description.split(None, 1)[0]
        yield FastxRecord(identifier, description, sequence, None)


def iter_fastq(handle: TextIO, path: str) -> Iterator[FastxRecord]:
    line_number = 0
    while True:
        header = handle.readline()
        if not header:
            return
        line_number += 1
        header = header.rstrip("\r\n")
        if not header.startswith("@"):
            raise ExtractionError(f"Invalid FASTQ header at {path}:{line_number}")
        description = header[1:]
        if not description:
            raise ExtractionError(f"Blank FASTQ header at {path}:{line_number}")

        sequence_lines: list[str] = []
        while True:
            line = handle.readline()
            if not line:
                raise ExtractionError(f"Truncated FASTQ sequence for {description!r} in {path}")
            line_number += 1
            line = line.rstrip("\r\n")
            if line.startswith("+"):
                break
            sequence_lines.append(line)
        sequence = "".join(sequence_lines)

        quality_lines: list[str] = []
        quality_length = 0
        while quality_length < len(sequence):
            line = handle.readline()
            if not line:
                raise ExtractionError(f"Truncated FASTQ quality for {description!r} in {path}")
            line_number += 1
            quality_line = line.rstrip("\r\n")
            quality_lines.append(quality_line)
            quality_length += len(quality_line)
        quality = "".join(quality_lines)
        if len(quality) != len(sequence):
            raise ExtractionError(
                f"Sequence and quality lengths differ for {description!r} in {path}"
            )
        identifier = description.split(None, 1)[0]
        yield FastxRecord(identifier, description, sequence, quality)


@contextlib.contextmanager
def fastx_records(path: str, file_format: str, use_pigz: bool) -> Iterator[Iterator[FastxRecord]]:
    with InputHandle(path, use_pigz) as handle:
        iterator = iter_fastq(handle, path) if file_format == "fastq" else iter_fasta(handle, path)
        yield iterator


class OutputHandle:
    def __init__(self, path: Path, gzip_output: bool, threads: int, use_pigz: bool):
        self.path = path
        self.gzip_output = gzip_output
        self.threads = max(1, threads)
        self.use_pigz = use_pigz
        self.file_handle: BinaryIO | None = None
        self.process: subprocess.Popen[bytes] | None = None
        self.handle: BinaryIO | None = None

    def __enter__(self) -> BinaryIO:
        self.file_handle = open(self.path, "wb")
        if self.gzip_output and self.use_pigz and shutil.which("pigz") is not None:
            self.process = subprocess.Popen(
                ["pigz", "-n", "-p", str(self.threads), "-c"],
                stdin=subprocess.PIPE,
                stdout=self.file_handle,
                stderr=subprocess.PIPE,
            )
            assert self.process.stdin is not None
            self.handle = self.process.stdin
        elif self.gzip_output:
            self.handle = gzip.GzipFile(fileobj=self.file_handle, mode="wb", mtime=0)
        else:
            self.handle = self.file_handle
        return self.handle

    def __exit__(self, exc_type, exc, traceback) -> None:
        if self.handle is not None and self.handle is not self.file_handle:
            self.handle.close()
        if self.process is not None:
            stderr = self.process.stderr.read().decode("utf-8", errors="replace") if self.process.stderr else ""
            return_code = self.process.wait()
            if return_code != 0 and exc_type is None:
                raise ExtractionError(f"pigz failed while writing {self.path}: {stderr.strip()}")
        if self.file_handle is not None:
            self.file_handle.close()


def render_record(record: FastxRecord, output_format: str) -> bytes:
    if output_format == "fastq":
        if record.quality is None:
            raise ExtractionError("FASTQ output requires FASTQ input")
        text = f"@{record.description}\n{record.sequence}\n+\n{record.quality}\n"
    else:
        wrapped = "\n".join(textwrap.wrap(record.sequence, width=60))
        text = f">{record.description}\n{wrapped}\n"
    return text.encode("utf-8")


def should_write(taxid: int, selected_taxids: set[int], exclude: bool) -> bool:
    matched = taxid in selected_taxids
    return not matched if exclude else matched


@contextlib.contextmanager
def staged_outputs(
    output_paths: list[str], gzip_output: bool, threads: int, input_workers: int, append: bool
) -> Iterator[list[BinaryIO]]:
    first_parent = Path(output_paths[0]).resolve().parent
    first_parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="kraken-extract-", dir=first_parent) as tmpdir:
        staged = [Path(tmpdir) / f"output-{index}" for index in range(len(output_paths))]
        available_output_threads = max(0, max(1, threads) - 1 - input_workers)
        use_pigz = available_output_threads >= len(output_paths)
        per_output_threads = max(1, available_output_threads // len(output_paths))
        with contextlib.ExitStack() as stack:
            handles = [
                stack.enter_context(OutputHandle(path, gzip_output, per_output_threads, use_pigz))
                for path in staged
            ]
            yield handles
        commit_outputs(staged, output_paths, append)


def commit_outputs(staged: list[Path], destinations: list[str], append: bool) -> None:
    replacement_paths: list[tuple[Path, Path]] = []
    committed: list[tuple[Path, Path | None]] = []
    try:
        for index, (source, destination_value) in enumerate(zip(staged, destinations)):
            destination = Path(destination_value)
            destination.parent.mkdir(parents=True, exist_ok=True)
            if append:
                replacement = destination.with_name(f".{destination.name}.new-{os.getpid()}-{index}")
                with open(replacement, "wb") as output_handle:
                    if destination.exists():
                        with open(destination, "rb") as existing_handle:
                            shutil.copyfileobj(existing_handle, output_handle)
                    with open(source, "rb") as source_handle:
                        shutil.copyfileobj(source_handle, output_handle)
            else:
                # The staged directory is deliberately created beside the outputs.
                # Rename the completed file into place instead of copying every byte
                # a second time. os.replace remains atomic on that filesystem.
                replacement = source
            replacement_paths.append((replacement, destination))
        try:
            for index, (replacement, destination) in enumerate(replacement_paths):
                backup = destination.with_name(f".{destination.name}.old-{os.getpid()}-{index}")
                if destination.exists():
                    os.replace(destination, backup)
                else:
                    backup = None
                try:
                    os.replace(replacement, destination)
                except BaseException:
                    if backup is not None:
                        os.replace(backup, destination)
                    raise
                committed.append((destination, backup))
        except BaseException:
            for destination, backup in reversed(committed):
                with contextlib.suppress(FileNotFoundError):
                    destination.unlink()
                if backup is not None:
                    os.replace(backup, destination)
            raise
        else:
            for _, backup in committed:
                if backup is not None:
                    backup.unlink()
    finally:
        # Append-mode replacements are outside the temporary staging directory.
        # Non-append staged files are either renamed or removed with that directory.
        if append:
            for replacement, _ in replacement_paths:
                with contextlib.suppress(FileNotFoundError):
                    replacement.unlink()


def validate_pair(record1: FastxRecord, record2: FastxRecord) -> str:
    id1 = normalized_id(record1.identifier)
    id2 = normalized_id(record2.identifier)
    if id1 != id2:
        raise ExtractionError(
            f"Paired input IDs do not match: {record1.identifier!r} and {record2.identifier!r}"
        )
    return id1


def ordered_extract(
    args: argparse.Namespace,
    selected_taxids: set[int],
    input_format: str,
    output_format: str,
    destinations: list[str],
) -> ExtractionStats:
    stats = ExtractionStats(mode="ordered")
    input_count = 2 if args.seq_file2 else 1
    output_count = len(destinations)
    use_pigz = args.threads >= 1 + input_count + output_count
    input_workers = input_count if use_pigz else 0
    with fastx_records(args.seq_file1, input_format, use_pigz) as records1:
        second_context = (
            fastx_records(args.seq_file2, input_format, use_pigz)
            if args.seq_file2
            else contextlib.nullcontext(None)
        )
        with second_context as records2:
            with staged_outputs(
                destinations, args.gzip_output, args.threads, input_workers, args.append
            ) as outputs:
                for assignment in iter_assignments(args.kraken_file):
                    try:
                        record1 = next(records1)
                    except StopIteration as exc:
                        raise OrderMismatch("Kraken output contains more assignments than sequence input") from exc
                    record2 = None
                    record_id = normalized_id(record1.identifier)
                    if records2 is not None:
                        try:
                            record2 = next(records2)
                        except StopIteration as exc:
                            raise ExtractionError("Second paired input ended before the first") from exc
                        record_id = validate_pair(record1, record2)
                    if normalized_id(assignment.identifier) != record_id:
                        raise OrderMismatch(
                            f"Assignment {assignment.identifier!r} does not match sequence {record1.identifier!r}"
                        )

                    stats.assignments_processed += 1
                    stats.records_processed += 1
                    if should_write(assignment.taxid, selected_taxids, args.exclude):
                        outputs[0].write(render_record(record1, output_format))
                        if record2 is not None:
                            outputs[1].write(render_record(record2, output_format))
                        stats.records_written += 1
                        if stats.records_written >= args.max_reads:
                            break
                else:
                    try:
                        next(records1)
                    except StopIteration:
                        pass
                    else:
                        raise OrderMismatch("Sequence input contains records without Kraken assignments")
                    if records2 is not None:
                        try:
                            next(records2)
                        except StopIteration:
                            pass
                        else:
                            raise ExtractionError("Second paired input contains extra records")
    return stats


def build_disk_index(
    database_path: Path,
    args: argparse.Namespace,
    selected_taxids: set[int],
) -> tuple[sqlite3.Connection, int, int]:
    connection = sqlite3.connect(database_path)
    connection.execute("PRAGMA journal_mode=OFF")
    connection.execute("PRAGMA synchronous=OFF")
    connection.execute("CREATE TABLE selected_reads (id TEXT PRIMARY KEY) WITHOUT ROWID")
    assignments_processed = 0
    selected_count = 0
    for assignment in iter_assignments(args.kraken_file):
        assignments_processed += 1
        if not should_write(assignment.taxid, selected_taxids, args.exclude):
            continue
        cursor = connection.execute(
            "INSERT OR IGNORE INTO selected_reads(id) VALUES (?)",
            (normalized_id(assignment.identifier),),
        )
        if cursor.rowcount:
            selected_count += 1
            if selected_count >= args.max_reads:
                break
    connection.commit()
    return connection, assignments_processed, selected_count


def fallback_extract(
    args: argparse.Namespace,
    selected_taxids: set[int],
    input_format: str,
    output_format: str,
    destinations: list[str],
) -> ExtractionStats:
    stats = ExtractionStats(mode="disk-index")
    input_count = 2 if args.seq_file2 else 1
    output_count = len(destinations)
    use_pigz = args.threads >= 1 + input_count + output_count
    input_workers = input_count if use_pigz else 0
    first_parent = Path(destinations[0]).resolve().parent
    with tempfile.TemporaryDirectory(prefix="kraken-index-", dir=first_parent) as tmpdir:
        connection, stats.assignments_processed, selected_count = build_disk_index(
            Path(tmpdir) / "selected.sqlite3", args, selected_taxids
        )
        try:
            with fastx_records(args.seq_file1, input_format, use_pigz) as records1:
                second_context = (
                    fastx_records(args.seq_file2, input_format, use_pigz)
                    if args.seq_file2
                    else contextlib.nullcontext(None)
                )
                with second_context as records2:
                    with staged_outputs(
                        destinations, args.gzip_output, args.threads, input_workers, args.append
                    ) as outputs:
                        for record1 in records1:
                            record2 = None
                            record_id = normalized_id(record1.identifier)
                            if records2 is not None:
                                try:
                                    record2 = next(records2)
                                except StopIteration as exc:
                                    raise ExtractionError("Second paired input ended before the first") from exc
                                record_id = validate_pair(record1, record2)
                            stats.records_processed += 1
                            found = connection.execute(
                                "SELECT 1 FROM selected_reads WHERE id = ?", (record_id,)
                            ).fetchone()
                            if found:
                                outputs[0].write(render_record(record1, output_format))
                                if record2 is not None:
                                    outputs[1].write(render_record(record2, output_format))
                                stats.records_written += 1
                                if stats.records_written >= selected_count:
                                    break
                        if records2 is not None and stats.records_written < selected_count:
                            try:
                                next(records2)
                            except StopIteration:
                                pass
                            else:
                                raise ExtractionError("Second paired input contains extra records")
        finally:
            connection.close()
    return stats


def validate_args(args: argparse.Namespace) -> None:
    if args.threads < 1:
        raise ExtractionError("--threads must be at least 1")
    if args.max_reads < 0:
        raise ExtractionError("--max must not be negative")
    if args.seq_file2 and not args.output_file2:
        raise ExtractionError("Must specify -o2/--output2 for paired input")
    if args.output_file2 and not args.seq_file2:
        raise ExtractionError("Second output was supplied without a second sequence input")
    if args.max_reads == 0:
        args.max_reads = 0


def write_status(path: str, stats: ExtractionStats) -> None:
    if not path:
        return
    payload = {
        "version": VERSION,
        "mode": stats.mode,
        "assignments_processed": stats.assignments_processed,
        "records_processed": stats.records_processed,
        "records_written": stats.records_written,
        "reads_present": stats.records_written > 0,
        "elapsed_seconds": round(stats.elapsed_seconds, 6),
    }
    status_path = Path(path)
    temporary = status_path.with_name(f".{status_path.name}.tmp-{os.getpid()}")
    temporary.write_text(json.dumps(payload, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, status_path)


def run(args: argparse.Namespace) -> ExtractionStats:
    validate_args(args)
    started = time.monotonic()
    requested_taxids = {parse_taxid(value) for value in args.taxid}
    selected_taxids = expanded_taxids(
        requested_taxids, args.report_file, args.parents, args.children
    )
    input_format = detect_format(args.seq_file1, args.threads > 1)
    if args.seq_file2:
        second_format = detect_format(args.seq_file2, args.threads > 1)
        if second_format != input_format:
            raise ExtractionError("Paired inputs use different sequence formats")
    if args.fastq_output and input_format != "fastq":
        raise ExtractionError("--fastq-output requires FASTQ input")
    output_format = "fastq" if args.fastq_output else "fasta"
    destinations = [args.output_file] + ([args.output_file2] if args.seq_file2 else [])

    if args.max_reads == 0:
        with staged_outputs(
            destinations, args.gzip_output, args.threads, 0, args.append
        ):
            pass
        stats = ExtractionStats(mode="ordered")
    else:
        try:
            stats = ordered_extract(
                args, selected_taxids, input_format, output_format, destinations
            )
        except OrderMismatch as exc:
            print(
                f"WARNING: {exc}; restarting with disk-backed compatibility mode",
                file=sys.stderr,
            )
            stats = fallback_extract(
                args, selected_taxids, input_format, output_format, destinations
            )
    stats.elapsed_seconds = time.monotonic() - started
    return stats


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    try:
        stats = run(args)
        write_status(args.status_file, stats)
    except (ExtractionError, OSError, sqlite3.Error) as exc:
        parser.exit(1, f"ERROR: {exc}\n")
    print(
        f"Extracted {stats.records_written} records from {stats.records_processed} "
        f"using {stats.mode} mode in {stats.elapsed_seconds:.2f}s"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

