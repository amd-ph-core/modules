#!/usr/bin/env python3
"""Check that a minimap2 index (MMI2) matches a reference FASTA.

Compares every stored reference name, length and sequence with the FASTA, in
order, and checks the index k-mer, window and flags against a minimap2 preset.
Exits non-zero with a message if anything differs.

"""

import hashlib
import os
import platform
import struct
import sys
from pathlib import Path

NT4 = bytes(
    "ACGT".find(chr(value).upper().replace("U", "T")) if chr(value).upper().replace("U", "T") in "ACGT" else 4
    for value in range(256)
)
LOW = bytes(value & 15 for value in range(256))
HIGH = bytes(value >> 4 for value in range(256))
PRESET_INDEX_SETTINGS = {
    # MMI2 headers store these fields as (window, k-mer, flags).
    "sr": (11, 21, 0),
    "asm20": (10, 19, 0),
}


def reference_fingerprints(path):
    references = {}
    name, length, digest = None, 0, None
    with open(path, "rb") as handle:
        for line in handle:
            if line.startswith(b">"):
                if name is not None:
                    references[name] = (length, digest.digest())
                name = line[1:].split()[0]
                if name in references:
                    raise ValueError(f"Duplicate reference name: {name!r}")
                length, digest = 0, hashlib.sha256()
            elif line.strip():
                if name is None:
                    raise ValueError("Expected a FASTA header")
                sequence = b"".join(line.split()).translate(NT4)
                length += len(sequence)
                digest.update(sequence)
    if name is not None:
        references[name] = (length, digest.digest())
    if not references:
        raise ValueError("Reference FASTA is empty")
    return references


def validate_index(index, fasta, preset=None):
    expected = reference_fingerprints(fasta)
    observed = []
    size = Path(index).stat().st_size
    settings = None
    with open(index, "rb") as handle:

        def read(count):
            value = handle.read(count)
            if len(value) != count:
                raise ValueError("Truncated minimap2 index")
            return value

        def skip(count):
            if handle.tell() + count > size:
                raise ValueError("Truncated minimap2 index")
            handle.seek(count, 1)

        while handle.tell() < size:
            if read(4) != b"MMI\\x02":
                raise ValueError("Expected an MMI2 index; rebuild with the pipeline's minimap2 version")
            window, kmer, buckets, count, flags = struct.unpack("<5I", read(20))
            if flags & 6:
                raise ValueError("Supplied index must retain reference names and sequences")
            if buckets > 30 or (1 << buckets) * 8 > size - handle.tell():
                raise ValueError("Invalid minimap2 index bucket count")
            current_settings = (window, kmer, flags)
            if settings is not None and settings != current_settings:
                raise ValueError("Inconsistent settings across multipart minimap2 index")
            settings = current_settings
            sequences = []
            for _ in range(count):
                name = read(read(1)[0])
                (length,) = struct.unpack("<I", read(4))
                if name not in expected or expected[name][0] != length:
                    raise ValueError(f"Index/reference name or length mismatch: {name!r}")
                sequences.append((name, length))
            # MMI2 stores positions and hash entries before its packed reference bases.
            for _ in range(1 << buckets):
                (positions,) = struct.unpack("<I", read(4))
                skip(8 * positions)
                (hashes,) = struct.unpack("<I", read(4))
                skip(16 * hashes)
            start = handle.tell()
            bases = sum(length for _, length in sequences)
            skip(((bases + 7) // 8) * 4)
            end = handle.tell()
            offset = 0
            for name, length in sequences:
                handle.seek(start + offset // 2)
                remaining, first = length, offset % 2
                digest = hashlib.sha256()
                while remaining:
                    packed = read(min(65536, (remaining + first + 1) // 2))
                    decoded = bytearray(len(packed) * 2)
                    decoded[::2], decoded[1::2] = packed.translate(LOW), packed.translate(HIGH)
                    sequence = decoded[first : first + remaining]
                    digest.update(sequence)
                    remaining -= len(sequence)
                    first = 0
                if digest.digest() != expected[name][1]:
                    raise ValueError(f"Index/reference sequence mismatch: {name!r}")
                observed.append(name)
                offset += length
            handle.seek(end)
    if observed != list(expected):
        raise ValueError("Index/reference sequence membership or order mismatch")
    if preset is not None:
        try:
            expected_settings = PRESET_INDEX_SETTINGS[preset]
        except KeyError as error:
            raise ValueError(f"Unsupported minimap2 preset: {preset!r}") from error
        if settings != expected_settings:
            expected_window, expected_kmer, expected_flags = expected_settings
            window, kmer, flags = settings
            raise ValueError(
                "Minimap2 index parameters do not match "
                f"preset {preset!r}: expected k={expected_kmer}, w={expected_window}, "
                f"flags={expected_flags}; found k={kmer}, w={window}, flags={flags}"
            )
    return settings


if __name__ == "__main__":
    index = "${index}"
    fasta = "${fasta}"
    preset = "${preset}"
    task_process = "${task.process}"

    try:
        window, kmer, flags = validate_index(index, fasta, preset)
    except ValueError as error:
        sys.exit(f"ERROR: {index} failed validation against {fasta}: {error}")
    print(f"Validated reference sequences: k={kmer}, w={window}, flags={flags}")

    # Pass the index on unchanged under validated/.
    os.makedirs("validated", exist_ok=True)
    os.symlink(os.path.join("..", index), os.path.join("validated", os.path.basename(index)))

    with open("versions.yml", "w", encoding="utf-8") as handle:
        handle.write(f'"{task_process}":\\n    python: "{platform.python_version()}"\\n')
