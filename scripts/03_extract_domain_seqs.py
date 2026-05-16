#!/usr/bin/env python3
"""
Step 3: Extend PLA2-like domain boundaries and extract truncated sequences.

For each protein in the boundary CSV, the aligned domain coordinates
(start, end — 1-based, inclusive) are extended by --flank residues on
both the N- and C-terminal sides, clamped to the actual sequence length.
The resulting subsequence is written to a FASTA file.

Performance notes vs the original Perl version
-----------------------------------------------
* The full-length FASTA is loaded once into a dict keyed by protein ID —
  lookups are O(1) regardless of file size.
* Multi-line sequences are joined with str.join (no per-char list appends).
* Whitespace stripping uses str.translate with a precomputed table — faster
  than str.replace or regex for large sequences.
* Output lines are accumulated in a list and written in one pass, minimising
  system-call overhead.

Inputs:
    data/processed/pfam_viralpla2_domain_ext.csv   (protein_id, start, end)
    data/processed/pfam_viralpla2_full.fasta

Output:
    results/pfam_viralpla2_full_domain_ext<N>.fasta
"""

import argparse
import logging
import sys
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

DEFAULT_CSV   = Path("data/processed/pfam_viralpla2_domain_ext.csv")
DEFAULT_FASTA = Path("data/processed/pfam_viralpla2_full.fasta")
DEFAULT_FLANK = 100

# Precomputed translation table to strip all whitespace from sequences.
_STRIP_WS = str.maketrans("", "", " \t\r\n")


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def load_boundaries(csv_path: Path) -> dict[str, tuple[int, int]]:
    """
    Parse a 3-column CSV (protein_id, start, end) with 1-based coordinates.
    Returns {protein_id: (start_1based, end_1based)}.
    Skips blank lines and lines beginning with '#'.
    An optional header row is detected and skipped automatically.
    """
    import csv
    boundaries: dict[str, tuple[int, int]] = {}

    with open(csv_path, newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh)
        for lineno, row in enumerate(reader, 1):
            if not row or row[0].startswith("#"):
                continue
            if lineno == 1 and not row[0].lstrip("-").isdigit() and not row[0][0].isalpha() or \
               (lineno == 1 and row[0].lower() in ("protein_id", "id", "accession", "protein")):
                continue  # skip header
            if len(row) < 3:
                log.warning(f"Line {lineno}: expected 3 columns, got {len(row)} — skipping.")
                continue
            try:
                boundaries[row[0].strip()] = (int(row[1].strip()), int(row[2].strip()))
            except ValueError:
                log.warning(f"Line {lineno}: non-integer coordinates — skipping.")

    log.info(f"Loaded boundaries for {len(boundaries)} proteins from {csv_path}")
    return boundaries


def load_fasta(fasta_path: Path) -> dict[str, str]:
    """
    Read a FASTA file into a dict {protein_id: sequence}.

    Only the first whitespace-delimited token of each header line is used
    as the key, matching the convention in load_boundaries.

    Sequences may be multi-line; all whitespace is stripped before storage.
    """
    sequences: dict[str, str] = {}
    current_id: str | None = None
    parts: list[str] = []

    with open(fasta_path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(parts).translate(_STRIP_WS)
                current_id = line[1:].split()[0]   # first token after '>'
                parts = []
            else:
                parts.append(line)

    # flush the final record
    if current_id is not None:
        sequences[current_id] = "".join(parts).translate(_STRIP_WS)

    log.info(f"Loaded {len(sequences)} sequences from {fasta_path}")
    return sequences


# ---------------------------------------------------------------------------
# Core extraction
# ---------------------------------------------------------------------------

def extract_domain(
    sequence: str,
    start_1: int,
    end_1: int,
    flank: int,
) -> tuple[str, int, int]:
    """
    Return the extended domain subsequence and its actual 1-based coordinates.

    Parameters
    ----------
    sequence : full-length protein sequence
    start_1  : 1-based domain start (inclusive)
    end_1    : 1-based domain end (inclusive)
    flank    : residues to add on each side (clamped to sequence boundaries)

    Returns
    -------
    (subsequence, actual_start_1based, actual_end_1based)
    """
    seq_len = len(sequence)
    ext_start = max(0, start_1 - 1 - flank)         # 0-based, inclusive
    ext_end   = min(seq_len, end_1 + flank)          # 0-based, exclusive
    return sequence[ext_start:ext_end], ext_start + 1, ext_end


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "-b", "--boundaries",
        type=Path, default=DEFAULT_CSV,
        help=f"Domain boundary CSV (default: {DEFAULT_CSV})",
    )
    p.add_argument(
        "-f", "--fasta",
        type=Path, default=DEFAULT_FASTA,
        help=f"Full-length sequence FASTA (default: {DEFAULT_FASTA})",
    )
    p.add_argument(
        "-n", "--flank",
        type=int, default=DEFAULT_FLANK,
        help=f"Residues to extend on each side (default: {DEFAULT_FLANK})",
    )
    p.add_argument(
        "-o", "--output",
        type=Path, default=None,
        help="Output FASTA path (auto-generated from input name if omitted)",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    # Auto-generate output path if not given
    if args.output is None:
        args.output = Path("results") / f"{args.fasta.stem}_domain_ext{args.flank}.fasta"

    # ---- Load inputs ----
    boundaries = load_boundaries(args.boundaries)
    if not boundaries:
        log.error("Boundary file is empty or could not be parsed.")
        sys.exit(1)

    sequences = load_fasta(args.fasta)
    if not sequences:
        log.error("FASTA file is empty or could not be parsed.")
        sys.exit(1)

    # ---- Extract and collect output lines ----
    args.output.parent.mkdir(parents=True, exist_ok=True)

    output_lines: list[str] = []
    written = skipped_no_seq = 0

    for protein_id, (start_1, end_1) in boundaries.items():
        if protein_id not in sequences:
            log.warning(f"  {protein_id}: in boundary file but not in FASTA — skipping.")
            skipped_no_seq += 1
            continue

        subseq, actual_start, actual_end = extract_domain(
            sequences[protein_id], start_1, end_1, args.flank
        )
        output_lines.append(
            f">{protein_id} "
            f"ext_start={actual_start} ext_end={actual_end} "
            f"flank={args.flank}\n"
            f"{subseq}\n"
        )
        written += 1

    # FASTA entries present but not in boundary file (informational)
    skipped_no_boundary = len(sequences) - written - skipped_no_seq

    # ---- Write output in one pass ----
    with open(args.output, "w", encoding="utf-8") as out_fh:
        out_fh.writelines(output_lines)

    log.info(
        f"Done.\n"
        f"  Sequences written  : {written}\n"
        f"  No FASTA match     : {skipped_no_seq}   (in CSV but missing from FASTA)\n"
        f"  No boundary entry  : {skipped_no_boundary}   (in FASTA but not in CSV)\n"
        f"  Output             : {args.output}"
    )


if __name__ == "__main__":
    main()
