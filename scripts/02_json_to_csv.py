#!/usr/bin/env python3
"""
Step 2: Flatten the raw InterPro JSON and export two CSVs and a FASTA.

Outputs
-------
pfam_sequences_flattened.csv
    One row per protein; every nested JSON field as a top-level column.

pfam_viralpla2_domain_ext.csv
    Three-column boundary file (protein_id, start, end) consumed by Step 3.
    Coordinates are 1-based, inclusive, taken from the first PF08398 domain
    fragment reported under entry_protein_locations for each protein.

pfam_viralpla2_full.fasta
    Full-length protein sequences in FASTA format, consumed by Step 3.
    Requires Step 1 to have been run with ?extra_fields=sequence in the URL.

All three files are written to the same output directory.

Input:
    data/raw/viralexport.json

Outputs:
    data/processed/pfam_sequences_flattened.csv
    data/processed/pfam_viralpla2_domain_ext.csv
    data/processed/pfam_viralpla2_full.fasta
"""

import json
import csv
import argparse
import logging
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

DEFAULT_INPUT  = Path("data/raw/viralexport.json")
DEFAULT_OUTPUT = Path("data/processed/pfam_sequences_flattened.csv")


# ---------------------------------------------------------------------------
# Flattening
# ---------------------------------------------------------------------------

def flatten(obj: dict, prefix: str = "") -> dict:
    """
    Recursively flatten a nested dict/list structure into a single-level dict.

    - Nested dicts  →  keys joined with '_'
    - Lists of dicts  →  each index appended to the key prefix
    - Lists of scalars  →  joined into a comma-separated string
    """
    out = {}
    for key, value in obj.items():
        full_key = f"{prefix}{key}" if prefix else key

        if isinstance(value, dict):
            out.update(flatten(value, full_key + "_"))

        elif isinstance(value, list):
            if value and isinstance(value[0], dict):
                for idx, item in enumerate(value):
                    out.update(flatten(item, f"{full_key}{idx}_"))
            else:
                out[full_key] = ", ".join(map(str, value))

        else:
            out[full_key] = value

    return out


# ---------------------------------------------------------------------------
# Boundary extraction
# ---------------------------------------------------------------------------

def extract_boundaries(entry: dict) -> tuple[str, int, int] | None:
    """
    Extract (protein_id, start, end) for the first PF08398 domain fragment
    from a single API result object.

    The InterPro protein-centric response embeds domain locations under:
        entry["entries"][i]["entry_protein_locations"][j]["fragments"][k]
            {"start": <int>, "end": <int>, ...}

    Only the first fragment of the first matching entry is used; proteins with
    no parseable location are skipped with a warning.
    """
    protein_id: str = (
        entry.get("metadata", {}).get("accession", "")
        or entry.get("accession", "")
    )
    if not protein_id:
        return None

    for entry_hit in entry.get("entries", []):
        locations = entry_hit.get("entry_protein_locations", [])
        for loc in locations:
            fragments = loc.get("fragments", [])
            if fragments:
                frag = fragments[0]
                try:
                    return protein_id, int(frag["start"]), int(frag["end"])
                except (KeyError, TypeError, ValueError):
                    continue

    log.warning(f"{protein_id}: no domain coordinates found — excluded from boundary CSV.")
    return None


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("-i", "--input",  type=Path, default=DEFAULT_INPUT,
                   help=f"Raw JSON from Step 1 (default: {DEFAULT_INPUT})")
    p.add_argument("-o", "--output", type=Path, default=DEFAULT_OUTPUT,
                   help=f"Flattened CSV output (default: {DEFAULT_OUTPUT})")
    return p.parse_args()


def main():
    args = parse_args()

    log.info(f"Loading {args.input} …")
    with open(args.input, encoding="utf-8") as fh:
        raw = json.load(fh)

    # Support both {"results": [...]} wrapper and a bare list
    entries = raw.get("results", raw) if isinstance(raw, dict) else raw

    if not entries:
        log.error("No entries found in the JSON file.")
        raise SystemExit(1)

    log.info(f"Processing {len(entries)} entries …")

    flat_rows  = []
    boundaries = []
    fasta_records = []   # list of (protein_id, sequence)

    for entry in entries:
        flat_rows.append(flatten(entry))

        coords = extract_boundaries(entry)
        if coords:
            boundaries.append(coords)

        protein_id = (
            entry.get("metadata", {}).get("accession", "")
            or entry.get("accession", "")
        )
        sequence = entry.get("extra_fields", {}).get("sequence", "")
        if protein_id and sequence:
            fasta_records.append((protein_id, sequence))

    # ---- Write flattened CSV ----
    all_keys = list(dict.fromkeys(k for row in flat_rows for k in row))
    args.output.parent.mkdir(parents=True, exist_ok=True)

    with open(args.output, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=all_keys, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(flat_rows)
    log.info(f"Flattened CSV: {len(flat_rows)} rows → {args.output}")

    # ---- Write boundary CSV ----
    boundary_path = args.output.parent / "pfam_viralpla2_domain_ext.csv"
    with open(boundary_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["protein_id", "start", "end"])
        writer.writerows(boundaries)
    log.info(f"Boundary CSV:  {len(boundaries)} proteins → {boundary_path}")

    # ---- Write full-length FASTA ----
    fasta_path = args.output.parent / "pfam_viralpla2_full.fasta"
    with open(fasta_path, "w", encoding="utf-8") as fh:
        for protein_id, sequence in fasta_records:
            fh.write(f">{protein_id}\n{sequence}\n")
    log.info(f"FASTA:         {len(fasta_records)} sequences → {fasta_path}")

    if not fasta_records:
        log.warning(
            "No sequences were written to FASTA. "
            "Make sure Step 1 was run with ?extra_fields=sequence in the URL."
        )
    if len(boundaries) < len(flat_rows):
        log.warning(
            f"{len(flat_rows) - len(boundaries)} proteins had no parseable domain "
            f"coordinates and were excluded from the boundary CSV."
        )


if __name__ == "__main__":
    main()
