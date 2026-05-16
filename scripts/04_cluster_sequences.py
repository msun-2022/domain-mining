#!/usr/bin/env python3
"""
Step 4: Sequence clustering with CD-HIT.

Two-pass strategy
-----------------
Pass 1  Remove near-identical sequences (default identity ≥ 0.99).
Pass 2  Cluster the de-duplicated set at a lower threshold (default ≥ 0.70)
        to group related sequences and pick one representative per cluster.

Requirements
------------
    cd-hit  must be on $PATH  (conda install -c bioconda cd-hit)

Input:
    results/pfam_viralpla2_full_domain_ext100.fasta  (from Step 3)

Outputs (all under results/clustering/):
    dedup_0.99.*          — CD-HIT output after deduplication
    clusters_0.70.*       — CD-HIT output after clustering
    representative_seqs.fasta   — one representative per cluster
    cluster_summary.csv         — cluster ID, representative, size, all members
"""

import argparse
import logging
import subprocess
import sys
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

DEFAULT_INPUT    = Path("results/pfam_viralpla2_full_domain_ext100.fasta")
DEFAULT_OUTDIR   = Path("results/clustering")
DEDUP_IDENTITY   = 0.99
CLUSTER_IDENTITY = 0.70


# ---------------------------------------------------------------------------
# CD-HIT helpers
# ---------------------------------------------------------------------------

def check_cdhit():
    """Abort early if cd-hit is not available."""
    result = subprocess.run(["which", "cd-hit"], capture_output=True)
    if result.returncode != 0:
        log.error(
            "cd-hit not found on PATH.\n"
            "Install with:  conda install -c bioconda cd-hit"
        )
        sys.exit(1)
    log.info(f"cd-hit found at: {result.stdout.decode().strip()}")


def word_size_for(identity: float) -> int:
    """
    CD-HIT recommended word size:
        ≥ 0.7 → 5 | ≥ 0.6 → 4 | ≥ 0.5 → 3 | < 0.5 → 2
    """
    if identity >= 0.7:
        return 5
    elif identity >= 0.6:
        return 4
    elif identity >= 0.5:
        return 3
    return 2


def run_cdhit(
    input_fasta: Path,
    output_prefix: Path,
    identity: float,
    threads: int = 4,
    memory_mb: int = 4000,
) -> Path:
    """
    Run cd-hit and return the path to the output FASTA (representatives).
    """
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "cd-hit",
        "-i", str(input_fasta),
        "-o", str(output_prefix),
        "-c", str(identity),
        "-n", str(word_size_for(identity)),
        "-T", str(threads),
        "-M", str(memory_mb),
        "-d", "0",        # full sequence description in .clstr
        "-g", "1",        # accurate (slower) clustering mode
    ]
    log.info(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log.error(f"cd-hit failed:\n{result.stderr}")
        sys.exit(result.returncode)
    log.info(result.stdout[-500:] if len(result.stdout) > 500 else result.stdout)
    return output_prefix   # cd-hit writes <prefix> (fasta) and <prefix>.clstr


# ---------------------------------------------------------------------------
# Parse .clstr file
# ---------------------------------------------------------------------------

def parse_clstr(clstr_path: Path) -> list[dict]:
    """
    Parse a CD-HIT .clstr file.

    Returns a list of dicts:
        {
          "cluster_id": int,
          "representative": str,    # sequence ID of the * member
          "size": int,
          "members": [str, ...]
        }
    """
    clusters = []
    current: dict | None = None

    with open(clstr_path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">Cluster"):
                if current:
                    clusters.append(current)
                cluster_id = int(line.split()[1])
                current = {"cluster_id": cluster_id, "representative": None,
                           "size": 0, "members": []}
            elif current is not None and line:
                # Example member line:
                # 0    142aa, >A0A023GPI8... *
                # 1    145aa, >B3A1T2...     at 97.18%
                parts = line.split(">")
                if len(parts) < 2:
                    continue
                seq_id = parts[1].split("...")[0].strip()
                current["members"].append(seq_id)
                current["size"] += 1
                if line.endswith("*"):
                    current["representative"] = seq_id

    if current:
        clusters.append(current)

    return clusters


def write_cluster_summary(clusters: list[dict], out_path: Path):
    import csv
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["cluster_id", "representative", "size", "members"])
        for c in clusters:
            writer.writerow([c["cluster_id"], c["representative"], c["size"],
                             "; ".join(c["members"])])
    log.info(f"Cluster summary written to {out_path}")


def extract_representatives(
    fasta_path: Path,
    representatives: set[str],
    out_path: Path,
):
    """Copy only representative sequences to out_path."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    keep = False
    with open(fasta_path, encoding="utf-8") as in_fh, \
         open(out_path, "w", encoding="utf-8") as out_fh:
        for line in in_fh:
            if line.startswith(">"):
                seq_id = line[1:].split()[0].strip()
                keep = seq_id in representatives
            if keep:
                out_fh.write(line)
                if line.startswith(">"):
                    written += 1
    log.info(f"{written} representative sequences written to {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "-i", "--input",
        type=Path, default=DEFAULT_INPUT,
        help=f"Input FASTA from Step 3 (default: {DEFAULT_INPUT})",
    )
    p.add_argument(
        "--outdir",
        type=Path, default=DEFAULT_OUTDIR,
        help=f"Output directory (default: {DEFAULT_OUTDIR})",
    )
    p.add_argument(
        "--dedup-id",
        type=float, default=DEDUP_IDENTITY,
        help=f"Identity threshold for deduplication (default: {DEDUP_IDENTITY})",
    )
    p.add_argument(
        "--cluster-id",
        type=float, default=CLUSTER_IDENTITY,
        help=f"Identity threshold for clustering (default: {CLUSTER_IDENTITY})",
    )
    p.add_argument(
        "-T", "--threads",
        type=int, default=4,
        help="Number of threads for cd-hit (default: 4)",
    )
    p.add_argument(
        "-M", "--memory",
        type=int, default=4000,
        help="Memory limit in MB for cd-hit (default: 4000)",
    )
    return p.parse_args()


def main():
    args = parse_args()
    check_cdhit()

    # ---- Pass 1: deduplication ----
    dedup_id_str  = str(args.dedup_id).replace(".", "")
    dedup_prefix  = args.outdir / f"dedup_{args.dedup_id}"
    run_cdhit(args.input, dedup_prefix, args.dedup_id,
              threads=args.threads, memory_mb=args.memory)

    # ---- Pass 2: clustering ----
    cluster_id_str  = str(args.cluster_id).replace(".", "")
    cluster_prefix  = args.outdir / f"clusters_{args.cluster_id}"
    run_cdhit(dedup_prefix, cluster_prefix, args.cluster_id,
              threads=args.threads, memory_mb=args.memory)

    # ---- Summarise clusters ----
    clstr_path = Path(str(cluster_prefix) + ".clstr")
    clusters   = parse_clstr(clstr_path)
    log.info(f"Total clusters: {len(clusters)}")

    summary_path = args.outdir / "cluster_summary.csv"
    write_cluster_summary(clusters, summary_path)

    # ---- Extract representatives from dedup FASTA (more sequences to pick from) ----
    reps = {c["representative"] for c in clusters if c["representative"]}
    rep_fasta = args.outdir / "representative_seqs.fasta"
    extract_representatives(dedup_prefix, reps, rep_fasta)

    log.info(
        f"\nSummary\n"
        f"  Input sequences     : {args.input}\n"
        f"  After dedup ({args.dedup_id})  : {dedup_prefix}\n"
        f"  Clusters ({args.cluster_id})   : {len(clusters)}\n"
        f"  Representatives     : {rep_fasta}\n"
        f"  Cluster summary     : {summary_path}"
    )


if __name__ == "__main__":
    main()
