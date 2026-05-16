# Functional Domain Mining Pipeline

A reproducible pipeline for retrieving, characterizing, and clustering short functional domain sequences from public protein databases. Originally developed as a case study on viral phospholipase A2-like (PLA2-like) domains, but designed to be directly reusable for any Pfam domain and any taxonomic group.

---

## Background & Motivation

Annotated domain boundaries in databases like Pfam/InterPro are derived from profile HMMs trained on known homologs. For small functional domains, these boundaries are a useful starting point but rarely the final word. Depending on the goal, you may need to:

- **extend the boundaries** to capture flanking structural elements (helices, loops) that contribute to fold stability or function — as done here with the `--flank` parameter
- **trim the boundaries** to focus on a catalytic core or a minimal binding motif
- **validate and refine** boundaries through structural alignment, mutagenesis data, or co-evolution analysis

This pipeline automates the retrieval-to-clustering leg of that workflow, producing a clean, non-redundant sequence set ready for MSA, phylogenetic analysis, or protein engineering.

### Case study: viral PLA2-like domains (PF08398)

Viral phospholipase A2-like domains are small (~150 aa) membrane-active domains found across RNA and DNA viruses, implicated in membrane disruption during viral entry or replication. Their compactness and functional potency make them interesting scaffolds for engineering membrane-permeabilising peptides and lipid-binding modules. This pipeline was first used to retrieve all annotated viral homologs (3,117 proteins from InterPro), extract and extend the domain regions, and produce a diverse representative set for MSA.

### Other domains this pipeline applies to directly

The only change needed to retarget the pipeline is the Pfam accession in `01_fetch_sequences.py` and the taxonomy filter (TaxID). Domains where this workflow is particularly relevant include:

- **Other small enzymatic domains** — lipases, proteases, nucleases where the annotated boundary excludes a catalytic loop
- **Antimicrobial peptide precursor domains** — e.g. defensins (PF00323), where the mature peptide boundaries differ from the Pfam match
- **Toxin/effector domains** — bacterial or viral effectors where flanking regions affect secretion or delivery
- **Binding modules** — SH2, SH3, PDZ, WW domains where the flanking linkers influence specificity
- **Coiled-coil and dimerisation domains** — where a few extra residues determine oligomeric state
- **Repeat domains** — ankyrins, TPRs, WD40s, where repeat boundaries inform engineering of consensus sequences

### Protein engineering applications

A common use case is generating **truncation libraries**: systematically varying the N- and C-terminal extension around an annotated domain to identify the minimal functional unit or to improve soluble expression. The `--flank` parameter in Step 3 directly supports this — running the pipeline multiple times with different flank values (e.g. 0, 25, 50, 100, 150) produces a set of truncation variants across the full sequence diversity, which can be used to:

- design expression constructs covering a range of boundary hypotheses
- identify positions that are conserved at the boundary (informative for where to cut)
- train or benchmark boundary-prediction tools against experimental expression/activity data

---

## Pipeline Overview

```
[EBI InterPro API]
       │
       ▼
01_fetch_sequences.py   →  data/raw/viralexport.json
       │
       ▼
02_json_to_csv.py       →  data/processed/pfam_sequences_flattened.csv
                        →  data/processed/pfam_viralpla2_domain_ext.csv
                        →  data/processed/pfam_viralpla2_full.fasta
       │
       ▼
03_extract_domain_seqs.py  →  results/pfam_viralpla2_full_domain_ext<N>.fasta
       │
       ▼
04_cluster_sequences.py    →  results/clustering/
                               ├── dedup_0.99.*
                               ├── clusters_0.70.*
                               ├── representative_seqs.fasta
                               └── cluster_summary.csv
```

| Step | Script | Description |
|------|--------|-------------|
| 1 | `01_fetch_sequences.py` | Paginate EBI InterPro API; collect all hits for a given Pfam + taxonomy |
| 2 | `02_json_to_csv.py` | Flatten JSON → metadata CSV; extract domain coordinates → boundary CSV; write full-length FASTA |
| 3 | `03_extract_domain_seqs.py` | Extend domain boundaries by ±N aa; write truncated FASTA |
| 4 | `04_cluster_sequences.py` | Two-pass CD-HIT: dedup at 0.99, cluster at 0.70 |

All intermediate files are generated automatically — no manual preparation needed between steps.

---

## Directory Structure

```
functional_domain_pipeline/
├── scripts/
│   ├── 01_fetch_sequences.py
│   ├── 02_json_to_csv.py
│   ├── 03_extract_domain_seqs.py
│   ├── 04_cluster_sequences.py
│   └── probe_json.py          # debug utility
├── data/
│   ├── raw/                   # API output JSON
│   └── processed/             # CSVs + full-length FASTA
├── results/
│   └── clustering/            # CD-HIT output + summaries
├── docs/
├── run_pipeline.sh
├── requirements.txt
└── README.md
```

> **Note:** `data/` and `results/` are listed in `.gitignore`. Commit only scripts and docs; large sequence files should be stored separately (e.g. Zenodo, OSF, or a private data store).

---

## Requirements

- Python ≥ 3.10 (standard library only — no `pip install` needed)
- [CD-HIT](https://github.com/weizhongli/cdhit) for Step 4

```bash
conda install -c bioconda cd-hit
```

---

## Quick Start

### Run the full pipeline

```bash
chmod +x run_pipeline.sh
./run_pipeline.sh
```

### Customise parameters

```bash
# Extend domain boundaries by 150 aa, cluster at 75% identity, use 8 threads
./run_pipeline.sh --flank 150 --cluster-id 0.75 --threads 8

# Dry-run: print commands without executing
./run_pipeline.sh --dry-run
```

### Retarget to a different domain / taxonomy

Edit the `BASE_URL` in `scripts/01_fetch_sequences.py`:

```python
# Example: human SH3 domains (PF00018), Homo sapiens TaxID 9606
BASE_URL = (
    "https://www.ebi.ac.uk/interpro/api/protein/UniProt"
    "/entry/pfam/PF00018/taxonomy/uniprot/9606/"
    "?extra_fields=sequence&page_size=200"
)
```

No other changes are needed — the rest of the pipeline is domain-agnostic.

### Run individual steps

```bash
# Step 1 — fetch
python scripts/01_fetch_sequences.py --output data/raw/viralexport.json

# Step 2 — process JSON → CSV + FASTA
python scripts/02_json_to_csv.py \
    --input  data/raw/viralexport.json \
    --output data/processed/pfam_sequences_flattened.csv

# Step 3 — extract extended domain sequences (flank = 100 aa)
python scripts/03_extract_domain_seqs.py \
    --boundaries data/processed/pfam_viralpla2_domain_ext.csv \
    --fasta      data/processed/pfam_viralpla2_full.fasta \
    --flank      100

# Step 4 — cluster
python scripts/04_cluster_sequences.py \
    --input      results/pfam_viralpla2_full_domain_ext100.fasta \
    --dedup-id   0.99 \
    --cluster-id 0.70 \
    --threads    4
```

---

## Outputs

| File | Description |
|------|-------------|
| `data/raw/viralexport.json` | Raw paginated API response |
| `data/processed/pfam_sequences_flattened.csv` | Full metadata table, one row per protein |
| `data/processed/pfam_viralpla2_domain_ext.csv` | Domain boundary file (`protein_id, start, end`) |
| `data/processed/pfam_viralpla2_full.fasta` | Full-length protein sequences |
| `results/pfam_viralpla2_full_domain_ext<N>.fasta` | Extended domain sequences (truncates) |
| `results/clustering/dedup_0.99.*` | CD-HIT output after deduplication at 99% identity |
| `results/clustering/clusters_0.70.*` | CD-HIT output after clustering at 70% identity |
| `results/clustering/representative_seqs.fasta` | One representative per cluster |
| `results/clustering/cluster_summary.csv` | Cluster ID, representative, size, all members |

---

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--flank` | `100` | Residues to extend on each side of the annotated domain boundary |
| `--dedup-id` | `0.99` | Identity threshold for removing near-identical sequences |
| `--cluster-id` | `0.70` | Identity threshold for diversity clustering |
| `--threads` | `4` | CD-HIT threads |

The `--flank` value is the most experiment-specific parameter. As a rough guide: 0 gives the raw annotated domain; 25–50 captures immediate secondary structure context; 100–150 is appropriate when the domain is part of a larger functional unit or when downstream structural modelling is planned.

---

## Acknowledgement

Ann Gregory (ex-Aera) for sharing the initial ideas of the pipeline. 
