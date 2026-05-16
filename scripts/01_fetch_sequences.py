#!/usr/bin/env python3
"""
Step 1: Retrieve PLA2-like domain sequences (PF08398) from viral genomes via EBI InterPro API.

Output:
    data/raw/viralexport.json  — raw paginated API results
"""

import sys
import json
import ssl
import argparse
import logging
from urllib import request
from urllib.error import HTTPError
from time import sleep
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# InterPro REST API — proteins (UniProt) that carry PF08398, filtered to
# viral taxonomy (TaxID 10239).
#
# ?extra_fields=sequence  — embed the full AA sequence in each result object,
#                           so no separate FASTA download is needed downstream.
# ?page_size=200          — maximum allowed page size; reduces round-trips.
BASE_URL = (
    "https://www.ebi.ac.uk/interpro/api/protein/UniProt"
    "/entry/pfam/PF08398/taxonomy/uniprot/10239/"
    "?extra_fields=sequence&page_size=200"
)

DEFAULT_OUTPUT = Path("data/raw/viralexport.json")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Core fetch logic
# ---------------------------------------------------------------------------

def fetch_all_results(base_url: str, max_retries: int = 3, retry_wait: int = 61) -> list:
    """
    Page through the InterPro API and collect all result objects.

    Returns a list of dicts, one per protein entry.
    """
    context = ssl._create_unverified_context()
    next_url = base_url
    all_results = []
    page = 0

    while next_url:
        page += 1
        log.info(f"Fetching page {page} …")
        attempts = 0

        while True:  # retry loop for transient errors
            try:
                req = request.Request(next_url, headers={"Accept": "application/json"})
                res = request.urlopen(req, context=context)

                if res.status == 408:
                    log.warning("Timeout (408) — waiting 61 s before retry …")
                    sleep(retry_wait)
                    continue

                if res.status == 204:
                    log.info("No content (204) — end of results.")
                    return all_results

                payload = json.loads(res.read().decode())
                break  # success

            except HTTPError as exc:
                if exc.code == 408 or attempts < max_retries:
                    attempts += 1
                    log.warning(f"HTTP {exc.code} — retry {attempts}/{max_retries} …")
                    sleep(retry_wait)
                else:
                    log.error(f"Fatal HTTP error at URL: {next_url}")
                    raise

        all_results.extend(payload.get("results", []))
        log.info(f"  … {len(all_results)} entries collected so far.")

        next_url = payload.get("next")
        if next_url:
            sleep(1)  # be polite to the server

    return all_results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "-o", "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Output JSON file (default: {DEFAULT_OUTPUT})",
    )
    p.add_argument(
        "--url",
        default=BASE_URL,
        help="Override the InterPro API base URL.",
    )
    return p.parse_args()


def main():
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    log.info(f"Starting retrieval from: {args.url}")
    results = fetch_all_results(args.url)

    if not results:
        log.warning("No results returned — check the URL or your network.")
        sys.exit(1)

    with open(args.output, "w", encoding="utf-8") as fh:
        json.dump({"results": results}, fh, indent=2)

    log.info(f"Done. {len(results)} entries written to {args.output}")


if __name__ == "__main__":
    main()
