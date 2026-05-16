#!/usr/bin/env python3
"""
Utility: inspect the raw API JSON to locate the sequence field.
Run once to debug, then discard.

Usage:
    python scripts/probe_json.py data/raw/viralexport.json
"""

import json
import sys


def find_key(obj, target: str, path: str = "root", max_depth: int = 8):
    """Recursively find all keys whose name contains `target`."""
    if max_depth == 0:
        return
    if isinstance(obj, dict):
        for k, v in obj.items():
            here = f"{path}.{k}"
            if target.lower() in k.lower():
                if isinstance(v, str):
                    preview = v[:120] + ("…" if len(v) > 120 else "")
                elif isinstance(v, (dict, list)):
                    preview = f"<{type(v).__name__}, len={len(v)}>"
                else:
                    preview = repr(v)
                print(f"  {here}\n    = {preview}\n")
            find_key(v, target, here, max_depth - 1)
    elif isinstance(obj, list):
        for i, item in enumerate(obj[:2]):
            find_key(item, target, f"{path}[{i}]", max_depth - 1)


def main():
    if len(sys.argv) < 2:
        print("Usage: python scripts/probe_json.py data/raw/viralexport.json")
        sys.exit(1)

    with open(sys.argv[1], encoding="utf-8") as fh:
        raw = json.load(fh)

    entries = raw.get("results", raw) if isinstance(raw, dict) else raw
    entry = entries[0]

    print(f"Total entries in file: {len(entries)}\n")
    print(f"Top-level keys of first entry:\n  {list(entry.keys())}\n")
    print(f"metadata keys:\n  {list(entry.get('metadata', {}).keys())}\n")

    print("=" * 60)
    print("All keys/paths containing 'seq':")
    print("=" * 60)
    find_key(entry, "seq")

    print("=" * 60)
    print("All keys/paths containing 'extra':")
    print("=" * 60)
    find_key(entry, "extra")


if __name__ == "__main__":
    main()
