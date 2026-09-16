#!/usr/bin/env python3
"""Summarize newline-delimited cfd_solvers benchmark JSON files as Markdown."""
from __future__ import annotations

import argparse
import json
import statistics
from collections import defaultdict
from pathlib import Path


def load(paths: list[Path]) -> list[dict]:
    rows: list[dict] = []
    for path in paths:
        for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            row = json.loads(line)
            if row.get("schema") != "cfd_solvers.benchmark.v1":
                raise ValueError(f"{path}:{number}: unsupported schema")
            rows.append(row)
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="+", type=Path)
    args = parser.parse_args()
    rows = load(args.files)
    groups: dict[tuple, list[dict]] = defaultdict(list)
    for row in rows:
        key = (row["lattice"], row["backend"], row["streaming"], row["device"], row["nx"], row["ny"], row["nz"])
        groups[key].append(row)
    print("| Lattice | Backend | Streaming | Device | Grid | Samples | Median MLUPS | Best MLUPS | Median GB/s |")
    print("|---|---|---|---|---:|---:|---:|---:|---:|")
    for key in sorted(groups):
        lattice, backend, streaming, device, nx, ny, nz = key
        items = groups[key]
        mlups = [float(r["mlups"]) for r in items]
        gbps = [float(r["estimated_ddf_gbps"]) for r in items]
        print(f"| {lattice} | {backend} | {streaming} | {device} | {nx}x{ny}x{nz} | {len(items)} | "
              f"{statistics.median(mlups):.3f} | {max(mlups):.3f} | {statistics.median(gbps):.3f} |")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
