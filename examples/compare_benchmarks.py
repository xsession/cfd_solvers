#!/usr/bin/env python3
"""Compare two cfd_solvers practical benchmark-suite JSON files."""
from __future__ import annotations

import argparse
import json
import math
import pathlib
import statistics


def load(path: pathlib.Path):
    data = json.loads(path.read_text())
    groups: dict[tuple[str, str, str], list[dict]] = {}
    for record in data["results"]:
        key = (record["feature_set"], record["case"], record["backend"])
        groups.setdefault(key, []).append(record)
    return groups, data


def med(records: list[dict], key: str) -> float:
    return statistics.median(float(r[key]) for r in records)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("baseline", type=pathlib.Path)
    ap.add_argument("candidate", type=pathlib.Path)
    ap.add_argument("--fail-on-workload-change", action="store_true")
    args = ap.parse_args()
    a, ad = load(args.baseline)
    b, bd = load(args.candidate)
    print(f"baseline:  {ad['timestamp_utc']}  {ad['host']['cpu']}  mode={ad.get('mode')}")
    print(f"candidate: {bd['timestamp_utc']}  {bd['host']['cpu']}  mode={bd.get('mode')}")
    if ad.get("mode") != bd.get("mode"):
        print("WARNING: benchmark modes differ; timings may not be directly comparable")
    if ad.get("build") != bd.get("build"):
        print("WARNING: build metadata differs; inspect compiler/backend flags before attributing speedups to code changes")
    print(f"{'case':38s} {'old ms':>11s} {'new ms':>11s} {'speedup':>9s} {'delta':>9s} {'workload':>10s}")
    print('-' * 96)
    changed = False
    for key in sorted(set(a) & set(b)):
        ar, br = a[key], b[key]
        old, new = med(ar, "simulation_ms"), med(br, "simulation_ms")
        speed = old / new if new else float("inf")
        delta = (new / old - 1.0) * 100.0 if old else 0.0
        same_shape = ar[0]["problem_size"] == br[0]["problem_size"] and ar[0]["steps"] == br[0]["steps"]
        old_sum, new_sum = med(ar, "checksum"), med(br, "checksum")
        scale = max(1.0, abs(old_sum), abs(new_sum))
        same_checksum = math.isfinite(old_sum) and math.isfinite(new_sum) and abs(old_sum - new_sum) <= 1.0e-9 * scale
        workload = "OK" if same_shape and same_checksum else "CHANGED"
        changed |= workload != "OK"
        print(f"{key[1]:38s} {old:11.4f} {new:11.4f} {speed:8.3f}x {delta:8.2f}% {workload:>10s}")
    missing = sorted(set(a) ^ set(b))
    if missing:
        changed = True
        print(f"WARNING: {len(missing)} benchmark case(s) exist in only one input")
    if args.fail_on_workload_change and changed:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
