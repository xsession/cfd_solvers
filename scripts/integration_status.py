#!/usr/bin/env python3
"""Report completion of the Markdown integration tracker and pinned upstreams."""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
phase_re = re.compile(r"^## Phase\s+(\d+[A-Za-z]?)\s*-\s*(.+)$")
box_re = re.compile(r"^\s*- \[([ xX])\]\s+(.+)$")


def parse_tracker(path: Path) -> dict[str, dict[str, object]]:
    phases: dict[str, dict[str, object]] = {}
    phase: str | None = None
    for line in path.read_text(encoding="utf-8").splitlines():
        m = phase_re.match(line)
        if m:
            phase = m.group(1)
            phases.setdefault(phase, {"title": m.group(2).strip(), "done": 0, "total": 0})
            continue
        m = box_re.match(line)
        if m and phase is not None:
            phases[phase]["total"] = int(phases[phase]["total"]) + 1
            if m.group(1).lower() == "x":
                phases[phase]["done"] = int(phases[phase]["done"]) + 1
    return phases


def pct(done: int, total: int) -> float:
    return 100.0 * done / total if total else 0.0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tracker", type=Path, default=ROOT / "docs" / "INTEGRATION_TRACKER.md")
    parser.add_argument("--manifest", type=Path, default=ROOT / "docs" / "upstreams.json")
    parser.add_argument("--json", action="store_true", dest="as_json")
    args = parser.parse_args()

    phases = parse_tracker(args.tracker)
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    overall_done = sum(int(v["done"]) for v in phases.values())
    overall_total = sum(int(v["total"]) for v in phases.values())

    report = {
        "overall": {"done": overall_done, "total": overall_total, "percent": pct(overall_done, overall_total)},
        "phases": [],
        "upstreams": [],
    }
    def phase_sort_key(value: str) -> tuple[int, str]:
        match = re.match(r"(\d+)([A-Za-z]?)$", value)
        return (int(match.group(1)), match.group(2)) if match else (10**9, value)

    for number in sorted(phases, key=phase_sort_key):
        item = phases[number]
        done, total = int(item["done"]), int(item["total"])
        report["phases"].append({
            "phase": number, "title": item["title"], "done": done, "total": total, "percent": pct(done, total)
        })
    for upstream in manifest["upstreams"]:
        raw_number = upstream.get("tracker_phase")
        number = str(raw_number)
        item = phases.get(number, {"done": 0, "total": 0})
        done, total = int(item["done"]), int(item["total"])
        pin = upstream.get("commit") or upstream.get("release") or "unversioned"
        report["upstreams"].append({
            "name": upstream["name"], "phase": number, "done": done, "total": total,
            "percent": pct(done, total), "pin": pin,
        })

    if args.as_json:
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0

    print("Integration tracker status")
    print(f"overall: {overall_done}/{overall_total} ({pct(overall_done, overall_total):.1f}%)")
    print("\nBy phase:")
    for item in report["phases"]:
        print(f"  Phase {item['phase']}: {item['done']:3d}/{item['total']:3d} ({item['percent']:5.1f}%)  {item['title']}")
    print("\nBy upstream family:")
    for item in report["upstreams"]:
        print(f"  {item['name']:<22} Phase {item['phase']}: {item['done']:3d}/{item['total']:3d} "
              f"({item['percent']:5.1f}%) pin={str(item['pin'])[:12]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
