#!/usr/bin/env python3
"""Run the cfd_solvers performance matrix and persist comparable history.

The executable emits one CFD_PERF JSON record per kernel. This runner adds the
experimental matrix, host/build provenance, repeated-run statistics, and an
optional baseline regression gate. It deliberately does not mix results from
different build flags or machines silently.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import pathlib
import platform
import statistics
import subprocess
import sys
from typing import Any


def executable_path(build: pathlib.Path) -> pathlib.Path:
    names = ["cfd-performance-bench"]
    if os.name == "nt":
        names.append("cfd-performance-bench.exe")
    for config in ("", "Release", "RelWithDebInfo", "MinSizeRel", "Debug"):
        for name in names:
            candidate = build / config / name if config else build / name
            if candidate.exists():
                if os.name != "nt" and not os.access(candidate, os.X_OK):
                    try:
                        candidate.chmod(candidate.stat().st_mode | 0o111)
                    except OSError:
                        pass
                return candidate
    raise FileNotFoundError(f"cfd-performance-bench not found below {build}")


def parse_size(value: str) -> tuple[int, int, int]:
    parts = value.lower().split("x")
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("size must be NXxNYxNZ, for example 32x32x32")
    try:
        result = tuple(int(part) for part in parts)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("size dimensions must be integers") from exc
    if any(part < 1 for part in result):
        raise argparse.ArgumentTypeError("size dimensions must be positive")
    return result  # type: ignore[return-value]


def cache_entries(build: pathlib.Path) -> dict[str, str]:
    cache = build / "CMakeCache.txt"
    if not cache.exists():
        return {}
    result: dict[str, str] = {}
    for line in cache.read_text(errors="ignore").splitlines():
        if not line or line.startswith("//") or line.startswith("#") or "=" not in line or ":" not in line:
            continue
        lhs, value = line.split("=", 1)
        result[lhs.split(":", 1)[0]] = value.strip()
    return result


def compiler_version(compiler: str) -> str:
    if not compiler or compiler == "unknown":
        return "unknown"
    try:
        result = subprocess.run([compiler, "--version"], text=True, capture_output=True, check=False, timeout=5)
        return (result.stdout or result.stderr).splitlines()[0].strip()
    except (OSError, subprocess.SubprocessError, IndexError):
        return "unknown"


def cpu_affinity() -> list[int] | str:
    try:
        return sorted(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        return "unknown"


def git_info(repo: pathlib.Path) -> dict[str, Any]:
    def run(*args: str) -> str:
        try:
            result = subprocess.run(["git", *args], cwd=repo, text=True, capture_output=True, check=False, timeout=5)
            return result.stdout.strip()
        except (OSError, subprocess.SubprocessError):
            return "unknown"

    return {
        "commit": run("rev-parse", "HEAD"),
        "branch": run("branch", "--show-current"),
        "dirty": bool(run("status", "--porcelain")),
    }


def run_one(executable: pathlib.Path, size: tuple[int, int, int], threads: int,
            warmup: int, steps: int, quick: bool) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    nx, ny, nz = size
    command = [str(executable), "--nx", str(nx), "--ny", str(ny), "--nz", str(nz),
               "--warmup", str(warmup), "--steps", str(steps), "--threads", str(threads)]
    if quick:
        command.append("--quick")
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if result.returncode != 0:
        sys.stderr.write(result.stdout + result.stderr)
        raise RuntimeError(f"benchmark failed with exit code {result.returncode}: {' '.join(command)}")
    records: list[dict[str, Any]] = []
    setup: dict[str, Any] = {}
    for line in result.stdout.splitlines():
        if line.startswith("CFD_PERF "):
            record = json.loads(line[len("CFD_PERF "):])
            record["requested_size"] = "x".join(str(value) for value in size)
            record["requested_threads"] = threads
            records.append(record)
        elif line.startswith("CFD_PERF_SETUP "):
            setup = json.loads(line[len("CFD_PERF_SETUP "):])
    if not records or not setup:
        raise RuntimeError("benchmark did not emit the required CFD_PERF and CFD_PERF_SETUP records")
    return records, setup


def summary_key(record: dict[str, Any]) -> tuple[Any, ...]:
    return (record["benchmark"], record["implementation"], record["nx"], record["ny"], record["nz"], record["threads"])


def summarize(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = {}
    for record in records:
        grouped.setdefault(summary_key(record), []).append(record)
    result: list[dict[str, Any]] = []
    for key, group in sorted(grouped.items(), key=lambda item: str(item[0])):
        elapsed = [float(item["elapsed_ms"]) for item in group]
        throughput = [float(item["cells_per_second"]) for item in group]
        gbps = [float(item["estimated_gbps"]) for item in group]
        working_set = [int(item["working_set_bytes"]) for item in group]
        errors = [float(item["oracle_max_abs_error"]) for item in group]
        residuals = [float(item["residual_l2"]) for item in group]
        result.append({
            "benchmark": key[0],
            "implementation": key[1],
            "nx": key[2], "ny": key[3], "nz": key[4], "threads": key[5],
            "runs": len(group),
            "elapsed_ms_median": statistics.median(elapsed),
            "elapsed_ms_min": min(elapsed),
            "elapsed_ms_max": max(elapsed),
            "elapsed_ms_stdev": statistics.stdev(elapsed) if len(elapsed) > 1 else 0.0,
            "cells_per_second_median": statistics.median(throughput),
            "estimated_gbps_median": statistics.median(gbps),
            "working_set_bytes_max": max(working_set),
            "oracle_max_abs_error_max": max(errors),
            "residual_l2_max": max(residuals),
        })
    return result


def compare_baseline(current: list[dict[str, Any]], baseline: list[dict[str, Any]], threshold: float) -> list[dict[str, Any]]:
    old = {summary_key(item): item for item in baseline}
    comparisons: list[dict[str, Any]] = []
    for item in current:
        key = summary_key(item)
        previous = old.get(key)
        if previous is None:
            comparisons.append({"key": key, "status": "new"})
            continue
        before = float(previous["elapsed_ms_median"])
        after = float(item["elapsed_ms_median"])
        change = (after / before - 1.0) * 100.0 if before > 0.0 else 0.0
        comparisons.append({
            "key": key,
            "status": "regression" if change > threshold else "pass",
            "baseline_ms": before,
            "current_ms": after,
            "change_percent": change,
            "speedup": before / after if after > 0.0 else 0.0,
        })
    return comparisons


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", required=True, type=pathlib.Path)
    parser.add_argument("--output", required=True, type=pathlib.Path)
    parser.add_argument("--history", type=pathlib.Path)
    parser.add_argument("--size", action="append", type=parse_size)
    parser.add_argument("--threads", action="append", type=int)
    parser.add_argument("--warmup", type=int, default=5)
    parser.add_argument("--steps", type=int, default=30)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--baseline", type=pathlib.Path)
    parser.add_argument("--max-regression-percent", type=float, default=5.0)
    parser.add_argument("--check", action="store_true", help="return non-zero when a baseline regression exceeds the threshold")
    args = parser.parse_args()
    if args.warmup < 0 or args.steps < 1 or args.repeats < 1:
        parser.error("warmup must be non-negative; steps and repeats must be positive")
    if args.max_regression_percent < 0.0:
        parser.error("--max-regression-percent must be non-negative")

    sizes = args.size or ([(16, 16, 16), (32, 32, 32)] if args.quick else [(32, 32, 32), (64, 64, 64)])
    threads = args.threads or [1, 2, 4]
    executable = executable_path(args.build_dir)
    raw_records: list[dict[str, Any]] = []
    setup_records: list[dict[str, Any]] = []
    for size in sizes:
        for thread_count in threads:
            for repeat in range(args.repeats):
                records, setup = run_one(executable, size, thread_count, args.warmup, args.steps, args.quick)
                for record in records:
                    record["repeat"] = repeat
                    raw_records.append(record)
                setup["repeat"] = repeat
                setup["requested_size"] = "x".join(str(value) for value in size)
                setup["requested_threads"] = thread_count
                setup_records.append(setup)
                print(f"{size[0]}x{size[1]}x{size[2]} threads={thread_count} repeat={repeat + 1}/{args.repeats} records={len(records)}")

    summary = summarize(raw_records)
    baseline_comparison: list[dict[str, Any]] = []
    if args.baseline:
        baseline_payload = json.loads(args.baseline.read_text())
        baseline_comparison = compare_baseline(summary, baseline_payload.get("summary", []), args.max_regression_percent)

    cache = cache_entries(args.build_dir)
    compiler = cache.get("CMAKE_CXX_COMPILER", "unknown")
    payload: dict[str, Any] = {
        "schema": "cfd_solvers.performance_suite.v1",
        "timestamp_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "git": git_info(pathlib.Path(__file__).resolve().parents[1]),
        "host": {
            "os": platform.platform(), "machine": platform.machine(),
            "processor": platform.processor(), "logical_cpus": os.cpu_count(),
            "cpu_affinity": cpu_affinity(), "python": platform.python_version(),
            "compiler": compiler, "compiler_version": compiler_version(compiler),
        },
        "build": {
            "type": cache.get("CMAKE_BUILD_TYPE", "unknown"),
            "feature_flags": {key: cache.get(key, "unknown") for key in (
                "CFD_ENABLE_OPENMP", "CFD_ENABLE_SYCL", "CFD_ENABLE_MPI",
                "CFD_ENABLE_NATIVE_ARCH", "CFD_FAST_MATH", "CFD_ENABLE_HDF5")},
        },
        "matrix": {"sizes": ["x".join(str(value) for value in size) for size in sizes],
                   "threads": threads, "warmup": args.warmup, "steps": args.steps,
                   "repeats": args.repeats, "quick": args.quick},
        "setup_summary": setup_records,
        "summary": summary,
        "baseline_comparison": baseline_comparison,
        "raw_records": raw_records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    if args.history:
        args.history.parent.mkdir(parents=True, exist_ok=True)
        with args.history.open("a", encoding="utf-8") as stream:
            stream.write(json.dumps(payload, separators=(",", ":")) + "\n")

    regressions = [item for item in baseline_comparison if item.get("status") == "regression"]
    print(f"wrote {args.output}; {len(summary)} summaries, {len(raw_records)} raw records")
    if regressions:
        print(f"baseline regressions above {args.max_regression_percent:.2f}%: {len(regressions)}")
        for item in regressions:
            print(f"  {item['key']}: {item['change_percent']:.2f}% slower")
    return 1 if args.check and regressions else 0


if __name__ == "__main__":
    raise SystemExit(main())
