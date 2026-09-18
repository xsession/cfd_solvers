#!/usr/bin/env python3
"""Run practical cfd_solvers examples and persist comparable timing records."""
from __future__ import annotations
import argparse, datetime as dt, json, os, pathlib, platform, re, subprocess, sys

TARGETS = [
    "cfd-example-phase01-hpc", "cfd-example-phase02-lbm", "cfd-example-phase03-fvm",
    "cfd-example-phase04-fem", "cfd-example-phase05-fdtd", "cfd-example-phase06-optics",
    "cfd-example-phase07-chemistry-corrosion", "cfd-example-phase08-multiphysics",
    "cfd-example-phase10-rf", "cfd-example-phase11-spice", "cfd-example-phase12-em-pic",
    "cfd-example-phase13-dem", "cfd-example-phase14-acoustics", "cfd-example-phase15-tcad",
    "cfd-example-phase16a-battery",
]

def cpu_model() -> str:
    if sys.platform.startswith("linux"):
        try:
            for line in pathlib.Path("/proc/cpuinfo").read_text(errors="ignore").splitlines():
                if line.lower().startswith("model name"):
                    return line.split(":", 1)[1].strip()
        except OSError:
            pass
    return platform.processor() or "unknown"

def cache_entries(build: pathlib.Path) -> dict[str, str]:
    cache = build / "CMakeCache.txt"
    if not cache.exists():
        return {}
    out: dict[str, str] = {}
    for line in cache.read_text(errors="ignore").splitlines():
        if not line or line.startswith("//") or line.startswith("#") or "=" not in line or ":" not in line:
            continue
        lhs, value = line.split("=", 1)
        name = lhs.split(":", 1)[0]
        out[name] = value.strip()
    return out

def compiler_version(compiler: str) -> str:
    if not compiler or compiler == "unknown":
        return "unknown"
    try:
        p = subprocess.run([compiler, "--version"], text=True, capture_output=True, check=False, timeout=5)
        line = (p.stdout or p.stderr).splitlines()
        return line[0].strip() if line else "unknown"
    except (OSError, subprocess.SubprocessError):
        return "unknown"

def cpu_affinity() -> list[int] | str:
    try:
        return sorted(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        return "unknown"

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--build-dir", required=True, type=pathlib.Path)
    ap.add_argument("--output", type=pathlib.Path)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--scale", type=int, default=1)
    ap.add_argument("--repeat", type=int, default=1)
    args = ap.parse_args()
    results=[]
    for repeat in range(args.repeat):
        for target in TARGETS:
            exe=args.build_dir/target
            cmd=[str(exe)]
            if args.quick: cmd.append("--quick")
            if args.scale != 1: cmd += ["--scale", str(args.scale)]
            proc=subprocess.run(cmd, text=True, capture_output=True, check=False)
            if proc.returncode != 0:
                sys.stderr.write(proc.stdout+proc.stderr)
                return proc.returncode
            lines=[line[len("CFD_BENCH "):] for line in proc.stdout.splitlines() if line.startswith("CFD_BENCH ")]
            if len(lines)!=1:
                raise RuntimeError(f"{target}: expected one CFD_BENCH record, got {len(lines)}")
            rec=json.loads(lines[0]); rec["repeat"]=repeat; rec["target"]=target; results.append(rec)
            print(f"{target:42s} {rec['simulation_ms']:12.4f} ms  {rec['throughput']:12.3g} {rec['work_unit_name']}/s")
    cache = cache_entries(args.build_dir)
    compiler = cache.get("CMAKE_CXX_COMPILER", "unknown")
    feature_flags = {k: cache.get(k, "unknown") for k in (
        "CFD_ENABLE_OPENMP", "CFD_ENABLE_SYCL", "CFD_ENABLE_MPI",
        "CFD_ENABLE_NATIVE_ARCH", "CFD_FAST_MATH", "CFD_ENABLE_HDF5")
    }
    payload={
        "schema":2,
        "timestamp_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
        "host":{
            "os":platform.platform(), "machine":platform.machine(), "cpu":cpu_model(),
            "logical_cpus":os.cpu_count(), "cpu_affinity":cpu_affinity(),
            "python":platform.python_version(), "compiler":compiler,
            "compiler_version":compiler_version(compiler),
        },
        "build":{"type":cache.get("CMAKE_BUILD_TYPE", "unknown"), "feature_flags":feature_flags},
        "environment":{"OMP_NUM_THREADS":os.environ.get("OMP_NUM_THREADS"), "CFD_BENCH_THREADS":os.environ.get("CFD_BENCH_THREADS")},
        "mode":{"quick":args.quick,"scale":args.scale,"repeat":args.repeat},
        "results":results,
    }
    out=args.output or pathlib.Path(f"benchmark_{dt.datetime.now().strftime('%Y%m%d_%H%M%S')}.json")
    out.parent.mkdir(parents=True,exist_ok=True); out.write_text(json.dumps(payload,indent=2)+"\n")
    print(f"wrote {out}")
    return 0
if __name__=="__main__": raise SystemExit(main())
