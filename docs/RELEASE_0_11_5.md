# Release 0.11.5 - multi-server supervision hardening

This release adds the next deployment/supervision layer on top of v0.11.4 supervised multi-server execution.

## Added

- Multi-server supervision plans that derive dashboard/operator artifacts from supervised remote jobs.
- Local/SSH access probes, Docker daemon probes and campaign-root writability checks.
- Retry launch wrappers with configurable attempt count and retry delay.
- Stdout/stderr log-tail command plans for each supervised remote job.
- Sensitive command-display redaction for inline `TOKEN=`, `PASSWORD=`, `SECRET=`, `KEY=`, `AWS_` and `GITHUB_` style values.
- Dashboard JSON and supervision TSV summary writers.
- CLI smoke case `particle-multiserver-supervision`.
- Focused regression target `cfd-v0115-multiserver-supervision-tests`.

## Validation

- Full CPU CTest matrix: 107/107 passed.
- Focused ASan + UBSan + leak detection for the v0.11.5 multi-server supervision slice: passed.

## Limitations

This remains an offline, deterministic supervision contract. It does not yet manage real SSH keys, store secrets, maintain persistent remote agents, stream logs through a web API, or validate live Docker daemons in the release environment.
