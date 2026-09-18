#!/usr/bin/env sh
set -eu

repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root"
git config core.hooksPath .githooks
echo "Configured Git hooks path: .githooks"
echo "The pre-commit hook will run clang-format on staged C/C++ files."
