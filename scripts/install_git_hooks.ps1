Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$repoRoot = git rev-parse --show-toplevel
if (-not $repoRoot) {
    throw "Not inside a Git repository."
}

Push-Location $repoRoot
try {
    git config core.hooksPath .githooks
    Write-Output "Configured Git hooks path: .githooks"
    Write-Output "The pre-commit hook will run clang-format on staged C/C++ files."
}
finally {
    Pop-Location
}
