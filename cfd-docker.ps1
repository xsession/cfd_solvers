# cfd_solvers Docker deployment helper (PowerShell).
#
# All logic lives here; the companion cfd-docker.bat is a thin launcher that
# invokes this script. Usage:
#   .\cfd-docker.bat build
#   .\cfd-docker.bat run lbm-d3q19-cpu
#   .\cfd-docker.bat run --threads 8 battery-dfn
#   .\cfd-docker.bat test
#   .\cfd-docker.bat mpi --ranks 4
#   .\cfd-docker.bat list
#   .\cfd-docker.bat ps          # list running cfd containers
#   .\cfd-docker.bat clean       # stop + remove containers
#
# Set $Env:CFD_TAG to override the image tag (default cfd_solvers:dev).
# Set $Env:CFD_MPI = "1" to build/run the MPI-enabled image.

param(
    [Parameter(Position = 0)][string]$Action = "help",
    [Parameter(Position = 1, ValueFromRemainingArguments = $true)][string[]]$Args
)

$ErrorActionPreference = "Stop"
$tag   = if ($Env:CFD_TAG) { $Env:CFD_TAG } else { "cfd_solvers:dev" }
$mpi   = if ($Env:CFD_MPI) { "ON" } else { "OFF" }
$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path

# Parse flags out of $Args.
$threads = 4
$ranks   = 4
$rest = @()
for ($i = 0; $i -lt $Args.Length; $i++) {
    switch ($Args[$i]) {
        "--threads" { $threads = [int]$Args[++$i]; continue }
        "--ranks"   { $ranks   = [int]$Args[++$i]; continue }
        default     { $rest += $Args[$i] }
    }
}

function Invoke-Docker([string[]]$DockerArgs) {
    & docker @DockerArgs
    if ($LASTEXITCODE -ne 0) { throw "docker $($DockerArgs -join ' ') failed with exit $LASTEXITCODE" }
}

function Build-Image([string]$MpiFlag) {
    Write-Host "Building $tag (MPI=$MpiFlag) ..."
    # args are key=value pairs so no quoting surprises.
    $dockerArgs = @("build","-f","Dockerfile","-t","$tag",
        "--build-arg","CMAKE_BUILD_TYPE=Release",
        "--build-arg","CFD_ENABLE_OPENMP=ON",
        "--build-arg","CFD_ENABLE_MPI=$MpiFlag",
        "--build-arg","CFD_ENABLE_NATIVE_ARCH=ON",
        "--build-arg","CFD_ENABLE_HDF5=ON",
        ".")
    Invoke-Docker $dockerArgs
    Write-Host "Built $tag"
}

switch ($Action) {
    "build" { Build-Image $mpi }

    "run" {
        $case = if ($rest.Count -ge 1) { $rest[0] } else { "--list" }
        Write-Host "Running cfd-solve $case (threads=$threads) ..."
        Invoke-Docker @("run","--rm","-e","OMP_NUM_THREADS=$threads","--entrypoint","/cfd_solvers/build/cfd-solve","$tag",$case)
    }

    "bench" {
        $rest = if ($rest.Count -ge 1) { $rest } else { @("--lattice","d3q19","--nx","64","--ny","64","--nz","64","--warmup","10","--steps","100","--streaming","both") }
        Write-Host "Running cfd-bench ..."
        Invoke-Docker @("run","--rm","-e","OMP_NUM_THREADS=$threads","--entrypoint","/cfd_solvers/build/cfd-bench","$tag") + $rest
    }

    "test" {
        Write-Host "Running full CTest suite ..."
        Invoke-Docker @("run","--rm","-e","OMP_NUM_THREADS=$threads","--entrypoint","bash","$tag","-c","ctest --test-dir /cfd_solvers/build --output-on-failure")
    }

    "mpi" {
        Build-Image "ON"
        Write-Host "Running distributed d3q19 (ranks=$ranks) ..."
        Invoke-Docker @("run","--rm","-e","OMP_NUM_THREADS=$threads","--entrypoint","mpirun","$tag",
            "--allow-run-as-root","--oversubscribe","-np","$ranks",
            "/cfd_solvers/build/cfd-distributed","--backend","cpu","--lattice","d3q19","--nx","192","--ny","128","--nz","96","--steps","100")
    }

    "list" {
        Invoke-Docker @("run","--rm","--entrypoint","/cfd_solvers/build/cfd-solve","$tag","--list")
    }

    "ps" { Invoke-Docker @("ps","--filter","ancestor=$tag") }
    "clean" { Invoke-Docker @("rm","-f", (& docker ps -aq --filter "ancestor=$tag") | Where-Object { $_ }) }
    "help" {
        Write-Host @"
cfd_solvers Docker helper

  .\cfd-docker.bat build                 Build the runtime image ($tag)
  .\cfd-docker.bat run [case]            Run a cfd-solve case (default --list)
  .\cfd-docker.bat run --threads 8 <case>
  .\cfd-docker.bat bench [args]          Run cfd-bench (defaults to d3q19 64^3)
  .\cfd-docker.bat test                   Run the full CTest suite in a container
  .\cfd-docker.bat mpi [--ranks 4]       Build+run the MPI distributed d3q19 case
  .\cfd-docker.bat list                   List available cfd-solve cases
  .\cfd-docker.bat ps                     List running cfd containers
  .\cfd-docker.bat clean                  Remove running cfd containers

Environment:
  CFD_TAG    image tag          (default cfd_solvers:dev)
  CFD_MPI=1  build MPI image    (default off)
"@
    }
    default { Write-Error "Unknown action '$Action'. Use 'help'." }
}
