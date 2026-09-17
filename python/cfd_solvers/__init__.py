"""Dependency-free ctypes bindings for the cfd_solvers C ABI."""
from .bindings import SolverLibrary, load
__all__ = ["SolverLibrary", "load"]
