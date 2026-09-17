import ctypes
import os
from pathlib import Path

class SolverLibrary:
    def __init__(self, path):
        self._lib = ctypes.CDLL(os.fspath(path))
        self._lib.cfd_solvers_version.restype = ctypes.c_char_p
        self._lib.cfd_write_xdmf_tri3.restype = ctypes.c_int
        self._lib.cfd_write_xdmf_tri3.argtypes = [ctypes.c_char_p,ctypes.POINTER(ctypes.c_double),ctypes.c_size_t,ctypes.POINTER(ctypes.c_size_t),ctypes.c_size_t,ctypes.POINTER(ctypes.c_double),ctypes.c_char_p,ctypes.c_char_p,ctypes.c_size_t]
    @property
    def version(self): return self._lib.cfd_solvers_version().decode("utf-8")
    def write_xdmf_tri3(self,path,points,triangles,scalar=None,name="field"):
        xy=[float(v) for point in points for v in point]
        cells=[int(v) for cell in triangles for v in cell]
        if len(xy)%2 or len(cells)%3: raise ValueError("points must be XY pairs and triangles must be index triples")
        xy_array=(ctypes.c_double*len(xy))(*xy);cell_array=(ctypes.c_size_t*len(cells))(*cells)
        scalar_array=None
        if scalar is not None:
            values=[float(v) for v in scalar]
            if len(values)!=len(xy)//2: raise ValueError("scalar length must equal point count")
            scalar_array=(ctypes.c_double*len(values))(*values)
        error=ctypes.create_string_buffer(1024)
        rc=self._lib.cfd_write_xdmf_tri3(os.fsencode(path),xy_array,len(xy)//2,cell_array,len(cells)//3,scalar_array,name.encode(),error,len(error))
        if rc: raise RuntimeError(error.value.decode("utf-8",errors="replace"))

def load(path=None):
    path=path or os.environ.get("CFD_SOLVERS_LIBRARY")
    if not path: raise RuntimeError("set CFD_SOLVERS_LIBRARY or pass the C API shared-library path")
    return SolverLibrary(Path(path))
