import pathlib
import sys
import tempfile

library, package_root = sys.argv[1:]
sys.path.insert(0, package_root)
import cfd_solvers

api = cfd_solvers.load(library)
assert api.version == "0.19.2"
with tempfile.TemporaryDirectory() as tmp:
    path = pathlib.Path(tmp) / "python-field.xdmf"
    api.write_xdmf_tri3(path, [(0, 0), (1, 0), (0, 1)], [(0, 1, 2)], [1, 2, 3], "temperature")
    text = path.read_text(encoding="utf-8")
    assert 'TopologyType="Triangle"' in text
    assert 'Name="temperature"' in text
print("v0.19.2 Python binding smoke passed")
