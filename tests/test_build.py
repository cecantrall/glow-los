import sys
from unittest.mock import patch

import pytest

from glow_los import build


def test_build_stops_before_compilation_without_compiler(monkeypatch):
    monkeypatch.setenv("FC", "missing-fortran")
    with patch.object(build.shutil, "which", return_value=None), patch.object(build.subprocess, "run") as run:
        with pytest.raises(RuntimeError, match="compiler.*was not found"):
            build.build()
        run.assert_not_called()


def test_build_uses_active_python_and_preserves_serious_ieee_diagnostics(monkeypatch, tmp_path):
    (tmp_path / "netcdf.mod").touch()
    (tmp_path / "libnetcdff.so").touch()
    monkeypatch.setenv("FC", "/toolchain/gfortran")
    monkeypatch.setenv("NETCDFF_INCLUDE", str(tmp_path))
    monkeypatch.setenv("NETCDFF_LIB", str(tmp_path))
    monkeypatch.delenv("FFLAGS", raising=False)
    with patch.object(build.shutil, "which", side_effect=lambda name: name), patch.object(build.subprocess, "run") as run:
        build.build()
    assert run.call_count == 3
    for call in run.call_args_list[:2]:
        command = call.args[0]
        assert command[:3] == [sys.executable, "-m", "numpy.f2py"]
        assert command[command.index("--backend") + 1] == "meson"
        assert call.kwargs["env"]["FC"] == "/toolchain/gfortran"
    redister = run.call_args_list[2].args[0]
    assert redister[0] == "/toolchain/gfortran"
    assert "-ffpe-summary=invalid,zero,overflow" in redister
    assert "-ffast-math" not in redister


def test_build_rejects_incomplete_netcdf_install_before_compiling(monkeypatch, tmp_path):
    monkeypatch.setenv("NETCDFF_INCLUDE", str(tmp_path))
    monkeypatch.setenv("NETCDFF_LIB", str(tmp_path))
    with patch.object(build.shutil, "which", side_effect=lambda name: name), patch.object(build.subprocess, "run") as run:
        with pytest.raises(RuntimeError, match="NETCDFF_INCLUDE"):
            build.build()
        run.assert_not_called()
