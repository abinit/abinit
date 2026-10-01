"""Tests for the read-only ROBODoc-to-Doxygen input filter."""

from __future__ import annotations

import hashlib
import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

ABINIT_ROOT = Path(__file__).resolve().parents[2]
FILTER = Path(__file__).with_name("robodoc_filter.py")
SPEC = importlib.util.spec_from_file_location("robodoc_filter", FILTER)
assert SPEC is not None
assert SPEC.loader is not None
ROBODOC_FILTER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ROBODOC_FILTER)
transform_text = ROBODOC_FILTER.transform_text


def test_transform_parameters_and_preserve_lines() -> None:
    """ROBODoc sections become Doxygen commands without shifting lines."""
    source = """!!****f* module/example
!! FUNCTION
!!  Demonstrate the conversion.
!! INPUTS
!!  [count]=Number of items
!!  workspace=Internal workspace, not a dummy argument
!! OUTPUT
!!  values(:)=Computed values
!! SIDE EFFECTS
!!  state<state_t>=State updated by the routine
!! SOURCE

subroutine example(count, values, state)
!!***
"""

    transformed = transform_text(source)

    assert transformed.count("\n") == source.count("\n")
    assert "!! @brief" in transformed
    assert "!! @param[in] count Number of items" in transformed
    assert "!!  workspace=Internal workspace, not a dummy argument" in transformed
    assert "!! @param[out] values Computed values" in transformed
    assert "!! @param[inout] state State updated by the routine" in transformed
    assert "!!\n!!\nsubroutine example" in transformed
    assert transformed.startswith("!>\n")


@pytest.mark.parametrize(
    "relative_path",
    [
        "shared/libpaw/src/m_pawrad.F90",
        "src/41_geometry/m_crystal.F90",
    ],
)
def test_filter_does_not_modify_source(relative_path: str) -> None:
    """Executing the filter leaves representative source files byte-identical."""
    source_path = ABINIT_ROOT / relative_path
    digest_before = hashlib.sha256(source_path.read_bytes()).digest()

    result = subprocess.run(
        [sys.executable, FILTER, source_path],
        check=True,
        stdout=subprocess.PIPE,
    )

    assert result.stdout
    assert hashlib.sha256(source_path.read_bytes()).digest() == digest_before
    assert result.stdout.count(b"\n") == source_path.read_bytes().count(b"\n")
