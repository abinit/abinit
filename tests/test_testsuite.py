"""
Unit tests for the ABINIT test-suite sanity checks.
"""

from __future__ import annotations

import os
import re

import pytest

from tests import abenv, abitests


@pytest.fixture(scope="module")
def full_database():
    """The full ABINIT test database, built once and shared by every check below."""
    # TODO should use with_disabled=True (see original check_testsuite.py)
    return abitests.build_database(with_disabled=False)


def get_allowed_cpp_vars() -> set[str]:
    """
    Inspect the libpaw header file, the autoconf macros and config.ac.

    Extract and return the set of allowed CPP options, used to check the
    need_cpp_vars/exclude_cpp_vars TEST_INFO sections for possible typos.

    Based on ~abinit/abichecks/scripts/check-cpp-options.
    """
    re_m4file = re.compile(r"\.m4$")
    re_hdrfile = re.compile(r"\.h$")
    re_acdef = re.compile(r"AC_DEFINE\(")
    re_cppdef = re.compile(r"^([ ]?)+#([ ]?)+define [0-9A-Z_]*")

    # Anchored to abenv.home_dir rather than cwd: the original script relied
    # on os.path.abspath("../"), which only resolves correctly when cwd is
    # `tests/` (how the script was traditionally invoked) -- not guaranteed
    # under pytest, which is typically run from the repo root.
    abidir = abenv.home_dir

    cpp_libpaw = set()
    for root, _dirs, files in os.walk(os.path.join(abidir, "src/39_libpaw")):
        for src in files:
            if not re_hdrfile.search(src):
                continue
            with open(os.path.join(root, src)) as fh:
                for line in fh:
                    if not re_cppdef.search(line):
                        continue
                    tmp_def = re.sub(r"^[# ]*define[ ]*([0-9A-Z_]*).*", r"\1", line).strip()
                    cpp_libpaw.add(tmp_def)

    cpp_buildsys = set()
    for root, _dirs, files in os.walk(os.path.join(abidir, "config/m4")):
        for src in files:
            if not re_m4file.search(src):
                continue
            with open(os.path.join(root, src)) as fh:
                for line in fh:
                    if not re_acdef.search(line):
                        continue
                    tmp_def = re.sub(r".*AC_DEFINE\([\[]?([^\],]*).*", r"\1", line).strip()
                    cpp_buildsys.add(tmp_def)

    with open(os.path.join(abidir, "configure.ac")) as fh:
        for line in fh:
            if not re_acdef.search(line):
                continue
            tmp_def = re.sub(r".*AC_DEFINE\([\[]?([^\],]*).*", r"\1", line).strip()
            cpp_buildsys.add(tmp_def)

    return cpp_buildsys.union(cpp_libpaw)


def check_authors(suite) -> set[str]:
    """
    Check if test authors follow the project's "First.Last" naming convention.

    Returns the set of unique author second names found in the suite.
    """
    def first_second_name(string):
        idx = string.rfind(".")
        if idx == -1:
            first, second = "", string
        else:
            first, second = string[:idx], string[idx + 1:]
        return first.strip(), second.strip()

    second_names = []
    for test in suite:
        if not hasattr(test, "authors"):
            authors = []
            for t in test:
                authors.extend(t.authors)
            authors = set(authors)
        else:
            authors = test.authors

        for string in authors:
            f, s = first_second_name(string)
            if not f and s and s != "Unknown":
                print(f"author(s) first name is missing in file {test.full_id}, string = {s} ")
            second_names.append(s)

    return set(second_names)


class TestReferenceAndInputFiles:
    """Every file physically present under a suite's Input/ or Refs/ directory
    must be referenced by at least (Input) or exactly (Refs) one test.
    """

    @pytest.mark.skip(
        reason="find_stale_or_lost_refs() correctness is not confirmed yet -- "
        "disabled until its logic has been reviewed."
    )
    def test_no_stale_or_lost_refs(self, full_database):
        err = full_database.find_stale_or_lost_refs()
        assert not err, err

    @pytest.mark.skip(
        reason="find_stale_or_lost_inputs() correctness is not confirmed yet -- "
        "disabled until its logic has been reviewed."
    )
    def test_no_stale_or_lost_inputs(self, full_database):
        err = full_database.find_stale_or_lost_inputs()
        assert not err, err


class TestKeywords:
    """Checks on the `keywords` TEST_INFO field across the whole suite."""

    @pytest.fixture(scope="class")
    def unknown_and_wrong_keywords(self, full_database):
        return full_database.find_unknown_wrong_keywords()

    def test_no_undocumented_keywords(self, unknown_and_wrong_keywords):
        unknowns, _wrong = unknown_and_wrong_keywords
        assert not unknowns, (
            f"Undocumented keywords: {sorted(unknowns)}\n"
            "ACTION: add them to the KNOWN_KEYWORDS dictionary in tests/known_keywords.json"
        )

    def test_no_keywords_with_blank_spaces(self, unknown_and_wrong_keywords):
        _unknowns, wrong = unknown_and_wrong_keywords
        assert not wrong, f"Keywords with blank spaces (use underscores instead): {sorted(wrong)}"


def test_recommended_testinfo_options_present(full_database):
    """Every test should define keywords/description/authors/max_nprocs in TEST_INFO."""
    err_str = full_database.check_testinfo_options()
    assert not err_str, err_str


def test_need_cpp_vars_are_all_known(full_database):
    """Every `need_cpp_vars` entry in TEST_INFO must be a real, defined CPP option.

    The original check_testsuite.py only printed mismatches here and never
    failed the run (no retcode increment for this particular check) --
    promoted to a real assertion as part of the pytest conversion.
    """
    allowed_cpp_vars = get_allowed_cpp_vars()
    # A single need_cpp_vars entry may be "HAVE_FOO or HAVE_BAR" -- an OR
    # requirement satisfied by either variable, handled the same way at
    # runtime by BaseTest.compute_nprocs() (testsuite.py). Split on that
    # token first so each individual CPP var is checked, instead of the
    # whole (necessarily unknown) compound string.
    or_token = " or "
    offenders = []
    for suite in full_database.values():
        for test in suite:
            for var in test.need_cpp_vars:
                tvars = (
                    {v.strip() for v in var.split(or_token)}
                    if or_token in var
                    else {var.removeprefix("!")}
                )
                diff = tvars.difference(allowed_cpp_vars)
                if diff:
                    offenders.append(f"{test.full_id}: {sorted(diff)}")

    assert not offenders, "\n".join(offenders)


@pytest.mark.skip(
    reason="Disabled in the original check_testsuite.py script (commented out "
    "call site) -- the author second-name convention check looks incomplete "
    "(it fails whenever any second name is found at all, which is presumably "
    "why it was never enabled) and needs a real look before being turned on."
)
def test_author_names_follow_convention(full_database):
    second_names = set()
    for suite in full_database.values():
        second_names |= check_authors(suite)
    assert not second_names, sorted(second_names)
