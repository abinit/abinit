"""Unit tests for the parametrized tests (tests/plans/parametrized_tests_plan.md, section 11.1)."""
from __future__ import annotations

import ast
import os
import pickle
import textwrap
from pathlib import Path

import pytest

pytest.importorskip("pydantic")
pytest.importorskip("yaml")

from tests import abenv  # noqa: E402

from ..testsuite import AbinitTestInfoParser  # noqa: E402
from .expand import (  # noqa: E402
    ExpansionError,
    build_instances,
    combinations,
    draw_nproc,
    expand_template,
    find_conflicts,
    insertion_lines,
    scan_datasets,
)
from .models import REACHABLE_NPROCS  # noqa: E402
from .paramtests import (  # noqa: E402
    BrokenTemplateTest, ParamFileToTest, ParametrizedAbinitTest, ParametrizedChain, ParametrizedTest)
from .template import TemplateError, load_template  # noqa: E402
from .yamltext import compare_yaml_text, extract_yaml_lines, fldiff_text, group_lines, merge_overrides  # noqa: E402

TESTS_DIR = Path(abenv.tests_dir)
T01 = TESTS_DIR / "parametrized" / "Input" / "t01.tmpl.abi"
T07 = TESTS_DIR / "parametrized" / "Input" / "t07.tmpl.abi"

BODY = """\
acell 3*10.18
ecut 8.0
pp_dirpath "$ABI_PSPDIR"
pseudos "PseudosTM_pwteter/14si.pspnc"
"""


def write_template(tmp_path: Path, yaml_block: str, body: str = BODY, name: str = "t07.tmpl.abi") -> Path:
    inp = tmp_path / "Input"
    inp.mkdir(exist_ok=True)
    path = inp / name
    path.write_text("%YAML\n" + textwrap.dedent(yaml_block) + "%ENDYAML\n" + body)
    return path


MINIMAL = """\
test_info:
  setup: {executable: abinit}
  files:
    files_to_test:
      - {name: "{id}.abo", mode: yaml_docs, tolabs: 1.1e-8, tolrel: 1.0e-8}
  extra_info: {keywords: [NC], authors: [J. Doe], description: "test"}
axes:
  - {type: grid, parameters: {useylm: [0, 1]}}
"""


# -------------------------------------------------------------------------------------------------
# Models and template loading
# -------------------------------------------------------------------------------------------------

class TestTemplateValidation:

    def test_t01_is_valid(self):
        tmpl = load_template(str(T01))
        assert tmpl.number == "01"
        assert tmpl.spec.test_info.files_to_test[0].mode == "yaml_docs"

    def test_minimal(self, tmp_path):
        tmpl = load_template(str(write_template(tmp_path, MINIMAL)))
        assert len(combinations(tmpl.spec)) == 2

    @pytest.mark.parametrize("old, new, match", [
        ("mode: yaml_docs,", "mode: yaml,", "obsolete"),
        ("executable: abinit}", "executable: abinit, psp_files: [a.psp8]}", "not allowed"),
        ("executable: abinit}", "executable: abinit, test_chain: [t11.abi]}", "must be the file name of a template"),
        ("executable: abinit}", "executable: abinit, test_chain: [t11.tmpl.abi, t11.tmpl.abi]}", "repeated"),
        ("executable: abinit}", "executable: abinit, max_nprocs: 2}", "belongs to section"),
        ("executable: abinit}", "executable: abinit, foo: 1}", "unknown TEST_INFO option"),
        ("tolrel: 1.0e-8}", "tolrel: 1.0e-8, overrides: {'Bad.path.x': ignore}}", "invalid override path"),
        ('{name: "{id}.abo", mode: yaml_docs, tolabs: 1.1e-8, tolrel: 1.0e-8}',
         '{name: "{id}.abo", tolabs: 1.1e-8, overrides: {DatasetInfo: ignore}}', "only with mode: yaml_docs"),
        ('{name: "{id}.abo", ', "{", "name"),                                        # the name is mandatory
        ('{name: "{id}.abo", ', '{name: "out.abo", ', "must contain {id}"),
        ('{name: "{id}.abo", ', '{name: "{stem}.abo", ', "only placeholder"),
        ('      - {name: "{id}.abo", mode: yaml_docs, tolabs: 1.1e-8, tolrel: 1.0e-8}\n',
         '      - {name: "{id}.abo", mode: yaml_docs, tolabs: 1.1e-8, tolrel: 1.0e-8}\n'
         '      - {name: "{id}.abo", tolabs: 1.0}\n', "same file more than once"),
        ("keywords: [NC]", "keywords: [NOT_A_KEYWORD]", "not registered"),
        ("parameters: {useylm: [0, 1]}", "parameters: {pseudos: [a, b]}", "cannot be swept"),
        ("parameters: {useylm: [0, 1]}", "parameters: {nproc: [1, 3]}", "reachable by testbot"),
        ("  - {type: grid, parameters: {useylm: [0, 1]}}\n",
         "  - {type: grid, parameters: {useylm: [0, 1]}}\n  - {type: grid, parameters: {useylm: [1]}}\n",
         "appears in axes"),
    ])
    def test_invalid(self, tmp_path, old, new, match):
        assert old in MINIMAL
        path = write_template(tmp_path, MINIMAL.replace(old, new))
        with pytest.raises(TemplateError, match=match):
            load_template(str(path), known_keywords={"NC": ""})

    def test_yaml_test_section_rejected(self, tmp_path):
        path = write_template(tmp_path, MINIMAL.replace("  extra_info:", "  yaml_test: {file: x.yaml}\n  extra_info:"))
        with pytest.raises(TemplateError, match="yaml_test"):
            load_template(str(path))

    def test_duplicate_keys_rejected(self, tmp_path):
        path = write_template(tmp_path, MINIMAL.replace("tolabs: 1.1e-8,", "tolabs: 1.1e-8, tolabs: 2.0,"))
        with pytest.raises(TemplateError, match="duplicate keys"):
            load_template(str(path))

    def test_exclude_unknown_key(self, tmp_path):
        path = write_template(tmp_path, MINIMAL + "exclude:\n  - {fftalg: 312}\n")
        with pytest.raises(TemplateError, match="exclude uses parameters"):
            load_template(str(path))

    def test_ini_block_in_body_rejected(self, tmp_path):
        path = write_template(tmp_path, MINIMAL, body=BODY + "#%%<BEGIN TEST_INFO>\n")
        with pytest.raises(TemplateError, match="must not contain"):
            load_template(str(path))


# -------------------------------------------------------------------------------------------------
# Expansion
# -------------------------------------------------------------------------------------------------

class TestExpansion:

    def test_product_zip_explicit_exclude(self, tmp_path):
        block = MINIMAL.replace("  - {type: grid, parameters: {useylm: [0, 1]}}\n", textwrap.dedent("""\
              - {type: grid, mode: zip, parameters: {useylm: [0, 1], fftalg: [112, 401]}}
              - type: explicit
                values:
                  - {nproc: 1, npband: 1}
                  - {nproc: 2, npband: 2}
            """)) + "exclude:\n  - {useylm: 1, nproc: 2}\n"
        tmpl = load_template(str(write_template(tmp_path, block)))
        combos = combinations(tmpl.spec)
        assert combos == [
            {"useylm": 0, "fftalg": 112, "nproc": 1, "npband": 1},
            {"useylm": 0, "fftalg": 112, "nproc": 2, "npband": 2},
            {"useylm": 1, "fftalg": 401, "nproc": 1, "npband": 1},
        ]
        insts = build_instances(tmpl)
        assert [i.test_id for i in insts] == ["ti0000_07", "ti0001_07_MPI2", "ti0002_07"]
        # nproc is never inserted; the other swept variables are.
        assert "nproc" not in insts[1].input_text.split("#%%")[0].split("inserted")[-1]
        assert "npband 2" in insts[1].input_text

    def test_parallel_test_info(self, tmp_path):
        block = MINIMAL.replace("parameters: {useylm: [0, 1]}", "parameters: {nproc: [1, 4]}")
        tmpl = load_template(str(write_template(tmp_path, block)))
        inst = build_instances(tmpl)[1]
        assert "[NCPU_4]" in inst.input_text and "nprocs_to_test = 4" in inst.input_text
        path = tmp_path / "inst.abi"
        path.write_text(inst.input_text)
        info = AbinitTestInfoParser(str(path)).generate_testinfo_nprocs(4)
        assert info.nprocs_to_test == [4] and info.max_nprocs == 4
        assert info.files_to_test[0].name == "ti0001_07_MPI4.abo"
        assert "PARAMETRIZED" in info.keywords

    def test_parallel_instance_without_reference(self, tmp_path):
        """--update-ref expands with no_check: parallel instances still need their [NCPU_N] section."""
        block = MINIMAL.replace("parameters: {useylm: [0, 1]}", "parameters: {nproc: [1, 2]}")
        tests = expand_template(str(write_template(tmp_path, block)), abenv, no_check=True)
        assert [type(t).__name__ for t in tests] == ["ParametrizedAbinitTest"] * 2
        assert tests[1].id == "ti0001_07_MPI2" and tests[1].nprocs_to_test == [2] and tests[1].files_to_test == []

    def test_max_nprocs_draw_is_deterministic_and_reachable(self):
        values = {draw_nproc("t05.tmpl.abi", 8) for _ in range(5)}
        assert len(values) == 1 and values.pop() in REACHABLE_NPROCS

    @pytest.mark.parametrize("body_line", ["useylm 1", "useylm3 0", "useylm: 0", "useylm+ 1", "USEYLM 1"])
    def test_conflicts(self, body_line):
        assert find_conflicts(BODY + body_line + "\n", ["useylm"]) == ["useylm"]

    def test_no_conflict_in_comments_or_prefixes(self):
        assert find_conflicts(BODY + "# useylm 1\nuseylmx 1\n", ["useylm"]) == []

    def test_datasets(self):
        assert scan_datasets(BODY) == (1, None)
        assert scan_datasets("ndtset 3\njdtset 1 2 5\n") == (3, [1, 2, 5])
        with pytest.raises(ExpansionError, match="udtset"):
            scan_datasets("udtset 2 2\n")

    def test_insertion_lines(self):
        params = {"useylm": 1, "nproc": 2, "ngkpt": [2, 2, 2]}
        assert insertion_lines(params, None, 1, None) == ["useylm 1", "ngkpt 2 2 2"]
        assert insertion_lines(params, 2, 2, None) == ["useylm2 1", "ngkpt2 2 2 2"]
        assert insertion_lines(params, "all", 2, None) == ["useylm 1", "ngkpt 2 2 2"]
        with pytest.raises(ExpansionError, match="must set `dataset`"):
            insertion_lines(params, None, 2, None)
        with pytest.raises(ExpansionError, match="not one of the datasets"):
            insertion_lines(params, 3, 3, [1, 2, 5])


# -------------------------------------------------------------------------------------------------
# Text-based YAML comparison (yamltext)
# -------------------------------------------------------------------------------------------------

REF_ABO = """\
 some text 1.0
--- !DatasetInfo
dimensions: {natom: 2, mpw: 302, }
...
--- !ResultsGS
etotal    :  -8.86517716E+00
convergence: {deltae: -1.4E-13, residm:  9.284E-09, }
cartesian_stress_tensor: # hartree/bohr^3
- [ -2.75462403E-05,   0.00000000E+00,   0.00000000E+00, ]
...
"""


def compare(tmp_path, out_text, rules=None, base=None):
    ref, out = tmp_path / "ref.abo", tmp_path / "out.abo"
    ref.write_text(REF_ABO)
    out.write_text(out_text)
    base = base or {"tolnlines": 0, "tolabs": 1.1e-8, "tolrel": 1.0e-8, "fld_options": "-medium"}
    return compare_yaml_text(str(ref), str(out), base, merge_overrides(rules), str(tmp_path / "work"))


class TestYamlText:

    def test_extraction_and_groups(self, tmp_path):
        (tmp_path / "ref.abo").write_text(REF_ABO)
        lines = extract_yaml_lines(str(tmp_path / "ref.abo"))
        assert {(yl.tag, yl.key) for yl in lines} >= {("ResultsGS", "etotal"), ("ResultsGS", "cartesian_stress_tensor")}
        groups = group_lines(lines, merge_overrides({"ResultsGS.convergence": {"fld_options": "-easy"}}))
        assert not any("mpw" in l for g in groups.values() for l in g)        # DatasetInfo ignored by default
        assert all("convergence" in l for l in groups["ResultsGS.convergence"])
        assert all(l.startswith(" ") for l in fldiff_text(groups[None]))      # first-column pitfall

    def test_identical(self, tmp_path):
        assert compare(tmp_path, REF_ABO).status == "succeeded"

    def test_run_description_is_ignored(self, tmp_path):
        assert compare(tmp_path, REF_ABO.replace("mpw: 302", "mpw: 151")).isok

    def test_real_change_fails(self, tmp_path):
        assert not compare(tmp_path, REF_ABO.replace("-8.86517716E+00", "-8.86517816E+00")).isok

    def test_array_row_change_fails(self, tmp_path):
        # Array rows start with "-": without the blank prefix fldiff would ignore them.
        assert not compare(tmp_path, REF_ABO.replace("[ -2.75462403E-05", "[ -2.85462403E-05")).isok

    def test_override_loosens_only_its_group(self, tmp_path):
        out = REF_ABO.replace("residm:  9.284E-09", "residm:  1.092E-12")
        rules = {"ResultsGS.convergence": {"fld_options": "-easy"}}
        assert compare(tmp_path, out, rules).isok
        out_bad = out.replace("-8.86517716E+00", "-8.86517816E+00")
        assert not compare(tmp_path, out_bad, rules).isok

    def test_ignore_rule(self, tmp_path):
        out = REF_ABO.replace("-8.86517716E+00", "-8.86517816E+00")
        assert compare(tmp_path, out, {"ResultsGS.etotal": "ignore"}).isok

    def test_empty_extract_is_an_error(self, tmp_path):
        res = compare(tmp_path, "no yaml at all\n")
        assert not res.isok


# -------------------------------------------------------------------------------------------------
# Test objects built by the hook
# -------------------------------------------------------------------------------------------------

class TestExpandTemplate:

    def test_t01_objects(self):
        tests = expand_template(str(T01), abenv)
        assert [t.id for t in tests] == ["ti0000_01", "ti0001_01", "ti0002_01", "ti0003_01"]
        for t in tests:
            assert isinstance(t, ParametrizedAbinitTest)
            assert t.suite_name == "parametrized"
            assert t.inp_fname == str(T01)                       # a real file (docs build, stale-input check)
            assert str(T01) in t.inputs_used
            assert "PARAMETRIZED" in t.keywords
            f = t.files_to_test[0]
            assert isinstance(f, ParamFileToTest) and f.param_mode == "yaml_docs" and f.ref_name == "t01.abo"
        restored = pickle.loads(pickle.dumps(tests))
        assert [t.id for t in restored] == [t.id for t in tests]

    def test_every_input_file_is_used(self):
        """Mimic find_stale_or_lost_inputs for the parametrized suite."""
        import importlib.util
        suite_dir = TESTS_DIR / "parametrized"
        spec = importlib.util.spec_from_file_location("param_suite", suite_dir / "__init__.py")
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        used = set()
        for name in mod.inp_files:
            for t in expand_template(str(suite_dir / "Input" / name.lstrip("-")), abenv):
                used.update(os.path.abspath(p) for p in t.inputs_used)
        present = {str(p) for p in (suite_dir / "Input").iterdir() if not p.name.startswith(".")}
        assert present <= used

    def test_broken_template_becomes_a_failing_test(self, tmp_path):
        path = write_template(tmp_path, MINIMAL, body=BODY + "useylm1 0\n")
        tests = expand_template(str(path), abenv)
        assert len(tests) == 1 and isinstance(tests[0], BrokenTemplateTest)
        assert "swept variables also appear in the body" in tests[0].error_message


# -------------------------------------------------------------------------------------------------
# Chains of prerequisites (plan 8.2)
# -------------------------------------------------------------------------------------------------

PREREQ = """\
test_info:
  setup: {executable: mrgddb}
  extra_info: {keywords: [NC], authors: [J. Doe], description: "prerequisite"}
"""


CHAIN = """\
test_info:
  setup: {executable: abinit, test_chain: [t04.tmpl.abi, t05.tmpl.abi]}
  files:
    files_to_test:
      - {name: "{id}.abo", tolabs: 1.0e-3}
  extra_info: {keywords: [NC], authors: [J. Doe], description: "chain"}
axes:
  - type: explicit
    values:
      - {nproc: 1, eph_np_pqbks: [1, 1, 1, 1, 1]}
      - {nproc: 2, eph_np_pqbks: [2, 1, 1, 1, 1]}
      - {nproc: 4, eph_np_pqbks: [2, 2, 1, 1, 1]}
      - {nproc: 4, eph_np_pqbks: [1, 2, 2, 1, 1]}
"""


def write_chain(tmp_path: Path) -> Path:
    """A swept template t07 with two prerequisites: t04 (mrgddb, stdin input) and t05 (abinit)."""
    write_template(tmp_path, PREREQ, body="MgO_DDB\nMgO phonons\n0\n", name="t04.tmpl.abi")
    write_template(tmp_path, PREREQ.replace("executable: mrgddb}", "executable: abinit, output_prefix: nscfo}"),
                   name="t05.tmpl.abi")
    return write_template(tmp_path, CHAIN)


class TestChains:

    def test_one_chain_per_nproc(self, tmp_path):
        chains = expand_template(str(write_chain(tmp_path)), abenv)
        assert all(isinstance(c, ParametrizedChain) for c in chains)
        assert [(c.id, c.nproc) for c in chains] == [("tc_07", 1), ("tc_07_MPI2", 2), ("tc_07_MPI4", 4)]
        c4 = chains[2]
        assert [t.id for t in c4.prerequisites] == ["tp1_07_MPI4", "tp2_07_MPI4"]
        assert [t.id for t in c4.instances] == ["ti0002_07_MPI4", "ti0003_07_MPI4"]
        assert [t.executable for t in c4.prerequisites] == ["mrgddb", "abinit"]
        assert [os.path.basename(t.inp_fname) for t in c4.prerequisites] == ["t04.tmpl.abi", "t05.tmpl.abi"]
        for t in c4.prerequisites:
            assert t.instance_index is None and t.no_check and not t.files_to_test
            assert t.nprocs_to_test == [4] and t.max_nprocs == 1
        for t in c4.instances:
            assert t.nprocs_to_test == [4] and t.max_nprocs == 4
        restored = pickle.loads(pickle.dumps(chains))
        assert [c.id for c in restored] == [c.id for c in chains]

    def test_prerequisites_run_on_one_process_in_dynamic_mode(self, tmp_path):
        """compute_nprocs: dynamic runs a prerequisite of the 4-process chain on 1 process, static skips it (P12)."""
        c4 = expand_template(str(write_chain(tmp_path)), abenv)[2]
        prereq, inst = c4.prerequisites[0], c4.instances[0]

        class Env:  # compute_nprocs only reads the CPP variables of the build
            defined_cppvars = set()
        assert prereq.compute_nprocs(Env(), 4, "dynamic") == (1, "")
        assert prereq.compute_nprocs(Env(), 4, "static")[0] == 0
        assert inst.compute_nprocs(Env(), 4, "dynamic") == (4, "")
        assert prereq.compute_nprocs(Env(), 2, "dynamic")[0] == 0   # selected only by the passes with 4 processes

    def test_stdin_input_starts_with_the_body(self, tmp_path):
        """mrgddb reads the name of the output file from the first line: the header goes at the end."""
        prereq = expand_template(str(write_chain(tmp_path)), abenv)[0].prerequisites[0]
        assert isinstance(prereq, ParametrizedTest)
        lines = prereq.input_text.splitlines()
        assert lines[0] == "MgO_DDB"
        assert any(l.startswith("# parametrized: prerequisite 1 (t04.tmpl.abi)") for l in lines)

    def test_missing_prerequisite(self, tmp_path):
        path = write_template(tmp_path, MINIMAL.replace("executable: abinit}", "executable: abinit, test_chain: [t01.tmpl.abi]}"))
        tests = expand_template(str(path), abenv)
        assert isinstance(tests[0], BrokenTemplateTest) and "not found" in tests[0].error_message

    @pytest.mark.parametrize("extra, match", [
        ("axes:\n  - {type: grid, parameters: {useylm: [0, 1]}}\n", "cannot have axes"),
        ("", "cannot have its own setup.test_chain"),
    ])
    def test_invalid_prerequisite(self, tmp_path, extra, match):
        prereq = PREREQ + extra
        if not extra:
            prereq = prereq.replace("executable: mrgddb}", "executable: mrgddb, test_chain: [t02.tmpl.abi]}")
        write_template(tmp_path, prereq, body="out_DDB\n", name="t01.tmpl.abi")
        path = write_template(tmp_path, MINIMAL.replace("executable: abinit}", "executable: abinit, test_chain: [t01.tmpl.abi]}"))
        tests = expand_template(str(path), abenv)
        assert isinstance(tests[0], BrokenTemplateTest) and match in tests[0].error_message

    def test_prerequisite_no_check_by_default(self, tmp_path):
        prereq = write_template(tmp_path, PREREQ, body="out_DDB\n", name="t01.tmpl.abi")
        assert load_template(str(prereq), prerequisite=True).spec.test_info.no_check
        with pytest.raises(TemplateError, match="at least one file"):
            load_template(str(prereq))   # the same file is not a valid swept template


# -------------------------------------------------------------------------------------------------
# Fixtures (plan 8.4)
# -------------------------------------------------------------------------------------------------

FIXTURE = """\
test_info:
  setup: {executable: mrgddb}
fixture: {outputs: [out_DDB]}
"""


def _build_counting(args):
    """Worker of test_cache_is_process_safe: every process asks for the same entry."""
    root, counter = args
    from . import fixture_cache as fc

    def build(tmp):
        with open(counter, "a") as fh:     # one line per actual build
            fh.write(f"{os.getpid()}\n")
        import time
        time.sleep(0.3)                    # keep the lock long enough for the others to wait on it
        Path(tmp, "out").write_text("data")
        return "ok", ""

    entry, st = fc.ensure(root, "fx", "k" * 64, build)
    return st["status"], Path(entry, "out").read_text()


class TestFixtures:

    def test_t07_fixtures(self):
        tests = expand_template(str(T07), abenv)
        assert [f.name for f in tests[0].fixtures] == ["mgo_ddb_444q", "mgo_dvdb_444q", "mgo_wfk_444k_70b"]
        ddb = tests[0].fixtures[0]
        assert ddb.outputs == ["MgO_444q_DDB"] and ddb.test.executable == "mrgddb"
        assert ddb.test.input_text.splitlines()[0] == "MgO_444q_DDB"      # stdin input: header at the end
        wfk = tests[0].fixtures[2].test
        assert wfk.executable == "abinit" and wfk.no_check and wfk.max_nprocs == 4   # fixture.nproc
        for t in tests:
            assert not isinstance(t, ParametrizedChain)
            used = [os.path.basename(p) for p in t.inputs_used]
            assert {"t07.tmpl.abi", "mgo_ddb_444q.fixture.abi", "mgo_dvdb_444q.fixture.abi", "mgo_wfk_444k_70b.fixture.abi"} <= set(used)
        restored = pickle.loads(pickle.dumps(tests))
        assert [f.name for f in restored[2].fixtures] == ["mgo_ddb_444q", "mgo_dvdb_444q", "mgo_wfk_444k_70b"]

    def test_templates_share_fixtures(self, tmp_path):
        """Two templates with the same DDB fixture and different WFK fixtures: one DDB entry, two WFK entries."""
        write_template(tmp_path, FIXTURE, body="out_DDB\nddb\n0\n", name="ddb.fixture.abi")
        wfk = PREREQ.replace("executable: mrgddb}", "executable: abinit, output_prefix: wfko}")
        write_template(tmp_path, wfk + "fixture: {outputs: [wfko_WFK]}\n", name="wfk_nb30.fixture.abi")
        write_template(tmp_path, wfk + "fixture: {outputs: [wfko_WFK]}\n", body=BODY + "nband 40\n",
                       name="wfk_nb40.fixture.abi")
        t1 = expand_template(str(write_template(tmp_path, MINIMAL + "fixtures: [ddb, wfk_nb30]\n", name="t01.tmpl.abi")), abenv)
        t2 = expand_template(str(write_template(tmp_path, MINIMAL + "fixtures: [ddb, wfk_nb40]\n", name="t02.tmpl.abi")), abenv)

        class Env:  # only path_of_bin is used
            def path_of_bin(self, name):
                return "/nonexistent/" + name
        (ddb1, wfk1), (ddb2, wfk2) = t1[0].fixtures, t2[0].fixtures
        assert ddb1.cache_key(Env()) == ddb2.cache_key(Env())          # shared cache entry
        assert wfk1.cache_key(Env()) != wfk2.cache_key(Env())          # different inputs, different entries

    def test_effective_nproc(self):
        from types import SimpleNamespace
        from .paramtests import Fixture
        fx = Fixture("wfk", "wfk.fixture.abi", ["WFK"], [], test=None, nproc=4)
        mpi, seq = SimpleNamespace(has_mpirun=True), SimpleNamespace(has_mpirun=False)
        assert fx.effective_nproc(mpi, 8) == 4
        assert fx.effective_nproc(mpi, 2) == 2            # CPU budget of the run
        assert fx.effective_nproc(seq, 8) == 1            # no MPI launcher
        assert fx.effective_nproc(None, None) == 1
        fx.nproc_strict = True
        assert fx.effective_nproc(mpi, 8) == 4
        assert fx.effective_nproc(mpi, 2) is None         # strict: never fewer processes
        assert Fixture("ddb", "p", ["DDB"], [], test=None).effective_nproc(None, None) == 1

    def test_nproc_of_t07_wfk_fixture(self):
        wfk = expand_template(str(T07), abenv)[0].fixtures[2]
        assert wfk.name == "mgo_wfk_444k_70b" and wfk.nproc == 4 and not wfk.nproc_strict
        assert wfk.test.max_nprocs == 4                   # the static compute_nprocs lets it run on 4

    def test_nproc_is_part_of_the_key_and_strict_skip_is_not_cached(self, tmp_path):
        from types import SimpleNamespace

        class Env:  # only path_of_bin is used
            def path_of_bin(self, name):
                return "/nonexistent/" + name
        wfk = expand_template(str(T07), abenv)[0].fixtures[2]
        key4 = wfk.cache_key(Env())
        wfk.nproc = 2
        assert wfk.cache_key(Env()) != key4
        wfk.nproc, wfk.nproc_strict = 4, True
        entry, st = wfk.ensure(Env(), str(tmp_path / "fixtures"), runner=SimpleNamespace(has_mpirun=False))
        assert st["status"] == "skipped" and "nproc_strict" in st["message"]
        assert not os.path.exists(entry)                  # another run, with MPI, can still build it

    def test_prebuild_selects_the_tests_that_run(self, monkeypatch):
        from . import prebuild
        tests = expand_template(str(T07), abenv)
        calls = []
        monkeypatch.setattr(prebuild, "prebuild_fixtures",
                            lambda sel, env, runner, root, max_cpus: calls.append(([t.id for t in sel], root)))

        class Env:
            defined_cppvars = set()
        prebuild.prebuild_for_suite(tests, Env(), None, "/tmp/ws", mpi_nprocs=4, runmode="dynamic", max_cpus=8)
        ids, root = calls[0]
        assert ids == ["ti0000_07", "ti0004_07_MPI4", "ti0005_07_MPI4", "ti0006_07_MPI4"] and root == "/tmp/ws/fixtures"
        assert [f.name for f in prebuild.fixtures_of(tests)] == ["mgo_ddb_444q", "mgo_dvdb_444q", "mgo_wfk_444k_70b"]

    def test_needs_cycle_and_outputs(self, tmp_path):
        write_template(tmp_path, FIXTURE.replace("outputs: [out_DDB]}", "outputs: [a], needs: [fb]}"), name="fa.fixture.abi")
        write_template(tmp_path, FIXTURE.replace("outputs: [out_DDB]}", "outputs: [b], needs: [fa]}"), name="fb.fixture.abi")
        path = write_template(tmp_path, MINIMAL + "fixtures: [fa]\n")
        tests = expand_template(str(path), abenv)
        assert isinstance(tests[0], BrokenTemplateTest) and "cycle in fixture.needs: fa -> fb -> fa" in tests[0].error_message

        write_template(tmp_path, FIXTURE, name="fa.fixture.abi")
        write_template(tmp_path, FIXTURE, name="fb.fixture.abi")
        path = write_template(tmp_path, MINIMAL + "fixtures: [fa, fb]\n")
        tests = expand_template(str(path), abenv)
        assert isinstance(tests[0], BrokenTemplateTest) and "both declare the output out_DDB" in tests[0].error_message

    def test_needs_are_resolved(self, tmp_path):
        write_template(tmp_path, FIXTURE.replace("outputs: [out_DDB]}", "outputs: [a]}"), name="fa.fixture.abi")
        write_template(tmp_path, FIXTURE.replace("outputs: [out_DDB]}", "outputs: [b], needs: [fa]}"), name="fb.fixture.abi")
        path = write_template(tmp_path, MINIMAL + "fixtures: [fb]\n")
        test = expand_template(str(path), abenv)[0]
        fb = test.fixtures[0]
        assert fb.name == "fb" and [n.name for n in fb.needs] == ["fa"]
        assert [f.name for f in fb.all_fixtures()] == ["fa", "fb"]
        assert any(p.endswith("fa.fixture.abi") for p in test.inputs_used)

    @pytest.mark.parametrize("yaml_block, match", [
        (FIXTURE.replace("fixture: {outputs: [out_DDB]}", "fixture: {outputs: []}"), "at least one file"),
        (FIXTURE.replace("[out_DDB]", "[dir/out_DDB]"), "plain file name"),
        (FIXTURE.replace("{executable: mrgddb}", "{executable: mrgddb, no_check: no}"), "not compared"),
        (FIXTURE.replace("{executable: mrgddb}", "{executable: mrgddb}\n  paral_info: {max_nprocs: 2}"), "fixture.nproc"),
        (FIXTURE.replace("{outputs: [out_DDB]}", "{outputs: [out_DDB], nproc: 4}"), "runs on 1 process only"),
        (FIXTURE.replace("{outputs: [out_DDB]}", "{outputs: [out_DDB], nproc: 0}"), "greater than 0"),
        (FIXTURE + "axes:\n  - {type: grid, parameters: {useylm: [0, 1]}}\n", "axes"),
        (FIXTURE.replace("{outputs: [out_DDB]}", "{outputs: [out_DDB], needs: [bad-name]}"), "invalid fixture name"),
        (FIXTURE.replace("{outputs: [out_DDB]}", "{outputs: [out_DDB], needs: [Mgo_ddb]}"), "lowercase"),
    ])
    def test_invalid_fixture(self, tmp_path, yaml_block, match):
        from .template import load_fixture
        path = write_template(tmp_path, yaml_block, name="fa.fixture.abi")
        with pytest.raises(TemplateError, match=match):
            load_fixture(str(path))

    def test_names_are_lowercase(self, tmp_path):
        from .template import load_fixture
        path = write_template(tmp_path, FIXTURE, name="MgO_ddb.fixture.abi")
        with pytest.raises(TemplateError, match="lowercase"):
            load_fixture(str(path))
        path = write_template(tmp_path, MINIMAL + "fixtures: [MgO_ddb]\n")
        with pytest.raises(TemplateError, match="lowercase"):
            load_template(str(path))

    def test_input_names_are_unique_ignoring_case(self):
        """Two files differing only by case cannot coexist on a case-insensitive filesystem."""
        names = [p.name for p in (TESTS_DIR / "parametrized" / "Input").iterdir()]
        assert len({n.lower() for n in names}) == len(names)

    def test_fixtures_and_test_chain_exclusive(self, tmp_path):
        text = MINIMAL.replace("executable: abinit}", "executable: abinit, test_chain: [t01.tmpl.abi]}") + "fixtures: [fa]\n"
        path = write_template(tmp_path, text)
        with pytest.raises(TemplateError, match="cannot be used together"):
            load_template(str(path))

    def test_missing_fixture(self, tmp_path):
        path = write_template(tmp_path, MINIMAL + "fixtures: [nope]\n")
        tests = expand_template(str(path), abenv)
        assert isinstance(tests[0], BrokenTemplateTest) and "nope.fixture.abi not found" in tests[0].error_message

    def test_cache_is_process_safe(self, tmp_path):
        """8 processes ask for the same entry at the same time: it is built once, and all of them get it."""
        import multiprocessing
        root, counter = str(tmp_path / "fixtures"), str(tmp_path / "builds.txt")
        with multiprocessing.get_context("spawn").Pool(8) as pool:
            results = pool.map(_build_counting, [(root, counter)] * 8)
        assert results == [("ok", "data")] * 8
        assert len(Path(counter).read_text().split()) == 1
        entries = sorted(p.name for p in Path(root).iterdir())
        assert entries == ["fx-kkkkkkkkkkkkkkkk", "fx-kkkkkkkkkkkkkkkk.lock"]   # no temporary directory left
        assert not os.access(Path(root, entries[0], "out"), os.W_OK)          # outputs are read-only

    def test_cache_records_failures(self, tmp_path):
        from . import fixture_cache as fc
        calls = []

        def build(tmp):
            calls.append(1)
            raise RuntimeError("boom")
        for _ in range(2):
            _, st = fc.ensure(str(tmp_path), "fx", "a" * 64, build)
            assert st["status"] == "failed" and "boom" in st["message"]
        assert len(calls) == 1                                             # not retried with the same key
        fc.remove_all(str(tmp_path))
        assert not tmp_path.exists()


# -------------------------------------------------------------------------------------------------
# Number of MPI processes: cross-check with the output, and runtests' results round trip
# -------------------------------------------------------------------------------------------------

class TestMpiNprocs:

    def test_parse_output(self, tmp_path):
        from .paramtests import mpi_nprocs_from_output
        path = tmp_path / "x.abo"
        path.write_text("header\n-   mpi_nproc: 4, omp_nthreads: 1 (-1 if OMP is not activated)\n")
        assert mpi_nprocs_from_output(str(path)) == 4
        assert mpi_nprocs_from_output(str(tmp_path / "missing.abo")) is None

    def test_mismatch_is_an_error(self, tmp_path):
        from .execute import _outcome_from_test
        test = [t for t in expand_template(str(T07), abenv) if t.id == "ti0001_07_MPI2"][0]
        (tmp_path / "ti0001_07_MPI2.abo").write_text("-   mpi_nproc: 1, omp_nthreads: 1\n")
        test.exceptions, test.exec_error = [], False
        out = _outcome_from_test(test, str(tmp_path), compared=False)
        assert out.status == "error" and out.nproc_run == 1 and "expected 2" in out.message
        (tmp_path / "ti0001_07_MPI2.abo").write_text("-   mpi_nproc: 2, omp_nthreads: 1\n")
        out = _outcome_from_test(test, str(tmp_path), compared=False)
        assert out.status == "ok" and out.nproc_run == 2

    def test_results_round_trip_keeps_mpi_nprocs(self):
        """runtests -j N: the parent copy of a test gets mpi_nprocs from the worker (it printed np=1 before)."""
        worker, parent = expand_template(str(T07), abenv)[2], expand_template(str(T07), abenv)[2]
        worker.mpi_nprocs, worker._executed, worker.exec_error = 4, True, False
        worker.stdout_fname, worker.workdir, worker.tot_etime, worker.run_etime = "x.stdout", "/tmp/x", 1.0, 1.0
        worker._status, worker._isok = "succeeded", True
        parent.load_results(worker.get_results())
        assert parent.mpi_nprocs == 4 and parent.full_id.endswith("[np=4]")


# -------------------------------------------------------------------------------------------------
# run_params.py -p/--patch
# -------------------------------------------------------------------------------------------------

class TestPatchOption:
    """run_params.py -p has the statuses, syntax and file selection of runtests.py -p."""

    def test_statuses(self):
        from tests.run_params import parse_patch_statuses
        assert parse_patch_statuses("failed") == {"failed"}
        assert parse_patch_statuses("failed+passed") == {"failed", "passed"}
        assert parse_patch_statuses("all") == {"failed", "passed", "succeeded"}
        assert parse_patch_statuses("not_succeeded") == {"failed", "passed", "skipped", "disabled"}
        with pytest.raises(SystemExit):
            parse_patch_statuses("differs")          # a report status, not a runtests.py status

    def test_pairs(self, tmp_path):
        from tests.run_params import patch_pairs, test_status
        from .execute import FileOutcome, InstanceOutcome
        tests = expand_template(str(T01), abenv)
        outcomes = []
        for t, status, fld in zip(tests, ["ok", "ok", "differs", "error"], ["succeeded", "passed", "failed", "failed"]):
            wd = tmp_path / t.id
            wd.mkdir()
            (wd / f"{t.id}.abo").write_text("x")
            outcomes.append(InstanceOutcome(t.id, t.instance_index, t.params, 1, str(wd), status,
                                            files=[FileOutcome(f"{t.id}.abo", "yaml_docs", fld, "")]))
        assert [test_status(o) for o in outcomes] == ["succeeded", "passed", "failed", "failed"]
        ref = str(T01.parent.parent / "Refs" / "t01.abo")
        names = lambda pairs: [(os.path.basename(o), r) for o, r in pairs]
        assert names(patch_pairs(tests, outcomes, {"failed"})) == [("ti0002_01.abo", ref), ("ti0003_01.abo", ref)]
        assert names(patch_pairs(tests, outcomes, {"passed"})) == [("ti0001_01.abo", ref)]
        assert len(patch_pairs(tests, outcomes, {"failed", "passed", "succeeded"})) == 4   # -p all: every file


# -------------------------------------------------------------------------------------------------
# run_params.py fixtures
# -------------------------------------------------------------------------------------------------

class TestFixturesCommand:

    def test_catalog_of_the_suite(self):
        from .docs import fixture_catalog
        catalog = {e["name"]: e for e in fixture_catalog(str(TESTS_DIR / "parametrized" / "Input"))}
        assert set(catalog) == {"mgo_ddb_444q", "mgo_dvdb_444q", "mgo_wfk_444k_70b"}
        wfk = catalog["mgo_wfk_444k_70b"]
        assert wfk["error"] == "" and wfk["executable"] == "abinit" and wfk["nproc"] == 4
        assert wfk["outputs"] == ["MgO_444k_70b_WFK"] and wfk["used_by"] == ["t07"]
        assert wfk["submodule"] == "MgO_eph_zpr" and "70 bands" in wfk["description"]

    def test_needs_invalid_and_unused(self, tmp_path):
        import io
        from .docs import describe_fixtures, fixture_catalog
        write_template(tmp_path, FIXTURE.replace("outputs: [out_DDB]}", "outputs: [a]}"), name="fa.fixture.abi")
        write_template(tmp_path, FIXTURE.replace("outputs: [out_DDB]}", "outputs: [b], needs: [fa]}"), name="fb.fixture.abi")
        write_template(tmp_path, FIXTURE.replace("outputs: [out_DDB]}", "outputs: [c]}"), name="fc.fixture.abi")
        write_template(tmp_path, FIXTURE.replace("outputs: [out_DDB]}", "outputs: []}"), name="bad.fixture.abi")
        write_template(tmp_path, MINIMAL + "fixtures: [fb]\n", name="t01.tmpl.abi")
        catalog = {e["name"]: e for e in fixture_catalog(str(tmp_path / "Input"))}
        assert catalog["fa"]["used_by"] == ["t01"]          # through fb.needs
        assert catalog["fb"]["used_by"] == ["t01"] and catalog["fc"]["used_by"] == []
        assert "at least one file" in catalog["bad"]["error"]
        fh = io.StringIO()
        assert describe_fixtures(str(tmp_path / "Input"), fh) == 1   # one invalid fixture
        text = fh.getvalue()
        assert "INVALID" in text and "(no template)" in text and "needs:      fa" in text

    def test_format_documentation(self, tmp_path):
        import io
        from .docs import FIXTURE_EXAMPLE, doc_fixture_format
        from .template import load_fixture
        fh = io.StringIO()
        doc_fixture_format(fh)
        text = fh.getvalue()
        for field in ("outputs:", "needs:", "nproc: positive integer", "nproc_strict:", "test_info:"):
            assert f"  {field}" in text, field
        path = tmp_path / "si_wfk.fixture.abi"
        path.write_text(FIXTURE_EXAMPLE)
        assert load_fixture(str(path)).spec.fixture.nproc == 4


# -------------------------------------------------------------------------------------------------
# tests/parametrized/README.md
# -------------------------------------------------------------------------------------------------

class TestReadme:
    """The examples of the README are valid, and its table lists the templates of the suite."""

    README = TESTS_DIR / "parametrized" / "README.md"

    def _yaml_blocks(self):
        import re as _re
        return [b for b in _re.findall(r"```text\n(.*?)```", self.README.read_text(), _re.S) if b.startswith("%YAML")]

    def test_examples_are_valid(self, tmp_path):
        from .template import load_fixture
        blocks = self._yaml_blocks()
        assert len(blocks) == 3
        (tmp_path / "Input").mkdir()
        nprocs = []
        for i, block in enumerate(blocks):
            if "\nfixture:" in block:
                path = tmp_path / "Input" / "si_wfk.fixture.abi"
                path.write_text(block)
                assert load_fixture(str(path)).spec.fixture.nproc == 4
            else:
                path = tmp_path / "Input" / f"t0{i + 1}.tmpl.abi"
                path.write_text(block)
                nprocs.append([inst.nproc for inst in build_instances(load_template(str(path)))])
        assert nprocs == [[1, 1], [1, 2, 4, 4]]     # the useylm example, and the MPI example

    def test_examples_table_lists_the_templates(self):
        import importlib.util
        spec = importlib.util.spec_from_file_location("param_suite", TESTS_DIR / "parametrized" / "__init__.py")
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        text = self.README.read_text()
        for name in mod.inp_files:
            stem = name.lstrip("-").split(".")[0]
            assert f"| `{stem}` |" in text, f"{stem} is missing from the Examples table of the README"


# -------------------------------------------------------------------------------------------------
# Documentation generated from the models (run_params.py doc)
# -------------------------------------------------------------------------------------------------

class TestDocs:

    def test_every_model_and_field_is_documented(self):
        import inspect as _inspect
        from pydantic import BaseModel
        from . import models
        from .docs import DOCUMENTED_MODELS
        documented = {m for _, m in DOCUMENTED_MODELS}
        defined = {obj for _, obj in _inspect.getmembers(models, _inspect.isclass)
                   if issubclass(obj, BaseModel) and obj is not BaseModel and obj.__module__ == models.__name__}
        assert defined == documented
        for model in documented:
            assert model.__doc__, model.__name__
            for name, field in model.model_fields.items():
                assert field.description, f"{model.__name__}.{name} has no description"

    def test_output_covers_fields_and_test_info_options(self):
        import io
        from ..testsuite import TESTCNF_KEYWORDS
        from .docs import DOCUMENTED_MODELS, doc_template_format
        from .models import FORBIDDEN_KEYS, TEMPLATE_KEY_NOTES
        fh = io.StringIO()
        doc_template_format(fh)
        text = fh.getvalue()
        for _, model in DOCUMENTED_MODELS:
            for name, field in model.model_fields.items():
                if not field.exclude:
                    assert f"  {name}: " in text, name
        for key in TESTCNF_KEYWORDS:
            if TESTCNF_KEYWORDS[key][2] != "yaml_test":
                assert f"    {key}" in text, key
        for key, why in FORBIDDEN_KEYS.items():
            assert f"{key}: NOT ALLOWED ({why})" in text
        assert set(TEMPLATE_KEY_NOTES) <= set(TESTCNF_KEYWORDS)
        assert "yaml_test:  NOT SUPPORTED" in text

    def test_example_is_a_valid_template(self, tmp_path):
        from .docs import EXAMPLE
        path = tmp_path / "Input" / "t01.tmpl.abi"
        path.parent.mkdir()
        path.write_text(EXAMPLE + BODY)
        assert len(build_instances(load_template(str(path)))) == 4


# -------------------------------------------------------------------------------------------------
# Import guard (plan section 1)
# -------------------------------------------------------------------------------------------------

def _imports(path: Path) -> list[tuple[str, bool, str]]:
    """Return (module, is_module_level, enclosing function) for the imports of `path`."""
    tree = ast.parse(path.read_text())
    result = []

    def visit(node, func, top):
        for child in ast.iter_child_nodes(node):
            if isinstance(child, (ast.Import, ast.ImportFrom)):
                names = [a.name for a in child.names] if isinstance(child, ast.Import) else [
                    ("." * child.level) + (child.module or "")]
                for n in names:
                    result.append((n, top, func))
            new_func = child.name if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)) else func
            new_top = top and not isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))
            visit(child, new_func, new_top)

    visit(tree, None, True)
    return result


class TestImportGuard:

    def test_existing_code_does_not_import_parametrized_at_module_level(self):
        files = [TESTS_DIR / "__init__.py", TESTS_DIR / "runtests.py", TESTS_DIR / "testbot.py"]
        files += [p for p in (TESTS_DIR / "pymods").glob("*.py")]
        for path in files:
            for mod, top, func in _imports(path):
                if "parametrized" not in mod:
                    continue
                assert not top, f"{path.name}: module-level import of {mod}"
                # The two hooks: template expansion, and the fixtures built before the tests (plan 8.4).
                assert (path.name, func) in {("testsuite.py", "make_abitests_from_inputs"), ("testsuite.py", "run_tests")}, \
                    f"{path.name}:{func}: only the hooks in make_abitests_from_inputs and run_tests may import {mod}"

    @pytest.mark.parametrize("name", ["paramtests.py", "yamltext.py", "execute.py", "report.py", "fixture_cache.py",
                                      "prebuild.py", "htmlreport.py"])
    def test_pickled_modules_do_not_need_optional_packages(self, name):
        """Modules reached by unpickling test objects must not import pydantic or yaml (plan section 9)."""
        banned = ("pydantic", "yaml", ".models", ".template", ".expand")
        for mod, _top, _func in _imports(TESTS_DIR / "pymods" / "parametrized" / name):
            assert not any(mod == b or mod.startswith(b + ".") for b in banned), f"{name} imports {mod}"


# -------------------------------------------------------------------------------------------------
# compare and reports (no ABINIT build needed: outputs are copies of the reference)
# -------------------------------------------------------------------------------------------------

class TestCompareAndReport:

    @staticmethod
    def make_outputs(tmp_path):
        """Workdirs whose outputs are copies of Refs/t01.abo: every instance must be ok."""
        import shutil
        tests = expand_template(str(T01), abenv)
        ref = TESTS_DIR / "parametrized" / "Refs" / "t01.abo"
        for t in tests:
            (tmp_path / t.id).mkdir()
            shutil.copy(ref, tmp_path / t.id / f"{t.id}.abo")
        return tests

    def test_compare_ok_differs_error(self, tmp_path):
        import re
        from .execute import compare_instances, is_failure
        tests = self.make_outputs(tmp_path)
        outcomes = compare_instances(tests, str(tmp_path))
        assert [o.status for o in outcomes] == ["ok"] * 4 and not is_failure(outcomes)

        out = tmp_path / "ti0002_01" / "ti0002_01.abo"
        out.write_text(re.sub(r"^(etotal\s*:\s*)(-?\d\.\d+)", lambda m: m.group(1) + f"{float(m.group(2)) + 1e-6:.8E}",
                              out.read_text(), count=1, flags=re.M))
        (tmp_path / "ti0003_01" / "ti0003_01.abo").unlink()
        outcomes = compare_instances(expand_template(str(T01), abenv), str(tmp_path))
        assert [o.status for o in outcomes] == ["ok", "ok", "differs", "error"] and is_failure(outcomes)
        assert "output not found" in outcomes[3].message
        assert os.path.exists(outcomes[2].files[0].fldiff_path)

    def test_compare_skipped_and_broken(self, tmp_path):
        from .execute import compare_instances
        tests = self.make_outputs(tmp_path)
        broken = expand_template(str(write_template(tmp_path, MINIMAL, body=BODY + "useylm 1\n")), abenv)
        outcomes = compare_instances(tests[:1] + broken, str(tmp_path), max_cpus=0)
        assert [o.status for o in outcomes] == ["skipped", "error"]

    @staticmethod
    def assert_balanced(text):
        """The tags of an HTML page are balanced (void elements excepted)."""
        import html.parser

        class Balance(html.parser.HTMLParser):
            def __init__(self):
                super().__init__()
                self.stack, self.bad = [], 0

            def handle_starttag(self, tag, attrs):
                if tag not in ("meta", "input", "br", "link"):
                    self.stack.append(tag)

            def handle_endtag(self, tag):
                if self.stack and self.stack[-1] == tag:
                    self.stack.pop()
                else:
                    self.bad += 1

        p = Balance()
        p.feed(text)
        assert p.bad == 0 and p.stack == []

    def test_reports(self, tmp_path):
        """comparison_report.json, and the pages of runtests.py written by run_params.py (template and instances)."""
        import json
        from .execute import compare_instances, write_html_pages
        from .report import build_report, write_reports
        tests = self.make_outputs(tmp_path)
        (tmp_path / "ti0001_01" / "ti0001_01.abo").unlink()
        outcomes = compare_instances(tests, str(tmp_path))
        report = build_report(str(T01), None, str(tmp_path), "compare", outcomes)
        json_path = write_reports(report, str(tmp_path / "rep" / "comparison_report.json"))

        data = json.loads(Path(json_path).read_text())
        assert data["summary"]["ok"] == 3 and data["summary"]["error"] == 1
        assert {"test_id", "params", "status", "files", "input_file", "workdir"} <= set(data["runs"][0])

        page = write_html_pages(str(T01), tests, outcomes, str(tmp_path), command="run_params.py compare t01")
        assert page == str(tmp_path / "index.html")
        text = Path(page).read_text()
        assert text.count('<tr data-status="') == 4 and text.count('<tr data-status="failed"') == 1
        assert text.count('<tr data-status="succeeded"') == 3
        assert "http" not in text                                          # self-contained
        assert 'href="ti0000_01/index.html"' in text and "reference" in text  # links to the instance pages
        assert "<th>useylm</th>" in text and "Differences" in text
        self.assert_balanced(text)

        # The page of an instance: the generic page of runtests.py plus the card of the parametrized test.
        inst = (tmp_path / "ti0002_01" / "index.html").read_text()
        assert "Instance 0002 of t01" in inst and 'href="../index.html"' in inst
        assert 'href="ti0002_01.abi"' in inst                               # input link: no .stdin file
        assert "yaml_docs group" in inst and "yaml_docs -medium" in inst
        self.assert_balanced(inst)

    def test_yaml_docs_diff_shows_the_compared_lines(self, tmp_path):
        """The text diff of a yaml_docs file is the diff of the extract (yamltext.write_extracts)."""
        import re
        from .execute import compare_instances, write_html_pages
        tests = self.make_outputs(tmp_path)
        out = tmp_path / "ti0002_01" / "ti0002_01.abo"
        out.write_text(re.sub(r"^(etotal\s*:\s*)(-?\d\.\d+)", lambda m: m.group(1) + f"{float(m.group(2)) + 1e-6:.8E}",
                              out.read_text(), count=1, flags=re.M))
        outcomes = compare_instances(tests, str(tmp_path))
        write_html_pages(str(T01), tests, outcomes, str(tmp_path))
        extract = tmp_path / "ti0002_01" / "ti0002_01.abo.yaml_docs"
        assert (extract / "output.txt").read_text().startswith("# yaml_docs[default]")
        diff = (tmp_path / "ti0002_01" / "ti0002_01.abo.diff").read_text()
        assert "output.txt" in diff and "reference.txt" in diff and "etotal" in diff
        assert "ti0002_01.abo.diff" in (tmp_path / "ti0002_01" / "index.html").read_text()

    def test_suite_pages(self, tmp_path):
        """runtests.py: one page per template in <suite workdir>/<suite>_<tNN>, linked from the index."""
        from .execute import compare_instances
        from .htmlreport import write_suite_pages
        tests = self.make_outputs(tmp_path)
        for t in tests:
            t.workdir = str(tmp_path / t.id)
        compare_instances(tests, str(tmp_path))
        paths, section = write_suite_pages(tests, str(tmp_path))
        assert paths == [str(tmp_path / "parametrized_t01" / "index.html")]
        assert 'href="parametrized_t01/index.html"' in section
        text = Path(paths[0]).read_text()
        assert 'href="../index.html"' in text              # back to the index of the suite
        self.assert_balanced(text)
