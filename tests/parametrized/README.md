# Parametrized tests

A parametrized test runs **the same calculation with several settings that must not change the results** (algorithms,
MPI layouts, CPU or GPU) and checks that all of them agree with **one** reference.

## Motivations

In the standard test suite, a test is one input file with one reference file per configuration, and each run is
compared only with its own frozen output. This has limitations:

- **One configuration, one input.** To check that two algorithms give the same results (e.g. `useylm 0` and `1`,
  two eigensolvers), one writes several datasets or several inputs that differ by one variable (`v9/t133` spells out
  the same GS + DFPT twice), each with its own reference. Nothing checks that the results agree with each other:
  each reference can drift on its own, and a bug that affects only one path is frozen in that path's reference when
  the references are updated. Moreover, the datasets of an input run one after the other in a single process.
- **MPI parallelism needs one reference file per number of processes.** In the `paral` and `mpiio` suites, a test
  lists its process counts in `nprocs_to_test`, and each count has its own `[NCPU_N]` section and its own reference:
  `mpiio/t01` has `t01_MPI1.abo`, `t01_MPI2.abo` and `t01_MPI4.abo`; `paral/Refs` holds 292 files, and 54 tests have
  several per-count references (13 for `paral/t81`). A parallel run is never compared with the sequential one, every
  reference must be regenerated and reviewed separately, and the distribution of the work (`npkpt`, `npband`, `npfft`,
  `eph_np_pqbks`, ...) is fixed in the input: testing another distribution means another input and more references.
- **The GPU versions are tested against GPU references.** The GPU suites (`gpu_omp`, ...) compare a GPU run with a
  reference produced by a GPU run: they detect changes of the GPU results, not differences between the GPU and the
  CPU implementations.

A parametrized test addresses these points with a single template and a single set of references:

- the swept settings are declared in the template and expanded into several tests, the *instances*;
- **the instances are independent tests, so they run concurrently** (`runtests.py -j N` distributes them over the workers),
  whereas the datasets of an input run sequentially in one process;
- **every instance is compared with the same golden reference**, produced by the first instance (`0000`, usually the
  simplest configuration: sequential, default algorithm, CPU), so the test checks that the code paths agree with each other;
- **the number of MPI processes is a parameter like the others** (`nproc`): one template runs the same input on
  1, 2, 4, ... processes, with different distributions of k-points, bands, FFT or perturbations, and checks every
  parallel run against the sequential reference (see [Changing the number of MPI processes](#changing-the-number-of-mpi-processes));
- **CPU versus GPU**: sweeping `gpu_option` (e.g. `[0, 2]`) compares the GPU implementation with the CPU one on the same
  input, instead of a GPU run with an old GPU output. Such a template requires a GPU build: list the CPP variables of the
  GPU suites in `test_info.setup.need_cpp_vars` (e.g. `[HAVE_GPU, HAVE_OPENMP_OFFLOAD]`), and compare with `yaml_docs`
  (the echo of `gpu_option` differs between the instances). No such template is in the suite yet.

## What a parametrized test is

An ABINIT input *template* is expanded into several tests, the *instances*, that run the same calculation with
different settings **that must not change the results**, for example:

- algorithms: `useylm 0/1`, eigensolvers (`wfoptalg`), `use_gemm_nonlop`, `rmm_diis`, `istwfk`;
- MPI layouts: `nproc` with `npkpt`/`npband`/`npfft`/`npspinor`, `eph_np_pqbks`;
- the GPU implementation: `gpu_option`.

The suite is **optional**: it needs `pyyaml` and `pydantic`, and runs only on request
(`runtests.py --with-parametrized`, or `with_parametrized` in `testbot.json`).
The design and its rationale are in [`../plans/parametrized_tests_plan.md`](../plans/parametrized_tests_plan.md);
this file describes how to use and write the tests.

## Layout

| Path | Content |
|---|---|
| `Input/tNN.tmpl.abi` | Templates. Listed in `__init__.py` (`inp_files`); numbered like the other suites, no information in the name. |
| `Input/<name>.fixture.abi` | Fixtures: steps that prepare files for the templates (merged DDB/DVDB, WFK, ...). Descriptive names, not listed in `inp_files`. |
| `Refs/tNN.abo`, `Refs/tNN<...>` | Golden references, one per compared file and template, produced by instance `0000`. |
| `../run_params.py` | Command line tool: `doc`, `fixtures`, `expand`, `run`, `compare`. |
| `../pymods/parametrized/` | Implementation (pydantic models, expansion, comparison, fixtures, reports) and its unit tests. |

The instances are generated on the fly and are **never committed**: they exist only in the workdirs.

## Running

```sh
cd _build/tests                                              # any build directory
../../tests/runtests.py parametrized --with-parametrized -j4          # 1-process instances
../../tests/runtests.py parametrized --with-parametrized -n 2 -j4     # plus the 2-process instances
../../tests/runtests.py "parametrized[5]" --with-parametrized         # the instances of t05 only
```

An instance with `nproc: N` runs only in the passes with exactly `N` MPI processes (`-n N`); in the other passes it is
reported as skipped. `run_params.py` runs every instance of one template in one call:

```sh
../../tests/run_params.py run ../../tests/parametrized/Input/t05.tmpl.abi --max-cpus 4
../../tests/run_params.py compare ../../tests/parametrized/Input/t05.tmpl.abi          # compare existing outputs again
../../tests/run_params.py compare ../../tests/parametrized/Input/t05.tmpl.abi -p failed  # open the diffs in $PATCHER
```

The exit status is non-zero if an instance differs or fails.

### HTML reports

`runtests.py` and `run_params.py` write the same pages, in the style of the other test reports:

- **One page per template**, with all its instances: status, MPI processes requested and actually used (`ran N`,
  read from the output), one column per swept parameter, the result of each compared file, the fixtures,
  and the fldiff excerpts of the failed instances with links to their diffs. The source of the template is included.
  `runtests.py` writes it in `Test_suite/parametrized_tNN/index.html` and lists the templates at the top of
  `Test_suite/index.html`; `run_params.py` writes it in the workdir root (default `run_params_workdir/index.html`).
- **One page per instance** (`<workdir>/index.html`): the page of a standard test plus a "Parametrized test" card
  with the parameters, the role (instance `0000` is the reference instance), the golden references, the fixtures,
  the result of each `yaml_docs` group and the generated input.
- **Diffs of `yaml_docs` files** (`.diff`, and `.diff.html` with `--make-html-diff`) show the compared lines,
  i.e. the YAML documents without the ignored entries, not the whole output files.

`run_params.py` also writes `comparison_report.json`, the machine-readable result. After `compare`, the pages are
rewritten with the new results (the run times are then 0: nothing ran).

## Writing a template

A template is an ABINIT input preceded by a `%YAML ... %ENDYAML` block:

```text
%YAML
test_info:
  setup: {executable: abinit}
  files:
    files_to_test:
      - {name: "{id}.abo", mode: yaml_docs, tolnlines: 0, tolabs: 1.1e-8, tolrel: 1.0e-8, fld_options: -medium}
  extra_info: {keywords: [NC], authors: [J. Doe],
               description: "Si diamond, GS: results must not depend on useylm. Motivation: ..."}
axes:
  - {type: grid, parameters: {useylm: [0, 1]}}
%ENDYAML
# Diamond Si at the ideal positions
acell 3*10.18
rprim  0.0 0.5 0.5   0.5 0.0 0.5   0.5 0.5 0.0
ntypat 1  znucl 14  natom 2  typat 1 1
xred  0.0 0.0 0.0   0.25 0.25 0.25
ecut 8.0
nband 4
kptopt 1  ngkpt 4 4 4  nshiftk 1  shiftk 0.0 0.0 0.0
nstep 50  tolvrs 1.0d-10
pp_dirpath "$ABI_PSPDIR"
pseudos "PseudosTM_pwteter/14si.pspnc"
```

- `test_info` is the TEST_INFO block of the ordinary tests written in YAML: same sections, keys and meaning.
  The expansion generates the INI block of each instance; `nprocs_to_test`, `psp_files` and `[NCPU_N]` are not allowed.
- `files_to_test`: every entry names its file, with `{id}` for the instance id (`{id}.abo` is the main output,
  `{id}o_DS2_DDB` a DDB). The reference is the name with `{id}` replaced by `tNN`: `Refs/tNN.abo`, `Refs/tNNo_DS2_DDB`.
  - `mode: yaml_docs` compares only the YAML documents of the file (energies, forces, stress, ...), as text, with fldiff.
    Preferred for the main output. `overrides: {ResultsGS.convergence: ignore}` ignores or loosens one document or key.
  - `mode: fldiff` (the default) compares the whole file, as in the ordinary tests. Use it for DDB files and for outputs
    without YAML results (EPH).
- `axes`: `grid` (lists of values, product or `zip`) or `explicit` (a list of combinations, for variables that change
  together). The instances are the product of the axes, minus `exclude`. `nproc` is reserved: the number of MPI
  processes of the instance, never written in the input.
- `dataset`: required when the body has `ndtset > 1` (see [Multi-dataset templates](#multi-dataset-templates)).
- `fixtures: [name, ...]`: files prepared once and linked into the workdir of every instance (see below).
- A swept variable must not appear in the body; the pseudopotentials (`pseudos`, `pp_dirpath`) cannot be swept.
- Put the **motivation** of the test in `description`: which code path it covers and why it is not covered elsewhere.

`../run_params.py doc` prints the complete reference of the `%YAML` block, generated from the pydantic models.

### Multi-dataset templates

A template can have several datasets (`ndtset`), e.g. a GS followed by DFPT or NSCF calculations. The `dataset` key
decides **at which level** the parametrization applies:

- `dataset: all`: the swept variables apply to **every dataset**. The expansion appends them without suffix.
  Use it for settings that must be the same throughout the calculation, e.g. `useylm` across GS and DFPT (`t03`, `t04`).
- `dataset: k`: the swept variables apply to **dataset `k` only**. The expansion appends them with the suffix `k`;
  the other datasets are identical in all the instances. Use it to sweep one step of a workflow, e.g. the
  eigensolver of an NSCF dataset, or the DFPT part only.

For example, with the three datasets of `t04` (GS, DDK, phonons + electric field) and `axes: [{type: grid, parameters: {useylm: [0, 1]}}]`,
the instance with `useylm 1` gets:

| `dataset` | Lines appended to the body |
|---|---|
| `all` | `useylm 1` (GS, DDK and phonons) |
| `3` | `useylm3 1` (phonons only) |

Rules:

- `dataset` is required when `ndtset > 1`, and not allowed otherwise; `k` must be one of the datasets of the body
  (`jdtset` is taken into account). `udtset` is not supported.
- A swept variable must not appear in the body in any form (plain, with a dataset suffix, or in series syntax).
- The YAML documents of **all** the datasets are compared, and so are the files of every dataset listed in `files_to_test`.
- Every instance runs all the datasets: the datasets that are not swept are recomputed in each instance.
  If they are expensive, compute them once in a [fixture](#fixtures) instead.

### Recipe

1. Choose a small system (a few seconds per instance) and the settings to sweep. Start from the closest example below.
2. `run_params.py expand Input/tNN.tmpl.abi` validates the template and lists the instances (`-o DIR` writes the inputs).
3. `run_params.py run Input/tNN.tmpl.abi --update-ref --max-cpus N` creates the references from instance `0000`
   and compares all the instances with them.
4. Look at the differences (`run_params_workdir/index.html`, `-p failed`) and set the tolerances: as tight as the instances allow.
   Prefer converging better (`tolvrs`, `tolwfr`) to loosening the tolerances.
5. **Sensitivity check**: change a value in one output (e.g. a force by 1e-7) and check that `compare` reports `differs`.
6. Add the template to `inp_files` in `__init__.py`, and run
   `runtests.py parametrized --with-parametrized -n N` for each process count used.

### Examples

| Template | What it shows |
|---|---|
| `t01` | Grid axis (`useylm` × `paral_kgb`) on 1 process; `overrides` for a solver-dependent document. |
| `t02` | Explicit axis with `nproc` and `npkpt`/`npband`/`npfft` (1 to 8 processes). |
| `t03` | Multi-dataset input (`dataset: all`), DFPT at Γ and at q ≠ 0, DDB files compared with fldiff. |
| `t04` | NC+SOC, `useylm` over GS, DDK and phonons + electric field (Born charges, dielectric tensor). |
| `t05` | Eigensolvers (CG, LOBPCG variants, RMM-DIIS, ChebFI) on NC+SOC, with the motivation in the description. |
| `t06` | PAW version of `t02`: `nproc` with `npkpt`/`npband`/`npfft`, FFT grids fixed by hand (`ngfft`, `ngfftdg`). |
| `t07` | EPH (ZPR) with fixtures from a git submodule, `eph_np_pqbks` on 1 to 8 processes, fldiff on the main output. |

## Changing the number of MPI processes

`nproc` in an axis sets the number of MPI processes of the instance. It is not an ABINIT variable: it goes in the
generated TEST_INFO, and the variables of the distribution go in the same axis entry, so they always match `nproc`.
Typical applications: check that a new parallelization level, or a change in the distribution of the work, does not
change the results; test the MPI paths of an algorithm that the sequential tests do not reach.

Example: GS of diamond Si with `paral_kgb 1`, on 1, 2 and 4 processes, distributing the k-points, the bands and the FFT:

```text
%YAML
test_info:
  setup: {executable: abinit}
  files:
    files_to_test:
      - {name: "{id}.abo", mode: yaml_docs, tolnlines: 0, tolabs: 1.1e-8, tolrel: 1.0e-8, fld_options: -medium}
  extra_info: {keywords: [NC], authors: [J. Doe],
               description: "Si diamond, GS with paral_kgb 1: results must not depend on the MPI distribution"}
axes:
  # The product npkpt * npband * npfft must be nproc.
  - type: explicit
    values:
      - {nproc: 1, npkpt: 1, npband: 1, npfft: 1}   # instance 0000: the sequential reference
      - {nproc: 2, npkpt: 1, npband: 2, npfft: 1}   # bands
      - {nproc: 4, npkpt: 2, npband: 2, npfft: 1}   # k-points and bands
      - {nproc: 4, npkpt: 1, npband: 2, npfft: 2}   # bands and FFT
%ENDYAML
# Diamond Si at the ideal positions
acell 3*10.18
rprim  0.0 0.5 0.5   0.5 0.0 0.5   0.5 0.5 0.0
ntypat 1  znucl 14  natom 2  typat 1 1
xred  0.0 0.0 0.0   0.25 0.25 0.25
ecut 8.0
nband 8
kptopt 1  ngkpt 4 4 4  nshiftk 1  shiftk 0.0 0.0 0.0
nstep 50  tolvrs 1.0d-10
paral_kgb 1
pp_dirpath "$ABI_PSPDIR"
pseudos "PseudosTM_pwteter/14si.pspnc"
```

- Each instance gets its own id (`ti0001_NN_MPI2`, `ti0002_NN_MPI4`, ...) and runs **only in the passes with exactly
  `nproc` processes**: `runtests.py -n 2` runs the 1- and 2-process instances, `-n 4` the 1- and 4-process ones,
  and the others are reported as skipped. testbot runs each instance once, in the pass with its number of processes.
  `run_params.py run --max-cpus 4` runs all of them in one call.
- `nproc` must be one of 1, 2, 4, 8, 10, 24, 64: the counts that testbot runs. Keep them small (1, 2, 4).
- After each run, the number of processes reported by ABINIT (`mpi_nproc` in the output) must be `nproc`, otherwise
  the instance fails: an instance cannot silently test another configuration than the one it claims.
- The same works for other drivers: `t07` sweeps `eph_np_pqbks` (perturbations, q-points, bands) on 1 to 8 processes.
- MPI and algorithms can be combined, e.g. an eigensolver axis crossed with an `nproc` axis, or `nproc` added to some
  entries of an explicit axis only.

## Comparing the outputs

Each file listed in `files_to_test` is compared with its reference with one of two modes. Both use fldiff, with
`tolnlines`/`tolabs`/`tolrel`/`fld_options` as in the ordinary tests; they differ in **what** is compared.

| Mode | Compares | Use it when |
|---|---|---|
| `fldiff` (default) | the whole file, line by line | the structure of the file does not depend on the swept settings |
| `yaml_docs` | only the YAML documents of the file | the swept settings change the rest of the output |

**How fldiff reads a file.** The first character of a line selects its treatment. Lines starting with `-` are ignored
and may even be missing in one of the files; lines starting with `P` are ignored by default; the other lines are
compared, numbers with a tolerance. ABINIT uses these marks for what depends on the parallel setup (memory estimates,
`mpi_nproc`, the distribution variables, e.g. `-    eph_np_pqbks  2 2 2 1 1`). The two files must have **the same
number of compared lines**: otherwise fldiff stops with a fatal error ("file 2 has more significant lines than file 1"),
whatever the tolerances.

**Whole-file `fldiff` is perfectly fine when the length of the output does not depend on the configuration.**
The typical case is a sweep over MPI distributions: the lines that change are marked with `-` or `P`, and everything
else, results included, must be identical. `t07` compares the whole main output of an EPH ZPR calculation over
8 `eph_np_pqbks` distributions this way. DDB files are also compared with `fldiff`: their structure depends only on
the perturbations (`t03`, `t04`); a header line such as `lmnmax`, which depends on `useylm`, counts as one differing line.

**Whole-file `fldiff` breaks when a swept setting changes what is printed**, for example:

- ABINIT echoes the input variables that differ from their default: with `useylm 1`, two lines (input echo and
  final echo) that the `useylm 0` instance does not have;
- eigensolvers and other algorithms converge in a different number of SCF iterations, and the main output prints
  one line per iteration.

**`yaml_docs` avoids the problem** by comparing only the YAML documents of the file (`ResultsGS`, `EnergyTerms`, ...),
which have the same structure in every configuration. The documents are extracted as text and compared with fldiff,
with two refinements:

- `DatasetInfo` and `BeginCycle.solver` are ignored by default (they describe the configuration, not the results);
- `overrides` ignore or loosen a single document or key, e.g. `{ResultsGS.convergence: ignore}`: the residuals
  depend on the solver even when the results agree.

**Files without YAML documents.** `yaml_docs` needs something to compare: if a file has no YAML document (apart from
the ignored ones), the comparison **fails** ("no YAML documents to compare"), it never passes silently. Today the GS writes its results as
YAML documents, but DFPT and EPH do not (they will, gradually). So:

| Output | Swept settings | Compare with |
|---|---|---|
| GS main output | anything (algorithms, solvers, MPI) | `yaml_docs` |
| DFPT | anything | the DDB files, `fldiff` (`t03`, `t04`); the main output with `yaml_docs` for the GS datasets |
| EPH, other drivers without YAML results | MPI distribution only | the main output, whole-file `fldiff` (`t07`) |
| EPH, other drivers without YAML results | settings that change what is printed | no reliable option yet: compare a file whose structure does not change, or wait for the YAML results |

With whole-file `fldiff`, the precision is the one of the printed numbers: the EPH self-energy table, for example,
has 3 decimals (meV), so the tolerances of `t07` accept a change of the last digit on at most 2 lines.

## Fixtures

A fixture is a step whose outputs several instances (or templates) need: merging the DDB and DFPT potentials of a
git submodule, an NSCF run for a WFK file, ... It is written in `Input/<name>.fixture.abi`:

```text
%YAML
test_info:
  setup: {executable: abinit, output_prefix: Si_444k_8b}
  extra_info: {keywords: [NC], authors: [J. Doe], description: "NSCF WFK of Si with 8 bands on the 4x4x4 k-mesh"}
fixture: {outputs: [Si_444k_8b_WFK], nproc: 4}
%ENDYAML
iscf -2  tolwfr 1.0d-18  nband 8
...
```

- The name of a fixture (the file name without `.fixture.abi`) uses **lowercase** letters, digits and `_`, and starts
  with a letter: on a case-insensitive filesystem (the default on macOS) names differing only by case would clash.
- `outputs` are literal file names, linked read-only into the workdir of each instance. They may describe the content
  (`MgO_444q_DDB`: 4×4×4 q-mesh) for readability; they are never parsed.
- `needs: [other_fixture]` for a fixture that reads the outputs of another one.
- `nproc` (default 1) is the number of MPI processes of the fixture. If the run cannot provide them, the fixture
  runs on fewer processes with a warning (`nproc_strict: yes` skips it instead). `mrgddb` and `mrgdv` need `nproc: 1`.
- Fixtures are built **once per workdir root**, before the tests (`Test_suite/fixtures` for runtests.py,
  `<workdir-root>/fixtures` for run_params.py), and shared by all the templates that use them.
  The cache is safe with `-j N`. `run_params.py run --rebuild-fixtures` removes it.
- `run_params.py fixtures` lists the available fixtures with their description, executable, `nproc`, outputs and the
  templates that use them; `run_params.py fixtures --doc` also prints the format of the fixture files.
- Chains of prerequisite templates (`setup.test_chain`) are also supported, but fixtures are preferred:
  they run once instead of once per process count, and testbot's static runmode does not skip them.

## Rules and pitfalls

- **Sweep only result-invariant settings.** Never physics or convergence parameters (`ecut`, `ngkpt`, `tsmear`, `nline`, ...),
  and not `mixprec` (single-precision FFTs change the results).
- **Converge what you compare.** With `tolvrs` the empty bands are converged differently by each eigensolver
  (`fermie` differs at 1e-8): to compare solvers, converge every band with `tolwfr` (`t05`).
- **fldiff flags a line only if both its absolute and relative differences exceed the threshold** of `fld_options`
  (`-medium`: 1.01e-8, `-easy`: 1.01e-5). `-easy` can hide real changes; `tolabs`/`tolrel` then bound the flagged lines.
- **Choose the comparison mode with care** (see [Comparing the outputs](#comparing-the-outputs)): whole-file `fldiff` stops
  with a fatal error if the swept settings change the number of printed lines; `yaml_docs` fails if the file
  has no YAML documents (DFPT and EPH today).
- **ChebFI** (`wfoptalg 111`) needs extra bands and `nbdbuf`; with NC pseudos and `nspinor 2` it also needs `cprj_in_memory 0`.
- **FFT distribution** (`npfft` > 1): the FFT grids must be divisible by `npfft`, so ABINIT may choose other grids than
  on 1 process (with PAW, a 30^3 instead of a 27^3 dense grid in `t06`) and the results change. Fix `ngfft` (and `ngfftdg`
  with PAW) in the body when you sweep `npfft`.
- **DFPT at Γ**: set `nqpt`/`qpt` in the response datasets only; set globally, they also apply to the GS dataset.
- **Process counts**: `nproc` must be one of 1, 2, 4, 8, 10, 24, 64 (the counts testbot runs). The report checks that
  ABINIT used the requested number (`mpi_nproc` in the output) and fails the instance otherwise.
- **Header lines that depend on the swept variable** (e.g. `lmnmax` in the DDB header with `useylm`) count as differing
  lines: allow them with `tolnlines`.

## Not supported yet

- Chains with more than 1 process in testbot's static runmode (the serial prerequisites are skipped): use fixtures.
- `run_params.py submit` (batch jobs) and the `pending` status.
- Whether to keep the `paral_info.max_nprocs` shorthand (one fixed process count for all the instances) is undecided:
  prefer an explicit `nproc` axis.
