"""
HTML reports of the parametrized tests, in the style of the pages written by runtests.py (testsuite._HTML_REPORT_CSS).

- instance_section(): the "Parametrized test" card of the page of an instance (<workdir>/index.html),
  added by BaseTest.write_html_report through the html_extra_sections hook;
- template page: one page per template, <suite workdir>/<suite>_<tNN>/index.html for runtests.py
  and <workdir root>/index.html for run_params.py, with all the instances, their parameters and statuses,
  the fixtures and the differences of the failed instances;
- write_suite_pages(): the template pages of a suite, and the card that links them from Test_suite/index.html.

The pages are built from the test objects (status, files_to_test, params), so the same code serves
runtests.py (after the workers returned their results, see ParametrizedMixin.get_results) and run_params.py.
Standard library only.
"""
from __future__ import annotations

import glob
import os
import re
import time
from html import escape
from typing import Any

from ..testsuite import _HTML_REPORT_CSS, _MY_NAME, ChainOfTests, sec2str, status2html
from . import fixture_cache
from .paramtests import BrokenTemplateTest, all_fixtures

#: Statuses of runtests.py, in the order of the summary tables.
STATUSES = ("failed", "passed", "succeeded", "skipped", "disabled")

#: Maximum number of lines of fldiff details shown for a failed file.
MAX_DETAIL_LINES = 80

#: Number of failed instances whose differences are expanded in the template page.
MAX_OPEN_FAILURES = 3

TEMPLATE_SUFFIX = ".tmpl.abi"

_FILTER_JS = """
<script>
(function () {
  const rows = Array.from(document.querySelectorAll("#instances tbody tr"));
  let active = null;
  document.querySelectorAll(".filters button").forEach(b => b.addEventListener("click", () => {
    active = active === b.dataset.status ? null : b.dataset.status;
    document.querySelectorAll(".filters button").forEach(x => x.classList.toggle("active", x.dataset.status === active));
    rows.forEach(r => { r.style.display = !active || r.dataset.status === active ? "" : "none"; });
  }));
})();
</script>
"""


def template_stem(path: str) -> str:
    """t03 for .../Input/t03.tmpl.abi."""
    name = os.path.basename(path)
    return name[:-len(TEMPLATE_SUFFIX)] if name.endswith(TEMPLATE_SUFFIX) else os.path.splitext(name)[0]


def template_of(test) -> str:
    """Path of the template of an instance, a prerequisite or a template that could not be expanded."""
    return getattr(test, "template_path", None) or test.inp_fname


def page_dir_name(suite_name: str, template_path: str) -> str:
    """Directory of the template page in the workdir of the suite (runtests.py), next to the test workdirs."""
    return f"{suite_name}_{template_stem(template_path)}"


def file_options(f) -> str:
    """Comparison mode and tolerances of a file to test."""
    opts = [getattr(f, "param_mode", "fldiff")] + list(f.fld_options)
    opts += [f"tolnlines={f.tolnlines}", f"tolabs={f.tolabs:g}", f"tolrel={f.tolrel:g}"]
    return " ".join(opts)


def _link(label: str, href: str) -> str:
    return f'<a href="{escape(href, quote=True)}">{escape(label)}</a>'


def _pills(params: dict[str, Any]) -> str:
    return "".join(f'<span class="pill">{escape(str(k))} = {escape(str(v))}</span>' for k, v in params.items())


def status_of(test) -> str:
    """The runtests.py status of a test; a test that did not run is skipped."""
    if getattr(test, "_status", None):
        return test._status
    if not getattr(test, "_executed", False):
        return "skipped"
    return test.status


def requested_nproc(test) -> int:
    """Number of MPI processes the test must run on (prerequisites and fixtures: the processes of the run)."""
    if getattr(test, "is_prerequisite", False):
        return getattr(test, "mpi_nprocs", 1)
    return test.nprocs_to_test[0] if getattr(test, "nprocs_to_test", None) else 1


def _nproc_html(test) -> str:
    want, ran = requested_nproc(test), getattr(test, "mpi_nprocs_run", None)
    if ran is None:
        return str(want)
    if ran == want:
        return f'{want} <span class="muted" title="mpi_nproc reported in the main output">(ran {ran})</span>'
    return f'{want} <span class="bad" title="mpi_nproc reported in the main output">(ran {ran})</span>'


def _role_pill(test) -> str:
    if getattr(test, "is_prerequisite", False):
        return '<span class="pill">prerequisite</span>'
    if getattr(test, "instance_index", None) == 0:
        return '<span class="pill pill-ref" title="Its outputs are the golden references (--update-ref)">reference</span>'
    return ""


def _read_head(path: str, nlines: int = MAX_DETAIL_LINES) -> str:
    try:
        with open(path, errors="replace") as fh:
            lines = fh.readlines()
    except OSError:
        return ""
    text = "".join(lines[:nlines])
    if len(lines) > nlines:
        text += f"... ({len(lines) - nlines} more lines in {os.path.basename(path)})\n"
    return text


# ----------------------------------------------------------------------------------------------------------------
# Page of an instance
# ----------------------------------------------------------------------------------------------------------------

def instance_section(test) -> str:
    """The "Parametrized test" card of the page of an instance (or of a prerequisite)."""
    template = template_of(test)
    stem = template_stem(template)
    href = test.template_page_href or f"../{page_dir_name(test.suite_name, template)}/index.html"

    if getattr(test, "is_prerequisite", False):
        title, role = f"Prerequisite of {stem}", "prerequisite (setup.test_chain): runs before the instances"
    else:
        title = f"Instance {test.instance_index:04d} of {stem}"
        role = ("reference instance: its outputs are the golden references of all the instances"
                if test.instance_index == 0 else "instance, compared with the golden references of instance 0000")

    items = [("Template", f'<span class="mono">{escape(template)}</span>'), ("Role", escape(role))]
    if test.params:
        items.append(("Parameters", _pills(test.params)))
    ran = getattr(test, "mpi_nprocs_run", None)
    nproc = f"requested {requested_nproc(test)}"
    if ran is not None:
        nproc += f", reported by ABINIT (mpi_nproc) {ran}"
    items.append(("MPI processes", nproc))
    refs = [_link(f.name, f"refs/{f.name}") + f' <span class="muted">&rarr; {escape(getattr(f, "ref_name", f.name))}</span>'
            for f in test.files_to_test]
    if refs:
        items.append(("Golden references", "<br>".join(refs)))
    fixtures = all_fixtures(list(getattr(test, "fixtures", ()) or ()))
    if fixtures:
        items.append(("Fixtures", " ".join(f'<span class="pill">{escape(fx.name)} (np={fx.nproc})</span>'
                                           for fx in fixtures)))
    rows = "".join(f"<dt>{k}</dt><dd>{v}</dd>" for k, v in items)

    groups = []
    for f in test.files_to_test:
        for g in getattr(f, "yaml_groups", None) or []:
            tol = " ".join(f"{k}={v}" for k, v in g["tolerances"].items())
            groups.append(f'<tr><td class="mono">{escape(f.name)}</td><td class="mono">{escape(g["path"])}</td>'
                          f'<td>{status2html(g["status"])}</td><td class="mono">{escape(tol)}</td>'
                          f'<td class="mono">{escape(g["msg"])}</td></tr>')
    groups_html = ""
    if groups:
        groups_html = ('<div class="table-wrap"><table><thead><tr><th>File</th><th>yaml_docs group</th><th>Status</th>'
                       '<th>Tolerances</th><th>Result</th></tr></thead><tbody>' + "".join(groups) +
                       '</tbody></table></div>'
                       '<div class="card-body muted">The text and HTML diffs of yaml_docs files show the compared '
                       'lines (the YAML documents without the ignored entries), not the whole files.</div>')

    return f"""
    <section class="card">
     <header class="card-header">
      <div><p class="eyebrow">Parametrized test</p><h2 class="report-title">{escape(title)}</h2></div>
      {_link(f"All the instances of {stem} →", href)}
     </header>
     <div class="card-body"><dl class="meta-list">{rows}</dl></div>
     {groups_html}
     <details><summary>Generated input</summary><pre class="log">{escape(test.input_text)}</pre></details>
    </section>
    """


# ----------------------------------------------------------------------------------------------------------------
# Page of a template
# ----------------------------------------------------------------------------------------------------------------

def _members(tests) -> list:
    out = []
    for t in tests:
        out.extend(list(t) if isinstance(t, ChainOfTests) else [t])
    return out


def fixture_statuses(root: str, fixtures) -> list[dict[str, Any]]:
    """The most recent entry of each fixture in the cache `root` (status.json), or a placeholder."""
    out = []
    for fx in fixtures:
        entries = sorted(glob.glob(os.path.join(root, fixture_cache.entry_name(fx.name, "*"))), key=os.path.getmtime)
        entries = [e for e in entries if os.path.isdir(e) and ".tmp-" not in e]
        st = fixture_cache.read_status(entries[-1]) if entries else None
        out.append({"name": fx.name, "nproc": fx.nproc, "entry": entries[-1] if entries else "",
                    "status": (st or {}).get("status", "not built"), "nproc_run": (st or {}).get("nproc"),
                    "build_time": (st or {}).get("build_time"), "message": (st or {}).get("message", "")})
    return out


def _description(test) -> str:
    """The description of the template: the instances append their parameters, e.g. `... [useylm=1]`."""
    text = (getattr(test, "description", "") or "").strip()
    if getattr(test, "params", None):
        text = re.sub(r"\s*\[[^\]]*\]$", "", text)
    return text


def render_template_page(template: str, tests, page_dir: str, back_href: str | None = None,
                         fixtures_root: str | None = None, command: str = "") -> str:
    """
    The page of the template `template` from its test objects `tests` (instances, chains, broken template).
    Links are relative to `page_dir`, the directory where the page is written.
    """
    members = _members(tests)
    stem = template_stem(template)
    counts = {s: sum(1 for m in members if status_of(m) == s) for s in STATUSES}

    # One column per parameter, in the order of appearance.
    keys: list[str] = []
    for m in members:
        for k in getattr(m, "params", {}) or {}:
            if k not in keys and k != "nproc":   # nproc: the np column
                keys.append(k)

    rows, failures = [], []
    for m in members:
        status = status_of(m)
        page = os.path.join(m.workdir, "index.html") if getattr(m, "workdir", None) else ""
        ident = _link(m.id, os.path.relpath(page, page_dir)) if page and os.path.isfile(page) else escape(m.id)
        params = getattr(m, "params", {}) or {}
        cells = "".join(f'<td class="mono">{escape(str(params[k])) if k in params else ""}</td>' for k in keys)
        files = " ".join(f'<span class="pill" title="{escape(f.fld_msg, quote=True)}">{escape(f.name)}: '
                         f'{escape(f.fld_status)}</span>' for f in m.files_to_test) if status != "skipped" else ""
        if isinstance(m, BrokenTemplateTest):
            files = f'<span class="mono">{escape(m.error_message)}</span>'
        etime = getattr(m, "run_etime", None)
        rows.append(f'<tr data-status="{status}"><td>{ident} {_role_pill(m)}</td><td>{status2html(status)}</td>'
                    f'<td class="num">{_nproc_html(m)}</td>{cells}<td>{files}</td>'
                    f'<td class="num">{sec2str(etime) if status != "skipped" and etime is not None else ""}</td></tr>')
        if status == "failed":
            failures.append(_failure_details(m, page_dir, is_open=len(failures) < MAX_OPEN_FAILURES))

    heads = "".join(f"<th>{escape(k)}</th>" for k in keys)
    first = next((m for m in members if not isinstance(m, BrokenTemplateTest)), members[0])
    buttons = "".join(f'<button data-status="{s}">{status2html(s)} {n}</button>' for s, n in counts.items() if n)

    meta = [("Description", escape(_description(first))),
            ("Template", f'<span class="mono">{escape(template)}</span>'),
            ("Authors", escape(", ".join(getattr(first, "authors", []) or []))),
            ("Keywords", escape(", ".join(getattr(first, "keywords", []) or [])))]
    refs = [f'<span class="mono">{escape(os.path.join(first.suite_refs_dir, getattr(f, "ref_name", f.name)))}</span>'
            for f in getattr(first, "files_to_test", []) if getattr(first, "suite_refs_dir", None)]
    if refs:
        meta.append(("Golden references", "<br>".join(refs) +
                     '<br><span class="muted">produced by instance 0000 (run_params.py run --update-ref)</span>'))
    if command:
        meta.append(("Command", f'<span class="mono">{escape(command)}</span>'))
    meta_html = "".join(f"<dt>{k}</dt><dd>{v}</dd>" for k, v in meta)

    try:
        with open(template, errors="replace") as fh:
            source = fh.read()
    except OSError:
        source = ""

    fixtures = all_fixtures([fx for m in members for fx in (getattr(m, "fixtures", ()) or ())])
    fixtures_html = ""
    if fixtures and fixtures_root:
        frows = "".join(_fixture_row(s) for s in fixture_statuses(fixtures_root, fixtures))
        fixtures_html = f"""
        <section class="card">
         <header class="card-header"><h2 class="report-title">Fixtures</h2></header>
         <div class="table-wrap"><table><thead><tr><th>Fixture</th><th>Status</th><th>np</th><th>Build time (s)</th>
         <th>Message / cache entry</th></tr></thead><tbody>{frows}</tbody></table></div>
        </section>"""

    failures_html = ""
    if failures:
        failures_html = f"""
        <section class="card">
         <header class="card-header"><h2 class="report-title">Differences</h2></header>
         {"".join(failures)}
        </section>"""

    back = f'<a class="back-link" href="{escape(back_href)}">&#8592; Back</a>' if back_href else ""
    etimes = [m.run_etime for m in members if status_of(m) != "skipped" and getattr(m, "run_etime", None) is not None]
    nprocs = sorted({requested_nproc(m) for m in members if not getattr(m, "is_prerequisite", False)})
    return f"""<!doctype html>
<html lang="en">
<head>
 <meta charset="utf-8">
 <meta name="viewport" content="width=device-width, initial-scale=1">
 <title>{escape(stem)} · Parametrized test</title>
 <style>{_HTML_REPORT_CSS}</style>
</head>
<body>
<main class="report">
 {back}
 <section class="card">
  <header class="card-header"><div><p class="eyebrow">Parametrized test</p><h1 class="report-title">{escape(stem)}</h1></div>
  {status2html("failed" if counts["failed"] else "succeeded" if counts["succeeded"] + counts["passed"] else "skipped")}</header>
  <div class="table-wrap"><table>
   <tr>{"".join(f"<th>{status2html(s)}</th>" for s in STATUSES)}</tr>
   <tr>{"".join(f"<td>{counts[s]}</td>" for s in STATUSES)}</tr>
  </table></div>
  <div class="card-body"><dl class="metrics">
   <div class="metric"><dt>Instances</dt><dd>{sum(1 for m in members if not getattr(m, "is_prerequisite", False))}</dd></div>
   <div class="metric"><dt>MPI processes</dt><dd>{", ".join(map(str, nprocs))}</dd></div>
   <div class="metric"><dt>Run time</dt><dd>{sec2str(sum(etimes))} s</dd></div>
  </dl></div>
 </section>
 <section class="card">
  <header class="card-header"><h2 class="report-title">Template</h2></header>
  <div class="card-body"><dl class="meta-list">{meta_html}</dl></div>
  <details><summary>Source of {escape(os.path.basename(template))}</summary><pre class="log">{escape(source)}</pre></details>
 </section>
 <section class="card">
  <header class="card-header"><h2 class="report-title">Instances</h2></header>
  <div class="filters">{buttons}</div>
  <div class="table-wrap"><table id="instances">
   <thead><tr><th>Instance</th><th>Status</th><th title="requested (reported by ABINIT)">np</th>{heads}<th>Files</th><th>Run time (s)</th></tr></thead>
   <tbody>{"".join(rows)}</tbody>
  </table></div>
 </section>
 {fixtures_html}
 {failures_html}
 <footer class="report-footer">Automatically generated by {escape(_MY_NAME)} on {time.asctime()}</footer>
 <a class="back-to-top" href="#" aria-label="Back to top">↑ Top</a>
</main>
{_FILTER_JS}
</body>
</html>
"""


def _fixture_row(s: dict[str, Any]) -> str:
    status = {"ok": "succeeded", "failed": "failed"}.get(s["status"], "skipped")
    nproc = str(s["nproc"])
    if s["nproc_run"] not in (None, s["nproc"]):
        nproc += f' <span class="bad">(ran {s["nproc_run"]})</span>'   # fallback to fewer processes
    build_time = "" if s["build_time"] is None else sec2str(s["build_time"])
    return (f'<tr><td class="mono">{escape(s["name"])}</td><td>{status2html(status)}</td><td class="num">{nproc}</td>'
            f'<td class="num">{build_time}</td>'
            f'<td class="mono">{escape(s["message"]) + "<br>" if s["message"] else ""}'
            f'<span class="muted">{escape(s["entry"])}</span></td></tr>')


def _failure_details(test, page_dir: str, is_open: bool = True) -> str:
    """The messages and the fldiff details of a failed instance, with links to its diffs."""
    parts = []
    errors = getattr(test, "result_errors", None) or [str(e) for e in getattr(test, "exceptions", [])]
    if isinstance(test, BrokenTemplateTest):
        errors = [test.error_message]
    for e in errors:
        parts.append(f'<div class="card-body"><div class="message">{escape(e)}</div></div>')
    workdir = getattr(test, "workdir", "") or ""
    for f in test.files_to_test:
        if f.fld_status != "failed" or not workdir:
            continue
        links = [_link(label, os.path.relpath(os.path.join(workdir, f.name + ext), page_dir))
                 for label, ext in (("fldiff", ".fldiff"), ("diff", ".diff"), ("HTML diff", ".diff.html"))
                 if os.path.isfile(os.path.join(workdir, f.name + ext))]
        details = _read_head(os.path.join(workdir, f.name + ".fldiff"))
        parts.append(f'<div class="card-body"><b class="mono">{escape(f.name)}</b> '
                     f'<span class="muted">({escape(file_options(f))})</span>: {escape(f.fld_msg)}'
                     f'<nav class="actions">{" ".join(links)}</nav></div>'
                     + (f'<pre class="log">{escape(details)}</pre>' if details else ""))
    return (f'<details{" open" if is_open else ""}><summary>{escape(test.id)} <span class="muted">{escape(_param_text(test))}</span></summary>'
            + "".join(parts) + "</details>")


def _param_text(test) -> str:
    return ", ".join(f"{k}={v}" for k, v in (getattr(test, "params", {}) or {}).items())


def write_template_page(path: str, template: str, tests, back_href: str | None = None,
                        fixtures_root: str | None = None, command: str = "") -> str:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(render_template_page(template, tests, os.path.dirname(os.path.abspath(path)),
                                      back_href=back_href, fixtures_root=fixtures_root, command=command))
    return path


# ----------------------------------------------------------------------------------------------------------------
# Suite (runtests.py)
# ----------------------------------------------------------------------------------------------------------------

def group_by_template(tests) -> dict[str, list]:
    """{template path: [tests]} for the parametrized tests (and chains) of `tests`, in order of appearance."""
    groups: dict[str, list] = {}
    for t in tests:
        members = list(t) if isinstance(t, ChainOfTests) else [t]
        param = [m for m in members if getattr(m, "template_path", None) or isinstance(m, BrokenTemplateTest)]
        if param:
            groups.setdefault(template_of(param[0]), []).append(t)
    return groups


def write_suite_pages(suite, workdir: str) -> tuple[list[str], str]:
    """
    Write the page of each parametrized template of `suite` in <workdir>/<suite>_<tNN>/index.html.
    Return the paths of the pages (for the tarball) and the card of Test_suite/index.html that links them.
    """
    paths, rows = [], []
    for template, tests in group_by_template(suite).items():
        members = _members(tests)
        name = page_dir_name(members[0].suite_name, template)
        path = write_template_page(os.path.join(workdir, name, "index.html"), template, tests, back_href="../index.html",
                                   fixtures_root=os.path.join(workdir, "fixtures"))
        paths.append(path)
        counts = {s: sum(1 for m in members if status_of(m) == s) for s in STATUSES}
        rows.append(f'<tr><td>{_link(template_stem(template), f"{name}/index.html")}</td>'
                    f'<td>{escape(_description(members[0]))}</td>'
                    + "".join(f'<td class="num">{counts[s]}</td>' for s in STATUSES) + "</tr>")
    if not rows:
        return paths, ""
    heads = "".join(f"<th>{status2html(s)}</th>" for s in STATUSES)
    return paths, f"""
    <section class="card">
     <header class="card-header"><div><p class="eyebrow">Parametrized tests</p><h2 class="report-title">Templates</h2></div></header>
     <div class="table-wrap"><table><thead><tr><th>Template</th><th>Description</th>{heads}</tr></thead>
     <tbody>{"".join(rows)}</tbody></table></div>
    </section>
    """


__all__ = ["file_options", "instance_section", "render_template_page", "write_suite_pages", "write_template_page"]
