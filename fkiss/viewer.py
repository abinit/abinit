"""Panel dashboard."""
import os
from fkiss.termcolor import cprint

try:
    import param
    import panel as pn
except ImportError as exc:
    cprint("Use `conda install panel` or `pip install panel` to install the python package.", "red")
    raise exc

def _df(df):
    return pn.widgets.DataFrame(df, disabled=True)


# Possible approach to display big SVG files: https://github.com/ariutta/svg-pan-zoom
#<script>
#document.getElementById('my-embed').addEventListener('load', function(){
#  // Will get called after embed element was loaded
#  svgPanZoom(document.getElementById('my-embed'));
#})
#</script>


class ProjectViewer(param.Parameterized):
    """
    A Dashboard to browse the source code, visualize connections among directories,
    files and procedures inside the Abinit project.
    Panel can can be executed either inside a jupyter notebook or as a standalone bokeh app.
    """

    engine = pn.widgets.Select(value="dot",
        options=['dot', 'neato', 'twopi', 'circo', 'fdp', 'sfdp', 'patchwork', 'osage'])

    def __init__(self, proj, **params):
        super().__init__(**params)
        self.proj = proj
        self._layout()

    def _layout(self):
        self.dir2files = self.proj.groupby_dirname()
        self.dirname2path = {os.path.basename(p): p for p in self.dir2files}
        self.dir_select = pn.widgets.Select(name="Directory", options=list(self.dirname2path.keys()))
        self.all_pubs = self.proj.get_all_public_procedures()

        width = 200
        self.file_select = pn.widgets.Select(name="Fortran File", width=width)
        self.pubproc_select = pn.widgets.Select(name="Public Procedure", width=width)
        self.datatype_select = pn.widgets.Select(name="Datatype", width=width)

        self.find_proc = pn.widgets.AutocompleteInput(name='Find Procedure', options=list(self.all_pubs.keys()),
                                                      placeholder='Enter procedure name', width=width)
        self.find_proc_btn = pn.widgets.Button(name="Find Procedure", button_type='primary', width=width)
        self.find_proc_btn.on_click(self.on_find_proc_btn)

        self.all_datatypes_and_fortfiles = self.proj.get_all_datatypes_and_fortfile()
        self.find_dtype = pn.widgets.AutocompleteInput(name='Find Datatype',
                                                       options=list(self.all_datatypes_and_fortfiles.keys()),
                                                       placeholder='Enter datatype name', width=width)
        self.find_dtype_btn = pn.widgets.Button(name="Find DataType", button_type='primary', width=width)
        self.find_dtype_btn.on_click(self.on_find_dtype_btn)

        self.tabs = pn.Tabs(
            ("Directory", self.view_dirname),
            ("File", self.view_fort_file),
            ("Procedure", self.view_pubproc),
            ("Datatype", self.view_datatype),
        )

        controllers = pn.Row(
                pn.Column(self.dir_select, self.file_select),
                pn.Column(self.pubproc_select, self.datatype_select),
                pn.Column(self.find_proc, self.find_proc_btn),
                pn.Column(self.find_dtype, self.find_dtype_btn),
                #pn.Column(self.engine),
                #, self.rerun_btn),
                #sizing_mode='scale_width'
        )

        self.panel = pn.Column(controllers, self.tabs, sizing_mode="scale_width")

    def _find_fort_file(self, dirpath):
        for fort_file in self.dir2files[dirpath]:
            if fort_file.name == self.file_select.value: return fort_file
        else:
            raise ValueError("Cannot find fortran file with name: `%s` in `%s`" % (
                             self.file_select.value, dirpath))

    @param.depends('dir_select.value')
    def view_dirname(self):
        dirpath = self.dirname2path[self.dir_select.value]
        # Update widgets.
        self.file_select.options = [f.name for f in self.dir2files[dirpath]]
        if hasattr(self, "tabs"): self.tabs.active = 0

        return pn.Row(_df(self.proj.get_stats_dir(dirpath)),
                      self.proj.get_graphviz_dir(dirpath, engine=self.engine.value),
                      sizing_mode="scale_width")

    @param.depends('file_select.value')
    def view_fort_file(self):
        dirpath = self.dirname2path[self.dir_select.value]
        fort_file = self._find_fort_file(dirpath)

        # Update widgets.
        self.pubproc_select.options = list(fort_file.all_public_procedures.keys())
        self.datatype_select.options = list(fort_file.all_datatypes.keys())
        if hasattr(self, "tabs"): self.tabs.active = 1

        return pn.Row(_df(fort_file.get_stats()),
                      fort_file.get_graphviz(engine=self.engine.value),
                      sizing_mode="scale_width")

    @param.depends('pubproc_select.value')
    def view_pubproc(self):
        pubname = self.pubproc_select.value
        if pubname is None: return
        obj = self.proj.find_public_entity(pubname)
        graph = self.proj.get_graphviz_pubname(pubname, engine=self.engine.value)
        if hasattr(self, "tabs"): self.tabs.active = 2
        return pn.Row(obj, graph, sizing_mode="scale_width")

    @param.depends('datatype_select.value')
    def view_datatype(self):
        typename = self.datatype_select.value
        if typename is None: return
        dirpath = self.dirname2path[self.dir_select.value]
        fort_file = self._find_fort_file(dirpath)
        #print("fort_file.path", fort_file.path)
        dtype = fort_file.all_datatypes[typename]
        if hasattr(self, "tabs"): self.tabs.active = 3
        return pn.Row(dtype)

    def on_find_proc_btn(self, event):
        pubname = self.find_proc.value
        if pubname is None: return # or pubname not in self.all_pubs: return
        proc = self.all_pubs[pubname]
        dirpath = proc.dirpath
        file_basename = os.path.basename(proc.path)
        fort_file = self.proj.fort_files[file_basename]

        # Update widgets.
        self.dir_select.value = os.path.basename(dirpath)
        self.file_select.options = [f.name for f in self.dir2files[dirpath]]
        self.pubproc_select.options = list(fort_file.all_public_procedures.keys())
        self.datatype_select.options = list(fort_file.all_datatypes.keys())

        self.pubproc_select.value = pubname
        if hasattr(self, "tabs"): self.tabs.active = 2

    def on_find_dtype_btn(self, event):
        dname = self.find_dtype.value
        dtype, fort_file = self.all_datatypes_and_fortfiles[dname]

        # Update widgets.
        self.dir_select.value = os.path.basename(fort_file.dirname)
        self.file_select.value = fort_file.basename
        self.datatype_select.value = dname
        if hasattr(self, "tabs"): self.tabs.active = 3
