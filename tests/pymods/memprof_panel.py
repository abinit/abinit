"""Panel dashboard to analyze the data stored in the mocc files (memory allocation info)."""
import os
from fkiss.termcolor import cprint

try:
    import param
    import panel as pn
except ImportError as exc:
    cprint("Use `conda install panel` or `pip install panel` to install the python package.", "red")
    raise exc

import bokeh.models.widgets as bkw

def _df(df):
    return pn.widgets.DataFrame(df, disabled=True)


class MoccViewer(param.Parameterized):
    """
    A Dashboard to browse and visualize the data stored in the .mocc file.
    It can can be executed either inside a jupyter notebook or as a standalone bokeh app.
    """

    #engine = pn.widgets.Select(value="dot",
    #    options=['dot', 'neato', 'twopi', 'circo', 'fdp', 'sfdp', 'patchwork', 'osage'])

    def __init__(self, mocc, **params):
        super().__init__(**params)
        self.mocc = mocc

    #@param.depends('dir_select.value')
    #def view_dirname(self):
    #    dirpath = self.dirname2path[self.dir_select.value]
    #    # Update widgets.
    #    self.file_select.options = [f.name for f in self.dir2files[dirpath]]
    #    if hasattr(self, "tabs"): self.tabs.active = 0

    #    return pn.Row(_df(self.proj.get_stats_dir(dirpath)),
    #                  self.proj.get_graphviz_dir(dirpath, engine=self.engine.value),
    #                  sizing_mode="scale_width")

    def get_panel(self):
        """Return tabs with widgets to interact with the DDB file."""
        tabs = pn.Tabs(); app = tabs.append
        mocc = self.mocc

        gspec = pn.GridSpec() #sizing_mode='scale_width')
        gspec[0, 0] = mocc.plot_memory_usage(show=False)
        gspec[1, 0] = mocc.plot_peaks(maxlen=10, show=False)
        app(("Plots", gspec))

        #app(("DataFrame", _df(mocc.dataframe)))
        #app(("Peaks", _df(mocc.get_peaks())))
        app(("Hotspots", _df(mocc.get_hotspots_dataframe())))
        app(("Intense", _df(mocc.get_intense_dataframe())))
        #app(("Memleaks", _df(mocc.find_memleaks())))

        row = pn.Row(bkw.PreText(text=mocc.to_string(verbose=0), sizing_mode="scale_both"))
        app(("Summary", row))
        #app(("Ph-bands", pn.Row(
        #    pn.Column("# PH-bands options",
        #              *[self.param[k] for k in ("nqsmall", "ndivsm", "asr", "chneut", "dipdip", "lo_to_splitting")],
        #              self.temp_range, self.plot_phbands_btn),
        #    self.plot_phbands_and_phdos)
        #))
        #app(("BECs", pn.Row(
        #    pn.Column("# Born effective charges options",
        #             *[self.param[k] for k in ("asr", "chneut", "dipdip", "gamma_ev")], self.get_epsinf_btn),
        #    self.get_epsinf)
        #))
        #app(("Global", pn.Column("# Global parameters", *[self.param[k] for k in ("units", "mpi_procs", "verbose")])))

        return tabs
