"""Panel dashboard to analyze the data stored in the mocc files (memory allocation info)."""
from fkiss.termcolor import cprint

try:
    import panel as pn
    import param
except ImportError as exc:
    cprint("Use `conda install panel` or `pip install panel` to install the python package.", "red")
    raise exc


def _df(df):
    """
    Create a Panel DataFrame widget.

    Args:
        df (pd.DataFrame): The dataframe to display.

    Returns:
        pn.widgets.DataFrame: The initialized widget.
    """
    return pn.widgets.DataFrame(df, disabled=True)


class MoccViewer(param.Parameterized):
    """
    A Dashboard to browse and visualize the data stored in the .mocc file.
    It can can be executed either inside a jupyter notebook or as a standalone bokeh app.
    """

    #engine = pn.widgets.Select(value="dot",
    #    options=['dot', 'neato', 'twopi', 'circo', 'fdp', 'sfdp', 'patchwork', 'osage'])

    def __init__(self, mocc, **params):
        """
        Initialize the MoccViewer.

        Args:
            mocc (AbimemFile): The memory profile data object.
            **params: Additional parameters for Parameterized.
        """
        super().__init__(**params)
        self.mocc = mocc

    def get_panel(self):
        """
        Build the Panel dashboard.

        Returns:
            pn.Tabs: The dashboard container.
        """
        tabs = pn.Tabs(); app = tabs.append
        mocc = self.mocc

        gspec = pn.GridSpec(sizing_mode="scale_width")
        gspec[0, 0] = mocc.plot_memory_usage(show=False, title="Memory usage")
        gspec[0, 1] = mocc.plot_hist(show=False, title="Allocation histogram")
        gspec[1, :] = mocc.plot_peaks(maxlen=10, title="Memory peaks", show=False)

        maxlen = 50
        df = mocc.get_peaks(maxlen=maxlen, as_dataframe=True)
        df.drop(columns=["locus", "line", "action", "ptr"], inplace=True)
        col = pn.Column(gspec,
                        f"## DataFrame with the first {maxlen} peaks",
                        _df(df),
                        sizing_mode="scale_width")
        app(("Plots", col))

        #app(("DataFrame", _df(mocc.dataframe)))
        hotdf = mocc.get_hotspots_dataframe()
        ax = hotdf.plot.pie(y="malloc_mb")
        import matplotlib
        fig = matplotlib.pyplot.gcf()
        app(("Hotspots",
             pn.Column(
                "### DataFrame with total memory allocated per Fortran file.",
                fig,
                _df(hotdf),
                sizing_mode="scale_width",
        )))
        app(("Intense",
             pn.Column(
                "### DataFrame with variables that are allocated/freed many times.",
                _df(mocc.get_intense_dataframe()),
                sizing_mode="scale_width")
        ))

        #retcode = memfile.find_memleaks()
        #app(("Memleaks", _df(mocc.find_memleaks())))

        return tabs
