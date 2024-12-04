# coding: utf-8

def get_ax_fig_plt(ax=None, **kwargs):
    """
    Helper function used in plot functions supporting an optional Axes argument.
    If ax is None, we build the `matplotlib` figure and create the Axes else
    we return the current active figure.

    Args:
        kwargs: keyword arguments are passed to plt.figure if ax is not None.

    Returns:
        ax: :class:`Axes` object
        figure: matplotlib figure
        plt: matplotlib pyplot module.
    """
    import matplotlib.pyplot as plt
    if ax is None:
        fig = plt.figure(**kwargs)
        ax = fig.add_subplot(1, 1, 1)
    else:
        fig = plt.gcf()

    return ax, fig, plt


def add_fig_kwargs(func):
    """
    Decorator that adds keyword arguments for functions returning matplotlib
    figures.

    The function should return either a matplotlib figure or None to signal
    some sort of error/unexpected event.
    See doc string below for the list of supported options.
    """
    from functools import wraps

    @wraps(func)
    def wrapper(*args, **kwargs):
        # pop the kwds used by the decorator.
        title = kwargs.pop("title", None)
        size_kwargs = kwargs.pop("size_kwargs", None)
        show = kwargs.pop("show", True)
        savefig = kwargs.pop("savefig", None)
        tight_layout = kwargs.pop("tight_layout", False)
        ax_grid = kwargs.pop("ax_grid", None)
        ax_annotate = kwargs.pop("ax_annotate", None)

        # Call func and return immediately if None is returned.
        fig = func(*args, **kwargs)
        if fig is None:
            return fig

        # Operate on matplotlib figure.
        if title is not None:
            fig.suptitle(title)

        if size_kwargs is not None:
            fig.set_size_inches(size_kwargs.pop("w"), size_kwargs.pop("h"),
                                **size_kwargs)

        if ax_grid is not None:
            for ax in fig.axes:
                ax.grid(bool(ax_grid))

        if ax_annotate:
            from string import ascii_letters
            tags = ascii_letters
            if len(fig.axes) > len(tags):
                tags = (1 + len(ascii_letters) // len(fig.axes)) * ascii_letters
            for ax, tag in zip(fig.axes, tags):
                ax.annotate("(%s)" % tag, xy=(0.05, 0.95), xycoords="axes fraction")

        if tight_layout:
            try:
                fig.tight_layout()
            except Exception as exc:
                # For some unknown reason, this problem shows up only on travis.
                # https://stackoverflow.com/questions/22708888/valueerror-when-using-matplotlib-tight-layout
                print("Ignoring Exception raised by fig.tight_layout\n", str(exc))

        if savefig:
            fig.savefig(savefig)

        if show:
            import matplotlib.pyplot as plt
            plt.show()

        return fig

    # Add docstring to the decorated method.
    s = "\n\n" + """\
        Keyword arguments controlling the display of the figure:

        ================  ====================================================
        kwargs            Meaning
        ================  ====================================================
        title             Title of the plot (Default: None).
        show              True to show the figure (default: True).
        savefig           "abc.png" or "abc.eps" to save the figure to a file.
        size_kwargs       Dictionary with options passed to fig.set_size_inches
                          e.g. size_kwargs=dict(w=3, h=4)
        tight_layout      True to call fig.tight_layout (default: False)
        ax_grid           True (False) to add (remove) grid from all axes in fig.
                          Default: None i.e. fig is left unchanged.
        ax_annotate       Add labels to  subplots e.g. (a), (b).
                          Default: False
        ================  ====================================================

"""

    if wrapper.__doc__ is not None:
        # Add s at the end of the docstring.
        wrapper.__doc__ += "\n" + s
    else:
        # Use s
        wrapper.__doc__ = s

    return wrapper


class MplExpose(object): # pragma: no cover
    """
    Example:

        with MplExpose() as e:
            e(obj.plot1(show=False))
            e(obj.plot2(show=False))
    """
    def __init__(self, slide_mode=False, slide_timeout=None, verbose=1):
        """
        Args:
            slide_mode: If true, iterate over figures. Default: Expose all figures at once.
            slide_timeout: Close figure after slide-timeout seconds Block if None.
            verbose: verbosity level
        """
        self.figures = []
        self.slide_mode = bool(slide_mode)
        self.timeout_ms = slide_timeout
        self.verbose = verbose
        if self.timeout_ms is not None:
            self.timeout_ms = int(self.timeout_ms * 1000)
            assert self.timeout_ms >= 0

        if self.verbose:
            if self.slide_mode:
                print("\nSliding matplotlib figures with slide timeout: %s [s]" % slide_timeout)
            else:
                print("\nLoading all matplotlib figures before showing them. It may take some time...")

        self.start_time = time.time()

    def __call__(self, obj):
        """
        Add an object to MplExpose. Support mpl figure, list of figures or
        generator yielding figures.
        """
        import types
        if isinstance(obj, (types.GeneratorType, list, tuple)):
            for fig in obj:
                self.add_fig(fig)
        else:
            self.add_fig(obj)

    def add_fig(self, fig):
        """Add a matplotlib figure."""
        if fig is None: return

        if not self.slide_mode:
            self.figures.append(fig)
        else:
            #print("Printing and closing", fig)
            import matplotlib.pyplot as plt
            if self.timeout_ms is not None:
                # Creating a timer object
                # timer calls plt.close after interval milliseconds to close the window.
                timer = fig.canvas.new_timer(interval=self.timeout_ms)
                timer.add_callback(plt.close, fig)
                timer.start()

            plt.show()
            fig.clear()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Activated at the end of the with statement. """
        self.expose()

    def expose(self):
        """Show all figures. Clear figures if needed."""
        if not self.slide_mode:
            print("All figures in memory, elapsed time: %.3f s" % (time.time() - self.start_time))
            import matplotlib.pyplot as plt
            plt.show()
            for fig in self.figures:
                fig.clear()

