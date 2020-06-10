from __future__ import print_function, division, unicode_literals

import time

from pprint import pprint
from itertools import groupby
from collections import namedtuple, deque
# OrderedDict was added in 2.7. ibm6 still uses python2.6
try:
    from collections import OrderedDict
except ImportError:
    from .ordereddict import OrderedDict

from .plotting import add_fig_kwargs, get_ax_fig_plt
from .tools import lazy_property


class Entry(namedtuple("Entry", "vname, ptr, action, size, file, line, tot_memory")):
    """
    vname: Variable name.
    prt: Address of variable.
    action: "A" for allocation, "D" for deallocation.
    size: Size of allocation in bits.
    file: Name of Fortran file in which allocation/deallocation is performed.
    line: Line number in file.
    tot_memory: Total memory in bits allocated so far.
    """

    @classmethod
    def from_line(cls, line):
        """Build entry from line."""
        vname = line[:59].strip().replace(" ", "")
        args = [vname] + line[59:].split()
        return cls(*args)

    def __new__(cls, *args):
        """Extends the base class adding type conversion of arguments."""
        # write(logunt,'(a,t60,a,1x,2(i0,1x),2(a,1x),2(i0,1x))')&
        # trim(vname), trim(act), addr, isize, trim(abimem_basename(file)), line, memtot_abi%memory
        return super(cls, Entry).__new__(cls,
        	vname=args[0],
        	action=args[1],
        	ptr=int(args[2]),
        	size=int(args[3]),
        	file=args[4],
                line=int(args[5]),
                tot_memory=int(args[6]),
        )

    def __repr__(self):
        return self.to_repr(with_addr=True)

    def to_repr(self, with_addr=True):
        if with_addr:
            return "<var=%s, %s@%s:%s, addr=%s, size_mb=%.3f>" % (
              self.vname, self.action, self.file, self.line, hex(self.ptr), self.size_mb)
        else:
            return "<var=%s, %s@%s:%s, size_mb=%.3f>" %  (
              self.vname, self.action, self.file, self.line, self.size_mb)

    @lazy_property
    def size_mb(self):
        """Size in Megabytes."""
        sign = {"A": +1, "D": -1}[self.action]
        return sign * self.size / (8 * 1024 ** 2)

    @lazy_property
    def tot_memory_mb(self):
        """Total memory in Mb."""
        return self.tot_memory / (8 * 1024 ** 2)

    @lazy_property
    def isalloc(self):
        """True if entry represents an allocation."""
        return self.action == "A"

    @lazy_property
    def isfree(self):
        """True if entry represents a deallocation."""
        return self.action == "D"

    @lazy_property
    def iszerosized(self):
        """True if this is a zero-sized alloc/free."""
        return self.size == 0

    @lazy_property
    def locus(self):
        """Location of the entry. This is (hopefully) unique."""
        #return "%s:%s" % (self.file, self.line)
        return "%s:%s@%s:%s" % (self.action, self.vname, self.file, self.line)

    def __hash__(self):
        return hash(self.locus, self.size)

    def __eq__(self, other):
        return self.locus == other.locus and self.size == other.size

    def __neq__(self, other):
        return not (self == other)

    def frees_onheap(self, other):
        if (not self.isfree) or other.isalloc: return False
        if self.size + other.size != 0: return False
        return True

    def frees_onstack(self, other):
        if (not self.isfree) or other.isalloc: return False
        if self.size + other.size != 0: return False
        if self.locus != other.locus: return False
        return True



def entries_to_dataframe(entries):
    """
    Convert list of entries to pandas DataFrame.
    """
    import pandas as pd
    rows, index = [], []
    for e in entries:
        rows.append(OrderedDict([
            ("locus", e.locus),
            ("vname", e.vname),
            ("file", e.file),
            ("line", e.line),
            ("action", e.action),
            ("size_mb", e.size_mb),
            ("tot_memory_mb", e.tot_memory_mb),
            ("ptr", e.ptr),
        ]))
        index.append(e.locus)

    return pd.DataFrame(rows, index=index, columns=list(rows[0].keys()))


class AbimemFile(object):
    def __init__(self, path):
        self.path = path

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        lines = []
        app = lines.append
        df = self.get_intense_dataframe()
        app(df.to_string())
        return "\n".join(lines)

    def find_small_allocs(self, nbytes=160):
        """Zero sized allocations are not counted."""
        smalles = []
        for e in self.all_entries:
            if not e.isalloc: continue
            if 0 < e.size <= nbytes: smalles.append(e)

        pprint(smalles)
        return smalles

    def get_intense_dataframe(self):
        """
        Return DataFrame with intensive spots i.e. variables that are allocated/freed many times.
        """
        df = self.dataframe
        index, rows = [], []
        for locus, g in self.dataframe.groupby(by="locus"):
            this_action = g.action.values[0]
            assert all(g.action.values == this_action)
            malloc_mb = g.size_mb.sum()
            rows.append(OrderedDict([
                ("ncalls", len(g)),
                ("malloc_mb", malloc_mb),
                ("mem_per_call_mb", malloc_mb / len(g)),
            ]))
            index.append(locus)

        import pandas as pd
        df = pd.DataFrame(rows, index=index, columns=list(rows[0].keys()))
        return df.sort_values(by="ncalls", ascending=False)

    def find_zerosized(self, as_dataframe=False):
        """
        Find zero-sized allocations.

        Args:
            as_dataframe: True to return a pandas dataframe instead of a deque.
        """
        elist = []
        eapp = elist.append
        for e in self.all_entries:
            if e.size == 0: eapp(e)

        return entries_to_dataframe(elist) if as_dataframe else elist

    def find_weird_ptrs(self):
        """Find negative or zero pointers."""
        elist = []
        eapp = elist.append
        for e in self.all_entries:
            if e.ptr <= 0: eapp(e)

        if elist:
            print("Found %d weird entries:" % len(elist))
            pprint(elist)
        else:
            print("No weird entries found")
        return elist

    @lazy_property
    def all_entries(self):
        """Parse file and create list of Entries."""
        all_entries = []
        app = all_entries.append
        with open(self.path, "rt") as fh:
            for lineno, line in enumerate(fh):
                # skip header line of abimem files
                if line.startswith("#"): continue
                try:
                    app(Entry.from_line(line))
                except Exception as exc:
                    print("Error while parsing lineno %d, line:\n%s" % (lineno, line))
                    raise exc

        return all_entries

    def get_peaks(self, maxlen=30, as_dataframe=False):
        """
        Find peaks in the allocation with the corresponding variable.

        Args:
            maxlen: Maximum number of peaks
            as_dataframe: True to return a pandas dataframe instead of a deque.
        """
        # The deque is bounded to the specified maximum length. Once a bounded length deque is full,
        # when new items are added, a corresponding number of items are discarded from the opposite end.
        peaks = deque(maxlen=maxlen)

        visited = set()
        for e in self.all_entries:
            # Avoid redundant entries:
            #k = (e.locus, e.size)
            k = e.locus
            if e.size == 0 or not e.isalloc or k in visited: continue
            if len(peaks) == 0:
                peaks.append(e)
                visited.add(k)
                continue

            if e.size > peaks[0].size:
                peaks.append(e)
                visited.add(k)
                # Keep peaks sorted
                peaks = deque(sorted(peaks, key=lambda x: x.size), maxlen=maxlen)

        peaks = deque(sorted(peaks, key=lambda x: x.size, reverse=True), maxlen=maxlen)
        return entries_to_dataframe(peaks) if as_dataframe else peaks

    @lazy_property
    def dataframe(self):
        """
        Return a |pandas-DataFrame| with **all** entries.
        """
        return entries_to_dataframe(self.all_entries)

    @add_fig_kwargs
    def plot_memory_usage(self, ax=None, **kwargs):
        """
        Plot total allocated memory in Mb on axis `ax`.
        """
        memory = [e.tot_memory_mb for e in self.all_entries]
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ax.plot(memory)
        ax.grid(True)
        ax.set_ylabel("Total Memory (Mb)")
        return fig

    @add_fig_kwargs
    def plot_peaks(self, ax=None, maxlen=20, fontsize=4, rotation=25, **kwargs):
        """
        Plot memory peaks as vertical bars.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            maxlen: Maximum number of peaks
            fontsize: fontsize for legends and titles
            rotation: Rotation angle for xticklabels.

        Returns: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        peaks = self.get_peaks(maxlen=maxlen)
        data = [e.size_mb for e in peaks]
        names = ["%s\n%s" % (e.vname, e.locus) for e in peaks]
        xs = list(range(len(data)))
        ax.bar(xs, data)
        ax.grid(True)
        ax.set_xticks(xs)
        ax.set_xticklabels(names, fontsize=fontsize, rotation=rotation)
        ax.set_ylabel("Memory (Mb)")
        return fig

    @add_fig_kwargs
    def plot_hist(self, ax=None, **kwargs):
        """
        Plot histogram with the number of arrays allocated for a given size

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.

        Returns: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        data = [e.size_mb for e in self.all_entries]
        ax.hist(data) #, bins=n_bins)
        ax.grid(True)
        ax.set_ylabel("Number of arrays")
        ax.set_xlabel("Memory (Mb)")
        return fig

    def get_hotspots_dataframe(self):
        """
        Return DataFrame with total memory allocated per Fortran file.
        """
        index, rows = [], []
        for filename, g in self.dataframe.groupby(by="file"):
            malloc_mb = g[g["action"] == "A"].size_mb.sum()
            free_mb = g[g["action"] == "D"].size_mb.sum()
            nalloc = len(g["action"] == "A")
            nfree = len(g["action"] == "D")
            rows.append(OrderedDict([
                ("malloc_mb", malloc_mb),
                ("free_mb", free_mb),
                #("diff_mb", malloc_mb + free_mb),
                ("nalloc", nalloc),
                ("nfree", nfree),
                #("npall", nalloc - nfree),
            ]))
            index.append(filename)

        import pandas as pd
        df = pd.DataFrame(rows, index=index, columns=list(rows[0].keys()))
        return df.sort_values(by="malloc_mb", ascending=False)

    def expose(self, slide_mode=False, slide_timeout=None, **kwargs):
        """
        Shows a predefined list of matplotlib figures with minimal input from the user.
        """
        #from abipy.tools.plotting import MplExpose
        with MplExpose(slide_mode=slide_mode, slide_timeout=slide_mode, verbose=1) as e:
            e(self.plot_memory_usage(show=False))
            e(self.plot_peaks(show=False))
            e(self.plot_hist(show=False))

    def find_memleaks(self, verbose=0):
        """
        Try to find memory leaks using the address of the arrays and the action performed (allocation/free).
        """
        heap, stack = Heap(), Stack()
        reallocs = []

        for newe in self.all_entries:
            p = newe.ptr
            if newe.size == 0: continue
            # Store new entry in list if the ptr is not in d
            # else we check if there's an allocation that matches a previous allocation
            # (zero-sized arrays are not included)
            # else there's a possible memory leak or some undected problems.
            if p not in heap:
                if newe.isalloc:
                    heap[p] = [newe]
                # isfree found but ptr has not been allocated:
                else:
                    # Likely comes from a reallocation
                    reallocs.append(newe)

            else:
                if newe.isfree and len(heap[p]) == 1 and heap[p][0].size + newe.size == 0:
                    heap.pop(p)
                else:
                    # In principle this should never happen but there are exceptions:
                    #
                    # 1) The compiler may decide to put the allocatable on the stack
                    #    In this case the ptr reported by gfortran is 0.
                    #
                    # 2) The allocatable variable is "reallocated" by the compiler (F2003).
                    #    Example:
                    #
                    #    allocate(foo(2,1))           ! p0 = &foo
                    #    foo = reshape([0,0], [2,1])  ! p1 = &foo. Reallocation of the LHS.
                    #                                 ! Use foo(:) to avoid that
                    #    deallocate(foo)              ! p2 = &foo
                    #
                    #    In this case, p2 != p0
                    if verbose:
                        print("WARN:", newe.ptr, newe, "ptr already on the heap ", len(heap[p]), \
                              " sizes = ", heap[p][0].size, newe.size)
                    #print("HEAP:", heap[newe.ptr])

                    locus = newe.locus
                    if locus not in stack:
                        stack[locus] = [newe]
                    else:
                        #if newe.ptr != 0: print(newe)
                        stack_loc = stack[locus]
                        ifind = -1
                        for i, olde in enumerate(stack_loc):
                            if newe.frees_onstack(olde):
                                ifind = i
                                break

                        if ifind != -1:
                            stack_loc.pop(ifind)
                        #else:
                        #    print(newe)

                    #if p == 0:
                    #    stack[p] = newe
                    #else:
                    #    print("varname", newe.vname, "in heap with size ",newe.size)
                    #    for weirde in heap[p]:
                    #        print("\tweird entry:", weirde)
                    #    heap[p].append(newe)

        if False and heap:
            # Possible memory leaks.
            count = -1
            keyfunc = lambda e: abs(e.size)
            for a, entries in heap.items():
                count += 1
                entries = [e for e in entries if e.size != 0]

                entries = sorted(entries, key=keyfunc)
                #if any(int(e.size) != 0 for e in l):

                #msizes = []
                for key, group in groupby(entries, keyfunc):
                    group = list(group)
                    #print([e.name for e in g])
                    pos_size = [e for e in group if e.size >0]
                    neg_size = [e for e in group if e.size <0]
                    if len(pos_size) != len(neg_size):
                        print("key", key)
                        for e in group:
                            print(e)
                        #print(list(g))

                #for i, e in enumerate(entries):
                #    print("\t[%d]" % i, e)
                #print("Count=%d" % count, 60 * "=")

        if heap: heap.show()
        if stack: stack.show()
        if verbose and reallocs:
            print("Possible reallocations:")
            pprint(reallocs)

        return len(heap) + len(stack) + len(reallocs)

    def get_panel(self):
        """
        Build panel with widgets to interact with the memocc file either in a notebook or in panel app.
        """
        from .memprof_panel import MoccViewer
        return MoccViewer(self).get_panel()


class Heap(dict):

    def show(self):
        print("=== HEAP OF LEN %s ===" % len(self))
        if not self: return
        # for p, elist in self.items():
        pprint(self, indent=4)
        print("")

    def pop_alloc(self, entry):
        if not entry.isfree: return 0
        elist = self.get[entry.ptr]
        if elist is None: return 0
        for i, olde in elist:
            if entry.size + olde.size != 0:
                elist.pop(i)
                return 1
        return 0


class Stack(dict):

    def show(self):
        print("=== STACK OF LEN %s ===" % len(self))
        if not self: return
        pprint(self)
        print("")


# Copied  from abipy.tools.plotting
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

