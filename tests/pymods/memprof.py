from __future__ import print_function, division, unicode_literals

from pprint import pprint
from itertools import groupby
from functools import wraps
from collections import namedtuple, deque
# OrderedDict was added in 2.7. ibm6 still uses python2.6
try:
    from collections import OrderedDict
except ImportError:
    from .ordereddict import OrderedDict


def group_entries_bylocus(entries):
    d = {}
    for e in entries:
        if e.locus not in d:
            d[e.locus] = [e]
        else:
            d[e.locus].append(e)
    return d


class Entry(namedtuple("Entry", "vname, ptr, action, size, file, func, line, tot_memory, sidx")):

    @classmethod
    def from_line(cls, line, sidx):
        args = line.split()
        args.append(sidx)
        return cls(*args)

    def __new__(cls, *args):
        """Extends the base class adding type conversion of arguments."""
        # write(logunt,'(a,t60,a,1x,2(i0,1x),2(a,1x),2(i0,1x))')&
        # trim(vname), trim(act), addr, isize, trim(basename(file)), trim(func), line, memtot_abi%memory
        return super(cls, Entry).__new__(cls,
        	vname=args[0],
        	action=args[1],
        	ptr=int(args[2]),
        	size=int(args[3]),
        	file=args[4],
        	func=args[5],
        	line=int(args[6]),
        	tot_memory=int(args[7]),
        	sidx=args[8],
        )

    def __repr__(self):
        return self.as_repr(with_addr=True)

    def as_repr(self, with_addr=True):
        if with_addr:
            return "<var=%s, %s@%s:%s:%s, addr=%s, size=%d, idx=%d>" % (
              self.vname, self.action, self.file, self.func, self.line, hex(self.ptr), self.size, self.sidx)
        else:
            return "<var=%s, %s@%s:%s:%s, size=%d, idx=%d>" %  (
              self.vname, self.action, self.file, self.func, self.line, self.size, self.sidx)

    @property
    def basename(self):
        return self.vname.split("%")[-1]

    @property
    def isalloc(self):
        """True if entry represents an allocation."""
        return self.action == "A"

    @property
    def isfree(self):
        """True if entry represents a deallocation."""
        return self.action == "D"

    @property
    def iszerosized(self):
        """True if this is a zero-sized alloc/free."""
        return self.size == 0

    @property
    def locus(self):
        """This is almost unique"""
        return self.func + "@" + self.file

    def frees_onheap(self, other):
        if (not self.isfree) or other.isalloc: return False
        if self.size + other.size != 0: return False
        return True

    def frees_onstack(self, other):
        if (not self.isfree) or other.isalloc: return False
        if self.size + other.size != 0: return False
        if self.locus != other.locus: return False
        return True


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
        print("=== STACK OF LEN %s ===)" % len(self))
        if not self: return
        pprint(self)
        print("")


def catchall(method):
    @wraps(method)
    def wrapper(*args, **kwargs):
        self = args[0]
        try:
            return method(*args, **kwargs)
        except Exception as exc:
            # Add info on file and re-raise.
            msg = "Exception while parsing file: %s\n" % self.path
            raise exc.__class__(msg + str(exc))
    return wrapper


class AbimemParser(object):
    def __init__(self, path):
        self.path = path

    #def __str__(self):
    #    lines = []
    #    app = lines.append
    #    return "\n".join(lines)

    @catchall
    def summarize(self):
        with open(self.path, "rt") as fh:
            l = fh.read()
            print(l)

    @catchall
    def find_small_allocs(self, nbytes=160):
        """Zero sized allocations are not counted."""
        smalles = []
        with open(self.path, "rt") as fh:
            for lineno, line in enumerate(fh):
                if lineno == 0: continue
                e = Entry.from_line(line, lineno)
                if not e.isalloc: continue
                if 0 < e.size <= nbytes: smalles.append(e)

        pprint(smalles)

        return smalles

    @catchall
    def find_intensive(self, threshold=2000):
        d = {}
        with open(self.path, "rt") as fh:
            for lineno, line in enumerate(fh):
                if lineno == 0: continue
                e = Entry.from_line(line, lineno)
                loc = e.locus
                if loc not in d:
                    d[loc] = [e]
                else:
                    d[loc].append(e)

        # Remove entries below the threshold and perform DSU sort
        dsu_list = [(elist, len(elist)) for _, elist in d.items() if len(elist) >= threshold]
        intensive = [t[0] for t in sorted(dsu_list, key=lambda x: x[1], reverse=True)]

        for elist in intensive:
            loc = elist[0].locus
            # assert all(e.locus == loc for e in elist)
            print("[%s] has %s allocations/frees" % (loc, len(elist)))

        return intensive

    #def show_peaks(self):

    @catchall
    def find_zerosized(self):
        elist = []
        eapp = elist.append
        for e in self.yield_all_entries():
            if e.size == 0: eapp(e)

        if elist:
            print("Found %d zero-sized entries:" % len(elist))
            pprint(elist)
        else:
            print("No zero-sized found")
        return elist

    @catchall
    def find_weird_ptrs(self):
        elist = []
        eapp = elist.append
        for e in self.yield_all_entries():
            if e.ptr <= 0: eapp(e)

        if elist:
            print("Found %d weird entries:" % len(elist))
            pprint(elist)
        else:
            print("No weird entries found")
        return elist

    def yield_all_entries(self):
        with open(self.path, "rt") as fh:
            for lineno, line in enumerate(fh):
                if lineno == 0: continue # skip header line of abimem files
                yield Entry.from_line(line, lineno)

    @catchall
    def find_peaks(self, maxlen=20):
        # the deque is bounded to the specified maximum length. Once a bounded length deque is full,
        # when new items are added, a corresponding number of items are discarded from the opposite end.
        peaks = deque(maxlen=maxlen)

        for e in self.yield_all_entries():
            size = e.size
            if size == 0 or not e.isalloc: continue

            if len(peaks) == 0:
                peaks.append(e); continue

            # TODO: Should remove redundant entries.
            if size > peaks[0].size:
                peaks.append(e)
                peaks = deque(sorted(peaks, key=lambda x: x.size), maxlen=maxlen)

        peaks = deque(sorted(peaks, key=lambda x: x.size, reverse=True), maxlen=maxlen)

        for peak in peaks:
            print(peak)

        return peaks

    @catchall
    def plot_memory_usage(self, show=True):
        memory = [e.tot_memory for e in self.yield_all_entries()]
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(memory)
        if show: plt.show()
        return fig

    #def get_dataframe(self):
    #    import pandas as pd
    #    frame = pd.DataFrame()
    #    return frame

    @catchall
    def find_memleaks(self):
        heap, stack = Heap(), Stack()
        reallocs = []

        with open(self.path, "rt") as fh:
            for lineno, line in enumerate(fh):
                if lineno == 0: continue
                newe = Entry.from_line(line, lineno)

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
                        # 1) The compiler could decide to put the allocatable on the stack
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
        if reallocs:
            print("Possible reallocations:")
            pprint(reallocs)

        return len(heap) + len(stack) + len(reallocs)
