"""
Main home of the test engine. Define the class used to traverse the data tree
while checking constraints as well as a few utility classes linked to that process.
"""
from .common import BaseDictWrapper, basestring, string


def short_repr(thing):
    """
    Shorten representation of things when the default one is too long.
    """
    s = str(thing)
    if len(s) > 30:
        if hasattr(thing, "short_str"):
            return thing.short_str()
        return f"<{type(thing).__name__} instance>"
    return s


class Issue:
    """
    Represent the result of a test.
    """
    def __init__(self, conf, msg):
        """
        Initialize the Issue object.

        Args:
            conf (DriverTestConf): The configuration driver.
            msg (str): The issue message.
        """
        self.path = conf.path or ("top level",)
        self.state = conf.current_state
        self.message = msg

    def __repr__(self):
        spath = ".".join(self.path)
        sstate = ", ".join("{}={}".format(*it) for it in self.state.items())
        return f"At {spath}({sstate}): {self.message}"

    def is_fail(self):
        return False


class Failure(Issue):
    """
    Represent the fail of a test.
    """
    def __init__(self, conf, msg, ref=None, tested=None):
        """
        Initialize the Failure object.

        Args:
            conf (DriverTestConf): The configuration driver.
            msg (str): The failure message.
            ref (optional): The reference value.
            tested (optional): The tested value.
        """
        self.ref = ref
        self.tested = tested
        Issue.__init__(self, conf, msg)

    def __repr__(self):
        if self.ref is not None:
            return Issue.__repr__(self) + f"\nref: {short_repr(self.ref)}\ntested: {short_repr(self.tested)}"
        return Issue.__repr__(self)

    def is_fail(self):
        return True


class DetailedFailure(Failure):
    """
    Represent the fail of a test when more info are available.
    """
    def __init__(self, conf, msg, details):
        """
        Initialize the DetailedFailure object.

        Args:
            conf (DriverTestConf): The configuration driver.
            msg (str): The failure message.
            details: Additional details about the failure.
        """
        self.details = details
        Issue.__init__(self, conf, msg)

    def __repr__(self):
        return Issue.__repr__(self) + f"\n{self.details}"


class Success(Issue):
    """Represent the success of a test."""


class Tester:
    """
    Drive the testing process.
    """
    def __init__(self, reference_docs, tested_docs, config):
        """
        Initialize the Tester object.

        Args:
            reference_docs (dict): Map of document IDs to reference Document objects.
            tested_docs (dict): Map of document IDs to tested Document objects.
            config (DriverTestConf): Configuration driver.
        """
        self.ref = reference_docs
        self.tested = tested_docs
        self.conf = config
        self.issues = []

    def check_this(self, name, ref, tested):
        """
        Check constraints applying to the 'tested' node of name 'name'
        against 'ref' and recursively apply on children.
        """
        #print("name:", name, "\nref:", ref, "\ntested:",tested)

        if ref is None or tested is None:
            if ref is not tested:
                msg = f"Expecting two None objects but got {ref!r} and {tested!r}"
                self.issues.append(Failure(self.conf, msg, ref, tested))
            else:
                #print("Got None None for name:", name)
                return

        def analyze(success, cons):
            """Analyse the result of a constraint check."""
            if success:
                msg = f"{cons.name} ok"
                self.issues.append(Success(self.conf, msg))

            elif hasattr(success, "details"):
                msg = f"{cons.name} ({short_repr(cons.value)}) failed"
                self.issues.append(DetailedFailure(self.conf, msg, success.details))

            else:
                msg = f"{cons.name} ({short_repr(cons.value)}) failed"
                self.issues.append(Failure(self.conf, msg, ref, tested))

        # we want to detect only dictionaries, not classes that inherit from it
        if type(ref) is dict:
            ref = BaseDictWrapper(ref)
        if type(tested) is dict:
            tested = BaseDictWrapper(tested)

        with self.conf.go_down(name):
            constraints = self.conf.get_constraints_for(ref)

            if self.conf.debug:
                # on debug mode, don't catch errors to keep the backtrace
                for cons in constraints:
                    success = cons.check(ref, tested, self.conf)
                    analyze(success, cons)
            else:
                for cons in constraints:
                    try:
                        success = cons.check(ref, tested, self.conf)
                    except Exception as e:
                        msg = (f"Exception while checking {cons.name} ({short_repr(ref)}/{short_repr(tested)}):\n"
                               f"{type(e).__name__}: {e!s}")
                        self.issues.append(Failure(self.conf, msg))
                    else:
                        # no exceptions
                        analyze(success, cons)

            if getattr(ref, "is_dict_like", False):
                # have children
                for child in ref:
                    if child not in tested:
                        msg = f"{child} was not present"
                        #raise RuntimeError(msg)
                        self.issues.append(Failure(self.conf, msg))
                    else:
                        self.check_this(child, ref[child], tested[child])

            elif hasattr(ref, "get_children"):
                # user made browsable
                try:
                    dref = ref.get_children()
                    dtest = tested.get_children()
                except Exception as e:
                    msg = (f"Tried to get a dict of item from {name} but failed:\n"
                           f"{type(e).__name__}: {e!s}")
                    self.issues.append(Failure(self.conf, msg))
                else:
                    for child in dref:
                        if child not in dtest:
                            msg = f"{child} was not present"
                            self.issues.append(Failure(self.conf, msg))
                        else:
                            self.check_this(child, dref[child], dtest[child])

            elif (hasattr(ref, "__iter__")
                  and not getattr(ref, "has_no_child", False)
                  and not isinstance(ref, basestring)):
                # strings have __iter__ but should not be browsed
                # user may want to neutralize the browsing of an iterable
                # using has_no_child (for example BaseArray)
                for index, (vref, vtest) in enumerate(zip(ref, tested)):
                    self.check_this(string(index), vref, vtest)

    def run(self):
        """
        Main entry point for testing.
        """
        top_cons = self.conf.get_top_level_constraints()
        for cons in top_cons:
            if top_cons[cons].apply_to("this"):
                # FIXME How to define and use top level constraints applying on
                # this like equations ? Is it worth it ?
                raise NotImplementedError("Top level constraints are not yet implemented")

        # browse document from the reference file
        for doc_id, ref_doc in self.ref.items():
            if doc_id not in self.tested:
                msg = f"Document ({doc_id}) is not present in the tested file."
                self.issues.append(Failure(self.conf, msg))

            # enter filters and start checking the docuement
            with self.conf.use_filter(ref_doc.iterators):
                self.check_this(ref_doc.tag, ref_doc.obj, self.tested[doc_id].obj)

        return self.issues
