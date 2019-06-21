from __future__ import print_function, division, unicode_literals
from .common import BaseDictWrapper, string, basestring


def short_repr(thing):
    '''
    Shorten representation of things when the default one is too long.
    '''
    s = str(thing)
    if len(s) > 30:
        if hasattr(thing, 'short_str'):
            return thing.short_str()
        else:
            return '<{} instance>'.format(type(thing).__name__)
    else:
        return s


class Issue(object):
    '''
        Represent the result of a test.
    '''
    def __init__(self, conf, msg):
        self.path = conf.path if conf.path else ('top level',)
        self.state = conf.current_state
        self.message = msg

    def __repr__(self):
        spath = '.'.join(self.path)
        sstate = ', '.join('{}={}'.format(*it) for it in self.state.items())
        return 'At {}({}): {}'.format(spath, sstate, self.message)

    def is_fail(self):
        return False


class Failure(Issue):
    '''
        Represent the fail of a test.
    '''
    def __init__(self, conf, msg, ref=None, tested=None):
        self.ref = ref
        self.tested = tested
        Issue.__init__(self, conf, msg)

    def __repr__(self):
        if self.ref is not None:
            return Issue.__repr__(self) + '\nref: {}\ntested: {}'.format(
                short_repr(self.ref), short_repr(self.tested)
            )
        else:
            return Issue.__repr__(self)

    def is_fail(self):
        return True


class DetailedFailure(Failure):
    def __init__(self, conf, msg, details):
        self.details = details
        Issue.__init__(self, conf, msg)

    def __repr__(self):
        return Issue.__repr__(self) + '\n{}'.format(self.details)


class Success(Issue):
    pass


class Tester(object):
    '''
        Drive the testing process.
    '''
    def __init__(self, reference_docs, tested_docs, config):
        self.ref = reference_docs
        self.tested = tested_docs
        self.conf = config
        self.issues = []

    def check_this(self, name, ref, tested):
        '''
            Check constraints applying to the 'tested' node of name 'name'
            against 'ref'
        '''

        def analyze(success, cons):
            if success:
                msg = '{} ok'.format(cons.name)
                self.issues.append(Success(self.conf, msg))
            elif hasattr(success, 'details'):
                msg = '{} ({}) failed'.format(cons.name, short_repr(cons.value))
                self.issues.append(DetailedFailure(self.conf, msg,
                                                   success.details))
            else:
                msg = '{} ({}) failed'.format(cons.name, short_repr(cons.value))
                self.issues.append(Failure(self.conf, msg,
                                           ref, tested))

        # we want to detect only dictionaries, not classes that inherit from it
        if type(ref) is dict:
            ref = BaseDictWrapper(ref)
        if type(tested) is dict:
            tested = BaseDictWrapper(tested)

        with self.conf.go_down(name):
            constraints = self.conf.get_constraints_for(ref)

            if self.conf.debug:
                for cons in constraints:
                    success = cons.check(ref, tested, self.conf)
                    analyze(success, cons)
            else:
                for cons in constraints:
                    try:
                        success = cons.check(ref, tested, self.conf)
                    except Exception as e:
                        msg = ('Exception while checking {} ({}/{}):\n'
                               '{}: {}').format(cons.name, short_repr(ref),
                                                short_repr(tested),
                                                type(e).__name__, str(e))
                        self.issues.append(Failure(self.conf, msg))
                    else:  # no exceptions
                        analyze(success, cons)

            if getattr(ref, 'is_dict_like', False):  # have children
                for child in ref:
                    if child not in tested:
                        msg = '{} was not present'.format(child)
                        self.issues.append(Failure(self.conf, msg))
                    else:
                        self.check_this(child, ref[child], tested[child])

            elif hasattr(ref, 'get_children'):  # user made browsable
                try:
                    dref = ref.get_children()
                    dtest = tested.get_children()
                except Exception as e:
                    msg = ('Tried to get a dict of item from {} but failed:\n'
                           '{}: {}').format(name, type(e).__name__, str(e))
                    self.issues.append(Failure(self.conf, msg))
                else:
                    for child in dref:
                        if child not in dtest:
                            msg = '{} was not present'.format(child)
                            self.issues.append(Failure(self.conf, msg))
                        else:
                            self.check_this(child, dref[child], dtest[child])

            elif (hasattr(ref, '__iter__')
                  and not getattr(ref, 'has_no_child', False)
                  and not isinstance(ref, basestring)):
                for index, (vref, vtest) in enumerate(zip(ref, tested)):
                    self.check_this(string(index), vref, vtest)

    def run(self):
        '''
            Main entry point for testing.
        '''
        top_cons = self.conf.get_top_level_constraints()
        for cons in top_cons:
            if top_cons[cons].apply_to('this'):
                # FIXME How to define and use top level constraints applying on
                # this like equations ?
                raise NotImplementedError('Top level constraints are not yet'
                                          ' implemented')
            # else:  if it does not apply on this it apply on a deeper level
            # so it is managed later

        for doc_id, ref_doc in self.ref.items():
            if doc_id not in self.tested:
                msg = ('Document ({}) is not present in the tested file.'
                       .format(doc_id))
                self.issues.append(Failure(self.conf, msg))
            with self.conf.use_filter(ref_doc.iterators):
                self.check_this(ref_doc.tag, ref_doc.obj,
                                self.tested[doc_id].obj)

        return self.issues
