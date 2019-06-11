from __future__ import print_function, division, unicode_literals
from .tricks import string
from .register_tag import BaseDictWrapper


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
        if self.ref:
            return Issue.__repr__(self) + '\nref: {}\ntested: {}'.format(
                self.ref, self.tested
            )
        else:
            return Issue.__repr__(self)

    def is_fail(self):
        return True


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
        self.issues = []

    def check_this(self, name, ref, tested):
        '''
            Check constraints applying to the 'tested' node (of name 'name')
            against 'ref'
        '''
        # we want to detect only dictionaries, not classes that inherit from it
        if type(ref) is dict:
            ref = BaseDictWrapper(ref)
        if type(tested) is dict:
            tested = BaseDictWrapper(tested)

        with self.conf.go_down(name):
            constraints = self.conf.get_constraints_for(ref)

            for cons in constraints:
                try:
                    success = cons.check(ref, tested, self.conf)
                except Exception as e:
                    msg = ('Exception while checking {} ({}/{}):\n'
                           '{}: {}').format(cons.name, ref, tested,
                                            e.__class__.__name__, str(e))
                    self.issues.append(Failure(self.conf, msg))
                else:
                    if success:
                        msg = '{} ok'.format(cons.name)
                        self.issues.append(Success(self.conf, msg))
                    else:
                        msg = '{} ({}) failed'.format(cons.name, cons.value)
                        self.issues.append(Failure(self.conf, msg,
                                                   ref, tested))

            if getattr(ref, '_is_dict_like', False):  # have children
                for child in ref:
                    if child not in tested:
                        msg = '{} was not present'.format(child)
                        self.issues.append(Failure(self.conf, msg))
                    else:
                        self.check_this(child, ref[child], tested[child])

            elif hasattr(ref, '__iter__') \
                    and not getattr(ref, '_has_no_child', False) \
                    and not isinstance(ref, string):
                for index, (vref, vtest) in enumerate(zip(ref, tested)):
                    self.check_this(index, vref, vtest)

    def run(self):
        '''
            Main entry point for testing.
        '''
        top_cons = self.conf.get_top_level_constraints()
        for cons in top_cons:
            if top_cons[cons].apply_to('this'):
                # ref = [doc.obj for doc in self.ref]
                # tested = [doc.obj for doc in self.tested]
                # success = top_cons[cons].check(ref, tested, self.conf)
                # if success:
                #     msg = '{} ok'.format(cons.name)
                #     self.issues.append(Success(self.conf, msg))
                # else:
                #     msg = '{} ({}) failed'.format(cons.name, cons.value)
                #     self.issues.append(Failure(self.conf, msg))
                # FIXME How to define and use top level constraints applying on
                # this like equations ?
                raise NotImplementedError('Top level constraints are not yet'
                                          ' implemented')

        if len(self.ref) != len(self.tested):
            msg = 'there is not the same number of documents in both side'
            self.issues.append(Failure(self.conf, msg))
        else:
            # FIXME Use a non linear matching of documents ?
            # ref_doc['iterators'] could be of some use here
            for ref_doc, tested_doc in zip(self.ref, self.tested):
                with self.conf.use_filter(ref_doc.iterators):
                    self.check_this(ref_doc.obj['label'], ref_doc.obj,
                                    tested_doc.obj)

        return self.issues
