from __future__ import print_function, division, unicode_literals
from .register_tag import BaseDictWrapper


class Issue(object):
    def __init__(self, conf, msg):
        self.path = conf.path if conf.path else ('top level',)
        self.state = conf.current_state
        self.message = msg

    def __repr__(self):
        spath = '.'.join(self.path)
        sstate = ', '.join('='.join(item) for item in self.state.items())
        return 'At {}({}): {}'.format(spath, sstate, self.message)


class Failure(Issue):
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


class Success(Issue):
    pass


class Tester(object):
    def __init__(self, reference_docs, tested_docs, config):
        self.ref = reference_docs
        self.tested = tested_docs
        self.conf = config
        self.failures = []
        self.success = []

    def check_this(self, name, ref, tested):
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
                    msg = ('Exception while checking {}:\n'
                           '{}: {}').format(cons.name, e.__class__.__name__,
                                            str(e))
                    self.failures.append(Failure(self.conf, msg))
                else:
                    if success:
                        msg = '{} ok'.format(cons.name)
                        self.success.append(Success(self.conf, msg))
                    else:
                        msg = '{} failed'.format(cons.name)
                        self.failures.append(Failure(self.conf, msg,
                                                     ref, tested))

            if isinstance(ref, BaseDictWrapper):  # have children
                for child in ref:
                    if child not in tested:
                        msg = '{} was not present'.format(child)
                        self.failures.append(Failure(self.conf, msg))
                    else:
                        self.check_this(child, ref[child], tested[child])

    def run(self):
        top_cons = self.conf.get_top_level_constraints()
        for cons in top_cons:
            if top_cons[cons].apply_to == 'this':
                ref = [doc['obj'] for doc in self.ref]
                tested = [doc['obj'] for doc in self.tested]
                success = top_cons[cons].check(ref, tested, self.conf)
                if success:
                    msg = '{} ok'.format(cons.name)
                    self.success.append(Success(self.conf, msg))
                else:
                    msg = '{} failed'.format(cons.name)
                    self.failures.append(Failure(self.conf, msg))

        if len(self.ref) != len(self.tested):
            msg = 'there is not the same number of documents in both side'
            self.failures.append(Failure(self.conf, msg))
        else:
            # FIXME Use a non linear matching of documents ?
            # ref_doc['iterators'] could be of some use here
            for ref_doc, tested_doc in zip(self.ref, self.tested):
                with self.conf.use_filter(ref_doc['iterators']):
                    self.check_this(ref_doc['obj']['label'], ref_doc['obj'],
                                    tested_doc['obj'])

        return self.failures, self.success
