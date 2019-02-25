from __future__ import print_function, division, unicode_literals


class Issue(object):
    def __init__(self, path, msg):
        self.path = path
        self.message = msg

    def __repr__(self):
        spath = '.'.join(self.path)
        return 'At {}: {}'.format(spath, self.message)


class Failure(Issue):
    pass


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
        with self.conf.go_down(name):
            constraints = self.conf.get_constraints_for(ref)

            for cons in constraints:
                success = cons.check(ref, tested, self.conf)
                if success:
                    msg = '{} ok'.format(cons.name)
                    self.success.append(Success(self.conf.path, msg))
                else:
                    msg = '{} failed'.format(cons.name)
                    self.failures.append(Failure(self.conf.path, msg))

            if hasattr(ref, '__getitem__'):
                for child in ref:
                    if child not in ref:
                        msg = '{} was not present'.format(child)
                        self.failures.append(Failure(self.conf.path, msg))
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
                    self.success.append(Success(('top level',), msg))
                else:
                    msg = '{} failed'.format(cons.name)
                    self.failures.append(Failure(('top level',), msg))

        if len(self.ref) != len(self.tested):
            msg = 'there is not the same number of documents in both side'
            self.failures.append(Failure(('top level',), msg))
        else:
            # FIXME Use a non linear matching of documents ?
            # ref_doc['iterator'] could be of some use here
            for ref_doc, tested_doc in zip(self.ref, self.tested):
                if self.conf.has_doc(ref_doc['obj'].label):
                    self.check_this(ref_doc['obj'].label, ref_doc['obj'],
                                    tested_doc['obj'])
                else:
                    # There is not specialization for this document
                    pass

        return self.failures, self.success


class Result(object):
    '''
        Analyse results and create a report.
    '''
    def __init__(self, failures, success=None, conf=None):
        self.failures = failures
        self.success = success
        self.conf = conf

    def report(self):
        if self.failures:
            return '\n'.join(repr(fail) for fail in self.failures)
        else:
            return 'success'

    def passed(self):
        if self.failures:
            return False
        return True

    def dump_details(self, f):
        f.write(self.report())
