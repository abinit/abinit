'''\
This is the explore_test shell. This tool let you inspect and explore YAML
files defining a test for Abinit. It also provide documentation about the
constraints and parameters available in test config files.
'''

from __future__ import print_function, division, unicode_literals
import os
import glob
import cmd
from subprocess import call
from .tester_conf import TesterConf
from .conf_parser import conf_parser
from .errors import ConfigError

intro = '''\
Welcome to the explore_test shell.
This tool let you inspect and explore a YAML file defining a test for Abinit.
You can also browse informations about parameters and constraints used to
define tests with the command show.

If your Python platform has readline support you can use TAB to auto-complete
command and some arguments.

Type help or ? to get the list of commands.
'''


try:
    import readline
    # remove / and ~ from delimiters to avoid problem with file path completion
    # readline.set_completer_delims(' \t\n`!@#$%^&*()-=+[{]}\\|;:\'",<>?')
    readline.set_completer_delims(' \n\t')
except ImportError:
    pass


def show_cons(cons, used=False):
    print('Constraint', cons.name)
    print(' Value type:', cons.type)
    print(' Inherited:', 'yes' if cons.inherited else 'no')
    print(' Parameters:')
    for param in cons.use_params:
        print(' -', param)
    print(' Exclude:')
    for name in cons.exclude:
        print(' -', name)
    print(' Current value:', cons.value if used else 'not defined here')
    print(' Description:')
    print(cons.test.__doc__)


def show_param(name, dic, used=False):
    print('Parameter', name)
    print(' Value type:', dic['type'])
    print(' Inherited:', 'yes' if dic['inherited'] else 'no')
    print(' Default value:', dic['default'])
    print(' Current value:', dic['value'] if used else 'default value')


def short_path(path):
    if len('.'.join(path)) > 50:
        shpath = path[:2] + '...' + '.'.join(path[-2:])
    else:
        shpath = '.'.join(path)
    return shpath


def print_iter(it):
    i = 0
    for elem in it:
        i += 1
        print(elem, end=', ')
        if i % 8 == 0:
            print()
    if i % 8 != 0:
        print()


class ExtendedTesterConf(TesterConf):
    def get_spec(self):
        '''
            Return the list of the specializations known at the current path.
        '''
        return self.tree.get_spec_at(self.current_path).keys()

    def get_spec_at(self, path):
        '''
            Return the list of the specializations known at the current path.
        '''
        new_path = list(self.current_path)
        for sp in path:
            if sp == 'TOP':
                new_path = []
            elif sp == 'UP':
                new_path.pop()
            else:
                new_path.append(sp)
        return self.tree.get_spec_at(new_path).keys()

    def go_root(self):
        '''
            Go back to root, flushing all stacks and reset path.
        '''
        self.current_path = []
        self.param_stack = []
        self.constraints_stack = []

    def get_all_constraints_here(self):
        '''
            return a dict of the constraints in the current scope.
        '''
        constraints = {}

        cursor = 0
        top = len(self.param_stack) - 1

        # traverse from top to bottom to overwrite if needed
        while cursor <= top:
            for name, cons in self.constraints_stack[cursor].items():
                if cursor < top and not cons.inherited:
                    # pass if we are no longer on
                    # the same level as the caller
                    continue
                else:
                    for name in cons.exclude:
                        if name in constraints:
                            del constraints[name]
                    constraints[cons.name] = cons
            cursor += 1

        return constraints

    def get_all_parameters_here(self):
        '''
            return a dict of the parameters in the current scope.
        '''
        parameters = {}

        cursor = 0
        top = len(self.param_stack) - 1

        # traverse from top to bottom to overwrite if needed
        while cursor <= top:
            for name, value in self.param_stack[cursor].items():
                if cursor < top and not self.known_params[name]['inherited']:
                    # pass if we are no longer on
                    # the same level as the caller
                    continue
                else:
                    parameters[name] = self.known_params[name].copy()
                    parameters[name]['value'] = value
            cursor += 1

        return parameters

    def get_known_constraints(self):
        return conf_parser.constraints.copy()

    def get_known_parameters(self):
        return self.known_params.copy()

    @property
    def path(self):
        return tuple(['TOP'] + self.current_path)

    def dump_state(self):
        return {
            'cons': self.constraints_stack.copy(),
            'param': self.param_stack.copy(),
            'path': self.current_path,
        }

    def restore_state(self, state):
        self.constraints_stack = state['con']
        self.param_stack = state['param']
        self.current_path = state['path']


class Explorer(cmd.Cmd):
    intro = intro

    debug = False
    prompt = '() '
    tree = None
    filename = ''
    full_path = ''

    def update_prompt(self):
        template = '{filename}: {path}> '
        if self.tree is None:
            self.prompt = '[no file loaded] > '
        else:
            self.prompt = template.format(
                filename=self.filename,
                path=short_path(self.tree.path)
            )

    # special hooks
    def emptyline(self):
        pass

    def preloop(self):
        self.update_prompt()

    def postcmd(self, stop, line):
        self.update_prompt()
        return stop

    def precmd(self, line):
        # most command should not be called if tree is not loaded
        if line == 'EOF':
            line = 'exit'
        if not self.tree and line and line.split(' ')[0] not in [
            'load', 'edit',
            'shell', '!', 'debug',
            'exit',
            'help', '?',
        ]:
            print('That command can only be used with an opened file.')
            return ''
        else:
            return line

    # commands
    def do_load(self, arg):
        '''
            Usage: load FILE
            Load a config file.
        '''
        filename = os.path.realpath(os.path.expanduser(arg))
        try:
            self.tree = ExtendedTesterConf.from_file(filename)
        except IOError:
            print('File not found.')
        except ConfigError as e:
            print('Invalid yaml config:')
            print(e)
        else:
            self.filename = os.path.basename(filename)
            self.full_path = filename
            print(arg, 'successfully loaded.')

    def do_up(self, arg):
        '''
            Usage: up
            Go up of one level.
        '''
        self.tree.go_up()

    def do_cd(self, arg):
        '''
            Usage: cd PATH
            Move to PATH relative to the current path the tree.  PATH is of the
            form name1.name2... name can be either a specialization or TOP to
            go back to root or UP to go up one level.
        '''
        if not arg.isspace():
            path = arg.replace('"', '').replace('\'', '').split('.')
            for spec in path:
                if spec == 'TOP':
                    self.tree.go_root()
                elif spec == 'UP':
                    self.tree.go_up()
                elif spec:
                    self.tree.go_down(spec)

    def complete_file_path(self, text):
        '''
            Autocompletion for file path arguments.
        '''
        return glob.glob(os.path.expanduser(text)+'*')

    def complete_load(self, text, line, begi, endi):
        return self.complete_file_path(text)

    def complete_edit(self, text, line, begi, endi):
        return self.complete_file_path(text)

    def complete_rel_path(self, text):
        '''
            Autocompletion for path arguments.
        '''

        specs = ['TOP', 'UP'] + list(self.tree.get_spec())
        path = text.replace('"', '').replace('\'', '').split('.')

        if not path:
            return specs
        elif len(path) > 1:
            partial = path[-1]
            path = path[:-1]
            return ['.'.join(path) + '.' + spec
                    for spec in self.tree.get_spec_at(path)
                    if spec.startswith(partial)]
        elif len(path) == 1:
            partial = path[0]
            return [spec for spec in specs if spec.startswith(partial)]

    def complete_ls(self, text, line, begi, endi):
        return self.complete_rel_path(text)

    def complete_cd(self, text, line, begi, endi):
        return self.complete_rel_path(text)

    def do_path(self, arg):
        '''
            Usage: path
            Print the current path.
        '''
        spath = '.'.join(self.tree.path)
        if spath:
            print(spath)
        else:
            print('top level')

    def do_ls(self, arg):
        '''
            Usage: ls PATH
            List nodes under the given PATH. See also cd.
        '''
        if not arg:
            print_iter(spec for spec in self.tree.get_spec())
        else:
            path = arg.replace('"', '').replace('\'', '').split('.')
            print_iter(spec for spec in self.tree.get_spec_at(path))

    def do_show(self, arg):
        '''
            Usage: show [ARG | *]
            If no argument is given, list all parameters and constraints
            visible from the current level.  If argument is *, list all
            parameters and constraints known.  If argument is ARG, show all
            informations about ARG.
        '''
        name = arg
        cons_scope = self.tree.get_all_constraints_here()
        param_scope = self.tree.get_all_parameters_here()
        cons_known = self.tree.get_known_constraints()
        param_known = self.tree.get_known_parameters()
        if not name:
            print('Current scope')
            print('Constraints:')
            print_iter(cons_scope)
            print('\nParameters:')
            print_iter(param_scope)
            print()
        elif name == '*':
            print('Everything known')
            print('Constraints:')
            print_iter(cons_known)
            print('\nParameters:')
            print_iter(param_known)
            print()
        else:
            if name in cons_scope:
                show_cons(cons_scope[name], used=True)
            elif name in cons_known:
                show_cons(cons_known[name], used=False)
            elif name in param_scope:
                show_param(name, param_scope[name], used=True)
            elif name in param_known:
                show_param(name, param_known[name], used=False)
            else:
                print(arg, 'is neither a known parameter,'
                      ' nor a known constraint')

    def complete_show(self, text, line, begi, endi):
        return [
            name for name in self.tree.get_known_constraints()
            if name.startswith(text)
        ] + [
            name for name in self.tree.get_known_parameters()
            if name.startswith(text)
        ]

    def do_parameters(self, arg):
        '''
            Usage: parameters
            List constraints applying at the current level.
        '''
        print_iter(self.tree.get_all_parameters_here())

    def do_constraints(self, arg):
        '''
            Usage: constraints
            List constraints applying at the current level.
        '''
        print_iter(self.tree.get_all_constraints_here())

    def do_tree(self, arg):
        '''
            Usage: tree
            Show the tree defined by the configuration.
        '''
        def show_rec(specs, indent=[]):
            for i, sp in enumerate(specs):
                print(''.join('  ' if last else '| ' for last in indent)
                      + '`--', sp, sep='')
                with self.tree.go_down(sp):
                    nspecs = self.tree.get_spec()
                    show_rec(nspecs, indent + [i + 1 == len(specs)])

        toplvl = self.tree.get_spec()
        print('.'.join(self.tree.path))
        show_rec(toplvl, [True])

    def do_shell(self, arg):
        '''
            Usage: shell CMD ARG1 ARG2...
            Pass command to the system shell.
        '''
        sh = os.environ.get('SHELL', '/bin/sh')
        try:
            call([sh, '-c', arg])
        except IOError:
            print('The shell command {} cannot be found.'.format(sh),
                  'You may want to set your SHELL envrionment variable to',
                  'select a different command.')

    def do_edit(self, arg):
        '''
            Usage: edit [FILE]
            Open a file in an editor and load it once done. If no path is
            provided the current file is used.  If it doesn't work try running
            the following in a shell:
            $ export EDITOR='nano'
            Replace nano with any editor you like. Then restart
            testconf_explorer.
        '''
        ed = os.environ.get('EDITOR', 'nano')
        if arg:
            filepath = arg
            path = ''
        else:
            if not self.tree:
                print('No file have been loaded yet.',
                      'Provide a file name to edit.')
                return
            filepath = self.full_path
            path = '.'.join(self.tree.path)

        try:
            call([ed, filepath])
        except OSError:
            print('The editor command {} cannot be found.'.format(ed),
                  'You may want to set your EDITOR envrionment variable to',
                  'select a different command.')
        else:
            self.do_load(filepath)
            self.do_cd(path)

    def do_cat(self, arg):
        '''
            Usage: cat
            Print the current config file verbatim.
        '''
        with open(self.full_path) as f:
            print(f.read())

    def do_exit(self, arg):
        '''
            Usage: exit
            Exit. CTRL-D is equivalent.
        '''
        print('Bye.')
        return True


class DebugExplorer(Explorer):
    def do_debug(self, arg):
        '''
            Usage debug PYTHON_EXPR
            Print the result of any arbitrary python expression.
            self give access to the main Explorer instance.
        '''
        print(eval(arg))
