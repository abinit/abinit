#!/usr/bin/env python
from __future__ import print_function, division, unicode_literals
import argparse
from os import path
from sys import version_info
from pymods.fldiff import Differ
from pymods.yaml_tools import is_available as has_yaml

if __name__ != '__main__':
    raise ImportError('testcli.py is not an importable module.')

PY3 = version_info.major >= 3

if has_yaml:
    import pymods.yaml_tools.explore_test as explore_test


class ArgParser(object):

    cmds = {
        'fldiff': {
            'help': 'make a diff between two output file like the test bot',
            'aliases': ['diff'],
        },
        'explore': {
            'help': 'explore a YAML test config and browse the documentation',
            'description': explore_test.__doc__,
            'aliases': ['exp', 'dig', 'sh'],
        }
    }

    def __init__(self):
        parser = argparse.ArgumentParser(
            description='tool box for Abinit test'
        )
        sub = parser.add_subparsers(dest='cmd')
        parser.set_defaults(cmd='not a command')

        for cmd in self.cmds:
            if not PY3:  # unfortunatly aliases are only available in Python 3
                del self.cmds[cmd]['aliases']
            cmd_parser = sub.add_parser(cmd, **self.cmds[cmd])
            getattr(self, 'parse_'+cmd)(cmd_parser)

        # Run
        args = parser.parse_args()

        self.unalias(args)

        if args.cmd == 'not a command':
            parser.parse_args(['--help'])
        else:
            getattr(self, args.cmd)(args)

    def alias(self, cmd):
        return self.aliases.get(cmd, {'aliases': []})['aliases']

    def unalias(self, args):
        if args.cmd in self.cmds:
            return
        for cmd, opts in self.cmds.items():
            if args.cmd in opts['aliases']:
                args.cmd = cmd
                return

    def parse_fldiff(self, parser):
        '''
            Create command line argument parser for the diff subcommand
        '''
        parser.add_argument('ref_file', metavar='REF', help='File reference')
        parser.add_argument('test_file', metavar='TESTED', help='File to be'
                            ' compared')
        parser.add_argument('-t', '--tolerance', metavar='TOL', type=float,
                            default=1.01e-10, help='Tolerance used to detect'
                            ' differences')
        parser.add_argument('-c', '--yaml-conf', metavar='YAML_CONF',
                            help='Provide a YAML file to configure the'
                            ' comparison')
        parser.add_argument('--include', help='Use "," meta character as "+"'
                            ' (default is "-")', action='store_true')
        parser.add_argument('--includeP', help='Use "P" meta character as "+"'
                            ' (default is "-")', action='store_true')
        parser.add_argument('--no-yaml', help='Disable YAML based comparison',
                            action='store_true')
        parser.add_argument('--no-fl', help='Disable legacy fldiff comparison',
                            action='store_true')
        parser.add_argument('--debug', '-d', action='store_true',
                            help='Enable debugging mode, no longer catch'
                            ' errors so you can read the tracebacks')

    def fldiff(self, args):
        '''
            diff subcommand, manual access to pymods.fldiff.Differ
        '''
        opts = {
            'tolerance': args.tolerance,
            'ignore': not args.include,
            'ignoreP': not args.includeP,
            'use_yaml': not args.no_yaml,
            'use_fl': not args.no_fl,
            'debug': args.debug,
        }

        yaml_test = {
            'file': args.yaml_conf
        }

        ref_name, test_name = args.ref_file, args.test_file

        differ = Differ(yaml_test, **opts)
        try:
            fld_result = differ.diff(ref_name, test_name)
            print(fld_result.dump_details())
            isok, _, _ = fld_result.passed_within_tols(0, 0.0, 0.0)
            ecode = 0 if isok else 1
            exit(ecode)
        except Exception as e:
            print('Something went wrong with this test:\n{}\n'.format(str(e)))
            exit(1)

    def parse_explore(self, parser):
        '''
            Create command line argument parser for the explore subcommand
        '''
        parser.add_argument('file', metavar='FILE', nargs='?', default=None,
                            help='YAML file defining a test.')
        parser.add_argument('-D', '--debug', action='store_true',
                            help='Enable debug command.')
        parser.add_argument('-d', '--default', metavar='DEFAULT_CONF', nargs=1,
                            help='YAML alternative default config for debug.',
                            default=None)

    def explore(self, args):
        '''
            explore subcommand, explore test config files and browse
            constraints and parameters documentation.
        '''
        if not has_yaml:
            print('YAML or Numpy support is not available. That command cannot'
                  ' run.')
            exit(1)
        if args.debug:
            if args.default:
                explore_test.ExtendedTestConf.default_conf = \
                        path.expanduser(args.default[0])
            explorer = explore_test.DebugExplorer()
        else:
            explorer = explore_test.Explorer()

        if args.file:
            filename = args.file
            explorer.do_load(filename)

        explorer.cmdloop()
        exit()


ArgParser()
