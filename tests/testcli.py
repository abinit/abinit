#!/usr/bin/env python
from __future__ import print_function, division, unicode_literals
import argparse
from pymods.fldiff import Differ
from pymods.yaml_tools import is_available as has_yaml

if __name__ != '__main__':
    raise ImportError('testcli.py is not an importable module.')

if has_yaml:
    import pymods.explore_test as explore_test


class ArgParser(object):

    def __init__(self):
        parser = argparse.ArgumentParser(
            description='tool box for Abinit test'
        )
        parser.set_defaults(cmd='not a command')
        sub = parser.add_subparsers()

        # Diff
        diff_parser = sub.add_parser('diff', help='make a diff between two'
                                     ' output file like the test bot')
        self.mk_diff(diff_parser)

        # Shell
        shell_parser = sub.add_parser('shell', help='explore a YAML test'
                                      ' config and browse the documentation',
                                      description=explore_test.__doc__)
        self.mk_shell(shell_parser)

        # Run
        args = parser.parse_args()

        if args.cmd == 'not a command':
            parser.parse_args(['--help'])
        else:
            getattr(self, args.cmd)(args)

    def mk_diff(self, parser):
        '''
            Create command line argument parser for the diff subcommand
        '''
        parser.set_defaults(cmd='diff')

        parser.add_argument('ref_file', metavar='REF', help='File reference')
        parser.add_argument('test_file', metavar='TESTED',
                            help='File to be compared')
        parser.add_argument('-t', '--tolerance', metavar='TOL', type=float,
                            default=1.01e-10)
        parser.add_argument('--include', action='store_true')
        parser.add_argument('--includeP', action='store_true')
        parser.add_argument('--no_yaml', action='store_true')
        parser.add_argument('--no_fl', action='store_true')
        parser.add_argument('--yaml-conf', metavar='YAML_CONF', nargs='?')

    def diff(self, args):
        '''
            diff subcommand, manual access to pymods.fldiff.Differ
        '''
        opts = {
            'tolerance': args.tolerance,
            'ignore': False if args.include else True,
            'ignoreP': False if args.includeP else True,
            'use_yaml': False if args.no_yaml else True,
            'use_fl': False if args.no_fl else True,
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

    def mk_shell(self, parser):
        '''
            Create command line argument parser for the shell subcommand
        '''
        parser.set_defaults(cmd='shell')

        parser.add_argument('file', metavar='FILE', nargs='?', default=None,
                            help='YAML file defining a test.')
        parser.add_argument('--debug', '-d', action='store_true',
                            default=False, help='Enable debug command.')


    def shell(self, args):
        '''
            shell subcommand, explore test config files and browse constraints
            and parameters documentation.
        '''
        if not has_yaml:
            print('YAML or Numpy support is not available. That command cannot'
                  ' run.')
            exit(1)
        if args.debug:
            explorer = explore_test.DebugExplorer()
        else:
            explorer = explore_test.Explorer()

        if args.file:
            filename = args.file
            explorer.do_load(filename)

        explorer.cmdloop()
        exit()


ArgParser()
