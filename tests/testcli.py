#!/usr/bin/env python
from __future__ import print_function, division, unicode_literals
import argparse
import sys
from pymods.fldiff import Differ
from pymods.yaml_tools import is_available as has_yaml

if __name__ != '__main__':
    raise ImportError('testcli.py is not an importable module.')

if has_yaml:
    import pymods.explore_test as explore_test


class ArgParser(object):

    def __init__(self):
        parser = argparse.ArgumentParser(
            description='tool box for Abinit test',
            usage='''testtools <command> [<args>]

Available commands are
   diff      make a diff between two output file like the test bot
   shell     explore a YAML test config and browse the documentation
''')
        parser.add_argument('command', help='Subcommand to run')
        # parse_args defaults to [1:] for args, but you need to
        # exclude the rest of the args too, or validation will fail
        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            print('Unrecognized command')
            parser.print_help()
            exit(1)
        # use dispatch pattern to invoke method with same name
        getattr(self, args.command)()

    def diff(self):
        parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter
        )

        # Minimal command line interface for debugging

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

        args = parser.parse_args()

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
        except Exception as e:
            print('Something went wrong with this test:\n{}\n'.format(str(e)))

    def shell(self):
        if not has_yaml:
            print('YAML or Numpy support is not available. That command cannot'
                  ' run.')
            return
        parser = argparse.ArgumentParser(description=explore_test.__doc__)
        parser.add_argument('file', metavar='FILE', nargs='?', default=None,
                            help='YAML file defining a test.')
        parser.add_argument('--debug', '-d', action='store_true',
                            default=False, help='Enable debug command.')
        args = parser.parse_args()

        if args.debug:
            explorer = explore_test.DebugExplorer()
        else:
            explorer = explore_test.Explorer()

        if args.file:
            filename = args.file
            explorer.do_load(filename)

        explorer.cmdloop()


ArgParser()
