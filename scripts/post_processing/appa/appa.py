#!/usr/bin/env python

"""
Abinit Post Process Application
author: Martin Alexandre
last edited: April 2013
"""

import sys

from PyQt4 import QtGui


def main():

    import gui.main_frame as mf  # GUI

    app = QtGui.QApplication(sys.argv)
    ex = mf.MainFrame()
    try:
        sys.exit(app.exec_())
    except SystemExit:
        import os
        import signal
        os.kill(os.getpid(), signal.SIGQUIT)


if __name__ == "__main__":
    try:
        main()
    except:
        print("To use APPA, you first need to compile a FORTRAN file:")
        print("   cd ~appa/fortran")
        print("f2py -c math.f90 --fcompiler='gfortran' -m math")
        print("f2py -c topology.f90 --fcompiler='gfortran' -m topology")
        print("f2py -c positions.f90 --fcompiler='gfortran' -m positions")
