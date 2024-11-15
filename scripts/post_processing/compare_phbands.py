#!/usr/bin/env python3
"""
Script to plot the phonon band structure and compare with the LWF phonon band structure.
 COPYRIGHT
 Copyright (C) 1999-2024 ABINIT group (HeXu)
 This file is distributed under the terms of the
 GNU General Public Licence, see ~abinit/COPYING
 or http://www.gnu.org/copyleft/gpl.txt .
 For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
"""

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt

Ha_cmm1 = 219474.6313705


def plot_file(fname, ax, **kwargs):
    d = np.loadtxt(fname)
    nband = d.shape[1]-1
    for i in range(1, nband+1):
        ax.plot(d[:, 0], d[:, i]*Ha_cmm1, **kwargs)
    ax.set_xlim(d[0, 0], d[-1, 0])
    return ax


def plot_compare(prefix, ax):
    fname1 = prefix+"_PHFRQ"
    fname2 = prefix+"_lwf_PHFRQ"
    if os.path.exists(fname1):
        plot_file(fname1, ax=ax, color='blue', alpha=0.9)
    else:
        raise IOError("Cannot find the file "+fname1+"Please check.")
    if os.path.exists(fname2):
        plot_file(fname2, ax=ax, marker='o', color='green', alpha=0.7)
    ax.axhline(linestyle='--', color='gray')
    ax.set_ylabel("Frequency ($cm^{-1}$)")
    return ax


def main():
    parser = argparse.ArgumentParser(
        description="This is a script to plot the phonon band structure and compare to the lattice Wannier function phonon band structure by comparing the PHFRQ files. If there are two files {prefix}_lwf_PHFRQ, and {prefix}_PHFRQ, the band structures will be compared. If only the  {prefix}_PHFRQ exists, it will plot the full phonon band structure.")
    parser.add_argument("prefix", type=str,
                        help="The prefix of the PHFRQ files, e.g. with the prefix 'run', the run_PHFRQ file contains the phonon bands. ", default="fun")
    parser.add_argument("--output", "-o", type=str,
                        help="The name of the output figure file, should end with usual image file names, like .png, .pdf, etc. If not specified, no file will be saved", default=None)
    parser.add_argument("--show", "-s",
                        action="store_true",
                        help="whether to show the band structure to screen.",
                        default=True)
    args = parser.parse_args()
    _fig, ax = plt.subplots()
    plot_compare(prefix=args.prefix, ax=ax)
    if args.output is not None:
        plt.savefig(args.output)
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
