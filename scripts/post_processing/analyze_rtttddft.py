#!/usr/bin/env python3
"""
Script to compute the optical conductivity and the dielectric
function from the time-dependent current density obtained
after a RT-TDDFT run with ABINIT
----------------
Partially based on the script in the tutorials of the exciting code
----------------
F. Brieuc
"""

import argparse
import sys

import numpy as np
from scipy.fftpack import fft, fftfreq
from scipy.integrate import trapz

# ------------------ Constants ----------------------
hb = 1.054571817e-34   # J.s
me = 9.1093837139e-31  # kg
a0 = 5.29177210544e-11 # m
qe = 1.602176634e-19   # C
Eh = hb**2/(me*a0**2)  # J

# ------------------ Conversions ----------------------
ev2j = 1.602176634e-19 # J -> eV
au2s = hb/Eh           # au -> second
au2as = au2s * 1e17     # au -> attosecond
au2fs = au2s * 1e15     # au -> femtosecond
au2Vpm = Eh/(qe*a0)      # au -> Volt/meter
au2ev = Eh / ev2j       # au -> eV
au2Ohm = me*a0**2/(qe**2*au2s) # au -> Ohm
au2Ohmcm1 = 1/(au2Ohm*a0*1e2)  # au -> 1/(Ohm*cm)

# ------------------- Functions -----------------------


def fourier_direct(time, signal, wcut, nfft):
    """
    Direct Fourier transform - from time to angular frequency
    Input:
     - time: time (array)
     - signal: signal (array) f(t) to be transformed
     - wcut: cut-off frequency (in Ha) for the exponential window (if > 0)
     - nfft: number of points to use in the FFT
             if greater than the length of signal then zero-padding is used
    Output:
      - w, F: Tuple containing angular frequencies (in Ha)
              and the fourier transform of signal*filter
    ------
    If wcut > 0 then the signal is multiplied in time by exp(-wcut*t)
    which leads to a convolution with a Lorentzian function in frequency space.
    Indeed the Fourier transform of exp(-wcut*t) is w_cut/(w_cut^2+w^2).
    Thus this window function generates a Lorentzian broadening in frequency
    with a FWHM given by 2*w_cut.
    If wcut == -1 use a third order polynomial damping function as filter.
    From Yabana et al., Phys. Stat. Sol. (B) 243,1121 (2006)
    """
    ntime = len(time)

    if len(signal) != ntime:
        msg = "! Error (in fourier_direct)"
        msg += " signal and time do not have the same length"
        sys.exit(msg)

    dt = time[1] - time[0]

    ## Apply filter if required

    # Exponential window: exp(-wcut*t)
    if wcut > 0.0:
        filter = np.exp(np.multiply(time,-wcut))
        f = np.multiply(filter,signal)

    # Damping function: 1-3(t/tmax)^2+2(t/tmax)^3
    elif wcut == -1:
        tmax = time[-1]
        filter = np.ones(len(time)) \
                 - np.multiply(time**2,3.0/tmax**2) \
                 + np.multiply(time**3,2.0/tmax**3)
        f = np.multiply(filter,signal)
    else:
        f = signal

    ## Zero-padding if required
    dn = nfft - ntime
    if dn > 0:
        ff = np.zeros(nfft)
        ff[0:ntime] = f
    else:
        ff = f

    ## Fourier transform
    # Define angular frequencies
    w = 2*np.pi*fftfreq(nfft,dt)
    ft = np.conj(np.multiply(fft(ff),dt))
    # Sort to get the right order in frequency
    w,ft = zip(*sorted(zip(w,ft)))

    return np.array(w), np.array(ft)


def search_keyword(lines, keyword):
    """
    Search for a specific keyword in list of strings (lines)
    """
    for line in lines:
        if keyword in line:
            splitted_line = line.split()
            for il, ls in enumerate(splitted_line):
                if keyword in ls:
                    return splitted_line[il+1]


def parse_output(file):
    """
    Extract required infos from abinit output file
    """
    # read file
    with open(file) as f:
        lines = f.readlines()

    tmp = search_keyword(lines,"nelect")
    nele = float(tmp[:-1])
    tmp = search_keyword(lines,"ucvol")
    vol = float(tmp)

    return nele, vol


def parse_command_line():
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser()

    help_str = "Name of the TDCURRENT file"
    parser.add_argument("-c", "--current", help=help_str, required=True)

    help_str = "Name of the TDEFIELD file or dirac if you used an impulse"
    help_str += " electric field (Dirac pulse) see also the ezero parameter"
    help_str += " in that case"
    parser.add_argument("-e", "--efield", help=help_str, required=False)

    help_str = "Amplitude of the electric field - Only used if -e dirac"
    parser.add_argument("-ez", "--ezero",  help=help_str, required=False,
                        type=float, default=1.0)

    help_str = "Name of the Abinit output file (often *.abo)"
    parser.add_argument("-o", "--outfile", help=help_str, required=False)

    help_str = "Direction (x, y or z) of electric field to consider"
    parser.add_argument("-d", "--dir", help=help_str, required=False,
                        default="x")

    help_str = "Cutoff freq. of exponential window [exp(-wcut*t)] (in Ha)"
    parser.add_argument("-wc", "--wcut", help=help_str, required=False,
                        type=float, default=0.0)

    help_str = "Use a third order polynomial damping function"
    help_str += "[1-3(t/tmax)^2+2(t/tmax)^3]"
    parser.add_argument("-damp", "--damp", help=help_str, required=False,
                        action="store_true", default=False)

    help_str = "Remove the first tshift time of current density"
    parser.add_argument("-ts", "--tshift", help=help_str, required=False,
                        type=float, default=0.0)

    help_str = "Stride time step (Default: 1)"
    parser.add_argument("-s", "--stride",  help=help_str,  required=False,
                        type=int, default=1)

    help_str = "If p>1 then zero-padding is used.\n"
    help_str += "The signal is considered to be of length p*ntime (Default: 1)"
    parser.add_argument("-p", "--padding", help=help_str, required=False,
                        type=int, default=1)

    help_str = "More output mostly for testing"
    parser.add_argument("-v", "--verbose", help=help_str, required=False,
                        action="store_true", default=False)

    args = parser.parse_args()

    return args


def check_input_params(args):
    """
    Check input parameters values
    """
    if args.stride < 0:
        sys.exit("Wrong value of stride! It should be larger than zero!")

    if args.padding < 1:
        sys.exit("Wrong value of padding! Should be larger or equal to one!")

    if args.wcut <= 0.:
        sys.exit("Wrong value of wcut! It should be larger than zero!")
    elif args.damp:
        msg = "Conflicting arguments: You need to set either wcut to use "
        msg += "an exponential window or damp to use a damping function!"
        sys.exit(msg)

    if args.tshift < 0:
        sys.exit("Wrong value of tshift! It should be larger than zero!")

    if args.dir != "x" and args.dir != "y" and args.dir != "z":
        sys.exit("Wrong value of dir! Should be x, y or z!")

    if args.ezero < 0.0:
        sys.exit("Wrong value of ezero! It should be larger than zero!")

    # Print summary of parameters
    print()
    print("# Summary of parameters")
    print("Current filename =", args.current)
    if args.efield.strip().lower() == "dirac":
        print("Dirac pulse electric field with amplitude", args.ezero, "au")
    else:
        print("Electric field filename", args.efield)
    print("Electric field direction used :", args.dir)
    print("Abinit output filename =", args.outfile)
    if args.tshift > 0.0:
        print("Time shift applied to the current ts = ", args.tshift)
    if args.stride > 1:
        print("Time stride applied to the current.")
        print("Reading every", args.tstride, "steps")
    if args.padding > 1:
        print("Use zero-padding for Fourier transform with p=", args.padding)
    if args.wcut > 0.0:
        print("Apply exponential window [exp(-wcut*t)] to the current density\
     with wcut = ", args.wcut, "Ha")
    elif args.damp:
        print("Apply damping function [1+3(t/tmax)^2-2(t/tmax)^3]\
     to the current density")

    calc_conducti = True
    dirac_pulse = False
    if args.efield is None:
        calc_conducti = False
        print()
        print("Conductivity and dielectric function will not be computed\
 since the electric field was not provided.")
    elif args.efield.strip().lower() == "dirac":
        dirac_pulse = True

    return calc_conducti, dirac_pulse


# -------------------- Main ------------------------

# Read command line arguments
args = parse_command_line()

# Check input parameters
calc_conducti, dirac_pulse = check_input_params(args)

# Read current
data = np.loadtxt(args.current.strip())
time = data[:,1]
current = data[:,2:5]
n0 = len(current)
dt0 = time[1] - time[0]

# Read electric field
if calc_conducti and not dirac_pulse:
    data = np.loadtxt(args.efield.strip())
    if args.dir == "x":
        efield = data[:,2]
    if args.dir == "y":
        efield = data[:,3]
    if args.dir == "z":
        efield = data[:,4]
    if np.all(np.abs(efield) < 1e-10):
        sys.exit("It seems that the electric field is zero at all times!\n\
If that is because you used an impulse electric field (Dirac pulse)\n\
then you should use the option -e dirac with -ez.")

# Remove tshift
current = current[time >= args.tshift]
if calc_conducti and not dirac_pulse:
    efield = efield[time >= args.tshift]
time = time[time >= args.tshift]

# Apply time stride
current = current[::args.stride]
dt = args.stride*dt0
n = len(current)
time = np.arange(0,n)*dt

# For zero-padding
nfft = n*args.padding

# Print summary
print()
print("Remove timesteps before tshift and apply time striding if required.")
print()
print("# Time-dependent current")
print("Original number of steps n0 =", n0)
print("Actual number of steps used n =", n)
print("Original timestep dt0 =", dt0, "au", "=", dt0*au2as, "as")
print("Actual timestep used dt =", dt, "au", "=", dt*au2as, "as")
print("Original length of trajectory tmax0 =",n0*dt0, "au",
      "=", n0*dt0*au2fs, "fs")
print("Actual length of trajectory tmax =",n*dt, "au", "=", n*dt*au2fs, "fs")

print()
print("# Fourier transform")
print("Number of frequency point used for FFT nfft =", nfft)
print("Minimum angular frequency wmin = dw =", 2*np.pi/(nfft*dt), "au",
                                         " =", 2*np.pi*au2ev/(nfft*dt), "eV")
print("Maximum angular frequency wmax = nfft*dw =", np.pi/dt, "au",
                                              " =", np.pi*au2ev/dt, "eV")

# Perform Fourier transform of current density
w, current_ft_x = fourier_direct(time, current[:,0], args.wcut, nfft)
w, current_ft_y = fourier_direct(time, current[:,1], args.wcut, nfft)
w, current_ft_z = fourier_direct(time, current[:,2], args.wcut, nfft)

# write out Fourier transform of current density
if args.verbose:
    header = "Input current density\n"
    header += "all quantities are in Hartree atomic units.\n"
    header += "time, J_x(t), J_y(t), J_z(t)"
    np.savetxt("current.dat",
               np.vstack([time,current[:,0],current[:,1],current[:,2]]).T,
               header=header)
if args.verbose or not calc_conducti:
    header = "FFT of current density\n"
    header += "all quantities are in Hartree atomic units.\n"
    header += "w(ang. freq.), Re[J_x(w)], Im[J_x(w)], Re[J_y(w)], Im[J_y(w)],"
    header += " Re[J_z(w)], Im[J_z(w)]"
    np.savetxt("current_ft.dat",
               np.vstack([w,np.real(current_ft_x),np.imag(current_ft_x),
                            np.real(current_ft_y),np.imag(current_ft_y),
                            np.real(current_ft_z),np.imag(current_ft_z)]).T,
               header=header)

# Perform Fourier transform of electric field
if calc_conducti:
    if dirac_pulse:
        nw = len(w)
        efield_ft = args.ezero*np.ones(nw)+1j*np.zeros(nw)
    else:
        w, efield_ft = fourier_direct(time,efield,0.0,nfft)
        nw = len(w)
        # check that we only get a real part
        if np.any(np.abs(np.imag(efield)) > 1e-10):
            print("Warning: FFT of the electric field seems to have a non \
zero imaginary part.")
        # write out electric field
        if args.verbose:
            header = "Input electric field\n"
            header += "all quantities are in Hartree atomic units.\n"
            header += "time, E(t)"
            np.savetxt("efield.dat", np.vstack([time,efield]).T, header=header)

    # write out Fourier transform of electric field
    if args.verbose:
        header = "FFT of electric field\n"
        header += "all quantities are in Hartree atomic units.\n"
        header += "w, Re[E(w)], Im[E(w)]"
        np.savetxt("efield_ft.dat",
                   np.vstack([w,np.real(efield_ft),np.imag(efield_ft)]).T,
                   header=header)

    # Compute optical conductivity
    nmin = int(nw/2)+1; nmax = len(w)
    w = w[nmin:nmax]
    sigma_x = np.divide(current_ft_x[nmin:nmax],
                        np.real(efield_ft[nmin:nmax]),
                        where=np.real(efield_ft[nmin:nmax]) != 0)
    sigma_y = np.divide(current_ft_y[nmin:nmax],
                        np.real(efield_ft[nmin:nmax]),
                        where=np.real(efield_ft[nmin:nmax]) != 0)
    sigma_z = np.divide(current_ft_z[nmin:nmax],
                        np.real(efield_ft[nmin:nmax]),
                        where=np.real(efield_ft[nmin:nmax]) != 0)

    # f-sum rule test
    if args.outfile is not None:
        nele, vol = parse_output(args.outfile)
        print()
        print("# f-sum rule (Thomas-Reiche-Kuhn) for conductivity:")
        wp2 = 4*np.pi*nele/vol
        pref = 8./wp2
        # x-direction
        I_sigma = trapz(np.real(sigma_x),x=w)
        if (I_sigma > 1e-6):
            print("x-direction:")
            print("Without high freq. correction", I_sigma*pref)
            I_sigma = trapz(np.real(sigma_x)-np.real(sigma_x)[-1],x=w)
            print("With high freq. correction", I_sigma*pref)
        # y-direction
        I_sigma = trapz(np.real(sigma_y),x=w)
        if (I_sigma > 1e-6):
            print("y-direction:")
            print("Without high freq. correction", I_sigma*pref)
            I_sigma = trapz(np.real(sigma_y)-np.real(sigma_y)[-1],x=w)
            print("With high freq. correction", I_sigma*pref)
        # y-direction
        I_sigma = trapz(np.real(sigma_z),x=w)
        if (I_sigma > 1e-6):
            print("z-direction:")
            print("Without high freq. correction", I_sigma*pref)
            I_sigma = trapz(np.real(sigma_z)-np.real(sigma_z)[-1],x=w)
            print("With high freq. correction", I_sigma*pref)
    else:
        print()
        msg = "Note that the f-sum rule (Thomas-Reiche-Kuhn) for conductivity"
        msg += "cannot be checked if you don't provide the abinit output file."
        print(msg)

    # write out conductivity
    header = "Optical conductivity sigma(w)\n w [Ha], w [eV], "
    header += "Re[sigma_x] [au], Im[sigma_x] [au], "
    header += "Re[sigma_y] [au], Im[sigma_y] [au], "
    header += "Re[sigma_z] [au], Im[sigma_z] [au], "
    header += "Re[sigma_x] [(Ohm.cm)^-1], Im[sigma_x] [(Ohm.cm)^-1], "
    header += "Re[sigma_y] [(Ohm.cm)^-1], Im[sigma_y] [(Ohm.cm)^-1], "
    header += "Re[sigma_z] [(Ohm.cm)^-1], Im[sigma_z] [(Ohm.cm)^-1]"
    np.savetxt("conductivity.dat",
               np.vstack([w,w*au2ev,
                          np.real(sigma_x), np.imag(sigma_x),
                          np.real(sigma_y), np.imag(sigma_y),
                          np.real(sigma_z), np.imag(sigma_z),
                          np.real(sigma_x)*au2Ohmcm1,
                          np.imag(sigma_x)*au2Ohmcm1,
                          np.real(sigma_y)*au2Ohmcm1,
                          np.imag(sigma_y)*au2Ohmcm1,
                          np.real(sigma_z)*au2Ohmcm1,
                          np.imag(sigma_z)*au2Ohmcm1]).T,
                          header=header)

    # Compute dielectric tensor
    eps_x = 4.0*np.pi*1j*np.divide(sigma_x, w)
    eps_y = 4.0*np.pi*1j*np.divide(sigma_y, w)
    eps_z = 4.0*np.pi*1j*np.divide(sigma_z, w)
    if args.dir == "x":
        eps_x = 1.0 + eps_x
    if args.dir == "y":
        eps_y = 1.0 + eps_y
    if args.dir == "z":
        eps_z = 1.0 + eps_z

    # write out dielectric tensor
    header = "Dielectric function epsilon(w) (unitless ie epsilon/epsilon_0)\n"
    header += "w [Ha], w [eV], "
    header += "Re[eps_x], Im[eps_x] Re[eps_y] Im[eps_y], Re[eps_z], Im[eps_z]"
    np.savetxt("dielectric.dat",
               np.vstack([w,w*au2ev,
                          np.real(eps_x), np.imag(eps_x),
                          np.real(eps_y), np.imag(eps_y),
                          np.real(eps_z), np.imag(eps_z)]).T,
               header=header)
