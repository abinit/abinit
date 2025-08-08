#!/usr/bin/env python3
"""
Script to compute the optical conductivity and the dielectric
function from the time-dependent current density obtained
after a RT-TDDFT run with ABINIT
----------------
Partially based on the script in the tutorials of the
exciting code
----------------
F. Brieuc
"""

import sys
import argparse
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
ev2j   = 1.602176634e-19 # J -> eV
au2s   = hb/Eh           # au -> second
au2as  = au2s * 1e17     # au -> attosecond
au2fs  = au2s * 1e15     # au -> femtosecond
au2Vpm = Eh/(qe*a0)      # au -> Volt/meter 
au2ev  = Eh / ev2j       # au -> eV
au2Ohm = me*a0**2/(qe**2*au2s) # au -> Ohm
au2Ohmcm1 = 1/(au2Ohm*a0*1e2)  # au -> 1/(Ohm*cm)

# ------------------- Functions -----------------------
def fourier_direct(time,signal,wcut,nfft):
    """
    Direct Fourier transform - from time to angular frequency
    Input:
     - time: time (array)
     - signal: signal (array) f(t) to be transformed
     - wcut: cut-off frequency (in Ha) for the exponential window (if > 0)
     - nfft: number of points to use in the FFT
             if greater than the length of signal then zero-padding is used
    Output:
      -w, F: Tuple containing angular frequencies (in Ha) and the
        fourier transform of signal*filter
    ------
    If wcut > 0 then the signal is multiplied in time by exp(-wcut*t)
    which leads to a convolution with a Lorentzian function in frequency space.
    Indeed the Fourier transform of exp(-wcut*t) is w_cut/(w_cut^2+w^2).
    Thus this window function generates a Lorentzian broadening in frequency
    with a FWHM given by 2*w_cut.
    If wcut == -1 use a third order polynomial damping function as filter.
    From Yabana et al., Phys. Stat. Sol. (B) 243,1121 (2006)
    """
    ntime =len(time)
    if len(signal) != ntime:
        sys.exit("!Error - in fourier_direct - signal and time do not have the same length")
    dt = time[1] - time[0]
    # Apply filter if required
    # Exponential windaow (exp(-wcut*t))
    if wcut > 0:
        filter = np.exp(np.multiply(time,-wcut))
        f = np.multiply(filter,signal)
    # Damping function (1-3(t/tmax)^2+2(t/tmax)^3)
    elif wcut == -1:
        tmax = time[-1]
        filter = np.ones(len(time)) - np.multiply(time**2,3.0/tmax**2) + np.multiply(time**3,2.0/tmax**3)
        f = np.multiply(filter,signal)
    else:
        f = signal
    # Zero-padding if required
    dn = nfft - ntime
    if dn > 0:
        ff = np.zeros(nfft)
        ff[0:ntime] = f
    else:
        ff = f
    # Define angular frequencies
    w = 2*np.pi*fftfreq(nfft,dt)
    # Fourier transform
    ft = np.conj(np.multiply(fft(ff),dt))
    # Sort to get the right order in frequency
    w,ft = zip(*sorted(zip(w,ft)))

    return np.array(w), np.array(ft)

def search_keyword(lines, keyword):
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
    
    tmp = search_keyword(lines,'nelect')
    nele = float(tmp[:-1])
    tmp = search_keyword(lines,'ucvol') 
    vol = float(tmp)

    return nele, vol

# -------------------- Main ------------------------

# Read command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--current', help='Name of the TDCURRENT file', required=True)
parser.add_argument('-e', '--efield',  help='Name of the TDEFIELD file\
                                             or dirac if you used an impulse electric field (Dirac pulse)\
                                             see also the ezero parameter in that case', required=False)
parser.add_argument('-ez', '--ezero',  help='Amplitude of the electric field - Only used if -e dirac', required=False,
                                       type=float, default=1.0)
parser.add_argument('-o', '--outfile', help='Name of the Abinit output file (often *.abo)', required=False)
parser.add_argument('-d', '--dir',     help='Direction (x, y or z) of electric field to consider (divide by E_dir)', required=False,
                                       type=str, default='x')
parser.add_argument('-wc', '--wcut',   help='Cutoff angular frequency exponential window [exp(-wcut*t)] (in Ha)', required=False,
                                       type=float, default=0.0)
parser.add_argument('-damp', '--damp', help='Use a third order polynomial damping function [1-3(t/tmax)^2+2(t/tmax)^3]', required=False,
                                       action='store_true', default=False)
parser.add_argument('-ts', '--tshift', help='Remove the first tshift time of current density', required=False,
                                       type=float, default=0.0)
parser.add_argument('-s', '--stride',  help='Stride time step (Default: 1)', required=False, type=int, default=1)
parser.add_argument('-p', '--padding', help='If p>1 then zero-padding is used.\n\
The signal is considered to be of length p*ntime (Default: 1)',
                                       required=False, type=int, default=1)
parser.add_argument('-v', '--verbose', help='More output mostly for testing', required=False,
                                       action='store_true', default=False)
args = parser.parse_args()

# Check input parameters
if args.stride < 0:
    sys.exit("Wrong value of stride! It should be larger than zero!")
else:
    stride=int(args.stride)

if args.padding < 1:
    sys.exit("Wrong value of padding! It should be larger or equal to one!")
else:
    padding=int(args.padding)

if args.wcut < 0:
    sys.exit("Wrong value of wcut! It should be larger than zero!")
else:
    wcut=float(args.wcut)

if (args.damp and wcut > 0):
    sys.exit("Conflicting arguments: You need to set either wcut to use an exponential window or damp to use a damping function")

if args.damp:
    wcut = -1

if args.tshift < 0:
    sys.exit("Wrong value of tshift! It should be larger than zero!")
else:
    tshift=float(args.tshift)

if args.outfile is None:
    outfile = False
else:
    outfile = True

if args.dir == 'x':
    dir=0
elif args.dir == 'y':
    dir=1
elif args.dir == 'z':
    dir=2
else:
    sys.exit("Wrong value of dir! Should be x, y or z!")

# Print summary of parameters
print("")
print("# Summary of parameters")
print("Current filename =", args.current)
if args.efield == "dirac":
    print("Dirac pulse electric field with amplitude", args.ezero, "au")
else:
    print("Electric field filename", args.efield)
    print("Electric field direction used :", args.dir)
print("Abinit output filename =", args.outfile)
print("Time shift applied to the current ts = ", args.tshift)
if args.stride > 1:
    print("Time stride applied to the current. Reading every", args.tstride, "steps")
if args.padding > 1:
    print("Use zero-padding for Fourier transform with p=", args.padding)
if wcut > 0:
    print("Apply exponential window [exp(-wcut*t)] to the current density before Fourier transform with wcut = ", wcut, "Ha")
elif wcut == -1:
    print("Apply damping function [1+3(t/tmax)^2-2(t/tmax)^3] to the current density before Fourier transform")

# Read current
data = np.loadtxt(args.current.strip())
time = data[:,1]
current = data[:,2:5]
n0 = len(current)
dt0 = time[1] - time[0]

# Read electric field
calc_conducti = True
dirac_pulse = False
if args.efield is None:
    calc_conducti = False
    print("")
    print("Conductivity and dielectric function will not be computed since the electric field was not provided.")
elif args.efield.strip().lower() == 'dirac':
    dirac_pulse = True
    if args.ezero is None:
        sys.exit("You should give the amplitude of the electric field -ez <E_0> if you use a Dirac pulse.")
    elif args.ezero > 0:
        ezero = args.ezero
    else:
        sys.exit("Wrong value of ezero! It should be larger than zero!")
else:
    data = np.loadtxt(args.efield.strip())
    efield = data[:,2+dir]
    if np.all(np.abs(efield) < 1e-10):
        sys.exit("It seems that the electric field is zero at all times!\n\
If that is because you used an impulse electric field (Dirac pulse)\n\
then you should use the option -e dirac with -ez.")

# Remove tshift
current = current[time>=tshift]
if calc_conducti and not dirac_pulse:
    efield = efield[time>=tshift]
time = time[time>=tshift]

# Apply time stride
current = current[::stride]
dt = stride*dt0
n = len(current)
time = np.arange(0,n)*dt

# For zero-padding
nfft = n*padding

# Print summary
print("")
print("Removing timesteps before tshift and applying time striding if required.")
print("")
print("# Time-dependent current")
print("Original number of steps n0 =", n0)
print("Actual number of steps used n =", n)
print("Original timestep dt0 =", dt0, "au", "=", dt0*au2as, "as")
print("Actual timestep used dt =", dt, "au", "=", dt*au2as, "as")
print("Original length of trajectory tmax0 =",n0*dt0, "au", "=", n0*dt0*au2fs, "fs")
print("Actual length of trajectory tmax =",n*dt, "au", "=", n*dt*au2fs, "fs")

print("")
print("# Fourier transform")
print("Number of frequency point used for FFT nfft =", nfft)
print("Minimum angular frequency wmin = dw = ", 2*np.pi/(nfft*dt), "au", "=", 2*np.pi*au2ev/(nfft*dt), "eV")
print("Maximum angular frequency wmax = nfft*dw = ", np.pi/dt, "au", "=", np.pi*au2ev/dt, "eV")

# Perform Fourier transform of current density
w, current_ft_x = fourier_direct(time,current[:,0],wcut,nfft)
w, current_ft_y = fourier_direct(time,current[:,1],wcut,nfft)
w, current_ft_z = fourier_direct(time,current[:,2],wcut,nfft)

# write out Fourier transform of current density
if (args.verbose):
    header = "Input current density\nall quantities are in Hartree atomic units.\ntime, J_x(t), J_y(t), J_z(t)"
    np.savetxt("current.dat",np.vstack([time,current[:,0],current[:,1],current[:,2]]).T,header=header)
    header = "FFT of current density\nall quantities are in Hartree atomic units.\n\
w(ang. freq.), Re[J_x(w)], Im[J_x(w)], Re[J_y(w)], Im[J_y(w)], Re[J_z(w)], Im[J_z(w)]"
    np.savetxt("current_ft.dat", np.vstack([w,np.real(current_ft_x),np.imag(current_ft_x),
                                              np.real(current_ft_y),np.imag(current_ft_y),
                                              np.real(current_ft_z),np.imag(current_ft_z)]).T, header=header)

# Perform Fourier transform of electric field
if calc_conducti:
    if dirac_pulse:
        nw = len(w)
        efield_ft = ezero*np.ones(nw)+1j*np.zeros(nw)
    else:
        w, efield_ft = fourier_direct(time,efield,0.0,nfft)
        nw = len(w)
        # check that we get a real part only
        if np.any(np.abs(np.imag(efield))>1e-10):
            print("Warning: FFT of the electric field seems to have a non zero imaginary part.")
        # write out electric field
        if (args.verbose):
            header = "Input electric field\nall quantities are in Hartree atomic units.\ntime, E(t)"
            np.savetxt("efield.dat",np.vstack([time,efield]).T,header=header)

    # write out Fourier transform of electric field
    if (args.verbose):
        header = "FFT of electric field\nall quantities are in Hartree atomic units.\nw, Re[E(w)], Im[E(w)]"
        np.savetxt("efield_ft.dat",np.vstack([w,np.real(efield_ft),np.imag(efield_ft)]).T,header=header)

    # Compute optical conductivity
    nmin = int(nw/2)+1; nmax = len(w)
    w = w[nmin:nmax]
    sigma_x = np.divide(current_ft_x[nmin:nmax], np.real(efield_ft[nmin:nmax]), where=np.real(efield_ft[nmin:nmax])!=0)
    sigma_y = np.divide(current_ft_y[nmin:nmax], np.real(efield_ft[nmin:nmax]), where=np.real(efield_ft[nmin:nmax])!=0)
    sigma_z = np.divide(current_ft_z[nmin:nmax], np.real(efield_ft[nmin:nmax]), where=np.real(efield_ft[nmin:nmax])!=0)

    # f-sum rule test
    if outfile:
        nele, vol = parse_output(args.outfile)
        print('')
        print("# f-sum rule (Thomas-Reiche-Kuhn) for conductivity:")
        wp2 = 4*np.pi*nele/vol 
        # x-direction
        I_sigma = trapz(np.real(sigma_x),x=w)
        if (I_sigma > 1e-6):
            print("x-direction:")
            print("Without high freq. correction", I_sigma*8/wp2)
            I_sigma = trapz(np.real(sigma_x)-np.real(sigma_x)[-1],x=w)
            print("With high freq. correction", I_sigma*8/wp2)
        # y-direction
        I_sigma = trapz(np.real(sigma_y),x=w)
        if (I_sigma > 1e-6):
            print("y-direction:")
            print("Without high freq. correction", I_sigma*8/wp2)
            I_sigma = trapz(np.real(sigma_y)-np.real(sigma_y)[-1],x=w)
            print("With high freq. correction", I_sigma*8/wp2)
        # y-direction
        I_sigma = trapz(np.real(sigma_z),x=w)
        if (I_sigma > 1e-6):
            print("z-direction:")
            print("Without high freq. correction", I_sigma*8/wp2)
            I_sigma = trapz(np.real(sigma_z)-np.real(sigma_z)[-1],x=w)
            print("With high freq. correction", I_sigma*8/wp2)
    else:
        print('')
        print("Note that the f-sum rule (Thomas-Reiche-Kuhn) for conductivity cannot be checked if you don't provide the abinit output file.")

    # write out conductivity
    header = "Optical conductivity sigma(w)\n w [Ha], w [eV], \
Re[sigma_x] [au], Im[sigma_x] [au], Re[sigma_y] [au], Im[sigma_y] [au], Re[sigma_z] [au], Im[sigma_z] [au], \
Re[sigma_x] [(Ohm.cm)^-1], Im[sigma_x] [(Ohm.cm)^-1], Re[sigma_y] [(Ohm.cm)^-1], Im[sigma_y] [(Ohm.cm)^-1], \
Re[sigma_z] [(Ohm.cm)^-1], Im[sigma_z(w)] [(Ohm.cm)^-1]" 
    np.savetxt("conductivity.dat",np.vstack([w,w*au2ev,\
                                  np.real(sigma_x), np.imag(sigma_x),\
                                  np.real(sigma_y), np.imag(sigma_y),\
                                  np.real(sigma_z), np.imag(sigma_z),\
                                  np.real(sigma_x)*au2Ohmcm1, np.imag(sigma_x)*au2Ohmcm1,\
                                  np.real(sigma_y)*au2Ohmcm1, np.imag(sigma_y)*au2Ohmcm1,\
                                  np.real(sigma_z)*au2Ohmcm1, np.imag(sigma_z)*au2Ohmcm1]).T,
                                  header=header)

    # Compute dielectric tensor
    eps_x = 4.0*np.pi*1j*np.divide(sigma_x, w)
    eps_y = 4.0*np.pi*1j*np.divide(sigma_y, w)
    eps_z = 4.0*np.pi*1j*np.divide(sigma_z, w)
    if dir==0:
        eps_x = 1.0 + eps_x
    if dir==1:
        eps_y = 1.0 + eps_y
    if dir==2:
        eps_z = 1.0 + eps_z

    # write out dielectric tensor
    header = "Dielectric function epsilon(w) (unitless ie epsilon/epsilon_0)\n w [Ha], w [eV], \
Re[eps_x], Im[eps_x], Re[eps_y], Im[eps_y], Re[eps_z], Im[eps_z]" 
    np.savetxt("dielectric.dat",np.vstack([w,w*au2ev,\
                                np.real(eps_x), np.imag(eps_x), np.real(eps_y), np.imag(eps_y),\
                                np.real(eps_z), np.imag(eps_z)]).T, header=header)
