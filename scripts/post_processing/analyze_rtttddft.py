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
import numpy as np
import argparse
from scipy.fftpack import fft, fftfreq

# -------------------- Units --------------------------
au2as = 2.41888432650478
au2fs = 0.0241888432650478
au2Vpm = 514220624373.482
au2ev = 27.2113838565563
au2Ohmcm = 45934.859262851380

# ------------------- Functions -----------------------
def fourier_direct(time,signal,wcut,nfft):
    """
    Direct Fourier transform - from time to angular frequency
    Input:
     - time: time (array)
     - signal: signal (array) f(t) to be transformed
     - wcut: cut-off frequency (in Ha) for the exponential window
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
    """
    ntime =len(time)
    if len(signal) != ntime:
        sys.exit("!Error - in fourier_direct - signal and time do not have the same length")
    dt = time[1] - time[0]
    # Apply filter if required
    if wcut > 0:
        filter = np.exp(np.multiply(time,-wcut))
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
    # Sort to get the righot order in frequency
    w,ft = zip(*sorted(zip(w,ft)))

    return np.array(w), np.array(ft)

# -------------------- Main ------------------------

# Read command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--current', help='Name of the TDCURRENT file', required=True)
parser.add_argument('-e', '--efield',  help='Name of the TDEFIELD file\
                                             or dirac if you used an impulse electric field (Dirac pulse)\
                                             see also the ezero parameter in that case', required=False)
parser.add_argument('-ez', '--ezero',  help='Amplitude of the electric field - Only used if -e dirac', required=False,
                                       type=float, default=1.0)
parser.add_argument('-d', '--dir',     help='Direction (x, y or z) of electric field to consider (divide by E_dir)', required=False,
                                       type=str, default='x')
parser.add_argument('-wc', '--wcut',   help='Cutoff angular frequency exponential window [exp(-wcut*t)] (in Ha)', required=False,
                                       type=float, default=0.04)
parser.add_argument('-ts', '--tshift', help='Remove the first tshift time of current density', required=False,
                                       type=float, default=0.0)
parser.add_argument('-s', '--stride',  help='Stride time step (Default: 1)', required=False, type=int, default=1)
parser.add_argument('-p', '--padding', help='If p>1 then zero-padding is used.\n\
The signal is considered to be of length p*ntime (Default: 1)',
                                       required=False, type=int, default=1)
parser.add_argument('-v', '--verbose',  help='More output mostly for testing', required=False,
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

if args.tshift < 0:
    sys.exit("Wrong value of tshift! It should be larger than zero!")
else:
    tshift=int(args.tshift)

if args.dir == 'x':
    dir=0
elif args.dir == 'y':
    dir=1
elif args.dir == 'z':
    dir=2
else:
    sys.exit("Wrong value of dir! Should be x, y or z!")

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
print("Maximum angular frequency wax = nfft*dw = ", np.pi/dt, "au", "=", np.pi*au2ev/dt, "eV")

# Perform Fourier transform of current density
w, current_ft_x = fourier_direct(time,current[:,0],wcut,nfft)
w, current_ft_y = fourier_direct(time,current[:,1],wcut,nfft)
w, current_ft_z = fourier_direct(time,current[:,2],wcut,nfft)

# write out Fourier transform of current density
nw = len(w)
if (args.verbose):
    header = "Input current density\nall quantities are in Hartree atomic units.\ntime, J_x(t), J_y(t), J_z(t)"
    np.savetxt("current.dat",np.vstack([time,current[:,0],current[:,1],current[:,2]]).T,header=header)
header = "FFT of current density\nall quantities are in Hartree atomic units.\n\
w(ang. freq.), Re[J_x(w)], Im[J_x(w)], Re[J_y(w)], Im[J_y(w)], Re[J_z(w)], Im[J_z(w)]"
np.savetxt("current_ft.dat", np.vstack([w,np.real(current_ft_x[0:nw]),np.imag(current_ft_x[0:nw]),
                                          np.real(current_ft_y[0:nw]),np.imag(current_ft_y[0:nw]),
                                          np.real(current_ft_z[0:nw]),np.imag(current_ft_z[0:nw])]).T, header=header)

# Perform Fourier transform of electric field
if calc_conducti:
    if dirac_pulse:
        efield_ft = ezero*np.ones(nw)+1j*np.zeros(nw)
    else:
        w, efield_ft = fourier_direct(time,efield,wcut,nfft)
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
    sigma_x = np.divide(current_ft_x, np.real(efield_ft), where=np.real(efield_ft)!=0)
    sigma_y = np.divide(current_ft_y, np.real(efield_ft), where=np.real(efield_ft)!=0)
    sigma_z = np.divide(current_ft_z, np.real(efield_ft), where=np.real(efield_ft)!=0)

    # write out conductivity
    header = "Optical conductivity\n w [Ha], w [eV], \
Re[sigma_x] [au], Im[sigma_x] [au], Re[sigma_y] [au], Im[sigma_y] [au], Re[sigma_z] [au], Im[sigma_z(w)] [au], \
Re[sigma_x] [Ohm.cm], Im[sigma_x] [Ohm.cm], Re[sigma_y] [Ohm.cm], Im[sigma_y] [Ohm.cm], Re[sigma_z] [Ohm.cm], Im[sigma_z(w)] [Ohm.cm]"
    np.savetxt("conductivity.dat",np.vstack([w[int(nw/2)+1:nw],w[int(nw/2)+1:nw]*au2ev,\
                                  np.real(sigma_x[int(nw/2)+1:nw]), np.imag(sigma_x[int(nw/2)+1:nw]),\
                                  np.real(sigma_y[int(nw/2)+1:nw]), np.imag(sigma_y[int(nw/2)+1:nw]),\
                                  np.real(sigma_z[int(nw/2)+1:nw]), np.imag(sigma_z[int(nw/2)+1:nw]),\
                                  np.real(sigma_x[int(nw/2)+1:nw])*au2Ohmcm, np.imag(sigma_x[int(nw/2)+1:nw])*au2Ohmcm,\
                                  np.real(sigma_y[int(nw/2)+1:nw])*au2Ohmcm, np.imag(sigma_y[int(nw/2)+1:nw])*au2Ohmcm,\
                                  np.real(sigma_z[int(nw/2)+1:nw])*au2Ohmcm, np.imag(sigma_z[int(nw/2)+1:nw])*au2Ohmcm]).T,
                                  header=header)

    # Compute dielectric tensor
    eps_x = 4.0*np.pi*1j*np.divide(sigma_x[int(nw/2)+1:nw], w[int(nw/2)+1:nw])
    eps_y = 4.0*np.pi*1j*np.divide(sigma_y[int(nw/2)+1:nw], w[int(nw/2)+1:nw])
    eps_z = 4.0*np.pi*1j*np.divide(sigma_z[int(nw/2)+1:nw], w[int(nw/2)+1:nw])
    if dir==0:
        eps_x = 1.0 + eps_x
    if dir==1:
        eps_y = 1.0 + eps_y
    if dir==2:
        eps_z = 1.0 + eps_z

    # write out dielectric tensor
    header = "Dielectric function\n w [Ha], w [eV], \
Re[eps_x] [au], Im[eps_x] [au], Re[eps_y] [au], Im[eps_y] [au], Re[eps_z] [au], Im[eps_z(w)] [au]"
    np.savetxt("dielectric.dat",np.vstack([w[int(nw/2)+1:nw],w[int(nw/2)+1:nw]*au2ev,\
                                np.real(eps_x), np.imag(eps_x), np.real(eps_y), np.imag(eps_y),\
                                np.real(eps_z), np.imag(eps_z)]).T, header=header)
