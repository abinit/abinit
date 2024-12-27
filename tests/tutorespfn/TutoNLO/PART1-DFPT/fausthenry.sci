// Scilab script to compute the Faust-Henry coefficient in
// zinc-blende semiconductors
// See Eqs. (2.90) - (2.94) of Light Scattering in Solids II.

clear

// ----------------------------------------------------------------------------
// ------- Parameters taken from the output files of ABINIT and ANADDB --------
// ----------------------------------------------------------------------------

// Name of the semiconductor
sc = 'AlAs'

// Unit cell volume in bohr^3
ucvol =  3.0109659E+02

// Atomic masses
m1 = 0.26981539000000D+02
m2 = 0.74921590000000D+02

// Born effective charge
zeff = 2.105999E+00

// Electronic dielectric constant
epsel = 9.94846084

// Non-linear optical susceptibilities (d_123 in pm/V)
d0mks = 32.772254

// First-order change in the electronic dielectric susceptibility
// induced by a displacement of atom 1 (Bohr^-1)
dchidtau = -0.099889084

// Energy (hartree) of TO and LO modes
wto = 1.641481E-03
wlo = 1.791368E-03

// ----------------------------------------------------------------------------
// -------------------------- Fundamental constants ---------------------------
// ----------------------------------------------------------------------------

e_Cb = 1.602176462e-19;
amu_emass = 1.66053873e-27/9.10938188e-31;
eps0  =8.854187817e-12;
Ha_cmm1 = 219474.6313710;
Bohr_Ang = 0.5291772083;

// Factor to convert the nlo susceptibilities from a.u. to pm/V
fac = 4*%pi*(Bohr_Ang**2)*1.0e-8*4*%pi*eps0/e_Cb

// Unit for output
unix('if [ -f fausthenry.out ]; then rm -f fausthenry.out; fi')
unit = file('open','fausthenry.out');

// ----------------------------------------------------------------------------
// ----- Computation of the Faust-Henry coefficient and related quantities ----
// ----------------------------------------------------------------------------

mes = ' Raman scattering in ' + sc
write(unit,mes)
mes = ' ========================'
write(unit,mes)
write(unit,'')

mes = ' General parameters:'
write(unit,mes)
mes = '                       epsel = ' + string(epsel)
write(unit,mes)
mes = '                d_123 (pm/V) = ' + string(d0mks)
write(unit,mes)
mes = '                        zeff = ' + string(zeff)
write(unit,mes)
mes = '                 wto (cm^-1) = ' + string(wto*Ha_cmm1)
write(unit,mes)
mes = '                 wlo (cm^-1) = ' + string(wlo*Ha_cmm1)
write(unit,mes)
write(unit,'')

d0au = d0mks/fac;
mu = amu_emass*m1*m2/(m1 + m2);

cfausth = zeff*dchidtau/(4*d0au*mu*wto**2);
mes = ' Faust-Henry coefficient:'
write(unit,mes)
mes = '                           C = ' + string(cfausth)
write(unit,mes)
write(unit,'')

ato = ucvol*dchidtau*Bohr_Ang**2;
alo = ato*(1-(wlo**2 - wto**2)/(cfausth*wto**2))
mes = ' Raman polarizability (A^2):' 
write(unit,mes)
mes = '                         ato = ' + string(ato)
write(unit,mes)
mes = '                         alo = ' + string(alo)
write(unit,mes)
write(unit,'')

dto = dchidtau;
dlo = dchidtau*(1-(wlo**2 - wto**2)/(cfausth*wto**2));
mes = ' Raman tensor of TO and LO modes:'
write(unit,mes)
mes = '                         dto = ' + string(dto)
write(unit,mes)
mes = '                         dlo = ' + string(dlo)
write(unit,mes)
mes = '                     dlo/dto = ' + string(dlo/dto)
write(unit,mes)
write(unit,'')

enhanc = (wto/wlo)*(dlo/dto)**2
mes = ' Enhancement factor:'
write(unit,mes)
mes = '      (wto/wlo)*(dlo/dto)**2 = ' + string(enhanc)
write(unit,mes)
write(unit,'')

rijk = -4*d0mks*(1+cfausth)/epsel**2,
mes = ' Electrooptic tensor (pm/V):'
write(unit,mes)
mes = '                        rijk = ' + string(rijk)
write(unit,mes)

// ----------------------------------------------------------------------------
file('close',unit)
quit
