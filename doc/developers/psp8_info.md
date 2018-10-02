## Information on the format 8 for pseudopotentials.

(Note: the implementation of format 8 was done by D. R. Hamann).

The format 8 for ABINIT norm-conserving separable pseudopotentials 
([[cite:Hamann1979]], [[cite:Hamann1989]] and [[cite:Kleinman1982]])
is designed to allow users
who wish to experiment with pseudopotentials, possibly with non-standard
features, to have great flexibility in doing so. It does not correspond 
to any publicly available tabulation. The open-source ONCVPSP package available at 
[www.mat-simresearch.com](http://www.mat-simresearch.com) produces pseudopotentials in this format.  

An annotated example is presented in detail below, followed by a second 
example incorporating spin-orbit coupling.  An extended discussion follows 
below the examples, with a separate section discussing spin-orbit as produced by ONCVPSP.

When first producing a pseudopotential in this format, it would be a very
good idea to do a well-converged calculation of an isolated atom in a big
box to confirm that you are correctly recreating the atomic valence levels you intend.

The format is best explained by an example, which is presented in
detail below.  (The example is excerpted from `Psps_for_tests/20ca_sic.drh`.
Another example of pspcod=8 is `Psps_for_tests/8o_sic.drh`.)
All but the last data block of this file are read in subroutine psp8in.f.
The "comment" (#) lines below are not permitted in the actual file. Nor are blank lines.

    Header

    line 1: Title
    line 2: atomic number, pseudoion charge, date
    line 3: pspcod=8 identifies this format
            pspxc=2 identifies the exchange-correlation functional (see
            ixc list in the input variables documentation)
            lmax=2 gives the largest angular momentum potential present
            lloc=4 (which is >lmax) indicates that the local potential
            is independent of that for any angular momentum.  In other
            cases, it can be <=lmax, and the local potential can be
            that of a particular angular momentum channel.  (The choice
            of 4 is arbitrary, but the index of the corresponding data block
            below must be consistent with the value assigned in the header.)
            mmax=313, giving the number of radial grid points (arbitrary).
            r2well=0 is an historical header holdover, and is not used
     line 4: rchrg=6.18.. is the radius beyond which the model core
             charge (if used) is zero or negligible.  The maximum
             radial mesh point must equal or exceed this value.
             fchrg=1.0 signals that a model core charge is present.  Any
             value >0.0 will do, and the value is not used otherwise.
             qchrg=0.0 is a historical holdover, not used.
     line 5: number of Bloechl-Kleinman-Bylander projectors nproj for each
             angular momentum l=0 to l=lmax.  The value of nproj for lloc 
             must be 0 (if lloc <= lmax).
     line 6: extension_switch=0, not initially used, but values >0 may
             be utilized to signal the presence of additonal lines of header
             information for special-purpose use in the future.
             The value 2 is now used to indicate the presence of spin-orbit
             projectors (see SPIN-ORBIT section below).

**EXPERIMENTAL** 
self-interaction-corrected psp for calcium (D. R. Hamann):

    20.0000      2.0000    040701      zatom,zion,pspd
    8     2     2     4   313     0    pspcod,pspxc,lmax,lloc,mmax,r2well
    6.18295050  1.00000000  0.00000000 rchrg fchrg qchrg
    2     1     1     0     0          nproj
    0                                  extension_switch

First data block

first line: angular momentum l=0, ekb(ii) for ii=1,2(=nproj(1))
        following lines labeled 1,2,...,mmax -

*   1st column: radial grid index
*   2nd column: radial grid mesh point (LINEAR mesh starting at 0.0)
*   3rd column: first  BKB projector for l=0
*   4th column: second BKB projector for l=0

The projectors should go smoothly to zero, or decay to negligible
values within the range of the radial grid.
In general there are nproj ekb and projector columns for each l, as
many as you want.
If nproj is zero, AND lloc/= the current l, this block is omitted.
While some may be absent, l blocks must be in ascending order.

    0                       1.8442549429363D+01 -1.5467656752626D+01
    1  0.0000000000000D+00  0.0000000000000D+00  0.0000000000000D+00
    2  2.0000000000000D-02  3.2695476840515D-03 -7.7803072439537D-05
    3  4.0000000000000D-02  6.5413292698932D-03 -1.5573364986856D-04
    .....
    .....

Second data block

first line - angular momentum l=1, ekb(ii) for ii=1,1(=nproj(2))
    following lines labeled 1,2,...,mmax -

* 1st column: radial grid index
* 2nd column: radial grid mesh point
* 3rd column: first (only)  BKB projector for l=1

``` 
    1                       1.0250215833500D+01
    1  0.0000000000000D+00  0.0000000000000D+00
    2  2.0000000000000D-02  4.1259561971920D-05
    3  4.0000000000000D-02  1.6501618805883D-04
    .....
    .....
```

Third data block as above, for l=2

    2                      -4.9318560092991D-01
    1  0.0000000000000D+00  0.0000000000000D+00
    2  2.0000000000000D-02 -2.8426374194219D-04
    3  4.0000000000000D-02 -2.2672423014226D-03
    .....
    .....

Fourth data block

first line - index corresponding to lloc in header.  In this case,
      the local potential does not correspond to that for any angular
      momentum.  lloc can be arbitrary in this case, but MUST be >lmax,
      and this block must follow the blocks for projectors with l <=lmax.
      following lines labeled 1,2,...,mmax -

*       1st column: radial grid index
*       2nd column: radial grid mesh point
*       3rd column: local potential

If lloc<lmax, the corresponding block must occur in its proper l order
among the blocks of nonlocal projectors.
The last line should, ideally, be equal to -zion/rad(mmax), and the
numerical data should approach -zion/r as continuously as possible,
since the Fourier transform is extended to infinity analytically
assumming this functional form.  This is not precisely satisfied for
the experimental potential in the example, and the result will be a
small-amplitude Gibbs oscillation as a function of the Q of the
Fourier-transformed local potential with a period of 2*Pi/rad(mmax).

    4
    1  0.0000000000000D+00 -9.3285752771727D-01
    2  2.0000000000000D-02 -9.3284556214441D-01
    3  4.0000000000000D-02 -9.3280966542583D-01
    .....
    .....
    313  6.2400000000000D+00 -3.3275722650983D-01

Fifth (last) data block (read in subroutine psp8cc)

lines labeled 1,2,...,mmax -

*       1st column: radial grid index
*       2nd column: radial grid mesh point (LINEAR mesh starting at 0.0)
*       3rd column: rhoc = 4*Pi*model core charge density
*       4th column: d rhoc / dr
*       5th column: d^2 rhoc / dr^2
*       6th column: d^3 rhoc / dr^3
*       7th column: d^4 rhoc / dr^4

This block is only read if the header variable fchrg >0.0.

     1  0.0000000000000D+00  1.8263604328455D+00 -1.3822038662342D-06 -3.4248477651502D-05  2.2949926081983D-05 -1.3721984594931D+01
     2  2.0000000000000D-02  1.8263581757935D+00 -1.8290946429655D-05 -2.7430630291457D-03 -2.7416158741933D-01 -1.3679142685104D+01
     3  4.0000000000000D-02  1.8263568047539D+00 -1.4618869292620D-04 -1.0954897796014D-02 -5.4658870398836D-01 -1.3549201500559D+01
    .....
    .....

## SPIN-ORBIT

Spin-orbit coupling is present in the file `Psps_for_tests/78_Pt_r.oncvpsp.psp8`.
This pseudopotential treats the 5s and 5p core states as valence.  
More details are given in the spin-orbit portion of the discussion section 

    Header: lines 1-4 as above
    line 5: number of Vanderbilt-Kleinman-Bylander projectors nproj for 
            the scalar-relativistic non-local potential for each
            angular momentum l=0 to l=lmax. The value of nproj for lloc 
            must be 0 (if lloc <= lmax).
    line 6: extension_switch=2 indicating spin-orbit projectors are present
    line 7: number of Vanderbilt-Kleinman-Bylander projectors nprojso for 
            the spin-orbit non-local potential for each angular momentum 
            l=1 to l=lmax.

    Pt    NOPTPSP  r_core=  2.21  2.52  2.40  3.02
         78.0000     18.0000      140202    zatom,zion,pspd
         8     2     3     4   500     0    pspcod,pspxc,lmax,lloc,mmax,r2well
      4.99000000  0.00000000  0.00000000    rchrg fchrg qchrg
         2     4     3     3    nproj
         2                 extension_switch
         4     4     3    nprojso

    Data blocks 1-4 are as-above for the scalar-relativistic non-local projector
    coefficients ("ekb"), grid indices, grid points, and projectors.

    Data block 5 is the local potential lloc=4 as above.

    Data blocks 6-8 are as-above for the spin-orbit non-local projector
    coefficients ("ekb"), grid indices, grid points, and projectors.

    Non-linear core corrections are not used in this example, but if they
    were, the model core charge and its derivatives would be the last data block

## DISCUSSION

This introduces a norm-conserving pseudopotential input file format
designed by D. R. Hamann which offers additional flexibility over
previously available Abinit psp formats, as well as possible improvements
in performance in certain respects. Its new features are:

- Numerical psp's, which offer the ability to experiment beyond existing
  tabulated collections.  Available psp's given as basis function sums 
  are easily convertible to this format.

- Linear radial grids, which have the following advantages:

    (a) Computational efficiency compared to log grids because the
        core region of those grids is highly redundant for psp's.

    (b) Increased accuracy at large ecut values.  The spacing of
        mesh points for log grids at larger radii, where there are
	usually still significant contributions to the psp's, becomes
	too large and introduces aliasing noise into the Fourier-
	transform psp's.  In fact, this code is written to double,
	triple, etc. the linear grid and interpolate prior to performing
	the transform to ensure accuracy, based on ecut.

    (c) Storage economy (probably not important in modern computers, but
        the files are easier to read).

    (d) While redundant, repeating the grid in the data block with each
        function facilitates inspection and plotting of the various inputs.

- Allowance for any number of nonlocal projectors for each anugular
momentum channel, previously available for several formats but without
the current flexibility.

- Direct input of the Bloechl-Kleinman-Bylander projectors rather
  than wave functions and semi-local potenials.  This allows greater
  flexibility important for some exploratory applications such
  as self-interaction-corrected psp's.  The projectors are applied
  to the wave functions as the following operators
  ```
                   Sum    |fkb(n,l),m> ekb(n,l) <fkb(n,l),m|
                   (nlm)
  ```
  The conventional Kleinman-Bylander [2] choice for a single projector
  (single n) would be
  ```
                   fkb(r,l) = (V(r,l) - V_loc(r))u(r,l),
  ```
  where V(r,l) is the semi-local norm-conserving pseudopotential for
  angular momentum l, V_loc(r) is the local potential, and u(r,l) is
  the radial Schroedinger equation solution (bound or scattering).[1,2]
  Note that u(r) has the conventional `r` factor, u(r)=r psi(r).
  In this case, the "energies" ekb (which have dimensions 1/energy
  for the above definition of fkb) are
  ```
                   ekb(l) = 1 / [Integral fkb(r,l)*u(r,l) dr],
                                  (0,inf)
  ```
  (no 4*Pi factor).  For the Bloechl generalization, see [[cite:Bloechl1990]].
  Other possible generalizations include the treatment of self-
  interaction corrections as additional projectors, see [[cite:Vogel1995]].</br>
  Finally, creating numerical projectors as part of the psp generation
  process (rather than within Abinit) allows one to solve the radial
  Schrodinger equation in its fully nonlocal form and compare logarithmic
  derivatives of the pseudo and all-electron wave functions over a wide
  energy range. This is a superior test of transferability, and will
  locate shallow "ghost resonances" above the reference valence
  states, which would be missed in the standard ghost tests 
  [[cite:Gonze1990a]], [[cite:Gonze1991]], but
  can lead to poor results.

- Ability to use an arbitrary local potential rather than that for
  one particular angular momentum channel.  This can be helpful in
  eliminating ghost states [[cite:Gonze1990a]], [[cite:Gonze1991]].</br>
  In the example above, a calcium pseudopotential, the "conventional"
  choice of the d-channel (lloc=2) for the local potential is particuarly
  bad because the incipient d shell represented by a shallow d
  resonance leads to a very deep d potential.  With this potential as
  local, it is nearly impossible to avoid ghost states.  Using p as
  the local potential would avoid this, but would place a spurious
  core-orthogonalization repulsive term in channels for higher angular
  momenta (f, g, h, ...) which could lead to errors in, for example,
  lattice constants.  The local potential here was formed by a
  smooth (polynomial) extrapolation of the tail-region all-electron
  potential to the origin.  While some other Abinit psp formats
  include the ability to form the local potential as a simple average
  over semi-local potentials, the d potential would likely still cause
  problems in this case.

- Spherical bessel functions for Fourier transforms in the `psp8*.f`
  routines are computed by a highly accurate recursion method in sbf8.f,
  and all angular momentum channels are treated simultaneously,
  increasing accuracy and efficiency compared to sin/cos formulations,
  especially for large and small arguments.

- The ability to experiment with various model core charges for the
  non-linear core correction [[cite:Louie1982]] is advantageous.  Functions that are
  not sufficiently smooth can lead to potentially severe convergence
  errors, especially in response function calculations.  For other
  formats, Abinit internally generates higher derivatives of the model
  core charge up to the 4th derivative, and piecewise-constructed
  model functions that LOOK smooth can cause havoc.  The requirement
  that all the derivatives to be present in the input file permits
  early detection of such problems, and aids their remediation.

- Release 7.x introduces the ablity to include spin-orbit coupling.  The
  Pt pseudopotential used as an example above was initially generated
  with two non-local projectors each from Dirac wave functions with j=l+1/2 
  and j=l-1/2 using Vanderbilt generalized norm conservation, which
  gives much greater accuracy than the Bloechl construction discussed
  above.[[cite:Hamann2013]]  Weighted sums and differences of these non-local operators
  were separately re-diagonalized to produce efficient orthonormal 
  projectors for the scalar-relativistic and spin-orbit terms which are
  required internally by abinit.  For a detailed description of the
  projectors, see Ref. 7.  The operation of SR and SO terms on the
  spinor wave functions is a straightforward generalization of
  ```
                   Sum    |fkb(n,l),m> ekb(n,l) <fkb(n,l),m|
                   (nlm)
  ```
  where "fkb" and "ekb" are appropriately reinterpreted.  Since the
  "fkb" for the Pt example are orthonormal, the "ekb" have the
  dimensions energy.  Projectors with negligibly small energies
  (<2E-5 Ha at present) are neglected.  There is full flexibility 
  to use other methods to generate the SR and SO terms.
