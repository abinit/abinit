---
authors: XG
---

# ABINIT package: context of development

The Corning code for electronic structure calculations began
to be developed in the late eighties by D.C. Allan, Corning Inc.
The fundamental algorithm of this code, the band-by-band
conjugate gradient algorithm, was proposed by M.P. Teter, M. Payne
and D.C.Allan, and described thoroughly in the review paper [[cite:Payne1992]]
where other references can be found.
It was written in Fortran 77. Some technical features
were: a plane-wave representation of wavefunctions, use of
pseudopotentials, local-density approximation (LDA) within
the Density-Functional Theory (DFT).

In July 1990, X. Gonze joined Cornell University, where
M.P.Teter was professor. Building upon the Corning code, the Respfn
(Response Function code) was written. See [[cite:Gonze1997a]]
for a complete description of the algorithms. While the Corning and
Respfn code were consistent at the end of X. Gonze stay
in Cornell University (in september 1992),
they began to diverge afterwards. Corning Inc. had
agreed to allow further development based on the version 920813
in Louvain-la-Neuve. Some features of Respfn were
developed by J.-C. Charlier and C. Lee both before the separation
of codes as well as after it.

In 1993 Corning Inc. agreed with Biosym Inc. that the Corning code,
renamed *Plane_Wave*, would be commercialized by Biosym Inc.
In 1995, Biosym Inc. merged with MSI (Molecular Simulation Inc.),
and shortly afterwards, the decision was taken by MSI not to
continue the development of the *Plane_Wave* code.

In 1996, D.C. Allan and X. Gonze explored the possibility to write
a new code, that should not be commercialized, but rendered
available to the scientific community as a freeware. It would be written
in Fortran 90, would include parallel features, would be based
on a SCF (Self-Consistency Field) algorithm. Its writing would
be facilitated by the availability of the
Respfn codes and Corning 920813 codes, as well as later
developments of the Corning code (then renamed *Plane_Wave*),
excluding the features developed by Biosym computer scientists,
or tightly bound with Biosym proprietary features.

The project was first named DFT2000, but the name ABINIT was definitely
adopted in September 1998. Version 1.9 of ABINIT was made available
in March 1999 to beta testers (for some it was done a bit earlier),
outside Louvain-la-Neuve or Corning. It had only "ground-state" features.

Version 2.0 was released in July 1999. The possibility to compute
response-functions (phonons, Born effective charges, dielectric constant)
was available.

Version 3.0, the first version to be proposed
on the Web under GNU GPL licence,
has features that considerably improve upon the existing
Respfn and the 920813 version of Corning  (or even the latest available
version of *Plane_Wave*), i.e. much faster, better parallelisation ...
ABINIT v3.0 was made available in December 2000.

An advisory commitee was set up in June 2000, and served
for the whole lifetime of the version 3.

Version 4.0 was delivered after the first international ABINIT
developer workshop, in January 2003, with a renewed advisory committee,
and the aim to bridge the speed gap with respect to other plane-wave
based first-principles codes: PAW, parallel FFT, better geometry
optimization, better molecular dynamics.

The second international ABINIT workshop, in May 2005 was held in Paris.

Version 5.0 was delivered in Autumn 2005, with a completely renewed
build system.
