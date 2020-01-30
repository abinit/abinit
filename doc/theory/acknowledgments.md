---
description: Suggested acknowledgments and references
---

# Acknowledgments  

This file provides a description of suggested acknowledgments and references
to be inserted in scientific papers whose results have been obtained thanks to
ABINIT. It discusses also briefly the problem of co-authorship.

## Introduction

In the next section, you will find several references we suggest you to cite in your papers that have benefited
from ABINIT. However, we wish first to clarify the spirit in which the present document has been written.
The users of the code have no formal obligations with respect to the ABINIT group (within the limits of the
GNU General Public License). However, it is common practice in the scientific
literature, to acknowledge the efforts of people that have made the research possible.
Please note the following:

1. The ABINIT project, in order to be viable, should be known as a robust tool, 
that has been tested, and that has allowed good scientific research.
This will be facilitated if the ABINIT project is mentioned properly in research papers. 

2. Some recent ideas and algorithms are coded, and it would be fair to cite these.

3. You might also register on the [ABINIT forum](https://forum.abinit.org). 
Indirectly, this also helps the ABINIT developer group, because the total number of registered people 
on the forum is often cited as an indicator of the user community size. 
In agreement with the GNU General Public License, there is no request for co-authorship 
of articles whose scientific results have been obtained thanks to ABINIT, by any ABINIT developer. 
This applies even for recently implemented features, as their availability in a public version 
is governed by the GNU GPL license.
If you think your work could benefit from collaboration with ABINIT developers, 
you can contact the ABINIT group for a possible arrangement, in which case co-authorship should be discussed. 
(Of course, the ABINIT developers also have the right to decline giving assistance to users).

## List of suggestions

The 2020 general ABINIT article [Gen.1](#g1),
first in the list of suggestions below, should be cited in papers that have benefited from the
ABINIT project, irrespective of their content.
There are four other ABINIT articles, [Gen.2](#g2), [Gen.3](#g3), [Gen.4](#b4), and [Gen.5](#g5) that might as well be
considered, irrespective of ther content, because these publications are
quite general as well, although they are older (2016, 2009, 2005 and 2002).  

There are also many articles that are more focused: they describe some
specific capability of ABINIT. The following list is actually not complete...
More references will be proposed by ABINIT itself (see the end of the output
file), see the [[varset:allvars|database]] of information on the input variables, 
in the [[topic:index|topics files]]. as well as in the references of the five general papers.

<a id="g1"></a>
- **Gen.1** At least, the most recent article [[cite:Gonze2020]] that describes the ABINIT project
should be mentioned in the bibliography section of your paper.
A version of this paper, that is not formatted for Computer Phys. Comm. is available
[here](https://www.abinit.org/sites/default/files/ABINIT20.pdf).
The licence allows the authors to put it on the Web.

<a id="g2"></a>
- **Gen.2** A version of the 2016 ABINIT article [[cite:Gonze2016]],
that is not formatted for Computer Phys. Comm. is available 
[here](https://www.abinit.org/sites/default/files/ABINIT16.pdf).
The licence allows the authors to put it on the Web. 

<a id="g3"></a>
- **Gen.3**. The 2009 ABINIT article [[cite:Gonze2009]] is especially detailed. A version of this paper, 
that is not formatted for Computer Phys. Comm. is available 
[here](https://www.abinit.org/sites/default/files/about/ABINIT_CPC_v10.pdf). 
The licence allows the authors to put it on the Web. 

<a id="g4"></a>
- **Gen.4**. The 2005 ABINIT article [[cite:Gonze2005]] is quite short. 
The .pdf of the latter paper is available [here](https://www.abinit.org/sites/default/files/zfk_0505-06_558-562.pdf). 
Note that it should not redistributed (Copyright by Oldenburg Wissenshaftverlag, 
the licence allows the authors to put it on the Web).

<a id="g5"></a>
- **Gen.5**. The very first paper on the ABINIT project [[cite:Gonze2002]] might also be considered for citation.

<a id="s1"></a>
- **Spe.1.** If the Projector-Augmented Wave method as implemented in ABINIT is used [[cite:Torrent2008]] should be mentioned.

<a id="s2"></a>
- **Spe.2.** Many ingredients needed for the calculations of responses to atomic displacements 
or homogeneous electric fields (dynamical matrices, effective charges and dielectric constants), 
as well as the Fourier interpolation implemented in the 'anaddb' code are described in [[cite:Gonze1997]] and [[cite:Gonze1997a]]. 

<a id="s3"></a>
- **Spe.3.** The methods used for the calculation of responses to homogeneous strain 
(elastic tensors, piezoelectric tensors, and internal force-response tensors) are described in [[cite:Hamann2005]].

<a id="s4"></a>
- **Spe.4.** If the "static" non-linear capabilities of ABINIT are used (Raman efficiencies, electro-optic coefficients ... ), 
cite [[cite:Veithen2005]] 

<a id="s5"></a>
- **Spe.5.** If the integration over the phonon degrees of freedom is used ([[anaddb:thmflag]]), cite [[cite:Lee1995]]. 

<a id="s6"></a>
- **Spe.6.** If the self-consistent capabilities of ABINIT beyond DFT are used (GW, COHSEX, HF, etc), 
cite [[cite:Bruneval2006]].

<a id="s7"></a>
- **Spe.7.** If the completeness relationship is used to speed up the convergence with respect to the number 
of bands in a GW calculation (input variables [[gwcomp]] and [[gwencomp]]), cite [[cite:Bruneval2008]] 

<a id="s8"></a>
- **Spe.8.** If the massive parallelism of ABINIT (coupled band/FFT or even coupled band/FFT/k-points ) is used, 
cite [[cite:Bottin2008]] ([available on Arxiv.org](https://arxiv.org/abs/0707.3405)). 

<a id="s9"></a>
- **Spe.9.** If the LDA+U method as implemented in ABINIT is used, cite [[cite:Amadon2008]]. 

<a id="s10"></a>
- **Spe.10.** If the ONCVPSP pseudopotentials are used, cite [[cite:Hamann2013]]. 

<a id="s11"></a>
- **Spe.11.** If the Van Der Waals DFT-D (Grimme) functionals are used, cite [[cite:Vantroeye2016]]. 

<a id="s12"></a>
- **Spe.12.** If the temperature dependence of the electronic structure or the zero-point motion effect 
on the electronic structure are computed, cite [[cite:Ponce2014]] and [[cite:Ponce2015]].

<a id="s13"></a>
- **Spe.13.** If the direct (DFPT) computation of effective masses is used, see [[topic:EffectiveMass]], cite [[cite:Laflamme2016]].
