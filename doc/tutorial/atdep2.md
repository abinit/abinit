---
authors: FB & JB
---

# Second tutorial on aTDEP

![INP SIMaP](https://simap.grenoble-inp.fr/medias/photo/thumbnail4-b_1521620102788-jpg)

## The 3$^{rd}$ and 4$^{th}$ order **effective** Interatomic Force Constants (IFC)

This tutorial shows how to build an anharmonic **Temperature Dependent Effective Potential** (TDEP) using the ABINIT package

In practice, this means to obtain the $3^{rd}$ and 4$^{th}$ order effective IFC.  Many quantities (Gr\"uneisen parameter, thermal expansion...) can be derived therefrom.

You will learn:

1. how to launch aTDEP just after an ABINIT simulation, 
2. the meaning and effects of the main input variables, and 
3. how to exploit the data provided in the output files.

You are supposed to have performed the [1$^{st}$ aTDEP tutorial](/tutorial/atdep1) and strongly encouraged to read the following documents:

* User guide: [[pdf:aTDEP_Guide| aTDEP guide]]  
* Theory: [[pdf:aTDEP_Paper|aTDEP paper]] corresponding to the article [[cite:Bottin2020]]

This tutorial should take about 1.5 hour.

[TUTORIAL_README]

## 1. Summary of the aTDEP method

In the previous tutorial, we have considered that the potential energy of a crystal can be rewritten using a Taylor expansion around the equilibrium. If this expansion is truncated at the 4$^{th}$ order, we obtain:
$$
U=   U_0 + 
\sum_{p\ge 1}^4 \frac{1}{p\ !} \sum_{\substack{\alpha_1...\alpha_p \\ i_1...i_p}}\overset{(p)}{\Theta}\vphantom{\Theta}_{i_1...i_p}^{\alpha_1...\alpha_p}\prod_{k=1}^p
 u_{i_k}^{\alpha_k}
$$

In the same way as previously, it is then possible to obtain the 3$^{rd}$ and 4$^{th}$ order **effective** IFC $\overset{(3)}{\Theta}\vphantom{\Theta}$ and $\overset{(4)}{\Theta}\vphantom{\Theta}$ by using a least squares method. These **effective** IFC are no longer constant and become temperature dependent by taking into account in an **effective** way all the terms neglected (above the truncation). The anharmonicity comes from both the presence of 3$^{rd}$ and 4$^{th}$ order **effective** IFC and their temperature dependancy.

## 2. Negative thermal expansion : Si-d

This calculation is similar to the one performed in the following article [[cite:Bottin2020]].

###	2.0 NetCDF input files

###	2.1 Convergence with respect to Rcut

###	2.2 Etot/FcatMDvsTDEP

###	2.3 The mode Gr\"uneisen parameters and the thermal expansion

###	2.4 Another question?

## 3. Temperature effect on the Gr\"uneisen parameters : MgO.

This calculation is similar to the one performed in the following article [[cite:Bouchet2019]].

###	3.1 Convergence with respect to Rcut

###	3.2 The LO-TO splitting

###	3.3 Another question?

## 4. Description of an alloy : UMo-$\gamma$

This calculation is similar to the one performed in the following article [[cite:Castellano2020]].

###	4.1 Convergence with respect to Rcut

###	4.2 The VCA and the SIFC approaches

###	4.3 Another question?

## 5. HOWTO use `qtdep` ?

![Harmonic Potential](https://upload.wikimedia.org/wikipedia/commons/6/60/Potential_approximation.png)
