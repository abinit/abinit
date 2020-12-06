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

You are supposed to have performed the [1$^{st}$ aTDEP tutorial](atdep1) and strongly encouraged to read the following documents:

* User guide: [[pdf:aTDEP_Guide| aTDEP guide]]  
* Theory: [[pdf:aTDEP_Paper|aTDEP paper]] corresponding to the article [[cite:Bottin2020]]

This tutorial should take about 1.5 hour.

[TUTORIAL_READMEV9]

## 1. Summary of the aTDEP method

In the previous tutorial, we have considered that the potential energy of a crystal can be rewritten using a Taylor expansion around the equilibrium. If this expansion is truncated at the 4$^{th}$ order, we obtain:
$$
U=   U_0 + 
\sum_{p\ge 1}^4 \frac{1}{p\ !} \sum_{\substack{\alpha_1...\alpha_p \\ i_1...i_p}}\overset{(p)}{\Theta}\vphantom{\Theta}_{i_1...i_p}^{\alpha_1...\alpha_p}\prod_{k=1}^p
 u_{i_k}^{\alpha_k}
$$

In the same way as previously, it is then possible to obtain the 3$^{rd}$ and 4$^{th}$ order **effective** IFC $\overset{(3)}{\Theta}\vphantom{\Theta}$ and $\overset{(4)}{\Theta}\vphantom{\Theta}$ by using a least squares method. These **effective** IFC are no longer constant and become temperature dependent by taking into account in an **effective** way all the terms neglected (above the truncation). The anharmonicity comes from both the presence of 3$^{rd}$ and 4$^{th}$ order **effective** IFC and their temperature dependancy.

## 2. Temperature effect on Gr\"uneisen parameter : MgO.

This calculation is similar to the one performed in the following article [[cite:Bouchet2019]].

###	2.1 Convergence with respect to Rcut

###	2.2 Another question?

## 3. A negative thermal expansion : Si-d

This calculation is similar to the one performed in the following article [[cite:Bottin2020]].

###	3.1 Convergence with respect to Rcut

###	3.2 Effect of the number of configurations

###	3.3 Another question?



## 4. Temperature effects on the sound speed : Fe-$\epsilon$

This calculation is similar to the one performed in the following article [[cite:Dewaele2021]].

###	4.1 Convergence with respect to Rcut

###	4.2 Another question?

