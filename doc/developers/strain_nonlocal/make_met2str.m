(*
   Mathematica program to create code for the program metstr.f
   for the strain derivative of the nonlocal pseudopotential
   operator.  The operator for angular momentum l has the form

     O_l = V_l(kp*kp) P_l(kp*kp, kp*k, k*k) V_l(k*k)

   where k and kp are the wavevectors of the input and output
   wavefunctions.  and P_l are Legendre polynomials.  The "*" sign
   represents a dot product.  The V_l and d V_l/d(k*k) are multiplied
   in a separate routine and are not explicit in this code.  The three
   sections of Mathematica code corresponding to iterm = 1, 2, 3 represent
   strain (s) derivatives of the 3 terms in the O_l product.  (The angular
   momentum "l" is represented by "rank" following Abinit convention in
   the routine nonlop which calls metstr.)

   For the straight application of O_l to a wavefunction, the rank of
   the input k tensor equals the rank of the output kp tensor.  This is
   also true for the strain derivative of P_l (iterm = 2).  However,
   the strain derivatives of the V_l terms  introduce two additional
   powers of k or kp, and thus couple to rank+2 tensors on either
   input or output.
*)

MakeMet2str[iterm_, outfile_] :=

Module[{tnk,tnkp,k,kp,dt,ks,kps,Plegendre,rank,
  rankin,rankout,limitin,limitout,poly,term,ii,jj},


(* symmetric tensors of k and kp components, ranks 1 through 7
   as 2-dimensional Mathematica arrays tnk[[rank]][[element]] and tnkp
   number of elements is (rank+1)(rank+2)/2                    *)

tnk = {
{1},  (* 0, 1 *)

{k1,  (* 1, 1 *)
 k2,  (* 1, 2 *)
 k3}, (* 1, 3 *)

{k1 k1,  (* 2, 1 *)
 k2 k2,  (* 2, 2 *)
 k3 k3,  (* 2, 3 *)
 k3 k2,  (* 2, 4 *)
 k3 k1,  (* 2, 5 *)
 k2 k1}, (* 2, 6 *)

{k1 k1 k1,  (* 3, 1 *)
 k2 k2 k1,  (* 3, 2 *)
 k3 k3 k1,  (* 3, 3 *)
 k3 k2 k1,  (* 3, 4 *)
 k3 k1 k1,  (* 3, 5 *)
 k2 k1 k1,  (* 3, 6 *)
 k2 k2 k2,  (* 3, 7 *)
 k3 k3 k2,  (* 3, 8 *)
 k3 k2 k2,  (* 3, 9 *)
 k3 k3 k3}, (* 3, 10*)

{k1 k1 k1 k1,  (* 4, 1 *)
 k2 k2 k1 k1,  (* 4, 2 *)
 k3 k3 k1 k1,  (* 4, 3 *)
 k3 k2 k1 k1,  (* 4, 4 *)
 k3 k1 k1 k1,  (* 4, 5 *)
 k2 k1 k1 k1,  (* 4, 6 *)
 k2 k2 k2 k1,  (* 4, 7 *)
 k3 k3 k2 k1,  (* 4, 8 *)
 k3 k2 k2 k1,  (* 4, 9 *)
 k3 k3 k3 k1,  (* 4, 10 *)
 k2 k2 k2 k2,  (* 4, 11 *)
 k3 k3 k2 k2,  (* 4, 12 *)
 k3 k2 k2 k2,  (* 4, 13 *)
 k3 k3 k3 k2,  (* 4, 14 *)
 k3 k3 k3 k3}, (* 4, 15 *)


{k1 k1 k1 k1 k1,  (* 5, 1 *)
 k2 k2 k1 k1 k1,  (* 5, 2 *)
 k3 k3 k1 k1 k1,  (* 5, 3 *)
 k3 k2 k1 k1 k1,  (* 5, 4 *)
 k3 k1 k1 k1 k1,  (* 5, 5 *)
 k2 k1 k1 k1 k1,  (* 5, 6 *)
 k2 k2 k2 k1 k1,  (* 5, 7 *)
 k3 k3 k2 k1 k1,  (* 5, 8 *)
 k3 k2 k2 k1 k1,  (* 5, 9 *)
 k3 k3 k3 k1 k1,  (* 5, 10 *)
 k2 k2 k2 k2 k1,  (* 5, 11 *)
 k3 k3 k2 k2 k1,  (* 5, 12 *)
 k3 k2 k2 k2 k1,  (* 5, 13 *)
 k3 k3 k3 k2 k1,  (* 5, 14 *)
 k3 k3 k3 k3 k1,  (* 5, 15 *)
 k2 k2 k2 k2 k2,  (* 5, 16 *)
 k3 k3 k2 k2 k2,  (* 5, 17 *)
 k3 k2 k2 k2 k2,  (* 5, 18 *)
 k3 k3 k3 k2 k2,  (* 5, 19 *)
 k3 k3 k3 k3 k2,  (* 5, 20 *)
 k3 k3 k3 k3 k3}, (* 5, 21 *)


{k1 k1 k1 k1 k1 k1,  (* 6, 1 *)
 k2 k2 k1 k1 k1 k1,  (* 6, 2 *)
 k3 k3 k1 k1 k1 k1,  (* 6, 3 *)
 k3 k2 k1 k1 k1 k1,  (* 6, 4 *)
 k3 k1 k1 k1 k1 k1,  (* 6, 5 *)
 k2 k1 k1 k1 k1 k1,  (* 6, 6 *)
 k2 k2 k2 k1 k1 k1,  (* 6, 7 *)
 k3 k3 k2 k1 k1 k1,  (* 6, 8 *)
 k3 k2 k2 k1 k1 k1,  (* 6, 9 *)
 k3 k3 k3 k1 k1 k1,  (* 6, 10 *)
 k2 k2 k2 k2 k1 k1,  (* 6, 11 *)
 k3 k3 k2 k2 k1 k1,  (* 6, 12 *)
 k3 k2 k2 k2 k1 k1,  (* 6, 13 *)
 k3 k3 k3 k2 k1 k1,  (* 6, 14 *)
 k3 k3 k3 k3 k1 k1,  (* 6, 15 *)
 k2 k2 k2 k2 k2 k1,  (* 6, 16 *)
 k3 k3 k2 k2 k2 k1,  (* 6, 17 *)
 k3 k2 k2 k2 k2 k1,  (* 6, 18 *)
 k3 k3 k3 k2 k2 k1,  (* 6, 19 *)
 k3 k3 k3 k3 k2 k1,  (* 6, 20 *)
 k3 k3 k3 k3 k3 k1,  (* 6, 21 *)
 k2 k2 k2 k2 k2 k2,  (* 6, 22 *)
 k3 k3 k2 k2 k2 k2,  (* 6, 23 *)
 k3 k2 k2 k2 k2 k2,  (* 6, 24 *)
 k3 k3 k3 k2 k2 k2,  (* 6, 25 *)
 k3 k3 k3 k3 k2 k2,  (* 6, 26 *)
 k3 k3 k3 k3 k3 k2,  (* 6, 27 *)
 k3 k3 k3 k3 k3 k3}, (* 6, 28 *)


{k1 k1 k1 k1 k1 k1 k1,  (* 7, 1 *)
 k2 k2 k1 k1 k1 k1 k1,  (* 7, 2 *)
 k3 k3 k1 k1 k1 k1 k1,  (* 7, 3 *)
 k3 k2 k1 k1 k1 k1 k1,  (* 7, 4 *)
 k3 k1 k1 k1 k1 k1 k1,  (* 7, 5 *)
 k2 k1 k1 k1 k1 k1 k1,  (* 7, 6 *)
 k2 k2 k2 k1 k1 k1 k1,  (* 7, 7 *)
 k3 k3 k2 k1 k1 k1 k1,  (* 7, 8 *)
 k3 k2 k2 k1 k1 k1 k1,  (* 7, 9 *)
 k3 k3 k3 k1 k1 k1 k1,  (* 7, 10 *)
 k2 k2 k2 k2 k1 k1 k1,  (* 7, 11 *)
 k3 k3 k2 k2 k1 k1 k1,  (* 7, 12 *)
 k3 k2 k2 k2 k1 k1 k1,  (* 7, 13 *)
 k3 k3 k3 k2 k1 k1 k1,  (* 7, 14 *)
 k3 k3 k3 k3 k1 k1 k1,  (* 7, 15 *)
 k2 k2 k2 k2 k2 k1 k1,  (* 7, 16 *)
 k3 k3 k2 k2 k2 k1 k1,  (* 7, 17 *)
 k3 k2 k2 k2 k2 k1 k1,  (* 7, 18 *)
 k3 k3 k3 k2 k2 k1 k1,  (* 7, 19 *)
 k3 k3 k3 k3 k2 k1 k1,  (* 7, 20 *)
 k3 k3 k3 k3 k3 k1 k1,  (* 7, 21 *)
 k2 k2 k2 k2 k2 k2 k1,  (* 7, 22 *)
 k3 k3 k2 k2 k2 k2 k1,  (* 7, 23 *)
 k3 k2 k2 k2 k2 k2 k1,  (* 7, 24 *)
 k3 k3 k3 k2 k2 k2 k1,  (* 7, 25 *)
 k3 k3 k3 k3 k2 k2 k1,  (* 7, 26 *)
 k3 k3 k3 k3 k3 k2 k1,  (* 7, 27 *)
 k3 k3 k3 k3 k3 k3 k1,  (* 7, 28 *)
 k2 k2 k2 k2 k2 k2 k2,  (* 7, 29 *)
 k3 k3 k2 k2 k2 k2 k2,  (* 7, 30 *)
 k3 k2 k2 k2 k2 k2 k2,  (* 7, 31 *)
 k3 k3 k3 k2 k2 k2 k2,  (* 7, 32 *)
 k3 k3 k3 k3 k2 k2 k2,  (* 7, 33 *)
 k3 k3 k3 k3 k3 k2 k2,  (* 7, 34 *)
 k3 k3 k3 k3 k3 k3 k2,  (* 7, 25 *)
 k3 k3 k3 k3 k3 k3 k3}};(* 7, 36 *)


tnkp = {
{1},  (* 0, 1 *)

{kp1,  (* 1, 1 *)
 kp2,  (* 1, 2 *)
 kp3}, (* 1, 3 *)

{kp1 kp1,  (* 2, 1 *)
 kp2 kp2,  (* 2, 2 *)
 kp3 kp3,  (* 2, 3 *)
 kp3 kp2,  (* 2, 4 *)
 kp3 kp1,  (* 2, 5 *)
 kp2 kp1}, (* 2, 6 *)

{kp1 kp1 kp1,  (* 3, 1 *)
 kp2 kp2 kp1,  (* 3, 2 *)
 kp3 kp3 kp1,  (* 3, 3 *)
 kp3 kp2 kp1,  (* 3, 4 *)
 kp3 kp1 kp1,  (* 3, 5 *)
 kp2 kp1 kp1,  (* 3, 6 *)
 kp2 kp2 kp2,  (* 3, 7 *)
 kp3 kp3 kp2,  (* 3, 8 *)
 kp3 kp2 kp2,  (* 3, 9 *)
 kp3 kp3 kp3}, (* 3, 10*)

{kp1 kp1 kp1 kp1,  (* 4, 1 *)
 kp2 kp2 kp1 kp1,  (* 4, 2 *)
 kp3 kp3 kp1 kp1,  (* 4, 3 *)
 kp3 kp2 kp1 kp1,  (* 4, 4 *)
 kp3 kp1 kp1 kp1,  (* 4, 5 *)
 kp2 kp1 kp1 kp1,  (* 4, 6 *)
 kp2 kp2 kp2 kp1,  (* 4, 7 *)
 kp3 kp3 kp2 kp1,  (* 4, 8 *)
 kp3 kp2 kp2 kp1,  (* 4, 9 *)
 kp3 kp3 kp3 kp1,  (* 4, 10 *)
 kp2 kp2 kp2 kp2,  (* 4, 11 *)
 kp3 kp3 kp2 kp2,  (* 4, 12 *)
 kp3 kp2 kp2 kp2,  (* 4, 13 *)
 kp3 kp3 kp3 kp2,  (* 4, 14 *)
 kp3 kp3 kp3 kp3}, (* 4, 15 *)


{kp1 kp1 kp1 kp1 kp1,  (* 5, 1 *)
 kp2 kp2 kp1 kp1 kp1,  (* 5, 2 *)
 kp3 kp3 kp1 kp1 kp1,  (* 5, 3 *)
 kp3 kp2 kp1 kp1 kp1,  (* 5, 4 *)
 kp3 kp1 kp1 kp1 kp1,  (* 5, 5 *)
 kp2 kp1 kp1 kp1 kp1,  (* 5, 6 *)
 kp2 kp2 kp2 kp1 kp1,  (* 5, 7 *)
 kp3 kp3 kp2 kp1 kp1,  (* 5, 8 *)
 kp3 kp2 kp2 kp1 kp1,  (* 5, 9 *)
 kp3 kp3 kp3 kp1 kp1,  (* 5, 10 *)
 kp2 kp2 kp2 kp2 kp1,  (* 5, 11 *)
 kp3 kp3 kp2 kp2 kp1,  (* 5, 12 *)
 kp3 kp2 kp2 kp2 kp1,  (* 5, 13 *)
 kp3 kp3 kp3 kp2 kp1,  (* 5, 14 *)
 kp3 kp3 kp3 kp3 kp1,  (* 5, 15 *)
 kp2 kp2 kp2 kp2 kp2,  (* 5, 16 *)
 kp3 kp3 kp2 kp2 kp2,  (* 5, 17 *)
 kp3 kp2 kp2 kp2 kp2,  (* 5, 18 *)
 kp3 kp3 kp3 kp2 kp2,  (* 5, 19 *)
 kp3 kp3 kp3 kp3 kp2,  (* 5, 20 *)
 kp3 kp3 kp3 kp3 kp3}, (* 5, 21 *)


{kp1 kp1 kp1 kp1 kp1 kp1,  (* 6, 1 *)
 kp2 kp2 kp1 kp1 kp1 kp1,  (* 6, 2 *)
 kp3 kp3 kp1 kp1 kp1 kp1,  (* 6, 3 *)
 kp3 kp2 kp1 kp1 kp1 kp1,  (* 6, 4 *)
 kp3 kp1 kp1 kp1 kp1 kp1,  (* 6, 5 *)
 kp2 kp1 kp1 kp1 kp1 kp1,  (* 6, 6 *)
 kp2 kp2 kp2 kp1 kp1 kp1,  (* 6, 7 *)
 kp3 kp3 kp2 kp1 kp1 kp1,  (* 6, 8 *)
 kp3 kp2 kp2 kp1 kp1 kp1,  (* 6, 9 *)
 kp3 kp3 kp3 kp1 kp1 kp1,  (* 6, 10 *)
 kp2 kp2 kp2 kp2 kp1 kp1,  (* 6, 11 *)
 kp3 kp3 kp2 kp2 kp1 kp1,  (* 6, 12 *)
 kp3 kp2 kp2 kp2 kp1 kp1,  (* 6, 13 *)
 kp3 kp3 kp3 kp2 kp1 kp1,  (* 6, 14 *)
 kp3 kp3 kp3 kp3 kp1 kp1,  (* 6, 15 *)
 kp2 kp2 kp2 kp2 kp2 kp1,  (* 6, 16 *)
 kp3 kp3 kp2 kp2 kp2 kp1,  (* 6, 17 *)
 kp3 kp2 kp2 kp2 kp2 kp1,  (* 6, 18 *)
 kp3 kp3 kp3 kp2 kp2 kp1,  (* 6, 19 *)
 kp3 kp3 kp3 kp3 kp2 kp1,  (* 6, 20 *)
 kp3 kp3 kp3 kp3 kp3 kp1,  (* 6, 21 *)
 kp2 kp2 kp2 kp2 kp2 kp2,  (* 6, 22 *)
 kp3 kp3 kp2 kp2 kp2 kp2,  (* 6, 23 *)
 kp3 kp2 kp2 kp2 kp2 kp2,  (* 6, 24 *)
 kp3 kp3 kp3 kp2 kp2 kp2,  (* 6, 25 *)
 kp3 kp3 kp3 kp3 kp2 kp2,  (* 6, 26 *)
 kp3 kp3 kp3 kp3 kp3 kp2,  (* 6, 27 *)
 kp3 kp3 kp3 kp3 kp3 kp3}, (* 6, 28 *)


{kp1 kp1 kp1 kp1 kp1 kp1 kp1,  (* 7, 1 *)
 kp2 kp2 kp1 kp1 kp1 kp1 kp1,  (* 7, 2 *)
 kp3 kp3 kp1 kp1 kp1 kp1 kp1,  (* 7, 3 *)
 kp3 kp2 kp1 kp1 kp1 kp1 kp1,  (* 7, 4 *)
 kp3 kp1 kp1 kp1 kp1 kp1 kp1,  (* 7, 5 *)
 kp2 kp1 kp1 kp1 kp1 kp1 kp1,  (* 7, 6 *)
 kp2 kp2 kp2 kp1 kp1 kp1 kp1,  (* 7, 7 *)
 kp3 kp3 kp2 kp1 kp1 kp1 kp1,  (* 7, 8 *)
 kp3 kp2 kp2 kp1 kp1 kp1 kp1,  (* 7, 9 *)
 kp3 kp3 kp3 kp1 kp1 kp1 kp1,  (* 7, 10 *)
 kp2 kp2 kp2 kp2 kp1 kp1 kp1,  (* 7, 11 *)
 kp3 kp3 kp2 kp2 kp1 kp1 kp1,  (* 7, 12 *)
 kp3 kp2 kp2 kp2 kp1 kp1 kp1,  (* 7, 13 *)
 kp3 kp3 kp3 kp2 kp1 kp1 kp1,  (* 7, 14 *)
 kp3 kp3 kp3 kp3 kp1 kp1 kp1,  (* 7, 15 *)
 kp2 kp2 kp2 kp2 kp2 kp1 kp1,  (* 7, 16 *)
 kp3 kp3 kp2 kp2 kp2 kp1 kp1,  (* 7, 17 *)
 kp3 kp2 kp2 kp2 kp2 kp1 kp1,  (* 7, 18 *)
 kp3 kp3 kp3 kp2 kp2 kp1 kp1,  (* 7, 19 *)
 kp3 kp3 kp3 kp3 kp2 kp1 kp1,  (* 7, 20 *)
 kp3 kp3 kp3 kp3 kp3 kp1 kp1,  (* 7, 21 *)
 kp2 kp2 kp2 kp2 kp2 kp2 kp1,  (* 7, 22 *)
 kp3 kp3 kp2 kp2 kp2 kp2 kp1,  (* 7, 23 *)
 kp3 kp2 kp2 kp2 kp2 kp2 kp1,  (* 7, 24 *)
 kp3 kp3 kp3 kp2 kp2 kp2 kp1,  (* 7, 25 *)
 kp3 kp3 kp3 kp3 kp2 kp2 kp1,  (* 7, 26 *)
 kp3 kp3 kp3 kp3 kp3 kp2 kp1,  (* 7, 27 *)
 kp3 kp3 kp3 kp3 kp3 kp3 kp1,  (* 7, 28 *)
 kp2 kp2 kp2 kp2 kp2 kp2 kp2,  (* 7, 29 *)
 kp3 kp3 kp2 kp2 kp2 kp2 kp2,  (* 7, 30 *)
 kp3 kp2 kp2 kp2 kp2 kp2 kp2,  (* 7, 31 *)
 kp3 kp3 kp3 kp2 kp2 kp2 kp2,  (* 7, 32 *)
 kp3 kp3 kp3 kp3 kp2 kp2 kp2,  (* 7, 33 *)
 kp3 kp3 kp3 kp3 kp3 kp2 kp2,  (* 7, 34 *)
 kp3 kp3 kp3 kp3 kp3 kp3 kp2,  (* 7, 25 *)
 kp3 kp3 kp3 kp3 kp3 kp3 kp3}};(* 7, 36 *) ;


(* quantities needed to set up wavevector polynomials
   k  = input wavefunction wavevector
   kp = output wavefunction wavevector
   m  = the metric tensor gmet *)

k = {k1, k2, k3};

kp = {kp1, kp2, kp3};

(* symmetry of metric tensor is built in here
   only symmetric strains are to be permitted so derivatives will
   retain this symmetry *)

m = {{m11[s1,s2], m12[s1,s2], m13[s1,s2]}, {m12[s1,s2], m22[s1,s2], m23[s1,s2]},
        {m13[s1,s2], m23[s1,s2], m33[s1,s2]}};

dt = kp.m.k;
ks = k.m.k;
kps = kp.m.kp;

(* The Legendre polynomials P_l are modified so that each term is
   of of order k^l * kp^l.  We only treat l = 1, 2, and 3. *)

Plegendre={1, dt, 1.5 dt^2 - 0.5 kps ks, 2.5 dt^3 - 1.5 kps dt ks};

Do[

poly = {D[kps,s2] D[kps,s1] Plegendre[[rank+1]],
        D[ks, s2] D[ks, s1] Plegendre[[rank+1]],
       (D[D[kps,s2],s1] Plegendre[[rank+1]]
      + D[kps,s1] D[Plegendre[[rank+1]],s2]
      + D[kps,s2] D[Plegendre[[rank+1]],s1]),
       (D[D[ks, s2],s1] Plegendre[[rank+1]]
      + D[ks, s1] D[Plegendre[[rank+1]],s2]
      + D[ks, s2] D[Plegendre[[rank+1]],s1]),
       (D[kps,s1] D[ks,s2] + D[kps,s2] D[ks,s1]) Plegendre[[rank+1]],
        D[D[Plegendre[[rank+1]],s2],s1]};

rankinarr  = {rank,rank+4,rank,rank+2,rank+2,rank};
rankoutarr = {rank+4,rank,rank+2,rank,rank+2,rank};

rankin = rankinarr[[iterm]];
rankout = rankoutarr[[iterm]];

Write[outfile,TextForm["elseif(rank=="],FortranForm[rank],
  TextForm[") then#"]];

limitin  = (rankin+1)(rankin+2)/2;
limitout  = (rankout+1)(rankout+2)/2;

(* this is the only real computational step -- extracting the tensor
components *)

Do[term = Simplify[
      Coefficient[poly[[iterm]], (tnkp[[rankout+1]][[jj]] *
        tnk[[rankin+1]][[ii]])]];

    Write[outfile,
      TextForm["cm("],
      FortranForm[ii],
      TextForm[","],
      FortranForm[jj],
      TextForm[")="],
      FortranForm[term],
      TextForm["#"]],

{jj,1,limitout}, {ii,1,limitin}],

{rank,0,3}];

(* {rank,1,3}]; *)
(*HACK: drop rank = 0 for term 6; result is zero and causes trouble *)

Close[outfile];

]
