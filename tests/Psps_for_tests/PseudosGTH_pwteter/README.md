# Information about the pseudopotentials in this directory

These pseudopotentials have the same characteristics:

* Same pseudopotential generator (GTH, see Phys. Rev. B 54, 1703 (1996))

* Same pseudopotential format (pspcod 2)

* Same exchange-correlation functional (pspxc 1 , Perdew-Wang parameterized by Teter)

They can also be found at the page https://www.abinit.org/sites/default/files/PrevAtomicData/psp-links/hgh.html
WARNING : BigDFT try to read an additional line giving rcutoff and rloc, present in some pseudopotentials.
If such a line exist after the value defined by lmax, but does NOT contain proper rcutoff and rloc (e.g. 0.0), BigDFT might fail by hanging
or SEGFAULT, without giving proper error message. This is difficult to debug.

