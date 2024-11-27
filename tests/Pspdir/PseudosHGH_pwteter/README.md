# Information about the pseudopotentials in this directory

These pseudopotentials have the same characteristics:

* Same pseudopotential generator (HGH, see Phys. Rev. B 58, 3641 (1998))

* Same pseudopotential format (pspcod 3)

* Same exchange-correlation functional (pspxc 1 , Perdew-Wang parameterized by Teter)

They can also be found at the page https://www.abinit.org/sites/default/files/PrevAtomicData/psp-links/hgh.html
WARNING : BigDFT try to read an additional line giving rcutoff and rloc, present in some pseudopotentials.
If such a line exist after the value defined by lmax, but does NOT contain proper rcutoff and rloc (e.g. 0.0), BigDFT might fail by hanging
or SEGFAULT, without giving proper error message. This is difficult to debug.

