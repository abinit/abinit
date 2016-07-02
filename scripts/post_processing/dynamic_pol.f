      program dynamic_pol
c
c     this small utility computes the dynamic polarisability
c     from a TDDFT calculation output

      implicit real*8(a-h,p-z)
      character*12 nomin,nomout
	character*5 string,check
      integer nomega,nexcit,omegamax


      real*8 eha(40),eev(40),ostr(40)

	omegamax=5.0d0
	nomega=100

      check = '  Osc'
	write(*,*)'Abinit output filename :'
	read(*,*) nomin
      write(*,*)'Polarisability output filename :'
	read(*,*) nomout

	open(5,file=nomin,status='OLD',err=201)
      string='rrrrr'
      do while (string.ne.check)
      read(5,10) string
	write(*,*) string
	enddo
      read(5,*)
      index=1
      
400	read(5,*,err=300)i,eha(index),eev(index),ostr(index)
      write(*,*)i,eha(index),eev(index),ostr(index)
      index=index+1
	goto 400

300   index=index-1
      write(*,*)'singlet excitations :',index

      string='rrrrr'
      do while (string.ne.check)
      read(5,10) string
	enddo
10    format(A5)
      read(5,*)
      
401	read(5,*,err=301)i,eha(index),eev(index),ostr(index)
      write(*,*)i,eha(index),eev(index),ostr(index)
      index=index+1
	goto 401

301   nexcit=index-1
      write(*,*)'singlet+triplet excitations :',nexcit

      open(6,file=nomout,status='UNKNOWN')
	write(6,*)'Omega,  Alpha'
	do iomega=1,nomega
	  om=omegamax*(iomega*1.0d0/(nomega*1.0d0))
        pola=0.0d0
	  do iexcit=1,nexcit
         pola=pola+ostr(iexcit)/(eev(iexcit)*eev(iexcit)-om*om)
	  enddo
      write(6,50)om,pola
	enddo
50    format(D13.5,' ,  ',D13.5)








      goto 500
201   write(*,*)'file ',nomin,' does not exist'
      goto 500
202   write(*,*)'cannot find string   Oscillator '
      goto 500
500   end