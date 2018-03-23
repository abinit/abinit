AC_INIT([Wannier90],[1.2],[yann.pouillon@ehu.es],[wannier90])
AC_REVISION([Fake configure script for Wannier90 1.2])

cat >make.sys <<EOF
F90 = ${F90}
FCOPTS = ${CPPFLAGS} ${F90FLAGS}
LDOPTS = ${CPPFLAGS} ${F90FLAGS} ${LDFLAGS}
LIBS = ${LIBS}

install:
	cp libwannier.a ${prefix}/lib
	cp wannier90.x ${prefix}/bin

EOF
