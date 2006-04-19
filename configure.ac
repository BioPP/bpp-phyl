AC_INIT(Bio++ PhylLib, 1.1.0, Julien.Dutheil@univ-montp2.fr, Bpp-Phyl)
AC_CONFIG_SRCDIR(Phyl/SubstitutionModel.h)
AM_INIT_AUTOMAKE(Bpp-Phyl,1.1.0)
AC_LANG_CPLUSPLUS
AC_CHECK_LIB(bpp-utils,   main, , [AC_MSG_ERROR([Library bpp-utils not found.])])
AC_CHECK_LIB(bpp-numcalc, main, , [AC_MSG_ERROR([Library bpp-numcalc not found.])])
AC_CHECK_LIB(bpp-seq,     main, , [AC_MSG_ERROR([Library bpp-seq not found.])])
AC_PROG_CC
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LIBTOOL
AC_ARG_WITH(
	[virtual-cov],
	[  --with-virtual-cov      Enable covariant return type with
	                  virtual inheritance (GCC >= 3.4)],
	[$with_virtual_cov=yes]
)
if test "$with_virtual_cov" == "yes"
then
	AC_MSG_WARN([Covariant return type with virtual inheritance enabled.]);
	AC_DEFINE(VIRTUAL_COV, yes)
fi
AC_OUTPUT(Makefile Phyl/Makefile)
