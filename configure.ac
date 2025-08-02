##
## -------------------------------------------------------------------------
## scitao - Scilab/Scicos toolbox for Adaptive Optics 
## Copyright (C) 2006  scitao
## Some inspiration came from SIP's configure.ac (Ricardo Fabbri)
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
## See the GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software Foundation,
## Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
## -------------------------------------------------------------------------
##

AC_PREREQ(2.57)
AC_INIT([scitao], [0.2.5], [goldollarch@gmail.com])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_SRCDIR([src/scao.c])
AM_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CC
if test x"$ac_cv_prog_cc_stdc" = x'no'; then
   AC_MSG_WARN([[found the $CC compiler but it is not ANSI-C 
                compliant.]])
fi

AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

AC_LIBTOOL_WIN32_DLL
AC_DISABLE_STATIC
AC_PROG_LIBTOOL

# Check to see if building shared libraries
libtool_build_shared_libs='no'
if ./libtool --config | grep 'build_libtool_libs=yes' > /dev/null
then
  libtool_build_shared_libs='yes'
fi

#
# Tests for windows
#
native_win32_build='no'
cygwin_build='no'
case $host_os in
  cygwin*) 
      cygwin_build='yes'
      ;;
  mingw* )
      native_win32_build='yes' 
      ;;
esac

# checks for xsltproc (to buld HTML help pages)
AC_CHECK_PROG([has_xsltproc],[xsltproc],[yes],[no], , )
if test x$has_xsltproc = xno; then
   # checks for sabcmd as an alternative
   AC_CHECK_PROG([has_sablotron],[sabcmd],[yes],[no], , )
   if test x$has_sablotron = xno; then
      AC_MSG_ERROR([[xsltproc or sabcmd programs were not found but \
at least one of them is required to build the scitao manual pages.
Download xslproc (from libxml) at:
   http://www.xmlsoft.org 
or sabcmd (from sablotron) at:
   http://www.gingerall.com/charlie/ga/xml/p_sab.xml]])
   fi
fi

# test if SCI is defined. If not,
# test if scilab is in path, then get SCI from scilab -nw.
if test -z "$SCI"; then
   if test $native_win32_build = yes; then
      dnl TODO take SCIDIR from the registry ??
      AC_MSG_WARN([[hardcoding SCI to /c/scilab]])
      SCIDIR="/c/scilab"
   else
      AC_CHECK_PROG([has_scilab],[scilab],[yes],[no], , )
      if test x$has_scilab = xno; then
         AC_MSG_ERROR([[neither Scilab was found in your PATH nor a SCI variable has been defined]])
      else
         cmd='F=mopen("path.incl","w");
              mfprintf(F,SCI);
              mclose(F);'
         echo "$cmd" | scilab -nw >/dev/null
         SCIDIR=`cat path.incl`
      fi
   fi
else
   SCIDIR="$SCI"
fi
CPPFLAGS="$CPPFLAGS -I${SCIDIR}/routines -I${SCIDIR}/routines/scicos"

AC_PREFIX_DEFAULT([SCI/contrib/scitao])

if test "x$prefix"  = xNONE; then 
   prefix="${SCIDIR}/contrib/scitao"
fi
TOOLBOXDIR=$prefix

# Set the default installation
AC_PREFIX_DEFAULT([$TOOLBOXDIR])
AC_SUBST(TOOLBOXDIR)
AC_SUBST(SCIDIR)

########
#
# Set defines required to build DLLs and modules using MinGW
# thx to Bob F.
#
LIBRARY_EXTRA_CPPFLAGS=''
LIBRARY_EXTRA_LIBS= ''

if test "${native_win32_build}" = 'yes'
then
  if test "${libtool_build_shared_libs}" = 'yes'
  then
    SHARED_LIB_SUFFIX=dll
    CPPFLAGS="$CPPFLAGS -DFORDLL"
    AC_MSG_WARN([[Don't forget to rename SCI/bin/LibScilab* to lowercase]])
    LDFLAGS="$LDFLAGS -L${SCIDIR}/bin"
    LIBRARY_EXTRA_CPPFLAGS="-D_SIMPLELIB"
    LIBRARY_EXTRA_LIBS="-lscilab"
  else
    CPPFLAGS="$CPPFLAGS -D_LIB"
  fi
else
   SHARED_LIB_SUFFIX=so
fi
AC_SUBST(SHARED_LIB_SUFFIX)

AC_SUBST(LIBRARY_EXTRA_CPPFLAGS)
AC_SUBST(LIBRARY_EXTRA_LIBS)

# Checks for libraries.
AC_CHECK_LIB([fftw3], [fftw_execute], [true])
AC_CHECK_LIB([fftw3f], [fftwf_execute], [true])
AC_CHECK_LIB([cfitsio], [ffinit], [true])
AC_CHECK_LIB([lapack], [sgbbrd_], [true])
AC_CHECK_LIB([blas], [isamax_], [true])
AC_CHECK_LIB([g2c], [derf_], [true])
AC_CHECK_LIB(m, sqrt, [true])
AC_CHECK_LIB([pthread], [pthread_attr_destroy], [true])

# Checks for header files.
AC_HEADER_STDC
if test x"$ac_cv_header_stdc" = x'no'; then
   AC_MSG_WARN([[this package uses ANSI-C headers but none were found.]])
fi

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_MODE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STAT
AC_CHECK_FUNCS([floor pow sqrt strchr strstr])

# FFTW libraries
FFTW_LIBS="-lfftw3 -lfftw3f"
AC_SUBST([FFTW_LIBS])

# Cfitsio library
CFITSIO_LIBS="-lcfitsio"
AC_SUBST([CFITSIO_LIBS])

# lapack libraries
case ${MACHTYPE} in
#  *x86_64*)
#    LAPACK_LIBS="-lacml -lg2c";;
  *powerpc*)
     LAPACK_LIBS="-framework vecLib -dynamic";;
  *)
    LAPACK_LIBS="-llapack -lblas";;
# Uncomment this to use Goto's BLAS:
# http://www.cs.utexas.edu/users/kgoto/signup_first.html
# LAPACK_LIBS="-llapack -lgoto_p4_512-r0.94 -lg2c";;
esac
AC_SUBST([LAPACK_LIBS])

case ${OSTYPE} in
  *cygwin*)
    LIBPATH="-L/usr/local/lib";
    LIBOBJECTS=`(cd config; ./extract_libobjnames.sh ../lib/Makefile.am)`;
    static_only=1;;
  *)
    LIBPATH="";
    LIBOBJECTS="";
    static_only=0;;
esac
AC_SUBST([LIBPATH])
AC_SUBST([LIBOBJECTS])
AM_CONDITIONAL([STATIC_ONLY], [test $static_only -ne 0])

LDFLAGS="$LDFLAGS -L${SCIDIR}/bin"
LIBS="$LIBS $CFITSIO_LIBS $FFTW_LIBS $LAPACK_LIBS -lstdc++"

AC_CONFIG_FILES([
     Makefile 
     src/Makefile 
     src/arroyo/Makefile 
     src/arroyo/win-util/Makefile 
     src/lightPipes/Makefile 
     src/wrapper/Makefile 
     src/interf/Makefile 
     src/scicos/Makefile 
     man/eng/optics/Makefile 
     man/eng/optics_scicos/Makefile 
     man/Makefile 
     man/eng/Makefile 
     macros/Makefile 
     demos/Makefile 
     examples/Makefile])
	 
AC_OUTPUT

echo 
echo "Final Compilation Flags:"
echo "  CC = $CC"
echo "  CFLAGS = $CFLAGS"
echo "  CPPFLAGS = $CPPFLAGS"
echo "  LDFLAGS = $LDFLAGS"
echo "  LIBS = $LIBS"
echo 
echo ===========================================================================
echo "This toolbox will be installed in $TOOLBOXDIR."
echo 'Now type "make" then "make install".'
echo "After that, if you want the toolbox to be automatically loaded, type"
echo '"make autoload" (as normal user). This will append the proper command'
echo "to your .scilab file."
echo ===========================================================================
echo 
