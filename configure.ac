AC_INIT(pyinit, 1.0.2)

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi

RBIN="${R_HOME}/bin/R"
CC=`"${RBIN}" CMD config CC`
CPP=`"${RBIN}" CMD config CPP`

AC_PROG_CC
AC_PROG_CPP

####################################################################################################
##
## Check for necessary support by C compiler
##
####################################################################################################
AC_LANG(C)
AC_CHECK_HEADERS(stdint.h inttypes.h limits.h float.h)

#check for -restrict CFLAGS support
oldCflags="$CFLAGS"
CFLAGS="$CFLAGS -restrict"
AC_MSG_CHECKING([whether C supports -restrict])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([])],
    [AC_MSG_RESULT([yes])]
    [RESTRICT_CFLAGS=-restrict],
    [AC_MSG_RESULT([no])]
)
CFLAGS="$oldCflags"
AC_SUBST([RESTRICT_CFLAGS])
LIBS="$oldLIBS"

AC_CONFIG_HEADERS([src/autoconfig.h])
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
