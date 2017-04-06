# Init

AC_INIT([SINA],[1.3.1-a1],[elmar@pruesse.net])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src/sina.cpp])
AM_INIT_AUTOMAKE([parallel-tests foreign subdir-objects])
AC_PREREQ([2.50])

# build rules for silent compile
AM_SILENT_RULES([yes]) 

# Check for programs

# only run AR check if AM_PROG_AR defined 
# (automake > 1.11 needs it somehow, older ones don't have it)
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
AC_PROG_CXX
AC_PROG_LIBTOOL
LT_INIT

# Check for C++ and required features

AC_SUBST(CPPFLAGS)
AC_SUBST(CXXFLAGS)
AC_SUBST(LDFLAGS)
AC_LANG(C++)
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_HEADER_STDBOOL
AX_CXX_COMPILE_STDCXX([11],[noext],[mandatory])

# Check general functions

AC_CHECK_FUNCS([strdup])
AC_CHECK_FUNCS([gettimeofday])
AC_CHECK_FUNCS([getpid]) 

# Check headers 

AC_CHECK_HEADERS([sys/time.h])
AC_CHECK_HEADER(cstdio)


# Check for libraries

AC_CHECK_LIB(z, compress, [Z_LIB=-lz], [AC_MSG_ERROR([libz.so not found!])])
AC_SUBST(Z_LIB)
AC_CHECK_LIB(rt, clock_gettime)

# Define something for PSTREAMS (to get rid of warnings)
AH_TEMPLATE(REDI_EVISCERATE_PSTREAMS, 
            [Define interpreted by PSTREAMS header lib, only
             set to supress Wundef warnings])
AC_DEFINE(REDI_EVISCERATE_PSTREAMS, [0])

# Allways show warnings
CXXFLAGS="$CXXFLAGS -W -Wall -Wundef -Wsign-promo -Wfloat-equal -fno-common"

# Export all symbols to get helpful backtraces
# (don't on Darwin, ld doesn't support these options there)
if test `uname` = Linux; then 
  CXXFLAGS="$CXXFLAGS -rdynamic"
  LDFLAGS="$LDFLAGS -export-dynamic"
fi

# Add debug option
# automatically sets enable_debug=yes/no
AC_ARG_ENABLE([debug],
  AS_HELP_STRING([--enable-debug]),
  CPPFLAGS="$CPPFLAGS -DDEBUG"
  CXXFLAGS="$CXXFLAGS -O0 -ggdb3",
  #else
  CPPFLAGS="$CPPFLAGS -O3 -DNDEBUG")

# Compile with address sanitizer if supported and --enable-debug
AS_IF([test "x$enable_debug" = "xyes"],[
  AX_CHECK_COMPILE_FLAG([-fsanitize=address],[ 
    CXXFLAGS="$CXXFLAGS -fsanitize=address -fno-omit-frame-pointer"
  ])
  AX_CHECK_LINK_FLAG([-fsanitize=address], [
    LDFLAGS="$LDFLAGS -fsanitize=address"
  ])
])

# Add profiling option
AC_ARG_ENABLE([profiling],
  AS_HELP_STRING([--enable-profiling]),
  CXXFLAGS="$CXXFLAGS -pg"
  LDFLAGS="$LDFLAGS -pg",)

# Add buildno option
AC_ARG_WITH([buildinfo],
  AS_HELP_STRING([--with-buildinfo=<NO>]),
  AC_DEFINE_UNQUOTED([PACKAGE_BUILDINFO], ["${withval}"], [Add extra buildinfo]))

####  Check for ARB and dependencies

# Get ARBHOME
AC_ARG_WITH([arbhome], AS_HELP_STRING([--with-arbhome=<PATH>]), 
  ARBHOME=${withval}
)
AC_SUBST(ARBHOME)

# Add ARB include and lib dirs to search paths
CPPFLAGS="$CPPFLAGS -I$ARBHOME/INCLUDE"
LDFLAGS="$LDFLAGS -L$ARBHOME/lib"

# Require GLib-2.0
PKG_CHECK_MODULES([GLIB], [glib-2.0 > 2.2])
LIBS="${GLIB_LIBS} $LIBS"
CPPFLAGS="$CPPFLAGS $GLIB_CFLAGS"

# Check for dynamic libraries
AC_SEARCH_LIBS([_Z14GB_print_errorv], [CORE])
AC_SEARCH_LIBS([_Z9GB_getcwdv], [ARBDB])
AC_REPLACE_FUNCS([GBT_find_sequence])

# Dynamic libs that need to be copied
if test `uname` = Darwin; then
  # need to copy glib for darwin (only present if provided by macports)
  GLIB_LIBPATH=`$PKG_CONFIG --variable=libdir glib-2.0`
  ARB_FILES="$ARB_FILES $GLIB_LIBPATH/libglib-2.0.0.dylib"
else
  ARB_FILES="$ARB_FILES $ARBHOME/lib/libCORE.so"
  ARB_FILES="$ARB_FILES $ARBHOME/lib/libARBDB.so"
fi


### Check for static PT server interface lib
AC_MSG_CHECKING([for ARB PROBE_COM static lib])
SAVED_LIBS="$LIBS"

# create test program
AC_LANG_CONFTEST([AC_LANG_PROGRAM(
        [[
#include "PT_com.h"
#include "client.h"
        ]],[[
T_PT_MAIN main;
GB_ERROR err;
aisc_open("", main, AISC_MAGIC_NUMBER, &err);
        ]]
)])

# test if common.a is needed
ARB_SLIBS=
for common in '' "$ARBHOME/PROBE_COM/common.a"; do
  LIBS="$ARBHOME/PROBE_COM/client.a $common $SAVED_LIBS"
  AC_LINK_IFELSE([],
    [ARB_SLIBS="$ARBHOME/PROBE_COM/client.a $common $ARB_SLIBS"
     break])
done
AS_IF([ test -n "$ARB_SLIBS"],
  AS_IF([test -z "$common"],
    AC_MSG_RESULT([yes]),
    AC_MSG_RESULT([yes (needs common.a)])
  ),
  AC_MSG_FAILURE([not found])
)

LIBS="$SAVED_LIBS"

# Check for static HELIX interface lib
AC_MSG_CHECKING([for ARB HELIX static lib])
SAVED_LIBS="$LIBS"
LIBS="$ARBHOME/SL/HELIX/HELIX.a $LIBS"
AC_LINK_IFELSE(
    [AC_LANG_PROGRAM(
        [#include "BI_helix.hxx"],
        [BI_helix h; h.init((GBDATA*)0,(const char*)0)];
    )],
    AC_MSG_RESULT(yes)
    ARB_SLIBS="$ARBHOME/SL/HELIX/HELIX.a $ARB_SLIBS"        
    ,
    AC_MSG_FAILURE([not found])
    )
LIBS="$SAVED_LIBS"

AC_SUBST(ARB_FILES)
AC_SUBST(ARB_LDFLAGS)
AC_SUBST(ARB_CPPFLAGS)
AC_SUBST(ARB_LIBS)
AC_SUBST(ARB_SLIBS)

# check boost libs (slow...)
BOOST_REQUIRE([1.34.1])
BOOST_THREADS([mt-s])
BOOST_PROGRAM_OPTIONS([mt-s])
BOOST_SERIALIZATION([mt-s])
BOOST_FOREACH([mt-s])
BOOST_SMART_PTR([mt-s])
BOOST_STATICASSERT([mt-s])
BOOST_STRING_ALGO([mt-s])
BOOST_SYSTEM([mt-s])
BOOST_TEST([mt])
BOOST_TUPLE([mt-s])
BOOST_VARIANT([mt-s])
BOOST_STATIC()


# Generate output

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
  Makefile
  src/sina
], [chmod +x src/sina])
AC_OUTPUT
