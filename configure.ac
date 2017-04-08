# Init

AC_INIT([SINA],[1.3.1-a3],[elmar@pruesse.net],,[http://github.com/epruesse/SINA])
AC_COPYRIGHT([Copyright (c) 2005-2017 Elmar Pruesse])
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
AC_PROG_SED
AC_PROG_LN_S
AC_PROG_INSTALL
AC_PROG_MKDIR_P
LT_INIT

# Check for C++ and required features

AC_SUBST(CPPFLAGS)
AC_SUBST(CXXFLAGS)
AC_SUBST(LDFLAGS)
AC_LANG([C++])
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
CXXFLAGS="$CXXFLAGS -w"

# Export all symbols to get helpful backtraces
# (don't on Darwin, ld doesn't support these options there)
#if test `uname` = Linux; then 
#  CXXFLAGS="$CXXFLAGS -rdynamic"
#  LDFLAGS="$LDFLAGS -export-dynamic"
#fi

# Add debug option
AX_ARG_ENABLE([debug], [disable optimizations and add debug symbols],
[
  CPPFLAGS="$CPPFLAGS -DDEBUG"
  CXXFLAGS="$CXXFLAGS -O0 -ggdb3"
  AX_CHECK_COMPILE_FLAG([-fsanitize=address],[ 
    CXXFLAGS="$CXXFLAGS -fsanitize=address -fno-omit-frame-pointer"
  ])
  AX_CHECK_LINK_FLAG([-fsanitize=address], [
    LDFLAGS="$LDFLAGS -fsanitize=address"
  ])
],[
  CPPFLAGS="$CPPFLAGS -DNDEBUG"
])

# Add option to disable building docs (needs latex...)
AX_ARG_DISABLE([docs], [don't build documentation])

# Add profiling option
AX_ARG_ENABLE([profiling], [add profiling information to build],
[
  CXXFLAGS="$CXXFLAGS -pg"
  LDFLAGS="$LDFLAGS -pg"
])

# Add coverage option
AX_CODE_COVERAGE

# Add buildno option
AC_ARG_WITH([buildinfo], AS_HELP_STRING([--with-buildinfo=<NO>]),
  AC_DEFINE_UNQUOTED([PACKAGE_BUILDINFO], ["${withval}"], [Add extra buildinfo])
  PACKAGE_BUILDINFO="$withval"
)

# ARB dependencies
PKG_CHECK_MODULES([GLIB], [glib-2.0 > 2.2])
if test x"$GLIB_LIBS" = x""; then
  missing_libs="$missing_libs\n glib-2.0 (libglib)"
fi

AX_LIB_ARBDB
if test x"$ARB_LIBS" = x""; then
  missing_libs="$missing_libs\n ARB libraries (libARBDB)"
fi

AX_LIB_ARB_PROBE
if test x"$ARB_PROBE_LIBS" = x""; then
  missing_libs="$missing_libs\n ARB PROBE library (PROBE_COM/client.a)"
fi

AX_LIB_ARB_HELIX
if test x"$ARB_HELIX_LIBS" = x""; then
  missing_libs="$missing_libs\n ARB HELIX library (SL/HELIX/HELIX.a)"
fi


case "$host_os" in
  darwin*)
#    LIBEXT=.dylib
    OSNAME=macos
    ;;
  *)
#    LIBEXT=.so
    OSNAME=linux
    ;;
esac
AC_SUBST(LIBEXT)
AC_SUBST(OSNAME)


# Dynamic libs that need to be copied
if test `uname` = Darwin; then
  # need to copy glib for darwin (only present if provided by macports)
  GLIB_LIBPATH=`$PKG_CONFIG --variable=libdir glib-2.0`
  ARB_FILES="$ARB_FILES $GLIB_LIBPATH/libglib-2.0.0.dylib"
fi

for lib in CORE ARBDB; do 
  if test -e "$ARBHOME/lib/lib${lib}${LIBEXT}"; then
    ARB_FILES="$ARB_FILES $ARBHOME/lib/lib${lib}${LIBEXT}"
  fi
done
AC_SUBST(ARB_FILES)


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
], [
  chmod +x src/sina
])



###### print some summaries

AC_MSG_NOTICE([---------------------------------------])
AC_MSG_NOTICE([  SINA summary])
AC_MSG_NOTICE([---------------------------------------])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([  prefix: $prefix])
AC_MSG_NOTICE([  build:  $build])
AC_MSG_NOTICE([  host:   $host])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([  debug:     $enable_debug])
AC_MSG_NOTICE([  profiling: $enable_profiling])
AC_MSG_NOTICE([  coverage:  $enable_coverage])
AC_MSG_NOTICE([  buildinfo:  $PACKAGE_BUILDINFO])
AC_MSG_NOTICE([  docs:      $enable_docs])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([  CFLAGS:   $CFLAGS])
AC_MSG_NOTICE([  CPPFLAGS: $CPPFLAGS])
AC_MSG_NOTICE([  CXXFLAGS: $CXXFLAGS])
AC_MSG_NOTICE([  LIBS:     $LIBS])
AC_MSG_NOTICE([  LDFLAGS:  $LDFLAGS])
AC_MSG_NOTICE([  CC:       $CC])
AC_MSG_NOTICE([  CXX:      $CXX])
AC_MSG_NOTICE([  LD:       $LD])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([  ARB:            ${ARBHOME}])
AC_MSG_NOTICE([    CPPFLAGS:     ${ARB_CPPFLAGS}])
AC_MSG_NOTICE([    LDFLAGS:      ${ARB_LDFLAGS}])
AC_MSG_NOTICE([    LIBS:         ${ARB_LIBS}])
AC_MSG_NOTICE([    PROBE LIBS:   ${ARB_PROBE_LIBS}])
AC_MSG_NOTICE([    HELIX LIBS:   ${ARB_HELIX_LIBS}])
AC_MSG_NOTICE([  boost:])
AC_MSG_NOTICE([    ldpath:     ${BOOST_LDPATH}])
AC_MSG_NOTICE([    cppflags:   ${BOOST_CPPFLAGS}])
AC_MSG_NOTICE([  zlib:         ${Z_LIB}])

### check that all deps are satisfied

if test x"$missing_libs" != x""; then
  AC_MSG_FAILURE([ Required libraries found missing:
  $missing_libs
  ])
fi

AC_OUTPUT
