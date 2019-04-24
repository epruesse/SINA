# Init

AC_INIT([SINA],[1.6.0-rc.1],[elmar@pruesse.net],,[http://github.com/epruesse/SINA])
AC_COPYRIGHT([Copyright (c) 2005-2018 Elmar Pruesse])
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

# AC_PROG_CXX sets -O2, which we don't want if debug enabled, so we store this
# as default flags for later decision
saved_CXXFLAGS="$CXXFLAGS"
AC_PROG_CXX
default_CXXFLAGS="$CXXFLAGS"
CXXFLAGS="$saved_CXXFLAGS"
AC_PROG_LIBTOOL
AC_PROG_SED
AC_PROG_LN_S
AC_PROG_INSTALL
AC_PROG_MKDIR_P
AC_CHECK_TOOL([BC], [bc], [:])
AC_CHECK_TOOL([TEE], [tee], [:])
AC_CHECK_TOOL([MKTEMP], [mktemp], [:])
AC_CHECK_TOOL([SPHINX_BUILD], [sphinx-build], [:])
AC_PROG_GREP

LT_INIT([shared static])

AC_REQUIRE_AUX_FILE([tap-driver.sh])

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

# Need dladdr to locate libARBDB and thereby ARBHOME
AC_SEARCH_LIBS([dladdr], [dl dld], [], [
  AC_MSG_ERROR([unable to find the dladdr() function])
])

# Define something for PSTREAMS (to get rid of warnings)
AH_TEMPLATE(REDI_EVISCERATE_PSTREAMS, 
            [Define interpreted by PSTREAMS header lib, only
             set to supress Wundef warnings])
AC_DEFINE(REDI_EVISCERATE_PSTREAMS, [0])

# Add coverage option
AX_CODE_COVERAGE
AM_COND_IF([CODE_COVERAGE_ENABLED], [
  AC_CACHE_CHECK([for coverage libs], [_cv_coverage_libs], [
    _cv_coverage_libs=no
    save_LIBS=$LIBS
    res=
    AC_LANG_CONFTEST([AC_LANG_CALL([], [main])])
    for lib in '' "$CODE_COVERAGE_LIBS"; do
      if test -z "$lib"; then
        res="none required"
      else
        res=$CODE_COVERAGE_LIBS
        LIBS="-l$lib $save_LIBS"
      fi
      AC_LINK_IFELSE([], [_cv_coverage_libs=$res])
      AS_IF([test -n "$_cv_coverage_libs"], [break])
    done
    LIBS=$save_LIBS
  ])
  AC_MSG_NOTICE([here - $_cv_coverage_libs])
  AS_IF([test x"$_cv_coverage_libs" != x"no"], [
    AS_IF([test x"$_cv_coverage_libs" == x"none required"], [
      CODE_COVERAGE_LIBS=
    ],[
      CODE_COVERAGE_LIBS=$_cv_coverage_libs
    ])
  ], [
    missing_libs="$missing_libs
    coverage libs ($CODE_COVERAGE_LIBS)"
  ])
])

# Export all symbols to get helpful backtraces
# (don't on Darwin, ld doesn't support these options there)
#if test `uname` = Linux; then 
#  CXXFLAGS="$CXXFLAGS -rdynamic"
#  LDFLAGS="$LDFLAGS -export-dynamic"
#fi

# Add debug option
AH_TEMPLATE([DEBUG])
AH_TEMPLATE([NDEBUG])
AX_ARG_ENABLE([debug], [disable optimizations and add debug symbols],
[
  AM_COND_IF([CODE_COVERAGE_ENABLED], [
    AC_MSG_WARN([Code Coverage enabled: using -DNDEBUG instead of -DDEBUG])
    AC_DEFINE([NDEBUG])
  ],[
    AC_DEFINE([DEBUG])
  ])
  CXXFLAGS="$(echo $CXXFLAGS | sed 's/-O@<:@@<:@:alnum:@:>@@:>@*//g')"
  CXXFLAGS="$CXXFLAGS -Og -ggdb3"
  AX_CHECK_COMPILE_FLAG([-fno-omit-frame-pointer],[
    CXXFLAGS="$CXXFLAGS -fno-omit-frame-pointer"
  ])
],[
  CPPFLAGS="$CPPFLAGS -DNDEBUG"
  CXXFLAGS="$default_CXXFLAGS $CXXFLAGS"
])

# Add address sanitizer option
AX_ARG_ENABLE([asan], [enable AddressSanitizer], [
  AX_CHECK_COMPILE_FLAG([-fsanitize=address],[
    CXXFLAGS="$CXXFLAGS -fsanitize=address"
  ])
  AX_CHECK_LINK_FLAG([-fsanitize=address],[
    LDFLAGS="$LDFLAGS -fsanitize=address"
  ])
])


# Add option to disable building docs (needs latex...)
AX_ARG_DISABLE([docs], [don't build documentation])

AX_ARG_ENABLE([fat_tar], [try to build a fat tar containing deps])

# Add profiling option
AX_ARG_ENABLE([profiling], [add profiling information to build],
[
  CXXFLAGS="$CXXFLAGS -pg"
  LDFLAGS="$LDFLAGS -pg"
])


# Add buildno option
AC_ARG_WITH([buildinfo], AS_HELP_STRING([--with-buildinfo=<NO>]),
  AC_DEFINE_UNQUOTED([PACKAGE_BUILDINFO], ["${withval}"], [Add extra buildinfo])
  PACKAGE_BUILDINFO="$withval"
)


# ARB dependencies
PKG_CHECK_MODULES([GLIB], [glib-2.0 > 2.2])
if test x"$GLIB_LIBS" = x""; then
  missing_libs="$missing_libs
  glib-2.0 (libglib)"
fi

AX_LIB_ARBDB
if test x"$ARB_LIBS" = x""; then
  missing_libs="$missing_libs
  ARB libraries (libARBDB)"
fi

AX_LIB_ARB_PROBE
if test x"$ARB_PROBE_LIBS" = x""; then
  missing_libs="$missing_libs
  ARB PROBE library (PROBE_COM/client.a)"
fi

AX_LIB_ARB_HELIX
if test x"$ARB_HELIX_LIBS" = x""; then
  missing_libs="$missing_libs
  ARB HELIX library (SL/HELIX/HELIX.a)"
fi

AX_ARB_CHECK_FUNC([GBT_FIND_SEQUENCE], [GBT_find_sequence(NULL, NULL)])
AX_ARB_STATUS_RETURN_TYPE

case "$host_os" in
  darwin*)
    LIBEXT=.dylib
    OSNAME=macos
    ;;
  *)
    LIBEXT=.so
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

AX_LIB_TBB
if test x"$TBB_LIB" = x""; then
  missing_libs="$missing_libs
  Intel Threading Building Blocks (libtbb)"
fi

AX_LIB_TBB_MALLOC

AX_BOOST_BASE([1.62])
AX_BOOST_THREAD
if test x"$BOOST_THREAD_LIB" = x""; then
  missing_libs="$missing_libs
  BOOST Thread (libboost_thread)"
fi
AX_BOOST_PROGRAM_OPTIONS
if test x"$BOOST_PROGRAM_OPTIONS_LIB" = x""; then
  missing_libs="$missing_libs
  BOOST Program Options (libboost_program_options)"
fi
AX_BOOST_IOSTREAMS
if test x"$BOOST_IOSTREAMS_LIB" = x""; then
  missing_libs="$missing_libs
  BOOST Program Options (libboost_iostreams)"
fi
AX_BOOST_FILESYSTEM
if test x"$BOOST_FILESYSTEM_LIB" = x""; then
  missing_libs="$missing_libs
  BOOST Program Options (libboost_filesystem)"
fi

AX_BOOST_SYSTEM
if test x"$BOOST_SYSTEM_LIB" = x""; then
  missing_libs="$missing_libs
  BOOST System (libboost_system)"
fi
AX_BOOST_UNIT_TEST_FRAMEWORK
#if test x"$BOOST_UNIT_TEST_FRAMEWORK_LIBS" = x""; then
#  missing_libs="$missing_libs\n BOOST Unit Test Framework (libboost_unit_test_framework)"
#fi

# Allways show warnings
CXXFLAGS="$CXXFLAGS -W"

# Generate output

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
  Makefile
  tests/test_helper.sh
])


###### print some summaries

AC_MSG_NOTICE([---------------------------------------])
AC_MSG_NOTICE([  AC_PACKAGE_STRING summary])
AC_MSG_NOTICE([---------------------------------------])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([  prefix: $prefix])
AC_MSG_NOTICE([  build:  $build])
AC_MSG_NOTICE([  host:   $host])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([  debug:     $enable_debug])
AC_MSG_NOTICE([  profiling: $enable_profiling])
AC_MSG_NOTICE([  coverage:  $CODE_COVERAGE_ENABLED])
AC_MSG_NOTICE([  buildinfo: $PACKAGE_BUILDINFO])
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
AM_COND_IF([CODE_COVERAGE_ENABLED], [
AC_MSG_NOTICE([  CODE COVERAGE:])
AC_MSG_NOTICE([    CPPFLAGS:        $CODE_COVERAGE_CPPFLAGS])
AC_MSG_NOTICE([    CFLAGS:          $CODE_COVERAGE_CFLAGS])
AC_MSG_NOTICE([    CXXFLAGS:        $CODE_COVERAGE_CXXFLAGS])
AC_MSG_NOTICE([    LIBS:            $CODE_COVERAGE_LIBS])
AC_MSG_NOTICE([])
])
AC_MSG_NOTICE([  ARB:])
AC_MSG_NOTICE([    ARBHOME:         ${ARBHOME}])
AC_MSG_NOTICE([    CPPFLAGS:        ${ARB_CPPFLAGS}])
AC_MSG_NOTICE([    LDFLAGS:         ${ARB_LDFLAGS}])
AC_MSG_NOTICE([    LIBS:            ${ARB_LIBS}])
AC_MSG_NOTICE([    PROBE LIBS:      ${ARB_PROBE_LIBS}])
AC_MSG_NOTICE([    HELIX LIBS:      ${ARB_HELIX_LIBS}])
AC_MSG_NOTICE([  BOOST:])
AC_MSG_NOTICE([    LDFLAGS:         ${BOOST_LDFLAGS}])
AC_MSG_NOTICE([    CPPFLAGS:        ${BOOST_CPPFLAGS}])
AC_MSG_NOTICE([    PROGRAM_OPTIONS: ${BOOST_PROGRAM_OPTIONS_LIB}])
AC_MSG_NOTICE([    THREAD:          ${BOOST_THREAD_LIB}])
AC_MSG_NOTICE([    IOSTREAMS:       ${BOOST_IOSTREAMS_LIB}])
AC_MSG_NOTICE([    SYSTEM:          ${BOOST_SYSTEM_LIB}])
AC_MSG_NOTICE([    FILESYSTEM:      ${BOOST_FILESYSTEM_LIB}])
AC_MSG_NOTICE([    UNIT_TEST:       ${BOOST_UNIT_TEST_FRAMEWORK_LIB}])
AC_MSG_NOTICE([  zlib:])
AC_MSG_NOTICE([    LIBS:            ${Z_LIB}])
AC_MSG_NOTICE([  tbb:])
AC_MSG_NOTICE([    CPPFLAGS:        ${TBB_CPPFLAGS}])
AC_MSG_NOTICE([    LDFLAGS:         ${TBB_LDFLAGS}])
AC_MSG_NOTICE([    LIBS:            ${TBB_LIB}])
AC_MSG_NOTICE([  tbb malloc:])
AC_MSG_NOTICE([    CPPFLAGS:        ${TBB_MALLOC_CPPFLAGS}])
AC_MSG_NOTICE([    LDFLAGS:         ${TBB_MALLOC_LDFLAGS}])
AC_MSG_NOTICE([    LIBS:            ${TBB_MALLOC_LIB}])


### check that all deps are satisfied

if test x"$missing_libs" != x""; then
  AC_MSG_FAILURE(Required libraries found missing:  $missing_libs )
fi

AC_OUTPUT
