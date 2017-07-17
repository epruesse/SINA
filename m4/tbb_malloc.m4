#
# SYNOPSIS
#
#   AX_LIB_TBB_MALLOC
#
# DESCRIPTION
#
#   Test for Intel Threaded Building Blocks libraries.
#   Tries to use _debug libs if $enable_debug=yes
#
#   This macro calls:
#
#     AC_SUBST(TBB_MALLOC_LIB)
#     AC_SUBST(TBB_MALLOC_LDFLAGS)
#     AC_SUBST(TBB_MALLOC_CPPFLAGS)
#
#   And sets:
#
#     HAVE_TBB_MALLOC
#     TBB_USE_DEBUG
#
# LICENSE
#
#   Copyright (c) 2017 Elmar Pruesse <elmar@pruesse.net>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

AC_DEFUN([AX_LIB_TBB_MALLOC],
[
    AH_TEMPLATE([HAVE_TBB_MALLOC], [Define to 1 if you have (and want to use) Intel TBB])
    AH_TEMPLATE([TBB_USE_DEBUG], [Define to 1 to use TBB debugging])

    AC_ARG_WITH([tbb-malloc],
        [AC_HELP_STRING([--with-tbb-malloc@<:@=ARG@:>@],
	    [use TBB Malloc library from standard location (ARG=yes),
	     from the specified location (ARG=<path>),
	     or disable it (ARG=no)
	     @<:@ARG=yes@: >@ ])],
        [
	if test "$withval" = "no"; then
            want_tbb_malloc="no"
	elif test "$withval" = "yes"; then
	    want_tbb_malloc="yes"
        else
            want_tbb_malloc="yes"
	    ax_tbb_malloc_path="$withval"
	fi
	],[want_tbb_malloc="yes"])

    TBB_MALLOC_LIB=
    TBB_MALLOC_CPPFLAGS=
    TBB_MALLOC_LDFLAGS=

    if test "$ax_tbb_malloc_path" != ""; then
        TBB_MALLOC_CPPFLAGS="-I$ax_tbb_malloc_path/include"
	TBB_MALLOC_LDFLAGS="-L$ax_tbb_malloc_path/lib"
    fi

    if test "$want_tbb_malloc" = "yes"; then
        saved_CPPFLAGS="$CPPFLAGS"
        saved_LDFLAGS="$LDFLAGS"
	CPPFLAGS="$CPPFLAGS $TBB_MALLOC_CPPFLAGS"
	LDFLAGS="$LDFLAGS $TBB_MALLOC_LDFLAGS"

        AC_CHECK_HEADER([tbb/tbb.h],[
            if test x"$enable_debug" = x"yes"; then
                AC_CHECK_LIB([tbbmalloc_debug], [main], [
                    TBB_MALLOC_LIB=-ltbbmalloc_debug
		    AC_DEFINE([TBB_USE_DEBUG])
                ])
            fi
            if test -z "$TBB_MALLOC_LIB"; then
                AC_CHECK_LIB([tbbmalloc], [main], [
                    TBB_MALLOC_LIB=-ltbbmalloc
		])
	    fi
        ])
    fi

    if test -n "$TBB_MALLOC_LIB"; then
        AC_DEFINE([HAVE_TBB_MALLOC])
    fi

    AC_SUBST([TBB_MALLOC_LIB])
    AC_SUBST([TBB_MALLOC_LDFLAGS])
    AC_SUBST([TBB_MALLOC_CFLAGS])
])
