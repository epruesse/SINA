#
# SYNOPSIS
#
#   AX_LIB_TBB
#
# DESCRIPTION
#
#   Test for Intel Threaded Building Blocks libraries.
#   Tries to use _debug libs if $enable_debug=yes
#
#   This macro calls:
#
#     AC_SUBST(TBB_LIB)
#     AC_SUBST(TBB_LDFLAGS)
#     AC_SUBST(TBB_CPPFLAGS)
#
#   And sets:
#
#     HAVE_TBB
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

AC_DEFUN([AX_LIB_TBB],
[
    AH_TEMPLATE([HAVE_TBB], [Define to 1 if you have (and want to use) Intel TBB])
    AH_TEMPLATE([TBB_USE_DEBUG], [Define to 1 to use TBB debugging])

    AC_ARG_WITH([tbb],
        [AC_HELP_STRING([--with-tbb@<:@=ARG@:>@],
	    [use TBB library from standard location (ARG=yes),
	     from the specified location (ARG=<path>),
	     or disable it (ARG=no)
	     @<:@ARG=yes@: >@ ])],
        [
	if test "$withval" = "no"; then
            want_tbb="no"
	elif test "$withval" = "yes"; then
	    want_tbb="yes"
	else
	    want_tbb="yes"
	    ax_tbb_path="$withval"
	fi
	],[want_tbb="yes"])

    TBB_LIB=
    TBB_CPPFLAGS=
    TBB_LDFLAGS=

    if test "$ax_tbb_path" != ""; then
        TBB_CPPFLAGS="-I$ax_tbb_path/include"
	TBB_LDFLAGS="-L$ax_tbb_path/lib"
    fi

    if test "$want_tbb" = "yes"; then
        saved_CPPFLAGS="$CPPFLAGS"
        saved_LDFLAGS="$LDFLAGS"
	CPPFLAGS="$CPPFLAGS $TBB_CPPFLAGS"
	LDFLAGS="$LDFLAGS $TBB_LDFLAGS"

        AC_CHECK_HEADER([tbb/tbb.h],[
            if test x"$enable_debug" = x"yes"; then
                AC_CHECK_LIB([tbb_debug], [main], [
                    TBB_LIB=-ltbb_debug
		    AC_DEFINE([TBB_USE_DEBUG])
                ])
            fi
            if test -z "$TBB_LIB"; then
                AC_CHECK_LIB([tbb], [main], [
                    TBB_LIB=-ltbb
		])
	    fi
        ])
    fi

    if test -n "$TBB_LIB"; then
        AC_DEFINE([HAVE_TBB])
    fi

    AC_SUBST([TBB_LIB])
    AC_SUBST([TBB_LDFLAGS])
    AC_SUBST([TBB_CFLAGS])
])
