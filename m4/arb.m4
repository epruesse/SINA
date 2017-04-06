AC_DEFUN([AX_LIB_ARBDB],
[
    AH_TEMPLATE([HAVE_ARB], [], [Defined to 1 if ARB libraries are present])
    # Check for dynamic libraries

    AC_ARG_WITH([arbhome],
        AC_HELP_STRING(
            [--with-arbhome=<PATH>],
            [point to ARBHOME]
        ),
        [
        if test x"$withval" = x"no"; then
            WANT_ARB="no"
        elif test x"$withval" = "yes"; then
            WANT_ARB="yes"
            ax_arb_path=""
        else
            WANT_ARB="yes"
            ax_arb_path="$withval"
        fi
        ],
        [WANT_ARB="yes"]
    )

    ARB_CPPFLAGS=""
    ARB_LDFLAGS=""
    ARB_LIBS=""

    if test x"$WANT_ARB" = x"yes"; then
        AC_MSG_CHECKING([for arbdb.h])
        for ax_arb_path_tmp in $ax_arb_path $ARBHOME /usr/lib/arb; do
            if test -f "$ax_arb_path_tmp/INCLUDE/arbdb.h" \
            && test -r "$ax_arb_path_tmp/INCLUDE/arbdb.h"; then
                ax_arb_path="$ax_arb_path_tmp"
                break
            fi
        done
        if test -n "$ax_arb_path"; then
            AC_MSG_RESULT([yes ($ax_arb_path/INCLUDE/arbdb.h)])
            success="yes"
        else
            AC_MSG_RESULT([no])
            success="no"
        fi
    fi
    
    if test x"$success" = x"yes"; then
        if test x"$GLIB_LIBS" = x""; then
            PKG_CHECK_MODULES([GLIB], [glib-2.0 > 2.2])
        fi

        AC_MSG_CHECKING([for libARBDB])

        ax_arb_ldflags="-L$ax_arb_path/lib"
        ax_arb_libs="-lARBDB -lCORE $GLIB_LIBS"
        ax_arb_cppflags="-I$ax_arb_path/INCLUDE $GLIB_CFLAGS"

        AC_LANG_PUSH(C++)
	
        saved_CPPFLAGS="$CPPFLAGS"
        CPPFLAGS="$CPPFLAGS $ax_arb_cppflags"
	saved_LIBS="$LIBS"
	LIBS="$LIBS $ax_arb_libs"
	saved_LDFLAGS="$LDFLAGS"
	LDFLAGS="$LDFLAGS $ax_arb_ldflags"

        AC_LINK_IFELSE([
            AC_LANG_PROGRAM([[
                #include <arbdb.h>
            ]], [[
                GB_open("","");
            ]])
        ],[
            AC_MSG_RESULT([yes])
            success="yes"
        ],[
            AC_MSG_RESULT([not found])
            success="no"
        ])
        CPPFLAGS="$saved_CPPFLAGS"
        AC_LANG_POP(C++)
    fi

    if test x"$success" == x"yes"; then
        ARB_CPPFLAGS="$ax_arb_cppflags"
        ARB_LDFLAGS="$ax_arb_ldflags"
	ARB_LIBS="$ax_arb_libs"
        ARBHOME="$ax_arb_path"
        AC_SUBST(ARB_CPPFLAGS)
        AC_SUBST(ARB_LDFLAGS)
	AC_SUBST(ARB_LIBS)
        AC_SUBST(ARBHOME)
        AC_DEFINE(HAVE_ARB)
    fi
])

### Check for static PT server interface lib
AC_DEFUN([AX_LIB_ARB_PROBE],
[
    AC_REQUIRE([AX_LIB_ARBDB])
    AH_TEMPLATE([HAVE_ARB_PROBE], [], [Defined to 1 if ARB PROBE_COM is present])
    
    AC_MSG_CHECKING([for ARB PROBE_COM static lib])
    if test -f "$ARBHOME/PROBE_COM/client.a" \
    && test -r "$ARBHOME"; then
        AC_MSG_RESULT([yes])
        success="yes"
    else
        AC_MSG_RESULT([no])
        success="no"
    fi

    if test x"$success" = x"yes"; then
        AC_MSG_CHECKING([whether ARB PROBE_COM needs common.a])
        AC_LANG_PUSH(C++)

        AC_LANG_CONFTEST([AC_LANG_PROGRAM([[
            #include "PT_com.h"
            #include "client.h"
        ]],[[
            T_PT_MAIN main;
            GB_ERROR err;
            aisc_open("", main, AISC_MAGIC_NUMBER, &err);
        ]])])

	saved_CPPFLAGS="$CPPFLAGS"
        CPPFLAGS="$CPPFLAGS $ARB_CPPFLAGS"

        saved_LIBS="$LIBS"
        ax_arb_probe_libs=""
        for common in '' "$ARBHOME/PROBE_COM/common.a"; do
            ax_arb_probe_libs_tmp="$ARBHOME/PROBE_COM/client.a $common"
            LIBS="$saved_LIBS $ax_arb_probe_libs_tmp"
            AC_LINK_IFELSE([], [
                ax_arb_probe_libs="$ax_arb_probe_libs_tmp"
                break
            ])
        done
        LIBS="$saved_LIBS"
	CPPFLAGS="$saved_CPPFLAGS"

        AC_LANG_POP(C++)

        success="yes"
        if test -z "$ax_arb_probe_libs"; then
            AC_MSG_RESULT([failed to link!])
            success="no"
        elif test -n "$common"; then
            AC_MSG_RESULT([yes])
        else
            AC_MSG_RESULT([no])
        fi
    fi
    if test x"$success" = x"yes"; then
        ARB_PROBE_LIBS="$ax_arb_probe_libs"
        AC_SUBST(ARB_PROBE_LIBS)
        AC_DEFINE(HAVE_ARB_PROBE)
	AC_REPLACE_FUNCS([GBT_find_sequence])
    fi
])

# Check for static HELIX interface lib
AC_DEFUN([AX_LIB_ARB_HELIX],
[
    AC_REQUIRE([AX_LIB_ARBDB])

    AH_TEMPLATE([HAVE_ARB_HELIX], [], [Defined to 1 if ARB SL HELIX is present])
    AC_MSG_CHECKING([for ARB HELIX static lib])

    saved_CPPFLAGS="$CPPFLAGS"
    CPPFLAGS="$CPPFLAGS $ARB_CPPFLAGS"
    saved_LIBS="$LIBS"
    LIBS="$LIBS $ARBHOME/SL/HELIX/HELIX.a"

    AC_LINK_IFELSE([
        AC_LANG_PROGRAM([[
            #include "BI_helix.hxx"
        ]],[[
            BI_helix h;
            h.init((GBDATA*)0,(const char*)0);
        ]])
    ],[
        AC_MSG_RESULT([yes])
        ARB_HELIX_LIBS="$ARBHOME/SL/HELIX/HELIX.a"
        AC_SUBST(ARB_HELIX_LIBS)
        AC_DEFINE(HAVE_ARB_HELIX)
    ],[
        AC_MSG_RESULT([no])
    ])
    
    LIBS="$saved_LIBS"
    CPPFLAGS="$saved_CPPFLAGS"
])

