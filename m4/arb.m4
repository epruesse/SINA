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
    success="no"
    if test x"$WANT_ARB" = x"yes"; then
        if test x"$GLIB_LIBS" = x""; then
            PKG_CHECK_MODULES([GLIB], [glib-2.0 > 2.2])
        fi
	if test x"$ax_arb_path" == x""; then
            AC_MSG_CHECKING([for location of ARB headers])
            for ax_arb_path_tmp in $ax_arb_path $ARBHOME /usr/lib/arb; do
                if test -f "$ax_arb_path_tmp/INCLUDE/arbdb.h" \
                && test -r "$ax_arb_path_tmp/INCLUDE/arbdb.h"; then
                    ax_arb_path="$ax_arb_path_tmp"
                    break
                fi
            done
            if test -n "$ax_arb_path"; then
                AC_MSG_RESULT([$ax_arb_path])
                success="yes"
            else
                AC_MSG_RESULT([not found])
                success="no"
            fi
	fi
	if test x"ax_arb_path" != x""; then
            saved_CPPFLAGS="$CPPFLAGS"
            CPPFLAGS="$CPPFLAGS -I$ax_arb_path/INCLUDE $GLIB_CFLAGS"
            AC_CHECK_HEADER([arbdb.h], [
                ax_arb_cppflags="-I$ax_arb_path/INCLUDE $GLIB_CFLAGS"
		success="yes"
            ],[
	        success="no"
	    ])
            CPPFLAGS="$saved_CPPFLAGS"
	fi
    fi

    if test x"$success" = x"yes"; then
        AC_MSG_CHECKING([for libARBDB location])

	ax_arb_lib_path=
	for libext in so dylib; do
            for ax_arb_lib_path_tmp in $ax_arb_path/lib /usr/lib/arb/lib; do
	        if test -f "$ax_arb_lib_path_tmp/libARBDB.$libext" \
                && test -r "$ax_arb_lib_path_tmp/libARBDB.$libext"; then
	            ax_arb_lib_path="$ax_arb_lib_path_tmp"
		    break
	        fi
                if test -f "$ax_arb_lib_path_tmp/libARBDB.$libext" \
                && test -r "$ax_arb_lib_path_tmp/libARBDB.$libext"; then
	            ax_arb_lib_path="$ax_arb_lib_path_tmp"
                    break
                fi
	    done
	done
	if test -n "$ax_arb_lib_path"; then
            AC_MSG_RESULT([$ax_arb_lib_path])
            success="yes"
        else
            AC_MSG_RESULT([not found])
            success="no"
        fi
    fi

    if test x"$success" = x"yes"; then
        AC_MSG_CHECKING([for GB_open in -lARBDB])
        saved_CPPFLAGS="$CPPFLAGS"
        saved_LIBS="$LIBS"
        saved_LDFLAGS="$LDFLAGS"

	ax_arb_ldflags="-L$ax_arb_lib_path -Wl,-rpath -Wl,$ax_arb_lib_path"
        ax_arb_libs="-lARBDB -lCORE $GLIB_LIBS"

        CPPFLAGS="$CPPFLAGS $ax_arb_cppflags"
        LIBS="$LIBS $ax_arb_libs"
        LDFLAGS="$LDFLAGS $ax_arb_ldflags"

        AC_LANG_PUSH(C++)
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
        AC_LANG_POP(C++)

        LIBS="$saved_LIBS"
        CPPFLAGS="$saved_CPPFLAGS"
        LDFLAGS="$saved_LDFLAGS"
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


### Check for a function in ARB
# Since ARB exports as C++, the full synopsis is needed.
# => Provide HAVE_xxx name in first argument and function in second argument
AC_DEFUN([AX_ARB_CHECK_FUNC],
[
    AC_REQUIRE([AX_LIB_ARBDB])
    AH_TEMPLATE([HAVE_$1])

    AC_MSG_CHECKING([for $1 in libARBDB])

    saved_CPPFLAGS="$CPPFLAGS"
    saved_LIBS="$LIBS"
    saved_LDFLAGS="$LDFLAGS"

    CPPFLAGS="$CPPFLAGS $ARB_CPPFLAGS"
    LIBS="$LIBS $ARB_LIBS"
    LDFLAGS="$LDFLAGS $ARB_LDFLAGS"

    AC_LANG_PUSH(C++)
    AC_LINK_IFELSE([
        AC_LANG_PROGRAM([[
            #include <arbdbt.h>
        ]], [[
            $2;
        ]])
    ],[
        AC_MSG_RESULT([yes])
        AC_DEFINE(HAVE_$1)
    ],[
        AC_MSG_RESULT([no])
    ])
    AC_LANG_POP(C++)

    LIBS="$saved_LIBS"
    CPPFLAGS="$saved_CPPFLAGS"
    LDFLAGS="$saved_LDFLAGS"
])

### Check for static PT server interface lib
AC_DEFUN([AX_LIB_ARB_PROBE],
[
    AC_REQUIRE([AX_LIB_ARBDB])
    AH_TEMPLATE([HAVE_ARB_PROBE], [], [Defined to 1 if ARB PROBE_COM is present])

    AC_MSG_CHECKING([for ARB PROBE_COM])
    if test -f "$ARBHOME/PROBE_COM/client.a" \
    && test -r "$ARBHOME"; then
        AC_MSG_RESULT([found])
        success="yes"
    else
        AC_MSG_RESULT([not found])
        success="no"
    fi

    if test x"$success" = x"yes"; then
        AC_MSG_CHECKING([for aisc_open in ARB PROBE_COM])
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
        saved_LDFLAGS="$LDFLAGS"
        saved_LIBS="$LIBS"

        CPPFLAGS="$CPPFLAGS $ARB_CPPFLAGS"
        LDFLAGS="$LDFLAGS $ARB_LDFLAGS"

        ax_arb_probe_libs=""
        for common in '' "$ARBHOME/PROBE_COM/common.a"; do
            ax_arb_probe_libs_tmp="$ARBHOME/PROBE_COM/client.a $common"
            LIBS="$saved_LIBS $ax_arb_probe_libs_tmp $ARB_LIBS"
            AC_LINK_IFELSE([], [
                ax_arb_probe_libs="$ax_arb_probe_libs_tmp"
                break
            ])
        done

        LIBS="$saved_LIBS"
        CPPFLAGS="$saved_CPPFLAGS"
        LDFLAGS="$saved_LDFLAGS"

        AC_LANG_POP(C++)

        success="yes"
        if test -z "$ax_arb_probe_libs"; then
            AC_MSG_RESULT([no])
            success="no"
        elif test -n "$common"; then
            AC_MSG_RESULT([yes])
        else
            AC_MSG_RESULT([yes (common.a not needed)])
        fi
    fi
    if test x"$success" = x"yes"; then
        ARB_PROBE_LIBS="$ax_arb_probe_libs"
        AC_SUBST(ARB_PROBE_LIBS)
        AC_DEFINE(HAVE_ARB_PROBE)
    fi
])

# Check for static HELIX interface lib
AC_DEFUN([AX_LIB_ARB_HELIX],
[
    AC_REQUIRE([AX_LIB_ARBDB])

    AH_TEMPLATE([HAVE_ARB_HELIX], [], [Defined to 1 if ARB SL HELIX is present])
    AC_MSG_CHECKING([for ARB HELIX static lib])

    saved_CPPFLAGS="$CPPFLAGS"
    saved_LDFLAGS="$LDFLAGS"
    saved_LIBS="$LIBS"

    CPPFLAGS="$CPPFLAGS $ARB_CPPFLAGS"
    LDFLAGS="$LDFLAGS $ARB_LDFLAGS"
    LIBS="$LIBS  $ARBHOME/SL/HELIX/HELIX.a $ARB_LIBS"

    AC_LANG_PUSH(C++)
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
    AC_LANG_POP(C++)

    LIBS="$saved_LIBS"
    CPPFLAGS="$saved_CPPFLAGS"
    LDFLAGS="$saved_LDFLAGS"
])

