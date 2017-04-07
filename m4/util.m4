# version of AX_ARG_ENABLE with less boiler plate

AC_DEFUN([AX_ARG_ENABLE],
[
  _AX_ARG_ENABLE([$1],[m4_translit([[$1]], [_], [-])],[$2],[$3],[$4],m4_default([$5], [no]))
])

AC_DEFUN([_AX_ARG_ENABLE],
[
      AC_ARG_ENABLE([$1], AS_HELP_STRING([m4_case([$6],
                                                  [no],[--enable-$2],
					          [yes],[--disable-$2])],
					 [$3]),
                    AS_CASE(["$enableval"],
		            [yes|no], [enable_$1="$enableval"],
			    AC_MSG_ERROR([bad value ${enableval} for $2])
			    ),
	            [enable_$1="$6"]
      )
      AS_CASE([$enable_$1], [yes], [$4], [no], [$5])
      AM_CONDITIONAL(m4_translit([[enable_$1]],[a-z],[A-Z]), [test x"$enable_$1" = "xyes"])
])

AC_DEFUN([AX_ARG_DISABLE],[AX_ARG_ENABLE([$1],[$2],[$3],[$4],[yes])])
