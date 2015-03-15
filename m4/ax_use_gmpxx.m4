
AC_DEFUN([AX_USE_GMPXX],
[
	AC_ARG_WITH([gmp-support],
        AS_HELP_STRING([--with-gmp-support], [add runtime support for the GMP big number library]),
        [
        if test "$withval" = "no"; then
            want_gmp="no"
        elif test "$withval" = "yes"; then
            want_gmp="yes"
            ax_gmp_lib=""
        else
            want_boost="gmp"
            ax_gmp_lib="$withval"
            fi
        ],
        [want_gmp="no"]
        )

LIBGMP=
          AS_IF([test "x$want_gmp" != xno],
            [AC_CHECK_LIB([gmpxx], [abs],
              [AC_SUBST([LIBGMP], ["-lgmp -lgmpxx"])
               AC_DEFINE([HAVE_GMP], [], [Enables GMP])
              ],
              [AC_MSG_FAILURE(
                 [--with-gmp-support was given, but test for gmp failed])],
              [])])


])


          

