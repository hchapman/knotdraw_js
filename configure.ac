#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.



AC_PREREQ(2.57)

# Don't forget to change the LIB_VERSION (and read the note if you are
# wondering what the numbers mean.
AC_INIT(planarmap, 1.3.0, hchaps@gmail.com)
#AC_CONFIG_AUX_DIR([.])
AM_INIT_AUTOMAKE()
LT_INIT

# Output
AC_CONFIG_FILES([Makefile
                planarmap/Makefile])

AC_OUTPUT

echo \
" ----------------------------------------------------------------------

  ${PACKAGE_NAME} Version ${PACKAGE_VERSION}

  Installing to prefix: ${prefix}

  Compiler: '${CC} ${CFLAGS} ${CPPFLAGS}'

  Now type 'make @<:@<target>@:>@'

    where the optional <target> is:

    all					- build all binaries
    install				- install everything
    check 				- run self-tests


  -----------------------------------------------------------------------"
