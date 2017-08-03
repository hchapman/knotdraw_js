#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# Autoconf 2.50 wants to disallow AM_ names.  We explicitly allow
# the ones we care about.
ifdef([m4_pattern_allow],
      [m4_pattern_allow([^AM_(CONFIG_HEADER|PROG_LEX)$])])dnl


AC_PREREQ(2.57)
AC_COPYRIGHT(
------------------------------- planarmap ---------------------------------
)

# Don't forget to change the LIB_VERSION (and read the note if you are
# wondering what the numbers mean.
AC_INIT(planarmap, 1.3.0, hchaps@gmail.com)
#AC_CONFIG_AUX_DIR([.])
AM_INIT_AUTOMAKE()

#
# The version number below has a very specific meaning.  It is in the form
# a:b:c.  The first number, a, is the interface version number.  Any time we
# add or remove a function, or change the order of parameters or change the way
# that the library reads or writes files, this number needs to be incremented.
# The second number, b, is the implementation number *of this interface*.  When
# a changes, b goes back to zero.  This is the number that needs to change
# every time we put out a new library distribution.  The third number, c, is
# the interface backward-compatibility count.  Suppose that we simply add a
# needed function to the library.  The new interface still provides everything
# that the old version provided and in the same way.  So the
# backward-compatibility count gets incremented.  Of course it gets reset to 0
# when a new interface is not binary-compatible with the old one.
#
# Due to various broken ideas, we have to have two versions of this number, one
# with colons below and one with periods above.  Please keep them in sync!
#
AC_SUBST(LIB_VERSION, 1:3:0)
# AC_CONFIG_SRCDIR([planarmap])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

#AC_GNU_SOURCE
LT_INIT

# Checks for programs.
#AC_PROG_AWK
#AC_PROG_YACC
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_CPP

#AM_PROG_LEX
AC_PROG_LN_S

AC_PROG_MAKE_SET
AC_PROG_LIBTOOL

# Checks for libraries. The syntax here is
#
# AC_SEARCH_LIBS([a function in the library],[the name of the library])
#
# if you find it, it will be automatically added to paths and the link commands
# for the project.
AC_SEARCH_LIBS([sqrt], [m])

AC_HEADER_STDC
AC_HEADER_ASSERT
AC_CHECK_HEADERS([stdlib.h math.h stdarg.h ctype.h stdio.h assert.h \
                  string.h stdbool.h float.h wait.h sys/wait.h unistd.h \
                  time.h sys/times.h sys/types.h fcntl.h signal.h complex.h \
		              limits.h strings.h stdint.h strings.h ])


# Signal to AUTOMAKE if we have argtable2

# AM_CONDITIONAL([HAVE_ARGTABLE2],[test "x${have_argtable2}" = xtrue])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_MODE_T
AC_TYPE_SIZE_T
AC_STRUCT_TM

# Checks for library functions.
#AC_FUNC_MALLOC
#AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_FUNC_STRFTIME

AC_CHECK_FUNCS([strcasestr mkstemp strcasecmp time floor sqrt ceil srandomdev sranddev random srandom])
AC_CHECK_FUNCS([bzero memset mkdir regcomp strchr strcspn strerror strrchr strstr strtol getline fstat64 memcpy])

AC_CHECK_FUNCS([bzero strchr strcspn strrchr strtol])
AC_CHECK_FUNC(getopt_long, SYS_GETOPTLONG=1, SYS_GETOPTLONG=0)
AC_CHECK_FUNC(regcomp,     SYS_REGEX=1,      SYS_REGEX=0)
AC_CHECK_FUNC(strptime,    SYS_STRPTIME=1,   SYS_STRPTIME=0)

# Define automake conditionals
AM_CONDITIONAL(USE_SYS_GETOPTLONG, test "$SYS_GETOPTLONG" = "1")
AM_CONDITIONAL(USE_ARGREX,         test "$SYS_REGEX" = "1")
AM_CONDITIONAL(USE_ARGDATE,        test "$SYS_STRPTIME" = "1")

# Can we inline and/or const?
AC_C_INLINE
AC_C_CONST

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
