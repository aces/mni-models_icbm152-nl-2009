#
#   Copyright (C) 2004  Andrew Janke
#
#   This macro is in the public domain -- do as you wish with it
#
# AC_PROG_MAKE_MODEL
#
# Test for make_model
# and set $make_model to the correct value.
#
#
dnl @synopsis AC_PROG_MAKE_MODEL
dnl
dnl This macro tests if make_model is installed. If it 
dnl is installed, it set $make_model to the right value
dnl
dnl @version 1.0
dnl @author Andrew Janke <rotor@bic.mni.mcgill.ca>
dnl
AC_DEFUN([AC_PROG_MAKE_MODEL],[
AC_CHECK_PROGS(make_model,[make_model],no)
export make_model;
if test $make_model = "no" ;
then
	AC_MSG_ERROR([Unable to find make_model, do you have mni_autoreg installed?]);
fi
AC_SUBST(make_model)
])
