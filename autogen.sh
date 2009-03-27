#! /bin/sh

set -e
 
aclocal -I local_m4
autoheader
automake --add-missing
autoreconf
