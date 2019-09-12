#!/bin/sh

autoreconf -fvi || exit $?

patch -p 1 --no-backup-if-mismatch --forward -r - -i contrib/autogen_patches/libtool.m4.arg_spaces.patch -d .
# The program 'patch' exits with exitcode=1 if the patch has already been applied.
# Consider this as a normal scenario:
exitcode=$?; test $exitcode -ne 0 && test $exitcode -ne 1 && exit $exitcode

patch -p 1 --no-backup-if-mismatch --forward -r - -i contrib/autogen_patches/libtool.m4.nag_wrapper.patch -d .
exitcode=$?; test $exitcode -ne 0 && test $exitcode -ne 1 && exit $exitcode

autoconf -f || exit $?

# Reset libtool.m4 timestamps to avoid confusing make
touch -r m4/ltversion.m4 m4/libtool.m4

patch -p 1 --no-backup-if-mismatch --forward -r - -i contrib/autogen_patches/ltmain_nag_pthread.patch -d .
# The program 'patch' exits with exitcode=1 if the patch has already been applied.
# Consider this as a normal scenario:
exitcode=$?; if test $exitcode -ne 0 && test $exitcode -ne 1; then exit $exitcode; fi
