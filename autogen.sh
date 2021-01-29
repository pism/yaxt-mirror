#!/bin/sh

autoreconf -fvi || exit $?

for patch_file in contrib/autogen_patches/00libtool.m4.arg_spaces.patch \
                  contrib/autogen_patches/01libtool.m4.mpi_wrappers.patch \
                  contrib/autogen_patches/02libtool.m4.nag_convenience.patch \
                  contrib/autogen_patches/03ltmain.sh.nag_pthread.patch; do
  patch -p 1 --no-backup-if-mismatch --forward -r - -i $patch_file
  # The program 'patch' exits with exitcode=1 if the patch has already been applied.
  # Consider this as a normal scenario:
  exitcode=$?; test $exitcode -ne 0 && test $exitcode -ne 1 && exit $exitcode
done

autoconf -f || exit $?

# Reset libtool.m4 timestamps to avoid confusing make
touch -r m4/ltversion.m4 m4/libtool.m4

