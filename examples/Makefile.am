noinst_PROGRAMS = 						\
	rrobin

rrobin_SOUES = = rrobin.c

AM_CFLAGS = -I$(top_srcdir)/src $(MPI_C_INCLUDE)
AM_FCFLAGS = -I$(top_srcdir)/src -I../src $(MPI_FC_INCLUDE) $(MPI_FC_MOD)
LDADD = ../src/libyaxt.a $(MPI_C_LIB)
