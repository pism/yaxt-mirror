noinst_PROGRAMS = 						\
	rrobin							\
	idxlist_examples

rrobin_SOURCES = rrobin.c
idxlist_examples_SOURCES = idxlist_examples.c

AM_CFLAGS = -I../inst_headers $(MPI_C_INCLUDE)
AM_FCFLAGS = -I$(top_srcdir)/src -I../src $(MPI_FC_INCLUDE) $(MPI_FC_MOD)
LDADD = ../src/libyaxt.la $(MPI_C_LIB)
