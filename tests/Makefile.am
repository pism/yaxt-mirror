noinst_PROGRAMS = 						\
	test_handles						\
	test_idxvec						\
	test_idxlist_collection					\
	test_idxsection						\
	test_idxstripes						\
	test_redist_p2p						\
	test_redist_p2p_parallel				\
	test_xmap						\
	test_xmap_parallel					\
	test_ut							\
	test_perf


test_idxvec_SOURCES = test_idxvec.c tests.h
test_idxlist_collection_SOURCES = test_idxlist_collection.c tests.h
test_idxsection_SOURCES = test_idxsection.c tests.h
test_idxstripes_SOURCES = test_idxstripes.c tests.h
test_redist_p2p_SOURCES = test_redist_p2p.c test.h
test_redist_p2p_parallel_SOURCES = test_redist_p2p_parallel.c test.h
test_xmap_SOURCES = test_xmap.c tests.h
test_xmap_parallel_SOURCES = test_xmap_parallel.c tests.h
test_handles_SOURCES = test_handles.c tests.h
test_ut_SOURCES = test_ut.f90
test_ut_LDADD = $(LDADD) $(MPI_FC_LIB)
test_perf_SOURCES = test_perf.f90

AM_CFLAGS = -I$(top_srcdir)/src
AM_FCFLAGS = -I$(top_srcdir)/src -I../src $(MPI_FC_INCLUDE)
LDADD = ../src/libyaxt.a

TESTS = 							\
	test_handles_run					\
	test_idxvec_run						\
	test_idxlist_collection_run				\
	test_idxsection_run					\
	test_idxstripes_run					\
	test_redist_p2p_run					\
	test_redist_p2p_parallel_run				\
	test_xmap_run						\
	test_xmap_parallel_run					\
	test_ut_run						\
	test_perf_run
