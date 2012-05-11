noinst_PROGRAMS = test_idxvec					\
	test_idxlist_collection					\
	test_idxsection						\
	test_redist_p2p						\
	test_redist_p2p_parallel				\
	test_xmap


test_idxvec_SOURCES = test_idxvec.c tests.h
test_idxlist_collection_SOURCES = test_idxlist_collection.c tests.h
test_redist_p2p_SOURCES = test_redist_p2p.c test.h
test_redist_p2p_parallel_SOURCES = test_redist_p2p_parallel.c test.h
test_xmap_SOURCES = test_xmap.c tests.h

AM_CFLAGS = -I$(top_srcdir)/src
LDADD = ../src/libyaxt.a

TESTS = test_idxvec 						\
	test_idxlist_collection					\
	test_idxsection						\
	test_redist_p2p						\
	test_redist_p2p_parallel_run				\
	test_xmap 
