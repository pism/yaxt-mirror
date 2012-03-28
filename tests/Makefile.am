noinst_PROGRAMS = test_idxlist					\
	test_idxvec_collection					\
	test_xmap


test_idxlist_SOURCES = test_idxlist.c tests.h
test_idxvec_collection_SOURCES = test_idxvec_collection.c tests.h
test_xmap_SOURCES = test_xmap.c tests.h

AM_CFLAGS = -I$(top_srcdir)/src
LDADD = ../src/libyaxt.a

TESTS = $(noinst_PROGRAMS)
