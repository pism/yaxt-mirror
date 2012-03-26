noinst_PROGRAMS = test_idxlist					\
	test_idxset_collection					\
	test_xmap


test_idxlist_SOURCES = test_idxlist.c tests.h
test_idxset_collection_SOURCES = test_idxset_collection.c tests.h
test_xmap_SOURCES = test_xmap.c tests.h

AM_CFLAGS = -I$(top_srcdir)/src
LDADD = ../src/libyaxt.a

TESTS = $(noinst_PROGRAMS)
