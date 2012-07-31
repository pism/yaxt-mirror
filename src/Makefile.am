AM_CPPFLAGS = -I../include -I$(srcdir)/../include
AM_FCFLAGS  = $(MPI_FC_INCLUDE)
lib_LIBRARIES = libyaxt.a 

libyaxt_a_SOURCES =						\
	xt_handles.c						\
	xt_handles.h						\
	xt_idxlist.c						\
	xt_idxlist.h						\
	xt_idxlist_collection.c					\
	xt_idxlist_collection.h					\
	xt_idxlist_unpack.c 					\
	xt_idxsection.c						\
	xt_idxsection.h						\
	xt_idxstripes.c						\
	xt_idxstripes.h						\
	xt_idxvec.c 						\
	xt_idxvec.h 						\
	xt_mpi.c						\
	xt_mpi.h						\
	xt_redist.c						\
	xt_redist.h						\
	xt_redist_p2p.c						\
	xt_redist_p2p.h						\
	xt_xmap.c						\
	xt_xmap.h						\
	xt_xmap_all2all.c					\
	xt_xmap_all2all.h					\
	xt_init.c						\
	xt_ut.f90						\
	xt_ut_c.c						\
	xt_ut_c.h						\
	core/core.c						\
	core/core.h						\
	core/ppm_xfuncs.h					\
	core/symprefix.h					\
	core/xmalloc.c						\
	core/xstdio.c						\
	ensure_array_size.c					\
	ensure_array_size.h					\
	quicksort.c						\
	quicksort.h						\
	mergesort.c						\
	mergesort.h						\
	sort_common.h						\
	xt_sort.c						\
	xt_sort.h						\
	yaxt_f2c.c						\
	yaxt.f90						\
	yaxt.h

clean-local:
	-rm -f *.$(FCMODEXT) i.*.L

