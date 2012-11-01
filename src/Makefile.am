AM_CPPFLAGS = -I../include -I$(srcdir)/../include $(MPI_C_INCLUDE)
AM_FCFLAGS  = $(MPI_FC_INCLUDE)
lib_LTLIBRARIES = libyaxt.la

# headers from subdir xt will be installed (after substituting
# definitions from config.h). These must be listed in a sequence in
# which they can be combined in yaxt.h.
libyaxt_la_SOURCES =						\
	xt/xt_core.h						\
	xt/xt_handles.h						\
	xt/sort_common.h					\
	xt/mergesort.h						\
	xt/quicksort.h						\
	xt/xt_sort.h						\
	xt/xt_idxlist.h						\
	xt/xt_idxlist_collection.h				\
	xt/xt_idxsection.h					\
	xt/xt_idxstripes.h					\
	xt/xt_idxvec.h 						\
	xt/xt_mpi.h						\
	xt/xt_redist.h						\
	xt/xt_redist_collection.h				\
	xt/xt_redist_collection_static.h			\
	xt/xt_redist_p2p.h					\
	xt/xt_xmap.h						\
	xt/xt_xmap_all2all.h					\
	xt/xt_xmap_dist_dir.h					\
	xt/xt_xmap_intersection.h				\
	xt/xt_idxmod.h						\
	xt_handles.c						\
	xt_idxlist.c						\
	xt_idxlist_collection.c					\
	xt_idxlist_unpack.c 					\
	xt_idxmod.c						\
	xt_idxsection.c						\
	xt_idxstripes.c						\
	xt_idxvec.c 						\
	xt_mpi.c						\
	xt_redist.c						\
	xt_redist_collection.c					\
	xt_redist_collection_static.c				\
	xt_redist_p2p.c						\
	xt_xmap.c						\
	xt_xmap_all2all.c					\
	xt_xmap_dist_dir.c					\
	xt_xmap_intersection.c					\
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
	mergesort.c						\
	xt_sort.c						\
	yaxt_f2c.c						\
	yaxt.f90

xtinclude_Makefile.inc: Makefile
	headers=`echo $(libyaxt_la_SOURCES) | sed 's/[ 	][ 	]*/\n/g' | grep 'xt/.*\.h$$'` ; \
	  echo "xtinclude_HEADERS = \\" >$@ ; \
	  for header in $$headers ; do echo "	$$header \\" ; done \
	    | sed '$$s/\\$$//' >>$@

all-local: xtinclude_Makefile.inc

clean-local:
	-rm -f *.$(FCMODEXT) i.*.L

install-exec-local:
	$(mkinstalldirs) "$(DESTDIR)$(libdir)/pkgconfig"
	$(install_sh_DATA) pkgconfig/yaxt.pc \
		"$(DESTDIR)$(libdir)/pkgconfig/yaxt.pc"

uninstall-local:
	rm -f "$(DESTDIR)$(libdir)/pkgconfig/yaxt.pc"
	-rmdir "$(DESTDIR)$(libdir)/pkgconfig"
