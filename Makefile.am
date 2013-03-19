SUBDIRS = . src inst_headers tests examples perf

ACLOCAL_AMFLAGS=-I m4

EXTRA_DIST = util/mpi_launch_poe util/sunf95preproc-wrapper \
	util/sxpreproc-wrapper util/xlfpreproc-wrapper \
	doc/Doxyfile doc/README.build doc/cfortran.doc \
	doc/cfortran.html doc/cfortran.txt doc/styleguide.txt \
	doc/images/decomp_idxlist_collection.png \
	doc/images/decomp_idxlist_collection.svg \
	doc/images/decomp_idxmod.png \
	doc/images/decomp_idxsection.png \
	doc/images/decomp_idxsection.svg \
	doc/images/decomp_idxstripes.png \
	doc/images/decomp_idxstripes.svg \
	doc/images/decomp_idxvec.png \
	doc/images/decomp_idxvec.svg \
	doc/images/neg_size_idxsection_nnnn.png \
	doc/images/neg_size_idxsection_nnnp.png \
	doc/images/neg_size_idxsection_nnpn.png \
	doc/images/neg_size_idxsection_nnpp.png \
	doc/images/neg_size_idxsection_npnn.png \
	doc/images/neg_size_idxsection_npnp.png \
	doc/images/neg_size_idxsection_nppn.png \
	doc/images/neg_size_idxsection_nppp.png \
	doc/images/neg_size_idxsection_pnnn.png \
	doc/images/neg_size_idxsection_pnnp.png \
	doc/images/neg_size_idxsection_pnpn.png \
	doc/images/neg_size_idxsection_pnpp.png \
	doc/images/neg_size_idxsection_ppnn.png \
	doc/images/neg_size_idxsection_ppnp.png \
	doc/images/neg_size_idxsection_pppn.png \
	doc/images/neg_size_idxsection_pppp.png \
	doc/images/stride_idxsection_x.png \
	doc/images/stride_idxsection_xy.png \
	doc/src/idxlists.dox \
	doc/src/idxsection.dox \
	doc/src/mainpage.dox \
	doc/src/redist_colls.dox \
	doc/src/rrobin.dox \
	doc/src/xmaps.dox \
	doc/src/yaxt_concept.dox \
	include/cfortran.h \
	m4/tj_find_type.m4

all-local:
	for subdir in $(SUBDIRS) ; do \
          test "$$subdir/Makefile" = . -o "$$subdir/Makefile"  -nt "$(top_srcdir)/$$subdir/Makefile.in" \
	    || ./config.status $$subdir/Makefile ; done

clean-local:
	for subdir in $(SUBDIRS) ; do \
          test "$$subdir/Makefile" = . -o "$$subdir/Makefile"  -nt "$(top_srcdir)/$$subdir/Makefile.in" \
	    || ./config.status $$subdir/Makefile ; done
