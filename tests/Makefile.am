#
# Keywords:
# Maintainer: Jörg Behrens <behrens@dkrz.de>
#             Moritz Hanke <hanke@dkrz.de>
#             Thomas Jahns <jahns@dkrz.de>
# URL: https://redmine.dkrz.de/doc/yaxt/html/index.html
# 
# Redistribution and use in source and binary forms, with or without
# modification, are  permitted provided that the following conditions are
# met:
# 
# Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
# 
# Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
# 
# Neither the name of the DKRZ GmbH nor the names of its contributors
# may be used to endorse or promote products derived from this software
# without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

noinst_PROGRAMS = 						\
	test_mpi_generate_datatype				\
	test_handles						\
	test_idxempty						\
	test_idxvec						\
	test_idxlist_collection					\
	test_idxmod						\
	test_idxsection						\
	test_idxstripes						\
	test_redist_collection					\
	test_redist_collection_parallel				\
	test_redist_collection_static				\
	test_redist_collection_static_parallel			\
	test_redist_p2p						\
	test_redist_p2p_parallel				\
	test_xmap_all2all					\
	test_xmap_all2all_parallel				\
	test_xmap_dist_dir					\
	test_xmap_dist_dir_parallel				\
	test_ut							\
	test_perf_idxsection_get_positions_of_indices		\
	test_perf						\
	test_perf_stripes					\
	test_sort						\
	test_yaxt


test_idxempty_SOURCES = test_idxempty.c tests.h test_idxlist_utils.h
test_idxvec_SOURCES = test_idxvec.c tests.h test_idxlist_utils.h
test_idxlist_collection_SOURCES = test_idxlist_collection.c tests.h \
	test_idxlist_utils.h
test_idxmod_SOURCES = test_idxmod.c tests.h
test_idxsection_SOURCES = test_idxsection.c tests.h test_idxlist_utils.h
test_idxstripes_SOURCES = test_idxstripes.c tests.h test_idxlist_utils.h
test_redist_collection_SOURCES = test_redist_collection.c tests.h
test_redist_collection_parallel_SOURCES = test_redist_collection_parallel.c \
	tests.h
test_redist_collection_static_SOURCES = test_redist_collection_static.c tests.h
test_redist_collection_static_parallel_SOURCES = \
	test_redist_collection_static_parallel.c tests.h
test_redist_p2p_SOURCES = test_redist_p2p.c tests.h
test_redist_p2p_parallel_SOURCES = test_redist_p2p_parallel.c tests.h
test_xmap_all2all_SOURCES = test_xmap_all2all.c tests.h
test_xmap_all2all_parallel_SOURCES = test_xmap_all2all_parallel.c tests.h
test_xmap_dist_dir_SOURCES = test_xmap_dist_dir.c tests.h
test_xmap_dist_dir_parallel_SOURCES = test_xmap_dist_dir_parallel.c tests.h
test_handles_SOURCES = test_handles.c tests.h
test_ut_SOURCES = test_ut.f90
test_ut_LDADD = ../src/libyaxt.la $(MPI_FC_LIB)
test_perf_idxsection_get_positions_of_indices_SOURCES = \
	test_perf_idxsection_get_positions_of_indices.c
test_perf_SOURCES = test_perf.f90
test_perf_LDADD = ../src/libyaxt.la $(MPI_FC_LIB)
test_perf_stripes_SOURCES = test_perf_stripes.f90
test_perf_stripes.$(OBJEXT): ../src/yaxt.$(OBJEXT)
test_perf_stripes_LDADD = ../src/libyaxt.la $(MPI_FC_LIB)
test_sort_SOURCES = test_sort.c tests.h
test_yaxt_SOURCES = test_yaxt.f90
test_yaxt_LDADD = ../src/libyaxt.la $(MPI_FC_LIB)
test_mpi_generate_datatype_SOURCES = test_mpi_generate_datatype.c tests.h

AM_CFLAGS = -I../inst_headers $(MPI_C_INCLUDE)
AM_FCFLAGS = -I$(top_srcdir)/src -I../src $(MPI_FC_INCLUDE) $(MPI_FC_MOD)
LDADD = ../src/libyaxt.la $(MPI_C_LIB)

TESTS = 							\
	test_mpi_generate_datatype_run				\
	test_handles_run					\
	test_idxempty_run					\
	test_idxvec_run						\
	test_idxlist_collection_run				\
	test_idxmod_run						\
	test_idxsection_run					\
	test_idxstripes_run					\
	test_redist_p2p_run					\
	test_redist_collection_run				\
	test_redist_collection_parallel_run			\
	test_redist_collection_static_run			\
	test_redist_collection_static_parallel_run		\
	test_redist_p2p_parallel_run				\
	test_xmap_all2all_run					\
	test_xmap_all2all_parallel_run				\
	test_xmap_dist_dir_run					\
	test_xmap_dist_dir_parallel_run				\
	test_ut_run						\
	test_perf_idxsection_get_positions_of_indices_run	\
	test_perf_run						\
	test_perf_stripes_run					\
	test_sort_run						\
	test_yaxt_run


AUTOMAKE_OPTIONS = color-tests
