#include <stdlib.h>
#include <assert.h>

#include <mpi.h>

#include <yaxt.h>

#include "tests.h"
#include "test_idxlist_utils.h"
#include "core/ppm_xfuncs.h"

int main(void) {

  // init:
  mpi_err_handler(MPI_Init(NULL, NULL), MPI_COMM_WORLD);
  xt_initialize(MPI_COMM_WORLD);

  { // idxvec modifier

    Xt_idx g_src_idx[] = {1,2,3,4,5,6,7,8,9};
    int g_src_num = sizeof(g_src_idx) / sizeof(g_src_idx[0]);
    Xt_idxlist g_src_idxlist = xt_idxvec_new(g_src_idx, g_src_num);

    Xt_idx g_dst_idx[] = {9,8,7,6,5,4,3,2,1};
    int g_dst_num = sizeof(g_dst_idx) / sizeof(g_dst_idx[0]);
    Xt_idxlist g_dst_idxlist = xt_idxvec_new(g_dst_idx, g_dst_num);

    struct Xt_modifier mod = {.extract=g_src_idxlist, .subst=g_dst_idxlist, .mask=0};

    Xt_idx patch_idx[] = {3,4,4,4,7,7,8};
    int patch_num = sizeof(patch_idx) / sizeof(patch_idx[0]);
    Xt_idxlist patch_idxlist = xt_idxvec_new(patch_idx, patch_num);

    // idx:{3,4,4,4,7,7,8} -> pos:{2,3,3,3,6,6,7} => idx:{7,6,6,6,3,3,2}
    Xt_idx ref_mpatch_idx[] = {7,6,6,6,3,3,2};

    Xt_idxlist mpatch_idxlist = xt_idxmod_new(patch_idxlist, &mod, 1, NULL);

    check_idxlist(mpatch_idxlist, ref_mpatch_idx, patch_num);

    xt_idxlist_delete(mpatch_idxlist);
    xt_idxlist_delete(patch_idxlist);
    xt_idxlist_delete(g_dst_idxlist);
    xt_idxlist_delete(g_src_idxlist);
  }

  { // idxstripes modifier

    struct Xt_stripe g_src_stripe[] = { {.start=1, .nstrides=20, .stride=1} };
    int g_src_num = sizeof(g_src_stripe) / sizeof(g_src_stripe[0]);
    Xt_idxlist g_src_idxlist = xt_idxstripes_new(g_src_stripe, g_src_num);

    struct Xt_stripe g_dst_stripe[] = { {.start=100, .nstrides=20, .stride=-1} };
    int g_dst_num = sizeof(g_dst_stripe) / sizeof(g_dst_stripe[0]);
    Xt_idxlist g_dst_idxlist = xt_idxstripes_new(g_dst_stripe, g_dst_num);

    struct Xt_modifier mod = {.extract=g_src_idxlist, .subst=g_dst_idxlist, .mask=32};

    Xt_idx patch_idx[] = {0,1,3,3,5,50,100,150};
    int patch_num = sizeof(patch_idx) / sizeof(patch_idx[0]);
    Xt_idxlist patch_idxlist = xt_idxvec_new(patch_idx, patch_num);

    // inter:{1,3,3,5} => extract_pos:{0,2,2,4} => subst_idx:{100,98,98,96}, patch_pos:{1,2,3,4} = > mpatch:{0,100,98,98,96,50,100,150}
    Xt_idx ref_mpatch_idx[] = {0,100,98,98,96,50,100,150};

    int mstate[]   = {1, 2, 3, 4, 5, 6, 7, 8};
    int ref_mstate[] = {1, 2|32, 3|32, 4|32, 5|32, 6, 7, 8};

    Xt_idxlist mpatch_idxlist = xt_idxmod_new(patch_idxlist, &mod, 1, mstate);

    check_idxlist(mpatch_idxlist, ref_mpatch_idx, patch_num);

    // check mstate:
    for(int i=0; i< patch_num; i++) {
      if (mstate[i] != ref_mstate[i])
        PUT_ERR("(mstate[i] != ref_mstate[i])\n");
    }

    xt_idxlist_delete(mpatch_idxlist);
    xt_idxlist_delete(patch_idxlist);
    xt_idxlist_delete(g_dst_idxlist);
    xt_idxlist_delete(g_src_idxlist);

  }

  { // multiple modifier, track modifier usage

    Xt_idx g1_src_idx[] = {1,2,3,4,5,6,7,8,9};
    int g1_src_num = sizeof(g1_src_idx) / sizeof(g1_src_idx[0]);
    Xt_idxlist g1_src_idxlist = xt_idxvec_new(g1_src_idx, g1_src_num);

    Xt_idx g1_dst_idx[] = {9,8,7,6,5,4,3,2,1};
    int g1_dst_num = sizeof(g1_dst_idx) / sizeof(g1_dst_idx[0]);
    Xt_idxlist g1_dst_idxlist = xt_idxvec_new(g1_dst_idx, g1_dst_num);

    Xt_idx g2_src_idx[] = {1,2,8,9,10};
    int g2_src_num = sizeof(g2_src_idx) / sizeof(g2_src_idx[0]);
    Xt_idxlist g2_src_idxlist = xt_idxvec_new(g2_src_idx, g2_src_num);

    Xt_idx g2_dst_idx[] = {8,2,8,2,5};
    int g2_dst_num = sizeof(g2_dst_idx) / sizeof(g2_dst_idx[0]);
    Xt_idxlist g2_dst_idxlist = xt_idxvec_new(g2_dst_idx, g2_dst_num);

    struct Xt_modifier mod[2] = { {.extract=g1_src_idxlist, .subst=g1_dst_idxlist, .mask=1},
                                  {.extract=g2_src_idxlist, .subst=g2_dst_idxlist, .mask=2} };

    Xt_idx patch_idx[] = {6,7,25,8,9,10};
    int patch_num = sizeof(patch_idx) / sizeof(patch_idx[0]);
    Xt_idxlist patch_idxlist = xt_idxvec_new(patch_idx, patch_num);
    int ref_mstate[] = {1|0, 1|0, 0|0, 1|2, 1|2, 0|2};
    int mstate[patch_num];
    // reset mstate:
    for (int i=0; i<patch_num; i++) {
      mstate[i] = 0;
    }

    // mod1: idx:{6,7,25,8,9,10} -> pos:{5,6,nil,7,8,nil} => idx:{4,3,25,2,1,10}
    // mod2: idx:{4,3,25,2,1,10} -> pos:{nil,nil,nil,1,0,4} => idx:{4,3,25,2,8,5}
    Xt_idx ref_mpatch_idx[] = {4,3,25,2,8,5};

    Xt_idxlist mpatch_idxlist = xt_idxmod_new(patch_idxlist, mod, 2, mstate);

    check_idxlist(mpatch_idxlist, ref_mpatch_idx, patch_num);

    // check mstate:
    for(int i=0; i< patch_num; i++) {
      if (mstate[i] != ref_mstate[i])
        PUT_ERR("(mstate[i] != ref_mstate[i])\n");
    }

    xt_idxlist_delete(mpatch_idxlist);
    xt_idxlist_delete(patch_idxlist);
    xt_idxlist_delete(g2_dst_idxlist);
    xt_idxlist_delete(g2_src_idxlist);
    xt_idxlist_delete(g1_dst_idxlist);
    xt_idxlist_delete(g1_src_idxlist);
  }

  // finalize:
  MPI_Finalize();
  return TEST_EXIT_CODE;
}
