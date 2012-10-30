#include <stdlib.h>
#include <stdio.h>

#include "yaxt.h"
#include "xt/xt_idxlist.h"
#include "xt_idxvec.h"
#include "xt/xt_idxstripes.h"
#include "xt/xt_idxsection.h"
#include "xt/xt_idxlist_collection.h"

static void print_index_list(Xt_idxlist idxlist);

int main(void) {

  {
    //! [index list example]
    Xt_idx indices[] = {0,3,32,26,44,14,48};
    Xt_count num_indices = sizeof(indices)/sizeof(indices[0]);

    Xt_idxlist idxvec = xt_idxvec_new(indices, num_indices);
    //! [index list example]

    puts("index vector example:");
    print_index_list(idxvec);

    xt_idxlist_delete(idxvec);
  }

  {
    //! [index stripes example]
    struct Xt_stripe stripes[] = {{.start = 0, .nstrides = 17, .stride = 3}};
    Xt_count num_stripes = sizeof(stripes)/sizeof(stripes[0]);

    Xt_idxlist idxstripes = xt_idxstripes_new(stripes, num_stripes);
    //! [index stripes example]

    puts("index stripes example:");
    print_index_list(idxstripes);

    xt_idxlist_delete(idxstripes);
  }

  {
    //! [index section example]
    Xt_idx start = 0;
    unsigned num_dimensions = 2;
    unsigned global_size[2] = {5,10};
    unsigned local_size[2] = {3,4};
    unsigned local_start[2] = {2,3};

    Xt_idxlist idxsection
      = xt_idxsection_new(start, num_dimensions, global_size,
                          local_size, local_start);
    //! [index section example]

    puts("index section example:");
    print_index_list(idxsection);

    xt_idxlist_delete(idxsection);
  }

  {
    //! [index list collection example]
    Xt_idx start = 0;
    unsigned num_dimensions = 2;
    unsigned global_size[2] = {5,10};
    unsigned local_size[2][2] = {{2,2},{3,4}};
    unsigned local_start[2][2] = {{0,0},{2,6}};

    Xt_idxlist idxlists[]
      = {xt_idxsection_new(start, num_dimensions, global_size,
                           local_size[0], local_start[0]),
         xt_idxsection_new(start, num_dimensions, global_size,
                           local_size[1], local_start[1])};
    Xt_count num_idxlists = sizeof(idxlists) / sizeof(idxlists[0]);

    Xt_idxlist idxcollection
      = xt_idxlist_collection_new(idxlists, num_idxlists);
    //! [index list collection example]

    puts("index list collection example:");
    print_index_list(idxcollection);

    for (Xt_count i = 0; i < num_idxlists; ++i)
      xt_idxlist_delete(idxlists[i]);
    xt_idxlist_delete(idxcollection);
  }

  return 0;
}

static void print_index_list(Xt_idxlist idxlist) {

  for (Xt_count i = 0; i < xt_idxlist_get_num_indices(idxlist); ++i) {

    Xt_idx index;

    xt_idxlist_get_index_at_position(idxlist, i, &index);
    printf("%d ", index);
  }
  puts("");
}
