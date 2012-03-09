#include <stdlib.h>
#include <stdio.h>

#include "mpi.h"

#include "xt_idxset.h"
#include "xt_xmap.h"

struct xt_xmap {

   xt_idxset src_idxset, dst_idxset;
   MPI_Comm comm;
};

struct xt_xmap*
xt_xmap_new(xt_idxset src_idxset, xt_idxset dst_idxset, MPI_Comm comm) {

   struct xt_xmap * xmap;

   xmap = malloc(1 * sizeof(*xmap));
   if (xmap == NULL) {
      fprintf(stderr, "malloc failed for exchange map\n");
      exit(EXIT_FAILURE);
   }

   xmap->src_idxset = src_idxset;
   xmap->dst_idxset = dst_idxset;
   xmap->comm = comm;

   return xmap;
}

void xt_xmap_delete(xt_xmap xmap) {

   free(xmap);
}
