#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "mpi.h"

#include "xt_idxset.h"
#include "xt_xmap.h"

struct xt_xmap {

   xt_idxset* src_idxsets_local, * dst_idxsets_local;
   unsigned num_src_idxsets_local, num_dst_idxsets_local;

   xt_idxset** dst_idxsets_all;
   unsigned * num_dst_idxsets_per_process;

   MPI_Comm comm;
};

struct xt_xmap*
xt_xmap_new_single(xt_idxset src_idxset, xt_idxset dst_idxset, MPI_Comm comm) {

   xt_xmap_new(&src_idxset, 1, &dst_idxset, 1, comm);
}

struct xt_xmap*
xt_xmap_new(xt_idxset * src_idxsets, unsigned num_src_idxsets,
            xt_idxset * dst_idxsets, unsigned num_dst_idxsets,
            MPI_Comm comm) {

   struct xt_xmap * xmap;

   xmap = malloc(1 * sizeof(*xmap));
   if (xmap == NULL) {
      fprintf(stderr, "malloc failed for exchange map\n");
      exit(EXIT_FAILURE);
   }

   if (num_src_idxsets > 0) {
   
      xmap->src_idxsets_local = malloc(num_src_idxsets * sizeof(*(xmap->src_idxsets_local)));

      if (xmap->src_idxsets_local == NULL) {
         fprintf(stderr, "malloc failed for local index sets\n");
         exit(EXIT_FAILURE);
      }

      memcpy(xmap->src_idxsets_local, src_idxsets, num_src_idxsets * sizeof(*(xmap->src_idxsets_local)));

   } else {
      xmap->src_idxsets_local = NULL;
   }

   if (num_dst_idxsets > 0) {
   
      xmap->dst_idxsets_local = malloc(num_dst_idxsets * sizeof(*(xmap->dst_idxsets_local)));

      if (xmap->src_idxsets_local == NULL) {
         fprintf(stderr, "malloc failed for local index sets\n");
         exit(EXIT_FAILURE);
      }

      memcpy(xmap->dst_idxsets_local, dst_idxsets, num_dst_idxsets * sizeof(*(xmap->dst_idxsets_local)));

   } else {
      xmap->dst_idxsets_local = NULL;
   }

   xmap->num_src_idxsets_local = num_src_idxsets;
   xmap->num_dst_idxsets_local = num_dst_idxsets;
   xmap->comm = comm;

   xmap->dst_idxsets_all = NULL;
   xmap->num_dst_idxsets_per_process = NULL;

   return xmap;
}

void xt_xmap_delete(xt_xmap xmap) {

   free(xmap->src_idxsets_local);
   free(xmap->dst_idxsets_local);

   free(xmap);
}
