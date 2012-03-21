#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "mpi.h"

#include "xt_idxset.h"
#include "xt_xmap.h"
#include "xt_mpi.h"

struct xt_xmap {

   xt_idxset src_idxset_local, dst_idxset_local;

   xt_idxset *src_idxsets_all, *dst_idxsets_all;

   MPI_Comm comm;
   int comm_size;
};

static void exchange_idxsets(struct xt_xmap * xmap) {

   unsigned i, j;

   // compute size of local index sets
   size_t send_buffer_size;
   send_buffer_size  = xt_idxset_get_pack_size(xmap->src_idxset_local, xmap->comm);
   send_buffer_size += xt_idxset_get_pack_size(xmap->dst_idxset_local, xmap->comm);

   void * send_buffer;

   send_buffer = malloc(send_buffer_size);
   if (send_buffer == NULL) {
      fprintf(stderr, "malloc failed for send buffer of index set exchange\n");
      exit(EXIT_FAILURE);
   }

   // pack local index sets
   int position = 0;
   xt_idxset_pack(xmap->src_idxset_local, send_buffer, send_buffer_size, &position,
                  xmap->comm);
   xt_idxset_pack(xmap->dst_idxset_local, send_buffer, send_buffer_size, &position,
                  xmap->comm);

   // exchange buffer sizes
   int send_buffer_header = send_buffer_size;
   int * recv_buffer_header;

   recv_buffer_header = malloc(xmap->comm_size * sizeof(*recv_buffer_header));
   if (recv_buffer_header == NULL) {
      fprintf(stderr, "malloc failed for receive buffer header of index set exchange\n");
      exit(EXIT_FAILURE);
   }

   mpi_err_handler(MPI_Allgather(&send_buffer_header, 1, MPI_INT, recv_buffer_header,
                                 xmap->comm_size, MPI_INT, xmap->comm), xmap->comm);

   void * recv_buffer;
   size_t recv_buffer_size = 0;
   int * displ;

   displ = malloc(xmap->comm_size * sizeof(*displ));
   if (displ == NULL) {
      fprintf(stderr, "malloc failed for displacment of index set exchange\n");
      exit(EXIT_FAILURE);
   }

   for (i = 0; i < xmap->comm_size; ++i) {

      displ[i] = i ? (displ[i-1] + recv_buffer_header[i]) : 0;
      recv_buffer_size += recv_buffer_header[i];
   }

   recv_buffer = malloc(recv_buffer_size);
   if (recv_buffer == NULL) {
      fprintf(stderr, "malloc failed for receive buffer of index set exchange\n");
      exit(EXIT_FAILURE);
   }

   // exchange buffers
   mpi_err_handler(MPI_Allgatherv(send_buffer, send_buffer_size, MPI_PACKED,
                                  recv_buffer, recv_buffer_header, displ, MPI_PACKED,
                                  xmap->comm), xmap->comm);

   xmap->src_idxsets_all = malloc(xmap->comm_size * sizeof(*(xmap->src_idxsets_all)));
   if (xmap->src_idxsets_all == NULL) {
      fprintf(stderr, "malloc failed for source index sets all of index set exchange\n");
      exit(EXIT_FAILURE);
   }

   xmap->dst_idxsets_all = malloc(xmap->comm_size * sizeof(*(xmap->dst_idxsets_all)));
   if (xmap->dst_idxsets_all == NULL) {
      fprintf(stderr, "malloc failed for destination index sets all of index set exchange\n");
      exit(EXIT_FAILURE);
   }

   // unpack buffers
   for (i = 0, position = 0; i < xmap->comm_size; ++i) {

      xmap->src_idxsets_all[i] = xt_idxset_unpack(recv_buffer, recv_buffer_size,
                                                  &position, xmap->comm);
      xmap->dst_idxsets_all[i] = xt_idxset_unpack(recv_buffer, recv_buffer_size,
                                                  &position, xmap->comm);
   }

   // clean up
   free(send_buffer);
   free(recv_buffer_header);
   free(recv_buffer);
   free(displ);
}

struct xt_xmap*
xt_xmap_new(xt_idxset src_idxset, xt_idxset dst_idxset, MPI_Comm comm) {

   struct xt_xmap * xmap;

   xmap = malloc(1 * sizeof(*xmap));
   if (xmap == NULL) {
      fprintf(stderr, "malloc failed for exchange map\n");
      exit(EXIT_FAILURE);
   }

   xmap->src_idxset_local = xt_idxset_copy(src_idxset);
   xmap->dst_idxset_local = xt_idxset_copy(dst_idxset);

   mpi_err_handler(MPI_Comm_dup(comm, &(xmap->comm)), comm);
   mpi_err_handler(MPI_Comm_size(xmap->comm, &(xmap->comm_size)), xmap->comm);

   // exchange index sets between all processes in comm
   exchange_idxsets(xmap);

   return xmap;
}

void xt_xmap_delete(xt_xmap xmap) {

   unsigned i;

   for (i = 0; i < xmap->comm_size; ++i) {

      xt_idxset_delete(xmap->src_idxsets_all[i]);
      xt_idxset_delete(xmap->dst_idxsets_all[i]);
   }

   free(xmap->src_idxsets_all);
   free(xmap->dst_idxsets_all);

   mpi_err_handler(MPI_Comm_free(&(xmap->comm)), MPI_COMM_WORLD);
   free(xmap);
}
