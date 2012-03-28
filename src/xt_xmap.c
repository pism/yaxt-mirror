#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "mpi.h"

#include "xt_idxvec.h"
#include "xt_xmap.h"
#include "xt_mpi.h"

struct xt_xmap {

   xt_idxvec src_idxvec_local, dst_idxvec_local;

   xt_idxvec *src_idxvecs_all, *dst_idxvecs_all;

   MPI_Comm comm;
   int comm_size;
};

static void exchange_idxvecs(struct xt_xmap * xmap) {

   unsigned i, j;

   // compute size of local index vectors
   size_t send_buffer_size;
   send_buffer_size  = xt_idxvec_get_pack_size(xmap->src_idxvec_local, xmap->comm);
   send_buffer_size += xt_idxvec_get_pack_size(xmap->dst_idxvec_local, xmap->comm);

   void * send_buffer;

   send_buffer = malloc(send_buffer_size);
   if (send_buffer == NULL) {
      fprintf(stderr, "malloc failed for send buffer of index vector exchange\n");
      exit(EXIT_FAILURE);
   }

   // pack local index vectors
   int position = 0;
   xt_idxvec_pack(xmap->src_idxvec_local, send_buffer, send_buffer_size, &position,
                  xmap->comm);
   xt_idxvec_pack(xmap->dst_idxvec_local, send_buffer, send_buffer_size, &position,
                  xmap->comm);

   // exchange buffer sizes
   int send_buffer_header = send_buffer_size;
   int * recv_buffer_header;

   recv_buffer_header = malloc(xmap->comm_size * sizeof(*recv_buffer_header));
   if (recv_buffer_header == NULL) {
      fprintf(stderr, "malloc failed for receive buffer header of index vector exchange\n");
      exit(EXIT_FAILURE);
   }

   mpi_err_handler(MPI_Allgather(&send_buffer_header, 1, MPI_INT, recv_buffer_header,
                                 xmap->comm_size, MPI_INT, xmap->comm), xmap->comm);

   void * recv_buffer;
   size_t recv_buffer_size = 0;
   int * displ;

   displ = malloc(xmap->comm_size * sizeof(*displ));
   if (displ == NULL) {
      fprintf(stderr, "malloc failed for displacment of index vector exchange\n");
      exit(EXIT_FAILURE);
   }

   for (i = 0; i < xmap->comm_size; ++i) {

      displ[i] = i ? (displ[i-1] + recv_buffer_header[i]) : 0;
      recv_buffer_size += recv_buffer_header[i];
   }

   recv_buffer = malloc(recv_buffer_size);
   if (recv_buffer == NULL) {
      fprintf(stderr, "malloc failed for receive buffer of index vector exchange\n");
      exit(EXIT_FAILURE);
   }

   // exchange buffers
   mpi_err_handler(MPI_Allgatherv(send_buffer, send_buffer_size, MPI_PACKED,
                                  recv_buffer, recv_buffer_header, displ, MPI_PACKED,
                                  xmap->comm), xmap->comm);

   xmap->src_idxvecs_all = malloc(xmap->comm_size * sizeof(*(xmap->src_idxvecs_all)));
   if (xmap->src_idxvecs_all == NULL) {
      fprintf(stderr, "malloc failed for source index vectors all of index vector exchange\n");
      exit(EXIT_FAILURE);
   }

   xmap->dst_idxvecs_all = malloc(xmap->comm_size * sizeof(*(xmap->dst_idxvecs_all)));
   if (xmap->dst_idxvecs_all == NULL) {
      fprintf(stderr, "malloc failed for destination index vectors all of index vector exchange\n");
      exit(EXIT_FAILURE);
   }

   // unpack buffers
   for (i = 0, position = 0; i < xmap->comm_size; ++i) {

      xmap->src_idxvecs_all[i] = xt_idxvec_unpack(recv_buffer, recv_buffer_size,
                                                  &position, xmap->comm);
      xmap->dst_idxvecs_all[i] = xt_idxvec_unpack(recv_buffer, recv_buffer_size,
                                                  &position, xmap->comm);
   }

   // clean up
   free(send_buffer);
   free(recv_buffer_header);
   free(recv_buffer);
   free(displ);
}

struct xt_xmap*
xt_xmap_new(xt_idxvec src_idxvec, xt_idxvec dst_idxvec, MPI_Comm comm) {

   struct xt_xmap * xmap;

   xmap = malloc(1 * sizeof(*xmap));
   if (xmap == NULL) {
      fprintf(stderr, "malloc failed for exchange map\n");
      exit(EXIT_FAILURE);
   }

   xmap->src_idxvec_local = xt_idxvec_copy(src_idxvec);
   xmap->dst_idxvec_local = xt_idxvec_copy(dst_idxvec);

   mpi_err_handler(MPI_Comm_dup(comm, &(xmap->comm)), comm);
   mpi_err_handler(MPI_Comm_size(xmap->comm, &(xmap->comm_size)), xmap->comm);

   // exchange index vectors between all processes in comm
   exchange_idxvecs(xmap);

   return xmap;
}

void xt_xmap_delete(xt_xmap xmap) {

   unsigned i;

   for (i = 0; i < xmap->comm_size; ++i) {

      xt_idxvec_delete(xmap->src_idxvecs_all[i]);
      xt_idxvec_delete(xmap->dst_idxvecs_all[i]);
   }

   free(xmap->src_idxvecs_all);
   free(xmap->dst_idxvecs_all);

   mpi_err_handler(MPI_Comm_free(&(xmap->comm)), MPI_COMM_WORLD);
   free(xmap);
}
