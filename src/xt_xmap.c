#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "mpi.h"

#include "xt_idxvec.h"
#include "xt_xmap.h"
#include "xt_mpi.h"

struct xt_xmap {
  struct {
    // list of global indices to send or receive from partner rank
    Xt_idxvec intersect;
    int rank;
  } *dst_msg, *src_msg;
  int ndst, nsrc;
  MPI_Comm comm;
};

static void exchange_idxvecs(struct xt_xmap * xmap,
                             xt_idxvec src_idxvec_local,
                             xt_idxvec dst_idxvec_local)
{

   xt_idxvec *src_idxvecs_all, *dst_idxvecs_all;

   unsigned i, j;
   int comm_size, comm_rank;

   mpi_err_handler(MPI_Comm_size(xmap->comm, &comm_size), xmap->comm);

   // compute size of local index vectors
   size_t send_buffer_size;
   send_buffer_size  = xt_idxvec_get_pack_size(src_idxvec_local, xmap->comm);
   send_buffer_size += xt_idxvec_get_pack_size(dst_idxvec_local, xmap->comm);

   void * send_buffer;

   send_buffer = malloc(send_buffer_size);
   if (send_buffer == NULL) {
     fputs("malloc failed for send buffer of index vector exchange\n", stderr);
      exit(EXIT_FAILURE);
   }

   // pack local index vectors
   int position = 0;
   xt_idxvec_pack(src_idxvec_local, send_buffer, send_buffer_size, &position,
                  xmap->comm);
   xt_idxvec_pack(dst_idxvec_local, send_buffer, send_buffer_size, &position,
                  xmap->comm);

   // exchange buffer sizes
   int send_buffer_header = send_buffer_size;
   int * recv_buffer_header;

   recv_buffer_header = malloc(comm_size * sizeof(*recv_buffer_header));
   if (recv_buffer_header == NULL) {
     fputs("malloc failed for receive buffer header"
           " of index vector exchange\n", stderr);
     exit(EXIT_FAILURE);
   }

   mpi_err_handler(MPI_Allgather(&send_buffer_header, 1, MPI_INT,
                                 recv_buffer_header,
                                 1, MPI_INT, xmap->comm), xmap->comm);

   void * recv_buffer;
   int * displ;

   displ = malloc(comm_size * sizeof(*displ));
   if (displ == NULL) {
     fputs("malloc failed for displacment of index vector exchange\n", stderr);
     exit(EXIT_FAILURE);
   }

   displ[0] = 0;
   size_t recv_buffer_size = recv_buffer_header[0];
   for (i = 1; i < comm_size; ++i) {

      displ[i] = displ[i-1] + recv_buffer_header[i-1];
      recv_buffer_size += recv_buffer_header[i];
   }

   recv_buffer = malloc(recv_buffer_size);
   if (recv_buffer == NULL) {
     fputs("malloc failed for receive buffer of index vector exchange\n",
           stderr);
     exit(EXIT_FAILURE);
   }

   // exchange buffers
   mpi_err_handler(MPI_Allgatherv(send_buffer, send_buffer_size, MPI_PACKED,
                                  recv_buffer, recv_buffer_header, displ,
                                  MPI_PACKED,
                                  xmap->comm), xmap->comm);

   src_idxvecs_all = malloc(comm_size * sizeof(*src_idxvecs_all));
   if (src_idxvecs_all == NULL) {
     fputs("malloc failed for source index vectors"
           " all of index vector exchange\n", stderr);
     exit(EXIT_FAILURE);
   }

   dst_idxvecs_all = malloc(comm_size * sizeof(*dst_idxvecs_all));
   if (dst_idxvecs_all == NULL) {
     fputs("malloc failed for destination index vectors"
           " all of index vector exchange\n", stderr);
     exit(EXIT_FAILURE);
   }

   // unpack buffers
   for (i = 0, position = 0; i < comm_size; ++i) {

      src_idxvecs_all[i] = xt_idxvec_unpack(recv_buffer, recv_buffer_size,
                                            &position, xmap->comm);
      dst_idxvecs_all[i] = xt_idxvec_unpack(recv_buffer, recv_buffer_size,
                                            &position, xmap->comm);
   }

   // allocate memory for intersections
   xmap->dst_msg = malloc(comm_size * sizeof (*xmap->dst_msg));
   xmap->ndst = 0;
   xmap->src_msg = malloc(comm_size * sizeof (*xmap->src_msg));
   xmap->nsrc = 0;
   if (!xmap->dst_msg || !xmap->src_msg)
     abort();


   // compute intersections
   for (i = 0; i < comm_size; ++i)
   {
     Xt_idxvec intersect
       = xt_idxvec_get_intersection(src_idxvecs_all[i], dst_idxvec_local);
     if (xt_idxvec_get_num_indices(intersect) > 0)
     {
       xmap->dst_msg[xmap->ndst].intersect = intersect;
       xmap->dst_msg[xmap->ndst].rank = i;
       xmap->ndst++;
     }
     else
       xt_idxvec_delete(intersect);
     intersect = xt_idxvec_get_intersection(src_idxvec_local,
                                            dst_idxvecs_all[i]);
     if (xt_idxvec_get_num_indices(intersect) > 0)
     {
       xmap->src_msg[xmap->nsrc].rank = i;
       xmap->src_msg[xmap->nsrc].intersect = intersect;
       xmap->nsrc++;
     }
     else
       xt_idxvec_delete(intersect);
   }

   for (i = 0; i < comm_size; ++i) {

      xt_idxvec_delete(src_idxvecs_all[i]);
      xt_idxvec_delete(dst_idxvecs_all[i]);
   }

   free(src_idxvecs_all);
   free(dst_idxvecs_all);

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


   mpi_err_handler(MPI_Comm_dup(comm, &(xmap->comm)), comm);

   // exchange index vectors between all processes in comm
   exchange_idxvecs(xmap, src_idxvec, dst_idxvec);

   return xmap;
}

void xt_xmap_delete(xt_xmap xmap) {

   unsigned i;


   mpi_err_handler(MPI_Comm_free(&(xmap->comm)), MPI_COMM_WORLD);
   free(xmap);
}
