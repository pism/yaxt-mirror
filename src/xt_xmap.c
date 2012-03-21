#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "mpi.h"

#include "xt_idxset.h"
#include "xt_xmap.h"
#include "xt_mpi.h"

struct xt_xmap {

   xt_idxset* src_idxsets_local, * dst_idxsets_local;
   unsigned num_src_idxsets_local, num_dst_idxsets_local;

   xt_idxset** src_idxsets_all, ** dst_idxsets_all;
   unsigned* num_src_idxsets_per_process, * num_dst_idxsets_per_process;

   MPI_Comm comm;
   int comm_size;
};

static void exchange_idxsets(struct xt_xmap * xmap) {

   unsigned i, j;

   // compute size of local index sets
   size_t send_buffer_size = 0;
   for (i = 0; i < xmap->num_src_idxsets_local; ++i)
      send_buffer_size += xt_idxset_get_pack_size(xmap->src_idxsets_local[i], xmap->comm);
   for (i = 0; i < xmap->num_dst_idxsets_local; ++i)
      send_buffer_size += xt_idxset_get_pack_size(xmap->dst_idxsets_local[i], xmap->comm);

   void * send_buffer;

   send_buffer = malloc(send_buffer_size);
   if (send_buffer == NULL) {
      fprintf(stderr, "malloc failed for send buffer of index set exchange\n");
      exit(EXIT_FAILURE);
   }

   // pack local index sets
   int position = 0;
   for (i = 0; i < xmap->num_src_idxsets_local; ++i)
      xt_idxset_pack(xmap->src_idxsets_local[i], send_buffer, send_buffer_size,
                     &position, xmap->comm);
   for (i = 0; i < xmap->num_dst_idxsets_local; ++i)
      xt_idxset_pack(xmap->dst_idxsets_local[i], send_buffer, send_buffer_size,
                     &position, xmap->comm);

   // exchange buffer sizes and number of source and destination index sets
   int send_buffer_header[3] = {send_buffer_size,
                                xmap->num_src_idxsets_local,
                                xmap->num_dst_idxsets_local};
   int * recv_buffer_header;

   recv_buffer_header = malloc(xmap->comm_size * 3 * sizeof(*recv_buffer_header));
   if (recv_buffer_header == NULL) {
      fprintf(stderr, "malloc failed for receive buffer header of index set exchange\n");
      exit(EXIT_FAILURE);
   }

   mpi_err_handler(MPI_Allgather(send_buffer_header, 3, MPI_INT, recv_buffer_header,
                                 3 * xmap->comm_size, MPI_INT, xmap->comm),
                   xmap->comm);

   xmap->num_src_idxsets_per_process = malloc(xmap->comm_size * sizeof(*(xmap->num_src_idxsets_per_process)));
   xmap->num_dst_idxsets_per_process = malloc(xmap->comm_size * sizeof(*(xmap->num_dst_idxsets_per_process)));

   void * recv_buffer;
   size_t recv_buffer_size = 0;
   int * displ;
   int * recvcount;

   displ = malloc(xmap->comm_size * sizeof(*displ));
   if (displ == NULL) {
      fprintf(stderr, "malloc failed for displacment of index set exchange\n");
      exit(EXIT_FAILURE);
   }
   recvcount = malloc(xmap->comm_size * sizeof(*recvcount));
   if (recvcount == NULL) {
      fprintf(stderr, "malloc failed for receive count of index set exchange\n");
      exit(EXIT_FAILURE);
   }

   for (i = 0; i < xmap->comm_size; ++i) {

      displ[i] = i ? (displ[i-1] + recv_buffer_header[3*i+0]) : 0;
      recv_buffer_size +=                    recv_buffer_header[3*i+0];
      recvcount[i] =                         recv_buffer_header[3*i+0];
      xmap->num_src_idxsets_per_process[i] = recv_buffer_header[3*i+1];
      xmap->num_dst_idxsets_per_process[i] = recv_buffer_header[3*i+2];
   }

   recv_buffer = malloc(recv_buffer_size);
   if (recv_buffer == NULL) {
      fprintf(stderr, "malloc failed for receive buffer of index set exchange\n");
      exit(EXIT_FAILURE);
   }

   // exchange buffers
   mpi_err_handler(MPI_Allgatherv(send_buffer, send_buffer_size, MPI_PACKED,
                                  recv_buffer, recvcount, displ, MPI_PACKED,
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

      if (xmap->num_src_idxsets_per_process[i]) {
         xmap->src_idxsets_all[i] = malloc(xmap->num_src_idxsets_per_process[i] *
                                           sizeof(*(xmap->src_idxsets_all[i])));
         if (xmap->src_idxsets_all[i] == NULL) {
            fprintf(stderr, "malloc failed for source index sets [i] all of index set exchange\n");
            exit(EXIT_FAILURE);
         }
      } else
         xmap->src_idxsets_all[i] = NULL;

      if (xmap->num_dst_idxsets_per_process[i]) {
         xmap->dst_idxsets_all[i] = malloc(xmap->num_dst_idxsets_per_process[i] *
                                           sizeof(*(xmap->dst_idxsets_all[i])));
         if (xmap->dst_idxsets_all[i] == NULL) {
            fprintf(stderr, "malloc failed for destination index sets [i] all of index set exchange\n");
            exit(EXIT_FAILURE);
         }
      } else
         xmap->dst_idxsets_all[i] = NULL;

      for (j = 0; j < xmap->num_src_idxsets_per_process[i]; ++j)
         xmap->src_idxsets_all[i][j] = xt_idxset_unpack(recv_buffer, recv_buffer_size,
                                                        &position, xmap->comm);
      for (j = 0; j < xmap->num_dst_idxsets_per_process[i]; ++j)
         xmap->dst_idxsets_all[i][j] = xt_idxset_unpack(recv_buffer, recv_buffer_size,
                                                        &position, xmap->comm);
   }

   // clean up
   free(send_buffer);
   free(recv_buffer_header);
   free(recv_buffer);
   free(recvcount);
   free(displ);
}

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
   mpi_err_handler(MPI_Comm_dup(comm, &(xmap->comm)), comm);
   mpi_err_handler(MPI_Comm_size(xmap->comm, &(xmap->comm_size)), xmap->comm);

   // exchange index sets between all processes in comm
   exchange_idxsets(xmap);

   return xmap;
}

void xt_xmap_delete(xt_xmap xmap) {

   unsigned i;

   free(xmap->src_idxsets_local);
   free(xmap->dst_idxsets_local);

   for (i = 0; i < xmap->comm_size; ++i) {

      free(xmap->src_idxsets_all[i]);
      free(xmap->dst_idxsets_all[i]);
   }

   free(xmap->src_idxsets_all);
   free(xmap->dst_idxsets_all);
   free(xmap->num_src_idxsets_per_process);
   free(xmap->num_dst_idxsets_per_process);
   mpi_err_handler(MPI_Comm_free(&(xmap->comm)), MPI_COMM_WORLD);
   free(xmap);
}
