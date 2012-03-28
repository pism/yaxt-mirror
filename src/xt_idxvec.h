#ifndef XT_IDXVEC_H
#define XT_IDXVEC_H

#include <stdlib.h>
#include "yaxt.h"

struct xt_idxvec_vtable {

   void  (*delete)(xt_idxvec);

   size_t (*get_pack_size)(xt_idxvec, MPI_Comm);
   void (*pack)(xt_idxvec, void*, int, int*, MPI_Comm);
   xt_idxvec (*get_intersection)(xt_idxvec,xt_idxvec);
   xt_idxvec (*copy)(xt_idxvec);
   int (*get_num_indices)(xt_idxvec);
   void (*get_indices)(xt_idxvec, int*);
};

enum xt_idxvec_type {
   LIST       = 1,
   COLLECTION = 2,
};

struct xt_idxvec {
   struct xt_idxvec_vtable *vtable;
};

struct xt_idxvec_unpack_routine_table {
   enum xt_idxvec_type type;
   xt_idxvec (*xt_idxvec_unpack)(void*,int,int*,MPI_Comm);
};

/**
 * \file xt_idxvec.h
 * \brief index vector declaration
 *
 * contains declaration for the index vector data type
 */
 
/**
 * computes the buffer size in byte that is required to pack the given
 * index list using MPI_Pack
 *
 * @param[in] idxvec index vector for which the packed size needs to be computed
 * @param[in] comm   MPI communicator that will be used to send the buffer for
 *                   which the size is computed
 * @return           returns the size in bytes of the packed index vector
 *
 * @remark You need to provide the same buffer for all related calls of
 *         xt_idxvec_get_pack_size , \ref xt_idxvec_pack and \ref xt_idxvec_unpack
 *
 * @see xt_idxvec_pack
 * @see xt_idxvec_unpack
 */
size_t xt_idxvec_get_pack_size(xt_idxvec idxvec, MPI_Comm comm);

/**
 * packs a index vector into previously allocated buffer using MPI_Pack
 *
 * @param[in]     idxvec      index vector, which needs be packed into the buffer
 * @param[in,out] buffer      previously allocated buffer into which the index
 *                            list will be packed
 * @param[in]     buffer_size size of buffer in bytes
 * @param[in,out] position    position in buffer at which the data of idxvec will be copied
 *                            (position is automatically updated to the first byte after
 *                            the packed data)
 * @param[in]     comm        MPI communicator that will be used to send the buffer
 *
 * @remark The provided buffer needs to have a size that is sufficient for the index vector. You can determin the required size using \ref xt_idxvec_get_pack_size .
 * @remark You need to provide the same buffer for all related calls of
 *         \ref xt_idxvec_get_pack_size , xt_idxvec_pack and \ref xt_idxvec_unpack
 *
 * @see xt_idxvec_get_pack_size
 * @see xt_idxvec_unpack
 *
 */
void xt_idxvec_pack(xt_idxvec idxvec, void* buffer, int buffer_size,
                    int* position, MPI_Comm comm);

/**
 * unpacks a index vector that received before
 *
 * @param[in]     buffer      buffer that contains the packed index vector
 * @param[in]     buffer_size size of the buffer in bytes
 * @param[in,out] position    position in the buffer at which the unpacking should start
 *                            (will automatically be set to the position after the unpacked
 *                            data)
 * @param[in]     comm        MPI communicator that was used to receive the buffer
 *
 * @remark You need to provide the same buffer for all related calls of
 *         \ref xt_idxvec_get_pack_size , \ref xt_idxvec_pack and xt_idxvec_unpack
 *
 * @see xt_idxvec_pack
 */
xt_idxvec xt_idxvec_unpack(void* buffer, int buffer_size, int* position, MPI_Comm comm);

/**
 * computes the intersection between two index vectors
 *
 * @param[in] idxvec_a index vector a
 * @param[in] idxvec_b index vector b
 * @return             return the intersection of idxvec_a and idxvec_b
 */
xt_idxvec xt_idxvec_get_intersection(xt_idxvec idxvec_a, xt_idxvec idxvec_b);

/**
 * generates a copy of a given index vector
 *
 * @param[in] idxvec index vector that is to be copied
 * @return           returns a copy of the given index vector
 */
xt_idxvec xt_idxvec_copy(xt_idxvec idxvec);

/**
 * returns the number of indices stored in the given index vector
 *
 * @param[in] idxvec index vector for which the number of indices is required
 * @return           number of indices in the given index vector
 */
int xt_idxvec_get_num_indices(xt_idxvec idxvec);
 
 /**
  * gets the indices stored in an index vector
  *
  * @param[in]     idxvec  index vector for which the indices are to be returned
  * @param[in,out] indices array into which the indices are written
  */
void xt_idxvec_get_indices(xt_idxvec idxvec, int *indices);

#endif // XT_IDXVEC_H
