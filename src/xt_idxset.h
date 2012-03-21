#ifndef XT_IDXSET_H
#define XT_IDXSET_H

#include "yaxt.h"

struct xt_idxset_vtable {

   void  (*delete)(xt_idxset);

   size_t (*get_pack_size)(xt_idxset, MPI_Comm);
   void (*pack)(xt_idxset, void*, int, int*, MPI_Comm);
   xt_idxset (*get_intersection)(xt_idxset,xt_idxset);
   xt_idxset (*copy)(xt_idxset);
   int (*get_num_indices)(xt_idxset);
   void (*get_indices)(xt_idxset, int*);
};

enum xt_idxset_type {
   LIST  = 1,
   MULTI = 2,
};

struct xt_idxset {
   struct xt_idxset_vtable *vtable;
};

/**
 * \file xt_idxset.h
 * \brief index set declaration
 *
 * contains declaration for the index set data type
 */
 
/**
 * computes the buffer size in byte that is required to pack the given
 * index list using MPI_Pack
 *
 * @param[in] idxset index set for which the packed size needs to be computed
 * @param[in] comm   MPI communicator that will be used to send the buffer for
 *                   which the size is computed
 * @return           returns the size in bytes of the packed index set
 *
 * @remark You need to provide the same buffer for all related calls of
 *         xt_idxset_get_pack_size , \ref xt_idxset_pack and \ref xt_idxset_unpack
 *
 * @see xt_idxset_pack
 * @see xt_idxset_unpack
 */
size_t xt_idxset_get_pack_size(xt_idxset idxset, MPI_Comm comm);

/**
 * packs a index set into previously allocated buffer using MPI_Pack
 *
 * @param[in]     idxset      index set, which needs be packed into the buffer
 * @param[in,out] buffer      previously allocated buffer into which the index
 *                            list will be packed
 * @param[in]     buffer_size size of buffer in bytes
 * @param[in,out] position    position in buffer at which the data of idxset will be copied
 *                            (position is automatically updated to the first byte after
 *                            the packed data)
 * @param[in]     comm        MPI communicator that will be used to send the buffer
 *
 * @remark The provided buffer needs to have a size that is sufficient for the index set. You can determin the required size using \ref xt_idxset_get_pack_size .
 * @remark You need to provide the same buffer for all related calls of
 *         \ref xt_idxset_get_pack_size , xt_idxset_pack and \ref xt_idxset_unpack
 *
 * @see xt_idxset_get_pack_size
 * @see xt_idxset_unpack
 *
 */
void xt_idxset_pack(xt_idxset idxset, void* buffer, int buffer_size,
                    int* position, MPI_Comm comm);

/**
 * unpacks a index set that received before
 *
 * @param[in]     buffer      buffer that contains the packed index set
 * @param[in]     buffer_size size of the buffer in bytes
 * @param[in,out] position    position in the buffer at which the unpacking should start
 *                            (will automatically be set to the position after the unpacked
 *                            data)
 * @param[in]     comm        MPI communicator that was used to receive the buffer
 *
 * @remark You need to provide the same buffer for all related calls of
 *         \ref xt_idxset_get_pack_size , \ref xt_idxset_pack and xt_idxset_unpack
 *
 * @see xt_idxset_pack
 */
xt_idxset xt_idxset_unpack(void* buffer, int buffer_size, int* position, MPI_Comm comm);

/**
 * computes the intersection between two index sets
 *
 * @param[in] idxset_a index set a
 * @param[in] idxset_b index set b
 * @return             return the intersection of idxset_a and idxset_b
 */
xt_idxset xt_idxset_get_intersection(xt_idxset idxset_a, xt_idxset idxset_b);

/**
 * generates a copy of a given index set
 *
 * @param[in] idxset index set that is to be copied
 * @return           returns a copy of the given index set
 */
xt_idxset xt_idxset_copy(xt_idxset idxset);

/**
 * returns the number of indices stored in the given index set
 *
 * @param[in] idxset index set for which that number of index sets is required
 * @return           number of indices in the given index set
 */
int xt_idxset_get_num_indices(xt_idxset idxset);
 
 /**
  * gets the indices stored in an index sets
  *
  * @param[in]     idxset  index set for which the indices are to be returned
  * @param[in,out] indices array into which the indices are written
  */
void xt_idxset_get_indices(xt_idxset idxset, int *indices);

#endif // XT_IDXSET_H
