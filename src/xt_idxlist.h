#ifndef XT_IDXLIST_H
#define XT_IDXLIST_H

#include <stdlib.h>
#include "yaxt.h"

struct Xt_stripe {

  Xt_idx start;
  unsigned extent;
  int stride;
};

struct xt_idxlist_vtable {

   void  (*delete)(xt_idxlist);

   size_t (*get_pack_size)(xt_idxlist, MPI_Comm);
   void (*pack)(xt_idxlist, void*, int, int*, MPI_Comm);
   xt_idxlist (*get_intersection)(xt_idxlist,xt_idxlist);
   xt_idxlist (*copy)(xt_idxlist);
   size_t (*get_num_indices)(xt_idxlist);
   void (*get_indices)(xt_idxlist, Xt_idx *indices);
   void (*get_index_stripes)(xt_idxlist, struct Xt_stripe **, unsigned *);
   int (*get_index)(Xt_idxlist, int, Xt_idx *);
   int (*get_position_of_index)(Xt_idxlist, Xt_idx, int *);
   int (*get_position_of_indices)(Xt_idxlist, Xt_idx *, unsigned, int *, int);
   int (*get_position_of_index_off)(Xt_idxlist, Xt_idx, int *, int);
   int (*get_position_of_indices_off)(Xt_idxlist, Xt_idx *, unsigned,
                                      int *, int *);
};

enum xt_idxvecype {
   VECTOR     = 1,
   COLLECTION = 2,
   SECTION    = 3,
   STRIPES    = 4,
};

struct xt_idxlist {
   struct xt_idxlist_vtable *vtable;
};

struct xt_idxlist_unpack_routine_table {
   enum xt_idxvecype type;
   xt_idxlist (*xt_idxlist_unpack)(void*,int,int*,MPI_Comm);
};

/**
 * \file xt_idxlist.h
 * \brief index list declaration
 *
 * contains declaration for the index list datatype
 */

/**
 * computes the buffer size in byte that is required to pack the given
 * index list using MPI_Pack
 *
 * @param[in] idxlist index list for which the packed size needs to be computed
 * @param[in] comm    MPI communicator that will be used to send the buffer for
 *                    which the size is computed
 * @return            returns the size in bytes of the packed index list
 *
 * @remark You need to provide the same buffer for all related calls of
 *         xt_idxlist_get_pack_size , \ref xt_idxlist_pack and \ref xt_idxlist_unpack
 *
 * @see xt_idxlist_pack
 * @see xt_idxlist_unpack
 */
size_t xt_idxlist_get_pack_size(xt_idxlist idxlist, MPI_Comm comm);

/**
 * packs an index list into previously allocated buffer using MPI_Pack
 *
 * @param[in]     idxlist     index list, which needs be packed into the buffer
 * @param[in,out] buffer      previously allocated buffer into which the index
 *                            list will be packed
 * @param[in]     buffer_size size of buffer in bytes
 * @param[in,out] position    position in buffer at which the data of idxlist will be copied
 *                            (position is automatically updated to the first byte after
 *                            the packed data)
 * @param[in]     comm        MPI communicator that will be used to send the buffer
 *
 * @remark The provided buffer needs to have a size that is sufficient for the
 *         index list. You can determin the required size using
 *         \ref xt_idxlist_get_pack_size .
 * @remark You need to provide the same buffer for all related calls of
 *         \ref xt_idxlist_get_pack_size , xt_idxlist_pack and \ref xt_idxlist_unpack
 *
 * @see xt_idxlist_get_pack_size
 * @see xt_idxlist_unpack
 *
 */
void xt_idxlist_pack(xt_idxlist idxlist, void* buffer, int buffer_size,
                     int* position, MPI_Comm comm);

/**
 * unpacks a index list that was received before
 *
 * @param[in]     buffer      buffer that contains the packed index list
 * @param[in]     buffer_size size of the buffer in bytes
 * @param[in,out] position    position in the buffer at which the unpacking should start
 *                            (will automatically be set to the position after the unpacked
 *                            data)
 * @param[in]     comm        MPI communicator that was used to receive the buffer
 *
 * @remark You need to provide the same buffer for all related calls of
 *         \ref xt_idxlist_get_pack_size , \ref xt_idxlist_pack and xt_idxlist_unpack
 *
 * @see xt_idxlist_pack
 */
xt_idxlist xt_idxlist_unpack(void* buffer, int buffer_size,
                             int* position, MPI_Comm comm);

/**
 * computes the intersection between two index list
 *
 * @param[in] idxlist_src index list for sender
 * @param[in] idxlist_dst index list for receiver
 * @return                return the intersection of idxlist_a and idxlist_b
 * @remark duplicate elements i in the receiver list will result in
 *         duplicate elements in the intersection if idxlist_src contains i
 * @remark the elements in the resulting index list are sorted
 */
xt_idxlist
xt_idxlist_get_intersection(xt_idxlist idxlist_src, xt_idxlist idxlist_dst);

/**
 * generates a copy of a given index list
 *
 * @param[in] idxlist index list that is to be copied
 * @return            returns a copy of the given index list
 */
xt_idxlist xt_idxlist_copy(xt_idxlist idxlist);

/**
 * returns the number of indices stored in the given index list
 *
 * @param[in] idxlist index list for which the number of indices is required
 * @return            number of indices in the given index list
 */
int xt_idxlist_get_num_indices(xt_idxlist idxlist);

/**
 * gets the indices stored in an index list
 *
 * @param[in]     idxlist index list for which the indices are to be returned
 * @param[in,out] indices array into which the indices are written
 */
void xt_idxlist_get_indices(xt_idxlist idxlist, Xt_idx *indices);

/**
 * returns the indices stored in an index list
 * the indices are returned in form of stripes (with stride 1)
 *
 * @param[in]  idxlist     index list for which the indices are to be returned
 * @param[out] stripes     array containing the stripes (user is responsible
 *                         for freeing the memory associated with the stripes)
 * @param[out] num_stripes number of stripes in array "stripes"
 */
void xt_idxlist_get_index_stripes(xt_idxlist idxlist,
                                  struct Xt_stripe ** stripes,
                                  unsigned * num_stripes);

/**
 * gets the index stored at the specified position in an index list
 *
 * @param[in]  idxlist  index list
 * @param[in]  position position of the index that is to be retrieved
 * @param[out] index    index stored in the index list at the given position
 * @return              0 if everything was okay
 *                      1 if the position is out of range
 */

int xt_idxlist_get_index(Xt_idxlist idxlist, int position, Xt_idx * index);

/**
 * gets the position of the first occurrence of the given index
 * in the stored index list
 *
 * @param[in]  idxlist  index list in which the given index is to be searched
 * @param[in]  index    index that is to be searched
 * @param[out] position position of the first occurrence of index in idxlist
 * @return              0 if everything was okay
 *                      1 if the index could not be found in idxlist
 */
int xt_idxlist_get_position_of_index(Xt_idxlist idxlist, Xt_idx index, int * position);

/**
 * gets the positions of the first occurrence of the given indices
 * in the stored index list
 *
 * @param[in]  idxlist     index list in which the given index is to be searched
 * @param[in]  indices     indices that are to be searched
 * @param[in]  num_indices number of indices in that array "indices"
 * @param[out] positions   positions of the first occurrence of the indices
 *                         in idxlist
 * @param[in]  single_match_only if true then do not consider a previous matching position again
 *                         in the next index search - this is required on the target side in order
 *                         to avoid multiple writes to the same memory position
 * @return                 0 if everything was okay
 *                         1 if the index could not be found in idxlist
 */
int xt_idxlist_get_position_of_indices(Xt_idxlist idxlist, Xt_idx * indices,
                                       unsigned num_indices, int * positions,
                                       int single_match_only);

/**
 * gets the position of the first occurrence of the given index
 * from the given offset in the stored index list idxlist
 *
 * @param[in]  idxlist  index list in which the given index is to be searched
 * @param[in]  index    index that is to be searched
 * @param[out] position position of the first occurrence of index in idxlist
 *                      with a position >= offset
 * @param[in]  offset   offset in the index list from which on the search
 *                      is to be conducted
 * @return              0 if everything was okay
 *                      1 if the index could not be found in idxlist
 */
int xt_idxlist_get_position_of_index_off(Xt_idxlist idxlist, Xt_idx index,
                                         int * position, int offset);

/**
 * gets the positions of the first occurrence of the given indices
 * from the given offsets in the stored index list
 *
 * @param[in]  idxlist     index list in which the given index is to be searched
 * @param[in]  indices     indices that are to be searched
 * @param[in]  num_indices number of indices in that array "indices"
 * @param[out] positions   positions of the first occurrence of the indices in idxlist
 *                         with a positions[i] >= offsets[i]
 * @param[in]  offsets     offsets in the index list from which on the search
 *                         is to be conducted
 * @return                 0 if everything was okay
 *                         1 if the index could not be found in idxlist
 *
 * \remarks a call to xt_idxlist_get_position_of_index_off_a with offset = NULL
 *          is the same as a call to xt_idxlist_get_position_of_indices
 */
int xt_idxlist_get_position_of_indices_off(Xt_idxlist idxlist, Xt_idx * indices,
                                           unsigned num_indices, int * positions,
                                           int * offsets);

#endif // XT_IDXLIST_H
