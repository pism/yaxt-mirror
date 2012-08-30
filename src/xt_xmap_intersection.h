#ifndef XT_XMAP_INTERSECTION_H
#define XT_XMAP_INTERSECTION_H

#include "yaxt.h"

/**
 * constructor for an exchange map \n
 * this operation is collective over all processes in comm \n
 * it uses the provided intersection information to generate the exchange map
 *
 * @param[in] src_intersections     array containing the source intersections with
 *                                  the other processes
 * @param[in] src_ranks             array containing for each source intersection
 *                                  the ranks of the other processes
 * @param[in] num_src_intersections number of source intersections
 * @param[in] dst_intersections     array containing the destination intersections with
 *                                  the other processes
 * @param[in] dst_ranks             array containing for each destination intersection
 *                                  the ranks of the other processes
 * @param[in] num_dst_intersections number of destination intersections
 * @param[in] src_idxlist           source index list
 * @param[in] dst_idxlist           destination index list
 * @param[in] comm                  MPI communicator that contains all processes
 *                                  that part in the exchange
 */
Xt_xmap xt_xmap_intersection_new(Xt_idxlist * src_intersections, int * src_ranks,
                                 Xt_count num_src_intersections,
                                 Xt_idxlist * dst_intersections, int * dst_ranks,
                                 Xt_count num_dst_intersections,
                                 Xt_idxlist src_idxlist, Xt_idxlist dst_idxlist,
                                 MPI_Comm comm);

#endif // XT_XMAP_INTERSECTION_H
