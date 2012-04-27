#ifndef XT_XMAP_H
#define XT_XMAP_H

#include "yaxt.h"

/**
 * \file xt_xmap.h
 * \brief exchange map declaration
 *
 * contains declaration for the exchange map data structure
 */

/**
 * gets the destination exchange list for a given rank
 *
 * @param[in] xmap exchange map
 * @param[in] rank rank of process for which the exchange list is to be retrieved
 * @return         destination exchange list for given rank
 *                 NULL in case there is no exchange with the given rank
 */
Xt_idxvec xt_xmap_get_dst_xlist(Xt_xmap xmap, int rank);

/**
 * gets the source exchange list for a given rank
 *
 * @param[in] xmap exchange map
 * @param[in] rank rank of process for which the exchange list is to be retrieved
 * @return         source exchange list for given rank
 *                 NULL in case there is no exchange with the given rank
 */
Xt_idxvec xt_xmap_get_src_xlist(Xt_xmap xmap, int rank);

#endif // XT_XMAP_H
