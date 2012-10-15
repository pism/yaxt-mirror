#ifndef XT_IDXMOD_H
#define XT_IDXMOD_H

#include "yaxt.h"
#include "xt_idxlist.h"

/** \example test_idxmod.c
 */

/**
 * A modifier is a pair of compact idxlists that describes a index-to-index stream.
 * Additionally there is a mask value that can be used to track a sequence of modifications (see \ref xt_idxmod_new).
 *
 */

struct Xt_modifier {
  struct Xt_idxlist* extract; ///<  idx values
  struct Xt_idxlist* subst; ///<  idx values, must describe at least as many indices as the extract idxlist
  int mask; ///< modifier mask, can be used to identify modifier usage in the result
};

/**
 * \brief generates a new index list based on an index list and a sequence of modifiers
 *
 * In the new index list all occurences of indices listed in the \ref Xt_modifier "modifier extract list"
 * are replaced by respective indices of the modifier substitute list.
 * Multiple modifiers are applied consecutively to the input index list. If a modification state array is given then it will be updated using the
 * mask values in the used modifiers (OR operation). An update only affects indices that are modified by the respective modifier.
 *
 * @param[in]     patch_idxlist  input index list
 * @param[in]     modifier       array of modifiers
 * @param[in]     modifier_num   number of modifiers
 * @param[in,out] mstate         modification state
 * @return returns a Xt_idxlist
 * @remark For performance considerations imagine patch_idxlist as a (small) task-local index list whose transform into an vector scales ideally
 * and the modifier idxlist-components as compact idxlists of global extent which should never downgrade to a elementwise index description.
 * @remark mstate is not initilaized
 * @remark mstate == NULL is valid (disables the modification state handling)
 * @remark if mstate != NULL then it is required to be at least as big as the input index list (patch_idxlist)
 */

struct Xt_idxlist *xt_idxmod_new(struct Xt_idxlist* patch_idxlist,  struct Xt_modifier *modifier, int modifier_num, int *mstate);

#endif
