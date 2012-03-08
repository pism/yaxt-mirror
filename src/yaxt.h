#ifndef YAXT_H
#define YAXT_H

/*
  \todo doxify this
*/

/*
  copy/delete a local selection of global indices
*/

xt_idxset *xt_idxset_new(int *idxlist);
void xt_idxset_delete(xt_idxset *idxset);


/*
  create/delete a communication map between index sets
*/

xt_xmap *xt_xmap_new(xt_idxset *src_idxset, xt_idxset *dst_idxset);
void xt_xmap_delete(xt_xmap *xmap);

/*
  inquiry function to check if the inverse is defined
*/
int xt_xmap_is_bijective(xt_xmap *xmap);

#endif // YAXT_H