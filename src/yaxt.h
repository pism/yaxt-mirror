/*
  TODO: doxygenify this
*/

/*
  copy/delete a local selection of global indices
*/

xt_idxset *xt_idxset_new(int *idxset);
void xt_idxset_delete(int *idxset);


/*
  create/delete a communication map between index sets
*/

xt_idxmap *xt_idxmap_new(xt_idxset *src_idxset, xt_idxset *dst_idxset);
void xt_idxmap_delete(xt_idxset *src_idxset, xt_idxset *dst_idxset);

/*
  inquiry function to check if the inverse is defined
*/
int xt_idxmap_is_bijective(xt_idxmap *idxmap);

