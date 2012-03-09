#include "tests.h"

#include "yaxt.h"

int main(void) {

   int index_list[] = {1,2,3,4,5,6,7};
   unsigned num_indices = sizeof(index_list) / sizeof(index_list[0]);

   xt_idxset idxset;

   idxset = xt_idxset_new(index_list, num_indices);

   xt_idxset_delete(idxset);

   return TEST_EXIT_CODE;
}
