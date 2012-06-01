#include "tests.h"

#include "xt_mpi.h"
#include "yaxt.h"
#include "xt_handles.h"

int main(void) {
  const int hsize = 32;
  int data[hsize];
  int handle[hsize];

  mpi_err_handler(MPI_Init(NULL, NULL), MPI_COMM_WORLD);
  Xt_initialize(MPI_COMM_WORLD);

  struct xt_handle_set_type *hset = xt_handle_set_new(4);
  for (int i = 0; i<hsize; i++) {
    data[i]=i;
    handle[i] = xt_handle_new(&hset, &data[i]);
  }

  for (int i = 0; i<hsize; i++) {
    if (!xt_handle_is_valid(hset,i)) PUT_ERR("unexpected invalid handle\n");
    int *p = xt_handle2pointer(hset, i);
    if (*p != data[i]) PUT_ERR("data lookup failed\n");
  }

  xt_handle_set_delete(&hset);

  MPI_Finalize();
  return TEST_EXIT_CODE;

}
