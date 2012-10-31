#include <mpi.h>
#include <yaxt.h>

#include "tests.h"

int main(void) {
  const int hsize = 32;
  int data[hsize];
  int handle[hsize];
  int state[hsize];

  xt_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);
  xt_initialize(MPI_COMM_WORLD);

  for (int i = 0; i<hsize; i++) {
    state[i] = 0;
  }

  struct xt_handle_set_type *hset = xt_handle_set_new(4);

  // new handles
  for (int i = 0; i<hsize; i++) {
    if (i % 2 == 0) {
      data[i]=i;
      handle[i] = xt_handle_new(hset, &data[i]);
      if (!xt_handle_is_valid(hset,handle[i])) PUT_ERR("unexpected invalid handle\n");
      state[i] = 1;
    }
  }

  // delete some handles
  for (int i = 1; i<hsize; i++) {
    if (state[i]) {
      if (i % 3 == 0) {
        if (!xt_handle_is_valid(hset,handle[i])) PUT_ERR("unexpected invalid handle\n");
        xt_handle_delete(hset,handle[i]);
        if (xt_handle_is_valid(hset,handle[i])) PUT_ERR("unexpected valid handle\n");
        state[i] = 0;
      }
    }
  }

  // more new handles
  for (int i = 0; i<hsize; i++) {
    if (!state[i]) {
      data[i]=i;
      handle[i] = xt_handle_new(hset, &data[i]);
      if (!xt_handle_is_valid(hset,handle[i])) PUT_ERR("unexpected invalid handle\n");
      state[i] = 1;
    }
  }

  // delete some handles
  for (int i = 0; i<hsize; i++) {
    if (state[i]) {
      if (i % 5 == 0) {
        int h = handle[i];
        if (!xt_handle_is_valid(hset,h)) PUT_ERR("unexpected invalid handle\n");
        xt_handle_delete(hset,h);
        if (xt_handle_is_valid(hset,h)) PUT_ERR("unexpected valid handle\n");
        state[i] = 0;
      }
    }
  }

  // check content
  for (int i = 0; i<hsize; i++) {
    if (state[i]) {
      int h = handle[i];
      if (!xt_handle_is_valid(hset,h)) PUT_ERR("unexpected invalid handle\n");
      int *p = xt_handle2pointer(hset, h);
      if (*p != data[i]) PUT_ERR("data lookup failed\n");
    }
  }

  xt_handle_set_delete(hset);

  MPI_Finalize();
  return TEST_EXIT_CODE;
}
