#include <stdlib.h>

#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"
#include "xt/xt_handles.h"


struct xt_handle_set_type {
  int num; // number of used handles
  int cap; // capacity
  int hnf; // hint for next free entry
  void **p; // vector of pointers to user data, p[handle] == pointer to user data
};

static const int default_handle_set_cap = 16;

static void extend_handle_set(struct xt_handle_set_type *hset);

int xt_handle_is_valid(struct xt_handle_set_type *hset, int ih) {
  if(!hset) return 0;
  if (ih < 0 || ih >= hset->cap) return 0;
  if (hset->p[ih]) return 1;
  return 0;
}

void *xt_handle2pointer(struct xt_handle_set_type *hset, int ih) {
  if (!xt_handle_is_valid(hset, ih)) return NULL;
  return hset->p[ih];
}

struct xt_handle_set_type *xt_handle_set_new(int cap) {
  struct xt_handle_set_type *hset;
  hset = malloc( sizeof(struct xt_handle_set_type) );
  hset->p = malloc(cap * sizeof(void*));
  hset->num = 0;
  hset->cap = cap;
  if (hset->cap < 1) hset->cap = default_handle_set_cap;
  hset->hnf = 0;
  for (int i=0; i<cap; i++) {
    hset->p[i] = NULL;
  }

  return hset;
}


static void extend_handle_set(struct xt_handle_set_type *hset) {
  if (!hset) {
    hset = xt_handle_set_new(default_handle_set_cap);
    return;
  }

  int old_cap, new_cap;
  old_cap = hset->cap;
  if (old_cap < 2) {
    new_cap = default_handle_set_cap;
  } else {
    new_cap = old_cap*2;
  }

  hset->p = xrealloc(hset->p, new_cap * sizeof(void*));

  for (int i=old_cap; i<new_cap; i++) {
    hset->p[i] = NULL;
  }
  hset->hnf = old_cap;
  hset->cap = new_cap;
}


void xt_handle_set_delete(struct xt_handle_set_type *hset) {
  free(hset->p);
  free(hset);
  hset=NULL;
}


int xt_handle_new(struct xt_handle_set_type *hset, void *p) {
  if (p == NULL) return -1;

  if (hset->num >= hset->cap) extend_handle_set(hset);
  for (int j=0; j<hset->cap; j++) {
    int i = (hset->hnf + j) % hset->cap;
    if (!hset->p[i]) {
      hset->p[i] = p;
      hset->hnf = (i+1) % hset->cap;
      hset->num++;
      return i;
    }
  }

  die("internal error");

  return 0;
}

void xt_handle_delete(struct xt_handle_set_type *hset, int handle) {
  if (!xt_handle_is_valid(hset, handle)) die("delete_handle: invalid handle");
  hset->p[handle]=NULL;
  hset->num--;
  if (hset->num < 0) die("delete_handle: internal error - negative occupation count");
  if (hset->p[hset->hnf]) hset->hnf = handle;
}

