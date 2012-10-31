#ifndef ENSURE_ARRAY_SIZE_H
#define ENSURE_ARRAY_SIZE_H

void
xt_realloc_array(void **array, size_t elem_size, size_t *curr_array_size,
                 size_t requested_size);

#define ENSURE_ARRAY_SIZE(arrayp, curr_array_size, req_size)            \
  do {                                                                  \
    if ((req_size) > (curr_array_size))                                 \
    {                                                                   \
      size_t casize = (curr_array_size);                                \
                                                                        \
      xt_realloc_array((void **)&(arrayp), sizeof(*(arrayp)), &casize,  \
                    (req_size));                                        \
      (curr_array_size) = casize;                                       \
    }                                                                   \
  }                                                                     \
  while(0)

#endif
