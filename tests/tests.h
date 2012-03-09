#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>

unsigned err_count__ = 0;

#ifdef VERBOSE
#define PUT_ERR(string) err_count__++, fputs((string), stderr);
#else
#define PUT_ERR(string) err_count__++;
#endif

#define INC_ERR (err_count__++)

#define TEST_EXIT_CODE ((!err_count__)?EXIT_SUCCESS:EXIT_FAILURE)

