/**
 * @copyright Copyright(C) 2023
 * @file bre2pcre.c
 * Original program by https://stackoverflow.com/users/9952196/shawn Shawn
 * provided under CC BY-SA 4.0
 * https://creativecommons.org/licenses/by-sa/4.0/
 * Adapted from:
 * https://stackoverflow.com/questions/76587209/how-to-convert-bre-to-perl-re/76589961#76589961
 *
 * Compile with:
 *    gcc -o bre2perl -O2 -march=native -Wall -Wextra bre2perl.c \
 *           $(pcre2-config --cflags --libs8)
 */

#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PCRE2_CODE_UNIT_WIDTH 8
#include <pcre2.h>

int main(void) {
  char *line = NULL;
  size_t linelength = 0;
  ssize_t len;

  while ((len = getline(&line, &linelength, stdin)) > 1) {
    if (line[len - 1] == '\n') { len -= 1; }
    PCRE2_UCHAR *perl_re = NULL;
    PCRE2_SIZE re_len = 0;
    int rc = pcre2_pattern_convert((PCRE2_SPTR)line, len,
                                   PCRE2_CONVERT_POSIX_BASIC,
                                   &perl_re, &re_len, NULL);
    if (rc == 0) {
      // There's a PCRE2-specific prefix to patterns that should be removed
      printf("%s\n", (const char *)perl_re
              + (strncmp((const char *)perl_re, "(*NUL)", 6) == 0 ? 6 : 0));
      pcre2_converted_pattern_free(perl_re);
    } else {
      PCRE2_UCHAR errmsg[1024];
      pcre2_get_error_message(rc, errmsg, sizeof errmsg);
      fprintf(stderr, "error converting pattern '%s': %s\n",
              line, (char *)errmsg);
      return EXIT_FAILURE;
    }
  }
  free(line);
  return 0;
}

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
