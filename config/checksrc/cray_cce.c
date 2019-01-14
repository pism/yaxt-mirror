/**
 * @file cray_cce.c
 * @brief demonstrates a problem in the Cray Compilation Environment
 *
 * @copyright Copyright  (C)  2019 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Thomas Jahns <jahns@dkrz.de>
 * URL: https://doc.redmine.dkrz.de/yaxt/html/
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are  permitted provided that the following conditions are
 * met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the DKRZ GmbH nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#include <stdlib.h>

/**
 * usage:
cc cray_cce.c && a.out

 * Program is expected to complete successfully. The Cray Compilation
 * Environment with versions before 8.7 and 9.0 may produce a segmentation
 * vault.
*/

#define COUNT (10)
#define BUG

void init(int * p_int) {

#ifdef BUG
  if (p_int == NULL) {
    p_int = malloc(COUNT * sizeof(p_int[0]));
    exit(EXIT_FAILURE);
  }
#endif

  for(int i = 0; i < COUNT; i++) p_int[i] = i;
}

int main(void) {

  int * p_int = malloc(2 * COUNT * sizeof(int));
  init(p_int + COUNT);
  free(p_int);

  return EXIT_SUCCESS;
}
