/*
!**********************************************************************!
!Copyright (c) 2017, Harald Braun, science + computing ag (s+c) / Atos S.A.
!                                  <harald.braun@atos.net>
!All rights reserved.
!
!Redistribution and use in source and binary forms, with or without
!modification, are permitted provided that the following conditions are met:
!    * Redistributions of source code must retain the above copyright
!      notice, this list of conditions and the following disclaimer.
!    * Redistributions in binary form must reproduce the above copyright
!      notice, this list of conditions and the following disclaimer in the
!      documentation and/or other materials provided with the distribution.
!    * Neither the name of the <organization> nor the
!      names of its contributors may be used to endorse or promote products
!      derived from this software without specific prior written permission.
!
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
!ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
!WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!DISCLAIMED. IN NO EVENT SHALL Harald Braun, or s+c BE LIABLE FOR ANY
!DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
!(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
!ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!**********************************************************************!

!----------------------------------------------------------------------

!Purpose and usage:
!
!
! include this header in your source and call
!
!    analyze_MPI_datatype(const MPI_Datatype DDT, const char * RECONST, const char * DDTFn_NAME)
!
! provide a 
!  1) an MPI_Datatype handle to analyze >DDT<
!  2) a char-string containing a filename prefix ('.c' or '.h' will be appended) for the 
!     reconstruction code >RECONST<.
!     If empty, no reconstructor code is generated.
!     Otherwise, two files are generated: RECONST.c and RECONST.h
!     To use the reconstructed datatyp compile <reconst>.c using an MPI-enabled compiler and
!     put the following ito your code:
!     "#include <RECONST.h>
!  3) a char-string containing a specific name for the MPI_Datatype reconstruction function >DDTFn_NAME<.
!     If empty, all reconstruction functions will have the very same name!
! 
!  To use the reconstructed MPI_Datatype compile <reconst>.c using an MPI-enabled compiler and
!  put the following ito your code:
!
!  #include <RECONST.h>
!  ...
!  MPI_Datatype my_ddt = DDTFn_NAME();
!  ...
!  ...commit and use the datatype ...
!  ...
!  ... do not forget to free it, once it is not used anymore ;)
!  ...
*/

#ifndef HAVE_HB_ANALYZE_MPI_DDT
#define HAVE_HB_ANALYZE_MPI_DDT
#include <sys/file.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>

static char const COMBINER_NAMED_STRING[]            = "MPI_COMBINER_NAMED: predefined datatype";
static char const COMBINER_DUP_STRING[]              = "MPI_COMBINER_DUP: MPI_Type_dup()";
static char const COMBINER_CONTIGUOUS_STRING[]       = "MPI_COMBINER_CONTIGUOUS: MPI_Type_contiguous()";
static char const COMBINER_VECTOR_STRING[]           = "MPI_COMBINER_VECTOR: MPI_Type_vector()";
static char const COMBINER_HVECTOR_INTEGER_STRING[]  = "MPI_COMBINER_HVECTOR_INTEGER: MPI_Type_hvector() from FORTRAN";
static char const COMBINER_HVECTOR_STRING[]          = "MPI_COMBINER_HVECTOR: MPI_Type_hvector() from C/C++ or MPI_Type_create_hvector()";
static char const COMBINER_INDEXED_STRING[]          = "MPI_COMBINER_INDEXED: MPI_Type_indexed()";
static char const COMBINER_HINDEXED_INTEGER_STRING[] = "MPI_COMBINER_HINDEXED_INTEGER: MPI_Type_hindexed() from FORTRAN";
static char const COMBINER_HINDEXED_STRING[]         = "MPI_COMBINER_HINDEXED: MPI_Type_hindexed() from C/C++ or MPI_Type_create_hindexed()";
static char const COMBINER_INDEXED_BLOCK_STRING[]    = "MPI_COMBINER_INDEXED_BLOCK: MPI_Type_create_indexed_block()";
static char const COMBINER_STRUCT_INTEGER_STRING[]   = "MPI_COMBINER_STRUCT_INTEGER: MPI_Type_struct() from FORTRAN";
static char const COMBINER_STRUCT_STRING[]           = "MPI_COMBINER_STRUCT: MPI_Type_struct() from C/C++ or MPI_Type_create_struct()";
static char const COMBINER_SUBARRAY_STRING[]         = "MPI_COMBINER_SUBARRAY: MPI_Type_create_subarray()";
static char const COMBINER_DARRAY_STRING[]           = "MPI_COMBINER_DARRAY: MPI_Type_create_darray()";
static char const COMBINER_F90_REAL_STRING[]         = "MPI_COMBINER_F90_REAL: MPI_Type_create_f90_real()";
static char const COMBINER_F90_COMPLEX_STRING[]      = "MPI_COMBINER_F90_COMPLEX: MPI_Type_create_f90_complex()";
static char const COMBINER_F90_INTEGER_STRING[]      = "MPI_COMBINER_F90_INTEGER: MPI_Type_create_f90_integer()";
static char const COMBINER_RESIZED_STRING[]          = "MPI_COMBINER_RESIZED: MPI_Type_create_resized()";
static char const COMBINER_UNKOWN_STRING[]           = "MPI_COMBINER_UNKOWN: This is an error!";

struct mdt_info;

struct mdt_info {
  char const *combiner_string;
  char name_string[MPI_MAX_OBJECT_NAME];
  MPI_Datatype mdt;
  MPI_Count pack_size;
  MPI_Count lb;
  MPI_Count extent;
  MPI_Count true_lb;
  MPI_Count true_extent;
  MPI_Aint *mdt_addrargs;
  MPI_Datatype *mdt_dtargs;
  int *mdt_intargs;
  struct mdt_info *next;
  struct mdt_info *previous;
  struct mdt_info *parent;
  struct mdt_info *first_child;
  struct mdt_info *last_child;
  int nest;
  int env_num_ints;
  int env_num_addr;
  int env_num_mdt;
  int env_combiner_id;
};

static void wmpi(const int mpi_err, const char * func_name, const int icall)
{
  if (mpi_err != MPI_SUCCESS) {
    fprintf(stderr, "got MPI_ERROR in %s call %d", func_name, icall);
    fflush(stderr);
    int aberr = MPI_Abort(MPI_COMM_WORLD, mpi_err);
    if(aberr != MPI_SUCCESS)
      abort();
  }
  return;
};

static char const * get_combiner_name(const int combiner_id)
{
  switch (combiner_id)
  {
    case MPI_COMBINER_NAMED:            return COMBINER_NAMED_STRING;
    case MPI_COMBINER_DUP:              return COMBINER_DUP_STRING;
    case MPI_COMBINER_CONTIGUOUS:       return COMBINER_CONTIGUOUS_STRING;
    case MPI_COMBINER_VECTOR:           return COMBINER_VECTOR_STRING;
    case MPI_COMBINER_HVECTOR_INTEGER:  return COMBINER_HVECTOR_INTEGER_STRING;
    case MPI_COMBINER_HVECTOR:          return COMBINER_HVECTOR_STRING;
    case MPI_COMBINER_INDEXED:          return COMBINER_INDEXED_STRING;
    case MPI_COMBINER_HINDEXED_INTEGER: return COMBINER_HINDEXED_INTEGER_STRING;
    case MPI_COMBINER_HINDEXED:         return COMBINER_HINDEXED_STRING;
    case MPI_COMBINER_INDEXED_BLOCK:    return COMBINER_INDEXED_BLOCK_STRING;
    case MPI_COMBINER_STRUCT_INTEGER:   return COMBINER_STRUCT_INTEGER_STRING;
    case MPI_COMBINER_STRUCT:           return COMBINER_STRUCT_STRING;
    case MPI_COMBINER_SUBARRAY:         return COMBINER_SUBARRAY_STRING;
    case MPI_COMBINER_DARRAY:           return COMBINER_DARRAY_STRING;
    case MPI_COMBINER_F90_REAL:         return COMBINER_F90_REAL_STRING;
    case MPI_COMBINER_F90_COMPLEX:      return COMBINER_F90_COMPLEX_STRING;
    case MPI_COMBINER_F90_INTEGER:      return COMBINER_F90_INTEGER_STRING;
    case MPI_COMBINER_RESIZED:          return COMBINER_RESIZED_STRING;
    default:                            return COMBINER_UNKOWN_STRING;
  };
};

static void print_basic(const struct mdt_info *info, const char *line_prefix)
{
  if (strlen(info->name_string) > 0)
  {
    fprintf(stdout, "%s ++ >%s<:: pack_size %lld; lb %lld, ext %lld; true: lb %lld, ext %lld\n", line_prefix,
            info->name_string, (long long)info->pack_size, (long long)info->lb, (long long)info->extent, (long long)info->true_lb, (long long)info->true_extent);
  } else {
    fprintf(stdout, "%s ++ pack_size %lld; lb %lld, ext %lld; true: lb %lld, ext %lld\n", line_prefix,
            (long long)info->pack_size, (long long)info->lb, (long long)info->extent, (long long)info->true_lb, (long long)info->true_extent);
  }
  return;
};

static void print_named(const struct mdt_info *info, const char *line_prefix)
{
  print_basic(info, line_prefix);
  fprintf(stdout, "%s ++ elemental datatype\n", line_prefix);
  return;
};

static void print_dup(const struct mdt_info *info, const char *line_prefix)
{
  print_basic(info, line_prefix);
  fprintf(stdout, "%s ++ copy of child\n", line_prefix);
  return;
};

static void print_contiguous(const struct mdt_info *info, const char *line_prefix)
{
  print_basic(info, line_prefix);
  fprintf(stdout, "%s ++ %d adjacent replicas of child\n", line_prefix, info->mdt_intargs[0]);
  return;
};

static void print_vector(const struct mdt_info *info, const char *line_prefix)
{
  print_basic(info, line_prefix);
  fprintf(stdout, "%s ++ %d blocks of %d replicas of child with stride %d\n", line_prefix, 
          info->mdt_intargs[0], info->mdt_intargs[1], info->mdt_intargs[2]);
  return;
};

static void print_hvector_int(const struct mdt_info *info, const char *line_prefix)
{
  print_basic(info, line_prefix);
  fprintf(stdout, "%s ++ %d blocks of %d replicas of child with byte-stride %d (Fortran)\n", line_prefix, 
          info->mdt_intargs[0], info->mdt_intargs[1], info->mdt_intargs[2]);
  return;
};

static void print_hvector(const struct mdt_info *info, const char *line_prefix)
{
  print_basic(info, line_prefix);
  fprintf(stdout, "%s ++ %d blocks of %d replicas of child with byte-stride %lld\n", line_prefix,
          info->mdt_intargs[0], info->mdt_intargs[1], (long long)info->mdt_addrargs[0]);
  return;
};

static void print_indexed(const struct mdt_info *info, const char *line_prefix)
{
  print_basic(info, line_prefix);
  fprintf(stdout, "%s ++ %d blocks of replicas of child with displacements\n", line_prefix,
          info->mdt_intargs[0]);
  return;
};

static void print_hindexed_int(const struct mdt_info *info, const char *line_prefix)
{
  print_basic(info, line_prefix);
  fprintf(stdout, "%s ++ %d blocks of replicas of child with byte-displacements (Fortran)\n", line_prefix,
          info->mdt_intargs[0]);
  return;
};

static void print_hindexed(const struct mdt_info *info, const char *line_prefix)
{
  print_basic(info, line_prefix);
  fprintf(stdout, "%s ++ %d blocks of replicas of child with byte-displacements\n", line_prefix,
          info->mdt_intargs[0]);
  return;
};

static void print_indexed_block(const struct mdt_info *info, const char *line_prefix)
{
  print_basic(info, line_prefix);
  fprintf(stdout, "%s ++ %d blocks of %d replicas of child with displacements\n", line_prefix,
          info->mdt_intargs[0], info->mdt_intargs[1]);
  return;
};

static void print_struct_int(const struct mdt_info *info, const char *line_prefix)
{
  print_basic(info, line_prefix);
  fprintf(stdout, "%s ++ %d blocks of replicas of child types with byte-displacements (Fortran)\n", line_prefix,
          info->mdt_intargs[0]);
  return;
};

static void print_struct(const struct mdt_info *info, const char *line_prefix)
{
  print_basic(info, line_prefix);
  fprintf(stdout, "%s ++ %d blocks of replicas of child types with byte-displacements\n", line_prefix,
          info->mdt_intargs[0]);
  return;
};

static char *  print_int_array(int size, int * array) {

  size_t tmp_size = 3 + size * (1 + 3 * sizeof(int)) + 1;
  char * tmp = malloc(tmp_size);

  snprintf(tmp, tmp_size, "{");
  for (int i = 0; i < size; i++) {
    snprintf(tmp, tmp_size, "%s,%d", tmp, array[i]);
  }
  snprintf(tmp, tmp_size, "}");
  return tmp;
}

static void print_subarray(const struct mdt_info *info, const char *line_prefix)
{
  print_basic(info, line_prefix);
  fprintf(stdout, "%s ++ %d dimensional subarray in %s ordering\n", line_prefix,
          info->mdt_intargs[0],
          (info->mdt_intargs[info->mdt_intargs[0]*3+1] == MPI_ORDER_C ? "C" : "Fortran"));
  char * full = print_int_array(info->mdt_intargs[0], info->mdt_intargs + 1);
  char * sub = print_int_array(info->mdt_intargs[0],
                               info->mdt_intargs + 1 + info->mdt_intargs[0]);
  char * start = print_int_array(info->mdt_intargs[0],
                                 info->mdt_intargs + 1 + 2*info->mdt_intargs[0]);
  fprintf(stdout, "%s ++ full %s sub %s start %s\n", line_prefix, full, sub, start);
  free(full);
  free(sub);
  free(start);
  return;
};

static void print_darray(const struct mdt_info *info, const char *line_prefix)
{
  print_basic(info, line_prefix);
  fprintf(stdout, "%s ++ %d dimensional distributed array\n", line_prefix,
          info->mdt_intargs[2]);
  return;
};

static void print_resized(const struct mdt_info *info, const char *line_prefix)
{
  print_basic(info, line_prefix);
  fprintf(stdout, "%s ++ child resized by lb %lld extent %lld\n", line_prefix,
          (long long)info->mdt_addrargs[0], (long long)info->mdt_addrargs[1]);
  return;
};

static void print_f90_real(const struct mdt_info *info, const char *line_prefix)
{ 
  print_basic(info, line_prefix);
  fprintf(stdout, "%s ++ f90 REAL %d digits %d exponent range\n", line_prefix,
          info->mdt_intargs[0], info->mdt_intargs[1]);
  return;
};

static void print_f90_complex(const struct mdt_info *info, const char *line_prefix)
{
  print_basic(info, line_prefix);
  fprintf(stdout, "%s ++ f90 COMPLEX %d digits %d exponent range\n", line_prefix,
          info->mdt_intargs[0], info->mdt_intargs[1]);
  return;
};

static void print_f90_integer(const struct mdt_info *info, const char *line_prefix)
{
  print_basic(info, line_prefix);
  fprintf(stdout, "%s ++ f90 INTEGER %d digits\n", line_prefix,
          info->mdt_intargs[0]);
  return;
};


static void print_combiner_info(const int combiner_id, const struct mdt_info *info, const char *line_prefix)
{
  switch (combiner_id)
  {
    case MPI_COMBINER_NAMED:
      print_named(info, line_prefix);
      break;
    case MPI_COMBINER_DUP:
      print_dup(info, line_prefix);
      break;
    case MPI_COMBINER_CONTIGUOUS:
      print_contiguous(info, line_prefix);
      break;
    case MPI_COMBINER_VECTOR:
      print_vector(info, line_prefix);
      break;
    case MPI_COMBINER_HVECTOR_INTEGER:
      print_hvector_int(info, line_prefix);
      break;
    case MPI_COMBINER_HVECTOR:
      print_hvector(info, line_prefix);
      break;
    case MPI_COMBINER_INDEXED:
      print_indexed(info, line_prefix);
      break;
    case MPI_COMBINER_HINDEXED_INTEGER:
      print_hindexed_int(info, line_prefix);
      break;
    case MPI_COMBINER_HINDEXED:
      print_hindexed(info, line_prefix);
      break;
    case MPI_COMBINER_INDEXED_BLOCK:
      print_indexed_block(info, line_prefix);
      break;
    case MPI_COMBINER_STRUCT_INTEGER:
      print_struct_int(info, line_prefix);
      break;
    case MPI_COMBINER_STRUCT:
      print_struct(info, line_prefix);
      break;
    case MPI_COMBINER_SUBARRAY:
      print_subarray(info, line_prefix);
      break;
    case MPI_COMBINER_DARRAY:
      print_darray(info, line_prefix);
      break;
    case MPI_COMBINER_F90_REAL:
      print_f90_real(info, line_prefix);
      break;
    case MPI_COMBINER_F90_COMPLEX:
      print_f90_complex(info, line_prefix);
      break;
    case MPI_COMBINER_F90_INTEGER:
      print_f90_integer(info, line_prefix);
      break;
    case MPI_COMBINER_RESIZED:
      print_resized(info, line_prefix);
      break;
    default:
      break;
  };
  return;
};


static void cleanup_mdt_infos(struct mdt_info *info);

static void cleanup_mdt_infos(struct mdt_info *info)
{
  while (info != NULL) {
    if (info->first_child != NULL) {
      cleanup_mdt_infos(info->first_child);
    }
    struct mdt_info * myempty = info;
    if (info->mdt_intargs != NULL)
      free(info->mdt_intargs);
    if (info->mdt_addrargs != NULL)
      free(info->mdt_addrargs);
    if (info->mdt_dtargs != NULL)
      free(info->mdt_dtargs);
    info = info->next;
    free(myempty);
  }
  return;
};

static void analyze_MPI_datatype_nest(const MPI_Datatype mpidt, const int nest, struct mdt_info *info);

static void analyze_MPI_datatype_nest(const MPI_Datatype mpidt, const int nest, struct mdt_info *info)
{
  if (mpidt == MPI_DATATYPE_NULL) return;

  static char const mynameis_format[] = "analyze_MPI_datatype_nest l%d";
  enum { mynameis_size = sizeof(mynameis_format) + 3 * sizeof(nest) + 1 };
  char mynameis[mynameis_size];
  snprintf(mynameis, mynameis_size, mynameis_format, nest);
  int icall = 0, dummy;
  info->mdt = mpidt;
  info->nest = nest;
  wmpi(MPI_Type_size_x(mpidt, &(info->pack_size)), mynameis, icall++);
  wmpi(MPI_Type_get_extent_x(mpidt, &(info->lb), &(info->extent)), mynameis, icall++);
  wmpi(MPI_Type_get_true_extent_x(mpidt, &(info->true_lb), &(info->true_extent)), mynameis, icall++);
  wmpi(MPI_Type_get_name(mpidt, info->name_string,  &(dummy)), mynameis, icall++);
  wmpi(MPI_Type_get_envelope(mpidt, &(info->env_num_ints), &(info->env_num_addr),
                             &(info->env_num_mdt), &(info->env_combiner_id)), mynameis, icall++);
  info->combiner_string = get_combiner_name(info->env_combiner_id);
  info->mdt_intargs = (int *) malloc(sizeof(int)*info->env_num_ints);
  info->mdt_addrargs = (MPI_Aint *) malloc(sizeof(MPI_Aint)*info->env_num_addr);
  info->mdt_dtargs = (MPI_Datatype *) malloc(sizeof(MPI_Datatype)*info->env_num_mdt);
  if (info->env_combiner_id != MPI_COMBINER_NAMED)
  {
    wmpi(MPI_Type_get_contents(mpidt, info->env_num_ints, info->env_num_addr, info->env_num_mdt,
                               info->mdt_intargs, info->mdt_addrargs, info->mdt_dtargs), mynameis, icall++);
  } else {
    info->first_child = NULL;
    info->last_child = NULL;
  }
  if (info->env_num_mdt > 0) {
    int i;
    struct mdt_info * last_info = NULL;
    for(i = 0; i < info->env_num_mdt; i++) {
      struct mdt_info * this_info = (struct mdt_info *) malloc(sizeof(struct mdt_info));
      if (i == 0) {
        info->first_child = this_info;
      }
      if (i == info->env_num_mdt - 1) {
        info->last_child = this_info;
      }
      this_info->previous = last_info;
      this_info->next = NULL;
      if (last_info != NULL) {
        last_info->next = this_info;
      }
      this_info->parent = info;
      analyze_MPI_datatype_nest(info->mdt_dtargs[i], info->nest + 1, this_info);
      last_info = this_info;
    }
  }
  return;
};

static void print_MPI_datatype_tree(struct mdt_info * const info, const char *line_prefix);

static void print_MPI_datatype_tree(struct mdt_info * const info, const char *line_prefix)
{
  static char const nest_prefix_format[] = "%s| ";
  size_t nest_prefix_size = sizeof(nest_prefix_format) + strlen(line_prefix) + 1;
  char * nest_prefix = malloc(nest_prefix_size);
  snprintf(nest_prefix, nest_prefix_size, nest_prefix_format, line_prefix);
  struct mdt_info *tinfo = info;
  while (tinfo != NULL) {
    fprintf(stdout, "%s-->%s\n", nest_prefix, tinfo->combiner_string);
    print_combiner_info(tinfo->env_combiner_id, tinfo, nest_prefix);
    fflush(stdout);
    if (tinfo->first_child != NULL) {
      print_MPI_datatype_tree(tinfo->first_child, nest_prefix);
      fprintf(stdout, "%s----#\n", nest_prefix);
      fflush(stdout);
    }
    tinfo = tinfo->next;
  }
  free(nest_prefix);
  return;
};

static char * get_first_child_name(
  const struct mdt_info * const info, const char * name_prefix) {

  char * first_child_name;
  if(info->first_child->env_combiner_id == MPI_COMBINER_NAMED) {
    static char const first_child_name_format[] = "%s";
    size_t first_child_name_size =
      sizeof(first_child_name_format) + strlen(info->first_child->name_string) + 1;
    first_child_name = malloc(first_child_name_size);
    snprintf(first_child_name, first_child_name_size, first_child_name_format,
             info->first_child->name_string);
  } else {
    static char const first_child_name_format[] = "%s_0";
    size_t first_child_name_size =
      sizeof(first_child_name_format) + strlen(name_prefix) + 1;
    first_child_name = malloc(first_child_name_size);
    snprintf(first_child_name, first_child_name_size, first_child_name_format,
             name_prefix);
  }
  return first_child_name;
}

static void write_null(
  const struct mdt_info * const info, FILE * reconst_file, const char * name_prefix)
{
  fprintf(
    reconst_file, "  MPI_Datatype %s = MPI_DATATYPE_NULL;\n", name_prefix);
  return;
};

static void write_dup(
  const struct mdt_info * const info, FILE * reconst_file, const char * name_prefix)
{
  char * first_child_name = get_first_child_name(info, name_prefix);
  fprintf(reconst_file, "  MPI_Datatype %s;\n  MPI_Type_dup(%s, &%s);\n",
          name_prefix, first_child_name, name_prefix);
  free(first_child_name);
  return;
};

static void write_contiguous(
  const struct mdt_info * const info, FILE * reconst_file, const char * name_prefix)
{
  char * first_child_name = get_first_child_name(info, name_prefix);
  fprintf(reconst_file, "  MPI_Datatype %s;\n  MPI_Type_contiguous(%d, %s, &%s);\n", 
          name_prefix, info->mdt_intargs[0], first_child_name, name_prefix);
  free(first_child_name);
  return;
};

static void write_vector(
  const struct mdt_info * const info, FILE * reconst_file, const char * name_prefix)
{
  char * first_child_name = get_first_child_name(info, name_prefix);
  fprintf(reconst_file, "  MPI_Datatype %s;\n  MPI_Type_vector(%d, %d, %d, %s, &%s);\n",
          name_prefix, info->mdt_intargs[0], info->mdt_intargs[1], info->mdt_intargs[2],
          first_child_name, name_prefix);
  free(first_child_name);
  return;
};

static void write_hvector(
  const struct mdt_info * const info, FILE * reconst_file, const char * name_prefix)
{
  char * first_child_name = get_first_child_name(info, name_prefix);
  fprintf(reconst_file, "  MPI_Datatype %s;\n  MPI_Type_create_hvector(%d, %d, %lld, %s, &%s);\n",
          name_prefix, info->mdt_intargs[0], info->mdt_intargs[1], (long long)info->mdt_addrargs[0],
          first_child_name, name_prefix);
  free(first_child_name);
  return;
};

static void write_hvector_int(
  const struct mdt_info * const info, FILE * reconst_file, const char * name_prefix)
{
  char * first_child_name = get_first_child_name(info, name_prefix);
  fprintf(reconst_file, "# this is not exactly MPI_COMBINER_HVECTOR_INTEGER ...!\n");
  fprintf(reconst_file, "  MPI_Datatype %s;\n  MPI_Type_create_hvector(%d, %d, (MPI_Aint)%d, %s, &%s);\n",
          name_prefix, info->mdt_intargs[0], info->mdt_intargs[1], info->mdt_intargs[2],
          first_child_name, name_prefix);
  free(first_child_name);
  return;
};

static void write_set_arrayvar_int(FILE * reconst_file, const char * arrayname, const char * prefix,
                            const int narray, const int * arrayval)
{
  fprintf(reconst_file, "  int %s%s[%d];\n", prefix, arrayname, narray);
  int l = 512;
  char line[l];
  char line2[l];
  int i;
  snprintf(line, l, "");
  for(i = 0; i < narray; i++) {
    snprintf(line2, l, "%s", line);
    snprintf(line, l, "%s %s%s[%d] = %d;", line2, prefix, arrayname, i, arrayval[i]);
    if(strlen(line) > 160) {
      fprintf(reconst_file, " %s\n", line);
      snprintf(line, l, "");
    }
  }
  if(strlen(line) > 0)
    fprintf(reconst_file, " %s\n", line);
};

static void write_set_arrayvar_addr(FILE * reconst_file, const char * arrayname, const char * prefix,
                             const int narray, const MPI_Aint * arrayval)
{
  fprintf(reconst_file, "  MPI_Aint %s%s[%d];\n", prefix, arrayname, narray);
  int i;
  for(i = 0; i < narray; i++) {
    fprintf(reconst_file, "%s%s[%d] = %lld;\n", prefix, arrayname, i, (long long)arrayval[i]);
  }
};

static void write_set_arrayvar_ddt(FILE * reconst_file, const char * arrayname, const char * prefix,
                             const int narray, const struct mdt_info * info)
{
  fprintf(reconst_file, "  MPI_Datatype %s%s[%d];\n", prefix, arrayname, narray);
  int i;
  struct mdt_info * cinfo = info->first_child;
  for(i = 0; i < narray; i++) {
    if(cinfo->env_combiner_id == MPI_COMBINER_NAMED)
    {
      fprintf(reconst_file, "%s%s[%d] = %s;\n", prefix, arrayname, i, cinfo->name_string);
    } else {
      fprintf(reconst_file, "%s%s[%d] = %s_%d;\n", prefix, arrayname, i, prefix, i);
    }
    cinfo = cinfo->next;
  }
};


static void write_indexed(
  const struct mdt_info * const info, FILE * reconst_file, const char * name_prefix)
{
  char * first_child_name = get_first_child_name(info, name_prefix);
  write_set_arrayvar_int(reconst_file, "_bkl", name_prefix, info->mdt_intargs[0], &(info->mdt_intargs[1]));
  write_set_arrayvar_int(reconst_file, "_dpl", name_prefix, info->mdt_intargs[0], &(info->mdt_intargs[1+info->mdt_intargs[0]]));
  fprintf(reconst_file, "  MPI_Datatype %s;\n  MPI_Type_indexed(%d, %s_bkl, %s_dpl, %s, &%s);\n",
          name_prefix, info->mdt_intargs[0], name_prefix, name_prefix,
          first_child_name, name_prefix);
  free(first_child_name);
  return;
};

static void write_hindexed(
  const struct mdt_info * const info, FILE * reconst_file, const char * name_prefix)
{
  char * first_child_name = get_first_child_name(info, name_prefix);
  write_set_arrayvar_int(reconst_file, "_bkl", name_prefix, info->mdt_intargs[0], &(info->mdt_intargs[1]));
  write_set_arrayvar_addr(reconst_file, "_dpl", name_prefix, info->mdt_intargs[0], info->mdt_addrargs);
  fprintf(reconst_file, "  MPI_Datatype %s;\n  MPI_Type_create_hindexed(%d, %s_bkl, %s_dpl, %s, &%s);\n",
          name_prefix, info->mdt_intargs[0], name_prefix, name_prefix,
          first_child_name, name_prefix);
  free(first_child_name);
  return;
};

static void write_hindexed_int(
  const struct mdt_info * const info, FILE * reconst_file, const char * name_prefix)
{
  char * first_child_name = get_first_child_name(info, name_prefix);
  write_set_arrayvar_int(reconst_file, "_bkl", name_prefix, info->mdt_intargs[0], &(info->mdt_intargs[1]));
  write_set_arrayvar_int(reconst_file, "_dpl", name_prefix, info->mdt_intargs[0], &(info->mdt_intargs[1+info->mdt_intargs[0]]));
  fprintf(reconst_file, "# this is not exactly MPI_COMBINER_HINDEXED_INTEGER ...!\n");
  fprintf(reconst_file, "  MPI_Datatype %s;\n  MPI_Type_create_hindexed(%d, %s_bkl, %s_dpl, %s, &%s);\n",
          name_prefix, info->mdt_intargs[0], name_prefix, name_prefix,
          first_child_name, name_prefix);
  free(first_child_name);
  return;
};

static void write_indexed_block(
  const struct mdt_info * const info, FILE * reconst_file, const char * name_prefix)
{
  char * first_child_name = get_first_child_name(info, name_prefix);
  write_set_arrayvar_int(reconst_file, "_dpl", name_prefix, info->mdt_intargs[0], &(info->mdt_intargs[2]));
  fprintf(reconst_file, "  MPI_Datatype %s;\n  MPI_Type_create_indexed_block(%d, %d, %s_dpl, %s, &%s);\n",
          name_prefix, info->mdt_intargs[0], info->mdt_intargs[1], name_prefix,
          first_child_name, name_prefix);
  free(first_child_name);
  return;
};

static void write_resized(
  const struct mdt_info * const info, FILE * reconst_file, const char * name_prefix)
{
  char * first_child_name = get_first_child_name(info, name_prefix);
  fprintf(reconst_file, "  MPI_Datatype %s;\n  MPI_Type_create_resized(%s, %lld, %lld, &%s);\n",
          name_prefix, first_child_name, (long long)info->mdt_addrargs[0], (long long)info->mdt_addrargs[1],
          name_prefix);
  free(first_child_name);
  return;
};

static void write_f90_real(
  const struct mdt_info * const info, FILE * reconst_file, const char * name_prefix)
{ 
  fprintf(reconst_file, "  MPI_Datatype %s;\n  MPI_Type_create_f90_real(%d, %d, &%s);\n",
          name_prefix, info->mdt_intargs[0], info->mdt_intargs[1], name_prefix);
  return;
};

static void write_f90_complex(
  const struct mdt_info * const info, FILE * reconst_file, const char * name_prefix)
{
  fprintf(reconst_file, "  MPI_Datatype %s;\n  MPI_Type_create_f90_complex(%d, %d, &%s);\n",
          name_prefix, info->mdt_intargs[0], info->mdt_intargs[1], name_prefix);
  return;
};

static void write_f90_integer(
  const struct mdt_info * const info, FILE * reconst_file, const char * name_prefix)
{
  fprintf(reconst_file, "  MPI_Datatype %s;\n  MPI_Type_create_f90_integer(%d, &%s);\n",
          name_prefix, info->mdt_intargs[0], name_prefix);
  return;
};

static void write_struct(
  const struct mdt_info * const info, FILE * reconst_file, const char * name_prefix)
{
  write_set_arrayvar_int(reconst_file, "_bkl", name_prefix, info->mdt_intargs[0], &(info->mdt_intargs[1]));
  write_set_arrayvar_addr(reconst_file, "_dpl", name_prefix, info->mdt_intargs[0], info->mdt_addrargs);
  write_set_arrayvar_ddt(reconst_file, "_dts", name_prefix, info->mdt_intargs[0], info);
  fprintf(reconst_file, "  MPI_Datatype %s;\n  MPI_Type_create_struct(%d, %s_bkl, %s_dpl, %s_dts, &%s);\n",
          name_prefix, info->mdt_intargs[0], name_prefix, name_prefix, name_prefix, name_prefix);
  return;
};

static void write_struct_int(
  const struct mdt_info * const info, FILE * reconst_file, const char * name_prefix)
{
  write_set_arrayvar_int(reconst_file, "_bkl", name_prefix, info->mdt_intargs[0], &(info->mdt_intargs[1]));
  write_set_arrayvar_addr(reconst_file, "_dpl", name_prefix, info->mdt_intargs[0], info->mdt_addrargs);
  write_set_arrayvar_ddt(reconst_file, "_dts", name_prefix, info->mdt_intargs[0], info);
  fprintf(reconst_file, "// this is not exactly MPI_COMBINER_STRUCT_INTEGER ...!\n");
  fprintf(reconst_file, "  MPI_Datatype %s;\n  MPI_Type_create_struct(%d, %s_bkl, %s_dpl, %s_dts, &%s);\n",
          name_prefix, info->mdt_intargs[0], name_prefix, name_prefix, name_prefix, name_prefix);
  return;
};

static void write_subarray(
  const struct mdt_info * const info, FILE * reconst_file, const char * name_prefix)

{
  char * first_child_name = get_first_child_name(info, name_prefix);
  fprintf(reconst_file, "int %s_full[%d], %s_sub[%d], %s_start[%d];\n", name_prefix, info->mdt_intargs[0] ,
          name_prefix, info->mdt_intargs[0], name_prefix, info->mdt_intargs[0]);
  int ict;
  for (ict = 0; ict < info->mdt_intargs[0]; ict++)
  {
    fprintf(reconst_file, "  %s_full[%d] = %d; %s_sub[%d] = %d; %s_start[%d] = %d;\n", 
            name_prefix, ict, info->mdt_intargs[1 + ict], 
            name_prefix, ict, info->mdt_intargs[1 + info->mdt_intargs[0] + ict],
            name_prefix, ict, info->mdt_intargs[1 + 2*info->mdt_intargs[0] + ict]);
  }
  fprintf(reconst_file, "  MPI_Datatype %s;\n  MPI_Type_create_subarray(%d, %s_full, %s_sub, %s_start, %d, %s_0, &%s);\n",
          name_prefix, info->mdt_intargs[0], name_prefix, name_prefix, name_prefix,
          info->mdt_intargs[1 + 3 * info->mdt_intargs[0]], first_child_name, name_prefix);
  free(first_child_name);
  return;
};

static void write_darray(
  const struct mdt_info * const info, FILE * reconst_file, const char * name_prefix)
{
  char * first_child_name = get_first_child_name(info, name_prefix);
  fprintf(reconst_file, "  int %s_glob[%d], %s_distr[%d], %s_darg[%d], %s_psiz[%d];\n", name_prefix, info->mdt_intargs[2] ,
          name_prefix, info->mdt_intargs[2], name_prefix, info->mdt_intargs[2], name_prefix, info->mdt_intargs[2]);
  int ict;
  for (ict = 0; ict < info->mdt_intargs[2]; ict++)
  {
    fprintf(reconst_file, "  %s_glob[%d] = %d; %s_distr[%d] = %d; %s_darg[%d] = %d; %s_psiz[%d] = %d;\n", 
            name_prefix, ict, info->mdt_intargs[3 + ict], 
            name_prefix, ict, info->mdt_intargs[3 + info->mdt_intargs[2] + ict],
            name_prefix, ict, info->mdt_intargs[3 + 2*info->mdt_intargs[2] + ict],
            name_prefix, ict, info->mdt_intargs[3 + 3*info->mdt_intargs[2] + ict]);
  }
  fprintf(reconst_file, "  MPI_Datatype %s;\n  MPI_Type_create_darray(%d, %d, %d, %s_glob, %s_disr, %s_darg, %s_psiz, %d, %s_0, &%s);\n",
          name_prefix, info->mdt_intargs[0], info->mdt_intargs[1], info->mdt_intargs[2], name_prefix, name_prefix,
          name_prefix, name_prefix, info->mdt_intargs[3 + 4 * info->mdt_intargs[2]],
          first_child_name, name_prefix);
  free(first_child_name);
  return;
};

static void write_combiner_code(
  const int combiner_id, const struct mdt_info * const info,
  FILE * reconst_file, const char * name_prefix)
{
  switch (combiner_id)
  {
    case MPI_COMBINER_NAMED:
      if (info->mdt == MPI_DATATYPE_NULL)
        write_null(info, reconst_file, name_prefix);
      break;
    case MPI_COMBINER_DUP:
      write_dup(info, reconst_file, name_prefix);
      break;
    case MPI_COMBINER_CONTIGUOUS:
      write_contiguous(info, reconst_file, name_prefix);
      break;
    case MPI_COMBINER_VECTOR:
      write_vector(info, reconst_file, name_prefix);
      break;
    case MPI_COMBINER_HVECTOR_INTEGER:
      write_hvector_int(info, reconst_file, name_prefix);
      break;
    case MPI_COMBINER_HVECTOR:
      write_hvector(info, reconst_file, name_prefix);
      break;
    case MPI_COMBINER_INDEXED:
      write_indexed(info, reconst_file, name_prefix);
      break;
    case MPI_COMBINER_HINDEXED_INTEGER:
      write_hindexed_int(info, reconst_file, name_prefix);
      break;
    case MPI_COMBINER_HINDEXED:
      write_hindexed(info, reconst_file, name_prefix);
      break;
    case MPI_COMBINER_INDEXED_BLOCK:
      write_indexed_block(info, reconst_file, name_prefix);
      break;
    case MPI_COMBINER_STRUCT_INTEGER:
      write_struct_int(info, reconst_file, name_prefix);
      break;
    case MPI_COMBINER_STRUCT:
      write_struct(info, reconst_file, name_prefix);
      break;
    case MPI_COMBINER_SUBARRAY:
      write_subarray(info, reconst_file, name_prefix);
      break;
    case MPI_COMBINER_DARRAY:
      write_darray(info, reconst_file, name_prefix);
      break;
    case MPI_COMBINER_F90_REAL:
      write_f90_real(info, reconst_file, name_prefix);
      break;
    case MPI_COMBINER_F90_COMPLEX:
      write_f90_complex(info, reconst_file, name_prefix);
      break;
    case MPI_COMBINER_F90_INTEGER:
      write_f90_integer(info, reconst_file, name_prefix);
      break;
    case MPI_COMBINER_RESIZED:
      write_resized(info, reconst_file, name_prefix);
      break;
    default:
      break;
  };
  return;
};

static void reconstruct_MPI_datatype_tree(const struct mdt_info * const info,
                                   FILE *reconst_file,
                                   const char *name_prefix);

static void reconstruct_MPI_datatype_tree(const struct mdt_info * const info,
                                   FILE *reconst_file,
                                   const char *name_prefix)
{
  const struct mdt_info *tinfo = info;
  int icombiner = 0;
  static char const curr_name_prefix_format[] = "%s_%d";
  size_t curr_name_prefix_size =
    sizeof(curr_name_prefix_format) + strlen(name_prefix) +
    3 * sizeof(icombiner) + 1;
  char * curr_name_prefix = malloc(curr_name_prefix_size);
  while (tinfo != NULL) {
    snprintf(curr_name_prefix, curr_name_prefix_size,
             curr_name_prefix_format, name_prefix, icombiner);
    if (tinfo->first_child != NULL)
      reconstruct_MPI_datatype_tree(
        tinfo->first_child, reconst_file, curr_name_prefix);
    write_combiner_code(
      tinfo->env_combiner_id, tinfo, reconst_file, curr_name_prefix);
    icombiner++;
    tinfo = tinfo->next;
  }
  free(curr_name_prefix);
  return;
};

static void write_cleanup_tree(const struct mdt_info * const info,
                        FILE *reconst_file,
                        const char *name_prefix,
                        int skip);

static void write_cleanup_tree(const struct mdt_info * const info,
                        FILE *reconst_file,
                        const char *name_prefix,
                        int skip)
{
  const struct mdt_info * tinfo = info;
  int icombiner = 0;
  static char const curr_name_prefix_format[] = "%s_%d";
  size_t curr_name_prefix_size =
    sizeof(curr_name_prefix_format) + strlen(name_prefix) +
    3 * sizeof(icombiner) + 1;
  char * curr_name_prefix = malloc(curr_name_prefix_size);
  while (tinfo != NULL) {
    snprintf(curr_name_prefix, curr_name_prefix_size,
             curr_name_prefix_format, name_prefix, icombiner);
    if (tinfo->first_child != NULL)
      write_cleanup_tree(
        tinfo->first_child, reconst_file, curr_name_prefix, 0);
    if (skip == 0)
      if (info->env_combiner_id != MPI_COMBINER_NAMED)
        fprintf(reconst_file, "  MPI_Type_free(&%s);\n", curr_name_prefix);
    icombiner++;
    tinfo = tinfo->next;
  }
  free(curr_name_prefix);
  return;
};


static void reconstruct_MPI_datatype(struct mdt_info * const info,
                              const char *reconst_file_name_in,
                              const char * datatype_name_in)
{
  size_t datatype_name_len = strlen(datatype_name_in);
  size_t reconst_file_name_len = strlen(reconst_file_name_in);
  static char const default_datatype_name[] = "MyReconstMPI_Datatype";
  char const * datatype_name =
    datatype_name_len > 0 ? datatype_name_in : default_datatype_name;
  char  * lock_file_name = malloc(reconst_file_name_len + 6 + 1);
  char  * reconst_cfile_name = malloc(reconst_file_name_len + 3 + 1);
  char  * reconst_hfile_name = malloc(reconst_file_name_len + 3 + 1);
  char dummy[64];
  snprintf(lock_file_name,     reconst_file_name_len + 6 + 1, "%s_lock", reconst_file_name_in);
  snprintf(reconst_cfile_name, reconst_file_name_len + 3 + 1, "%s.c",    reconst_file_name_in);
  snprintf(reconst_hfile_name, reconst_file_name_len + 3 + 1, "%s.h",    reconst_file_name_in);

  FILE * lock_file = fopen(lock_file_name, "w+");
  flock(fileno(lock_file), LOCK_EX);

  FILE * reconst_cfile = fopen(reconst_cfile_name, "a");
  FILE * reconst_hfile = fopen(reconst_hfile_name, "a");
  fprintf(reconst_hfile, "\n#ifndef MPIH_DEFINED\n#include <mpi.h>\n#define MPIH_DEFINED\n#endif\n\n");
  fprintf(reconst_cfile, "\n#ifndef MPIH_DEFINED\n#include <mpi.h>\n#define MPIH_DEFINED\n#endif\n\n");
  fprintf(reconst_cfile, "/* automatically generated reconstruction of MPI_Datatype\n");
  fprintf(reconst_cfile, "   see source file >analyze_MPI_Datatype.h< for details */\n");
  fprintf(reconst_cfile, "MPI_Datatype get_%s(void)\n{\n", 
          ((datatype_name_len > 0) ? datatype_name : "reconstuct_rename_me"));
  fprintf(reconst_cfile, "// packed size - - - - - - : %lld\n// lower bound - - - - - - : %lld\n",
          (long long)info->pack_size, (long long)info->lb);
  fprintf(reconst_cfile, "// extent  - - - - - - - - : %lld\n// min. allocation packed  : %lld\n",
          (long long)info->extent, (long long)info->pack_size);
  fprintf(reconst_cfile, "// min. allocation unpacked: %lld\n",
          (long long)(info->extent+info->lb));
  fprintf(reconst_hfile, "MPI_Datatype get_%s(void);\n\n", 
          ((datatype_name_len > 0) ? datatype_name : "reconstuct_rename_me"));
  reconstruct_MPI_datatype_tree(info, reconst_cfile, datatype_name);
  write_cleanup_tree(info, reconst_cfile, datatype_name, 1);
  fprintf(reconst_cfile, "  return %s_0;\n", datatype_name);
  fprintf(reconst_cfile, "};\n");
  fclose(reconst_cfile);
  fclose(reconst_hfile);
  flock(fileno(lock_file), LOCK_UN);
  fclose(lock_file);
  free(reconst_hfile_name);
  free(reconst_cfile_name);
  free(lock_file_name);
  return;
};

static struct mdt_info * mdt_info_new() {

  struct mdt_info *info =
    (struct mdt_info *) malloc(sizeof(struct mdt_info));
  static char const null_combiner[] = "";
  static char const null_name_string[] = "MPI_DATATYPE_NULL";

  info->combiner_string = null_combiner;
  snprintf(info->name_string, MPI_MAX_OBJECT_NAME, "%s", null_name_string);
  info->mdt = MPI_DATATYPE_NULL;
  info->pack_size = 0;
  info->lb = 0;
  info->extent = 0;
  info->true_lb = 0;
  info->true_extent = 0;
  info->mdt_addrargs = NULL;
  info->mdt_dtargs = NULL;
  info->mdt_intargs = NULL;
  info->next = NULL;
  info->previous = NULL;
  info->parent = NULL;
  info->first_child = NULL;
  info->last_child = NULL;
  info->nest = 0;
  info->env_num_ints = 0;
  info->env_num_addr = 0;
  info->env_num_mdt = 0;
  info->env_combiner_id = MPI_COMBINER_NAMED;
  return info;
}

static void analyze_MPI_datatype(
  const MPI_Datatype mpi_type, const char * reconst, const char * datatype_name) {
  int rank;
  wmpi(MPI_Comm_rank(MPI_COMM_WORLD, &rank), "analyze_MPI_datatype", 0);
  printf("rank %d analyzing provided MPI_Datatype... \n", rank);
  struct mdt_info *info = mdt_info_new();
  analyze_MPI_datatype_nest(mpi_type, 0, info);
  // char line_prefix[1024];
  // snprintf(line_prefix, 1024, "%-6d:: ", rank);
  // print_MPI_datatype_tree(info, line_prefix);
  if (strlen(reconst) > 0)
  {
    reconstruct_MPI_datatype(info, reconst, datatype_name);
  }
  cleanup_mdt_infos(info);
  return;
};

#endif
