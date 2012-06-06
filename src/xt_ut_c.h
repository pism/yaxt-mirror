#ifndef XT_UT_C_H
#define XT_UT_C_H

/**
 * \file xt_ut_c.h
 * \brief supportes the unitrans interfaces/wrappers in xt_ut.f90 with yaxt functionality.
 *
 */

/**
 * Initialization of support data for unitrans (start sizes of allocated data).
 * We use the unitrans names. Unlike unitrans, the allocated data can grow.
 *
 * @param[in] decomp_size       number of decompositions
 * @param[in] comm_tmpl_size    number of communication templates
 * @param[in] comm_size         number of transpositions
 * @param[in] debug_lvl         debug level
 * @param[in] mode              default communication mode
 * @param[in] debug_unit        fortran file unit for debug output (not used)
 */
void xt_ut_init(int decomp_size, int comm_tmpl_size, int comm_size, int debug_lvl, int mode, int debug_unit);

/**
 * creates new decomposition and returns access handle
 *
 * @param[in] idx_vec           local index vector
 * @param[in] idx_vec_n         size of idx_vec
 * @param[in] domain_size       number of transpositions
 * @return returns handle
 */
int xt_ut_init_decomposition_1d(int *idx_vec, int idx_vec_n, int domain_size);

/**
 * creates new unitrans compatible oneway transposition decomposition and returns access handle
 *
 * @param[in] decomp_handle_in   access handle for input decomposition
 * @param[in] decomp_handle_out  access handle for output decomposition
 * @param[in] mpi_world          mpi communicator
 * @param[in] icheck_unique      enable/disbale check of uniqueness
 * @return returns handle
 */
int xt_ut_init_oneway_transposition_template(int decomp_handle_in, int decomp_handle_out,
                                             int mpi_world, int icheck_unique);
/**
 * destroys transposition template
 *
 * @param[in] id                 transposition handle
 */
void xt_ut_destroy_transposition_template(int id);

#endif
