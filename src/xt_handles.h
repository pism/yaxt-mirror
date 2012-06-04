#ifndef XT_HANDLES_H
#define XT_HANDLES_H

/**
 * \file xt_handles.h
 * \brief maps integers to data pointers
 *
 * contains utility routines for handle sets
 */


typedef struct xt_handle_set_type **xt_handle_set_type;

/**
 * checks if a given handle is valid within a handle set
 *
 * @param[in] hset       handle set
 * @param[in] handle     handle within handle set
 */
int xt_handle_is_valid(struct xt_handle_set_type *hset, int handle);

/**
 * constructor for handle sets
 *
 * @param[in] cap        start capacity (size of handle space)
 * @return returns an empty handle set of capacity cap
 */
struct xt_handle_set_type *xt_handle_set_new(int cap);

/**
 * destructor
 *
 * @param[in,out] hset handle set
 */
void xt_handle_set_delete(struct xt_handle_set_type *hset);

/**
 * registers user pointer with a new handle
 *
 * @param[in,out] hset     handle set
 * @param[in] p            pointer to user data
 * @return returns integer handle
 */
int xt_handle_new(struct xt_handle_set_type *hset, void *p);

/**
 * unregisters handle
 *
 * @param[in,out] hset     handle set
 * @param[in] handle       handle
 */
void xt_handle_delete(struct xt_handle_set_type *hset, int handle);

/**
 * inquires about user pointer connected with handle
 *
 * @param[in] hset     handle set
 * @param[in] handle       handle
 */
void *xt_handle2pointer(struct xt_handle_set_type *hset, int handle);

#endif
