include(forloop2.m4)dnl
define({interface_gen_nd},
{forloop({dim},{1},{7},{{    MODULE PROCEDURE $2_$1_}dim{d}
})})dnl
dnl xt_redist_s_exchange_impl_gen(typenametag,type)
define({xt_redist_s_exchange_impl_gen},{
  ! see @ref xt_redist_s_exchange
  SUBROUTINE xt_redist_s_exchange_}{$1}{_1d_as(redist, src_size, src_data, &
       dst_size, dst_data)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER, INTENT(in) :: src_size, dst_size
    $2, TARGET, INTENT(in) :: src_data(src_size)
    $2, TARGET, INTENT(inout) :: dst_data(dst_size)
    TYPE(c_ptr) :: src_data_cptr, dst_data_cptr
ifelse($2,LOGICAL,
{#ifdef HAVE_FC_LOGICAL_INTEROP
})dnl
    src_data_cptr = C_LOC(src_data)
    dst_data_cptr = C_LOC(dst_data)
ifelse($2,{LOGICAL},
{#else
    CALL xt_slice_c_loc(src_data, src_data_cptr)
    CALL xt_slice_c_loc(dst_data, dst_data_cptr)
#endif
})dnl
    CALL xt_redist_s_exchange1(redist, src_data_cptr, dst_data_cptr)
  END SUBROUTINE xt_redist_s_exchange_}{$1}{_1d_as
forloop({dim},{1},{7},{
  ! see @ref xt_redist_s_exchange
  SUBROUTINE {xt_redist_s_exchange_}$1{_}dim{d(redist, src_data, dst_data)
    TYPE(xt_redist), INTENT(in) :: redist
    $2, TARGET, INTENT(in) :: src_data(:}forloop({pdim},{2},dim,{,:}){)
    $2, TARGET, INTENT(inout) :: dst_data(:}forloop({pdim},{2},dim,{,:}){)

    $2, POINTER :: src_p(:}forloop({pdim},{2},dim,{,:}){), dst_p(:}forloop({pdim},{2},dim,{,:}){)
    $2, TARGET :: dummy(1}forloop({pdim},{2},dim,{,1}){)
    INTEGER :: src_size, dst_size
    src_size = SIZE(src_data)
    dst_size = SIZE(dst_data)
    IF (src_size > 0) THEN
      src_p => src_data
    ELSE
      src_p => dummy
      src_size = 1
    END IF
    IF (dst_size > 0) THEN
      dst_p => dst_data
    ELSE
      dst_p => dummy
      dst_size = 1
    END IF
    CALL xt_redist_s_exchange_}{$1}{_1d_as(redist, src_size, src_p, dst_size, dst_p)
  END SUBROUTINE xt_redist_s_exchange_}{$1}{_}dim{d
}})})dnl
dnl
define({xt_redist_a_exchange_impl_gen},{
  ! see @ref xt_redist_a_exchange
  SUBROUTINE xt_redist_a_exchange_}{$1}{_1d_as(redist, src_size, src_data, &
       dst_size, dst_data, request)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER, INTENT(in) :: src_size, dst_size
    $2, TARGET, INTENT(in) :: src_data(src_size)
    $2, TARGET, INTENT(inout) :: dst_data(dst_size)
    TYPE(xt_request), INTENT(out) :: request

    $2, TARGET :: dummy(1)
    TYPE(c_ptr) :: src_data_cptr, dst_data_cptr
    IF (src_size > 0) THEN
ifelse($2,LOGICAL,
{#ifdef HAVE_FC_LOGICAL_INTEROP
})dnl
      src_data_cptr = C_LOC(src_data)
ifelse($2,{LOGICAL},
{#else
      CALL xt_slice_c_loc(src_data, src_data_cptr)
#endif
})dnl
    ELSE
ifelse($2,LOGICAL,
{#ifdef HAVE_FC_LOGICAL_INTEROP
})dnl
      src_data_cptr = C_LOC(dummy)
ifelse($2,{LOGICAL},
{#else
      CALL xt_slice_c_loc(dummy, src_data_cptr)
#endif
})dnl
    END IF
    IF (dst_size > 0) THEN
ifelse($2,LOGICAL,
{#ifdef HAVE_FC_LOGICAL_INTEROP
})dnl
      dst_data_cptr = C_LOC(dst_data)
ifelse($2,{LOGICAL},
{#else
      CALL xt_slice_c_loc(dst_data, dst_data_cptr)
#endif
})dnl
    ELSE
ifelse($2,LOGICAL,
{#ifdef HAVE_FC_LOGICAL_INTEROP
})dnl
      dst_data_cptr = C_LOC(dummy)
ifelse($2,{LOGICAL},
{#else
      CALL xt_slice_c_loc(dummy, dst_data_cptr)
#endif
})dnl
    END IF
    CALL xt_redist_a_exchange1(redist, src_data_cptr, dst_data_cptr, request)
  END SUBROUTINE xt_redist_a_exchange_}{$1}{_1d_as
forloop({dim},{1},{7},{
  ! see @ref xt_redist_a_exchange
  SUBROUTINE {xt_redist_a_exchange_}$1{_}dim{d(redist, src_data, dst_data, &
       request)
    TYPE(xt_redist), INTENT(in) :: redist
    $2, TARGET, INTENT(in) :: src_data(:}forloop({pdim},{2},dim,{,:}){)
    $2, TARGET, INTENT(inout) :: dst_data(:}forloop({pdim},{2},dim,{,:}){)
    TYPE(xt_request), INTENT(out) :: request

    INTEGER :: src_size, dst_size
    src_size = SIZE(src_data)
    dst_size = SIZE(dst_data)
#ifdef HAVE_FC_IS_CONTIGUOUS
    IF (.NOT. (IS_CONTIGUOUS(src_data) .AND. IS_CONTIGUOUS(dst_data))) &
      CALL xt_abort('arguments to xt_redist_a_exchange must be contiguous!',&
      filename, __LINE__)
#endif
    CALL xt_redist_a_exchange_}{$1}{_1d_as(redist, src_size, src_data, dst_size, &
         dst_data, request)
  END SUBROUTINE xt_redist_a_exchange_}{$1}{_}dim{d
}})})dnl
dnl
dnl
dnl
dnl Local Variables:
dnl f90-continuation-indent: 5
dnl coding: utf-8
dnl mode: f90
dnl indent-tabs-mode: nil
dnl show-trailing-whitespace: t
dnl require-trailing-newline: t
dnl End:
dnl
