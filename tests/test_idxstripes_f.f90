!>
!! @file test_idxstripes_f.f90
!!
!! @copyright Copyright  (C)  2013 Jörg Behrens <behrens@dkrz.de>
!!                                 Moritz Hanke <hanke@dkrz.de>
!!                                 Thomas Jahns <jahns@dkrz.de>
!!
!! @author Jörg Behrens <behrens@dkrz.de>
!!         Moritz Hanke <hanke@dkrz.de>
!!         Thomas Jahns <jahns@dkrz.de>
!!

!
! Keywords:
! Maintainer: Jörg Behrens <behrens@dkrz.de>
!             Moritz Hanke <hanke@dkrz.de>
!             Thomas Jahns <jahns@dkrz.de>
! URL: https://doc.redmine.dkrz.de/yaxt/html/
!
! Redistribution and use in source and binary forms, with or without
! modification, are  permitted provided that the following conditions are
! met:
!
! Redistributions of source code must retain the above copyright notice,
! this list of conditions and the following disclaimer.
!
! Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the following disclaimer in the
! documentation and/or other materials provided with the distribution.
!
! Neither the name of the DKRZ GmbH nor the names of its contributors
! may be used to endorse or promote products derived from this software
! without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
! IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
! OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
PROGRAM test_idxstripes_f
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort
  USE mpi
  USE test_idxlist_utils, ONLY: check_idxlist, test_err_count, &
       idxlist_pack_unpack_copy
  USE yaxt, ONLY: xt_initialize, xt_finalize, xt_bounds, xt_pos_ext, &
       xt_stripe, xt_idxlist, xt_idxlist_delete, xt_idxstripes_new, &
       xt_idxvec_from_stripes_new, xt_int_kind, xt_idxlist_copy, &
       xt_idxvec_from_stripes_new, &
       xt_idxlist_get_index_stripes, xt_idxlist_get_intersection, &
       xt_idxlist_get_index_at_position, xt_idxlist_get_indices_at_positions, &
       xt_idxlist_get_bounding_box, OPERATOR(/=), &
       xt_idxlist_get_pos_exts_of_index_stripes, &
       xt_idxlist_get_num_indices, xt_idxvec_new, &
       xt_idxstripes_from_idxlist_new
  USE iso_c_binding, ONLY: c_int
  IMPLICIT NONE
  INTEGER, PARAMETER :: xi = xt_int_kind

  CALL init_mpi
  CALL xt_initialize(mpi_comm_world)
  CALL stripe_test_general1
  CALL stripe_test_general2
  CALL stripe_test_general3
  CALL stripe_test_general4
  CALL stripe_test_general5
  CALL test_intersection1
  CALL test_intersection2
  CALL test_intersection3
  CALL test_intersection4
  CALL test_intersection5
  CALL test_intersection6
  CALL test_intersection7
  CALL test_intersection8
  CALL test_intersection9
  CALL test_intersection10
  CALL test_intersection11
  CALL test_intersection12
  CALL test_intersection13
  CALL test_intersection14
  CALL test_intersection15
  CALL test_intersection_stripe2vec
  CALL test_idxlist_stripes_pos_ext1
  CALL test_idxlist_stripes_pos_ext2
  CALL test_idxlist_stripes_pos_ext3
  CALL test_idxlist_stripes_pos_ext_randomized1
  CALL test_get_pos1
  CALL test_get_pos2
  CALL test_get_pos3
  CALL test_get_pos4
  CALL test_stripe_overlap
  CALL test_stripe_bb1
  CALL test_stripe_bb2
  CALL check_pos_ext1
  CALL check_pos_ext2
  CALL check_pos_ext3
  CALL check_pos_ext4
  CALL check_pos_ext5
  CALL check_pos_ext6
  CALL check_pos_ext7
  CALL check_pos_ext8
  IF (test_err_count() /= 0) &
       CALL test_abort("non-zero error count!", &
       __FILE__, &
       __LINE__)
  CALL xt_finalize
  CALL finish_mpi

CONTAINS
  SUBROUTINE stripe_test_general(stripes, ref_indices)
    TYPE(xt_stripe), INTENT(in) :: stripes(:)
    INTEGER(xt_int_kind), INTENT(in) :: ref_indices(:)

    TYPE(xt_idxlist) :: idxstripes, idxvec
    INTEGER :: num_ext, num_unmatched, num_pos, i
    INTEGER(c_int) :: ext_size
    TYPE(xt_pos_ext), ALLOCATABLE :: pos_ext(:)

    idxstripes = xt_idxstripes_new(stripes, SIZE(stripes))
    CALL do_tests(idxstripes, ref_indices)

    num_unmatched = xt_idxlist_get_pos_exts_of_index_stripes(idxstripes, &
         stripes, pos_ext, .TRUE.)
    IF (num_unmatched /= 0) &
         CALL test_abort("stripes not found", &
         __FILE__, &
         __LINE__)

    num_pos = 0
    num_ext = SIZE(pos_ext)
    DO i = 1, num_ext
      ext_size = pos_ext(i)%size
      IF (num_pos /= pos_ext(i)%start) &
           CALL test_abort("position/start mismatch", &
           __FILE__, &
           __LINE__)
      num_pos = num_pos + ext_size
    END DO
    IF (num_pos /= xt_idxlist_get_num_indices(idxstripes)) &
         CALL test_abort("index list length/positions overlap mismatch", &
         __FILE__, &
         __LINE__)

    DEALLOCATE(pos_ext)
    CALL xt_idxlist_delete(idxstripes)

    ! test recreation of stripes from reference vector
    idxvec = xt_idxvec_new(ref_indices)
    idxstripes = xt_idxstripes_from_idxlist_new(idxvec)
    CALL check_idxlist(idxstripes, ref_indices)
    CALL xt_idxlist_delete(idxvec)
    CALL xt_idxlist_delete(idxstripes)
  END SUBROUTINE stripe_test_general

  SUBROUTINE stripe_test_general1
    TYPE(xt_stripe), PARAMETER :: stripes(3) = (/ xt_stripe(0, 1, 5), &
         xt_stripe(10, 1, 5), xt_stripe(20, 1, 5) /);
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(15) &
         = (/ 0_xi,  1_xi,  2_xi,  3_xi,  4_xi, &
         &   10_xi, 11_xi, 12_xi, 13_xi, 14_xi, &
         &   20_xi, 21_xi, 22_xi, 23_xi, 24_xi /)
    CALL stripe_test_general(stripes, ref_indices)
  END SUBROUTINE stripe_test_general1

  SUBROUTINE stripe_test_general2
    TYPE(xt_stripe), PARAMETER :: stripes(3) = (/ xt_stripe(0, 1, 5), &
         xt_stripe(10, 2, 5), xt_stripe(20, 3, 5) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(15) &
         = (/ 0_xi,  1_xi,  2_xi,  3_xi,  4_xi, &
         &   10_xi, 12_xi, 14_xi, 16_xi, 18_xi, &
         &   20_xi, 23_xi, 26_xi, 29_xi, 32_xi /)
    CALL stripe_test_general(stripes, ref_indices)
  END SUBROUTINE stripe_test_general2

  SUBROUTINE stripe_test_general3
    TYPE(xt_stripe), PARAMETER :: stripes(2) = (/ xt_stripe(0, 6, 5), &
         xt_stripe(1, 3, 5) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(10) &
         = (/ 0_xi, 6_xi, 12_xi, 18_xi, 24_xi, &
         &    1_xi, 4_xi,  7_xi, 10_xi, 13_xi /)
    CALL stripe_test_general(stripes, ref_indices)
  END SUBROUTINE stripe_test_general3

  SUBROUTINE stripe_test_general4
    TYPE(xt_stripe), PARAMETER :: stripes(2) = (/ xt_stripe(0, -1, 5), &
         xt_stripe(1, 1, 5) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(10) &
         = (/ 0_xi, -1_xi, -2_xi, -3_xi, -4_xi, &
         &    1_xi,  2_xi,  3_xi,  4_xi,  5_xi /)
    CALL stripe_test_general(stripes, ref_indices)
  END SUBROUTINE stripe_test_general4

  SUBROUTINE stripe_test_general5
    TYPE(xt_stripe), PARAMETER :: stripes(2) = (/ xt_stripe(9, -2, 5), &
         xt_stripe(0, 2, 5) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(10) &
         = (/ 9_xi,  7_xi,  5_xi,  3_xi,  1_xi, &
         &    0_xi,  2_xi,  4_xi,  6_xi,  8_xi /)
    CALL stripe_test_general(stripes, ref_indices)
  END SUBROUTINE stripe_test_general5

  SUBROUTINE test_intersection(stripes_a, stripes_b, ref_indices_a, ref_indices_b)
    TYPE(xt_stripe), INTENT(in) :: stripes_a(:), stripes_b(:)
    INTEGER(xt_int_kind), INTENT(in) :: ref_indices_a(:)
    INTEGER(xt_int_kind), OPTIONAL, INTENT(in) :: ref_indices_b(:)
    TYPE(xt_idxlist) :: idxstripes_a, idxstripes_b, intersection(2)

    idxstripes_a = xt_idxstripes_new(stripes_a)
    idxstripes_b = xt_idxstripes_new(stripes_b)
    intersection(1) = xt_idxlist_get_intersection(idxstripes_a, idxstripes_b)
    intersection(2) = xt_idxlist_get_intersection(idxstripes_b, idxstripes_a)
    CALL do_tests(intersection(1), ref_indices_a)
    IF (PRESENT(ref_indices_b)) THEN
      CALL do_tests(intersection(2), ref_indices_b)
    ELSE
      CALL do_tests(intersection(2), ref_indices_a)
    END IF
    CALL xt_idxlist_delete(intersection(2))
    CALL xt_idxlist_delete(intersection(1))
    CALL xt_idxlist_delete(idxstripes_a)
    CALL xt_idxlist_delete(idxstripes_b)
  END SUBROUTINE test_intersection

  SUBROUTINE test_intersection1
    TYPE(xt_stripe), PARAMETER :: stripes_a(2) = (/ xt_stripe(0, 1, 4), &
         xt_stripe(6, 1, 4) /), &
         stripes_b(1) = (/ xt_stripe(1, 1, 8) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(6) &
         = (/ 1_xi, 2_xi, 3_xi, 6_xi, 7_xi, 8_xi /)
    CALL test_intersection(stripes_a, stripes_b, ref_indices)
  END SUBROUTINE test_intersection1

  SUBROUTINE test_intersection2
    TYPE(xt_stripe), PARAMETER :: stripes_a(3) = (/ xt_stripe(0, 1, 4), &
         xt_stripe(6, 1, 4), xt_stripe(11, 1, 4) /), &
         stripes_b(2) = (/ xt_stripe(1, 1, 7), xt_stripe(9, 1, 5) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(9) &
         = (/ 1_xi, 2_xi, 3_xi, 6_xi, 7_xi, 9_xi, 11_xi, 12_xi, 13_xi /)
    CALL test_intersection(stripes_a, stripes_b, ref_indices)
  END SUBROUTINE test_intersection2

  SUBROUTINE test_intersection3
    TYPE(xt_stripe), PARAMETER :: stripes_a(2) = (/ xt_stripe(0, 1, 3), &
         xt_stripe(8, 1, 3) /), &
         stripes_b(2) = (/ xt_stripe(3, 1, 5), xt_stripe(11, 1, 3) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(1) = (/ -1_xi /)
    CALL test_intersection(stripes_a, stripes_b, ref_indices(1:0))
  END SUBROUTINE test_intersection3

  SUBROUTINE test_intersection4
    TYPE(xt_stripe), PARAMETER :: stripes_a(1) = (/ xt_stripe(0, 1, 10) /), &
         stripes_b(2) = (/ xt_stripe(0, 2, 5), xt_stripe(9, -2, 5) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(10) &
         = (/ 0_xi, 1_xi, 2_xi, 3_xi, 4_xi, 5_xi, 6_xi, 7_xi, 8_xi, 9_xi /)
    CALL test_intersection(stripes_a, stripes_b, ref_indices)
  END SUBROUTINE test_intersection4

  SUBROUTINE test_intersection5
    TYPE(xt_stripe), PARAMETER :: stripes_a(2) = (/ xt_stripe(0, 3, 5), &
         xt_stripe(1, 7, 5) /), &
         stripes_b(2) = (/ xt_stripe(0, 2, 7), xt_stripe(24, -1, 10) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(6) &
         = (/ 0_xi, 6_xi, 8_xi, 12_xi, 15_xi, 22_xi /)
    CALL test_intersection(stripes_a, stripes_b, ref_indices)
  END SUBROUTINE test_intersection5

  SUBROUTINE test_intersection6
    TYPE(xt_stripe), PARAMETER :: stripes_a(1) = (/ xt_stripe(0, 1, 10) /), &
         stripes_b(2) = (/ xt_stripe(5, 1, 5), xt_stripe(4, -1, 5) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(10) &
         = (/ 0_xi, 1_xi, 2_xi, 3_xi, 4_xi, 5_xi, 6_xi, 7_xi, 8_xi, 9_xi /)
    CALL test_intersection(stripes_a, stripes_b, ref_indices)
  END SUBROUTINE test_intersection6

  SUBROUTINE test_intersection7
    TYPE(xt_stripe), PARAMETER :: stripes_a(2) = (/ xt_stripe(0, 1, 10) , &
            xt_stripe(20, 1, 5) /), &
         stripes_b(2) = (/ xt_stripe(3, 1, 5), xt_stripe(17, 1, 5) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(7) &
         = (/ 3_xi, 4_xi, 5_xi, 6_xi, 7_xi, 20_xi, 21_xi /)
    CALL test_intersection(stripes_a, stripes_b, ref_indices)
  END SUBROUTINE test_intersection7

  SUBROUTINE test_intersection8
    TYPE(xt_stripe), PARAMETER :: stripes_a(10) = (/ xt_stripe(0, 1, 2), &
         xt_stripe(3, 1, 2), xt_stripe(5, 1, 2), xt_stripe(8, 1, 2), &
         xt_stripe(10, 1, 2), xt_stripe(14, 1, 2), xt_stripe(17, 1, 2), &
         xt_stripe(20, 1, 2), xt_stripe(23, 1, 2), xt_stripe(25, 1, 2) /), &
         stripes_b(5) = (/ xt_stripe(5, 1, 3), xt_stripe(8, 1, 2), &
         xt_stripe(19, 1, 1), xt_stripe(20, 1, 2), xt_stripe(30, 1, 2) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(6) &
         = (/ 5_xi, 6_xi, 8_xi, 9_xi, 20_xi, 21_xi /)
    CALL test_intersection(stripes_a, stripes_b, ref_indices)
  END SUBROUTINE test_intersection8

  SUBROUTINE test_intersection9
    TYPE(xt_stripe), PARAMETER :: stripes_a(3) = (/ xt_stripe(0, 1, 5), &
         xt_stripe(1, 1, 5), xt_stripe(2, 1, 5) /), &
         stripes_b(1) = (/ xt_stripe(-2, 1, 10) /)
    INTEGER(xt_int_kind) :: i
    INTEGER(xt_int_kind), PARAMETER :: ref_indices_a(7) &
         = (/ (i, i=0_xi,6_xi) /), &
         ref_indices_b(15) = (/ 0_xi, 1_xi, 1_xi, 2_xi, 2_xi, 2_xi, 3_xi, &
         &                      3_xi, 3_xi, 4_xi, 4_xi, 4_xi, 5_xi, 5_xi, 6_xi /)
    CALL test_intersection(stripes_a, stripes_b, ref_indices_a, ref_indices_b)
  END SUBROUTINE test_intersection9

  SUBROUTINE test_intersection10
    TYPE(xt_stripe), PARAMETER :: stripes_a(1) = (/ xt_stripe(0, 2, 5) /), &
         stripes_b(1) = (/ xt_stripe(1, 2, 5) /)
    INTEGER(xt_int_kind), PARAMETER :: dummy(1) = (/ -1_xi /)
    CALL test_intersection(stripes_a, stripes_b, dummy(1:0))
  END SUBROUTINE test_intersection10

  SUBROUTINE test_intersection11
    TYPE(xt_stripe), PARAMETER :: stripes_a(1) = (/ xt_stripe(0, 5, 20) /), &
         stripes_b(1) = (/ xt_stripe(1, 7, 15) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(3) = (/ 15_xi, 50_xi, 85_xi /)
    CALL test_intersection(stripes_a, stripes_b, ref_indices)
  END SUBROUTINE test_intersection11

  ! both ranges overlap in range but have no
  ! indices in common because of stride
  SUBROUTINE test_intersection12
    TYPE(xt_stripe), PARAMETER :: stripes_a(1) = (/ xt_stripe(34, 29, 12) /), &
         stripes_b(1) = (/ xt_stripe(36, 7, 2) /)
    INTEGER(xt_int_kind), PARAMETER :: dummy(1) = (/ -1_xi /)

    CALL test_intersection(stripes_a, stripes_b, dummy(1:0))
  END SUBROUTINE test_intersection12

  ! same as test_intersection12 but with negative stride
  SUBROUTINE test_intersection13
    TYPE(xt_stripe), PARAMETER :: &
         stripes_a(1) = (/ xt_stripe(353, -29, 12) /), &
         stripes_b(1) = (/ xt_stripe(36, 7, 2) /)
    INTEGER(xt_int_kind), PARAMETER :: dummy(1) = (/ -1_xi /)

    CALL test_intersection(stripes_a, stripes_b, dummy(1:0))
  END SUBROUTINE test_intersection13

  SUBROUTINE test_intersection14
    TYPE(xt_stripe), PARAMETER :: &
         stripes_a(1) = (/ xt_stripe(95, -29, 2) /), &
         stripes_b(1) = (/ xt_stripe(81, 14, 2) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(1) = (/ 95_xi /)

    CALL test_intersection(stripes_a, stripes_b, ref_indices)
  END SUBROUTINE test_intersection14

  SUBROUTINE test_intersection15
    TYPE(xt_stripe), PARAMETER :: &
         stripes_a(1) = (/ xt_stripe(546, 14, 2) /), &
         stripes_b(1) = (/ xt_stripe(354, 206, 2) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(1) = (/ 560_xi /)

    CALL test_intersection(stripes_a, stripes_b, ref_indices)
  END SUBROUTINE test_intersection15

  SUBROUTINE test_intersection_stripe2vec
    INTEGER, PARAMETER :: num_stripes = 3
    TYPE(xt_stripe), PARAMETER :: stripes(num_stripes) &
         = (/ xt_stripe(4, 1, 1), xt_stripe(5, 1, 1), xt_stripe(10, -10, 2) /)
    TYPE(xt_idxlist) :: idxvec_a, idxvec_b, intersection
    INTEGER(xt_int_kind), PARAMETER :: index_vector(1) = (/ 5_xi /)
    INTEGER(xt_int_kind) :: intersection_idx
    LOGICAL :: not_found
    idxvec_a = xt_idxvec_from_stripes_new(stripes)
    idxvec_b = xt_idxvec_new(index_vector)
    intersection = xt_idxlist_get_intersection(idxvec_a, idxvec_b)
    IF (xt_idxlist_get_num_indices(intersection) /= 1) &
         CALL test_abort("unexpected number of indices in intersection!", &
         __FILE__, &
         __LINE__)
    not_found = xt_idxlist_get_index_at_position(intersection, 0, &
         intersection_idx)
    IF (intersection_idx /= index_vector(1)) &
         CALL test_abort("unexpected index in intersection!", &
         __FILE__, &
         __LINE__)
    CALL xt_idxlist_delete(intersection)
    CALL xt_idxlist_delete(idxvec_a)
    CALL xt_idxlist_delete(idxvec_b)
  END SUBROUTINE test_intersection_stripe2vec

  SUBROUTINE test_idxlist_stripes_pos_ext1
    INTEGER, PARAMETER :: num_indices = 223
    INTEGER(xt_int_kind), PARAMETER :: index_vector(num_indices) = (/ &
         3375_xi, 3376_xi, 3379_xi, 3380_xi, 3381_xi, 3387_xi, 3388_xi, &
         3389_xi, 3390_xi, 3391_xi, 3392_xi, 3393_xi, 3421_xi, 3422_xi, &
         3423_xi, 3424_xi, 3425_xi, 3426_xi, 3427_xi, 3444_xi, 3458_xi, &
         3459_xi, 3461_xi, 3462_xi, 3463_xi, 3464_xi, 3465_xi, 3466_xi, &
         3467_xi, 3468_xi, 3469_xi, 3470_xi, 3471_xi, 3472_xi, 3473_xi, &
         3474_xi, 3475_xi, 3476_xi, 3477_xi, 3478_xi, 3479_xi, 3480_xi, &
         3529_xi, 3606_xi, 3607_xi, 3608_xi, 3611_xi, 3612_xi, 3613_xi, &
         3614_xi, 3617_xi, 3620_xi, 3621_xi, 3622_xi, 3623_xi, 3624_xi, &
         3625_xi, 3626_xi, 3627_xi, 3628_xi, 3629_xi, 3630_xi, 3631_xi, &
         3684_xi, 3685_xi, 3686_xi, 3687_xi, 3688_xi, 3689_xi, 3690_xi, &
         3691_xi, 3692_xi, 3693_xi, 3694_xi, 3695_xi, 3696_xi, 3697_xi, &
         3698_xi, 3699_xi, 3700_xi, 3701_xi, 3702_xi, 3703_xi, 3704_xi, &
         3705_xi, 3706_xi, 3707_xi, 3708_xi, 3709_xi, 3713_xi, 3714_xi, &
         3715_xi, 3716_xi, 3717_xi, 3718_xi, 3719_xi, 3720_xi, 3721_xi, &
         3722_xi, 3723_xi, 3724_xi, 3725_xi, 3726_xi, 3727_xi, 3728_xi, &
         3729_xi, 3730_xi, 3731_xi, 3741_xi, 3742_xi, 3931_xi, 3932_xi, &
         3374_xi, 3382_xi, 3385_xi, 3394_xi, 3404_xi, 3408_xi, 3412_xi, &
         3440_xi, 3443_xi, 3457_xi, 3481_xi, 3483_xi, 3527_xi, 3619_xi, &
         3735_xi, 3743_xi, 3925_xi, 3930_xi, 3377_xi, 3378_xi, 3383_xi, &
         3384_xi, 3386_xi, 3395_xi, 3397_xi, 3398_xi, 3400_xi, 3402_xi, &
         3403_xi, 3407_xi, 3409_xi, 3410_xi, 3413_xi, 3420_xi, 3441_xi, &
         3442_xi, 3445_xi, 3448_xi, 3449_xi, 3451_xi, 3460_xi, 3482_xi, &
         3519_xi, 3520_xi, 3526_xi, 3528_xi, 3530_xi, 3592_xi, 3593_xi, &
         3595_xi, 3596_xi, 3597_xi, 3609_xi, 3610_xi, 3615_xi, 3616_xi, &
         3618_xi, 3644_xi, 3710_xi, 3711_xi, 3712_xi, 3732_xi, 3733_xi, &
         3736_xi, 3737_xi, 3748_xi, 3749_xi, 3753_xi, 3754_xi, 3759_xi, &
         3760_xi, 3766_xi, 3767_xi, 3919_xi, 3920_xi, 3924_xi, 3926_xi, &
         3933_xi, 3934_xi, 2589_xi, 2602_xi, 2680_xi, 3326_xi, 3340_xi, &
         3341_xi, 3396_xi, 3401_xi, 3411_xi, 3414_xi, 3418_xi, 3446_xi, &
         3447_xi, 3450_xi, 3515_xi, 3521_xi, 3525_xi, 3582_xi, 3590_xi, &
         3591_xi, 3594_xi, 3642_xi, 3734_xi, 3738_xi, 3747_xi, 3750_xi, &
         3761_xi, 3765_xi, 3865_xi, 3918_xi, 3923_xi, 3935_xi /)
    INTEGER, PARAMETER :: num_stripes = 26
    TYPE(xt_stripe), PARAMETER :: stripes(num_stripes) = (/ &
         xt_stripe(3326, 14,  2), xt_stripe(3341, 33,  1), &
         xt_stripe(3374,  1, 25), xt_stripe(3400,  1,  5), &
         xt_stripe(3407,  1,  8), xt_stripe(3418,  2,  1), &
         xt_stripe(3420,  1,  8), xt_stripe(3440,  1, 12), &
         xt_stripe(3457,  1, 27), xt_stripe(3515,  4,  1), &
         xt_stripe(3519,  1,  3), xt_stripe(3525,  1,  6), &
         xt_stripe(3582,  8,  1), xt_stripe(3590,  1,  8), &
         xt_stripe(3606,  1, 26), xt_stripe(3642,  2,  2), &
         xt_stripe(3684,  1, 55), xt_stripe(3741,  1,  3), &
         xt_stripe(3747,  1,  4), xt_stripe(3753,  1,  2), &
         xt_stripe(3759,  1,  3), xt_stripe(3765,  1,  3), &
         xt_stripe(3865, 53,  1), xt_stripe(3918,  1,  3), &
         xt_stripe(3923,  1,  4), xt_stripe(3930,  1,  6) /)
    TYPE(xt_idxlist) :: idxlist

    idxlist = xt_idxvec_new(index_vector, num_indices)
    CALL check_idxlist_stripes_pos_ext(idxlist, stripes)
    CALL xt_idxlist_delete(idxlist)
  END SUBROUTINE test_idxlist_stripes_pos_ext1

  SUBROUTINE test_idxlist_stripes_pos_ext2
    INTEGER, PARAMETER :: num_indices = 201
    INTEGER(xt_int_kind), PARAMETER :: index_vector(num_indices) = (/ &
         &  178_xi,  179_xi,  180_xi,  181_xi,  182_xi,  183_xi,  184_xi, &
         &  186_xi,  187_xi,  188_xi,  189_xi,  190_xi,  194_xi,  195_xi, &
         &  196_xi,  197_xi,  198_xi,  199_xi,  200_xi,  201_xi,  202_xi, &
         &  203_xi,  204_xi,  205_xi,  206_xi,  207_xi,  208_xi,  209_xi, &
         &  210_xi,  211_xi,  212_xi,  217_xi,  223_xi,  426_xi,  428_xi, &
         &  429_xi,  430_xi,  434_xi,  435_xi,  436_xi,  437_xi,  438_xi, &
         &  439_xi,  440_xi,  442_xi,  443_xi,  444_xi,  445_xi,  446_xi, &
         &  447_xi,  448_xi,  449_xi,  450_xi,  451_xi,  452_xi,  453_xi, &
         &  454_xi,  455_xi,  456_xi,  457_xi,  458_xi,  670_xi,  671_xi, &
         &  672_xi,  673_xi,  674_xi,  675_xi,  676_xi,  677_xi,  682_xi, &
         &  684_xi,  685_xi,  686_xi,  687_xi,  688_xi,  689_xi,  690_xi, &
         &  692_xi,  695_xi,  703_xi,  704_xi,  705_xi,  706_xi,  707_xi, &
         &  894_xi,  895_xi,  896_xi,  897_xi,  898_xi,  899_xi,  900_xi, &
         &  901_xi,  906_xi,  907_xi,  908_xi,  913_xi,  915_xi,  921_xi, &
         &  922_xi,  923_xi,  924_xi,  925_xi,  926_xi,  927_xi, 1096_xi, &
         & 1097_xi, 1098_xi, 1099_xi, 1100_xi, 1101_xi, 1102_xi, 1103_xi, &
         & 1107_xi, 1108_xi, 1109_xi, 1110_xi, 1111_xi, 1113_xi, 1114_xi, &
         & 1119_xi, 1120_xi, 1121_xi, 2095_xi, 2096_xi, 2097_xi, 2098_xi, &
         & 2100_xi, 2102_xi, 2103_xi, 2104_xi, 2105_xi, 2107_xi, 2108_xi, &
         & 2109_xi, 2110_xi, 2112_xi, 2118_xi, 2120_xi, 2121_xi, 2122_xi, &
         & 2123_xi, 2124_xi, 2125_xi, 2127_xi, 2128_xi, 2129_xi, 2130_xi, &
         & 2134_xi, 2140_xi, 2141_xi, 2142_xi, 2143_xi, 2145_xi, 2148_xi, &
         & 2149_xi, 2151_xi, 2152_xi, 2153_xi, 2154_xi, 2155_xi, 2156_xi, &
         &  683_xi,  691_xi,  903_xi,  914_xi, 1105_xi, 1115_xi, 2099_xi, &
         & 2106_xi, 2111_xi, 2115_xi, 2126_xi, 2132_xi, 2139_xi, 2144_xi, &
         & 2147_xi, 2150_xi, 2305_xi,  427_xi,  465_xi,  466_xi,  678_xi, &
         &  693_xi,  902_xi,  909_xi, 1104_xi, 1112_xi, 2101_xi, 2113_xi, &
         & 2114_xi, 2116_xi, 2117_xi, 2119_xi, 2131_xi, 2136_xi, 2138_xi, &
         & 2146_xi, 2297_xi, 2302_xi, 2304_xi, 2307_xi /)
    INTEGER, PARAMETER :: num_stripes = 8
    TYPE(xt_stripe), PARAMETER :: stripes(num_stripes) = (/ &
      xt_stripe(670, 1,  9), xt_stripe(682, 1, 12), &
      xt_stripe(695, 8,  1), xt_stripe(703, 1,  5), &
      xt_stripe(894, 1, 10), xt_stripe(906, 1,  4), &
      xt_stripe(913, 1,  3), xt_stripe(921, 1,  7) /)
    TYPE(xt_idxlist) :: idxlist

    idxlist = xt_idxvec_new(index_vector)
    CALL check_idxlist_stripes_pos_ext(idxlist, stripes)
    CALL xt_idxlist_delete(idxlist)
  END SUBROUTINE test_idxlist_stripes_pos_ext2

  SUBROUTINE test_idxlist_stripes_pos_ext3
    INTEGER, PARAMETER :: num_indices = 1144
    INTEGER(xt_int_kind), PARAMETER :: index_vector(num_indices) = (/ &
      2055, 2056, 2060, 2193, 2199, 2203, 2211, 2212, 2278, 2281, 2311, &
      2312, 2316, 2317, 2322, 2332, 2447, 2448, 2452, 2585, 2591, 2595, &
      2603, 2604, 2670, 2673, 2703, 2704, 2708, 2709, 2714, 2724, 2839, &
      2840, 2844, 2977, 2983, 2987, 2995, 2996, 3062, 3065, 3095, 3096, &
      3100, 3101, 3106, 3116, 3231, 3232, 3236, 3369, 3375, 3379, 3387, &
      3388, 3454, 3457, 3487, 3488, 3492, 3493, 3498, 3508, 3623, 3624, &
      3628, 3761, 3767, 3771, 3779, 3780, 3846, 3849, 3879, 3880, 3884, &
      3885, 3890, 3900, 3997, 4001, 4002, 4053, 4057, 4084, 4085, 4092, &
      4102, 4188, 4192, 4201, 4373, 4377, 4378, 4429, 4433, 4460, 4461, &
      4468, 4478, 4564, 4568, 4577, 4749, 4753, 4754, 4805, 4809, 4836, &
      4837, 4844, 4854, 4945, 4953, 5125, 5129, 5130, 5181, 5185, 5212, &
      5213, 5220, 5230, 5321, 5329, 5501, 5505, 5506, 5557, 5561, 5588, &
      5589, 5596, 5606, 5697, 5705,  162,  163,  166,  168,  171,  172, &
       173,  177,  181,  362,  363,  367,  369,  375,  378,  382,  383, &
       386,  570,  571,  574,  576,  579,  580,  581,  585,  589,  758, &
       759,  763,  765,  769,  774,  775,  778,  962,  963,  966,  968, &
       971,  972,  973,  977,  981, 1150, 1151, 1155, 1157, 1161, 1166, &
      1167, 1170, 1354, 1355, 1358, 1360, 1363, 1364, 1365, 1369, 1373, &
      1542, 1543, 1547, 1549, 1553, 1558, 1559, 1562, 1746, 1747, 1750, &
      1752, 1755, 1756, 1757, 1761, 1918, 1919, 1923, 1925, 1929, 1934, &
      1935, 1938, 1988, 1989, 2024, 2025, 2032, 2033, 2036, 2038, 2039, &
      2048, 2049, 2053, 2054, 2057, 2058, 2061, 2076, 2077, 2091, 2092, &
      2093, 2095, 2097, 2126, 2127, 2144, 2145, 2149, 2150, 2156, 2198, &
      2204, 2205, 2207, 2245, 2253, 2254, 2256, 2268, 2269, 2277, 2279, &
      2280, 2283, 2287, 2298, 2299, 2307, 2308, 2309, 2310, 2333, 2334, &
      2380, 2381, 2416, 2417, 2424, 2425, 2428, 2430, 2431, 2440, 2441, &
      2445, 2446, 2449, 2450, 2453, 2468, 2469, 2483, 2484, 2485, 2487, &
      2489, 2518, 2519, 2536, 2537, 2541, 2542, 2548, 2590, 2596, 2597, &
      2599, 2637, 2645, 2646, 2648, 2660, 2661, 2669, 2671, 2672, 2675, &
      2679, 2690, 2691, 2699, 2700, 2701, 2702, 2725, 2726, 2772, 2773, &
      2808, 2809, 2816, 2817, 2820, 2822, 2823, 2832, 2833, 2837, 2838, &
      2841, 2842, 2845, 2860, 2861, 2875, 2876, 2877, 2879, 2881, 2910, &
      2911, 2928, 2929, 2933, 2934, 2940, 2982, 2988, 2989, 2991, 3029, &
      3037, 3038, 3040, 3052, 3053, 3061, 3063, 3064, 3067, 3071, 3082, &
      3083, 3091, 3092, 3093, 3094, 3117, 3118, 3164, 3165, 3200, 3201, &
      3208, 3209, 3212, 3214, 3215, 3224, 3225, 3229, 3230, 3233, 3234, &
      3237, 3252, 3253, 3267, 3268, 3269, 3271, 3273, 3302, 3303, 3320, &
      3321, 3325, 3326, 3332, 3374, 3380, 3381, 3383, 3421, 3429, 3430, &
      3432, 3444, 3445, 3453, 3455, 3456, 3459, 3463, 3474, 3475, 3483, &
      3484, 3485, 3486, 3509, 3510, 3556, 3557, 3592, 3593, 3600, 3601, &
      3604, 3606, 3607, 3616, 3617, 3621, 3622, 3625, 3626, 3629, 3644, &
      3645, 3659, 3660, 3661, 3663, 3665, 3694, 3695, 3712, 3713, 3717, &
      3718, 3724, 3766, 3772, 3773, 3775, 3813, 3821, 3822, 3824, 3836, &
      3837, 3845, 3847, 3848, 3851, 3855, 3866, 3867, 3875, 3876, 3877, &
      3878, 3901, 3902, 3948, 3949, 3984, 3985, 3992, 3993, 3996, 3998, &
      3999, 4008, 4009, 4013, 4014, 4017, 4018, 4021, 4036, 4037, 4051, &
      4052, 4054, 4055, 4058, 4090, 4091, 4093, 4108, 4109, 4112, 4113, &
      4114, 4158, 4164, 4165, 4193, 4199, 4200, 4212, 4213, 4222, 4223, &
      4225, 4227, 4231, 4242, 4243, 4250, 4251, 4271, 4272, 4274, 4324, &
      4325, 4360, 4361, 4368, 4369, 4372, 4374, 4375, 4384, 4385, 4389, &
      4390, 4393, 4394, 4397, 4412, 4413, 4427, 4428, 4430, 4431, 4434, &
      4466, 4467, 4469, 4484, 4485, 4488, 4489, 4490, 4534, 4540, 4541, &
      4569, 4575, 4576, 4588, 4589, 4598, 4599, 4601, 4603, 4607, 4618, &
      4619, 4626, 4627, 4647, 4648, 4650, 4700, 4701, 4736, 4737, 4744, &
      4745, 4748, 4750, 4751, 4760, 4761, 4765, 4766, 4769, 4770, 4773, &
      4788, 4789, 4803, 4804, 4806, 4807, 4810, 4842, 4843, 4845, 4860, &
      4861, 4864, 4865, 4866, 4910, 4916, 4917, 4951, 4952, 4964, 4965, &
      4974, 4975, 4977, 4979, 4983, 4994, 4995, 5002, 5003, 5023, 5024, &
      5026, 5076, 5077, 5112, 5113, 5120, 5121, 5124, 5126, 5127, 5136, &
      5137, 5141, 5142, 5145, 5146, 5149, 5164, 5165, 5179, 5180, 5182, &
      5183, 5186, 5218, 5219, 5221, 5236, 5237, 5240, 5241, 5242, 5286, &
      5292, 5293, 5327, 5328, 5340, 5341, 5350, 5351, 5353, 5355, 5359, &
      5370, 5371, 5378, 5379, 5399, 5400, 5402, 5452, 5453, 5488, 5489, &
      5496, 5497, 5500, 5502, 5503, 5512, 5513, 5517, 5518, 5521, 5522, &
      5525, 5540, 5541, 5555, 5556, 5558, 5559, 5562, 5594, 5595, 5597, &
      5612, 5613, 5616, 5617, 5618, 5662, 5668, 5669, 5703, 5704, 5716, &
      5717, 5726, 5727, 5729, 5731, 5735, 5746, 5747, 5754, 5755, 5775, &
      5776, 5778, 5958, 5959, 5962, 5964, 5967, 5968, 5971, 5973, 6154, &
      6155, 6159, 6161, 6167, 6170, 6172, 6173, 6350, 6351, 6354, 6356, &
      6359, 6360, 6363, 6530, 6531, 6535, 6537, 6543, 6546, 6548, 6549, &
      6726, 6727, 6730, 6732, 6735, 6736, 6739, 6906, 6907, 6911, 6913, &
      6919, 6922, 6924, 6925, 7102, 7103, 7106, 7108, 7111, 7112, 7115, &
      7282, 7283, 7287, 7289, 7295, 7298, 7300, 7301, 7478, 7479, 7482, &
      7484, 7487, 7488, 7491, 7646, 7647, 7651, 7653, 7657, 7660, 7661, &
       130,  161,  169,  170,  336,  361,  366,  384,  538,  569,  577, &
       578,  736,  757,  762,  776,  930,  961,  969,  970, 1128, 1149, &
      1154, 1168, 1322, 1353, 1361, 1362, 1520, 1541, 1546, 1560, 1714, &
      1745, 1753, 1754, 1896, 1917, 1922, 1936, 1985, 2019, 2031, 2035, &
      2040, 2044, 2052, 2059, 2062, 2071, 2087, 2090, 2094, 2140, 2148, &
      2153, 2157, 2206, 2257, 2263, 2267, 2284, 2288, 2293, 2295, 2305, &
      2306, 2377, 2411, 2423, 2427, 2432, 2436, 2444, 2451, 2454, 2463, &
      2479, 2482, 2486, 2532, 2540, 2545, 2549, 2598, 2649, 2655, 2659, &
      2676, 2680, 2685, 2687, 2697, 2698, 2769, 2803, 2815, 2819, 2824, &
      2828, 2836, 2843, 2846, 2855, 2871, 2874, 2878, 2924, 2932, 2937, &
      2941, 2990, 3041, 3047, 3051, 3068, 3072, 3077, 3079, 3089, 3090, &
      3161, 3195, 3207, 3211, 3216, 3220, 3228, 3235, 3238, 3247, 3263, &
      3266, 3270, 3316, 3324, 3329, 3333, 3382, 3433, 3439, 3443, 3460, &
      3464, 3469, 3471, 3481, 3482, 3553, 3587, 3599, 3603, 3608, 3612, &
      3620, 3627, 3630, 3639, 3655, 3658, 3662, 3708, 3716, 3721, 3725, &
      3774, 3825, 3831, 3835, 3852, 3856, 3861, 3863, 3873, 3874, 3945, &
      3979, 3991, 3995, 4000, 4004, 4012, 4019, 4022, 4031, 4033, 4047, &
      4050, 4104, 4106, 4115, 4207, 4221, 4228, 4232, 4237, 4249, 4252, &
      4321, 4355, 4367, 4371, 4376, 4380, 4388, 4395, 4398, 4407, 4409, &
      4423, 4426, 4480, 4482, 4491, 4583, 4597, 4604, 4608, 4613, 4625, &
      4628, 4697, 4731, 4743, 4747, 4752, 4756, 4764, 4771, 4774, 4783, &
      4785, 4799, 4802, 4856, 4858, 4867, 4959, 4973, 4980, 4984, 4989, &
      5001, 5004, 5073, 5107, 5119, 5123, 5128, 5132, 5140, 5147, 5150, &
      5159, 5161, 5175, 5178, 5232, 5234, 5243, 5335, 5349, 5356, 5360, &
      5365, 5377, 5380, 5449, 5483, 5495, 5499, 5504, 5508, 5516, 5523, &
      5526, 5535, 5537, 5551, 5554, 5608, 5610, 5619, 5711, 5725, 5732, &
      5736, 5741, 5753, 5756, 5930, 5957, 5965, 5966, 6128, 6153, 6158, &
      6174, 6322, 6349, 6357, 6358, 6504, 6529, 6534, 6550, 6698, 6725, &
      6733, 6734, 6880, 6905, 6910, 6926, 7074, 7101, 7109, 7110, 7256, &
      7281, 7286, 7302, 7450, 7477, 7485, 7486, 7624, 7645, 7650, 7662 /)
    INTEGER, PARAMETER :: num_stripes = 187
    TYPE(xt_stripe), PARAMETER :: stripes(num_stripes) = (/ &
      xt_stripe(173, 408, 2), xt_stripe(973, 392, 3), xt_stripe(1985, 4, 2), &
      xt_stripe(2044, 4, 2), xt_stripe(2049, 3, 1), xt_stripe(2052, 1, 9), &
      xt_stripe(2062, 131, 2), xt_stripe(2198, 1, 2), xt_stripe(2203, 1, 5), &
      xt_stripe(2211, 1, 2), xt_stripe(2263, 4, 1), xt_stripe(2267, 1, 3), &
      xt_stripe(2277, 1, 5), xt_stripe(2283, 1, 2), xt_stripe(2287, 1, 2), &
      xt_stripe(2293, 2, 2), xt_stripe(2298, 1, 2), xt_stripe(2305, 1, 8), &
      xt_stripe(2316, 1, 2), xt_stripe(2322, 10, 1), xt_stripe(2332, 1, 3), &
      xt_stripe(2377, 4, 2), xt_stripe(2436, 4, 2), xt_stripe(2441, 3, 1), &
      xt_stripe(2444, 1, 9), xt_stripe(2454, 131, 2), xt_stripe(2590, 1, 2), &
      xt_stripe(2595, 1, 5), xt_stripe(2603, 1, 2), xt_stripe(2655, 4, 1), &
      xt_stripe(2659, 1, 3), xt_stripe(2669, 1, 5), xt_stripe(2675, 1, 2), &
      xt_stripe(2679, 1, 2), xt_stripe(2685, 2, 2), xt_stripe(2690, 1, 2), &
      xt_stripe(2697, 1, 8), xt_stripe(2708, 1, 2), xt_stripe(2714, 10, 1), &
      xt_stripe(2724, 1, 3), xt_stripe(2769, 4, 2), xt_stripe(2828, 4, 2), &
      xt_stripe(2833, 3, 1), xt_stripe(2836, 1, 9), xt_stripe(2846, 131, 2), &
      xt_stripe(2982, 1, 2), xt_stripe(2987, 1, 5), xt_stripe(2995, 1, 2), &
      xt_stripe(3047, 4, 1), xt_stripe(3051, 1, 3), xt_stripe(3061, 1, 5), &
      xt_stripe(3067, 1, 2), xt_stripe(3071, 1, 2), xt_stripe(3077, 2, 2), &
      xt_stripe(3082, 1, 2), xt_stripe(3089, 1, 8), xt_stripe(3100, 1, 2), &
      xt_stripe(3106, 10, 1), xt_stripe(3116, 1, 3), xt_stripe(3161, 4, 2), &
      xt_stripe(3220, 4, 2), xt_stripe(3225, 3, 1), xt_stripe(3228, 1, 9), &
      xt_stripe(3238, 131, 2), xt_stripe(3374, 1, 2), xt_stripe(3379, 1, 5), &
      xt_stripe(3387, 1, 2), xt_stripe(3439, 4, 1), xt_stripe(3443, 1, 3), &
      xt_stripe(3453, 1, 5), xt_stripe(3459, 1, 2), xt_stripe(3463, 1, 2), &
      xt_stripe(3469, 2, 2), xt_stripe(3474, 1, 2), xt_stripe(3481, 1, 8), &
      xt_stripe(3492, 1, 2), xt_stripe(3498, 10, 1), xt_stripe(3508, 1, 3), &
      xt_stripe(3553, 4, 2), xt_stripe(3612, 4, 2), xt_stripe(3617, 3, 1), &
      xt_stripe(3620, 1, 9), xt_stripe(3630, 131, 2), xt_stripe(3766, 1, 2), &
      xt_stripe(3771, 1, 5), xt_stripe(3779, 1, 2), xt_stripe(3831, 4, 1), &
      xt_stripe(3835, 1, 3), xt_stripe(3845, 1, 5), xt_stripe(3851, 1, 2), &
      xt_stripe(3855, 1, 2), xt_stripe(3861, 2, 2), xt_stripe(3866, 1, 2), &
      xt_stripe(3873, 1, 8), xt_stripe(3884, 1, 2), xt_stripe(3890, 10, 1), &
      xt_stripe(3900, 1, 3), xt_stripe(3945, 3, 2), xt_stripe(3979, 5, 2), &
      xt_stripe(3985, 6, 1), xt_stripe(3991, 1, 3), xt_stripe(3995, 2, 1), &
      xt_stripe(3997, 1, 6), xt_stripe(4031, 2, 2), xt_stripe(4036, 1, 2), &
      xt_stripe(4047, 3, 1), xt_stripe(4050, 1, 6), xt_stripe(4057, 1, 2), &
      xt_stripe(4084, 1, 2), xt_stripe(4090, 1, 4), xt_stripe(4102, 2, 4), &
      xt_stripe(4109, 3, 1), xt_stripe(4112, 1, 4), xt_stripe(4188, 4, 2), &
      xt_stripe(4193, 6, 1), xt_stripe(4199, 1, 3), xt_stripe(4321, 3, 2), &
      xt_stripe(4355, 5, 2), xt_stripe(4361, 6, 1), xt_stripe(4367, 1, 3), &
      xt_stripe(4371, 2, 1), xt_stripe(4373, 1, 6), xt_stripe(4407, 2, 2), &
      xt_stripe(4412, 1, 2), xt_stripe(4423, 3, 1), xt_stripe(4426, 1, 6), &
      xt_stripe(4433, 1, 2), xt_stripe(4460, 1, 2), xt_stripe(4466, 1, 4), &
      xt_stripe(4478, 2, 4), xt_stripe(4485, 3, 1), xt_stripe(4488, 1, 4), &
      xt_stripe(4564, 4, 2), xt_stripe(4569, 6, 1), xt_stripe(4575, 1, 3), &
      xt_stripe(4697, 3, 2), xt_stripe(4731, 5, 2), xt_stripe(4737, 6, 1), &
      xt_stripe(4743, 1, 3), xt_stripe(4747, 2, 1), xt_stripe(4749, 1, 6), &
      xt_stripe(4783, 2, 2), xt_stripe(4788, 1, 2), xt_stripe(4799, 3, 1), &
      xt_stripe(4802, 1, 6), xt_stripe(4809, 1, 2), xt_stripe(4836, 1, 2), &
      xt_stripe(4842, 1, 4), xt_stripe(4854, 2, 4), xt_stripe(4861, 3, 1), &
      xt_stripe(4864, 1, 4), xt_stripe(4945, 6, 1), xt_stripe(4951, 1, 3), &
      xt_stripe(5107, 5, 2), xt_stripe(5113, 6, 1), xt_stripe(5119, 1, 3), &
      xt_stripe(5123, 2, 1), xt_stripe(5125, 1, 6), xt_stripe(5159, 2, 2), &
      xt_stripe(5164, 1, 2), xt_stripe(5175, 3, 1), xt_stripe(5178, 1, 6), &
      xt_stripe(5185, 1, 2), xt_stripe(5212, 1, 2), xt_stripe(5218, 1, 4), &
      xt_stripe(5230, 2, 4), xt_stripe(5237, 3, 1), xt_stripe(5240, 1, 4), &
      xt_stripe(5321, 6, 1), xt_stripe(5327, 1, 3), xt_stripe(5483, 5, 2), &
      xt_stripe(5489, 6, 1), xt_stripe(5495, 1, 3), xt_stripe(5499, 2, 1), &
      xt_stripe(5501, 1, 6), xt_stripe(5535, 2, 2), xt_stripe(5540, 1, 2), &
      xt_stripe(5551, 3, 1), xt_stripe(5554, 1, 6), xt_stripe(5561, 1, 2), &
      xt_stripe(5588, 1, 2), xt_stripe(5594, 1, 4), xt_stripe(5606, 2, 4), &
      xt_stripe(5613, 3, 1), xt_stripe(5616, 1, 4), xt_stripe(5697, 6, 1), &
      xt_stripe(5703, 1, 3) /)
    TYPE(xt_idxlist) :: idxlist

    idxlist = xt_idxvec_new(index_vector)
    CALL check_idxlist_stripes_pos_ext(idxlist, stripes)
    CALL xt_idxlist_delete(idxlist)
  END SUBROUTINE test_idxlist_stripes_pos_ext3

  SUBROUTINE test_idxlist_stripes_pos_ext_randomized1
    INTEGER, PARAMETER :: num_iterations=128, &
         max_num_indices=1024, max_index=1024

    INTEGER, ALLOCATABLE :: rseed(:)
    INTEGER :: rseed_size, i, iteration, num_indices, num_stripes
    INTEGER(xt_int_kind), ALLOCATABLE :: indices(:)
    REAL, ALLOCATABLE :: rvals(:)
    TYPE(xt_idxlist) :: idxlist
    TYPE(xt_stripe), ALLOCATABLE :: stripes(:)

    CALL random_seed(size=rseed_size)
    ALLOCATE(rseed(rseed_size))
    DO i = 1, rseed_size
      rseed(i) = 4711
    END DO
    CALL random_seed(put=rseed)
    ALLOCATE(indices(max_num_indices), rvals(max_num_indices))
    DO iteration = 1, num_iterations
      CALL random_number(rvals(1))
      num_indices = NINT(rvals(1) * REAL(max_num_indices))

      CALL random_number(rvals(1:num_indices))
      DO i = 1, num_indices
        indices(i) = NINT(rvals(i)*REAL((2*max_index)-max_index), xt_int_kind)
      END DO
      idxlist = xt_idxvec_new(indices, num_indices)

      CALL xt_idxlist_get_index_stripes(idxlist, stripes)
      num_stripes = SIZE(stripes)

      CALL check_idxlist_stripes_pos_ext(idxlist, stripes)

      CALL xt_idxlist_delete(idxlist)
    END DO
  END SUBROUTINE test_idxlist_stripes_pos_ext_randomized1

  SUBROUTINE check_idxlist_stripes_pos_ext(idxlist, stripes)
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    TYPE(xt_stripe), INTENT(in) :: stripes(:)

    TYPE(xt_pos_ext), ALLOCATABLE :: pos_ext(:)
    INTEGER :: num_stripes, num_ext, num_unmatched
    INTEGER :: abs_pos_ext_size, jsign, i, j, k, send_pos
    INTEGER(xt_int_kind) :: intersection_index, orig_index
    LOGICAL, PARAMETER :: single_match_only = .TRUE.
    LOGICAL :: unmatched_in_intersection, unmatched_in_idxlist
    TYPE(xt_idxlist) :: intersection
    num_stripes = SIZE(stripes)

    num_unmatched = xt_idxlist_get_pos_exts_of_index_stripes( &
         idxlist, num_stripes, stripes, num_ext, pos_ext, single_match_only)

    ! testing of results
    IF (num_unmatched /= 0) &
         CALL test_abort("error in xt_idxlist_get_pos_exts_of_index_stripes", &
         __FILE__, &
         __LINE__)
    intersection = xt_idxvec_from_stripes_new(stripes)
    k = 0
    DO i = 1, num_ext
      abs_pos_ext_size = INT(ABS(pos_ext(i)%size))
      jsign = MERGE(1, -1, pos_ext(i)%size >= 0)
      DO j = 0, abs_pos_ext_size-1
        unmatched_in_intersection &
             = xt_idxlist_get_index_at_position(intersection, k, &
             intersection_index)
        send_pos = pos_ext(i)%start + jsign * j
        unmatched_in_idxlist &
             = xt_idxlist_get_index_at_position(idxlist, send_pos, orig_index)
        IF (unmatched_in_intersection .OR. unmatched_in_idxlist &
             .OR. intersection_index /= orig_index) THEN
          WRITE (0, '(4(a,i0))') "intersection pos ", k, &
               " index ", intersection_index, &
               " orig pos ", send_pos, &
               " index ", orig_index
          CALL test_abort("error in xt_idxlist_get_pos_exts_of_index_stripes", &
               __FILE__, &
               __LINE__)
        END IF
        k = k + 1
      END DO
    END DO
    CALL xt_idxlist_delete(intersection)
  END SUBROUTINE check_idxlist_stripes_pos_ext

  SUBROUTINE test_get_pos(stripes, pos)
    TYPE(xt_stripe), INTENT(in) :: stripes(:)
    INTEGER, INTENT(in) :: pos(:)
    INTEGER(xt_int_kind), PARAMETER :: dummy = 1_xi
    INTEGER(xt_int_kind) :: ref_sel_idx(SIZE(pos)), sel_idx(SIZE(pos))
    INTEGER(xt_int_kind), PARAMETER :: undef_idx = -HUGE(dummy)
    INTEGER :: num_pos, ip, p, ref_undef_count, undef_count
    TYPE(xt_idxlist) :: idxlist
    idxlist = xt_idxstripes_new(stripes)
    num_pos = SIZE(pos)
    ref_undef_count = 0
    DO ip = 1, num_pos
      p = pos(ip)
      IF (xt_idxlist_get_index_at_position(idxlist, p, ref_sel_idx(ip))) THEN
        ref_sel_idx(ip) = undef_idx
        ref_undef_count = ref_undef_count + 1
      END IF
    END DO
    undef_count = xt_idxlist_get_indices_at_positions(idxlist, pos, sel_idx, &
         undef_idx)
    IF (undef_count /= ref_undef_count) &
         CALL test_abort("inequal undef count!", &
         __FILE__, &
         __LINE__)
    IF (ANY(sel_idx /= ref_sel_idx)) &
         CALL test_abort("incorrect index returned for position!", &
         __FILE__, &
         __LINE__)
    CALL xt_idxlist_delete(idxlist)
  END SUBROUTINE test_get_pos

  SUBROUTINE test_get_pos1
    TYPE(xt_stripe), PARAMETER :: stripes(3) = (/ xt_stripe(0, 1, 5), &
         xt_stripe(10, 1, 5), xt_stripe(20, -1, 5) /)
    INTEGER, PARAMETER :: pos(13) = &
         (/   0,   2,   7,   9,  11, &
         &  100,  11, 200,   9, 300, &
         &   18, 400,   5 /)
    CALL test_get_pos(stripes, pos)
  END SUBROUTINE test_get_pos1

  SUBROUTINE test_get_pos2
    TYPE(xt_stripe), PARAMETER :: stripes(4) = (/ xt_stripe(0, 1, 3), &
         xt_stripe(10, 1, 2), xt_stripe(20, -1, 6), xt_stripe(30, -1, 7) /)
    INTEGER, PARAMETER :: pos(19) = &
         (/   -1,    0,    1,    2,    3,    4,   23,    5,    6,    7, &
         &     8,    9,   10,   11,   12,    0,    2,  100, 2000 /)
    CALL test_get_pos(stripes, pos)
  END SUBROUTINE test_get_pos2

  SUBROUTINE test_get_pos3
    TYPE(xt_stripe), PARAMETER :: stripes(4) = (/ xt_stripe(0, 1, 3), &
         xt_stripe(10, 1, 2), xt_stripe(20, -1, 6), xt_stripe(30, -1, 7) /)
    INTEGER, PARAMETER :: pos(13) = &
         (/    4,    7,    2,    5,    9,    0,   10,    6,   11,    8, &
         &    12,    1,    3 /)
    CALL test_get_pos(stripes, pos)
  END SUBROUTINE test_get_pos3

  SUBROUTINE test_get_pos4
    TYPE(xt_stripe), PARAMETER :: stripes(3) = (/ xt_stripe(0, 1, 5), &
         xt_stripe(10, 1, 5), xt_stripe(20, -1, 5) /)
    INTEGER, PARAMETER :: pos(7) = &
         (/  -10,  200,  700,   90,   90,   18,  141 /)
    CALL test_get_pos(stripes, pos)
  END SUBROUTINE test_get_pos4

  SUBROUTINE test_stripe_overlap
    TYPE(xt_stripe), PARAMETER :: stripes(2) = (/ xt_stripe(0, 1, 5), &
         xt_stripe(1, 1, 5) /)
    INTEGER(xt_int_kind) :: i, j
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(10) &
         = (/ ((i + j, i=0,4), j = 0, 1) /)
    CALL stripe_test_general(stripes, ref_indices)
  END SUBROUTINE test_stripe_overlap

  SUBROUTINE test_stripe_bb(stripes, global_size, global_start_index, bounds_ref)
    TYPE(xt_stripe), INTENT(in) :: stripes(:)
    INTEGER(xt_int_kind), INTENT(in) :: global_size(:), global_start_index
    TYPE(xt_bounds), INTENT(in) :: bounds_ref(:)

    TYPE(xt_bounds) :: bounds(SIZE(global_size))
    TYPE(xt_idxlist) :: idxstripes

    IF (SIZE(global_size) /= SIZE(bounds_ref)) &
         CALL test_abort("size mismatch for bounding-box", &
         __FILE__, &
         __LINE__)
    idxstripes = xt_idxstripes_new(stripes, SIZE(stripes))

    bounds = xt_idxlist_get_bounding_box(idxstripes, global_size, &
         global_start_index)
    IF (ANY(bounds /= bounds_ref)) &
         CALL test_abort("boundary box doesn't match reference", &
         __FILE__, &
         __LINE__)
    CALL xt_idxlist_delete(idxstripes)
  END SUBROUTINE test_stripe_bb

  SUBROUTINE test_stripe_bb1
    TYPE(xt_stripe), PARAMETER :: stripes(1) = (/ xt_stripe(-1, -1, -1) /)
    INTEGER(xt_int_kind), PARAMETER :: global_size(3) = 4_xi, &
         global_start_index = 0
    TYPE(xt_bounds), PARAMETER :: bounds_ref(3) = xt_bounds(0, 0)
    CALL test_stripe_bb(stripes(1:0), global_size, global_start_index, bounds_ref)
  END SUBROUTINE test_stripe_bb1

  SUBROUTINE test_stripe_bb2
    TYPE(xt_stripe), PARAMETER :: stripes(3) = (/ xt_stripe(47, -12, 2), &
         xt_stripe(32, 12, 2), xt_stripe(36, 12, 2) /)
    INTEGER(xt_int_kind), PARAMETER :: global_size(3) = (/ 5_xi, 4_xi, 3_xi /), &
         global_start_index = 1
    TYPE(xt_bounds), PARAMETER :: bounds_ref(3) = (/ xt_bounds(2, 2), &
         xt_bounds(2, 2), xt_bounds(1, 2) /)
    CALL test_stripe_bb(stripes, global_size, global_start_index, bounds_ref)
  END SUBROUTINE test_stripe_bb2

  SUBROUTINE do_tests(idxlist, ref_indices)
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    INTEGER(xt_int_kind), INTENT(in) :: ref_indices(:)

    TYPE(xt_stripe), ALLOCATABLE :: stripes(:)
    TYPE(xt_stripe), PARAMETER :: dummy(1) = (/ xt_stripe(0,0,0) /)
    INTEGER :: num_stripes
    TYPE(xt_idxlist) :: temp_idxlist, idxlist_copy

    CALL check_idxlist(idxlist, ref_indices)
    CALL xt_idxlist_get_index_stripes(idxlist, stripes)
    IF (ALLOCATED(stripes)) THEN
      num_stripes = SIZE(stripes)
      temp_idxlist = xt_idxvec_from_stripes_new(stripes, num_stripes)
    ELSE
      num_stripes = 0
      temp_idxlist = xt_idxvec_from_stripes_new(dummy, num_stripes)
    END IF
    CALL check_idxlist(temp_idxlist, ref_indices)

    CALL xt_idxlist_delete(temp_idxlist)

    IF (ALLOCATED(stripes)) DEALLOCATE(stripes)

    ! test packing and unpacking
    idxlist_copy = idxlist_pack_unpack_copy(idxlist)

    ! check copy
    CALL check_idxlist(idxlist_copy, ref_indices)

    CALL xt_idxlist_delete(idxlist_copy)

    ! test copying
    idxlist_copy = xt_idxlist_copy(idxlist)

    ! check copy
    CALL check_idxlist(idxlist_copy, ref_indices)

    ! clean up
    CALL xt_idxlist_delete(idxlist_copy)
  END SUBROUTINE do_tests

  SUBROUTINE check_pos_ext(stripes, search_stripes, ref_pos_ext, &
       single_match_only, ref_unmatched, test_desc)
    TYPE(xt_stripe), INTENT(in) :: stripes(:), search_stripes(:)
    TYPE(xt_pos_ext), intent(in) :: ref_pos_ext(:)
    LOGICAL, INTENT(in) :: single_match_only
    INTEGER, INTENT(in) :: ref_unmatched
    CHARACTER(len=*) :: test_desc

    INTEGER :: num_stripes, num_search_stripes, num_ref_pos_ext, num_ext, &
         unmatched
    TYPE(xt_idxlist) :: idxstripes
    TYPE(xt_pos_ext), ALLOCATABLE :: pos_ext(:)

    num_stripes = SIZE(stripes)
    num_search_stripes = SIZE(search_stripes)
    num_ref_pos_ext = SIZE(ref_pos_ext)

    idxstripes = xt_idxstripes_new(stripes)
    unmatched = xt_idxlist_get_pos_exts_of_index_stripes(idxstripes, &
         num_search_stripes, search_stripes, &
         num_ext, pos_ext, single_match_only)
    IF (unmatched /= ref_unmatched) &
         CALL test_abort("error in number of unmatched indices for " &
         // test_desc, &
         __FILE__, &
         __LINE__)
    IF (num_ext < 0 .OR. num_ext /= num_ref_pos_ext) &
         CALL test_abort("error finding " // test_desc, &
         __FILE__, &
         __LINE__)
    IF (ANY(pos_ext /= ref_pos_ext)) &
         CALL test_abort("incorrect position extent length found in "&
         // test_desc, &
         __FILE__, &
         __LINE__)
    DEALLOCATE(pos_ext)
    CALL xt_idxlist_delete(idxstripes)
  END SUBROUTINE check_pos_ext

  SUBROUTINE check_pos_ext1
    INTEGER, PARAMETER :: num_stripes = 1, num_ref_pos_ext = 1, &
         num_ref_unmatched = 0

    TYPE(Xt_stripe), PARAMETER :: stripes(num_stripes) &
         = (/ xt_stripe(1_xi, 1_xi, 10) /), &
         search_stripes(1) = (/ xt_stripe(10_xi, -1_xi, 5) /)

    TYPE(xt_pos_ext), PARAMETER :: ref_pos_ext(num_ref_pos_ext) &
         = (/ xt_pos_ext(9, -5) /)

    CALL check_pos_ext(stripes, search_stripes, ref_pos_ext, .TRUE., &
         num_ref_unmatched, "simple inverted stripe")
  END SUBROUTINE check_pos_ext1

  SUBROUTINE check_pos_ext2
    INTEGER, PARAMETER :: num_stripes = 1, num_ref_pos_ext = 1, &
         num_ref_unmatched = 5

    TYPE(Xt_stripe), PARAMETER :: stripes(num_stripes) &
         = (/ xt_stripe(1_xi, 1_xi, 10) /), &
         search_stripes(2) = xt_stripe(10_xi, -1_xi, 5)

    TYPE(xt_pos_ext), PARAMETER :: ref_pos_ext(num_ref_pos_ext) &
         = (/ xt_pos_ext(9, -5) /)

    CALL check_pos_ext(stripes, search_stripes, ref_pos_ext, .TRUE., &
         num_ref_unmatched, "simple inverted stripe")
  END SUBROUTINE check_pos_ext2

  SUBROUTINE check_pos_ext3
    INTEGER, PARAMETER :: num_stripes = 2, num_ref_pos_ext = 1, &
         num_ref_unmatched = 4

    TYPE(Xt_stripe), PARAMETER :: stripes(num_stripes) &
         = (/ xt_stripe(1_xi, 1_xi, 10), xt_stripe(15_xi, 1_xi, 10) /), &
         search_stripes(1) = xt_stripe(10_xi, 1_xi, 6)

    TYPE(xt_pos_ext), PARAMETER :: ref_pos_ext(num_ref_pos_ext) &
         = (/ xt_pos_ext(9, 2) /)

    CALL check_pos_ext(stripes, search_stripes, ref_pos_ext, .TRUE., &
         num_ref_unmatched, "search inc stripe over inc gap")
  END SUBROUTINE check_pos_ext3

  SUBROUTINE check_pos_ext4
    INTEGER, PARAMETER :: num_stripes = 2, num_ref_pos_ext = 1, &
         num_ref_unmatched = 4

    TYPE(Xt_stripe), PARAMETER :: stripes(num_stripes) &
         = (/ xt_stripe(25_xi, -1_xi, 11), xt_stripe(10_xi, -1_xi, 10) /), &
         search_stripes(1) = xt_stripe(10_xi, 1_xi, 6)

    TYPE(xt_pos_ext), PARAMETER :: ref_pos_ext(num_ref_pos_ext) &
         = (/ xt_pos_ext(11, -2) /)

    CALL check_pos_ext(stripes, search_stripes, ref_pos_ext, .TRUE., &
         num_ref_unmatched, "search inc stripe over dec gap")
  END SUBROUTINE check_pos_ext4

  SUBROUTINE check_pos_ext5
    INTEGER, PARAMETER :: num_stripes = 2, num_ref_pos_ext = 1, &
         num_ref_unmatched = 4

    TYPE(Xt_stripe), PARAMETER :: stripes(num_stripes) &
         = (/ xt_stripe(25_xi, -1_xi, 11), xt_stripe(10_xi, -1_xi, 10) /), &
         search_stripes(1) = xt_stripe(15_xi, -1_xi, 6)

    TYPE(xt_pos_ext), PARAMETER :: ref_pos_ext(num_ref_pos_ext) &
         = (/ xt_pos_ext(10, 2) /)

    CALL check_pos_ext(stripes, search_stripes, ref_pos_ext, .TRUE., &
         num_ref_unmatched, "search dec stripe over dec gap")
  END SUBROUTINE check_pos_ext5

  SUBROUTINE check_pos_ext6
    INTEGER, PARAMETER :: num_stripes = 2, num_ref_pos_ext = 1, &
         num_ref_unmatched = 4

    TYPE(Xt_stripe), PARAMETER :: stripes(num_stripes) &
         = (/ xt_stripe(1_xi, 1_xi, 10), xt_stripe(15_xi, 1_xi, 10) /), &
         search_stripes(1) = xt_stripe(15_xi, -1_xi, 6)

    TYPE(xt_pos_ext), PARAMETER :: ref_pos_ext(num_ref_pos_ext) &
         = (/ xt_pos_ext(10, -2) /)

    CALL check_pos_ext(stripes, search_stripes, ref_pos_ext, .TRUE., &
         num_ref_unmatched, "search dec stripe over inc gap")
  END SUBROUTINE check_pos_ext6

  SUBROUTINE check_pos_ext7
    INTEGER, PARAMETER :: num_stripes = 3, num_ref_pos_ext = 1, &
         num_ref_unmatched = 8

    TYPE(Xt_stripe), PARAMETER :: stripes(num_stripes) &
         = (/ xt_stripe(1_xi, 1_xi, 10), xt_stripe(15_xi, 1_xi, 10), &
         &    xt_stripe(29_xi, 1_xi, 10) /), &
         search_stripes(1) = xt_stripe(32_xi, -1_xi, 30)

    TYPE(xt_pos_ext), PARAMETER :: ref_pos_ext(num_ref_pos_ext) &
         = (/ xt_pos_ext(23, -22) /)

    CALL check_pos_ext(stripes, search_stripes, ref_pos_ext, .TRUE., &
         num_ref_unmatched, "search dec stripe over 2 inc gap")
  END SUBROUTINE check_pos_ext7

  SUBROUTINE check_pos_ext8
    INTEGER, PARAMETER :: num_stripes = 5, num_ref_pos_ext = 5, &
         num_ref_unmatched = 0

    TYPE(Xt_stripe), PARAMETER :: stripes(num_stripes) &
         = (/ xt_stripe(1_xi, 1_xi, 10), xt_stripe(15_xi, 1_xi, 10), &
         &    xt_stripe(29_xi, 1_xi, 10), xt_stripe(14_xi, -1_xi, 4), &
         &    xt_stripe(28_xi, -1_xi, 4) /), &
         search_stripes(1) = xt_stripe(32_xi, -1_xi, 30)

    TYPE(xt_pos_ext), PARAMETER :: ref_pos_ext(num_ref_pos_ext) &
         = (/ xt_pos_ext(23, -4), xt_pos_ext(34, 4), xt_pos_ext(19, -10), &
         &    xt_pos_ext(30, 4), xt_pos_ext(9, -8) /)

    CALL check_pos_ext(stripes, search_stripes, ref_pos_ext, .TRUE., &
         num_ref_unmatched, "search dec stripe over jumbled stripes")
  END SUBROUTINE check_pos_ext8

END PROGRAM test_idxstripes_f
