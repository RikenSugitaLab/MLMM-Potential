!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_energy_angles_mod
!> @brief   calculate angle energy
!! @authors Yuji Sugita (YS), Takaharu Mori (TM), Jaewoon Jung (JJ)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_energy_angles_mod

  use at_enefunc_str_mod
  use timers_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public :: compute_energy_angle
  public :: compute_energy_angle_g96
  public :: compute_energy_urey

contains
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_angle
  !> @brief        calculate angle energy
  !! @authors      YS, TM
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[inout] eangle  : ene_angle angle energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_angle(enefunc, coord, force, virial, eangle)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: eangle

    ! local variables
    real(wp)                 :: d12(1:3), r12_2, inv_r12_2
    real(wp)                 :: d32(1:3), r32_2, inv_r32_2
    real(wp)                 :: r12r32, inv_r12r32
    real(wp)                 :: cos_t, sin_t, t123, t_dif, vtmp
    real(wp)                 :: cc_frc, cc_frc2, cc_frc3
    integer                  :: i, j, k, istart, iend

    real(wp),        pointer :: fc(:), theta0(:)
    real(wp),        pointer :: work(:,:)
    integer,         pointer :: list(:,:)


    call timer(TimerAngle, TimerOn)

    istart = enefunc%istart_angle
    iend   = enefunc%iend_angle

    list   => enefunc%angl_list
    fc     => enefunc%angl_force_const
    theta0 => enefunc%angl_theta_min
    work   => enefunc%work


    ! calculation of angle energy and gradient
    !
    !$omp parallel do default(none)                                    &
    !$omp private(i, j, k, d12, d32, r12_2, r32_2, r12r32, inv_r12r32, &
    !$omp         inv_r12_2, inv_r32_2, cos_t, sin_t, t123, t_dif,     &
    !$omp         cc_frc, cc_frc2, cc_frc3, vtmp)                      &
    !$omp shared(istart, iend, fc, theta0, work, list, coord)          &
    !$omp reduction(+:eangle) reduction(+:virial)
    !
    do i = istart, iend

      ! angle energy: E=K[t-t0]^2
      !
      d12(1:3) = coord(1:3,list(1,i)) - coord(1:3,list(2,i))
      d32(1:3) = coord(1:3,list(3,i)) - coord(1:3,list(2,i))
      r12_2    = d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3)
      r32_2    = d32(1)*d32(1) + d32(2)*d32(2) + d32(3)*d32(3)
      r12r32   = sqrt( r12_2*r32_2 )
      inv_r12r32 = 1.0_wp / r12r32
      inv_r12_2  = 1.0_wp / r12_2
      inv_r32_2  = 1.0_wp / r32_2

      cos_t  = ( d12(1)*d32(1) + d12(2)*d32(2) + d12(3)*d32(3) ) * inv_r12r32
      cos_t  = min(  1.0_wp, cos_t )
      cos_t  = max( -1.0_wp, cos_t )
      t123   = acos(cos_t)
      t_dif  = t123 - theta0(i)
      eangle = eangle + fc(i) * t_dif * t_dif

      ! gradient: dE/dX
      !
      sin_t  = sin(t123)
      sin_t  = max( EPS, sin_t )
      cc_frc = - ( 2.0_wp * fc(i) * t_dif ) / sin_t
      cc_frc2 = cos_t * inv_r12_2
      cc_frc3 = cos_t * inv_r32_2
      work(1:3,i) = cc_frc * ( d32(1:3) * inv_r12r32 - d12(1:3) * cc_frc2 )
      work(4:6,i) = cc_frc * ( d12(1:3) * inv_r12r32 - d32(1:3) * cc_frc3 )

      ! virial from angle
      !
      do j = 1, 3
        do k = j+1, 3
          vtmp = d12(k)*work(j,i) + d32(k)*work(j+3,i)
          virial(k,j) = virial(k,j) - vtmp
          virial(j,k) = virial(j,k) - vtmp
        end do
        vtmp = d12(j)*work(j,i) + d32(j)*work(j+3,i)
        virial(j,j) = virial(j,j) - vtmp
      end do

    end do
    !$omp end parallel do

    ! store force
    ! 
    do i = istart, iend
      force(1:3,list(1,i)) = force(1:3,list(1,i)) - work(1:3,i)
      force(1:3,list(2,i)) = force(1:3,list(2,i)) + work(1:3,i) + work(4:6,i)
      force(1:3,list(3,i)) = force(1:3,list(3,i)) - work(4:6,i)
    end do

    call timer(TimerAngle, TimerOff)
 
    return
  
  end subroutine compute_energy_angle

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_angle_g96
  !> @brief        calculate angle energy
  !! @authors      CK
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[inout] eangle  : ene_angle angle energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_angle_g96(enefunc, coord, force, virial, eangle)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: eangle

    ! local variables
    real(wp)                 :: d12(1:3), r12_2, inv_r12_2
    real(wp)                 :: d32(1:3), r32_2, inv_r32_2
    real(wp)                 :: r12r32, inv_r12r32
    real(wp)                 :: cos_t, dif, cc_frc, cc_frc2, cc_frc3, vtmp
    integer                  :: i, j, k, istart, iend

    real(wp),        pointer :: fc(:), theta0(:)
    real(wp),        pointer :: work(:,:)
    integer,         pointer :: list(:,:)


    call timer(TimerAngle, TimerOn)

    istart = enefunc%istart_angle
    iend   = enefunc%iend_angle

    list   => enefunc%angl_list
    fc     => enefunc%angl_force_const
    theta0 => enefunc%angl_theta_min
    work   => enefunc%work


    ! calculation of angle energy and gradient
    !
    !$omp parallel do default(none)                                         &
    !$omp private(i, j, k, d12, d32, r12_2, inv_r12_2, r32_2, inv_r32_2,    &
    !$omp         r12r32, inv_r12r32, cos_t, dif, cc_frc, cc_frc2, cc_frc3, &
    !$omp         vtmp)                                                     &
    !$omp shared(istart, iend, fc, theta0, work, list, coord)               &
    !$omp reduction(+:eangle) reduction(+:virial)
    !
    do i = istart, iend

      ! angle energy: E=K[t-t0]^2
      !
      d12(1:3) = coord(1:3,list(1,i)) - coord(1:3,list(2,i))
      d32(1:3) = coord(1:3,list(3,i)) - coord(1:3,list(2,i))

      r12_2   = d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3)
      r32_2   = d32(1)*d32(1) + d32(2)*d32(2) + d32(3)*d32(3)
      r12r32  = sqrt( r12_2*r32_2 )

      inv_r12r32 = 1.0_wp / r12r32
      inv_r12_2  = 1.0_wp / r12_2
      inv_r32_2  = 1.0_wp / r32_2

      cos_t   = (d12(1)*d32(1) + d12(2)*d32(2) + d12(3)*d32(3)) * inv_r12r32

      dif     = cos_t - cos(theta0(i))
      cc_frc  = fc(i) * dif
      eangle  = eangle + ( cc_frc * dif )

      ! gradient: dE/dX
      !
      cc_frc  = 2.0_wp * cc_frc
      cc_frc2 = cos_t * inv_r12_2
      cc_frc3 = cos_t * inv_r32_2
      work(1:3,i) = cc_frc * ( d32(1:3) * inv_r12r32 - d12(1:3) * cc_frc2 )
      work(4:6,i) = cc_frc * ( d12(1:3) * inv_r12r32 - d32(1:3) * cc_frc3 )

      ! virial from angle
      !
      do j = 1, 3
        do k = j+1, 3
          vtmp = d12(k)*work(j,i) + d32(k)*work(j+3,i)
          virial(k,j) = virial(k,j) - vtmp
          virial(j,k) = virial(j,k) - vtmp
        end do
        vtmp = d12(j)*work(j,i) + d32(j)*work(j+3,i)
        virial(j,j) = virial(j,j) - vtmp
      end do

    end do
    !$omp end parallel do

    ! store force
    !
    do i = istart, iend
      force(1:3,list(1,i)) = force(1:3,list(1,i)) - work(1:3,i)
      force(1:3,list(2,i)) = force(1:3,list(2,i)) + work(1:3,i) + work(4:6,i)
      force(1:3,list(3,i)) = force(1:3,list(3,i)) - work(4:6,i)
    end do

    call timer(TimerAngle, TimerOff)

    return

  end subroutine compute_energy_angle_g96

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_urey
  !> @brief        calculate Urey-Bradley energy
  !! @authors      YS, TM
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[inout] eurey   : urey-b energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_urey(enefunc, coord, force, virial, eurey)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: eurey

    ! local variables
    real(wp)                 :: d13(1:3), r13, ub_dif, cc_frc, vtmp
    integer                  :: i, j, k, istart, iend

    real(wp),        pointer :: fc_ub(:), r0_ub(:)
    real(wp),        pointer :: work(:,:)
    integer,         pointer :: list(:,:)


    call timer(TimerAngle, TimerOn)

    istart = enefunc%istart_urey
    iend   = enefunc%iend_urey

    list   => enefunc%urey_list
    fc_ub  => enefunc%urey_force_const
    r0_ub  => enefunc%urey_rmin
    work   => enefunc%work


    ! calculation of angle energy and gradient
    !
    !$omp parallel do default(none)                             &
    !$omp private(i, j, k, d13, r13, ub_dif, cc_frc, vtmp)      &
    !$omp shared(istart, iend, fc_ub, r0_ub, work, coord, list) &
    !$omp reduction(+:eurey) reduction(+:virial)
    !
    do i = istart, iend

      ! urey-bradley energy: E=K[ub-ub0]^2
      !
      d13(1:3) = coord(1:3,list(1,i)) - coord(1:3,list(2,i))
      r13    = sqrt( d13(1)*d13(1) + d13(2)*d13(2) + d13(3)*d13(3) )
      ub_dif = r13 - r0_ub(i)
      eurey  = eurey + fc_ub(i) * ub_dif * ub_dif

      ! gradient: dE/dx
      !
      cc_frc = ( 2.0_wp * fc_ub(i) * ub_dif ) / r13
      work(1:3,i) = cc_frc * d13(1:3)

      ! virial from urey-bradley
      !
      do j = 1, 3
        do k = j+1, 3
          vtmp = d13(k)*work(j,i) 
          virial(k,j) = virial(k,j) - vtmp
          virial(j,k) = virial(j,k) - vtmp
        end do
        vtmp = d13(j)*work(j,i) 
        virial(j,j) = virial(j,j) - vtmp
      end do

    end do
    !$omp end parallel do


    ! store force
    !
    do i = istart, iend
      force(1:3,list(1,i)) = force(1:3,list(1,i)) - work(1:3,i)
      force(1:3,list(2,i)) = force(1:3,list(2,i)) + work(1:3,i)
    end do


    call timer(TimerAngle, TimerOff)
 
    return
  
  end subroutine compute_energy_urey
  
end module at_energy_angles_mod
