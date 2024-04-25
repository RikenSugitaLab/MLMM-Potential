!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_energy_bonds_mod
!> @brief   calculate bond energy
!! @authors Yuji Sugita (YS), Takaharu Mori (TM), Jaewoon Jung (JJ)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
! 
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_energy_bonds_mod

  use at_enefunc_str_mod
  use timers_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public :: compute_energy_bond

contains
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_bond
  !> @brief        calculate bond energy
  !! @authors      YS, TM, JJ
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[inout] ebond   : bond energy of target systems
  !! @note         Fujitsu suggested that force reduction in OpenMP is too slow.
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_bond(enefunc, coord, force, virial, ebond)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: ebond

    ! local variables
    real(wp)                 :: d12(1:3), r12, r_dif, cc_frc, vtmp
    integer                  :: i, j, k, istart, iend

    real(wp),        pointer :: fc(:), r0(:), work(:,:)
    integer,         pointer :: list(:,:)


    call timer(TimerBond, TimerOn)

    ! use pointers
    !
    istart = enefunc%istart_bond
    iend   = enefunc%iend_bond
    list   => enefunc%bond_list
    fc     => enefunc%bond_force_const
    r0     => enefunc%bond_dist_min
    work   => enefunc%work


    ! calculation of bond energy and gradient
    !
    !$omp parallel do default(none)                          &
    !$omp private(i, j, k, d12, r12, r_dif, cc_frc, vtmp)    &
    !$omp shared(istart, iend, fc, r0, work, coord, list)    &
    !$omp reduction(+:ebond) reduction(+:virial) 
    !
    do i = istart, iend
    
      ! bond energy: E=K[b-b0]^2
      !
      d12(1:3) = coord(1:3,list(1,i)) - coord(1:3,list(2,i))
      r12   = sqrt( d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3) )
      r_dif = r12 - r0(i)
      ebond = ebond + fc(i) * r_dif * r_dif

      ! gradient: dE/dX
      !
      cc_frc    = (2.0_wp * fc(i) * r_dif) / r12
      work(1:3,i) = cc_frc * d12(1:3)

      ! virial coefficient
      !
      do j = 1, 3
        do k = j+1, 3
          vtmp = d12(k) * work(j, i) 
          virial(k,j) = virial(k,j) - vtmp
          virial(j,k) = virial(j,k) - vtmp
        end do
        vtmp = d12(j) * work(j, i) 
        virial(j,j) = virial(j,j) - vtmp
      end do

    end do
    !$omp end parallel do

    ! store force: F=-dE/dX
    !
    do i = istart, iend
      force(1:3,list(1,i)) = force(1:3,list(1,i)) - work(1:3,i)
      force(1:3,list(2,i)) = force(1:3,list(2,i)) + work(1:3,i)
    end do

    call timer(TimerBond, TimerOff)

    return

  end subroutine compute_energy_bond

end module at_energy_bonds_mod
