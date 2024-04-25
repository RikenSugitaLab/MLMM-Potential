!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_energy_go_mod
!> @brief   calculate contact and non-contact energy (Go potential)
!! @authors Takaharu Mori (TM), Chigusa Kobayashi (CK)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
! 
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_energy_go_mod

  use at_pairlist_str_mod
  use at_enefunc_str_mod
  use molecules_str_mod
  use timers_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public :: compute_energy_contact_126
  public :: compute_energy_contact_1210
  public :: compute_energy_contact_12106
  public :: compute_energy_noncontact_nobc

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_contact_126
  !> @brief        calculate contact energy for all-atom Go model
  !! @authors      TM, CK
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] econt    : contact energy of target systems
  !! @note         P.C.Whitford et al., Proteins, 75, 430-441 (2009)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_contact_126(enefunc, coord, force, virial, econt)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: econt

    ! local variables
    real(wp)                 :: dij(3)
    real(wp)                 :: rij2, inv_rij2, inv_rij6, inv_rij12
    real(wp)                 :: lj6, lj12, coef
    real(wp)                 :: work(3)
    integer                  :: i, l, istart, iend
    integer                  :: a1, a2

    real(wp),        pointer :: contact_lj6(:), contact_lj12(:)
    integer,         pointer :: list(:,:)


    call timer(TimerNonBond, TimerOn)
    ! use pointers
    !
    istart       =  enefunc%istart_contact
    iend         =  enefunc%iend_contact
    list         => enefunc%contact_list
    contact_lj6  => enefunc%contact_lj6
    contact_lj12 => enefunc%contact_lj12

    ! calculate energy and gradient
    !
    !$omp parallel do default(none)                                       &
    !$omp private(i, l, dij, rij2, inv_rij2, inv_rij6, inv_rij12, coef,   &
    !$omp         work, lj6, lj12,a1, a2)                                 &
    !$omp shared(istart, iend, list, coord, contact_lj6, contact_lj12) &
    !$omp reduction(+:force) reduction(+:virial) reduction(+:econt)
    !
    do i = istart, iend

      ! compute distance
      !
      a1 = list(1,i)
      a2 = list(2,i)
      dij(1:3) = coord(1:3,a1) - coord(1:3,a2)
      rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

      lj6  = contact_lj6(i)
      lj12 = contact_lj12(i)

      ! lj energy
      !
      inv_rij2  = 1.0_wp / rij2
      inv_rij6  = inv_rij2 * inv_rij2 * inv_rij2
      inv_rij12 = inv_rij6 * inv_rij6
      econt = econt + (lj12 * inv_rij12 - lj6 * inv_rij6)
 
      ! gradient
      !
      coef = - inv_rij2 * (12.0_wp*lj12*inv_rij12-6.0_wp*lj6*inv_rij6)
      work(1:3) = coef*dij(1:3)

      ! store force
      !
      force(1:3,list(1,i)) = force(1:3,list(1,i)) - work(1:3)
      force(1:3,list(2,i)) = force(1:3,list(2,i)) + work(1:3)
        
      ! virial
      !
      do l = 1, 3
        virial(1:3,l) = virial(1:3,l) - dij(1:3)*work(l)
      end do
    end do
    !$omp end parallel do
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_contact_126

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_contact_1210
  !> @brief        calculate contact energy for C alpha Go model
  !! @authors      TM, CK
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] econt    : contact energy of target systems
  !! @note         C.Clementi et al., J.Mol.Biol., 298, 937-953 (2000)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_contact_1210(enefunc, coord, force, virial, econt)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: econt

    ! local variables
    real(wp)                 :: dij(3), rij2
    real(wp)                 :: inv_rij2, inv_rij6, inv_rij10, inv_rij12
    real(wp)                 :: lj10, lj12, coef
    real(wp)                 :: work(3)
    integer                  :: i, l, istart, iend
    integer                  :: a1, a2

    real(wp),        pointer :: contact_lj10(:), contact_lj12(:)
    integer,         pointer :: list(:,:)


    call timer(TimerNonBond, TimerOn)
    ! use pointers
    !
    istart       =  enefunc%istart_contact
    iend         =  enefunc%iend_contact
    list         => enefunc%contact_list
    contact_lj10 => enefunc%contact_lj10
    contact_lj12 => enefunc%contact_lj12

    ! calculate energy and gradient
    !
    !$omp parallel do default(none)                                          &
    !$omp private(i, l, dij, rij2, inv_rij2, inv_rij6, inv_rij10, inv_rij12, &
    !$omp         coef, work, lj10, lj12, a1, a2)                            &
    !$omp shared(istart, iend, list, coord, contact_lj10, contact_lj12)      &
    !$omp reduction(+:force) reduction(+:virial) reduction(+:econt)
    !
    do i = istart, iend

      ! compute distance
      !
      a1 = list(1,i)
      a2 = list(2,i)
      dij(1:3) = coord(1:3,a1) - coord(1:3,a2)
      rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

      lj10 = contact_lj10(i)
      lj12 = contact_lj12(i)

      ! lj energy
      !
      inv_rij2  = 1.0_wp / rij2
      inv_rij6  = inv_rij2 * inv_rij2 * inv_rij2
      inv_rij10 = inv_rij6 * inv_rij2 * inv_rij2
      inv_rij12 = inv_rij6 * inv_rij6
      econt = econt + (lj12 * inv_rij12 - lj10 * inv_rij10)
 
      ! gradient
      !
      coef = - inv_rij2 * (12.0_wp*lj12*inv_rij12-10.0_wp*lj10*inv_rij10)
      work(1:3) = coef*dij(1:3)

      ! store force
      !
      force(1:3,list(1,i)) = force(1:3,list(1,i)) - work(1:3)
      force(1:3,list(2,i)) = force(1:3,list(2,i)) + work(1:3)
        
      ! virial
      !
      do l = 1, 3
        virial(1:3,l) = virial(1:3,l) - dij(1:3)*work(l)
      end do
    end do
    !$omp end parallel do
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_contact_1210

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_contact_12106
  !> @brief        calculate contact energy for KBGO
  !! @authors      TM, CK
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] econt    : contact energy of target systems
  !! @note         Karanicolas & Brooks, Protein Sci., 11, 2351-2361 (2002)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_contact_12106(enefunc, coord, force, virial, econt)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: econt

    ! local variables
    real(wp)                 :: dij(3), rij2
    real(wp)                 :: inv_rij2, inv_rij6, inv_rij10, inv_rij12
    real(wp)                 :: lj6, lj10, lj12, coef
    real(wp)                 :: work(3)
    integer                  :: i, l, istart, iend
    integer                  :: a1, a2

    real(wp),        pointer :: contact_lj10(:), contact_lj12(:)
    real(wp),        pointer :: contact_lj6(:)
    integer,         pointer :: list(:,:)


    call timer(TimerNonBond, TimerOn)
    ! use pointers
    !
    istart       =  enefunc%istart_contact
    iend         =  enefunc%iend_contact
    list         => enefunc%contact_list
    contact_lj6  => enefunc%contact_lj6
    contact_lj10 => enefunc%contact_lj10
    contact_lj12 => enefunc%contact_lj12

    ! calculate energy and gradient
    !
    !$omp parallel do default(none)                                          &
    !$omp private(i, l, dij, rij2, inv_rij2, inv_rij6, inv_rij10, inv_rij12, &
    !$omp         coef, work, lj6, lj10, lj12, a1, a2)                       &
    !$omp shared(istart, iend, list, coord, contact_lj6,                     &
    !$omp       contact_lj10, contact_lj12)                                  &
    !$omp reduction(+:force) reduction(+:virial) reduction(+:econt)
    !
    do i = istart, iend

      ! compute distance
      !
      a1 = list(1,i)
      a2 = list(2,i)
      dij(1:3) = coord(1:3,a1) - coord(1:3,a2)
      rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

      lj6  = contact_lj6(i)
      lj10 = contact_lj10(i)
      lj12 = contact_lj12(i)

      ! lj energy
      !
      inv_rij2  = 1.0_wp / rij2
      inv_rij6  = inv_rij2 * inv_rij2 * inv_rij2
      inv_rij10 = inv_rij6 * inv_rij2 * inv_rij2
      inv_rij12 = inv_rij6 * inv_rij6
      econt = econt + (lj12 * inv_rij12 - lj10 * inv_rij10 + lj6 * inv_rij6)
 
      ! gradient
      !
      coef = - inv_rij2 * (12.0_wp*lj12*inv_rij12-10.0_wp*lj10*inv_rij10 &
                          + 6.0_wp * lj6 * inv_rij6)
      work(1:3) = coef*dij(1:3)

      ! store force
      !
      force(1:3,list(1,i)) = force(1:3,list(1,i)) - work(1:3)
      force(1:3,list(2,i)) = force(1:3,list(2,i)) + work(1:3)
        
      ! virial
      !
      do l = 1, 3
        virial(1:3,l) = virial(1:3,l) - dij(1:3)*work(l)
      end do
    end do
    !$omp end parallel do
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_contact_12106

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_noncontact_nobc
  !> @brief        calculate non-contact energy with pairlist (NOBC)
  !! @authors      TM, CK
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] encont   : contact energy of target systems
  !! @note         C.Clementi et al., J.Mol.Biol., 298, 937-953 (2000)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_noncontact_nobc(enefunc, molecule, pairlist, &
                                            coord, force, virial, encont)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: encont

    ! local variables
    real(wp)                  :: dij(3), rij2
    real(wp)                  :: inv_rij2, inv_rij6, inv_rij12, term_lj12
    real(wp)                  :: lj12, coef, noncontact_lj12
    real(wp)                  :: work(3)
    real(wp)                  :: cutoff2, cutoff
    integer                   :: i, j, k, l, natom, id
    integer                   :: num_nb15, ini_nb15, fin_nb15
    integer                   ::  nthread, my_id
    integer                   :: omp_get_num_threads, omp_get_thread_num
    integer                   :: forcefield
         
    integer, pointer          :: num_nb15_calc(:,:), nb15_calc_list(:,:)
    integer, pointer          :: atmcls(:)
    real(wp), pointer         :: nonb_lj12(:,:)


    call timer(TimerNonBond, TimerOn)
    ! use pointers
    !
    natom          =  molecule%num_atoms
    num_nb15_calc  => pairlist%num_nb15_calc
    nb15_calc_list => pairlist%nb15_calc_list
    num_nb15       = 0
    atmcls         => molecule%atom_cls_no

    noncontact_lj12=  enefunc%noncontact_lj12
    forcefield     =  enefunc%forcefield
    cutoff         =  enefunc%cutoffdist
    cutoff2        =  cutoff * cutoff
    nonb_lj12      => enefunc%nonb_lj12

    ! calculate energy and gradient
    !
    !$omp parallel default(none)                                    &
    !$omp firstprivate(num_nb15)                                    & 
    !$omp private(ini_nb15, fin_nb15, i, k, j, l, dij, coef, work,  &
    !$omp         rij2, inv_rij2, inv_rij6, inv_rij12, term_lj12,   &
    !$omp         my_id, lj12, id)                                  &
    !$omp shared(natom, num_nb15_calc, nb15_calc_list, cutoff2,     &
    !$omp        coord, nthread,  my_city_rank, nproc_city, atmcls, &
    !$omp        forcefield, noncontact_lj12, nonb_lj12)            &
    !$omp reduction(+:force) reduction(+:virial) reduction(+:encont) 
    !
#ifdef OMP
    id      = omp_get_thread_num()
    nthread = omp_get_num_threads()
#else
    id      = 0
    nthread = 1
#endif
    my_id   = my_city_rank * nthread + id
!    if (forcefield /= ForcefieldKBGO) lj12 = noncontact_lj12

    do i = 1, natom-1

      ini_nb15 = num_nb15 + 1
      fin_nb15 = num_nb15 + num_nb15_calc(i,id+1)
      num_nb15 = fin_nb15

      if (mod(i-1,nproc_city*nthread) /= my_id) ini_nb15 = fin_nb15 + 1

      do k = ini_nb15, fin_nb15
        
        j = nb15_calc_list(k,id+1)

        ! compute distance
        !
        dij(1:3) = coord(1:3,i) - coord(1:3,j)
        rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
!        if (forcefield == ForcefieldKBGO)  &
           lj12 = nonb_lj12(atmcls(i),atmcls(j))
        ! cutoff
        !
        if (rij2 >= cutoff2) cycle

        ! non-contact energy
        !
        inv_rij2  = 1.0_wp / rij2
        inv_rij6  = inv_rij2 * inv_rij2 * inv_rij2
        inv_rij12 = inv_rij6 * inv_rij6
        term_lj12 = lj12 * inv_rij12
        
        encont = encont + term_lj12
        
        ! gradient
        !
        coef = - inv_rij2 * (12.0_wp*term_lj12)
        work(1:3) = coef*dij(1:3)
        
        ! store force
        !
        force(1:3,i) = force(1:3,i) - work(1:3)
        force(1:3,j) = force(1:3,j) + work(1:3)
        
        ! virial
        !
        do l = 1, 3
          virial(1:3,l) = virial(1:3,l) - dij(1:3)*work(l)
        end do

      end do

    end do
    !$omp end parallel
    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_noncontact_nobc

end module at_energy_go_mod
