!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_energy_table_linear_bondcorr_mod
!> @brief   calculate bond correction with linear interpolation table
!! @authors Jaewoon Jung(JJ)
!  
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_energy_table_linear_bondcorr_mod

  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  !
  public   :: pme_bond_corr_linear
  private  :: pme_bond_corr_linear_general
  private  :: pme_bond_corr_linear_gro_amber
  private  :: pme_bond_corr_linear_general_check
  private  :: pme_bond_corr_linear_gro_amber_check
  ! FEP
  public   :: pme_bond_corr_linear_fep
  private  :: pme_bond_corr_linear_general_fep
  private  :: pme_bond_corr_linear_general_check_fep
  private  :: pme_bond_corr_linear_gro_amber_fep
  private  :: pme_bond_corr_linear_gro_amber_check_fep

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_linear
  !> @brief        calculate pme bond correction with linear lookup table
  !  @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_linear(domain, enefunc, force, eelec)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)

    if (enefunc%nonb_limiter) then
      call pme_bond_corr_linear_general_check(domain, enefunc, force, eelec)
     
      if (enefunc%forcefield /= ForcefieldCHARMM) then
     
        call pme_bond_corr_linear_gro_amber_check(domain, enefunc, force, eelec)
     
      end if
    else
    
      call pme_bond_corr_linear_general(domain, enefunc, force, eelec)
     
      if (enefunc%forcefield /= ForcefieldCHARMM) then
     
        call pme_bond_corr_linear_gro_amber(domain, enefunc, force, eelec)
     
      end if
    endif

    return

   end subroutine pme_bond_corr_linear

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_linear_general
  !> @brief        calculate bond correction term in PME (general)
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[inout] force   : forces of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_linear_general(domain, enefunc, force, eelec)

    ! formal arguments 
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: dij(3), rij2, work(3), cutoff2
    real(wp)                 :: R, term_elec, ccr, coef
    integer                  :: i, ix, iy, j, k, ij, L, iwater
    integer                  :: num_excl, ini_excl, fin_excl
    integer                  :: id, omp_get_thread_num
    integer                  :: ncell_local

    real(dp),        pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: table_ecor(:), table_decor(:)
    real(wp),        pointer :: density, cutoff
    integer,         pointer :: cell_pair(:,:)
    integer,         pointer :: natom(:), nwater(:), water_list(:,:,:)
    integer,         pointer :: num_nonb_excl1(:,:), nonb_excl_list1(:,:)
    integer,         pointer :: num_nonb_excl(:,:), nonb_excl_list(:,:)


    cell_pair       => domain%cell_pairlist1
    natom           => domain%num_atom
    nwater          => domain%num_water
    water_list      => domain%water_list
    coord           => domain%coord
    charge          => domain%charge

    table_ecor      => enefunc%table%table_ecor
    table_decor     => enefunc%table%table_decor
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    num_nonb_excl1  => enefunc%num_nonb_excl1
    num_nonb_excl   => enefunc%num_nonb_excl
    nonb_excl_list1 => enefunc%nonb_excl_list1
    nonb_excl_list  => enefunc%nonb_excl_list

    ncell_local = domain%num_cell_local
    cutoff2     = cutoff*cutoff

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, j, k, ij, ix, iy, ini_excl, fin_excl, num_excl,       &
    !$omp         dij, rij2, L, R, coef, work, term_elec,  ccr)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    ! calculate energy and gradient
    !
    do i = id+1, ncell_local, nthread

      num_excl = 0

      do ix = 1, natom(i) - 1

        ini_excl = num_excl + 1
        fin_excl = num_excl + num_nonb_excl1(ix,i)
        num_excl = fin_excl

        ! skip for zero charge
        !
        if (abs(charge(ix,i)) < EPS) cycle

        do k = ini_excl, fin_excl

          iy = nonb_excl_list1(k,i)

          ! skip for zero charge
          !
          if (abs(charge(iy,i)) < EPS) cycle

          ! compute distance
          !
          dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,i)
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          ! energy and gradient
          !
          rij2 = cutoff2*density/rij2
          L    = int(rij2)
          R    = rij2 - L

          term_elec = table_ecor(L) + R*(table_ecor(L+1)-table_ecor(L))
          ccr       = charge(ix,i) * charge(iy,i) * term_elec
          eelec(id+1) = eelec(id+1) + ccr

          term_elec = table_decor(L) + R*(table_decor(L+1)-table_decor(L))
          coef      = charge(ix,i) * charge(iy,i) * term_elec

          work(1:3) = coef*dij(1:3)

          force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
          force(1:3,iy,i,id+1) = force(1:3,iy,i,id+1) + work(1:3)

        end do
  
      end do

    end do

    do ij = id+1, maxcell_near, nthread

      i = cell_pair(1,ij)
      j = cell_pair(2,ij)

      num_excl = 0

      do ix = 1, natom(i)
        ini_excl = num_excl + 1
        fin_excl = num_excl + num_nonb_excl(ix,ij)
        num_excl = fin_excl

        ! skip for zero charge
        !
        if (abs(charge(ix,i)) < EPS) cycle

        do k = ini_excl, fin_excl

          iy = nonb_excl_list(k,ij)

          ! skip for zero charge
          !
          if (abs(charge(iy,j)) < EPS) cycle

          ! compute distance
          !
          dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,j)
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          ! energy and gradient
          !
          rij2 = cutoff2*density/rij2
          L    = int(rij2)
          R    = rij2 - L

          term_elec = table_ecor(L) + R*(table_ecor(L+1)-table_ecor(L))

          ccr       = charge(ix,i) * charge(iy,j) * term_elec
          eelec(id+1) = eelec(id+1) + ccr

          term_elec = table_decor(L) + R*(table_decor(L+1)-table_decor(L))
          coef      = charge(ix,i) * charge(iy,j) * term_elec

          work(1:3) = coef*dij(1:3)
          force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
          force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + work(1:3)

        end do

      end do
    end do

    !$omp end parallel

    return

  end subroutine pme_bond_corr_linear_general

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_linear_gro_amber
  !> @brief        calculate bodn correction relating with 14 scaling
  !  @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_linear_gro_amber(domain, enefunc, &
                                          force, eelec)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: dij(3), rij2, work(3), cutoff2
    real(wp)                 :: R, term_elec, coef, cc, qq_scale
    integer                  :: i, ix, iy, j, k, ij, L
    integer                  :: num_nb14, ini_nb14, fin_nb14
    integer                  :: ncell_local, id, omp_get_thread_num

    real(dp),        pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: table_ene(:), table_grad(:)
    real(wp),        pointer :: density, cutoff
    integer,         pointer :: cell_pair(:,:), atmcls(:,:), natom(:)
    integer,         pointer :: num_nb14_calc1(:,:), nb14_calc_list1(:,:)
    integer,         pointer :: num_nb14_calc(:,:), nb14_calc_list(:,:)


    coord           => domain%coord
    charge          => domain%charge
    cell_pair       => domain%cell_pairlist1
    natom           => domain%num_atom

    table_ene       => enefunc%table%table_ecor
    table_grad      => enefunc%table%table_decor
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    num_nb14_calc1  => enefunc%num_nb14_calc1
    num_nb14_calc   => enefunc%num_nb14_calc
    nb14_calc_list1 => enefunc%nb14_calc_list1
    nb14_calc_list  => enefunc%nb14_calc_list

    cutoff2         =  cutoff * cutoff
    ncell_local     =  domain%num_cell_local

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, j, k, ij, ix, iy, ini_nb14, fin_nb14, num_nb14,       &
    !$omp         dij, rij2, L, R, coef, work, term_elec, cc, qq_scale)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    ! calculate energy and gradient
    !
    do i = id+1, ncell_local, nthread
      num_nb14 = 0

      do ix = 1, natom(i) - 1

        ini_nb14 = num_nb14 + 1
        fin_nb14 = num_nb14 + num_nb14_calc1(ix,i)
        num_nb14 = fin_nb14

        ! skip for zero charge
        !
        if (abs(charge(ix,i)) < EPS) cycle

        do k = ini_nb14, fin_nb14

          iy = nb14_calc_list1(k,i)

          ! skip for zero charge
          !
          if (abs(charge(iy,i)) < EPS) cycle

          ! compute distance
          !
          dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,i)
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          ! energy and gradient
          !
          rij2 = cutoff2*density/rij2
          L    = int(rij2)
          R    = rij2 - L

          qq_scale  = -enefunc%nb14_qq_scale1(k,i)+1.0_wp
          term_elec = table_ene(L) + R*(table_ene(L+1)-table_ene(L))
          cc        = charge(ix,i)*charge(iy,i)*qq_scale
          eelec(id+1) = eelec(id+1) + term_elec*cc

          term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
          coef      = cc*term_elec

          work(1:3) = coef*dij(1:3)

          force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
          force(1:3,iy,i,id+1) = force(1:3,iy,i,id+1) + work(1:3)

        end do

      end do
    end do

    do ij = id+1, maxcell_near, nthread

      i = cell_pair(1,ij)
      j = cell_pair(2,ij)

      num_nb14 = 0

      do ix = 1, natom(i)

        ini_nb14 = num_nb14 + 1
        fin_nb14 = num_nb14 + num_nb14_calc(ix,ij)
        num_nb14 = fin_nb14

        ! skip for zero charge
        !
        if (abs(charge(ix,i)) < EPS) cycle

        do k = ini_nb14, fin_nb14

          iy = nb14_calc_list(k,ij)

          ! skip for zero charge
          !
          if (abs(charge(iy,j)) < EPS) cycle

          ! compute distance
          !
          dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,j)
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          ! energy and gradient
          !
          rij2 = cutoff2*density/rij2
          L    = int(rij2)
          R    = rij2 - L

          qq_scale  = -enefunc%nb14_qq_scale(k,ij)+1.0_wp
          cc        = charge(ix,i)*charge(iy,j)*qq_scale
          term_elec = table_ene(L) + R*(table_ene(L+1)-table_ene(L))
          eelec(id+1) = eelec(id+1) + term_elec*cc

          term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
          coef      = cc*term_elec

          work(1:3) = coef*dij(1:3)

          force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
          force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + work(1:3)

        end do

      end do

    end do

    !$omp end parallel

    return

  end subroutine pme_bond_corr_linear_gro_amber

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_linear_general_check
  !> @brief        calculate bond correction term in PME (general)
  !! @authors      JJ, CK
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[inout] force   : forces of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_linear_general_check(domain, enefunc, force, eelec)

    ! formal arguments 
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: dij(3), rij2, work(3), cutoff2
    real(wp)                 :: R, term_elec, ccr, coef
    real(wp)                 :: minimum_contact
    integer                  :: i, ix, iy, j, k, ij, L, iwater
    integer                  :: num_excl, ini_excl, fin_excl
    integer                  :: id, omp_get_thread_num
    integer                  :: ncell_local

    real(dp),        pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: table_ecor(:), table_decor(:)
    real(wp),        pointer :: density, cutoff
    integer,         pointer :: cell_pair(:,:)
    integer,         pointer :: natom(:), nwater(:), water_list(:,:,:)
    integer,         pointer :: num_nonb_excl1(:,:), nonb_excl_list1(:,:)
    integer,         pointer :: num_nonb_excl(:,:), nonb_excl_list(:,:)


    cell_pair       => domain%cell_pairlist1
    natom           => domain%num_atom
    nwater          => domain%num_water
    water_list      => domain%water_list
    coord           => domain%coord
    charge          => domain%charge

    table_ecor      => enefunc%table%table_ecor
    table_decor     => enefunc%table%table_decor
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    num_nonb_excl1  => enefunc%num_nonb_excl1
    num_nonb_excl   => enefunc%num_nonb_excl
    nonb_excl_list1 => enefunc%nonb_excl_list1
    nonb_excl_list  => enefunc%nonb_excl_list

    ncell_local = domain%num_cell_local
    cutoff2     = cutoff*cutoff
    minimum_contact =  enefunc%minimum_contact

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, j, k, ij, ix, iy, ini_excl, fin_excl, num_excl,       &
    !$omp         dij, rij2, L, R, coef, work, term_elec,  ccr)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    ! calculate energy and gradient
    !
    do i = id+1, ncell_local, nthread

      num_excl = 0

      do ix = 1, natom(i) - 1

        ini_excl = num_excl + 1
        fin_excl = num_excl + num_nonb_excl1(ix,i)
        num_excl = fin_excl

        ! skip for zero charge
        !
        if (abs(charge(ix,i)) < EPS) cycle

        do k = ini_excl, fin_excl

          iy = nonb_excl_list1(k,i)

          ! skip for zero charge
          !
          if (abs(charge(iy,i)) < EPS) cycle

          ! compute distance
          !
          dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,i)
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          rij2     = max(rij2, minimum_contact)

          ! energy and gradient
          !
          rij2 = cutoff2*density/rij2
          L    = int(rij2)
          R    = rij2 - L

          term_elec = table_ecor(L) + R*(table_ecor(L+1)-table_ecor(L))
          ccr       = charge(ix,i) * charge(iy,i) * term_elec
          eelec(id+1) = eelec(id+1) + ccr

          term_elec = table_decor(L) + R*(table_decor(L+1)-table_decor(L))
          coef      = charge(ix,i) * charge(iy,i) * term_elec

          work(1:3) = coef*dij(1:3)

          force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
          force(1:3,iy,i,id+1) = force(1:3,iy,i,id+1) + work(1:3)

        end do
  
      end do

    end do

    do ij = id+1, maxcell_near, nthread

      i = cell_pair(1,ij)
      j = cell_pair(2,ij)

      num_excl = 0

      do ix = 1, natom(i)
        ini_excl = num_excl + 1
        fin_excl = num_excl + num_nonb_excl(ix,ij)
        num_excl = fin_excl

        ! skip for zero charge
        !
        if (abs(charge(ix,i)) < EPS) cycle

        do k = ini_excl, fin_excl

          iy = nonb_excl_list(k,ij)

          ! skip for zero charge
          !
          if (abs(charge(iy,j)) < EPS) cycle

          ! compute distance
          !
          dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,j)
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          rij2     = max(rij2, minimum_contact)

          ! energy and gradient
          !
          rij2 = cutoff2*density/rij2
          L    = int(rij2)
          R    = rij2 - L

          term_elec = table_ecor(L) + R*(table_ecor(L+1)-table_ecor(L))

          ccr       = charge(ix,i) * charge(iy,j) * term_elec
          eelec(id+1) = eelec(id+1) + ccr

          term_elec = table_decor(L) + R*(table_decor(L+1)-table_decor(L))
          coef      = charge(ix,i) * charge(iy,j) * term_elec

          work(1:3) = coef*dij(1:3)
          force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
          force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + work(1:3)

        end do

      end do
    end do

    !$omp end parallel

    return

  end subroutine pme_bond_corr_linear_general_check

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_linear_gro_amber_check
  !> @brief        calculate bodn correction relating with 14 scaling
  !  @authors      JJ, CK
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_linear_gro_amber_check(domain, enefunc, &
                                          force, eelec)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: dij(3), rij2, work(3), cutoff2
    real(wp)                 :: R, term_elec, coef, cc, qq_scale
    real(wp)                 :: minimum_contact
    integer                  :: i, ix, iy, j, k, ij, L
    integer                  :: num_nb14, ini_nb14, fin_nb14
    integer                  :: ncell_local, id, omp_get_thread_num

    real(dp),        pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: table_ene(:), table_grad(:)
    real(wp),        pointer :: density, cutoff
    integer,         pointer :: cell_pair(:,:), atmcls(:,:), natom(:)
    integer,         pointer :: num_nb14_calc1(:,:), nb14_calc_list1(:,:)
    integer,         pointer :: num_nb14_calc(:,:), nb14_calc_list(:,:)


    coord           => domain%coord
    charge          => domain%charge
    cell_pair       => domain%cell_pairlist1
    natom           => domain%num_atom

    table_ene       => enefunc%table%table_ecor
    table_grad      => enefunc%table%table_decor
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    num_nb14_calc1  => enefunc%num_nb14_calc1
    num_nb14_calc   => enefunc%num_nb14_calc
    nb14_calc_list1 => enefunc%nb14_calc_list1
    nb14_calc_list  => enefunc%nb14_calc_list

    cutoff2         =  cutoff * cutoff
    ncell_local     =  domain%num_cell_local
    minimum_contact =  enefunc%minimum_contact

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, j, k, ij, ix, iy, ini_nb14, fin_nb14, num_nb14,       &
    !$omp         dij, rij2, L, R, coef, work, term_elec, cc, qq_scale)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    ! calculate energy and gradient
    !
    do i = id+1, ncell_local, nthread
      num_nb14 = 0

      do ix = 1, natom(i) - 1

        ini_nb14 = num_nb14 + 1
        fin_nb14 = num_nb14 + num_nb14_calc1(ix,i)
        num_nb14 = fin_nb14

        ! skip for zero charge
        !
        if (abs(charge(ix,i)) < EPS) cycle

        do k = ini_nb14, fin_nb14

          iy = nb14_calc_list1(k,i)

          ! skip for zero charge
          !
          if (abs(charge(iy,i)) < EPS) cycle

          ! compute distance
          !
          dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,i)
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          rij2     = max(rij2, minimum_contact)

          ! energy and gradient
          !
          rij2 = cutoff2*density/rij2
          L    = int(rij2)
          R    = rij2 - L

          qq_scale  = -enefunc%nb14_qq_scale1(k,i)+1.0_wp
          term_elec = table_ene(L) + R*(table_ene(L+1)-table_ene(L))
          cc        = charge(ix,i)*charge(iy,i)*qq_scale
          eelec(id+1) = eelec(id+1) + term_elec*cc

          term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
          coef      = cc*term_elec

          work(1:3) = coef*dij(1:3)

          force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
          force(1:3,iy,i,id+1) = force(1:3,iy,i,id+1) + work(1:3)

        end do

      end do
    end do

    do ij = id+1, maxcell_near, nthread

      i = cell_pair(1,ij)
      j = cell_pair(2,ij)

      num_nb14 = 0

      do ix = 1, natom(i)

        ini_nb14 = num_nb14 + 1
        fin_nb14 = num_nb14 + num_nb14_calc(ix,ij)
        num_nb14 = fin_nb14

        ! skip for zero charge
        !
        if (abs(charge(ix,i)) < EPS) cycle

        do k = ini_nb14, fin_nb14

          iy = nb14_calc_list(k,ij)

          ! skip for zero charge
          !
          if (abs(charge(iy,j)) < EPS) cycle

          ! compute distance
          !
          dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,j)
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          rij2     = max(rij2, minimum_contact)

          ! energy and gradient
          !
          rij2 = cutoff2*density/rij2
          L    = int(rij2)
          R    = rij2 - L

          qq_scale  = -enefunc%nb14_qq_scale(k,ij)+1.0_wp
          cc        = charge(ix,i)*charge(iy,j)*qq_scale
          term_elec = table_ene(L) + R*(table_ene(L+1)-table_ene(L))
          eelec(id+1) = eelec(id+1) + term_elec*cc

          term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
          coef      = cc*term_elec

          work(1:3) = coef*dij(1:3)

          force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
          force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + work(1:3)

        end do

      end do

    end do

    !$omp end parallel

    return

  end subroutine pme_bond_corr_linear_gro_amber_check

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_linear_fep
  !> @brief        calculate pme bond correction with linear lookup table
  !                for FEP calculation
  !  @authors      NK
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_linear_fep(domain, enefunc, force, eelec)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)

    if (enefunc%nonb_limiter) then

      call pme_bond_corr_linear_general_check_fep(domain, enefunc, &
                             force, eelec)

      if (enefunc%forcefield /= ForcefieldCHARMM) then

        call pme_bond_corr_linear_gro_amber_check_fep(domain, enefunc, &
                                          force, eelec)
      end if
    else

      call pme_bond_corr_linear_general_fep(domain, enefunc, &
                                            force, eelec)

      if (enefunc%forcefield /= ForcefieldCHARMM) then
        call pme_bond_corr_linear_gro_amber_fep(domain, enefunc, &
                                                force, eelec)
      end if
    endif

    return

   end subroutine pme_bond_corr_linear_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_linear_general_fep
  !> @brief        calculate bond correction term in PME (general)
  !                for FEP calculation
  !! @authors      NK
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[inout] force   : forces of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_linear_general_fep(domain, enefunc, &
                                              force, eelec)

    ! formal arguments 
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: dij(3), rij2, work(3), cutoff2
    real(wp)                 :: R, term_elec, ccr, coef
    integer                  :: i, ix, iy, j, k, ij, L, iwater
    integer                  :: num_excl, ini_excl, fin_excl
    integer                  :: id, omp_get_thread_num
    integer                  :: ncell_local

    real(dp),        pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: table_ecor(:), table_decor(:)
    real(wp),        pointer :: density, cutoff
    integer,         pointer :: cell_pair(:,:)
    integer,         pointer :: natom(:), nwater(:), water_list(:,:,:)
    integer,         pointer :: num_nonb_excl1(:,:), nonb_excl_list1(:,:)
    integer,         pointer :: num_nonb_excl(:,:), nonb_excl_list(:,:)

    ! FEP
    real(wp)                 :: lambel
    real(wp)                 :: eelec_tmp
    integer                  :: fg1, fg2
    integer,         pointer :: fepgrp(:,:)
    real(wp),        pointer :: table_nonb_lambda(:,:,:)

    cell_pair       => domain%cell_pairlist1
    natom           => domain%num_atom
    nwater          => domain%num_water
    water_list      => domain%water_list
    coord           => domain%coord
    charge          => domain%charge

    table_ecor      => enefunc%table%table_ecor
    table_decor     => enefunc%table%table_decor
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist

    num_nonb_excl1  => enefunc%num_nonb_excl1
    num_nonb_excl   => enefunc%num_nonb_excl
    nonb_excl_list1 => enefunc%nonb_excl_list1
    nonb_excl_list  => enefunc%nonb_excl_list

    ncell_local = domain%num_cell_local
    cutoff2     = cutoff*cutoff

    ! FEP
    fepgrp            => domain%fepgrp
    table_nonb_lambda => enefunc%table_nonb_lambda

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, j, k, ij, ix, iy, ini_excl, fin_excl, num_excl,       &
    !$omp         dij, rij2, L, R, coef, work, term_elec,  ccr, lambel,        &
    !$omp         fg1, fg2, eelec_tmp)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    ! calculate energy and gradient
    !
    do i = id+1, ncell_local, nthread

      num_excl = 0

      do ix = 1, natom(i) - 1

        ini_excl = num_excl + 1
        fin_excl = num_excl + num_nonb_excl1(ix,i)
        num_excl = fin_excl

        ! skip for zero charge
        !
        if (abs(charge(ix,i)) < EPS) cycle

        do k = ini_excl, fin_excl

          iy = nonb_excl_list1(k,i)

          ! skip for zero charge
          !
          if (abs(charge(iy,i)) < EPS) cycle

          ! compute distance
          !
          dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,i)
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          ! FEP: Determine lamblj, lambel, and softcore
          fg1 = fepgrp(ix,i)
          fg2 = fepgrp(iy,i)
          lambel = table_nonb_lambda(2,fg1,fg2)

          ! energy and gradient
          !
          rij2 = cutoff2*density/rij2
          L    = int(rij2)
          R    = rij2 - L

          term_elec = table_ecor(L) + R*(table_ecor(L+1)-table_ecor(L))
          ccr       = charge(ix,i) * charge(iy,i) * term_elec
          eelec_tmp = lambel*ccr
          eelec(id+1) = eelec(id+1) + eelec_tmp

          term_elec = table_decor(L) + R*(table_decor(L+1)-table_decor(L))
          coef      = charge(ix,i) * charge(iy,i) * term_elec

          work(1:3) = lambel*coef*dij(1:3)

          force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
          force(1:3,iy,i,id+1) = force(1:3,iy,i,id+1) + work(1:3)

        end do

      end do

    end do

    do ij = id+1, maxcell_near, nthread

      i = cell_pair(1,ij)
      j = cell_pair(2,ij)

      num_excl = 0

      do ix = 1, natom(i)
        ini_excl = num_excl + 1
        fin_excl = num_excl + num_nonb_excl(ix,ij)
        num_excl = fin_excl

        ! skip for zero charge
        !
        if (abs(charge(ix,i)) < EPS) cycle

        do k = ini_excl, fin_excl

          iy = nonb_excl_list(k,ij)

          ! skip for zero charge
          !
          if (abs(charge(iy,j)) < EPS) cycle

          ! compute distance
          !
          dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,j)
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          ! FEP: Determine lamblj, lambel, and softcore
          fg1 = fepgrp(ix,i)
          fg2 = fepgrp(iy,j)
          lambel = table_nonb_lambda(2,fg1,fg2)

          ! energy and gradient
          !
          rij2 = cutoff2*density/rij2
          L    = int(rij2)
          R    = rij2 - L

          term_elec = table_ecor(L) + R*(table_ecor(L+1)-table_ecor(L))

          ccr       = charge(ix,i) * charge(iy,j) * term_elec
          eelec_tmp = lambel*ccr
          eelec(id+1) = eelec(id+1) + eelec_tmp

          term_elec = table_decor(L) + R*(table_decor(L+1)-table_decor(L))
          coef      = charge(ix,i) * charge(iy,j) * term_elec

          work(1:3) = lambel*coef*dij(1:3)
          force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
          force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + work(1:3)

        end do

      end do
    end do

    !$omp end parallel

    return

  end subroutine pme_bond_corr_linear_general_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_linear_general_check_fep
  !> @brief        calculate bond correction term in PME (general)
  !                for FEP calculations
  !! @authors      NK
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[inout] force   : forces of target systems
  !! @param[inout] eelec   : electrostatic energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_linear_general_check_fep(domain, enefunc, &
                                              force, eelec)

    ! formal arguments 
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: dij(3), rij2, work(3), cutoff2
    real(wp)                 :: R, term_elec, ccr, coef
    real(wp)                 :: minimum_contact
    integer                  :: i, ix, iy, j, k, ij, L, iwater
    integer                  :: num_excl, ini_excl, fin_excl
    integer                  :: id, omp_get_thread_num
    integer                  :: ncell_local

    real(dp),        pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: table_ecor(:), table_decor(:)
    real(wp),        pointer :: density, cutoff
    integer,         pointer :: cell_pair(:,:)
    integer,         pointer :: natom(:), nwater(:), water_list(:,:,:)
    integer,         pointer :: num_nonb_excl1(:,:), nonb_excl_list1(:,:)
    integer,         pointer :: num_nonb_excl(:,:), nonb_excl_list(:,:)

    ! FEP
    real(wp)                 :: lambel
    real(wp)                 :: eelec_tmp
    integer                  :: fg1, fg2
    integer,         pointer :: fepgrp(:,:)
    real(wp),        pointer :: table_nonb_lambda(:,:,:)

    cell_pair       => domain%cell_pairlist1
    natom           => domain%num_atom
    nwater          => domain%num_water
    water_list      => domain%water_list
    coord           => domain%coord
    charge          => domain%charge

    table_ecor      => enefunc%table%table_ecor
    table_decor     => enefunc%table%table_decor
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist

    num_nonb_excl1  => enefunc%num_nonb_excl1
    num_nonb_excl   => enefunc%num_nonb_excl
    nonb_excl_list1 => enefunc%nonb_excl_list1
    nonb_excl_list  => enefunc%nonb_excl_list

    ncell_local     = domain%num_cell_local
    cutoff2         = cutoff*cutoff
    minimum_contact =  enefunc%minimum_contact

    ! FEP
    fepgrp            => domain%fepgrp
    table_nonb_lambda => enefunc%table_nonb_lambda

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, j, k, ij, ix, iy, ini_excl, fin_excl, num_excl,       &
    !$omp         dij, rij2, L, R, coef, work, term_elec,  ccr, lambel,        &
    !$omp         eelec_tmp, fg1, fg2)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    ! calculate energy and gradient
    !
    do i = id+1, ncell_local, nthread

      num_excl = 0

      do ix = 1, natom(i) - 1

        ini_excl = num_excl + 1
        fin_excl = num_excl + num_nonb_excl1(ix,i)
        num_excl = fin_excl

        ! skip for zero charge
        !
        if (abs(charge(ix,i)) < EPS) cycle

        do k = ini_excl, fin_excl

          iy = nonb_excl_list1(k,i)

          ! skip for zero charge
          !
          if (abs(charge(iy,i)) < EPS) cycle

          ! compute distance
          !
          dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,i)
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          rij2     = max(rij2, minimum_contact)

          ! FEP: Determine lamblj, lambel, and softcore
          fg1 = fepgrp(ix,i)
          fg2 = fepgrp(iy,i)
          lambel = table_nonb_lambda(2,fg1,fg2)

          ! energy and gradient
          !
          rij2 = cutoff2*density/rij2
          L    = int(rij2)
          R    = rij2 - L

          term_elec = table_ecor(L) + R*(table_ecor(L+1)-table_ecor(L))
          ccr       = charge(ix,i) * charge(iy,i) * term_elec
          eelec_tmp = lambel*ccr
          eelec(id+1) = eelec(id+1) + eelec_tmp

          term_elec = table_decor(L) + R*(table_decor(L+1)-table_decor(L))
          coef      = charge(ix,i) * charge(iy,i) * term_elec

          work(1:3) = lambel*coef*dij(1:3)

          force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
          force(1:3,iy,i,id+1) = force(1:3,iy,i,id+1) + work(1:3)

        end do
  
      end do

    end do

    do ij = id+1, maxcell_near, nthread

      i = cell_pair(1,ij)
      j = cell_pair(2,ij)

      num_excl = 0

      do ix = 1, natom(i)
        ini_excl = num_excl + 1
        fin_excl = num_excl + num_nonb_excl(ix,ij)
        num_excl = fin_excl

        ! skip for zero charge
        !
        if (abs(charge(ix,i)) < EPS) cycle

        do k = ini_excl, fin_excl

          iy = nonb_excl_list(k,ij)

          ! skip for zero charge
          !
          if (abs(charge(iy,j)) < EPS) cycle

          ! compute distance
          !
          dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,j)
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          rij2     = max(rij2, minimum_contact)

          ! FEP: Determine lamblj, lambel, and softcore
          fg1 = fepgrp(ix,i)
          fg2 = fepgrp(iy,j)
          lambel = table_nonb_lambda(2,fg1,fg2)

          ! energy and gradient
          !
          rij2 = cutoff2*density/rij2
          L    = int(rij2)
          R    = rij2 - L

          term_elec = table_ecor(L) + R*(table_ecor(L+1)-table_ecor(L))

          ccr       = charge(ix,i) * charge(iy,j) * term_elec
          eelec_tmp = lambel*ccr
          eelec(id+1) = eelec(id+1) + eelec_tmp

          term_elec = table_decor(L) + R*(table_decor(L+1)-table_decor(L))
          coef      = charge(ix,i) * charge(iy,j) * term_elec

          work(1:3) = lambel*coef*dij(1:3)
          force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
          force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + work(1:3)

        end do

      end do
    end do

    !$omp end parallel

    return

  end subroutine pme_bond_corr_linear_general_check_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_linear_gro_amber_fep
  !> @brief        calculate bodn correction relating with 14 scaling
  !                for FEP calculations
  !  @authors      NK, HO
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_linear_gro_amber_fep(domain, enefunc, &
                                              force, eelec)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: dij(3), rij2, work(3), cutoff2
    real(wp)                 :: R, term_elec, coef, cc, qq_scale
    integer                  :: i, ix, iy, j, k, ij, L
    integer                  :: num_nb14, ini_nb14, fin_nb14
    integer                  :: ncell_local, id, omp_get_thread_num

    real(dp),        pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: table_ene(:), table_grad(:)
    real(wp),        pointer :: density, cutoff
    real(wp),        pointer :: qq_scale_tmp(:,:), qq_scale1_tmp(:,:)
    integer,         pointer :: cell_pair(:,:), atmcls(:,:), natom(:)
    integer,         pointer :: num_nb14_calc1(:,:), nb14_calc_list1(:,:)
    integer,         pointer :: num_nb14_calc(:,:), nb14_calc_list(:,:)

    ! FEP
    real(wp)                 :: lambel
    real(wp)                 :: eelec_tmp
    integer                  :: fg1, fg2
    integer,         pointer :: fepgrp(:,:)
    real(wp),        pointer :: table_nonb_lambda(:,:,:)

    coord           => domain%coord
    charge          => domain%charge
    cell_pair       => domain%cell_pairlist1
    natom           => domain%num_atom

    table_ene       => enefunc%table%table_ecor
    table_grad      => enefunc%table%table_decor
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist

    num_nb14_calc1  => enefunc%num_nb14_calc1
    num_nb14_calc   => enefunc%num_nb14_calc
    nb14_calc_list1 => enefunc%nb14_calc_list1
    nb14_calc_list  => enefunc%nb14_calc_list
    qq_scale1_tmp   => enefunc%nb14_qq_scale1
    qq_scale_tmp    => enefunc%nb14_qq_scale

    cutoff2         =  cutoff * cutoff
    ncell_local     =  domain%num_cell_local

    ! FEP
    fepgrp            => domain%fepgrp
    table_nonb_lambda => enefunc%table_nonb_lambda

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, j, k, ij, ix, iy, ini_nb14, fin_nb14, num_nb14,       &
    !$omp         dij, rij2, L, R, coef, work, term_elec, cc, qq_scale,        &
    !$omp         lambel, eelec_tmp, fg1, fg2)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    ! calculate energy and gradient
    !
    do i = id+1, ncell_local, nthread

      num_nb14 = 0

      do ix = 1, natom(i) - 1

        ini_nb14 = num_nb14 + 1
        fin_nb14 = num_nb14 + num_nb14_calc1(ix,i)
        num_nb14 = fin_nb14

        ! skip for zero charge
        !
        if (abs(charge(ix,i)) < EPS) cycle

        do k = ini_nb14, fin_nb14

          iy = nb14_calc_list1(k,i)

          ! skip for zero charge
          !
          if (abs(charge(iy,i)) < EPS) cycle

          ! compute distance
          !
          dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,i)
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          ! FEP: Determine lamblj, lambel, and softcore
          fg1 = fepgrp(ix,i)
          fg2 = fepgrp(iy,i)
          lambel = table_nonb_lambda(2,fg1,fg2)

          ! energy and gradient
          !
          rij2 = cutoff2*density/rij2
          L    = int(rij2)
          R    = rij2 - L

          qq_scale  = -qq_scale1_tmp(k,i)+1.0_wp
          term_elec = table_ene(L) + R*(table_ene(L+1)-table_ene(L))
          cc        = charge(ix,i)*charge(iy,i)*qq_scale
          eelec_tmp = lambel*term_elec*cc
          eelec(id+1) = eelec(id+1) + eelec_tmp

          term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
          coef      = cc*term_elec

          work(1:3) = lambel*coef*dij(1:3)

          force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
          force(1:3,iy,i,id+1) = force(1:3,iy,i,id+1) + work(1:3)

        end do

      end do
    end do

    do ij = id+1, maxcell_near, nthread

      i = cell_pair(1,ij)
      j = cell_pair(2,ij)

      num_nb14 = 0

      do ix = 1, natom(i)
        ini_nb14 = num_nb14 + 1
        fin_nb14 = num_nb14 + num_nb14_calc(ix,ij)
        num_nb14 = fin_nb14

        ! skip for zero charge
        !
        if (abs(charge(ix,i)) < EPS) cycle

        do k = ini_nb14, fin_nb14

          iy = nb14_calc_list(k,ij)

          ! skip for zero charge
          !
          if (abs(charge(iy,j)) < EPS) cycle

          ! compute distance
          !
          dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,j)
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          ! FEP: Determine lamblj, lambel, and softcore
          fg1 = fepgrp(ix,i)
          fg2 = fepgrp(iy,j)
          lambel = table_nonb_lambda(2,fg1,fg2)

          ! energy and gradient
          !
          rij2 = cutoff2*density/rij2
          L    = int(rij2)
          R    = rij2 - L

          qq_scale  = -qq_scale_tmp(k,ij)+1.0_wp
          cc        = charge(ix,i)*charge(iy,j)*qq_scale
          term_elec = table_ene(L) + R*(table_ene(L+1)-table_ene(L))
          eelec_tmp = lambel*term_elec*cc
          eelec(id+1) = eelec(id+1) + eelec_tmp

          term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
          coef      = cc*term_elec

          work(1:3) = lambel*coef*dij(1:3)

          force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
          force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + work(1:3)

        end do

      end do

    end do

    !$omp end parallel

    return

  end subroutine pme_bond_corr_linear_gro_amber_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_bond_corr_linear_gro_amber_check_fep
  !> @brief        calculate bodn correction relating with 14 scaling
  !                for FEP calculations
  !  @authors      NK
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_bond_corr_linear_gro_amber_check_fep(domain, enefunc, &
                                          force, eelec)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: dij(3), rij2, work(3), cutoff2
    real(wp)                 :: R, term_elec, coef, cc, qq_scale
    real(wp)                 :: minimum_contact
    integer                  :: i, ix, iy, j, k, ij, L
    integer                  :: num_nb14, ini_nb14, fin_nb14
    integer                  :: ncell_local, id, omp_get_thread_num

    real(dp),        pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: table_ene(:), table_grad(:)
    real(wp),        pointer :: density, cutoff
    real(wp),        pointer :: qq_scale_tmp(:,:), qq_scale1_tmp(:,:)
    integer,         pointer :: cell_pair(:,:), atmcls(:,:), natom(:)
    integer,         pointer :: num_nb14_calc1(:,:), nb14_calc_list1(:,:)
    integer,         pointer :: num_nb14_calc(:,:), nb14_calc_list(:,:)

    ! FEP
    real(wp)                 :: lambel
    real(wp)                 :: eelec_tmp
    integer                  :: fg1, fg2
    integer,         pointer :: fepgrp(:,:)
    real(wp),        pointer :: table_nonb_lambda(:,:,:)

    coord           => domain%coord
    charge          => domain%charge
    cell_pair       => domain%cell_pairlist1
    natom           => domain%num_atom

    table_ene       => enefunc%table%table_ecor
    table_grad      => enefunc%table%table_decor
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist

    num_nb14_calc1  => enefunc%num_nb14_calc1
    num_nb14_calc   => enefunc%num_nb14_calc
    nb14_calc_list1 => enefunc%nb14_calc_list1
    nb14_calc_list  => enefunc%nb14_calc_list
    qq_scale1_tmp   => enefunc%nb14_qq_scale1
    qq_scale_tmp    => enefunc%nb14_qq_scale

    cutoff2         =  cutoff * cutoff
    ncell_local     =  domain%num_cell_local
    minimum_contact =  enefunc%minimum_contact

    ! FEP
    fepgrp            => domain%fepgrp
    table_nonb_lambda => enefunc%table_nonb_lambda

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, j, k, ij, ix, iy, ini_nb14, fin_nb14, num_nb14,       &
    !$omp         dij, rij2, L, R, coef, work, term_elec, cc, qq_scale,        &
    !$omp         lambel, eelec_tmp, fg1, fg2)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    ! calculate energy and gradient
    !
    do i = id+1, ncell_local, nthread
      num_nb14 = 0

      do ix = 1, natom(i) - 1

        ini_nb14 = num_nb14 + 1
        fin_nb14 = num_nb14 + num_nb14_calc1(ix,i)
        num_nb14 = fin_nb14

        ! skip for zero charge
        !
        if (abs(charge(ix,i)) < EPS) cycle

        do k = ini_nb14, fin_nb14

          iy = nb14_calc_list1(k,i)

          ! skip for zero charge
          !
          if (abs(charge(iy,i)) < EPS) cycle

          ! compute distance
          !
          dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,i)
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          rij2     = max(rij2, minimum_contact)

          ! FEP: Determine lamblj, lambel, and softcore
          fg1 = fepgrp(ix,i)
          fg2 = fepgrp(iy,i)
          lambel = table_nonb_lambda(2,fg1,fg2)

          ! energy and gradient
          !
          rij2 = cutoff2*density/rij2
          L    = int(rij2)
          R    = rij2 - L

          qq_scale  = -qq_scale1_tmp(k,i)+1.0_wp
          term_elec = table_ene(L) + R*(table_ene(L+1)-table_ene(L))
          cc        = charge(ix,i)*charge(iy,i)*qq_scale
          eelec_tmp = lambel*term_elec*cc
          eelec(id+1) = eelec(id+1) + eelec_tmp

          term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
          coef      = cc*term_elec

          work(1:3) = lambel*coef*dij(1:3)

          force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
          force(1:3,iy,i,id+1) = force(1:3,iy,i,id+1) + work(1:3)

        end do

      end do
    end do

    do ij = id+1, maxcell_near, nthread

      i = cell_pair(1,ij)
      j = cell_pair(2,ij)

      num_nb14 = 0

      do ix = 1, natom(i)

        ini_nb14 = num_nb14 + 1
        fin_nb14 = num_nb14 + num_nb14_calc(ix,ij)
        num_nb14 = fin_nb14

        ! skip for zero charge
        !
        if (abs(charge(ix,i)) < EPS) cycle

        do k = ini_nb14, fin_nb14

          iy = nb14_calc_list(k,ij)

          ! skip for zero charge
          !
          if (abs(charge(iy,j)) < EPS) cycle

          ! compute distance
          !
          dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,j)
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          rij2     = max(rij2, minimum_contact)

          ! FEP: Determine lamblj, lambel, and softcore
          fg1 = fepgrp(ix,i)
          fg2 = fepgrp(iy,j)
          lambel = table_nonb_lambda(2,fg1,fg2)

          ! energy and gradient
          !
          rij2 = cutoff2*density/rij2
          L    = int(rij2)
          R    = rij2 - L

          qq_scale  = -qq_scale_tmp(k,ij)+1.0_wp
          cc        = charge(ix,i)*charge(iy,j)*qq_scale
          term_elec = table_ene(L) + R*(table_ene(L+1)-table_ene(L))
          eelec_tmp = lambel*term_elec*cc
          eelec(id+1) = eelec(id+1) + eelec_tmp

          term_elec = table_grad(L) + R*(table_grad(L+1)-table_grad(L))
          coef      = cc*term_elec

          work(1:3) = lambel*coef*dij(1:3)

          force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
          force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + work(1:3)

        end do

      end do

    end do

    !$omp end parallel

    return

  end subroutine pme_bond_corr_linear_gro_amber_check_fep

end module sp_energy_table_linear_bondcorr_mod
