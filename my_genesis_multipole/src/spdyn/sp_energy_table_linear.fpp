!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_energy_table_linear_mod
!> @brief   calculate nonbonded energy with table and with linear interpolation
!! @authors Jaewoon Jung(JJ)
!  
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_energy_table_linear_mod

  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use mpi_parallel_mod
  use constants_mod
  use timers_mod

  implicit none
  private

  ! subroutines
  !
  public  :: compute_energy_nonbond14_table_linear
  private :: compute_energy_nonbond14_table_linear_charmm
  private :: compute_energy_nonbond14_table_linear_gro_amber
  private :: compute_energy_nonbond14_table_linear_charmm_check
  private :: compute_energy_nonbond14_table_linear_gro_amber_check
  ! FEP
  public  :: compute_energy_nonbond14_table_linear_fep
  private :: compute_energy_nonbond14_table_linear_charmm_fep
  private :: compute_energy_nonbond14_table_linear_charmm_check_fep
  private :: compute_energy_nonbond14_table_linear_gro_amber_fep
  private :: compute_energy_nonbond14_table_linear_gro_amber_check_fep

contains
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_table_linear
  !> @brief        calculate nonbonded14 energy with PME and lookup table
  !!               (linear)
  !  @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_table_linear(domain, enefunc, &
                                               force, eelec, evdw)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw(nthread)

    ! ==> Type 6/7
    if (enefunc%forcefield == ForcefieldCHARMM) then

      if (enefunc%nonb_limiter) then
        call compute_energy_nonbond14_table_linear_charmm_check(domain, &
                                            enefunc, force, eelec, evdw)
      else
        call compute_energy_nonbond14_table_linear_charmm(domain, enefunc, &
                                            force, eelec, evdw)
      endif

    ! ==> Type 12/13
    else ! ForcefieldAMBER, ForcefieldGROAMBER, ForcefieldGROMARTINI
    
      if (enefunc%nonb_limiter) then
        call compute_energy_nonbond14_table_linear_gro_amber_check(domain,  &
                                            enefunc, force, eelec, evdw)
      else
        call compute_energy_nonbond14_table_linear_gro_amber(domain, enefunc, &
                                            force, eelec, evdw)
      endif
    
    end if

    return
  end subroutine compute_energy_nonbond14_table_linear

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_table_linear_charmm
  !> @brief        calculate nonbonded14 energy with PME and lookup table
  !!               (linear)
  !  @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_table_linear_charmm(domain, enefunc, &
                                               force, eelec, evdw)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                 :: dij(1:3), rij2
    real(wp)                 :: R, table(6)
    real(wp)                 :: term_lj12, term_lj6, term_elec
    real(wp)                 :: cutoff2, grad_coef, lj12, lj6
    real(wp)                 :: work(1:3)
    real(wp)                 :: force_local(1:3)
    integer                  :: i, ix, iy, j, k, ij, m, L
    integer                  :: num_nb14, ini_nb14, fin_nb14
    integer                  :: ncell_local, id, omp_get_thread_num

    real(dp),        pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: nb14_lj12(:,:), nb14_lj6(:,:)
    real(wp),        pointer :: table_ene(:), table_grad(:)
    real(wp),        pointer :: density, cutoff
    integer,         pointer :: cell_pair(:,:), atmcls(:,:), natom(:)
    integer,         pointer :: num_nb14_calc1(:,:), nb14_calc_list1(:,:)
    integer,         pointer :: num_nb14_calc(:,:), nb14_calc_list(:,:)


    coord           => domain%coord
    charge          => domain%charge
    cell_pair       => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom

    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nb14_lj12       => enefunc%nb14_lj12
    nb14_lj6        => enefunc%nb14_lj6
    num_nb14_calc1  => enefunc%num_nb14_calc1
    num_nb14_calc   => enefunc%num_nb14_calc
    nb14_calc_list1 => enefunc%nb14_calc_list1
    nb14_calc_list  => enefunc%nb14_calc_list

    cutoff2         =  cutoff * cutoff
    ncell_local     =  domain%num_cell_local

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, j, k, ij, ix, iy, m, ini_nb14, fin_nb14, num_nb14,    &
    !$omp         dij, rij2, L, R, term_lj12, term_lj6, term_elec, grad_coef,  &
    !$omp         work, table, lj12, lj6)
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

        do k = ini_nb14, fin_nb14

          iy = nb14_calc_list1(k,i)

          ! compute distance
          !
          dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,i)
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          ! energy and gradient
          !
          rij2  = cutoff2*density/rij2
          L     = int(rij2)
          R     = rij2 - L
          lj6   = nb14_lj6 (atmcls(ix,i),atmcls(iy,i))
          lj12  = nb14_lj12(atmcls(ix,i),atmcls(iy,i))

          table(1:6)  = table_ene(3*L-2:3*L+3)
          term_lj12   = table(1) + R*(table(4)-table(1))
          term_lj6    = table(2) + R*(table(5)-table(2))
          term_elec   = table(3) + R*(table(6)-table(3))
          evdw(id+1)  = evdw(id+1) + term_lj12*lj12 - term_lj6*lj6
          eelec(id+1) = eelec(id+1)+ charge(ix,i)*charge(iy,i)*term_elec

          table(1:6)  = table_grad(3*L-2:3*L+3)
          term_lj12   = table(1) + R*(table(4)-table(1))
          term_lj6    = table(2) + R*(table(5)-table(2))
          term_elec   = table(3) + R*(table(6)-table(3))
          grad_coef   = term_lj12*lj12 - term_lj6*lj6 +     &
                        charge(ix,i)*charge(iy,i)*term_elec

          work(1:3) = grad_coef*dij(1:3)

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

        do k = ini_nb14, fin_nb14

          iy = nb14_calc_list(k,ij)

          ! compute distance
          !
          dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,j)
          rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          ! energy and gradient
          !
          rij2  = cutoff2*density/rij2
          L     = int(rij2)
          R     = rij2 - L
          lj6   = nb14_lj6(atmcls(ix,i),atmcls(iy,j))
          lj12  = nb14_lj12(atmcls(ix,i),atmcls(iy,j))

          table(1:6)  = table_ene(3*L-2:3*L+3)
          term_lj12   = table(1) + R*(table(4)-table(1))
          term_lj6    = table(2) + R*(table(5)-table(2))
          term_elec   = table(3) + R*(table(6)-table(3))
          evdw(id+1)  = evdw(id+1) + term_lj12*lj12 - term_lj6*lj6
          eelec(id+1) = eelec(id+1) + charge(ix,i)*charge(iy,j)*term_elec

          table(1:6)  = table_grad(3*L-2:3*L+3)
          term_lj12   = table(1) + R*(table(4)-table(1))
          term_lj6    = table(2) + R*(table(5)-table(2))
          term_elec   = table(3) + R*(table(6)-table(3))
          grad_coef   = term_lj12*lj12 - term_lj6*lj6 + &
                        charge(ix,i)*charge(iy,j)*term_elec

          work(1:3) = grad_coef*dij(1:3)

          force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
          force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + work(1:3)

        end do

      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond14_table_linear_charmm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_table_linear_gro_amber
  !> @brief        calculate nonbonded14 energy with PME and lookup table
  !!               (linear)
  !  @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_table_linear_gro_amber(domain, enefunc, &
                                                             force, eelec, evdw)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                 :: dij(1:3), rij2
    real(wp)                 :: R, table(6)
    real(wp)                 :: term_lj12, term_lj6, term_elec
    real(wp)                 :: cutoff2, grad_coef, lj12, lj6
    real(wp)                 :: work(1:3)
    real(wp)                 :: lj_scale, qq_scale, cc
    real(wp)                 :: inv_r12, inv_r121, inv_r123, inv_r126, inv_r1212
    integer                  :: i, ix, iy, j, k, ij, m, L
    integer                  :: num_nb14, ini_nb14, fin_nb14
    integer                  :: ncell_local, id, omp_get_thread_num

    real(dp),        pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: nb14_lj12(:,:), nb14_lj6(:,:)
    real(wp),        pointer :: table_ene(:), table_grad(:)
    real(wp),        pointer :: density, cutoff
    integer,         pointer :: cell_pair(:,:), atmcls(:,:), natom(:)
    integer,         pointer :: num_nb14_calc1(:,:), nb14_calc_list1(:,:)
    integer,         pointer :: num_nb14_calc(:,:), nb14_calc_list(:,:)


    coord           => domain%coord
    charge          => domain%charge
    cell_pair       => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom

    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nb14_lj12       => enefunc%nb14_lj12
    nb14_lj6        => enefunc%nb14_lj6
    num_nb14_calc1  => enefunc%num_nb14_calc1
    num_nb14_calc   => enefunc%num_nb14_calc
    nb14_calc_list1 => enefunc%nb14_calc_list1
    nb14_calc_list  => enefunc%nb14_calc_list

    cutoff2         =  cutoff * cutoff
    ncell_local     =  domain%num_cell_local


    !$omp parallel default(shared)                                             &
    !$omp private(id, i, j, k, ij, ix, iy, m, ini_nb14, fin_nb14, num_nb14,    &
    !$omp         dij, rij2, L, R, term_lj12, term_lj6, term_elec, grad_coef,  &
    !$omp         work, table, lj12, lj6, qq_scale, lj_scale, cc, inv_r12,     &
    !$omp         inv_r121, inv_r123, inv_r126, inv_r1212)
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

        do k = ini_nb14, fin_nb14

          iy = nb14_calc_list1(k,i)

          ! compute distance
          !
          dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,i)
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          qq_scale = enefunc%nb14_qq_scale1(k,i)
          lj_scale = enefunc%nb14_lj_scale1(k,i)
          inv_r12  = 1.0_wp / rij2
          inv_r121 = sqrt(inv_r12)
          inv_r123 = inv_r12 * inv_r121
          inv_r126 = inv_r123 * inv_r123
          inv_r1212 = inv_r126 * inv_r126

          ! energy and gradient
          !
          rij2  = cutoff2*density/rij2
          L     = int(rij2)
          R     = rij2 - L
          lj6   = nb14_lj6 (atmcls(ix,i),atmcls(iy,i))
          lj12  = nb14_lj12(atmcls(ix,i),atmcls(iy,i))
          cc    = charge(ix,i)*charge(iy,i)*qq_scale

          term_lj12   = inv_r1212 
          term_lj6    = inv_r126
          term_elec   = table_ene(3*L) + R*(table_ene(3*L+3)-table_ene(3*L))
          evdw(id+1)  = evdw(id+1) + (term_lj12*lj12 - term_lj6*lj6)*lj_scale
          eelec(id+1) = eelec(id+1) + cc*term_elec

          term_lj12   = -12.0_wp * inv_r1212 * inv_r12
          term_lj6    = -6.0_wp * inv_r126 * inv_r12
          term_elec   = table_grad(3*L) + R*(table_grad(3*L+3)-table_grad(3*L))
          grad_coef   = (term_lj12*lj12 - term_lj6*lj6)*lj_scale + cc*term_elec

          work(1:3) = grad_coef*dij(1:3)

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

        do k = ini_nb14, fin_nb14

          iy = nb14_calc_list(k,ij)

          ! compute distance
          !
          dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,j)
          rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          inv_r12 = 1.0_wp / rij2
          inv_r121 = sqrt(inv_r12)
          inv_r123 = inv_r12 * inv_r121
          inv_r126 = inv_r123 * inv_r123
          inv_r1212 = inv_r126 * inv_r126

          qq_scale = enefunc%nb14_qq_scale(k,ij)
          lj_scale = enefunc%nb14_lj_scale(k,ij)

          ! energy and gradient
          !
          rij2  = cutoff2*density/rij2
          L     = int(rij2)
          R     = rij2 - L
          lj6   = nb14_lj6(atmcls(ix,i),atmcls(iy,j))
          lj12  = nb14_lj12(atmcls(ix,i),atmcls(iy,j))
          cc    = charge(ix,i)*charge(iy,j)*qq_scale

          term_lj12   = inv_r1212
          term_lj6    = inv_r126
          term_elec   = table_ene(3*L) + R*(table_ene(3*L+3)-table_ene(3*L))
          evdw(id+1)  = evdw(id+1) + (term_lj12*lj12 - term_lj6*lj6)*lj_scale
          eelec(id+1) = eelec(id+1) + cc*term_elec

          term_lj12   = -12.0_wp * inv_r1212 * inv_r12
          term_lj6    = -6.0_wp * inv_r126 * inv_r12
          term_elec   = table_grad(3*L) + R*(table_grad(3*L+3)-table_grad(3*L))
          grad_coef   = (term_lj12*lj12 - term_lj6*lj6)*lj_scale + cc*term_elec

          work(1:3) = grad_coef*dij(1:3)

          force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
          force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + work(1:3)

        end do

      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond14_table_linear_gro_amber

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_table_linear_charmm_check
  !> @brief        calculate nonbonded14 energy with PME and lookup table
  !!               (linear)
  !  @authors      JJ, CK
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_table_linear_charmm_check(domain, &
                                               enefunc, force, eelec, evdw)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                 :: dij(1:3), rij2
    real(wp)                 :: R, table(6)
    real(wp)                 :: term_lj12, term_lj6, term_elec
    real(wp)                 :: cutoff2, grad_coef, lj12, lj6
    real(wp)                 :: work(1:3)
    real(wp)                 :: force_local(1:3)
    real(wp)                 :: minimum_contact
    integer                  :: i, ix, iy, j, k, ij, m, L
    integer                  :: num_nb14, ini_nb14, fin_nb14
    integer                  :: ncell_local, id, omp_get_thread_num

    real(dp),        pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: nb14_lj12(:,:), nb14_lj6(:,:)
    real(wp),        pointer :: table_ene(:), table_grad(:)
    real(wp),        pointer :: density, cutoff
    integer,         pointer :: cell_pair(:,:), atmcls(:,:), natom(:)
    integer,         pointer :: num_nb14_calc1(:,:), nb14_calc_list1(:,:)
    integer,         pointer :: num_nb14_calc(:,:), nb14_calc_list(:,:)


    coord           => domain%coord
    charge          => domain%charge
    cell_pair       => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom

    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nb14_lj12       => enefunc%nb14_lj12
    nb14_lj6        => enefunc%nb14_lj6
    num_nb14_calc1  => enefunc%num_nb14_calc1
    num_nb14_calc   => enefunc%num_nb14_calc
    nb14_calc_list1 => enefunc%nb14_calc_list1
    nb14_calc_list  => enefunc%nb14_calc_list

    cutoff2         =  cutoff * cutoff
    ncell_local     =  domain%num_cell_local
    minimum_contact =  enefunc%minimum_contact

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, j, k, ij, ix, iy, m, ini_nb14, fin_nb14, num_nb14,    &
    !$omp         dij, rij2, L, R, term_lj12, term_lj6, term_elec, grad_coef,  &
    !$omp         work, table, lj12, lj6)
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

        do k = ini_nb14, fin_nb14

          iy = nb14_calc_list1(k,i)

          ! compute distance
          !
          dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,i)
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          rij2     = max(rij2, minimum_contact)

          ! energy and gradient
          !
          rij2  = cutoff2*density/rij2
          L     = int(rij2)
          R     = rij2 - L
          lj6   = nb14_lj6 (atmcls(ix,i),atmcls(iy,i))
          lj12  = nb14_lj12(atmcls(ix,i),atmcls(iy,i))

          table(1:6)  = table_ene(3*L-2:3*L+3)
          term_lj12   = table(1) + R*(table(4)-table(1))
          term_lj6    = table(2) + R*(table(5)-table(2))
          term_elec   = table(3) + R*(table(6)-table(3))
          evdw(id+1)  = evdw(id+1) + term_lj12*lj12 - term_lj6*lj6
          eelec(id+1) = eelec(id+1)+ charge(ix,i)*charge(iy,i)*term_elec

          table(1:6)  = table_grad(3*L-2:3*L+3)
          term_lj12   = table(1) + R*(table(4)-table(1))
          term_lj6    = table(2) + R*(table(5)-table(2))
          term_elec   = table(3) + R*(table(6)-table(3))
          grad_coef   = term_lj12*lj12 - term_lj6*lj6 +     &
                        charge(ix,i)*charge(iy,i)*term_elec

          work(1:3) = grad_coef*dij(1:3)

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

        do k = ini_nb14, fin_nb14

          iy = nb14_calc_list(k,ij)

          ! compute distance
          !
          dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,j)
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          rij2     = max(rij2, minimum_contact)

          ! energy and gradient
          !
          rij2  = cutoff2*density/rij2
          L     = int(rij2)
          R     = rij2 - L
          lj6   = nb14_lj6(atmcls(ix,i),atmcls(iy,j))
          lj12  = nb14_lj12(atmcls(ix,i),atmcls(iy,j))

          table(1:6)  = table_ene(3*L-2:3*L+3)
          term_lj12   = table(1) + R*(table(4)-table(1))
          term_lj6    = table(2) + R*(table(5)-table(2))
          term_elec   = table(3) + R*(table(6)-table(3))
          evdw(id+1)  = evdw(id+1) + term_lj12*lj12 - term_lj6*lj6
          eelec(id+1) = eelec(id+1) + charge(ix,i)*charge(iy,j)*term_elec

          table(1:6)  = table_grad(3*L-2:3*L+3)
          term_lj12   = table(1) + R*(table(4)-table(1))
          term_lj6    = table(2) + R*(table(5)-table(2))
          term_elec   = table(3) + R*(table(6)-table(3))
          grad_coef   = term_lj12*lj12 - term_lj6*lj6 + &
                        charge(ix,i)*charge(iy,j)*term_elec

          work(1:3) = grad_coef*dij(1:3)

          force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
          force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + work(1:3)

        end do

      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond14_table_linear_charmm_check

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_table_linear_gro_amber_check
  !> @brief        calculate nonbonded14 energy with PME and lookup table
  !!               (linear)
  !  @authors      JJ, CK
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_table_linear_gro_amber_check(domain,  &
                                               enefunc, force, eelec, evdw)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                 :: dij(1:3), rij2
    real(wp)                 :: R, table(6)
    real(wp)                 :: term_lj12, term_lj6, term_elec
    real(wp)                 :: cutoff2, grad_coef, lj12, lj6
    real(wp)                 :: work(1:3)
    real(wp)                 :: lj_scale, qq_scale, cc
    real(wp)                 :: inv_r12, inv_r121, inv_r123, inv_r126, inv_r1212
    real(wp)                 :: minimum_contact
    integer                  :: i, ix, iy, j, k, ij, m, L
    integer                  :: num_nb14, ini_nb14, fin_nb14
    integer                  :: ncell_local, id, omp_get_thread_num

    real(dp),        pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: nb14_lj12(:,:), nb14_lj6(:,:)
    real(wp),        pointer :: table_ene(:), table_grad(:)
    real(wp),        pointer :: density, cutoff
    integer,         pointer :: cell_pair(:,:), atmcls(:,:), natom(:)
    integer,         pointer :: num_nb14_calc1(:,:), nb14_calc_list1(:,:)
    integer,         pointer :: num_nb14_calc(:,:), nb14_calc_list(:,:)


    coord           => domain%coord
    charge          => domain%charge
    cell_pair       => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom

    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nb14_lj12       => enefunc%nb14_lj12
    nb14_lj6        => enefunc%nb14_lj6
    num_nb14_calc1  => enefunc%num_nb14_calc1
    num_nb14_calc   => enefunc%num_nb14_calc
    nb14_calc_list1 => enefunc%nb14_calc_list1
    nb14_calc_list  => enefunc%nb14_calc_list

    cutoff2         =  cutoff * cutoff
    ncell_local     =  domain%num_cell_local
    minimum_contact =  enefunc%minimum_contact


    !$omp parallel default(shared)                                             &
    !$omp private(id, i, j, k, ij, ix, iy, m, ini_nb14, fin_nb14, num_nb14,    &
    !$omp         dij, rij2, L, R, term_lj12, term_lj6, term_elec, grad_coef,  &
    !$omp         work, table, lj12, lj6, qq_scale, lj_scale, cc, inv_r12,     &
    !$omp         inv_r121, inv_r123, inv_r126, inv_r1212)
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

        do k = ini_nb14, fin_nb14

          iy = nb14_calc_list1(k,i)

          ! compute distance
          !
          dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,i)
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          rij2     = max(rij2, minimum_contact)
          qq_scale = enefunc%nb14_qq_scale1(k,i)
          lj_scale = enefunc%nb14_lj_scale1(k,i)
          inv_r12  = 1.0_wp / rij2
          inv_r121 = sqrt(inv_r12)
          inv_r123 = inv_r12 * inv_r121
          inv_r126 = inv_r123 * inv_r123
          inv_r1212 = inv_r126 * inv_r126

          ! energy and gradient
          !
          rij2  = cutoff2*density/rij2
          L     = int(rij2)
          R     = rij2 - L
          lj6   = nb14_lj6 (atmcls(ix,i),atmcls(iy,i))
          lj12  = nb14_lj12(atmcls(ix,i),atmcls(iy,i))
          cc    = charge(ix,i)*charge(iy,i)*qq_scale

          term_lj12   = inv_r1212 
          term_lj6    = inv_r126
          term_elec   = table_ene(3*L) + R*(table_ene(3*L+3)-table_ene(3*L))
          evdw(id+1)  = evdw(id+1) + (term_lj12*lj12 - term_lj6*lj6)*lj_scale
          eelec(id+1) = eelec(id+1) + cc*term_elec

          term_lj12   = -12.0_wp * inv_r1212 * inv_r12
          term_lj6    = -6.0_wp * inv_r126 * inv_r12
          term_elec   = table_grad(3*L) + R*(table_grad(3*L+3)-table_grad(3*L))
          grad_coef   = (term_lj12*lj12 - term_lj6*lj6)*lj_scale + cc*term_elec

          work(1:3) = grad_coef*dij(1:3)

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

        do k = ini_nb14, fin_nb14

          iy = nb14_calc_list(k,ij)

          ! compute distance
          !
          dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,j)
          rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          rij2   = max(rij2, minimum_contact)
          inv_r12 = 1.0_wp / rij2
          inv_r121 = sqrt(inv_r12)
          inv_r123 = inv_r12 * inv_r121
          inv_r126 = inv_r123 * inv_r123
          inv_r1212 = inv_r126 * inv_r126

          qq_scale = enefunc%nb14_qq_scale(k,ij)
          lj_scale = enefunc%nb14_lj_scale(k,ij)

          ! energy and gradient
          !
          rij2  = cutoff2*density/rij2
          L     = int(rij2)
          R     = rij2 - L
          lj6   = nb14_lj6(atmcls(ix,i),atmcls(iy,j))
          lj12  = nb14_lj12(atmcls(ix,i),atmcls(iy,j))
          cc    = charge(ix,i)*charge(iy,j)*qq_scale

          term_lj12   = inv_r1212
          term_lj6    = inv_r126
          term_elec   = table_ene(3*L) + R*(table_ene(3*L+3)-table_ene(3*L))
          evdw(id+1)  = evdw(id+1) + (term_lj12*lj12 - term_lj6*lj6)*lj_scale
          eelec(id+1) = eelec(id+1) + cc*term_elec

          term_lj12   = -12.0_wp * inv_r1212 * inv_r12
          term_lj6    = -6.0_wp * inv_r126 * inv_r12
          term_elec   = table_grad(3*L) + R*(table_grad(3*L+3)-table_grad(3*L))
          grad_coef   = (term_lj12*lj12 - term_lj6*lj6)*lj_scale + cc*term_elec

          work(1:3) = grad_coef*dij(1:3)

          force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
          force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + work(1:3)

        end do

      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond14_table_linear_gro_amber_check

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_table_linear_fep
  !> @brief        calculate nonbonded14 energy with PME and lookup table
  !!               (linear) for FEP calculation
  !  @authors      NK
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_table_linear_fep(domain, enefunc, &
                                               force, eelec, evdw)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw(nthread)

    ! ==> Type 6/7
    if (enefunc%forcefield == ForcefieldCHARMM) then

      if (enefunc%nonb_limiter) then
        call compute_energy_nonbond14_table_linear_charmm_check_fep(domain, &
                                            enefunc, force, &
                                            eelec, evdw)

      else
        call compute_energy_nonbond14_table_linear_charmm_fep(domain, enefunc, &
                                            force, eelec, evdw)

      endif

    ! ==> Type 12/13
    else ! ForcefieldAMBER, ForcefieldGROAMBER
    
      if (enefunc%nonb_limiter) then
        call compute_energy_nonbond14_table_linear_gro_amber_check_fep(domain, &
                                            enefunc, force, &
                                            eelec, evdw)

      else
        call compute_energy_nonbond14_table_linear_gro_amber_fep(domain, enefunc, &
                                            force, eelec, evdw)

      endif
    
    end if

    return
  end subroutine compute_energy_nonbond14_table_linear_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_table_linear_charmm_fep
  !> @brief        calculate nonbonded14 energy with PME, lookup table (linear),
  !!               and soft core for FEP calculation
  !  @authors      NK, HO
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] eelec0  : electrostatic energy of reference systems
  !! @param[inout] evdw0   : van der Waals energy of reference systems
  !! @param[inout] eelec1  : electrostatic energy of target systems
  !! @param[inout] evdw1   : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_table_linear_charmm_fep(domain, enefunc, &
                                               force, eelec, evdw)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                 :: dij(1:3)
    real(wp)                 :: rij2
    real(wp)                 :: R, table(6)
    real(wp)                 :: evdw_tmp, eelec_tmp, rij2_tmp
    real(wp)                 :: term_lj12, term_lj6, term_elec
    real(wp)                 :: cutoff2, grad_coef, lj12, lj6
    real(wp)                 :: work(1:3)
    real(wp)                 :: force_local(1:3)
    integer                  :: i, ix, iy, j, k, ij, m, L
    integer                  :: num_nb14, ini_nb14, fin_nb14
    integer                  :: ncell_local, id, omp_get_thread_num

    real(dp),        pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: nb14_lj12(:,:), nb14_lj6(:,:)
    real(wp),        pointer :: table_ene(:), table_grad(:)
    real(wp),        pointer :: density, cutoff
    integer,         pointer :: cell_pair(:,:), atmcls(:,:), natom(:)
    integer,         pointer :: num_nb14_calc1(:,:), nb14_calc_list1(:,:)
    integer,         pointer :: num_nb14_calc(:,:), nb14_calc_list(:,:)

    ! FEP
    real(wp)                 :: rij2_sc_lj, rij2_sc_el
    real(wp)                 :: lamblj, lambel
    integer                  :: fg1, fg2
    integer,         pointer :: fepgrp(:,:)
    real(wp),        pointer :: table_nonb_lambda(:,:,:)

    natom           => domain%num_atom
    coord           => domain%coord
    charge          => domain%charge
    cell_pair       => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nb14_lj12       => enefunc%nb14_lj12
    nb14_lj6        => enefunc%nb14_lj6

    num_nb14_calc1  => enefunc%num_nb14_calc1
    num_nb14_calc   => enefunc%num_nb14_calc
    nb14_calc_list1 => enefunc%nb14_calc_list1
    nb14_calc_list  => enefunc%nb14_calc_list

    cutoff2         =  cutoff * cutoff
    ncell_local     =  domain%num_cell_local

    ! FEP
    fepgrp            => domain%fepgrp
    table_nonb_lambda => enefunc%table_nonb_lambda

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, j, k, ij, ix, iy, m, ini_nb14, fin_nb14, num_nb14,    &
    !$omp         dij, rij2, L, R, term_lj12, term_lj6, term_elec, grad_coef,  &
    !$omp         work, table, lj12, lj6, evdw_tmp, eelec_tmp, rij2_tmp,       &
    !$omp         rij2_sc_lj, rij2_sc_el, lamblj, lambel, fg1, fg2)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    ! calculate energy and gradient
    !
    do i = id+1, ncell_local, nthread
    if (natom(i) == 0) cycle

      num_nb14 = 0

      do ix = 1, natom(i) - 1

        ini_nb14 = num_nb14 + 1
        fin_nb14 = num_nb14 + num_nb14_calc1(ix,i)
        num_nb14 = fin_nb14

        do k = ini_nb14, fin_nb14

          iy = nb14_calc_list1(k,i)

          ! compute distance
          !
          dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,i)
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          ! FEP: Determine lamblj, lambel, and softcore
          fg1 = fepgrp(ix,i)
          fg2 = fepgrp(iy,i)
          lamblj = table_nonb_lambda(1,fg1,fg2)
          lambel = table_nonb_lambda(2,fg1,fg2)
          rij2_sc_lj = rij2 + table_nonb_lambda(3,fg1,fg2)
          rij2_sc_el = rij2 + table_nonb_lambda(4,fg1,fg2)

          ! FEP: LJ with soft core for reference
          rij2_tmp  = cutoff2*density/rij2_sc_lj
          L     = int(rij2_tmp)
          R     = rij2_tmp - L
          lj6   = nb14_lj6 (atmcls(ix,i),atmcls(iy,i))
          lj12  = nb14_lj12(atmcls(ix,i),atmcls(iy,i))
          table(1:6)  = table_ene(3*L-2:3*L+3)
          term_lj12   = table(1) + R*(table(4)-table(1))
          term_lj6    = table(2) + R*(table(5)-table(2))
          evdw_tmp = lamblj*(term_lj12*lj12 - term_lj6*lj6)
          evdw(id+1) = evdw(id+1) + evdw_tmp
          table(1:6)  = table_grad(3*L-2:3*L+3)
          term_lj12   = table(1) + R*(table(4)-table(1))
          term_lj6    = table(2) + R*(table(5)-table(2))

          ! FEP: EL with soft core for reference
          rij2_tmp  = cutoff2*density/rij2_sc_el
          L     = int(rij2_tmp)
          R     = rij2_tmp - L
          table(1:6)  = table_ene(3*L-2:3*L+3)
          term_elec   = table(3) + R*(table(6)-table(3))
          eelec_tmp = lambel * charge(ix,i) * charge(iy,i)*term_elec
          eelec(id+1) = eelec(id+1) + eelec_tmp
          table(1:6)  = table_grad(3*L-2:3*L+3)
          term_elec   = table(3) + R*(table(6)-table(3))

          ! gradient
          grad_coef   = lamblj*(term_lj12*lj12 - term_lj6*lj6) + &
                        lambel*charge(ix,i)*charge(iy,i)*term_elec

          work(1:3) = grad_coef*dij(1:3)

          force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
          force(1:3,iy,i,id+1) = force(1:3,iy,i,id+1) + work(1:3)

        end do

      end do

    end do

    do ij = id+1, maxcell_near, nthread

      i = cell_pair(1,ij)
      j = cell_pair(2,ij)
      if (natom(i) == 0) cycle
      if (natom(j) == 0) cycle

      num_nb14 = 0

      do ix = 1, natom(i)

        ini_nb14 = num_nb14 + 1
        fin_nb14 = num_nb14 + num_nb14_calc(ix,ij)
        num_nb14 = fin_nb14

        do k = ini_nb14, fin_nb14

          iy = nb14_calc_list(k,ij)

          ! compute distance
          !
          dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,j)
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          ! FEP: Determine lamblj, lambel, and softcore
          fg1 = fepgrp(ix,i)
          fg2 = fepgrp(iy,j)
          lamblj = table_nonb_lambda(1,fg1,fg2)
          lambel = table_nonb_lambda(2,fg1,fg2)
          rij2_sc_lj = rij2 + table_nonb_lambda(3,fg1,fg2)
          rij2_sc_el = rij2 + table_nonb_lambda(4,fg1,fg2)

          ! FEP: LJ with soft core for reference
          rij2_tmp  = cutoff2*density/rij2_sc_lj
          L     = int(rij2_tmp)
          R     = rij2_tmp - L
          lj6   = nb14_lj6(atmcls(ix,i),atmcls(iy,j))
          lj12  = nb14_lj12(atmcls(ix,i),atmcls(iy,j))
          table(1:6)  = table_ene(3*L-2:3*L+3)
          term_lj12   = table(1) + R*(table(4)-table(1))
          term_lj6    = table(2) + R*(table(5)-table(2))
          evdw_tmp = lamblj*(term_lj12*lj12 - term_lj6*lj6)
          evdw(id+1) = evdw(id+1) + evdw_tmp
          table(1:6)  = table_grad(3*L-2:3*L+3)
          term_lj12   = table(1) + R*(table(4)-table(1))
          term_lj6    = table(2) + R*(table(5)-table(2))

          ! FEP: EL with soft core for reference
          rij2_tmp  = cutoff2*density/rij2_sc_el
          L     = int(rij2_tmp)
          R     = rij2_tmp - L
          table(1:6)  = table_ene(3*L-2:3*L+3)
          term_elec   = table(3) + R*(table(6)-table(3))
          eelec_tmp = lambel * charge(ix,i) * charge(iy,j)*term_elec
          eelec(id+1) = eelec(id+1) + eelec_tmp
          table(1:6)  = table_grad(3*L-2:3*L+3)
          term_elec   = table(3) + R*(table(6)-table(3))

          ! gradient
          grad_coef   = lamblj*(term_lj12*lj12 - term_lj6*lj6) + &
                        lambel*charge(ix,i)*charge(iy,j)*term_elec

          work(1:3) = grad_coef*dij(1:3)

          force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
          force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + work(1:3)

        end do

      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond14_table_linear_charmm_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_table_linear_charmm_check_fep
  !> @brief        calculate nonbonded14 energy with PME, lookup table,
  !!               and soft core for FEP calculations
  !  @authors      NK, HO
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] eelec0  : electrostatic energy of reference systems
  !! @param[inout] evdw0   : van der Waals energy of reference systems
  !! @param[inout] eelec1  : electrostatic energy of target systems
  !! @param[inout] evdw1   : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_table_linear_charmm_check_fep(domain, enefunc, &
                                               force, eelec, evdw)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                 :: dij(1:3)
    real(wp)                 :: rij2
    real(wp)                 :: R, table(6)
    real(wp)                 :: evdw_tmp, eelec_tmp, rij2_tmp
    real(wp)                 :: term_lj12, term_lj6, term_elec
    real(wp)                 :: cutoff2, grad_coef, lj12, lj6
    real(wp)                 :: work(1:3)
    real(wp)                 :: force_local(1:3)
    real(wp)                 :: minimum_contact
    integer                  :: i, ix, iy, j, k, ij, m, L
    integer                  :: num_nb14, ini_nb14, fin_nb14
    integer                  :: ncell_local, id, omp_get_thread_num

    real(dp),        pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: nb14_lj12(:,:), nb14_lj6(:,:)
    real(wp),        pointer :: table_ene(:), table_grad(:)
    real(wp),        pointer :: density, cutoff
    integer,         pointer :: cell_pair(:,:), atmcls(:,:), natom(:)
    integer,         pointer :: num_nb14_calc1(:,:), nb14_calc_list1(:,:)
    integer,         pointer :: num_nb14_calc(:,:), nb14_calc_list(:,:)

    ! FEP
    real(wp)                 :: rij2_sc_lj, rij2_sc_el
    real(wp)                 :: lamblj, lambel
    integer                  :: fg1, fg2
    integer,         pointer :: fepgrp(:,:)
    real(wp),        pointer :: table_nonb_lambda(:,:,:)

    natom           => domain%num_atom
    coord           => domain%coord
    charge          => domain%charge
    cell_pair       => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nb14_lj12       => enefunc%nb14_lj12
    nb14_lj6        => enefunc%nb14_lj6

    num_nb14_calc1  => enefunc%num_nb14_calc1
    num_nb14_calc   => enefunc%num_nb14_calc
    nb14_calc_list1 => enefunc%nb14_calc_list1
    nb14_calc_list  => enefunc%nb14_calc_list

    cutoff2         =  cutoff * cutoff
    ncell_local     =  domain%num_cell_local
    minimum_contact =  enefunc%minimum_contact

    ! FEP
    fepgrp            => domain%fepgrp
    table_nonb_lambda => enefunc%table_nonb_lambda

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, j, k, ij, ix, iy, m, ini_nb14, fin_nb14, num_nb14,    &
    !$omp         dij, rij2, L, R, term_lj12, term_lj6, term_elec, grad_coef,  &
    !$omp         work, table, lj12, lj6, evdw_tmp, eelec_tmp, rij2_tmp,       &
    !$omp         rij2_sc_lj, rij2_sc_el, lamblj, lambel, fg1, fg2)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    ! calculate energy and gradient
    !
    do i = id+1, ncell_local, nthread
    if (natom(i) == 0) cycle

      num_nb14 = 0

      do ix = 1, natom(i) - 1

        ini_nb14 = num_nb14 + 1
        fin_nb14 = num_nb14 + num_nb14_calc1(ix,i)
        num_nb14 = fin_nb14

        do k = ini_nb14, fin_nb14

          iy = nb14_calc_list1(k,i)

          ! compute distance
          !
          dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,i)
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          rij2     = max(rij2, minimum_contact)

          ! FEP: Determine lamblj, lambel, and softcore
          fg1 = fepgrp(ix,i)
          fg2 = fepgrp(iy,i)
          lamblj = table_nonb_lambda(1,fg1,fg2)
          lambel = table_nonb_lambda(2,fg1,fg2)
          rij2_sc_lj = rij2 + table_nonb_lambda(3,fg1,fg2)
          rij2_sc_el = rij2 + table_nonb_lambda(4,fg1,fg2)

          ! FEP: LJ with soft core for reference
          rij2_tmp  = cutoff2*density/rij2_sc_lj
          L     = int(rij2_tmp)
          R     = rij2_tmp - L
          lj6   = nb14_lj6 (atmcls(ix,i),atmcls(iy,i))
          lj12  = nb14_lj12(atmcls(ix,i),atmcls(iy,i))
          table(1:6)  = table_ene(3*L-2:3*L+3)
          term_lj12   = table(1) + R*(table(4)-table(1))
          term_lj6    = table(2) + R*(table(5)-table(2))
          evdw_tmp = lamblj*(term_lj12*lj12 - term_lj6*lj6)
          evdw(id+1) = evdw(id+1) + evdw_tmp
          table(1:6)  = table_grad(3*L-2:3*L+3)
          term_lj12   = table(1) + R*(table(4)-table(1))
          term_lj6    = table(2) + R*(table(5)-table(2))

          ! FEP: EL with soft core for reference
          rij2_tmp  = cutoff2*density/rij2_sc_el
          L     = int(rij2_tmp)
          R     = rij2_tmp - L
          table(1:6)  = table_ene(3*L-2:3*L+3)
          term_elec   = table(3) + R*(table(6)-table(3))
          eelec_tmp = lambel * charge(ix,i) * charge(iy,i)*term_elec
          eelec(id+1) = eelec(id+1) + eelec_tmp
          table(1:6)  = table_grad(3*L-2:3*L+3)
          term_elec   = table(3) + R*(table(6)-table(3))

          ! gradient
          grad_coef   = lamblj*(term_lj12*lj12 - term_lj6*lj6) + &
                        lambel*charge(ix,i)*charge(iy,i)*term_elec

          work(1:3) = grad_coef*dij(1:3)

          force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
          force(1:3,iy,i,id+1) = force(1:3,iy,i,id+1) + work(1:3)

        end do

      end do

    end do

    do ij = id+1, maxcell_near, nthread

      i = cell_pair(1,ij)
      j = cell_pair(2,ij)
      if (natom(i) == 0) cycle
      if (natom(j) == 0) cycle

      num_nb14 = 0

      do ix = 1, natom(i)

        ini_nb14 = num_nb14 + 1
        fin_nb14 = num_nb14 + num_nb14_calc(ix,ij)
        num_nb14 = fin_nb14

        do k = ini_nb14, fin_nb14

          iy = nb14_calc_list(k,ij)

          ! compute distance
          !
          dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,j)
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          rij2     = max(rij2, minimum_contact)

          ! FEP: Determine lamblj, lambel, and softcore
          fg1 = fepgrp(ix,i)
          fg2 = fepgrp(iy,j)
          lamblj = table_nonb_lambda(1,fg1,fg2)
          lambel = table_nonb_lambda(2,fg1,fg2)
          rij2_sc_lj = rij2 + table_nonb_lambda(3,fg1,fg2)
          rij2_sc_el = rij2 + table_nonb_lambda(4,fg1,fg2)

          ! FEP: LJ with soft core for reference
          rij2_tmp  = cutoff2*density/rij2_sc_lj
          L     = int(rij2_tmp)
          R     = rij2_tmp - L
          lj6   = nb14_lj6(atmcls(ix,i),atmcls(iy,j))
          lj12  = nb14_lj12(atmcls(ix,i),atmcls(iy,j))
          table(1:6)  = table_ene(3*L-2:3*L+3)
          term_lj12   = table(1) + R*(table(4)-table(1))
          term_lj6    = table(2) + R*(table(5)-table(2))
          evdw_tmp = lamblj*(term_lj12*lj12 - term_lj6*lj6)
          evdw(id+1) = evdw(id+1) + evdw_tmp
          table(1:6)  = table_grad(3*L-2:3*L+3)
          term_lj12   = table(1) + R*(table(4)-table(1))
          term_lj6    = table(2) + R*(table(5)-table(2))

          ! FEP: EL with soft core for reference
          rij2_tmp  = cutoff2*density/rij2_sc_el
          L     = int(rij2_tmp)
          R     = rij2_tmp - L
          table(1:6)  = table_ene(3*L-2:3*L+3)
          term_elec   = table(3) + R*(table(6)-table(3))
          eelec_tmp = lambel * charge(ix,i) * charge(iy,j)*term_elec
          eelec(id+1) = eelec(id+1) + eelec_tmp
          table(1:6)  = table_grad(3*L-2:3*L+3)
          term_elec   = table(3) + R*(table(6)-table(3))

          ! gradient
          grad_coef   = lamblj*(term_lj12*lj12 - term_lj6*lj6) + &
                        lambel*charge(ix,i)*charge(iy,j)*term_elec

          work(1:3) = grad_coef*dij(1:3)

          force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
          force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + work(1:3)

        end do

      end do

    end do


    !$omp end parallel

    return

  end subroutine compute_energy_nonbond14_table_linear_charmm_check_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_table_linear_gro_amber_fep
  !> @brief        calculate nonbonded14 energy with PME and lookup table
  !!               (linear) for FEP calculations
  !  @authors      NK, HO
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_table_linear_gro_amber_fep(domain, enefunc, &
                                               force, eelec, evdw)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                 :: dij(1:3), rij2
    real(wp)                 :: R, table(6)
    real(wp)                 :: evdw_tmp, eelec_tmp
    real(wp)                 :: term_lj12, term_lj6, term_elec
    real(wp)                 :: cutoff2, grad_coef, lj12, lj6
    real(wp)                 :: work(1:3)
    real(wp)                 :: lj_scale, qq_scale, cc
    real(wp)                 :: inv_r12, inv_r121, inv_r123, inv_r126, inv_r1212
    integer                  :: i, ix, iy, j, k, ij, m, L
    integer                  :: num_nb14, ini_nb14, fin_nb14
    integer                  :: ncell_local, id, omp_get_thread_num

    real(dp),        pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: nb14_lj12(:,:), nb14_lj6(:,:)
    real(wp),        pointer :: table_ene(:), table_grad(:)
    real(wp),        pointer :: density, cutoff
    integer,         pointer :: cell_pair(:,:), atmcls(:,:), natom(:)
    integer,         pointer :: num_nb14_calc1(:,:), nb14_calc_list1(:,:)
    integer,         pointer :: num_nb14_calc(:,:), nb14_calc_list(:,:)
    real(wp),        pointer :: nb14_lj_scale1(:,:), nb14_lj_scale(:,:)
    real(wp),        pointer :: nb14_qq_scale1(:,:), nb14_qq_scale(:,:)

    ! FEP
    real(wp)                 :: rij2_sc_lj, rij2_sc_el
    real(wp)                 :: lamblj, lambel
    integer                  :: fg1, fg2
    integer,         pointer :: fepgrp(:,:)
    real(wp),        pointer :: table_nonb_lambda(:,:,:)

    natom           => domain%num_atom
    coord           => domain%coord
    charge          => domain%charge
    cell_pair       => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no

    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nb14_lj12       => enefunc%nb14_lj12
    nb14_lj6        => enefunc%nb14_lj6

    num_nb14_calc1  => enefunc%num_nb14_calc1
    num_nb14_calc   => enefunc%num_nb14_calc
    nb14_calc_list1 => enefunc%nb14_calc_list1
    nb14_calc_list  => enefunc%nb14_calc_list
    nb14_lj_scale1  => enefunc%nb14_lj_scale1
    nb14_lj_scale   => enefunc%nb14_lj_scale
    nb14_qq_scale1  => enefunc%nb14_qq_scale1
    nb14_qq_scale   => enefunc%nb14_qq_scale

    cutoff2         =  cutoff * cutoff
    ncell_local     =  domain%num_cell_local

    ! FEP
    fepgrp            => domain%fepgrp
    table_nonb_lambda => enefunc%table_nonb_lambda

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, j, k, ij, ix, iy, m, ini_nb14, fin_nb14, num_nb14,    &
    !$omp         dij, rij2, L, R, term_lj12, term_lj6, term_elec, grad_coef,  &
    !$omp         work, table, lj12, lj6, qq_scale, lj_scale, cc, inv_r12,     &
    !$omp         inv_r121, inv_r123, inv_r126, inv_r1212,                     &
    !$omp         evdw_tmp, eelec_tmp, rij2_sc_lj, rij2_sc_el,                 &
    !$omp         lamblj, lambel, fg1, fg2)
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

        do k = ini_nb14, fin_nb14

          iy = nb14_calc_list1(k,i)

          lj_scale = nb14_lj_scale1(k,i)
          qq_scale = nb14_qq_scale1(k,i)
          lj6   = nb14_lj6 (atmcls(ix,i),atmcls(iy,i)) * lj_scale
          lj12  = nb14_lj12(atmcls(ix,i),atmcls(iy,i)) * lj_scale
          cc    = charge(ix,i)*charge(iy,i)*qq_scale

          ! compute distance
          !
          dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,i)
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          ! FEP: Determine lamblj, lambel, and softcore
          fg1 = fepgrp(ix,i)
          fg2 = fepgrp(iy,i)
          lamblj = table_nonb_lambda(1,fg1,fg2)
          lambel = table_nonb_lambda(2,fg1,fg2)
          rij2_sc_lj = rij2 + table_nonb_lambda(3,fg1,fg2)
          rij2_sc_el = rij2 + table_nonb_lambda(4,fg1,fg2)

          ! FEP: LJ with soft core for reference
          inv_r12    = 1.0_wp / rij2_sc_lj
          inv_r121   = sqrt(inv_r12)
          inv_r123   = inv_r12 * inv_r121
          inv_r126   = inv_r123 * inv_r123
          inv_r1212  = inv_r126 * inv_r126
          term_lj12  = inv_r1212
          term_lj6   = inv_r126
          evdw_tmp = lamblj*(term_lj12*lj12 - term_lj6*lj6)
          evdw(id+1) = evdw(id+1) + evdw_tmp
          term_lj12  = -12.0_wp * inv_r1212 * inv_r12
          term_lj6   = -6.0_wp * inv_r126 * inv_r12

          ! FEP: EL with soft core for reference
          rij2  = cutoff2*density/rij2_sc_el
          L     = int(rij2)
          R     = rij2 - L
          term_elec   = table_ene(3*L) + R*(table_ene(3*L+3)-table_ene(3*L))
          eelec_tmp = lambel * cc*term_elec
          eelec(id+1) = eelec(id+1) + eelec_tmp
          term_elec   = table_grad(3*L) + R*(table_grad(3*L+3)-table_grad(3*L))

          ! gradient
          grad_coef   = lamblj*(term_lj12*lj12 - term_lj6*lj6) + &
                        lambel*cc*term_elec

          work(1:3) = grad_coef*dij(1:3)

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

        do k = ini_nb14, fin_nb14

          iy = nb14_calc_list(k,ij)

          lj_scale = nb14_lj_scale(k,ij)
          qq_scale = nb14_qq_scale(k,ij)
          lj6   = nb14_lj6(atmcls(ix,i),atmcls(iy,j)) * lj_scale
          lj12  = nb14_lj12(atmcls(ix,i),atmcls(iy,j)) * lj_scale
          cc    = charge(ix,i)*charge(iy,j)*qq_scale

          ! compute distance
          !
          dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,j)
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          ! FEP: Determine lamblj, lambel, and softcore
          fg1 = fepgrp(ix,i)
          fg2 = fepgrp(iy,j)
          lamblj = table_nonb_lambda(1,fg1,fg2)
          lambel = table_nonb_lambda(2,fg1,fg2)
          rij2_sc_lj = rij2 + table_nonb_lambda(3,fg1,fg2)
          rij2_sc_el = rij2 + table_nonb_lambda(4,fg1,fg2)

          ! FEP: LJ with soft core for reference
          inv_r12    = 1.0_wp / rij2_sc_lj
          inv_r121   = sqrt(inv_r12)
          inv_r123   = inv_r12 * inv_r121
          inv_r126   = inv_r123 * inv_r123
          inv_r1212  = inv_r126 * inv_r126
          term_lj12  = inv_r1212
          term_lj6   = inv_r126
          evdw_tmp = lamblj*(term_lj12*lj12 - term_lj6*lj6)
          evdw(id+1) = evdw(id+1) + evdw_tmp
          term_lj12  = -12.0_wp * inv_r1212 * inv_r12
          term_lj6   = -6.0_wp * inv_r126 * inv_r12

          ! FEP: EL with soft core for reference
          rij2  = cutoff2*density/rij2_sc_el
          L     = int(rij2)
          R     = rij2 - L
          term_elec   = table_ene(3*L) + R*(table_ene(3*L+3)-table_ene(3*L))
          eelec_tmp = lambel * cc*term_elec
          eelec(id+1) = eelec(id+1) + eelec_tmp
          term_elec   = table_grad(3*L) + R*(table_grad(3*L+3)-table_grad(3*L))

          ! gradient
          grad_coef   = lamblj*(term_lj12*lj12 - term_lj6*lj6) + &
                        lambel*cc*term_elec

          work(1:3) = grad_coef*dij(1:3)

          force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
          force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + work(1:3)

        end do

      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond14_table_linear_gro_amber_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_table_linear_gro_amber_check_fep
  !> @brief        calculate nonbonded14 energy with PME and lookup table
  !!               (linear) for FEP calculations
  !  @authors      NK, HO
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions
  !! @param[inout] force   : forces for each cell
  !! @param[inout] eelec   : electrostatic energy of target systems
  !! @param[inout] evdw    : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_table_linear_gro_amber_check_fep(domain, &
                                               enefunc, force, &
                                               eelec, evdw)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                 :: dij(1:3), rij2
    real(wp)                 :: R, table(6)
    real(wp)                 :: evdw_tmp, eelec_tmp
    real(wp)                 :: term_lj12, term_lj6, term_elec
    real(wp)                 :: cutoff2, grad_coef, lj12, lj6
    real(wp)                 :: work(1:3)
    real(wp)                 :: lj_scale, qq_scale, cc
    real(wp)                 :: inv_r12, inv_r121, inv_r123, inv_r126, inv_r1212
    real(wp)                 :: minimum_contact
    integer                  :: i, ix, iy, j, k, ij, m, L
    integer                  :: num_nb14, ini_nb14, fin_nb14
    integer                  :: ncell_local, id, omp_get_thread_num

    real(dp),        pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    real(wp),        pointer :: nb14_lj12(:,:), nb14_lj6(:,:)
    real(wp),        pointer :: table_ene(:), table_grad(:)
    real(wp),        pointer :: density, cutoff
    integer,         pointer :: cell_pair(:,:), atmcls(:,:), natom(:)
    integer,         pointer :: num_nb14_calc1(:,:), nb14_calc_list1(:,:)
    integer,         pointer :: num_nb14_calc(:,:), nb14_calc_list(:,:)
    real(wp),        pointer :: nb14_lj_scale1(:,:), nb14_lj_scale(:,:)
    real(wp),        pointer :: nb14_qq_scale1(:,:), nb14_qq_scale(:,:)

    ! FEP
    real(wp)                 :: rij2_sc_lj, rij2_sc_el
    real(wp)                 :: lamblj, lambel
    integer                  :: fg1, fg2
    integer,         pointer :: fepgrp(:,:)
    real(wp),        pointer :: table_nonb_lambda(:,:,:)

    natom           => domain%num_atom
    coord           => domain%coord
    charge          => domain%charge
    cell_pair       => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no

    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nb14_lj12       => enefunc%nb14_lj12
    nb14_lj6        => enefunc%nb14_lj6

    num_nb14_calc1  => enefunc%num_nb14_calc1
    num_nb14_calc   => enefunc%num_nb14_calc
    nb14_calc_list1 => enefunc%nb14_calc_list1
    nb14_calc_list  => enefunc%nb14_calc_list
    nb14_lj_scale1  => enefunc%nb14_lj_scale1
    nb14_lj_scale   => enefunc%nb14_lj_scale
    nb14_qq_scale1  => enefunc%nb14_qq_scale1
    nb14_qq_scale   => enefunc%nb14_qq_scale

    cutoff2         =  cutoff * cutoff
    ncell_local     =  domain%num_cell_local
    minimum_contact =  enefunc%minimum_contact

    ! FEP
    fepgrp            => domain%fepgrp
    table_nonb_lambda => enefunc%table_nonb_lambda

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, j, k, ij, ix, iy, m, ini_nb14, fin_nb14, num_nb14,    &
    !$omp         dij, rij2, L, R, term_lj12, term_lj6, term_elec, grad_coef,  &
    !$omp         work, table, lj12, lj6, qq_scale, lj_scale, cc, inv_r12,     &
    !$omp         inv_r121, inv_r123, inv_r126, inv_r1212,                     &
    !$omp         evdw_tmp, eelec_tmp, rij2_sc_lj, rij2_sc_el,                 &
    !$omp         lamblj, lambel, fg1, fg2)
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

        do k = ini_nb14, fin_nb14

          iy = nb14_calc_list1(k,i)

          lj_scale = nb14_lj_scale1(k,i)
          qq_scale = nb14_qq_scale1(k,i)
          lj6   = nb14_lj6 (atmcls(ix,i),atmcls(iy,i)) * lj_scale
          lj12  = nb14_lj12(atmcls(ix,i),atmcls(iy,i)) * lj_scale
          cc    = charge(ix,i)*charge(iy,i)*qq_scale

          ! compute distance
          !
          dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,i)
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          rij2     = max(rij2, minimum_contact)

          ! FEP: Determine lamblj, lambel, and softcore
          fg1 = fepgrp(ix,i)
          fg2 = fepgrp(iy,i)
          lamblj = table_nonb_lambda(1,fg1,fg2)
          lambel = table_nonb_lambda(2,fg1,fg2)
          rij2_sc_lj = rij2 + table_nonb_lambda(3,fg1,fg2)
          rij2_sc_el = rij2 + table_nonb_lambda(4,fg1,fg2)

          ! FEP: LJ with soft core for reference
          inv_r12    = 1.0_wp / rij2_sc_lj
          inv_r121   = sqrt(inv_r12)
          inv_r123   = inv_r12 * inv_r121
          inv_r126   = inv_r123 * inv_r123
          inv_r1212  = inv_r126 * inv_r126
          term_lj12  = inv_r1212
          term_lj6   = inv_r126
          evdw_tmp = lamblj*(term_lj12*lj12 - term_lj6*lj6)
          evdw(id+1) = evdw(id+1) + evdw_tmp
          term_lj12  = -12.0_wp * inv_r1212 * inv_r12
          term_lj6   = -6.0_wp * inv_r126 * inv_r12

          ! FEP: EL with soft core for reference
          rij2  = cutoff2*density/rij2_sc_el
          L     = int(rij2)
          R     = rij2 - L
          term_elec   = table_ene(3*L) + R*(table_ene(3*L+3)-table_ene(3*L))
          eelec_tmp = lambel * cc*term_elec
          eelec(id+1) = eelec(id+1) + eelec_tmp
          term_elec   = table_grad(3*L) + R*(table_grad(3*L+3)-table_grad(3*L))

          ! gradient
          grad_coef   = lamblj*(term_lj12*lj12 - term_lj6*lj6) + &
                        lambel*cc*term_elec

          work(1:3) = grad_coef*dij(1:3)

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

        do k = ini_nb14, fin_nb14

          iy = nb14_calc_list(k,ij)

          lj_scale = nb14_lj_scale(k,ij)
          qq_scale = nb14_qq_scale(k,ij)
          lj6   = nb14_lj6(atmcls(ix,i),atmcls(iy,j)) * lj_scale
          lj12  = nb14_lj12(atmcls(ix,i),atmcls(iy,j)) * lj_scale
          cc    = charge(ix,i)*charge(iy,j)*qq_scale

          ! compute distance
          !
          dij(1:3) = coord(1:3,ix,i) - coord(1:3,iy,j)
          rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          rij2     = max(rij2, minimum_contact)

          ! FEP: Determine lamblj, lambel, and softcore
          fg1 = fepgrp(ix,i)
          fg2 = fepgrp(iy,j)
          lamblj = table_nonb_lambda(1,fg1,fg2)
          lambel = table_nonb_lambda(2,fg1,fg2)
          rij2_sc_lj = rij2 + table_nonb_lambda(3,fg1,fg2)
          rij2_sc_el = rij2 + table_nonb_lambda(4,fg1,fg2)

          ! FEP: LJ with soft core for reference
          inv_r12    = 1.0_wp / rij2_sc_lj
          inv_r121   = sqrt(inv_r12)
          inv_r123   = inv_r12 * inv_r121
          inv_r126   = inv_r123 * inv_r123
          inv_r1212  = inv_r126 * inv_r126
          term_lj12  = inv_r1212
          term_lj6   = inv_r126
          evdw_tmp = lamblj*(term_lj12*lj12 - term_lj6*lj6)
          evdw(id+1) = evdw(id+1) + evdw_tmp
          term_lj12  = -12.0_wp * inv_r1212 * inv_r12
          term_lj6   = -6.0_wp * inv_r126 * inv_r12

          ! FEP: EL with soft core for reference
          rij2  = cutoff2*density/rij2_sc_el
          L     = int(rij2)
          R     = rij2 - L
          term_elec   = table_ene(3*L) + R*(table_ene(3*L+3)-table_ene(3*L))
          eelec_tmp = lambel * cc*term_elec
          eelec(id+1) = eelec(id+1) + eelec_tmp
          term_elec   = table_grad(3*L) + R*(table_grad(3*L+3)-table_grad(3*L))

          ! gradient
          grad_coef   = lamblj*(term_lj12*lj12 - term_lj6*lj6) + &
                        lambel*cc*term_elec

          work(1:3) = grad_coef*dij(1:3)

          force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) - work(1:3)
          force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + work(1:3)

        end do

      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond14_table_linear_gro_amber_check_fep

end module sp_energy_table_linear_mod
