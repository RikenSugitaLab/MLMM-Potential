!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_energy_table_linear_nowater_mod
!> @brief   calculate nonbonded energy with table and with linear interpolation
!! @authors Jaewoon Jung(JJ)
!  
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_energy_table_linear_nowater_mod

  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use mpi_parallel_mod
  use timers_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  !
  public  :: compute_energy_nonbond_table_solute_linear
  public  :: compute_force_nonbond_table_solute_linear
#ifndef USE_GPU
  public  :: compute_energy_nonbond_table_linear_check
  public  :: compute_energy_nonbond_table_linear
  public  :: compute_force_nonbond_table_linear
  ! FEP
  public  :: compute_energy_nonbond_table_linear_check_fep
  public  :: compute_energy_nonbond_table_linear_fep
  public  :: compute_force_nonbond_table_linear_fep
#else
  public  :: compute_energy_nonbond_table_linear_gpu
  public  :: compute_force_nonbond_table_linear_gpu
  private :: cpu_compute_force_intra_cell_univ
  ! FEP
  public  :: compute_energy_nonbond_table_linear_gpu_fep
  public  :: compute_force_nonbond_table_linear_gpu_fep
  private :: cpu_compute_force_intra_cell_univ_fep
#endif


contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_table_solute_linear
  !> @brief        calculate nonbonded energy without solvents with lookup table
  !! @authors      JJ 
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] virial   : virial term of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : van der Waals energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_table_solute_linear(   &
                                                 domain, enefunc, &
                                                 pairlist, coord_pbc, force, &
                                                 virial, eelec, evdw)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:)
    real(dp),                 intent(inout) :: eelec(nthread)
    real(dp),                 intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: lj6, lj12
    real(wp)                  :: R
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: trans_x, trans_y, trans_z
    real(wp)                  :: rtmp(1:3), qtmp
    real(wp)                  :: elec_temp, evdw_temp
    real(wp)                  :: dij_list(3,MaxAtom)
    real(wp)                  :: force_local(3), rij2_list(MaxAtom)
    real(wp)                  :: force_localj(3,MaxAtom)
    real(wp)                  :: ieps, jeps, eps, irmin, jrmin, rmin
    integer                   :: i, ix, iy, j, k, ij, iix, iwater, list, L, L1
    integer                   :: num_nb15, ini_nb15, fin_nb15, num_count
    integer                   :: id, omp_get_thread_num, j_list(MaxAtom)
    integer                   :: ncell, ncell_local
    integer                   :: iatmcls, jatmcls
    integer                   :: check_virial

    real(dp),         pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: cell_move(:,:,:), system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: density, cutoff
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer,          pointer :: cell_pairlist(:,:), atmcls(:,:)
    integer,          pointer :: natom(:), nsolute(:), nwater(:)
    integer,          pointer :: nb15_list(:,:), nb15_listw(:,:)
    integer,          pointer :: nb15_cell(:), nb15_cellw(:)
    integer,          pointer :: num_nb15_calc1(:,:), num_nb15_calc(:,:)
    integer,          pointer :: num_nb15_calcw(:,:)
    integer,          pointer :: nb15_calc_list1(:,:), nb15_calc_list(:,:)
    integer,          pointer :: virial_check(:,:)


    cell_pairlist   => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    nsolute         => domain%num_solute
    nwater          => domain%num_water
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6

    nb15_list       => pairlist%nb15_list
    nb15_listw      => pairlist%nb15_listw
    nb15_cell       => pairlist%nb15_cell
    nb15_cellw      => pairlist%nb15_cellw
    num_nb15_calc1  => pairlist%num_nb15_calc1
    num_nb15_calc   => pairlist%num_nb15_calc
    num_nb15_calcw  => pairlist%num_nb15_calcw
    nb15_calc_list1 => pairlist%nb15_calc_list1
    nb15_calc_list  => pairlist%nb15_calc_list

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff


    ! calculate energy and gradient
    !

    !$omp parallel default(shared)                                          &
    !$omp private(id, i, ix, num_nb15, rtmp, qtmp, ini_nb15, fin_nb15,      &
    !$omp         k, iy, dij, rij2, L, R, term_lj12, term_lj6, grad_coef,   &
    !$omp         work, term_elec, iwater, ij, j, trans_x,                  &
    !$omp         trans_y, trans_z, iix, list, num_count, j_list, dij_list, &
    !$omp         rij2_list, force_local, force_localj, lj6, lj12,  L1,     &
    !$omp         elec_temp, evdw_temp, iatmcls, jatmcls, ieps, jeps, eps,  &
    !$omp         irmin, jrmin, rmin, check_virial)             
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell, nthread
      do ix = 1, natom(i)
        coord_pbc(1:3,ix,i) = coord(1:3,ix,i) + trans1(1:3,ix,i)
      end do
    end do

    !$omp barrier

    ! energy within a cell(solute-solute)
    !
    do i = id+1, ncell_local, nthread

      num_nb15 = 0
      do ix = 1, nsolute(i) - 1

        rtmp(1:3) = coord_pbc(1:3,ix,i)
        qtmp = charge(ix,i)

        ini_nb15  = num_nb15 + 1
        fin_nb15  = num_nb15 + num_nb15_calc1(ix,i)
        num_nb15  = fin_nb15
        num_count = 0

        force_local(1:3) = 0.0_wp
        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp

!ocl norecurrence(force)
!ocl swp
        do k = ini_nb15, fin_nb15

          iy = nb15_calc_list1(k,i)

          ! compute distance
          !
          dij(1:3) = rtmp(1:3) - coord_pbc(1:3,iy,i)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1:3,num_count)  = dij(1:3)
            rij2_list(num_count) = rij2
            j_list(num_count)   = iy
          end if
        end do

        do k = 1, num_count
          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)
          iy   = j_list(k)
          rij2  = cutoff2*density/rij2_list(k)
          lj6  = nonb_lj6(atmcls(ix,i),atmcls(iy,i))
          lj12 = nonb_lj12(atmcls(ix,i),atmcls(iy,i))

          L    = int(rij2)
          R    = rij2 - L

          L1   = 3*L - 2
          term_lj12 = table_ene(L1)   + R*(table_ene(L1+3)-table_ene(L1))
          term_lj6  = table_ene(L1+1) + R*(table_ene(L1+4)-table_ene(L1+1))
          term_elec = table_ene(L1+2) + R*(table_ene(L1+5)-table_ene(L1+2))
          evdw_temp = evdw_temp + term_lj12*lj12 - term_lj6*lj6
          elec_temp = elec_temp + qtmp*charge(iy,i)*term_elec

          term_lj12 = table_grad(L1)   + R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6 &
                    + qtmp*charge(iy,i)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force_localj(1,k)  = work(1)
          force_localj(2,k)  = work(2)
          force_localj(3,k)  = work(3)

        end do

        do k = 1, num_count
          iy = j_list(k)
          force(1:3,iy,i,id+1) = force(1:3,iy,i,id+1) + force_localj(1:3,k)
        end do
        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw(id+1)  = evdw(id+1)  + evdw_temp
      end do

      do ix = 1, nsolute(i)

        rtmp(1:3) = coord_pbc(1:3,ix,i)
        qtmp = charge(ix,i)
        num_count = 0
        force_local(1:3) = 0.0_wp
        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp

!ocl norecurrence(force)
!ocl swp
        do iy = nsolute(i)+1, natom(i)
          ! compute distance
          !
          dij(1) = rtmp(1) - coord_pbc(1,iy,i)
          dij(2) = rtmp(2) - coord_pbc(2,iy,i)
          dij(3) = rtmp(3) - coord_pbc(3,iy,i)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1:3,num_count)  = dij(1:3)
            rij2_list(num_count) = rij2
            j_list(num_count)   = iy
          end if
        end do

        do k = 1, num_count
          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)
          iy   = j_list(k)
          rij2  = cutoff2*density/rij2_list(k)
          lj6  = nonb_lj6 (atmcls(ix,i),atmcls(iy,i))
          lj12 = nonb_lj12(atmcls(ix,i),atmcls(iy,i))

          L    = int(rij2)
          R    = rij2 - L

          L1   = 3*L - 2
          term_lj12 = table_ene(L1) + R*(table_ene(L1+3)-table_ene(L1))
          term_lj6  = table_ene(L1+1) + R*(table_ene(L1+4)-table_ene(L1+1))
          term_elec = table_ene(L1+2) + R*(table_ene(L1+5)-table_ene(L1+2))
          evdw_temp = evdw_temp + term_lj12*lj12 - term_lj6*lj6
          elec_temp = elec_temp + qtmp*charge(iy,i)*term_elec

          term_lj12 = table_grad(L1) + R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6  &
                    + qtmp*charge(iy,i)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force_localj(1,k)  = work(1)
          force_localj(2,k)  = work(2)
          force_localj(3,k)  = work(3)

        end do

        do k = 1, num_count
          iy = j_list(k)
          force(1:3,iy,i,id+1) = force(1:3,iy,i,id+1) + force_localj(1:3,k)
        end do
        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw(id+1) = evdw(id+1) + evdw_temp
      end do
    end do


    ! interaction between different cells
    !
    do ij = id+1, maxcell, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)

      num_nb15 = 0

      trans_x = cell_move(1,j,i) * system_size(1)
      trans_y = cell_move(2,j,i) * system_size(2)
      trans_z = cell_move(3,j,i) * system_size(3)
      check_virial = virial_check(j,i)

      do iix = 1, nb15_cell(ij)

        ix   = nb15_list(iix,ij)
        rtmp(1:3) = coord_pbc(1:3,ix,i)
        qtmp = charge(ix,i)

        ini_nb15 = num_nb15 + 1
        fin_nb15 = num_nb15 + num_nb15_calc(ix,ij)
        num_nb15 = fin_nb15

        num_count = 0
        force_local(1:3) = 0.0_wp
        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp
!ocl norecurrence(force)
!ocl swp
        do k = ini_nb15, fin_nb15
          iy   = nb15_calc_list(k,ij)
          dij(1) = rtmp(1) - coord_pbc(1,iy,j) + trans_x
          dij(2) = rtmp(2) - coord_pbc(2,iy,j) + trans_y
          dij(3) = rtmp(3) - coord_pbc(3,iy,j) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1:3,num_count)  = dij(1:3)
            rij2_list(num_count) = rij2
            j_list(num_count)   = iy
          end if
        end do

        do k = 1, num_count
          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)
          iy   = j_list(k)
          rij2  = cutoff2*density/rij2_list(k)
          lj6  = nonb_lj6 (atmcls(ix,i),atmcls(iy,j))
          lj12 = nonb_lj12(atmcls(ix,i),atmcls(iy,j))

          L    = int(rij2)
          R    = rij2 - L

          L1 = 3*L - 2
          term_lj12 = table_ene(L1) + R*(table_ene(L1+3)-table_ene(L1))
          term_lj6  = table_ene(L1+1) + R*(table_ene(L1+4)-table_ene(L1+1))
          term_elec = table_ene(L1+2) + R*(table_ene(L1+5)-table_ene(L1+2))
          evdw_temp = evdw_temp + term_lj12*lj12 - term_lj6*lj6
          elec_temp = elec_temp + qtmp*charge(iy,j)*term_elec

          term_lj12 = table_grad(L1) + R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6  &
                    + qtmp*charge(iy,j)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force_localj(1,k) = work(1)
          force_localj(2,k) = work(2)
          force_localj(3,k) = work(3)

        end do

        do k = 1, num_count
          iy = j_list(k)
          force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + force_localj(1:3,k)
        end do
        if (check_virial == 1) &
          virial(1:3,ij) = virial(1:3,ij) - force_local(1:3)
        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw(id+1) = evdw(id+1) + evdw_temp

      end do

      do iix = 1, nb15_cellw(ij)

        iwater   = nb15_listw(iix,ij)
        ini_nb15 = num_nb15 + 1
        fin_nb15 = num_nb15 + num_nb15_calcw(iwater,ij)
        num_nb15 = fin_nb15

        do list = 1, 3
          ix = nsolute(i) + 3*(iwater-1) + list

          num_count = 0
          force_local(1:3) = 0.0_wp
          elec_temp = 0.0_wp
          evdw_temp = 0.0_wp

!ocl norecurrence(force)
!ocl swp
          do k = ini_nb15, fin_nb15
            iy = nb15_calc_list(k,ij)

            ! compute distance
            !
            dij(1) = coord_pbc(1,ix,i) - coord_pbc(1,iy,j) + trans_x
            dij(2) = coord_pbc(2,ix,i) - coord_pbc(2,iy,j) + trans_y
            dij(3) = coord_pbc(3,ix,i) - coord_pbc(3,iy,j) + trans_z
            rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
            if (rij2 < cutoff2) then
              num_count = num_count + 1
              dij_list(1:3,num_count)  = dij(1:3)
              rij2_list(num_count) = rij2
              j_list(num_count)   = iy
            end if
          end do

          do k = 1, num_count
            dij(1) = dij_list(1,k)
            dij(2) = dij_list(2,k)
            dij(3) = dij_list(3,k)
            iy   = j_list(k)
            rij2  = cutoff2*density/rij2_list(k)
            lj6  = nonb_lj6 (atmcls(ix,i),atmcls(iy,j))
            lj12 = nonb_lj12(atmcls(ix,i),atmcls(iy,j))

            L    = int(rij2)
            R    = rij2 - L

            L1   = 3*L - 2
            term_lj12 = table_ene(L1) + R*(table_ene(L1+3)-table_ene(L1))
            term_lj6  = table_ene(L1+1) + R*(table_ene(L1+4)-table_ene(L1+1))
            term_elec = table_ene(L1+2) + R*(table_ene(L1+5)-table_ene(L1+2))
            evdw_temp = evdw_temp + term_lj12*lj12 - term_lj6*lj6
            elec_temp = elec_temp + charge(ix,i)*charge(iy,j)*term_elec

            term_lj12 = table_grad(L1) + R*(table_grad(L1+3)-table_grad(L1))
            term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
            term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
            grad_coef = term_lj12*lj12 - term_lj6*lj6  &
                      + charge(ix,i)*charge(iy,j)*term_elec

            work(1) = grad_coef*dij(1)
            work(2) = grad_coef*dij(2)
            work(3) = grad_coef*dij(3)

            ! store force
            !
            !
            force_local(1) = force_local(1) - work(1)
            force_local(2) = force_local(2) - work(2)
            force_local(3) = force_local(3) - work(3)
            force_localj(1,k) = work(1)
            force_localj(2,k) = work(2)
            force_localj(3,k) = work(3)

          end do

          do k = 1, num_count
            iy = j_list(k)
            force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + force_localj(1:3,k)
          end do
          if (check_virial == 1) &
            virial(1:3,ij) = virial(1:3,ij) - force_local(1:3)
          force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
          eelec(id+1) = eelec(id+1) + elec_temp
          evdw(id+1) = evdw(id+1) + evdw_temp

        end do
      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond_table_solute_linear

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_table_solute_linear
  !> @brief        calculate nonbonded force without solvents with lookup table
  !! @authors      JJ 
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] virial   : virial term of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

#ifndef KCOMP
  subroutine compute_force_nonbond_table_solute_linear( &
                                                 domain, enefunc, &
                                                 pairlist, coord_pbc, &
                                                 force, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:)

    ! local variables
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: R, lj12, lj6
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: trans_x, trans_y, trans_z
    real(wp)                  :: rtmp(1:3), qtmp
    real(wp)                  :: dij_list(4,MaxAtom)
    real(wp)                  :: force_local(3)
    real(wp)                  :: force_localj(3,MaxAtom)
    real(wp)                  :: ieps, jeps, eps, irmin, jrmin, rmin
    integer                   :: i, ix, iy, j, k, ij, iix, iwater, list, L, L1
    integer                   :: num_nb15, ini_nb15, fin_nb15, num_count
    integer                   :: id, omp_get_thread_num, j_list(MaxAtom)
    integer                   :: ncell, ncell_local
    integer                   :: iatmcls,jatmcls
    integer                   :: check_virial

    real(dp),         pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: cell_move(:,:,:), system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: density, cutoff
    real(wp),         pointer :: table_grad(:)
    integer,          pointer :: cell_pairlist(:,:), atmcls(:,:)
    integer,          pointer :: natom(:), nsolute(:), nwater(:)
    integer,          pointer :: nb15_list(:,:), nb15_listw(:,:)
    integer,          pointer :: nb15_cell(:), nb15_cellw(:)
    integer,          pointer :: num_nb15_calc1(:,:), num_nb15_calc(:,:)
    integer,          pointer :: num_nb15_calcw(:,:)
    integer,          pointer :: nb15_calc_list1(:,:), nb15_calc_list(:,:)
    integer,          pointer :: virial_check(:,:)


    cell_pairlist   => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    nsolute         => domain%num_solute
    nwater          => domain%num_water
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6 

    nb15_list       => pairlist%nb15_list
    nb15_listw      => pairlist%nb15_listw
    nb15_cell       => pairlist%nb15_cell
    nb15_cellw      => pairlist%nb15_cellw
    num_nb15_calc1  => pairlist%num_nb15_calc1
    num_nb15_calc   => pairlist%num_nb15_calc 
    num_nb15_calcw  => pairlist%num_nb15_calcw
    nb15_calc_list1 => pairlist%nb15_calc_list1
    nb15_calc_list  => pairlist%nb15_calc_list

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff

    ! calculate energy and gradient
    !

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, ix, num_nb15, rtmp, qtmp, ini_nb15, fin_nb15,         &
    !$omp         k, iy, rij2, L, R, term_lj12, term_lj6, grad_coef,           &
    !$omp         work, term_elec, iwater, ij, j, trans_x,                     &
    !$omp         trans_y, trans_z, iix, list, num_count, j_list, dij_list,    &
    !$omp         force_local, force_localj, lj12, lj6, L1, iatmcls, jatmcls,  &
    !$omp         ieps, jeps, eps, irmin, jrmin, rmin, dij, check_virial)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell, nthread
      do ix = 1, natom(i)
        coord_pbc(1,ix,i) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(2,ix,i) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(3,ix,i) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do

    !$omp barrier

    do i = id+1, ncell_local, nthread
      num_nb15 = 0

      do ix = 1, nsolute(i) - 1

        rtmp(1:3) = coord_pbc(1:3,ix,i)
        qtmp = charge(ix,i)
        ini_nb15 = num_nb15 + 1
        fin_nb15 = num_nb15 + num_nb15_calc1(ix,i)
        num_nb15 = fin_nb15

        num_count = 0
        force_local(1:3) = 0.0_wp

!ocl norecurrence(force)
!ocl swp
        do k = ini_nb15, fin_nb15

          iy   = nb15_calc_list1(k,i)
          dij(1) = rtmp(1) - coord_pbc(1,iy,i)
          dij(2) = rtmp(2) - coord_pbc(2,iy,i)
          dij(3) = rtmp(3) - coord_pbc(3,iy,i)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1,num_count) = dij(1)
            dij_list(2,num_count) = dij(2)
            dij_list(3,num_count) = dij(3)
            dij_list(4,num_count) = rij2
            j_list(num_count) = iy
          end if

        end do

        do k = 1, num_count

          iy = j_list(k)
          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)
          rij2  = cutoff2 * density / dij_list(4,k)
          lj12 = nonb_lj12(atmcls(ix,i),atmcls(iy,i))
          lj6  = nonb_lj6(atmcls(ix,i),atmcls(iy,i))

          L    = int(rij2)
          R    = rij2 - L

          L1 = 3*L - 2
          term_lj12 = table_grad(L1)   + R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6  &
                    + qtmp*charge(iy,i)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force_localj(1,k) = work(1)
          force_localj(2,k) = work(2)
          force_localj(3,k) = work(3)

        end do

        do k = 1, num_count
          iy = j_list(k)
          force(1:3,iy,i,id+1) = force(1:3,iy,i,id+1) + force_localj(1:3,k)
        end do
        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
      end do

      do ix = 1, nsolute(i) 

        rtmp(1:3) = coord_pbc(1:3,ix,i)
        qtmp = charge(ix,i)

        num_count = 0
        force_local(1:3) = 0.0_wp
!ocl norecurrence(force)
!ocl swp
        do iy = nsolute(i)+1, natom(i)

          ! compute distance
          !
          dij(1) = rtmp(1) - coord_pbc(1,iy,i)
          dij(2) = rtmp(2) - coord_pbc(2,iy,i)
          dij(3) = rtmp(3) - coord_pbc(3,iy,i)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1,num_count) = dij(1)
            dij_list(2,num_count) = dij(2)
            dij_list(3,num_count) = dij(3)
            dij_list(4,num_count) = rij2
            j_list(num_count) = iy
          end if
        end do

        do k = 1, num_count
          iy   = j_list(k)
          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)
          rij2  = cutoff2 * density / dij_list(4,k)
          lj12  = nonb_lj12(atmcls(ix,i),atmcls(iy,i))
          lj6   = nonb_lj6(atmcls(ix,i),atmcls(iy,i))

          L     = int(rij2)
          R     = rij2 - L
          L1    = 3*L - 2
          term_lj12 = table_grad(L1)   + R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6  &
                    + qtmp*charge(iy,i)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force_localj(1,k) = work(1)
          force_localj(2,k) = work(2)
          force_localj(3,k) = work(3)

        end do 

        do k = 1, num_count
          iy = j_list(k)
          force(1:3,iy,i,id+1) = force(1:3,iy,i,id+1) + force_localj(1:3,k)
        end do
        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
      end do
    end do

    ! interaction between different cells
    !
    do ij = id+1, maxcell, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)

      num_nb15 = 0

      trans_x = cell_move(1,j,i) * system_size(1)
      trans_y = cell_move(2,j,i) * system_size(2)
      trans_z = cell_move(3,j,i) * system_size(3)

      check_virial = virial_check(j,i)

      do iix = 1, nb15_cell(ij)

        ix = nb15_list(iix,ij)
        rtmp(1:3) = coord_pbc(1:3,ix,i)
        qtmp = charge(ix,i)

        ini_nb15  = num_nb15 + 1
        fin_nb15  = num_nb15 + num_nb15_calc(ix,ij)
        num_nb15  = fin_nb15

        num_count = 0
        force_local(1:3) = 0.0_wp

!ocl norecurrence(force)
!ocl swp
        do k = ini_nb15, fin_nb15
          iy = nb15_calc_list(k,ij)

          dij(1) = rtmp(1) - coord_pbc(1,iy,j) + trans_x
          dij(2) = rtmp(2) - coord_pbc(2,iy,j) + trans_y
          dij(3) = rtmp(3) - coord_pbc(3,iy,j) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1,num_count) = dij(1)
            dij_list(2,num_count) = dij(2)
            dij_list(3,num_count) = dij(3)
            dij_list(4,num_count) = rij2
            j_list(num_count)   = iy
          end if
        end do

        do k = 1, num_count

          iy   = j_list(k)
          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)
          rij2  = cutoff2 * density / dij_list(4,k)
          lj12 = nonb_lj12(atmcls(ix,i),atmcls(iy,j))
          lj6  = nonb_lj6(atmcls(ix,i),atmcls(iy,j))
          L    = int(rij2)
          R    = rij2 - L
          L1   = 3*L - 2
          term_lj12 = table_grad(L1)   + R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6  &
                    + qtmp*charge(iy,j)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force_localj(1,k) = work(1)
          force_localj(2,k) = work(2)
          force_localj(3,k) = work(3)

        end do

        if (check_virial == 1) &
          virial(1:3,ij) = virial(1:3,ij) - force_local(1:3)
        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
        do k = 1, num_count
          iy = j_list(k)
          force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + force_localj(1:3,k)
        end do

      end do

      do iix = 1, nb15_cellw(ij)

        iwater = nb15_listw(iix,ij)

        ini_nb15 = num_nb15 + 1
        fin_nb15 = num_nb15 + num_nb15_calcw(iwater,ij)
        num_nb15 = fin_nb15

        do list = 1, 3

          ix   = nsolute(i) + 3*(iwater-1) + list
          rtmp(1:3) = coord_pbc(1:3,ix,i)
          qtmp = charge(ix,i)
          num_count = 0
          force_local(1:3) = 0.0_wp
!ocl norecurrence(force)
!ocl swp
          do k = ini_nb15, fin_nb15
            iy = nb15_calc_list(k,ij)

            ! compute distance
            !
            dij(1) = rtmp(1) - coord_pbc(1,iy,j) + trans_x
            dij(2) = rtmp(2) - coord_pbc(2,iy,j) + trans_y
            dij(3) = rtmp(3) - coord_pbc(3,iy,j) + trans_z
            rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            if (rij2 < cutoff2) then
              num_count = num_count + 1
              dij_list(1,num_count) = dij(1)
              dij_list(2,num_count) = dij(2)
              dij_list(3,num_count) = dij(3)
              dij_list(4,num_count) = rij2
              j_list(num_count) = iy
            end if
          end do

          do k = 1, num_count
            iy = j_list(k)
            dij(1) = dij_list(1,k)
            dij(2) = dij_list(2,k)
            dij(3) = dij_list(3,k)
            iy = j_list(k)
            rij2 = cutoff2 * density / dij_list(4,k)
            lj12 = nonb_lj12(atmcls(ix,i),atmcls(iy,j))
            lj6  = nonb_lj6(atmcls(ix,i),atmcls(iy,j))

            L    = int(rij2)
            R    = rij2 - L
            L1   = 3*L - 2
            term_lj12 = table_grad(L1)   + R*(table_grad(L1+3)-table_grad(L1))
            term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
            term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
            grad_coef = term_lj12*lj12 - term_lj6*lj6   &
                      + qtmp*charge(iy,j)*term_elec

            work(1) = grad_coef*dij(1)
            work(2) = grad_coef*dij(2)
            work(3) = grad_coef*dij(3)

            ! store force
            !
            force_local(1) = force_local(1) - work(1)
            force_local(2) = force_local(2) - work(2)
            force_local(3) = force_local(3) - work(3)
            force_localj(1,k) = work(1)
            force_localj(2,k) = work(2)
            force_localj(3,k) = work(3)

          end do

          do k = 1, num_count
            iy = j_list(k)
            force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + force_localj(1:3,k)
          end do
          if (check_virial == 1) &
            virial(1:3,ij) = virial(1:3,ij) - force_local(1:3)
          force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)

        end do

      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_force_nonbond_table_solute_linear

#else

  subroutine compute_force_nonbond_table_solute_linear( &
                                                 domain, enefunc, &
                                                 pairlist, coord_pbc, &
                                                 force, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:)

    ! local variables
    real(wp)                  :: dij(1:3), rij2, inv_r2
    real(wp)                  :: R, lj12, lj6
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: trans_x, trans_y, trans_z
    real(wp)                  :: rtmp(1:3), qtmp
    real(wp)                  :: force_local(3)
    real(wp)                  :: force_localj(3,MaxAtom)
    integer                   :: i, ix, iy, j, k, ij, iix, iwater, list, L, L1
    integer                   :: num_nb15, ini_nb15, fin_nb15, num_count
    integer                   :: id, omp_get_thread_num
    integer                   :: ncell, ncell_local
    integer                   :: iatmcls,jatmcls
    integer                   :: check_virial

    real(dp),         pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: cell_move(:,:,:), system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: density, cutoff
    real(wp),         pointer :: table_grad(:)
    integer,          pointer :: cell_pairlist(:,:), atmcls(:,:)
    integer,          pointer :: natom(:), nsolute(:), nwater(:)
    integer,          pointer :: nb15_list(:,:), nb15_listw(:,:)
    integer,          pointer :: nb15_cell(:), nb15_cellw(:)
    integer,          pointer :: num_nb15_calc1(:,:), num_nb15_calc(:,:)
    integer,          pointer :: num_nb15_calcw(:,:)
    integer,          pointer :: nb15_calc_list1(:,:), nb15_calc_list(:,:)
    integer,          pointer :: virial_check(:,:)


    cell_pairlist   => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    nsolute         => domain%num_solute
    nwater          => domain%num_water
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6

    nb15_list       => pairlist%nb15_list
    nb15_listw      => pairlist%nb15_listw
    nb15_cell       => pairlist%nb15_cell
    nb15_cellw      => pairlist%nb15_cellw
    num_nb15_calc1  => pairlist%num_nb15_calc1
    num_nb15_calc   => pairlist%num_nb15_calc
    num_nb15_calcw  => pairlist%num_nb15_calcw
    nb15_calc_list1 => pairlist%nb15_calc_list1
    nb15_calc_list  => pairlist%nb15_calc_list

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff

    ! calculate energy and gradient
    !

    !$omp parallel default(shared)                                     &
    !$omp private(id, i, ix, num_nb15, rtmp, qtmp, ini_nb15, fin_nb15, &
    !$omp         k, iy, rij2, L, R, term_lj12, term_lj6, grad_coef,   &
    !$omp         work, term_elec, iwater, ij, j, trans_x,             &
    !$omp         trans_y, trans_z, iix, list, num_count,              &
    !$omp         force_local, force_localj, lj12, lj6, L1,            &
    !$omp         iatmcls, jatmcls, dij, inv_r2, check_virial)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell, nthread
      do ix = 1, natom(i)
        coord_pbc(1,ix,i) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(2,ix,i) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(3,ix,i) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do

    !$omp barrier

    do i = id+1, ncell_local, nthread
      num_nb15 = 0

      do ix = 1, nsolute(i) - 1

        rtmp(1:3) = coord_pbc(1:3,ix,i)
        qtmp = charge(ix,i)
        iatmcls = atmcls(ix,i)

        ini_nb15 = num_nb15 + 1
        fin_nb15 = num_nb15 + num_nb15_calc1(ix,i)
        num_nb15 = fin_nb15

        force_local(1:3) = 0.0_wp
!ocl norecurrence(force)
!ocl swp
        do k = ini_nb15, fin_nb15

          iy   = nb15_calc_list1(k,i)
          dij(1) = rtmp(1) - coord_pbc(1,iy,i)
          dij(2) = rtmp(2) - coord_pbc(2,iy,i)
          dij(3) = rtmp(3) - coord_pbc(3,iy,i)
          rij2   = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          inv_r2 = 1.0_wp / rij2

          rij2   = cutoff2 * density * inv_r2

          jatmcls = atmcls(iy,i)
          lj12   = nonb_lj12(iatmcls,jatmcls)
          lj6    = nonb_lj6(iatmcls,jatmcls)

          L     = int(rij2)
          R     = rij2 - L

          L1    = 3*L - 2
          term_lj12 = table_grad(L1  ) + R*(table_grad(L1+3)-table_grad(L1  ))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6   &
                    + qtmp*charge(iy,i)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force(1,iy,i,id+1) = force(1,iy,i,id+1) + work(1)
          force(2,iy,i,id+1) = force(2,iy,i,id+1) + work(2)
          force(3,iy,i,id+1) = force(3,iy,i,id+1) + work(3)

        end do
        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
      end do

      do ix = 1, nsolute(i)

        rtmp(1:3) = coord_pbc(1:3,ix,i)
        qtmp = charge(ix,i)
        iatmcls = atmcls(ix,i)

        force_local(1:3) = 0.0_wp
!ocl norecurrence(force)
!ocl swp
        do iy = nsolute(i)+1, natom(i)

          ! compute distance
          !
          dij(1) = rtmp(1) - coord_pbc(1,iy,i)
          dij(2) = rtmp(2) - coord_pbc(2,iy,i)
          dij(3) = rtmp(3) - coord_pbc(3,iy,i)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          inv_r2 = 1.0_wp / rij2
          rij2  = cutoff2 * density *inv_r2

          jatmcls = atmcls(iy,i)
          lj12    = nonb_lj12(iatmcls,jatmcls)
          lj6     = nonb_lj6(iatmcls,jatmcls)

          L     = int(rij2)
          R     = rij2 - L
          L1    = 3*L - 2
          term_lj12 = table_grad(L1  ) + R*(table_grad(L1+3)-table_grad(L1  ))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6   &
                    + qtmp*charge(iy,i)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force(1,iy,i,id+1) = force(1,iy,i,id+1) + work(1)
          force(2,iy,i,id+1) = force(2,iy,i,id+1) + work(2)
          force(3,iy,i,id+1) = force(3,iy,i,id+1) + work(3)

        end do

        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
      end do
    end do

    ! interaction between different cells
    !
    do ij = id+1, maxcell, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)

      num_nb15 = 0

      trans_x = cell_move(1,j,i) * system_size(1)
      trans_y = cell_move(2,j,i) * system_size(2)
      trans_z = cell_move(3,j,i) * system_size(3)

      check_virial = virial_check(j,i)

      do iix = 1, nb15_cell(ij)

        ix = nb15_list(iix,ij)
        rtmp(1:3) = coord_pbc(1:3,ix,i)
        qtmp = charge(ix,i)
        iatmcls = atmcls(ix,i)

        ini_nb15  = num_nb15 + 1
        fin_nb15  = num_nb15 + num_nb15_calc(ix,ij)
        num_nb15  = fin_nb15

        num_count = 0
        force_local(1:3) = 0.0_wp

!ocl norecurrence(force)
!ocl swp
        do k = ini_nb15, fin_nb15
          iy = nb15_calc_list(k,ij)

          dij(1) = rtmp(1) - coord_pbc(1,iy,j) + trans_x
          dij(2) = rtmp(2) - coord_pbc(2,iy,j) + trans_y
          dij(3) = rtmp(3) - coord_pbc(3,iy,j) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          inv_r2 = 1.0_wp / rij2
          rij2 = cutoff2 * density * inv_r2

          jatmcls = atmcls(iy,j)
          lj12    = nonb_lj12(iatmcls,jatmcls)
          lj6     = nonb_lj6(iatmcls,jatmcls)

          L    = int(rij2)
          R    = rij2 - L
          L1   = 3*L - 2
          term_lj12 = table_grad(L1  ) + R*(table_grad(L1+3)-table_grad(L1  ))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6   &
                    + qtmp*charge(iy,j)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force(1,iy,j,id+1) = force(1,iy,j,id+1) + work(1)
          force(2,iy,j,id+1) = force(2,iy,j,id+1) + work(2)
          force(3,iy,j,id+1) = force(3,iy,j,id+1) + work(3)

        end do

        if (check_virial == 1) &
          virial(1:3,ij) = virial(1:3,ij) - force_local(1:3)
        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)

      end do

      do iix = 1, nb15_cellw(ij)

        iwater = nb15_listw(iix,ij)

        ini_nb15 = num_nb15 + 1
        fin_nb15 = num_nb15 + num_nb15_calcw(iwater,ij)
        num_nb15 = fin_nb15

        do list = 1, 3

          ix   = nsolute(i) + 3*(iwater-1) + list
          rtmp(1:3) = coord_pbc(1:3,ix,i)
          qtmp = charge(ix,i)
          iatmcls = atmcls(ix,i)

          force_local(1:3) = 0.0_wp
!ocl norecurrence(force)
!ocl swp
          do k = ini_nb15, fin_nb15
            iy = nb15_calc_list(k,ij)

            ! compute distance
            !
            dij(1) = rtmp(1) - coord_pbc(1,iy,j) + trans_x
            dij(2) = rtmp(2) - coord_pbc(2,iy,j) + trans_y
            dij(3) = rtmp(3) - coord_pbc(3,iy,j) + trans_z
            rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            inv_r2 = 1.0_wp / rij2
            rij2 = cutoff2 * density *inv_r2
            jatmcls = atmcls(iy,j)
            lj12    = nonb_lj12(iatmcls,jatmcls)
            lj6     = nonb_lj6(iatmcls,jatmcls)

            L    = int(rij2)
            R    = rij2 - L
            L1   = 3*L - 2
            term_lj12 = table_grad(L1  ) + R*(table_grad(L1+3)-table_grad(L1  ))
            term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
            term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
            grad_coef = term_lj12*lj12 - term_lj6*lj6   &
                      + qtmp*charge(iy,j)*term_elec

            work(1) = grad_coef*dij(1)
            work(2) = grad_coef*dij(2)
            work(3) = grad_coef*dij(3)

            ! store force
            !
            force_local(1) = force_local(1) - work(1)
            force_local(2) = force_local(2) - work(2)
            force_local(3) = force_local(3) - work(3)
            force(1,iy,j,id+1) = force(1,iy,j,id+1) + work(1)
            force(2,iy,j,id+1) = force(2,iy,j,id+1) + work(2)
            force(3,iy,j,id+1) = force(3,iy,j,id+1) + work(3)

          end do

          if (check_virial == 1) &
            virial(1:3,ij) = virial(1:3,ij) - force_local(1:3)
          force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)

        end do

      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_force_nonbond_table_solute_linear

#endif

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_table_linear
  !> @brief        calculate nonbonded energy with lookup table
  !! @authors      JJ 
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] virial   : virial term of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : van der Waals energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_table_linear( &
                                                 domain, enefunc, pairlist, &
                                                 coord_pbc, force, virial,  &
                                                 eelec, evdw)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:)
    real(dp),                 intent(inout) :: eelec(nthread)
    real(dp),                 intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: lj6, lj12
    real(wp)                  :: R
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, term_lj, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: trans_x, trans_y, trans_z
    real(wp)                  :: rtmp(1:3), qtmp
    real(wp)                  :: elec_temp, evdw_temp
    real(wp)                  :: dij_list(3,MaxAtom)
    real(wp)                  :: force_local(3), rij2_list(MaxAtom)
    real(wp)                  :: force_localj(3,MaxAtom)
    integer                   :: i, ix, iy, j, k, ij, iix, iwater, list, L, L1
    integer                   :: num_nb15, ini_nb15, fin_nb15, num_count
    integer                   :: id, omp_get_thread_num, j_list(MaxAtom)
    integer                   :: ncell, ncell_local
    integer                   :: check_virial

    real(dp),         pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: cell_move(:,:,:), system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: density, cutoff
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer,          pointer :: cell_pairlist(:,:), atmcls(:,:)
    integer,          pointer :: natom(:), nsolute(:), nwater(:)
    integer,          pointer :: nb15_list(:,:)
    integer,          pointer :: nb15_cell(:)
    integer,          pointer :: num_nb15_calc1(:,:), num_nb15_calc(:,:)
    integer,          pointer :: nb15_calc_list1(:,:), nb15_calc_list(:,:)
    integer,          pointer :: virial_check(:,:)


    cell_pairlist   => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    nsolute         => domain%num_solute
    nwater          => domain%num_water
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6

    nb15_list       => pairlist%nb15_list
    nb15_cell       => pairlist%nb15_cell
    num_nb15_calc1  => pairlist%num_nb15_calc1
    num_nb15_calc   => pairlist%num_nb15_calc
    nb15_calc_list1 => pairlist%nb15_calc_list1
    nb15_calc_list  => pairlist%nb15_calc_list

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff

    ! calculate energy and gradient
    !

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, ix, num_nb15, rtmp, qtmp, ini_nb15, fin_nb15, k, iy,  &
    !$omp         rij2, L, R, term_lj12, term_lj6, term_lj, grad_coef, work,   &
    !$omp         term_elec, iwater, ij, j, trans_x, trans_y,                  &
    !$omp         trans_z, iix, list, num_count, j_list, dij_list, rij2_list,  &
    !$omp         force_local, force_localj, lj6, lj12, L1, elec_temp, &
    !$omp         evdw_temp, dij, check_virial)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell, nthread
      do ix = 1, natom(i)
        coord_pbc(1,ix,i) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(2,ix,i) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(3,ix,i) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do

    !$omp barrier

    ! energy within a cell
    !
    do i = id+1, ncell_local, nthread

      num_nb15 = 0
      do ix = 1, natom(i) - 1

        rtmp(1:3) = coord_pbc(1:3,ix,i)
        qtmp = charge(ix,i)
        ini_nb15  = num_nb15 + 1
        fin_nb15  = num_nb15 + num_nb15_calc1(ix,i)
        num_nb15  = fin_nb15
        num_count = 0

        do k = ini_nb15, fin_nb15

          iy = nb15_calc_list1(k,i)

          ! compute distance
          !
          dij(1) = rtmp(1) - coord_pbc(1,iy,i)
          dij(2) = rtmp(2) - coord_pbc(2,iy,i)
          dij(3) = rtmp(3) - coord_pbc(3,iy,i)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1:3,num_count)  = dij(1:3)
            rij2_list(num_count) = rij2
            j_list(num_count)   = iy
          end if
        end do

        force_local(1:3) = 0.0_wp
        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp

!ocl norecurrence(force)
!ocl swp
        do k = 1, num_count
          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)
          iy   = j_list(k)
          rij2  = cutoff2*density/rij2_list(k)
          lj6  = nonb_lj6(atmcls(ix,i),atmcls(iy,i))
          lj12 = nonb_lj12(atmcls(ix,i),atmcls(iy,i))

          L    = int(rij2)
          R    = rij2 - L

          L1   = 3*L - 2
          term_lj12 = table_ene(L1) + R*(table_ene(L1+3)-table_ene(L1))
          term_lj6  = table_ene(L1+1) + R*(table_ene(L1+4)-table_ene(L1+1))
          term_elec = table_ene(L1+2) + R*(table_ene(L1+5)-table_ene(L1+2))
          evdw_temp = evdw_temp + term_lj12*lj12 - term_lj6*lj6
          elec_temp = elec_temp + qtmp*charge(iy,i)*term_elec

          term_lj12 = table_grad(L1) + R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6                           &
                     + qtmp*charge(iy,i)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force_localj(1,k)  = work(1)
          force_localj(2,k)  = work(2)
          force_localj(3,k)  = work(3)

        end do

        do k = 1, num_count
          iy = j_list(k)
          force(1:3,iy,i,id+1) = force(1:3,iy,i,id+1) + force_localj(1:3,k)
        end do
        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw(id+1) = evdw(id+1) + evdw_temp

      end do

    end do


    ! interaction between different cells
    !
    do ij = id+1, maxcell, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)

      num_nb15 = 0

      trans_x = cell_move(1,j,i) * system_size(1)
      trans_y = cell_move(2,j,i) * system_size(2)
      trans_z = cell_move(3,j,i) * system_size(3)
      check_virial = virial_check(j,i)

      do iix = 1, nb15_cell(ij)

        ix   = nb15_list(iix,ij)
        rtmp(1:3) = coord_pbc(1:3,ix,i)
        qtmp = charge(ix,i)

        ini_nb15 = num_nb15 + 1
        fin_nb15 = num_nb15 + num_nb15_calc(ix,ij)
        num_nb15 = fin_nb15

        num_count = 0
        do k = ini_nb15, fin_nb15
          iy   = nb15_calc_list(k,ij)
          dij(1) = rtmp(1) - coord_pbc(1,iy,j) + trans_x
          dij(2) = rtmp(2) - coord_pbc(2,iy,j) + trans_y
          dij(3) = rtmp(3) - coord_pbc(3,iy,j) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1:3,num_count)  = dij(1:3)
            rij2_list(num_count) = rij2
            j_list(num_count)   = iy
          end if
        end do

        force_local(1:3) = 0.0_wp
        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp

!ocl norecurrence(force)
!ocl swp
        do k = 1, num_count
          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)
          iy   = j_list(k)
          rij2  = cutoff2*density/rij2_list(k)
          lj6  = nonb_lj6 (atmcls(ix,i),atmcls(iy,j))
          lj12 = nonb_lj12(atmcls(ix,i),atmcls(iy,j))

          L    = int(rij2)
          R    = rij2 - L

          L1   = 3*L - 2
          term_lj12 = table_ene(L1) + R*(table_ene(L1+3)-table_ene(L1))
          term_lj6  = table_ene(L1+1) + R*(table_ene(L1+4)-table_ene(L1+1))
          term_elec = table_ene(L1+2) + R*(table_ene(L1+5)-table_ene(L1+2))
          evdw_temp = evdw_temp + term_lj12*lj12 - term_lj6*lj6
          elec_temp = elec_temp + qtmp*charge(iy,j)*term_elec

          term_lj12 = table_grad(L1) + R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6                           &
                     + qtmp*charge(iy,j)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)

          force_localj(1,k) = work(1)
          force_localj(2,k) = work(2)
          force_localj(3,k) = work(3)

        end do

        do k = 1, num_count
          iy = j_list(k)
          force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + force_localj(1:3,k)
        end do
        if (check_virial == 1) &
          virial(1:3,ij) = virial(1:3,ij) - force_local(1:3)
        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw(id+1) = evdw(id+1) + evdw_temp

      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond_table_linear

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_table_linear_check
  !> @brief        calculate nonbonded energy with lookup table
  !! @authors      JJ , CK
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] virial   : virial term of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : van der Waals energy of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_table_linear_check( &
                                                 domain, enefunc, pairlist, &
                                                 coord_pbc, force, virial,  &
                                                 eelec, evdw)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:)
    real(dp),                 intent(inout) :: eelec(nthread)
    real(dp),                 intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: lj6, lj12
    real(wp)                  :: R
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, term_lj, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: trans_x, trans_y, trans_z
    real(wp)                  :: rtmp(1:3), qtmp
    real(wp)                  :: elec_temp, evdw_temp
    real(wp)                  :: dij_list(3,MaxAtom)
    real(wp)                  :: force_local(3), rij2_list(MaxAtom)
    real(wp)                  :: force_localj(3,MaxAtom)
    real(wp)                  :: minimum_contact
    integer                   :: i, ix, iy, j, k, ij, iix, iwater, list, L, L1
    integer                   :: num_nb15, ini_nb15, fin_nb15, num_count
    integer                   :: id, omp_get_thread_num, j_list(MaxAtom)
    integer                   :: ncell, ncell_local
    integer                   :: check_virial

    real(dp),         pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: cell_move(:,:,:), system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: density, cutoff
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer,          pointer :: cell_pairlist(:,:), atmcls(:,:)
    integer,          pointer :: natom(:), nsolute(:), nwater(:)
    integer,          pointer :: nb15_list(:,:)
    integer,          pointer :: nb15_cell(:)
    integer,          pointer :: num_nb15_calc1(:,:), num_nb15_calc(:,:)
    integer,          pointer :: nb15_calc_list1(:,:), nb15_calc_list(:,:)
    integer,          pointer :: virial_check(:,:)


    cell_pairlist   => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    nsolute         => domain%num_solute
    nwater          => domain%num_water
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6

    nb15_list       => pairlist%nb15_list
    nb15_cell       => pairlist%nb15_cell
    num_nb15_calc1  => pairlist%num_nb15_calc1
    num_nb15_calc   => pairlist%num_nb15_calc
    nb15_calc_list1 => pairlist%nb15_calc_list1
    nb15_calc_list  => pairlist%nb15_calc_list

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff
    minimum_contact =  enefunc%minimum_contact

    ! calculate energy and gradient
    !

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, ix, num_nb15, rtmp, qtmp, ini_nb15, fin_nb15, k, iy,  &
    !$omp         rij2, L, R, term_lj12, term_lj6, term_lj, grad_coef, work,   &
    !$omp         term_elec, iwater, ij, j, trans_x, trans_y,                  &
    !$omp         trans_z, iix, list, num_count, j_list, dij_list, rij2_list,  &
    !$omp         force_local, force_localj, lj6, lj12, L1, elec_temp, &
    !$omp         evdw_temp, dij, check_virial)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell, nthread
      do ix = 1, natom(i)
        coord_pbc(1,ix,i) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(2,ix,i) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(3,ix,i) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do

    !$omp barrier

    ! energy within a cell
    !
    do i = id+1, ncell_local, nthread

      num_nb15 = 0
      do ix = 1, natom(i) - 1

        rtmp(1:3) = coord_pbc(1:3,ix,i)
        qtmp = charge(ix,i)
        ini_nb15  = num_nb15 + 1
        fin_nb15  = num_nb15 + num_nb15_calc1(ix,i)
        num_nb15  = fin_nb15
        num_count = 0

        do k = ini_nb15, fin_nb15

          iy = nb15_calc_list1(k,i)

          ! compute distance
          !
          dij(1) = rtmp(1) - coord_pbc(1,iy,i)
          dij(2) = rtmp(2) - coord_pbc(2,iy,i)
          dij(3) = rtmp(3) - coord_pbc(3,iy,i)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1:3,num_count)  = dij(1:3)
            rij2_list(num_count) = max(rij2, minimum_contact)
            j_list(num_count)   = iy
          end if
        end do

        force_local(1:3) = 0.0_wp
        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp

!ocl norecurrence(force)
!ocl swp
        do k = 1, num_count
          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)
          iy   = j_list(k)
          rij2  = cutoff2*density/rij2_list(k)
          lj6  = nonb_lj6(atmcls(ix,i),atmcls(iy,i))
          lj12 = nonb_lj12(atmcls(ix,i),atmcls(iy,i))

          L    = int(rij2)
          R    = rij2 - L

          L1   = 3*L - 2
          term_lj12 = table_ene(L1) + R*(table_ene(L1+3)-table_ene(L1))
          term_lj6  = table_ene(L1+1) + R*(table_ene(L1+4)-table_ene(L1+1))
          term_elec = table_ene(L1+2) + R*(table_ene(L1+5)-table_ene(L1+2))
          evdw_temp = evdw_temp + term_lj12*lj12 - term_lj6*lj6
          elec_temp = elec_temp + qtmp*charge(iy,i)*term_elec

          term_lj12 = table_grad(L1) + R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6                           &
                     + qtmp*charge(iy,i)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force_localj(1,k)  = work(1)
          force_localj(2,k)  = work(2)
          force_localj(3,k)  = work(3)

        end do

        do k = 1, num_count
          iy = j_list(k)
          force(1:3,iy,i,id+1) = force(1:3,iy,i,id+1) + force_localj(1:3,k)
        end do
        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw(id+1) = evdw(id+1) + evdw_temp

      end do

    end do


    ! interaction between different cells
    !
    do ij = id+1, maxcell, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)

      num_nb15 = 0

      trans_x = cell_move(1,j,i) * system_size(1)
      trans_y = cell_move(2,j,i) * system_size(2)
      trans_z = cell_move(3,j,i) * system_size(3)
      check_virial = virial_check(j,i)

      do iix = 1, nb15_cell(ij)

        ix   = nb15_list(iix,ij)
        rtmp(1:3) = coord_pbc(1:3,ix,i)
        qtmp = charge(ix,i)

        ini_nb15 = num_nb15 + 1
        fin_nb15 = num_nb15 + num_nb15_calc(ix,ij)
        num_nb15 = fin_nb15

        num_count = 0
        do k = ini_nb15, fin_nb15
          iy   = nb15_calc_list(k,ij)
          dij(1) = rtmp(1) - coord_pbc(1,iy,j) + trans_x
          dij(2) = rtmp(2) - coord_pbc(2,iy,j) + trans_y
          dij(3) = rtmp(3) - coord_pbc(3,iy,j) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1:3,num_count)  = dij(1:3)
            rij2_list(num_count) = max(rij2, minimum_contact)
            j_list(num_count)   = iy
          end if
        end do

        force_local(1:3) = 0.0_wp
        elec_temp = 0.0_wp
        evdw_temp = 0.0_wp

!ocl norecurrence(force)
!ocl swp
        do k = 1, num_count
          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)
          iy   = j_list(k)
          rij2  = cutoff2*density/rij2_list(k)
          lj6  = nonb_lj6 (atmcls(ix,i),atmcls(iy,j))
          lj12 = nonb_lj12(atmcls(ix,i),atmcls(iy,j))

          L    = int(rij2)
          R    = rij2 - L

          L1   = 3*L - 2
          term_lj12 = table_ene(L1) + R*(table_ene(L1+3)-table_ene(L1))
          term_lj6  = table_ene(L1+1) + R*(table_ene(L1+4)-table_ene(L1+1))
          term_elec = table_ene(L1+2) + R*(table_ene(L1+5)-table_ene(L1+2))
          evdw_temp = evdw_temp + term_lj12*lj12 - term_lj6*lj6
          elec_temp = elec_temp + qtmp*charge(iy,j)*term_elec

          term_lj12 = table_grad(L1) + R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6                           &
                     + qtmp*charge(iy,j)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)

          force_localj(1,k) = work(1)
          force_localj(2,k) = work(2)
          force_localj(3,k) = work(3)

        end do

        do k = 1, num_count
          iy = j_list(k)
          force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + force_localj(1:3,k)
        end do
        if (check_virial == 1) &
          virial(1:3,ij) = virial(1:3,ij) - force_local(1:3)
        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw(id+1) = evdw(id+1) + evdw_temp

      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond_table_linear_check

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_table_linear
  !> @brief        calculate nonbonded force without solvents with lookup table
  !! @authors      JJ 
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] virial   : virial term of target systems
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

#ifndef KCOMP

  subroutine compute_force_nonbond_table_linear( &
                                                 domain, enefunc, pairlist, &
                                                 coord_pbc, force, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:)

    ! local variables
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: R, lj12, lj6
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: trans_x, trans_y, trans_z
    real(wp)                  :: rtmp(1:3), qtmp
    real(wp)                  :: force_local(3), force_local_iy(3,MaxAtom)
    real(wp)                  :: dij_list(4,MaxAtom)
    real(wp)                  :: ieps, jeps, eps, irmin, jrmin, rmin
    integer                   :: i, ix, iy, j, k, ij, iix, iwater, list, L, L1
    integer                   :: num_nb15, ini_nb15, fin_nb15, num_count
    integer                   :: id, omp_get_thread_num, j_list(MaxAtom)
    integer                   :: iatmcls,jatmcls
    integer                   :: ncell, ncell_local
    integer                   :: check_virial

    real(dp),         pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: cell_move(:,:,:), system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: density, cutoff
    real(wp),         pointer :: table_grad(:)
    integer,          pointer :: cell_pairlist(:,:), atmcls(:,:)
    integer,          pointer :: natom(:), nsolute(:), nwater(:)
    integer,          pointer :: nb15_list(:,:)
    integer,          pointer :: nb15_cell(:)
    integer,          pointer :: num_nb15_calc1(:,:), num_nb15_calc(:,:)
    integer,          pointer :: nb15_calc_list1(:,:), nb15_calc_list(:,:)
    integer,          pointer :: virial_check(:,:)


    cell_pairlist   => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    nsolute         => domain%num_solute
    nwater          => domain%num_water
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6 

    nb15_list       => pairlist%nb15_list
    nb15_cell       => pairlist%nb15_cell
    num_nb15_calc1  => pairlist%num_nb15_calc1
    num_nb15_calc   => pairlist%num_nb15_calc 
    nb15_calc_list1 => pairlist%nb15_calc_list1
    nb15_calc_list  => pairlist%nb15_calc_list

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff

    ! calculate energy and gradient
    !

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, ix, num_nb15, rtmp, qtmp, ini_nb15, fin_nb15, k, iy,  &
    !$omp         rij2, L, R, term_lj12, term_lj6, grad_coef, work, term_elec, &
    !$omp         iwater, ij, j, trans_x, trans_y, trans_z,                    &
    !$omp         iix, list, num_count, j_list, force_local, force_local_iy,   &
    !$omp         lj12, lj6, L1, dij_list, iatmcls, jatmcls, ieps, jeps, eps,  &
    !$omp         irmin, jrmin, rmin, dij, check_virial)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i = id+1, ncell, nthread
      do ix = 1, natom(i)
        coord_pbc(1,ix,i) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(2,ix,i) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(3,ix,i) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do

    !$omp barrier

    ! energy within a cell(solute-solute)
    !
    do i = id+1, ncell_local, nthread

      num_nb15 = 0

      do ix = 1, natom(i) - 1

        rtmp(1:3) = coord_pbc(1:3,ix,i)
        qtmp = charge(ix,i)
        iatmcls = atmcls(ix,i)
        ieps    = lj_coef(1,iatmcls)
        irmin   = lj_coef(2,iatmcls)

        ini_nb15  = num_nb15 + 1
        fin_nb15  = num_nb15 + num_nb15_calc1(ix,i)
        num_nb15  = fin_nb15
        force_local(1:3) = 0.0_wp
        num_count = 0

!ocl norecurrence(force)
!ocl swp
        do k = ini_nb15, fin_nb15

          iy = nb15_calc_list1(k,i)
          dij(1) = rtmp(1) - coord_pbc(1,iy,i)
          dij(2) = rtmp(2) - coord_pbc(2,iy,i)
          dij(3) = rtmp(3) - coord_pbc(3,iy,i)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1,num_count) = dij(1)
            dij_list(2,num_count) = dij(2)
            dij_list(3,num_count) = dij(3)
            dij_list(4,num_count) = rij2
            j_list(num_count) = iy
          end if
        end do
        force_local(1:3) = 0.0_wp

        do k = 1, num_count
          iy   = j_list(k)
          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)
          rij2  = cutoff2 * density / dij_list(4,k)

          lj12 = nonb_lj12(atmcls(ix,i),atmcls(iy,i))
          lj6  = nonb_lj6(atmcls(ix,i),atmcls(iy,i))

          L    = int(rij2)
          R    = rij2 - L

          L1 = 3*L - 2
          term_lj12 = table_grad(L1) + R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6  &
                    + qtmp*charge(iy,i)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force_local_iy(1,k) = work(1)
          force_local_iy(2,k) = work(2)
          force_local_iy(3,k) = work(3)

        end do

        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
        do k = 1, num_count
          iy = j_list(k) 
          force(1:3,iy,i,id+1) = force(1:3,iy,i,id+1) + force_local_iy(1:3,k)
        end do
      end do

    end do

    ! interaction between different cells
    !
    do ij = id+1, maxcell, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)

      num_nb15 = 0

      trans_x = cell_move(1,j,i) * system_size(1)
      trans_y = cell_move(2,j,i) * system_size(2)
      trans_z = cell_move(3,j,i) * system_size(3)
      check_virial = virial_check(j,i)

      do iix = 1, nb15_cell(ij)

        ix = nb15_list(iix,ij)
        rtmp(1:3) = coord_pbc(1:3,ix,i)
        qtmp = charge(ix,i)
        iatmcls = atmcls(ix,i)
        ieps    = lj_coef(1,iatmcls)
        irmin   = lj_coef(2,iatmcls)

        ini_nb15  = num_nb15 + 1
        fin_nb15  = num_nb15 + num_nb15_calc(ix,ij)
        num_nb15  = fin_nb15

        force_local(1:3) = 0.0_wp
        num_count = 0

!ocl norecurrence(force)
!ocl swp
        do k = ini_nb15, fin_nb15
          iy = nb15_calc_list(k,ij)

          dij(1) = rtmp(1) - coord_pbc(1,iy,j) + trans_x
          dij(2) = rtmp(2) - coord_pbc(2,iy,j) + trans_y
          dij(3) = rtmp(3) - coord_pbc(3,iy,j) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1,num_count) = dij(1)
            dij_list(2,num_count) = dij(2)
            dij_list(3,num_count) = dij(3)
            dij_list(4,num_count) = rij2
            j_list(num_count) = iy
          end if
        end do

        force_local(1:3) = 0.0_wp
        do k = 1, num_count
          iy = j_list(k)
          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)
          rij2  = cutoff2 * density / dij_list(4,k)
          lj12 = nonb_lj12(atmcls(ix,i),atmcls(iy,j))
          lj6  = nonb_lj6(atmcls(ix,i),atmcls(iy,j))

          L  = int(rij2)
          R  = rij2 - L
          L1 = 3*L - 2
          term_lj12 = table_grad(L1) + R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6   &
                    + qtmp*charge(iy,j)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force_local_iy(1,k) = work(1)
          force_local_iy(2,k) = work(2)
          force_local_iy(3,k) = work(3)

        end do

        do k = 1, num_count
          iy = j_list(k)
          force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + force_local_iy(1:3,k)
        end do
        if (check_virial == 1) &
          virial(1:3,ij) = virial(1:3,ij) - force_local(1:3)
        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_force_nonbond_table_linear

#else

  subroutine compute_force_nonbond_table_linear( &
                                                 domain, enefunc, pairlist, &
                                                 coord_pbc, force, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:)

    ! local variables
    real(wp)                  :: dij(1:3), rij2, inv_r2
    real(wp)                  :: R, lj12, lj6
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: trans_x, trans_y, trans_z
    real(wp)                  :: rtmp(1:3), qtmp
    real(wp)                  :: force_local(3)
    integer                   :: i, ix, iy, j, k, ij, iix, iwater, list, L, L1
    integer                   :: num_nb15, ini_nb15, fin_nb15, num_count
    integer                   :: id, omp_get_thread_num
    integer                   :: iatmcls,jatmcls
    integer                   :: ncell, ncell_local
    integer                   :: check_virial

    real(dp),         pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: cell_move(:,:,:), system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: density, cutoff
    real(wp),         pointer :: table_grad(:)
    integer,          pointer :: cell_pairlist(:,:), atmcls(:,:)
    integer,          pointer :: natom(:), nsolute(:), nwater(:)
    integer,          pointer :: nb15_list(:,:)
    integer,          pointer :: nb15_cell(:)
    integer,          pointer :: num_nb15_calc1(:,:), num_nb15_calc(:,:)
    integer,          pointer :: nb15_calc_list1(:,:), nb15_calc_list(:,:)
    integer,          pointer :: virial_check(:,:)


    cell_pairlist   => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    nsolute         => domain%num_solute
    nwater          => domain%num_water
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6

    nb15_list       => pairlist%nb15_list
    nb15_cell       => pairlist%nb15_cell
    num_nb15_calc1  => pairlist%num_nb15_calc1
    num_nb15_calc   => pairlist%num_nb15_calc
    nb15_calc_list1 => pairlist%nb15_calc_list1
    nb15_calc_list  => pairlist%nb15_calc_list

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff

    ! calculate energy and gradient
    !

    !$omp parallel default(shared)                                      &
    !$omp private(id, i, ix, num_nb15, rtmp, qtmp, ini_nb15, fin_nb15,  & 
    !$omp         k, iy, iwater, ij, j, iix, list, num_count, L, L1,    &
    !$omp         iatmcls, jatmcls, rij2, R, term_lj12, term_lj6,       &
    !$omp         dij, lj12, lj6, grad_coef, work, term_elec, inv_r2,   &
    !$omp         trans_x, trans_y, trans_z, force_local, check_virial)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i = id+1, ncell, nthread
      do ix = 1, natom(i)
        coord_pbc(1:3,ix,i) = coord(1:3,ix,i) + trans1(1:3,ix,i)
      end do
    end do

    !$omp barrier

    ! energy within a cell(solute-solute)
    !
    do i = id+1, ncell_local, nthread

      num_nb15 = 0

      do ix = 1, natom(i) - 1

        rtmp(1:3) = coord_pbc(1:3,ix,i)
        qtmp      = charge(ix,i)
        iatmcls   = atmcls(ix,i)

        ini_nb15  = num_nb15 + 1
        fin_nb15  = num_nb15 + num_nb15_calc1(ix,i)
        num_nb15  = fin_nb15

        force_local(1:3) = 0.0_wp

!ocl norecurrence(force)
!ocl swp
        do k = ini_nb15, fin_nb15

          iy = nb15_calc_list1(k,i)

          dij(1) = rtmp(1) - coord_pbc(1,iy,i)
          dij(2) = rtmp(2) - coord_pbc(2,iy,i)
          dij(3) = rtmp(3) - coord_pbc(3,iy,i)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          inv_r2 = 1.0_wp / rij2
          rij2  = cutoff2 * density * inv_r2

          jatmcls = atmcls(iy,i)
          lj12 = nonb_lj12(iatmcls,jatmcls)
          lj6  = nonb_lj6(iatmcls,jatmcls)

          L  = int(rij2)
          R  = rij2 - L
          L1 = 3*L - 2

          term_lj12 = table_grad(L1  ) + R*(table_grad(L1+3)-table_grad(L1  ))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6   &
                    + qtmp*charge(iy,i)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force(1,iy,i,id+1) = force(1,iy,i,id+1) + work(1)
          force(2,iy,i,id+1) = force(2,iy,i,id+1) + work(2)
          force(3,iy,i,id+1) = force(3,iy,i,id+1) + work(3)

        end do

        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
      end do

    end do

    ! interaction between different cells
    !
    do ij = id+1, maxcell, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)

      num_nb15 = 0

      trans_x = cell_move(1,j,i) * system_size(1)
      trans_y = cell_move(2,j,i) * system_size(2)
      trans_z = cell_move(3,j,i) * system_size(3)
      check_virial = virial_check(j,i)

      do iix = 1, nb15_cell(ij)

        ix = nb15_list(iix,ij)
        rtmp(1:3) = coord_pbc(1:3,ix,i)
        qtmp = charge(ix,i)
        iatmcls = atmcls(ix,i)

        ini_nb15  = num_nb15 + 1
        fin_nb15  = num_nb15 + num_nb15_calc(ix,ij)
        num_nb15  = fin_nb15

        force_local(1:3) = 0.0_wp

!ocl norecurrence(force)
!ocl swp
        do k = ini_nb15, fin_nb15
          iy = nb15_calc_list(k,ij)

          dij(1) = rtmp(1) - coord_pbc(1,iy,j) + trans_x
          dij(2) = rtmp(2) - coord_pbc(2,iy,j) + trans_y
          dij(3) = rtmp(3) - coord_pbc(3,iy,j) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          inv_r2 = 1.0_wp / rij2
          rij2 = cutoff2 * density * inv_r2

          jatmcls = atmcls(iy,j)
          lj12 = nonb_lj12(iatmcls,jatmcls)
          lj6  = nonb_lj6(iatmcls,jatmcls)

          L  = int(rij2)
          R  = rij2 - L
          L1 = 3*L - 2
          term_lj12 = table_grad(L1  ) + R*(table_grad(L1+3)-table_grad(L1  ))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6   &
                    + qtmp*charge(iy,j)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force(1,iy,j,id+1) = force(1,iy,j,id+1) + work(1)
          force(2,iy,j,id+1) = force(2,iy,j,id+1) + work(2)
          force(3,iy,j,id+1) = force(3,iy,j,id+1) + work(3)

        end do

        if (check_virial == 1) &
          virial(1:3,ij) = virial(1:3,ij) - force_local(1:3)
        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_force_nonbond_table_linear

#endif

#ifdef USE_GPU
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_table_linear_gpu
  !> @brief        calculate nonbonded energy with gpu
  !! @authors      JJ
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[in]    npt      : logical to check virial calculation is required
  !! @param[inout] coord_pbc: pbc oriented coordinates
  !! @param[inout] force    : forces for each cell
  !! @param[inout] virial   : virial term of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : van der Waals energy of target systems
  !! @param[inout] ene_virial : energy&virial outputs
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_table_linear_gpu(domain, enefunc, pairlist,&
                                                     npt, coord_pbc, force,    &
                                                     virial, eelec, evdw,      &
                                                     ene_virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    logical,                  intent(in)    :: npt
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)
    real(dp),                 intent(inout) :: eelec(nthread)
    real(dp),                 intent(inout) :: evdw(nthread)
    real(dp),                 intent(inout) :: ene_virial(:)

    ! local variables
    real(wp)                  :: cutoff2, pairlistdist2
    integer                   :: ncell, ncell_local, ncell_bound

    real(dp),         pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: cell_move(:,:,:), system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: density, cutoff, pairlistdist
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer,          pointer :: atmcls(:,:)
    integer,          pointer :: natom(:)
    integer,          pointer :: virial_check(:,:)

    ! for GPU
    integer,          pointer :: univ_cell_pairlist1(:,:)
    integer,          pointer :: univ_cell_pairlist2(:,:)
    integer,          pointer :: univ_ix_natom(:)
    integer(1),       pointer :: univ_ix_list(:,:)
    integer,          pointer :: univ_iy_natom(:)
    integer(1),       pointer :: univ_iy_list(:,:)
    integer,          pointer :: univ_ij_sort_list(:)
    integer(1),       pointer :: univ_mask2(:,:)
    integer,          pointer :: univ_mask2_size
    integer                   :: univ_ncell_nonzero
    integer                   :: univ_update
      
    integer                   :: ncell_max, ij, i, j
    integer                   :: cutoff_int
    integer                   :: num_atom_cls
    integer                   :: ret_shape(1)
    integer                   :: univ_ij, univ_ij0, iiy, idx, base_idx

    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6
    pairlistdist    => enefunc%pairlistdist

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    ncell_bound     =  domain%num_cell_boundary
    cutoff2         =  cutoff * cutoff
    pairlistdist2   =  pairlistdist * pairlistdist

    ret_shape = shape(natom)
    ncell_max = ret_shape(1)
    cutoff_int = enefunc%table%cutoff_int
    num_atom_cls = enefunc%num_atom_cls

    !
    ! univ
    univ_cell_pairlist1 => domain%univ_cell_pairlist1
    univ_cell_pairlist2 => domain%univ_cell_pairlist2
    univ_ix_natom       => pairlist%univ_ix_natom
    univ_ix_list        => pairlist%univ_ix_list
    univ_iy_natom       => pairlist%univ_iy_natom
    univ_iy_list        => pairlist%univ_iy_list
    univ_ij_sort_list   => pairlist%univ_ij_sort_list
    univ_mask2          => pairlist%univ_mask2
    univ_mask2_size     => pairlist%univ_mask2_size
    univ_ncell_nonzero  =  pairlist%univ_ncell_nonzero
    univ_update         =  pairlist%univ_update

    !
    ! launch GPU kernels
    call gpu_launch_compute_energy_nonbond_table_linear_univ( &
         coord_pbc, force(:,:,:,1), ene_virial,               &
         coord, trans1, cell_move, charge, atmcls, natom,     &
         nonb_lj12, nonb_lj6, table_ene, table_grad,          &
         univ_cell_pairlist1, univ_mask2,                     &
         univ_ix_natom, univ_ix_list, univ_iy_natom,          &
         univ_iy_list, univ_ij_sort_list, virial_check,       &
         MaxAtom, MaxAtomCls, num_atom_cls,                   &
         ncell_local, ncell_bound, ncell_max,                 &
         cutoff_int, univ_maxcell, univ_maxcell1,             &
         univ_ncell_nonzero, univ_ncell_near, univ_update,    &
         univ_mask2_size, univ_natom_max, maxcell, density,   &
         cutoff2,                                             &
         system_size(1), system_size(2), system_size(3) )

    return

  end subroutine compute_energy_nonbond_table_linear_gpu

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_table_linear_gpu
  !> @brief        calculate nonbonded force with gpu
  !! @authors      JJ
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[in]    npt      : logical to check virial calculation is required
  !! @param[inout] coord_pbc: pbc oriented coordinates
  !! @param[inout] force    : forces for each cell
  !! @param[inout] virial   : virial term of target systems
  !! @param[inout] ene_virial : energy&virial outputs
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_table_linear_gpu(domain, enefunc, pairlist, &
                                                    npt, cpu_calc, coord_pbc,  &
                                                    force_omp, force, virial,  &
                                                    ene_virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    logical,                  intent(in)    :: npt
    logical,                  intent(in)    :: cpu_calc
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force_omp(:,:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)
    real(dp),                 intent(inout) :: ene_virial(:)

    ! local variables
    real(wp)                  :: cutoff2, pairlistdist2
    integer                   :: ncell, ncell_local, ncell_bound
    ! integer                   :: check_virial

    real(dp),         pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: cell_move(:,:,:), system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: density, cutoff, pairlistdist
    real(wp),         pointer :: table_grad(:)
    integer,          pointer :: cell_pairlist(:,:), atmcls(:,:)
    integer,          pointer :: natom(:)
    integer,          pointer :: virial_check(:,:)

    ! for GPU
    integer,          pointer :: univ_cell_pairlist1(:,:)
    integer,          pointer :: univ_cell_pairlist2(:,:)
    integer,          pointer :: univ_ix_natom(:)
    integer(1),       pointer :: univ_ix_list(:,:)
    integer,          pointer :: univ_iy_natom(:)
    integer(1),       pointer :: univ_iy_list(:,:)
    integer,          pointer :: univ_ij_sort_list(:)
    integer(1),       pointer :: univ_mask2(:,:)
    integer,          pointer :: univ_mask2_size
    integer                   :: univ_ncell_nonzero, univ_gpu_start_index
    integer                   :: univ_update
    integer                   :: ncell_max
    integer                   :: cutoff_int
    integer                   :: num_atom_cls
    integer                   :: ret_shape(1)
    integer                   :: univ_ij, univ_ij0, iiy, idx, base_idx

    cell_pairlist   => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6
    pairlistdist    => enefunc%pairlistdist

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    ncell_bound     =  domain%num_cell_boundary
    cutoff2         =  cutoff * cutoff
    pairlistdist2   =  pairlistdist * pairlistdist

    ret_shape = shape(natom)
    ncell_max = ret_shape(1)
    cutoff_int = enefunc%table%cutoff_int
    num_atom_cls = enefunc%num_atom_cls

    !
    ! univ
    univ_cell_pairlist1 => domain%univ_cell_pairlist1
    univ_cell_pairlist2 => domain%univ_cell_pairlist2
    univ_ix_natom       => pairlist%univ_ix_natom
    univ_ix_list        => pairlist%univ_ix_list
    univ_iy_natom       => pairlist%univ_iy_natom
    univ_iy_list        => pairlist%univ_iy_list
    univ_ij_sort_list   => pairlist%univ_ij_sort_list
    univ_mask2          => pairlist%univ_mask2
    univ_mask2_size     => pairlist%univ_mask2_size
    univ_ncell_nonzero  =  pairlist%univ_ncell_nonzero
    univ_update         =  pairlist%univ_update

    !  launch GPU kernels (data transfer from CPU to GPU as well)
    !
    if (cpu_calc) then
      univ_gpu_start_index = univ_gpu_start
    else
      univ_gpu_start = ncell_local
    end if
    call gpu_launch_compute_force_nonbond_table_linear_univ( &
         coord_pbc, force(:,:,:,1), ene_virial,              &  ! argments
         coord, trans1, cell_move, charge, atmcls, natom,    &  ! arrays
         nonb_lj12, nonb_lj6, table_grad,                    &  ! arrays
         univ_cell_pairlist1, univ_mask2, univ_ix_natom,     &
         univ_ix_list, univ_iy_natom, univ_iy_list,          &
         univ_ij_sort_list, virial_check,                    &
         MaxAtom, MaxAtomCls, num_atom_cls,                  &  ! size
         ncell_local, ncell_bound, ncell_max, cutoff_int,    &  ! size
         univ_maxcell, univ_maxcell1, univ_ncell_nonzero,    &
         univ_ncell_near, univ_update, univ_mask2_size,      &
         univ_natom_max, npt, cpu_calc, density, cutoff2,    &
         pairlistdist2, univ_gpu_start,                      &
         system_size(1), system_size(2), system_size(3))

    if (cpu_calc) then
      call timer(TimerTest10, TimerOn)
      call cpu_compute_force_intra_cell_univ(                &
         coord, trans1, cell_move, system_size, charge,      &
         atmcls, natom, nonb_lj12, nonb_lj6, table_grad,     &  
         univ_cell_pairlist1, univ_mask2, univ_ix_natom,     &
         univ_ix_list, univ_iy_natom, univ_iy_list,          &
         univ_ij_sort_list,                                  &
         ncell_local, ncell_bound, univ_ncell_nonzero,       &
         univ_update, density, cutoff2,                      &
         coord_pbc, force_omp)
      call timer(TimerTest10, TimerOff)
    end if

    return

  end subroutine compute_force_nonbond_table_linear_gpu

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    cpu_compute_force_intra_cell_univ    
  !> @brief        calculate nonbonded force within each cell
  !! @authors      JJ 
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cpu_compute_force_intra_cell_univ(coord, trans1, cell_move,   &
             system_size, charge, atmcls, natom, nonb_lj12, nonb_lj6,      &
             table_grad, univ_cell_pairlist1, univ_mask2,                  &
             univ_ix_natom, univ_ix_list, univ_iy_natom, univ_iy_list,     &
             univ_ij_sort_list, ncell_local, ncell_bound,                  &
             univ_ncell_nonzero, univ_update, density, cutoff2,            &
             coord_pbc, force)

    ! formal arguments
    real(dp),                 intent(in)    :: coord(:,:,:)
    real(wp),                 intent(in)    :: trans1(:,:,:)
    real(wp),                 intent(in)    :: cell_move(:,:,:)
    real(wp),                 intent(in)    :: system_size(:)
    real(wp),                 intent(in)    :: charge(:,:)
    integer,                  intent(in)    :: atmcls(:,:)
    integer,                  intent(in)    :: natom(:)
    real(wp),                 intent(in)    :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),                 intent(in)    :: table_grad(:)
    integer,                  intent(in)    :: univ_cell_pairlist1(:,:)
    integer(1),               intent(in)    :: univ_mask2(:,:)
    integer,                  intent(in)    :: univ_ix_natom(:)
    integer(1),               intent(in)    :: univ_ix_list(:,:)
    integer,                  intent(in)    :: univ_iy_natom(:)
    integer(1),               intent(in)    :: univ_iy_list(:,:)
    integer,                  intent(in)    :: univ_ij_sort_list(:)
    integer,                  intent(in)    :: ncell_local, ncell_bound
    integer,                  intent(in)    :: univ_ncell_nonzero
    integer,                  intent(in)    :: univ_update 
    real(wp),                 intent(in)    :: cutoff2, density
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)

    ! local variables
    real(wp)                  :: dij_list(1:3,MaxAtom), rij2_list(MaxAtom)
    real(wp)                  :: cell_move_local(1:3), dij(1:3), rij2, inv_r2
    real(wp)                  :: R, lj12, lj6
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: work(1:3), grad_coef
    real(wp)                  :: rtmp(1:3), qtmp
    real(wp)                  :: force_local(3), force_j(1:3,MaxAtom)
    integer                   :: j_list(MaxAtom)
    integer                   :: index, univ_ij
    integer                   :: i, j, ix, iy, iix, iiy, k, L, L1, idx
    integer                   :: ix_natom, iy_natom
    integer                   :: id, omp_get_thread_num
    integer                   :: iatmcls,jatmcls
    integer                   :: num_count

    !$omp parallel default(shared)                                            &
    !$omp private(id, i, j, ix, iy, iix, iiy, k, univ_ij, index, idx, j_list, &
    !$omp         ix_natom, iy_natom, num_count, L, L1, iatmcls, jatmcls,     &
    !$omp         cell_move_local, rtmp, qtmp, rij2, R, term_lj12, term_lj6,  &
    !$omp         term_elec, grad_coef, work, dij, inv_r2, lj12, lj6,         &
    !$omp         dij_list, rij2_list, force_local, force_j)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell_local+ncell_bound, nthread
      do ix = 1, natom(i)
        coord_pbc(1:3,ix,i) = coord(1:3,ix,i) + trans1(1:3,ix,i)
      end do
    end do

    !$omp barrier

    do i = id+1, ncell_local, nthread

      do ix = 1, natom(i) - 1

        rtmp(1:3) = coord_pbc(1:3,ix,i)
        qtmp      = charge(ix,i)
        iatmcls   = atmcls(ix,i)

        force_local(1:3) = 0.0_wp
        num_count        = 0

!ocl norecurrence(force)
!ocl swp
        do iy = ix + 1, natom(i)

          idx = iy + (ix-1)*univ_natom_max

          if (univ_mask2(idx,i) == 1) then

            dij(1:3) = rtmp(1:3) - coord_pbc(1:3,iy,i)
            rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            if (rij2 < cutoff2) then
              num_count = num_count + 1
              dij_list(1:3,num_count)  = dij(1:3)
              rij2_list(num_count) = rij2
              j_list(num_count)   = iy
            end if
          end if
        end do

        do k = 1, num_count

          dij(1:3)  = dij_list(1:3,k)
          iy      = j_list(k)
          rij2    = cutoff2*density/rij2_list(k)
          jatmcls = atmcls(iy,i)
          lj6     = nonb_lj6(iatmcls,jatmcls)
          lj12    = nonb_lj12(iatmcls,jatmcls)

          L  = int(rij2)
          R  = rij2 - L
          L1 = 3*L - 2

          term_lj12 = table_grad(L1  ) + R*(table_grad(L1+3)-table_grad(L1  ))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6   &
                    + qtmp*charge(iy,i)*term_elec

          work(1:3) = grad_coef*dij(1:3)

          ! store force
          !
          force_local(1:3) = force_local(1:3) - work(1:3)
          force_j(1:3,k) = work(1:3)

        end do

        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
        do k = 1, num_count
          iy = j_list(k)
          force(1:3,iy,i,id+1) = force(1:3,iy,i,id+1) + force_j(1:3,k)
        end do

      end do

    end do

    do index = ncell_local+id+1, univ_gpu_start, nthread

      univ_ij = univ_ij_sort_list(index)
      i = univ_cell_pairlist1(1,univ_ij)
      j = univ_cell_pairlist1(2,univ_ij)
!     ix_natom = univ_ix_natom(univ_ij)
!     iy_natom = univ_iy_natom(univ_ij)
      cell_move_local(1:3) = cell_move(1:3,j,i) * system_size(1:3)

      if (univ_ij > univ_ncell_near) then

!       do iix = 1, ix_natom
        do ix = 1, natom(i)

!         ix = univ_ix_list(iix,univ_ij)
          rtmp(1:3) = coord_pbc(1:3,ix,i) + cell_move_local(1:3)
          qtmp      = charge(ix,i)
          iatmcls   = atmcls(ix,i)

          force_local(1:3) = 0.0_wp
          num_count = 0

!         do iiy = 1, iy_natom
          do iy = 1, natom(j)
!           iy = univ_iy_list(iiy,univ_ij)
            dij(1:3) = rtmp(1:3) - coord_pbc(1:3,iy,j)  
            rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
            if (rij2 < cutoff2) then
              num_count = num_count + 1
              dij_list(1:3,num_count)  = dij(1:3)
              rij2_list(num_count) = rij2
              j_list(num_count)   = iy
            end if
          end do

          do k = 1, num_count

            dij(1)  = dij_list(1,k)
            dij(2)  = dij_list(2,k)
            dij(3)  = dij_list(3,k)
            iy      = j_list(k)
            rij2    = cutoff2*density/rij2_list(k)
            jatmcls = atmcls(iy,j)
            lj6     = nonb_lj6(iatmcls,jatmcls)
            lj12    = nonb_lj12(iatmcls,jatmcls)

            L  = int(rij2)
            R  = rij2 - L
            L1 = 3*L - 2

            term_lj12 = table_grad(L1  ) + R*(table_grad(L1+3)-table_grad(L1))
            term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
            term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
            grad_coef = term_lj12*lj12 - term_lj6*lj6   &
                      + qtmp*charge(iy,j)*term_elec
            work(1:3) = grad_coef*dij(1:3)

            ! store force
            !
            force_local(1:3) = force_local(1:3) - work(1:3)
            force_j(1:3,k)   = work(1:3)

          end do

          force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
          do k = 1, num_count
            iy = j_list(k)
            force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + force_j(1:3,k)
          end do

        end do

      else

!       do iix = 1, ix_natom
        do ix = 1, natom(i)

!         ix = univ_ix_list(iix,univ_ij)
          rtmp(1:3) = coord_pbc(1:3,ix,i) + cell_move_local(1:3)
          qtmp      = charge(ix,i)
          iatmcls   = atmcls(ix,i)

          force_local(1:3) = 0.0_wp
          num_count = 0

!         do iiy = 1, iy_natom
          do iy = 1, natom(j)
!           iy = univ_iy_list(iiy,univ_ij)
            idx = iy + (ix-1)*univ_natom_max
            if (univ_mask2(idx,univ_ij) == 1) then
              dij(1:3) = rtmp(1:3) - coord_pbc(1:3,iy,j)
              rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
              if (rij2 < cutoff2) then
                num_count = num_count + 1
                dij_list(1:3,num_count)  = dij(1:3)
                rij2_list(num_count) = rij2
                j_list(num_count)   = iy
              end if
            end if
          end do

          do k = 1, num_count

            dij(1)  = dij_list(1,k)
            dij(2)  = dij_list(2,k)
            dij(3)  = dij_list(3,k)
            iy      = j_list(k)
            rij2    = cutoff2*density/rij2_list(k)
            jatmcls = atmcls(iy,j)
            lj6     = nonb_lj6(iatmcls,jatmcls)
            lj12    = nonb_lj12(iatmcls,jatmcls)

            L  = int(rij2)
            R  = rij2 - L
            L1 = 3*L - 2

            term_lj12 = table_grad(L1  ) + R*(table_grad(L1+3)-table_grad(L1))
            term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
            term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
            grad_coef = term_lj12*lj12 - term_lj6*lj6   &
                      + qtmp*charge(iy,j)*term_elec
            work(1:3) = grad_coef*dij(1:3)

            ! store force
            !
            force_local(1:3) = force_local(1:3) - work(1:3)
            force_j(1:3,k)   = work(1:3)

          end do

          force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
          do k = 1, num_count
            iy = j_list(k)
            force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + force_j(1:3,k)
          end do

        end do

      end if

    end do

    !$omp end parallel 

    return

  end subroutine cpu_compute_force_intra_cell_univ

#endif

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_table_linear_check_fep
  !> @brief        calculate nonbonded energy with lookup table and soft core
  !                for FEP calculations
  !! @authors      NK, HO
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] virial   : virial term of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_table_linear_check_fep( &
                                                 domain, enefunc, pairlist, &
                                                 coord_pbc, force, virial,  &
                                                 eelec, evdw)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:)
    real(dp),                 intent(inout) :: eelec(nthread)
    real(dp),                 intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: lj6, lj12
    real(wp)                  :: R
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, term_lj, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: trans_x, trans_y, trans_z
    real(wp)                  :: rtmp(1:3), qtmp
    real(wp)                  :: dij_list(3,MaxAtom)
    real(wp)                  :: force_local(3), rij2_list(MaxAtom)
    real(wp)                  :: force_localj(3,MaxAtom)
    real(wp)                  :: minimum_contact
    integer                   :: i, ix, iy, j, k, ij, iix, iwater, list, L, L1
    integer                   :: num_nb15, ini_nb15, fin_nb15, num_count
    integer                   :: id, omp_get_thread_num, j_list(MaxAtom)
    integer                   :: ncell, ncell_local
    integer                   :: check_virial

    real(dp),         pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: cell_move(:,:,:), system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: density, cutoff
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer,          pointer :: cell_pairlist(:,:), atmcls(:,:)
    integer,          pointer :: natom(:), nsolute(:), nwater(:)
    integer,          pointer :: nb15_list(:,:)
    integer,          pointer :: nb15_cell(:)
    integer,          pointer :: num_nb15_calc1(:,:), num_nb15_calc(:,:)
    integer,          pointer :: nb15_calc_list1(:,:), nb15_calc_list(:,:)
    integer,          pointer :: virial_check(:,:)

    ! FEP
    real(wp)                  :: evdw_temp, elec_temp
    real(wP)                  :: lamblj, lambel
    integer                   :: fg1, fg2
    integer,          pointer :: nb15_list_fep(:,:)
    integer,          pointer :: nb15_cell_fep(:)
    integer,          pointer :: num_nb15_calc1_fep(:,:), num_nb15_calc_fep(:,:)
    integer,          pointer :: nb15_calc_list1_fep(:,:), nb15_calc_list_fep(:,:)
    integer,          pointer :: fepgrp(:,:)
    real(wp),         pointer :: table_nonb_lambda(:,:,:)

    cell_pairlist   => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    nsolute         => domain%num_solute
    nwater          => domain%num_water
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6

    nb15_list       => pairlist%nb15_list
    nb15_cell       => pairlist%nb15_cell
    num_nb15_calc1  => pairlist%num_nb15_calc1
    num_nb15_calc   => pairlist%num_nb15_calc
    nb15_calc_list1 => pairlist%nb15_calc_list1
    nb15_calc_list  => pairlist%nb15_calc_list

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff
    minimum_contact =  enefunc%minimum_contact

    ! FEP
    nb15_list_fep       => pairlist%nb15_list_fep
    nb15_cell_fep       => pairlist%nb15_cell_fep
    num_nb15_calc1_fep  => pairlist%num_nb15_calc1_fep
    num_nb15_calc_fep   => pairlist%num_nb15_calc_fep
    nb15_calc_list1_fep => pairlist%nb15_calc_list1_fep
    nb15_calc_list_fep  => pairlist%nb15_calc_list_fep
    fepgrp              => domain%fepgrp
    table_nonb_lambda   => enefunc%table_nonb_lambda

    ! calculate energy and gradient
    !

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, ix, num_nb15, rtmp, qtmp, ini_nb15, fin_nb15, k, iy,  &
    !$omp         rij2, L, R, term_lj12, term_lj6, term_lj, grad_coef, work,   &
    !$omp         term_elec, iwater, ij, j, trans_x, trans_y,                  &
    !$omp         trans_z, iix, list, num_count, j_list, dij_list,             &
    !$omp         force_local, force_localj, lj6, lj12, L1, dij, check_virial, &
    !$omp         elec_temp, evdw_temp, rij2_list, fg1, fg2, lamblj, lambel)

#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell, nthread
      do ix = 1, natom(i)
        coord_pbc(1,ix,i) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(2,ix,i) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(3,ix,i) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do

    !$omp barrier

    ! Unperturbed region
    !

    ! energy within a cell
    !
    do i = id+1, ncell_local, nthread

      num_nb15 = 0
      do ix = 1, natom(i) - 1

        rtmp(1:3) = coord_pbc(1:3,ix,i)
        qtmp = charge(ix,i)
        ini_nb15  = num_nb15 + 1
        fin_nb15  = num_nb15 + num_nb15_calc1(ix,i)
        num_nb15  = fin_nb15
        num_count = 0

        do k = ini_nb15, fin_nb15

          iy = nb15_calc_list1(k,i)

          ! compute distance
          !
          dij(1) = rtmp(1) - coord_pbc(1,iy,i)
          dij(2) = rtmp(2) - coord_pbc(2,iy,i)
          dij(3) = rtmp(3) - coord_pbc(3,iy,i)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1:3,num_count)  = dij(1:3)
            rij2_list(num_count) = max(rij2, minimum_contact)
            j_list(num_count)   = iy
          end if
        end do

        force_local(1:3) = 0.0_wp
        evdw_temp   = 0.0_wp
        elec_temp   = 0.0_wp

!ocl norecurrence(force)
!ocl swp
        do k = 1, num_count
          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)
          iy   = j_list(k)
          rij2  = cutoff2*density/rij2_list(k)
          lj6  = nonb_lj6(atmcls(ix,i),atmcls(iy,i))
          lj12 = nonb_lj12(atmcls(ix,i),atmcls(iy,i))

          L    = int(rij2)
          R    = rij2 - L

          L1   = 3*L - 2
          term_lj12 = table_ene(L1) + R*(table_ene(L1+3)-table_ene(L1))
          term_lj6  = table_ene(L1+1) + R*(table_ene(L1+4)-table_ene(L1+1))
          term_elec = table_ene(L1+2) + R*(table_ene(L1+5)-table_ene(L1+2))
          evdw_temp = evdw_temp + (term_lj12*lj12 - term_lj6*lj6)
          elec_temp = elec_temp + qtmp*charge(iy,i)*term_elec

          term_lj12 = table_grad(L1) + R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6                           &
                     + qtmp*charge(iy,i)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force_localj(1,k)  = work(1)
          force_localj(2,k)  = work(2)
          force_localj(3,k)  = work(3)

        end do

        do k = 1, num_count
          iy = j_list(k)
          force(1:3,iy,i,id+1) = force(1:3,iy,i,id+1) + force_localj(1:3,k)
        end do
        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw(id+1) = evdw(id+1) + evdw_temp

      end do

    end do


    ! interaction between different cells
    !
    do ij = id+1, maxcell, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)

      num_nb15 = 0

      trans_x = cell_move(1,j,i) * system_size(1)
      trans_y = cell_move(2,j,i) * system_size(2)
      trans_z = cell_move(3,j,i) * system_size(3)
      check_virial = virial_check(j,i)

      do iix = 1, nb15_cell(ij)

        ix   = nb15_list(iix,ij)
        rtmp(1:3) = coord_pbc(1:3,ix,i)
        qtmp = charge(ix,i)

        ini_nb15 = num_nb15 + 1
        fin_nb15 = num_nb15 + num_nb15_calc(ix,ij)
        num_nb15 = fin_nb15

        num_count = 0
        do k = ini_nb15, fin_nb15
          iy   = nb15_calc_list(k,ij)
          dij(1) = rtmp(1) - coord_pbc(1,iy,j) + trans_x
          dij(2) = rtmp(2) - coord_pbc(2,iy,j) + trans_y
          dij(3) = rtmp(3) - coord_pbc(3,iy,j) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1:3,num_count)  = dij(1:3)
            rij2_list(num_count) = max(rij2, minimum_contact)
            j_list(num_count)   = iy
          end if
        end do

        force_local(1:3) = 0.0_wp
        evdw_temp  = 0.0_wp
        elec_temp  = 0.0_wp

!ocl norecurrence(force)
!ocl swp
        do k = 1, num_count
          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)
          iy   = j_list(k)
          rij2  = cutoff2*density/rij2_list(k)
          lj6  = nonb_lj6 (atmcls(ix,i),atmcls(iy,j))
          lj12 = nonb_lj12(atmcls(ix,i),atmcls(iy,j))

          L    = int(rij2)
          R    = rij2 - L

          L1   = 3*L - 2
          term_lj12 = table_ene(L1) + R*(table_ene(L1+3)-table_ene(L1))
          term_lj6  = table_ene(L1+1) + R*(table_ene(L1+4)-table_ene(L1+1))
          term_elec = table_ene(L1+2) + R*(table_ene(L1+5)-table_ene(L1+2))
          evdw_temp = evdw_temp + (term_lj12*lj12 - term_lj6*lj6)
          elec_temp = elec_temp + qtmp*charge(iy,j)*term_elec

          term_lj12 = table_grad(L1) + R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6                           &
                     + qtmp*charge(iy,j)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)

          force_localj(1,k) = work(1)
          force_localj(2,k) = work(2)
          force_localj(3,k) = work(3)

        end do

        do k = 1, num_count
          iy = j_list(k)
          force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + force_localj(1:3,k)
        end do
        if (check_virial == 1) &
          virial(1:3,ij) = virial(1:3,ij) - force_local(1:3)
        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
        evdw(id+1)  = evdw(id+1)  + evdw_temp
        eelec(id+1) = eelec(id+1) + elec_temp

      end do

    end do


    ! Perturbed region including
    ! singleA, singleB, dualA, and dualB
    !

    ! energy within a cell
    !
    do i = id+1, ncell_local, nthread

      num_nb15 = 0
      do ix = 1, natom(i) - 1

        rtmp(1:3) = coord_pbc(1:3,ix,i)
        qtmp = charge(ix,i)
        ini_nb15  = num_nb15 + 1
        fin_nb15  = num_nb15 + num_nb15_calc1_fep(ix,i)
        num_nb15  = fin_nb15
        num_count = 0

        do k = ini_nb15, fin_nb15

          iy = nb15_calc_list1_fep(k,i)

          ! compute distance
          !
          dij(1) = rtmp(1) - coord_pbc(1,iy,i)
          dij(2) = rtmp(2) - coord_pbc(2,iy,i)
          dij(3) = rtmp(3) - coord_pbc(3,iy,i)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1:3,num_count)  = dij(1:3)
            rij2_list(num_count) = max(rij2, minimum_contact)
            j_list(num_count)   = iy
          end if
        end do

        force_local(1:3) = 0.0_wp
        evdw_temp   = 0.0_wp
        elec_temp   = 0.0_wp

!ocl norecurrence(force)
!ocl swp
        do k = 1, num_count
          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)
          iy   = j_list(k)
          lj6  = nonb_lj6(atmcls(ix,i),atmcls(iy,i))
          lj12 = nonb_lj12(atmcls(ix,i),atmcls(iy,i))

          ! FEP: Determine lamblj and lambel
          fg1 = fepgrp(ix,i)
          fg2 = fepgrp(iy,i)
          lamblj = table_nonb_lambda(1,fg1,fg2)
          lambel = table_nonb_lambda(2,fg1,fg2)

          ! FEP: LJ with soft core for reference
          rij2      = cutoff2*density/(rij2_list(k)+table_nonb_lambda(3,fg1,fg2))
          L         = int(rij2)
          R         = rij2 - L
          L1        = 3*L - 2
          term_lj12 = table_ene(L1) + R*(table_ene(L1+3)-table_ene(L1))
          term_lj6  = table_ene(L1+1) + R*(table_ene(L1+4)-table_ene(L1+1))
          evdw_temp = evdw_temp + lamblj*(term_lj12*lj12 - term_lj6*lj6)
          term_lj12 = table_grad(L1)+ R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1)+ R*(table_grad(L1+4)-table_grad(L1+1))

          ! FEP: EL with soft core for reference
          rij2      = cutoff2*density/(rij2_list(k)+table_nonb_lambda(4,fg1,fg2))
          L         = int(rij2)
          R         = rij2 - L
          L1        = 3*L - 2
          term_elec = table_ene(L1+2) + R*(table_ene(L1+5)-table_ene(L1+2))
          elec_temp = elec_temp + lambel*qtmp*charge(iy,i)*term_elec
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))

          ! gradient
          grad_coef = lamblj*(term_lj12*lj12 - term_lj6*lj6) &
            + lambel*qtmp*charge(iy,i)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force_localj(1,k)  = work(1)
          force_localj(2,k)  = work(2)
          force_localj(3,k)  = work(3)
        end do

        do k = 1, num_count
          iy = j_list(k)
          force(1:3,iy,i,id+1) = force(1:3,iy,i,id+1) + force_localj(1:3,k)
        end do
        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
        evdw(id+1)  = evdw(id+1)  + evdw_temp
        eelec(id+1) = eelec(id+1) + elec_temp

      end do

    end do


    ! interaction between different cells
    !
    do ij = id+1, maxcell, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)

      num_nb15 = 0

      trans_x = cell_move(1,j,i) * system_size(1)
      trans_y = cell_move(2,j,i) * system_size(2)
      trans_z = cell_move(3,j,i) * system_size(3)
      check_virial = virial_check(j,i)

      do iix = 1, nb15_cell_fep(ij)

        ix   = nb15_list_fep(iix,ij)
        rtmp(1:3) = coord_pbc(1:3,ix,i)
        qtmp = charge(ix,i)

        ini_nb15 = num_nb15 + 1
        fin_nb15 = num_nb15 + num_nb15_calc_fep(ix,ij)
        num_nb15 = fin_nb15

        num_count = 0
        do k = ini_nb15, fin_nb15
          iy   = nb15_calc_list_fep(k,ij)
          dij(1) = rtmp(1) - coord_pbc(1,iy,j) + trans_x
          dij(2) = rtmp(2) - coord_pbc(2,iy,j) + trans_y
          dij(3) = rtmp(3) - coord_pbc(3,iy,j) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1:3,num_count)  = dij(1:3)
            rij2_list(num_count) = max(rij2, minimum_contact)
            j_list(num_count)   = iy
          end if
        end do

        force_local(1:3) = 0.0_wp
        evdw_temp  = 0.0_wp
        elec_temp  = 0.0_wp

!ocl norecurrence(force)
!ocl swp
        do k = 1, num_count
          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)
          iy   = j_list(k)
          lj6  = nonb_lj6 (atmcls(ix,i),atmcls(iy,j))
          lj12 = nonb_lj12(atmcls(ix,i),atmcls(iy,j))

          ! FEP: Determine lamblj and lambel
          fg1 = fepgrp(ix,i)
          fg2 = fepgrp(iy,j)
          lamblj = table_nonb_lambda(1,fg1,fg2)
          lambel = table_nonb_lambda(2,fg1,fg2)

          ! FEP: LJ with soft core for reference
          rij2      = cutoff2*density/(rij2_list(k)+table_nonb_lambda(3,fg1,fg2))
          L         = int(rij2)
          R         = rij2 - L
          L1        = 3*L - 2
          term_lj12 = table_ene(L1) + R*(table_ene(L1+3)-table_ene(L1))
          term_lj6  = table_ene(L1+1) + R*(table_ene(L1+4)-table_ene(L1+1))
          evdw_temp = evdw_temp + lamblj*(term_lj12*lj12 - term_lj6*lj6)
          term_lj12 = table_grad(L1)+ R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1)+ R*(table_grad(L1+4)-table_grad(L1+1))

          ! FEP: EL with soft core for reference
          rij2      = cutoff2*density/(rij2_list(k)+table_nonb_lambda(4,fg1,fg2))
          L         = int(rij2)
          R         = rij2 - L
          L1        = 3*L - 2
          term_elec = table_ene(L1+2) + R*(table_ene(L1+5)-table_ene(L1+2))
          elec_temp = elec_temp + lambel*qtmp*charge(iy,j)*term_elec
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))

          ! gradient
          grad_coef = lamblj*(term_lj12*lj12 - term_lj6*lj6) &
            + lambel*qtmp*charge(iy,j)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)

          force_localj(1,k) = work(1)
          force_localj(2,k) = work(2)
          force_localj(3,k) = work(3)

        end do

        do k = 1, num_count
          iy = j_list(k)
          force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + force_localj(1:3,k)
        end do
        if (check_virial .eq. 1) &
          virial(1:3,ij) = virial(1:3,ij) - force_local(1:3)
        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
        evdw(id+1)  = evdw(id+1)  + evdw_temp
        eelec(id+1) = eelec(id+1) + elec_temp

      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond_table_linear_check_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_table_linear_fep
  !> @brief        calculate nonbonded energy with lookup table and soft core
  !!               for FEP calculation
  !! @authors      NK, HO
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] virial   : virial term of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_table_linear_fep( &
                                                 domain, enefunc, pairlist, &
                                                 coord_pbc, force, virial,  &
                                                 eelec, evdw)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:)
    real(dp),                 intent(inout) :: eelec(nthread)
    real(dp),                 intent(inout) :: evdw(nthread)

    ! local variables
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: lj6, lj12
    real(wp)                  :: R
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, term_lj, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: trans_x, trans_y, trans_z
    real(wp)                  :: rtmp(1:3), qtmp
    real(wp)                  :: dij_list(3,MaxAtom)
    real(wp)                  :: force_local(3), rij2_list(MaxAtom)
    real(wp)                  :: force_localj(3,MaxAtom)
    integer                   :: i, ix, iy, j, k, ij, iix, iwater, list, L, L1
    integer                   :: num_nb15, ini_nb15, fin_nb15, num_count
    integer                   :: id, omp_get_thread_num, j_list(MaxAtom)
    integer                   :: ncell, ncell_local
    integer                   :: check_virial

    real(dp),         pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: cell_move(:,:,:), system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: density, cutoff
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer,          pointer :: cell_pairlist(:,:), atmcls(:,:)
    integer,          pointer :: natom(:), nsolute(:), nwater(:)
    integer,          pointer :: nb15_list(:,:)
    integer,          pointer :: nb15_cell(:)
    integer,          pointer :: num_nb15_calc1(:,:), num_nb15_calc(:,:)
    integer,          pointer :: nb15_calc_list1(:,:), nb15_calc_list(:,:)
    integer,          pointer :: virial_check(:,:)

    ! FEP
    real(wp)                  :: evdw_temp, elec_temp
    real(wP)                  :: lamblj, lambel
    integer                   :: fg1, fg2
    integer,          pointer :: nb15_list_fep(:,:)
    integer,          pointer :: nb15_cell_fep(:)
    integer,          pointer :: num_nb15_calc1_fep(:,:), num_nb15_calc_fep(:,:)
    integer,          pointer :: nb15_calc_list1_fep(:,:), nb15_calc_list_fep(:,:)
    integer,          pointer :: fepgrp(:,:)
    real(wp),         pointer :: table_nonb_lambda(:,:,:)

    cell_pairlist   => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    nsolute         => domain%num_solute
    nwater          => domain%num_water
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6

    nb15_list       => pairlist%nb15_list
    nb15_cell       => pairlist%nb15_cell
    num_nb15_calc1  => pairlist%num_nb15_calc1
    num_nb15_calc   => pairlist%num_nb15_calc
    nb15_calc_list1 => pairlist%nb15_calc_list1
    nb15_calc_list  => pairlist%nb15_calc_list

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff

    ! FEP
    nb15_list_fep       => pairlist%nb15_list_fep
    nb15_cell_fep       => pairlist%nb15_cell_fep
    num_nb15_calc1_fep  => pairlist%num_nb15_calc1_fep
    num_nb15_calc_fep   => pairlist%num_nb15_calc_fep
    nb15_calc_list1_fep => pairlist%nb15_calc_list1_fep
    nb15_calc_list_fep  => pairlist%nb15_calc_list_fep
    fepgrp              => domain%fepgrp
    table_nonb_lambda   => enefunc%table_nonb_lambda

    ! calculate energy and gradient
    !

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, ix, num_nb15, rtmp, qtmp, ini_nb15, fin_nb15, k, iy,  &
    !$omp         rij2, L, R, term_lj12, term_lj6, term_lj, grad_coef, work,   &
    !$omp         term_elec, iwater, ij, j, trans_x, trans_y,                  &
    !$omp         trans_z, iix, list, num_count, j_list, dij_list,             &
    !$omp         force_local, force_localj, lj6, lj12, L1, dij, check_virial, &
    !$omp         elec_temp, evdw_temp, rij2_list, fg1, fg2, lamblj, lambel)

#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell, nthread
      do ix = 1, natom(i)
        coord_pbc(1,ix,i) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(2,ix,i) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(3,ix,i) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do

    !$omp barrier

    ! Unperturbed region
    !

    ! energy within a cell
    !
    do i = id+1, ncell_local, nthread

      num_nb15 = 0
      do ix = 1, natom(i) - 1

        rtmp(1:3) = coord_pbc(1:3,ix,i)
        qtmp = charge(ix,i)
        ini_nb15  = num_nb15 + 1
        fin_nb15  = num_nb15 + num_nb15_calc1(ix,i)
        num_nb15  = fin_nb15
        num_count = 0

        do k = ini_nb15, fin_nb15

          iy = nb15_calc_list1(k,i)

          ! compute distance
          !
          dij(1) = rtmp(1) - coord_pbc(1,iy,i)
          dij(2) = rtmp(2) - coord_pbc(2,iy,i)
          dij(3) = rtmp(3) - coord_pbc(3,iy,i)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1:3,num_count)  = dij(1:3)
            rij2_list(num_count) = rij2
            j_list(num_count)   = iy
          end if
        end do

        force_local(1:3) = 0.0_wp
        evdw_temp   = 0.0_wp
        elec_temp   = 0.0_wp

!ocl norecurrence(force)
!ocl swp
        do k = 1, num_count
          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)
          iy   = j_list(k)
          rij2  = cutoff2*density/rij2_list(k)
          lj6  = nonb_lj6(atmcls(ix,i),atmcls(iy,i))
          lj12 = nonb_lj12(atmcls(ix,i),atmcls(iy,i))

          L    = int(rij2)
          R    = rij2 - L

          L1   = 3*L - 2
          term_lj12 = table_ene(L1) + R*(table_ene(L1+3)-table_ene(L1))
          term_lj6  = table_ene(L1+1) + R*(table_ene(L1+4)-table_ene(L1+1))
          term_elec = table_ene(L1+2) + R*(table_ene(L1+5)-table_ene(L1+2))
          evdw_temp = evdw_temp + (term_lj12*lj12 - term_lj6*lj6)
          elec_temp = elec_temp + qtmp*charge(iy,i)*term_elec

          term_lj12 = table_grad(L1) + R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6                           &
                     + qtmp*charge(iy,i)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force_localj(1,k)  = work(1)
          force_localj(2,k)  = work(2)
          force_localj(3,k)  = work(3)

        end do

        do k = 1, num_count
          iy = j_list(k)
          force(1:3,iy,i,id+1) = force(1:3,iy,i,id+1) + force_localj(1:3,k)
        end do
        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
        eelec(id+1) = eelec(id+1) + elec_temp
        evdw(id+1) = evdw(id+1) + evdw_temp

      end do

    end do


    ! interaction between different cells
    !
    do ij = id+1, maxcell, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)

      num_nb15 = 0

      trans_x = cell_move(1,j,i) * system_size(1)
      trans_y = cell_move(2,j,i) * system_size(2)
      trans_z = cell_move(3,j,i) * system_size(3)
      check_virial = virial_check(j,i)

      do iix = 1, nb15_cell(ij)

        ix   = nb15_list(iix,ij)
        rtmp(1:3) = coord_pbc(1:3,ix,i)
        qtmp = charge(ix,i)

        ini_nb15 = num_nb15 + 1
        fin_nb15 = num_nb15 + num_nb15_calc(ix,ij)
        num_nb15 = fin_nb15

        num_count = 0
        do k = ini_nb15, fin_nb15
          iy   = nb15_calc_list(k,ij)
          dij(1) = rtmp(1) - coord_pbc(1,iy,j) + trans_x
          dij(2) = rtmp(2) - coord_pbc(2,iy,j) + trans_y
          dij(3) = rtmp(3) - coord_pbc(3,iy,j) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1:3,num_count)  = dij(1:3)
            rij2_list(num_count) = rij2
            j_list(num_count)   = iy
          end if
        end do

        force_local(1:3) = 0.0_wp
        evdw_temp  = 0.0_wp
        elec_temp  = 0.0_wp

!ocl norecurrence(force)
!ocl swp
        do k = 1, num_count
          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)
          iy   = j_list(k)
          rij2  = cutoff2*density/rij2_list(k)
          lj6  = nonb_lj6 (atmcls(ix,i),atmcls(iy,j))
          lj12 = nonb_lj12(atmcls(ix,i),atmcls(iy,j))

          L    = int(rij2)
          R    = rij2 - L

          L1   = 3*L - 2
          term_lj12 = table_ene(L1) + R*(table_ene(L1+3)-table_ene(L1))
          term_lj6  = table_ene(L1+1) + R*(table_ene(L1+4)-table_ene(L1+1))
          term_elec = table_ene(L1+2) + R*(table_ene(L1+5)-table_ene(L1+2))
          evdw_temp = evdw_temp + (term_lj12*lj12 - term_lj6*lj6)
          elec_temp = elec_temp + qtmp*charge(iy,j)*term_elec

          term_lj12 = table_grad(L1) + R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6                           &
                     + qtmp*charge(iy,j)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)

          force_localj(1,k) = work(1)
          force_localj(2,k) = work(2)
          force_localj(3,k) = work(3)

        end do

        do k = 1, num_count
          iy = j_list(k)
          force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + force_localj(1:3,k)
        end do
        if (check_virial == 1) &
          virial(1:3,ij) = virial(1:3,ij) - force_local(1:3)
        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
        evdw(id+1)  = evdw(id+1)  + evdw_temp
        eelec(id+1) = eelec(id+1) + elec_temp

      end do

    end do


    ! Perturbed region including
    ! singleA, singleB, dualA, and dualB
    !

    ! energy within a cell
    !
    do i = id+1, ncell_local, nthread

      num_nb15 = 0
      do ix = 1, natom(i) - 1

        rtmp(1:3) = coord_pbc(1:3,ix,i)
        qtmp = charge(ix,i)
        ini_nb15  = num_nb15 + 1
        fin_nb15  = num_nb15 + num_nb15_calc1_fep(ix,i)
        num_nb15  = fin_nb15
        num_count = 0

        do k = ini_nb15, fin_nb15

          iy = nb15_calc_list1_fep(k,i)

          ! compute distance
          !
          dij(1) = rtmp(1) - coord_pbc(1,iy,i)
          dij(2) = rtmp(2) - coord_pbc(2,iy,i)
          dij(3) = rtmp(3) - coord_pbc(3,iy,i)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1:3,num_count)  = dij(1:3)
            rij2_list(num_count) = rij2
            j_list(num_count)   = iy
          end if
        end do

        force_local(1:3) = 0.0_wp
        evdw_temp   = 0.0_wp
        elec_temp   = 0.0_wp

!ocl norecurrence(force)
!ocl swp
        do k = 1, num_count
          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)
          iy   = j_list(k)
          lj6  = nonb_lj6(atmcls(ix,i),atmcls(iy,i))
          lj12 = nonb_lj12(atmcls(ix,i),atmcls(iy,i))

          ! FEP: Determine lamblj and lambel
          fg1 = fepgrp(ix,i)
          fg2 = fepgrp(iy,i)
          lamblj = table_nonb_lambda(1,fg1,fg2)
          lambel = table_nonb_lambda(2,fg1,fg2)

          ! FEP: LJ with soft core for reference
          rij2      = cutoff2*density/(rij2_list(k)+table_nonb_lambda(3,fg1,fg2))
          L         = int(rij2)
          R         = rij2 - L
          L1        = 3*L - 2
          term_lj12 = table_ene(L1) + R*(table_ene(L1+3)-table_ene(L1))
          term_lj6  = table_ene(L1+1) + R*(table_ene(L1+4)-table_ene(L1+1))
          evdw_temp = evdw_temp + lamblj*(term_lj12*lj12 - term_lj6*lj6)
          term_lj12 = table_grad(L1)+ R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1)+ R*(table_grad(L1+4)-table_grad(L1+1))

          ! FEP: EL with soft core for reference
          rij2      = cutoff2*density/(rij2_list(k)+table_nonb_lambda(4,fg1,fg2))
          L         = int(rij2)
          R         = rij2 - L
          L1        = 3*L - 2
          term_elec = table_ene(L1+2) + R*(table_ene(L1+5)-table_ene(L1+2))
          elec_temp = elec_temp + lambel*qtmp*charge(iy,i)*term_elec
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))

          ! gradient
          grad_coef = lamblj*(term_lj12*lj12 - term_lj6*lj6) &
            + lambel*qtmp*charge(iy,i)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force_localj(1,k)  = work(1)
          force_localj(2,k)  = work(2)
          force_localj(3,k)  = work(3)
        end do

        do k = 1, num_count
          iy = j_list(k)
          force(1:3,iy,i,id+1) = force(1:3,iy,i,id+1) + force_localj(1:3,k)
        end do
        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
        evdw(id+1)  = evdw(id+1)  + evdw_temp
        eelec(id+1) = eelec(id+1) + elec_temp

      end do

    end do


    ! interaction between different cells
    !
    do ij = id+1, maxcell, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)

      num_nb15 = 0

      trans_x = cell_move(1,j,i) * system_size(1)
      trans_y = cell_move(2,j,i) * system_size(2)
      trans_z = cell_move(3,j,i) * system_size(3)
      check_virial = virial_check(j,i)

      do iix = 1, nb15_cell_fep(ij)

        ix   = nb15_list_fep(iix,ij)
        rtmp(1:3) = coord_pbc(1:3,ix,i)
        qtmp = charge(ix,i)

        ini_nb15 = num_nb15 + 1
        fin_nb15 = num_nb15 + num_nb15_calc_fep(ix,ij)
        num_nb15 = fin_nb15

        num_count = 0
        do k = ini_nb15, fin_nb15
          iy   = nb15_calc_list_fep(k,ij)
          dij(1) = rtmp(1) - coord_pbc(1,iy,j) + trans_x
          dij(2) = rtmp(2) - coord_pbc(2,iy,j) + trans_y
          dij(3) = rtmp(3) - coord_pbc(3,iy,j) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1:3,num_count)  = dij(1:3)
            rij2_list(num_count) = rij2
            j_list(num_count)   = iy
          end if
        end do

        force_local(1:3) = 0.0_wp
        evdw_temp  = 0.0_wp
        elec_temp  = 0.0_wp

!ocl norecurrence(force)
!ocl swp
        do k = 1, num_count
          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)
          iy   = j_list(k)
          lj6  = nonb_lj6 (atmcls(ix,i),atmcls(iy,j))
          lj12 = nonb_lj12(atmcls(ix,i),atmcls(iy,j))

          ! FEP: Determine lamblj and lambel
          fg1 = fepgrp(ix,i)
          fg2 = fepgrp(iy,j)
          lamblj = table_nonb_lambda(1,fg1,fg2)
          lambel = table_nonb_lambda(2,fg1,fg2)

          ! FEP: LJ with soft core for reference
          rij2      = cutoff2*density/(rij2_list(k)+table_nonb_lambda(3,fg1,fg2))
          L         = int(rij2)
          R         = rij2 - L
          L1        = 3*L - 2
          term_lj12 = table_ene(L1) + R*(table_ene(L1+3)-table_ene(L1))
          term_lj6  = table_ene(L1+1) + R*(table_ene(L1+4)-table_ene(L1+1))
          evdw_temp = evdw_temp + lamblj*(term_lj12*lj12 - term_lj6*lj6)
          term_lj12 = table_grad(L1)+ R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1)+ R*(table_grad(L1+4)-table_grad(L1+1))

          ! FEP: EL with soft core for reference
          rij2      = cutoff2*density/(rij2_list(k)+table_nonb_lambda(4,fg1,fg2))
          L         = int(rij2)
          R         = rij2 - L
          L1        = 3*L - 2
          term_elec = table_ene(L1+2) + R*(table_ene(L1+5)-table_ene(L1+2))
          elec_temp = elec_temp + lambel*qtmp*charge(iy,j)*term_elec
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))

          ! gradient
          grad_coef = lamblj*(term_lj12*lj12 - term_lj6*lj6) &
            + lambel*qtmp*charge(iy,j)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)

          force_localj(1,k) = work(1)
          force_localj(2,k) = work(2)
          force_localj(3,k) = work(3)

        end do

        do k = 1, num_count
          iy = j_list(k)
          force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + force_localj(1:3,k)
        end do
        if (check_virial .eq. 1) &
          virial(1:3,ij) = virial(1:3,ij) - force_local(1:3)
        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
        evdw(id+1)  = evdw(id+1)  + evdw_temp
        eelec(id+1) = eelec(id+1) + elec_temp

      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_nonbond_table_linear_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_force_nonbond_table_linear_fep
  !> @brief        calculate nonbonded force without solvents with lookup table
  !!               and soft core for FEP calculations
  !! @authors      NK, HO
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] virial   : virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_table_linear_fep( &
                                                 domain, enefunc, pairlist, &
                                                 coord_pbc, force, virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:)

    ! local variables
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: R, lj12, lj6
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: trans_x, trans_y, trans_z
    real(wp)                  :: rtmp(1:3), qtmp
    real(wp)                  :: force_local(3), force_local_iy(3,MaxAtom)
    real(wp)                  :: dij_list(4,MaxAtom)
    real(wp)                  :: ieps, jeps, eps, irmin, jrmin, rmin
    integer                   :: i, ix, iy, j, k, ij, iix, iwater, list, L, L1
    integer                   :: num_nb15, ini_nb15, fin_nb15, num_count
    integer                   :: id, omp_get_thread_num, j_list(MaxAtom)
    integer                   :: iatmcls,jatmcls
    integer                   :: ncell, ncell_local
    integer                   :: check_virial

    real(dp),         pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: cell_move(:,:,:), system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: density, cutoff
    real(wp),         pointer :: table_grad(:)
    integer,          pointer :: cell_pairlist(:,:), atmcls(:,:)
    integer,          pointer :: natom(:), nsolute(:), nwater(:)
    integer,          pointer :: nb15_list(:,:)
    integer,          pointer :: nb15_cell(:)
    integer,          pointer :: num_nb15_calc1(:,:), num_nb15_calc(:,:)
    integer,          pointer :: nb15_calc_list1(:,:), nb15_calc_list(:,:)
    integer,          pointer :: virial_check(:,:)

    ! FEP
    real(wP)                  :: lamblj, lambel
    integer                   :: fg1, fg2
    integer,          pointer :: nb15_list_fep(:,:)
    integer,          pointer :: nb15_cell_fep(:)
    integer,          pointer :: num_nb15_calc1_fep(:,:), num_nb15_calc_fep(:,:)
    integer,          pointer :: nb15_calc_list1_fep(:,:), nb15_calc_list_fep(:,:)
    integer,          pointer :: fepgrp(:,:)
    real(wp),         pointer :: table_nonb_lambda(:,:,:)

    cell_pairlist   => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    nsolute         => domain%num_solute
    nwater          => domain%num_water
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6 

    nb15_list       => pairlist%nb15_list
    nb15_cell       => pairlist%nb15_cell
    num_nb15_calc1  => pairlist%num_nb15_calc1
    num_nb15_calc   => pairlist%num_nb15_calc 
    nb15_calc_list1 => pairlist%nb15_calc_list1
    nb15_calc_list  => pairlist%nb15_calc_list

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    cutoff2         =  cutoff * cutoff

    ! FEP
    nb15_list_fep       => pairlist%nb15_list_fep
    nb15_cell_fep       => pairlist%nb15_cell_fep
    num_nb15_calc1_fep  => pairlist%num_nb15_calc1_fep
    num_nb15_calc_fep   => pairlist%num_nb15_calc_fep
    nb15_calc_list1_fep => pairlist%nb15_calc_list1_fep
    nb15_calc_list_fep  => pairlist%nb15_calc_list_fep
    fepgrp              => domain%fepgrp
    table_nonb_lambda   => enefunc%table_nonb_lambda

    ! calculate energy and gradient
    !

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, ix, num_nb15, rtmp, qtmp, ini_nb15, fin_nb15, k, iy,  &
    !$omp         rij2, L, R, term_lj12, term_lj6, grad_coef, work, term_elec, &
    !$omp         iwater, ij, j, trans_x, trans_y, trans_z,                    &
    !$omp         iix, list, num_count, j_list, force_local, force_local_iy,   &
    !$omp         lj12, lj6, L1, dij_list, iatmcls, jatmcls, ieps, jeps, eps,  &
    !$omp         irmin, jrmin, rmin, dij, check_virial,                       &
    !$omp         lamblj, lambel, fg1, fg2)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do i = id+1, ncell, nthread
      do ix = 1, natom(i)
        coord_pbc(1,ix,i) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(2,ix,i) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(3,ix,i) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do

    !$omp barrier

    ! Unperturbed region
    !

    ! energy within a cell(solute-solute)
    !
    do i = id+1, ncell_local, nthread

      num_nb15 = 0

      do ix = 1, natom(i) - 1

        rtmp(1:3) = coord_pbc(1:3,ix,i)
        qtmp = charge(ix,i)
        iatmcls = atmcls(ix,i)
        ieps    = lj_coef(1,iatmcls)
        irmin   = lj_coef(2,iatmcls)

        ini_nb15  = num_nb15 + 1
        fin_nb15  = num_nb15 + num_nb15_calc1(ix,i)
        num_nb15  = fin_nb15
        num_count = 0

!ocl norecurrence(force)
!ocl swp
        do k = ini_nb15, fin_nb15
          iy = nb15_calc_list1(k,i)
          dij(1) = rtmp(1) - coord_pbc(1,iy,i)
          dij(2) = rtmp(2) - coord_pbc(2,iy,i)
          dij(3) = rtmp(3) - coord_pbc(3,iy,i)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1,num_count) = dij(1)
            dij_list(2,num_count) = dij(2)
            dij_list(3,num_count) = dij(3)
            dij_list(4,num_count) = rij2
            j_list(num_count) = iy
          end if
        end do

        force_local(1:3) = 0.0_wp
        do k = 1, num_count
          iy   = j_list(k)
          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)
          rij2  = cutoff2 * density / dij_list(4,k)

          lj12 = nonb_lj12(atmcls(ix,i),atmcls(iy,i))
          lj6  = nonb_lj6(atmcls(ix,i),atmcls(iy,i))

          L    = int(rij2)
          R    = rij2 - L

          L1 = 3*L - 2
          term_lj12 = table_grad(L1) + R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6  &
                    + qtmp*charge(iy,i)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force_local_iy(1,k) = work(1)
          force_local_iy(2,k) = work(2)
          force_local_iy(3,k) = work(3)

        end do

        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
        do k = 1, num_count
          iy = j_list(k) 
          force(1:3,iy,i,id+1) = force(1:3,iy,i,id+1) + force_local_iy(1:3,k)
        end do
      end do

    end do

    ! interaction between different cells
    !
    do ij = id+1, maxcell, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)

      num_nb15 = 0

      trans_x = cell_move(1,j,i) * system_size(1)
      trans_y = cell_move(2,j,i) * system_size(2)
      trans_z = cell_move(3,j,i) * system_size(3)
      check_virial = virial_check(j,i)

      do iix = 1, nb15_cell(ij)

        ix = nb15_list(iix,ij)
        rtmp(1:3) = coord_pbc(1:3,ix,i)
        qtmp = charge(ix,i)
        iatmcls = atmcls(ix,i)
        ieps    = lj_coef(1,iatmcls)
        irmin   = lj_coef(2,iatmcls)

        ini_nb15  = num_nb15 + 1
        fin_nb15  = num_nb15 + num_nb15_calc(ix,ij)
        num_nb15  = fin_nb15

        num_count = 0

!ocl norecurrence(force)
!ocl swp
        do k = ini_nb15, fin_nb15
          iy = nb15_calc_list(k,ij)

          dij(1) = rtmp(1) - coord_pbc(1,iy,j) + trans_x
          dij(2) = rtmp(2) - coord_pbc(2,iy,j) + trans_y
          dij(3) = rtmp(3) - coord_pbc(3,iy,j) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1,num_count) = dij(1)
            dij_list(2,num_count) = dij(2)
            dij_list(3,num_count) = dij(3)
            dij_list(4,num_count) = rij2
            j_list(num_count) = iy
          end if
        end do

        force_local(1:3) = 0.0_wp
        do k = 1, num_count
          iy = j_list(k)
          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)
          rij2  = cutoff2 * density / dij_list(4,k)
          lj12 = nonb_lj12(atmcls(ix,i),atmcls(iy,j))
          lj6  = nonb_lj6(atmcls(ix,i),atmcls(iy,j))

          L  = int(rij2)
          R  = rij2 - L
          L1 = 3*L - 2
          term_lj12 = table_grad(L1) + R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6   &
                    + qtmp*charge(iy,j)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force_local_iy(1,k) = work(1)
          force_local_iy(2,k) = work(2)
          force_local_iy(3,k) = work(3)

        end do

        do k = 1, num_count
          iy = j_list(k)
          force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + force_local_iy(1:3,k)
        end do
        if (check_virial == 1) &
          virial(1:3,ij) = virial(1:3,ij) - force_local(1:3)
        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
      end do

    end do

    ! Perturbed region including
    ! singleA, singleB, dualA, and dualB
    !

    ! energy within a cell(solute-solute)
    !
    do i = id+1, ncell_local, nthread

      num_nb15 = 0

      do ix = 1, natom(i) - 1

        rtmp(1:3) = coord_pbc(1:3,ix,i)
        qtmp = charge(ix,i)
        iatmcls = atmcls(ix,i)
        ieps    = lj_coef(1,iatmcls)
        irmin   = lj_coef(2,iatmcls)

        ini_nb15  = num_nb15 + 1
        fin_nb15  = num_nb15 + num_nb15_calc1_fep(ix,i)
        num_nb15  = fin_nb15
        num_count = 0

!ocl norecurrence(force)
!ocl swp
        do k = ini_nb15, fin_nb15
          iy = nb15_calc_list1_fep(k,i)
          dij(1) = rtmp(1) - coord_pbc(1,iy,i)
          dij(2) = rtmp(2) - coord_pbc(2,iy,i)
          dij(3) = rtmp(3) - coord_pbc(3,iy,i)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1,num_count) = dij(1)
            dij_list(2,num_count) = dij(2)
            dij_list(3,num_count) = dij(3)
            dij_list(4,num_count) = rij2
            j_list(num_count) = iy
          end if
        end do

        force_local(1:3) = 0.0_wp
        do k = 1, num_count
          iy   = j_list(k)
          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)

          lj12 = nonb_lj12(atmcls(ix,i),atmcls(iy,i))
          lj6  = nonb_lj6(atmcls(ix,i),atmcls(iy,i))

          ! FEP: Determine lamblj and lambel
          fg1 = fepgrp(ix,i)
          fg2 = fepgrp(iy,i)
          lamblj = table_nonb_lambda(1,fg1,fg2)
          lambel = table_nonb_lambda(2,fg1,fg2)

          ! FEP: LJ with soft-core for reference
          rij2      = cutoff2 * density / (dij_list(4,k) + table_nonb_lambda(3,fg1,fg2))
          L         = int(rij2)
          R         = rij2 - L
          L1        = 3*L - 2
          term_lj12 = table_grad(L1)+ R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1)+ R*(table_grad(L1+4)-table_grad(L1+1))

          ! FEP: EL with soft-core for reference
          rij2      = cutoff2 * density / (dij_list(4,k) + table_nonb_lambda(4,fg1,fg2))
          L         = int(rij2)
          R         = rij2 - L
          L1        = 3*L - 2
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))

          ! Gradient for reference
          grad_coef = lamblj * (term_lj12*lj12 - term_lj6*lj6)  &
                  + lambel * qtmp*charge(iy,i)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force_local_iy(1,k) = work(1)
          force_local_iy(2,k) = work(2)
          force_local_iy(3,k) = work(3)

        end do

        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
        do k = 1, num_count
          iy = j_list(k) 
          force(1:3,iy,i,id+1) = force(1:3,iy,i,id+1) + force_local_iy(1:3,k)
        end do
      end do

    end do

    ! interaction between different cells
    !
    do ij = id+1, maxcell, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)

      num_nb15 = 0

      trans_x = cell_move(1,j,i) * system_size(1)
      trans_y = cell_move(2,j,i) * system_size(2)
      trans_z = cell_move(3,j,i) * system_size(3)
      check_virial = virial_check(j,i)

      do iix = 1, nb15_cell_fep(ij)

        ix = nb15_list_fep(iix,ij)
        rtmp(1:3) = coord_pbc(1:3,ix,i)
        qtmp = charge(ix,i)
        iatmcls = atmcls(ix,i)
        ieps    = lj_coef(1,iatmcls)
        irmin   = lj_coef(2,iatmcls)

        ini_nb15  = num_nb15 + 1
        fin_nb15  = num_nb15 + num_nb15_calc_fep(ix,ij)
        num_nb15  = fin_nb15

        num_count = 0

!ocl norecurrence(force)
!ocl swp
        do k = ini_nb15, fin_nb15
          iy = nb15_calc_list_fep(k,ij)

          dij(1) = rtmp(1) - coord_pbc(1,iy,j) + trans_x
          dij(2) = rtmp(2) - coord_pbc(2,iy,j) + trans_y
          dij(3) = rtmp(3) - coord_pbc(3,iy,j) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          if (rij2 < cutoff2) then
            num_count = num_count + 1
            dij_list(1,num_count) = dij(1)
            dij_list(2,num_count) = dij(2)
            dij_list(3,num_count) = dij(3)
            dij_list(4,num_count) = rij2
            j_list(num_count) = iy
          end if
        end do

        force_local(1:3) = 0.0_wp
        do k = 1, num_count
          iy = j_list(k)
          dij(1) = dij_list(1,k)
          dij(2) = dij_list(2,k)
          dij(3) = dij_list(3,k)
          lj12 = nonb_lj12(atmcls(ix,i),atmcls(iy,j))
          lj6  = nonb_lj6(atmcls(ix,i),atmcls(iy,j))

          ! FEP: Determine lamblj and lambel
          fg1 = fepgrp(ix,i)
          fg2 = fepgrp(iy,j)
          lamblj = table_nonb_lambda(1,fg1,fg2)
          lambel = table_nonb_lambda(2,fg1,fg2)

          ! FEP: LJ with soft-core for reference
          rij2      = cutoff2 * density / (dij_list(4,k) + table_nonb_lambda(3,fg1,fg2))
          L         = int(rij2)
          R         = rij2 - L
          L1        = 3*L - 2
          term_lj12 = table_grad(L1)+ R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1)+ R*(table_grad(L1+4)-table_grad(L1+1))

          ! FEP: EL with soft-core for reference
          rij2      = cutoff2 * density / (dij_list(4,k) + table_nonb_lambda(4,fg1,fg2))
          L         = int(rij2)
          R         = rij2 - L
          L1        = 3*L - 2
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))

          ! Gradient for reference
          grad_coef = lamblj * (term_lj12*lj12 - term_lj6*lj6)   &
                  + lambel * qtmp*charge(iy,j)*term_elec

          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force_local_iy(1,k) = work(1)
          force_local_iy(2,k) = work(2)
          force_local_iy(3,k) = work(3)

        end do

        do k = 1, num_count
          iy = j_list(k)
          force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + force_local_iy(1:3,k)
        end do
        if (check_virial == 1) &
          virial(1:3,ij) = virial(1:3,ij) - force_local(1:3)
        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_force_nonbond_table_linear_fep

#ifdef USE_GPU
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_table_linear_gpu_fep(domain, enefunc, pairlist,&
                                                     npt, coord_pbc, force,    &
                                                     virial, eelec, evdw,      &
                                                     ene_virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    logical,                  intent(in)    :: npt
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)
    real(dp),                 intent(inout) :: eelec(nthread)
    real(dp),                 intent(inout) :: evdw(nthread)
    real(dp),                 intent(inout) :: ene_virial(:)

    ! local variables
    real(wp)                  :: cutoff2, pairlistdist2
    integer                   :: ncell, ncell_local, ncell_bound

    real(dp),         pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: cell_move(:,:,:), system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: density, cutoff, pairlistdist
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer,          pointer :: atmcls(:,:)
    integer,          pointer :: natom(:)
    integer,          pointer :: virial_check(:,:)

    ! for GPU
    integer,          pointer :: univ_cell_pairlist1(:,:)
    integer,          pointer :: univ_cell_pairlist2(:,:)
    integer,          pointer :: univ_ix_natom(:)
    integer(1),       pointer :: univ_ix_list(:,:)
    integer,          pointer :: univ_iy_natom(:)
    integer(1),       pointer :: univ_iy_list(:,:)
    integer,          pointer :: univ_ij_sort_list(:)
    integer(1),       pointer :: univ_mask2(:,:)
    integer,          pointer :: univ_mask2_size
    integer                   :: univ_ncell_nonzero
    integer                   :: univ_update
      
    integer                   :: ncell_max, ij, i, j
    integer                   :: cutoff_int
    integer                   :: num_atom_cls
    integer                   :: ret_shape(1)
    integer                   :: univ_ij, univ_ij0, iiy, idx, base_idx

    ! FEP
    integer,          pointer :: fepgrp(:,:)
    real(wp),         pointer :: table_nonb_lambda(:,:,:)

    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_ene       => enefunc%table%table_ene
    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6
    pairlistdist    => enefunc%pairlistdist

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    ncell_bound     =  domain%num_cell_boundary
    cutoff2         =  cutoff * cutoff
    pairlistdist2   =  pairlistdist * pairlistdist

    ret_shape = shape(natom)
    ncell_max = ret_shape(1)
    cutoff_int = enefunc%table%cutoff_int
    num_atom_cls = enefunc%num_atom_cls

    !
    ! univ
    univ_cell_pairlist1 => domain%univ_cell_pairlist1
    univ_cell_pairlist2 => domain%univ_cell_pairlist2
    univ_ix_natom       => pairlist%univ_ix_natom
    univ_ix_list        => pairlist%univ_ix_list
    univ_iy_natom       => pairlist%univ_iy_natom
    univ_iy_list        => pairlist%univ_iy_list
    univ_ij_sort_list   => pairlist%univ_ij_sort_list
    univ_mask2          => pairlist%univ_mask2
    univ_mask2_size     => pairlist%univ_mask2_size
    univ_ncell_nonzero  =  pairlist%univ_ncell_nonzero
    univ_update         =  pairlist%univ_update

    ! FEP
    fepgrp              => domain%fepgrp
    table_nonb_lambda   => enefunc%table_nonb_lambda

    !
    ! launch GPU kernels
    call gpu_launch_compute_energy_nonbond_table_linear_univ_fep( &
         coord_pbc, force(:,:,:,1), ene_virial,               &
         coord, trans1, cell_move, charge, atmcls, fepgrp, natom,       &
         nonb_lj12, nonb_lj6, table_ene, table_grad, table_nonb_lambda, &
         univ_cell_pairlist1, univ_mask2,                     &
         univ_ix_natom, univ_ix_list, univ_iy_natom,          &
         univ_iy_list, univ_ij_sort_list, virial_check,       &
         MaxAtom, MaxAtomCls, num_atom_cls,                   &
         ncell_local, ncell_bound, ncell_max,                 &
         cutoff_int, univ_maxcell, univ_maxcell1,             &
         univ_ncell_nonzero, univ_ncell_near, univ_update,    &
         univ_mask2_size, univ_natom_max, maxcell, density,   &
         cutoff2,                                             &
         system_size(1), system_size(2), system_size(3) )


    return

  end subroutine compute_energy_nonbond_table_linear_gpu_fep

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_force_nonbond_table_linear_gpu_fep(domain, enefunc, pairlist, &
                                                    npt, cpu_calc, coord_pbc,  &
                                                    force_omp, force, virial,  &
                                                    ene_virial)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    logical,                  intent(in)    :: npt
    logical,                  intent(in)    :: cpu_calc
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force_omp(:,:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: virial(:,:,:)
    real(dp),                 intent(inout) :: ene_virial(:)

    ! local variables
    real(wp)                  :: cutoff2, pairlistdist2
    integer                   :: ncell, ncell_local, ncell_bound
    ! integer                   :: check_virial

    real(dp),         pointer :: coord(:,:,:)
    real(wp),         pointer :: trans1(:,:,:), charge(:,:)
    real(wp),         pointer :: cell_move(:,:,:), system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: density, cutoff, pairlistdist
    real(wp),         pointer :: table_grad(:)
    integer,          pointer :: cell_pairlist(:,:), atmcls(:,:)
    integer,          pointer :: natom(:)
    integer,          pointer :: virial_check(:,:)

    ! for GPU
    integer,          pointer :: univ_cell_pairlist1(:,:)
    integer,          pointer :: univ_cell_pairlist2(:,:)
    integer,          pointer :: univ_ix_natom(:)
    integer(1),       pointer :: univ_ix_list(:,:)
    integer,          pointer :: univ_iy_natom(:)
    integer(1),       pointer :: univ_iy_list(:,:)
    integer,          pointer :: univ_ij_sort_list(:)
    integer(1),       pointer :: univ_mask2(:,:)
    integer,          pointer :: univ_mask2_size
    integer                   :: univ_ncell_nonzero, univ_gpu_start_index
    integer                   :: univ_update
    integer                   :: ncell_max
    integer                   :: cutoff_int
    integer                   :: num_atom_cls
    integer                   :: ret_shape(1)
    integer                   :: univ_ij, univ_ij0, iiy, idx, base_idx

    ! FEP
    integer,          pointer :: fepgrp(:,:)
    real(wp),         pointer :: table_nonb_lambda(:,:,:)

    cell_pairlist   => domain%cell_pairlist1
    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    charge          => domain%charge
    cell_move       => domain%cell_move
    system_size     => domain%system_size
    virial_check    => domain%virial_check

    table_grad      => enefunc%table%table_grad
    density         => enefunc%table%density
    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12
    nonb_lj6        => enefunc%nonb_lj6
    pairlistdist    => enefunc%pairlistdist

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    ncell_local     =  domain%num_cell_local
    ncell_bound     =  domain%num_cell_boundary
    cutoff2         =  cutoff * cutoff
    pairlistdist2   =  pairlistdist * pairlistdist

    ret_shape = shape(natom)
    ncell_max = ret_shape(1)
    cutoff_int = enefunc%table%cutoff_int
    num_atom_cls = enefunc%num_atom_cls

    ! FEP
    fepgrp            => domain%fepgrp
    table_nonb_lambda => enefunc%table_nonb_lambda

    !
    ! univ
    univ_cell_pairlist1 => domain%univ_cell_pairlist1
    univ_cell_pairlist2 => domain%univ_cell_pairlist2
    univ_ix_natom       => pairlist%univ_ix_natom
    univ_ix_list        => pairlist%univ_ix_list
    univ_iy_natom       => pairlist%univ_iy_natom
    univ_iy_list        => pairlist%univ_iy_list
    univ_ij_sort_list   => pairlist%univ_ij_sort_list
    univ_mask2          => pairlist%univ_mask2
    univ_mask2_size     => pairlist%univ_mask2_size
    univ_ncell_nonzero  =  pairlist%univ_ncell_nonzero
    univ_update         =  pairlist%univ_update

    !  launch GPU kernels (data transfer from CPU to GPU as well)
    !
    if (cpu_calc) then
      univ_gpu_start_index = univ_gpu_start
    else
      univ_gpu_start = ncell_local
    end if
    call gpu_launch_compute_force_nonbond_table_linear_univ_fep( &
         coord_pbc, force(:,:,:,1), ene_virial,                  &  ! argments
         coord, trans1, cell_move, charge, atmcls, fepgrp, natom,&  ! arrays
         nonb_lj12, nonb_lj6, table_grad, table_nonb_lambda,     &  ! arrays
         univ_cell_pairlist1, univ_mask2, univ_ix_natom,     &
         univ_ix_list, univ_iy_natom, univ_iy_list,          &
         univ_ij_sort_list, virial_check,                    &
         MaxAtom, MaxAtomCls, num_atom_cls,                  &  ! size
         ncell_local, ncell_bound, ncell_max, cutoff_int,    &  ! size
         univ_maxcell, univ_maxcell1, univ_ncell_nonzero,    &
         univ_ncell_near, univ_update, univ_mask2_size,      &
         univ_natom_max, npt, cpu_calc, density, cutoff2,    &
         pairlistdist2, univ_gpu_start,                      &
         system_size(1), system_size(2), system_size(3))

    if (cpu_calc) then
      call timer(TimerTest10, TimerOn)
      call cpu_compute_force_intra_cell_univ_fep(            &
         coord, trans1, cell_move, system_size, charge,      &
         atmcls, fepgrp, natom, nonb_lj12, nonb_lj6, table_grad, table_nonb_lambda,    &  
         univ_cell_pairlist1, univ_mask2, univ_ix_natom,     &
         univ_ix_list, univ_iy_natom, univ_iy_list,          &
         univ_ij_sort_list,                                  &
         ncell_local, ncell_bound, univ_ncell_nonzero,       &
         univ_update, density, cutoff2,                      &
         coord_pbc, force_omp)
      call timer(TimerTest10, TimerOff)
    end if

    return

  end subroutine compute_force_nonbond_table_linear_gpu_fep

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine cpu_compute_force_intra_cell_univ_fep(coord, trans1, cell_move,   &
             system_size, charge, atmcls, fepgrp, natom, nonb_lj12, nonb_lj6,      &
             table_grad, table_nonb_lambda, univ_cell_pairlist1, univ_mask2,                  &
             univ_ix_natom, univ_ix_list, univ_iy_natom, univ_iy_list,     &
             univ_ij_sort_list, ncell_local, ncell_bound,                  &
             univ_ncell_nonzero, univ_update, density, cutoff2,            &
             coord_pbc, force)

    ! formal arguments
    real(dp),                 intent(in)    :: coord(:,:,:)
    real(wp),                 intent(in)    :: trans1(:,:,:)
    real(wp),                 intent(in)    :: cell_move(:,:,:)
    real(wp),                 intent(in)    :: system_size(:)
    real(wp),                 intent(in)    :: charge(:,:)
    integer,                  intent(in)    :: atmcls(:,:)
    integer,                  intent(in)    :: fepgrp(:,:)
    integer,                  intent(in)    :: natom(:)
    real(wp),                 intent(in)    :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),                 intent(in)    :: table_grad(:)
    real(wp),                 intent(in)    :: table_nonb_lambda(:,:,:)
    integer,                  intent(in)    :: univ_cell_pairlist1(:,:)
    integer(1),               intent(in)    :: univ_mask2(:,:)
    integer,                  intent(in)    :: univ_ix_natom(:)
    integer(1),               intent(in)    :: univ_ix_list(:,:)
    integer,                  intent(in)    :: univ_iy_natom(:)
    integer(1),               intent(in)    :: univ_iy_list(:,:)
    integer,                  intent(in)    :: univ_ij_sort_list(:)
    integer,                  intent(in)    :: ncell_local, ncell_bound
    integer,                  intent(in)    :: univ_ncell_nonzero
    integer,                  intent(in)    :: univ_update 
    real(wp),                 intent(in)    :: cutoff2, density
    real(wp),                 intent(inout) :: coord_pbc(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)


    ! local variables
    real(wp)                  :: dij_list(1:3,MaxAtom), rij2_list(MaxAtom)
    real(wp)                  :: cell_move_local(1:3), dij(1:3), rij2, inv_r2
    real(wp)                  :: R, lj12, lj6
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: work(1:3), grad_coef
    real(wp)                  :: rtmp(1:3), qtmp
    real(wp)                  :: force_local(3), force_j(1:3,MaxAtom)
    integer                   :: j_list(MaxAtom)
    integer                   :: index, univ_ij
    integer                   :: i, j, ix, iy, iix, iiy, k, L, L1, idx
    integer                   :: ix_natom, iy_natom
    integer                   :: id, omp_get_thread_num
    integer                   :: iatmcls,jatmcls
    integer                   :: num_count

    ! FEP
    integer                   :: fg1, fg2
    real(wp)                  :: lamblj, lambel

    !$omp parallel default(shared)                                            &
    !$omp private(id, i, j, ix, iy, iix, iiy, k, univ_ij, index, idx, j_list, &
    !$omp         ix_natom, iy_natom, num_count, L, L1, iatmcls, jatmcls,     &
    !$omp         cell_move_local, rtmp, qtmp, rij2, R, term_lj12, term_lj6,  &
    !$omp         term_elec, grad_coef, work, dij, inv_r2, lj12, lj6,         &
    !$omp         dij_list, rij2_list, force_local, force_j, fg1, fg2,        &
    !$omp         lamblj, lambel)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell_local+ncell_bound, nthread
      do ix = 1, natom(i)
        coord_pbc(1:3,ix,i) = coord(1:3,ix,i) + trans1(1:3,ix,i)
      end do
    end do

    !$omp barrier

    do i = id+1, ncell_local, nthread

      do ix = 1, natom(i) - 1

        rtmp(1:3) = coord_pbc(1:3,ix,i)
        qtmp      = charge(ix,i)
        iatmcls   = atmcls(ix,i)

        force_local(1:3) = 0.0_wp
        num_count        = 0

        ! FEP
        fg1 = fepgrp(ix,i)
!ocl norecurrence(force)
!ocl swp
        do iy = ix + 1, natom(i)

          ! FEP
          fg2 = fepgrp(iy,i)

          idx = iy + (ix-1)*univ_natom_max

          ! FEP
          if ((univ_mask2(idx,i) == 1) .and. &
            (int(table_nonb_lambda(5,fg1,fg2)) == 1)) then

            dij(1:3) = rtmp(1:3) - coord_pbc(1:3,iy,i)
            rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            if (rij2 < cutoff2) then
              num_count = num_count + 1
              dij_list(1:3,num_count)  = dij(1:3)
              rij2_list(num_count) = rij2
              j_list(num_count)   = iy
            end if
          end if
        end do

        do k = 1, num_count

          dij(1:3)  = dij_list(1:3,k)
          iy      = j_list(k)
          jatmcls = atmcls(iy,i)
          lj6     = nonb_lj6(iatmcls,jatmcls)
          lj12    = nonb_lj12(iatmcls,jatmcls)

          ! FEP: Determine lamblj and lambel
          fg2 = fepgrp(iy,i)
          lamblj = table_nonb_lambda(1,fg1,fg2)
          lambel = table_nonb_lambda(2,fg1,fg2)

          rij2    = cutoff2*density/(rij2_list(k)+table_nonb_lambda(3,fg1,fg2))
          L  = int(rij2)
          R  = rij2 - L
          L1 = 3*L - 2
          term_lj12 = table_grad(L1  ) + R*(table_grad(L1+3)-table_grad(L1  ))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))

          rij2    = cutoff2*density/(rij2_list(k)+table_nonb_lambda(4,fg1,fg2))
          L  = int(rij2)
          R  = rij2 - L
          L1 = 3*L - 2
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))

          grad_coef = lamblj*(term_lj12*lj12 - term_lj6*lj6)   &
                    + lambel*qtmp*charge(iy,i)*term_elec

          work(1:3) = grad_coef*dij(1:3)

          ! store force
          !
          force_local(1:3) = force_local(1:3) - work(1:3)
          force_j(1:3,k) = work(1:3)

        end do

        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
        do k = 1, num_count
          iy = j_list(k)
          force(1:3,iy,i,id+1) = force(1:3,iy,i,id+1) + force_j(1:3,k)
        end do

      end do

    end do

    do index = ncell_local+id+1, univ_gpu_start, nthread

      univ_ij = univ_ij_sort_list(index)
      i = univ_cell_pairlist1(1,univ_ij)
      j = univ_cell_pairlist1(2,univ_ij)
!     ix_natom = univ_ix_natom(univ_ij)
!     iy_natom = univ_iy_natom(univ_ij)
      cell_move_local(1:3) = cell_move(1:3,j,i) * system_size(1:3)

      if (univ_ij > univ_ncell_near) then

!       do iix = 1, ix_natom
        do ix = 1, natom(i)

!         ix = univ_ix_list(iix,univ_ij)
          rtmp(1:3) = coord_pbc(1:3,ix,i) + cell_move_local(1:3)
          qtmp      = charge(ix,i)
          iatmcls   = atmcls(ix,i)

          force_local(1:3) = 0.0_wp
          num_count = 0

          ! FEP
          fg1 = fepgrp(ix,i)

!         do iiy = 1, iy_natom
          do iy = 1, natom(j)
!           iy = univ_iy_list(iiy,univ_ij)
            ! FEP
            fg2 = fepgrp(iy,j)
            if (int(table_nonb_lambda(5,fg1,fg2)) == 1) then
              dij(1:3) = rtmp(1:3) - coord_pbc(1:3,iy,j)  
              rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
              if (rij2 < cutoff2) then
                num_count = num_count + 1
                dij_list(1:3,num_count)  = dij(1:3)
                rij2_list(num_count) = rij2
                j_list(num_count)   = iy
              end if
            end if
          end do

          do k = 1, num_count

            dij(1)  = dij_list(1,k)
            dij(2)  = dij_list(2,k)
            dij(3)  = dij_list(3,k)
            iy      = j_list(k)
            jatmcls = atmcls(iy,j)
            lj6     = nonb_lj6(iatmcls,jatmcls)
            lj12    = nonb_lj12(iatmcls,jatmcls)

            ! FEP: Determine lamblj and lambel
            fg2 = fepgrp(iy,j)
            lamblj = table_nonb_lambda(1,fg1,fg2)
            lambel = table_nonb_lambda(2,fg1,fg2)

            rij2    = cutoff2*density/(rij2_list(k)+table_nonb_lambda(3,fg1,fg2))
            L  = int(rij2)
            R  = rij2 - L
            L1 = 3*L - 2
            term_lj12 = table_grad(L1  ) + R*(table_grad(L1+3)-table_grad(L1))
            term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))

            rij2    = cutoff2*density/(rij2_list(k)+table_nonb_lambda(4,fg1,fg2))
            L  = int(rij2)
            R  = rij2 - L
            L1 = 3*L - 2
            term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))

            grad_coef = lamblj*(term_lj12*lj12 - term_lj6*lj6)   &
                      + lambel*qtmp*charge(iy,j)*term_elec
            work(1:3) = grad_coef*dij(1:3)

            ! store force
            !
            force_local(1:3) = force_local(1:3) - work(1:3)
            force_j(1:3,k)   = work(1:3)

          end do

          force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
          do k = 1, num_count
            iy = j_list(k)
            force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + force_j(1:3,k)
          end do

        end do

      else

!       do iix = 1, ix_natom
        do ix = 1, natom(i)

!         ix = univ_ix_list(iix,univ_ij)
          rtmp(1:3) = coord_pbc(1:3,ix,i) + cell_move_local(1:3)
          qtmp      = charge(ix,i)
          iatmcls   = atmcls(ix,i)

          force_local(1:3) = 0.0_wp
          num_count = 0

          ! FEP
          fg1 = fepgrp(ix,i)

!         do iiy = 1, iy_natom
          do iy = 1, natom(j)
!           iy = univ_iy_list(iiy,univ_ij)
            ! FEP
            fg2 = fepgrp(iy,j)
            idx = iy + (ix-1)*univ_natom_max
            if ((univ_mask2(idx,univ_ij) == 1) .and. &
              (int(table_nonb_lambda(5,fg1,fg2)) == 1)) then
              dij(1:3) = rtmp(1:3) - coord_pbc(1:3,iy,j)
              rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
              if (rij2 < cutoff2) then
                num_count = num_count + 1
                dij_list(1:3,num_count)  = dij(1:3)
                rij2_list(num_count) = rij2
                j_list(num_count)   = iy
              end if
            end if
          end do

          do k = 1, num_count

            dij(1)  = dij_list(1,k)
            dij(2)  = dij_list(2,k)
            dij(3)  = dij_list(3,k)
            iy      = j_list(k)
            jatmcls = atmcls(iy,j)
            lj6     = nonb_lj6(iatmcls,jatmcls)
            lj12    = nonb_lj12(iatmcls,jatmcls)

            ! FEP
            ! Determine lamblj and lambel
            fg2 = fepgrp(iy,j)
            lamblj = table_nonb_lambda(1,fg1,fg2)
            lambel = table_nonb_lambda(2,fg1,fg2)

            rij2    = cutoff2*density/(rij2_list(k)+table_nonb_lambda(3,fg1,fg2))
            L  = int(rij2)
            R  = rij2 - L
            L1 = 3*L - 2
            term_lj12 = table_grad(L1  ) + R*(table_grad(L1+3)-table_grad(L1))
            term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))

            rij2    = cutoff2*density/(rij2_list(k)+table_nonb_lambda(4,fg1,fg2))
            L  = int(rij2)
            R  = rij2 - L
            L1 = 3*L - 2
            term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))

            grad_coef = lamblj*(term_lj12*lj12 - term_lj6*lj6)   &
                      + lambel*qtmp*charge(iy,j)*term_elec
            work(1:3) = grad_coef*dij(1:3)

            ! store force
            !
            force_local(1:3) = force_local(1:3) - work(1:3)
            force_j(1:3,k)   = work(1:3)

          end do

          force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
          do k = 1, num_count
            iy = j_list(k)
            force(1:3,iy,j,id+1) = force(1:3,iy,j,id+1) + force_j(1:3,k)
          end do

        end do

      end if

    end do

    !$omp end parallel 

    return

  end subroutine cpu_compute_force_intra_cell_univ_fep

#endif

end module sp_energy_table_linear_nowater_mod
