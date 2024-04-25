!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_fep_energy_mod
!> @brief   calculate energy differnces between adjacent states in FEP 
!! @authors Hiraku Oshima (HO)
!
!  (c) Copyright 2019 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_fep_energy_mod

  use sp_energy_pme_mod
  use sp_energy_mod          
  use sp_dynvars_str_mod
  use sp_boundary_str_mod
  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_energy_str_mod
  use sp_ensemble_str_mod
  use sp_remd_str_mod
  use sp_alchemy_str_mod
  use sp_fep_utils_mod
  use sp_ensemble_str_mod
  use sp_domain_str_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use math_libs_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: compute_fep_energy
  public  :: assign_condition_feprest

  private :: compute_energy_ref
  private :: compute_energy_tgt
  private :: assign_condition_fep
  private :: assign_condition_feprest_angle_ub
  private :: assign_condition_feprest_lj
  private :: assign_condition_feprest_internal

contains

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_fep_energy(domain, enefunc, dynvars, &
                               pairlist, ensemble, boundary,      &
                               alchemy, remd)

    ! formal arguments
    type(s_domain),    target, intent(inout) :: domain
    type(s_enefunc),           intent(inout) :: enefunc
    type(s_dynvars),   target, intent(inout) :: dynvars
    type(s_pairlist),          intent(in)    :: pairlist
    type(s_ensemble),          intent(in)    :: ensemble
    type(s_boundary),          intent(inout) :: boundary
    type(s_alchemy),   target, intent(inout) :: alchemy
    type(s_remd),              intent(inout) :: remd

    ! local variables
    integer,         pointer :: lambid
    integer                  :: lambid_forward, lambid_backward
    integer                  :: replica_id
    integer                  :: param_id, param_id_forward, param_id_backward
    integer                  :: i, ix, j, k
    type(s_energy)           :: energy_ref
    type(s_energy)           :: energy_tgt

    lambid    => alchemy%lambid

    if (alchemy%fep_md_type == FEP_Parallel) then
      ! FEP-REMD case

      ! energy calculation by changing lambdas
      replica_id = my_country_no + 1
      param_id = remd%repid2parmsetid(replica_id)

      call compute_energy_ref(domain, enefunc, boundary, dynvars%energy, energy_ref)
      dynvars%energy%deltU_fep(1) = energy_ref%total

      ! evalutate total energy for target state
      ! and calculate the difference from reference state

      ! shift to the backward state
      param_id_backward = param_id - 1
      if (param_id_backward > 0) then
        call assign_condition_fep(param_id_backward, alchemy, domain,&
                                  ensemble, enefunc, remd)
        call compute_energy_tgt(domain, enefunc, pairlist, boundary, energy_tgt)
        dynvars%energy%deltU_fep(2)=energy_tgt%total-energy_ref%total
      end if

      ! shift to the forward state
      param_id_forward  = param_id + 1
      if (param_id_forward <= remd%nreplicas(1)) then
        call assign_condition_fep(param_id_forward, alchemy, domain, &
                                  ensemble, enefunc, remd)
        call compute_energy_tgt(domain, enefunc, pairlist, boundary, energy_tgt)
        dynvars%energy%deltU_fep(3)=energy_tgt%total-energy_ref%total
      end if

      ! return to the original state
      call assign_condition_fep(param_id, alchemy, domain, ensemble, &
                                enefunc, remd)
    else
      ! FEP-MD case

      ! evaluate total energy for reference state
      call compute_energy_ref(domain, enefunc, boundary, dynvars%energy, energy_ref)
      dynvars%energy%deltU_fep(1) = energy_ref%total

      ! evalutate total energy for target state
      ! and calculate the difference from reference state
      if (alchemy%fep_direction == FEP_Bothsides) then

        if (lambid == 1) then
          lambid_backward = lambid
          lambid_forward  = lambid + 1
        else if (lambid == alchemy%num_fep_windows) then
          lambid_backward = lambid - 1
          lambid_forward  = lambid
        else
          lambid_backward = lambid - 1
          lambid_forward  = lambid + 1
        end if

        enefunc%lambljA      = alchemy%lambljA(lambid_backward)
        enefunc%lambljB      = alchemy%lambljB(lambid_backward)
        enefunc%lambelA      = alchemy%lambelA(lambid_backward)
        enefunc%lambelB      = alchemy%lambelB(lambid_backward)
        enefunc%lambbondA    = alchemy%lambbondA(lambid_backward)
        enefunc%lambbondB    = alchemy%lambbondB(lambid_backward)
        enefunc%lambrest     = alchemy%lambrest(lambid_backward)
        ! Make table of lambda and softcore in FEP
        call set_lambda_table_fep(enefunc)
        call compute_energy_tgt(domain, enefunc, pairlist, boundary, energy_tgt)
        dynvars%energy%deltU_fep(2)=energy_tgt%total-energy_ref%total

        enefunc%lambljA      = alchemy%lambljA(lambid_forward)
        enefunc%lambljB      = alchemy%lambljB(lambid_forward)
        enefunc%lambelA      = alchemy%lambelA(lambid_forward)
        enefunc%lambelB      = alchemy%lambelB(lambid_forward)
        enefunc%lambbondA    = alchemy%lambbondA(lambid_forward)
        enefunc%lambbondB    = alchemy%lambbondB(lambid_forward)
        enefunc%lambrest     = alchemy%lambrest(lambid_forward)
        ! Make table of lambda and softcore in FEP
        call set_lambda_table_fep(enefunc)
        call compute_energy_tgt(domain, enefunc, pairlist, boundary, energy_tgt)
        dynvars%energy%deltU_fep(3)=energy_tgt%total-energy_ref%total

      else if (alchemy%fep_direction == FEP_Forward) then

        lambid_forward       = lambid + 1
        enefunc%lambljA      = alchemy%lambljA(lambid_forward)
        enefunc%lambljB      = alchemy%lambljB(lambid_forward)
        enefunc%lambelA      = alchemy%lambelA(lambid_forward)
        enefunc%lambelB      = alchemy%lambelB(lambid_forward)
        enefunc%lambbondA    = alchemy%lambbondA(lambid_forward)
        enefunc%lambbondB    = alchemy%lambbondB(lambid_forward)
        enefunc%lambrest     = alchemy%lambrest(lambid_forward)
        ! Make table of lambda and softcore in FEP
        call set_lambda_table_fep(enefunc)
        call compute_energy_tgt(domain, enefunc, pairlist, boundary, energy_tgt)
        dynvars%energy%deltU_fep(2)=energy_tgt%total-energy_ref%total

      else if (alchemy%fep_direction == FEP_Reverse) then

        lambid_backward      = lambid - 1
        enefunc%lambljA      = alchemy%lambljA(lambid_backward)
        enefunc%lambljB      = alchemy%lambljB(lambid_backward)
        enefunc%lambelA      = alchemy%lambelA(lambid_backward)
        enefunc%lambelB      = alchemy%lambelB(lambid_backward)
        enefunc%lambbondA    = alchemy%lambbondA(lambid_backward)
        enefunc%lambbondB    = alchemy%lambbondB(lambid_backward)
        enefunc%lambrest     = alchemy%lambrest(lambid_backward)
        ! Make table of lambda and softcore in FEP
        call set_lambda_table_fep(enefunc)
        call compute_energy_tgt(domain, enefunc, pairlist, boundary, energy_tgt)
        dynvars%energy%deltU_fep(2)=energy_tgt%total-energy_ref%total

      end if

      ! return to the original state
      enefunc%lambljA      = alchemy%lambljA(lambid)
      enefunc%lambljB      = alchemy%lambljB(lambid)
      enefunc%lambelA      = alchemy%lambelA(lambid)
      enefunc%lambelB      = alchemy%lambelB(lambid)
      enefunc%lambbondA    = alchemy%lambbondA(lambid)
      enefunc%lambbondB    = alchemy%lambbondB(lambid)
      enefunc%lambrest     = alchemy%lambrest(lambid)
      ! Make table of lambda and softcore in FEP
      call set_lambda_table_fep(enefunc)

    end if

  end subroutine compute_fep_energy

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_ref(domain, enefunc, boundary, energy, energy_ref)

    ! formal arguments
    type(s_domain),   target, intent(inout) :: domain
    type(s_enefunc),  target, intent(inout) :: enefunc
    type(s_boundary),         intent(in)    :: boundary
    type(s_energy),           intent(inout) :: energy
    type(s_energy),           intent(inout) :: energy_ref

    ! local variables
    real(dp)                 :: volume

    call mpi_allreduce(energy%total, energy_ref%total, 1, mpi_real8,  &
                       mpi_sum, mpi_comm_country, ierror)

    energy_ref%total = energy_ref%total + &
               energy%restraint_distance + &
               energy%restraint_rmsd     + &
               energy%restraint_emfit    + &
               energy%contact            + &
               energy%noncontact

    ! Dispersion correction
    if (enefunc%dispersion_corr /= Disp_corr_NONE) then
      volume =  boundary%box_size_x_ref * &
                boundary%box_size_y_ref * &
                boundary%box_size_z_ref
      energy_ref%disp_corr_energy = enefunc%dispersion_energy_preserve &
        + enefunc%lambljA*enefunc%dispersion_energy_vanish &
        + enefunc%lambljB*enefunc%dispersion_energy_appear
      energy_ref%disp_corr_energy = energy_ref%disp_corr_energy / volume
      energy_ref%total = energy_ref%total + energy_ref%disp_corr_energy
    end if

    return

  end subroutine compute_energy_ref

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_tgt(domain, enefunc, pairlist, boundary, energy)

    ! formal arguments
    type(s_domain),   target, intent(inout) :: domain
    type(s_enefunc),  target, intent(inout) :: enefunc
    type(s_pairlist),         intent(in)    :: pairlist
    type(s_boundary),         intent(in)    :: boundary
    type(s_energy),           intent(inout) :: energy

    integer                           :: nc, nv
    real(dp),       allocatable, save :: force(:,:,:)
    real(dp),       allocatable, save :: force_long(:,:,:)
    real(wp),       allocatable, save :: force_omp(:,:,:,:)
    real(wp),       allocatable, save :: force_pbc(:,:,:,:)
    real(dp),       allocatable, save :: virial_cell(:,:)
    real(dp)                          :: virial(3,3)
    real(dp)                          :: virial_long(3,3)
    real(dp)                          :: virial_extern(3,3)

    ! REST allocation
    if ( .not. allocated(force) ) then
      nc = size(domain%force(3,MaxAtom,:))
      nv = size(domain%virial_cellpair(3,:))
      allocate(force(3,MaxAtom,nc), &
               force_long(3,MaxAtom,nc), &
               force_omp(3,MaxAtom,nc,nthread), &
               force_pbc(3,MaxAtom,nc,nthread), &
               virial_cell(3,nv))
    end if

    call compute_energy_fep(domain, enefunc, pairlist, boundary, &
                            domain%coord, .false., &
                            .true., .true., .true., .false., &
                            energy, domain%translated, force, force_long, &
                            force_omp, force_pbc, virial_cell, virial, &
                            virial_long, virial_extern)

    energy%total = energy%total + &
               energy%restraint_distance + &
               energy%restraint_rmsd     + &
               energy%restraint_emfit    + &
               energy%contact            + &
               energy%noncontact

    if (enefunc%dispersion_corr /= Disp_Corr_NONE) then
      energy%total = energy%total + energy%disp_corr_energy
    end if

    deallocate(force, force_long, force_omp, force_pbc, virial_cell)

    return

  end subroutine compute_energy_tgt

  !=============================================================================

  subroutine assign_condition_fep(parmsetid, alchemy, domain, ensemble, &
                                  enefunc, remd)

    ! formal arguments
    integer,                 intent(in)    :: parmsetid
    type(s_alchemy),         intent(in)    :: alchemy
    type(s_domain),          intent(inout) :: domain
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_remd),  optional, intent(inout) :: remd

    ! local variables
    integer                  :: i, j, k
    integer                  :: umbrid, funcid, atmid
    integer                  :: lambid
    real(wp)                 :: rest_param

    do i = 1, remd%dimension
      lambid = remd%iparameters(i,remd%parmidsets(parmsetid,i))
      enefunc%lambljA   = remd%dlambljA(lambid)
      enefunc%lambljB   = remd%dlambljB(lambid)
      enefunc%lambelA   = remd%dlambelA(lambid)
      enefunc%lambelB   = remd%dlambelB(lambid)
      enefunc%lambbondA = remd%dlambbondA(lambid)
      enefunc%lambbondB = remd%dlambbondB(lambid)
      enefunc%lambrest  = remd%dlambrest(lambid)
      call set_lambda_table_fep(enefunc)

      if (remd%types(i) == RemdAlchemyRest) then

        if (ensemble%ensemble == EnsembleNVE ) then
          call error_msg('Fep_Rest> Fep_Rest is not available for NVE!')
        end if

        rest_param = remd%dparameters(i,remd%parmidsets(parmsetid,i))
        call assign_condition_feprest( remd%fep_rest(i), &
                                                rest_param, domain, &
                                                ensemble, enefunc )

      end if

    end do

    return

  end subroutine assign_condition_fep

  !=============================================================================

  subroutine assign_condition_feprest( soltemp, rest_param, &
                                                domain, ensemble, enefunc )

    ! formal arguments
    type(s_soltemp),         intent(inout) :: soltemp
    real(wp),                intent(in)    :: rest_param
    type(s_domain),          intent(inout) :: domain
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_enefunc),         intent(inout) :: enefunc

    integer  :: alist(8), i, j, ix, num, tgt, n, org, ncell, counter
    real(wp) :: coeff_full, coeff_half, wgt, tmp1, tmp2
    real(wp) :: el_fact, alpha

    coeff_full = ensemble%temperature / rest_param
    coeff_half = sqrt( coeff_full )

    if ( soltemp%done_setup ) then
      tmp1 = coeff_full
      tmp2 = coeff_half

      coeff_full = coeff_full / soltemp%rest_param_full
      coeff_half = coeff_half / soltemp%rest_param_half

      ! remember coeffs for next exchange
      !
      soltemp%rest_param_full = tmp1
      soltemp%rest_param_half = tmp2
    else
      soltemp%rest_param_full = coeff_full
      soltemp%rest_param_half = coeff_half
      soltemp%done_setup = .true.
    end if

    ncell = domain%num_cell_local + domain%num_cell_boundary

    if ( soltemp%sw_charge ) then
      do i = 1, ncell
        ! charge
        do ix = 1, domain%num_atom(i)
          if ( soltemp%is_solute(domain%id_l2g(ix,i)) > 0 ) then
            domain%charge(ix,i) = coeff_half * domain%charge(ix,i)
          end if
        end do
      end do
#ifdef USE_GPU
      call gpu_upload_charge( domain%charge )
#endif /* USE_GPU */

      ! Calculating self energy
      u_self_preserve = 0.0_dp
      u_self_appear   = 0.0_dp
      u_self_vanish   = 0.0_dp
      el_fact = ELECOEF / enefunc%dielec_const
      alpha = enefunc%pme_alpha
      do i = 1, domain%num_cell_local
        do ix = 1, domain%num_atom_preserve(i)
          u_self_preserve = u_self_preserve &
            + domain%charge(domain%pmelist_preserve(ix,i),i)**2
        end do
        do ix = 1, domain%num_atom_appear_gr(i)
          u_self_appear = u_self_appear &
            + domain%charge(domain%pmelist_appear_gr(ix,i),i)**2
        end do
        do ix = 1, domain%num_atom_vanish_gr(i)
          u_self_vanish = u_self_vanish &
            + domain%charge(domain%pmelist_vanish_gr(ix,i),i)**2
        end do
      end do
      u_self_preserve = - u_self_preserve * el_fact * alpha/sqrt(PI)
      u_self_appear   = - u_self_appear * el_fact * alpha/sqrt(PI)
      u_self_vanish   = - u_self_vanish * el_fact * alpha/sqrt(PI)
    end if

    do i = 1, domain%num_cell_local

      ! bond
      if ( soltemp%sw_bonds ) then
        call assign_condition_feprest_internal( soltemp, &
          2, coeff_full, &
          enefunc%num_bond(i), &
          enefunc%bond_list(:,:,i), &
          enefunc%bond_force_const(:,i) )
      end if

      ! do not add sw_angles flags here!
      call assign_condition_feprest_angle_ub( soltemp, coeff_full, &
        enefunc%num_angle(i), &
        enefunc%angle_list(:,:,i), &
        enefunc%angle_force_const(:,i), &
        enefunc%urey_force_const(:,i) )

      ! dihedral
      if ( soltemp%sw_dihedrals ) then
        call assign_condition_feprest_internal( soltemp, &
          4, coeff_full, &
          enefunc%num_dihedral(i), &
          enefunc%dihe_list(:,:,i), &
          enefunc%dihe_force_const(:,i) )
      end if

      ! impropers
      if ( soltemp%sw_impropers ) then
        call assign_condition_feprest_internal( soltemp, &
          4, coeff_full, &
          enefunc%num_improper(i), &
          enefunc%impr_list(:,:,i), &
          enefunc%impr_force_const(:,i) )
      end if

    end do

    ! lj
    if ( soltemp%sw_lj ) then
      n = enefunc%num_atom_cls
      if ( allocated(enefunc%nb14_lj6) ) then
        call assign_condition_feprest_lj( &
                 n, soltemp%rest_param_half, soltemp%rest_param_full, &
                 soltemp, enefunc%nb14_lj6 )
      end if
      if ( allocated(enefunc%nb14_lj12) ) then
        call assign_condition_feprest_lj( &
                 n, soltemp%rest_param_half, soltemp%rest_param_full, &
                 soltemp, enefunc%nb14_lj12 )
      end if
      if ( allocated(enefunc%nonb_lj6) ) then
        call assign_condition_feprest_lj( &
                 n, soltemp%rest_param_half, soltemp%rest_param_full, &
                 soltemp, enefunc%nonb_lj6 )
      end if
      if ( allocated(enefunc%nonb_lj12) ) then
        call assign_condition_feprest_lj( &
                 n, soltemp%rest_param_half, soltemp%rest_param_full, &
                 soltemp, enefunc%nonb_lj12 )
      end if
#ifdef USE_GPU
      if ( allocated(enefunc%nonb_lj12) .or. allocated(enefunc%nonb_lj6) ) then
        call gpu_upload_lj_coeffs( n, enefunc%nonb_lj12, enefunc%nonb_lj6 )
      end if
#endif /* USE_GPU */
    end if

    ! cmap
    if ( soltemp%num_cmap_type > 0 .and. soltemp%sw_cmaps ) then
      do i = 1, soltemp%num_cmap_type
        tgt = i + soltemp%istart_cmap_type - 1
        org = soltemp%cmap_type_org(i)
        wgt = soltemp%rest_param_full ** soltemp%cmap_weight(i)
        enefunc%cmap_coef(1:4,1:4,1:24,1:24,tgt) = &
          wgt * enefunc%cmap_coef(1:4,1:4,1:24,1:24,org)
      end do
    end if

    return

  end subroutine assign_condition_feprest

  !=============================================================================

  subroutine assign_condition_feprest_angle_ub( soltemp, &
                                  coeff_full, nangles, aindex, fc, fc_ub )

    ! formal arguments
    type(s_soltemp),         intent(in)    :: soltemp
    real(wp),                intent(in)    :: coeff_full
    integer,                 intent(in)    :: nangles
    integer,                 intent(in)    :: aindex(:,:)
    real(wp),                intent(inout) :: fc(:)
    real(wp),                intent(inout) :: fc_ub(:)

    ! local variables
    integer  :: alist(3), ix, j, num
    real(wp) :: wgt

    ! angle and urey
    do ix = 1, nangles
      alist(1:3) = aindex(1:3,ix)
      num = 0
      do j = 1, 3
        if (soltemp%is_solute(alist(j)) > 0) num = num + 1
      end do
      if ( num > 0 .and. soltemp%sw_angles ) then
        wgt = real(num,wp) / 3.0_wp
        wgt = coeff_full ** wgt
        fc(ix) = wgt * fc(ix)
      end if
      if ( fc_ub(ix) > EPS .and. soltemp%sw_ureys ) then
        num = 0
        if (soltemp%is_solute(alist(1)) > 0) num = num + 1
        if (soltemp%is_solute(alist(3)) > 0) num = num + 1
        if ( num > 0 ) then
          wgt = real(num,wp) / 2.0_wp
          wgt = coeff_full ** wgt
          fc_ub(ix) = wgt * fc_ub(ix)
        end if
      end if
    end do

    return

  end subroutine assign_condition_feprest_angle_ub

  !=============================================================================

  subroutine assign_condition_feprest_lj( n, coeff_half, coeff_full, &
                                                   soltemp, nbcoeff )

    ! formal arguments
    integer,                 intent(in)    :: n
    real(wp),                intent(in)    :: coeff_half
    real(wp),                intent(in)    :: coeff_full
    type(s_soltemp),         intent(in)    :: soltemp
    real(wp),                intent(inout) :: nbcoeff(n,n)

    ! local variables
    integer :: i, j, oldcount, newcount, org, orgj

    oldcount = soltemp%istart_atom_cls - 1
    newcount = oldcount + soltemp%num_atom_cls

    do i = 1, soltemp%num_atom_cls
      org = soltemp%atom_cls_no_org(i)
      do j = 1, oldcount
        nbcoeff(j,i+oldcount) = coeff_half * nbcoeff(j,org)
        nbcoeff(i+oldcount,j) = coeff_half * nbcoeff(j,org)
      end do
      do j = oldcount + 1, newcount
        orgj = soltemp%atom_cls_no_org(j-oldcount)
        nbcoeff(j,i+oldcount) = coeff_full * nbcoeff(orgj,org)
        nbcoeff(i+oldcount,j) = coeff_full * nbcoeff(orgj,org)
      end do
    end do

    return

  end subroutine assign_condition_feprest_lj

  !=============================================================================

  subroutine assign_condition_feprest_internal( soltemp, n, &
                                  coeff_full, n_internal, aindex, fc )

    ! formal arguments
    type(s_soltemp),         intent(in)    :: soltemp
    integer,                 intent(in)    :: n
    real(wp),                intent(in)    :: coeff_full
    integer,                 intent(in)    :: n_internal
    integer,                 intent(in)    :: aindex(:,:)
    real(wp),                intent(inout) :: fc(:)

    ! local variables
    integer  :: alist(1:n), ix, j, num
    real(wp) :: wgt

    do ix = 1, n_internal
      alist(1:n) = aindex(1:n,ix)
      num = 0
      do j = 1, n
        if (soltemp%is_solute(alist(j)) > 0) num = num + 1
      end do
      if ( num > 0 ) then
        wgt = real(num,wp) / real(n,wp)
        wgt = coeff_full ** wgt
        fc(ix) = wgt * fc(ix)
      end if
    end do

    return

  end subroutine assign_condition_feprest_internal

end module sp_fep_energy_mod
