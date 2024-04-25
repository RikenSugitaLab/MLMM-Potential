!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_dynvars_mod
!> @brief   define dynamic variables
!! @authors Jaewoon Jung (JJ), Takaharu Mori (TM), Yuji Sugita (YS), 
!!          Chigusa Kobayashi (CK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_dynvars_mod

  use sp_output_str_mod
  use sp_dynamics_str_mod
  use sp_dynvars_str_mod
  use sp_ensemble_str_mod
  use sp_boundary_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use molecules_str_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use sp_alchemy_str_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif
  
  implicit none
  private

  ! variables
  integer, public, save :: DynvarsOut = 6

  logical,         save :: etitle = .true.

  ! subroutines
  public  :: setup_dynvars
  public  :: output_dynvars
  public  :: compute_dynvars
  private :: output_dynvars_genesis
  private :: output_dynvars_charmm
  private :: output_dynvars_namd
  private :: output_dynvars_gromacs
  private :: compute_pressure
  private :: reduce_property

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_dynvars
  !> @brief        copy molecule information into dynvars
  !! @authors      JJ
  !! @param[out]   dynvars  : dynamic variables
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_dynvars(dynvars, dynamics)

    ! formal arguments
    type(s_dynvars),            intent(inout) :: dynvars
    type(s_dynamics), optional, intent(in)    :: dynamics

    ! initialize variables in dynvars
    !
    call init_dynvars(dynvars)

    ! setup dynvars
    !
    if (present(dynamics)) then
      if (dynamics%restart .and. dynamics%integrator == IntegratorLEAP) then
        dynvars%step = 1
      else
        dynvars%step = 0
      endif
    end if

    return

  end subroutine setup_dynvars

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_dynvars
  !> @brief        output dynamical variables
  !! @authors      YS
  !! @param[in]    output  : output information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    dynvars : dynamical variables information
  !! @param[in]    ensemble : ensemble information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_dynvars(output, enefunc, dynvars, ensemble, boundary)

    ! formal arguments
    type(s_output),             intent(in) :: output
    type(s_enefunc),            intent(in) :: enefunc
    type(s_dynvars),            intent(in) :: dynvars
    type(s_ensemble), optional, intent(in) :: ensemble
    type(s_boundary), optional, intent(in) :: boundary

    ! If REMD, output is done on replica_main_rank
    ! If MD, output is done on main_rank
    !
    if (output%replica) then
      if (.not. replica_main_rank) return
    else if (output%rpath) then
      if (.not. replica_main_rank) return
    else
      if (.not. main_rank) return
    end if

    select case (enefunc%output_style)

    case (OutputStyleGENESIS)

      call output_dynvars_genesis(output, enefunc, dynvars, ensemble, boundary)

    case (OutputStyleCHARMM)

      call output_dynvars_charmm (enefunc, dynvars)

    case (OutputStyleNAMD)

      call output_dynvars_namd   (dynvars)

    case (OutputStyleGROMACS)

      call output_dynvars_gromacs(enefunc, dynvars)

    end select

    return

  end subroutine output_dynvars

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_dynvars
  !> @brief        calculate other dynamical variables with domain
  !! @authors      JJ
  !! @param[in]    dynamics : dynamics information
  !! @param[in]    boundary : boundary information
  !! @param[in]    ensemble : ensemble information
  !! @param[inout] domain   : domain information
  !! @param[inout] dynvars  : dynamical variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_dynvars(enefunc, dynamics, boundary, ensemble, domain, &
                             dynvars)

    ! formal arguments
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_boundary),        intent(in)    :: boundary
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_domain),  target, intent(inout) :: domain
    type(s_dynvars), target, intent(inout) :: dynvars

    ! local variables
    real(dp)                 :: kex, key, kez
    real(dp)                 :: ekin, ekin_old, ekin_new
    real(dp)                 :: vel_oldx, vel_oldy, vel_oldz
    real(dp)                 :: rmsg, viri_ke, viri_ke_ext
    real(dp)                 :: gamma_t, dt, scale_v
    real(dp)                 :: virial_sum(1:3,1:3)
    integer                  :: i, ix, k, l

    real(dp),        pointer :: mass(:,:)
    real(dp),        pointer :: vel(:,:,:), vel_ref(:,:,:), vel_full(:,:,:)
    real(dp),        pointer :: force(:,:,:)
    real(dp),        pointer :: virial(:,:), virial_long(:,:), virial_ext(:,:)
    real(dp),        pointer :: viri_const(:,:)
    integer,         pointer :: ncell, natom(:), nsolute(:), nwater(:)
    integer,         pointer :: water_list(:,:,:)


    ncell       => domain%num_cell_local
    natom       => domain%num_atom
    nsolute     => domain%num_solute
    nwater      => domain%num_water 
    water_list  => domain%water_list
    mass        => domain%mass
    vel         => domain%velocity
    vel_ref     => domain%velocity_ref
    vel_full    => domain%velocity_full
    force       => domain%force

    virial      => dynvars%virial
    virial_long => dynvars%virial_long
    virial_ext  => dynvars%virial_extern
    viri_const  => dynvars%virial_const

    gamma_t     =  ensemble%gamma_t*AKMA_PS
    dt          =  dynamics%timestep/AKMA_PS
    if (dynamics%xo_respa) dt =  dt * real(dynamics%thermo_period,dp)

    ! Compute velocities(t + dt)
    !   v(t+dt) = 0.5*(v(t+3/2dt) + v(t+1/2dt))
    !
    if (dynamics%integrator == IntegratorLEAP) then

      do i = 1, ncell
        do ix = 1, natom(i)
          vel_full(1,ix,i) = 0.5_dp*(vel(1,ix,i) + vel_ref(1,ix,i))
          vel_full(2,ix,i) = 0.5_dp*(vel(2,ix,i) + vel_ref(2,ix,i))
          vel_full(3,ix,i) = 0.5_dp*(vel(3,ix,i) + vel_ref(3,ix,i))
        end do
      end do

    else

      do i = 1, ncell
        do ix = 1, natom(i)
          vel_full(1,ix,i) = vel(1,ix,i)
          vel_full(2,ix,i) = vel(2,ix,i)
          vel_full(3,ix,i) = vel(3,ix,i)
        end do
      end do

    end if

    ! kinetic energy(t+dt) and temperature(t+dt)
    ! 
    kex = 0.0_dp
    key = 0.0_dp
    kez = 0.0_dp

    do i = 1, ncell
      do ix = 1, natom(i)
        ! FEP: skip singleB to avoid duplication
        if (domain%fep_use) then
          if (domain%fepgrp(ix,i) == 2) cycle
        end if
        kex = kex + mass(ix,i)*vel_full(1,ix,i)*vel_full(1,ix,i)
        key = key + mass(ix,i)*vel_full(2,ix,i)*vel_full(2,ix,i)
        kez = kez + mass(ix,i)*vel_full(3,ix,i)*vel_full(3,ix,i)
      end do
    end do

    ekin  = 0.5_dp * (kex + key + kez)

    ! other dynvars(t+dt)
    !
    ekin_old = 0.0_dp
    ekin_new = 0.0_dp
    rmsg     = 0.0_dp

    if (dynamics%integrator == IntegratorLEAP) then

      do i = 1, ncell
        do ix = 1, natom(i)

          ! FEP: skip singleB to avoid duplication
          if (domain%fep_use) then
            if (domain%fepgrp(ix,i) == 2) cycle
          end if
          vel_oldx = 2.0_dp*vel_full(1,ix,i)-vel(1,ix,i)
          vel_oldy = 2.0_dp*vel_full(2,ix,i)-vel(2,ix,i)
          vel_oldz = 2.0_dp*vel_full(3,ix,i)-vel(3,ix,i)
          ekin_old = ekin_old + &
                     (vel_oldx**2+vel_oldy**2+vel_oldz**2)*mass(ix,i)
          ekin_new = ekin_new + &
                     (vel(1,ix,i)**2+ vel(2,ix,i)**2+vel(3,ix,i)**2)*mass(ix,i)
          rmsg     = rmsg+ force(1,ix,i)**2+force(2,ix,i)**2+force(3,ix,i)**2

         end do
      end do

      ekin_old = 0.5_dp * ekin_old
      ekin_new = 0.5_dp * ekin_new

    else

      do i = 1, ncell
        do ix = 1, nsolute(i)
          ! FEP: skip singleB to avoid duplication
          if (domain%fep_use) then
            if (domain%fepgrp(ix,i) == 2) cycle
          end if
          rmsg = rmsg + force(1,ix,i)*force(1,ix,i) &
                      + force(2,ix,i)*force(2,ix,i) &
                      + force(3,ix,i)*force(3,ix,i)
        end do
        do l = 1, nwater(i)
          do k = 1, 3
            ix = water_list(k,l,i)
            ! FEP: skip singleB to avoid duplication
            if (domain%fep_use) then
              if (domain%fepgrp(ix,i) == 2) cycle
            end if
            rmsg = rmsg + force(1,ix,i)*force(1,ix,i) &
                        + force(2,ix,i)*force(2,ix,i) &
                        + force(3,ix,i)*force(3,ix,i)
          end do
        end do
      end do

    end if

    if (ensemble%tpcontrol == TpcontrolLangevin) then
      scale_v = 1.0_dp + 0.5_dp*gamma_t*dt
      kex  = kex * scale_v
      key  = key * scale_v
      kez  = kez * scale_v
      ekin = ekin * scale_v
      ekin_old= ekin_old * scale_v
      ekin_new = ekin_new * scale_v
    end if

    ! virial (virial + virial_const)
    !

    if (dynamics%integrator /= IntegratorLEAP)  then
      virial_sum(1:3,1:3) = virial(1:3,1:3) + viri_const(1:3,1:3)
    else
      virial_sum(1:3,1:3) = virial(1:3,1:3)
    end if

    if (dynamics%integrator == IntegratorVRES) &
      virial_sum(1:3,1:3) = virial_sum(1:3,1:3) + virial_long(1:3,1:3)

#ifdef HAVE_MPI_GENESIS
    call reduce_property(dynvars%energy%bond,               &
                         dynvars%energy%angle,              &
                         dynvars%energy%urey_bradley,       &
                         dynvars%energy%dihedral,           &
                         dynvars%energy%improper,           &
                         dynvars%energy%cmap,               &
                         dynvars%energy%electrostatic,      &
                         dynvars%energy%van_der_waals,      &
                         dynvars%energy%restraint_position, &
                         dynvars%energy%contact,            &
                         dynvars%energy%noncontact,         &
                         dynvars%energy%total,              &
                         virial_sum(1,1), virial_sum(2,2), virial_sum(3,3), &
                         virial_ext(1,1), virial_ext(2,2), virial_ext(3,3), &
                         kex, key, kez, ekin, ekin_old, ekin_new, rmsg)

#endif

    dynvars%total_kene  = ekin
    dynvars%temperature = 2.0_dp*ekin/(real(domain%num_deg_freedom,dp)*KBOLTZ)

    ! potential energy(t+dt) and total energy(t+dt)
    !
    dynvars%total_pene = &
               dynvars%energy%bond               + &
               dynvars%energy%angle              + &
               dynvars%energy%urey_bradley       + &
               dynvars%energy%dihedral           + &
               dynvars%energy%improper           + &
               dynvars%energy%cmap               + &
               dynvars%energy%electrostatic      + &
               dynvars%energy%van_der_waals      + &
               dynvars%energy%restraint_distance + &
               dynvars%energy%restraint_position + &
               dynvars%energy%restraint_rmsd     + &
               dynvars%energy%restraint_emfit    + &
               dynvars%energy%contact            + &
               dynvars%energy%noncontact

    if (enefunc%dispersion_corr /= Disp_Corr_NONE) then
      dynvars%total_pene = dynvars%total_pene + &
                           dynvars%energy%disp_corr_energy
    end if
    dynvars%total_energy = &
               dynvars%total_pene + dynvars%total_kene

    ! GaMD
    !
    if (enefunc%gamd_use) then
      if (enefunc%gamd%boost_dih) then
        dynvars%total_energy = dynvars%total_energy + &
                             dynvars%energy%dihedral_gamd
      else if (enefunc%gamd%boost_pot) then
        dynvars%total_energy = dynvars%total_energy + &
                             dynvars%energy%total_gamd
      else if (enefunc%gamd%boost_dual) then
        dynvars%total_energy = dynvars%total_energy + &
                             dynvars%energy%total_gamd + &
                             dynvars%energy%dihedral_gamd
      end if
    end if

    ! volume(t+dt) and pressure(t+dt)
    !
    dynvars%box_size_x = boundary%box_size_x_ref
    dynvars%box_size_y = boundary%box_size_y_ref
    dynvars%box_size_z = boundary%box_size_z_ref

    dynvars%volume =  boundary%box_size_x_ref &
                    * boundary%box_size_y_ref &
                    * boundary%box_size_z_ref

    call compute_pressure(domain%num_deg_freedom, dynvars%volume, &
                          dynvars%temperature, virial_sum,          &
                          dynvars%internal_pressure)

    call compute_pressure(domain%num_deg_freedom, dynvars%volume,   &
                          dynvars%temperature, virial_ext,          &
                          dynvars%external_pressure)

    dynvars%pressure_xx = P_ATMOS*(kex + virial_sum(1,1))/dynvars%volume
    dynvars%pressure_yy = P_ATMOS*(key + virial_sum(2,2))/dynvars%volume
    dynvars%pressure_zz = P_ATMOS*(kez + virial_sum(3,3))/dynvars%volume

    dynvars%hfc_kene     = 0.5_dp*(ekin_old + ekin_new)
    ! FEP: remove degrees of freedom of duplicate atoms
    if (domain%fep_use) then
      dynvars%rms_gradient = sqrt(rmsg/real(3*(domain%num_atom_all - &
        domain%num_atom_single_all),dp))
    else
      dynvars%rms_gradient = sqrt(rmsg/real(3*domain%num_atom_all,dp))
    end if

    viri_ke =  virial_sum(1,1) + virial_sum(2,2) + virial_sum(3,3)

    dynvars%internal_virial = viri_ke/3.0_dp
    dynvars%virial_kene     = -0.5_dp*viri_ke

    viri_ke_ext =virial_ext(1,1) + virial_ext(2,2) + virial_ext(3,3)
    dynvars%external_virial = viri_ke_ext/3.0_dp

    return

  end subroutine compute_dynvars

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_dynvars_genesis
  !> @brief        output dynvars in GENESIS style
  !! @authors      YS, TM
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    dynvars : dynamical variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_dynvars_genesis(output, enefunc, dynvars, ensemble, boundary)

    ! formal arguments
    type(s_output),             intent(in) :: output
    type(s_enefunc),            intent(in) :: enefunc
    type(s_dynvars),            intent(in) :: dynvars
    type(s_ensemble), optional, intent(in) :: ensemble
    type(s_boundary), optional, intent(in) :: boundary

    ! local variables
    integer,parameter        :: clength=16, flength=4
    integer                  :: i, ifm
    character(16)            :: title
    character(16)            :: category(999)
    character                :: frmt*5, frmt_res*10, rfrmt*7
    character                :: rfrmt_cont*9,frmt_cont*7
    real(dp)                 :: values(999)
    real(dp)                 :: ene_restraint


    title = 'INFO:       STEP'
    write(frmt,'(A2,I2,A)') '(A',clength,')'
    write(frmt_cont,'(A2,I2,A3)') '(A',clength,',$)'
    write(frmt_res,'(A2,I2,A6)') '(A',clength-3,',I3.3)'
    write(rfrmt,'(A2,I2,A1,I1,A1)') '(F',clength,'.',flength,')'
    write(rfrmt_cont,'(A2,I2,A1,I1,A3)') '(F',clength,'.',flength,',$)'


    ifm = 1


    if (present(ensemble)) then
      write(category(ifm),frmt) 'TIME'
      values(ifm) = dynvars%time
      ifm = ifm+1

      write(category(ifm),frmt) 'TOTAL_ENE'
      values(ifm) = dynvars%total_energy
      ifm = ifm+1
    endif

    write(category(ifm),frmt) 'POTENTIAL_ENE'
    values(ifm) = dynvars%total_pene
    ifm = ifm+1

    if (present(ensemble)) then
      write(category(ifm),frmt) 'KINETIC_ENE'
      values(ifm) = dynvars%total_kene
      ifm = ifm+1
    endif

    write(category(ifm),frmt) 'RMSG'
    values(ifm) = dynvars%rms_gradient
    ifm = ifm+1

    if (enefunc%num_bond_all > 0) then
      write(category(ifm),frmt) 'BOND'
      values(ifm) = dynvars%energy%bond
      ifm = ifm+1
    endif

    if (enefunc%num_angl_all > 0) then
      write(category(ifm),frmt) 'ANGLE'
      values(ifm) = dynvars%energy%angle
      ifm = ifm+1

      if (enefunc%forcefield == ForcefieldCHARMM) then
        write(category(ifm),frmt) 'UREY-BRADLEY'
        values(ifm) = dynvars%energy%urey_bradley
        ifm = ifm+1
      endif
    endif

    if (enefunc%num_dihe_all > 0 .or. enefunc%num_rb_dihe_all > 0) then
      write(category(ifm),frmt) 'DIHEDRAL'
      values(ifm) = dynvars%energy%dihedral
      ifm = ifm+1
    endif

    if (enefunc%num_impr_all > 0 ) then
      write(category(ifm),frmt) 'IMPROPER'
      values(ifm) = dynvars%energy%improper
      ifm = ifm+1
    endif

    if (enefunc%forcefield == ForcefieldCHARMM) then

      if (enefunc%num_cmap_all > 0 ) then
        write(category(ifm),frmt) 'CMAP'
        values(ifm) = dynvars%energy%cmap
        ifm = ifm+1
      endif
    endif

    if (enefunc%forcefield == ForcefieldAAGO) then

      write(category(ifm),frmt) 'NATIVE_CONTACT'
      values(ifm) = dynvars%energy%contact
      ifm = ifm+1

      write(category(ifm),frmt) 'NON-NATIVE_CONT'
      values(ifm) = dynvars%energy%noncontact
      ifm = ifm+1

    else

      write(category(ifm),frmt) 'VDWAALS'
      values(ifm) = dynvars%energy%van_der_waals
      ifm = ifm+1

      if (enefunc%dispersion_corr /= Disp_Corr_NONE) then
        write(category(ifm),frmt) 'DISP-CORR_ENE'
        values(ifm) = dynvars%energy%disp_corr_energy
        ifm = ifm+1
      endif

      write(category(ifm),frmt) 'ELECT'
      values(ifm) = dynvars%energy%electrostatic
      ifm = ifm+1

    end if

    if (enefunc%num_restraintfuncs > 0) then
      if (enefunc%restraint_rmsd) then
        write(category(ifm),frmt) 'RMSD'
        values(ifm) = dynvars%energy%rmsd
        ifm = ifm+1
        write(category(ifm),frmt) 'RMSD_target'
        values(ifm) = enefunc%rmsd_target
        ifm = ifm+1
      endif

      if (enefunc%restraint_emfit) then
        write(category(ifm),frmt) 'EMCORR'
        values(ifm) = dynvars%energy%emcorr
        ifm = ifm+1
      end if

      ene_restraint =   dynvars%energy%restraint_distance &
                      + dynvars%energy%restraint_position &
                      + dynvars%energy%restraint_rmsd     &
                      + dynvars%energy%restraint_emfit
      write(category(ifm),frmt) 'RESTRAINT_TOTAL'
      values(ifm) = ene_restraint
      ifm = ifm+1

    endif

    if (present(ensemble) .and. present(boundary)) then

      write(category(ifm),frmt) 'TEMPERATURE'
      values(ifm) = dynvars%temperature
      ifm = ifm+1
    
      write(category(ifm),frmt) 'VOLUME'
      if (boundary%type == BoundaryTypeNOBC) then
        values(ifm) = 0.0_dp
      else
        values(ifm) = dynvars%volume
      end if
      ifm = ifm+1

      if (output%verbose                    .or. &
          ensemble%ensemble == EnsembleNPT  .or. &
          ensemble%ensemble == EnsembleNPAT .or. &
          ensemble%ensemble == EnsembleNPgT ) then
     
        write(category(ifm),frmt) 'BOXX'
        values(ifm) = dynvars%box_size_x
        ifm = ifm+1
     
        write(category(ifm),frmt) 'BOXY'
        values(ifm) = dynvars%box_size_y
        ifm = ifm+1
     
        write(category(ifm),frmt) 'BOXZ'
        values(ifm) = dynvars%box_size_z
        ifm = ifm+1
     
        write(category(ifm),frmt) 'VIRIAL'
        values(ifm) = dynvars%internal_virial
        ifm = ifm+1
     
        if (enefunc%dispersion_corr == Disp_Corr_EPRESS) then
          write(category(ifm),frmt) 'DISP-CORR_VIR'
          values(ifm) = dynvars%energy%disp_corr_virial
          ifm = ifm+1
        endif
     
        write(category(ifm),frmt) 'PRESSURE'
        values(ifm) = dynvars%internal_pressure
        ifm = ifm+1
     
        write(category(ifm),frmt) 'PRESSXX'
        values(ifm) = dynvars%pressure_xx
        ifm = ifm+1
     
        write(category(ifm),frmt) 'PRESSYY'
        values(ifm) = dynvars%pressure_yy
        ifm = ifm+1
     
        write(category(ifm),frmt) 'PRESSZZ'
        values(ifm) = dynvars%pressure_zz
        ifm = ifm+1
     
      endif
    endif

    if (enefunc%gamd_use) then
        if (enefunc%gamd%boost_pot) then
            write(category(ifm),frmt) 'POTENTIAL_GAMD'
            values(ifm) = dynvars%energy%total_gamd
            ifm = ifm+1
        else if (enefunc%gamd%boost_dih) then
            write(category(ifm),frmt) 'DIHEDRAL_GAMD'
            values(ifm) = dynvars%energy%dihedral_gamd
            ifm = ifm+1
        else if (enefunc%gamd%boost_dual) then
            write(category(ifm),frmt) 'POTENTIAL_GAMD'
            values(ifm) = dynvars%energy%total_gamd
            ifm = ifm+1
            write(category(ifm),frmt) 'DIHEDRAL_GAMD'
            values(ifm) = dynvars%energy%dihedral_gamd
            ifm = ifm+1
        end if
    end if

    if (etitle) then

      write(DynvarsOut,'(A,$)') title

      do i = 1, ifm-1

        if (i == ifm-1) then
          write(DynvarsOut,frmt) category(i)
        else
          write(DynvarsOut,frmt_cont) category(i)
        endif
      end do

      write(DynvarsOut,'(A80)') ' --------------- ---------------'//&
                 ' --------------- --------------- ---------------'
      etitle = .false.
    end if

    write(DynvarsOut,'(A5,1x,I10,$)') 'INFO:', dynvars%step

    do i = 1, ifm-1
      if (i == ifm-1) then
        write(DynvarsOut,rfrmt) values(i)
      else
        write(DynvarsOut,rfrmt_cont) values(i)
      endif
    end do
    write(DynvarsOut,'(A)') ''

    return

  end subroutine output_dynvars_genesis

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_dynvars_charmm
  !> @brief        output dynvars in CHARMM style
  !! @authors      YS, CK
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    dynvars : dynamical variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_dynvars_charmm(enefunc, dynvars)

    ! formal arguments
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_dynvars), target, intent(in)    :: dynvars

    ! local variables
    real(dp)                 :: time, totener, totke, energy, temperature
    real(dp)                 :: grms, hfctote, hfcke, ehfcor, virke
    real(dp)                 ::                hbonds, asp, user
    real(dp)                 :: vire, viri, presse, pressi, volume
    real(dp)                 :: boxx, boxy, boxz, pixx, piyy, pizz
    real(dp)                 :: cdihe, cintcr, noe
    integer                  :: step

    real(dp),        pointer :: bonds, angles, urey_b, dihedrals, impropers
    real(dp),        pointer :: cmaps, disp_corr
    real(dp),        pointer :: vdwaals, elec
    real(dp),        pointer :: posicon, restdist


    step        = 0
    time        = 0.0_dp
    hbonds      = 0.0_dp
    totener     = 0.0_dp
    totke       = 0.0_dp
    energy      = 0.0_dp
    temperature = 0.0_dp
    asp         = 0.0_dp
    user        = 0.0_dp
    vire        = 0.0_dp
    viri        = 0.0_dp
    presse      = 0.0_dp
    pressi      = 0.0_dp
    volume      = 0.0_dp
    boxx        = 0.0_dp
    boxy        = 0.0_dp
    boxz        = 0.0_dp
    pixx        = 0.0_dp
    piyy        = 0.0_dp
    pizz        = 0.0_dp
    hfctote     = 0.0_dp
    ehfcor      = 0.0_dp
    volume      = 0.0_dp


    step        = dynvars%step
    time        = dynvars%time
    energy      = dynvars%total_pene
    totke       = dynvars%total_kene
    totener     = dynvars%total_energy
    temperature = dynvars%temperature
    grms        = dynvars%rms_gradient
    hfcke       = dynvars%hfc_kene
    virke       = dynvars%virial_kene

    cdihe       = 0.0_dp
    cintcr      = 0.0_dp
    noe         = 0.0_dp

    viri        = dynvars%internal_virial
    vire        = dynvars%external_virial
    pressi      = dynvars%internal_pressure
    presse      = dynvars%external_pressure
    volume      = dynvars%volume
    boxx        = dynvars%box_size_x
    boxy        = dynvars%box_size_y
    boxz        = dynvars%box_size_z
    pixx        = dynvars%pressure_xx 
    piyy        = dynvars%pressure_yy
    pizz        = dynvars%pressure_zz

    bonds      => dynvars%energy%bond
    angles     => dynvars%energy%angle
    urey_b     => dynvars%energy%urey_bradley
    dihedrals  => dynvars%energy%dihedral
    impropers  => dynvars%energy%improper
    elec       => dynvars%energy%electrostatic
    vdwaals    => dynvars%energy%van_der_waals
    cmaps      => dynvars%energy%cmap
    restdist   => dynvars%energy%restraint_distance
    posicon    => dynvars%energy%restraint_position
    disp_corr  => dynvars%energy%disp_corr_energy

    ! write title if necessary
    !
    if (etitle) then
      write(DynvarsOut,'(A)') 'Output_Energy> CHARMM_Style is used'
      write(DynvarsOut,'(A)') ' '
      write(DynvarsOut,'(A79)') 'DYNA DYN: Step         Time      TOTEner        TOTKe       ENERgy  TEMPerature'
      write(DynvarsOut,'(A79)') 'DYNA PROP:             GRMS      HFCTote        HFCKe       EHFCor        VIRKe'
      write(DynvarsOut,'(A79)') 'DYNA INTERN:          BONDs       ANGLes       UREY-b    DIHEdrals    IMPRopers'
      if (enefunc%num_cmap_all > 0) then
        write(DynvarsOut,'(A79)') 'DYNA CROSS:           CMAPs                                                    '
      end if
      if (enefunc%dispersion_corr /= Disp_Corr_NONE) then
        write(DynvarsOut,'(A79)') 'DYNA  DISP:       Disp-Corr                                                    '
      endif
      write(DynvarsOut,'(A79)') 'DYNA EXTERN:        VDWaals         ELEC       HBONds          ASP         USER'
      if (enefunc%restraint) then
        write(DynvarsOut,'(A79)') 'DYNA CONSTR:       HARMonic    CDIHedral        CIC     RESDistance       NOE'
      end if
      write(DynvarsOut,'(A79)') 'DYNA PRESS:            VIRE         VIRI       PRESSE       PRESSI       VOLUme'
      write(DynvarsOut,'(A79)') ' ----------       ---------    ---------    ---------    ---------    ---------'
      etitle = .false.
    end if


    ! write energy in CHARMM-style
    !
    write(DynvarsOut,'(A5,I9,5F13.5)') 'DYNA>', step, time, totener, totke, energy, temperature
    write(DynvarsOut,'(A14,5F13.5)')   'DYNA PROP>    ', grms, hfctote, hfcke, ehfcor, virke
    write(DynvarsOut,'(A14,5F13.5)')   'DYNA INTERN>  ', bonds, angles, urey_b, dihedrals, impropers

    if (enefunc%num_cmap_all > 0) then
      write(DynvarsOut,'(A14,F13.5)')    'DYNA CROSS>   ', cmaps
    end if
    if (enefunc%dispersion_corr /= Disp_Corr_NONE) then
      write(DynvarsOut,'(A14,F13.5)')    'DYNA DISP>    ', disp_corr
    end if

    write(DynvarsOut,'(A14,5F13.5)')   'DYNA EXTERN>  ', vdwaals, elec, hbonds, asp, user

    if (enefunc%restraint) then
      write(DynvarsOut,'(A14,5F13.5)')   'DYNA CONSTR>  ', posicon, cdihe, cintcr, restdist, noe
    end if

    write(DynvarsOut,'(A14,5F13.5)')   'DYNA PRESS>   ', vire, viri, presse, pressi, volume
    write(DynvarsOut,'(A79)') ' ----------       ---------    ---------    ---------    ---------    ---------'

    return

  end subroutine output_dynvars_charmm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_dynvars_namd
  !> @brief        output dynvars in NAMD style
  !! @authors      YS, CK
  !! @param[in]    dynvars : dynamical variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_dynvars_namd(dynvars)

    ! formal arguments
    type(s_dynvars), target, intent(in)    :: dynvars

    ! local variables
    real(dp)                 :: ebond, eangle
    real(dp)                 ::               misc, kinetic
    real(dp)                 :: total, temp, total2, total3, tempavg
    real(dp)                 :: pressure, gpressure, volume
    real(dp)                 :: pressavg, gpressavg
    integer                  :: nstep

    real(dp),        pointer :: edihed, eimprp
    real(dp),        pointer :: eelect, evdw
    real(dp),        pointer :: eboundary


    nstep     =  0
    ebond     =  dynvars%energy%bond  + dynvars%energy%restraint_distance
    eangle    =  dynvars%energy%angle + dynvars%energy%urey_bradley

    edihed    => dynvars%energy%dihedral
    eimprp    => dynvars%energy%improper
    eelect    => dynvars%energy%electrostatic
    evdw      => dynvars%energy%van_der_waals
    eboundary => dynvars%energy%restraint_position


    ! write title if necessary
    !
    if (etitle) then
      write(DynvarsOut,'(A)') 'Output_Energy> NAMD_Style is used'
      write(DynvarsOut,'(A)') ' '
      write(DynvarsOut,'(A75)') 'ETITLE:      TS           BOND          ANGLE          DIHED          IMPRP'
      write(DynvarsOut,'(A75)') '          ELECT            VDW       BOUNDARY           MISC        KINETIC'
      write(DynvarsOut,'(A75)') '          TOTAL           TEMP         TOTAL2         TOTAL3        TEMPAVG'
      write(DynvarsOut,'(A75)') '       PRESSURE      GPRESSURE         VOLUME       PRESSAVG      GPRESSAVG'
      write(DynvarsOut,'(A)') ' '
      etitle = .false.
    end if


    ! write energy in NAMD-style
    !
    write(DynvarsOut,'(A7,I8,4F15.4)') &
                                 'ENERGY:', nstep, ebond, eangle, edihed, eimprp
    write(DynvarsOut,'(5F15.4)') eelect, evdw, eboundary, misc, kinetic
    write(DynvarsOut,'(5F15.4)') total, temp, total2, total3, tempavg
    write(DynvarsOut,'(5F15.4)') pressure, gpressure, volume, pressavg,gpressavg
    write(DynvarsOut,'(A)') ' '


    return

  end subroutine output_dynvars_namd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_dynvars_gromacs
  !> @brief        output dynvars in GROMACS style
  !! @authors      YS, CK
  !! @param[in]    dynvars : dynamical variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_dynvars_gromacs(enefunc, dynvars)

    ! formal arguments
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_dynvars), target, intent(in)    :: dynvars

    ! local variables
    real(dp)                 :: time, bonds, angles, urey_b
    real(dp)                 :: dihedrals, impropers, vdwaals, elec
    real(dp)                 :: restdist, posicon
    real(dp)                 :: energy, totke, totener
    real(dp)                 :: temperature, pressi, presse, disp_corr
    real(dp)                 :: pres_dc
    integer                  :: i, step


    step        = dynvars%step
    time        = dynvars%time
    bonds       = dynvars%energy%bond               * CAL2JOU
    angles      = dynvars%energy%angle              * CAL2JOU
    urey_b      = dynvars%energy%urey_bradley       * CAL2JOU
    dihedrals   = dynvars%energy%dihedral           * CAL2JOU
    impropers   = dynvars%energy%improper           * CAL2JOU
    vdwaals     = dynvars%energy%van_der_waals      * CAL2JOU
    elec        = dynvars%energy%electrostatic      * CAL2JOU
    restdist    = dynvars%energy%restraint_distance * CAL2JOU
    posicon     = dynvars%energy%restraint_position * CAL2JOU
    energy      = dynvars%total_pene                * CAL2JOU
    totke       = dynvars%total_kene                * CAL2JOU
    totener     = dynvars%total_energy              * CAL2JOU
    temperature = dynvars%temperature
    pressi      = dynvars%internal_pressure         * ATM2BAR
    presse      = dynvars%external_pressure         * ATM2BAR
    disp_corr   = dynvars%energy%disp_corr_energy   * CAL2JOU
    pres_dc     = dynvars%energy%disp_corr_virial   * P_ATMOS/dynvars%volume &
                  * ATM2BAR


    write(DynvarsOut,'(3A15)') &
         'Step', 'Time', 'Lambda'
    write(DynvarsOut,'(I15,2F15.5)') &
         step, time, 0.0_dp
    write(DynvarsOut,'(A)') &
         ' '
    write(DynvarsOut,'(A)') &
         '   Energies (kJ/mol)'

    write(DynvarsOut,'(5A15)') &
         'Bond', 'Angle', 'Urey-bradley', 'Dihedral', 'Improper Dih.'
    write(DynvarsOut,'(5ES15.5E2)') &
         bonds,angles,urey_b,dihedrals,impropers

    if (enefunc%dispersion_corr == Disp_Corr_NONE) then
      write(DynvarsOut,'(5A15)') &
           'LJ (1-4,SR', ' Coulomb(1-4,SR', 'Position Rest.', 'Potential', 'Kinetic En.'
      write(DynvarsOut,'(5ES15.5E2)') &
           vdwaals,elec,posicon,energy,totke
    else
      write(DynvarsOut,'(5A15)') &
         'LJ (1-4,SR', ' Coulomb(1-4,SR', 'Disper. corr.', 'Position Rest.', &
         'Pres. DC (bar)'
      write(DynvarsOut,'(5ES15.5E2)') &
           vdwaals,elec,disp_corr,posicon,pres_dc
      write(DynvarsOut,'(3A15)') &
         'Position Rest.', 'Potential', 'Kinetic En.'
      write(DynvarsOut,'(3ES15.5E2)') &
          posicon,energy,totke

    endif

    write(DynvarsOut,'(5A15)') &
         'Total Energy', 'Temperature', 'Pressure(int.)', 'Pressure(ext.)'
    write(DynvarsOut,'(5ES15.5E2)') &
         totener,temperature,pressi,presse

    write(DynvarsOut,'(A)')  ' '

    return

  end subroutine output_dynvars_gromacs

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_pressure
  !> @brief        compute interal pressure
  !! @authors      TM
  !! @param[in]    num_degree  : number of degrees of freedom
  !! @param[in]    volume      : volume (A^3)
  !! @param[in]    temperature : temperature (K)
  !! @param[in]    virial      : internal virial
  !! @param[out]   pressure    : pressure (atm)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_pressure(num_degree, volume, temperature, virial, pressure)

    ! formal arguments
    integer,                 intent(in)    :: num_degree
    real(dp),                intent(in)    :: volume
    real(dp),                intent(in)    :: temperature
    real(dp),                intent(in)    :: virial(3,3)
    real(dp),                intent(out)   :: pressure

    ! local variables
    real(dp)                 :: viri_ke


    ! compute internal pressure
    !
    viri_ke  = virial(1,1) + virial(2,2) + virial(3,3)
    pressure = P_ATMOS / (3.0_dp * volume)                      &
              * (num_degree * KBOLTZ * temperature + viri_ke)


    return

  end subroutine compute_pressure

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    reduce_property
  !> @brief        reduce property
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine reduce_property(val1, val2, val3, val4, val5, val6, val7, val8,   &
                             val9, val10, val11, val12, val13, val14, val15,   &
                             val16, val17, val18, val19, val20, val21, val22,  &
                             val23, val24, val25)

    ! formal arguments
    real(dp),                intent(inout) ::  val1,  val2,  val3
    real(dp),                intent(inout) ::  val4,  val5,  val6
    real(dp),                intent(inout) ::  val7,  val8,  val9
    real(dp),                intent(inout) :: val10, val11, val12
    real(dp),                intent(inout) :: val13, val14, val15
    real(dp),                intent(inout) :: val16, val17, val18
    real(dp),                intent(inout) :: val19, val20, val21
    real(dp),                intent(inout) :: val22, val23, val24, val25

#ifdef HAVE_MPI_GENESIS
    ! local variables
    real(dp)                 :: before_reduce(25), after_reduce(25)


    before_reduce(1)  = val1
    before_reduce(2)  = val2
    before_reduce(3)  = val3
    before_reduce(4)  = val4
    before_reduce(5)  = val5
    before_reduce(6)  = val6
    before_reduce(7)  = val7
    before_reduce(8)  = val8
    before_reduce(9)  = val9
    before_reduce(10) = val10
    before_reduce(11) = val11
    before_reduce(12) = val12
    before_reduce(13) = val13
    before_reduce(14) = val14
    before_reduce(15) = val15
    before_reduce(16) = val16
    before_reduce(17) = val17
    before_reduce(18) = val18
    before_reduce(19) = val19
    before_reduce(20) = val20
    before_reduce(21) = val21
    before_reduce(22) = val22
    before_reduce(23) = val23
    before_reduce(24) = val24
    before_reduce(25) = val25

    call mpi_allreduce(before_reduce, after_reduce, 25, mpi_real8, &
                       mpi_sum, mpi_comm_country, ierror)

    val1  = after_reduce(1)
    val2  = after_reduce(2)
    val3  = after_reduce(3)
    val4  = after_reduce(4)
    val5  = after_reduce(5)
    val6  = after_reduce(6)
    val7  = after_reduce(7)
    val8  = after_reduce(8)
    val9  = after_reduce(9)
    val10 = after_reduce(10)
    val11 = after_reduce(11)
    val12 = after_reduce(12)
    val13 = after_reduce(13)
    val14 = after_reduce(14)
    val15 = after_reduce(15)
    val16 = after_reduce(16)
    val17 = after_reduce(17)
    val18 = after_reduce(18)
    val19 = after_reduce(19)
    val20 = after_reduce(20)
    val21 = after_reduce(21)
    val22 = after_reduce(22)
    val23 = after_reduce(23)
    val24 = after_reduce(24)
    val25 = after_reduce(25)

#endif
    return

  end subroutine reduce_property

end module sp_dynvars_mod
