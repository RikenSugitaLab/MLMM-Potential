!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_md_respa_mod
!> @brief   perform molecular dynamics simulation with velocity verlet
!!          and respa algorithm
!! @authors Jaewoon Jung (JJ)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_md_respa_mod

  use sp_output_mod
  use sp_update_domain_mod
  use sp_assign_velocity_mod
  use sp_dynvars_mod
  use sp_constraints_mod
  use sp_communicate_mod
  use sp_energy_str_mod
  use sp_energy_mod
  use sp_output_str_mod
  use sp_dynamics_str_mod
  use sp_dynvars_str_mod
  use sp_ensemble_str_mod
  use sp_constraints_str_mod
  use sp_constraints_mod
  use sp_boundary_str_mod
  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use sp_remd_str_mod    ! fep/rest
  use sp_alchemy_str_mod ! fep/rest
  use random_mod
  use math_libs_mod
  use messages_mod
  use timers_mod
  use mpi_parallel_mod
  use constants_mod
  use sp_fep_utils_mod
  use sp_fep_energy_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: vverlet_respa_dynamics
  private :: initial_vverlet
  private :: integrate_vv1
  private :: integrate_vv2
  private :: nve_vv1
  private :: nve_vv2
  private :: nosehoover_thermostat
  private :: mtk_barostat_vv1
  private :: mtk_thermostat
  private :: langevin_thermostat_vv1
  private :: langevin_thermostat_vv2
  private :: langevin_barostat_vv1
  private :: langevin_barostat_vv2
  private :: update_barostat
  ! FEP
  public  :: vverlet_respa_dynamics_fep
  private :: initial_vverlet_fep

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    vverlet_respa_dynamics
  !> @brief        velocity verlet integrator using respa
  !! @authors      JJ
  !! @param[inout] output      : output information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] dynvars     : dynamic variables information
  !! @param[inout] dynamics    : dynamics information
  !! @param[inout] pairlist    : non-bond pair list information
  !! @param[inout] boundary    : boundary conditions information
  !! @param[inout] constraints : bond constraint information
  !! @param[inout] ensemble    : ensemble information
  !! @param[inout] comm        : information of communication
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine vverlet_respa_dynamics(output, domain, enefunc, dynvars, &
                                    dynamics, pairlist, boundary,     &
                                    constraints, ensemble, comm, remd)

    ! formal arguments
    type(s_output),          intent(inout) :: output
    type(s_domain),  target, intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_dynvars), target, intent(inout) :: dynvars
    type(s_dynamics),target, intent(inout) :: dynamics
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_boundary),        intent(inout) :: boundary
    type(s_constraints),     intent(inout) :: constraints
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_comm),            intent(inout) :: comm
    type(s_remd),            intent(inout) :: remd

    ! local variables
    real(dp)                 :: simtim, temperature, factor
    real(dp)                 :: dt_short, dt_long
    real(dp)                 :: min_time
    integer                  :: i, ii, j, k, jx, nsteps
    integer                  :: istep, multistep
    integer                  :: iseed, istart, iend
    integer                  :: min_k, k_min, k_max
    logical                  :: npt

    real(dp),        pointer :: coord(:,:,:), coord_ref(:,:,:)
    real(wp),        pointer :: coord_pbc(:,:,:)
    real(dp),        pointer :: coord_old(:,:,:)
    real(dp),        pointer :: vel(:,:,:), vel_ref(:,:,:), mass(:,:)
    real(wp),        pointer :: force_pbc(:,:,:,:), force_omp(:,:,:,:)
    real(dp),        pointer :: force(:,:,:)
    real(dp),        pointer :: virial_cell(:,:), virial(:,:)
    real(dp),        pointer :: force_short(:,:,:), force_long(:,:,:)
    integer,         pointer :: ncell, natom(:)
    logical,         pointer :: XI_RESPA, XO_RESPA


    ncell       => domain%num_cell_local
    natom       => domain%num_atom
    mass        => domain%mass
    coord       => domain%coord
    coord_ref   => domain%coord_ref
    coord_old   => domain%coord_old
    coord_pbc   => domain%translated
    force       => domain%force    
    force_omp   => domain%force_omp
    force_short => domain%force_short
    force_long  => domain%force_long
    vel         => domain%velocity
    vel_ref     => domain%velocity_ref
    force_pbc   => domain%force_pbc
    virial_cell => domain%virial_cellpair
    virial      => dynvars%virial
    XI_RESPA    => dynamics%xi_respa
    XO_RESPA    => dynamics%xo_respa

    temperature =  ensemble%temperature
    nsteps      =  dynamics%nsteps
    istart      =  dynamics%istart_step
    iend        =  dynamics%iend_step
    multistep   =  dynamics%elec_long_period
    dt_short    =  dynamics%timestep/AKMA_PS
    dt_long     =  real(multistep,dp)*dt_short
    simtim      =  dynamics%initial_time
    iseed       =  dynamics%iseed_init_velocity
    npt         = ensemble%use_barostat

    if (abs(dynamics%initial_rmsd) < 0.001_dp)  &
      dynamics%initial_rmsd = dynvars%energy%rmsd
    if (dynamics%target_md) enefunc%rmsd_force = 1 / (dt_long*dt_long)

    ! Check restart
    !
    if (.not. dynamics%restart) then

      call initial_velocity(temperature,           &
                            domain%num_atom_all,   &
                            domain%num_cell_local, &
                            domain%id_g2l,         &
                            domain%mass,           &
                            iseed,                 &
                            domain%velocity)

      call stop_trans_rotation(domain%num_cell_local,         &
                               domain%num_atom,               &
                               dynamics%stop_com_translation, &
                               dynamics%stop_com_rotation,    &
                               domain%mass,                   &
                               domain%coord,                  &
                               domain%velocity)

      call initial_vverlet(npt, output, enefunc, dynamics,       &
                           pairlist, boundary, ensemble, constraints, &
                           domain, dynvars, comm)

    else

      ! After 2nd cycle of REMD simulation (istart /= 1), this is skipped
      !
      if (istart == 1) then

        if (constraints%tip4 .or. enefunc%table%tip4) &
        call decide_dummy(domain, constraints, coord)

        call communicate_coor(domain, comm)

        call compute_energy(domain, enefunc, pairlist, boundary, coord,     &
                            npt, .false., .true., &
                            .false., .false., dynvars%energy, coord_pbc,     &
                            force_short, force_long, force_omp,     &
                            force_pbc, virial_cell, dynvars%virial, &
                            dynvars%virial_long, dynvars%virial_extern)
  
        call communicate_force(domain, comm, force_short)
        call communicate_force(domain, comm, force_long)

        if (constraints%tip4 .or. enefunc%table%tip4) then
          call water_force_redistribution(constraints, domain, force_long,  &
                                          virial)
          call water_force_redistribution(constraints, domain, force_short, &
                                          virial)
        end if

      end if
    
    end if

    ! Main loop
    !
!    do i = 1, nsteps / multistep
    do i = int((istart-1)/multistep)+1, int(iend / multistep)

      simtim = simtim + dt_long * AKMA_PS
      dynvars%time = simtim
!      dynvars%step = i * multistep + istart - 1
      dynvars%step = i * multistep
      enefunc%rpath_sum_mf_flag = .false.

      ! Inner integrator
      !
      do j = 1, multistep

        ! decide the current steps
        !
        istep = (i-1) * multistep + j

        if (dynamics%target_md .or. dynamics%steered_md) &
          enefunc%rmsd_target = dynamics%initial_rmsd &
                              + (dynamics%final_rmsd-dynamics%initial_rmsd) &
                               *real(istep,dp)/real(nsteps,dp)

        ! save coordiantes for constraint
        !
        call timer(TimerIntegrator, TimerOn)
        call timer(TimerUpdate, TimerOn)
        do k = 1, ncell
          do jx = 1, natom(k)
            coord_ref(1:3,jx,k) = coord(1:3,jx,k)
            vel_ref  (1:3,jx,k) = vel  (1:3,jx,k)
          end do
        end do
        call timer(TimerUpdate, TimerOff)

        ! VV1
        !
        call integrate_vv1(dynamics, istep, j, dt_long, dt_short, ensemble, &
                           domain, constraints, boundary, dynvars)

        call timer(TimerIntegrator, TimerOff)

        ! cell migration and update cell pairlist
        !

        if (dynamics%nbupdate_period > 0 .and. &
            j == multistep .and. i > 0) then

          call domain_interaction_update_md(istep, dynamics, domain, enefunc, &
                                          pairlist, boundary, constraints, comm)

        end if

        ! short range forces
        !
        call timer(TimerIntegrator, TimerOn)

        if (constraints%tip4 .or. enefunc%table%tip4) &
        call decide_dummy(domain, constraints, coord)

        call communicate_coor(domain, comm)

        call timer(TimerIntegrator, TimerOff)

        if (j < multistep) then

          pairlist%univ_update = 0
          enefunc%rpath_sum_mf_flag = enefunc%rpath_flag
          call compute_energy_short(domain, enefunc, pairlist, boundary, coord, &
                                    npt, mod(istep,dynamics%eneout_period) == 0,&
                                    dynvars%energy,                             &
                                    coord_pbc,                                  &
                                    force_short,                                &
                                    force_omp,                                  &
                                    force_pbc,                                  &
                                    virial_cell,                                &
                                    dynvars%virial,                             &
                                    dynvars%virial_extern)

          call timer(TimerIntegrator, TimerOn)
          call timer(TimerComm2, TimerOn)
          call communicate_force(domain, comm, force_short)
          if (constraints%tip4 .or. enefunc%table%tip4) &
            call water_force_redistribution(constraints, domain, force_short, &
                                            virial)

          call timer(TimerComm2, TimerOff)
          call timer(TimerIntegrator, TimerOff)

        end if
        
        ! full step velocities with foce_short
        ! long range force for last inner step
        !
        if (j == multistep) then

          call compute_energy(domain, enefunc, pairlist, boundary, coord,        &
                              npt, .false., mod(istep,dynamics%eneout_period)==0,&
                              .false.,                                           &
                              .false.,                                           &
                              dynvars%energy,                                    &
                              coord_pbc,                                         &
                              force_short,                                       &
                              force_long,                                        &
                              force_omp,                                         &
                              force_pbc,                                         &
                              virial_cell,                                       &
                              dynvars%virial,                                    &
                              dynvars%virial_long,                               &
                              dynvars%virial_extern)

          call timer(TimerIntegrator, TimerOn)
          call timer(TimerComm2, TimerOn)
          call communicate_force(domain, comm, force_long)
          call communicate_force(domain, comm, force_short)
          if (constraints%tip4 .or. enefunc%table%tip4) then
            call water_force_redistribution(constraints, domain, force_long,  &
                                            virial)
            call water_force_redistribution(constraints, domain, force_short, &
                                            virial)
          end if
          call timer(TimerComm2, TimerOff)
          call timer(TimerIntegrator, TimerOff)

        end if

        ! VV2
        !
        call integrate_vv2(dynamics, istep, j, dt_long, dt_short, ensemble,    &
                           domain, constraints, boundary, dynvars)

        ! Remove translational and rotational motion about COM(t + dt)
        !
        if (dynamics%stoptr_period > 0) then
          if (mod(istep,dynamics%stoptr_period) == 0) then

              call stop_trans_rotation(domain%num_cell_local,         &
                                       domain%num_atom,               &
                                       dynamics%stop_com_translation, &
                                       dynamics%stop_com_rotation,    &
                                       domain%mass,                   &
                                       domain%coord,                  &
                                       domain%velocity)
          end if
        end if

      end do

      ! energy output
      !
      if (dynamics%eneout_period > 0) then
        if (mod(istep,dynamics%eneout_period) == 0 ) then
          do k = 1, ncell
            do jx = 1, natom(k)
              force(1:3,jx,k) = force_short(1:3,jx,k) + force_long(1:3,jx,k)
            end do
          end do
        end if
      end if

      call output_md(output, dynamics, boundary, pairlist, &
                     ensemble, dynvars, domain, enefunc, remd)

      ! output parallel I/O restart
      !
      call output_prst_md(output, enefunc, dynamics, boundary,                 &
                                  dynvars, domain, constraints)

    end do

    ! Close output files
    !
    !call close_output(output)

    return

  end subroutine vverlet_respa_dynamics

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    initial_vverlet
  !> @brief        compute the first step (0+dt)
  !! @authors      JJ
  !! @param[in]    npt         : flag for NPT or not
  !! @param[in]    output      : output information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    pairlist    : pairlist information
  !! @param[in]    boundary    : boundary information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] constraints : constraints information
  !! @param[inout] domain      : domain information
  !! @param[inout] dynvars     : dynamic variables information
  !! @param[inout] comm        : communication information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

   subroutine initial_vverlet(npt, output, enefunc, dynamics, pairlist, &
                              boundary, ensemble, constraints, domain,  &
                              dynvars, comm)

    ! formal arguments
    logical,                 intent(in)    :: npt
    type(s_output),          intent(in)    :: output
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_constraints),     intent(inout) :: constraints
    type(s_domain),  target, intent(inout) :: domain
    type(s_dynvars), target, intent(inout) :: dynvars
    type(s_comm),            intent(inout) :: comm

    ! local variables
    real(dp)                 :: factor, temperature, energy(18), temp(18)
    real(dp)                 :: imass, simtim, dt, friction
    integer                  :: i, ix, j, jx, ncell, k, l

    real(dp),        pointer :: coord(:,:,:), coord_ref(:,:,:)
    real(wp),        pointer :: coord_pbc(:,:,:)
    real(dp),        pointer :: vel(:,:,:), vel_ref(:,:,:), mass(:,:)
    real(dp),        pointer :: force(:,:,:)
    real(wp),        pointer :: force_omp(:,:,:,:)
    real(wp),        pointer :: force_pbc(:,:,:,:)
    real(dp),        pointer :: virial_cell(:,:), virial(:,:)
    real(dp),        pointer :: force_short(:,:,:), force_long(:,:,:)
    integer,         pointer :: natom(:), nsolute(:), nwater(:)
    integer,         pointer :: water_list(:,:,:)


    natom       => domain%num_atom
    nwater      => domain%num_water
    nsolute     => domain%num_solute
    water_list  => domain%water_list
    coord       => domain%coord
    coord_ref   => domain%coord_ref
    coord_pbc   => domain%translated
    force       => domain%force
    force_short => domain%force_short
    force_long  => domain%force_long 
    force_omp   => domain%force_omp
    vel         => domain%velocity
    vel_ref     => domain%velocity_ref
    mass        => domain%mass
    force_pbc   => domain%force_pbc
    virial_cell => domain%virial_cellpair
    virial      => dynvars%virial

    dt          =  dynamics%timestep/AKMA_PS
    simtim      =  dynamics%initial_time
    temperature =  ensemble%temperature
    friction    =  ensemble%gamma_t * AKMA_PS
    ncell       =  domain%num_cell_local

    dynvars%time = simtim
    dynvars%step = 0

    ! save coordinates(0) and velocities(0)
    ! if rigid-body on, update coordinates(0)
    !
    do i = 1, ncell
      do ix = 1, natom(i)
        coord_ref(1:3,ix,i) = coord(1:3,ix,i)
        vel_ref(1:3,ix,i)   = vel(1:3,ix,i)
      end do
    end do

    if (constraints%rigid_bond) then

      call compute_constraints(ConstraintModeLEAP, .false., dt, coord_ref, &
                               domain, constraints, coord, vel,            &
                               dynvars%virial_const)

      do i = 1, ncell
        do ix = 1, natom(i)
          coord_ref(1:3,ix,i) = coord(1:3,ix,i)
        end do
      end do
    end if

    ! calculate energy(0) and forces(0)
    !
    if (constraints%tip4 .or. enefunc%table%tip4) &
    call decide_dummy(domain, constraints, coord)

    call communicate_coor(domain, comm)

    call compute_energy(domain, enefunc, pairlist, boundary, coord,     &
                        npt, .false., .true.,                   &
                        .false., .false., dynvars%energy, coord_pbc,     &
                        force_short, force_long, force_omp,     &
                        force_pbc, virial_cell, dynvars%virial, &
                        dynvars%virial_long, dynvars%virial_extern)

    call communicate_force(domain, comm, force_short)
    call communicate_force(domain, comm, force_long)

    if (constraints%tip4 .or. enefunc%table%tip4) then
      call water_force_redistribution(constraints, domain, force_long,  &
                                      virial)
      call water_force_redistribution(constraints, domain, force_short, &
                                      virial)
    end if

    do i = 1, ncell
      do ix = 1, natom(i)
        force(1:3,ix,i) = force_short(1:3,ix,i) + force_long(1:3,ix,i)
      end do
    end do

    ! if rigid-body on, update velocity(0 + dt/2 and 0 - dt/2)
    !
    if (constraints%rigid_bond) then

      ! calculate velocities(0 + dt/2) and coordinates(0 + dt)
      ! update coordinates(0 + dt) and velocities(0 + dt/2)
      !
      do i = 1, ncell
        do ix = 1, nsolute(i)
          imass = 1.0_dp/mass(ix,i)
          vel(1:3,ix,i) = vel_ref(1:3,ix,i) + 0.5_dp*dt*force(1:3,ix,i)*imass
          coord(1:3,ix,i) = coord_ref(1:3,ix,i) + dt*vel(1:3,ix,i)
        end do
        do l = 1, nwater(i)
          do k = 1, 3
            ix = water_list(k,l,i)
            imass = 1.0_dp/mass(ix,i)
            vel(1:3,ix,i) = vel_ref(1:3,ix,i) + 0.5_dp*dt*force(1:3,ix,i)*imass
            coord(1:3,ix,i) = coord_ref(1:3,ix,i) + dt*vel(1:3,ix,i)
          end do
        end do
      end do

      call compute_constraints(ConstraintModeVVER1, .false., dt, coord_ref,    &
                               domain, constraints, coord, vel,                &
                               dynvars%virial_const)

      ! calculate velocities(0 - dt/2) and coordinates(0 - dt)
      ! update coordinates(0 - dt) and velocities(0 - dt/2)
      !
      do i = 1, ncell
        do ix = 1, nsolute(i)
          imass = 1.0_dp/mass(ix,i)
          vel_ref(1:3,ix,i) = vel_ref(1:3,ix,i)                                &
                            - 0.5_dp*dt*force(1:3,ix,i)*imass
          coord(1:3,ix,i) = coord_ref(1:3,ix,i) - dt*vel_ref(1:3,ix,i)
        end do
        do l = 1, nwater(i)
          do k = 1, 3
            ix = water_list(k,l,i)
            imass = 1.0_dp/mass(ix,i)
            vel_ref(1:3,ix,i) = vel_ref(1:3,ix,i)                                &
                              - 0.5_dp*dt*force(1:3,ix,i)*imass
            coord(1:3,ix,i) = coord_ref(1:3,ix,i) - dt*vel_ref(1:3,ix,i)
          end do
        end do
      end do

      call compute_constraints(ConstraintModeVVER1, .false., -dt, coord_ref,   &
                               domain, constraints, coord, vel_ref,            &
                               dynvars%virial_const)

      ! calculate velocity(0)
      !
      do i = 1, ncell
        do ix = 1, natom(i)
          vel_ref(1:3,ix,i) = 0.5_dp*(vel(1:3,ix,i) + vel_ref(1:3,ix,i))
        end do
      end do

      ! vel <= updated velocities (0) and coord <= constrained
      ! coordinates(0)
      !
      do i = 1, ncell
        do ix = 1, natom(i)
          vel(1:3,ix,i) = vel_ref(1:3,ix,i)
          coord(1:3,ix,i) = coord_ref(1:3,ix,i)
        end do
      end do

    end if

    ! output dynvars(0)
    !
    dynvars%time = 0.0_dp
    dynvars%step = 0

    call compute_dynvars(enefunc, dynamics, boundary, ensemble, domain, dynvars)

    call output_dynvars(output, enefunc, dynvars, ensemble, boundary)

    ! at this point
    !   coord, velocity, and force are at 0

    return

  end subroutine initial_vverlet

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    integrate_vv1
  !> @brief        VV1 with thermostat/barostat
  !! @authors      JJ
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    istep       : dynamics step        
  !! @param[in]    dt_long     : outer loop time step 
  !! @param[in]    dt_short    : inner loop time step 
  !! @param[inout] ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] boundary    : boundary information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine integrate_vv1(dynamics, istep, inner_step, dt_long, dt_short, &
                           ensemble, domain, constraints, boundary, dynvars)

    ! formal arguments
    type(s_dynamics),        intent(in)    :: dynamics
    integer,                 intent(in)    :: istep
    integer,                 intent(in)    :: inner_step
    real(dp),                intent(in)    :: dt_long
    real(dp),                intent(in)    :: dt_short
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_domain),          intent(inout) :: domain
    type(s_constraints),     intent(inout) :: constraints
    type(s_boundary),        intent(inout) :: boundary
    type(s_dynvars),         intent(inout) :: dynvars

    integer  :: alloc_stat, ncell, nh_chain_length

    if (ensemble%ensemble == EnsembleNVE) then

      call nve_vv1(dynamics, inner_step, dt_long, dt_short, domain, &
                   constraints, dynvars)

    else if (ensemble%ensemble == EnsembleNPT .and.  &
             ensemble%tpcontrol == TpcontrolMTK) then

      call mtk_barostat_vv1(dynamics, istep, dt_long, dt_short, ensemble, &
                            domain, constraints, boundary, dynvars)

    else if (ensemble%tpcontrol == TpcontrolLangevin) then

      if (istep == 1 .and. .not. allocated(ensemble%random_force)) then
        ncell = domain%num_cell_local
        allocate(ensemble%random_force(3,MaxAtom,ncell), stat=alloc_stat)
        if (alloc_stat /= 0) call error_msg_alloc
      end if

      if (ensemble%ensemble == EnsembleNVT) then

        call langevin_thermostat_vv1(dynamics, istep, inner_step, &
                                     dt_long, dt_short, ensemble, &
                                     domain, constraints, dynvars)

      else if (ensemble%ensemble == EnsembleNPT) then

        if (dynamics%xi_respa) then
          write(MsgOut,'(A)') 'XI_Respa is not available for NPT'
          call error_msg
        else
          call langevin_barostat_vv1(dynamics, istep, inner_step, dt_long,    &
                                     dt_short, ensemble, domain, constraints, &
                                     boundary, dynvars)
        end if

      end if

    else if (ensemble%ensemble == EnsembleNVT .and. &
             ensemble%tpcontrol == TpcontrolNoseHoover) then
     
      if (dynamics%xi_respa) &

        call nosehoover_thermostat(dynamics, dt_short, ensemble, domain,     &
             dynvars)

      call nve_vv1(dynamics, inner_step, dt_long, dt_short, domain,          &
                   constraints, dynvars)

    else if (ensemble%tpcontrol == TpcontrolBussi) then

      if (ensemble%ensemble == EnsembleNVT) then

        call nve_vv1(dynamics, inner_step, dt_long, dt_short, domain, &
                     constraints, dynvars)

      else if (ensemble%ensemble == EnsembleNPT) then

        call bussi_barostat_vv1(dynamics, istep, inner_step, dt_long,    &
                                dt_short, ensemble, domain, constraints, &
                                boundary, dynvars)

      end if

    end if

    return

  end subroutine integrate_vv1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    integrate_vv2
  !> @brief        VV2 with thermostat/barostat
  !! @authors      JJ
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    istep       : dynamics step        
  !! @param[in]    dt_long     : outer loop time step 
  !! @param[in]    dt_short    : inner loop time step 
  !! @param[inout] ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] boundary    : boundary information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine integrate_vv2(dynamics, istep, inner_step, dt_long, dt_short, ensemble, &
                           domain, constraints, boundary, dynvars)

    ! formal arguments
    type(s_dynamics),        intent(in)    :: dynamics
    integer,                 intent(in)    :: istep
    integer,                 intent(in)    :: inner_step
    real(dp),                intent(in)    :: dt_long 
    real(dp),                intent(in)    :: dt_short
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_domain),          intent(inout) :: domain
    type(s_constraints),     intent(inout) :: constraints
    type(s_boundary),        intent(inout) :: boundary
    type(s_dynvars),         intent(inout) :: dynvars

    integer  :: alloc_stat, ncell

    if (ensemble%ensemble == EnsembleNVE) then

      call timer(TimerIntegrator, TimerOn)
      call nve_vv2(dynamics, inner_step, dt_long, dt_short, domain, &
                   constraints, dynvars)
      call timer(TimerIntegrator, TimerOff)

    else if (ensemble%ensemble == EnsembleNPT .and.  &
             ensemble%tpcontrol == TpcontrolMTK) then

      call mtk_barostat_vv2(dynamics, dt_long, dt_short, ensemble, domain, &
                            constraints, boundary, dynvars)

    else if (ensemble%tpcontrol == TpcontrolLangevin) then

      if (ensemble%ensemble == EnsembleNVT) then

        call langevin_thermostat_vv2(dynamics, istep, inner_step, dt_long, &
                                     dt_short, ensemble, domain,           &
                                     constraints, dynvars)

      else if (ensemble%ensemble == EnsembleNPT) then

        if (dynamics%xi_respa) then
          write(MsgOut,'(A)') 'XI_Respa is not available for NPT'
          call error_msg
        else
          call mpi_barrier(mpi_comm_country, ierror)
          call langevin_barostat_vv2  (dynamics, istep, inner_step, dt_long, &
                                       dt_short, ensemble, domain,           &
                                       constraints, boundary, dynvars)
          call mpi_barrier(mpi_comm_country, ierror)
        end if

      end if

    else if (ensemble%ensemble == EnsembleNVT .and.  &
             ensemble%tpcontrol == TpcontrolNoseHoover) then

      call nve_vv2(dynamics, inner_step, dt_long, dt_short, domain,            &
                   constraints, dynvars)

      if (dynamics%xi_respa) &
        call nosehoover_thermostat(dynamics, dt_short, ensemble, domain,       &
                                   dynvars)

    else if (ensemble%tpcontrol == TpcontrolBussi) then

      if (ensemble%ensemble == EnsembleNVT) then

        call nve_vv2(dynamics, inner_step, dt_long, dt_short, domain,         &
                     constraints, dynvars)

        if (mod(istep,dynamics%thermo_period) == 0) &
          call bussi_thermostat(dynamics, dt_short, ensemble, domain, dynvars)

      else if (ensemble%ensemble == EnsembleNPT) then

        call bussi_barostat_vv2(dynamics, istep, inner_step, dt_long,         &
                                dt_short, ensemble, domain,                   &
                                constraints, boundary, dynvars)

      end if

    end if

    return

  end subroutine integrate_vv2

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    nve_vv1
  !> @brief        VV1 with NVE
  !! @authors      JJ
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    dt_long     : long time step        
  !! @param[in]    dt_short    : short time step        
  !! @param[inout] domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine nve_vv1(dynamics, inner_step, dt_long, dt_short, domain, &
                     constraints, dynvars)

    ! formal arguments
    type(s_dynamics),        intent(in)    :: dynamics
    integer,                 intent(in)    :: inner_step
    real(dp),                intent(in)    :: dt_long
    real(dp),                intent(in)    :: dt_short
    type(s_domain),  target, intent(inout) :: domain
    type(s_constraints),     intent(inout) :: constraints
    type(s_dynvars), target, intent(inout) :: dynvars

    ! local variables
    integer                  :: i, ix, k, l, ncell
    real(dp)                 :: h_dt_long, h_dt_short
    real(dp)                 :: factor

    integer,         pointer :: natom(:), nsolute(:), nwater(:)
    integer,         pointer :: water_list(:,:,:)
    real(dp),        pointer :: coord(:,:,:), coord_ref(:,:,:)
    real(dp),        pointer :: vel(:,:,:)
    real(dp),        pointer :: force_long(:,:,:), force_short(:,:,:)
    real(dp),        pointer :: mass(:,:)
    real(dp),        pointer :: viri_const(:,:)

    ! use pointers
    !
    ncell       =  domain%num_cell_local
    natom       => domain%num_atom
    nsolute     => domain%num_solute
    nwater      => domain%num_water
    water_list  => domain%water_list
    mass        => domain%mass
    coord       => domain%coord
    coord_ref   => domain%coord_ref
    vel         => domain%velocity
    force_long  => domain%force_long
    force_short => domain%force_short
    viri_const  => dynvars%virial_const

    h_dt_long   = dt_long  / 2.0_dp
    h_dt_short  = dt_short / 2.0_dp

    call timer(TimerUpdate, TimerOn)

    ! VV1
    if (inner_step == 1) then

      do i = 1, ncell
        do ix = 1, nsolute(i)
          factor = h_dt_long / mass(ix,i)
          vel(1:3,ix,i)   = vel(1:3,ix,i)   + factor*force_long(1:3,ix,i)
        end do
        do l = 1, nwater(i)
          do k = 1, 3
            ix = water_list(k,l,i)
            factor = h_dt_long / mass(ix,i)
            vel(1:3,ix,i)   = vel(1:3,ix,i)   + factor*force_long(1:3,ix,i)
          end do
        end do
      end do

    end if

    do i = 1, ncell
      do ix = 1, nsolute(i)
        factor = h_dt_short / mass(ix,i)
        vel(1:3,ix,i)   = vel(1:3,ix,i)   + factor*force_short(1:3,ix,i)
        coord(1:3,ix,i) = coord(1:3,ix,i) + dt_short*vel(1:3,ix,i)
      end do
      do l = 1, nwater(i)
        do k = 1, 3
          ix = water_list(k,l,i)
          factor = h_dt_short / mass(ix,i)
          vel(1:3,ix,i)   = vel(1:3,ix,i)   + factor*force_short(1:3,ix,i)
          coord(1:3,ix,i) = coord(1:3,ix,i) + dt_short*vel(1:3,ix,i)
        end do
      end do
    end do
    call timer(TimerUpdate, TimerOff)

    ! RATTLE VV1
    !
    if (constraints%rigid_bond) then
      call compute_constraints(ConstraintModeLEAP, .true., dt_short,  &
                               coord_ref, domain, constraints, coord, &
                               vel, viri_const)
      dynvars%virial_const(1:3,1:3) = 2.0_dp * dynvars%virial_const(1:3,1:3)

    end if

    return

  end subroutine nve_vv1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    nve_vv2
  !> @brief        VV2 with NVE
  !! @authors      JJ
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    dt_long     : long time step        
  !! @param[in]    dt_short    : short time step        
  !! @param[inout] domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine nve_vv2(dynamics, inner_step, dt_long, dt_short, domain, &
                     constraints, dynvars)

    ! formal arguments
    type(s_dynamics),        intent(in)    :: dynamics
    integer,                 intent(in)    :: inner_step
    real(dp),                intent(in)    :: dt_long
    real(dp),                intent(in)    :: dt_short
    type(s_domain),  target, intent(inout) :: domain
    type(s_constraints),     intent(inout) :: constraints
    type(s_dynvars), target, intent(inout) :: dynvars

    ! local variables
    integer                  :: i, ix, k, l, ncell
    real(dp)                 :: dt, h_dt_long, h_dt_short
    real(dp)                 :: factor

    integer,         pointer :: natom(:), nsolute(:), nwater(:)
    integer,         pointer :: water_list(:,:,:)
    real(dp),        pointer :: coord(:,:,:), coord_ref(:,:,:)
    real(dp),        pointer :: vel(:,:,:)
    real(dp),        pointer :: force_long(:,:,:), force_short(:,:,:)
    real(dp),        pointer :: mass(:,:)
    real(dp),        pointer :: viri_const(:,:)

    ! use pointers
    !
    ncell       =  domain%num_cell_local
    natom       => domain%num_atom
    nsolute     => domain%num_solute
    nwater      => domain%num_water 
    water_list  => domain%water_list
    mass        => domain%mass
    coord       => domain%coord
    coord_ref   => domain%coord_ref
    vel         => domain%velocity
    force_long  => domain%force_long
    force_short => domain%force_short
    viri_const  => dynvars%virial_const

    h_dt_long   =  dt_long  / 2.0_dp
    h_dt_short  =  dt_short / 2.0_dp

    call timer(TimerUpdate, TimerOn)

    ! VV2
    !
    if (inner_step == dynamics%elec_long_period) then

      do i = 1, ncell
        do ix = 1, nsolute(i)
          factor = h_dt_long / mass(ix,i)
          vel(1:3,ix,i)   = vel(1:3,ix,i)   + factor*force_long(1:3,ix,i)
        end do
        do l = 1, nwater(i)
          do k = 1, 3
            ix = water_list(k,l,i)
            factor = h_dt_long / mass(ix,i)
            vel(1:3,ix,i)   = vel(1:3,ix,i)   + factor*force_long(1:3,ix,i)
          end do
        end do 
      end do

    end if

    do i = 1, ncell
      do ix = 1, nsolute(i)
        factor = h_dt_short / mass(ix,i)
        vel(1:3,ix,i)   = vel(1:3,ix,i)   + factor*force_short(1:3,ix,i)
      end do
      do l = 1, nwater(i)
        do k = 1, 3
          ix = water_list(k,l,i)
          factor = h_dt_short / mass(ix,i)
          vel(1:3,ix,i)   = vel(1:3,ix,i)   + factor*force_short(1:3,ix,i)
        end do
      end do
    end do

    call timer(TimerUpdate, TimerOff)

    ! RATTLE VV2
    !
    if (constraints%rigid_bond) then
      call compute_constraints(ConstraintModeVVER2, .false., dt_short, &
                               coord_ref, domain, constraints, coord,  &
                               vel, viri_const)
    end if

    return

  end subroutine nve_vv2

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    nosehoover_thermostat
  !> @brief        control temperature using Nose-Hoover thermostat
  !! @authors      JJ
  !! @param[in]    dynamics : dynamics information
  !! @param[in]    ensemble : ensemble information
  !! @param[inout] domain   : domain information
  !! @param[inout] dynvars  : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine nosehoover_thermostat(dynamics, dt, ensemble, domain, dynvars)

    ! formal arguments
    type(s_dynamics),        intent(in)    :: dynamics
    real(dp),                intent(in)    :: dt
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_domain),  target, intent(inout) :: domain
    type(s_dynvars), target, intent(inout) :: dynvars

    ! local variables
    real(dp)                 :: kin(1:3), w(1:3)
    real(dp)                 :: temp0, tau_t, dt_small
    real(dp)                 :: dt_1, dt_2, dt_4, dt_8
    real(dp)                 :: ekf, tempf, factor1
    real(dp)                 :: KbT, degree
    real(dp)                 :: scale_kin, scale_vel
    integer                  :: i, j, k, ix
    integer                  :: num_degree, ncell
    integer                  :: nh_length, nh_step

    real(dp),        pointer :: nh_mass(:), nh_vel(:)
    real(dp),        pointer :: nh_force(:), nh_coef(:)
    real(dp),        pointer :: vel(:,:,:), mass(:,:)
    integer,         pointer :: natom(:)


    temp0       =  ensemble%temperature
    tau_t       =  ensemble%tau_t/AKMA_PS
    nh_length   =  ensemble%nhchain
    nh_step     =  ensemble%nhmultistep
    num_degree  =  domain%num_deg_freedom
    ncell       =  domain%num_cell_local
    KbT         =  KBOLTZ * temp0
    degree      =  real(num_degree, dp)

    natom       => domain%num_atom
    vel         => domain%velocity
    mass        => domain%mass
    nh_mass     => dynvars%nh_mass
    nh_vel      => dynvars%nh_velocity
    nh_force    => dynvars%nh_force
    nh_coef     => dynvars%nh_coef

    ! Nose-Hoover mass
    !
    nh_mass(2:nh_length) = KbT * (tau_t**2)
    nh_mass(1)           = degree * nh_mass(2)

    ! Yoshida coefficient
    !
    w(1) = 1.0_dp / (2.0_dp - 2.0_dp**(1.0_dp/3.0_dp))
    w(3) = w(1)
    w(2) = 1.0_dp - w(1) - w(3)

    ! temperature scale factor
    !
    scale_vel = 1.0_dp

    ! calculate kinetic energy
    !
    kin(1:3) = 0.0_dp
    do i = 1, ncell
      do ix = 1, natom(i)
        kin(1:3) = kin(1:3) + mass(ix,i)*vel(1:3,ix,i)*vel(1:3,ix,i)
      end do
    end do
    ekf = kin(1) + kin(2) + kin(3)

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, ekf, 1, mpi_real8, mpi_sum, &
                       mpi_comm_country, ierror)
#endif

    ! decide Nose-Hoover theremostat chain momentum
    !
    dt_small = dt / real(nh_step)

    do i = 1, nh_step
      do j = 1, 3

        dt_1 = w(j) * dt_small
        dt_2 = dt_1 * 0.5_dp
        dt_4 = dt_2 * 0.5_dp
        dt_8 = dt_4 * 0.5_dp

        nh_force(nh_length) = nh_mass(nh_length-1)*nh_vel(nh_length-1)**2-KbT
        nh_force(nh_length) = nh_force(nh_length) / nh_mass(nh_length)
        nh_vel(nh_length) = nh_vel(nh_length) + nh_force(nh_length)*dt_4

        do k = nh_length-1, 2, -1
          nh_force(k) = (nh_mass(k-1)*nh_vel(k-1)**2-KbT) / nh_mass(k)
          nh_coef(k)  = exp(-nh_vel(k+1)*dt_8)
          nh_vel(k)   = nh_vel(k) * nh_coef(k)
          nh_vel(k)   = nh_vel(k) + nh_force(k)*dt_4
          nh_vel(k)   = nh_vel(k) * nh_coef(k)
        end do

        nh_force(1) = (ekf - degree*KbT) / nh_mass(1)
        nh_coef(1)  = exp(-nh_vel(2)*dt_8)
        nh_vel(1)   = nh_vel(1) * nh_coef(1)
        nh_vel(1)   = nh_vel(1) + nh_force(1)*dt_4
        nh_vel(1)   = nh_vel(1) * nh_coef(1)

        scale_kin = exp(-nh_vel(1)*dt_1)
        scale_vel = scale_vel * exp(-nh_vel(1)*dt_2)
        ekf       = ekf * scale_kin

        nh_force(1) = (ekf - degree*KbT) / nh_mass(1)
        nh_vel(1)   = nh_vel(1) * nh_coef(1)
        nh_vel(1)   = nh_vel(1) + nh_force(1)*dt_4
        nh_vel(1)   = nh_vel(1) * nh_coef(1)

        do k = 2, nh_length-1
          nh_force(k) = (nh_mass(k-1)*nh_vel(k-1)**2-KbT) / nh_mass(k)
          nh_vel(k)   = nh_vel(k) * nh_coef(k)
          nh_vel(k)   = nh_vel(k) + nh_force(k)*dt_4
          nh_vel(k)   = nh_vel(k) * nh_coef(k)
        end do

        nh_force(nh_length) = nh_mass(nh_length-1)*nh_vel(nh_length-1)**2-KbT
        nh_force(nh_length) = nh_force(nh_length) / nh_mass(nh_length)
        nh_vel(nh_length) = nh_vel(nh_length) + nh_force(nh_length)*dt_4

      end do
    end do

    ! velocity scaling
    ! 
    do i = 1, ncell
      do ix = 1, natom(i)
        vel(1:3,ix,i) = scale_vel*vel(1:3,ix,i)
      end do
    end do

    return

  end subroutine nosehoover_thermostat


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    mtk_barostat_vv1
  !> @brief        control temperature and pressure using MTK barostat
  !! @authors      JJ
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    istep       : current md step
  !! @param[inout] ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] boundary    : boundary information
  !! @param[inout] dynvars     : dynamic variables information

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine mtk_barostat_vv1(dynamics, istep, dt_long, dt_short, ensemble, &
                              domain, constraints, boundary, dynvars)

    ! formal arguments
    type(s_dynamics),            intent(in)    :: dynamics
    integer,                     intent(in)    :: istep
    real(dp),                    intent(in)    :: dt_long
    real(dp),                    intent(in)    :: dt_short
    type(s_ensemble),    target, intent(inout) :: ensemble
    type(s_domain),      target, intent(inout) :: domain
    type(s_constraints), target, intent(inout) :: constraints
    type(s_boundary),    target, intent(inout) :: boundary
    type(s_dynvars),     target, intent(inout) :: dynvars

    ! local variables
    real(dp)                 :: dt, half_dt, quart_dt, inv_dt
    real(dp)                 :: temp0, press0, degree
    real(dp)                 :: tau_t, tau_p
    real(dp)                 :: KbT
    real(dp)                 :: kin(1:3), ekin
    real(dp)                 :: vel_change(1:3)
    real(dp)                 :: press(1:3), pressxy, pressxyz, pressz
    real(dp)                 :: volume
    real(dp)                 :: bmoment_ref(1:3), scale_b(1:3)
    real(dp)                 :: virial_constraint(1:3,1:3), virial_sum(1:3)
    real(dp)                 :: factor
    real(dp)                 :: size_scale(1:3), baro_force(1:3)
    integer                  :: num_degree, nh_length, nh_step
    integer                  :: ncell, nboundary
    integer                  :: maxiter
    integer                  :: i, j, jx, ij, k, l

    real(dp),        pointer :: mass(:,:)
    real(dp),        pointer :: coord(:,:,:), coord_ref(:,:,:)
    real(dp),        pointer :: vel(:,:,:), vel_ref(:,:,:), force(:,:,:)
    real(dp),        pointer :: temporary(:,:,:)
    real(dp),        pointer :: virial(:,:), viri_const(:,:)
    real(dp),        pointer :: pmass
    real(dp),        pointer :: bmoment(:)
    real(dp),        pointer :: nh_mass(:)
    integer,         pointer :: natom(:), nsolute(:), nwater(:)
    integer,         pointer :: water_list(:,:,:)

    ! use pointers
    !
    dt         =  dynamics%timestep/AKMA_PS
    inv_dt     =  1.0_dp/dt
    temp0      =  ensemble%temperature
    press0     =  ensemble%pressure * ATMOS_P
    tau_t      =  ensemble%tau_t / AKMA_PS
    tau_p      =  ensemble%tau_p / AKMA_PS
    nh_length  =  ensemble%nhchain
    nh_step    =  ensemble%nhmultistep
    num_degree =  domain%num_deg_freedom
    degree     =  real(num_degree, dp)
    ncell      =  domain%num_cell_local
    nboundary  =  domain%num_cell_boundary
    KbT        =  KBOLTZ * temp0

    natom      => domain%num_atom
    nsolute    => domain%num_solute
    nwater     => domain%num_water 
    water_list => domain%water_list
    mass       => domain%mass
    coord      => domain%coord
    coord_ref  => domain%coord_ref
    vel        => domain%velocity
    vel_ref    => domain%velocity_ref
    force      => domain%force
    temporary  => domain%coord_old
    pmass      => ensemble%pmass
    nh_mass    => dynvars%nh_mass
    virial     => dynvars%virial
    viri_const => dynvars%virial_const
    bmoment    => dynvars%barostat_momentum


    ! pressure mass
    !
    if (istep == 1) then
      nh_mass(2:nh_length) = KbT * (tau_t**2)
      nh_mass(1)           = degree * nh_mass(2)
      pmass                = (degree + 3.0_dp) * KbT * (tau_p**2)
    end if

    ! save box size(t+dt) and compute volume(t+dt)
    !
    boundary%box_size_x_ref = boundary%box_size_x
    boundary%box_size_y_ref = boundary%box_size_y
    boundary%box_size_z_ref = boundary%box_size_z
    volume = boundary%box_size_x * boundary%box_size_y * boundary%box_size_z

    ! initial constraint virial
    !
    virial_constraint(1:3,1:3) = 0.0d0

    ! time step
    !
    half_dt   = dt * 0.5_dp
    quart_dt  = dt *0.25_dp

    ! maximum iteration
    !
    if (constraints%rigid_bond) then
      maxiter = 4
    else
      maxiter = 1
    end if

    ! NHC thermostat
    !
    call mtk_thermostat(dynamics, dt_short, ensemble, domain, dynvars)

    ! reference barostat/thermostat coefficients are saved
    !
    bmoment_ref(1:3) = bmoment(1:3)
    do j = 1, ncell
      do jx = 1, natom(j)
        vel_ref(1:3,jx,j) = vel(1:3,jx,j)
      end do
    end do

    ! iteration of barostats
    !
    do i = 1, maxiter

      ! kinetic energy component
      !
      kin(1:3) = 0.0_dp
      do j = 1, ncell
        do jx = 1, natom(j)
          kin(1:3) = kin(1:3) + mass(jx,j)*vel_ref(1:3,jx,j)*vel_ref(1:3,jx,j)
        end do
      end do
      ekin = kin(1) + kin(2) + kin(3)

      ! virial + virial+constraint
      !
      virial_sum(1) = virial(1,1) + virial_constraint(1,1)
      virial_sum(2) = virial(2,2) + virial_constraint(2,2)
      virial_sum(3) = virial(3,3) + virial_constraint(3,3)

      call reduce_pres(kin, ekin, virial_sum)

      ! compute pressure
      !
      press(1:3) = (kin(1:3)+virial_sum(1:3)) / volume
      pressxyz   = (press(1)+press(2)+press(3)) / 3.0_dp
      pressxy    = (press(1)+press(2)) / 2.0_dp
      pressz     = press(3)

      ! update barostat and velocity scale
      !
      if (ensemble%isotropy == IsotropyISO) then
        baro_force(1:3) = 3.0_dp*(ekin/degree+volume*(pressxyz-press0))/pmass
      else if (ensemble%isotropy == IsotropyANISO) then
        baro_force(1:3) = 3.0_dp*(ekin/degree+volume*(press(1:3)-press0))/pmass
      else if (ensemble%isotropy == IsotropySEMI_ISO) then
        baro_force(1:2) = 3.0_dp*(ekin/degree+volume*(pressxy-press0))/pmass
        baro_force(3)   = 3.0_dp*(ekin/degree+volume*(press(3)-press0))/pmass
      else if (ensemble%isotropy == IsotropyXY_Fixed) then
        baro_force(1:2) = 0.0_dp
        baro_force(3) = 3.0_dp*(ekin/degree+volume*(press(3)-press0))/pmass
      end if
      bmoment(1:3) = bmoment_ref(1:3) + baro_force(1:3)*half_dt

      ! VV1 (velocity)
      !
      do j = 1, ncell
        do jx = 1, nsolute(j)
          factor = half_dt/mass(jx,j)
          vel(1:3,jx,j) = vel_ref(1:3,jx,j) + factor*force(1:3,jx,j)
        end do
        do l = 1, nwater(j)
          do k = 1, 3
            jx = water_list(k,l,j)
            factor = half_dt/mass(jx,j)
            vel(1:3,jx,j) = vel_ref(1:3,jx,j) + factor*force(1:3,jx,j)
          end do
        end do
      end do

      ! VV1 (coordinate)
      !
      size_scale(1:3) = exp(bmoment(1:3)*dt)
      do j = 1, ncell
        do jx = 1, natom(j)
          coord(1:3,jx,j) = size_scale(1:3)*coord_ref(1:3,jx,j) &
                          + vel(1:3,jx,j)*dt
        end do
      end do

      ! RATTLE VV1
      !
      if (constraints%rigid_bond) then
        do j = 1, ncell
          do jx = 1, natom(j)
            temporary(1:3,jx,j) = coord(1:3,jx,j)
          end do
        end do
        call compute_constraints(ConstraintModeLEAP, .false., dt_short, &
                                 coord_ref, domain, constraints, coord, &
                                 vel, viri_const)
        virial_constraint(1:3,1:3) = virial_constraint(1:3,1:3)                &
                                   + 2.0_dp * viri_const(1:3,1:3)

        do j = 1, ncell
          do jx = 1, natom(j)
            vel_change(1:3) = (coord(1:3,jx,j) - temporary(1:3,jx,j))*inv_dt
            force(1:3,jx,j) = force(1:3,jx,j) + mass(jx,j)*vel_change(1:3)/half_dt
          end do
        end do
      end if

    end do

    viri_const(1,1) = virial_constraint(1,1)
    viri_const(2,2) = virial_constraint(2,2)
    viri_const(3,3) = virial_constraint(3,3)

    ! compute box size(t+2dt)
    !   size(t+2dt) = exp[eta(t+3/2dt)*dt] * size(t+dt)
    !
    scale_b(1:3) = exp(bmoment(1:3)*dt)
    boundary%box_size_x = scale_b(1) * boundary%box_size_x_ref
    boundary%box_size_y = scale_b(2) * boundary%box_size_y_ref
    boundary%box_size_z = scale_b(3) * boundary%box_size_z_ref
    call bcast_boxsize(boundary%box_size_x, boundary%box_size_y,               &
                       boundary%box_size_z)
    boundary%cell_size_x = boundary%box_size_x / real(boundary%num_cells_x,dp)
    boundary%cell_size_y = boundary%box_size_y / real(boundary%num_cells_y,dp)
    boundary%cell_size_z = boundary%box_size_z / real(boundary%num_cells_z,dp)

    ! update boudary conditions
    !
    dynvars%barostat_momentum(1:3) = bmoment(1:3)
    do j = 1, ncell+nboundary
      do jx = 1, natom(j)
        domain%trans_vec(1:3,jx,j) = domain%trans_vec(1:3,jx,j) * scale_b(1:3)
      end do
    end do

    boundary%origin_x = boundary%origin_x * scale_b(1)
    boundary%origin_y = boundary%origin_y * scale_b(2)
    boundary%origin_z = boundary%origin_z * scale_b(3)

    domain%system_size(1) = boundary%box_size_x
    domain%system_size(2) = boundary%box_size_y
    domain%system_size(3) = boundary%box_size_z

    return

  end subroutine mtk_barostat_vv1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    mtk_barostat_vv2
  !> @brief        control temperature and pressure using MTK barostat
  !! @authors      JJ
  !! @param[in]    dynamics    : dynamics information
  !! @param[inout] ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] boundary    : boundary information
  !! @param[inout] dynvars     : dynamic variables information

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine mtk_barostat_vv2(dynamics, dt_long, dt_short, ensemble, domain, &
                              constraints, boundary, dynvars)

    ! formal arguments
    type(s_dynamics),            intent(in)    :: dynamics
    real(dp),                    intent(in)    :: dt_long
    real(dp),                    intent(in)    :: dt_short
    type(s_ensemble),    target, intent(inout) :: ensemble
    type(s_domain),      target, intent(inout) :: domain
    type(s_constraints), target, intent(inout) :: constraints
    type(s_boundary),    target, intent(inout) :: boundary
    type(s_dynvars),     target, intent(inout) :: dynvars

    ! local variables
    real(dp)                 :: dt, half_dt, quart_dt, inv_dt
    real(dp)                 :: temp0, press0, degree
    real(dp)                 :: tau_t, tau_p
    real(dp)                 :: KbT
    real(dp)                 :: kin(1:3), ekin
    real(dp)                 :: vel_change(1:3)
    real(dp)                 :: press(1:3), pressxy, pressxyz, pressz
    real(dp)                 :: volume
    real(dp)                 :: bmoment_ref(1:3), scale_b(1:3)
    real(dp)                 :: virial_constraint(1:3,1:3), virial_sum(1:3)
    real(dp)                 :: factor
    real(dp)                 :: size_scale(1:3), baro_force(1:3)
    integer                  :: num_degree, nh_length, nh_step
    integer                  :: ncell, nboundary
    integer                  :: maxiter
    integer                  :: i, j, jx, ij, k, l

    real(dp),        pointer :: mass(:,:)
    real(dp),        pointer :: coord(:,:,:), coord_ref(:,:,:)
    real(dp),        pointer :: vel(:,:,:), vel_ref(:,:,:)
    real(dp),        pointer :: force_short(:,:,:), force_long(:,:,:)
    real(dp),        pointer :: force_const(:,:,:)
    real(dp),        pointer :: virial(:,:), viri_const(:,:)
    real(dp),        pointer :: pmass
    real(dp),        pointer :: bmoment(:)
    real(dp),        pointer :: nh_mass(:)
    integer,         pointer :: natom(:), nsolute(:), nwater(:)
    integer,         pointer :: water_list(:,:,:)

    ! use pointers
    !
    dt         =  dynamics%timestep/AKMA_PS
    inv_dt     =  1.0_dp/dt
    temp0      =  ensemble%temperature
    press0     =  ensemble%pressure * ATMOS_P
    tau_t      =  ensemble%tau_t / AKMA_PS
    tau_p      =  ensemble%tau_p / AKMA_PS
    nh_length  =  ensemble%nhchain
    nh_step    =  ensemble%nhmultistep
    num_degree =  domain%num_deg_freedom
    degree     =  real(num_degree, dp)
    ncell      =  domain%num_cell_local
    nboundary  =  domain%num_cell_boundary
    KbT        =  KBOLTZ * temp0

    natom      => domain%num_atom
    nsolute    => domain%num_solute
    nwater     => domain%num_water 
    water_list => domain%water_list
    mass       => domain%mass
    coord      => domain%coord
    coord_ref  => domain%coord_ref
    vel        => domain%velocity
    vel_ref    => domain%velocity_ref
    force_long => domain%force_long
    force_short=> domain%force_short
    force_const=> domain%coord_old
    pmass      => ensemble%pmass
    nh_mass    => dynvars%nh_mass
    virial     => dynvars%virial
    viri_const => dynvars%virial_const
    bmoment    => dynvars%barostat_momentum

    ! save box size(t+dt) and compute volume(t+dt)
    !
    boundary%box_size_x_ref = boundary%box_size_x
    boundary%box_size_y_ref = boundary%box_size_y
    boundary%box_size_z_ref = boundary%box_size_z
    volume = boundary%box_size_x * boundary%box_size_y * boundary%box_size_z

    ! initial constraint virial
    !
    virial_constraint(1:3,1:3) = 0.0d0

    ! time step
    !
    half_dt   = dt_short * 0.5_dp
    quart_dt  = dt_short * 0.25_dp

    ! maximum iteration
    !
    if (constraints%rigid_bond) then
      maxiter = 4
    else
      maxiter = 1
    end if

    ! initial barostat
    !
    bmoment_ref(1:3) = bmoment(1:3)

    ! initial constraint force
    !
    do j = 1, ncell
      do jx = 1, natom(j)
        force_const(1:3,jx,j) = 0.0_dp
      end do
    end do

    ! current velocity
    !
    do j = 1, ncell
      do jx =1 ,natom(j)
        vel_ref(1:3,jx,j) = vel(1:3,jx,j)
      end do
    end do

    ! iteration of barostat
    !
    do i = 1, maxiter

      ! VV2
      !
      do j = 1, ncell
        do jx = 1, nsolute(j)
          factor = half_dt / mass(jx,j)
          vel(1:3,jx,j) = vel_ref(1:3,jx,j) + factor*force_short(1:3,jx,j)
          vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_const(1:3,jx,j)
        end do
        do l = 1, nwater(j)
          do k = 1, 3
            jx = water_list(k,l,j)
            factor = half_dt / mass(jx,j)
            vel(1:3,jx,j) = vel_ref(1:3,jx,j) + factor*force_short(1:3,jx,j)
            vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_const(1:3,jx,j)
          end do
        end do
      end do

      ! RATTLE VV2
      !
      if (constraints%rigid_bond) then
        do j = 1, ncell
          do jx = 1, natom(j)
            coord_ref(1:3,jx,j) = vel(1:3,jx,j)
          end do
        end do
        call compute_constraints(ConstraintModeVVER2, .false., dt, coord_ref, &
                                 domain, constraints, coord, vel, &
                                 viri_const)
        do j = 1, ncell
          do jx = 1, natom(j)
            vel_change(1:3) = vel(1:3,jx,j) - coord_ref(1:3,jx,j)
            force_const(1:3,jx,j) = force_const(1:3,jx,j)                  &
                                    + mass(jx,j)*vel_change(1:3)/half_dt
          end do
        end do

        ! constraint virial
        ! 
        if (i == maxiter) then
          do j = 1, ncell
            do jx = 1, natom(j)
              virial_constraint(1,1) = virial_constraint(1,1)             + &
                                       force_const(1,jx,j)*coord(1,jx,j)
              virial_constraint(2,2) = virial_constraint(2,2)             + &
                                       force_const(2,jx,j)*coord(2,jx,j)
              virial_constraint(3,3) = virial_constraint(3,3)             + &
                                       force_const(3,jx,j)*coord(3,jx,j)
            end do
          end do
        end if

      end if

      ! kinetic energy component
      !
      kin(1:3) = 0.0_dp
      do j = 1, ncell
        do jx = 1, natom(j)
          kin(1:3) = kin(1:3) + mass(jx,j)*vel(1:3,jx,j)*vel(1:3,jx,j)
        end do
      end do
      ekin = kin(1) + kin(2) + kin(3)

      ! virial + virial+constraint
      !
      virial_sum(1) = virial(1,1) + virial_constraint(1,1)
      virial_sum(2) = virial(2,2) + virial_constraint(2,2)
      virial_sum(3) = virial(3,3) + virial_constraint(3,3)

      call reduce_pres(kin, ekin, virial_sum)

      ! compute pressure
      !
      press(1:3) = (kin(1:3)+virial_sum(1:3)) / volume
      pressxyz   = (press(1)+press(2)+press(3)) / 3.0_dp
      pressxy    = (press(1)+press(2)) / 2.0_dp
      pressz     = press(3)
      press0     = press0 * ensemble%pressure_short

      ! update barostat and velocity scale
      !
      if (ensemble%isotropy == IsotropyISO) then
        baro_force(1:3) = 3.0_dp*(ekin/degree+volume*(pressxyz-press0))/pmass
      else if (ensemble%isotropy == IsotropyANISO) then
        baro_force(1:3) = 3.0_dp*(ekin/degree+volume*(press(1:3)-press0))/pmass
      else if (ensemble%isotropy == IsotropySEMI_ISO) then
        baro_force(1:2) = 3.0_dp*(ekin/degree+volume*(pressxy-press0))/pmass
        baro_force(3)   = 3.0_dp*(ekin/degree+volume*(press(3)-press0))/pmass
      else if (ensemble%isotropy == IsotropyXY_Fixed) then
        baro_force(1:2) = 0.0_dp
        baro_force(3) = 3.0_dp*(ekin/degree+volume*(press(3)-press0))/pmass
      end if
      bmoment(1:3) = bmoment_ref(1:3) + baro_force(1:3)*half_dt

    end do

    viri_const(1,1) = virial_constraint(1,1)
    viri_const(2,2) = virial_constraint(2,2)
    viri_const(3,3) = virial_constraint(3,3)

    ! thermostat
    !
    call mtk_thermostat(dynamics, dt_short, ensemble, domain, dynvars)

    return

  end subroutine mtk_barostat_vv2

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    mtk_thermostat
  !> @brief        control temperature of particles and barostats
  !! @authors      JJ
  !! @param[in]    dynamics : dynamics information
  !! @param[in]    dt       : time step
  !! @param[in]    ensemble : ensemble information
  !! @param[inout] domain   : domain information
  !! @param[inout] dynvars  : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine mtk_thermostat(dynamics, dt, ensemble, domain, dynvars)

    ! formal arguments
    type(s_dynamics),         intent(in)    :: dynamics
    real(dp),                 intent(in)    :: dt
    type(s_ensemble), target, intent(in)    :: ensemble
    type(s_domain),   target, intent(inout) :: domain
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    real(dp)                 :: kin(1:3), w(1:3)
    real(dp)                 :: temp0, tau_t, dt_small
    real(dp)                 :: dt_1, dt_2, dt_4, dt_8
    real(dp)                 :: ekin_ptl, ekin_baro
    real(dp)                 :: KbT, degree
    real(dp)                 :: scale_kin, scale_ptl, scale_baro
    integer                  :: i, j, k, ix
    integer                  :: num_degree, ncell
    integer                  :: nh_length, nh_step

    real(dp),        pointer :: nh_mass(:), nh_vel(:)
    real(dp),        pointer :: nh_force(:), nh_coef(:)
    real(dp),        pointer :: nh_baro_vel(:), nh_baro_coef(:)
    real(dp),        pointer :: nh_baro_force(:)
    real(dp),        pointer :: vel(:,:,:), mass(:,:)
    real(dp),        pointer :: pmass, bmoment(:)
    integer,         pointer :: natom(:)


    temp0         =  ensemble%temperature
    tau_t         =  ensemble%tau_t/AKMA_PS
    nh_length     =  ensemble%nhchain
    nh_step       =  ensemble%nhmultistep
    num_degree    =  domain%num_deg_freedom
    ncell         =  domain%num_cell_local
    KbT           =  KBOLTZ * temp0
    degree        =  real(num_degree, dp)

    pmass         => ensemble%pmass
    natom         => domain%num_atom
    vel           => domain%velocity
    mass          => domain%mass
    nh_mass       => dynvars%nh_mass
    nh_vel        => dynvars%nh_velocity
    nh_force      => dynvars%nh_force
    nh_coef       => dynvars%nh_coef
    nh_baro_vel   => dynvars%nh_baro_velocity
    nh_baro_force => dynvars%nh_baro_force
    nh_baro_coef  => dynvars%nh_baro_coef
    bmoment       => dynvars%barostat_momentum

    ! Yoshida coefficient
    !
    w(1) = 1.0_dp / (2.0_dp - 2.0_dp**(1.0_dp/3.0_dp))
    w(3) = w(1)
    w(2) = 1.0_dp - w(1) - w(3)

    ! temperature scale factor
    !
    scale_ptl  = 1.0_dp
    scale_baro = 1.0_dp

    ! calculate kinetic energy
    !
    kin(1:3) = 0.0_dp
    do i = 1, ncell
      do ix = 1, natom(i)
        kin(1:3) = kin(1:3) + mass(ix,i)*vel(1:3,ix,i)**2
      end do
    end do
    ekin_ptl = kin(1) + kin(2) + kin(3)

    ! barostat kinetic energy
    !
    kin(1:3) = 0.0_dp
    kin(1:3) = kin(1:3) + pmass*bmoment(1:3)**2
    ekin_baro = kin(1) + kin(2) + kin(3)
    ekin_baro = ekin_baro / 3.0_dp

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, ekin_ptl, 1, mpi_real8, mpi_sum, &
                       mpi_comm_country, ierror)
    call mpi_allreduce(mpi_in_place, ekin_baro, 1, mpi_real8, mpi_sum, &
                       mpi_comm_country, ierror)
#endif

    ! decide Nose-Hoover theremostat chain momentum
    !
    dt_small = dt / real(nh_step)

    do i = 1, nh_step
      do j = 1, 3

        dt_1 = w(j) * dt_small
        dt_2 = dt_1 * 0.5_dp
        dt_4 = dt_2 * 0.5_dp
        dt_8 = dt_4 * 0.5_dp

        ! Nose-Hoover chains for particles
        !
        nh_force(nh_length) = nh_mass(nh_length-1)*nh_vel(nh_length-1)**2-KbT
        nh_force(nh_length) = nh_force(nh_length) / nh_mass(nh_length)
        nh_vel(nh_length) = nh_vel(nh_length) + nh_force(nh_length)*dt_4

        do k = nh_length-1, 2, -1
          nh_force(k) = (nh_mass(k-1)*nh_vel(k-1)**2-KbT) / nh_mass(k)
          nh_coef(k)  = exp(-nh_vel(k+1)*dt_8)
          nh_vel(k)   = nh_vel(k) * nh_coef(k)
          nh_vel(k)   = nh_vel(k) + nh_force(k)*dt_4
          nh_vel(k)   = nh_vel(k) * nh_coef(k)
        end do

        nh_force(1) = (ekin_ptl - degree*KbT) / nh_mass(1)
        nh_coef(1)  = exp(-nh_vel(2)*dt_8)
        nh_vel(1)   = nh_vel(1) * nh_coef(1)
        nh_vel(1)   = nh_vel(1) + nh_force(1)*dt_4
        nh_vel(1)   = nh_vel(1) * nh_coef(1)

        ! Nose-Hoover chains for barostat
        !
        nh_baro_force(nh_length) = nh_mass(nh_length-1)*nh_baro_vel(nh_length-1)**2-KbT
        nh_baro_force(nh_length) = nh_baro_force(nh_length) / nh_mass(nh_length)
        nh_baro_vel(nh_length) = nh_baro_vel(nh_length) + nh_baro_force(nh_length)*dt_4

        do k = nh_length-1, 2, -1
          nh_baro_force(k) = (nh_mass(k-1)*nh_baro_vel(k-1)**2-KbT) / nh_mass(k)
          nh_baro_coef(k)  = exp(-nh_baro_vel(k+1)*dt_8)
          nh_baro_vel(k)   = nh_baro_vel(k) * nh_baro_coef(k)
          nh_baro_vel(k)   = nh_baro_vel(k) + nh_baro_force(k)*dt_4
          nh_baro_vel(k)   = nh_baro_vel(k) * nh_baro_coef(k)
        end do
        nh_baro_force(1) = (ekin_baro - KbT) / nh_mass(1)
        nh_baro_coef(1)  = exp(-nh_baro_vel(2)*dt_8)
        nh_baro_vel(1)   = nh_baro_vel(1) * nh_baro_coef(1)
        nh_baro_vel(1)   = nh_baro_vel(1) + nh_baro_force(1)*dt_4
        nh_baro_vel(1)   = nh_baro_vel(1) * nh_baro_coef(1)

        ! scale kinetic energy and barostat kinetic energy
        !
        scale_kin  = exp(-nh_vel(1)*dt_1)
        ekin_ptl   = ekin_ptl * scale_kin
        scale_ptl  = scale_ptl * exp(-nh_vel(1)*dt_2)

        scale_kin  = exp(-nh_baro_vel(1)*dt_1)
        ekin_baro  = ekin_baro * scale_kin
        scale_baro = scale_baro * exp(-nh_baro_vel(1)*dt_2)

        ! Nose-Hoover chains for particles 
        !
        nh_force(1) = (ekin_ptl - degree*KbT) / nh_mass(1)
        nh_vel(1)   = nh_vel(1) * nh_coef(1)
        nh_vel(1)   = nh_vel(1) + nh_force(1)*dt_4
        nh_vel(1)   = nh_vel(1) * nh_coef(1)

        do k = 2, nh_length-1
          nh_force(k) = (nh_mass(k-1)*nh_vel(k-1)**2-KbT) / nh_mass(k)
          nh_vel(k)   = nh_vel(k) * nh_coef(k)
          nh_vel(k)   = nh_vel(k) + nh_force(k)*dt_4
          nh_vel(k)   = nh_vel(k) * nh_coef(k)
        end do

        nh_force(nh_length) = nh_mass(nh_length-1)*nh_vel(nh_length-1)**2-KbT
        nh_force(nh_length) = nh_force(nh_length) / nh_mass(nh_length)
        nh_vel(nh_length) = nh_vel(nh_length) + nh_force(nh_length)*dt_4

        nh_baro_force(1) = (ekin_baro - KbT) / nh_mass(1)
        nh_baro_vel(1)   = nh_baro_vel(1) * nh_baro_coef(1)
        nh_baro_vel(1)   = nh_baro_vel(1) + nh_baro_force(1)*dt_4
        nh_baro_vel(1)   = nh_baro_vel(1) * nh_baro_coef(1)

        do k = 2, nh_length-1
          nh_baro_force(k) = (nh_mass(k-1)*nh_baro_vel(k-1)**2-KbT) / nh_mass(k)
          nh_baro_vel(k)   = nh_baro_vel(k) * nh_baro_coef(k)
          nh_baro_vel(k)   = nh_baro_vel(k) + nh_baro_force(k)*dt_4
          nh_baro_vel(k)   = nh_baro_vel(k) * nh_baro_coef(k)
        end do

        nh_baro_force(nh_length) = nh_mass(nh_length-1)*nh_baro_vel(nh_length-1)**2-KbT
        nh_baro_force(nh_length) = nh_baro_force(nh_length) / nh_mass(nh_length)
        nh_baro_vel(nh_length) = nh_baro_vel(nh_length) + nh_baro_force(nh_length)*dt_4

      end do
    end do

    ! velocity scaling
    ! 
    do i = 1, ncell
      do ix = 1, natom(i)
        vel(1:3,ix,i) = scale_ptl*vel(1:3,ix,i)
      end do
    end do
    bmoment(1:3) = scale_baro*bmoment(1:3)

    return

  end subroutine mtk_thermostat

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    langevin_thermostat_vv1
  !> @brief        control temperature using Langevin thermostat
  !! @authors      JJ
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine langevin_thermostat_vv1(dynamics, istep, inner_step, &
                                     dt_long, dt_short, ensemble, &
                                     domain, constraints, dynvars)

    ! formal arguments
    type(s_dynamics),         intent(in)    :: dynamics
    integer,                  intent(in)    :: istep
    integer,                  intent(in)    :: inner_step
    real(dp),                 intent(in)    :: dt_long, dt_short
    type(s_ensemble), target, intent(in)    :: ensemble
    type(s_domain),   target, intent(inout) :: domain
    type(s_constraints),      intent(inout) :: constraints
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    real(dp)                 :: imass, half_dt, inv_dt
    real(dp)                 :: dt_therm, half_dt_therm
    real(dp)                 :: temp0, gamma_t, scale_v
    real(dp)                 :: factor, sigma
    real(dp)                 :: rsq, v1, v2, grandom(1:3)
    real(dp)                 :: kBT, vel_tmp(1:3)
    integer                  :: j, jx, k, l, ncell

    real(dp),        pointer :: vel(:,:,:), vel_ref(:,:,:)
    real(dp),        pointer :: coord(:,:,:), coord_ref(:,:,:)
    real(dp),        pointer :: force_long(:,:,:), force_short(:,:,:)
    real(dp),        pointer :: random_f(:,:,:), mass(:,:)
    real(dp),        pointer :: temporary(:,:,:)
    real(dp),        pointer :: virial(:,:), viri_const(:,:)
    integer,         pointer :: natom(:), nsolute(:), nwater(:)
    integer,         pointer :: water_list(:,:,:)


    inv_dt        =  1.0_dp / dt_short
    dt_therm      =  dt_short * real(dynamics%thermo_period, dp)
    half_dt_therm =  dt_therm / 2.0_dp
    temp0         =  ensemble%temperature
    gamma_t       =  ensemble%gamma_t *AKMA_PS
    random_f      => ensemble%random_force
    ncell         =  domain%num_cell_local

    natom         => domain%num_atom
    nsolute       => domain%num_solute
    nwater        => domain%num_water 
    water_list    => domain%water_list
    mass          => domain%mass
    coord         => domain%coord
    coord_ref     => domain%coord_ref
    vel           => domain%velocity
    vel_ref       => domain%velocity_ref
    force_long    => domain%force_long
    force_short   => domain%force_short
    temporary     => domain%coord_old
    virial        => dynvars%virial
    viri_const    => dynvars%virial_const

    ! setup variables
    !
    kBT      = KBOLTZ * temp0

    ! scale factor for velocities
    !
    scale_v = exp(-gamma_t*half_dt_therm)

    ! random force
    !
    if (istep == 1) then

      factor  = 1.0_dp - scale_v*scale_v
      factor  = factor*KBOLTZ*temp0

      do j = 1, ncell
        do jx = 1, natom(j)

          if (abs(mass(jx,j)) < EPS) cycle

          sigma = sqrt(factor/mass(jx,j))
          rsq = 2.0_dp

          do while (rsq >= 1.0_dp)
            v1  = 2.0_dp*random_get() - 1.0_dp
            v2  = 2.0_dp*random_get() - 1.0_dp
            rsq = v1*v1 + v2*v2
          end do

          rsq   = sqrt(-2.0_dp * log(rsq) / rsq)
          grandom(1) = rsq * v1
          rsq = 2.0_dp

          do while (rsq >= 1.0_dp)
            v1  = 2.0_dp*random_get() - 1.0_dp
            v2  = 2.0_dp*random_get() - 1.0_dp
            rsq = v1*v1 + v2*v2
          end do

          rsq   = sqrt(-2.0_dp * log(rsq) / rsq)
          grandom(2) = rsq * v1
          rsq = 2.0_dp

          do while (rsq >= 1.0_dp)
            v1  = 2.0_dp*random_get() - 1.0_dp
            v2  = 2.0_dp*random_get() - 1.0_dp
            rsq = v1*v1 + v2*v2
          end do

          rsq   = sqrt(-2.0_dp * log(rsq) / rsq)
          grandom(3) = rsq * v1
          random_f(1:3,jx,j) = sigma*grandom(1:3)

        end do
      end do

    end if

    ! Thermostat before VV1
    !
    if (mod(istep-1, dynamics%thermo_period) == 0) then

      ! Thermostat
      !
      do j = 1, ncell
        do jx = 1, natom(j)
          vel(1:3,jx,j) = vel_ref(1:3,jx,j)*scale_v
          vel(1:3,jx,j) = vel(1:3,jx,j) + random_f(1:3,jx,j)
        end do
      end do

      ! FEP
      ! synchronize single B with single A
      if (domain%fep_use) call sync_single_fep(domain, vel)

    end if

    ! VV1 with long range force
    !
    if (inner_step == 1) then
      do j = 1, ncell
        do jx = 1, nsolute(j)
          factor = 0.5_dp * dt_long / mass(jx,j)
          vel  (1:3,jx,j) = vel(1:3,jx,j) + factor*force_long(1:3,jx,j)
        end do
        do l = 1, nwater(j)
          do k = 1, 3
            jx = water_list(k,l,j)
            factor = 0.5_dp * dt_long / mass(jx,j)
            vel  (1:3,jx,j) = vel(1:3,jx,j) + factor*force_long(1:3,jx,j)
          end do
        end do
      end do
    end if

    ! VV1
    !
    do j = 1, ncell
      do jx = 1, nsolute(j)
        factor = 0.5_dp * dt_short / mass(jx,j)
        vel  (1:3,jx,j) = vel(1:3,jx,j) + factor*force_short(1:3,jx,j)
        coord(1:3,jx,j) = coord_ref(1:3,jx,j)+vel(1:3,jx,j)*dt_short
      end do
      do l = 1, nwater(j)
        do k = 1, 3
          jx = water_list(k,l,j)
          factor = 0.5_dp * dt_short / mass(jx,j)
          vel  (1:3,jx,j) = vel(1:3,jx,j) + factor*force_short(1:3,jx,j)
          coord(1:3,jx,j) = coord_ref(1:3,jx,j)+vel(1:3,jx,j)*dt_short
        end do
      end do
    end do

    ! Coordinate constraint (RATTLE VV1)
    !
    if (constraints%rigid_bond) then
      do j = 1, ncell
        do jx = 1, natom(j)
          temporary(1,jx,j) = coord(1,jx,j)
          temporary(2,jx,j) = coord(2,jx,j)
          temporary(3,jx,j) = coord(3,jx,j)
        end do
      end do

      call compute_constraints(ConstraintModeLEAP, .false., dt_short, &
                               coord_ref, domain, constraints, coord, &
                               vel, viri_const)
      dynvars%virial_const(1:3,1:3) = 2.0_dp * dynvars%virial_const(1:3,1:3)

      do j = 1, ncell
        do jx = 1, natom(j)
          vel_tmp(1:3)  = (coord(1:3,jx,j) - temporary(1:3,jx,j))*inv_dt
          vel(1:3,jx,j) = vel(1:3,jx,j) + vel_tmp(1:3)
        end do
      end do
    end if

    return

  end subroutine langevin_thermostat_vv1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    langevin_thermostat_vv2
  !> @brief        Langevin thermostat and barostat
  !! @authors      JJ
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine langevin_thermostat_vv2(dynamics, istep, inner_step, dt_long,    &
                                     dt_short, ensemble, domain, constraints, &
                                     dynvars)

    ! formal arguments
    type(s_dynamics),         intent(in)    :: dynamics
    integer,                  intent(in)    :: istep
    integer,                  intent(in)    :: inner_step
    real(dp),                 intent(in)    :: dt_long
    real(dp),                 intent(in)    :: dt_short
    type(s_ensemble), target, intent(in)    :: ensemble
    type(s_domain),   target, intent(inout) :: domain
    type(s_constraints),      intent(inout) :: constraints
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    real(dp)                 :: inv_dt, temp0
    real(dp)                 :: scale_v, factor
    real(dp)                 :: gamma_t
    real(dp)                 :: sigma
    real(dp)                 :: v1, v2, rsq, grandom(1:3)
    real(dp)                 :: half_dt, quart_dt
    real(dp)                 :: dt_therm, half_dt_therm
    integer                  :: i, j, jx, k, l, ncell

    real(dp),        pointer :: mass(:,:), viri_const(:,:)
    real(dp),        pointer :: random_f(:,:,:)
    real(dp),        pointer :: vel(:,:,:), vel_ref(:,:,:)
    real(dp),        pointer :: force_long(:,:,:), force_short(:,:,:)
    real(dp),        pointer :: coord(:,:,:), coord_ref(:,:,:)
    integer,         pointer :: natom(:), nsolute(:), nwater(:)
    integer,         pointer :: water_list(:,:,:)


    inv_dt     =  1.0_dp/dt_short
    temp0      =  ensemble%temperature
    gamma_t    =  ensemble%gamma_t * AKMA_PS
    random_f   => ensemble%random_force
    ncell      =  domain%num_cell_local

    mass       => domain%mass
    natom      => domain%num_atom
    nsolute    => domain%num_solute
    nwater     => domain%num_water 
    water_list => domain%water_list
    coord      => domain%coord
    coord_ref  => domain%coord_ref
    vel        => domain%velocity
    vel_ref    => domain%velocity_ref
    force_long => domain%force_long
    force_short=> domain%force_short
    viri_const => dynvars%virial_const

    ! time step
    !
    dt_therm       = dt_short * real(dynamics%thermo_period,dp)
    half_dt_therm  = dt_therm / 2.0_dp

    ! VV2
    !
    do j = 1, ncell
      do jx = 1, nsolute(j)
        factor = 0.5_dp * dt_short/mass(jx,j)
        vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_short(1:3,jx,j)
      end do
      do l = 1, nwater(j)
        do k = 1, 3
          jx = water_list(k,l,j)
          factor = 0.5_dp * dt_short/mass(jx,j)
          vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_short(1:3,jx,j)
        end do
      end do
    end do

    if (inner_step == dynamics%elec_long_period) then

      ! VV2 (long range force)
      !
      do j = 1, ncell
        do jx = 1, nsolute(j)
          factor = 0.5_dp * dt_long / mass(jx,j)
          vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_long(1:3,jx,j)
        end do
        do l = 1, nwater(j)
          do k = 1, 3
            jx = water_list(k,l,j)
            factor = 0.5_dp * dt_long / mass(jx,j)
            vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_long(1:3,jx,j)
          end do
        end do
      end do

      if (mod(istep,  dynamics%thermo_period) == 0) then

        ! random force
        !
        scale_v = exp(-gamma_t*half_dt_therm)
        factor   = 1.0_dp - scale_v*scale_v
        factor   = factor*KBOLTZ*temp0/2.0_dp

        do j = 1, ncell
          do jx = 1, natom(j)

            if (abs(mass(jx,j)) < EPS)  cycle

            sigma = sqrt(factor/mass(jx,j))
            rsq = 2.0_dp

            do while (rsq >= 1.0_dp)
              v1  = 2.0_dp*random_get() - 1.0_dp
              v2  = 2.0_dp*random_get() - 1.0_dp
              rsq = v1*v1 + v2*v2
            end do
  
            rsq   = sqrt(-2.0_dp * log(rsq) / rsq)
            grandom(1) = rsq * v1
            rsq = 2.0_dp
  
            do while (rsq >= 1.0_dp)
              v1  = 2.0_dp*random_get() - 1.0_dp
              v2  = 2.0_dp*random_get() - 1.0_dp
              rsq = v1*v1 + v2*v2
            end do
  
            rsq   = sqrt(-2.0_dp * log(rsq) / rsq)
            grandom(2) = rsq * v1
            rsq = 2.0_dp
  
            do while (rsq >= 1.0_dp)
              v1  = 2.0_dp*random_get() - 1.0_dp
              v2  = 2.0_dp*random_get() - 1.0_dp
              rsq = v1*v1 + v2*v2
            end do
  
            rsq   = sqrt(-2.0_dp * log(rsq) / rsq)
            grandom(3) = rsq * v1
            random_f(1:3,jx,j) = sigma*grandom(1:3)
  
          end do
        end do

        ! Thermostat
        !
        do j = 1, ncell
          do jx = 1, natom(j)
            vel(1:3,jx,j) = vel(1:3,jx,j) * scale_v
            vel(1:3,jx,j) = vel(1:3,jx,j) + random_f(1:3,jx,j)
          end do
        end do

        ! FEP
        ! synchronize single B with single A
        if (domain%fep_use) call sync_single_fep(domain, vel)

      end if

    end if

    ! RATTLE VV2
    !
    if (constraints%rigid_bond) &
      call compute_constraints(ConstraintModeVVER2, .false., dt_short, &
                               coord_ref, domain, constraints, coord,  &
                               vel, viri_const)

    return

  end subroutine langevin_thermostat_vv2


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    langevin_barostat_vv1
  !> @brief        Langevin thermostat and barostat
  !! @authors      JJ
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] boundary    : boundary information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine langevin_barostat_vv1(dynamics, istep, inner_step, dt_long,    &
                                   dt_short, ensemble, domain, constraints, &
                                   boundary, dynvars)

    ! formal arguments
    type(s_dynamics),         intent(in)    :: dynamics
    integer,                  intent(in)    :: istep
    integer,                  intent(in)    :: inner_step
    real(dp),                 intent(in)    :: dt_long
    real(dp),                 intent(in)    :: dt_short
    type(s_ensemble), target, intent(inout) :: ensemble
    type(s_domain),   target, intent(inout) :: domain
    type(s_constraints),      intent(inout) :: constraints
    type(s_boundary),         intent(inout) :: boundary
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    real(dp)                 :: inv_dt, temp0, press0, d_ndegf
    real(dp)                 :: kin(1:3), kin_temp(1:3), ekin
    real(dp)                 :: delta_vel(1:3)
    real(dp)                 :: volume, press(1:3)
    real(dp)                 :: pressxy, pressxyz
    real(dp)                 :: crdx, crdy, crdz, factor
    real(dp)                 :: bmoment_ref(3), scale_b(1:3)
    real(dp)                 :: gamma_t, gamma_p
    real(dp)                 :: sigma
    real(dp)                 :: v1, v2, rsq, grandom(1:3)
    real(dp)                 :: half_dt, quart_dt, half_dt_short
    real(dp)                 :: dt_baro, half_dt_baro, quart_dt_baro
    real(dp)                 :: dt_therm, half_dt_therm
    real(dp)                 :: size_scale(1:3), vel_scale
    real(dp)                 :: virial_constraint(3,3), virial_sum(3)
    integer                  :: i, j, ij, jx, k, l, i_ndegf, maxiter
    integer                  :: ncell, nboundary

    real(dp),        pointer :: pmass, pforce(:), random_f(:,:,:)
    real(dp),        pointer :: mass(:,:)
    real(dp),        pointer :: coord(:,:,:), coord_ref(:,:,:)
    real(dp),        pointer :: temporary(:,:,:), temporary1(:,:,:)
    real(dp),        pointer :: vel(:,:,:), vel_ref(:,:,:)
    real(dp),        pointer :: force_long(:,:,:), force_short(:,:,:)
    real(dp),        pointer :: virial_short(:,:), virial_long(:,:)
    real(dp),        pointer :: viri_const(:,:), kin_ref(:)
    real(dp),        pointer :: bmoment(:)
    integer,         pointer :: natom(:), nsolute(:), nwater(:)
    integer,         pointer :: water_list(:,:,:)


    inv_dt       =  1.0_dp/dt_short
    temp0        =  ensemble%temperature
    press0       =  ensemble%pressure * ATMOS_P
    gamma_t      =  ensemble%gamma_t * AKMA_PS
    gamma_p      =  ensemble%gamma_p * AKMA_PS
    pmass        => ensemble%pmass
    pforce       => ensemble%pforce
    random_f     => ensemble%random_force
    i_ndegf      =  domain%num_deg_freedom
    d_ndegf      =  real(i_ndegf,dp)
    ncell        =  domain%num_cell_local
    nboundary    =  domain%num_cell_boundary

    natom        => domain%num_atom
    nsolute      => domain%num_solute
    nwater       => domain%num_water 
    water_list   => domain%water_list
    mass         => domain%mass
    coord        => domain%coord
    coord_ref    => domain%coord_ref
    vel          => domain%velocity
    vel_ref      => domain%velocity_ref
    force_long   => domain%force_long
    force_short  => domain%force_short
    temporary    => domain%coord_old
    temporary1   => domain%velocity_full
    virial_short => dynvars%virial
    virial_long  => dynvars%virial_long
    viri_const   => dynvars%virial_const
    bmoment      => dynvars%barostat_momentum
    kin_ref      => dynvars%kinetic_ref

    ! save box size(t+dt) and compute volume(t+dt)
    !
    boundary%box_size_x_ref = boundary%box_size_x
    boundary%box_size_y_ref = boundary%box_size_y
    boundary%box_size_z_ref = boundary%box_size_z
    volume = boundary%box_size_x * boundary%box_size_y * boundary%box_size_z

    ! time step
    !
    half_dt_short = dt_short / 2.0_dp
    dt_therm      = dt_short * real(dynamics%thermo_period, dp)
    half_dt_therm = dt_therm / 2.0_dp
    dt_baro       = dt_short * real(dynamics%baro_period, dp)
    half_dt_baro  = dt_baro / 2.0_dp
    quart_dt_baro = half_dt_baro / 2.0_dp

    ! scale factor for veloctiy rescaling
    !
    vel_scale = exp(-gamma_t*half_dt_therm)

    ! maximum iteration
    !
    if (mod(istep-1,dynamics%baro_period) == 0) then
      if (constraints%rigid_bond) then
        maxiter = 4
      else
        maxiter = 1
      end if
    else
      maxiter = 1
    end if

    ! constraint virial
    !
    virial_constraint(1:3,1:3) = 0.0_dp

    ! barostate coefficient
    !
    bmoment_ref(1:3) = bmoment(1:3)

    ! pmass and stochastic force (pforce)
    !
    if (istep == 1) then

      if (ensemble%isotropy == IsotropyISO) then

        pmass = (d_ndegf+3.0_dp)*KBOLTZ*temp0 / (2.0_dp*PI*gamma_p)**2
        sigma = sqrt(gamma_p*pmass*KBOLTZ*temp0/quart_dt_baro)
        pforce(1) = sigma * random_get_gauss()
        pforce(2) = pforce(1)
        pforce(3) = pforce(1)

      else if (ensemble%isotropy == IsotropySEMI_ISO) then

        pmass = (d_ndegf+3.0_dp)*KBOLTZ*temp0 / (3.0_dp*(2.0_dp*PI*gamma_p)**2)
        sigma = sqrt(gamma_p*pmass*KBOLTZ*temp0/quart_dt_baro)
        pforce(1) = sigma * random_get_gauss()
        pforce(2) = pforce(1)
        pforce(3) = sigma * random_get_gauss()

      else if (ensemble%isotropy == IsotropyANISO) then

        pmass = (d_ndegf+3.0_dp)*KBOLTZ*temp0 / (3.0_dp*(2.0_dp*PI*gamma_p)**2)
        sigma = sqrt(gamma_p*pmass*KBOLTZ*temp0/quart_dt_baro)
        pforce(1) = sigma * random_get_gauss()
        pforce(2) = sigma * random_get_gauss()
        pforce(3) = sigma * random_get_gauss()

      else if (ensemble%isotropy == IsotropyXY_Fixed) then

        pmass = (d_ndegf+1.0_dp)*KBOLTZ*temp0 / (3.0_dp*(2.0_dp*PI*gamma_p)**2)
        sigma = sqrt(gamma_p*pmass*KBOLTZ*temp0/quart_dt_baro)
        pforce(1) = 0.0_dp
        pforce(2) = 0.0_dp
        pforce(3) = sigma * random_get_gauss()

      end if

#ifdef HAVE_MPI_GENESIS
      call mpi_bcast(pforce, 3, mpi_real8, 0, mpi_comm_country, ierror)
#endif

      ! random force
      !
      factor = 1.0_dp - vel_scale*vel_scale
      factor = factor*KBOLTZ*temp0/2.0_dp

      do j = 1, ncell
        do jx = 1, natom(j)

          if (abs(mass(jx,j)) < EPS) cycle

          sigma = sqrt(factor/mass(jx,j))
          rsq = 2.0_dp

          do while (rsq >= 1.0_dp)
            v1  = 2.0_dp*random_get() - 1.0_dp
            v2  = 2.0_dp*random_get() - 1.0_dp
            rsq = v1*v1 + v2*v2
          end do

          rsq   = sqrt(-2.0_dp * log(rsq) / rsq)
          grandom(1) = rsq * v1
          rsq = 2.0_dp
          do while (rsq >= 1.0_dp)
            v1  = 2.0_dp*random_get() - 1.0_dp
            v2  = 2.0_dp*random_get() - 1.0_dp
            rsq = v1*v1 + v2*v2
          end do

          rsq   = sqrt(-2.0_dp * log(rsq) / rsq)
          grandom(2) = rsq * v1
          rsq = 2.0_dp

          do while (rsq >= 1.0_dp)
            v1  = 2.0_dp*random_get() - 1.0_dp
            v2  = 2.0_dp*random_get() - 1.0_dp
            rsq = v1*v1 + v2*v2
          end do

          rsq   = sqrt(-2.0_dp * log(rsq) / rsq)
          grandom(3) = rsq * v1
          random_f(1:3,jx,j) = sigma*grandom(1:3)

        end do
      end do

      kin_ref(1:3) = 0.0_dp
      do j = 1, ncell
        do jx = 1, natom(j)

          ! FEP
          ! skip singleB to avoid duplication
          if (domain%fep_use) then
            if (domain%fepgrp(jx,j) == 2) cycle
          end if

          kin_ref(1:3) = kin_ref(1:3) &
                       + mass(jx,j)*vel(1:3,jx,j)**2
        end do
      end do
      call mpi_allreduce(mpi_in_place, kin_ref, 3, mpi_real8, mpi_sum, &
                         mpi_comm_country, ierror)

    end if

    do j = 1, ncell
      do jx = 1, natom(j)
        vel_ref(1:3,jx,j) = vel(1:3,jx,j)
        temporary1(1:3,jx,j) = 0.0_dp
      end do
    end do

    if (mod(istep-1,dynamics%baro_period) == 0) then

      do i = 1, maxiter

        ! Barostat 1
        !
        kin_temp(1:3) = 0.0_dp
        do j = 1, ncell
          do jx = 1, natom(j)

            ! FEP
            ! skip singleB to avoid duplication
            if (domain%fep_use) then
              if (domain%fepgrp(jx,j) == 2) cycle
            end if

            kin_temp(1:3) = kin_temp(1:3)  &
                          + mass(jx,j)*vel(1:3,jx,j)*vel(1:3,jx,j)
          end do
        end do

        ! virial + virial_constraint
        !
        virial_sum(1) = virial_short(1,1) + virial_long(1,1)
        virial_sum(2) = virial_short(2,2) + virial_long(2,2)
        virial_sum(3) = virial_short(3,3) + virial_long(3,3)
        virial_sum(1) = virial_sum(1) &
                      + virial_constraint(1,1)
        virial_sum(2) = virial_sum(2) &
                      + virial_constraint(2,2)
        virial_sum(3) = virial_sum(3) &
                      + virial_constraint(3,3)

        call reduce_pres(kin_temp, ekin, virial_sum)
        kin(1:3) = 0.5_dp*(kin_ref(1:3)+kin_temp(1:3))
        ekin = 0.5_dp*(kin(1)+kin(2)+kin(3))

        press(1:3) = (kin(1:3) + virial_sum(1:3))/volume
        pressxyz = (press(1) + press(2) + press(3))/3.0_dp
        pressxy  = (press(1) + press(2))/2.0_dp

        ! update barostat
        !
        do j = 1, ncell
          do jx = 1, natom(j)
            vel(1:3,jx,j) = vel_ref(1:3,jx,j)
          end do
        end do
        call update_barostat(ensemble, boundary, bmoment_ref, pforce,   &
                             press, pressxyz, pressxy, press0, volume,  &
                             d_ndegf, pmass, gamma_p, ekin,             &
                             dt_baro, half_dt_baro, natom, ncell, vel,  &
                             bmoment)
        size_scale(1:3) = exp(bmoment(1:3)*dt_baro)

        do j = 1, ncell
          do jx = 1, natom(j)
            vel(1:3,jx,j) = vel(1:3,jx,j) * vel_scale
            vel(1:3,jx,j) = vel(1:3,jx,j) + random_f(1:3,jx,j)
          end do
        end do

        ! FEP
        ! synchronize single B with single A
        if (domain%fep_use) call sync_single_fep(domain, vel)

        do j = 1, ncell
          do jx = 1, nsolute(j)
            factor = half_dt_short / mass(jx,j)
            vel  (1:3,jx,j) = vel(1:3,jx,j) + factor*force_short(1:3,jx,j)
            vel  (1:3,jx,j) = vel(1:3,jx,j) + factor*temporary1(1:3,jx,j)
            vel  (1:3,jx,j) = vel(1:3,jx,j) + factor*force_long(1:3,jx,j)
            coord(1:3,jx,j) = size_scale(1:3)*coord_ref(1:3,jx,j) &
                            + vel(1:3,jx,j)*dt_short
          end do
          do l = 1, nwater(j)
            do k = 1, 3
              jx = water_list(k,l,j)
              factor = half_dt_short / mass(jx,j)
              vel  (1:3,jx,j) = vel(1:3,jx,j) + factor*force_short(1:3,jx,j)
              vel  (1:3,jx,j) = vel(1:3,jx,j) + factor*temporary1(1:3,jx,j)
              vel  (1:3,jx,j) = vel(1:3,jx,j) + factor*force_long(1:3,jx,j)
              coord(1:3,jx,j) = size_scale(1:3)*coord_ref(1:3,jx,j) &
                              + vel(1:3,jx,j)*dt_short
            end do
          end do
        end do

        ! RATTLE VV1
        !
        if (constraints%rigid_bond) then
          do j = 1, ncell
            do jx = 1, natom(j)
              temporary(1:3,jx,j) = coord(1:3,jx,j)
            end do
          end do
          call compute_constraints(ConstraintModeLEAP, .false., dt_short, &
                                   coord_ref, domain, constraints, coord, &
                                   vel, viri_const)

          virial_constraint(1:3,1:3) = virial_constraint(1:3,1:3) &
                                     + 2.0_dp * viri_const(1:3,1:3)

          do j = 1, ncell
            do jx = 1, natom(j)
              delta_vel(1:3) = (coord(1:3,jx,j) - temporary(1:3,jx,j))*inv_dt
              vel(1:3,jx,j) = vel(1:3,jx,j) + delta_vel(1:3)
              temporary1(1:3,jx,j) = temporary1(1:3,jx,j) &
                                   + mass(jx,j)*delta_vel(1:3)/half_dt_short
            end do
          end do
        end if

      end do

      ! thermostat
      !
      do j = 1, ncell
        do jx = 1, natom(j)
          vel(1:3,jx,j) = vel_ref(1:3,jx,j) * vel_scale
          vel(1:3,jx,j) = vel(1:3,jx,j) + random_f(1:3,jx,j)
        end do
      end do

      ! FEP
      ! synchronize single B with single A
      if (domain%fep_use) call sync_single_fep(domain, vel)

      ! VV1 (long range force)
      !
      do j = 1, ncell
        do jx = 1, nsolute(j)
          factor = 0.5_dp * dt_long / mass(jx,j)
          vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_long(1:3,jx,j)
        end do
        do l = 1, nwater(j)
          do k = 1, 3
            jx = water_list(k,l,j)
            factor = 0.5_dp * dt_long / mass(jx,j)
            vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_long(1:3,jx,j)
          end do
        end do
      end do

      do j = 1, ncell
        do jx = 1, nsolute(j)
          factor = half_dt_short / mass(jx,j)
          vel  (1:3,jx,j) = vel(1:3,jx,j) + factor*force_short(1:3,jx,j)
          coord(1:3,jx,j) = size_scale(1:3)*coord_ref(1:3,jx,j) &
                          + vel(1:3,jx,j)*dt_short
        end do
        do l = 1, nwater(j)
          do k = 1, 3
            jx = water_list(k,l,j)
            factor = half_dt_short / mass(jx,j)
            vel  (1:3,jx,j) = vel(1:3,jx,j) + factor*force_short(1:3,jx,j)
            coord(1:3,jx,j) = size_scale(1:3)*coord_ref(1:3,jx,j) &
                            + vel(1:3,jx,j)*dt_short
          end do
        end do
      end do

      ! RATTLE VV1
      !
      if (constraints%rigid_bond) then
        do j = 1, ncell
          do jx = 1, natom(j)
            temporary(1:3,jx,j) = coord(1:3,jx,j)
          end do
        end do
        call compute_constraints(ConstraintModeLEAP, .false., dt_short, &
                                 coord_ref, domain, constraints, coord, &
                                 vel, viri_const)

        do j = 1, ncell
          do jx = 1, natom(j)
            delta_vel(1:3) = (coord(1:3,jx,j) - temporary(1:3,jx,j))*inv_dt
            vel(1:3,jx,j) = vel(1:3,jx,j) + delta_vel(1:3)
          end do
        end do
      end if

    else

      do j = 1, ncell
        do jx = 1, natom(j)
          vel(1:3,jx,j) = vel_ref(1:3,jx,j)
        end do
      end do

      ! Thermostat
      !
      if (mod(istep-1,dynamics%thermo_period) == 0) then
        do j = 1, ncell
          do jx = 1, natom(j)
            vel(1:3,jx,j) = vel(1:3,jx,j) * vel_scale
            vel(1:3,jx,j) = vel(1:3,jx,j) + random_f(1:3,jx,j)
          end do
        end do
      end if

      ! FEP
      ! synchronize single B with single A
      if (domain%fep_use) call sync_single_fep(domain, vel)

      do j = 1, ncell
        do jx = 1, nsolute(j)
          factor = half_dt_short / mass(jx,j)
          vel  (1:3,jx,j) = vel(1:3,jx,j) + factor*force_short(1:3,jx,j)
          vel  (1:3,jx,j) = vel(1:3,jx,j) + factor*force_long(1:3,jx,j)
          coord(1:3,jx,j) = coord_ref(1:3,jx,j) &
                          + vel(1:3,jx,j)*dt_short
        end do
        do l = 1, nwater(j)
          do k = 1, 3
            jx = water_list(k,l,j)
            factor = half_dt_short / mass(jx,j)
            vel  (1:3,jx,j) = vel(1:3,jx,j) + factor*force_short(1:3,jx,j)
            vel  (1:3,jx,j) = vel(1:3,jx,j) + factor*force_long(1:3,jx,j)
            coord(1:3,jx,j) = coord_ref(1:3,jx,j) &
                            + vel(1:3,jx,j)*dt_short
          end do
        end do
      end do

      ! RATTLE VV1
      !
      if (constraints%rigid_bond) then
        do j = 1, ncell
          do jx = 1, natom(j)
            temporary(1:3,jx,j) = coord(1:3,jx,j)
          end do
        end do
        call compute_constraints(ConstraintModeLEAP, .false., dt_short, &
                                 coord_ref, domain, constraints, coord, &
                                 vel, viri_const)

        do j = 1, ncell
          do jx = 1, natom(j)
            delta_vel(1:3) = (coord(1:3,jx,j) - temporary(1:3,jx,j))*inv_dt
            vel(1:3,jx,j) = vel(1:3,jx,j) + delta_vel(1:3)
          end do
        end do
      end if

      kin_ref(1:3) = 0.0_dp
      do j = 1, ncell
        do jx = 1, natom(j)

          ! FEP
          ! skip singleB to avoid duplication
          if (domain%fep_use) then
            if (domain%fepgrp(jx,j) == 2) cycle
          end if

          kin_ref(1:3) = kin_ref(1:3) + mass(jx,j)*vel(1:3,jx,j)*vel(1:3,jx,j)
        end do
      end do
      call mpi_allreduce(mpi_in_place, kin_ref, 3, mpi_real8, mpi_sum, &
                         mpi_comm_country, ierror)


      do j = 1, ncell
        do jx = 1, natom(j)
          vel(1:3,jx,j) = vel_ref(1:3,jx,j)
        end do
      end do

      ! Thermostat
      !
      if (mod(istep-1,dynamics%thermo_period) == 0) then
        do j = 1, ncell
          do jx = 1, natom(j)
            vel(1:3,jx,j) = vel(1:3,jx,j) * vel_scale
            vel(1:3,jx,j) = vel(1:3,jx,j) + random_f(1:3,jx,j)
          end do
        end do
      end if

      ! FEP
      ! synchronize single B with single A
      if (domain%fep_use) call sync_single_fep(domain, vel)

      ! VV1 (long range force)
      !
      if (mod(istep-1,dynamics%elec_long_period) == 0) then
        do j = 1, ncell
          do jx = 1, nsolute(j)
            factor = 0.5_dp * dt_long / mass(jx,j)
            vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_long(1:3,jx,j)
          end do
          do l = 1, nwater(j)
            do k = 1, 3
              jx = water_list(k,l,j)
              factor = 0.5_dp * dt_long / mass(jx,j)
              vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_long(1:3,jx,j) 
            end do
          end do
        end do

      end if

      do j = 1, ncell
        do jx = 1, nsolute(j)
          factor = half_dt_short / mass(jx,j)
          vel  (1:3,jx,j) = vel(1:3,jx,j) + factor*force_short(1:3,jx,j)
          coord(1:3,jx,j) = coord_ref(1:3,jx,j) &
                          + vel(1:3,jx,j)*dt_short
        end do
        do l = 1, nwater(j)
          do k = 1, 3
            jx = water_list(k,l,j)
            factor = half_dt_short / mass(jx,j)
            vel  (1:3,jx,j) = vel(1:3,jx,j) + factor*force_short(1:3,jx,j)
            coord(1:3,jx,j) = coord_ref(1:3,jx,j) &
                            + vel(1:3,jx,j)*dt_short
          end do
        end do
      end do

      ! RATTLE VV1
      !
      if (constraints%rigid_bond) then
        do j = 1, ncell
          do jx = 1, natom(j)
            temporary(1:3,jx,j) = coord(1:3,jx,j)
          end do
        end do
        call compute_constraints(ConstraintModeLEAP, .false., dt_short, &
                                 coord_ref, domain, constraints, coord, &
                                 vel, viri_const)

        do j = 1, ncell
          do jx = 1, natom(j)
            delta_vel(1:3) = (coord(1:3,jx,j) - temporary(1:3,jx,j))*inv_dt
            vel(1:3,jx,j) = vel(1:3,jx,j) + delta_vel(1:3)
          end do
        end do
      end if

    end if

    if (mod(istep-1,dynamics%baro_period) == 0) then

      viri_const(1:3,1:3) = virial_constraint(1:3,1:3)

      ! compute box size(t+2dt)
      !   size(t+2dt) = exp[eta(t+3/2dt)*dt] * size(t+dt)
      !
      scale_b(1:3) = exp(bmoment(1:3)*dt_baro)
      boundary%box_size_x = scale_b(1) * boundary%box_size_x_ref
      boundary%box_size_y = scale_b(2) * boundary%box_size_y_ref
      boundary%box_size_z = scale_b(3) * boundary%box_size_z_ref

      call bcast_boxsize(boundary%box_size_x, boundary%box_size_y, &
                         boundary%box_size_z)

      boundary%cell_size_x = boundary%box_size_x / real(boundary%num_cells_x,dp)
      boundary%cell_size_y = boundary%box_size_y / real(boundary%num_cells_y,dp)
      boundary%cell_size_z = boundary%box_size_z / real(boundary%num_cells_z,dp)

      ! update boudary conditions
      !
      dynvars%barostat_momentum(1:3) = bmoment(1:3)
      do j = 1, ncell+nboundary
        do jx = 1, natom(j)
          domain%trans_vec(1:3,jx,j) = domain%trans_vec(1:3,jx,j) * scale_b(1:3)
        end do
      end do

      domain%system_size(1) = boundary%box_size_x
      domain%system_size(2) = boundary%box_size_y
      domain%system_size(3) = boundary%box_size_z
    
    end if

    return

  end subroutine langevin_barostat_vv1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    langevin_barostat_vv2
  !> @brief        Langevin thermostat and barostat
  !! @authors      JJ, TM
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] boundary    : boundary information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine langevin_barostat_vv2(dynamics, istep, inner_step, dt_long,    &
                                   dt_short, ensemble, domain, constraints, &
                                   boundary, dynvars)

    ! formal arguments
    type(s_dynamics),         intent(in)    :: dynamics
    integer,                  intent(in)    :: istep
    integer,                  intent(in)    :: inner_step
    real(dp),                 intent(in)    :: dt_long
    real(dp),                 intent(in)    :: dt_short
    type(s_ensemble), target, intent(inout) :: ensemble
    type(s_domain),   target, intent(inout) :: domain
    type(s_constraints),      intent(inout) :: constraints
    type(s_boundary),         intent(inout) :: boundary
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    real(dp)                 :: inv_dt, temp0, press0, d_ndegf
    real(dp)                 :: kin(1:3), ekin
    real(dp)                 :: velx, vely, velz
    real(dp)                 :: volume, press(1:3)
    real(dp)                 :: pressxy, pressxyz
    real(dp)                 :: factor
    real(dp)                 :: bmoment_ref(3)
    real(dp)                 :: gamma_t, gamma_p
    real(dp)                 :: sigma, virial_constraint(1:3,1:3)
    real(dp)                 :: virial_sum(1:3)
    real(dp)                 :: v1, v2, rsq, grandom(1:3)
    real(dp)                 :: half_dt, quart_dt, half_dt_baro, quart_dt_baro
    real(dp)                 :: half_dt_short
    real(dp)                 :: dt_therm, half_dt_therm
    real(dp)                 :: size_scale(1:3), vel_scale
    real(dp)                 :: vel_change(1:3)
    integer                  :: i, j, jx, k, l, i_ndegf, ncell, maxiter

    real(dp),        pointer :: pmass, pforce(:), random_f(:,:,:)
    real(dp),        pointer :: mass(:,:)
    real(dp),        pointer :: vel(:,:,:), vel_ref(:,:,:), force_add(:,:,:)
    real(dp),        pointer :: force_long(:,:,:), force_short(:,:,:)
    real(dp),        pointer :: coord(:,:,:), coord_ref(:,:,:)
    real(dp),        pointer :: coord_deri(:,:,:)
    real(dp),        pointer :: virial_short(:,:), virial_long(:,:)
    real(dp),        pointer :: viri_const(:,:)
    real(dp),        pointer :: bmoment(:)
    integer,         pointer :: natom(:), nsolute(:), nwater(:)
    integer,         pointer :: water_list(:,:,:)


    inv_dt       =  1.0_dp/dt_short
    temp0        =  ensemble%temperature
    press0       =  ensemble%pressure * ATMOS_P
    gamma_t      =  ensemble%gamma_t * AKMA_PS
    gamma_p      =  ensemble%gamma_p * AKMA_PS
    pmass        => ensemble%pmass
    pforce       => ensemble%pforce
    random_f     => ensemble%random_force
    i_ndegf      =  domain%num_deg_freedom
    d_ndegf      =  real(i_ndegf,dp)
    ncell        =  domain%num_cell_local

    mass         => domain%mass
    natom        => domain%num_atom
    nsolute      => domain%num_solute
    nwater       => domain%num_water 
    water_list   => domain%water_list
    coord        => domain%coord
    coord_ref    => domain%coord_ref
    coord_deri   => domain%velocity_full
    vel          => domain%velocity
    vel_ref      => domain%velocity_ref
    force_long   => domain%force_long
    force_short  => domain%force_short
    force_add    => domain%coord_old
    virial_short => dynvars%virial
    virial_long  => dynvars%virial_long
    viri_const   => dynvars%virial_const
    bmoment      => dynvars%barostat_momentum

    ! save box size(t+dt) and compute volume(t+dt)
    !
    boundary%box_size_x_ref = boundary%box_size_x
    boundary%box_size_y_ref = boundary%box_size_y
    boundary%box_size_z_ref = boundary%box_size_z
    volume = boundary%box_size_x * boundary%box_size_y * boundary%box_size_z

    ! time step
    !
    half_dt_short = dt_short / 2.0_dp
    dt_therm      = dt_short * real(dynamics%thermo_period, dp)
    half_dt_therm = dt_therm / 2.0_dp
    half_dt_baro  = dt_short * real(dynamics%baro_period, dp) / 2.0_dp
    quart_dt_baro = half_dt_baro / 2.0_dp

    ! initial constraint virial
    !
    virial_constraint(1:3,1:3) = 0.0_dp

    ! initial constraint force
    !
    do j = 1, ncell
      do jx = 1, natom(j)
        force_add(1:3,jx,j) = 0.0_dp
      end do
    end do

    ! Langevin piston force
    !
    if (mod(istep, dynamics%baro_period) == 0) then

      ! barostate coefficient
      !
      bmoment_ref(1:3) = bmoment(1:3)

      ! pmass and stochastic force (pforce)
      !
      if (ensemble%isotropy == IsotropyISO) then

        sigma = sqrt(2.0_dp*gamma_p*pmass*KBOLTZ*temp0/half_dt_baro)
        pforce(1) = sigma * random_get_gauss()
        pforce(2) = pforce(1)
        pforce(3) = pforce(1)
  
      else if (ensemble%isotropy == IsotropySEMI_ISO) then

        sigma = sqrt(gamma_p*pmass*KBOLTZ*temp0/half_dt_baro)
        pforce(1) = sigma * random_get_gauss()
        pforce(2) = pforce(1)
        pforce(3) = sigma * random_get_gauss()

      else if (ensemble%isotropy == IsotropyANISO) then

        sigma = sqrt(gamma_p*pmass*KBOLTZ*temp0/half_dt_baro)
        pforce(1) = sigma * random_get_gauss()
        pforce(2) = sigma * random_get_gauss()
        pforce(3) = sigma * random_get_gauss()
  
      else if (ensemble%isotropy == IsotropyXY_Fixed) then
  
        sigma = sqrt(gamma_p*pmass*KBOLTZ*temp0/half_dt_baro)
        pforce(1) = 0.0_dp
        pforce(2) = 0.0_dp
        pforce(3) = sigma * random_get_gauss()
  
      end if

#ifdef HAVE_MPI_GENESIS
      call mpi_bcast(pforce, 3, mpi_real8, 0, mpi_comm_country, ierror)
#endif
    end if

    ! Langevin random force
    !
    if (mod(istep, dynamics%thermo_period) == 0) then

      vel_scale = exp(-gamma_t*half_dt_therm)
      factor    = 1.0_dp - vel_scale*vel_scale
      factor    = factor*KBOLTZ*temp0/2.0_dp
      do j = 1, ncell
        do jx = 1, natom(j)

          if (abs(mass(jx,j)) < EPS) cycle

          sigma = sqrt(factor/mass(jx,j))
          rsq = 2.0_dp

          do while (rsq >= 1.0_dp)
            v1  = 2.0_dp*random_get() - 1.0_dp
            v2  = 2.0_dp*random_get() - 1.0_dp
            rsq = v1*v1 + v2*v2
          end do

          rsq   = sqrt(-2.0_dp * log(rsq) / rsq)
          grandom(1) = rsq * v1
          rsq = 2.0_dp
  
          do while (rsq >= 1.0_dp)
            v1  = 2.0_dp*random_get() - 1.0_dp
            v2  = 2.0_dp*random_get() - 1.0_dp
            rsq = v1*v1 + v2*v2
          end do
  
          rsq   = sqrt(-2.0_dp * log(rsq) / rsq)
          grandom(2) = rsq * v1
          rsq = 2.0_dp
  
          do while (rsq >= 1.0_dp)
            v1  = 2.0_dp*random_get() - 1.0_dp
            v2  = 2.0_dp*random_get() - 1.0_dp
            rsq = v1*v1 + v2*v2
          end do

          rsq   = sqrt(-2.0_dp * log(rsq) / rsq)
          grandom(3) = rsq * v1
          random_f(1:3,jx,j) = sigma*grandom(1:3)

        end do
      end do

    end if

    ! VV2 (short range force)
    !
    do j = 1, ncell
      do jx = 1, nsolute(j)
        factor = half_dt_short / mass(jx,j)
        vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_short(1:3,jx,j)
        vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_add(1:3,jx,j)
      end do
      do l = 1, nwater(j)
        do k = 1, 3
          jx = water_list(k,l,j)
          factor = half_dt_short / mass(jx,j)
          vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_short(1:3,jx,j)
          vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_add(1:3,jx,j)
        end do
      end do     
    end do

    ! VV2 (long range force)
    !
    if (inner_step == dynamics%elec_long_period) then
      do j = 1, ncell
        do jx = 1, nsolute(j)
          factor = 0.5_dp * dt_long / mass(jx,j)
          vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_long(1:3,jx,j)
        end do
        do l = 1, nwater(j)
          do k = 1, 3
            jx = water_list(k,l,j)
            factor = 0.5_dp * dt_long / mass(jx,j)
            vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_long(1:3,jx,j)
          end do
        end do
      end do
    end if

    ! Thermostat
    !
    if (mod(istep, dynamics%thermo_period) == 0) then
      do j = 1, ncell
        do jx = 1, natom(j)
          vel(1:3,jx,j) = vel(1:3,jx,j) * vel_scale
          vel(1:3,jx,j) = vel(1:3,jx,j) + random_f(1:3,jx,j)
        end do
      end do
    end if

    ! FEP
    ! synchronize single B with single A
    if (domain%fep_use) call sync_single_fep(domain, vel)

    ! RATTLE VV2
    !
    if (constraints%rigid_bond) then
      do j = 1, ncell
        do jx = 1, natom(j)
          coord_deri(1:3,jx,j) = vel(1:3,jx,j) &
                               + bmoment(1:3)*coord(1:3,jx,j)
          coord_ref(1:3,jx,j) = coord_deri(1:3,jx,j)
        end do
      end do
      call compute_constraints(ConstraintModeVVER2, .false., dt_short, &
                               coord_ref, domain, constraints, coord,  &
                               coord_deri, viri_const)

      if (mod(istep,dynamics%baro_period) == 0) then

        do j = 1, ncell
          do jx = 1, natom(j)
            vel_change(1:3) = coord_deri(1:3,jx,j) - coord_ref(1:3,jx,j)
            vel(1:3,jx,j) = vel(1:3,jx,j) + vel_change(1:3)
            force_add(1:3,jx,j) = force_add(1:3,jx,j)                &
                                + mass(jx,j)*vel_change(1:3)/half_dt_short
          end do
        end do

        ! constraint virial
        ! 
        do j = 1, ncell
          do jx = 1, natom(j)

            ! FEP
            ! skip singleB to avoid duplication
            if (domain%fep_use) then
              if (domain%fepgrp(jx,j) == 2) cycle
            end if

            virial_constraint(1,1) = virial_constraint(1,1)             + &
                                     force_add(1,jx,j)*coord(1,jx,j)
            virial_constraint(2,2) = virial_constraint(2,2)             + &
                                     force_add(2,jx,j)*coord(2,jx,j)
            virial_constraint(3,3) = virial_constraint(3,3)             + &
                                     force_add(3,jx,j)*coord(3,jx,j)
          end do
        end do
      end if

    end if

    if (mod(istep, dynamics%baro_period) == 0) &
      viri_const(1:3,1:3) = virial_constraint(1:3,1:3)

    return

  end subroutine langevin_barostat_vv2

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    bussi_thermostat
  !> @brief        control temperature using Bussi's stochastic re-scaling
  !! @authors      TA, JJ
  !! @param[in]    dynamics   : dynamics information
  !! @param[in]    dt_short   : short time step
  !! @param[in]    ensemble   : ensemble information
  !! @param[inout] domain     : domain information
  !! @param[inout] dynvars    : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine bussi_thermostat(dynamics, dt_short, ensemble, domain, dynvars)

    ! formal arguments
    type(s_dynamics),        intent(in)    :: dynamics
    real(dp),                intent(in)    :: dt_short
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_domain),  target, intent(inout) :: domain
    type(s_dynvars), target, intent(inout) :: dynvars

    ! local variables
    real(dp)                 :: temp0, tau_t, dt_therm
    real(dp)                 :: ekin
    real(dp)                 :: tempf, tempt, factor, rr, sigma
    real(dp)                 :: degree
    integer                  :: i, ix
    integer                  :: num_degree, ncell

    real(dp),        pointer :: vel(:,:,:), mass(:,:)
    integer,         pointer :: natom(:)


    ! use pointers
    !
    dt_therm    =  dt_short * real(dynamics%thermo_period, dp)
    temp0       =  ensemble%temperature
    tau_t       =  ensemble%tau_t/AKMA_PS
    num_degree  =  domain%num_deg_freedom
    degree      =  real(num_degree, dp)
    ncell       =  domain%num_cell_local
    natom       => domain%num_atom
    vel         => domain%velocity
    mass        => domain%mass

    ! calculate temperature
    !
    ekin = 0.0_dp
    do i = 1, ncell
      do ix = 1, natom(i)

        ! FEP
        ! skip singleB to avoid duplication
        if (domain%fep_use) then
          if (domain%fepgrp(ix,i) == 2) cycle
        end if

        ekin = ekin + &
             mass(ix,i)*(vel(1,ix,i)*vel(1,ix,i) + vel(2,ix,i)*vel(2,ix,i) + &
                         vel(3,ix,i)*vel(3,ix,i))
      end do
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, ekin, 1, mpi_real8, mpi_sum, &
                       mpi_comm_country, ierror)
#endif
    tempf = ekin/(degree*KBOLTZ)

    ! calculate scaling factor
    !
    factor = exp(-dt_therm/tau_t)
    rr = random_get_gauss()
    tempt = tempf*factor &
          + temp0/degree*(1.0_dp-factor)*(sum_gauss(num_degree-1)+rr*rr) &
          + 2.0_dp*sqrt(tempf*temp0/degree*(1.0_dp-factor)*factor)*rr

    factor = sqrt(tempt/tempf)

#ifdef HAVE_MPI_GENESIS
    call mpi_bcast(factor, 1, mpi_real8, 0, mpi_comm_country, ierror)
#endif

    ! scale velocities
    !
    do i = 1, ncell
      do ix = 1, natom(i)
        vel(1:3,ix,i) = factor*vel(1:3,ix,i)
      end do
    end do

    return

  end subroutine bussi_thermostat

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    bussi_barostat_vv1
  !> @brief        Bussi thermostat and barostat
  !! @authors      TA, JJ
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    ensemble    : ensemble information
  !! @param[in]    inner_step  : inner step
  !! @param[in]    dt_long     : long time step
  !! @param[in]    dt_short    : shot time step
  !! @param[inout] domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] boundary    : boundary information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine bussi_barostat_vv1(dynamics, istep, inner_step, dt_long,    &
                                dt_short, ensemble, domain, constraints, &
                                boundary, dynvars)

    ! formal arguments
    type(s_dynamics),         intent(in)    :: dynamics
    integer,                  intent(in)    :: istep
    integer,                  intent(in)    :: inner_step
    real(dp),                 intent(in)    :: dt_long
    real(dp),                 intent(in)    :: dt_short
    type(s_ensemble), target, intent(inout) :: ensemble
    type(s_domain),   target, intent(inout) :: domain
    type(s_constraints),      intent(inout) :: constraints
    type(s_boundary),         intent(inout) :: boundary
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    real(dp)                 :: dt, inv_dt, temp0, press0, d_ndegf, degree
    real(dp)                 :: dt_therm, half_dt_therm
    real(dp)                 :: dt_baro, half_dt_baro
    real(dp)                 :: kin(1:3), kin_temp(1:3), ekin
    real(dp)                 :: vel_tmp(1:3)
    real(dp)                 :: volume, press(1:3)
    real(dp)                 :: pressxy, pressxyz
    real(dp)                 :: factor
    real(dp)                 :: bmoment_ref(3), scale_b(1:3), vel_scale_2(1:3)
    real(dp)                 :: tau_t, tau_p, tempt, tempf, sigma
    real(dp)                 :: gr, ekin0, vel_scale
    real(dp)                 :: half_dt, half_dt_long, size_scale(1:3)
    real(dp)                 :: virial_constraint(3,3), virial_sum(3)
    integer                  :: i, j, ij, jx, k, l, i_ndegf, maxiter
    integer                  :: ncell, nboundary

    real(dp),        pointer :: pmass 
    real(dp),        pointer :: mass(:,:)
    real(dp),        pointer :: coord(:,:,:), coord_ref(:,:,:)
    real(dp),        pointer :: temporary(:,:,:), temporary1(:,:,:)
    real(dp),        pointer :: vel(:,:,:), vel_ref(:,:,:), kin_ref(:)
    real(dp),        pointer :: force_short(:,:,:), force_long(:,:,:)
    real(dp),        pointer :: virial(:,:), viri_const(:,:), virial_long(:,:)
    real(dp),        pointer :: bmoment(:)
    integer,         pointer :: natom(:), nsolute(:), nwater(:)
    integer,         pointer :: water_list(:,:,:)


    dt          =  dt_short
    inv_dt      =  1.0_dp/dt
    temp0       =  ensemble%temperature
    press0      =  ensemble%pressure * ATMOS_P
    tau_t       =  ensemble%tau_t / AKMA_PS
    tau_p       =  ensemble%tau_p / AKMA_PS
    pmass       => ensemble%pmass
    i_ndegf     =  domain%num_deg_freedom
    d_ndegf     =  real(i_ndegf,dp)
    ncell       =  domain%num_cell_local
    nboundary   =  domain%num_cell_boundary
    natom       => domain%num_atom
    nsolute     => domain%num_solute
    nwater      => domain%num_water 
    water_list  => domain%water_list
    mass        => domain%mass
    coord       => domain%coord
    coord_ref   => domain%coord_ref
    vel         => domain%velocity
    vel_ref     => domain%velocity_ref
    force_long  => domain%force_long
    force_short => domain%force_short
    temporary   => domain%coord_old
    temporary1  => domain%velocity_full
    virial      => dynvars%virial
    viri_const  => dynvars%virial_const
    virial_long => dynvars%virial_long
    bmoment     => dynvars%barostat_momentum
    kin_ref     => dynvars%kinetic_ref

    ! save box size(t+dt) and compute volume(t+dt)
    !
    boundary%box_size_x_ref = boundary%box_size_x
    boundary%box_size_y_ref = boundary%box_size_y
    boundary%box_size_z_ref = boundary%box_size_z
    volume = boundary%box_size_x * boundary%box_size_y * boundary%box_size_z

    ! time step
    !
    half_dt  = dt / 2.0_dp
    half_dt_long = dt_long / 2.0_dp
    dt_therm      = dt_short * real(dynamics%thermo_period, dp)
    half_dt_therm = dt_therm / 2.0_dp
    dt_baro       = dt_short * real(dynamics%baro_period, dp)
    half_dt_baro  = dt_baro / 2.0_dp

    ! maximum iteration
    !
    if (mod(istep-1,dynamics%baro_period) == 0) then
      if (constraints%rigid_bond) then
        maxiter = 4
      else
        maxiter = 1
      end if
    else
      maxiter = 1
    end if

    ! constraint virial
    !
    virial_constraint(1:3,1:3) = 0.0_dp

    ! barostate coefficient
    !
    bmoment_ref(1:3) = bmoment(1:3)

    ! pmass and ekin0
    !
    if (ensemble%isotropy == IsotropyISO) then
      degree = d_ndegf+3.0_dp
      pmass = degree*KBOLTZ*temp0 * tau_p**2
      ekin0 = 0.5_dp*KBOLTZ*temp0 * degree
    else if (ensemble%isotropy == IsotropySEMI_ISO) then
      degree = d_ndegf+3.0_dp
      pmass = degree*KBOLTZ*temp0 * tau_p**2 / 3.0_dp
      ekin0 = 0.5_dp*KBOLTZ*temp0 * degree
    else if (ensemble%isotropy == IsotropyANISO) then
      degree = d_ndegf+3.0_dp
      pmass = degree*KBOLTZ*temp0 * tau_p**2 / 3.0_dp
      ekin0 = 0.5_dp*KBOLTZ*temp0 * degree
    else if (ensemble%isotropy == IsotropyXY_Fixed) then
      degree = d_ndegf+1.0_dp
      pmass = degree*KBOLTZ*temp0 * tau_p**2
      ekin0 = 0.5_dp*KBOLTZ*temp0 * degree
    end if

    vel_scale = 1.0_dp

    if (istep == 1) then
      kin_ref(1:3) = 0.0_dp
      do j = 1, ncell
        do jx = 1, natom(j)

          ! FEP
          ! skip singleB to avoid duplication
          if (domain%fep_use) then
            if (domain%fepgrp(jx,j) == 2) cycle
          end if

          kin_ref(1:3) = kin_ref(1:3)  &
                       + mass(jx,j)*vel(1:3,jx,j)*vel(1:3,jx,j)
        end do
      end do
      call mpi_allreduce(mpi_in_place, kin_ref, 3, mpi_real8, mpi_sum, &
                         mpi_comm_country, ierror)
    end if

    ! thermostat
    !
    if (mod(istep-1,dynamics%thermo_period) == 0) then

      kin(1:3) = 0.0_dp
      do j = 1, ncell
        do jx = 1, natom(j)

          ! FEP
          ! skip singleB to avoid duplication
          if (domain%fep_use) then
            if (domain%fepgrp(jx,j) == 2) cycle
          end if

          kin(1:3)  = kin(1:3) + mass(jx,j)*vel_ref(1:3,jx,j)*vel_ref(1:3,jx,j)
        end do
      end do
      ekin   = 0.5_dp * (kin(1)+kin(2)+kin(3))

#ifdef HAVE_MPI_GENESIS
      call mpi_allreduce(mpi_in_place, ekin, 1, mpi_real8, mpi_sum, &
                         mpi_comm_country, ierror)
#endif

      ! compute stochastic scaling factor
      !
      factor = exp(-dt_therm/tau_t)
      ekin  = ekin + 0.5_dp*pmass*dot_product(bmoment_ref(1:3),bmoment_ref(1:3))
      gr = random_get_gauss()
      vel_scale = factor &
                + (1.0_dp-factor)*ekin0/(degree*ekin) &
                 *(sum_gauss(i_ndegf+2)+gr**2)        &
                + 2.0_dp*gr*sqrt(ekin0/(degree*ekin)*factor*(1.0_dp-factor))
      vel_scale = sqrt(vel_scale)
      vel_scale = sign(vel_scale, &
                       gr+sqrt(factor*degree*ekin/((1.0_dp-factor)*ekin0)))

#ifdef HAVE_MPI_GENESIS
      call mpi_bcast(vel_scale, 1, mpi_real8, 0, mpi_comm_country, ierror)
#endif

    end if

    !reference
    !
    do j = 1, ncell
      do jx = 1, natom(j)
        vel_ref(1:3,jx,j) = vel(1:3,jx,j)
        temporary1(1:3,jx,j) = 0.0_dp
      end do
    end do

    ! decide barostat momentum
    !
    if (mod(istep-1,dynamics%baro_period) == 0) then

      do i = 1, maxiter

        ! scale bmoment (eta) for barostat
        !
        bmoment(1:3) = bmoment_ref(1:3)*vel_scale

        ! compute kinetic energy
        !
        kin_temp(1:3) = 0.0_dp
        do j = 1, ncell
          do jx = 1, natom(j)

            ! FEP
            ! skip singleB to avoid duplication
            if (domain%fep_use) then
              if (domain%fepgrp(jx,j) == 2) cycle
            end if

            kin_temp(1:3)  = kin_temp(1:3) &
                           + mass(jx,j) * vel(1:3,jx,j)*vel(1:3,jx,j)
          end do
        end do

        ! virial + virial_constraint
        !
        virial_sum(1) = virial(1,1) + virial_long(1,1) &
                      + virial_constraint(1,1)
        virial_sum(2) = virial(2,2) + virial_long(2,2) &
                      + virial_constraint(2,2)
        virial_sum(3) = virial(3,3) + virial_long(3,3) &
                      + virial_constraint(3,3)

        call reduce_pres(kin_temp, ekin, virial_sum)
        kin(1:3) = (kin_temp(1:3)+kin_ref(1:3)) / 2.0_dp
        ekin = 0.5_dp*(kin(1)+kin(2)+kin(3))

        ! compute pressure in the unit of kcal/mol*A3
        !
        press(1:3) = (kin(1:3) + virial_sum(1:3))/volume
        pressxyz = (press(1)+press(2)+press(3))/3.0_dp
        pressxy  = (press(1)+press(2))/2.0_dp

        ! update barostat
        !
        call update_barostat_mtk(ensemble, press(1), press(2), press(3), &
                                 pressxyz, pressxy, press0, volume,      &
                                 d_ndegf, pmass, dt_baro, ekin, bmoment)

        ! scale velocities
        !
        do j = 1, ncell
          do jx = 1, natom(j)
            vel(1:3,jx,j) = vel_ref(1:3,jx,j)*vel_scale
          end do
        end do

        size_scale(1:3)  = exp(bmoment(1:3)*dt)
        vel_scale_2(1:3) = powersinh(bmoment(1:3)*dt)*dt
        do j = 1, ncell
          do jx = 1, nsolute(j)
            factor = half_dt / mass(jx,j)
            vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_short(1:3,jx,j)
            vel(1:3,jx,j) = vel(1:3,jx,j) + factor*temporary1(1:3,jx,j)
            vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_long(1:3,jx,j)
            coord(1:3,jx,j) = size_scale(1:3)*coord_ref(1:3,jx,j) &
                            + vel_scale_2*vel(1:3,jx,j)
            vel  (1:3,jx,j) = vel(1:3,jx,j)/size_scale(1:3)
          end do
          do l = 1, nwater(j)
            do k = 1, 3
              jx = water_list(k,l,j)
              factor = half_dt / mass(jx,j)
              vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_short(1:3,jx,j)
              vel(1:3,jx,j) = vel(1:3,jx,j) + factor*temporary1(1:3,jx,j)
              vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_long(1:3,jx,j)
              coord(1:3,jx,j) = size_scale(1:3)*coord_ref(1:3,jx,j) &
                              + vel_scale_2*vel(1:3,jx,j)
              vel  (1:3,jx,j) = vel(1:3,jx,j)/size_scale(1:3)
            end do
          end do
        end do

        ! RATTLE VV1
        !
        if (constraints%rigid_bond) then

          do j = 1, ncell
            do jx = 1, natom(j)
              temporary(1:3,jx,j) = coord(1:3,jx,j)
            end do
          end do

          call compute_constraints(ConstraintModeLEAP, .false., dt, coord_ref, &
                                   domain, constraints, coord, vel, &
                                   viri_const)

          virial_constraint(1:3,1:3) = virial_constraint(1:3,1:3) &
                                     + 2.0_dp * viri_const(1:3,1:3)

          do j = 1, ncell
            do jx = 1, natom(j)
              vel_tmp(1:3) = (coord(1:3,jx,j)-temporary(1:3,jx,j)) * inv_dt
              vel(1:3,jx,j) = vel(1:3,jx,j) + vel_tmp(1:3)
              temporary1(1:3,jx,j) = temporary1(1:3,jx,j) &
                                   + mass(jx,j)*vel_tmp(1:3)/half_dt
            end do
          end do
        end if

      end do

      ! thermostat
      !
      do j = 1, ncell
        do jx = 1, natom(j)
          vel(1:3,jx,j) = vel_ref(1:3,jx,j)*vel_scale
        end do
      end do

      ! VV1 (long range)
      !
      do j = 1, ncell
        do jx = 1, nsolute(j)
          factor = 0.5_dp*dt_long / mass(jx,j)
          vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_long(1:3,jx,j)
        end do
        do l = 1, nwater(j)
          do k = 1, 3
            jx = water_list(k,l,j)
            factor = 0.5_dp*dt_long / mass(jx,j)
            vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_long(1:3,jx,j)
          end do
        end do 
      end do

      ! VV1 (short range)
      !
      do j = 1, ncell
        do jx = 1, nsolute(j)
          factor = half_dt / mass(jx,j)
          vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_short(1:3,jx,j)
          coord(1:3,jx,j) = size_scale(1:3)*coord_ref(1:3,jx,j) &
                          + vel_scale_2*vel(1:3,jx,j)
          vel  (1:3,jx,j) = vel(1:3,jx,j)/size_scale(1:3)
        end do
        do l = 1, nwater(j)
          do k = 1, 3
            jx = water_list(k,l,j)
            factor = half_dt / mass(jx,j)
            vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_short(1:3,jx,j)
            coord(1:3,jx,j) = size_scale(1:3)*coord_ref(1:3,jx,j) &
                            + vel_scale_2*vel(1:3,jx,j)
            vel  (1:3,jx,j) = vel(1:3,jx,j)/size_scale(1:3)
          end do
        end do
      end do

      ! RATTLE VV1
      !
      if (constraints%rigid_bond) then
        do j = 1, ncell
          do jx = 1, natom(j)
            temporary(1:3,jx,j) = coord(1:3,jx,j)
          end do
        end do
        call compute_constraints(ConstraintModeLEAP, .false., dt_short, &
                                 coord_ref, domain, constraints, coord, &
                                 vel, viri_const)

        do j = 1, ncell
          do jx = 1, natom(j)
            vel_tmp(1:3) = (coord(1:3,jx,j) - temporary(1:3,jx,j))*inv_dt
            vel(1:3,jx,j) = vel(1:3,jx,j) + vel_tmp(1:3)
          end do
        end do
      end if

    else

      ! kinetic energy component calculation
      !
      do j = 1, ncell
        do jx = 1, natom(j)
          vel(1:3,jx,j) = vel_ref(1:3,jx,j)
        end do
      end do

      if (mod(istep-1,dynamics%thermo_period) == 0) then
        do j = 1, ncell
          do jx = 1, natom(j)
            vel(1:3,jx,j) = vel(1:3,jx,j) * vel_scale
          end do
        end do
      end if

      size_scale(1:3)  = exp(bmoment(1:3)*dt)
      vel_scale_2(1:3) = powersinh(bmoment(1:3)*dt)*dt

      do j = 1, ncell
        do jx = 1, nsolute(j)
          factor = half_dt / mass(jx,j)
          vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_long(1:3,jx,j)
          vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_short(1:3,jx,j)
          coord(1:3,jx,j) = size_scale(1:3)*coord_ref(1:3,jx,j) &
                          + vel_scale_2*vel(1:3,jx,j)
          vel  (1:3,jx,j) = vel(1:3,jx,j)/size_scale(1:3)
        end do
        do l = 1, nwater(j)
          do k = 1, 3
            jx = water_list(k,l,j)
            factor = half_dt / mass(jx,j)
            vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_long(1:3,jx,j)
            vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_short(1:3,jx,j)
            coord(1:3,jx,j) = size_scale(1:3)*coord_ref(1:3,jx,j) &
                            + vel_scale_2*vel(1:3,jx,j)
            vel  (1:3,jx,j) = vel(1:3,jx,j)/size_scale(1:3)
          end do
        end do
      end do

      ! RATTLE VV1
      !
      if (constraints%rigid_bond) then
        do j = 1, ncell
          do jx = 1, natom(j)
            temporary(1:3,jx,j) = coord(1:3,jx,j)
          end do
        end do
        call compute_constraints(ConstraintModeLEAP, .false., dt_short, &
                                 coord_ref, domain, constraints, coord, &
                                 vel, viri_const)
  
        do j = 1, ncell
          do jx = 1, natom(j)
            vel_tmp(1:3) = (coord(1:3,jx,j) - temporary(1:3,jx,j))*inv_dt
            vel(1:3,jx,j) = vel(1:3,jx,j) + vel_tmp(1:3)
          end do
        end do
      end if
 
      kin_ref(1:3) = 0.0_dp
      do j = 1, ncell
        do jx = 1, natom(j)

          ! FEP
          ! skip singleB to avoid duplication
          if (domain%fep_use) then
            if (domain%fepgrp(jx,j) == 2) cycle
          end if

          kin_ref(1:3) = kin_ref(1:3) + mass(jx,j)*vel(1:3,jx,j)*vel(1:3,jx,j)
        end do
      end do
      call mpi_allreduce(mpi_in_place, kin_ref, 3, mpi_real8, mpi_sum, &
                         mpi_comm_country, ierror)

      ! from here we start VV1
      !
      do j = 1, ncell
        do jx = 1, natom(j)
          vel(1:3,jx,j) = vel_ref(1:3,jx,j)
        end do
      end do

      ! thermostat
      !
      if (mod(istep-1,dynamics%thermo_period) == 0) then
        do j = 1, ncell
          do jx = 1, natom(j)
            vel(1:3,jx,j) = vel(1:3,jx,j) * vel_scale
          end do
        end do
      end if

      ! VV1 (long range force)
      !
      if (mod(istep-1,dynamics%elec_long_period) == 0) then
        do j = 1, ncell
          do jx = 1, nsolute(j)
            factor = 0.5_dp * dt_long / mass(jx,j)
            vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_long(1:3,jx,j)
          end do
          do l = 1, nwater(j)
            do k = 1, 3
              jx = water_list(k,l,j)
              factor = 0.5_dp * dt_long / mass(jx,j)
              vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_long(1:3,jx,j)
            end do
          end do
        end do

      end if

      ! VV1
      !
      do j = 1, ncell         
        do jx = 1, nsolute(j)
          factor = half_dt / mass(jx,j)
          vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_short(1:3,jx,j)
          coord(1:3,jx,j) = size_scale(1:3)*coord_ref(1:3,jx,j) &
                          + vel_scale_2*vel(1:3,jx,j)
          vel  (1:3,jx,j) = vel(1:3,jx,j)/size_scale(1:3)
        end do
        do l = 1, nwater(j)
          do k = 1, 3
            jx = water_list(k,l,j)
            factor = half_dt / mass(jx,j)
            vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_short(1:3,jx,j)
            coord(1:3,jx,j) = size_scale(1:3)*coord_ref(1:3,jx,j) &
                            + vel_scale_2*vel(1:3,jx,j)
            vel  (1:3,jx,j) = vel(1:3,jx,j)/size_scale(1:3)
          end do
        end do
      end do

      ! RATTLE VV1
      !
      if (constraints%rigid_bond) then
        do j = 1, ncell
          do jx = 1, natom(j)
            temporary(1:3,jx,j) = coord(1:3,jx,j)
          end do
        end do
        call compute_constraints(ConstraintModeLEAP, .false., dt_short, &
                                 coord_ref, domain, constraints, coord, &
                                 vel, viri_const)

        do j = 1, ncell
          do jx = 1, natom(j)
            vel_tmp(1:3) = (coord(1:3,jx,j) - temporary(1:3,jx,j))*inv_dt
            vel(1:3,jx,j) = vel(1:3,jx,j) + vel_tmp(1:3)
          end do
        end do
      end if

    end if

    ! update virial constraint
    !
    viri_const(1:3,1:3) = virial_constraint(1:3,1:3)

    ! update barostat momentum
    !
    dynvars%barostat_momentum(1:3) = bmoment(1:3)

    ! compute box size
    !
    scale_b(1:3) = exp(bmoment(1:3)*dt)
    boundary%box_size_x = scale_b(1) * boundary%box_size_x_ref
    boundary%box_size_y = scale_b(2) * boundary%box_size_y_ref
    boundary%box_size_z = scale_b(3) * boundary%box_size_z_ref

    call bcast_boxsize(boundary%box_size_x, boundary%box_size_y, &
                       boundary%box_size_z)

    boundary%cell_size_x = boundary%box_size_x / real(boundary%num_cells_x,dp)
    boundary%cell_size_y = boundary%box_size_y / real(boundary%num_cells_y,dp)
    boundary%cell_size_z = boundary%box_size_z / real(boundary%num_cells_z,dp)

    ! update boudary conditions
    !
    do j = 1, ncell+nboundary
      do jx = 1, natom(j)
        domain%trans_vec(1:3,jx,j) = domain%trans_vec(1:3,jx,j) * scale_b(1:3)
      end do
    end do

    domain%system_size(1) = boundary%box_size_x
    domain%system_size(2) = boundary%box_size_y
    domain%system_size(3) = boundary%box_size_z

    return

  end subroutine bussi_barostat_vv1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    bussi_barostat_vv2
  !> @brief        Bussi thermostat and barostat
  !! @authors      TA, JJ
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    ensemble    : ensemble information
  !! @param[in]    istep       : step number
  !! @param[in]    inner_step  : inner step
  !! @param[in]    dt_long     : long time step
  !! @param[in]    dt_short    : shot time step
  !! @param[inout] domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] boundary    : boundary information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine bussi_barostat_vv2(dynamics, istep, inner_step, dt_long,    &
                                dt_short, ensemble, domain, constraints, &
                                boundary, dynvars)

    ! formal arguments
    type(s_dynamics),         intent(in)    :: dynamics
    integer,                  intent(in)    :: istep
    integer,                  intent(in)    :: inner_step
    real(dp),                 intent(in)    :: dt_long
    real(dp),                 intent(in)    :: dt_short
    type(s_ensemble), target, intent(inout) :: ensemble
    type(s_domain),   target, intent(inout) :: domain
    type(s_constraints),      intent(inout) :: constraints
    type(s_boundary),         intent(inout) :: boundary
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    real(dp)                 :: dt, inv_dt, temp0, press0, d_ndegf, degree
    real(dp)                 :: dt_therm, half_dt_therm
    real(dp)                 :: dt_baro, half_dt_baro
    real(dp)                 :: kin(1:3), ekin
    real(dp)                 :: vel_tmp(1:3), cm(1:8)
    real(dp)                 :: volume, press(1:3)
    real(dp)                 :: pressxy, pressxyz
    real(dp)                 :: factor
    real(dp)                 :: bmoment_ref(3)
    real(dp)                 :: tau_t, tau_p, tempt, tempf, sigma
    real(dp)                 :: gr, ekin0
    real(dp)                 :: virial_constraint(1:3,1:3)
    real(dp)                 :: virial_sum(1:3)
    real(dp)                 :: half_dt, half_dt_long, size_scale(1:3), vel_scale
    real(dp)                 :: vel_change(1:3), vel_scal_2(1:3)
    integer                  :: i, j, jx, k, l, i_ndegf, ncell, maxiter

    real(dp),        pointer :: pmass, pforce(:)
    real(dp),        pointer :: mass(:,:)
    real(dp),        pointer :: vel(:,:,:), vel_ref(:,:,:)
    real(dp),        pointer :: force_add(:,:,:), coord_deri(:,:,:)
    real(dp),        pointer :: force_short(:,:,:), force_long(:,:,:)
    real(dp),        pointer :: coord(:,:,:), coord_ref(:,:,:)
    real(dp),        pointer :: virial(:,:), viri_const(:,:), virial_long(:,:)
    real(dp),        pointer :: bmoment(:)
    integer,         pointer :: natom(:), nsolute(:), nwater(:)
    integer,         pointer :: water_list(:,:,:)

    dt          =  dt_short
    inv_dt      =  1.0_dp/dt
    temp0       =  ensemble%temperature
    press0      =  ensemble%pressure * ATMOS_P
    tau_t       =  ensemble%tau_t / AKMA_PS
    tau_p       =  ensemble%tau_p / AKMA_PS
    pmass       => ensemble%pmass
    pforce      => ensemble%pforce
    i_ndegf     =  domain%num_deg_freedom
    d_ndegf     =  real(i_ndegf,dp)
    ncell       =  domain%num_cell_local
    mass        => domain%mass
    natom       => domain%num_atom
    nsolute     => domain%num_solute
    nwater      => domain%num_water 
    water_list  => domain%water_list
    coord       => domain%coord
    coord_ref   => domain%coord_ref
    vel         => domain%velocity
    vel_ref     => domain%velocity_ref
    force_short => domain%force_short
    force_long  => domain%force_long
    force_add   => domain%coord_old
    coord_deri  => domain%velocity_full
    virial      => dynvars%virial
    viri_const  => dynvars%virial_const
    virial_long => dynvars%virial_long
    bmoment     => dynvars%barostat_momentum

    ! save box size(t+dt) and compute volume(t+dt)
    !
    boundary%box_size_x_ref = boundary%box_size_x
    boundary%box_size_y_ref = boundary%box_size_y
    boundary%box_size_z_ref = boundary%box_size_z
    volume = boundary%box_size_x * boundary%box_size_y * boundary%box_size_z

    ! time step
    !
    half_dt = dt / 2.0_dp
    half_dt_long = dt_long / 2.0_dp
    dt_therm      = dt_short * real(dynamics%thermo_period, dp)
    half_dt_therm = dt_therm / 2.0_dp
    dt_baro       = dt_short * real(dynamics%baro_period, dp)
    half_dt_baro  = dt_baro / 2.0_dp

    ! initial constraint virial
    !
    virial_constraint(1:3,1:3) = 0.0_dp

    ! initial constraint force
    !
    do j = 1, ncell
      do jx = 1, natom(j)
        force_add(1:3,jx,j) = 0.0_dp
      end do
    end do

    ! shift velocity, force
    !
    cm(1:8) = 0.0_dp
    do j = 1, ncell
      do jx = 1, nsolute(j)

        ! FEP
        ! skip singleB to avoid duplication
        if (domain%fep_use) then
          if (domain%fepgrp(jx,j) == 2) cycle
        end if

        cm(1:3)  = cm(1:3) + mass(jx,j)*vel(1:3,jx,j)
        cm(4)    = cm(4)   + mass(jx,j)
        cm(5:7)  = cm(5:7) + force_long(1:3,jx,j)
        cm(8)    = cm(8)   + 1.0_dp
      end do
      do l = 1, nwater(j)
        do k = 1, 3
          jx = water_list(k,l,j)
          cm(1:3)  = cm(1:3) + mass(jx,j)*vel(1:3,jx,j)
          cm(4)    = cm(4)   + mass(jx,j)
          cm(5:7)  = cm(5:7) + force_long(1:3,jx,j)
          cm(8)    = cm(8)   + 1.0_dp
        end do
      end do
    end do
    call mpi_allreduce(mpi_in_place, cm, 8, mpi_real8, mpi_sum, &
                       mpi_comm_city, ierror)
    do j = 1, ncell
      do jx = 1, natom(j)
        vel(1:3,jx,j)        = vel(1:3,jx,j)        - cm(1:3)/cm(4)
        force_long(1:3,jx,j) = force_long(1:3,jx,j) - cm(5:7)/cm(8)
      end do
    end do

    ! pmass and ekin0
    !
    if (ensemble%isotropy == IsotropyISO) then
      degree = d_ndegf+3.0_dp
      pmass = degree*KBOLTZ*temp0 * tau_p**2
      ekin0 = 0.5_dp*KBOLTZ*temp0 * degree
    else if (ensemble%isotropy == IsotropySEMI_ISO) then
      degree = d_ndegf+3.0_dp
      pmass = degree*KBOLTZ*temp0 * tau_p**2 / 3.0_dp
      ekin0 = 0.5_dp*KBOLTZ*temp0 * degree
    else if (ensemble%isotropy == IsotropyANISO) then
      degree = d_ndegf+3.0_dp
      pmass = degree*KBOLTZ*temp0 * tau_p**2 / 3.0_dp
      ekin0 = 0.5_dp*KBOLTZ*temp0 * real(i_ndegf+3,dp)
    else if (ensemble%isotropy == IsotropyXY_Fixed) then
      degree = d_ndegf+1.0_dp
      pmass = degree*KBOLTZ*temp0 * tau_p**2
      ekin0 = 0.5_dp*KBOLTZ*temp0 * degree
    end if

    ! VV2 (long range force)
    !
    if (mod(istep, dynamics%elec_long_period) == 0 ) then
      do j = 1, ncell
        do jx = 1, nsolute(j)
          factor = half_dt_long / mass(jx,j)
          vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_long(1:3,jx,j)
        end do
        do l = 1, nwater(j)
          do k = 1, 3
            jx = water_list(k,l,j)
            factor = half_dt_long / mass(jx,j)
            vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_long(1:3,jx,j)
          end do
        end do
      end do
    end if

    ! VV2 (short range force)
    !
    do j = 1, ncell
      do jx = 1, nsolute(j)
        factor = half_dt / mass(jx,j)
        vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_short(1:3,jx,j)
      end do
      do l = 1, nwater(j)
        do k = 1, 3
          jx = water_list(k,l,j)
          factor = half_dt / mass(jx,j)
          vel(1:3,jx,j) = vel(1:3,jx,j) + factor*force_short(1:3,jx,j)
        end do
      end do
    end do

    ! RATTLE VV2
    !
    if (constraints%rigid_bond) then
      do j = 1, ncell
        do jx = 1, natom(j)
          coord_deri(1:3,jx,j) = vel(1:3,jx,j) &
                               + bmoment(1:3)*coord(1:3,jx,j)
          coord_ref(1:3,jx,j)  = coord_deri(1:3,jx,j)
        end do
      end do
      call compute_constraints(ConstraintModeVVER2, .false., dt, coord_ref, &
                               domain, constraints, coord, coord_deri,      &
                               viri_const)
      do j = 1, ncell
        do jx = 1, natom(j)
          vel_change(1:3) = coord_deri(1:3,jx,j) - coord_ref(1:3,jx,j)
          vel(1:3,jx,j) = vel(1:3,jx,j) + vel_change(1:3)
          force_add(1:3,jx,j) = mass(jx,j)*vel_change(1:3)/half_dt
        end do
      end do

      ! constraint virial
      !
      virial_constraint(1:3,1:3) = 0.0_dp
      do j = 1, ncell
        do jx = 1, natom(j)
          virial_constraint(1,1) = virial_constraint(1,1)             + &
                                   force_add(1,jx,j)*coord(1,jx,j)
          virial_constraint(2,2) = virial_constraint(2,2)             + &
                                   force_add(2,jx,j)*coord(2,jx,j)
          virial_constraint(3,3) = virial_constraint(3,3)             + &
                                   force_add(3,jx,j)*coord(3,jx,j)
        end do
      end do

    end if

    viri_const(1:3,1:3) = virial_constraint(1:3,1:3)

    return

  end subroutine bussi_barostat_vv2

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_barostat
  !> @brief        update barostat parameter bmoment (eta)
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_barostat(ensemble, boundary, bmoment_ref, pforce, &
                             press, pressxyz, pressxy, press0, &
                             volume, d_ndegf, pmass, gamma_p,  &
                             ekin, half_dt, quart_dt, natom,   &
                             ncell, vel, bmoment)

    ! formal arguments
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_boundary),        intent(in)    :: boundary
    real(dp),                intent(in)    :: bmoment_ref(:)
    real(dp),                intent(in)    :: pforce(:)
    real(dp),                intent(in)    :: ekin
    real(dp),                intent(in)    :: press(:)
    real(dp),                intent(in)    :: pressxyz
    real(dp),                intent(in)    :: pressxy
    real(dp),                intent(in)    :: press0
    real(dp),                intent(in)    :: volume
    real(dp),                intent(in)    :: d_ndegf
    real(dp),                intent(in)    :: pmass
    real(dp),                intent(in)    :: gamma_p
    real(dp),                intent(in)    :: quart_dt
    real(dp),                intent(in)    :: half_dt
    integer,                 intent(in)    :: natom(:)
    integer,                 intent(in)    :: ncell
    real(dp),                intent(inout) :: vel(:,:,:)
    real(dp),                intent(inout) :: bmoment(:)

    ! local variable
    real(dp)         :: gamma0, pressxy0, trace
    real(dp)         :: vel_scale(1:3)
    integer          :: i, ix

    gamma0 = ensemble%gamma*ATMOS_P*100.0_dp/1.01325_dp

    ! eta(t+1/4dt)
    !
    bmoment(1:3) = exp(-gamma_p*quart_dt/2.0_dp)*bmoment_ref(1:3)

    ! eta(t+1/4dt) is scaled according to pressure
    !
    if (ensemble%isotropy == IsotropyISO) then

      bmoment(1) = bmoment(1) + quart_dt*(3.0_dp*volume*(pressxyz - press0) &
                              + 6.0_dp*ekin/d_ndegf + pforce(1))/pmass
      bmoment(2) = bmoment(1)
      bmoment(3) = bmoment(1)

    else if (ensemble%isotropy == IsotropySEMI_ISO) then

      if (ensemble%ensemble == EnsembleNPT) then
        bmoment(1) = bmoment(1) + quart_dt*(volume*(pressxy - press0)   &
                                + 2.0_dp*ekin/d_ndegf + pforce(1))/pmass
        bmoment(2) = bmoment(2) + quart_dt*(volume*(pressxy - press0)   &
                                + 2.0_dp*ekin/d_ndegf + pforce(2))/pmass
        bmoment(3) = bmoment(3) + quart_dt*(volume*(press(3)  - press0) &
                                + 2.0_dp*ekin/d_ndegf + pforce(3))/pmass
      else if (ensemble%ensemble == EnsembleNPgT) then
        pressxy0 = press0 - gamma0 / boundary%box_size_z
        bmoment(1) = bmoment(1) + quart_dt*(volume*(pressxy - pressxy0) &
                                + 2.0_dp*ekin/d_ndegf + pforce(1))/pmass
        bmoment(2) = bmoment(2) + quart_dt*(volume*(pressxy - pressxy0) &
                                + 2.0_dp*ekin/d_ndegf + pforce(2))/pmass
        bmoment(3) = bmoment(3) + quart_dt*(volume*(press(3)  - press0) &
                                + 2.0_dp*ekin/d_ndegf + pforce(3))/pmass
      end if

    else if (ensemble%isotropy == IsotropyANISO) then

      if (ensemble%ensemble == EnsembleNPT) then
        bmoment(1) = bmoment(1) + quart_dt*(volume*(press(1) - press0)   &
                                    + 2.0_dp*ekin/d_ndegf + pforce(1))/pmass
        bmoment(2) = bmoment(2) + quart_dt*(volume*(press(2) - press0)   &
                                    + 2.0_dp*ekin/d_ndegf + pforce(2))/pmass
        bmoment(3) = bmoment(3) + quart_dt*(volume*(press(3) - press0)   &
                                    + 2.0_dp*ekin/d_ndegf + pforce(3))/pmass
      else if (ensemble%ensemble == EnsembleNPgT) then
        pressxy0 = press0 - gamma0 / boundary%box_size_z
        bmoment(1) = bmoment(1) + quart_dt*(volume*(press(1) - pressxy0) &
                                    + 2.0_dp*ekin/d_ndegf + pforce(1))/pmass
        bmoment(2) = bmoment(2) + quart_dt*(volume*(press(2) - pressxy0) &
                                    + 2.0_dp*ekin/d_ndegf + pforce(2))/pmass
        bmoment(3) = bmoment(3) + quart_dt*(volume*(press(3) - press0)   &
                                    + 2.0_dp*ekin/d_ndegf + pforce(3))/pmass
      end if

    else if (ensemble%isotropy == IsotropyXY_Fixed) then

      bmoment(1) = 0.0_dp
      bmoment(2) = 0.0_dp
      bmoment(3) = bmoment(3) + quart_dt*(volume*(press(3) - press0) &
                              + 2.0_dp*ekin/d_ndegf + pforce(3))/pmass

    end if

    bmoment(1:3)  = exp(-gamma_p*quart_dt/2.0_dp)*bmoment(1:3)

    ! velocities are scaled according to scaled eta(t+1/4dt)
    !
    trace = bmoment(1)+bmoment(2)+bmoment(3)
    vel_scale(1:3) = exp(-half_dt*(bmoment(1:3)+trace/d_ndegf))
    do i = 1, ncell
      do ix = 1, natom(i)
        vel(1:3,ix,i) = vel(1:3,ix,i) * vel_scale(1:3)
      end do
    end do

    ! eta(t+1/4dt) is scaled
    !
    bmoment(1:3) = exp(-gamma_p*quart_dt/2.0_dp)*bmoment(1:3)

    ! eta(t+1/2dt)
    !
    if (ensemble%isotropy == IsotropyISO) then

      bmoment(1) = bmoment(1) + quart_dt*(3.0_dp*volume*(pressxyz - press0) &
                              + 6.0_dp*ekin/d_ndegf + pforce(1))/pmass
      bmoment(2) = bmoment(1)
      bmoment(3) = bmoment(1)

    else if (ensemble%isotropy == IsotropySEMI_ISO) then

      if (ensemble%ensemble == EnsembleNPT) then
        bmoment(1) = bmoment(1) + quart_dt*(volume*(pressxy - press0)   &
                                + 2.0_dp*ekin/d_ndegf + pforce(1))/pmass
        bmoment(2) = bmoment(2) + quart_dt*(volume*(pressxy - press0)   &
                                + 2.0_dp*ekin/d_ndegf + pforce(2))/pmass
        bmoment(3) = bmoment(3) + quart_dt*(volume*(press(3)  - press0) &
                                + 2.0_dp*ekin/d_ndegf + pforce(3))/pmass
      else if (ensemble%ensemble == EnsembleNPgT) then
        pressxy0 = press0 - gamma0 / boundary%box_size_z
        bmoment(1) = bmoment(1) + quart_dt*(volume*(pressxy - pressxy0) &
                                + 2.0_dp*ekin/d_ndegf + pforce(1))/pmass
        bmoment(2) = bmoment(2) + quart_dt*(volume*(pressxy - pressxy0) &
                                + 2.0_dp*ekin/d_ndegf + pforce(2))/pmass
        bmoment(3) = bmoment(3) + quart_dt*(volume*(press(3)  - press0) &
                                + 2.0_dp*ekin/d_ndegf + pforce(3))/pmass
      end if

    else if (ensemble%isotropy == IsotropyANISO) then

      if (ensemble%ensemble == EnsembleNPT) then
        bmoment(1) = bmoment(1) + quart_dt*(volume*(press(1) - press0)   &
                                    + 2.0_dp*ekin/d_ndegf + pforce(1))/pmass
        bmoment(2) = bmoment(2) + quart_dt*(volume*(press(2) - press0)   &
                                    + 2.0_dp*ekin/d_ndegf + pforce(2))/pmass
        bmoment(3) = bmoment(3) + quart_dt*(volume*(press(3) - press0)   &
                                    + 2.0_dp*ekin/d_ndegf + pforce(3))/pmass
      else if (ensemble%ensemble == EnsembleNPgT) then
        pressxy0 = press0 - gamma0 / boundary%box_size_z
        bmoment(1) = bmoment(1) + quart_dt*(volume*(press(1) - pressxy0) &
                                    + 2.0_dp*ekin/d_ndegf + pforce(1))/pmass
        bmoment(2) = bmoment(2) + quart_dt*(volume*(press(2) - pressxy0) &
                                    + 2.0_dp*ekin/d_ndegf + pforce(2))/pmass
        bmoment(3) = bmoment(3) + quart_dt*(volume*(press(3) - press0)   &
                                    + 2.0_dp*ekin/d_ndegf + pforce(3))/pmass
      end if

    else if (ensemble%isotropy == IsotropyXY_Fixed) then

      bmoment(1) = 0.0_dp
      bmoment(2) = 0.0_dp
      bmoment(3) = bmoment(3) + quart_dt*(volume*(press(3) - press0) &
                              + 2.0_dp*ekin/d_ndegf + pforce(3))/pmass

    end if

    bmoment(1:3)  = exp(-gamma_p*quart_dt/2.0_dp)*bmoment(1:3)

    return

  end subroutine update_barostat

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_barostat_mtk
  !> @brief        update barostat parameter bmoment for MTK
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_barostat_mtk(ensemble, pressx, pressy, pressz, pressxyz, &
                                 pressxy, press0, volume, d_ndegf, pmass,    &
                                 half_dt, ekin, bmoment)

    ! formal arguments
    type(s_ensemble),        intent(in)    :: ensemble
    real(dp),                intent(in)    :: pressx
    real(dp),                intent(in)    :: pressy
    real(dp),                intent(in)    :: pressz
    real(dp),                intent(in)    :: pressxyz
    real(dp),                intent(in)    :: pressxy
    real(dp),                intent(in)    :: press0
    real(dp),                intent(in)    :: volume
    real(dp),                intent(in)    :: d_ndegf
    real(dp),                intent(in)    :: pmass
    real(dp),                intent(in)    :: half_dt
    real(dp),                intent(in)    :: ekin
    real(dp),                intent(inout) :: bmoment(:)

    if (ensemble%isotropy == IsotropyISO) then

      bmoment(1) = bmoment(1) + half_dt*(3.0_dp*volume*(pressxyz - press0) &
                              + 6.0_dp*ekin/d_ndegf)/pmass
      bmoment(2) = bmoment(1)
      bmoment(3) = bmoment(1)

    else if (ensemble%isotropy == IsotropySEMI_ISO) then

      bmoment(1) = bmoment(1) + half_dt*(volume*(pressxy - press0)   &
                              + 2.0_dp*ekin/d_ndegf)/pmass
      bmoment(2) = bmoment(2) + half_dt*(volume*(pressxy - press0)   &
                              + 2.0_dp*ekin/d_ndegf)/pmass
      bmoment(3) = bmoment(3) + half_dt*(volume*(pressz  - press0)   &
                              + 2.0_dp*ekin/d_ndegf)/pmass

    else if (ensemble%isotropy == IsotropyANISO) then

      bmoment(1) = bmoment(1) + half_dt*(volume*(pressx - press0)    &
                              + 2.0_dp*ekin/d_ndegf)/pmass
      bmoment(2) = bmoment(2) + half_dt*(volume*(pressy - press0)    &
                              + 2.0_dp*ekin/d_ndegf)/pmass
      bmoment(3) = bmoment(3) + half_dt*(volume*(pressz - press0)    &
                              + 2.0_dp*ekin/d_ndegf)/pmass

    else if (ensemble%isotropy == IsotropyXY_Fixed) then

      bmoment(1) = 0.0_dp
      bmoment(2) = 0.0_dp
      bmoment(3) = bmoment(3) + half_dt*(volume*(pressz - press0) &
                              + 2.0_dp*ekin/d_ndegf)/pmass

    end if

    return

  end subroutine update_barostat_mtk

#ifdef HAVE_MPI_GENESIS
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine reduce_pres
  !> @brief
  !! @authors   JJ
  !! @param[in]
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine reduce_pres(val1, val2, val3)

    ! formal arguments
    real(dp), intent(inout) :: val1(:)
    real(dp), intent(inout) :: val2
    real(dp), intent(inout) :: val3(:)

    ! local variables
    double precision              :: before_reduce(7), after_reduce(7)

    before_reduce(1:3) = val1(1:3)
    before_reduce(4)   = val2
    before_reduce(5:7) = val3(1:3)

    call mpi_allreduce(before_reduce, after_reduce, 7, mpi_double_precision,   &
                       mpi_sum, mpi_comm_country, ierror)

    val1(1:3) = after_reduce(1:3)
    val2      = after_reduce(4)
    val3(1:3) = after_reduce(5:7)

    return

  end subroutine reduce_pres
#endif

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine bcast_boxsize
  !> @brief
  !! @authors   JJ
  !! @param[in]
  !! @date    2012/09/10 (JJ)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine bcast_boxsize(val1, val2, val3)

    ! formal arguments
    real(dp), intent(inout) :: val1, val2, val3

    ! local variables
    double precision              :: list(3)

    list(1)  = val1
    list(2)  = val2
    list(3)  = val3

    call mpi_bcast(list, 3, mpi_real8, 0, mpi_comm_country, ierror)

    val1 = list(1)
    val2 = list(2)
    val3 = list(3)

    return

  end subroutine bcast_boxsize

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine vverlet_respa_dynamics_fep(output, domain, enefunc, dynvars, &
                                    dynamics, pairlist, boundary,         &
                                    constraints, ensemble, comm,          &
                                    remd, alchemy)

    ! formal arguments
    type(s_output),            intent(inout) :: output
    type(s_domain),    target, intent(inout) :: domain
    type(s_enefunc),           intent(inout) :: enefunc
    type(s_dynvars),   target, intent(inout) :: dynvars
    type(s_dynamics),  target, intent(inout) :: dynamics
    type(s_pairlist),          intent(inout) :: pairlist
    type(s_boundary),          intent(inout) :: boundary
    type(s_constraints),       intent(inout) :: constraints
    type(s_ensemble),          intent(inout) :: ensemble
    type(s_comm),              intent(inout) :: comm
    type(s_remd),              intent(inout) :: remd
    type(s_alchemy), optional, intent(inout) :: alchemy

    ! local variables
    real(dp)                 :: simtim, temperature, factor
    real(dp)                 :: dt_short, dt_long
    real(dp)                 :: min_time
    integer                  :: i, j, k, jx, nsteps
    integer                  :: istep, multistep
    integer                  :: iseed, istart, iend
    integer                  :: min_k, k_min, k_max
    logical                  :: npt

    real(dp),        pointer :: coord(:,:,:), coord_ref(:,:,:)
    real(wp),        pointer :: coord_pbc(:,:,:)
    real(dp),        pointer :: coord_old(:,:,:)
    real(dp),        pointer :: vel(:,:,:), vel_ref(:,:,:), mass(:,:)
    real(wp),        pointer :: force_pbc(:,:,:,:), force_omp(:,:,:,:)
    real(dp),        pointer :: force(:,:,:)
    real(dp),        pointer :: virial_cell(:,:), virial(:,:)
    real(dp),        pointer :: force_short(:,:,:), force_long(:,:,:)
    integer,         pointer :: ncell, natom(:)
    logical,         pointer :: XI_RESPA, XO_RESPA

    ! fep/rest
    integer :: parmsetid, parmsetid_forward, parmsetid_backward
    integer :: replicaid
    type(s_energy), save :: frenergy


    ncell       => domain%num_cell_local
    natom       => domain%num_atom
    mass        => domain%mass
    coord       => domain%coord
    coord_ref   => domain%coord_ref
    coord_old   => domain%coord_old
    coord_pbc   => domain%translated
    force       => domain%force    
    force_omp   => domain%force_omp
    force_short => domain%force_short
    force_long  => domain%force_long
    vel         => domain%velocity
    vel_ref     => domain%velocity_ref
    force_pbc   => domain%force_pbc
    virial_cell => domain%virial_cellpair
    virial      => dynvars%virial
    XI_RESPA    => dynamics%xi_respa
    XO_RESPA    => dynamics%xo_respa

    temperature =  ensemble%temperature
    nsteps      =  dynamics%nsteps
    istart      =  dynamics%istart_step
    iend        =  dynamics%iend_step
    multistep   =  dynamics%elec_long_period
    dt_short    =  dynamics%timestep/AKMA_PS
    dt_long     =  real(multistep,dp)*dt_short
    simtim      =  dynamics%initial_time
    iseed       =  dynamics%iseed_init_velocity
    npt         = ensemble%use_barostat

    ! Setup for FEP
    !
    ! Copy lambda of enefunc to domain, because they are needed for 
    ! calculation of virial in SHAKE.
    domain%lambljA   = enefunc%lambljA
    domain%lambljB   = enefunc%lambljB
    domain%lambelA   = enefunc%lambelA
    domain%lambelB   = enefunc%lambelB
    domain%lambbondA = enefunc%lambbondA
    domain%lambbondB = enefunc%lambbondB
    domain%lambrest  = enefunc%lambrest
    ! synchronize single B with single A
    call sync_single_fep(domain, coord)
    call sync_single_fep(domain, vel)

    if (abs(dynamics%initial_rmsd) < 0.001_dp)  &
      dynamics%initial_rmsd = dynvars%energy%rmsd
    if (dynamics%target_md) enefunc%rmsd_force = 1 / (dt_long*dt_long)

    ! Check restart
    !
    if (.not. dynamics%restart) then

      call initial_velocity(temperature,           &
                            domain%num_atom_all,   &
                            domain%num_cell_local, &
                            domain%id_g2l,         &
                            domain%mass,           &
                            iseed,                 &
                            domain%velocity)

      call stop_trans_rotation_fep(domain%num_cell_local,         &
                               domain%num_atom,               &
                               dynamics%stop_com_translation, &
                               dynamics%stop_com_rotation,    &
                               domain%mass,                   &
                               domain%coord,                  &
                               domain%velocity, domain)

      call initial_vverlet_fep(npt, output, enefunc, dynamics,       &
                           pairlist, boundary, ensemble, constraints, &
                           domain, dynvars, comm)

    else

      ! After 2nd cycle of REMD simulation (istart /= 1), this is skipped
      !
      if (istart == 1) then

        if (constraints%tip4 .or. enefunc%table%tip4) &
        call decide_dummy(domain, constraints, coord)

        call communicate_coor(domain, comm)

        call compute_energy_fep(domain, enefunc, pairlist, boundary, coord,     &
                            npt, .false., .true., &
                            .false., .false., dynvars%energy, coord_pbc,     &
                            force_short, force_long, force_omp,     &
                            force_pbc, virial_cell, dynvars%virial, &
                            dynvars%virial_long, dynvars%virial_extern)
  
        call communicate_force(domain, comm, force_short)
        call communicate_force(domain, comm, force_long)

        if (constraints%tip4 .or. enefunc%table%tip4) then
          call water_force_redistribution(constraints, domain, force_long,  &
                                          virial)
          call water_force_redistribution(constraints, domain, force_short, &
                                          virial)
        end if
        
        ! FEP
        ! merge forces of single A and single B
        call add_single_fep(domain, force_short)
        call add_single_fep(domain, force_long)

      end if
    
    end if

    ! Main loop
    !
!    do i = 1, nsteps / multistep
    do i = int((istart-1)/multistep)+1, int(iend / multistep)

      simtim = simtim + dt_long * AKMA_PS
      dynvars%time = simtim
!      dynvars%step = i * multistep + istart - 1
      dynvars%step = i * multistep
      enefunc%rpath_sum_mf_flag = .false.

      ! Inner integrator
      !
      do j = 1, multistep

        ! decide the current steps
        !
        istep = (i-1) * multistep + j

        if (dynamics%target_md .or. dynamics%steered_md) &
          enefunc%rmsd_target = dynamics%initial_rmsd &
                              + (dynamics%final_rmsd-dynamics%initial_rmsd) &
                               *real(istep,dp)/real(nsteps,dp)

        ! save coordiantes for constraint
        !
        call timer(TimerIntegrator, TimerOn)
        call timer(TimerUpdate, TimerOn)
        do k = 1, ncell
          do jx = 1, natom(k)
            coord_ref(1:3,jx,k) = coord(1:3,jx,k)
            vel_ref  (1:3,jx,k) = vel  (1:3,jx,k)
          end do
        end do
        call timer(TimerUpdate, TimerOff)

        ! VV1
        !
        call integrate_vv1(dynamics, istep, j, dt_long, dt_short, ensemble, &
                           domain, constraints, boundary, dynvars)

        call timer(TimerIntegrator, TimerOff)

        ! cell migration and update cell pairlist
        !

        if (dynamics%nbupdate_period > 0 .and. &
            j == multistep .and. i > 0) then

          call domain_interaction_update_md_fep(istep, dynamics, domain, enefunc, &
                                          pairlist, boundary, constraints, comm)

        end if

        ! short range forces
        !
        call timer(TimerIntegrator, TimerOn)

        if (constraints%tip4 .or. enefunc%table%tip4) &
        call decide_dummy(domain, constraints, coord)

        call communicate_coor(domain, comm)

        call timer(TimerIntegrator, TimerOff)

        if (j < multistep) then

          pairlist%univ_update = 0
          enefunc%rpath_sum_mf_flag = enefunc%rpath_flag
          call compute_energy_short_fep(domain, enefunc, pairlist, boundary, coord, &
                                    npt, mod(istep,dynamics%eneout_period) == 0,&
                                    dynvars%energy,                             &
                                    coord_pbc,                                  &
                                    force_short,                                &
                                    force_omp,                                  &
                                    force_pbc,                                  &
                                    virial_cell,                                &
                                    dynvars%virial,                             &
                                    dynvars%virial_extern)

          call timer(TimerIntegrator, TimerOn)
          call timer(TimerComm2, TimerOn)
          call communicate_force(domain, comm, force_short)
          if (constraints%tip4 .or. enefunc%table%tip4) &
            call water_force_redistribution(constraints, domain, force_short, &
                                            virial)

          call timer(TimerComm2, TimerOff)

          ! FEP
          ! merge forces of single A and single B
          call add_single_fep(domain, force_short)

          call timer(TimerIntegrator, TimerOff)

        end if
        
        ! full step velocities with foce_short
        ! long range force for last inner step
        !
        if (j == multistep) then

          call compute_energy_fep(domain, enefunc, pairlist, boundary, coord, &
                          npt, .false., mod(istep,dynamics%eneout_period)==0, &
                          .false.,                                            &
                          .false.,                                            &
                          dynvars%energy,                                     &
                          coord_pbc,                                          &
                          force_short,                                        &
                          force_long,                                         &
                          force_omp,                                          &
                          force_pbc,                                          &
                          virial_cell,                                        &
                          dynvars%virial,                                     &
                          dynvars%virial_long,                                &
                          dynvars%virial_extern)

          ! FEP: compute and output energy differnces between adjacent states
          if ((domain%fep_use).and.(present(alchemy))) then
            if (dynamics%fepout_period > 0) then
              if ((mod(istep,dynamics%fepout_period) == 0) .and. &
                (mod(istep,dynamics%eneout_period) == 0)) then
                if (istep > dynamics%equilsteps) then
                  call compute_fep_energy(domain, enefunc, dynvars, pairlist, &
                    ensemble, boundary, alchemy, remd)
                  call output_fep_energy(output, enefunc, dynvars)
                end if
              end if
            end if
          end if

          call timer(TimerIntegrator, TimerOn)
          call timer(TimerComm2, TimerOn)
          call communicate_force(domain, comm, force_long)
          call communicate_force(domain, comm, force_short)
          if (constraints%tip4 .or. enefunc%table%tip4) then
            call water_force_redistribution(constraints, domain, force_long,  &
                                            virial)
            call water_force_redistribution(constraints, domain, force_short, &
                                            virial)
          end if

          call timer(TimerComm2, TimerOff)

          ! FEP
          ! merge forces of single A and single B
          call add_single_fep(domain, force_long)
          call add_single_fep(domain, force_short)

          call timer(TimerIntegrator, TimerOff)

        end if

        ! VV2
        !
        call integrate_vv2(dynamics, istep, j, dt_long, dt_short, ensemble,    &
                           domain, constraints, boundary, dynvars)

        ! Remove translational and rotational motion about COM(t + dt)
        !
        if (dynamics%stoptr_period > 0) then
          if (mod(istep,dynamics%stoptr_period) == 0) then

            call stop_trans_rotation_fep(domain%num_cell_local,         &
                                       domain%num_atom,               &
                                       dynamics%stop_com_translation, &
                                       dynamics%stop_com_rotation,    &
                                       domain%mass,                   &
                                       domain%coord,                  &
                                       domain%velocity, domain)

          end if
        end if

      end do

      ! energy output
      !
      if (dynamics%eneout_period > 0) then
        if (mod(istep,dynamics%eneout_period) == 0 ) then
          do k = 1, ncell
            do jx = 1, natom(k)
              force(1:3,jx,k) = force_short(1:3,jx,k) + force_long(1:3,jx,k)
            end do
          end do
        end if
      end if

      call output_md(output, dynamics, boundary, pairlist, &
                     ensemble, dynvars, domain, enefunc, remd, alchemy)


      ! output parallel I/O restart
      !
      call output_prst_md(output, enefunc, dynamics, boundary,                 &
                                  dynvars, domain, constraints)

    end do

    ! Close output files
    !
!    call close_output(output)

    return

  end subroutine vverlet_respa_dynamics_fep

  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine initial_vverlet_fep(npt, output, enefunc, dynamics, pairlist, &
                              boundary, ensemble, constraints, domain,  &
                              dynvars, comm)

    ! formal arguments
    logical,                 intent(in)    :: npt
    type(s_output),          intent(in)    :: output
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_constraints),     intent(inout) :: constraints
    type(s_domain),  target, intent(inout) :: domain
    type(s_dynvars), target, intent(inout) :: dynvars
    type(s_comm),            intent(inout) :: comm

    ! local variables
    real(dp)                 :: factor, temperature, energy(18), temp(18)
    real(dp)                 :: imass, simtim, dt, friction
    integer                  :: i, ix, j, jx, ncell, k, l

    real(dp),        pointer :: coord(:,:,:), coord_ref(:,:,:)
    real(wp),        pointer :: coord_pbc(:,:,:)
    real(dp),        pointer :: vel(:,:,:), vel_ref(:,:,:), mass(:,:)
    real(dp),        pointer :: force(:,:,:)
    real(wp),        pointer :: force_omp(:,:,:,:)
    real(wp),        pointer :: force_pbc(:,:,:,:)
    real(dp),        pointer :: virial_cell(:,:), virial(:,:)
    real(dp),        pointer :: force_short(:,:,:), force_long(:,:,:)
    integer,         pointer :: natom(:), nsolute(:), nwater(:)
    integer,         pointer :: water_list(:,:,:)


    natom       => domain%num_atom
    nwater      => domain%num_water
    nsolute     => domain%num_solute
    water_list  => domain%water_list
    coord       => domain%coord
    coord_ref   => domain%coord_ref
    coord_pbc   => domain%translated
    force       => domain%force
    force_short => domain%force_short
    force_long  => domain%force_long 
    force_omp   => domain%force_omp
    vel         => domain%velocity
    vel_ref     => domain%velocity_ref
    mass        => domain%mass
    force_pbc   => domain%force_pbc
    virial_cell => domain%virial_cellpair
    virial      => dynvars%virial

    dt          =  dynamics%timestep/AKMA_PS
    simtim      =  dynamics%initial_time
    temperature =  ensemble%temperature
    friction    =  ensemble%gamma_t * AKMA_PS
    ncell       =  domain%num_cell_local

    dynvars%time = simtim
    dynvars%step = 0

    ! Setup for FEP
    ! synchronize single B with single A
    call sync_single_fep(domain, coord)
    call sync_single_fep(domain, vel)

    ! save coordinates(0) and velocities(0)
    ! if rigid-body on, update coordinates(0)
    !
    do i = 1, ncell
      do ix = 1, natom(i)
        coord_ref(1:3,ix,i) = coord(1:3,ix,i)
        vel_ref(1:3,ix,i)   = vel(1:3,ix,i)
      end do
    end do

    if (constraints%rigid_bond) then

      call compute_constraints(ConstraintModeLEAP, .false., dt, coord_ref, &
                               domain, constraints, coord, vel,            &
                               dynvars%virial_const)

      do i = 1, ncell
        do ix = 1, natom(i)
          coord_ref(1:3,ix,i) = coord(1:3,ix,i)
        end do
      end do
    end if

    ! calculate energy(0) and forces(0)
    !
    if (constraints%tip4 .or. enefunc%table%tip4) &
    call decide_dummy(domain, constraints, coord)

    call communicate_coor(domain, comm)

    call compute_energy_fep(domain, enefunc, pairlist, boundary, coord,     &
                        npt, .false., .true.,                   &
                        .false., .false., dynvars%energy, coord_pbc,     &
                        force_short, force_long, force_omp,     &
                        force_pbc, virial_cell, dynvars%virial, &
                        dynvars%virial_long, dynvars%virial_extern)

    call communicate_force(domain, comm, force_short)
    call communicate_force(domain, comm, force_long)

    if (constraints%tip4 .or. enefunc%table%tip4) then
      call water_force_redistribution(constraints, domain, force_long,  &
                                      virial)
      call water_force_redistribution(constraints, domain, force_short, &
                                      virial)
    end if

    ! FEP
    ! merge forces of single A and single B
    call add_single_fep(domain, force_long)
    call add_single_fep(domain, force_short)

    do i = 1, ncell
      do ix = 1, natom(i)
        force(1:3,ix,i) = force_short(1:3,ix,i) + force_long(1:3,ix,i)
      end do
    end do

    ! if rigid-body on, update velocity(0 + dt/2 and 0 - dt/2)
    !
    if (constraints%rigid_bond) then

      ! calculate velocities(0 + dt/2) and coordinates(0 + dt)
      ! update coordinates(0 + dt) and velocities(0 + dt/2)
      !
      do i = 1, ncell
        do ix = 1, nsolute(i)
          imass = 1.0_dp/mass(ix,i)
          vel(1:3,ix,i) = vel_ref(1:3,ix,i) + 0.5_dp*dt*force(1:3,ix,i)*imass
          coord(1:3,ix,i) = coord_ref(1:3,ix,i) + dt*vel(1:3,ix,i)
        end do
        do l = 1, nwater(i)
          do k = 1, 3
            ix = water_list(k,l,i)
            imass = 1.0_dp/mass(ix,i)
            vel(1:3,ix,i) = vel_ref(1:3,ix,i) + 0.5_dp*dt*force(1:3,ix,i)*imass
            coord(1:3,ix,i) = coord_ref(1:3,ix,i) + dt*vel(1:3,ix,i)
          end do
        end do
      end do

      call compute_constraints(ConstraintModeVVER1, .false., dt, coord_ref,    &
                               domain, constraints, coord, vel,                &
                               dynvars%virial_const)

      ! calculate velocities(0 - dt/2) and coordinates(0 - dt)
      ! update coordinates(0 - dt) and velocities(0 - dt/2)
      !
      do i = 1, ncell
        do ix = 1, nsolute(i)
          imass = 1.0_dp/mass(ix,i)
          vel_ref(1:3,ix,i) = vel_ref(1:3,ix,i)                                &
                            - 0.5_dp*dt*force(1:3,ix,i)*imass
          coord(1:3,ix,i) = coord_ref(1:3,ix,i) - dt*vel_ref(1:3,ix,i)
        end do
        do l = 1, nwater(i)
          do k = 1, 3
            ix = water_list(k,l,i)
            imass = 1.0_dp/mass(ix,i)
            vel_ref(1:3,ix,i) = vel_ref(1:3,ix,i)                                &
                              - 0.5_dp*dt*force(1:3,ix,i)*imass
            coord(1:3,ix,i) = coord_ref(1:3,ix,i) - dt*vel_ref(1:3,ix,i)
          end do
        end do
      end do

      call compute_constraints(ConstraintModeVVER1, .false., -dt, coord_ref,   &
                               domain, constraints, coord, vel_ref,            &
                               dynvars%virial_const)

      ! calculate velocity(0)
      !
      do i = 1, ncell
        do ix = 1, natom(i)
          vel_ref(1:3,ix,i) = 0.5_dp*(vel(1:3,ix,i) + vel_ref(1:3,ix,i))
        end do
      end do

      ! vel <= updated velocities (0) and coord <= constrained
      ! coordinates(0)
      !
      do i = 1, ncell
        do ix = 1, natom(i)
          vel(1:3,ix,i) = vel_ref(1:3,ix,i)
          coord(1:3,ix,i) = coord_ref(1:3,ix,i)
        end do
      end do

    end if

    ! output dynvars(0)
    !
    dynvars%time = 0.0_dp
    dynvars%step = 0

    call compute_dynvars(enefunc, dynamics, boundary, ensemble, domain, &
                             dynvars)

    call output_dynvars(output, enefunc, dynvars, ensemble, boundary)


    ! at this point
    !   coord, velocity, and force are at 0

    return

  end subroutine initial_vverlet_fep

end module sp_md_respa_mod
