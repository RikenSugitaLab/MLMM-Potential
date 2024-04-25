!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_md_leapfrog_mod
!> @brief   perform molecular dynamics simulation
!! @authors Jaewoon Jung (JJ), Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_md_leapfrog_mod

  use sp_output_mod
  use sp_update_domain_mod
  use sp_assign_velocity_mod
  use sp_dynvars_mod
  use sp_constraints_mod
  use sp_parallel_io_mod
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
  use sp_enefunc_gamd_mod
  use sp_remd_str_mod    ! fep/rest
  use sp_alchemy_str_mod ! fep/rest
  use sp_energy_pme_mod
  use random_mod
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
  public  :: leapfrog_dynamics
  private :: initial_leapfrog
  private :: control_temp_pres_leap
  private :: berendsen_leapfrog
  private :: nose_hoover_leapfrog
  private :: gaussian_leapfrog
  private :: langevin_leapfrog_nvt
  private :: langevin_leapfrog_npt
  private :: simulated_annealing_leapfrog
  private :: reduce_pres
  private :: bcast_boxsize

  ! FEP
  public  :: leapfrog_dynamics_fep
  private :: initial_leapfrog_fep

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    leapfrog_dynamics
  !> @brief        leapfrog integrator with domain decomposition
  !! @authors      JJ, TM
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

  subroutine leapfrog_dynamics(output, domain, enefunc, dynvars, dynamics, &
                               pairlist, boundary, constraints, ensemble,  &
                               comm, remd)

    ! formal arguments
    type(s_output),          intent(inout) :: output
    type(s_domain),  target, intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_dynvars), target, intent(inout) :: dynvars
    type(s_dynamics),        intent(inout) :: dynamics
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_boundary),        intent(inout) :: boundary
    type(s_constraints),     intent(inout) :: constraints
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_comm),            intent(inout) :: comm   
    type(s_remd),            intent(inout) :: remd

    ! local variables
    real(dp)                 :: simtim, dt, temperature, factor
    real(dp)                 :: energy(18), temp(18), temp_prev(18)
    integer                  :: i, j, k, jx, nsteps, ekin, iseed
    integer                  :: istart, iend
    logical                  :: npt

    real(dp),        pointer :: coord(:,:,:), coord_ref(:,:,:)
    real(dp),        pointer :: coord_old(:,:,:)
    real(wp),        pointer :: coord_pbc(:,:,:)
    real(dp),        pointer :: vel(:,:,:), vel_ref(:,:,:)
    real(dp),        pointer :: force(:,:,:), force_long(:,:,:)
    real(dp),        pointer :: mass(:,:)
    real(wp),        pointer :: force_omp(:,:,:,:)
    real(wp),        pointer :: force_pbc(:,:,:,:)
    real(dp),        pointer :: virial_cell(:,:)
    real(dp),        pointer :: virial(:,:)
    integer,         pointer :: ncell, nboundary, natom(:)


    ncell       => domain%num_cell_local
    nboundary   => domain%num_cell_boundary
    natom       => domain%num_atom
    mass        => domain%mass
    coord       => domain%coord
    coord_ref   => domain%coord_ref
    coord_old   => domain%coord_old
    coord_pbc   => domain%translated
    force       => domain%force
    force_long  => domain%force_long
    force_omp   => domain%force_omp
    force_pbc   => domain%force_pbc
    virial_cell => domain%virial_cellpair
    vel         => domain%velocity
    vel_ref     => domain%velocity_ref
    virial      => dynvars%virial

    temperature =  ensemble%temperature
    npt         =  ensemble%use_barostat
    nsteps      =  dynamics%nsteps
    istart      =  dynamics%istart_step
    iend        =  dynamics%iend_step
    dt          =  dynamics%timestep/AKMA_PS
    simtim      =  dynamics%initial_time
    iseed       =  dynamics%iseed_init_velocity


    if (abs(dynamics%initial_rmsd) < 0.001_wp)  &
      dynamics%initial_rmsd = dynvars%energy%rmsd
    if (dynamics%target_md) enefunc%rmsd_force = 1 / (dt*dt)

    ! Check restart
    !
    if (.not. dynamics%restart .and. .not. pio_restart) then

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

      call initial_leapfrog(npt, output, enefunc, dynamics,            &
                            pairlist, boundary, ensemble, constraints, &
                            domain, dynvars, comm)
    end if


    ! Main loop 
    !   coord is at 0 +  dt and vel is at 0 + 1/2dt, if restart off
    !   coord is at t + 2dt and vel is at t + 3/2dt, if restart on
    !
    do i = istart, iend

      call timer(TimerIntegrator, TimerOn)

      simtim = simtim + dynamics%timestep
      dynvars%time = simtim
      dynvars%step = i
      if (dynamics%target_md .or. dynamics%steered_md) &
        enefunc%rmsd_target = dynamics%initial_rmsd &
                            + (dynamics%final_rmsd-dynamics%initial_rmsd) &
                             *real(dynvars%step,wp)/real(nsteps,wp)
      enefunc%rpath_sum_mf_flag = enefunc%rpath_flag

      ! send the coordinate data
      !
      call timer(TimerComm1, TimerOn)

      if (constraints%tip4 .or. enefunc%table%tip4) &
      call decide_dummy(domain, constraints, coord)

      call communicate_coor(domain, comm)

      call timer(TimerComm1, TimerOff)

      ! Save coordinates(t + dt) and velocities(t + 1/2dt)
      !
      !$omp parallel do default(shared) private(j,jx)
      !
      do j = 1, ncell+nboundary
        do jx = 1, natom(j)
          coord_ref(1:3,jx,j) = coord(1:3,jx,j)
          vel_ref  (1:3,jx,j) = vel  (1:3,jx,j)
        end do
      end do
      !$omp end parallel do
      call timer(TimerIntegrator, TimerOff)

      ! Compute energy(t + dt), force(t + dt), and internal virial(t + dt)
      !
      call compute_energy(domain, enefunc, pairlist, boundary, coord, &
                          npt, .false., mod(i,dynamics%eneout_period)==0, &
                          .true.,                  &
                          enefunc%nonb_limiter,    &
                          dynvars%energy,          &
                          coord_pbc,               &
                          force,                   &
                          force_long,              &
                          force_omp,               &
                          force_pbc,               &
                          virial_cell,             &
                          dynvars%virial,          &
                          dynvars%virial_long, &
                          dynvars%virial_extern)

      call timer(TimerIntegrator, TimerOn)
      call timer(TimerComm2, TimerOn)
                          
      call communicate_force(domain, comm, force)
      if (constraints%tip4 .or. enefunc%table%tip4) &
        call water_force_redistribution(constraints, domain, force, virial)

      call timer(TimerComm2, TimerOff)

      if (ensemble%tpcontrol /= TpcontrolLangevin) then

        ! Newtonian dynamics
        !   v(t+3/2dt) = v(t+1/2dt) + dt*F(t+dt)/m
        !   r(t+2dt) = r(t+dt) + dt*v(t+3/2dt)
        !
        !$omp parallel do default(shared) private(j, jx, factor)
        do j = 1, ncell
          do jx = 1, natom(j)
            factor = dt/mass(jx,j)
            vel  (1:3,jx,j) = vel  (1:3,jx,j) + factor*force(1:3,jx,j)
            coord(1:3,jx,j) = coord(1:3,jx,j) + dt*vel(1:3,jx,j)
          end do
        end do
        !$omp end parallel do


        ! Bond constraint
        !   coord (unconstrained) is at t + 2dt, vel     is at t + 3/2dt
        !   coord_ref             is at t +  dt, vel_ref is at t + 1/2dt
        !   compute constrained coordinates(t+2dt) and constraint virial(t+dt)
        !   update velocities(t + 3/2dt):
        !     v(t+3/2dt) = (r(t+2dt) - r(t+dt))/dt
        !   add constraint virial(t + dt)
        !
        if (constraints%rigid_bond) then
          call compute_constraints(ConstraintModeLEAP, .true., dt, coord_ref, &
                                   domain, constraints, coord, vel,  &
                                   dynvars%virial_const)

          dynvars%virial(1:3,1:3) = dynvars%virial(1:3,1:3) &
                                  + dynvars%virial_const(1:3,1:3)
        end if

        ! Control temperature and pressure
        !   coord     is at t + 2dt, vel     is at t + 3/2dt
        !   coord_ref is at t +  dt, vel_ref is at t + 1/2dt
        !   scale velocities(t + 3/2dt) and coordinates(t + 2dt)
        !

        if (ensemble%ensemble /= EnsembleNVE) then

          call control_temp_pres_leap(dynamics, ensemble, domain, constraints, &
                                      boundary, dynvars)

        end if

      else

        ! Langevin dynamics
        !
        if (ensemble%ensemble == EnsembleNVT) then

          call langevin_leapfrog_nvt(dynamics, ensemble, domain, &
                                     constraints, dynvars)

        else if (ensemble%ensemble == EnsembleNPT  .or. &
                 ensemble%ensemble == EnsembleNPAT .or. &
                 ensemble%ensemble == EnsembleNPgT) then

          call langevin_leapfrog_npt(dynamics, ensemble, domain, &
                                     constraints, boundary, dynvars)

        end if

      end if

      call timer(TimerIntegrator, TimerOff)


      if (dynamics%stoptr_period > 0) then

        if (mod(i,dynamics%stoptr_period) == 0) then

          do j = 1, ncell
            do jx = 1, natom(j)
              coord_old(1:3,jx,j) = 0.5_dp*(coord(1:3,jx,j)+coord_ref(1:3,jx,j))
            end do
          end do
          call stop_trans_rotation(ncell, natom,                  &
                                   dynamics%stop_com_translation, &
                                   dynamics%stop_com_rotation,    &
                                   mass, coord_old, vel)
        end if

      end if

      ! ------------------------------------------------------------------------
      ! OUTPUT energy, trajectory, and restart file
      !   coord     is at t + 2dt, vel      is at t + 3/2dt
      !   coord_ref is at t +  dt, vel_ref  is at t +    dt
      !   volume    is at t +  dt, box_size is at t +   2dt

      call random_push_stock

      ! output dynvars(t + dt)
      !
      call output_md(output, dynamics, boundary, pairlist, &
                     ensemble, dynvars, domain, enefunc, remd)

      ! update interaction
      !
      if (dynamics%nbupdate_period > 0) &
        call domain_interaction_update_md(i, dynamics, domain, enefunc, &
                                          pairlist, boundary, constraints, comm)

      ! Simulated annealing
      !
      if (dynamics%anneal_period > 0) then
        if (mod(i,dynamics%anneal_period) == 0) then

          call simulated_annealing_leapfrog(dynamics, ensemble)

        end if
      end if

      ! Update GAMD
      !
      if (enefunc%gamd%update_period > 0) then
        if (mod(i,enefunc%gamd%update_period) == 0) then
          call update_gamd_leapfrog(output, enefunc, dynvars)
        end if
      end if

      ! output parallel I/O restart
      !
      call output_prst_md(output, enefunc, dynamics, boundary, &
                                  dynvars, domain, constraints)
      
    end do


    return

  end subroutine leapfrog_dynamics

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    initial_leapfrog
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

  subroutine initial_leapfrog(npt, output, enefunc, dynamics, pairlist, &
                              boundary, ensemble, constraints, domain,  &
                              dynvars, comm)

    ! formal arguments
    logical,                 intent(in)    :: npt
    type(s_output),          intent(in)    :: output
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_dynamics),        intent(inout) :: dynamics
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_constraints),     intent(inout) :: constraints
    type(s_domain),  target, intent(inout) :: domain
    type(s_dynvars), target, intent(inout) :: dynvars
    type(s_comm),            intent(inout) :: comm

    ! local variables
    real(dp)                 :: factor, temperature, energy(18)
    real(dp)                 :: temp(18),ekin1, ekin2
    real(dp)                 :: imass, simtim, dt, friction
    integer                  :: i, ix, j, jx, ncell, k, l

    real(dp),        pointer :: coord(:,:,:), coord_ref(:,:,:)
    real(wp),        pointer :: coord_pbc(:,:,:)
    real(dp),        pointer :: vel(:,:,:), vel_ref(:,:,:)
    real(dp),        pointer :: force(:,:,:), force_long(:,:,:)
    real(dp),        pointer :: mass(:,:)
    real(wp),        pointer :: force_omp(:,:,:,:)
    real(wp),        pointer :: force_pbc(:,:,:,:)
    real(dp),        pointer :: virial_cell(:,:), virial(:,:)
    integer,         pointer :: natom(:), nsolute(:), nwater(:)
    integer,         pointer :: water_list(:,:,:)


    natom        => domain%num_atom
    nsolute      => domain%num_solute
    nwater       => domain%num_water 
    water_list   => domain%water_list
    coord        => domain%coord
    coord_ref    => domain%coord_ref
    coord_pbc    => domain%translated
    force        => domain%force
    force_long   => domain%force_long
    force_omp    => domain%force_omp
    force_pbc    => domain%force_pbc
    virial_cell  => domain%virial_cellpair
    vel          => domain%velocity
    vel_ref      => domain%velocity_ref
    mass         => domain%mass
    virial       => dynvars%virial

    ncell        =  domain%num_cell_local
    dt           =  dynamics%timestep/AKMA_PS
    simtim       =  dynamics%initial_time
    temperature  =  ensemble%temperature
    friction     =  ensemble%gamma_t * AKMA_PS

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
                               domain, constraints, coord, vel,  &
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

    call compute_energy(domain, enefunc, pairlist, boundary, coord,  &
                        npt, .false., .true.,  &
                        .true.,                &
                        enefunc%nonb_limiter,  &
                        dynvars%energy,        &
                        coord_pbc,             &
                        force,                 &
                        force_long,            &
                        force_omp,             &
                        force_pbc,             &
                        virial_cell,           &
                        dynvars%virial,        &
                        dynvars%virial_long,   &
                        dynvars%virial_extern)

    ! get forces by communication
    !
    call communicate_force(domain, comm, force)

    if (constraints%tip4 .or. enefunc%table%tip4) &
      call water_force_redistribution(constraints, domain, force, virial)

    ! Calculate velocities(0 - dt/2)
    !
    if (.not. constraints%rigid_bond) then

     do i = 1, ncell
       do ix = 1, nsolute(i)
         imass = 1.0_wp/mass(ix,i)
         vel(1:3,ix,i) = vel(1:3,ix,i) - 0.5_wp*dt*force(1:3,ix,i)*imass
       end do
       do l = 1, nwater(i)
         do k = 1, 3
           ix = water_list(k,l,i)
           imass = 1.0_wp/mass(ix,i)
           vel(1:3,ix,i) = vel(1:3,ix,i) - 0.5_wp*dt*force(1:3,ix,i)*imass
         end do
       end do
     end do

    else

      ! calculate velocities(0 + dt/2) and coordinates(0 + dt)
      ! update coordinates(0 + dt) and velocities(0 + dt/2)
      !
      do i = 1, ncell
        do ix = 1, nsolute(i)
          imass = 1.0_wp/mass(ix,i)
          vel(1:3,ix,i) = vel_ref(1:3,ix,i) + 0.5_wp*dt*force(1:3,ix,i)*imass
          coord(1:3,ix,i) = coord_ref(1:3,ix,i) + dt*vel(1:3,ix,i)
        end do
        do l = 1, nwater(i)
          do k = 1, 3
            ix = water_list(k,l,i)
            imass = 1.0_wp/mass(ix,i)
            vel(1:3,ix,i) = vel_ref(1:3,ix,i) + 0.5_wp*dt*force(1:3,ix,i)*imass
            coord(1:3,ix,i) = coord_ref(1:3,ix,i) + dt*vel(1:3,ix,i)
          end do
        end do  
      end do

      call compute_constraints(ConstraintModeLEAP, .true., dt, coord_ref, &
                               domain, constraints, coord, vel,           &
                               dynvars%virial_const)

      ! calculate velocities(0 - dt/2) and coordinates(0 - dt)
      ! update coordinates(0 - dt) and velocities(0 - dt/2)
      !
      do i = 1, ncell
        do ix = 1, nsolute(i)
          imass = 1.0_wp/mass(ix,i)
          vel_ref(1:3,ix,i) = vel_ref(1:3,ix,i)                            &
                             - 0.5_wp*dt*force(1:3,ix,i)*imass
          coord(1:3,ix,i) = coord_ref(1:3,ix,i) - dt*vel_ref(1:3,ix,i)
        end do
        do l = 1, nwater(i)
          do k = 1, 3
            ix = water_list(k,l,i)
            imass = 1.0_wp/mass(ix,i)
            vel_ref(1:3,ix,i) = vel_ref(1:3,ix,i)                            &
                               - 0.5_wp*dt*force(1:3,ix,i)*imass
            coord(1:3,ix,i) = coord_ref(1:3,ix,i) - dt*vel_ref(1:3,ix,i)
          end do
        end do
      end do

      call compute_constraints(ConstraintModeLEAP, .true., -dt, coord_ref, &
                               domain, constraints, coord, vel_ref,        &
                               dynvars%virial_const)

      ! vel <= vel(0 - dt/2) and coord <= coord(0)
      !
      do i = 1, ncell
        do ix = 1, natom(i)
          vel(1:3,ix,i)   = vel_ref(1:3,ix,i)
          coord(1:3,ix,i) = coord_ref(1:3,ix,i)
        end do
      end do

    end if

    ! first step dynamics
    !   Calculate coordinates(0 + dt) and velocities(0 + 1/2dt)
    !   Note: this routine corresponds to i = 0 in the main loop
    !

    ! coord_ref <= coord(0) and vel_ref <= vel(0 - 1/2dt)
    !
    do i = 1, ncell
      do ix = 1, natom(i)
        coord_ref(1:3,ix,i) = coord(1:3,ix,i)
        vel_ref(1:3,ix,i)   = vel(1:3,ix,i)
      end do
    end do

    ! calculate energy(0) and forces(0)
    !
    if (constraints%tip4 .or. enefunc%table%tip4) &
    call decide_dummy(domain, constraints, coord)

    call communicate_coor(domain, comm)

    call compute_energy(domain, enefunc, pairlist, boundary, coord,  &
                        npt, .false., .true.,  &
                        .true.,                &
                        enefunc%nonb_limiter,  &
                        dynvars%energy,        &
                        coord_pbc,             &
                        force,                 &
                        force_long,            &
                        force_omp,             &
                        force_pbc,             &
                        virial_cell,           &
                        dynvars%virial,        &
                        dynvars%virial_long,   &
                        dynvars%virial_extern)

    call communicate_force(domain, comm, force)

    if (constraints%tip4 .or. enefunc%table%tip4) &
      call water_force_redistribution(constraints, domain, force, virial)

    ! Newtonian dynamics
    !   v(0+1/2dt) = v(0-1/2dt) + dt*F(0)/m
    !   r(0+dt) = r(0) + dt*v(0+1/2dt)
    !
    do i = 1, ncell
      do ix = 1, nsolute(i)
        factor = dt/mass(ix,i)
        vel(1:3,ix,i) = vel(1:3,ix,i) + factor*force(1:3,ix,i)
        coord(1:3,ix,i) = coord(1:3,ix,i) + dt*vel(1:3,ix,i)
      end do
      do l = 1, nwater(i)
        do k = 1, 3
          ix = water_list(k,l,i)
          factor = dt/mass(ix,i)
          vel(1:3,ix,i) = vel(1:3,ix,i) + factor*force(1:3,ix,i)
          coord(1:3,ix,i) = coord(1:3,ix,i) + dt*vel(1:3,ix,i)
        end do
      end do
    end do


    ! SHAKE and SETTLE
    !   calculate constrained coordinates(0 + dt) and constraint virial(0)
    !   update velocities(0 + 1/2dt):
    !     v(0+1/2dt) = (r(0+dt) - r(0))/dt
    !
    if (constraints%rigid_bond) then
      call compute_constraints(ConstraintModeLEAP, .true., dt, coord_ref, &
                               domain, constraints, coord, vel, &
                               dynvars%virial_const)

      dynvars%virial(1:3,1:3) = &
        dynvars%virial(1:3,1:3) + dynvars%virial_const(1:3,1:3)

    end if

    ! output dynvars(0)
    !
    call compute_dynvars(enefunc, dynamics, boundary, ensemble, domain, dynvars)

    call output_dynvars(output, enefunc, dynvars, ensemble, boundary)

    dynamics%restart = .true.

    ! at this point
    !   coord is at 0 + dt, and vel is at 0 + 1/2dt

    return

  end subroutine initial_leapfrog
 
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    control_temp_pres_leap
  !> @brief        driver to control temperature and pressure
  !! @authors      TM, JJ
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] boundary    : boundary information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine control_temp_pres_leap(dynamics, ensemble, domain, constraints, &
                                    boundary, dynvars)

    ! formal arguments
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_domain),          intent(inout) :: domain
    type(s_constraints),     intent(inout) :: constraints
    type(s_boundary),        intent(inout) :: boundary
    type(s_dynvars),         intent(inout) :: dynvars 


    select case (ensemble%tpcontrol)

    case (TpcontrolNoseHoover)

      call nose_hoover_leapfrog(dynamics, ensemble, domain, constraints, &
                                boundary, dynvars)

    case (TpcontrolBerendsen)

      call berendsen_leapfrog  (dynamics, ensemble, domain, constraints, &
                                boundary, dynvars)

    case (TpcontrolGauss)

      call gaussian_leapfrog   (dynamics, ensemble, domain, constraints, &
                                dynvars)

    end select

    return

  end subroutine control_temp_pres_leap

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    nose_hoover_leapfrog
  !> @brief        control temperature and pressure using Nose-Hoover method
  !! @authors      TM, JJ
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] boundary    : boundary information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine nose_hoover_leapfrog(dynamics, ensemble, domain, constraints, &
                                  boundary, dynvars)

    ! formal arguments
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_domain),  target, intent(inout) :: domain
    type(s_constraints),     intent(inout) :: constraints
    type(s_boundary),        intent(inout) :: boundary
    type(s_dynvars), target, intent(inout) :: dynvars

    ! local variables
    real(dp)                 :: dt, inv_dt, temp0, press0, tau_t
    real(dp)                 :: kin(1:3)
    real(dp)                 :: volume, press(1:3)
    real(dp)                 :: vel_tmp(1:3), ekin, tpress, tempt
    real(dp)                 :: scale_v, scale_c(3), factor
    real(dp)                 :: sigma, qmass
    integer                  :: i, ix, k, l, num_degree

    real(dp),        pointer :: fric_ref, fric
    real(dp),        pointer :: coord(:,:,:), coord_ref(:,:,:), force(:,:,:)
    real(dp),        pointer :: vel(:,:,:), vel_ref(:,:,:), virial(:,:)
    real(dp),        pointer :: mass(:,:)
    integer,         pointer :: ncell, natom(:), nsolute(:), nwater(:)
    integer,         pointer :: water_list(:,:,:)


    natom      => domain%num_atom
    nsolute    => domain%num_solute
    nwater     => domain%num_water 
    water_list => domain%water_list
    coord      => domain%coord
    coord_ref  => domain%coord_ref
    vel        => domain%velocity
    vel_ref    => domain%velocity_ref
    force      => domain%force
    mass       => domain%mass
    ncell      => domain%num_cell_local

    virial     => dynvars%virial
    fric_ref   => dynvars%thermostat_momentum_ref
    fric       => dynvars%thermostat_momentum

    dt         =  dynamics%timestep/AKMA_PS
    inv_dt     =  1.0_dp/dt
    num_degree =  domain%num_deg_freedom
    temp0      =  ensemble%temperature
    press0     =  ensemble%pressure
    tau_t      =  ensemble%tau_t/AKMA_PS


    ! save friction coeff(t+1/2dt)
    !
    dynvars%thermostat_momentum_ref = dynvars%thermostat_momentum

    ! compute kinetic energy(t+dt)
    !   v(t+dt) = 0.5*(v(t+3/2dt) + v(t+1/2dt))
    !
    kin(1:3) = 0.0_dp

    do i = 1, ncell
      do ix = 1, natom(i)
        vel_tmp(1:3) = 0.5_dp*(vel(1:3,ix,i) + vel_ref(1:3,ix,i))
        kin(1:3) = kin(1:3) + mass(ix,i) * vel_tmp(1:3)*vel_tmp(1:3)
      end do
    end do
    ekin   = 0.5_dp * (kin(1)+kin(2)+kin(3))

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, ekin, 1, mpi_real8, &
                       mpi_sum, mpi_comm_country, ierror)
#endif

    ! thermostat
    !   fric(t+3/2dt) = fric(t+1/2dt) + dt*[2Ekin(t+dt)-2s]/qmass
    !   fric(t+dt)    = [fric(t+3/2dt) + fric(t+1/2dt)]/2
    !   v(t+dt)    = [v(t+3/2dt) + v(t+1/2dt)]/2
    !   v(t+3/2dt) = v(t+1/2dt) + dt*[F(t+dt)/m - fric(t+dt)*v(t+dt)]
    !   r(t+2dt)   = r(t+dt) + v(t+3/2dt)*dt
    !
    sigma = 0.5_dp * num_degree * KBOLTZ * temp0
    qmass = 2.0_dp * sigma * tau_t**2
    fric  = fric_ref + dt*(2.0_dp*ekin - 2.0_dp*sigma)/qmass

    factor = 0.5_dp*(fric_ref + fric)
    do i = 1, ncell
      do ix = 1, nsolute(i)
        vel_tmp(1:3)   = 0.5_dp*(vel(1:3,ix,i) + vel_ref(1:3,ix,i))
        vel(1:3,ix,i)  = vel_ref(1:3,ix,i) +  &
                         (force(1:3,ix,i)/mass(ix,i)-factor*vel_tmp(1:3))*dt
        coord(1:3,ix,i) = coord_ref(1:3,ix,i) + vel(1:3,ix,i)*dt
      end do
      do l = 1, nwater(i)
        do k = 1, 3
          ix = water_list(k,l,i) 
          vel_tmp(1:3)   = 0.5_dp*(vel(1:3,ix,i) + vel_ref(1:3,ix,i))
          vel(1:3,ix,i)  = vel_ref(1:3,ix,i) +  &
                           (force(1:3,ix,i)/mass(ix,i)-factor*vel_tmp(1:3))*dt
          coord(1:3,ix,i) = coord_ref(1:3,ix,i) + vel(1:3,ix,i)*dt
        end do
      end do
    end do


    ! SHAKE and SETTLE
    !   compute constrained coordinates(t + 2dt)
    !   update velocities(t + 3/2dt):
    !     vel(t+3/2dt) = (r(t+2dt) - r(t+dt))/dt
    !
    if (constraints%rigid_bond) then
      call compute_constraints(ConstraintModeLEAP, .true., dt, coord_ref, &
                               domain,  constraints, coord, vel, &
                               dynvars%virial_const)

    end if

    return

  end subroutine nose_hoover_leapfrog

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    berendsen_leapfrog
  !> @brief        control temperature and pressure using Berendsen method
  !! @authors      TM, JJ
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] boundary    : boundary information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine berendsen_leapfrog(dynamics, ensemble, domain, constraints, &
                                boundary, dynvars)

    ! formal arguments
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_domain),  target, intent(inout) :: domain
    type(s_constraints),     intent(inout) :: constraints
    type(s_boundary),        intent(inout) :: boundary
    type(s_dynvars), target, intent(inout) :: dynvars

    ! local variables
    real(dp)                 :: dt, inv_dt, temp0, press0, gamma0, pressxy0
    real(dp)                 :: kin(1:3), tau_t, tau_p, compress
    real(dp)                 :: volume, press(1:3)
    real(dp)                 :: vel_tmp(1:3), ekin, tpress, tempt
    real(dp)                 :: scale_v, scale_c(3), factor
    real(dp)                 :: virial_sum(1:3)
    integer                  :: i, ix, j, k, l, ij, num_degree

    real(dp),        pointer :: coord(:,:,:), coord_ref(:,:,:), force(:,:,:)
    real(wp),        pointer :: trans(:,:,:)
    real(dp),        pointer :: vel(:,:,:), vel_ref(:,:,:)
    real(dp),        pointer :: mass(:,:)
    real(dp),        pointer :: virial(:,:)
    integer,         pointer :: ncell, nboundary
    integer,         pointer :: natom(:), nsolute(:), nwater(:)
    integer,         pointer :: water_list(:,:,:)


    ncell      => domain%num_cell_local
    nboundary  => domain%num_cell_boundary
    mass       => domain%mass
    natom      => domain%num_atom
    nsolute    => domain%num_solute
    nwater     => domain%num_water 
    water_list => domain%water_list
    coord      => domain%coord
    coord_ref  => domain%coord_ref
    vel        => domain%velocity
    vel_ref    => domain%velocity_ref
    force      => domain%force
    virial     => dynvars%virial
    trans      => domain%trans_vec

    dt         =  dynamics%timestep/AKMA_PS
    inv_dt     =  1.0_dp/dt
    num_degree =  domain%num_deg_freedom
    temp0      =  ensemble%temperature
    press0     =  ensemble%pressure 
    tau_t      =  ensemble%tau_t/AKMA_PS
    tau_p      =  ensemble%tau_p/AKMA_PS
    compress   =  ensemble%compressibility

    ! save box size(t+dt) and compute volume(t+dt)
    !
    boundary%box_size_x_ref = boundary%box_size_x
    boundary%box_size_y_ref = boundary%box_size_y
    boundary%box_size_z_ref = boundary%box_size_z
    volume = boundary%box_size_x * boundary%box_size_y * boundary%box_size_z

    ! compute temperature(t+dt)
    !   v(t+dt) = 0.5*(v(t+3/2dt) + v(t+1/2dt))
    !
    kin(1:3) = 0.0_dp

    do i = 1, ncell
      do ix = 1, natom(i)
        vel_tmp(1:3) = 0.5_dp*(vel(1:3,ix,i) + vel_ref(1:3,ix,i))
        kin(1:3) = kin(1:3) + mass(ix,i) * vel_tmp(1:3)*vel_tmp(1:3)
      end do
    end do

    ekin   = 0.5_dp * (kin(1)+kin(2)+kin(3))

    call reduce_pres(kin, ekin, virial(1,1), virial(2,2), &
                     virial(3,3), virial_sum)
    tempt  = 2.0_dp * ekin / (num_degree * KBOLTZ)


    ! thermostat
    !   scale_v    = sqrt(1+dt/tau_t(temp_ref/temp(t+dt) -1))
    !   v(t+3/2dt) = (v(t+1/2dt) + dtF(t+dt)/m) * scale_v
    !   r(t+2dt)   = r(t+dt) + v(t+3/2dt)*dt
    !
    scale_v = sqrt(1.0_dp + dt * (temp0/tempt - 1.0_dp)/tau_t)

    do i = 1, ncell
      do ix = 1, nsolute(i)
        factor   = dt/mass(ix,i)
        vel(1:3,ix,i) = (vel_ref(1:3,ix,i) + factor*force(1:3,ix,i))*scale_v
        coord(1:3,ix,i) = coord_ref(1:3,ix,i) + vel(1:3,ix,i)*dt
      end do
      do l = 1, nwater(i)
        do k = 1, 3
          ix = water_list(k,l,i)
          factor   = dt/mass(ix,i)
          vel(1:3,ix,i) = (vel_ref(1:3,ix,i) + factor*force(1:3,ix,i))*scale_v
          coord(1:3,ix,i) = coord_ref(1:3,ix,i) + vel(1:3,ix,i)*dt
        end do
      end do
    end do

    ! SHAKE and SETTLE
    !   compute constrained coordinates(t + 2dt)
    !   update velocities(t + 3/2dt):
    !     vel(t+3/2dt) = (r(t+2dt) - r(t+dt))/dt
    !
    if (constraints%rigid_bond) then
      call compute_constraints(ConstraintModeLEAP, .true., dt, coord_ref, &
                               domain, constraints, coord, vel, &
                               dynvars%virial_const)
    end if

    if (ensemble%use_barostat) then

      ! compute pressure(t+dt)
      !
      press(1:3) = P_ATMOS * (kin(1:3) + virial_sum(1:3)) / volume
      tpress = (press(1) + press(2) + press(3))/3.0_dp

      ! barostat
      !   scale coordinates(t+2dt) and box_size(t+2dt)
      !   scale_c  = (1- compress*dt*(Pext - P(t+dt)/tau_p))**1/3
      !   r(t+2dt)    = scale_c * r(t+2dt)
      !   size(t+2dt) = scale_c * size(t+dt)
      !
      if (ensemble%isotropy == IsotropyISO) then

        scale_c(1) = (1.0_dp - compress*dt*(press0 - tpress)/tau_p)**ONE_THIRD
        scale_c(2) = scale_c(1)
        scale_c(3) = scale_c(1)

      else if (ensemble%isotropy == IsotropySEMI_ISO) then

        scale_c(1) = (1.0_dp -                                                 &
             compress*dt*(press0 - 0.5_dp*(press(1)+press(2)))/tau_p)**ONE_THIRD
        scale_c(2) = (1.0_dp -                                                 &
             compress*dt*(press0 - 0.5_dp*(press(1)+press(2)))/tau_p)**ONE_THIRD
        scale_c(3) = (1.0_dp -                                                 &
             compress*dt*(press0 - press(3))/tau_p)**ONE_THIRD

      else if (ensemble%isotropy == IsotropyANISO) then

        scale_c(1:3) = (1.0_dp -                                               &
            compress*dt*(press0 - press(1:3))/tau_p)**ONE_THIRD

      else if (ensemble%isotropy == IsotropyXY_Fixed) then

        scale_c(1) = 1.0_dp
        scale_c(2) = 1.0_dp
        scale_c(3) = (1.0_dp -                                                 &
             compress*dt*(press0 - press(3))/tau_p)**ONE_THIRD

      end if

      do i = 1, ncell
        do ix = 1, natom(i)
          coord(1:3,ix,i) = scale_c(1:3) * coord(1:3,ix,i)
          trans(1:3,ix,i) = scale_c(1:3) * trans(1:3,ix,i)
        end do
      end do

      do i = ncell+1, ncell+nboundary
        do ix = 1, natom(i)
          trans(1:3,ix,i) = scale_c(1:3) * trans(1:3,ix,i)
        end do
      end do

      boundary%box_size_x = scale_c(1) * boundary%box_size_x_ref
      boundary%box_size_y = scale_c(2) * boundary%box_size_y_ref
      boundary%box_size_z = scale_c(3) * boundary%box_size_z_ref

      call bcast_boxsize(boundary%box_size_x, boundary%box_size_y, &
                         boundary%box_size_z)

      boundary%cell_size_x = boundary%box_size_x / real(boundary%num_cells_x,dp)
      boundary%cell_size_y = boundary%box_size_y / real(boundary%num_cells_y,dp)
      boundary%cell_size_z = boundary%box_size_z / real(boundary%num_cells_z,dp)

      domain%system_size(1) = boundary%box_size_x
      domain%system_size(2) = boundary%box_size_y
      domain%system_size(3) = boundary%box_size_z

    end if

    return

  end subroutine berendsen_leapfrog

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    gaussian_leapfrog
  !> @brief        control temperature using Gaussian thermostat
  !! @authors      TM, JJ
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine gaussian_leapfrog(dynamics, ensemble, domain, constraints, dynvars)

    ! formal arguments
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_domain),  target, intent(inout) :: domain
    type(s_constraints),     intent(inout) :: constraints
    type(s_dynvars),         intent(inout) :: dynvars

    ! local variables
    real(dp)                 :: dt, inv_dt, ekin, tempt, beta
    real(dp)                 :: factor1, factor2, factor
    real(dp)                 :: vel_tmp(1:3), temp0
    integer                  :: i, ix, j, k, l, num_degree

    real(dp),        pointer :: mass(:,:)
    real(dp),        pointer :: coord(:,:,:), coord_ref(:,:,:)
    real(dp),        pointer :: vel(:,:,:), vel_ref(:,:,:), force(:,:,:)
    integer,         pointer :: ncell, natom(:), nsolute(:), nwater(:)
    integer,         pointer :: water_list(:,:,:)


    ncell      => domain%num_cell_local
    natom      => domain%num_atom
    nsolute    => domain%num_solute
    nwater     => domain%num_water 
    water_list => domain%water_list
    coord      => domain%coord
    coord_ref  => domain%coord_ref
    vel        => domain%velocity
    vel_ref    => domain%velocity_ref
    force      => domain%force
    mass       => domain%mass

    dt         =  dynamics%timestep/AKMA_PS
    inv_dt     =  1.0_dp/dt
    temp0      =  ensemble%temperature
    num_degree =  domain%num_deg_freedom


    ! compute temperature(t+dt)
    !   v(t+dt) = 0.5*(v(t+3/2dt) + v(t+1/2dt))
    !
    ekin = 0.0_dp
    do i = 1, ncell
      do ix = 1, natom(i)
        vel_tmp(1:3) = 0.5_dp*(vel(1:3,ix,i) + vel_ref(1:3,ix,i))
        ekin = ekin + mass(ix,i)*(vel_tmp(1)*vel_tmp(1)+vel_tmp(2)*vel_tmp(2)+ &
                                  vel_tmp(3)*vel_tmp(3))
      end do
    end do
    ekin   = 0.5_dp * ekin

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, ekin, 1, mpi_real8, &
                       mpi_sum, mpi_comm_country, ierror)
#endif
    tempt  = 2.0_dp * ekin / (num_degree * KBOLTZ)

    ! calculate velocities(t + 3/2dt) and coordinates(t + 2dt)
    !   v(t+3/2dt) = v(t+1/2dt)*(2*beta-1) + beta*dt*F(t+dt)/m
    !   r(t+2dt) = r(t+dt) + dt*vel(t+3/2dt)
    !
    beta    = sqrt(temp0/tempt)
    factor  = beta * dt
    factor1 = 2.0_dp*beta - 1.0_dp

    do i = 1, ncell
      do ix = 1, natom(i)
        factor2  = factor/mass(ix,i)
        vel(1:3,ix,i) = factor1*vel_ref(1:3,ix,i) + factor2*force(1:3,ix,i)
        coord(1:3,ix,i) = coord_ref(1:3,ix,i) + dt*vel(1:3,ix,i)
      end do
      do l = 1, nwater(i)
        do k = 1, 3
          ix = water_list(k,l,i)
          factor2  = factor/mass(ix,i)
          vel(1:3,ix,i) = factor1*vel_ref(1:3,ix,i) + factor2*force(1:3,ix,i)
          coord(1:3,ix,i) = coord_ref(1:3,ix,i) + dt*vel(1:3,ix,i)
        end do
      end do
    end do


    ! SHAKE and SETTLE
    !   compute constrained coordinates(t + 2dt)
    !   update velocities(t + 3/2dt):
    !     vel(t+3/2dt) = (r(t+2dt) - r(t+dt))/dt
    !
    if (constraints%rigid_bond) then
      call compute_constraints(ConstraintModeLEAP, .true., dt, coord_ref, &
                               domain, constraints, coord, vel, &
                               dynvars%virial_const)

    end if

    return

  end subroutine gaussian_leapfrog

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    langevin_leapfrog_nvt
  !> @brief        control temperature using Langevin
  !! @authors      JJ, TM
  !! @param[in]    dynamics    : dynamics information
  !! @param[in]    ensemble    : ensemble information
  !! @param[inout] domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine langevin_leapfrog_nvt(dynamics, ensemble, domain, &
                                   constraints, dynvars)

    ! formal arguments
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_domain),  target, intent(inout) :: domain
    type(s_constraints),     intent(inout) :: constraints
    type(s_dynvars), target, intent(inout) :: dynvars

    ! local variables
    real(dp)                 :: dt, temp0, gamma_t, inv_dt
    real(dp)                 :: scale_v, scale_f, factor
    real(dp)                 :: sigma
    real(dp)                 :: v1, v2, rsq, grandom(3)
    integer                  :: j, jx, k, l, ncell

    real(dp),        pointer :: coord(:,:,:), force(:,:,:), coord_ref(:,:,:)
    real(dp),        pointer :: vel(:,:,:), viri(:,:), viri_const(:,:)
    real(dp),        pointer :: mass(:,:)
    integer,         pointer :: natom(:), nsolute(:), nwater(:)
    integer,         pointer :: water_list(:,:,:)


    natom      => domain%num_atom
    nsolute    => domain%num_solute
    nwater     => domain%num_water 
    water_list => domain%water_list
    mass       => domain%mass
    coord      => domain%coord
    coord_ref  => domain%coord_ref
    vel        => domain%velocity
    force      => domain%force
    viri       => dynvars%virial
    viri_const => dynvars%virial_const

    dt         =  dynamics%timestep/AKMA_PS
    inv_dt     =  1.0_dp/dt
    temp0      =  ensemble%temperature
    ncell      =  domain%num_cell_local
    gamma_t    =  ensemble%gamma_t*AKMA_PS


    ! add random force R(t+dt) to F(t+dt)
    !   R(t+dt) = sqrt(2gmKbT/dt)*Gauss(0,1)
    !
    factor = 2.0_dp*gamma_t*KBOLTZ*temp0/dt

    do j = 1, ncell
      do jx = 1, natom(j)

        sigma = sqrt(mass(jx,j) * factor)
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

        force(1:3,jx,j) = force(1:3,jx,j) + sigma * grandom(1:3)

      end do
    end do

    ! FEP
    ! synchronize single B with single A
    if (domain%fep_use) call sync_single_fep(domain, force)

    ! Langevin dynamics (Langevin thermostat)
    !   calculate velocities(t + 3/2dt)
    !   v(t+3/2dt) = scale_v*v(t+1/2dt) + scale_f*(F(t+dt)+R(t+dt))/m
    !   r(t+2dt)   = r(t+dt) + dt*v(t+3/2dt)
    !
    factor  = 1.0_dp + 0.5_dp * dt * gamma_t
    scale_v = (2.0_dp/factor) - 1.0_dp
    scale_f = dt/factor

    !$omp parallel do default(none)                                            &
    !$omp private(j, jx)                                                       &
    !$omp shared(ncell, natom, nsolute, nwater, water_list, dt, scale_v,       &
    !$omp        scale_f, vel, mass, force, coord)
    !
    do j = 1, ncell
      do jx = 1, nsolute(j)
        vel(1:3,jx,j) = scale_v*vel(1:3,jx,j)                                  &
                       + scale_f*force(1:3,jx,j)/mass(jx,j)
        coord(1:3,jx,j) = coord(1:3,jx,j) + dt*vel(1:3,jx,j)
      end do
      do l = 1, nwater(j)
        do k = 1, 3
          jx = water_list(k,l,j)
          vel(1:3,jx,j) = scale_v*vel(1:3,jx,j)                                  &
                         + scale_f*force(1:3,jx,j)/mass(jx,j)
          coord(1:3,jx,j) = coord(1:3,jx,j) + dt*vel(1:3,jx,j)
        end do
      end do
    end do
    !$omp end parallel do


    ! Bond constraint
    !   coord (unconstrained) is at t + 2dt, vel     is at t + 3/2dt
    !   coord_ref             is at t +  dt, vel_ref is at t + 1/2dt
    !   compute constrained coordinates(t+2dt) and constraint virial(t+dt)
    !   update velocities(t + 3/2dt):
    !     v(t+3/2dt) = (r(t+2dt) - r(t+dt))/dt
    !   add constraint virial(t + dt)
    !
    if (constraints%rigid_bond) then

      call compute_constraints(ConstraintModeLEAP, .true., dt, coord_ref, &
                               domain, constraints, coord, vel, &
                               dynvars%virial_const)

      !needed?
      viri(1:3,1:3) = viri(1:3,1:3) + viri_const(1:3,1:3)

    end if

    return

  end subroutine langevin_leapfrog_nvt

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    langevin_leapfrog_npt
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

  subroutine langevin_leapfrog_npt(dynamics, ensemble, domain, &
                                   constraints, boundary, dynvars)

    ! formal arguments
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_domain),  target, intent(inout) :: domain
    type(s_constraints),     intent(inout) :: constraints
    type(s_boundary),        intent(inout) :: boundary
    type(s_dynvars), target, intent(inout) :: dynvars

    ! local variables
    real(dp)                 :: dt, inv_dt, temp0, press0, gamma0, pressxy0
    real(dp)                 :: d_ndegf, kin(1:3), ekin
    real(dp)                 :: virial_sum(1:3)
    real(dp)                 :: volume, press(1:3)
    real(dp)                 :: pressxy, pressxyz
    real(dp)                 :: vel_tmp(1:3), coord_tmp(1:3)
    real(dp)                 :: scale_f(3), scale_v(3), fact(3), factor
    real(dp)                 :: bmoment_ref(3), bmoment2(3), scale_b(3)
    real(dp)                 :: gamma_t, gamma_p, pmass, pforce(3)
    real(dp)                 :: sigma
    real(dp)                 :: v1, v2, rsq, grandom(3)
    integer                  :: i, j, ij, jx, k, l, i_ndegf, ncell, nboundary

    real(dp),        pointer :: mass(:,:)
    real(dp),        pointer :: coord(:,:,:)
    real(dp),        pointer :: coord_ref(:,:,:), coord_old(:,:,:)
    real(dp),        pointer :: vel(:,:,:), vel_ref(:,:,:), force(:,:,:)
    real(dp),        pointer :: virial(:,:), viri_const(:,:)
    real(dp),        pointer :: bmoment(:)
    integer,         pointer :: natom(:), nsolute(:), nwater(:)
    integer,         pointer :: water_list(:,:,:)


    mass       => domain%mass
    natom      => domain%num_atom
    nsolute    => domain%num_solute
    nwater     => domain%num_water 
    water_list => domain%water_list
    coord      => domain%coord
    coord_ref  => domain%coord_ref
    vel        => domain%velocity
    vel_ref    => domain%velocity_ref
    force      => domain%force
    coord_old  => domain%coord_old
    virial     => dynvars%virial
    viri_const => dynvars%virial_const
    bmoment    => dynvars%barostat_momentum

    dt         =  dynamics%timestep/AKMA_PS
    inv_dt     =  1.0_dp/dt
    temp0      =  ensemble%temperature
    press0     =  ensemble%pressure * ATMOS_P
    gamma0     =  ensemble%gamma    * ATMOS_P*100.0_dp/1.01325_dp !TODO
    gamma_t    =  ensemble%gamma_t * AKMA_PS
    gamma_p    =  ensemble%gamma_p * AKMA_PS
    i_ndegf    =  domain%num_deg_freedom
    d_ndegf    =  real(i_ndegf,dp)
    ncell      =  domain%num_cell_local
    nboundary  =  domain%num_cell_boundary


    ! save box size(t+dt) and compute volume(t+dt)
    !
    boundary%box_size_x_ref = boundary%box_size_x
    boundary%box_size_y_ref = boundary%box_size_y
    boundary%box_size_z_ref = boundary%box_size_z
    volume = boundary%box_size_x * boundary%box_size_y * boundary%box_size_z

    bmoment_ref(1:3) = bmoment(1:3)

    ! Initial guess to calculate eta(t+dt)
    !   v(t+3/2dt) = v(t+1/2dt) + dt*F(t+dt)/m
    !   r(t+2dt)   = r(t+dt) + dt*v(t+3/2dt)
    !
    !$omp parallel do default(shared) private(j, jx)                      
    !
    do j = 1, ncell
      do jx = 1, nsolute(j)
        vel(1:3,jx,j)   = vel(1:3,jx,j) + dt*force(1:3,jx,j)/mass(jx,j)
        coord(1:3,jx,j) = coord(1:3,jx,j) + dt*vel(1:3,jx,j)
      end do
      do l = 1, nwater(j)
        do k = 1, 3
          jx = water_list(k,l,j)
          vel(1:3,jx,j)   = vel(1:3,jx,j) + dt*force(1:3,jx,j)/mass(jx,j)
          coord(1:3,jx,j) = coord(1:3,jx,j) + dt*vel(1:3,jx,j)
        end do
      end do
    end do
    !$omp end parallel do

    ! add random_force R(t+dt) to F(t+dt)
    !   R(t+dt) = sqrt(2gmKbT/dt)*Gauss(0,1)
    !
    factor   = 2.0_dp*gamma_t*KBOLTZ*temp0/dt

    do j = 1, ncell
      do jx = 1, natom(j)

        sigma = sqrt(mass(jx,j) * factor)
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

        force(1:3,jx,j) = force(1:3,jx,j) + sigma * grandom(1:3)

      end do
    end do

    ! FEP
    ! synchronize single B with single A
    if (domain%fep_use) call sync_single_fep(domain, force)

    ! calculate piston mass (pmass) and stochastic force (pforce)
    ! acting on the barostat
    !
    if (ensemble%isotropy == IsotropyISO) then

      pmass = real(i_ndegf+3,dp)*KBOLTZ*temp0 / (2.0_dp*PI*gamma_p)**2
      sigma = sqrt(2.0_dp*gamma_p*pmass*KBOLTZ*temp0/dt)

      pforce(1) = sigma * random_get_gauss()
      pforce(2) = pforce(1)
      pforce(3) = pforce(1)

    else if (ensemble%isotropy == IsotropySEMI_ISO) then

      pmass = real(i_ndegf+3,dp)*KBOLTZ*temp0 / (3.0_dp*(2.0_dp*PI*gamma_p)**2)
      sigma = sqrt(2.0_dp*gamma_p*pmass*KBOLTZ*temp0/dt)

      pforce(1) = sigma * random_get_gauss()
      pforce(2) = pforce(1)
      pforce(3) = sigma * random_get_gauss()

    else if (ensemble%isotropy == IsotropyANISO) then

      pmass = real(i_ndegf+3,dp)*KBOLTZ*temp0 / (3.0_dp*(2.0_dp*PI*gamma_p)**2)
      sigma = sqrt(2.0_dp*gamma_p*pmass*KBOLTZ*temp0/dt)

      pforce(1) = sigma * random_get_gauss()
      pforce(2) = sigma * random_get_gauss()
      pforce(3) = sigma * random_get_gauss()

    else if (ensemble%isotropy == IsotropyXY_FIXED) then

      pmass = real(i_ndegf+1,dp)*KBOLTZ*temp0 / (3.0_dp*(2.0_dp*PI*gamma_p)**2)
      sigma = sqrt(2.0_dp*gamma_p*pmass*KBOLTZ*temp0/dt)

      pforce(1) = 0.0_dp
      pforce(2) = 0.0_dp
      pforce(3) = sigma * random_get_gauss()

    end if

#ifdef HAVE_MPI_GENESIS
    call mpi_bcast(pforce, 3, mpi_real8, 0, mpi_comm_country, ierror)
#endif

    ! iteration
    !   Since leapfrog algorithm cannot directly compute vel(t+dt),
    !   iteration scheme is required to obtain T(t+dt) and P(t+dt).
    !
    do i = 1, 8

      ! compute kinetic energy(t+dt)
      !   v(t+dt) = 0.5*(v(t+3/2dt) + v(t+1/2dt))
      !
      kin(1:3) = 0.0_dp

      do j = 1, ncell
        do jx = 1, natom(j)

          ! FEP
          ! skip singleB to avoid duplication
          if (domain%fep_use) then
            if (domain%fepgrp(jx,j) == 2) cycle
          end if

          kin(1:3) = kin(1:3) + mass(jx,j) * vel(1:3,jx,j)*vel(1:3,jx,j)
          kin(1:3) = kin(1:3) + mass(jx,j)*vel_ref(1:3,jx,j)*vel_ref(1:3,jx,j)
        end do
      end do
      kin(1:3) = kin(1:3) / 2.0_dp

!     ! scale factor for harmonic oscillators (R.Pastor et al., Mol.Phys. 1988)
!     ! it is not used in the barostat (Jaewoon Jung)
!     if (i /= 1) then
!       kin(1:3) = kin(1:3) * (1.0_dp + 0.5_dp*gamma_t*dt)
!     end if

      ekin   = 0.5_dp * (kin(1) + kin(2) + kin(3))
      call reduce_pres(kin, ekin, virial(1,1), virial(2,2), virial(3,3),    &
                       virial_sum)

      ! compute pressure(t+dt) in the unit of kcal/mol*A3
      !
      press(1:3) = (kin(1:3) + virial_sum(1:3)) / volume
      pressxyz   = (press(1) + press(2) + press(3))/3.0_dp
      pressxy    = (press(1) + press(2))/2.0_dp

      ! compute thermostat and barostat parameters
      !   for isotropic systems
      !     eta(t+3/2dt) = eta(t+1/2dt) + dt[dim*V(t+dt)*(P(t+dt)-Pext)
      !                                     + dim*2Ekin(t+dt)/Nf +Rp]/pmass
      !     eta(t+3/2dt) = eps(-gamma_p*dt)*eta(t+3/2dt)
      !     eta(t+dt)    = 0.5*(eta(t+3/2dt) + eta(t+1/2dt))
      !     factor       = 1 + [gamma_t+(1+dim/Nf)*eta(t+dt)]dt/2
      !
      if (ensemble%isotropy == IsotropyISO) then

        bmoment(1) = bmoment_ref(1) + dt*(3.0_dp*volume*(pressxyz - press0)    &
                                    + 6.0_dp*ekin/d_ndegf + pforce(1))/pmass
        bmoment(2) = bmoment(1)
        bmoment(3) = bmoment(1)

        bmoment(1:3)  = exp(-gamma_p*dt)*bmoment(1:3)
        bmoment2(1:3) = 0.5_dp*(bmoment(1:3) + bmoment_ref(1:3))

        fact(1:3) = 1.0_dp+                                                    &
           (gamma_t+bmoment2(1:3)*(1.0_dp+3.0_dp/d_ndegf))*0.5_dp*dt

      else if (ensemble%isotropy == IsotropySEMI_ISO) then

        if (ensemble%ensemble == EnsembleNPT) then

          bmoment(1:2) = bmoment_ref(1:2) + dt*(volume*(pressxy - press0)      &
                                      + 2.0_dp*ekin/d_ndegf + pforce(1:2))/pmass
          bmoment(3)   = bmoment_ref(3) + dt*(volume*(press(3)  - press0)      &
                                      + 2.0_dp*ekin/d_ndegf + pforce(3))/pmass

        else if (ensemble%ensemble == EnsembleNPgT) then

          pressxy0 = press0 - gamma0 / boundary%box_size_z

          bmoment(1:2) = bmoment_ref(1:2) + dt*(volume*(pressxy - pressxy0)    &
                                      + 2.0_dp*ekin/d_ndegf + pforce(1:2))/pmass
          bmoment(3) = bmoment_ref(3) + dt*(volume*(press(3)  - press0)        &
                                      + 2.0_dp*ekin/d_ndegf + pforce(3))/pmass

        end if

        bmoment(1:3)  = exp(-gamma_p*dt)*bmoment(1:3)
        bmoment2(1:3) = 0.5_dp*(bmoment(1:3) + bmoment_ref(1:3))

        factor    = bmoment2(1) + bmoment2(2) + bmoment2(3)
        fact(1:3) = 1.0_dp+(gamma_t+bmoment2(1:3)+factor/d_ndegf)*0.5_dp*dt

      else if (ensemble%isotropy == IsotropyANISO) then

        if (ensemble%ensemble == EnsembleNPT) then

          bmoment(1:3) = bmoment_ref(1:3) + dt*(volume*(press(1:3) - press0)   &
                                      + 2.0_dp*ekin/d_ndegf + pforce(1:3))/pmass

        else if (ensemble%ensemble == EnsembleNPgT) then

          pressxy0 = press0 - gamma0 / boundary%box_size_z

          bmoment(1:2) = bmoment_ref(1:2) + dt*(volume*(press(1:2) - pressxy0) &
                                      + 2.0_dp*ekin/d_ndegf + pforce(1:2))/pmass
          bmoment(3) = bmoment_ref(3) + dt*(volume*(press(3) - press0)         &
                                      + 2.0_dp*ekin/d_ndegf + pforce(3))/pmass

        end if

        bmoment(1:3)  = exp(-gamma_p*dt)*bmoment(1:3)
        bmoment2(1:3) = 0.5_dp*(bmoment(1:3) + bmoment_ref(1:3))

        factor    = bmoment2(1) + bmoment2(2) + bmoment2(3)
        fact(1:3) = 1.0_dp+(gamma_t+bmoment2(1:3)+factor/d_ndegf)*0.5_dp*dt

      else if (ensemble%isotropy == IsotropyXY_Fixed) then

        bmoment(1) = 0.0_dp
        bmoment(2) = 0.0_dp
        bmoment(3) = bmoment_ref(3) + dt*(volume*(press(3) - press0)           &
                                    + 2.0_dp*ekin/d_ndegf + pforce(3))/pmass

        bmoment(1:3)  = exp(-gamma_p*dt)*bmoment(1:3)
        bmoment2(1:3) = 0.5_dp*(bmoment(1:3) + bmoment_ref(1:3))

        factor    = bmoment2(1) + bmoment2(2) + bmoment2(3)
        fact(1:3) = 1.0_dp+(gamma_t+bmoment2(1:3)+factor/d_ndegf)*0.5_dp*dt

      end if

      scale_v(1:3) = 2.0_dp/fact(1:3) - 1.0_dp
      scale_f(1:3) = dt/fact(1:3)


      ! Langevin dynamics
      !   calculate velocities(t + 3/2dt) and coordinates(t + 2dt)
      !   v(t+3/2dt) = scale_v*v(t+1/2dt) + scale_f*(F(t+dt)+R(t+dt))/m
      !   r(t+3/2dt) = 0.5*(r(t+2dt) + r(t+dt))
      !   r(t+2dt)   = r(t+dt) + dt*(v(t+3/2dt) + eta(t+3/2dt)*r(t+3/2dt)))
      !
      do j = 1, ncell
        do jx = 1, nsolute(j)
          vel(1:3,jx,j) = scale_v(1:3)*vel_ref(1:3,jx,j) + &
                          scale_f(1:3)*force(1:3,jx,j)/mass(jx,j)
          coord_tmp(1:3) = 0.5_dp*(coord_ref(1:3,jx,j) + coord(1:3,jx,j))
          coord(1:3,jx,j) = coord_ref(1:3,jx,j) + dt*(vel(1:3,jx,j) + &
                            bmoment(1:3)*coord_tmp(1:3))
        end do
        do l = 1, nwater(j)
          do k = 1, 3
            jx = water_list(k,l,j)
            vel(1:3,jx,j) = scale_v(1:3)*vel_ref(1:3,jx,j) + &
                            scale_f(1:3)*force(1:3,jx,j)/mass(jx,j)
            coord_tmp(1:3) = 0.5_dp*(coord_ref(1:3,jx,j) + coord(1:3,jx,j))
            coord(1:3,jx,j) = coord_ref(1:3,jx,j) + dt*(vel(1:3,jx,j) + &
                              bmoment(1:3)*coord_tmp(1:3)) 
          end do
        end do
      end do

      ! Bond constraint
      !   store unconstrained coord(t + 2dt)
      !   compute constrained coordinates(t + 2dt)
      !   add contribution of constraints to virial
      !   correct velocities(t + 3/2dt) and forces(t + dt)
      !
      if (constraints%rigid_bond) then

        do j = 1, ncell
          do jx = 1, natom(j)
            coord_old(1:3,jx,j) = coord(1:3,jx,j)
          end do
        end do

        call compute_constraints(ConstraintModeLEAP, .false., dt, coord_ref, &
                                 domain, constraints, coord, vel, &
                                 viri_const)

        virial(1:3,1:3) = virial(1:3,1:3) + viri_const(1:3,1:3)

        !$omp parallel do default(shared) private(j, jx, vel_tmp)      
        !
        do j = 1, ncell
          do jx =1, natom(j)
            vel_tmp(1:3)    = (coord(1:3,jx,j) - coord_old(1:3,jx,j))*inv_dt
            vel(1:3,jx,j)   = vel(1:3,jx,j) + vel_tmp(1:3)
            force(1:3,jx,j) = force(1:3,jx,j) + mass(jx,j)*vel_tmp(1:3)*inv_dt
          end do
        end do
        !$omp end parallel do

      end if

    end do


    ! compute box size(t+2dt)
    !   size(t+2dt) = exp[eta(t+3/2dt)*dt] * size(t+dt)
    !

    scale_b(1:3) = exp(bmoment(1:3)*dt)
    boundary%box_size_x = scale_b(1) * boundary%box_size_x_ref
    boundary%box_size_y = scale_b(2) * boundary%box_size_y_ref
    boundary%box_size_z = scale_b(3) * boundary%box_size_z_ref

    call bcast_boxsize(boundary%box_size_x, &
                       boundary%box_size_y, &
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

    return

  end subroutine langevin_leapfrog_npt

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    simulated_annealing_leapfrog
  !> @brief        change target temperature linearly
  !! @authors      TM
  !! @param[in]    dynamics: dynamics information
  !! @param[out]   ensemble: ensemble information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine simulated_annealing_leapfrog(dynamics, ensemble)
  
    ! formal arguments
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_ensemble),        intent(inout) :: ensemble
    
    ! local variable
    real(dp)                 :: old_temperature
    

    if (.not. dynamics%annealing) return

    old_temperature      = ensemble%temperature
    ensemble%temperature = ensemble%temperature + dynamics%dtemperature
    
    if (main_rank) then
      write(MsgOut,'(A,F10.3,A,F10.3)')                              &
            'Simulated_Annealing_Leapfrog> Anneal temperature from', &
            old_temperature, ' to ', ensemble%temperature

      if (ensemble%temperature < 0.0_dp) &
        call error_msg( &
        'Simulated_Annealing_Leapfrog> Error: Temperature is less than 0 K')

    end if

    return

  end subroutine simulated_annealing_leapfrog

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_gamd_leapfrog
  !> @brief        update GaMD parameters
  !! @authors      HO
  !! @param[inout] output      : output information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] dynvars     : dynamic variables information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
 
  subroutine update_gamd_leapfrog(output, enefunc, dynvars)
    ! formal arguments
    type(s_output),           intent(inout) :: output
    type(s_enefunc),  target, intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars

    ! local variables
    type(s_enefunc_gamd), pointer    :: gamd

    gamd => enefunc%gamd

    ! Compute and update gamd parameters
    !
    call setup_enefunc_gamd(gamd)

    ! output gamd parameters
    !
    call output_gamd(output, dynvars, enefunc)

    return

  end subroutine update_gamd_leapfrog

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    reduce_pres
  !> @brief        reduce pres
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine reduce_pres(val1, val2, val3, val4, val5, val6)

    ! formal arguments
    real(dp),                intent(inout) :: val1(1:3), val2
    real(dp),                intent(inout) :: val3, val4, val5
    real(dp),                intent(inout) :: val6(1:3)

#ifdef HAVE_MPI_GENESIS

    ! local variables
    real(dp)                 :: before_reduce(7), after_reduce(7)

    before_reduce(1:3)  = val1(1:3)
    before_reduce(4)    = val2
    before_reduce(5)    = val3
    before_reduce(6)    = val4
    before_reduce(7)    = val5

    call mpi_allreduce(before_reduce, after_reduce, 7, mpi_real8, &
                       mpi_sum, mpi_comm_country, ierror)

    val1(1:3) = after_reduce(1:3)
    val2      = after_reduce(4)
    val6(1:3) = after_reduce(5:7)

#else

    val6(1) = val3
    val6(2) = val4
    val6(3) = val5

#endif

    return

  end subroutine reduce_pres

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    bcast_boxsize
  !> @brief        bcast box size
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine bcast_boxsize(val1, val2, val3)

    ! formal arguments
    real(dp),                intent(inout) :: val1, val2, val3

#ifdef HAVE_MPI_GENESIS

    ! local variables
    real(dp)                 :: list(3)


    list(1) = val1
    list(2) = val2
    list(3) = val3

    call mpi_bcast(list, 3, mpi_real8, 0, mpi_comm_country, ierror)

    val1 = list(1)
    val2 = list(2)
    val3 = list(3)

#endif

    return

  end subroutine bcast_boxsize

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    leapfrog_dynamics_fep
  !> @brief        leapfrog integrator with domain decomposition
  !!               for FEP calculations
  !! @authors      NK
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

  subroutine leapfrog_dynamics_fep(output, domain, enefunc, dynvars, dynamics, &
                               pairlist, boundary, constraints, ensemble,      &
                               comm, remd, alchemy)

    ! formal arguments
    type(s_output),            intent(inout) :: output
    type(s_domain),    target, intent(inout) :: domain
    type(s_enefunc),           intent(inout) :: enefunc
    type(s_dynvars),   target, intent(inout) :: dynvars
    type(s_dynamics),          intent(inout) :: dynamics
    type(s_pairlist),          intent(inout) :: pairlist
    type(s_boundary),          intent(inout) :: boundary
    type(s_constraints),       intent(inout) :: constraints
    type(s_ensemble),          intent(inout) :: ensemble
    type(s_comm),              intent(inout) :: comm   
    type(s_remd),              intent(inout) :: remd
    type(s_alchemy), optional, intent(inout) :: alchemy

    ! local variables
    real(dp)                 :: simtim, dt, temperature, factor
    real(dp)                 :: energy(18), temp(18), temp_prev(18)
    integer                  :: i, j, k, jx, nsteps, ekin, iseed
    integer                  :: istart, iend
    logical                  :: npt

    real(dp),        pointer :: coord(:,:,:), coord_ref(:,:,:)
    real(dp),        pointer :: coord_old(:,:,:)
    real(wp),        pointer :: coord_pbc(:,:,:)
    real(dp),        pointer :: vel(:,:,:), vel_ref(:,:,:)
    real(dp),        pointer :: force(:,:,:), force_long(:,:,:)
    real(dp),        pointer :: mass(:,:)
    real(wp),        pointer :: force_omp(:,:,:,:)
    real(wp),        pointer :: force_pbc(:,:,:,:)
    real(dp),        pointer :: virial_cell(:,:)
    real(dp),        pointer :: virial(:,:)
    integer,         pointer :: ncell, nboundary, natom(:)

    ! fep/rest
    integer :: parmsetid, parmsetid_forward, parmsetid_backward
    integer :: replicaid
    type(s_energy), save :: frenergy


    ncell       => domain%num_cell_local
    nboundary   => domain%num_cell_boundary
    natom       => domain%num_atom
    mass        => domain%mass
    coord       => domain%coord
    coord_ref   => domain%coord_ref
    coord_old   => domain%coord_old
    coord_pbc   => domain%translated
    force       => domain%force
    force_long  => domain%force_long
    force_omp   => domain%force_omp
    force_pbc   => domain%force_pbc
    virial_cell => domain%virial_cellpair
    vel         => domain%velocity
    vel_ref     => domain%velocity_ref
    virial      => dynvars%virial

    temperature =  ensemble%temperature
    npt         =  ensemble%use_barostat
    nsteps      =  dynamics%nsteps
    istart      =  dynamics%istart_step
    iend        =  dynamics%iend_step
    dt          =  dynamics%timestep/AKMA_PS
    simtim      =  dynamics%initial_time
    iseed       =  dynamics%iseed_init_velocity

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

    if (abs(dynamics%initial_rmsd) < 0.001_wp)  &
      dynamics%initial_rmsd = dynvars%energy%rmsd
    if (dynamics%target_md) enefunc%rmsd_force = 1 / (dt*dt)

    ! Check restart
    !
    if (.not. dynamics%restart .and. .not. pio_restart) then

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

      call initial_leapfrog_fep(npt, output, enefunc, dynamics,            &
                            pairlist, boundary, ensemble, constraints, &
                            domain, dynvars, comm)
    end if


    ! Main loop 
    !   coord is at 0 +  dt and vel is at 0 + 1/2dt, if restart off
    !   coord is at t + 2dt and vel is at t + 3/2dt, if restart on
    !
    do i = istart, iend

      call timer(TimerIntegrator, TimerOn)

      simtim = simtim + dynamics%timestep
      dynvars%time = simtim
      dynvars%step = i
      if (dynamics%target_md .or. dynamics%steered_md) &
        enefunc%rmsd_target = dynamics%initial_rmsd &
                            + (dynamics%final_rmsd-dynamics%initial_rmsd) &
                             *real(dynvars%step,wp)/real(nsteps,wp)
      enefunc%rpath_sum_mf_flag = enefunc%rpath_flag

      ! send the coordinate data
      !
      call timer(TimerComm1, TimerOn)

      if (constraints%tip4 .or. enefunc%table%tip4) &
      call decide_dummy(domain, constraints, coord)

      call communicate_coor(domain, comm)

      call timer(TimerComm1, TimerOff)

      ! Save coordinates(t + dt) and velocities(t + 1/2dt)
      !
      !$omp parallel do default(shared) private(j,jx)
      !
      do j = 1, ncell+nboundary
        do jx = 1, natom(j)
          coord_ref(1:3,jx,j) = coord(1:3,jx,j)
          vel_ref  (1:3,jx,j) = vel  (1:3,jx,j)
        end do
      end do
      !$omp end parallel do
      call timer(TimerIntegrator, TimerOff)

      ! Compute energy(t + dt), force(t + dt), and internal virial(t + dt)
      !
      call compute_energy_fep(domain, enefunc, pairlist, boundary, coord, &
                          npt, .false., mod(i,dynamics%eneout_period)==0, &
                          .true.,                  &
                          enefunc%nonb_limiter,    &
                          dynvars%energy,          &
                          coord_pbc,               &
                          force,                   &
                          force_long,              &
                          force_omp,               &
                          force_pbc,               &
                          virial_cell,             &
                          dynvars%virial,          &
                          dynvars%virial_long, &
                          dynvars%virial_extern)

      ! FEP: compute and output energy differnces between adjacent states
      if ((domain%fep_use).and.(present(alchemy))) then
        if (dynamics%fepout_period > 0) then
          if ((mod(i,dynamics%fepout_period) == 0) .and. &
            (mod(i,dynamics%eneout_period) == 0)) then
            if (i > dynamics%equilsteps) then
              call compute_fep_energy(domain, enefunc, dynvars, pairlist, &
                ensemble, boundary, alchemy, remd)
              call output_fep_energy(output, enefunc, dynvars)
            end if
          end if
        end if
      end if

      call timer(TimerIntegrator, TimerOn)
      call timer(TimerComm2, TimerOn)
                          
      call communicate_force(domain, comm, force)
      if (constraints%tip4 .or. enefunc%table%tip4) &
        call water_force_redistribution(constraints, domain, force, virial)

      call timer(TimerComm2, TimerOff)

      ! FEP: merge forces of single A and single B
      call add_single_fep(domain, force)

      if (ensemble%tpcontrol /= TpcontrolLangevin) then

        ! Newtonian dynamics
        !   v(t+3/2dt) = v(t+1/2dt) + dt*F(t+dt)/m
        !   r(t+2dt) = r(t+dt) + dt*v(t+3/2dt)
        !
        !$omp parallel do default(shared) private(j, jx, factor)
        do j = 1, ncell
          do jx = 1, natom(j)
            factor = dt/mass(jx,j)
            vel  (1:3,jx,j) = vel  (1:3,jx,j) + factor*force(1:3,jx,j)
            coord(1:3,jx,j) = coord(1:3,jx,j) + dt*vel(1:3,jx,j)
          end do
        end do
        !$omp end parallel do


        ! Bond constraint
        !   coord (unconstrained) is at t + 2dt, vel     is at t + 3/2dt
        !   coord_ref             is at t +  dt, vel_ref is at t + 1/2dt
        !   compute constrained coordinates(t+2dt) and constraint virial(t+dt)
        !   update velocities(t + 3/2dt):
        !     v(t+3/2dt) = (r(t+2dt) - r(t+dt))/dt
        !   add constraint virial(t + dt)
        !
        if (constraints%rigid_bond) then
          call compute_constraints(ConstraintModeLEAP, .true., dt, coord_ref, &
                                   domain, constraints, coord, vel,  &
                                   dynvars%virial_const)
          dynvars%virial(1:3,1:3) = dynvars%virial(1:3,1:3) &
                                  + dynvars%virial_const(1:3,1:3)
        end if

        ! Control temperature and pressure
        !   coord     is at t + 2dt, vel     is at t + 3/2dt
        !   coord_ref is at t +  dt, vel_ref is at t + 1/2dt
        !   scale velocities(t + 3/2dt) and coordinates(t + 2dt)
        !

        if (ensemble%ensemble /= EnsembleNVE) then

          call error_msg('Leapfrog_dynamics_fep> Only Langevin is allowed in FEP')

        end if

      else

        ! Langevin dynamics
        !
        if (ensemble%ensemble == EnsembleNVT) then

          call langevin_leapfrog_nvt(dynamics, ensemble, domain, &
                                     constraints, dynvars)

        else if (ensemble%ensemble == EnsembleNPT  .or. &
                 ensemble%ensemble == EnsembleNPAT .or. &
                 ensemble%ensemble == EnsembleNPgT) then

          call langevin_leapfrog_npt(dynamics, ensemble, domain, &
                                     constraints, boundary, dynvars)

        end if

      end if

      call timer(TimerIntegrator, TimerOff)


      if (dynamics%stoptr_period > 0) then

        if (mod(i,dynamics%stoptr_period) == 0) then

          do j = 1, ncell
            do jx = 1, natom(j)
              coord_old(1:3,jx,j) = 0.5_dp*(coord(1:3,jx,j)+coord_ref(1:3,jx,j))
            end do
          end do
          call stop_trans_rotation_fep(ncell, natom,                  &
                                   dynamics%stop_com_translation, &
                                   dynamics%stop_com_rotation,    &
                                   mass, coord_old, vel, domain)

        end if

      end if

      ! ------------------------------------------------------------------------
      ! OUTPUT energy, trajectory, and restart file
      !   coord     is at t + 2dt, vel      is at t + 3/2dt
      !   coord_ref is at t +  dt, vel_ref  is at t +    dt
      !   volume    is at t +  dt, box_size is at t +   2dt

      call random_push_stock

      ! output dynvars(t + dt)
      !
      call output_md(output, dynamics, boundary, pairlist, &
                     ensemble, dynvars, domain, enefunc, remd, alchemy)

      ! update interaction
      !
      if (dynamics%nbupdate_period > 0) &
        call domain_interaction_update_md_fep(i, dynamics, domain, enefunc, &
                                          pairlist, boundary, constraints, comm)

      ! Simulated annealing
      !
      if (dynamics%anneal_period > 0) then
        if (mod(i,dynamics%anneal_period) == 0) then

          call simulated_annealing_leapfrog(dynamics, ensemble)

        end if
      end if

      ! output parallel I/O restart
      !
      call output_prst_md(output, enefunc, dynamics, boundary, &
                                  dynvars, domain, constraints)
      
    end do

    return

  end subroutine leapfrog_dynamics_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    initial_leapfrog_fep
  !> @brief        compute the first step (0+dt) for FEP calculations
  !! @authors      NK
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

  subroutine initial_leapfrog_fep(npt, output, enefunc, dynamics, pairlist, &
                              boundary, ensemble, constraints, domain,  &
                              dynvars, comm)

    ! formal arguments
    logical,                 intent(in)    :: npt
    type(s_output),          intent(in)    :: output
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_dynamics),        intent(inout) :: dynamics
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    type(s_ensemble),        intent(in)    :: ensemble
    type(s_constraints),     intent(inout) :: constraints
    type(s_domain),  target, intent(inout) :: domain
    type(s_dynvars), target, intent(inout) :: dynvars
    type(s_comm),            intent(inout) :: comm

    ! local variables
    real(dp)                 :: factor, temperature, energy(18)
    real(dp)                 :: temp(18),ekin1, ekin2
    real(dp)                 :: imass, simtim, dt, friction
    integer                  :: i, ix, j, jx, ncell, k, l

    real(dp),        pointer :: coord(:,:,:), coord_ref(:,:,:)
    real(wp),        pointer :: coord_pbc(:,:,:)
    real(dp),        pointer :: vel(:,:,:), vel_ref(:,:,:)
    real(dp),        pointer :: force(:,:,:), force_long(:,:,:)
    real(dp),        pointer :: mass(:,:)
    real(wp),        pointer :: force_omp(:,:,:,:)
    real(wp),        pointer :: force_pbc(:,:,:,:)
    real(dp),        pointer :: virial_cell(:,:), virial(:,:)
    integer,         pointer :: natom(:), nsolute(:), nwater(:)
    integer,         pointer :: water_list(:,:,:)


    natom        => domain%num_atom
    nsolute      => domain%num_solute
    nwater       => domain%num_water 
    water_list   => domain%water_list
    coord        => domain%coord
    coord_ref    => domain%coord_ref
    coord_pbc    => domain%translated
    force        => domain%force
    force_long   => domain%force_long
    force_omp    => domain%force_omp
    force_pbc    => domain%force_pbc
    virial_cell  => domain%virial_cellpair
    vel          => domain%velocity
    vel_ref      => domain%velocity_ref
    mass         => domain%mass
    virial       => dynvars%virial

    ncell        =  domain%num_cell_local
    dt           =  dynamics%timestep/AKMA_PS
    simtim       =  dynamics%initial_time
    temperature  =  ensemble%temperature
    friction     =  ensemble%gamma_t * AKMA_PS

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
                               domain, constraints, coord, vel,  &
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

    call compute_energy_fep(domain, enefunc, pairlist, boundary, coord,  &
                        npt, .false., .true.,  &
                        .true.,                &
                        enefunc%nonb_limiter,  &
                        dynvars%energy,        &
                        coord_pbc,             &
                        force,                 &
                        force_long,            &
                        force_omp,             &
                        force_pbc,             &
                        virial_cell,           &
                        dynvars%virial,        &
                        dynvars%virial_long,   &
                        dynvars%virial_extern)

    ! get forces by communication
    !
    call communicate_force(domain, comm, force)

    if (constraints%tip4 .or. enefunc%table%tip4) &
      call water_force_redistribution(constraints, domain, force, virial)

    ! FEP
    ! merge forces of single A and single B
    call add_single_fep(domain, force)

    ! Calculate velocities(0 - dt/2)
    !
    if (.not. constraints%rigid_bond) then

     do i = 1, ncell
       do ix = 1, nsolute(i)
         imass = 1.0_wp/mass(ix,i)
         vel(1:3,ix,i) = vel(1:3,ix,i) - 0.5_wp*dt*force(1:3,ix,i)*imass
       end do
       do l = 1, nwater(i)
         do k = 1, 3
           ix = water_list(k,l,i)
           imass = 1.0_wp/mass(ix,i)
           vel(1:3,ix,i) = vel(1:3,ix,i) - 0.5_wp*dt*force(1:3,ix,i)*imass
         end do
       end do
     end do

    else

      ! calculate velocities(0 + dt/2) and coordinates(0 + dt)
      ! update coordinates(0 + dt) and velocities(0 + dt/2)
      !
      do i = 1, ncell
        do ix = 1, nsolute(i)
          imass = 1.0_wp/mass(ix,i)
          vel(1:3,ix,i) = vel_ref(1:3,ix,i) + 0.5_wp*dt*force(1:3,ix,i)*imass
          coord(1:3,ix,i) = coord_ref(1:3,ix,i) + dt*vel(1:3,ix,i)
        end do
        do l = 1, nwater(i)
          do k = 1, 3
            ix = water_list(k,l,i)
            imass = 1.0_wp/mass(ix,i)
            vel(1:3,ix,i) = vel_ref(1:3,ix,i) + 0.5_wp*dt*force(1:3,ix,i)*imass
            coord(1:3,ix,i) = coord_ref(1:3,ix,i) + dt*vel(1:3,ix,i)
          end do
        end do  
      end do

      call compute_constraints(ConstraintModeLEAP, .true., dt, coord_ref, &
                               domain, constraints, coord, vel,           &
                               dynvars%virial_const)

      ! calculate velocities(0 - dt/2) and coordinates(0 - dt)
      ! update coordinates(0 - dt) and velocities(0 - dt/2)
      !
      do i = 1, ncell
        do ix = 1, nsolute(i)
          imass = 1.0_wp/mass(ix,i)
          vel_ref(1:3,ix,i) = vel_ref(1:3,ix,i)                            &
                             - 0.5_wp*dt*force(1:3,ix,i)*imass
          coord(1:3,ix,i) = coord_ref(1:3,ix,i) - dt*vel_ref(1:3,ix,i)
        end do
        do l = 1, nwater(i)
          do k = 1, 3
            ix = water_list(k,l,i)
            imass = 1.0_wp/mass(ix,i)
            vel_ref(1:3,ix,i) = vel_ref(1:3,ix,i)                            &
                               - 0.5_wp*dt*force(1:3,ix,i)*imass
            coord(1:3,ix,i) = coord_ref(1:3,ix,i) - dt*vel_ref(1:3,ix,i)
          end do
        end do
      end do

      call compute_constraints(ConstraintModeLEAP, .true., -dt, coord_ref, &
                               domain, constraints, coord, vel_ref,        &
                               dynvars%virial_const)

      ! vel <= vel(0 - dt/2) and coord <= coord(0)
      !
      do i = 1, ncell
        do ix = 1, natom(i)
          vel(1:3,ix,i)   = vel_ref(1:3,ix,i)
          coord(1:3,ix,i) = coord_ref(1:3,ix,i)
        end do
      end do

    end if

    ! first step dynamics
    !   Calculate coordinates(0 + dt) and velocities(0 + 1/2dt)
    !   Note: this routine corresponds to i = 0 in the main loop
    !

    ! coord_ref <= coord(0) and vel_ref <= vel(0 - 1/2dt)
    !
    do i = 1, ncell
      do ix = 1, natom(i)
        coord_ref(1:3,ix,i) = coord(1:3,ix,i)
        vel_ref(1:3,ix,i)   = vel(1:3,ix,i)
      end do
    end do

    ! calculate energy(0) and forces(0)
    !
    if (constraints%tip4 .or. enefunc%table%tip4) &
    call decide_dummy(domain, constraints, coord)

    call communicate_coor(domain, comm)

    call compute_energy_fep(domain, enefunc, pairlist, boundary, coord,  &
                        npt, .false., .true.,  &
                        .true.,                &
                        enefunc%nonb_limiter,  &
                        dynvars%energy,        &
                        coord_pbc,             &
                        force,                 &
                        force_long,            &
                        force_omp,             &
                        force_pbc,             &
                        virial_cell,           &
                        dynvars%virial,        &
                        dynvars%virial_long,   &
                        dynvars%virial_extern)

    call communicate_force(domain, comm, force)

    if (constraints%tip4 .or. enefunc%table%tip4) &
      call water_force_redistribution(constraints, domain, force, virial)

    ! FEP
    ! merge forces of single A and single B
    call add_single_fep(domain, force)

    ! Newtonian dynamics
    !   v(0+1/2dt) = v(0-1/2dt) + dt*F(0)/m
    !   r(0+dt) = r(0) + dt*v(0+1/2dt)
    !
    do i = 1, ncell
      do ix = 1, nsolute(i)
        factor = dt/mass(ix,i)
        vel(1:3,ix,i) = vel(1:3,ix,i) + factor*force(1:3,ix,i)
        coord(1:3,ix,i) = coord(1:3,ix,i) + dt*vel(1:3,ix,i)
      end do
      do l = 1, nwater(i)
        do k = 1, 3
          ix = water_list(k,l,i)
          factor = dt/mass(ix,i)
          vel(1:3,ix,i) = vel(1:3,ix,i) + factor*force(1:3,ix,i)
          coord(1:3,ix,i) = coord(1:3,ix,i) + dt*vel(1:3,ix,i)
        end do
      end do
    end do


    ! SHAKE and SETTLE
    !   calculate constrained coordinates(0 + dt) and constraint virial(0)
    !   update velocities(0 + 1/2dt):
    !     v(0+1/2dt) = (r(0+dt) - r(0))/dt
    !
    if (constraints%rigid_bond) then
      call compute_constraints(ConstraintModeLEAP, .true., dt, coord_ref, &
                               domain, constraints, coord, vel, &
                               dynvars%virial_const)

      dynvars%virial(1:3,1:3) = &
        dynvars%virial(1:3,1:3) + dynvars%virial_const(1:3,1:3)

    end if

    ! output dynvars(0)
    !
    call compute_dynvars(enefunc, dynamics, boundary, ensemble, domain, &
                             dynvars)

    call output_dynvars(output, enefunc, dynvars, ensemble, boundary)


    dynamics%restart = .true.

    ! at this point
    !   coord is at 0 + dt, and vel is at 0 + 1/2dt

    return

  end subroutine initial_leapfrog_fep

end module sp_md_leapfrog_mod
