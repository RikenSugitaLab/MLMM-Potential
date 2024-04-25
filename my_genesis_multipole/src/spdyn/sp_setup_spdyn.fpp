!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_setup_spdyn
!> @brief   setup variables and structures in MD (DD) simulaton
!! @authors Jaewoon Jung (JJ), Takaharu Mori (TM), Chigusa Kobayashi (CK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_setup_spdyn_mod

  use sp_control_mod
  use sp_restart_mod
  use sp_output_mod
  use sp_input_mod
  use sp_minimize_mod
  use sp_dynamics_mod
  use sp_dynvars_mod
  use sp_domain_mod
  use sp_ensemble_mod
  use sp_enefunc_table_mod
  use sp_restraints_mod
  use sp_constraints_mod
  use sp_boundary_mod
  use sp_pairlist_mod
  use sp_enefunc_mod
  use sp_energy_mod
  use sp_energy_pme_mod
  use sp_parallel_io_mod
  use sp_communicate_mod
  use sp_output_str_mod
  use sp_minimize_str_mod
  use sp_dynamics_str_mod
  use sp_dynvars_str_mod
  use sp_ensemble_str_mod
  use sp_restraints_str_mod
  use sp_constraints_str_mod
  use sp_boundary_str_mod
  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_enefunc_fit_mod
  use sp_energy_str_mod
  use sp_domain_str_mod
  use sp_remd_str_mod
  use sp_rpath_str_mod
  use sp_remd_mod
  use sp_rpath_mod
  use sp_experiments_mod
  use sp_gamd_mod
  use molecules_mod
  use molecules_str_mod
  use fitting_mod
  use fitting_str_mod
  use sp_fep_topology_mod
  use fileio_localres_mod
  use fileio_grocrd_mod
  use fileio_grotop_mod
  use fileio_ambcrd_mod
  use fileio_prmtop_mod
  use fileio_top_mod
  use fileio_par_mod
  use fileio_gpr_mod
  use fileio_psf_mod
  use fileio_pdb_mod
  use fileio_crd_mod
  use fileio_rst_mod
  use fileio_mode_mod
  use structure_check_mod
  use messages_mod
  use timers_mod
  use mpi_parallel_mod
  use constants_mod
  use sp_alchemy_mod
  use sp_alchemy_str_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: setup_spdyn_md
  public  :: setup_spdyn_min
  public  :: setup_spdyn_remd
  public  :: setup_spdyn_rpath
  public  :: setup_spdyn_md_pio
  public  :: setup_spdyn_min_pio
  private :: read_parallel_io_rst
  private :: save_parallel_io_t0
  !FEP
  public  :: setup_spdyn_md_fep
  public  :: setup_spdyn_min_fep
  public  :: setup_spdyn_remd_fep

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_spdyn_md
  !> @brief        setup variables and structures in MD simulation
  !! @authors      JJ
  !! @param[in]    ctrl_data   : information of control parameters
  !! @param[out]   output      : information of output
  !! @param[out]   molecule    : information of molecules
  !! @param[out]   enefunc     : information of energy function
  !! @param[out]   pairlist    : information of nonbonded pairlist
  !! @param[out]   dynvars     : information of dynamic variables
  !! @param[out]   dynamics    : information of molecular dynamics
  !! @param[out]   constraints : information of constraints
  !! @param[out]   ensemble    : information of ensemble
  !! @param[out]   boundary    : information of boundary condition
  !! @param[out]   domain      : information of each domain
  !! @param[out]   comm        : communicator for domain
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_spdyn_md(ctrl_data, output, molecule, enefunc, pairlist,    &
                           dynvars, dynamics, constraints, ensemble, boundary, &
                           domain, comm)

    ! formal arguments
    type(s_ctrl_data),       intent(inout) :: ctrl_data
    type(s_output),          intent(inout) :: output
    type(s_molecule),        intent(inout) :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_dynamics),        intent(inout) :: dynamics
    type(s_constraints),     intent(inout) :: constraints
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_boundary),        intent(inout) :: boundary
    type(s_domain),          intent(inout) :: domain
    type(s_comm),            intent(inout) :: comm

    ! local variables
    type(s_restraints)       :: restraints
    type(s_top)              :: top
    type(s_par)              :: par
    type(s_gpr)              :: gpr
    type(s_psf)              :: psf
    type(s_prmtop)           :: prmtop
    type(s_grotop)           :: grotop
    type(s_pdb)              :: pdb
    type(s_crd)              :: crd
    type(s_ambcrd)           :: ambcrd
    type(s_grocrd)           :: grocrd
    type(s_rst)              :: rst
    type(s_pdb)              :: ref
    type(s_pdb)              :: fit
    type(s_ambcrd)           :: ambref
    type(s_grocrd)           :: groref
    type(s_mode)             :: mode
    type(s_localres)         :: localres


    ! setup parallel I/O
    !

    if (pio_check_ranked_file(ctrl_data%inp_info%rstfile)) then

      call setup_spdyn_md_pio(ctrl_data, output, enefunc, pairlist,     &
                              dynvars, dynamics, constraints, ensemble, &
                              boundary, domain, comm)
      return

    end if


    ! read input files
    !
    call input_md(ctrl_data%inp_info, top, par, psf, prmtop, grotop,  &
                  pdb, crd, ambcrd, grocrd, rst, ref, ambref, groref, &
                  localres, mode)


    ! define molecules
    !
    call define_molecules(molecule, pdb, crd, top, par, gpr, psf, ref, fit, &
                          mode, prmtop, ambcrd, ambref, grotop, grocrd, groref)

    call dealloc_pdb_all(pdb)
    call dealloc_crd_all(crd)
    call dealloc_top_all(top)
    call dealloc_gpr_all(gpr)
    call dealloc_psf_all(psf)
    call dealloc_pdb_all(ref)
    call dealloc_pdb_all(fit)
    call dealloc_ambcrd_all(ambcrd)
    call dealloc_ambcrd_all(ambref)
    call dealloc_grocrd_all(grocrd)
    call dealloc_grocrd_all(groref)


    ! restart coordinates, velocity and boundary
    !
    if (rst%rstfile_type /= RstfileTypeUndef) then
      call setup_restart_pre(rst, molecule)
    end if


    ! set parameters for boundary condition
    !
    call setup_boundary(ctrl_data%bound_info,            &
                        ctrl_data%ene_info%table,        &
                        ctrl_data%ene_info%pairlistdist, &
                        ctrl_data%ene_info%water_model,  &
                        ctrl_data%ens_info%ensemble,     &
                        ctrl_data%cons_info%rigid_bond,  &
                        ctrl_data%ene_info%dsize_cg,     &
                        ctrl_data%ene_info%dmin_size_cg, &
                        molecule, rst, boundary)


    ! set parameters for domain 
    !
    call setup_domain(ctrl_data%ene_info,  &
                      ctrl_data%cons_info, &
                      boundary, molecule, enefunc, constraints, domain)


    ! set parameters for communication
    !
    call setup_communicate(boundary, domain, comm)
    call setup_communicate_size(domain, comm)


    ! set parameters for restraints
    !
    call setup_restraints(ctrl_data%res_info, &
                          ctrl_data%sel_info, &
                          molecule, restraints)


    ! setup enefunc in each domain
    !
    call define_enefunc(ctrl_data%ene_info,  &
                        par, prmtop, grotop, &
                        localres, molecule, constraints, restraints, &
                        domain, enefunc)

    call setup_fitting_spdyn(.false., ctrl_data%fit_info, ctrl_data%sel_info, &
                             domain, molecule, enefunc)

    ! setup experiments
    !
    call setup_experiments(ctrl_data%exp_info, molecule, restraints, &
                           enefunc)

    call dealloc_restraints_all(restraints)
    call dealloc_localres(localres, LocalRestraint)


    ! set parameters for pairlist
    !
    call setup_pairlist(enefunc, domain, pairlist)


    ! set parameters for PME
    !
    if (ctrl_data%ene_info%electrostatic == ElectrostaticPME) then
      call setup_pme(domain, boundary, enefunc)
      call pme_pre  (domain, boundary)
    end if


    ! set parameters for dynamics
    !
    call setup_dynamics(ctrl_data%dyn_info,   &
                        ctrl_data%bound_info, &
                        ctrl_data%res_info, molecule, dynamics)

    ! set parameters for dynamic variables
    !
    call setup_dynvars(dynvars, dynamics)


    ! set parameters for ensemble
    !
    call setup_ensemble(ctrl_data%ens_info, dynamics, ensemble)


    ! set parameters for constraints
    !
    call setup_constraints(ctrl_data%cons_info, &
                           par, prmtop, grotop, molecule, enefunc, constraints)

    call dealloc_molecules_all(molecule)
    call dealloc_par_all(par)
    call dealloc_prmtop_all(prmtop)
    call dealloc_grotop_all(grotop)


    ! set gamd
    !
    call setup_gamd(ctrl_data%gamd_info, dynamics, domain, enefunc)


    ! set output
    !
    call setup_output_md(ctrl_data%out_info, dynamics, output)


    ! restart other variables
    !
    if (rst%rstfile_type /= RstfileTypeUndef) then
      call setup_restart_post(rst, dynamics, dynvars)
      call dealloc_rst_all(rst)
    end if


    ! save t0 information ( case of normal start and parallel restart out )
    !
    if (pio_check_ranked_file(output%rstfile)) then
      call save_parallel_io_t0(ctrl_data, boundary, .true.)
    end if

    if (enefunc%nonb_limiter .and. main_rank) then
      write(MsgOut,'(A,F12.8)')  &
        'Setup_Spdyn_Md> nonb_limiter : minimim distance= ', &
          sqrt(enefunc%minimum_contact)
      write(MsgOut,'(A)') 
    endif

    domain%num_deg_freedom = molecule%num_deg_freedom

    return

  end subroutine setup_spdyn_md

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_spdyn_min
  !> @brief        setup variables and structures in minimization
  !! @authors      TM, JJ
  !! @param[in]    ctrl_data   : information of control parameters
  !! @param[out]   output      : information of output
  !! @param[out]   molecule    : information of molecules
  !! @param[out]   enefunc     : information of energy function
  !! @param[out]   pairlist    : information of nonbonded pairlist
  !! @param[out]   dynvars     : information of dynamic variables
  !! @param[out]   minimize    : information of minimize
  !! @param[out]   constraints : information of constraints
  !! @param[out]   boundary    : information of boundary condition
  !! @param[out]   domain      : information of each domain
  !! @param[out]   comm        : communicator for domain
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_spdyn_min(ctrl_data, output, molecule, enefunc, pairlist, &
                             dynvars, minimize, constraints, boundary,       &
                             domain, comm)

    ! formal arguments
    type(s_ctrl_data),       intent(inout) :: ctrl_data
    type(s_output),          intent(inout) :: output
    type(s_molecule),        intent(inout) :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_minimize),        intent(inout) :: minimize
    type(s_constraints),     intent(inout) :: constraints
    type(s_boundary),        intent(inout) :: boundary
    type(s_domain),          intent(inout) :: domain
    type(s_comm),            intent(inout) :: comm

    ! local variables
    type(s_restraints)       :: restraints
    type(s_top)              :: top
    type(s_par)              :: par
    type(s_gpr)              :: gpr
    type(s_psf)              :: psf
    type(s_prmtop)           :: prmtop
    type(s_grotop)           :: grotop
    type(s_pdb)              :: pdb
    type(s_crd)              :: crd
    type(s_ambcrd)           :: ambcrd
    type(s_grocrd)           :: grocrd
    type(s_rst)              :: rst
    type(s_pdb)              :: ref
    type(s_pdb)              :: fit
    type(s_ambcrd)           :: ambref
    type(s_grocrd)           :: groref
    type(s_mode)             :: mode
    type(s_localres)         :: localres
 
 
    ! setup parallel I/O
    !

    if (pio_check_ranked_file(ctrl_data%inp_info%rstfile)) then

      call setup_spdyn_min_pio(ctrl_data, output, enefunc, pairlist,     &
                               dynvars, minimize, constraints, boundary, &
                               domain, comm)
      return

    end if


    ! read input files
    !
    call input_min(ctrl_data%inp_info, top, par, psf, prmtop, grotop,  &
                   pdb, crd, ambcrd, grocrd, rst, ref, ambref, groref, &
                   localres, mode)


    ! define molecules
    !
    call define_molecules(molecule, pdb, crd, top, par, gpr, psf, ref, fit, &
                          mode, prmtop, ambcrd, ambref, grotop, grocrd, groref)

    call dealloc_pdb_all(pdb)
    call dealloc_crd_all(crd)
    call dealloc_top_all(top)
    call dealloc_psf_all(psf)
    call dealloc_pdb_all(ref)
    call dealloc_pdb_all(fit)
    call dealloc_ambcrd_all(ambcrd)
    call dealloc_ambcrd_all(ambref)
    call dealloc_grocrd_all(grocrd)
    call dealloc_grocrd_all(groref)


    ! restart coordinates, velocity and boundary
    !
    if (rst%rstfile_type /= RstfileTypeUndef) then
      call setup_restart_pre(rst, molecule)
    end if

    ! setup structure_check
    !
    call setup_structure_check(molecule,                                      &
                     ctrl_data%ene_info%forcefield_char, molecule%atom_coord, &
                     ctrl_data%min_info%check_structure,                      &
                     ctrl_data%min_info%fix_ring_error,                       &
                     ctrl_data%min_info%fix_chirality_error,                  &
                     ctrl_data%min_info%exclude_ring_grpid,                   &
                     ctrl_data%min_info%exclude_chiral_grpid)

    ! set parameters for boundary condition
    !
    call setup_boundary(ctrl_data%bound_info,            &
                        ctrl_data%ene_info%table,        &
                        ctrl_data%ene_info%pairlistdist, &
                        ctrl_data%ene_info%water_model,  &
                        ctrl_data%ens_info%ensemble,     &
                        ctrl_data%cons_info%rigid_bond,  &
                        ctrl_data%ene_info%dsize_cg,     &
                        ctrl_data%ene_info%dmin_size_cg, &
                        molecule, rst, boundary)


    ! set parameters for domain
    !
    if (ctrl_data%cons_info%rigid_bond) then
      ctrl_data%cons_info%rigid_bond = .false.
      if (main_rank) then
        write(Msgout, '(A)') 'Setup_Constraints> WARNING : &
                            & constraints are applied only for water'
        write(Msgout, '(A)') 
      end if
    end if
      
    call setup_domain(ctrl_data%ene_info,  &
                      ctrl_data%cons_info, &
                      boundary, molecule, enefunc, constraints, domain)


    ! set parameters for communication
    !
    call setup_communicate(boundary, domain, comm)
    call setup_communicate_size(domain, comm)


    ! set parameters for restraints
    !
    call setup_restraints(ctrl_data%res_info, &
                          ctrl_data%sel_info, &
                          molecule, restraints)


    ! setup enefunc in each domain
    !
    call define_enefunc(ctrl_data%ene_info,  &
                        par, prmtop, grotop, &
                        localres, molecule, constraints, restraints, &
                        domain, enefunc)

    call setup_fitting_spdyn(.false., ctrl_data%fit_info, ctrl_data%sel_info, &
                             domain, molecule, enefunc)

    call setup_constraints(ctrl_data%cons_info, par, prmtop, &
                           grotop, molecule, enefunc, constraints)

    ! setup experiments
    !
    call setup_experiments(ctrl_data%exp_info, molecule, restraints, &
                           enefunc)

    call dealloc_restraints_all(restraints)
    call dealloc_molecules_all(molecule)
    call dealloc_par_all(par)
    call dealloc_prmtop_all(prmtop)
    call dealloc_grotop_all(grotop)
    call dealloc_localres(localres, LocalRestraint)
    

    ! set parameters for pairlist
    !
    call setup_pairlist(enefunc, domain, pairlist)


    ! set parameters for PME
    !
    if (ctrl_data%ene_info%electrostatic == ElectrostaticPME) then
      call setup_pme(domain, boundary, enefunc)
      call pme_pre  (domain, boundary)
    end if


    ! set parameters for minimize
    !
    call setup_minimize(ctrl_data%min_info, minimize)


    ! set parameters for dynamic variables
    !
    call setup_dynvars(dynvars)

    ! set output
    !
    call setup_output_min(ctrl_data%out_info, minimize, output)


    ! save t0 information ( case of normal start and parallel restart out )
    !
    if (pio_check_ranked_file(output%rstfile)) then
      call save_parallel_io_t0(ctrl_data, boundary, .false.)
    end if

    domain%num_deg_freedom = molecule%num_deg_freedom

    if (enefunc%nonb_limiter .and. main_rank) then
      write(MsgOut,'(A,F12.8)')  &
        'Setup_Spdyn_Min> nonb_limiter : minimim distance= ', &
          sqrt(enefunc%minimum_contact)
      write(MsgOut,'(A)') 
    endif

    return

  end subroutine setup_spdyn_min

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_spdyn_remd
  !> @brief        setup variables and structures in REMD simulation
  !! @authors      TM
  !! @param[in]    ctrl_data   : information of control parameters
  !! @param[out]   output      : information of output
  !! @param[out]   molecule    : information of molecules
  !! @param[out]   enefunc     : information of energy function
  !! @param[out]   pairlist    : information of nonbonded pairlist
  !! @param[out]   dynvars     : information of dynamic variables
  !! @param[out]   dynamics    : information of molecular dynamics
  !! @param[out]   constraints : information of constraints
  !! @param[out]   ensemble    : information of ensemble
  !! @param[out]   boundary    : information of boundary condition
  !! @param[out]   domain      : information of each domain
  !! @param[out]   comm        : communicator for domain
  !! @param[out]   remd        : information of remd
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_spdyn_remd(ctrl_data, output, molecule, enefunc, pairlist,  &
                           dynvars, dynamics, constraints, ensemble, boundary, &
                           domain, comm, remd)

    ! formal arguments
    type(s_ctrl_data),       intent(inout) :: ctrl_data
    type(s_output),          intent(inout) :: output
    type(s_molecule),        intent(inout) :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_dynamics),        intent(inout) :: dynamics
    type(s_constraints),     intent(inout) :: constraints
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_boundary),        intent(inout) :: boundary
    type(s_domain),          intent(inout) :: domain
    type(s_comm),            intent(inout) :: comm
    type(s_remd),            intent(inout) :: remd

    ! local variables
    type(s_restraints)       :: restraints
    type(s_top)              :: top
    type(s_par)              :: par
    type(s_prmtop)           :: prmtop
    type(s_grotop)           :: grotop
    type(s_gpr)              :: gpr
    type(s_psf)              :: psf
    type(s_pdb)              :: pdb
    type(s_crd)              :: crd
    type(s_ambcrd)           :: ambcrd
    type(s_grocrd)           :: grocrd
    type(s_rst)              :: rst
    type(s_pdb)              :: ref
    type(s_pdb)              :: fit
    type(s_ambcrd)           :: ambref
    type(s_grocrd)           :: groref
    type(s_mode)             :: mode
    type(s_localres)         :: localres


    ! read input files
    !
    call input_remd(ctrl_data%inp_info, top, par, psf, prmtop, grotop, &
                    pdb, crd, ambcrd, grocrd, rst, ref, ambref, groref,&
                    localres, mode)

    ! define molecules
    !
    call define_molecules(molecule, pdb, crd, top, par, gpr, psf, ref, fit, &
                          mode, prmtop, ambcrd, ambref, grotop, grocrd, groref)

    call dealloc_pdb_all(pdb)
    call dealloc_crd_all(crd)
    call dealloc_top_all(top)
    call dealloc_gpr_all(gpr)
    call dealloc_psf_all(psf)
    call dealloc_pdb_all(ref)
    call dealloc_pdb_all(fit)
    call dealloc_ambcrd_all(ambcrd)
    call dealloc_ambcrd_all(ambref)
    call dealloc_grocrd_all(grocrd)
    call dealloc_grocrd_all(groref)

    ! restart coordinates, velocity and boundary
    !
    if (rst%rstfile_type /= RstfileTypeUndef) then
      call setup_restart_pre(rst, molecule)
    end if

    ! set parameters for boundary condition
    !
    call setup_boundary(ctrl_data%bound_info,            &
                        ctrl_data%ene_info%table,        &
                        ctrl_data%ene_info%pairlistdist, &
                        ctrl_data%ene_info%water_model,  &
                        ctrl_data%ens_info%ensemble,     &
                        ctrl_data%cons_info%rigid_bond,  &
                        ctrl_data%ene_info%dsize_cg,     &
                        ctrl_data%ene_info%dmin_size_cg, &
                        molecule, rst, boundary)

    ! set parameters for restraints
    !
    call setup_restraints(ctrl_data%res_info, &
                          ctrl_data%sel_info, &
                          molecule, restraints)

    call setup_solute_tempering(ctrl_data%rep_info, &
                                molecule, restraints, ctrl_data%cons_info)

    ! set parameters for domain 
    !
    call setup_domain(ctrl_data%ene_info,  &
                      ctrl_data%cons_info, &
                      boundary, molecule, enefunc, constraints, domain)

    ! set parameters for communication
    !
    call setup_communicate(boundary, domain, comm)
    call setup_communicate_size(domain, comm)

    ! setup enefunc in each domain
    !
    call define_enefunc(ctrl_data%ene_info, par, prmtop, grotop,     &
                        localres, molecule, constraints, restraints, &
                        domain, enefunc)

    call setup_fitting_spdyn(.false., ctrl_data%fit_info, ctrl_data%sel_info, &
                             domain, molecule, enefunc)

    ! setup experiments
    !
    call setup_experiments(ctrl_data%exp_info, molecule, restraints, &
                           enefunc)

    call dealloc_localres(localres, LocalRestraint)

    ! set parameters for pairlist
    !
    call setup_pairlist(enefunc, domain, pairlist)

    ! set parameters for PME
    !
    if (ctrl_data%ene_info%electrostatic == ElectrostaticPME) then
      call setup_pme(domain, boundary, enefunc)
      call pme_pre  (domain, boundary)
    end if

    ! set parameters for dynamics
    !
    call setup_dynamics(ctrl_data%dyn_info,   &
                        ctrl_data%bound_info, &
                        ctrl_data%res_info, molecule, dynamics)

    ! set parameters for dynamic variables
    !
    call setup_dynvars(dynvars, dynamics)

    ! set parameters for ensemble
    !
    call setup_ensemble(ctrl_data%ens_info, dynamics, ensemble)

    ! set parameters for constraints
    !
    call setup_constraints(ctrl_data%cons_info, &
                           par, prmtop, grotop, molecule, enefunc, constraints)
    call dealloc_par_all(par)
    call dealloc_prmtop_all(prmtop)
    call dealloc_grotop_all(grotop)

    ! setup remd
    !
    call setup_remd(ctrl_data%rep_info, rst, boundary, dynamics, molecule, &
                    domain, restraints, ensemble, enefunc, remd)
    call dealloc_restraints_all(restraints)
    call dealloc_molecules_all(molecule)

    ! set gamd
    !
    call setup_gamd(ctrl_data%gamd_info, dynamics, domain, enefunc, remd)

    ! set output
    !
    call setup_output_remd(ctrl_data%out_info, dynamics, remd, output)

    ! restart other variables
    !
    if (rst%rstfile_type /= RstfileTypeUndef) then
      call setup_restart_post(rst, dynamics, dynvars)
      call dealloc_rst_all(rst)
    end if

    domain%num_deg_freedom = molecule%num_deg_freedom

    if (enefunc%nonb_limiter .and. main_rank) then
      write(MsgOut,'(A,F12.8)')  &
        'Setup_Spdyn_Remd> nonb_limiter : minimim distance= ', &
          sqrt(enefunc%minimum_contact)
      write(MsgOut,'(A)') 
    endif

    return

  end subroutine setup_spdyn_remd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_spdyn_rpath
  !> @brief        setup variables and structures in RPATH simulation
  !! @authors      YK, YM
  !! @param[in]    ctrl_data   : information of control parameters
  !! @param[out]   output      : information of output
  !! @param[out]   molecule    : information of molecules
  !! @param[out]   enefunc     : information of energy function
  !! @param[out]   pairlist    : information of nonbonded pairlist
  !! @param[out]   dynvars     : information of dynamic variables
  !! @param[out]   dynamics    : information of molecular dynamics
  !! @param[out]   constraints : information of constraints
  !! @param[out]   ensemble    : information of ensemble
  !! @param[out]   boundary    : information of boundary condition
  !! @param[out]   domain      : information of each domain
  !! @param[out]   comm        : communicator for domain
  !! @param[out]   rpath       : information of rpath
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_spdyn_rpath(ctrl_data, output, molecule, enefunc, pairlist, &
                           dynvars, dynamics, constraints, ensemble, boundary, &
                           domain, comm, rpath)

    ! formal arguments
    type(s_ctrl_data),       intent(inout) :: ctrl_data
    type(s_output),          intent(inout) :: output
    type(s_molecule),        intent(inout) :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_dynamics),        intent(inout) :: dynamics
    type(s_constraints),     intent(inout) :: constraints
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_boundary),        intent(inout) :: boundary
    type(s_domain),          intent(inout) :: domain
    type(s_comm),            intent(inout) :: comm
    type(s_rpath),           intent(inout) :: rpath

    ! local variables
    type(s_restraints)       :: restraints
    type(s_top)              :: top
    type(s_par)              :: par
    type(s_prmtop)           :: prmtop
    type(s_grotop)           :: grotop
    type(s_gpr)              :: gpr
    type(s_psf)              :: psf
    type(s_pdb)              :: pdb
    type(s_crd)              :: crd
    type(s_ambcrd)           :: ambcrd
    type(s_grocrd)           :: grocrd
    type(s_rst)              :: rst
    type(s_pdb)              :: ref
    type(s_pdb)              :: fit
    type(s_ambcrd)           :: ambref
    type(s_grocrd)           :: groref
    type(s_mode)             :: mode
    type(s_localres)         :: localres


    ! define replica
    !
    call define_nreplica(ctrl_data%rpath_info, rpath)

    ! read input files
    !
    call input_rpath(ctrl_data%inp_info, top, par, psf, prmtop, grotop, &
                    pdb, crd, ambcrd, grocrd, rst, ref, fit, ambref, groref,&
                    localres, mode)

    ! define molecules
    !
    call define_molecules(molecule, pdb, crd, top, par, gpr, psf, ref, fit, &
                          mode, prmtop, ambcrd, ambref, grotop, grocrd, groref)

    call dealloc_pdb_all(pdb)
    call dealloc_crd_all(crd)
    call dealloc_top_all(top)
    call dealloc_gpr_all(gpr)
    call dealloc_psf_all(psf)
    call dealloc_pdb_all(ref)
    call dealloc_pdb_all(fit)
    call dealloc_ambcrd_all(ambcrd)
    call dealloc_ambcrd_all(ambref)
    call dealloc_grocrd_all(grocrd)
    call dealloc_grocrd_all(groref)

    ! restart coordinates, velocity and boundary
    !
    if (rst%rstfile_type /= RstfileTypeUndef) then
      call setup_restart_pre(rst, molecule)
    end if

    ! set parameters for boundary condition
    !
    call setup_boundary(ctrl_data%bound_info,            &
                        ctrl_data%ene_info%table,        &
                        ctrl_data%ene_info%pairlistdist, &
                        ctrl_data%ene_info%water_model,  &
                        ctrl_data%ens_info%ensemble,     &
                        ctrl_data%cons_info%rigid_bond,  &
                        ctrl_data%ene_info%dsize_cg,     &
                        ctrl_data%ene_info%dmin_size_cg, &
                        molecule, rst, boundary)

    ! set parameters for domain 
    !
    call setup_domain(ctrl_data%ene_info,  &
                      ctrl_data%cons_info, &
                      boundary, molecule, enefunc, constraints, domain)

    ! set parameters for communication
    !
    call setup_communicate(boundary, domain, comm)
    call setup_communicate_size(domain, comm)

    ! set parameters for restraints
    !
    call setup_restraints(ctrl_data%res_info, &
                          ctrl_data%sel_info, &
                          molecule, restraints)

    ! setup enefunc in each domain
    !
    call define_enefunc(ctrl_data%ene_info, par, prmtop, grotop,     &
                        localres, molecule, constraints, restraints, &
                        domain, enefunc)

    call setup_fitting_spdyn(.true., ctrl_data%fit_info, ctrl_data%sel_info, &
                             domain, molecule, enefunc)

    call dealloc_localres(localres, LocalRestraint)

    ! set parameters for pairlist
    !
    call setup_pairlist(enefunc, domain, pairlist)

    ! set parameters for PME
    !
    if (ctrl_data%ene_info%electrostatic == ElectrostaticPME) then
      call setup_pme(domain, boundary, enefunc)
      call pme_pre  (domain, boundary)
    end if

    ! set parameters for dynamics
    !
    call setup_dynamics(ctrl_data%dyn_info,   &
                        ctrl_data%bound_info, &
                        ctrl_data%res_info, molecule, dynamics)

    ! set parameters for dynamic variables
    !
    call setup_dynvars(dynvars, dynamics)

    ! set parameters for ensemble
    !
    call setup_ensemble(ctrl_data%ens_info, dynamics, ensemble)

    ! set parameters for constraints
    !
    call setup_constraints(ctrl_data%cons_info, &
                           par, prmtop, grotop, molecule, enefunc, constraints)
    call dealloc_par_all(par)
    call dealloc_prmtop_all(prmtop)
    call dealloc_grotop_all(grotop)

    ! setup rpath
    !
    call setup_rpath(ctrl_data%rpath_info, rst, boundary, dynamics, molecule, &
                    domain, restraints, ensemble, enefunc, rpath)
    call dealloc_restraints_all(restraints)
    call dealloc_molecules_all(molecule)

    ! set output
    !
    call setup_output_rpath(ctrl_data%out_info, dynamics, output)

    ! restart other variables
    !
    if (rst%rstfile_type /= RstfileTypeUndef) then
      call setup_restart_post(rst, dynamics, dynvars)
      call dealloc_rst_all(rst)
    end if

    domain%num_deg_freedom = molecule%num_deg_freedom

    if (enefunc%nonb_limiter .and. main_rank) then
      write(MsgOut,'(A,F12.8)')  &
        'Setup_Spdyn_Rpath> nonb_limiter : minimim distance= ', &
          sqrt(enefunc%minimum_contact)
      write(MsgOut,'(A)') 
    endif

    return

  end subroutine setup_spdyn_rpath

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_spdyn_md_pio
  !> @brief        setup variables and structures in MD simulation
  !! @authors      JJ
  !! @param[in]    ctrl_data   : information of control parameters
  !! @param[out]   output      : information of output
  !! @param[out]   enefunc     : information of energy function
  !! @param[out]   pairlist    : information of nonbonded pairlist
  !! @param[out]   dynvars     : information of dynamic variables
  !! @param[out]   dynamics    : information of molecular dynamics
  !! @param[out]   constraints : information of constraints
  !! @param[out]   ensemble    : information of ensemble
  !! @param[out]   boundary    : information of boundary condition
  !! @param[out]   domain      : information of each domain
  !! @param[out]   comm        : communicator for domain
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_spdyn_md_pio(ctrl_data, output, enefunc, pairlist,     &
                                dynvars, dynamics, constraints, ensemble, &
                                boundary, domain, comm)

    ! formal arguments
    type(s_ctrl_data),       intent(inout) :: ctrl_data
    type(s_output),          intent(inout) :: output
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_dynamics),        intent(inout) :: dynamics
    type(s_constraints),     intent(inout) :: constraints
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_boundary),        intent(inout) :: boundary
    type(s_domain),          intent(inout) :: domain
    type(s_comm),            intent(inout) :: comm

    ! local variables
    type(s_localres)         :: localres


    ! read domain restart file (Parallel I/O)
    !
    call read_parallel_io_rst(ctrl_data,   &
                              domain,      &
                              enefunc,     &
                              boundary,    &
                              constraints, &
                              dynvars,     &
                              dynamics)

    ! determine the number of degree of freedom
    !
    domain%num_deg_freedom = 3*domain%num_atom_all

    !
    ! determine restart or not

    dynamics%restart=pio_restart

    ! read input files
    !
    if (ctrl_data%inp_info%local_resfile /= '') &
    call input_localres(ctrl_data%inp_info%local_resfile, localres)


    ! set parameters for boundary condition
    !
    call setup_boundary_pio(ctrl_data%bound_info,            &
                            ctrl_data%ene_info%table,        &
                            ctrl_data%ene_info%pairlistdist, &
                            ctrl_data%ene_info%water_model,  &
                            ctrl_data%ens_info%ensemble,     &
                            ctrl_data%cons_info%rigid_bond,  &
                            ctrl_data%ene_info%dsize_cg,     &
                            ctrl_data%ene_info%dmin_size_cg, &
                            boundary)


    ! set parameters for domain 
    !
    enefunc%contact_check=ctrl_data%ene_info%contact_check
    enefunc%minimum_contact=ctrl_data%ene_info%minimum_contact
    call setup_domain_pio(ctrl_data%cons_info, &
                          boundary, enefunc, constraints, domain)


    ! set parameters for communication
    !
    call setup_communicate(boundary, domain, comm)
    call setup_communicate_size(domain, comm)


    ! setup enefunc in each domain
    !
    call define_enefunc_pio(ctrl_data%ene_info, &
                            localres, constraints, domain, enefunc)

    call dealloc_localres(localres, LocalRestraint)


    ! set pairlist
    !
    call setup_pairlist(enefunc, domain, pairlist)


    ! set parameters for PME
    !
    if (ctrl_data%ene_info%electrostatic == ElectrostaticPME) then
      call setup_pme(domain, boundary, enefunc)
      call pme_pre  (domain, boundary)
    end if


    ! set parameters for dynamics
    !
    call setup_dynamics_pio(ctrl_data%dyn_info,   &
                            ctrl_data%bound_info, &
                            domain, dynamics)


    ! set parameters for ensemble
    !
    call setup_ensemble(ctrl_data%ens_info, dynamics, ensemble)


    ! set parameters for constraints
    !
    call setup_constraints_pio(ctrl_data%cons_info, pio_restart, &
                               enefunc, constraints, domain)


    ! set output
    !
    call setup_output_md(ctrl_data%out_info, dynamics, output)

    ! change step number
    !
    if (dynamics%restart .and. dynamics%integrator == IntegratorLEAP) then
     dynvars%step = 1
    else
     dynvars%step = 0
    endif

    if (enefunc%nonb_limiter .and. main_rank) then
      write(MsgOut,'(A,F12.8)')  &
        'Setup_Spdyn_Md_Pio> nonb_limiter : minimim distance= ', &
          sqrt(enefunc%minimum_contact)
      write(MsgOut,'(A)') 
    endif

    return

  end subroutine setup_spdyn_md_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_spdyn_min_pio
  !> @brief        setup variables and structures in minimization
  !! @authors      JJ
  !! @param[in]    ctrl_data   : information of control parameters
  !! @param[out]   output      : information of output
  !! @param[out]   enefunc     : information of energy function
  !! @param[out]   pairlist    : information of nonbonded pairlist
  !! @param[out]   dynvars     : information of dynamic variables
  !! @param[out]   minimize    : information of minimize
  !! @param[out]   constraints : information of constraints
  !! @param[out]   boundary    : information of boundary condition
  !! @param[out]   domain      : information of each domain
  !! @param[out]   comm        : communicator for domain
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_spdyn_min_pio(ctrl_data, output, enefunc, pairlist,     &
                                 dynvars, minimize, constraints, boundary, &
                                 domain, comm)

    ! formal arguments
    type(s_ctrl_data),       intent(inout) :: ctrl_data
    type(s_output),          intent(inout) :: output
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_minimize),        intent(inout) :: minimize
    type(s_constraints),     intent(inout) :: constraints
    type(s_boundary),        intent(inout) :: boundary
    type(s_domain),          intent(inout) :: domain
    type(s_comm),            intent(inout) :: comm

    ! local variables
    type(s_localres)         :: localres
 
 
    ! read domain restart file (Parallel I/O)
    !
    call read_parallel_io_rst(ctrl_data, &
                              domain,    &
                              enefunc,   &
                              boundary)

    ! read input files
    !
    if (ctrl_data%inp_info%local_resfile /= '') &
    call input_localres(ctrl_data%inp_info%local_resfile, localres)


    ! set parameters for boundary condition
    !
    call setup_boundary_pio(ctrl_data%bound_info,            &
                            ctrl_data%ene_info%table,        &
                            ctrl_data%ene_info%pairlistdist, &
                            ctrl_data%ene_info%water_model,  &
                            ctrl_data%ens_info%ensemble,     &
                            ctrl_data%cons_info%rigid_bond,  &
                            ctrl_data%ene_info%dsize_cg,     &
                            ctrl_data%ene_info%dmin_size_cg, &
                            boundary)


    ! set parameters for domain 
    !
    enefunc%contact_check=ctrl_data%ene_info%contact_check
    enefunc%minimum_contact=ctrl_data%ene_info%minimum_contact
    call setup_domain_pio(ctrl_data%cons_info, &
                          boundary, enefunc, constraints, domain)


    ! set parameters for communication
    !
    call setup_communicate(boundary, domain, comm)
    call setup_communicate_size(domain, comm)


    ! setup enefunc in each domain
    !
    call define_enefunc_pio(ctrl_data%ene_info, &
                            localres, constraints, domain, enefunc)

    call dealloc_localres(localres, LocalRestraint)


    ! set pairlist
    !
    call setup_pairlist(enefunc, domain, pairlist)


    ! set parameters for PME
    !
    if (ctrl_data%ene_info%electrostatic == ElectrostaticPME) then
      call setup_pme(domain, boundary, enefunc)
      call pme_pre  (domain, boundary)
    end if


    ! set parameters for minimize
    !
    call setup_minimize(ctrl_data%min_info, minimize)


    ! set output
    !
    call setup_output_min(ctrl_data%out_info, minimize, output)

    ! change step number
    !
    dynvars%step = 0

    if (enefunc%nonb_limiter .and. main_rank) then
      write(MsgOut,'(A,F12.8)')  &
        'Setup_Spdyn_Min_Pio> nonb_limiter : minimim distance= ', &
          sqrt(enefunc%minimum_contact)
      write(MsgOut,'(A)') 
    endif

    return

  end subroutine setup_spdyn_min_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_parallel_io_rst
  !> @brief        read parallel I/O restart file
  !! @authors      JJ
  !! @param[in]    ctrl_data   : information of control parameters
  !! @param[out]   domain      : information of each domain
  !! @param[out]   enefunc     : information of energy function
  !! @param[out]   boundary    : information of boundary condition
  !! @param[out]   constraints : information of constraints
  !! @param[out]   dynvars     : information of dynamic variables
  !! @param[out]   dynamics    : information of molecular dynamics
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_parallel_io_rst(ctrl_data, domain, enefunc, boundary, &
                                  constraints, dynvars, dynamics)

    ! formal arguments
    type(s_ctrl_data),             intent(in)    :: ctrl_data
    type(s_domain),                intent(inout) :: domain
    type(s_enefunc),               intent(inout) :: enefunc
    type(s_boundary),              intent(inout) :: boundary
    type(s_constraints), optional, intent(inout) :: constraints
    type(s_dynvars),     optional, intent(inout) :: dynvars
    type(s_dynamics),    optional, intent(inout) :: dynamics

    ! local variables
    integer                                      :: i
    logical                                      :: fit_check

    if (main_rank) then
      write(MsgOut,'(A)') 'Read_Parallel_Io_Rst> '
      write(MsgOut,'(A)') ' '
    end if

    call pio_read_domain_rst( &
                          pio_get_ranked_filename(ctrl_data%inp_info%rstfile), &
                          boundary,    &
                          domain,      &
                          enefunc,     &
                          constraints, &
                          dynvars,     &
                          dynamics)

    if (.not. main_rank) &
      return

    write(MsgOut,'(A,A)') ' Parallel I/O file: ', &
                                             trim(ctrl_data%inp_info%rstfile)
    write(MsgOut,'(A,L)') '          Restart : ', pio_restart
    write(MsgOut,'(A)')   ' '

    call pio_check_compatible( &
                          ctrl_data%ene_info%pairlistdist, &
                          ctrl_data%ene_info%table,        &
                          ctrl_data%ene_info%water_model,  &
                          ctrl_data%cons_info%rigid_bond,  &
                          ctrl_data%cons_info%fast_water,  &
                          ctrl_data%cons_info%water_model, &
                          ctrl_data%ens_info%ensemble,     &
                          ctrl_data%bound_info%pbc_info%box_size_x, &
                          ctrl_data%bound_info%pbc_info%box_size_y, &
                          ctrl_data%bound_info%pbc_info%box_size_z, &
                          ctrl_data%bound_info%origin_x,   &
                          ctrl_data%bound_info%origin_y,   &
                          ctrl_data%bound_info%origin_z,   &
                          ctrl_data%bound_info%domain_x,   &
                          ctrl_data%bound_info%domain_y,   &
                          ctrl_data%bound_info%domain_z,   &
                          ctrl_data%sel_info%groups,       &
                          ctrl_data%res_info%function,     &
                          ctrl_data%res_info%constant,     &
                          ctrl_data%res_info%select_index)

    ! compatible check for fit

    if (ctrl_data%fit_info%fitting_method /= FittingMethodNO) then
      fit_check=.false.
      do i = 1,enefunc%num_restraintfuncs
        ! restaint posi is not fitted in pio
        if (enefunc%restraint_kind(i) == RestraintsFuncRMSD .or. &
            enefunc%restraint_kind(i) == RestraintsFuncRMSDCOM .or. &
            enefunc%restraint_kind(i) == RestraintsFuncPC .or. &
            enefunc%restraint_kind(i) == RestraintsFuncPCCOM) then
            fit_check = .true.
            exit
        endif
      end do
      if (fit_check) then
        call error_msg('Setup_Spdyn_Md_Pio> WARNING: Fitting is not allowed')
      else
        if (main_rank) then
          write(MsgOut,*) "Setup_Fitting_Spdyn> NO fitting is applied, skip"
          write(MsgOut,*) 
        endif
        enefunc%fitting_method=FittingMethodNO
      endif
    endif

    return

  end subroutine read_parallel_io_rst

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    save_parallel_io_t0
  !> @brief        save t0 status to parallel I/O system
  !! @authors      JJ
  !! @param[in]    ctrl_data : information of control parameters
  !! @param[in]    domain    : information of each domain
  !! @param[in]    boundary  : information of boundary condition
  !! @param[in]    md        : flag for MD or minimization
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine save_parallel_io_t0(ctrl_data, boundary, md)

    ! formal arguments
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_boundary),        intent(in)    :: boundary
    logical,                 intent(in)    :: md


    pio_t0_info%input_topfile            = ctrl_data%inp_info%topfile
    pio_t0_info%input_parfile            = ctrl_data%inp_info%parfile
    pio_t0_info%energy_pairlistdist      = ctrl_data%ene_info%pairlistdist
    pio_t0_info%energy_table             = ctrl_data%ene_info%table
    pio_t0_info%energy_watermodel        = ctrl_data%ene_info%water_model
    pio_t0_info%boundary_boxsizex        = boundary%box_size_x
    pio_t0_info%boundary_boxsizey        = boundary%box_size_y
    pio_t0_info%boundary_boxsizez        = boundary%box_size_z
    pio_t0_info%boundary_originx         = boundary%origin_x
    pio_t0_info%boundary_originy         = boundary%origin_y
    pio_t0_info%boundary_originz         = boundary%origin_z
    pio_t0_info%boundary_domain_xyz      = boundary%num_domain(1)* &
                                           boundary%num_domain(2)* &
                                           boundary%num_domain(3)

    allocate(pio_t0_info%selection_group(size(ctrl_data%sel_info%groups)))
    pio_t0_info%selection_group(:)       = ctrl_data%sel_info%groups(:)

    allocate(pio_t0_info%restraint_func (size(ctrl_data%res_info%function)))
    allocate(pio_t0_info%restraint_const(size(ctrl_data%res_info%function)))
    allocate(pio_t0_info%restraint_index(size(ctrl_data%res_info%function)))
    pio_t0_info%restraint_func(:)        = ctrl_data%res_info%function(:)
    pio_t0_info%restraint_const(:)       = ctrl_data%res_info%constant(:)
    pio_t0_info%restraint_index(:)       = ctrl_data%res_info%select_index(:)

    if (md) then

      pio_t0_info%constraint_rigidbond   = ctrl_data%cons_info%rigid_bond
      pio_t0_info%constraint_fastwater   = ctrl_data%cons_info%fast_water
      pio_t0_info%constraint_watermodel  = ctrl_data%cons_info%water_model
      pio_t0_info%ensemble_type          = ctrl_data%ens_info%ensemble

    else

      pio_t0_info%constraint_rigidbond   = .false.
      pio_t0_info%constraint_fastwater   = .false.
      pio_t0_info%constraint_watermodel  = 'NONE'
      pio_t0_info%ensemble_type          = 0

    end if

    return

  end subroutine save_parallel_io_t0

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_spdyn_md_fep
  !> @brief        setup variables and structures in MD simulation for FEP
  !! @authors      HO
  !! @param[in]    ctrl_data   : information of control parameters
  !! @param[out]   output      : information of output
  !! @param[out]   molecule    : information of molecules
  !! @param[out]   enefunc     : information of energy function
  !! @param[out]   pairlist    : information of nonbonded pairlist
  !! @param[out]   dynvars     : information of dynamic variables
  !! @param[out]   dynamics    : information of molecular dynamics
  !! @param[out]   constraints : information of constraints
  !! @param[out]   ensemble    : information of ensemble
  !! @param[out]   boundary    : information of boundary condition
  !! @param[out]   domain      : information of each domain
  !! @param[out]   comm        : communicator for domain
  !! @param[out]   alchemy     : information of alchemy
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_spdyn_md_fep(ctrl_data, output, molecule, enefunc, pairlist,    &
                           dynvars, dynamics, constraints, ensemble, boundary, &
                           domain, comm, alchemy)

    ! formal arguments
    type(s_ctrl_data),       intent(inout) :: ctrl_data
    type(s_output),          intent(inout) :: output
    type(s_molecule),        intent(inout) :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_dynamics),        intent(inout) :: dynamics
    type(s_constraints),     intent(inout) :: constraints
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_boundary),        intent(inout) :: boundary
    type(s_domain),          intent(inout) :: domain
    type(s_comm),            intent(inout) :: comm
    type(s_alchemy),         intent(inout) :: alchemy

    ! local variables
    type(s_restraints)       :: restraints
    type(s_top)              :: top
    type(s_par)              :: par
    type(s_gpr)              :: gpr
    type(s_psf)              :: psf
    type(s_prmtop)           :: prmtop
    type(s_grotop)           :: grotop
    type(s_pdb)              :: pdb
    type(s_crd)              :: crd
    type(s_ambcrd)           :: ambcrd
    type(s_grocrd)           :: grocrd
    type(s_rst)              :: rst
    type(s_pdb)              :: ref
    type(s_pdb)              :: fit
    type(s_ambcrd)           :: ambref
    type(s_grocrd)           :: groref
    type(s_mode)             :: mode
    type(s_localres)         :: localres

    ! setup parallel I/O
    !

    if (pio_check_ranked_file(ctrl_data%inp_info%rstfile)) then

      call setup_spdyn_md_pio(ctrl_data, output, enefunc, pairlist,     &
                              dynvars, dynamics, constraints, ensemble, &
                              boundary, domain, comm)
      return

    end if


    ! read input files
    !
    call input_md(ctrl_data%inp_info, top, par, psf, prmtop, grotop,  &
                  pdb, crd, ambcrd, grocrd, rst, ref, ambref, groref, &
                  localres, mode)


    ! define molecules
    !
    call define_molecules(molecule, pdb, crd, top, par, gpr, psf, ref, fit, &
                          mode, prmtop, ambcrd, ambref, grotop, grocrd, groref)

    ! FEP: define singleA, singleB, dualA, dualB, and preserved regions
    call define_fep_topology(molecule, par, prmtop, ctrl_data%sel_info, &
                                ctrl_data%alch_info)

    call dealloc_pdb_all(pdb)
    call dealloc_crd_all(crd)
    call dealloc_top_all(top)
    call dealloc_gpr_all(gpr)
    call dealloc_psf_all(psf)
    call dealloc_pdb_all(ref)
    call dealloc_pdb_all(fit)
    call dealloc_ambcrd_all(ambcrd)
    call dealloc_ambcrd_all(ambref)
    call dealloc_grocrd_all(grocrd)
    call dealloc_grocrd_all(groref)


    ! restart coordinates, velocity and boundary
    !
    if (rst%rstfile_type /= RstfileTypeUndef) then
      call setup_restart_pre(rst, molecule)
    end if


    ! set parameters for boundary condition
    !
    call setup_boundary(ctrl_data%bound_info,            &
                        ctrl_data%ene_info%table,        &
                        ctrl_data%ene_info%pairlistdist, &
                        ctrl_data%ene_info%water_model,  &
                        ctrl_data%ens_info%ensemble,     &
                        ctrl_data%cons_info%rigid_bond,  &
                        ctrl_data%ene_info%dsize_cg,     &
                        ctrl_data%ene_info%dmin_size_cg, &
                        molecule, rst, boundary)


    ! set parameters for domain 
    ! FEP
    call setup_domain_fep(ctrl_data%ene_info,  &
                      ctrl_data%cons_info, &
                      boundary, molecule, enefunc, constraints, domain)


    ! set parameters for communication
    !
    call setup_communicate(boundary, domain, comm)
    call setup_communicate_size(domain, comm)


    ! set parameters for restraints
    !
    call setup_restraints(ctrl_data%res_info, &
                          ctrl_data%sel_info, &
                          molecule, restraints)


    ! setup enefunc in each domain
    !
    call define_enefunc(ctrl_data%ene_info,  &
                        par, prmtop, grotop, &
                        localres, molecule, constraints, restraints, &
                        domain, enefunc)

    call setup_fitting_spdyn(.false., ctrl_data%fit_info, ctrl_data%sel_info, &
                             domain, molecule, enefunc)

    ! setup experiments
    !
    call setup_experiments(ctrl_data%exp_info, molecule, restraints, &
                           enefunc)

    call dealloc_restraints_all(restraints)
    call dealloc_localres(localres, LocalRestraint)


    ! setup alchemy
    ! FEP
    call setup_alchemy_md(ctrl_data%alch_info, enefunc, alchemy)


    ! set parameters for pairlist
    ! FEP
    call setup_pairlist_fep(enefunc, domain, pairlist)


    ! set parameters for PME
    !
    if (ctrl_data%ene_info%electrostatic == ElectrostaticPME) then
      call setup_pme(domain, boundary, enefunc)
      ! FEP
      call pme_pre_fep(domain, boundary)
    end if


    ! set parameters for dynamics
    ! FEP
    call setup_dynamics_fep(ctrl_data%dyn_info,   &
                        ctrl_data%bound_info, &
                        ctrl_data%res_info,   &
                        ctrl_data%alch_info,   &
                        molecule, dynamics)


    ! set parameters for dynamic variables
    !
    call setup_dynvars(dynvars, dynamics)


    ! set parameters for ensemble
    !
    call setup_ensemble(ctrl_data%ens_info, dynamics, ensemble)

    ! FEP: Some thermostat and barostat are not available in FEP
    if (ensemble%tpcontrol  == TpcontrolBerendsen) then
      call error_msg('Setup_Ensemble> Berendsen is not allowed in FEP')
    end if
    if (ensemble%tpcontrol  == TpcontrolNoseHoover) then
      call error_msg('Setup_Ensemble> NoseHoover is not allowed in FEP')
    end if
    if (ensemble%tpcontrol  == TpcontrolMTK) then
      call error_msg('Setup_Ensemble> MTK is not allowed in FEP')
    end if

    ! set parameters for constraints
    !
    call setup_constraints(ctrl_data%cons_info, &
                           par, prmtop, grotop, molecule, enefunc, constraints)

    ! FEP: count hbonds in single topology
    call count_hbonds_single_fep(molecule)

    ! FEP: remove degree of freedom of singleB
    if (constraints%rigid_bond) then
      call update_num_deg_freedom('After removing degrees of freedom &
        &of singleB in FEP',    &
        -3*molecule%num_atoms_fep(2) + molecule%num_hbonds_singleB, &
        molecule%num_deg_freedom)
    else
      call update_num_deg_freedom('After removing degrees of freedom &
        &of singleB in FEP',    &
        -3*molecule%num_atoms_fep(2), &
        molecule%num_deg_freedom)
    end if

    call dealloc_molecules_all(molecule)
    call dealloc_par_all(par)
    call dealloc_prmtop_all(prmtop)
    call dealloc_grotop_all(grotop)


    ! set gamd
    !
    call setup_gamd(ctrl_data%gamd_info, dynamics, domain, enefunc)


    ! set output
    !
    call setup_output_md(ctrl_data%out_info, dynamics, output)


    ! restart other variables
    !
    if (rst%rstfile_type /= RstfileTypeUndef) then
      call setup_restart_post(rst, dynamics, dynvars)
      call dealloc_rst_all(rst)
    end if


    ! save t0 information ( case of normal start and parallel restart out )
    !
    if (pio_check_ranked_file(output%rstfile)) then
      call save_parallel_io_t0(ctrl_data, boundary, .true.)
    end if

    if (enefunc%nonb_limiter .and. main_rank) then
      write(MsgOut,'(A,F12.8)')  &
        'Setup_Spdyn_Md> nonb_limiter : minimim distance= ', &
          sqrt(enefunc%minimum_contact)
      write(MsgOut,'(A)') 
    endif

    domain%num_deg_freedom = molecule%num_deg_freedom

    return

  end subroutine setup_spdyn_md_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_spdyn_min_fep
  !> @brief        setup variables and structures in minimization for FEP
  !! @authors      HO
  !! @param[in]    ctrl_data   : information of control parameters
  !! @param[out]   output      : information of output
  !! @param[out]   molecule    : information of molecules
  !! @param[out]   enefunc     : information of energy function
  !! @param[out]   pairlist    : information of nonbonded pairlist
  !! @param[out]   dynvars     : information of dynamic variables
  !! @param[out]   minimize    : information of minimize
  !! @param[out]   constraints : information of constraints
  !! @param[out]   boundary    : information of boundary condition
  !! @param[out]   domain      : information of each domain
  !! @param[out]   comm        : communicator for domain
  !! @param[out]   alchemy     : information of alchemy
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_spdyn_min_fep(ctrl_data, output, molecule, enefunc, pairlist, &
                             dynvars, minimize, constraints, boundary,       &
                             domain, comm, alchemy)

    ! formal arguments
    type(s_ctrl_data),       intent(inout) :: ctrl_data
    type(s_output),          intent(inout) :: output
    type(s_molecule),        intent(inout) :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_minimize),        intent(inout) :: minimize
    type(s_constraints),     intent(inout) :: constraints
    type(s_boundary),        intent(inout) :: boundary
    type(s_domain),          intent(inout) :: domain
    type(s_comm),            intent(inout) :: comm
    type(s_alchemy),         intent(inout) :: alchemy

    ! local variables
    type(s_restraints)       :: restraints
    type(s_top)              :: top
    type(s_par)              :: par
    type(s_gpr)              :: gpr
    type(s_psf)              :: psf
    type(s_prmtop)           :: prmtop
    type(s_grotop)           :: grotop
    type(s_pdb)              :: pdb
    type(s_crd)              :: crd
    type(s_ambcrd)           :: ambcrd
    type(s_grocrd)           :: grocrd
    type(s_rst)              :: rst
    type(s_pdb)              :: ref
    type(s_pdb)              :: fit
    type(s_ambcrd)           :: ambref
    type(s_grocrd)           :: groref
    type(s_mode)             :: mode
    type(s_localres)         :: localres


    ! setup parallel I/O
    !

    if (pio_check_ranked_file(ctrl_data%inp_info%rstfile)) then

      call setup_spdyn_min_pio(ctrl_data, output, enefunc, pairlist,     &
                               dynvars, minimize, constraints, boundary, &
                               domain, comm)
      return

    end if


    ! read input files
    !
    call input_min(ctrl_data%inp_info, top, par, psf, prmtop, grotop,  &
                   pdb, crd, ambcrd, grocrd, rst, ref, ambref, groref, &
                   localres, mode)

    ! define molecules
    !
    call define_molecules(molecule, pdb, crd, top, par, gpr, psf, ref, fit, &
                          mode, prmtop, ambcrd, ambref, grotop, grocrd, groref)

    ! FEP
    call define_fep_topology(molecule, par, prmtop, ctrl_data%sel_info, &
                                ctrl_data%alch_info)

    call dealloc_pdb_all(pdb)
    call dealloc_crd_all(crd)
    call dealloc_top_all(top)
    call dealloc_psf_all(psf)
    call dealloc_pdb_all(ref)
    call dealloc_pdb_all(fit)
    call dealloc_ambcrd_all(ambcrd)
    call dealloc_ambcrd_all(ambref)
    call dealloc_grocrd_all(grocrd)
    call dealloc_grocrd_all(groref)


    ! restart coordinates, velocity and boundary
    !
    if (rst%rstfile_type /= RstfileTypeUndef) then
      call setup_restart_pre(rst, molecule)
    end if


    ! set parameters for boundary condition
    !
    call setup_boundary(ctrl_data%bound_info,            &
                        ctrl_data%ene_info%table,        &
                        ctrl_data%ene_info%pairlistdist, &
                        ctrl_data%ene_info%water_model,  &
                        ctrl_data%ens_info%ensemble,     &
                        ctrl_data%cons_info%rigid_bond,  &
                        ctrl_data%ene_info%dsize_cg,     &
                        ctrl_data%ene_info%dmin_size_cg, &
                        molecule, rst, boundary)


    ! set parameters for domain
    !
    if (ctrl_data%cons_info%rigid_bond) then
      ctrl_data%cons_info%rigid_bond = .false.
      if (main_rank) then
        write(Msgout, '(A)') 'Setup_Constraints> WARNING : &
                            & constraints are applied only for water'
        write(Msgout, '(A)') 
      end if
    end if
      
    ! FEP
    call setup_domain_fep(ctrl_data%ene_info,  &
                      ctrl_data%cons_info, &
                      boundary, molecule, enefunc, constraints, domain)


    ! set parameters for communication
    !
    call setup_communicate(boundary, domain, comm)
    call setup_communicate_size(domain, comm)


    ! set parameters for restraints
    !
    call setup_restraints(ctrl_data%res_info, &
                          ctrl_data%sel_info, &
                          molecule, restraints)


    ! setup enefunc in each domain
    !
    call define_enefunc(ctrl_data%ene_info, par, prmtop, grotop,     &
                        localres, molecule, constraints, restraints, &
                        domain, enefunc)

    call setup_fitting_spdyn(.false., ctrl_data%fit_info, ctrl_data%sel_info, &
                             domain, molecule, enefunc)

    call setup_constraints(ctrl_data%cons_info, par, prmtop, &
                           grotop, molecule, enefunc, constraints)

    ! setup experiments
    !
    call setup_experiments(ctrl_data%exp_info, molecule, restraints, &
                           enefunc)

    ! FEP: count hbonds in single topology
    call count_hbonds_single_fep(molecule)

    ! FEP: remove degree of freedom of singleB
    if (constraints%rigid_bond) then
      call update_num_deg_freedom('After removing degrees of freedom &
        &of singleB in FEP',    &
        -3*molecule%num_atoms_fep(2) + molecule%num_hbonds_singleB, &
        molecule%num_deg_freedom)
    else
      call update_num_deg_freedom('After removing degrees of freedom &
        &of singleB in FEP',    &
        -3*molecule%num_atoms_fep(2), &
        molecule%num_deg_freedom)
    end if

    call dealloc_restraints_all(restraints)
    call dealloc_molecules_all(molecule)
    call dealloc_par_all(par)
    call dealloc_prmtop_all(prmtop)
    call dealloc_grotop_all(grotop)
    call dealloc_localres(localres, LocalRestraint)
    
    ! setup alchemy
    ! FEP
    call setup_alchemy_min(ctrl_data%alch_info, enefunc, alchemy)

    ! set parameters for pairlist
    ! FEP
    call setup_pairlist_fep(enefunc, domain, pairlist)

    ! set parameters for PME
    !
    if (ctrl_data%ene_info%electrostatic == ElectrostaticPME) then
      call setup_pme(domain, boundary, enefunc)
      ! FEP
      call pme_pre_fep(domain, boundary)
    end if

    ! set parameters for minimize
    !
    call setup_minimize(ctrl_data%min_info, minimize)

    ! set parameters for dynamic variables
    !
    call setup_dynvars(dynvars)

    ! set output
    !
    call setup_output_min(ctrl_data%out_info, minimize, output)


    ! save t0 information ( case of normal start and parallel restart out )
    !
    if (pio_check_ranked_file(output%rstfile)) then
      call save_parallel_io_t0(ctrl_data, boundary, .false.)
    end if

    domain%num_deg_freedom = molecule%num_deg_freedom

    if (enefunc%nonb_limiter .and. main_rank) then
      write(MsgOut,'(A,F12.8)')  &
        'Setup_Spdyn_Min> nonb_limiter : minimim distance= ', &
          sqrt(enefunc%minimum_contact)
      write(MsgOut,'(A)') 
    endif

    return

  end subroutine setup_spdyn_min_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_spdyn_remd_fep
  !> @brief        setup variables and structures in FEP simulation
  !! @authors      NK, HO
  !! @param[in]    ctrl_data   : information of control parameters
  !! @param[out]   output      : information of output
  !! @param[out]   molecule    : information of molecules
  !! @param[out]   enefunc     : information of energy function
  !! @param[out]   pairlist    : information of nonbonded pairlist
  !! @param[out]   dynvars     : information of dynamic variables
  !! @param[out]   dynamics    : information of molecular dynamics
  !! @param[out]   constraints : information of constraints
  !! @param[out]   ensemble    : information of ensemble
  !! @param[out]   boundary    : information of boundary condition
  !! @param[out]   domain      : information of each domain
  !! @param[out]   comm        : communicator for domain
  !! @param[out]   remd        : information of remd
  !! @param[out]   alchemy     : information of alchemy
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_spdyn_remd_fep(ctrl_data, output, molecule, enefunc, &
                           pairlist, dynvars, dynamics, constraints, ensemble, &
                           boundary, domain, comm, remd, alchemy)

    ! formal arguments
    type(s_ctrl_data),       intent(inout) :: ctrl_data
    type(s_output),          intent(inout) :: output
    type(s_molecule),        intent(inout) :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_dynvars),         intent(inout) :: dynvars
    type(s_dynamics),        intent(inout) :: dynamics
    type(s_constraints),     intent(inout) :: constraints
    type(s_ensemble),        intent(inout) :: ensemble
    type(s_boundary),        intent(inout) :: boundary
    type(s_domain),          intent(inout) :: domain
    type(s_comm),            intent(inout) :: comm
    type(s_remd),            intent(inout) :: remd
    type(s_alchemy),         intent(inout) :: alchemy

    ! local variables
    type(s_restraints)       :: restraints
    type(s_top)              :: top
    type(s_par)              :: par
    type(s_prmtop)           :: prmtop
    type(s_grotop)           :: grotop
    type(s_gpr)              :: gpr
    type(s_psf)              :: psf
    type(s_pdb)              :: pdb
    type(s_crd)              :: crd
    type(s_ambcrd)           :: ambcrd
    type(s_grocrd)           :: grocrd
    type(s_rst)              :: rst
    type(s_pdb)              :: ref
    type(s_pdb)              :: fit
    type(s_ambcrd)           :: ambref
    type(s_grocrd)           :: groref
    type(s_mode)             :: mode
    type(s_localres)         :: localres


    ! read input files
    !
    call input_remd(ctrl_data%inp_info, top, par, psf, prmtop, grotop, &
                    pdb, crd, ambcrd, grocrd, rst, ref, ambref, groref,&
                    localres, mode)

    ! define molecules
    !
    call define_molecules(molecule, pdb, crd, top, par, gpr, psf, ref, fit, &
                          mode, prmtop, ambcrd, ambref, grotop, grocrd, groref)

    ! FEP
    call define_fep_topology(molecule, par, prmtop, ctrl_data%sel_info, &
                                ctrl_data%alch_info)

    call dealloc_pdb_all(pdb)
    call dealloc_crd_all(crd)
    call dealloc_top_all(top)
    call dealloc_gpr_all(gpr)
    call dealloc_psf_all(psf)
    call dealloc_pdb_all(ref)
    call dealloc_pdb_all(fit)
    call dealloc_ambcrd_all(ambcrd)
    call dealloc_ambcrd_all(ambref)
    call dealloc_grocrd_all(grocrd)
    call dealloc_grocrd_all(groref)

    ! restart coordinates, velocity and boundary
    !
    if (rst%rstfile_type /= RstfileTypeUndef) then
      call setup_restart_pre(rst, molecule)
    end if

    ! set parameters for boundary condition
    !
    call setup_boundary(ctrl_data%bound_info,            &
                        ctrl_data%ene_info%table,        &
                        ctrl_data%ene_info%pairlistdist, &
                        ctrl_data%ene_info%water_model,  &
                        ctrl_data%ens_info%ensemble,     &
                        ctrl_data%cons_info%rigid_bond,  &
                        ctrl_data%ene_info%dsize_cg,     &
                        ctrl_data%ene_info%dmin_size_cg, &
                        molecule, rst, boundary)

    ! set parameters for restraints
    !
    call setup_restraints(ctrl_data%res_info, &
                          ctrl_data%sel_info, &
                          molecule, restraints)

    call setup_solute_tempering(ctrl_data%rep_info, &
                                molecule, restraints, ctrl_data%cons_info)

    ! set parameters for domain 
    ! FEP
    call setup_domain_fep(ctrl_data%ene_info,  &
                      ctrl_data%cons_info, &
                      boundary, molecule, enefunc, constraints, domain)

    ! set parameters for communication
    !
    call setup_communicate(boundary, domain, comm)
    call setup_communicate_size(domain, comm)

    ! setup enefunc in each domain
    !
    call define_enefunc(ctrl_data%ene_info, par, prmtop, grotop,     &
                        localres, molecule, constraints, restraints, &
                        domain, enefunc)

    call setup_fitting_spdyn(.false., ctrl_data%fit_info, ctrl_data%sel_info, &
                             domain, molecule, enefunc)

    ! setup experiments
    !
    call setup_experiments(ctrl_data%exp_info, molecule, restraints, &
                           enefunc)

    call dealloc_localres(localres, LocalRestraint)

    ! setup alchemy
    ! FEP
    call setup_alchemy_remd(ctrl_data%alch_info, enefunc, alchemy)

    ! set parameters for pairlist
    ! FEP
    call setup_pairlist_fep(enefunc, domain, pairlist)

    ! set parameters for PME
    !
    if (ctrl_data%ene_info%electrostatic == ElectrostaticPME) then
      call setup_pme(domain, boundary, enefunc)
      ! FEP
      call pme_pre_fep(domain, boundary)
    end if

    ! set parameters for dynamics
    ! FEP
    call setup_dynamics_fep(ctrl_data%dyn_info,   &
                        ctrl_data%bound_info, &
                        ctrl_data%res_info,   &
                        ctrl_data%alch_info,   &
                        molecule, dynamics)

    ! set parameters for dynamic variables
    !
    call setup_dynvars(dynvars, dynamics)

    ! set parameters for ensemble
    !
    call setup_ensemble(ctrl_data%ens_info, dynamics, ensemble)

    ! FEP: Some thermostat and barostat are not available in FEP
    if (ensemble%tpcontrol  == TpcontrolBerendsen) then
      call error_msg('Setup_Ensemble> Berendsen is not allowed in FEP')
    end if
    if (ensemble%tpcontrol  == TpcontrolNoseHoover) then
      call error_msg('Setup_Ensemble> NoseHoover is not allowed in FEP')
    end if
    if (ensemble%tpcontrol  == TpcontrolMTK) then
      call error_msg('Setup_Ensemble> MTK is not allowed in FEP')
    end if

    ! set parameters for constraints
    !
    call setup_constraints(ctrl_data%cons_info, &
                           par, prmtop, grotop, molecule, enefunc, constraints)
    call dealloc_par_all(par)
    call dealloc_prmtop_all(prmtop)
    call dealloc_grotop_all(grotop)

    ! setup remd
    !
    call setup_remd(ctrl_data%rep_info, rst, boundary, dynamics, molecule, &
                    domain, restraints, ensemble, enefunc, remd, alchemy)

    ! FEP: count hbonds in single topology
    call count_hbonds_single_fep(molecule)

    ! FEP:remove degree of freedom of singleB
    if (constraints%rigid_bond) then
      call update_num_deg_freedom('After removing degrees of freedom &
        &of singleB in FEP',    &
        -3*molecule%num_atoms_fep(2) + molecule%num_hbonds_singleB, &
        molecule%num_deg_freedom)
    else
      call update_num_deg_freedom('After removing degrees of freedom &
        &of singleB in FEP',    &
        -3*molecule%num_atoms_fep(2), &
        molecule%num_deg_freedom)
    end if

    call dealloc_restraints_all(restraints)
    call dealloc_molecules_all(molecule)

    ! set gamd
    !
    call setup_gamd(ctrl_data%gamd_info, dynamics, domain, enefunc, remd)

    ! set output
    !
    call setup_output_remd(ctrl_data%out_info, dynamics, remd, output)

    ! restart other variables
    !
    if (rst%rstfile_type /= RstfileTypeUndef) then
      call setup_restart_post(rst, dynamics, dynvars)
      call dealloc_rst_all(rst)
    end if

    domain%num_deg_freedom = molecule%num_deg_freedom

    if (enefunc%nonb_limiter .and. main_rank) then
      write(MsgOut,'(A,F12.8)')  &
        'Setup_Spdyn_Remd> nonb_limiter : minimim distance= ', &
          sqrt(enefunc%minimum_contact)
      write(MsgOut,'(A)') 
    endif

    return

  end subroutine setup_spdyn_remd_fep

end module sp_setup_spdyn_mod
