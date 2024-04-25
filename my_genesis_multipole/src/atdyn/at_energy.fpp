!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_energy_mod
!> @brief   compute energy
!! @authors Yuji Sugita (YS), Takaharu Mori (TM), Takashi Imai (TI),
!!          Chigusa Kobayashi (CK), Jaewoon Jung (JJ), Norio Takase (NT),
!!          Yoshinobu Akinaga (YA), Kiyoshi Yagi (KY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_energy_mod

  use at_energy_go_mod
  use at_energy_pme_mod
  use at_energy_eef1_mod
  use at_energy_gbsa_mod
  use at_energy_restraints_mod
  use at_energy_nonbonds_mod
  use at_energy_dihedrals_mod
  use at_energy_angles_mod
  use at_energy_bonds_mod
  use at_boundary_str_mod
  use at_boundary_mod
  use at_pairlist_str_mod
  use at_enefunc_str_mod
  use at_energy_str_mod
  use at_restraints_str_mod
  use at_constraints_str_mod
  use at_constraints_mod
  use at_qmmm_mod
  use at_energy_gamd_mod
  use molecules_str_mod
  use fileio_control_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use math_libs_mod
  use string_mod
  ! added by YKL
  use at_driver_mod
  ! END
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  integer, parameter      :: FakeDefault      = 0

  ! structures
  type, public :: s_ene_info
    integer               :: forcefield       = ForcefieldCHARMM
    character(20)         :: forcefield_char  = 'CHARMM'
    integer               :: electrostatic    = ElectrostaticPME
    real(wp)              :: switchdist       = 10.0_wp
    real(wp)              :: cutoffdist       = 12.0_wp
    real(wp)              :: pairlistdist     = 13.5_wp
    real(wp)              :: dielec_const     =  1.0_wp
    logical               :: vdw_force_switch = .false.
    logical               :: vdw_shift        = .false.
    logical               :: cmap_pspline     = .false.
    logical               :: contact_check    = .false.
    logical               :: nonb_limiter     = .false.
    real(wp)              :: pme_alpha ! do not initialize here
    real(wp)              :: pme_alpha_tol    = 1.0e-5
    integer               :: pme_ngrid_x      = FakeDefault
    integer               :: pme_ngrid_y      = FakeDefault
    integer               :: pme_ngrid_z      = FakeDefault
    integer               :: pme_nspline      = 4
    logical               :: pme_multiple     = .false.
    real(wp)              :: pme_mul_ratio    = 1.0_wp
    real(wp)              :: pme_max_spacing  = 1.2_wp
    logical               :: table            = .true.
    real(wp)              :: table_density    = 20.0_wp
    integer               :: table_order      = 1
    character(5)          :: water_model      = 'TIP3'
    integer               :: output_style     = OutputStyleGENESIS
    integer               :: dispersion_corr  = Disp_corr_NONE
    real(wp)              :: minimum_contact  = 0.5_wp
    integer               :: implicit_solvent = ImplicitSolventNONE
    real(wp)              :: gbsa_eps_solvent = 78.5_wp
    real(wp)              :: gbsa_eps_solute  = 1.0_wp
    real(wp)              :: gbsa_alpha       = 1.0_wp
    real(wp)              :: gbsa_beta        = 0.8_wp
    real(wp)              :: gbsa_gamma       = 4.85_wp
    real(wp)              :: gbsa_salt_cons   = 0.2_wp
    real(wp)              :: gbsa_vdw_offset  = 0.09_wp
    real(wp)              :: gbsa_surf_tens   = 0.005_wp
    real(wp)              :: imm1_memb_thick  = 27.0_wp
    integer               :: imm1_exponent_n  = 10
    real(wp)              :: imm1_factor_a    = 0.91_wp
    logical               :: imm1_make_pore   = .false.
    real(wp)              :: imm1_pore_radius = 5.0_wp
    real(wp)              :: imic_axis_a      = 18.0_wp
    real(wp)              :: imic_axis_b      = 18.0_wp
    real(wp)              :: imic_axis_c      = 18.0_wp
    real(wp)              :: imic_exponent_m1 = 1.0_wp
    real(wp)              :: imic_exponent_m2 = 1.0_wp
    real(wp)              :: imic_steepness   = 0.5_wp
    ! if ensemble section does not exist, these temperatures are used (hidden option).
    real(wp)              :: gbsa_temperature = 298.15_wp
    real(wp)              :: eef1_temperature = 298.15_wp
  end type s_ene_info

  ! varibles
  logical, save           :: etitle           = .true.
  logical, save           :: vervose          = .true. 

  ! subroutines
  public  :: show_ctrl_energy
  public  :: read_ctrl_energy
  public  :: setup_energy
  public  :: compute_energy
  public  :: output_energy
  private :: compute_energy_charmm
  private :: compute_energy_charmm19
  private :: compute_energy_go
  private :: compute_energy_amber
  private :: compute_energy_gro_amber
  private :: compute_energy_gro_martini
  public  :: compute_energy_qmmm
  public  :: compute_energy_qmmm_gamd
  private :: compute_energy_spot
  private :: output_energy_genesis
  private :: output_energy_charmm
  private :: output_energy_namd
  private :: output_energy_gromacs
  private :: compute_stats
#ifdef HAVE_MPI_GENESIS
  private :: allreduce_eneforce
#endif

#ifdef QSIMULATE
#include "../lib/qsimulate/interface.hpp"
#endif

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_energy
  !> @brief        show ENERGY section usage
  !! @authors      NT, TM
  !! @param[in]    show_all : show all usage or not
  !! @param[in]    run_mode : run mode string : "md", "min", "remd", "rpath"
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine show_ctrl_energy(show_all, run_mode)

    ! formal arguments
    logical,                 intent(in)    :: show_all
    character(*),            intent(in)    :: run_mode


    if (show_all) then

      select case (run_mode)

      case ('md', 'min', 'vib', 'remd', 'rpath')

        write(MsgOut,'(A)') '[ENERGY]'
        write(MsgOut,'(A)') 'forcefield       = CHARMM  # [CHARMM,AAGO,CAGO,KBGO,AMBER,GROAMBER,GROMARTINI]'
        write(MsgOut,'(A)') 'electrostatic    = PME     # [CUTOFF,PME]'
        write(MsgOut,'(A)') 'switchdist       = 10.0    # switch distance'
        write(MsgOut,'(A)') 'cutoffdist       = 12.0    # cutoff distance'
        write(MsgOut,'(A)') 'pairlistdist     = 13.5    # pair-list distance'
        write(MsgOut,'(A)') '# dielec_const     = 1.0     # dielectric constant'
        write(MsgOut,'(A)') '# vdw_force_switch = NO      # force switch option for van der Waals'
        write(MsgOut,'(A)') '# vdw_shift        = NO      # force switch option for van der Waals'
        write(MsgOut,'(A)') '# cmap_pspline     = NO      # periodic spline for CMAP'
        write(MsgOut,'(A)') '# pme_alpha        = auto    # width of Gaussian distribution in [PME]'
        write(MsgOut,'(A)') '# pme_alpha_tol    = 1.0e-5  # param for auto pme_alpha determination'
        write(MsgOut,'(A)') '# pme_nspline      = 4       # order of B-spline in [PME]'
        write(MsgOut,'(A)') '# pme_multiple     = NO      # separate real and reciprocal parts of [PME] in MPI run'
        write(MsgOut,'(A)') '# pme_mul_ratio    = 1.0     # Nodes ratio (# of real calc. node / # of reciprocal calc. node)'
        write(MsgOut,'(A)') '# table_density    = 20.0    # number of bins used for lookup table'
        write(MsgOut,'(A)') '# water_model      = TIP3    # water model'
        write(MsgOut,'(A)') '# dispersion_corr  = NONE    # dispersion correction [NONE,Energy,EPress]'
        write(MsgOut,'(A)') '# implicit_solvent = NONE    # [NONE,GBSA]'
        write(MsgOut,'(A)') '# gbsa_eps_solvent = 78.5    # solvent dielectric constant in GB'
        write(MsgOut,'(A)') '# gbsa_eps_solute  = 1.0     # solute dielectric constant in GB'
        write(MsgOut,'(A)') '# gbsa_alpha       = 1.0     # scaling factor alpha in GB'
        write(MsgOut,'(A)') '# gbsa_beta        = 0.8     # scaling factor beta in GB'
        write(MsgOut,'(A)') '# gbsa_gamma       = 4.85    # scaling factor gamma in GB'
        write(MsgOut,'(A)') '# gbsa_salt_cons   = 0.2     # salt concentration (mol/L) in GB'
        write(MsgOut,'(A)') '# gbsa_vdw_offset  = 0.09    # van der Waals radius offset in GB'
        write(MsgOut,'(A)') '# gbsa_surf_tens   = 0.005   # surface tension (kcal/mol/A^2) in SA'
        write(MsgOut,'(A)') '# output_style     = GENESIS # format of energy output [GENESIS,CHARMM,NAMD,GROMACS]'
        if (run_mode == 'min') then
          write(MsgOut,'(A)') '# contact_check   = YES     # check atomic clash'
          write(MsgOut,'(A)') '# nonb_limiter    = YES     # avoid failure due to atomic clash'
          write(MsgOut,'(A)') '# minimum_contact = 0.5     # definition of atomic clash distance'
        end if
        write(MsgOut,'(A)') ' '

      end select

    else

      select case (run_mode)

      case ('md', 'min', 'vib', 'remd', 'rpath')

        write(MsgOut,'(A)') '[ENERGY]'
        write(MsgOut,'(A)') 'forcefield       = CHARMM  # [CHARMM,AAGO,CAGO,KBGO,AMBER,GROAMBER,GROMARTINI]'
        write(MsgOut,'(A)') 'electrostatic    = PME     # [CUTOFF,PME]'
        write(MsgOut,'(A)') 'switchdist       = 10.0    # switch distance'
        write(MsgOut,'(A)') 'cutoffdist       = 12.0    # cutoff distance'
        write(MsgOut,'(A)') 'pairlistdist     = 13.5    # pair-list distance'
        write(MsgOut,'(A)') ' '


      end select

    end if

    return

  end subroutine show_ctrl_energy

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      read_ctrl_energy
  !> @brief        read ENERGY section in the control file
  !! @authors      YS, TI, TM, JJ
  !! @param[in]    handle   : unit number of control files
  !! @param[out]   ene_info : ENERGY section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_energy(handle, ene_info)

    ! parameters
    character(*),            parameter     :: Section = 'Energy'

    ! formal arguments
    integer,                 intent(in)    :: handle
    type(s_ene_info),        intent(inout) :: ene_info

    character(MaxLine)                     :: pme_alpha = "auto"


    ! read parameters from control file
    !
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_type   (handle, Section, 'forcefield',       &
                               ene_info%forcefield, ForceFieldTypes)
    call read_ctrlfile_string (handle, Section, 'forcefield',       &
                               ene_info%forcefield_char)
    call read_ctrlfile_type   (handle, Section, 'electrostatic',    &
                               ene_info%electrostatic, ElectrostaticTypes)
    call read_ctrlfile_real   (handle, Section, 'switchdist',       &
                               ene_info%switchdist)
    call read_ctrlfile_real   (handle, Section, 'cutoffdist',       &
                               ene_info%cutoffdist)
    call read_ctrlfile_real   (handle, Section, 'pairlistdist',     &
                               ene_info%pairlistdist)
    call read_ctrlfile_real   (handle, Section, 'dielec_const',     &
                               ene_info%dielec_const)
    call read_ctrlfile_logical(handle, Section, 'vdw_force_switch', &
                               ene_info%vdw_force_switch)
    call read_ctrlfile_logical(handle, Section, 'vdw_shift',        &
                               ene_info%vdw_shift)
    call read_ctrlfile_string (handle, Section, 'pme_alpha',        &
                               pme_alpha)
    call read_ctrlfile_real   (handle, Section, 'pme_alpha_tol',    &
                               ene_info%pme_alpha_tol)
    call read_ctrlfile_integer(handle, Section, 'pme_ngrid_x',      &
                               ene_info%pme_ngrid_x)
    call read_ctrlfile_integer(handle, Section, 'pme_ngrid_y',      &
                               ene_info%pme_ngrid_y)
    call read_ctrlfile_integer(handle, Section, 'pme_ngrid_z',      &
                               ene_info%pme_ngrid_z)
    call read_ctrlfile_integer(handle, Section, 'pme_nspline',      &
                               ene_info%pme_nspline)
    call read_ctrlfile_logical(handle, Section, 'pme_multiple',     &
                               ene_info%pme_multiple)
    call read_ctrlfile_real   (handle, Section, 'pme_mul_ratio',    &
                               ene_info%pme_mul_ratio)
    call read_ctrlfile_real   (handle, Section, 'pme_max_spacing',  &
                               ene_info%pme_max_spacing)
    call read_ctrlfile_real   (handle, Section, 'table_density',    &
                               ene_info%table_density)
    call read_ctrlfile_type   (handle, Section, 'implicit_solvent', &
                               ene_info%implicit_solvent, ImplicitSolventTypes)
    call read_ctrlfile_real   (handle, Section, 'eef1_temperature', &
                               ene_info%eef1_temperature)
    call read_ctrlfile_real   (handle, Section, 'imm1_memb_thick',  &
                               ene_info%imm1_memb_thick)
    call read_ctrlfile_integer(handle, Section, 'imm1_exponent_n',  &
                               ene_info%imm1_exponent_n)
    call read_ctrlfile_real   (handle, Section, 'imm1_factor_a',    &
                               ene_info%imm1_factor_a)
    call read_ctrlfile_logical(handle, Section, 'imm1_make_pore',   &
                               ene_info%imm1_make_pore)
    call read_ctrlfile_real   (handle, Section, 'imm1_pore_radius', &
                               ene_info%imm1_pore_radius)
    call read_ctrlfile_real   (handle, Section, 'imic_axis_a',      &
                               ene_info%imic_axis_a)
    call read_ctrlfile_real   (handle, Section, 'imic_axis_b',      &
                               ene_info%imic_axis_b)
    call read_ctrlfile_real   (handle, Section, 'imic_axis_c',      &
                               ene_info%imic_axis_c)
    call read_ctrlfile_real   (handle, Section, 'imic_exponent_m1', &
                               ene_info%imic_exponent_m1)
    call read_ctrlfile_real   (handle, Section, 'imic_exponent_m2', &
                               ene_info%imic_exponent_m2)
    call read_ctrlfile_real   (handle, Section, 'imic_steepness',   &
                               ene_info%imic_steepness)
    call read_ctrlfile_real   (handle, Section, 'gbsa_eps_solvent', &
                               ene_info%gbsa_eps_solvent)
    call read_ctrlfile_real   (handle, Section, 'gbsa_eps_solute',  &
                               ene_info%gbsa_eps_solute)
    call read_ctrlfile_real   (handle, Section, 'gbsa_alpha',       &
                               ene_info%gbsa_alpha)
    call read_ctrlfile_real   (handle, Section, 'gbsa_beta',        &
                               ene_info%gbsa_beta)
    call read_ctrlfile_real   (handle, Section, 'gbsa_gamma',       &
                               ene_info%gbsa_gamma)
    call read_ctrlfile_real   (handle, Section, 'gbsa_salt_cons',   &
                               ene_info%gbsa_salt_cons)
    call read_ctrlfile_real   (handle, Section, 'gbsa_vdw_offset',  &
                               ene_info%gbsa_vdw_offset)
    call read_ctrlfile_real   (handle, Section, 'gbsa_surf_tens',   &
                               ene_info%gbsa_surf_tens)
    call read_ctrlfile_real   (handle, Section, 'gbsa_temperature', &
                               ene_info%gbsa_temperature)
    call read_ctrlfile_string (handle, Section, 'water_model',      &
                               ene_info%water_model)
    call read_ctrlfile_type   (handle, Section, 'output_style',     &
                               ene_info%output_style, OutputStyleTypes)
    call read_ctrlfile_type   (handle, Section, 'dispersion_corr',  &
                               ene_info%dispersion_corr, Disp_corr_Types)
    call read_ctrlfile_logical(handle, Section, 'contact_check',    &
                               ene_info%contact_check)
    call read_ctrlfile_logical(handle, Section, 'nonb_limiter',     &
                               ene_info%nonb_limiter)
    call read_ctrlfile_real   (handle, Section, 'minimum_contact',  &
                               ene_info%minimum_contact)

    call end_ctrlfile_section(handle)


    ! check table
    !
    if (ene_info%forcefield == ForcefieldKBGO     .or. &
        ene_info%forcefield == ForcefieldAAGO   .or. &
        ene_info%forcefield == ForcefieldCAGO) then
        ene_info%table=.false.
    end if

    if (ene_info%table) then
      if (ene_info%electrostatic == ElectrostaticCutoff ) &
          ene_info%table_order = 3
      if (ene_info%electrostatic == ElectrostaticPME ) &
          ene_info%table_order = 1
    end if

    ! error check for inputs
    !
    if (ene_info%switchdist > ene_info%cutoffdist) then
      call error_msg( &
         'Read_Ctrl_Energy> switchdist must be less than cutoffdist')
    end if

    if (ene_info%cutoffdist >= ene_info%pairlistdist) then
      call error_msg( &
         'Read_Ctrl_Energy> cutoffdist must be less than pairlistdist')
    end if

    if (ene_info%pme_nspline < 3 ) then
      call error_msg( &
         'Read_Ctrl_Energy> "pme_nspline" is too small')

    else if (mod(ene_info%pme_nspline,2) == 1 ) then
      call error_msg( &
         'Read_Ctrl_Energy> "pme_nspline" should be even (in current version)')

    end if

    if (ene_info%water_model .ne. "NONE" .and. &
        ene_info%electrostatic == ElectrostaticCutoff) then
        if (main_rank)      &
        write(MsgOut,'(A)') &
          'Read_Ctrl_Energy>  WARNING: water_model is not available'

          ene_info%water_model = "NONE"

    end if

    ! MK: if pme_alpha is not specified, determine from cutoffdist
    call tolower(pme_alpha)
    if (trim(pme_alpha) .eq. "auto") then
      ene_info%pme_alpha = get_ewald_alpha(ene_info%cutoffdist, &
                                           ene_info%pme_alpha_tol)
    else
      read(pme_alpha,*) ene_info%pme_alpha
    end if

    ! error check for each FFs
    !
    if (ene_info%forcefield == ForcefieldAMBER) then

      if (ene_info%electrostatic == ElectrostaticPME) then
        if (ene_info%switchdist /= ene_info%cutoffdist) then
          if (main_rank)      &
          write(MsgOut,'(A)') &
            'Read_Ctrl_Energy>  WARNING: switchdist should be equal'//     &
            ' to cutoffdist if forcefield is AMBER'
            ene_info%switchdist = ene_info%cutoffdist
        end if

        if (ene_info%vdw_force_switch) then
          if (main_rank)      &
          write(MsgOut,'(A)') &
            'Read_Ctrl_Energy>  WARNING: vdW force switch should not be set'// &
            ' if forcefield is AMBER and electrostatic is PME'
            ene_info%vdw_force_switch = .false.
        end if

        if (ene_info%dispersion_corr /= Disp_corr_EPress) then
          if (main_rank)      &
          write(MsgOut,'(A)') &
            'Read_Ctrl_Energy>  WARNING: dispersion_corr should'//    &
            ' be set to EPRESS if forcefield is AMBER and electrostatic is PME'
          ene_info%dispersion_corr = Disp_corr_EPress
        end if

      else if (ene_info%electrostatic == ElectrostaticCUTOFF) then

        if (ene_info%dispersion_corr /= Disp_corr_NONE) then
          if (main_rank)      &
          write(MsgOut,'(A)') &
            'Read_Ctrl_Energy>  WARNING: dispersion_corr should'//    &
            ' not be set if forcefield is AMBER and electrostatic is cutoff'
          ene_info%dispersion_corr = Disp_corr_NONE
        end if

      end if

    end if

    if (ene_info%forcefield == ForcefieldCHARMM) then
      if (ene_info%dispersion_corr /= Disp_corr_NONE) then
        if (main_rank)      &
        write(MsgOut,'(A)') &
          'Read_Ctrl_Energy>  WARNING: dispersion correction can not'//   &
          ' be set ene_info%cutoffdist if forcefield is CHARMM'
          ene_info%dispersion_corr = Disp_corr_NONE
      end if

    end if

    if (ene_info%forcefield == ForcefieldGROMARTINI) then
      if (ene_info%water_model .ne. "NONE") then
        call error_msg( &
         'Read_Ctrl_Energy> water_model is not allowed, if forcefield'//  &
         ' is GROMARTINI ')
      end if
      if (ene_info%electrostatic .ne. ElectrostaticCUTOFF) then
        call error_msg( &
         'Read_Ctrl_Energy> PME is not allowed, if forcefield'//  &
         ' is GROMARTINI ')
      endif
      if (ene_info%electrostatic == ElectrostaticCUTOFF .and. &
          .not. ene_info%vdw_shift) then
        if (main_rank)      &
          write(MsgOut,'(A)') &
          'Read_Ctrl_Energy>  WARNING: vdW shift is automatically'//     &
          ' set if forcefield is GROMARTINI and CUTOFF'
          ene_info%vdw_shift = .true.
      endif
    end if

    if (ene_info%forcefield /= ForcefieldCHARMM .and. &
        ene_info%forcefield /= ForcefieldAMBER) then
      if (ene_info%vdw_force_switch) then
        if (main_rank)      &
        write(MsgOut,'(A)') &
          'Read_Ctrl_Energy>  WARNING: vdW_force_switch can be used'// &
          ' only for CHARMM'
          ene_info%vdw_force_switch = .false.
      end if
    end if

    if (ene_info%forcefield /= ForcefieldGROAMBER .and. &
        ene_info%forcefield /= ForcefieldGROMARTINI) then
      if (ene_info%vdw_shift) then
        if (main_rank)      &
        write(MsgOut,'(A)') &
          'Read_Ctrl_Energy>  WARNING: vdW_shift can be used only for'//     &
          ' GROAMBER and GROMARTINI'
          ene_info%vdw_shift = .false.
      end if
    end if

    if (ene_info%implicit_solvent == ImplicitSolventEEF1 .or. &
        ene_info%implicit_solvent == ImplicitSolventIMM1) then
      if (ene_info%forcefield /= ForcefieldCHARMM .and. &
          ene_info%forcefield /= ForcefieldCHARMM19) then
        call error_msg('Read_Ctrl_Energy> EEF1/IMM1 is available with CHARMM/CHARMM19')
      end if

      if (ene_info%electrostatic /= ElectrostaticCutoff) then
        call error_msg('Read_Ctrl_Energy> EEF1/IMM1 is available with CUTOFF')
      end if
    end if

    if (ene_info%implicit_solvent == ImplicitSolventGBSA) then
      if (ene_info%forcefield /= ForcefieldCHARMM) then
        call error_msg('Read_Ctrl_Energy> GBSA is available with CHARMM currently')
      end if

      if (ene_info%electrostatic /= ElectrostaticCutoff) then
        call error_msg('Read_Ctrl_Energy> GBSA is available with CUTOFF')
      end if
    end if

    ! Comment out the following lines if longer cutoff is needed.
    ! Note, however, that longer cutoff distance requires a larger memory.
    !
    if (ene_info%switchdist > 100.0_wp) then
      call error_msg( &
         'Read_Ctrl_Energy> switchdist must be less than 100.0')
    end if 

    if (ene_info%cutoffdist > 100.0_wp) then
      call error_msg( &
         'Read_Ctrl_Energy> cutoffdist must be less than 100.0')
    end if

    if (ene_info%pairlistdist > 100.0_wp) then
      call error_msg( &
         'Read_Ctrl_Energy> pairlistdist must be less than 100.0')
    end if 
    
    ! error check
    !
    if (ene_info%contact_check) then
      if (main_rank) &
        write(MsgOut,'(A)') &
      'Read_Ctrl_Energy>  WARNING: in atdyn, contact_check is not available '//&
      'use nonb_limiter=YES'
      ene_info%nonb_limiter = .true.
    endif
    if (.not. ene_info%table) then
      if (ene_info%nonb_limiter) then
        if (main_rank) &
          write(MsgOut,'(A)') &
       'Read_Ctrl_Energy>  WARNING: nonb_limiter is only in table'
        ene_info%nonb_limiter = .false.
      endif
    endif

    ! write parameters to MsgOut
    !
    if (main_rank) then

      write(MsgOut,'(A)')  'Read_Ctrl_Energy> Parameters of Energy Calculations'
      write(MsgOut,'(A20,A10)')                                   &
            '  forcefield      = ',                               &
                trim(ForceFieldTypes(ene_info%forcefield))      
                                                                
      write(MsgOut,'(A20,F10.3,A20,F10.3)')                       &
            '  switchdist      = ', ene_info%switchdist,          &
            '  cutoffdist      = ', ene_info%cutoffdist         

      if (ene_info%electrostatic == ElectrostaticCUTOFF .and. &
           ene_info%vdw_shift) &
        write(MsgOut,'(A20,F10.3)') '  switchdist(ele) = ', 0.0_wp

      write(MsgOut,'(A20,F10.3,A20,F10.3)')                       &
            '  pairlistdist    = ', ene_info%pairlistdist,        &
            '  dielec_const    = ', ene_info%dielec_const       
                                                                
      if (ene_info%vdw_force_switch) then                       
        write(MsgOut,'(A)') '  vdW force_switch=        yes'    
      else                                                      
        if (ene_info%forcefield == ForcefieldGROAMBER .or.        &
            ene_info%forcefield == ForcefieldGROMARTINI ) then
          if (ene_info%vdw_shift) then
            if (ene_info%forcefield == ForcefieldGROMARTINI .and. &
               ene_info%electrostatic == ElectrostaticCUTOFF) then
              write(MsgOut,'(A)') '  vdW/elec  shift =       yes'
            else
              write(MsgOut,'(A)') '  vdW shift       =       yes'
            end if
          elseif (ene_info%cutoffdist-ene_info%switchdist < EPS) then
            write(MsgOut,'(A)') '  vdW switch      =        no'
          else 
            write(MsgOut,'(A)') '  vdW switch      =       yes'
          end if
        else
          write(MsgOut,'(A)') ' vdW force_switch =         no'
        end if
      end if

      ! if PME
      if (ene_info%electrostatic == ElectrostaticPME) then

        write(MsgOut,'(A20,A10)')                            &
              '  electrostatic   = ',                        &
              trim(ElectrostaticTypes(ene_info%electrostatic))
        if (trim(pme_alpha) .eq. "auto") then
          write(MsgOut,'(A20,F10.5,A20,E10.3)')                &
                '  pme_alpha       = ', ene_info%pme_alpha,    &
                '  pme_alpha_tol   = ', ene_info%pme_alpha_tol
        else
          write(MsgOut,'(A20,F10.5)')                &
                '  pme_alpha       = ', ene_info%pme_alpha
        endif
        if (ene_info%pme_ngrid_x .ne. FakeDefault .or.  &
            ene_info%pme_ngrid_y .ne. FakeDefault .or.  &
            ene_info%pme_ngrid_z .ne. FakeDefault) then
          write(MsgOut,'(A20,3I10)')                           &
                '  pme_ngrid(x,y,z)= ', ene_info%pme_ngrid_x   &
                                      , ene_info%pme_ngrid_y   &
                                      , ene_info%pme_ngrid_z
        endif
        write(MsgOut,'(A20,I10)')                            &
              '  pme_nspline     = ', ene_info%pme_nspline

        if (ene_info%pme_multiple) then
          write(MsgOut,'(A)') '  pme_multiple    =        yes'
          write(MsgOut,'(A20,F10.3)')                        &
              '  pme_mul_ratio   = ', ene_info%pme_mul_ratio
        else
          write(MsgOut,'(A)') '  pme_multiple    =         no'
        end if

      ! if CUTOFF
      else if (ene_info%electrostatic == ElectrostaticCutoff) then

        write(MsgOut,'(A20,A7)')                             &
              '  electrostatic   = ',                        &
              trim(ElectrostaticTypes(ene_info%electrostatic))

      end if

      if (ene_info%implicit_solvent == ImplicitSolventNONE) then
        write(MsgOut,'(A)')   '  implicit_solvent=       none'
      else if (ene_info%implicit_solvent == ImplicitSolventEEF1) then
        write(MsgOut,'(A)')   '  implicit_solvent=       eef1'
      else if (ene_info%implicit_solvent == ImplicitSolventIMM1) then
        write(MsgOut,'(A)')   '  implicit_solvent=       imm1'
        write(MsgOut,'(A20,F10.5)')                &
                '  imm1_memb_thick = ', ene_info%imm1_memb_thick
        write(MsgOut,'(A20,I10)')                  &
                '  imm1_exponent_n = ', ene_info%imm1_exponent_n
        write(MsgOut,'(A20,F10.5)')                &
                '  imm1_factor_a   = ', ene_info%imm1_factor_a
        if (ene_info%imm1_make_pore) then
          write(MsgOut,'(A)') '  imm1_make_pore  =        yes'
          write(MsgOut,'(A20,F10.5)')                &
                '  imm1_pore_radius= ', ene_info%imm1_pore_radius
        else
          write(MsgOut,'(A)') '  imm1_make_pore  =         no'
        end if
      else if (ene_info%implicit_solvent == ImplicitSolventIMIC) then
        write(MsgOut,'(A)')   '  implicit_solvent=       imic'
        write(MsgOut,'(A20,F10.5)')                &
                '  imm1_memb_thick = ', ene_info%imm1_memb_thick
        write(MsgOut,'(A20,I10)')                  &
                '  imm1_exponent_n = ', ene_info%imm1_exponent_n
        write(MsgOut,'(A20,F10.5)')                &
                '  imm1_factor_a   = ', ene_info%imm1_factor_a
        write(MsgOut,'(A20,F10.5)')                &
                '  imic_axis_a     = ', ene_info%imic_axis_a
        write(MsgOut,'(A20,F10.5)')                &
                '  imic_axis_b     = ', ene_info%imic_axis_b
        write(MsgOut,'(A20,F10.5)')                &
                '  imic_axis_c     = ', ene_info%imic_axis_c
        write(MsgOut,'(A20,F10.5)')                &
                '  imic_exponent_m1= ', ene_info%imic_exponent_m1
        write(MsgOut,'(A20,F10.5)')                &
                '  imic_exponent_m2= ', ene_info%imic_exponent_m2
        write(MsgOut,'(A20,F10.5)')                &
                '  imic_steepness  = ', ene_info%imic_steepness
      else if (ene_info%implicit_solvent == ImplicitSolventGBSA) then
        write(MsgOut,'(A)')   '  implicit_solvent=       gbsa'
        write(MsgOut,'(2(A20,F10.5))')                &
                '  gbsa_eps_solvent= ', ene_info%gbsa_eps_solvent, &
                '  gbsa_eps_solute = ', ene_info%gbsa_eps_solute
        write(MsgOut,'(2(A20,F10.5))')                &
                '  gbsa_alpha      = ', ene_info%gbsa_alpha,       &
                '  gbsa_beta       = ', ene_info%gbsa_beta
        write(MsgOut,'(2(A20,F10.5))')                &
                '  gbsa_gamma      = ', ene_info%gbsa_gamma,       &
                '  gbsa_salt_cons  = ', ene_info%gbsa_salt_cons
        write(MsgOut,'(2(A20,F10.5))')                &
                '  gbsa_vdw_offset = ', ene_info%gbsa_vdw_offset,  &
                '  gbsa_surf_tens  = ', ene_info%gbsa_surf_tens
      end if

      if (ene_info%table) then
        write(MsgOut,'(A20,I10)') '  table_order     = ', &
             ene_info%table_order
        write(MsgOut,'(A20,F10.3)') '  table_density   = ', &
             ene_info%table_density
      else
        write(MsgOut,'(A)') '  table           =         no'
      end if

      write(MsgOut,'(A20,A10)') '  water_model     = ',     &
                    trim(ene_info%water_model)

      write(MsgOut,'(A20,A10)')                             &
            '  output_style    = ',                         &
            trim(OutputStyleTypes(ene_info%output_style))

      if (ene_info%dispersion_corr == Disp_corr_NONE) then
        write(MsgOut,'(A)') '  dispersion_corr =       none'
      else if (ene_info%dispersion_corr == Disp_corr_Energy) then
        write(MsgOut,'(A)') '  dispersion_corr =     energy'
      else if (ene_info%dispersion_corr == Disp_corr_EPress) then
        write(MsgOut,'(A)') '  dispersion_corr =     epress'
      end if

      if (ene_info%nonb_limiter) then
        write(MsgOut,'(A)') '  nonb_limiter    =     yes'
        write(MsgOut,'(A,F10.3)') ' minimum_contact  = ', &
             ene_info%minimum_contact
      else
        write(MsgOut,'(A)') '  nonb_limiter    =      no'
      endif

      write(MsgOut,'(A)') ' '

    end if
    ene_info%minimum_contact=ene_info%minimum_contact*ene_info%minimum_contact

    return

  end subroutine read_ctrl_energy

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_energy
  !> @brief        setup and cleanup energy information
  !! @authors      CK
  !! @param[in]    restraints : molecular information
  !! @param[out]   energy     : dynamic variables
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_energy(restraints, energy)

    ! formal arguments
    type(s_restraints),        intent(in)    :: restraints
    type(s_energy),            intent(inout) :: energy

    ! local variables
    integer                  :: nfuncs

    ! allocate dynvars
    !
    nfuncs = restraints%nfunctions
    call alloc_energy(energy, EneRestraints, nfuncs)

    return

  end subroutine setup_energy

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy
  !> @brief        compute potential energy
  !! @authors      TM, CK, KY
  !! @param[in]    molecule      : information of molecules
  !! @param[in]    enefunc       : information of potential functions
  !! @param[in]    pairlist      : information of nonbonded pair list
  !! @param[in]    boundary      : information of boundary
  !! @param[in]    nonb_ene      : flag for calculate nonbonded energy
  !! @param[in]    nonb_limiter  : flag for nonbond limiter
  !! @param[in]    coord         : coordinates of target systems
  !! @param[inout] energy        : energy information
  !! @param[inout] temporary     : temporary coordinates
  !!                               (used for MPI allreduce)
  !! @param[inout] force         : forces of target systems
  !! @param[inout] virial        : virial of target systems
  !! @param[inout] virial_ext    : extern virial of target systems
  !! @param[in]    constraints   : information of constraints (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy(molecule, enefunc, pairlist, boundary, nonb_ene, &
                            nonb_limiter,                                    &
                            coord, energy, temporary, force, virial,         &
                            virial_ext, constraints)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule
    type(s_enefunc), target, intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    logical,                 intent(in)    :: nonb_ene
    logical,                 intent(in)    :: nonb_limiter
    real(wp),                intent(inout) :: coord(:,:)
    type(s_energy),          intent(inout) :: energy
    real(wp),                intent(inout) :: temporary(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: virial_ext(3,3)
    type(s_constraints), optional, intent(in) :: constraints

    ! local variables
    real(wp)                               :: volume
    integer                                :: i, natom
    type(s_qmmm), pointer                  :: qmmm

    logical                                :: nonb_ene_flag


    call timer(TimerEnergy, TimerOn)

    qmmm => enefunc%qmmm
    natom = size(coord(1,:))

    if (qmmm%do_qmmm) then

      ! set the charges of QM atoms
      if (qmmm%qm_classical) then
        
        do i = 1, qmmm%qm_natoms
          molecule%charge(qmmm%qmatom_id(i)) = qmmm%qm_charge(i)
        end do
        energy%qm_ene = 0.0_wp

      else

        do i = 1, qmmm%qm_natoms
          molecule%charge(qmmm%qmatom_id(i)) = 0.0D+00
        end do

      end if

    end if

    if (enefunc%gamd_use) then 
      nonb_ene_flag = .true.
    else
      nonb_ene_flag = nonb_ene
    end if

    select case (enefunc%forcefield)

    case (ForcefieldCHARMM)

      call compute_energy_charmm(                                    &
                             molecule, enefunc, pairlist, boundary,  &
                             nonb_ene_flag, nonb_limiter,            &
                             coord, energy, temporary,               &
                             force, virial, virial_ext)


    case (ForcefieldCHARMM19)

      call compute_energy_charmm19(                                  &
                             molecule, enefunc, pairlist, boundary,  &
                             nonb_ene_flag, coord, energy, temporary,&
                             force, virial, virial_ext)

    case (ForcefieldAAGO, ForcefieldCAGO, ForcefieldKBGO)

      call compute_energy_go(molecule, enefunc, pairlist, boundary,  &
                             coord, energy, temporary,               &
                             force, virial, virial_ext)

    case (ForcefieldAMBER)

      call compute_energy_amber(                                    &
                             molecule, enefunc, pairlist, boundary, &
                             nonb_ene_flag, nonb_limiter,           &
                             coord, energy, temporary,              &
                             force, virial, virial_ext)

    case (ForcefieldGROAMBER)

      call compute_energy_gro_amber(                                &
                             molecule, enefunc, pairlist, boundary, &
                             nonb_ene_flag, nonb_limiter,           &
                             coord, energy, temporary,              &
                             force, virial, virial_ext)

    case (ForcefieldGROMARTINI)

      call compute_energy_gro_martini(                              &
                             molecule, enefunc, pairlist, boundary, &
                             nonb_ene_flag, nonb_limiter,           &
                             coord, energy, temporary,              &
                             force, virial, virial_ext)

    end select

    ! Dispersion correction
    if (enefunc%dispersion_corr /= Disp_corr_NONE) then
      volume = boundary%box_size_x_ref * &
               boundary%box_size_y_ref * &
               boundary%box_size_z_ref
      energy%disp_corr_energy = enefunc%dispersion_energy/volume
      if (enefunc%dispersion_corr == Disp_corr_EPress) then
        energy%disp_corr_virial = enefunc%dispersion_virial/volume
        virial(1,1) = virial(1,1) + energy%disp_corr_virial
        virial(2,2) = virial(2,2) + energy%disp_corr_virial
        virial(3,3) = virial(3,3) + energy%disp_corr_virial
      end if

    end if

    if (enefunc%rpath_sum_mf_flag) then
      call compute_stats(enefunc)
    end if

    ! QM/MM energy and gradient calculation, added by ykl to boost qm energy
    if (qmmm%do_qmmm) then
      if (enefunc%gamd_use) then !added by ykl
        call compute_energy_qmmm_gamd(enefunc, qmmm, molecule, natom, coord, force, &
              energy)
      else
        call compute_energy_qmmm(molecule, coord, qmmm, energy, force)
      end if
    end if
    
    if (enefunc%gamd_use) then !added by ykl
      call boost_gamd(enefunc, natom, energy, temporary, force, virial)
    end if

    ! remove forces of fixed atoms
    !
    if (present(constraints)) then
      if (constraints%num_fixatm > 0) &
        call clear_fixatm_component(constraints, molecule%num_atoms, force)
    end if

    !write(MsgOut,'(''debug> force'')')
    !do i = 1, molecule%num_atoms
    !  if (main_rank) then
    !    write(MsgOut,'(''debug> '',i8,3f20.6,f8.4)') &
    !                             i, force(1:3,i), molecule%charge(i)
    !  endif
    !end do
    !write(MsgOut,*)

    call timer(TimerEnergy, TimerOff)

    return

  end subroutine compute_energy

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_energy
  !> @brief        output energy
  !! @authors      YS
  !! @param[in]    step    : step
  !! @param[in]    enefunc : potential energy functions information
  !! @param[inout] energy  : energy information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_energy(step, enefunc, energy)

    ! formal arguments
    integer,                 intent(in)    :: step
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_energy),          intent(inout) :: energy


    if (.not. main_rank) return

    select case (enefunc%output_style)

    case (OutputStyleGENESIS)

      call output_energy_genesis(step, enefunc, energy)

    case (OutputStyleCHARMM)

      call output_energy_charmm (step, enefunc, energy)

    case (OutputStyleNAMD)

      call output_energy_namd   (step, enefunc, energy)

    case (OutputStyleGROMACS)

      call output_energy_gromacs(step, enefunc, energy)

    end select

    return

  end subroutine output_energy

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_charmm
  !> @brief        compute CHARMM potential energy
  !! @authors      YS, TI, JJ, TM, CK
  !! @param[in]    molecule     : information of molecules
  !! @param[inout] enefunc      : information of potential functions
  !! @param[in]    pairlist     : information of nonbonded pair list
  !! @param[in]    boundary     : information of boundary
  !! @param[in]    nonb_ene     : flag for calculate nonbonded energy
  !! @param[in]    nonb_limiter : flag for nonbond limimiter
  !! @param[in]    coord        : coordinates of target systems
  !! @param[inout] energy       : energy information
  !! @param[inout] temporary    : temporary coordinates (used for MPI allreduce)
  !! @param[inout] force        : forces of target systems
  !! @param[inout] virial       : virial of target systems
  !! @param[inout] virial_ext   : extended virial of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_charmm(molecule, enefunc, pairlist, boundary, &
                                   nonb_ene, nonb_limiter,                &
                                   coord, energy, temporary,    &
                                   force, virial, virial_ext)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    logical,                 intent(in)    :: nonb_ene
    logical,                 intent(in)    :: nonb_limiter
    real(wp),                intent(in)    :: coord(:,:)
    type(s_energy),          intent(inout) :: energy
    real(wp),                intent(inout) :: temporary(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: virial_ext(3,3)

    ! local variables
    integer                  :: i, j, k, natom
    real(wp)                 :: resene
    logical                  :: calc_force
    integer                  :: dimno_i, dimno_j
    integer                  :: atom_i, atom_j

    ! initialization of energy and forces
    !
    call init_energy(energy)

    natom = size(coord(1,:))

    force     (1:3,1:natom) = 0.0_wp
    virial    (1:3,1:3)     = 0.0_wp
    virial_ext(1:3,1:3)     = 0.0_wp

    if (real_calc) then

      ! bond energy
      !
      call compute_energy_bond(enefunc, coord, force, virial,   &
                            energy%bond)

      ! angle energy
      !
      call compute_energy_angle(enefunc, coord, force, virial,  &
                            energy%angle)

      call compute_energy_urey(enefunc, coord, force, virial,   &
                            energy%urey_bradley)

      ! dihedral and cmap energies
      !
      if (enefunc%gamd_use) then

        call compute_energy_dihed_gamd(enefunc, natom, coord, &
          force, virial, energy)

      else

        call compute_energy_dihed(enefunc, coord, force, virial,&
                            energy%dihedral)

        call compute_energy_cmap(enefunc, coord, force, virial, &
                            energy%cmap)

      end if

      ! improper energy
      !
      call compute_energy_improp(enefunc, coord, force, virial, &
                            energy%improper)

      ! restraint energy
      !
      if (enefunc%restraint_flag) then

        if (enefunc%gamd_use) then

          call compute_energy_restraints_gamd(enefunc, boundary, natom, coord, &
                            force, virial, virial_ext, energy)

        else

          call compute_energy_restraint(enefunc, boundary, coord, force,    &
                            virial, virial_ext, energy)

        end if

      end if

    end if


    ! nonbonded energy
    !
    select case (boundary%type)

    case (BoundaryTypeNOBC)

      if (enefunc%eef1_use) then
        call compute_energy_nonbond_eef1(enefunc, molecule, pairlist, nonb_ene, &
                              coord, force, virial, &
                              energy%electrostatic, &
                              energy%van_der_waals, &
                              energy%solvation)

      else

        call compute_energy_nonbond_nobc(enefunc, molecule, pairlist, nonb_ene, &
                              coord, force, virial, &
                              energy%electrostatic, &
                              energy%van_der_waals)

        if (enefunc%gbsa_use) then
          call compute_energy_gbsa(enefunc, molecule, pairlist, nonb_ene, &
                              coord, force, virial, &
                              energy%solvation)
        end if

      end if

      ! spherical boundary potential
      if (enefunc%spot_use) &
        call compute_energy_spot(molecule, boundary, coord, force, virial, energy)

    case (BoundaryTypePBC)

      if (enefunc%pme_use) then

        call compute_energy_nonbond_pme(enefunc, molecule, pairlist,    &
                            boundary, nonb_ene,   nonb_limiter,  &
                            coord, force, virial, &
                            energy%electrostatic, &
                            energy%van_der_waals)

      else

        call compute_energy_nonbond_cutoff(enefunc, molecule, pairlist, &
                            boundary, nonb_ene,   &
                            coord, force, virial, &
                            energy%electrostatic, &
                            energy%van_der_waals)

      end if

    case default

      call error_msg('Compute_Energy_CHARMM> Unknown boundary condition')

    end select


    ! allreduce energy, force, and virial
    !
#ifdef HAVE_MPI_GENESIS
    call allreduce_eneforce(enefunc, natom, energy, temporary, &
                            force, virial, virial_ext)
#endif

    ! total energy
    !
    resene = 0.0_wp
    if (enefunc%num_restraintfuncs > 0) then
      resene = sum(energy%restraint(1:enefunc%num_restraintfuncs))
    endif

    energy%total = energy%bond          + energy%angle                 &
                 + energy%urey_bradley  + energy%dihedral              &
                 + energy%cmap          + energy%improper              &
                 + energy%electrostatic + energy%van_der_waals         &
                 + energy%solvation     + energy%spot                  &
                 + resene

    ! GaMD boost and statistics
    !
    !if (enefunc%gamd_use) then 
    !  call boost_gamd(enefunc, natom, energy, temporary, force, virial)
    !end if

    return

  end subroutine compute_energy_charmm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_charmm19
  !> @brief        compute CHARMM19 potential energy
  !! @authors      NT
  !! @param[in]    molecule   : information of molecules
  !! @param[in]    enefunc    : information of potential functions
  !! @param[in]    pairlist   : information of nonbonded pair list
  !! @param[in]    boundary   : information of boundary
  !! @param[in]    nonb_ene   : flag for calculate nonbonded energy
  !! @param[in]    coord      : coordinates of target systems
  !! @param[inout] energy     : energy information
  !! @param[inout] temporary  : temporary coordinates (used for MPI allreduce)
  !! @param[inout] force      : forces of target systems
  !! @param[inout] virial     : virial of target systems
  !! @param[inout] virial_ext : extended virial of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_charmm19(molecule, enefunc, pairlist, boundary, &
                                     nonb_ene, coord, energy, temporary, force,&
                                     virial, virial_ext)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    logical,                 intent(in)    :: nonb_ene
    real(wp),                intent(in)    :: coord(:,:)
    type(s_energy),          intent(inout) :: energy
    real(wp),                intent(inout) :: temporary(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: virial_ext(3,3)

    ! local variables
    integer                  :: i, natom
    real(wp)                 :: resene


    ! initialization of energy and forces
    !
    call init_energy(energy)

    natom = size(coord(1,:))

    force(1:3,1:natom)  = 0.0_wp
    virial(1:3,1:3)     = 0.0_wp
    virial_ext(1:3,1:3) = 0.0_wp

    if (real_calc) then

      ! bond energy
      !
      call compute_energy_bond(enefunc, coord, force, virial,   &
                            energy%bond)

      ! angle energy
      !
      call compute_energy_angle(enefunc, coord, force, virial,  &
                            energy%angle)

      ! dihedral energy
      !
      call compute_energy_dihed(enefunc, coord, force, virial,  &
                            energy%dihedral)

      ! improper energy
      !
      call compute_energy_improp(enefunc, coord, force, virial, &
                            energy%improper)

      ! restraint energy
      !
      if (enefunc%restraint_flag) then
        call compute_energy_restraint(enefunc, boundary, coord, force,    &
                            virial, virial_ext, energy)
      end if

    end if


    ! nonbonded energy
    !

    select case (boundary%type)

    case (BoundaryTypeNOBC)

      if (enefunc%eef1_use) then

        call compute_energy_nonbond_eef1(enefunc, molecule, pairlist, nonb_ene, &
                              coord, force, virial, &
                              energy%electrostatic, &
                              energy%van_der_waals, &
                              energy%solvation)

      else

        call compute_energy_nonbond_nobc(enefunc, molecule, pairlist, nonb_ene, &
                              coord, force, virial, &
                              energy%electrostatic, &
                              energy%van_der_waals)

      end if

    case default

      call error_msg &
      ('Compute_Energy_CHARMM19> boundarytype = PBC : currently unavailable!')

    end select


    ! allreduce energy, force, and virial
    !
#ifdef HAVE_MPI_GENESIS
    call allreduce_eneforce(enefunc, natom, energy, temporary, &
                            force, virial, virial_ext)
#endif

    ! total energy
    !
    resene = 0.0_wp
    if (enefunc%num_restraintfuncs > 0) then
      resene = sum(energy%restraint(1:enefunc%num_restraintfuncs))
    endif
    energy%total = energy%bond          + energy%angle                 &
                 + energy%dihedral      + energy%improper              &
                 + energy%electrostatic + energy%van_der_waals         &
                 + energy%solvation                                    &
                 + resene

    return

  end subroutine compute_energy_charmm19

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_go
  !> @brief        compute all-atom Go potential energy
  !! @authors      TM, CK
  !! @param[in]    molecule   : information of molecules
  !! @param[in]    enefunc    : information of potential functions
  !! @param[in]    pairlist   : information of nonbonded pair list
  !! @param[in]    boundary   : information of boundary
  !! @param[in]    coord      : coordinates of target systems
  !! @param[inout] energy     : energy information
  !! @param[inout] temporary  : temporary coordinates (used for MPI allreduce)
  !! @param[inout] force      : forces of target systems
  !! @param[inout] virial     : virial of target systems
  !! @param[inout] virial_ext : extended virial of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_go(molecule, enefunc, pairlist, boundary, coord, &
                               energy, temporary, force, virial, virial_ext)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(wp),                intent(in)    :: coord(:,:)
    type(s_energy),          intent(inout) :: energy
    real(wp),                intent(inout) :: temporary(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: virial_ext(3,3)

    ! local variables
    integer                  :: i, natom
    real(wp)                 :: resene

    ! initialization of energy and forces
    !
    call init_energy(energy)

    natom = size(coord(1,:))

    force(1:3,1:natom)  = 0.0_wp
    virial(1:3,1:3)     = 0.0_wp
    virial_ext(1:3,1:3) = 0.0_wp

    ! bond energy
    !
    call compute_energy_bond(enefunc, coord, force, virial,     &
                            energy%bond)

    ! angle energy
    !
    call compute_energy_angle(enefunc, coord, force, virial,    &
                            energy%angle)

    ! dihedral energy
    !
    if (enefunc%gamd_use) then

      call compute_energy_dihed_gamd(enefunc, natom, coord, &
          force, virial, energy)

    else

      call compute_energy_dihed(enefunc, coord, force, virial,    &
                            energy%dihedral)

    end if

    ! improper energy
    !
    if (enefunc%forcefield == ForceFieldAAGO) then

      call compute_energy_improp(enefunc, coord, force, virial, &
                            energy%improper)

    end if

    ! contact energy
    !
    if (enefunc%forcefield == ForceFieldAAGO) then

      call compute_energy_contact_126(enefunc, coord, force, virial,  &
                            energy%contact)

    else if (enefunc%forcefield == ForceFieldCAGO) then

      call compute_energy_contact_1210(enefunc, coord, force, virial, &
                            energy%contact)

    else if (enefunc%forcefield == ForceFieldKBGO) then

      call compute_energy_contact_12106(enefunc, coord, force, virial, &
                            energy%contact)

    end if

    ! non-contact energy
    !   Note that 1-4 interactions are not calculated in the all-atom Go model
    !
    select case(boundary%type)

    case (BoundaryTypeNOBC)

      call compute_energy_noncontact_nobc(enefunc, molecule, pairlist, &
                               coord, force, virial, energy%noncontact)

    case default

      call error_msg('Compute_Energy_Go> Bad boundary condition for Go model')

    end select

    ! restraint energy
    !
    if (enefunc%restraint_flag) then
    end if

    ! restraint energy
    !
    if (enefunc%restraint_flag) then

      if (enefunc%gamd_use) then

        call compute_energy_restraints_gamd(enefunc, boundary, natom, coord, &
                          force, virial, virial_ext, energy)

      else

        call compute_energy_restraint(enefunc, boundary, coord, force, virial, &
                                    virial_ext, energy)

      end if

    end if

    ! allreduce energy, force, and virial
    !
#ifdef HAVE_MPI_GENESIS
    call allreduce_eneforce(enefunc, natom, energy, temporary, &
                            force, virial, virial_ext)
#endif

    ! total energy
    !
    resene = 0.0_wp
    if (enefunc%num_restraintfuncs > 0) then
      resene = sum(energy%restraint(1:enefunc%num_restraintfuncs))
    endif

    energy%total = energy%bond     + energy%angle                      &
                 + energy%dihedral + energy%improper                   &
                 + energy%contact  + energy%noncontact                 &
                 + resene

    ! GaMD boost and statistics
    !
    !if (enefunc%gamd_use) then 
    !  call boost_gamd(enefunc, natom, energy, temporary, force, virial)
    !end if

    return

  end subroutine compute_energy_go

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_amber
  !> @brief        compute AMBER99 potential energy
  !! @authors      NT
  !! @param[in]    molecule   : information of molecules
  !! @param[inout] enefunc    : information of potential functions
  !! @param[in]    pairlist   : information of nonbonded pair list
  !! @param[in]    boundary   : information of boundary
  !! @param[in]    nonb_ene   : flag for calculate nonbonded energy
  !! @param[in]    nonb_limiter : flag for nonbond limimiter
  !! @param[in]    coord      : coordinates of target systems
  !! @param[inout] energy     : energy information
  !! @param[inout] temporary  : temporary coordinates (used for MPI allreduce)
  !! @param[inout] force      : forces of target systems
  !! @param[inout] virial     : virial of target systems
  !! @param[inout] virial_ext : extended virial of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_amber(molecule, enefunc, pairlist, boundary, &
                                  nonb_ene, nonb_limiter,                &
                                  coord, energy, temporary,              &
                                  force, virial, virial_ext)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    logical,                 intent(in)    :: nonb_ene
    logical,                 intent(in)    :: nonb_limiter
    real(wp),                intent(inout) :: coord(:,:)
    type(s_energy),          intent(inout) :: energy
    real(wp),                intent(inout) :: temporary(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: virial_ext(3,3)

    ! local variables
    integer                  :: i, natom
    real(wp)                 :: resene

    if (.not. enefunc%pme_use .and. boundary%type == BoundaryTypePBC) then
      call error_msg('Compute_Energy_Amber> CUTOFF is only for NOBC ')
    endif

    ! initialization of energy and forces
    !
    call init_energy(energy)

    natom = size(coord(1,:))

    force     (1:3,1:natom) = 0.0_wp
    virial    (1:3,1:3)     = 0.0_wp
    virial_ext(1:3,1:3)     = 0.0_wp

    if (real_calc) then

      ! bond energy
      !
      call compute_energy_bond(enefunc, coord, force, virial,   &
                            energy%bond)

      ! angle energy
      !
      call compute_energy_angle(enefunc, coord, force, virial,  &
                            energy%angle)

      ! dihedral energy
      !
      if (enefunc%gamd_use) then

        call compute_energy_dihed_gamd(enefunc, natom, coord, &
            force, virial, energy)

      else

        call compute_energy_dihed(enefunc, coord, force, virial,  &
                              energy%dihedral)

      end if

      ! improper energy
      !
      call compute_energy_improp_cos(enefunc, coord, force, virial,  &
                            energy%improper)

      ! restraint energy
      !
      if (enefunc%restraint_flag) then

        if (enefunc%gamd_use) then

          call compute_energy_restraints_gamd(enefunc, boundary, natom, coord, &
                              force, virial, virial_ext, energy)

        else

          call compute_energy_restraint(enefunc, boundary, coord, force, virial, &
                              virial_ext, energy)

        end if

      end if

    end if


    ! nonbonded energy
    !
    select case (boundary%type)

    case (BoundaryTypeNOBC)

      call compute_energy_nonbond_nobc(enefunc, molecule, pairlist, nonb_ene, &
                            coord, force, virial, &
                            energy%electrostatic, &
                            energy%van_der_waals)

      ! spherical boundary potential
      if (enefunc%spot_use) &
        call compute_energy_spot(molecule, boundary, coord, force, virial, energy)

    case (BoundaryTypePBC)

      if (enefunc%pme_use) then

        call compute_energy_nonbond_pme(enefunc, molecule, pairlist, &
                            boundary, nonb_ene, nonb_limiter,   &
                            coord, force, virial, &
                            energy%electrostatic, &
                            energy%van_der_waals)

      else

        call compute_energy_nonbond_cutoff(enefunc, molecule, pairlist, &
                            boundary, nonb_ene,   &
                            coord, force, virial, &
                            energy%electrostatic, &
                            energy%van_der_waals)

      end if

    case default

      call error_msg('Compute_Energy_Amber> Unknown boundary condition')

    end select


    ! allreduce energy, force, and virial
    !
#ifdef HAVE_MPI_GENESIS
    call allreduce_eneforce(enefunc, natom, energy, temporary, force, &
                            virial, virial_ext)
#endif

    ! total energy
    !
    resene = 0.0_wp
    if (enefunc%num_restraintfuncs > 0) then
      resene = sum(energy%restraint(1:enefunc%num_restraintfuncs))
    endif

    energy%total = energy%bond          + energy%angle                 &
                 + energy%urey_bradley  + energy%dihedral              &
                 + energy%cmap          + energy%improper              &
                 + energy%electrostatic + energy%van_der_waals         &
                 + energy%spot          + resene

    ! GaMD boost and statistics
    !
    !if (enefunc%gamd_use) then 
    !  call boost_gamd(enefunc, natom, energy, temporary, force, virial)
    !end if

    return

  end subroutine compute_energy_amber


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_gro_amber
  !> @brief        compute GROMACS-AMBER potential energy
  !! @authors      NT
  !! @param[in]    molecule     : information of molecules
  !! @param[inout] enefunc      : information of potential functions
  !! @param[in]    pairlist     : information of nonbonded pair list
  !! @param[in]    boundary     : information of boundary
  !! @param[in]    nonb_ene     : flag for calculate nonbonded energy
  !! @param[in]    nonb_limiter : flag for nonbond limiter
  !! @param[in]    coord        : coordinates of target systems
  !! @param[inout] energy       : energy information
  !! @param[inout] temporary    : temporary coordinates (used for MPI allreduce)
  !! @param[inout] force        : forces of target systems
  !! @param[inout] virial       : virial of target systems
  !! @param[inout] virial_ext   : extended virial of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_gro_amber(molecule, enefunc, pairlist, boundary, &
                                      nonb_ene, nonb_limiter,                &
                                      coord, energy, temporary,              &
                                      force, virial, virial_ext)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    logical,                 intent(in)    :: nonb_ene
    logical,                 intent(in)    :: nonb_limiter
    real(wp),                intent(in)    :: coord(:,:)
    type(s_energy),          intent(inout) :: energy
    real(wp),                intent(inout) :: temporary(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: virial_ext(3,3)

    ! local variables
    integer                  :: i, natom
    real(wp)                 :: resene


    ! initialization of energy and forces
    !
    call init_energy(energy)

    natom = size(coord(1,:))

    force     (1:3,1:natom) = 0.0_wp
    virial    (1:3,1:3)     = 0.0_wp
    virial_ext(1:3,1:3)     = 0.0_wp

    if (real_calc) then

      ! bond energy
      !
      call compute_energy_bond(enefunc, coord, force, virial,   &
                            energy%bond)

      ! angle energy
      !
      call compute_energy_angle(enefunc, coord, force, virial,  &
                            energy%angle)

      ! dihedral energy
      !
      call compute_energy_dihed(enefunc, coord, force, virial,  &
                            energy%dihedral)

      ! Ryckaert-Bellemans dihedral energy
      !
      call compute_energy_rb_dihed(enefunc, coord, force, virial, &
                            energy%dihedral)

      ! restraint energy
      !
      if (enefunc%restraint_flag) then
        call compute_energy_restraint(enefunc, boundary, coord, force, virial, &
                            virial_ext, energy)
      end if

    end if


    ! nonbonded energy
    !
    select case (boundary%type)

    case (BoundaryTypeNOBC)

      call compute_energy_nonbond_nobc(enefunc, molecule, pairlist, nonb_ene, &
                            coord, force, virial, &
                            energy%electrostatic, &
                            energy%van_der_waals)

    case (BoundaryTypePBC)

      if (enefunc%pme_use) then

        call compute_energy_nonbond_pme(enefunc, molecule, pairlist, &
                            boundary, nonb_ene, nonb_limiter,   &
                            coord, force, virial, &
                            energy%electrostatic, &
                            energy%van_der_waals)

      else

        call compute_energy_nonbond_cutoff(enefunc, molecule, pairlist, &
                            boundary, nonb_ene,   &
                            coord, force, virial, &
                            energy%electrostatic, &
                            energy%van_der_waals)

      end if

    case default

      call error_msg('Compute_Energy_Gro_Amber> Unknown boundary condition')

    end select


    ! allreduce energy, force, and virial
    !
#ifdef HAVE_MPI_GENESIS
    call allreduce_eneforce(enefunc, natom, energy, temporary, force, &
                            virial, virial_ext)
#endif

    ! total energy
    !
    resene = 0.0_wp
    if (enefunc%num_restraintfuncs > 0) then
      resene = sum(energy%restraint(1:enefunc%num_restraintfuncs))
    endif

    energy%total = energy%bond          + energy%angle                 &
                 + energy%urey_bradley  + energy%dihedral              &
                 + energy%cmap          + energy%improper              &
                 + energy%electrostatic + energy%van_der_waals         &
                 + resene

    return

  end subroutine compute_energy_gro_amber

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_gro_martini
  !> @brief        compute GROMACS-MARTINI potential energy
  !! @authors      NT
  !! @param[in]    molecule     : information of molecules
  !! @param[inout] enefunc      : information of potential functions
  !! @param[in]    pairlist     : information of nonbonded pair list
  !! @param[in]    boundary     : information of boundary
  !! @param[in]    nonb_ene     : flag for calculate nonbonded energy
  !! @param[in]    nonb_limiter : flag for nonbond limiter
  !! @param[in]    coord        : coordinates of target systems
  !! @param[inout] energy       : energy information
  !! @param[inout] temporary    : temporary coordinates (used for MPI allreduce)
  !! @param[inout] force        : forces of target systems
  !! @param[inout] virial       : virial of target systems
  !! @param[inout] virial_ext   : extended virial of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_gro_martini(molecule, enefunc, pairlist, boundary, &
                                        nonb_ene, nonb_limiter,                &
                                        coord, energy, temporary,              &
                                        force, virial, virial_ext)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    logical,                 intent(in)    :: nonb_ene
    logical,                 intent(in)    :: nonb_limiter
    real(wp),                intent(in)    :: coord(:,:)
    type(s_energy),          intent(inout) :: energy
    real(wp),                intent(inout) :: temporary(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: virial_ext(3,3)

    ! local variables
    integer                  :: i, natom
    real(wp)                 :: resene


    ! initialization of energy and forces
    !
    call init_energy(energy)

    natom = size(coord(1,:))

    force     (1:3,1:natom) = 0.0_wp
    virial    (1:3,1:3)     = 0.0_wp
    virial_ext(1:3,1:3)     = 0.0_wp

    if (real_calc) then

      ! bond energy
      !
      call compute_energy_bond(enefunc, coord, force, virial,      &
                            energy%bond)

      ! angle energy
      !
      call compute_energy_angle_g96(enefunc, coord, force, virial, &
                            energy%angle)

      ! dihedral energy
      !
      call compute_energy_dihed(enefunc, coord, force, virial, &
                            energy%dihedral)

      ! restraint energy
      !
      if (enefunc%restraint_flag) then
        call compute_energy_restraint(enefunc, boundary, coord, force, virial, &
                            virial_ext, energy)
      end if

    end if


    ! nonbonded energy
    !
    select case (boundary%type)

    case (BoundaryTypeNOBC)

      call compute_energy_nonbond_nobc(enefunc, molecule, pairlist, nonb_ene, &
                            coord, force, virial, &
                            energy%electrostatic, &
                            energy%van_der_waals)

    case (BoundaryTypePBC)

      call compute_energy_nonbond_cutoff(enefunc, molecule, pairlist, &
                          boundary, nonb_ene,   &
                          coord, force, virial, &
                          energy%electrostatic, &
                          energy%van_der_waals)

    case default

      call error_msg('Compute_Energy_Gro_Martini> Unknown boundary condition')

    end select


    ! allreduce energy, force, and virial
    !
#ifdef HAVE_MPI_GENESIS
    call allreduce_eneforce(enefunc, natom, energy, temporary, force, &
                            virial, virial_ext)
#endif

    ! total energy
    !
    resene = 0.0_wp
    if (enefunc%num_restraintfuncs > 0) then
      resene = sum(energy%restraint(1:enefunc%num_restraintfuncs))
    endif

    energy%total = energy%bond          + energy%angle                 &
                 + energy%urey_bradley  + energy%dihedral              &
                 + energy%cmap          + energy%improper              &
                 + energy%electrostatic + energy%van_der_waals         &
                 + resene

    return

  end subroutine compute_energy_gro_martini

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_qmmm
  !> @brief        compute QMMM energy and gradients
  !! @authors      YA, KY
  !! @param[in]    molecule   : information of molecules
  !! @param[in]    coord      : coordinates of target systems
  !! @param[in]    qmmm       : QM/MM information
  !! @param[inout] energy     : energy information
  !! @param[inout] force      : force information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_qmmm(molecule, coord, qmmm, energy, force)
#ifdef QSIMULATE
    use iso_c_binding
#endif

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule
    real(wp),                intent(in)    :: coord(:,:)
    type(s_qmmm), target,    intent(inout) :: qmmm
    type(s_energy),          intent(inout) :: energy
    real(wp),                intent(inout) :: force(:,:)

    ! local variables
    character(len_exec)      :: exec
    character(100)           :: count, trial
    !character(len_exec)      :: host !added by YKL
    integer(4)               :: istat, rename
    integer                  :: i, k !, socket, inet !added by YKL
    logical, pointer         :: qmmm_error
    logical                  :: iftest !added by YKL

#ifdef QSIMULATE
    integer               :: iatom
    integer(c_int)        :: natoms
    integer               :: qm_atom, mm_atom
    integer(4)            :: unlink
    logical               :: retry
    real(wp)              :: rr01, r0(3), r1(3), r01(3)
    real(c_double)              :: bagel_energy
    integer(c_int),  allocatable, target :: bagel_atoms(:)
    real(c_double), allocatable, target :: bagel_coord(:, :)
    real(c_double), allocatable, target :: bagel_charges(:)
    real(c_double), allocatable, target :: bagel_force(:,:)
    real(c_double)              :: bagel_dipole(3)
    real(c_double), allocatable, target :: bagel_qmcharges(:)
#endif

    integer               :: n_atoms

    qmmm_error => qmmm%qmmm_error

    if (.not. qmmm%qm_classical) then
      qmmm%qmtrial = 0
      qmmm_error   = .true.
      do while (qmmm_error)

        call timer(TimerQMMM, TimerOn)

        ! prepare QM filenames
        !
        write(count,'(i0)') qmmm%qm_count
        qmmm%qmbasename = trim(qmmm%qmbase0)//count
        if (len(trim(qmmm%qmindex)) > 0) &
          qmmm%qmbasename = trim(qmmm%qmbasename)//'_'//trim(qmmm%qmindex)

        qmmm%qminp  = trim(qmmm%qmbasename)//'.inp'
        qmmm%qmout  = trim(qmmm%qmbasename)//'.log'
        qmmm%qminfo = trim(qmmm%qmbasename)//'.qminfo'

        if(mod(qmmm%qm_count, qmmm%qmsave_period) == 0) then
          qmmm%savefile = .true.
        else
          qmmm%savefile = .false.
        endif

#ifdef QSIMULATE
        if (qmmm%qmtyp == QMtypQSimulate) then
          call qm_generate_input(molecule, coord, qmmm)

          ! we create temporary here
          natoms = qmmm%qm_natoms + qmmm%num_qmmmbonds + qmmm % mm_natoms

          allocate(bagel_atoms(natoms), bagel_coord(3,natoms), bagel_charges(natoms), &
                   bagel_force(3, natoms), bagel_qmcharges(natoms))

          do i = 1, qmmm%qm_natoms
            bagel_atoms(i) = qmmm%qm_atomic_no(i)
            bagel_coord(:, i) = coord(:, qmmm%qmatom_id(i))
            bagel_charges(i) = 0.0_wp
          enddo
          do i = 1, qmmm%num_qmmmbonds
            iatom = qmmm%qm_natoms + i
            bagel_atoms(iatom) = -1
            bagel_coord(:, iatom) = qmmm%linkatom_coord(:, i)
            bagel_charges(iatom) = 0.0_wp
          enddo
          do i = 1, qmmm%mm_natoms
            iatom = qmmm%qm_natoms + qmmm%num_qmmmbonds + i
            bagel_atoms(iatom) = 0
            bagel_coord(:, iatom) = coord(:, qmmm%mmatom_id(i))
            bagel_charges(iatom) = molecule%charge(qmmm%mmatom_id(i))
          enddo

          qmmm%bagel_output = trim(qmmm%workdir)//'/'//trim(qmmm%qmout) // c_null_char
          retry = qmmm%qmtrial/=0

          call qsimulate_interface(mpi_comm_country, my_world_rank, &
                 qmmm%qm_input, qmmm%bagel_output, natoms, c_loc(bagel_atoms), &
                 c_loc(bagel_coord), c_loc(bagel_charges), &
                 bagel_energy, c_loc(bagel_force), bagel_dipole, &
                 c_loc(bagel_qmcharges), retry, qmmm_error)

          energy%qm_ene = bagel_energy
          qmmm%qm_dipole(1:3) = bagel_dipole(1:3)

          if (.not. qmmm%ene_only) then
            do i = 1, qmmm%mm_natoms
              iatom = qmmm%qm_natoms + qmmm%num_qmmmbonds + i
              qmmm%mm_force(:, i) = -bagel_force(:, iatom)
            enddo
            do i = 1, qmmm%num_qmmmbonds
              iatom = qmmm%qm_natoms + i
              qmmm%linkatom_force(:, i) = -bagel_force(:, iatom)
              qmmm%linkatom_charge(i) = bagel_qmcharges(iatom)
            enddo
            do i = 1, qmmm%qm_natoms
              qmmm%qm_force(:, i) = -bagel_force(:, i)
              qmmm%qm_charge(i) = bagel_qmcharges(i)
            enddo
            qmmm%is_qm_charge = .true.
          end if

          if (replica_main_rank) then
            if (.not. qmmm%savefile .and. .not. qmmm_error) then
              if(qmmm%qm_debug) write(MsgOut,'(''QMMM_debug> remove '',a)') &
                                                 trim(qmmm%workdir)//'/'//trim(qmmm%qmout)
              istat = unlink(trim(qmmm%workdir)//'/'//trim(qmmm%qmout))
              if (istat /= 0) then
                call error_msg ('Compute_Energy_QMMM> Error while removing '// &
                                 trim(qmmm%workdir)//'/'//trim(qmmm%qmout))
              end if
            end if

            ! changed by YKL to print QM/MM force and energy
            if (qmmm%savefile .and. qmmm%save_qminfo) then
              call write_qminfo(qmmm, energy%qm_ene)
            end if
            ! end
          end if

#ifdef HAVE_MPI_GENESIS
          ! broadcast QM information to other processes
          call mpi_bcast(energy%qm_ene, 1, mpi_wp_real, 0, mpi_comm_country, ierror)
          call mpi_bcast(qmmm%qm_force(1,1), qmmm%qm_natoms*3, mpi_wp_real, 0,   &
              mpi_comm_country, ierror)
          if (qmmm%mm_natoms /= 0) then
              call mpi_bcast(qmmm%mm_force(1,1), qmmm%mm_natoms*3, mpi_wp_real, 0,     &
                  mpi_comm_country, ierror)
          end if
          call mpi_bcast(qmmm%qm_charge(1) , qmmm%qm_natoms, mpi_wp_real, 0,       &
                         mpi_comm_country, ierror)
          call mpi_bcast(qmmm%is_qm_charge, 1, mpi_logical, 0,                &
                         mpi_comm_country, ierror)
          call mpi_bcast(qmmm%qm_dipole(1), 3, mpi_wp_real, 0, mpi_comm_country, ierror)

          if (qmmm%num_qmmmbonds > 0) then
            call mpi_bcast(qmmm%linkatom_force(1,1), qmmm%num_qmmmbonds*3, mpi_wp_real, 0, &
                           mpi_comm_country, ierror)
            call mpi_bcast(qmmm%linkatom_charge(1) , qmmm%num_qmmmbonds, mpi_wp_real, 0,   &
                           mpi_comm_country, ierror)
          end if
#endif

          if(.not. qmmm_error) call qm_post_process(coord, qmmm, energy, force)

          deallocate(bagel_atoms, bagel_coord, bagel_charges, bagel_force, &
                     bagel_qmcharges)

        else ! all the cases except for QSimulate
#endif   

          if (qmmm%qmtyp == QMtypML) then !added by YKL
            call qm_generate_input(molecule, coord, qmmm)
            n_atoms = qmmm%qm_natoms + qmmm%num_qmmmbonds

            !call compute_electrostatic_mlmm(molecule, coord, qmmm, energy, qmmm%temporary, 1)
            call compute_electrostatic_mlmm(molecule, coord, qmmm)

            if (replica_main_rank) then
              call ml_model( molecule, coord, qmmm, energy)
            end if
            qmmm%is_qm_charge = .true.

#ifdef HAVE_MPI_GENESIS
         ! broadcast QM information to other processes
            call mpi_bcast(energy%qm_ene, 1, mpi_wp_real, 0, mpi_comm_country, ierror)
            call mpi_bcast(qmmm%is_qm_charge, 1, mpi_logical, 0,                &
                           mpi_comm_country, ierror)
            call mpi_bcast(qmmm%qm_force(1,1), qmmm%qm_natoms*3, mpi_wp_real, 0,     &
                           mpi_comm_country, ierror)
            call mpi_bcast(qmmm%T0_grad(1) , n_atoms, mpi_wp_real, 0,       &
                           mpi_comm_country, ierror)
            call mpi_bcast(qmmm%T1_grad(1,1) , n_atoms*3, mpi_wp_real, 0,       &
                           mpi_comm_country, ierror)
            if (qmmm%max_order>1) then
              call mpi_bcast(qmmm%T2_grad(1,1,1) , n_atoms*9, mpi_wp_real, 0,       &
                             mpi_comm_country, ierror)
            endif

            if (qmmm%max_order>2) then
              call mpi_bcast(qmmm%T3_grad(1,1,1,1) , n_atoms*27, mpi_wp_real, 0,       &
                             mpi_comm_country, ierror)
            endif

            if (qmmm%num_qmmmbonds > 0) then
              call mpi_bcast(qmmm%linkatom_force(1,1), qmmm%num_qmmmbonds*3, mpi_wp_real, 0, &
                             mpi_comm_country, ierror)
              call mpi_bcast(qmmm%linkatom_charge(1) , qmmm%num_qmmmbonds, mpi_wp_real, 0,   &
                             mpi_comm_country, ierror)
            end if
#endif

            !call compute_electrostatic_mlmm(molecule, coord, qmmm, energy, qmmm%temporary, 2)
            call compute_qmmm_force(molecule, qmmm, qmmm%temporary)

            !qmmm%is_qm_charge = .true.
            if (replica_main_rank) then
              if (qmmm%savefile .and. qmmm%save_qminfo) then
            !    write(*,*) "is_qm_charge:", qmmm%is_qm_charge
                call write_qminfo(qmmm, energy%qm_ene)
              end if
            end if
            qmmm_error = .false.
            if(.not. qmmm_error) call qm_post_process(coord, qmmm, energy, force)

          else ! all the cases except for QSimulate and MLMM
            if (replica_main_rank) then

              exec = 'cd '//trim(qmmm%workdir)//'; '//trim(qmmm%qmexe)// &
                     ' '//trim(qmmm%qminp)//' '//trim(qmmm%qmout)

              exec          = trim(exec)//' '//trim(count)

              if(qmmm%qm_debug .and. replica_main_rank) then
                write(MsgOut,'(''QMMM_debug> launching QM program '',a)') &
                                                        trim(QMtypTypes(qmmm%qmtyp))
                write(MsgOut,'(''QMMM_debug>  QM command : '',a)') trim(exec)
                write(MsgOut,'(''QMMM_debug>  QMexe      : '',a)') trim(qmmm%qmexe)
                write(MsgOut,'(''QMMM_debug>  QMfolder   : '',a)') trim(qmmm%workdir)
                write(MsgOut,'(''QMMM_debug>  QMinp      : '',a)') trim(qmmm%qminp)
                write(MsgOut,'(''QMMM_debug>  QMout      : '',a)') trim(qmmm%qmout)
                write(MsgOut,'(''QMMM_debug>  SaveFile   : '',l)') qmmm%savefile
              endif

              ! generate QM input
              !
              call qm_generate_input(molecule, coord, qmmm)

              if(qmmm%dryrun) then
                qmmm_error = .false.
                exit
              end if

              ! execute QM program
              !
              call system(trim(exec))

              ! MPI spawn is currently not working
              !
              !    if (trim(qmmm%qm_exemode) .eq. 'system') then
              !      call system(trim(exec))
              !    else if (trim(qmmm%qm_exemode) .eq. 'spawn') then
              !#ifdef HAVE_MPI_GENESIS
              !      call mpi_comm_spawn(exec, MPI_ARGV_NULL, qmmm%qm_nprocs, MPI_INFO_NULL, &
              !              0, mpi_comm_country, qmmm%inter_comm, MPI_ERRCODES_IGNORE, ierror)
              !      if (replica_main_rank) call qm_monitor(qmmm)
              !      call mpi_barrier(mpi_comm_country, ierror)
              !#endif
              !    end if

            else
              if(qmmm%dryrun) then
                qmmm_error = .false.
                exit
              endif

            end if

            ! retrieve QM energy and gradients
            !
            call qm_read_output(molecule, qmmm, energy, qmmm_error)
            if(.not. qmmm_error) call qm_post_process(coord, qmmm, energy, force)

            ! changed by YKL to print QM/MM force and energy
            if (replica_main_rank) then
              if (qmmm%savefile .and. qmmm%save_qminfo) then
                call write_qminfo(qmmm, energy%qm_ene)
              end if
            end if
          ! end
          end if

#ifdef QSIMULATE
        end if
#endif

        call timer(TimerQMMM, TimerOff)

        if (qmmm%ignore_qm_error) exit

        if (qmmm_error) then
          if (replica_main_rank) then
            write(MsgOut,'(A47,I5)') &
             'Compute_Energy_QMMM> Failed in QM calculation. Retry ...  rank_no = ', &
              my_world_rank
            write(trial,'(i0)') qmmm%qmtrial
            istat = rename(trim(qmmm%workdir)//'/'//trim(qmmm%qmout), &
                           trim(qmmm%workdir)//'/'//trim(qmmm%qmout)//'_'//trim(trial))
            if (istat /= 0) then
              call error_msg ('Compute_Energy_QMMM> Error while saveing '// &
                              trim(qmmm%workdir)//'/'//trim(qmmm%qmout))
            end if
          end if

          qmmm%qmtrial = qmmm%qmtrial + 1
          if (qmmm%qmtrial > qmmm%qmmaxtrial) &
            call error_msg ('Compute_Energy_QMMM> maximum number of &
                            &QM retry reached')
        end if

      end do


      ! update count number
      if (nrep_per_proc == 1) then
        qmmm%qm_count = qmmm%qm_count + 1
      else
        if (mod(my_replica_no, nrep_per_proc) == 0) qmmm%qm_count = qmmm%qm_count + 1
      end if  

    end if

    ! output gradient for debug purpose
    if (qmmm%qm_debug .and. .not. qmmm%ene_only .and. replica_main_rank) then
       write(MsgOut,'(''QMMM_debug> QMMM force'')')
       do i = 1, molecule%num_atoms
         write(MsgOut,'(''QMMM_debug> '',i8,3f20.6,f8.4)') &
                                         i, force(1:3,i), molecule%charge(i)
       end do
    end if

    return

  end subroutine compute_energy_qmmm

 !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_qmmm_gamd
  !> @brief        calculate qmmm energies for GaMD
  !! @authors      YKL
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    natom    : number of atoms
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] energy   : energy information
  !! @detail       force and virial are saved to scale them after energy
  !!               evaluation
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_qmmm_gamd(enefunc, qmmm, molecule, natom, coord, &
    force, energy)

    ! formal arguments
    type(s_enefunc), target, intent(inout) :: enefunc
    integer,                 intent(in)    :: natom
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:)
    type(s_molecule),        intent(inout) :: molecule
    type(s_qmmm),            intent(inout) :: qmmm
    type(s_energy),          intent(inout) :: energy

    ! local variables
    integer                  :: id, i, j, k, ix, ic, jc
    real(wp), pointer :: f_qm(:,:)

    if (enefunc%gamd%boost_qm) then

      f_qm => enefunc%gamd%f_qm

      do i = 1, natom
        f_qm(1:3,i) = force(1:3,i)
      end do

    end if

    ! qm energy
    !
    call compute_energy_qmmm(molecule, coord, qmmm, energy, force)

    if (enefunc%gamd%boost_qm) then 

      do i = 1, natom
        f_qm(1:3,i) = force(1:3,i) - f_qm(1:3,i)
      end do

    end if

    return

  end subroutine compute_energy_qmmm_gamd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_spot
  !> @brief        calculate spherical boundary 
  !! @authors      KY
  !! @param[in]    molecule : molecule information
  !! @param[in]    boundary : boundary information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] energy   : energy information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_spot(molecule, boundary, coord, force, virial, energy)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_boundary),        intent(in)    :: boundary
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    type(s_energy),          intent(inout) :: energy

    ! local variables
    integer  :: natom, nfunc, nc
    integer  :: i, j, k, n, nsave
    integer  :: id, my_id
    real(wp) :: rn, rmin, din(3), rin
    real(wp) :: espot, force_local(3), virial_local(3,3), vtmp
    real(wp) :: const, fcoeff
    integer  :: expo
#ifdef OMP
    integer  :: omp_get_thread_num, omp_get_max_threads
#endif

#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif

    natom = molecule%num_atoms
    nfunc = boundary%nfunctions
    nc    = 0

    espot         = 0.0_wp
    force_local   = 0.0_wp
    virial_local  = 0.0_wp

    !$omp parallel                                        &
    !$omp private(i, j, k, n, nsave,                      &
    !$omp         force_local, vtmp,                      &
    !$omp         rn, rmin, din, rin,                     &
    !$omp         const, fcoeff, expo, id, my_id)         &
    !$omp reduction(+:espot, virial_local, nc) 

#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    my_id = my_city_rank * nthread + id

    do i = 1, natom
      if (mod(i-1,nproc_city*nthread) /= my_id) cycle

      ! skip fix atoms
      if (boundary%fixatm(i)) cycle

      ! find the nearest center
      rmin = -1.0_wp
      do n = 1, nfunc
        rn = boundary%radius(n)

        din = coord(:,i) - boundary%center(:,n)
        rin = sqrt(din(1)*din(1) + din(2)*din(2) + din(3)*din(3))

        if (rin < rn) then
          nsave = -1
          exit

        else
          if (rmin < 0.0_wp .or. (rin - rn) < rmin) then
            rmin  = rin - rn
            nsave = n
          end if

        end if

      end do

      ! skip if inside sphere
      if (nsave < 0) cycle

      ! calc ene, force, and virial
      nc = nc + 1
      n  = nsave
      rn    = boundary%radius(n)
      const = boundary%const(n)
      expo  = boundary%exponent(n)

      din = coord(:,i) - boundary%center(:,n)
      rin = sqrt(din(1)*din(1) + din(2)*din(2) + din(3)*din(3))
      din  = din/rin

      if (expo == 2) then
        espot  = espot + const*(rin - rn)*(rin - rn)
        fcoeff = - 2.0_wp*const*(rin - rn)

      else
        espot  = espot + const*(rin - rn)**(expo)
        fcoeff = - real(expo,wp)*const*((rin - rn)**(expo-1))

      end if
      force_local = fcoeff*din
      force(:,i)  = force(:,i)  + force_local

      do j = 1, 3
        do k = j+1, 3
          vtmp = - coord(k,i) * force_local(j)
          virial_local(k,j) = virial_local(k,j) + vtmp
          virial_local(j,k) = virial_local(j,k) + vtmp
        end do
        vtmp = - coord(j,i) * force_local(j)
        virial_local(j,j) = virial_local(j,j) + vtmp
      end do

    end do
    !$omp end parallel

    energy%spot  = espot
    virial       = virial + virial_local

  end subroutine compute_energy_spot

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_energy_genesis
  !> @brief        output energy in GENESIS style
  !! @authors      TM
  !! @param[in]    step    : step
  !! @param[in]    enefunc : information of potential functions
  !! @param[in]    energy  : energy information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_energy_genesis(step, enefunc, energy)

    ! formal arguments
    integer,                 intent(in)    :: step
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_energy),          intent(in)    :: energy

    ! local variables
    integer,parameter        :: clength=16, flength=4
    integer                  :: i, ifm
    character(16)            :: title
    character(16)            :: category(999)
    character                :: frmt*5, frmt_res*10, rfrmt*7
    character                :: rfrmt_cont*9,frmt_cont*7
    real(wp)                 :: values(999)

    write(title,'(A16)') 'STEP'
    write(frmt,'(A2,I2,A)') '(A',clength,')'
    write(frmt_cont,'(A2,I2,A3)') '(A',clength,',$)'
    write(frmt_res,'(A2,I2,A6)') '(A',clength-3,',I3.3)'
    write(rfrmt,'(A2,I2,A1,I1,A1)') '(F',clength,'.',flength,')'
    write(rfrmt_cont,'(A2,I2,A1,I1,A3)') '(F',clength,'.',flength,',$)'

    ifm = 1

    if (enefunc%num_bonds > 0) then
      write(category(ifm),frmt) 'BOND'
      values(ifm) = energy%bond
      ifm = ifm+1
    end if

    if (enefunc%num_angles > 0 .or. enefunc%num_ureys > 0) then
      write(category(ifm),frmt) 'ANGLE'
      values(ifm) = energy%angle
      ifm = ifm+1

      if (enefunc%forcefield == ForcefieldCHARMM) then
        write(category(ifm),frmt) 'UREY-BRADLEY'
        values(ifm) = energy%urey_bradley
        ifm = ifm+1
      end if
    end if

    if (enefunc%num_dihedrals > 0 .or. enefunc%num_rb_dihedrals > 0) then
      write(category(ifm),frmt) 'DIHEDRAL'
      values(ifm) = energy%dihedral
      ifm = ifm+1
    end if

    if (enefunc%num_impropers > 0 ) then
      write(category(ifm),frmt) 'IMPROPER'
      values(ifm) = energy%improper
      ifm = ifm+1
    end if

    if (enefunc%forcefield == ForcefieldCHARMM .or.  &
        enefunc%forcefield == ForcefieldCHARMM19) then

      if (enefunc%num_cmaps > 0 ) then

        write(category(ifm),frmt) 'CMAP'
        values(ifm) = energy%cmap
        ifm = ifm+1
      end if
    end if

    if (enefunc%forcefield == ForcefieldAAGO .or. &
        enefunc%forcefield == ForcefieldCAGO .or. &
        enefunc%forcefield == ForcefieldKBGO) then

      write(category(ifm),frmt) 'NATIVE_CONTACT'
      values(ifm) = energy%contact
      ifm = ifm+1

      write(category(ifm),frmt) 'NON-NATIVE_CONT'
      values(ifm) = energy%noncontact
      ifm = ifm+1
    else

      write(category(ifm),frmt) 'VDWAALS'
      values(ifm) = energy%van_der_waals
      ifm = ifm+1

      if (enefunc%dispersion_corr /= Disp_Corr_NONE) then
        write(category(ifm),frmt) 'DISP-CORR_ENE'
        values(ifm) = energy%disp_corr_energy
        ifm = ifm+1
      end if

      write(category(ifm),frmt) 'ELECT'
      values(ifm) = energy%electrostatic
      ifm = ifm+1

      if (enefunc%eef1_use .or. enefunc%gbsa_use) then
        write(category(ifm),frmt) 'SOLVATION'
        values(ifm) = energy%solvation
        ifm = ifm+1
      endif

    end if

    do i = 1, enefunc%num_restraintfuncs
      write(category(ifm),frmt_res) 'RESTRAINT', i
      values(ifm) = energy%restraint(i)
      ifm = ifm+1
    end do
    do i = 1, enefunc%num_restraintfuncs
      write(category(ifm),frmt_res) 'RESTR_CVS', i
      values(ifm) = energy%restraint_cv(i)
      ifm = ifm+1
    end do

    if (enefunc%steered_function > 0) then
      write(category(ifm),frmt) 'SMD_CV'
      values(ifm) = enefunc%target_value
      ifm = ifm+1
    end if
    if (enefunc%target_function > 0) then
      write(category(ifm),frmt) 'TMD_CV'
      values(ifm) = enefunc%target_value
      ifm = ifm+1
    end if

    if (enefunc%qmmm%do_qmmm) then
      write(category(ifm),frmt) 'QM'
      values(ifm) = energy%qm_ene  * CONV_UNIT_ENE
      ifm = ifm + 1
    end if

    if (enefunc%spot_use) then
      write(category(ifm),frmt) 'SPOT'
      values(ifm) = energy%spot
      ifm = ifm + 1
    end if

    if (enefunc%gamd_use) then
      if (enefunc%gamd%boost_pot) then
        write(category(ifm),frmt) 'POTENTIAL_GAMD'
        values(ifm) = energy%total_gamd
        ifm = ifm+1

        if (enefunc%gamd%boost_qm) then
          write(category(ifm),frmt) 'QM_GAMD' !added by ykl
          values(ifm) = energy%qm_gamd
          ifm = ifm+1
        end if
      else if (enefunc%gamd%boost_dih) then
        write(category(ifm),frmt) 'DIHEDRAL_GAMD'
        values(ifm) = energy%dihedral_gamd
        ifm = ifm+1
      else if (enefunc%gamd%boost_dual) then
        write(category(ifm),frmt) 'POTENTIAL_GAMD'
        values(ifm) = energy%total_gamd
        ifm = ifm+1
        write(category(ifm),frmt) 'DIHEDRAL_GAMD'
        values(ifm) = energy%dihedral_gamd
        ifm = ifm+1
        if (enefunc%gamd%boost_qm) then
          write(category(ifm),frmt) 'QM_GAMD'
          values(ifm) = energy%qm_gamd
          ifm = ifm+1
        end if
      else if (enefunc%gamd%boost_qm_only) then
        write(category(ifm),frmt) 'QM_GAMD'
        values(ifm) = energy%qm_gamd
        ifm = ifm+1
      end if
    end if

    if (etitle) then

      write(MsgOut,'(A,$)') title

      do i = 1, ifm-1

        if (i == ifm-1) then
          write(MsgOut,frmt) category(i)
        else
          write(MsgOut,frmt_cont) category(i)
        end if
      end do

      write(MsgOut,'(A80)') ' --------------- --------------- --------------- --------------- ---------------'
      etitle = .false.
    end if

    write(MsgOut,'(6x,I10,$)') step

    do i = 1, ifm-1
      if (i == ifm-1) then
        write(MsgOut,rfrmt) values(i)
      else
        write(MsgOut,rfrmt_cont) values(i)
      end if
    end do

    write(MsgOut,'(A)') ''

    return

  end subroutine output_energy_genesis

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_energy_charmm
  !> @brief        output energy in CHARMM style
  !! @authors      YS, CK
  !! @param[in]    step    : step
  !! @param[in]    enefunc : information of potential functions
  !! @param[in]    energy  : energy information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_energy_charmm(step, enefunc, energy)

    ! formal arguments
    integer,                 intent(in)    :: step
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_energy),  target, intent(in)    :: energy

    ! local variables
    real(wp)                 :: time, totener, totke, energy_, temperature
    real(wp)                 :: grms, hfctote, hfcke, ehfcor, virke
    real(wp)                 :: hbonds, asp, user
    real(wp)                 :: imnbvdw, imelec, imhbnd, rxnfield, extelec
    real(wp)                 :: ewksum, ewself, ewexcl, ewqcor, ewutil
    real(wp)                 :: vire, viri, presse, pressi, volume
    real(wp)                 :: cdihe, cintcr, noe
    real(wp)                 :: restdist, posicon
    integer                  :: i

    real(wp),        pointer :: bonds, angles, urey_b, dihedrals, impropers
    real(wp),        pointer :: vdwaals, elec
    real(wp),        pointer :: cmaps, contact, noncontact, disp_corr


    time        = 0.0_wp
    hbonds      = 0.0_wp
    totener     = 0.0_wp
    totke       = 0.0_wp
    energy_     = 0.0_wp
    temperature = 0.0_wp
    asp         = 0.0_wp
    user        = 0.0_wp
    imnbvdw     = 0.0_wp
    imelec      = 0.0_wp
    imhbnd      = 0.0_wp
    rxnfield    = 0.0_wp
    extelec     = 0.0_wp
    ewksum      = 0.0_wp
    ewself      = 0.0_wp
    ewexcl      = 0.0_wp
    ewqcor      = 0.0_wp
    ewutil      = 0.0_wp
    vire        = 0.0_wp
    viri        = 0.0_wp
    presse      = 0.0_wp
    pressi      = 0.0_wp
    volume      = 0.0_wp
    grms        = 0.0_wp
    hfctote     = 0.0_wp
    hfcke       = 0.0_wp
    ehfcor      = 0.0_wp
    virke       = 0.0_wp
    volume      = 0.0_wp
    noe         = 0.0_wp
    cdihe       = 0.0_wp
    cintcr      = 0.0_wp

    ! write title if necessary
    !
    if (etitle) then
      write(MsgOut,'(A)') 'Output_Energy> CHARMM_Style is used'
      write(MsgOut,'(A)') ' '
      write(MsgOut,'(A79)') 'DYNA DYN: Step         Time      TOTEner        TOTKe       ENERgy  TEMPerature'
      write(MsgOut,'(A79)') 'DYNA PROP:             GRMS      HFCTote        HFCKe       EHFCor        VIRKe'
      write(MsgOut,'(A79)') 'DYNA INTERN:          BONDs       ANGLes       UREY-b    DIHEdrals    IMPRopers'

      if (enefunc%num_cmaps > 0) then
        write(MsgOut,'(A79)') 'DYNA CROSS:           CMAPs                                                    '
      end if
      if (enefunc%dispersion_corr /= Disp_Corr_NONE) then
        write(MsgOut,'(A79)') 'DYNA  DISP:       Disp-Corr                                                    '
      end if

      if (enefunc%forcefield == ForcefieldAAGO .or. &
          enefunc%forcefield == ForcefieldCAGO .or. &
          enefunc%forcefield == ForcefieldKBGO) then
        write(MsgOut,'(A79)') 'DYNA GO:            CONTACT     NCONTACT                                       '
      else
        write(MsgOut,'(A79)') 'DYNA EXTERN:        VDWaals         ELEC       HBONds          ASP         USER'
      end if

      write(MsgOut,'(A79)') 'DYNA IMAGES:        IMNBvdw       IMELec       IMHBnd       RXNField    EXTElec'
      write(MsgOut,'(A79)') 'DYNA EWALD:          EWKSum       EWSElf       EWEXcl       EWQCor       EWUTil'

      if (enefunc%restraint_flag) then
        write(MsgOut,'(A79)') 'DYNA CONSTR:       HARMonic    CDIHedral          CIC     RESDistance       NOE'
      end if

      write(MsgOut,'(A79)') 'DYNA PRESS:            VIRE         VIRI       PRESSE       PRESSI       VOLUme'
      write(MsgOut,'(A79)') ' ----------       ---------    ---------    ---------    ---------    ---------'
      !if (nrep_per_proc == 1) vervose = .false.
      !if (.not. vervose) etitle = .false.
      !vervose = .false.

    end if

    ! use pointers
    !
    bonds      => energy%bond
    angles     => energy%angle
    urey_b     => energy%urey_bradley
    dihedrals  => energy%dihedral
    impropers  => energy%improper
    elec       => energy%electrostatic
    vdwaals    => energy%van_der_waals
    cmaps      => energy%cmap
    contact    => energy%contact
    noncontact => energy%noncontact
    disp_corr  => energy%disp_corr_energy

    ! Restraint energy
    !
    restdist = 0.0_wp
    posicon  = 0.0_wp
    do i = 1, enefunc%num_restraintfuncs
      select case (enefunc%restraint_kind(i))

      case (RestraintsFuncPOSI,                         &
            RestraintsFuncRMSD:RestraintsFuncRMSDCOM)

        posicon  = posicon  + energy%restraint(i)

      case (RestraintsFuncDIST:RestraintsFuncDISTCOM,   &
            RestraintsFuncANGLE:RestraintsFuncANGLECOM, &
            RestraintsFuncDIHED:RestraintsFuncDIHEDCOM, &
            RestraintsFuncREPUL:RestraintsFuncREPULCOM, &
            RestraintsFuncFB:RestraintsFuncFBCOM)
 
        restdist = restdist + energy%restraint(i)

      end select
    end do


    ! write energy in CHARMM-style
    !
    write(MsgOut,'(A5,I9,5F13.5)')   'DYNA>'         , step, time, totener, totke, energy_, temperature
    write(MsgOut,'(A14,5F13.5)')     'DYNA PROP>    ', grms, hfctote, hfcke, ehfcor, virke
    write(MsgOut,'(A14,5F13.5)')     'DYNA INTERN>  ', bonds, angles, urey_b, dihedrals, impropers

    if (enefunc%num_cmaps > 0) then
      write(MsgOut,'(A14, F13.5)')   'DYNA CROSS>   ', cmaps
    end if
    if (enefunc%dispersion_corr /= Disp_Corr_NONE) then
      write(MsgOut,'(A14,F13.5)')    'DYNA DISP>    ', disp_corr
    end if

    if (enefunc%forcefield == ForcefieldAAGO .or. &
        enefunc%forcefield == ForcefieldCAGO) then
      write(MsgOut,'(A14,2F13.5)')   'DYNA GO>      ', contact, noncontact
    else
      write(MsgOut,'(A14,5F13.5)')   'DYNA EXTERN>  ', vdwaals, elec, hbonds, asp, user
    end if

    write(MsgOut,'(A14,5F13.5)')     'DYNA IMAGES>  ', imnbvdw, imelec, imhbnd, rxnfield, extelec
    write(MsgOut,'(A14,5F13.5)')     'DYNA EWALD>   ', ewksum, ewself, ewexcl, ewqcor, ewutil

    if (enefunc%restraint_flag) then
      write(MsgOut,'(A14,5F13.5)')   'DYNA CONSTR>  ', posicon, cdihe, cintcr, restdist, noe
    end if

    write(MsgOut,'(A14,5F13.5)')     'DYNA PRESS>   ', vire, viri, presse, pressi, volume
    write(MsgOut,'(A79)')            ' ----------       ---------    ---------    ---------    ---------    ---------'
    write(MsgOut,'(A)')              ' '

    return

  end subroutine output_energy_charmm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_energy_namd
  !> @brief        output energy in NAMD style
  !! @authors      YS, CK
  !! @param[in]    step    : step
  !! @param[in]    enefunc : potential energy functions information
  !! @param[inout] energy  : energy information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_energy_namd(step, enefunc, energy)

    ! formal arguments
    integer,                 intent(in)    :: step
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_energy),  target, intent(inout) :: energy

    ! local variables
    real(wp)                 :: ebond, eangle
    real(wp)                 :: misc, kinetic
    real(wp)                 :: total, temp, total2, total3, tempavg
    real(wp)                 :: pressure, gpressure, volume
    real(wp)                 :: pressavg, gpressavg
    real(wp)                 :: restdist, posicon, eboundary
    integer                  :: i

    real(wp),        pointer :: edihed, eimprp
    real(wp),        pointer :: eelect, evdw


    ! write title if necessary
    !
    if (etitle) then
      write(MsgOut,'(A)') 'Output_Energy> NAMD_Style is used'
      write(MsgOut,'(A)') ' '
      write(MsgOut,'(A75)') 'ETITLE:      TS           BOND          ANGLE          DIHED          IMPRP'
      write(MsgOut,'(A75)') '          ELECT            VDW       BOUNDARY           MISC        KINETIC'
      write(MsgOut,'(A75)') '          TOTAL           TEMP         TOTAL2         TOTAL3        TEMPAVG'
      write(MsgOut,'(A75)') '       PRESSURE      GPRESSURE         VOLUME       PRESSAVG      GPRESSAVG'
      write(MsgOut,'(A)') ' '
      etitle = .false.
    end if

    ! restraint energy
    !
    restdist = 0.0_wp
    posicon  = 0.0_wp
    do i = 1, enefunc%num_restraintfuncs
      select case (enefunc%restraint_kind(i))

      case (RestraintsFuncPOSI,                         &
            RestraintsFuncRMSD:RestraintsFuncRMSDCOM)

        posicon = posicon + energy%restraint(i)

      case (RestraintsFuncDIST:RestraintsFuncDISTCOM,   &
            RestraintsFuncANGLE:RestraintsFuncANGLECOM, &
            RestraintsFuncDIHED:RestraintsFuncDIHEDCOM, &
            RestraintsFuncREPUL:RestraintsFuncREPULCOM, &
            RestraintsFuncFB:RestraintsFuncFBCOM)
 
        restdist = restdist + energy%restraint(i)

      end select
    end do

    ebond     =  energy%bond  + restdist
    eangle    =  energy%angle + energy%urey_bradley
    edihed    => energy%dihedral
    eimprp    => energy%improper
    eelect    => energy%electrostatic
    evdw      => energy%van_der_waals
    eboundary =  posicon

    ! write energy in NAMD-style
    !
    write(MsgOut,'(A7,I8,4F15.4)')'ENERGY:', step, ebond, eangle, edihed,eimprp
    write(MsgOut,'(5F15.4)') eelect, evdw, eboundary, misc, kinetic
    write(MsgOut,'(5F15.4)') total, temp, total2, total3, tempavg
    write(MsgOut,'(5F15.4)') pressure, gpressure, volume, pressavg, gpressavg
    write(MsgOut,'(A)') ' '

    return

  end subroutine output_energy_namd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_energy_gromacs
  !> @brief        output energy in GROMACS style
  !! @authors      NT
  !! @param[in]    step    : step
  !! @param[in]    enefunc : information of potential functions
  !! @param[in]    energy  : energy information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_energy_gromacs(step, enefunc, energy)

    ! formal arguments
    integer,                 intent(in)    :: step
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_energy),  target, intent(in)    :: energy

    ! local variables
    real(wp)                 :: time, bonds, angles, urey_b
    real(wp)                 :: dihedrals, impropers, vdwaals, elec
    real(wp)                 :: restdist, posicon
    real(wp)                 :: energy_, totke, totener
    real(wp)                 :: temperature, pressi, presse, disp_corr
    real(wp)                 :: noncontact, contact
    integer                  :: i


    time        = 0.0_wp
    bonds       = energy%bond          * CAL2JOU
    angles      = energy%angle         * CAL2JOU
    urey_b      = energy%urey_bradley  * CAL2JOU
    dihedrals   = energy%dihedral      * CAL2JOU
    impropers   = energy%improper      * CAL2JOU
    vdwaals     = energy%van_der_waals * CAL2JOU
    elec        = energy%electrostatic * CAL2JOU
    disp_corr   = energy%disp_corr_energy * CAL2JOU
    noncontact  = energy%noncontact    * CAL2JOU
    contact     = energy%contact       * CAL2JOU

    restdist    = 0.0_wp
    posicon     = 0.0_wp

    energy_     = 0.0_wp
    totke       = 0.0_wp
    totener     = 0.0_wp
    temperature = 0.0_wp
    pressi      = 0.0_wp
    presse      = 0.0_wp


    write(MsgOut,'(3A15)') &
         'Step', 'Time', 'Lambda'
    write(MsgOut,'(I15,2F15.5)') &
         step, time, 0.0_wp
    write(MsgOut,'(A)') &
         ' '
    write(MsgOut,'(A)') &
         '   Energies (kJ/mol)'

    write(MsgOut,'(5A15)') &
         'Bond', 'Angle', 'Urey-bradley', 'Dihedral', 'Improper Dih.'
    write(MsgOut,'(5ES15.5E2)') &
         bonds,angles,urey_b,dihedrals,impropers

    if (enefunc%forcefield == ForcefieldAAGO .or. &
             enefunc%forcefield == ForcefieldCAGO .or. &
             enefunc%forcefield == ForcefieldKBGO) then
      write(MsgOut,'(5A15)') &
           ' Contact  ', ' Noncontact    ', 'Position Rest.', 'Potential',  &
           'Kinetic En.'
      write(MsgOut,'(5ES15.5E2)') &
           contact,noncontact,posicon,energy_,totke
    else if (enefunc%dispersion_corr /= Disp_Corr_NONE) then
      write(MsgOut,'(5A15)') &
           'LJ (1-4,SR', ' Coulomb(1-4,SR', 'Position Rest.', 'Potential',  &
           'Kinetic En.'
      write(MsgOut,'(5ES15.5E2)') &
           vdwaals,elec,posicon,energy_,totke
    else
      write(MsgOut,'(4A15)') &
         'LJ (1-4,SR', ' Coulomb(1-4,SR', 'Disper. corr.', 'Position Rest.'
      write(MsgOut,'(4ES15.5E2)') &
           vdwaals,elec,disp_corr,posicon
      write(MsgOut,'(2A15)') &
         'Potential', 'Kinetic En.'
      write(MsgOut,'(2ES15.5E2)') &
          energy_,totke
    end if

    write(MsgOut,'(5A15)') &
         'Total Energy', 'Temperature', 'Pressure(int.)', 'Pressure(ext.)'
    write(MsgOut,'(5ES15.5E2)') &
         totener,temperature,pressi,presse

    write(MsgOut,'(A)')  ' '

    return

  end subroutine output_energy_gromacs

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_stats
  !> @brief        compute statistical quantities for RPATH
  !! @authors      YM
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_stats(enefunc)
    ! formal arguments
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k
    integer                  :: dimno_i, dimno_j
    integer                  :: atom_i, atom_j
    integer                  :: irep_num, irep_first
    real(dp)                 :: etmp, stmp
    real(dp)                 :: d(1:3)

!    if (replica_main_rank) then
    if (.not. replica_main_rank) return

    irep_first = 0
    i = enefunc%istart_restraint
    do while(i <= enefunc%iend_restraint .and. irep_first == 0) 
      if (enefunc%restraint_replica_index(i) > 0) irep_first = i
      i = i + 1
    end do
    if (irep_first == 0) return

    if(enefunc%restraint_kind(irep_first) == RestraintsFuncPOSI) then
      do i = 1, enefunc%stats_dimension
        enefunc%stats_force(i) = enefunc%stats_force(i) + &
           2.0_wp * enefunc%restraint_const(1,irep_first) * &
           enefunc%stats_delta(i)
      end do
    else
      irep_num = 0
      do i = 1, enefunc%iend_restraint
        if (enefunc%restraint_replica_index(i) > 0) then
          irep_num = irep_num + 1
          enefunc%stats_force(irep_num) = enefunc%stats_force(irep_num) + &
            2.0_wp * enefunc%restraint_const(1,i) * enefunc%stats_delta(irep_num)
        endif
      end do
    end if

    if (enefunc%restraint_kind(irep_first) == RestraintsFuncPOSI .or.  &
        enefunc%restraint_kind(irep_first) == RestraintsFuncPC .or.    &
        enefunc%restraint_kind(irep_first) == RestraintsFuncPCCOM) then
      do dimno_i = 1, enefunc%stats_dimension
        etmp = 0.0_dp
        do i = 1, enefunc%stats_natom
          d(1:3) = enefunc%stats_grad(1:3,i,dimno_i)
          do k = 1, 3
            etmp = etmp + (1.0_dp/enefunc%stats_mass(i,dimno_i))*d(k)*d(k)
          end do
        end do
        enefunc%stats_metric(dimno_i, dimno_i) =  &
           enefunc%stats_metric(dimno_i, dimno_i) + etmp
      end do
    else
      do dimno_i = 1, enefunc%stats_dimension
        do dimno_j = 1, enefunc%stats_dimension
          etmp = 0.0_dp
          do i = 1, enefunc%stats_natom
            atom_i = enefunc%stats_atom(i,dimno_i)
            stmp = (1.0_dp / enefunc%stats_mass(i,dimno_i)) 
            d(1:3) = enefunc%stats_grad(1:3,i,dimno_i)
            do j = 1, enefunc%stats_natom
              atom_j = enefunc%stats_atom(j,dimno_j)
              if (atom_i == atom_j) then
                do k = 1, 3
!                  enefunc%stats_metric(dimno_i,dimno_j) = &
!                    enefunc%stats_metric(dimno_i,dimno_j) &
!                    + (1.0_wp / enefunc%stats_mass(i,dimno_i)) &
!                    * enefunc%stats_grad(k,i,dimno_i) * enefunc%stats_grad(k,j,dimno_j)
                  etmp = etmp + stmp * d(k) * enefunc%stats_grad(k,j,dimno_j)
                end do
              end if
            end do
          end do
          enefunc%stats_metric(dimno_i, dimno_j) =  &
             enefunc%stats_metric(dimno_i, dimno_j) + etmp
        end do
      end do


    endif

    return

  end subroutine compute_stats

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    allreduce_eneforce
  !> @brief
  !! @authors      JJ, TM, CK
  !! @param[in]    enefunc    : potential energy functions information
  !! @param[in]    natom      : number of atoms
  !! @param[inout] energy     : energy information
  !! @param[inout] temporary  : temporary coordinates (used for MPI allreduce)
  !! @param[inout] force      : forces of target systems
  !! @param[inout] virial     : virial of target systems
  !! @param[inout] virial_ext : extended virial of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

#ifdef HAVE_MPI_GENESIS
  subroutine allreduce_eneforce(enefunc, natom, energy, temporary, force, &
                                virial, virial_ext)

    ! formal arguments
    type(s_enefunc),         intent(in)    :: enefunc
    integer,                 intent(in)    :: natom
    type(s_energy),          intent(inout) :: energy
    real(wp),                intent(inout) :: temporary(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: virial_ext(3,3)

    ! local variables
    real(wp)                 :: after_allreduce(999), before_allreduce(999)
    integer                  :: i, j, n
    integer                  :: ncycle, icycle, nlen, ixx


    ! Allreduce force
    !
    temporary(1:3,1:natom) = 0.0_wp
    do j = 1, natom
      temporary(1,j) = force(1,j)
      temporary(2,j) = force(2,j)
      temporary(3,j) = force(3,j)
    end do

    ncycle = (natom - 1) / mpi_drain + 1
    nlen   = mpi_drain
    ixx    = 1


    do icycle = 1, ncycle
      if (icycle == ncycle) nlen = natom - (ncycle-1) * mpi_drain
      call mpi_allreduce(temporary(1,ixx), force(1,ixx), 3*nlen, &
                      mpi_wp_real, mpi_sum, mpi_comm_country, ierror)
      ixx = ixx + nlen
    end do


    ! Allreduce virial and energy components
    !
    n = 0
    do i = 1, 3
      do j = 1, 3
        n = n + 1
        before_allreduce(n) = virial(i,j)
      end do
    end do

    do i = 1, 3
      do j = 1, 3
        n = n + 1
        before_allreduce(n) = virial_ext(i,j)
      end do
    end do

    before_allreduce(19) = energy%bond
    before_allreduce(20) = energy%angle
    before_allreduce(21) = energy%urey_bradley
    before_allreduce(22) = energy%dihedral
    before_allreduce(23) = energy%improper
    before_allreduce(24) = energy%cmap
    before_allreduce(25) = energy%electrostatic
    before_allreduce(26) = energy%van_der_waals
    before_allreduce(27) = energy%contact
    before_allreduce(28) = energy%noncontact
    before_allreduce(29) = energy%solvation
    before_allreduce(30) = energy%spot

    call mpi_allreduce(before_allreduce, after_allreduce, 30, &
                       mpi_wp_real,  mpi_sum,                 &
                       mpi_comm_country, ierror)

    n = 0
    do i = 1, 3
      do j = 1, 3
        n = n + 1
        virial(i,j) = after_allreduce(n)
      end do
    end do

    do i = 1, 3
      do j = 1, 3
        n = n + 1
        virial_ext(i,j) = after_allreduce(n)
      end do
    end do

    energy%bond               = after_allreduce(19)
    energy%angle              = after_allreduce(20)
    energy%urey_bradley       = after_allreduce(21)
    energy%dihedral           = after_allreduce(22)
    energy%improper           = after_allreduce(23)
    energy%cmap               = after_allreduce(24)
    energy%electrostatic      = after_allreduce(25)
    energy%van_der_waals      = after_allreduce(26)
    energy%contact            = after_allreduce(27)
    energy%noncontact         = after_allreduce(28)
    energy%solvation          = after_allreduce(29)
    energy%spot               = after_allreduce(30)

    return

  end subroutine allreduce_eneforce
#endif

end module at_energy_mod
