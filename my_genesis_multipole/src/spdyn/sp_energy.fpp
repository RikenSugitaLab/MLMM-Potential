!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_energy_mod
!> @brief   compute energy
!! @authors Jaewoon Jung(JJ), Yuji Sugita (YS), Takaharu Mori (TM),
!!          Chigusa Kobayashi (CK), Norio Takase (NT)
!  
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_energy_mod

  use sp_energy_table_linear_mod          
  use sp_energy_table_linear_bondcorr_mod 
  use sp_energy_nonbonds_mod
  use sp_energy_restraints_mod
  use sp_energy_dihedrals_mod
  use sp_energy_angles_mod
  use sp_energy_bonds_mod
  use sp_energy_go_mod
  use sp_energy_gamd_mod
  use sp_boundary_str_mod
  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_energy_str_mod
  use sp_restraints_str_mod
  use sp_domain_str_mod
  use sp_experiments_mod
  use sp_enefunc_gamd_mod
  use fileio_control_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use math_libs_mod
  use string_mod
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
    real(wp)              :: dielec_const     = 1.0_wp
    real(wp)              :: dmin_size_cg     = 20.0_wp
    logical               :: vdw_force_switch = .false.
    logical               :: vdw_shift        = .false.
    logical               :: cmap_pspline     = .false.
    logical               :: contact_check    = .false.
    logical               :: nonb_limiter     = .false.
    integer               :: structure_check  = StructureCheckNone
    logical               :: dsize_cg         = .false.
    real(wp)              :: pme_alpha ! do not initialize here
    real(wp)              :: pme_alpha_tol    = 1.0e-5
    integer               :: pme_ngrid_x      = FakeDefault
    integer               :: pme_ngrid_y      = FakeDefault
    integer               :: pme_ngrid_z      = FakeDefault
    integer               :: pme_nspline      = 4
    real(wp)              :: pme_max_spacing  = 1.2_wp
    integer               :: fft_scheme       = FFT_1dalltoall
    logical               :: table            = .true.
    integer               :: table_order      = 1
    real(wp)              :: table_density    = 20.0_wp
    character(5)          :: water_model      = 'NONE'
    integer               :: output_style     = OutputStyleGENESIS
    integer               :: dispersion_corr  = Disp_corr_NONE
    real(wp)              :: minimum_contact  = 0.5_wp
    real(wp)              :: err_minimum_contact = 0.3_wp
    logical               :: vacuum           = .false.
  end type s_ene_info

  ! varibles
  logical, save           :: etitle = .true.
  !FEP
  logical, save           :: etitle_fep = .true.

  ! subroutines
  public  :: show_ctrl_energy
  public  :: read_ctrl_energy
  public  :: compute_energy
  public  :: compute_energy_short
  public  :: compute_energy_long
  public  :: output_energy
  private :: compute_energy_charmm
  private :: compute_energy_amber
  private :: compute_energy_gro_amber
  private :: compute_energy_gro_martini
  private :: compute_energy_charmm_short
  private :: compute_energy_amber_short
  private :: compute_energy_gro_amber_short
  private :: compute_energy_general_long
  private :: output_energy_genesis
  private :: output_energy_charmm
  private :: output_energy_namd
  private :: output_energy_gromacs
  private :: reduce_ene
  private :: compute_stats
  !FEP
  public  :: compute_energy_fep
  public  :: output_energy_fep
  private :: compute_energy_charmm_fep
  private :: compute_energy_amber_fep
  private :: compute_energy_gro_amber_fep
  private :: output_energy_genesis_fep
  public  :: compute_energy_short_fep
  private :: compute_energy_charmm_short_fep
  private :: compute_energy_amber_short_fep
  private :: compute_energy_gro_amber_short_fep

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_energy
  !> @brief        show ENERGY section usage
  !! @authors      NT
  !! @param[in]    show_all : show all usage or not
  !! @param[in]    run_mode : run mode string : "md", "min"
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_energy(show_all, run_mode)

    ! formal arguments
    logical,                 intent(in)    :: show_all
    character(*),            intent(in)    :: run_mode


    if (show_all) then

      select case (run_mode)

      case ('md', 'min', 'remd', 'rpath')

        write(MsgOut,'(A)') '[ENERGY]'
        write(MsgOut,'(A)') 'forcefield    = CHARMM    # [CHARMM,AMBER,GROAMBER,GROMARTINI]'
        write(MsgOut,'(A)') 'electrostatic = PME       # [CUTOFF,PME]'
        write(MsgOut,'(A)') 'switchdist    = 10.0      # switch distance'
        write(MsgOut,'(A)') 'cutoffdist    = 12.0      # cutoff distance'
        write(MsgOut,'(A)') 'pairlistdist  = 13.5      # pair-list distance'
        write(MsgOut,'(A)') '# dielec_const  = 1.0       # dielectric constant'
        write(MsgOut,'(A)') '# vdw_force_switch = NO     # force switch option for van der Waals'
        write(MsgOut,'(A)') '# vdw_shift        = NO     # shift option for van der Waals'
        write(MsgOut,'(A)') '# pme_alpha     = auto      # width of Gaussian distribution in [PME]'
        write(MsgOut,'(A)') '# pme_alpha_tol = 1.0e-5    # param for auto pme_alpha determination'
        write(MsgOut,'(A)') '# pme_nspline   = 4         # order of B-spline in [PME]'
        write(MsgOut,'(A)') '# pme_max_spacing  = 1.2    # Max grid spacing allowed '
        write(MsgOut,'(A)') '# fft_scheme    = 1dalltoall# FFT parallel scheme [1dalltoall, 2dalltoall, 1dallgather]'
        write(MsgOut,'(A)') '# table_density = 20.0      # number of bins used for lookup table'
        write(MsgOut,'(A)') '# output_style  = GENESIS   # format of energy output [GENESIS,CHARMM,NAMD,GROMACS]'
        write(MsgOut,'(A)') '# dispersion_corr = NONE    # dispersion correction [NONE,Energy,EPress]'
        write(MsgOut,'(A)') '# vacuum = NO               # vacuum option'
        if (run_mode == 'min') then
          write(MsgOut,'(A)') '# contact_check   = YES     # check atomic clash'
          write(MsgOut,'(A)') '# nonb_limiter    = YES     # avoid failure due to atomic clash'
          write(MsgOut,'(A)') '# minimum_contact = 0.5     # definition of atomic clash distance'
        end if
        write(MsgOut,'(A)') ' '

      end select

    else

      select case (run_mode)

      case ('md', 'min', 'remd', 'rpath')

        write(MsgOut,'(A)') '[ENERGY]'
        write(MsgOut,'(A)') 'forcefield    = CHARMM    # [CHARMM,AMBER,GROAMBER,GROMARTINI]'
        write(MsgOut,'(A)') 'electrostatic = PME       # [CUTOFF,PME]'
        write(MsgOut,'(A)') 'switchdist    = 10.0      # switch distance'
        write(MsgOut,'(A)') 'cutoffdist    = 12.0      # cutoff distance'
        write(MsgOut,'(A)') 'pairlistdist  = 13.5      # pair-list distance'
        write(MsgOut,'(A)') ' '


      end select

    end if

    return

  end subroutine show_ctrl_energy

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      read_ctrl_energy
  !> @brief        read ENERGY section in the control file
  !! @authors      YS, TM, JJ
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
    real(wp)                               :: cutoff, cutoff2, mind
    real(wp)                               :: cutoff2_water
    integer                                :: cutoff_int2


    ! read parameters from control file
    !
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_type   (handle, Section, 'forcefield',    &
                               ene_info%forcefield, ForceFieldTypes)
    call read_ctrlfile_string (handle, Section, 'forcefield',    &
                               ene_info%forcefield_char)
    call read_ctrlfile_type   (handle, Section, 'electrostatic', &
                               ene_info%electrostatic, ElectrostaticTypes)
    call read_ctrlfile_real   (handle, Section, 'switchdist',    &
                               ene_info%switchdist)
    call read_ctrlfile_real   (handle, Section, 'cutoffdist',    &
                               ene_info%cutoffdist)
    call read_ctrlfile_real   (handle, Section, 'pairlistdist',  &
                               ene_info%pairlistdist)
    call read_ctrlfile_real   (handle, Section, 'dmin_size_cg',  &
                               ene_info%dmin_size_cg)
    call read_ctrlfile_real   (handle, Section, 'dielec_const',  &
                               ene_info%dielec_const)
    call read_ctrlfile_logical(handle, Section, 'vdw_force_switch', &
                               ene_info%vdw_force_switch)
    call read_ctrlfile_logical(handle, Section, 'vdw_shift',     &
                               ene_info%vdw_shift)
    call read_ctrlfile_string (handle, Section, 'pme_alpha',     &
                               pme_alpha)
    call read_ctrlfile_real   (handle, Section, 'pme_alpha_tol', &
                               ene_info%pme_alpha_tol)
    call read_ctrlfile_integer(handle, Section, 'pme_ngrid_x',   &
                               ene_info%pme_ngrid_x)
    call read_ctrlfile_integer(handle, Section, 'pme_ngrid_y',   &
                               ene_info%pme_ngrid_y)
    call read_ctrlfile_integer(handle, Section, 'pme_ngrid_z',   &
                               ene_info%pme_ngrid_z)
    call read_ctrlfile_integer(handle, Section, 'pme_nspline',   &
                               ene_info%pme_nspline)
    call read_ctrlfile_real   (handle, Section, 'pme_max_spacing', &
                               ene_info%pme_max_spacing)
    call read_ctrlfile_type   (handle, Section, 'FFT_scheme',    &
                               ene_info%fft_scheme, FFT_Types)
    call read_ctrlfile_real   (handle, Section, 'table_density', &
                               ene_info%table_density)
    call read_ctrlfile_string (handle, Section, 'water_model',   &
                               ene_info%water_model)
    call read_ctrlfile_type   (handle, Section, 'output_style',  &
                               ene_info%output_style, OutputStyleTypes)
    call read_ctrlfile_type   (handle, Section, 'dispersion_corr',  &
                               ene_info%dispersion_corr, Disp_corr_Types)
    call read_ctrlfile_type   (handle, Section, 'structure_check',  &
                               ene_info%structure_check, StructureCheckTypes)
    call read_ctrlfile_logical(handle, Section, 'contact_check',  &
                               ene_info%contact_check)
    if (ene_info%contact_check) then
      ene_info%nonb_limiter=.true.
    endif
    call read_ctrlfile_logical(handle, Section, 'nonb_limiter',  &
                               ene_info%nonb_limiter)
    call read_ctrlfile_real   (handle, Section, 'minimum_contact', &
                               ene_info%minimum_contact)
    call read_ctrlfile_logical(handle, Section, 'vacuum',     &
                               ene_info%vacuum)
    call end_ctrlfile_section(handle)


    ! check vacuum
    !
    if (ene_info%vacuum) then
      ! Tables are not used.
      ene_info%table = .false.
      ! Switch functions are not used.
      ene_info%vdw_force_switch = .false.
    end if

    ! check table
    !
    if (ene_info%table) then
      if (ene_info%electrostatic == ElectrostaticCutoff ) &
          ene_info%table_order = 3
      if (ene_info%electrostatic == ElectrostaticPME ) &
          ene_info%table_order = 1
    endif

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

    if (ene_info%water_model /= "NONE" .and. &
        ene_info%electrostatic == ElectrostaticCUTOFF) then
        if (main_rank)      &
        write(MsgOut,'(A)') &
          'Read_Ctrl_Energy>  WARNING: water_model is not available'

          ene_info%water_model = "NONE"

    endif

#ifdef USE_GPU
    if (ene_info%electrostatic == ElectrostaticCUTOFF) &
      call error_msg( &
         'Read_Ctrl_Energy> cutoff is not available with GPU')

    if (ene_info%water_model /= "NONE") ene_info%water_model = 'NONE'

    if (ene_info%nonb_limiter) &
      call error_msg( &
         'Read_Ctrl_Energy> nonb_limiter is not available with GPU')

    if (ene_info%structure_check /= StructureCheckNone) &
      call error_msg( &
         'Read_Ctrl_Energy> structure_check is not available with GPU')
#endif

    call tolower(pme_alpha)
    if (trim(pme_alpha) == "auto") then
      ene_info%pme_alpha = get_ewald_alpha(ene_info%cutoffdist, &
                                           ene_info%pme_alpha_tol)
    else
      read(pme_alpha,*) ene_info%pme_alpha
    end if

    ! error check for each FFs
    !
    if (ene_info%forcefield == ForcefieldAMBER) then
!      if (ene_info%electrostatic /= ElectrostaticPME) then
!        call error_msg('Read_Ctrl_Energy> CUTOFF is not allowed in amber')
!      endif

      if (ene_info%switchdist /= ene_info%cutoffdist) then
        if (main_rank)      &
        write(MsgOut,'(A)') &
          'Read_Ctrl_Energy>  WARNING: switchdist should be equal to '//  &
          'cutoffdist if forcefield is AMBER'
          ene_info%switchdist = ene_info%cutoffdist
      endif
      if (ene_info%dispersion_corr /= Disp_corr_EPress) then
        if (main_rank)      &
        write(MsgOut,'(A)') &
          'Read_Ctrl_Energy>  WARNING: dispersion correction should be '//&
          'set ene_info%cutoffdist if forcefield is AMBER'
          ene_info%dispersion_corr = Disp_corr_EPress
      endif
    endif

    if (ene_info%forcefield == ForcefieldCHARMM) then
      if (ene_info%dispersion_corr /= Disp_corr_NONE) then
        if (main_rank)      &
        write(MsgOut,'(A)') &
          'Read_Ctrl_Energy>  WARNING: dispersion correction can not be'//&
          ' set ene_info%cutoffdist if forcefield is CHARMM'
          ene_info%dispersion_corr = Disp_corr_NONE
      endif

    endif

    if (ene_info%forcefield == ForcefieldGROMARTINI) then
      ene_info%dsize_cg  = .true.
      if (ene_info%dmin_size_cg < ene_info%pairlistdist ) then
        if (main_rank)      &
          write(MsgOut,'(A)') &
          'Read_Ctrl_Energy>  WARNING: dmin_size_cg is automatically'//     &
          ' set to pairlistdist if forcefield is GROMARTINI and '//         &
          'dmin_size_cg < pairlistdist'
        ene_info%dmin_size_cg = ene_info%pairlistdist
      endif
      if (ene_info%water_model /= "NONE") then
        call error_msg( &
         'Read_Ctrl_Energy> water_model is not allowed, if forcefield'//  &
         ' is GROMARTINI ')
      endif
      if (ene_info%electrostatic /= ElectrostaticCUTOFF) then
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
    endif

    if (ene_info%forcefield /= ForcefieldCHARMM) then
      if (ene_info%vdw_force_switch) then
        if (main_rank)      &
        write(MsgOut,'(A)') &
          'Read_Ctrl_Energy>  WARNING: vdW force switch can not be'//     &
          ' set if forcefield is not CHARMM'
          ene_info%vdw_force_switch = .false.
      endif
    endif

    if (ene_info%forcefield /= ForcefieldGROAMBER .and. &
        ene_info%forcefield /= ForcefieldGROMARTINI) then
      if (ene_info%vdw_shift) then
        if (main_rank)      &
        write(MsgOut,'(A)') &
          'Read_Ctrl_Energy>  WARNING: vdW shift can not be set'//        &
          ' if forcefield is not GROAMBER'
          ene_info%vdw_shift = .false.
      endif
    endif

    ! error check
    !

    if (.not. ene_info%table) then
      if (ene_info%nonb_limiter ) then
        if (main_rank) &
          write(MsgOut,'(A)') &
       'Read_Ctrl_Energy>  WARNING: nonb_limiter is only in table'
        ene_info%nonb_limiter = .false.
      endif
    endif

    if (ene_info%nonb_limiter .and. .not. ene_info%contact_check) then
      if (main_rank)       &
      write(MsgOut,'(A)') &
       'Read_Ctrl_Energy>  WARNING: contact_check is available in '//&
       'nonb_limiter=YES'
      ene_info%contact_check = .true.
    endif
    if (ene_info%contact_check .and.                            &
        ene_info%structure_check == StructureCheckNone) then
      ene_info%structure_check = StructureCheckFirst
    endif

    ! error check for CUTOFF and vacuum
    !
    if (ene_info%vacuum) then
      if (ene_info%electrostatic /= ElectrostaticCutoff) then
        call error_msg( &
         'Read_Ctrl_Energy> vacuum is not allowed unless CUTOFF is used')
      end if
    end if

    ! write parameters to MsgOut
    !
    if (main_rank) then

      write(MsgOut,'(A)')  'Read_Ctrl_Energy> Parameters of Energy Calculations'
      write(MsgOut,'(A20,A10)')                            &
            '  forcefield      = ',                        &
                trim(ForceFieldTypes(ene_info%forcefield))

      write(MsgOut,'(A20,F10.3,A20,F10.3)')             &
            '  switchdist      = ', ene_info%switchdist,          &
            '  cutoffdist      = ', ene_info%cutoffdist         

      if (ene_info%electrostatic == ElectrostaticCUTOFF .and. &
           ene_info%vdw_shift) &
        write(MsgOut,'(A20,F10.3)') '  switchdist(ele) = ', 0.0_wp

      write(MsgOut,'(A20,F10.3,A20,F10.3)')                &
            '  pairlistdist    = ', ene_info%pairlistdist, &
            '  dielec_const    = ', ene_info%dielec_const
     if (ene_info%forcefield == ForcefieldGROMARTINI) then
       write(MsgOut,'(A20,F10.3)')                &
            '  dmin_size_cg    = ', ene_info%dmin_size_cg
     endif

      if (ene_info%vdw_force_switch) then
        write(MsgOut,'(A)') ' vdW force_switch =        yes'
      else
        if (ene_info%forcefield == ForcefieldGROAMBER .or. &
            ene_info%forcefield == ForcefieldGROMARTINI ) then
          if (ene_info%vdw_shift) then
            if (ene_info%forcefield == ForcefieldGROMARTINI .and. &
               ene_info%electrostatic == ElectrostaticCUTOFF) then
              write(MsgOut,'(A)') '  vdW/elec  shift =       yes'
            endif
          else 
            write(MsgOut,'(A)') '  vdW switch      =       yes'
          endif
        else
          write(MsgOut,'(A)') ' vdW force_switch =         no'
        endif
      end if

      ! if PME
      if (ene_info%electrostatic == ElectrostaticPME) then

        write(MsgOut,'(A20,A10)')                            &
              '  electrostatic   = ',                        &
              trim(ElectrostaticTypes(ene_info%electrostatic))
        if (ene_info%pme_ngrid_x /= FakeDefault .or.  &
            ene_info%pme_ngrid_y /= FakeDefault .or.  &
            ene_info%pme_ngrid_z /= FakeDefault) then
          write(MsgOut,'(A20,3I10)')                           &
                '  pme_ngrid(x,y,z)= ', ene_info%pme_ngrid_x   &
                                      , ene_info%pme_ngrid_y   &
                                      , ene_info%pme_ngrid_z
        endif
        write(MsgOut,'(A20,I10)')                            &
              '  pme_nspline     = ', ene_info%pme_nspline
        if (trim(pme_alpha) == "auto") then
          write(MsgOut,'(A20,F10.5,A20,E10.3)')                &
              '  pme_alpha       = ', ene_info%pme_alpha,    &
              '  pme_alpha_tol   = ', ene_info%pme_alpha_tol
        else
          write(MsgOut,'(A20,F10.5)')                &
              '  pme_alpha       = ', ene_info%pme_alpha
        endif
        write(MsgOut,'(A20,A10)')                            &
              '  fft_scheme      = ', trim(FFT_types(ene_info%fft_scheme))

      ! if CUTOFF
      else if (ene_info%electrostatic == ElectrostaticCutoff) then

        write(MsgOut,'(A20,A10)')                            &
              '  electrostatic   = ',                        &
              trim(ElectrostaticTypes(ene_info%electrostatic))

      end if

      if (ene_info%table) then
        write(MsgOut,'(A20,I10)')   '  table_order     = ', &
             ene_info%table_order
        write(MsgOut,'(A20,F10.3)') '  table_density   = ', &
             ene_info%table_density
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
      if (ene_info%contact_check) then
        write(MsgOut,'(A)') '  contact_check   =     yes'
        write(MsgOut,'(A,F10.3)') '  minimum_contact = ', &
             ene_info%minimum_contact
      else
        write(MsgOut,'(A)') '  contact_check   =      no'
      endif
      if (ene_info%structure_check /= StructureCheckNone) then
        write(MsgOut,'(A20,A6)')                             &
              '  structure_check = ',                        &
              trim(StructureCheckTypes(ene_info%structure_check))
      endif

      ! if vacuum
      if (ene_info%vacuum) then
        write(MsgOut,'(A)') '  vacuum          =     yes'
      else
        write(MsgOut,'(A)') '  vacuum          =      no'
      end if

      write(MsgOut,'(A)') ' '

    end if
    ene_info%minimum_contact=ene_info%minimum_contact*ene_info%minimum_contact

    if (ene_info%table_order == 1) then
      cutoff = ene_info%cutoffdist
      cutoff2 = cutoff*cutoff
      cutoff2_water  = (cutoff+4.5_wp) * (cutoff+4.5_wp)
      cutoff_int2    = int(cutoff2_water*ene_info%table_density)
      mind=cutoff2*ene_info%table_density/real(cutoff_int2,wp)+0.001_wp
      ene_info%err_minimum_contact=mind
    endif

    ene_info%err_minimum_contact=ene_info%err_minimum_contact* \
                                 ene_info%err_minimum_contact

    return
  
  end subroutine read_ctrl_energy

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy
  !> @brief        compute potential energy
  !! @authors      JJ
  !! @param[in]    domain        : domain information
  !! @param[in]    enefunc       : potential energy functions information
  !! @param[in]    pairlist      : pair list information
  !! @param[in]    boundary      : boundary information
  !! @param[in]    coord         : coordinates of target systems
  !! @param[in]    npt           : flag for NPT or not
  !! @param[in]    reduce        : flag for reduce energy and virial
  !! @param[in]    nonb_ene      : flag for calculate nonbonded energy
  !! @param[in]    nonb_limiter  : flag for nonbond limiter
  !! @param[inout] energy        : energy information
  !! @param[input] coord_pbc     : !TODO
  !! @param[inout] force         : forces of target systems
  !! @param[inout] force_omp     : temprary forces of target systems
  !! @param[inout] force_pbc     : !TODO
  !! @param[inout] virial_cell   : !TODO
  !! @param[inout] virial        : virial term of target systems
  !! @param[inout] virial_ext    : extern virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy(domain, enefunc, pairlist, boundary, coord,  &
                            npt, reduce, nonb_ene, merge_force,          &
                            nonb_limiter, energy,  &
                            coord_pbc, force, force_long, force_omp,     &
                            force_pbc, virial_cell, virial, virial_long, &
                            virial_ext)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(dp),                intent(in)    :: coord(:,:,:)
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: reduce
    logical,                 intent(in)    :: nonb_ene
    logical,                 intent(in)    :: merge_force
    logical,                 intent(in)    :: nonb_limiter
    type(s_energy),          intent(inout) :: energy
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(dp),                intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: force_long(:,:,:)
    real(wp),                intent(inout) :: force_omp(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3)
    real(dp),                intent(inout) :: virial_long(3,3)
    real(dp),                intent(inout) :: virial_ext(3,3)

    ! local variables
    real(dp)                 :: volume


    call timer(TimerEnergy, TimerOn)


    select case (enefunc%forcefield)

    case (ForcefieldCHARMM)

      if (enefunc%gamd_use) then 

        call compute_energy_charmm( &
                              domain, enefunc, pairlist, boundary, coord, &
                              npt, reduce, .true., merge_force,         &
                              nonb_limiter,                               &
                              energy, coord_pbc, force, force_long,       &
                              force_omp, force_pbc, virial_cell, virial,  &
                              virial_long, virial_ext)

      else

        call compute_energy_charmm( &
                              domain, enefunc, pairlist, boundary, coord, &
                              npt, reduce, nonb_ene, merge_force,         &
                              nonb_limiter,                               &
                              energy, coord_pbc, force, force_long,       &
                              force_omp, force_pbc, virial_cell, virial,  &
                              virial_long, virial_ext)

      end if

    case (ForcefieldAMBER)

      if (enefunc%gamd_use) then 

        call compute_energy_amber( &
                              domain, enefunc, pairlist, boundary, coord, &
                              npt, reduce, .true., merge_force,         &
                              nonb_limiter,                               &
                              energy, coord_pbc, force, force_long,       &
                              force_omp, force_pbc, virial_cell, virial,  &
                              virial_long, virial_ext)

      else

        call compute_energy_amber( &
                              domain, enefunc, pairlist, boundary, coord, &
                              npt, reduce, nonb_ene, merge_force,         &
                              nonb_limiter,                               &
                              energy, coord_pbc, force, force_long,       &
                              force_omp, force_pbc, virial_cell, virial,  &
                              virial_long, virial_ext)

      end if

    case (ForcefieldGROAMBER)

      call compute_energy_gro_amber( &
                              domain, enefunc, pairlist, boundary, coord, &
                              npt, reduce, nonb_ene, merge_force,         &
                              nonb_limiter,                               &
                              energy, coord_pbc, force, force_long,       &
                              force_omp, force_pbc, virial_cell, virial,  &
                              virial_long, virial_ext)


    case (ForcefieldGROMARTINI)

      call compute_energy_gro_martini( &
                              domain, enefunc, pairlist, boundary, coord, &
                              npt, reduce, nonb_ene, nonb_limiter,        &
                              energy, coord_pbc,                          &
                              force, force_long, force_omp, force_pbc,    &
                              virial_cell, virial, virial_ext)

    case (ForcefieldAAGO)

      call compute_energy_go(domain, enefunc, pairlist, boundary,  &
                             coord, reduce, nonb_ene, energy,      &
                             force, force_long, force_omp,         &
                             virial_long, virial, virial_ext)

    end select

    ! Dispersion correction
    if (enefunc%dispersion_corr /= Disp_corr_NONE) then
      volume =  boundary%box_size_x_ref * &
                boundary%box_size_y_ref * &
                boundary%box_size_z_ref
      energy%disp_corr_energy = enefunc%dispersion_energy/volume
      if (enefunc%dispersion_corr == Disp_corr_EPress) then
        energy%disp_corr_virial = enefunc%dispersion_virial/volume
        if (replica_main_rank .or. main_rank) then
          virial(1,1) = virial(1,1) + real(energy%disp_corr_virial,dp)
          virial(2,2) = virial(2,2) + real(energy%disp_corr_virial,dp)
          virial(3,3) = virial(3,3) + real(energy%disp_corr_virial,dp)
        endif
      end if

    end if

    if (enefunc%rpath_sum_mf_flag) then
      call compute_stats(enefunc)
    end if

    call timer(TimerEnergy, TimerOff)

    return

  end subroutine compute_energy

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_short
  !> @brief        compute potential energy with short range interaction
  !! @authors      JJ
  !! @param[in]    domain      : domain information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[in]    pairlist    : pair list information
  !! @param[in]    boundary    : boundary information
  !! @param[in]    coord       : coordinates of target systems
  !! @param[in]    nonb_ene    : flag for calculate nonbonded energy
  !! @param[inout] energy      : energy information
  !! @param[input] coord_pbc   : !TODO
  !! @param[inout] force       : forces of target systems
  !! @param[inout] force_omp   : temprary forces of target systems
  !! @param[inout] force_pbc   : !TODO
  !! @param[inout] virial_cell : !TODO
  !! @param[inout] virial      : virial term of target systems
  !! @param[inout] virial_ext  : extern virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_short(domain, enefunc, pairlist, boundary, coord, &
                                  npt, nonb_ene, energy, coord_pbc, &
                                  force, force_omp, force_pbc, &
                                  virial_cell, virial, virial_ext)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(dp),                intent(in)    :: coord(:,:,:)
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: nonb_ene
    type(s_energy),          intent(inout) :: energy
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(dp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: force_omp(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3)
    real(dp),                intent(inout) :: virial_ext(3,3)


    call timer(TimerEnergy, TimerOn)


    select case (enefunc%forcefield)

    case (ForcefieldCHARMM)

      call compute_energy_charmm_short( &
                              domain, enefunc, pairlist, boundary, coord, &
                              npt, nonb_ene, energy, coord_pbc,           &
                              force, force_omp, force_pbc,                &
                              virial_cell, virial, virial_ext)

    case (ForcefieldAMBER)

      call compute_energy_amber_short( &
                              domain, enefunc, pairlist, boundary, coord, &
                              npt, nonb_ene, energy, coord_pbc,           &
                              force, force_omp, force_pbc,                &
                              virial_cell, virial, virial_ext)

    case (ForcefieldGROAMBER)

      call compute_energy_gro_amber_short( &
                              domain, enefunc, pairlist, boundary, coord, &
                              npt, nonb_ene, energy, coord_pbc,           &
                              force, force_omp, force_pbc,                &
                              virial_cell, virial, virial_ext)

    case (ForcefieldGROMARTINI)

      !TODO

    end select


    call timer(TimerEnergy, TimerOff)

    return

  end subroutine compute_energy_short

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_long
  !> @brief        compute potential energy
  !! @authors      JJ
  !! @param[in]    domain    : domain information
  !! @param[in]    enefunc   : potential energy functions information
  !! @param[in]    pairlist  : pair list information
  !! @param[in]    boundary  : boundary information
  !! @param[in]    coord     : coordinates of target systems
  !! @param[in]    npt       : flag for NPT or not
  !! @param[inout] energy    : energy information
  !! @param[inout] force     : forces of target systems
  !! @param[inout] force_omp : temprary forces of target systems
  !! @param[inout] virial    : virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_long(domain, enefunc, pairlist, boundary, coord,  &
                                 npt, energy, force, force_omp, virial)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(dp),                intent(in)    :: coord(:,:,:)
    logical,                 intent(in)    :: npt
    type(s_energy),          intent(inout) :: energy
    real(dp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: force_omp(:,:,:,:)
    real(dp),                intent(inout) :: virial(3,3)


    call timer(TimerEnergy, TimerOn)

    call compute_energy_general_long( &
                              domain, enefunc, pairlist, boundary, coord,  &
                              npt, energy, force, virial)

    call timer(TimerEnergy, TimerOff)

    return

  end subroutine compute_energy_long

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_energy
  !> @brief        output energy
  !! @authors      YS
  !! @param[in]    step    : step 
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    energy  : energy information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_energy(step, enefunc, energy)

    ! formal arguments
    integer,                 intent(in)    :: step
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_energy),          intent(in)    :: energy


    if (.not. main_rank) return

    select case (enefunc%output_style)

    case (OutputStyleGENESIS)

      call output_energy_genesis(step, enefunc, energy)

    case (OutputStyleCHARMM)

      call output_energy_charmm(step, enefunc, energy)

    case (OutputStyleNAMD)

      call output_energy_namd(step, energy)

    case (OutputStyleGROMACS)

      call output_energy_gromacs(step, energy)

    end select

    return

  end subroutine output_energy

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_charmm
  !> @brief        compute potential energy with charmm force field
  !! @authors      JJ
  !! @param[in]    domain        : domain information
  !! @param[in]    enefunc       : potential energy functions information
  !! @param[in]    pairlist      : pair list information
  !! @param[in]    boundary      : boundary information
  !! @param[in]    coord         : coordinates of target systems
  !! @param[in]    npt           : flag for NPT or not
  !! @param[in]    reduce        : flag for reduce energy and virial
  !! @param[in]    nonb_ene      : flag for calculate nonbonded energy
  !! @param[in]    nonb_limiter  : flag for nonbond limiter
  !! @param[inout] energy        : energy information
  !! @param[input] coord_pbc     : !TODO
  !! @param[inout] force         : forces of target systems
  !! @param[inout] force_omp     : temprary forces of target systems
  !! @param[inout] force_pbc     : !TODO
  !! @param[inout] virial_cell   : !TODO
  !! @param[inout] virial        : virial term of target systems
  !! @param[inout] virial_ext    : extern virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_charmm(domain, enefunc, pairlist, boundary, coord, &
                                   npt, reduce, nonb_ene, merge_force,         &
                                   nonb_limiter,                               &
                                   energy, coord_pbc, force, force_long,       &
                                   force_omp, force_pbc, virial_cell, virial,  &
                                   virial_long, virial_ext)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(dp),                intent(in)    :: coord(:,:,:)
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: reduce
    logical,                 intent(in)    :: nonb_ene
    logical,                 intent(in)    :: merge_force
    logical,                 intent(in)    :: nonb_limiter
    type(s_energy),          intent(inout) :: energy
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(dp),                intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: force_long(:,:,:)
    real(wp),                intent(inout) :: force_omp(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3)
    real(dp),                intent(inout) :: virial_long(3,3)
    real(dp),                intent(inout) :: virial_ext(3,3)

    ! local variable
    real(dp)                 :: virial_omp(3,3,nthread)
    real(dp)                 :: virial_ext_omp(3,3,nthread)
    real(dp)                 :: elec_omp   (nthread)
    real(dp)                 :: evdw_omp   (nthread)
    real(dp)                 :: ebond_omp  (nthread)
    real(dp)                 :: eangle_omp (nthread)
    real(dp)                 :: eurey_omp  (nthread)
    real(dp)                 :: edihed_omp (nthread)
    real(dp)                 :: eimprop_omp(nthread)
    real(dp)                 :: ecmap_omp  (nthread)
    real(dp)                 :: eposi_omp  (nthread)
    real(dp)                 :: trans(1:3)
    real(dp)                 :: force_tmp(1:3)
    integer                  :: ncell, natom, id, i, j, k, ix, ic, jc
    integer                  :: omp_get_thread_num
    integer                  :: dimno_i, dimno_j
    integer                  :: atom_i, atom_j

    ! number of cells and atoms
    !
    ncell = domain%num_cell_local + domain%num_cell_boundary
    natom = domain%max_num_atom

    call timer(TimerTest1, TimerOn)
    ! initialization of energy and forces
    !
    call init_energy(energy)

    force      (1:3,1:natom,1:ncell) = 0.0_dp
    force_long (1:3,1:natom,1:ncell) = 0.0_dp
    virial     (1:3,1:3)             = 0.0_dp
    virial_long(1:3,1:3)             = 0.0_dp
    virial_ext (1:3,1:3)             = 0.0_dp

#ifndef USE_GPU
    !$omp parallel do
    do id = 1, nthread
      force_omp(1:3,1:natom,1:ncell,id) = 0.0_wp
      force_pbc(1:3,1:natom,1:ncell,id) = 0.0_wp
    end do
    !$omp end parallel do
#else
    !$omp parallel do
    do id = 1, nthread
      force_omp(1:3,1:natom,1:ncell,id) = 0.0_wp
    end do
    !$omp end parallel do
    !$omp parallel do
    do i = 1, ncell
      force_pbc(1:3,1:natom,i,1) = 0.0_wp
    end do
    !$omp end parallel do
#endif

    virial_omp    (1:3,1:3,1:nthread) = 0.0_dp
    virial_ext_omp(1:3,1:3,1:nthread) = 0.0_dp
    virial_cell   (1:3,1:maxcell) = 0.0_dp
    ebond_omp     (1:nthread) = 0.0_dp
    eangle_omp    (1:nthread) = 0.0_dp
    eurey_omp     (1:nthread) = 0.0_dp
    edihed_omp    (1:nthread) = 0.0_dp
    eimprop_omp   (1:nthread) = 0.0_dp
    ecmap_omp     (1:nthread) = 0.0_dp
    elec_omp      (1:nthread) = 0.0_dp
    evdw_omp      (1:nthread) = 0.0_dp
    eposi_omp     (1:nthread) = 0.0_dp
    call timer(TimerTest1, TimerOff)

    ! setup for emfit
    !
    if (enefunc%do_emfit) then
      call emfit_pre(domain, boundary, coord)
    end if

    ! nonbonded energy
    !
    select case (boundary%type)

    case (BoundaryTypePBC)

      if (enefunc%pme_use) then

        call compute_energy_nonbond_pme(domain, enefunc, pairlist, boundary, &
                              npt, nonb_ene, nonb_limiter, coord_pbc,        &
                              force_long,               &
                              force_omp, force_pbc,     &
                              virial_cell, virial_omp,  &
                              elec_omp, evdw_omp)
        if (.not.merge_force) then
          do id = 1, nthread
            do k = 1, 3
              virial_long(k,k) = virial_long(k,k) + virial_omp(k,k,id)
            end do
          end do
          virial_omp(1:3,1:3,1:nthread) = 0.0_dp
        end if

      else

        call compute_energy_nonbond_cutoff(domain, enefunc, pairlist, &
                              nonb_ene, force_pbc, virial_omp, &
                              elec_omp, evdw_omp)
      end if

    case default

      call error_msg('Compute_Energy_Charmm> Unknown boundary condition')

    end select
 
    if (real_calc) then
    
      ! bond energy
      !
      call compute_energy_bond(domain, enefunc, coord, &
                              force_omp, ebond_omp)

      ! angle energy
      !
      call compute_energy_angle(domain, enefunc, coord, &
                              force_omp, eangle_omp, eurey_omp)

      ! dihedral and cmap energies
      !
      if (enefunc%gamd_use) then

        if (enefunc%gamd%boost_dih .or. enefunc%gamd%boost_dual) then

          call compute_energy_dihed_gamd(domain, enefunc, npt, nonb_ene, &
            coord, force_omp, edihed_omp, ecmap_omp, eimprop_omp)

        else

          if (enefunc%local_restraint) then
            call compute_energy_dihed_localres(domain, enefunc, coord, &
                                  force_omp, edihed_omp)
          else
            call compute_energy_dihed(domain, enefunc, coord, &
                                  force_omp, edihed_omp)
          end if
      
          call compute_energy_cmap(domain, enefunc, coord, &
                                  force_omp, ecmap_omp)

        end if

      else

        if (enefunc%local_restraint) then
          call compute_energy_dihed_localres(domain, enefunc, coord, &
                                force_omp, edihed_omp)
        else
          call compute_energy_dihed(domain, enefunc, coord, &
                                force_omp, edihed_omp)
        end if
    
        call compute_energy_cmap(domain, enefunc, coord, &
                                force_omp, ecmap_omp)

      end if

      ! improper energy
      !
      call compute_energy_improp(domain, enefunc, coord, &
                              force_omp, eimprop_omp)
 
      ! 1-4 interaction 
      !
      call timer(TimerTest2, TimerOn)
      if (enefunc%pme_use) then
 
        call compute_energy_nonbond14_table_linear(domain, enefunc, force_omp, &
                              elec_omp, evdw_omp)
        call pme_bond_corr_linear(domain, enefunc, force_omp, elec_omp)
 
      end if
      call timer(TimerTest2, TimerOff)

      call timer(TimerTest3, TimerOn)

      ! virial for bonding
      !
      if (nonb_ene .or. npt) then

        !$omp parallel default(shared) private(id, i, ix, k) 
#ifdef OMP
        id = omp_get_thread_num()
#else
        id = 0
#endif
        do i = 1, ncell
          do ix = 1, domain%num_atom(i)
            do k = 1, 3
              virial_omp(k,k,id+1) = virial_omp(k,k,id+1) + &
                                     coord(k,ix,i)*force_omp(k,ix,i,id+1)
            end do
          end do
        end do
        !$omp end parallel

      end if
      call timer(TimerTest3, TimerOff)

      ! restraint energy
      !
      if (enefunc%restraint) then

        if (enefunc%gamd_use) then

          if (enefunc%gamd%boost_pot .or. enefunc%gamd%boost_dual) then

            call compute_energy_restraints_gamd(.true., .true., domain, &
                            boundary, enefunc, coord, &
                            force_omp, virial_omp, virial_ext_omp,           &
                            eposi_omp, energy)

          else

            call compute_energy_restraints(.true., .true., domain, boundary, &
                        enefunc, coord, &
                        force_omp, virial_omp, virial_ext_omp,           &
                        eposi_omp, energy%restraint_rmsd,                &
                        energy%rmsd, energy%restraint_distance,          &
                        energy%restraint_emfit, energy%emcorr)

          end if

        else

          call compute_energy_restraints(.true., .true., domain, boundary, &
                        enefunc, coord, &
                        force_omp, virial_omp, virial_ext_omp,           &
                        eposi_omp, energy%restraint_rmsd,                &
                        energy%rmsd, energy%restraint_distance,          &
                        energy%restraint_emfit, energy%emcorr)

        end if

      end if

    end if

    ! finish GPU
    !
#ifdef USE_GPU
    if (enefunc%pme_use) then
      call compute_energy_nonbond_pme_wait( &
           domain, enefunc, pairlist, boundary, &
           npt, nonb_ene, coord_pbc, force_omp, &
           force_pbc, virial_cell, virial_omp,  &
           elec_omp, evdw_omp)
    end if
#endif

    call timer(TimerTest4, TimerOn)

    ! virial with periodic boundary condition
    !
    if (enefunc%pme_use .and. (nonb_ene .or. npt)) then

      !$omp parallel default(shared) private(id, i, ix, ic, jc, trans, k)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
#ifndef USE_GPU
      do i = 1, ncell
        do ix = 1, domain%num_atom(i)
          do k = 1, 3
            virial_omp(k,k,id+1) = virial_omp(k,k,id+1) + &
                domain%translated(k,ix,i)*force_pbc(k,ix,i,id+1)
          end do
        end do
      end do
      do i = id+1, maxcell, nthread
        ic = domain%cell_pairlist1(1,i)
        jc = domain%cell_pairlist1(2,i)
        if (domain%virial_check(jc,ic) == 1) then
          trans(1:3) = domain%cell_move(1:3,jc,ic) * domain%system_size(1:3)
          do k = 1, 3
            virial_omp(k,k,id+1) = virial_omp(k,k,id+1)   &
                                   - trans(k)*virial_cell(k,i)
          end do
        end if
      end do
#else
      !$omp do schedule(dynamic,1)
      do i = 1, ncell
        do ix = 1, domain%num_atom(i)
          do k = 1, 3
            virial_omp(k,k,id+1) = virial_omp(k,k,id+1) + &
                                   domain%translated(k,ix,i)*force_pbc(k,ix,i,1)
          end do
        end do
      end do
#endif
      !$omp end parallel

    end if
    call timer(TimerTest4, TimerOff)

    call timer(TimerTest5, TimerOn)
    ! gather values
    !
    !$omp parallel do default(shared) private(id, i, ix, force_tmp) &
    !$omp schedule(dynamic,1)
    do i = 1, ncell
      do ix = 1, domain%num_atom(i)
        force_tmp(1:3) = force_omp(1:3,ix,i,1)
        force_tmp(1:3) = force_tmp(1:3) + force_pbc(1:3,ix,i,1)
        if (merge_force) force_tmp(1:3) = force_tmp(1:3) + force_long(1:3,ix,i)
#ifndef USE_GPU
        do id = 2, nthread
          force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,id) 
          force_tmp(1:3) = force_tmp(1:3) + force_pbc(1:3,ix,i,id) 
        end do
#else
        do id = 2, nthread
          force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,id)
        end do
#endif
        force(1:3,ix,i) = force_tmp(1:3)
      end do
    end do
    !$omp end parallel do

    do id = 1, nthread

      virial    (1:3,1:3) = virial    (1:3,1:3) + virial_omp    (1:3,1:3,id)
      virial_ext(1:3,1:3) = virial_ext(1:3,1:3) + virial_ext_omp(1:3,1:3,id)

      energy%bond               = energy%bond               + ebond_omp(id)
      energy%angle              = energy%angle              + eangle_omp(id)
      energy%urey_bradley       = energy%urey_bradley       + eurey_omp(id)
      energy%dihedral           = energy%dihedral           + edihed_omp(id)
      energy%improper           = energy%improper           + eimprop_omp(id)
      energy%cmap               = energy%cmap               + ecmap_omp(id)
      energy%electrostatic      = energy%electrostatic      + elec_omp(id)
      energy%van_der_waals      = energy%van_der_waals      + evdw_omp(id)
      energy%restraint_position = energy%restraint_position + eposi_omp(id)

    end do

    call timer(TimerTest5, TimerOff)
   
    ! total energy
    !
    energy%total = energy%bond          &
                 + energy%angle         &
                 + energy%urey_bradley  &
                 + energy%dihedral      &
                 + energy%cmap          &
                 + energy%improper      &
                 + energy%electrostatic &
                 + energy%van_der_waals &
                 + energy%restraint_position

    ! GaMD boost and statistics
    !
    if (enefunc%gamd_use) then
      call boost_gamd(domain, enefunc, energy, force, virial)
    end if

    if (reduce) &
      call reduce_ene(energy, virial)

    return

  end subroutine compute_energy_charmm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_amber
  !> @brief        compute potential energy with AMBER99 force field
  !! @authors      JJ
  !! @param[in]    domain        : domain information
  !! @param[in]    enefunc       : potential energy functions information
  !! @param[in]    pairlist      : pair list information
  !! @param[in]    boundary      : boundary information
  !! @param[in]    coord         : coordinates of target systems
  !! @param[in]    npt           : flag for NPT or not
  !! @param[in]    reduce        : flag for reduce energy and virial
  !! @param[in]    nonb_ene      : flag for calculate nonbonded energy
  !! @param[in]    nonb_limiter  : flag for nonbond limiter
  !! @param[inout] energy        : energy information
  !! @param[input] coord_pbc     : !TODO
  !! @param[inout] force         : forces of target systems
  !! @param[inout] force_omp     : temprary forces of target systems
  !! @param[inout] force_pbc     : !TODO
  !! @param[inout] virial_cell   : !TODO
  !! @param[inout] virial        : virial term of target systems
  !! @param[inout] virial_ext    : extern virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_amber(domain, enefunc, pairlist, boundary, coord,  &
                                  npt, reduce, nonb_ene, merge_force,          &
                                  nonb_limiter,                                &
                                  energy, coord_pbc, force, force_long,        &
                                  force_omp, force_pbc, virial_cell, virial,   &
                                  virial_long, virial_ext)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(dp),                intent(in)    :: coord(:,:,:)
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: reduce
    logical,                 intent(in)    :: nonb_ene
    logical,                 intent(in)    :: merge_force
    logical,                 intent(in)    :: nonb_limiter
    type(s_energy),          intent(inout) :: energy
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(dp),                intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: force_long(:,:,:)
    real(wp),                intent(inout) :: force_omp(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3)
    real(dp),                intent(inout) :: virial_long(3,3)
    real(dp),                intent(inout) :: virial_ext(3,3)

    ! local variable
    real(dp)                 :: virial_omp(3,3,nthread)
    real(dp)                 :: virial_ext_omp(3,3,nthread)
    real(dp)                 :: elec_omp   (nthread)
    real(dp)                 :: evdw_omp   (nthread)
    real(dp)                 :: ebond_omp  (nthread)
    real(dp)                 :: eangle_omp (nthread)
    real(dp)                 :: eurey_omp  (nthread)
    real(dp)                 :: edihed_omp (nthread)
    real(dp)                 :: eimprop_omp(nthread)
    real(dp)                 :: ecmap_omp  (nthread)
    real(dp)                 :: eposi_omp  (nthread)
    real(dp)                 :: trans(1:3)
    real(dp)                 :: force_tmp(1:3)
    integer                  :: ncell, natom, id, i, ix, ic, jc, k
    integer                  :: omp_get_thread_num

    ! number of cells and atoms
    !
    ncell = domain%num_cell_local + domain%num_cell_boundary
    natom = domain%max_num_atom

    ! initialization of energy and forces
    !
    call init_energy(energy)

    force      (1:3,1:natom,1:ncell) = 0.0_dp
    force_long (1:3,1:natom,1:ncell) = 0.0_dp
    virial     (1:3,1:3)             = 0.0_dp
    virial_long(1:3,1:3)             = 0.0_dp
    virial_ext (1:3,1:3)             = 0.0_dp

#ifndef USE_GPU
    force_omp     (1:3,1:natom,1:ncell,1:nthread) = 0.0_wp
    force_pbc     (1:3,1:natom,1:ncell,1:nthread) = 0.0_wp
#else
    force_omp     (1:3,1:natom,1:ncell,1:nthread) = 0.0_wp
    force_pbc     (1:3,1:natom,1:ncell,1) = 0.0_wp
#endif

    virial_omp    (1:3,1:3,1:nthread) = 0.0_dp
    virial_ext_omp(1:3,1:3,1:nthread) = 0.0_dp
    virial_cell   (1:3,1:maxcell) = 0.0_dp
    ebond_omp     (1:nthread) = 0.0_dp
    eangle_omp    (1:nthread) = 0.0_dp
    eurey_omp     (1:nthread) = 0.0_dp
    edihed_omp    (1:nthread) = 0.0_dp
    eimprop_omp   (1:nthread) = 0.0_dp
    ecmap_omp     (1:nthread) = 0.0_dp
    elec_omp      (1:nthread) = 0.0_dp
    evdw_omp      (1:nthread) = 0.0_dp
    eposi_omp     (1:nthread) = 0.0_dp

    ! setup for emfit
    !
    if (enefunc%do_emfit) then
      call emfit_pre(domain, boundary, coord)
    end if

    ! nonbonded energy
    !
    select case (boundary%type)

    case (BoundaryTypePBC)

      if (enefunc%pme_use) then

        call compute_energy_nonbond_pme(domain, enefunc, pairlist, boundary, &
                              npt, nonb_ene, nonb_limiter, coord_pbc, &
                              force_long,               &
                              force_omp, force_pbc,     &
                              virial_cell, virial_omp,  &
                              elec_omp, evdw_omp)
        if (.not.merge_force) then
          do id = 1, nthread
            do k = 1, 3
              virial_long(k,k) = virial_long(k,k) + virial_omp(k,k,id)
            end do
          end do
          virial_omp(1:3,1:3,1:nthread) = 0.0_dp
        end if

      else

        call compute_energy_nonbond_cutoff(domain, enefunc, pairlist, &
                              nonb_ene, force_pbc, virial_omp, &
                              elec_omp, evdw_omp)
      end if

    case default

      call error_msg('Compute_Energy_Amber> Unknown boundary condition')

    end select

    if (real_calc) then

      ! bond energy
      !
      call compute_energy_bond(domain, enefunc, coord, &
                              force_omp, ebond_omp)

      ! angle energy
      !
      call compute_energy_angle(domain, enefunc, coord, &
                              force_omp, eangle_omp, eurey_omp)

      ! dihedral energy
      !
      if (enefunc%gamd_use) then
        if (enefunc%gamd%boost_dih .or. enefunc%gamd%boost_dual) then

          call compute_energy_dihed_gamd(domain, enefunc, npt, nonb_ene, &
            coord, force_omp, edihed_omp, ecmap_omp, eimprop_omp)

        else

          if (enefunc%local_restraint) then
            call compute_energy_dihed_localres(domain, enefunc, coord, &
                                  force_omp, edihed_omp)
          else
            call compute_energy_dihed(domain, enefunc, coord, &
                                  force_omp, edihed_omp)
          end if

        end if

      else

        if (enefunc%local_restraint) then
          call compute_energy_dihed_localres(domain, enefunc, coord, &
                                force_omp, edihed_omp)
        else
          call compute_energy_dihed(domain, enefunc, coord, &
                                force_omp, edihed_omp)
        end if

      end if

      ! improper energy
      !
      call compute_energy_improp_cos(domain, enefunc, coord, &
                              force_omp, eimprop_omp)

      ! 1-4 interaction with linear table
      !
      if (enefunc%pme_use) then
 
        call compute_energy_nonbond14_table_linear(domain, enefunc, force_omp, &
                              elec_omp, evdw_omp)

        call pme_bond_corr_linear(domain, enefunc, force_omp, elec_omp)
 
      end if

      ! virial for bonding
      !
      if (nonb_ene .or. npt) then
        
        !$omp parallel default(shared) private(id, i, ix, k)
#ifdef OMP
        id = omp_get_thread_num()
#else
        id = 0
#endif
        do i = 1, ncell
          do ix = 1, domain%num_atom(i)
            do k = 1, 3
              virial_omp(k,k,id+1) = virial_omp(k,k,id+1) + &
                                     coord(k,ix,i)*force_omp(k,ix,i,id+1)
            end do
          end do
        end do
        !$omp end parallel

      end if

      ! restraint energy
      !
      if (enefunc%restraint) then

        if (enefunc%gamd_use) then

          if (enefunc%gamd%boost_pot .or. enefunc%gamd%boost_dual) then

            call compute_energy_restraints_gamd(.true., .true., domain, &
                         boundary, enefunc, coord, &
                         force_omp, virial_omp, virial_ext_omp,           &
                         eposi_omp, energy)

          else

            call compute_energy_restraints(.true., .true., domain, boundary, &
                          enefunc, coord, &
                          force_omp, virial_omp, virial_ext_omp,           &
                          eposi_omp, energy%restraint_rmsd,                &
                          energy%rmsd, energy%restraint_distance,          &
                          energy%restraint_emfit, energy%emcorr)

          end if

        else

          call compute_energy_restraints(.true., .true., domain, boundary, &
                          enefunc, coord, &
                          force_omp, virial_omp, virial_ext_omp,           &
                          eposi_omp, energy%restraint_rmsd,                &
                          energy%rmsd, energy%restraint_distance,          &
                          energy%restraint_emfit, energy%emcorr)

        end if
      end if
    end if

    ! finish GPU
    !
#ifdef USE_GPU
    if (enefunc%pme_use) then
      call compute_energy_nonbond_pme_wait( &
           domain, enefunc, pairlist, boundary, &
           npt, nonb_ene, coord_pbc, force_omp, &
           force_pbc, virial_cell, virial_omp,  &
           elec_omp, evdw_omp)
    end if
#endif

    ! virial with periodic boundary condition
    !
    if (enefunc%pme_use .and. (nonb_ene .or. npt)) then

      !$omp parallel default(shared) private(id, i, ix, ic, jc, trans, k)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
#ifndef USE_GPU
      do i = 1, ncell
        do ix = 1, domain%num_atom(i)
          do k = 1, 3
            virial_omp(k,k,id+1) = virial_omp(k,k,id+1) + &
                domain%translated(k,ix,i)*force_pbc(k,ix,i,id+1)
          end do
        end do
      end do
      do i = id+1, maxcell, nthread
        ic = domain%cell_pairlist1(1,i)
        jc = domain%cell_pairlist1(2,i)
        if (domain%virial_check(jc,ic) == 1) then
          trans(1:3) = domain%cell_move(1:3,jc,ic) * domain%system_size(1:3)
          do k = 1, 3
            virial_omp(k,k,id+1) = virial_omp(k,k,id+1) &
                                 - trans(k)*virial_cell(k,i)
          end do
        end if
      end do
#else
      !$omp do schedule(dynamic,1)
      do i = 1, ncell
        do ix = 1, domain%num_atom(i)
          do k = 1, 3
            virial_omp(k,k,id+1) = virial_omp(k,k,id+1) + &
                domain%translated(k,ix,i)*force_pbc(k,ix,i,1)
          end do
        end do
      end do
#endif
      !$omp end parallel

    end if

    ! gather values
    !
    !$omp parallel do default(shared) private(id, i, ix, force_tmp) &
    !$omp schedule(dynamic,1)
    do i = 1, ncell
      do ix = 1, domain%num_atom(i)
        force_tmp(1:3) = force_omp(1:3,ix,i,1)
        force_tmp(1:3) = force_tmp(1:3) + force_pbc(1:3,ix,i,1)
        if (merge_force) force_tmp(1:3) = force_tmp(1:3) + force_long(1:3,ix,i)
#ifndef USE_GPU
        do id = 2, nthread
          force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,id)
          force_tmp(1:3) = force_tmp(1:3) + force_pbc(1:3,ix,i,id)
        end do
#else
        do id = 2, nthread
          force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,id)
        end do
#endif
        force(1:3,ix,i) = force_tmp(1:3)
      end do
    end do
    !$omp end parallel do

    do id = 1, nthread

      virial    (1:3,1:3) = virial    (1:3,1:3) + virial_omp    (1:3,1:3,id)
      virial_ext(1:3,1:3) = virial_ext(1:3,1:3) + virial_ext_omp(1:3,1:3,id)

      energy%bond               = energy%bond               + ebond_omp(id)
      energy%angle              = energy%angle              + eangle_omp(id)
      energy%urey_bradley       = energy%urey_bradley       + eurey_omp(id)
      energy%dihedral           = energy%dihedral           + edihed_omp(id)
      energy%improper           = energy%improper           + eimprop_omp(id)
      energy%electrostatic      = energy%electrostatic      + elec_omp(id)
      energy%van_der_waals      = energy%van_der_waals      + evdw_omp(id)
      energy%restraint_position = energy%restraint_position + eposi_omp(id)

    end do

    ! total energy
    !
    energy%total = energy%bond          &
                 + energy%angle         &
                 + energy%urey_bradley  &
                 + energy%dihedral      &
                 + energy%improper      &
                 + energy%electrostatic &
                 + energy%van_der_waals &
                 + energy%restraint_position

    ! GaMD boost and statistics
    !
    if (enefunc%gamd_use) then 
      call boost_gamd(domain, enefunc, energy, force, virial)
    end if

    if (reduce) &
      call reduce_ene(energy, virial)

    return

  end subroutine compute_energy_amber

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_gro_amber
  !> @brief        compute potential energy with GROMACS-AMBER force field
  !! @authors      JJ
  !! @param[in]    domain      : domain information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[in]    pairlist    : pair list information
  !! @param[in]    boundary    : boundary information
  !! @param[in]    coord       : coordinates of target systems
  !! @param[in]    npt         : flag for NPT or not
  !! @param[in]    reduce      : flag for reduce energy and virial
  !! @param[in]    nonb_ene    : flag for calculate nonbonded energy
  !! @param[in]    nonb_limiter  : flag for nonbond limiter
  !! @param[inout] energy      : energy information
  !! @param[input] coord_pbc   : !TODO
  !! @param[inout] force       : forces of target systems
  !! @param[inout] force_omp   : temprary forces of target systems
  !! @param[inout] force_pbc   : !TODO
  !! @param[inout] virial_cell : !TODO
  !! @param[inout] virial      : virial term of target systems
  !! @param[inout] virial_ext  : extern virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_gro_amber(domain, enefunc, pairlist, boundary,     &
                                   coord, npt, reduce, nonb_ene, merge_force,  &
                                   nonb_limiter,                               &
                                   energy, coord_pbc, force, force_long,       &
                                   force_omp, force_pbc, virial_cell, virial,  &
                                   virial_long, virial_ext)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(dp),                intent(in)    :: coord(:,:,:)
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: reduce
    logical,                 intent(in)    :: nonb_ene
    logical,                 intent(in)    :: merge_force
    logical,                 intent(in)    :: nonb_limiter
    type(s_energy),          intent(inout) :: energy
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(dp),                intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: force_long(:,:,:)
    real(wp),                intent(inout) :: force_omp(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3)
    real(dp),                intent(inout) :: virial_long(3,3)
    real(dp),                intent(inout) :: virial_ext(3,3)

    ! local variable
    real(dp)                 :: virial_omp(3,3,nthread)
    real(dp)                 :: virial_ext_omp(3,3,nthread)
    real(dp)                 :: elec_omp    (nthread)
    real(dp)                 :: evdw_omp    (nthread)
    real(dp)                 :: ebond_omp   (nthread)
    real(dp)                 :: eangle_omp  (nthread)
    real(dp)                 :: eurey_omp   (nthread)
    real(dp)                 :: edihed_omp  (nthread)
    real(dp)                 :: erbdihed_omp(nthread)
    real(dp)                 :: eposi_omp   (nthread)
    real(dp)                 :: trans(1:3)
    real(dp)                 :: force_tmp(1:3)
    integer                  :: ncell, natom, id, i, ix, ic, jc, k
    integer                  :: omp_get_thread_num

    ! number of cells and atoms
    !
    ncell = domain%num_cell_local + domain%num_cell_boundary
    natom = domain%max_num_atom

    ! initialization of energy and forces
    !
    call init_energy(energy)

    force      (1:3,1:natom,1:ncell)   = 0.0_dp
    force_long (1:3,1:natom,1:ncell)   = 0.0_dp
    virial     (1:3,1:3)               = 0.0_dp
    virial_long(1:3,1:3)               = 0.0_dp
    virial_ext (1:3,1:3)               = 0.0_dp

#ifndef USE_GPU
    force_omp     (1:3,1:natom,1:ncell,1:nthread) = 0.0_wp
    force_pbc     (1:3,1:natom,1:ncell,1:nthread) = 0.0_wp
#else
    force_omp     (1:3,1:natom,1:ncell,1:nthread) = 0.0_wp
    force_pbc     (1:3,1:natom,1:ncell,1) = 0.0_wp
#endif
    virial_omp    (1:3,1:3,1:nthread) = 0.0_dp
    virial_ext_omp(1:3,1:3,1:nthread) = 0.0_dp
    virial_cell   (1:3,1:maxcell)     = 0.0_dp
    ebond_omp     (1:nthread)         = 0.0_dp
    eangle_omp    (1:nthread)         = 0.0_dp
    eurey_omp     (1:nthread)         = 0.0_dp
    edihed_omp    (1:nthread)         = 0.0_dp
    erbdihed_omp  (1:nthread)         = 0.0_dp
    elec_omp      (1:nthread)         = 0.0_dp
    evdw_omp      (1:nthread)         = 0.0_dp
    eposi_omp     (1:nthread)         = 0.0_dp

    ! setup for emfit
    !
    if (enefunc%do_emfit) then
      call emfit_pre(domain, boundary, coord)
    end if

    ! nonbonded energy
    !
    select case (boundary%type)

    case (BoundaryTypePBC)

      if (enefunc%pme_use) then

        call compute_energy_nonbond_pme(domain, enefunc, pairlist, boundary, &
                              npt, nonb_ene, nonb_limiter, coord_pbc,        &
                              force_long,               &
                              force_omp, force_pbc,     &
                              virial_cell, virial_omp,  &
                              elec_omp, evdw_omp)
        if (.not.merge_force) then
          do id = 1, nthread
            do k = 1, 3
              virial_long(k,k) = virial_long(k,k) + virial_omp(k,k,id)
            end do
          end do
          virial_omp(1:3,1:3,1:nthread) = 0.0_dp
        end if

      else

        call compute_energy_nonbond_cutoff(domain, enefunc, pairlist, &
                              nonb_ene, force_pbc, virial_omp, &
                              elec_omp, evdw_omp)
      end if

    case default

      call error_msg('Compute_Energy_Gro_Amber> Unknown boundary condition')

    end select

    if (real_calc) then
    
      ! bond energy
      !
      call compute_energy_bond(domain, enefunc, coord, &
                              force_omp, ebond_omp)

      ! angle energy
      !
      call compute_energy_angle(domain, enefunc, coord, &
                              force_omp, eangle_omp, eurey_omp)

      ! dihedral energy
      !
      call compute_energy_dihed(domain, enefunc, coord, &
                              force_omp, edihed_omp)
  
      ! Ryckaert-Bellemans dihedral energy
      !
      call compute_energy_rb_dihed(domain, enefunc, coord, &
                              force_omp, erbdihed_omp)

      if (enefunc%pme_use) then
 
        call compute_energy_nonbond14_table_linear(domain, enefunc, force_omp, &
                              elec_omp, evdw_omp)

        call pme_bond_corr_linear(domain, enefunc, force_omp, elec_omp)
 
      end if

      ! virial for bonding
      !
      if (nonb_ene .or. npt) then

        !$omp parallel default(shared) private(id, i, ix, k)
#ifdef OMP
        id = omp_get_thread_num()
#else
        id = 0
#endif
        do i = 1, ncell
          do ix = 1, domain%num_atom(i)
            do k = 1, 3
              virial_omp(k,k,id+1) = virial_omp(k,k,id+1) + &
                                     coord(k,ix,i)*force_omp(k,ix,i,id+1)
            end do
          end do
        end do
        !$omp end parallel

      end if

      ! restraint energy
      !
      if (enefunc%restraint) &
        call compute_energy_restraints(.true., .true., domain, boundary, enefunc, coord, &
                              force_omp, virial_omp, virial_ext_omp,           &
                              eposi_omp, energy%restraint_rmsd,                &
                              energy%rmsd, energy%restraint_distance,          &
                              energy%restraint_emfit, energy%emcorr)
    end if

    ! finish GPU
    !
#ifdef USE_GPU
    if (enefunc%pme_use) then
      call compute_energy_nonbond_pme_wait( &
           domain, enefunc, pairlist, boundary, &
           npt, nonb_ene, coord_pbc, force_omp, &
           force_pbc, virial_cell, virial_omp,  &
           elec_omp, evdw_omp)
    end if
#endif

    ! virial with periodic boundary condition
    !
    if (enefunc%pme_use .and. (nonb_ene .or. npt)) then

      !$omp parallel default(shared) private(id, i, ix, ic, jc, trans, k)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
#ifndef USE_GPU
      do i = 1, ncell
        do ix = 1, domain%num_atom(i)
          do k = 1, 3
            virial_omp(k,k,id+1) = virial_omp(k,k,id+1) + &
                domain%translated(k,ix,i)*force_pbc(k,ix,i,id+1)
          end do
        end do
      end do
      do i = id+1, maxcell, nthread
        ic = domain%cell_pairlist1(1,i)
        jc = domain%cell_pairlist1(2,i)
        if (domain%virial_check(jc,ic) == 1) then
          trans(1:3) = domain%cell_move(1:3,jc,ic) * domain%system_size(1:3)
          do k = 1, 3
            virial_omp(k,k,id+1) = virial_omp(k,k,id+1) - &
                                   trans(k)*virial_cell(k,i)
          end do
        end if
      end do
#else
      !$omp do schedule(dynamic,1)
      do i = 1, ncell
        do ix = 1, domain%num_atom(i)
          do k = 1, 3
            virial_omp(k,k,id+1) = virial_omp(k,k,id+1) + &
                                   domain%translated(k,ix,i)*force_pbc(k,ix,i,1)
          end do
        end do
      end do
#endif
      !$omp end parallel

    end if

    ! gather values
    !
    !$omp parallel do default(shared) private(id, i, ix, force_tmp) &
    !$omp schedule(dynamic,1)
    do i = 1, ncell
      do ix = 1, domain%num_atom(i)
        force_tmp(1:3) = force_omp(1:3,ix,i,1)
        force_tmp(1:3) = force_tmp(1:3) + force_pbc(1:3,ix,i,1)
        if (merge_force) force_tmp(1:3) = force_tmp(1:3) + force_long(1:3,ix,i)
#ifndef USE_GPU
        do id = 2, nthread
          force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,id)
          force_tmp(1:3) = force_tmp(1:3) + force_pbc(1:3,ix,i,id)
        end do
#else
        do id = 2, nthread
          force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,id)
        end do
#endif
        force(1:3,ix,i) = force_tmp(1:3)
      end do
    end do
    !$omp end parallel do

    do id = 1, nthread

      virial    (1:3,1:3) = virial    (1:3,1:3) + virial_omp    (1:3,1:3,id)
      virial_ext(1:3,1:3) = virial_ext(1:3,1:3) + virial_ext_omp(1:3,1:3,id)

      energy%bond               = energy%bond               + ebond_omp(id)
      energy%angle              = energy%angle              + eangle_omp(id)
      energy%urey_bradley       = energy%urey_bradley       + eurey_omp(id)
      energy%dihedral           = energy%dihedral           + edihed_omp(id)
      energy%dihedral           = energy%dihedral           + erbdihed_omp(id)
      energy%electrostatic      = energy%electrostatic      + elec_omp(id)
      energy%van_der_waals      = energy%van_der_waals      + evdw_omp(id)
      energy%restraint_position = energy%restraint_position + eposi_omp(id)

    end do

    ! total energy
    !
    energy%total = energy%bond          &
                 + energy%angle         &
                 + energy%urey_bradley  &
                 + energy%dihedral      &
                 + energy%improper      &
                 + energy%electrostatic &
                 + energy%van_der_waals &
                 + energy%restraint_position

    if (reduce) &
      call reduce_ene(energy, virial)

    return

  end subroutine compute_energy_gro_amber

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_gro_martini
  !> @brief        compute potential energy with GROMACS-MARTINI force field
  !! @authors      JJ
  !! @param[in]    domain      : domain information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[in]    pairlist    : pair list information
  !! @param[in]    boundary    : boundary information
  !! @param[in]    coord       : coordinates of target systems
  !! @param[in]    npt         : flag for NPT or not
  !! @param[in]    reduce      : flag for reduce energy and virial
  !! @param[in]    nonb_ene    : flag for calculate nonbonded energy
  !! @param[in]    nonb_limiter  : flag for nonbond limiter
  !! @param[inout] energy      : energy information
  !! @param[input] coord_pbc   : !TODO
  !! @param[inout] force       : forces of target systems
  !! @param[inout] force_omp   : temprary forces of target systems
  !! @param[inout] force_pbc   : !TODO
  !! @param[inout] virial_cell : !TODO
  !! @param[inout] virial      : virial term of target systems
  !! @param[inout] virial_ext  : extern virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_gro_martini(domain, enefunc, pairlist, boundary,   &
                              coord, npt, reduce, nonb_ene, nonb_limiter,      &
                              energy, coord_pbc, &
                              force, force_long, force_omp, force_pbc, &
                              virial_cell, virial, virial_ext)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(dp),                intent(in)    :: coord(:,:,:)
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: reduce
    logical,                 intent(in)    :: nonb_ene
    logical,                 intent(in)    :: nonb_limiter
    type(s_energy),          intent(inout) :: energy
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(dp),                intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: force_long(:,:,:)
    real(wp),                intent(inout) :: force_omp(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3)
    real(dp),                intent(inout) :: virial_ext(3,3)

    ! local variable
    real(dp)                 :: virial_omp(3,3,nthread)
    real(dp)                 :: virial_ext_omp(3,3,nthread)
    real(dp)                 :: elec_omp    (nthread)
    real(dp)                 :: evdw_omp    (nthread)
    real(dp)                 :: ebond_omp   (nthread)
    real(dp)                 :: eangle_omp  (nthread)
    real(dp)                 :: edihed_omp  (nthread)
    real(dp)                 :: eposi_omp   (nthread)
    integer                  :: ncell, natom, id, i, ix, k


    ! number of cells and atoms
    !
    ncell = domain%num_cell_local + domain%num_cell_boundary
    natom = domain%max_num_atom

    ! initialization of energy and forces
    !
    call init_energy(energy)

    force     (1:3,1:natom,1:ncell)   = 0.0_dp
    force_long(1:3,1:natom,1:ncell)   = 0.0_dp
    virial    (1:3,1:3)               = 0.0_dp
    virial_ext(1:3,1:3)               = 0.0_dp

    force_omp     (1:3,1:natom,1:ncell,1:nthread) = 0.0_wp
    virial_omp    (1:3,1:3,1:nthread) = 0.0_dp
    virial_ext_omp(1:3,1:3,1:nthread) = 0.0_dp
    ebond_omp     (1:nthread)         = 0.0_dp
    eangle_omp    (1:nthread)         = 0.0_dp
    edihed_omp    (1:nthread)         = 0.0_dp
    elec_omp      (1:nthread)         = 0.0_dp
    evdw_omp      (1:nthread)         = 0.0_dp
    eposi_omp     (1:nthread)         = 0.0_dp

    ! setup for emfit
    !
    if (enefunc%do_emfit) then
      call emfit_pre(domain, boundary, coord)
    end if

    if (real_calc) then

      ! bond energy
      !
      call compute_energy_bond(domain, enefunc, coord, &
                              force_omp, ebond_omp)

      ! angle energy
      !
      call compute_energy_angle_g96(domain, enefunc, coord, &
                              force_omp, eangle_omp)

      ! dihedral energy
      !
      call compute_energy_dihed(domain, enefunc, coord, &
                              force_omp, edihed_omp)

      ! virial for bonding
      !
      do id = 1, nthread
        do i = 1, ncell
          do ix = 1, domain%num_atom(i)
            do k = 1, 3
              virial(k,k) = virial(k,k) + coord(k,ix,i)*force_omp(k,ix,i,id)
            end do
          end do
        end do
      end do

      ! restraint energy
      !
      if (enefunc%restraint) &
        call compute_energy_restraints(.true., .true., domain, boundary, enefunc, coord, &
                              force_omp, virial_omp, virial_ext_omp,           &
                              eposi_omp, energy%restraint_rmsd,                &
                              energy%rmsd, energy%restraint_distance,          &
                              energy%restraint_emfit, energy%emcorr)
    end if


    ! nonbonded energy
    !
    select case (boundary%type)

    case (BoundaryTypePBC)

      call compute_energy_nonbond_cutoff(domain, enefunc, pairlist, &
                            nonb_ene, &
                            force_omp, virial_omp, elec_omp, evdw_omp)
    case default

      call error_msg('Compute_Energy_Gro_Martini> Unknown boundary condition') 

    end select


    ! gather values
    !
    do id = 1, nthread

      do i = 1, ncell
        do ix = 1, domain%num_atom(i)
          force(1,ix,i) = force(1,ix,i) + force_omp(1,ix,i,id)
          force(2,ix,i) = force(2,ix,i) + force_omp(2,ix,i,id)
          force(3,ix,i) = force(3,ix,i) + force_omp(3,ix,i,id)
        end do
      end do

      virial    (1:3,1:3) = virial    (1:3,1:3) + virial_omp    (1:3,1:3,id)
      virial_ext(1:3,1:3) = virial_ext(1:3,1:3) + virial_ext_omp(1:3,1:3,id)

      energy%bond               = energy%bond               + ebond_omp(id)
      energy%angle              = energy%angle              + eangle_omp(id)
      energy%dihedral           = energy%dihedral           + edihed_omp(id)
      energy%electrostatic      = energy%electrostatic      + elec_omp(id)
      energy%van_der_waals      = energy%van_der_waals      + evdw_omp(id)
      energy%restraint_position = energy%restraint_position + eposi_omp(id)

    end do

    ! total energy
    !
    energy%total = energy%bond          &
                 + energy%angle         &
                 + energy%dihedral      &
                 + energy%improper      &
                 + energy%electrostatic &
                 + energy%van_der_waals &
                 + energy%restraint_position

    if (reduce) &
      call reduce_ene(energy, virial)

    return

  end subroutine compute_energy_gro_martini

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_go
  !> @brief        compute potential energy with GROMACS-AMBER force field
  !! @authors      JJ
  !! @param[in]    domain        : domain information
  !! @param[in]    enefunc       : potential energy functions information
  !! @param[in]    pairlist      : pair list information
  !! @param[in]    boundary      : boundary information
  !! @param[in]    coord         : coordinates of target systems
  !! @param[in]    reduce        : flag for reduce energy and virial
  !! @param[in]    nonb_ene      : flag for calculate nonbonded energy
  !! @param[inout] energy        : energy information
  !! @param[inout] force         : forces of target systems
  !! @param[inout] virial        : virial term of target systems
  !! @param[inout] virial_ext    : extern virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_go(domain, enefunc, pairlist, boundary, coord,   &
                               reduce, nonb_ene, energy, force, force_long,  &
                               force_omp, virial_long, virial, virial_ext)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(dp),                intent(in)    :: coord(:,:,:)
    logical,                 intent(in)    :: reduce
    logical,                 intent(in)    :: nonb_ene
    type(s_energy),          intent(inout) :: energy
    real(dp),                intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: force_long(:,:,:)
    real(wp),                intent(inout) :: force_omp(:,:,:,:)
    real(dp),                intent(inout) :: virial_long(3,3)
    real(dp),                intent(inout) :: virial(3,3)
    real(dp),                intent(inout) :: virial_ext(3,3)

    ! local variable
    real(dp)                 :: virial_omp(3,3,nthread)
    real(dp)                 :: virial_ext_omp(3,3,nthread)
    real(dp)                 :: elec_omp    (nthread)
    real(dp)                 :: evdw_omp    (nthread)
    real(dp)                 :: ebond_omp   (nthread)
    real(dp)                 :: eangle_omp  (nthread)
    real(dp)                 :: eurey_omp   (nthread)
    real(dp)                 :: edihed_omp  (nthread)
    real(dp)                 :: eimprop_omp (nthread)
    real(dp)                 :: eposi_omp   (nthread)
    real(dp)                 :: econtact_omp(nthread)
    real(dp)                 :: enoncontact_omp (nthread)
    real(wp)                 :: force_tmp(1:3)
    integer                  :: ncell, natom, id, i, ix, ic, jc, k
    integer                  :: omp_get_thread_num

    ! number of cells and atoms
    !
    ncell = domain%num_cell_local + domain%num_cell_boundary
    natom = domain%max_num_atom

    ! initialization of energy and forces
    !
    call init_energy(energy)

    force      (1:3,1:natom,1:ncell)   = 0.0_wp
    virial     (1:3,1:3)               = 0.0_dp
    virial_ext (1:3,1:3)               = 0.0_dp

    force_omp     (1:3,1:natom,1:ncell,1:nthread) = 0.0_wp
    virial_omp    (1:3,1:3,1:nthread) = 0.0_dp
    virial_ext_omp(1:3,1:3,1:nthread) = 0.0_dp
    ebond_omp     (1:nthread)         = 0.0_dp
    eangle_omp    (1:nthread)         = 0.0_dp
    eurey_omp     (1:nthread)         = 0.0_dp
    edihed_omp    (1:nthread)         = 0.0_dp
    eimprop_omp   (1:nthread)         = 0.0_dp
    eposi_omp     (1:nthread)         = 0.0_dp
    econtact_omp  (1:nthread)         = 0.0_dp
    enoncontact_omp(1:nthread)        = 0.0_dp

    ! setup for emfit
    !
    if (enefunc%do_emfit) then
      call emfit_pre(domain, boundary, coord)
    end if


    ! bond energy
    !
    call compute_energy_bond(domain, enefunc, coord, &
                             force_omp, ebond_omp)

    ! angle energy
    !
    call compute_energy_angle(domain, enefunc, coord, &
                              force_omp, eangle_omp, eurey_omp)

    ! dihedral energy
    !
    call compute_energy_dihed(domain, enefunc, coord, &
                              force_omp, edihed_omp)

    ! improper energy
    !
    if (enefunc%forcefield == ForceFieldAAGO)          &
      call compute_energy_improp(domain, enefunc, coord, &
                                 force_omp, eimprop_omp)

    ! contact energy
    !
    call timer(TimerNonBond, TimerOn)
    if (enefunc%forcefield == ForceFieldAAGO) then
      call compute_energy_contact_126(domain, enefunc, coord,  &
                                      force_omp, econtact_omp, &
                                      enoncontact_omp)
    end if

    ! non-contact energy
    !
    call compute_energy_noncontact_nobc(domain, enefunc, pairlist, coord, &
                                        force_omp, enoncontact_omp)
    call timer(TimerNonBond, TimerOff)

    ! restraint energy
    !
    if (enefunc%restraint) &
      call compute_energy_restraints(.true., .true., domain, boundary, enefunc, coord, &
                            force_omp, virial_omp, virial_ext_omp,           &
                            eposi_omp, energy%restraint_rmsd,                &
                            energy%rmsd, energy%restraint_distance,          &
                            energy%restraint_emfit, energy%emcorr)

    ! gather values
    !
    !$omp parallel do default(shared) private(id, i, ix, force_tmp) &
    !$omp schedule(dynamic,1)
    do i = 1, ncell
      do ix = 1, domain%num_atom(i)
        force_tmp(1:3) = force_omp(1:3,ix,i,1)
        do id = 2, nthread
          force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,id)
        end do
        force(1:3,ix,i) = force_tmp(1:3)
      end do
    end do
    !$omp end parallel do

    do id = 1, nthread

      virial    (1:3,1:3) = virial    (1:3,1:3) + virial_omp    (1:3,1:3,id)
      virial_ext(1:3,1:3) = virial_ext(1:3,1:3) + virial_ext_omp(1:3,1:3,id)

      energy%bond               = energy%bond               + ebond_omp(id)
      energy%angle              = energy%angle              + eangle_omp(id)
      energy%dihedral           = energy%dihedral           + edihed_omp(id)
      energy%improper           = energy%improper           + eimprop_omp(id)
      energy%contact            = energy%contact            + econtact_omp(id)
      energy%noncontact         = energy%noncontact         + enoncontact_omp(id)
      energy%restraint_position = energy%restraint_position + eposi_omp(id)

    end do

    ! total energy
    !
    energy%total = energy%bond          &
                 + energy%angle         &
                 + energy%dihedral      &
                 + energy%improper      &
                 + energy%contact       &
                 + energy%noncontact    &
                 + energy%restraint_position

    if (reduce) &
      call reduce_ene_go(energy)

    return

  end subroutine compute_energy_go

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_charmm_short
  !> @brief        compute potential energy with charmm force field
  !! @authors      JJ
  !! @param[in]    domain      : domain information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[in]    pairlist    : pair list information
  !! @param[in]    boundary    : boundary information
  !! @param[in]    coord       : coordinates of target systems
  !! @param[in]    nonb_ene    : flag for calculate nonbonded energy
  !! @param[inout] energy      : energy information
  !! @param[input] coord_pbc   : !TODO
  !! @param[inout] force       : forces of target systems
  !! @param[inout] force_omp   : temprary forces of target systems
  !! @param[inout] force_pbc   : !TODO
  !! @param[inout] virial_cell : !TODO
  !! @param[inout] virial      : virial term of target systems
  !! @param[inout] virial_ext  : extern virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_charmm_short(domain, enefunc, pairlist, boundary, &
                                   coord, npt, nonb_ene, energy, coord_pbc, &
                                   force, force_omp, force_pbc, &
                                   virial_cell, virial, virial_ext)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(dp),                intent(in)    :: coord(:,:,:)
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: nonb_ene
    type(s_energy),          intent(inout) :: energy
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(dp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: force_omp(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3)
    real(dp),                intent(inout) :: virial_ext(3,3)

    ! local variable
    real(dp)                 :: virial_omp(3,3,nthread)
    real(dp)                 :: virial_ext_omp(3,3,nthread)
    real(dp)                 :: elec_omp   (nthread)
    real(dp)                 :: evdw_omp   (nthread)
    real(dp)                 :: ebond_omp  (nthread)
    real(dp)                 :: eangle_omp (nthread)
    real(dp)                 :: eurey_omp  (nthread)
    real(dp)                 :: edihed_omp (nthread)
    real(dp)                 :: eimprop_omp(nthread)
    real(dp)                 :: ecmap_omp  (nthread)
    real(dp)                 :: eposi_omp  (nthread)
    real(dp)                 :: trans(1:3)
    real(dp)                 :: force_tmp(1:3)
    integer                  :: ncell, natom, id, i, k, ix, ic, jc
    integer                  :: omp_get_thread_num

    ! number of cells and atoms
    !
    ncell = domain%num_cell_local + domain%num_cell_boundary
    natom = domain%max_num_atom

    ! initialization of energy and forces
    !
    call init_energy(energy)

    force     (1:3,1:natom,1:ncell) = 0.0_dp
    virial    (1:3,1:3)             = 0.0_dp
    virial_ext(1:3,1:3)             = 0.0_dp

#ifndef USE_GPU
    force_omp     (1:3,1:natom,1:ncell,1:nthread) = 0.0_wp
    force_pbc     (1:3,1:natom,1:ncell,1:nthread) = 0.0_wp
#else
    force_omp     (1:3,1:natom,1:ncell,1:nthread) = 0.0_wp
    force_pbc     (1:3,1:natom,1:ncell,1) = 0.0_wp
#endif
    virial_omp    (1:3,1:3,1:nthread) = 0.0_dp
    virial_ext_omp(1:3,1:3,1:nthread) = 0.0_dp
    virial_cell   (1:3,1:maxcell) = 0.0_dp
    ebond_omp     (1:nthread) = 0.0_dp
    eangle_omp    (1:nthread) = 0.0_dp
    eurey_omp     (1:nthread) = 0.0_dp
    edihed_omp    (1:nthread) = 0.0_dp
    eimprop_omp   (1:nthread) = 0.0_dp
    ecmap_omp     (1:nthread) = 0.0_dp
    elec_omp      (1:nthread) = 0.0_dp
    evdw_omp      (1:nthread) = 0.0_dp
    eposi_omp     (1:nthread) = 0.0_dp

    ! nonbonded energy
    !
    select case (boundary%type)

    case (BoundaryTypePBC)

      if (enefunc%pme_use) then

        call compute_energy_nonbond_pme_short(domain, enefunc, pairlist, &
                              npt, nonb_ene, coord_pbc,  &
                              force_omp, force_pbc, &
                              virial_cell, virial_omp, elec_omp, evdw_omp)

      else

        call error_msg( &
        'Compute_Energy_Charmm_Short>  Multi-step is available only with PME')

      end if

    case default

      call error_msg('Compute_Energy_Charmm_Short> Unknown boundary condition')

    end select
 
    ! bond energy
    !
    call compute_energy_bond(domain, enefunc, coord, &
                              force_omp, ebond_omp)

    ! angle energy
    !
    call compute_energy_angle(domain, enefunc, coord, &
                              force_omp, eangle_omp, eurey_omp)

    ! dihedral energy
    !
    if (enefunc%local_restraint) then
      call compute_energy_dihed_localres(domain, enefunc, coord, &
                              force_omp, edihed_omp)
    else
      call compute_energy_dihed(domain, enefunc, coord, &
                              force_omp, edihed_omp)
    end if

    ! improper energy
    !
    call compute_energy_improp(domain, enefunc, coord, &
                              force_omp, eimprop_omp)

    ! cmap energy
    ! 
    call compute_energy_cmap(domain, enefunc, coord, &
                              force_omp, ecmap_omp)

    ! 1-4 interaction with linear table
    !
    call compute_energy_nonbond14_table_linear(domain, enefunc, force_omp,   &
                             elec_omp, evdw_omp)
    call pme_bond_corr_linear(domain, enefunc, force_omp, elec_omp)

    ! restraint energy
    !
    if (enefunc%restraint) &
      call compute_energy_restraints(.true., .true., domain, boundary, enefunc, coord, &
                              force_omp, virial_omp, virial_ext_omp,         &
                              eposi_omp, energy%restraint_rmsd,              &
                              energy%rmsd, energy%restraint_distance,        &
                              energy%restraint_emfit, energy%emcorr)

    ! finish GPU
    !
#ifdef USE_GPU
    if (enefunc%pme_use) then
      call compute_energy_nonbond_pme_wait( &
           domain, enefunc, pairlist, boundary, &
           npt, nonb_ene, coord_pbc, force_omp, &
           force_pbc, virial_cell, virial_omp,  &
           elec_omp, evdw_omp)
    end if
#endif

    ! gather values
    !
    !$omp parallel do default(shared) private(id, i, ix, force_tmp) &
    !$omp schedule(dynamic,1)
    do i = 1, ncell
      do ix = 1, domain%num_atom(i)
        force_tmp(1:3) = force_omp(1:3,ix,i,1)
        force_tmp(1:3) = force_tmp(1:3) + force_pbc(1:3,ix,i,1)
#ifndef USE_GPU
        do id = 2, nthread
          force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,id)
          force_tmp(1:3) = force_tmp(1:3) + force_pbc(1:3,ix,i,id)
        end do
#else
        do id = 2, nthread
          force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,id)
        end do
#endif
        force(1:3,ix,i) = force_tmp(1:3)
      end do
    end do
    !$omp end parallel do

    return

  end subroutine compute_energy_charmm_short

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_amber_short
  !> @brief        compute potential energy with AMBER99 force field
  !! @authors      JJ
  !! @param[in]    domain      : domain information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[in]    pairlist    : pair list information
  !! @param[in]    boundary    : boundary information
  !! @param[in]    coord       : coordinates of target systems
  !! @param[in]    nonb_ene    : flag for calculate nonbonded energy
  !! @param[inout] energy      : energy information
  !! @param[input] coord_pbc   : !TODO
  !! @param[inout] force       : forces of target systems
  !! @param[inout] force_omp   : temprary forces of target systems
  !! @param[inout] force_pbc   : !TODO
  !! @param[inout] virial_cell : !TODO
  !! @param[inout] virial      : virial term of target systems
  !! @param[inout] virial_ext  : extern virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_amber_short(domain, enefunc, pairlist, boundary, &
                              coord, npt, nonb_ene, energy, coord_pbc, &
                              force, force_omp, force_pbc, &
                              virial_cell, virial, virial_ext)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(dp),                intent(in)    :: coord(:,:,:)
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: nonb_ene
    type(s_energy),          intent(inout) :: energy
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(dp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: force_omp(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3)
    real(dp),                intent(inout) :: virial_ext(3,3)

    ! local variable
    real(dp)                 :: virial_omp(3,3,nthread)
    real(dp)                 :: virial_ext_omp(3,3,nthread)
    real(dp)                 :: elec_omp   (nthread)
    real(dp)                 :: evdw_omp   (nthread)
    real(dp)                 :: ebond_omp  (nthread)
    real(dp)                 :: eangle_omp (nthread)
    real(dp)                 :: eurey_omp  (nthread)
    real(dp)                 :: edihed_omp (nthread)
    real(dp)                 :: eimprop_omp(nthread)
    real(dp)                 :: ecmap_omp  (nthread)
    real(dp)                 :: eposi_omp  (nthread)
    real(dp)                 :: trans(1:3)
    real(dp)                 :: force_tmp(1:3)
    integer                  :: ncell, natom, id, i, ix, ic, jc, k
    integer                  :: omp_get_thread_num

    ! number of cells and atoms
    !
    ncell = domain%num_cell_local + domain%num_cell_boundary
    natom = domain%max_num_atom

    ! initialization of energy and forces
    !
    call init_energy(energy)

    force     (1:3,1:natom,1:ncell) = 0.0_dp
    virial    (1:3,1:3)             = 0.0_dp
    virial_ext(1:3,1:3)             = 0.0_dp

#ifndef USE_GPU
    force_omp     (1:3,1:natom,1:ncell,1:nthread) = 0.0_wp
    force_pbc     (1:3,1:natom,1:ncell,1:nthread) = 0.0_wp
#else
    force_omp     (1:3,1:natom,1:ncell,1:nthread) = 0.0_wp
    force_pbc     (1:3,1:natom,1:ncell,1) = 0.0_wp
#endif
    virial_omp    (1:3,1:3,1:nthread) = 0.0_dp
    virial_ext_omp(1:3,1:3,1:nthread) = 0.0_dp
    virial_cell   (1:3,1:maxcell) = 0.0_dp
    ebond_omp     (1:nthread) = 0.0_dp
    eangle_omp    (1:nthread) = 0.0_dp
    eurey_omp     (1:nthread) = 0.0_dp
    edihed_omp    (1:nthread) = 0.0_dp
    eimprop_omp   (1:nthread) = 0.0_dp
    ecmap_omp     (1:nthread) = 0.0_dp
    elec_omp      (1:nthread) = 0.0_dp
    evdw_omp      (1:nthread) = 0.0_dp
    eposi_omp     (1:nthread) = 0.0_dp

    ! nonbonded energy
    !
    select case (boundary%type)

    case (BoundaryTypePBC)

      if (enefunc%pme_use) then

        call compute_energy_nonbond_pme_short(domain, enefunc, pairlist, &
                              npt, nonb_ene, coord_pbc,  &
                              force_omp, force_pbc, &
                              virial_cell, virial_omp, elec_omp, evdw_omp)

      else

        call error_msg( &
        'Compute_Energy_Amber_Short>  Multi-step is available only with PME')

      end if

    case default

      call error_msg('Compute_Energy_Amber_Short> Unknown boundary condition')

    end select

    ! bond energy
    !
    call compute_energy_bond(domain, enefunc, coord, &
                              force_omp, ebond_omp)

    ! angle energy
    !
    call compute_energy_angle(domain, enefunc, coord, &
                              force_omp, eangle_omp, eurey_omp)

    ! dihedral energy
    !
    if (enefunc%local_restraint) then
      call compute_energy_dihed_localres(domain, enefunc, coord, &
                              force_omp, edihed_omp)
    else
      call compute_energy_dihed(domain, enefunc, coord, &
                              force_omp, edihed_omp)
    end if

    ! improper energy
    !
    call compute_energy_improp_cos(domain, enefunc, coord, &
                              force_omp, eimprop_omp)

    ! 1-4 interaction with linear table
    !
    call compute_energy_nonbond14_table_linear(domain, enefunc, force_omp, &
                              elec_omp, evdw_omp)
    call pme_bond_corr_linear(domain, enefunc, force_omp, elec_omp)

    ! restraint energy
    !
    if (enefunc%restraint) &
      call compute_energy_restraints(.true., .true., domain, boundary, enefunc, coord, &
                              force_omp, virial_omp, virial_ext_omp,         &
                              eposi_omp, energy%restraint_rmsd,              &
                              energy%rmsd, energy%restraint_distance,        &
                              energy%restraint_emfit, energy%emcorr)

    ! finish GPU
    !
#ifdef USE_GPU
    if (enefunc%pme_use) then
      call compute_energy_nonbond_pme_wait( &
           domain, enefunc, pairlist, boundary, &
           npt, nonb_ene, coord_pbc, force_omp, &
           force_pbc, virial_cell, virial_omp,  &
           elec_omp, evdw_omp)
    end if
#endif

    ! gather values
    !
    !$omp parallel do default(shared) private(id, i, ix, force_tmp) &
    !$omp schedule(dynamic,1)
    do i = 1, ncell
      do ix = 1, domain%num_atom(i)
        force_tmp(1:3) = force_omp(1:3,ix,i,1)
        force_tmp(1:3) = force_tmp(1:3) + force_pbc(1:3,ix,i,1)
#ifndef USE_GPU
        do id = 2, nthread
          force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,id)
          force_tmp(1:3) = force_tmp(1:3) + force_pbc(1:3,ix,i,id)
        end do
#else
        do id = 2, nthread
          force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,id)
        end do
#endif
        force(1:3,ix,i) = force_tmp(1:3)
      end do
    end do
    !$omp end parallel do

    return

  end subroutine compute_energy_amber_short

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_gro_amber_short
  !> @brief        compute potential energy with GROMACS-AMBER force field
  !! @authors      JJ
  !! @param[in]    domain      : domain information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[in]    pairlist    : pair list information
  !! @param[in]    boundary    : boundary information
  !! @param[in]    coord       : coordinates of target systems
  !! @param[in]    nonb_ene    : flag for calculate nonbonded energy
  !! @param[inout] energy      : energy information
  !! @param[input] coord_pbc   : !TODO
  !! @param[inout] force       : forces of target systems
  !! @param[inout] force_omp   : temprary forces of target systems
  !! @param[inout] force_pbc   : !TODO
  !! @param[inout] virial_cell : !TODO
  !! @param[inout] virial      : virial term of target systems
  !! @param[inout] virial_ext  : extern virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_gro_amber_short(domain, enefunc, pairlist,         &
                            boundary, coord, npt, nonb_ene, energy, coord_pbc, &
                            force, force_omp, force_pbc,                       &
                            virial_cell, virial, virial_ext)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(dp),                intent(in)    :: coord(:,:,:)
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: nonb_ene
    type(s_energy),          intent(inout) :: energy
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(dp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: force_omp(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3)
    real(dp),                intent(inout) :: virial_ext(3,3)

    ! local variable
    real(dp)                 :: virial_omp(3,3,nthread)
    real(dp)                 :: virial_ext_omp(3,3,nthread)
    real(dp)                 :: elec_omp    (nthread)
    real(dp)                 :: evdw_omp    (nthread)
    real(dp)                 :: ebond_omp   (nthread)
    real(dp)                 :: eangle_omp  (nthread)
    real(dp)                 :: eurey_omp   (nthread)
    real(dp)                 :: edihed_omp  (nthread)
    real(dp)                 :: erbdihed_omp(nthread)
    real(dp)                 :: eposi_omp   (nthread)
    real(dp)                 :: trans(1:3)
    real(dp)                 :: force_tmp(1:3)
    integer                  :: ncell, natom, id, i, ix, ic, jc, k
    integer                  :: omp_get_thread_num

    ! number of cells and atoms
    !
    ncell = domain%num_cell_local + domain%num_cell_boundary
    natom = domain%max_num_atom

    ! initialization of energy and forces
    !
    call init_energy(energy)

    force     (1:3,1:natom,1:ncell) = 0.0_dp
    virial    (1:3,1:3)             = 0.0_dp
    virial_ext(1:3,1:3)             = 0.0_dp

#ifndef USE_GPU
    force_omp     (1:3,1:natom,1:ncell,1:nthread) = 0.0_wp
    force_pbc     (1:3,1:natom,1:ncell,1:nthread) = 0.0_wp
#else
    force_omp     (1:3,1:natom,1:ncell,1:nthread) = 0.0_wp
    force_pbc     (1:3,1:natom,1:ncell,1) = 0.0_wp
#endif
    virial_omp    (1:3,1:3,1:nthread) = 0.0_dp
    virial_ext_omp(1:3,1:3,1:nthread) = 0.0_dp
    virial_cell   (1:3,1:maxcell) = 0.0_dp
    ebond_omp     (1:nthread) = 0.0_dp
    eangle_omp    (1:nthread) = 0.0_dp
    eurey_omp     (1:nthread) = 0.0_dp
    edihed_omp    (1:nthread) = 0.0_dp
    erbdihed_omp  (1:nthread) = 0.0_dp
    elec_omp      (1:nthread) = 0.0_dp
    evdw_omp      (1:nthread) = 0.0_dp
    eposi_omp     (1:nthread) = 0.0_dp

    ! nonbonded energy
    !
    select case (boundary%type)

    case (BoundaryTypePBC)

      if (enefunc%pme_use) then

        call compute_energy_nonbond_pme_short(domain, enefunc, pairlist, &
                              npt, nonb_ene, coord_pbc,  &
                              force_omp, force_pbc, &
                              virial_cell, virial_omp, elec_omp, evdw_omp)

      else

        call error_msg( &
        'Compute_Energy_Gro_Amber_Short>  Multi-step is available only with PME')

      end if

    case default

      call error_msg('Compute_Energy_Gro_Amber_Short> Unknown boundary condition')

    end select

    ! bond energy
    !
    call compute_energy_bond(domain, enefunc, coord, &
                             force_omp, ebond_omp)

    ! angle energy
    !
    call compute_energy_angle(domain, enefunc, coord, &
                              force_omp, eangle_omp, eurey_omp)

    ! dihedral energy
    !
    call compute_energy_dihed(domain, enefunc, coord, &
                              force_omp, edihed_omp)
  
    ! Ryckaert-Bellemans dihedral energy
    !
    call compute_energy_rb_dihed(domain, enefunc, coord, &
                                 force_omp, erbdihed_omp)

    ! 1-4 interaction with linear table
    !
    call compute_energy_nonbond14_table_linear(domain, enefunc, force_omp, &
                                               elec_omp, evdw_omp)
    call pme_bond_corr_linear(domain, enefunc, force_omp, elec_omp)

    ! restraint energy
    !
    if (enefunc%restraint) &
      call compute_energy_restraints(.true., .true., domain, boundary, enefunc, coord, &
                              force_omp, virial_omp, virial_ext_omp,         &
                              eposi_omp, energy%restraint_rmsd,              &
                              energy%rmsd, energy%restraint_distance,        &
                              energy%restraint_emfit, energy%emcorr)

    ! finish GPU
    !
#ifdef USE_GPU
    if (enefunc%pme_use) then
      call compute_energy_nonbond_pme_wait( &
           domain, enefunc, pairlist, boundary, &
           npt, nonb_ene, coord_pbc, force_omp, &
           force_pbc, virial_cell, virial_omp,  &
           elec_omp, evdw_omp)
    end if
#endif

    ! gather values
    !
    !$omp parallel do default(shared) private(id, i, ix, force_tmp) &
    !$omp schedule(dynamic,1)
    do i = 1, ncell
      do ix = 1, domain%num_atom(i)
        force_tmp(1:3) = force_omp(1:3,ix,i,1)
        force_tmp(1:3) = force_tmp(1:3) + force_pbc(1:3,ix,i,1)
#ifndef USE_GPU
        do id = 2, nthread
          force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,id)
          force_tmp(1:3) = force_tmp(1:3) + force_pbc(1:3,ix,i,id)
        end do
#else
        do id = 2, nthread
          force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,id)
        end do
#endif
        force(1:3,ix,i) = force_tmp(1:3)
      end do
    end do
    !$omp end parallel do

    return

  end subroutine compute_energy_gro_amber_short

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_general_long
  !> @brief        compute long range interaction 
  !! @authors      JJ
  !! @param[in]    domain     : domain information
  !! @param[in]    enefunc    : potential energy functions information
  !! @param[in]    pairlist   : pair list information
  !! @param[in]    boundary   : boundary information
  !! @param[in]    coord      : coordinates of target systems
  !! @param[in]    npt        : flag for NPT or not
  !! @param[inout] energy     : energy information
  !! @param[inout] force      : forces of target systems
  !! @param[inout] force_omp  : temprary forces of target systems
  !! @param[inout] virial     : virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_general_long(domain, enefunc, pairlist, boundary, &
                                         coord,  npt, energy, &
                                         force, virial)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(dp),                intent(in)    :: coord(:,:,:)
    logical,                 intent(in)    :: npt
    type(s_energy),          intent(inout) :: energy
    real(dp),                intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: virial(3,3)

    ! local variable
    real(dp)                 :: virial_omp(3,3,nthread)
    real(dp)                 :: elec_long
    real(dp)                 :: elec_omp(nthread)
    integer                  :: ncell, natom, id, i, ix, k


    ! number of cells and atoms
    !
    ncell = domain%num_cell_local + domain%num_cell_boundary
    natom = domain%max_num_atom

    ! initialization of energy and forces
    !
    force(1:3,1:natom,1:ncell) = 0.0_dp
    virial(1:3,1:3)            = 0.0_dp
    elec_long                  = 0.0_dp

    virial_omp(1:3,1:3,1:nthread) = 0.0_dp
    elec_omp  (1:nthread) = 0.0_dp

    if (enefunc%pme_use) &
      call compute_energy_nonbond_pme_long(domain, enefunc, boundary, &
                              npt, force, virial_omp, elec_omp)

    do id = 1, nthread

      virial(1:3,1:3) = virial(1:3,1:3) + virial_omp(1:3,1:3,id)

      elec_long = elec_long + elec_omp(id)

    end do

    energy%electrostatic = energy%electrostatic + elec_long

    ! total energy
    energy%total = energy%total + elec_long

    return

  end subroutine compute_energy_general_long


#ifdef USE_GPU
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_pme_wait
  !> @brief        get the information from GPU to CPU
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_pme_wait(domain, enefunc, pairlist,        &
                                             boundary, npt, nonb_ene,          &
                                             coord_pbc, force, force_pbc,      &
                                             virial_cell, virial, eelec, evdw)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: nonb_ene
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3,nthread)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw(nthread)

    ! local
    real(dp) :: ene_virial(1:5)

    !
    ! GPU is used and overlapping is enabled
    !
    if (nonb_ene) then
      call gpu_wait_compute_energy_nonbond_table_linear_univ(   &
           coord_pbc, force_pbc(:,:,:,1), ene_virial)
      eelec(1) = eelec(1) + ene_virial(1)
      evdw (1) = evdw(1)  + ene_virial(2)
      virial(1,1,1) = virial(1,1,1) + ene_virial(3)
      virial(2,2,1) = virial(2,2,1) + ene_virial(4)
      virial(3,3,1) = virial(3,3,1) + ene_virial(5)
    else
      call gpu_wait_compute_force_nonbond_table_linear_univ( &
           coord_pbc, force_pbc(:,:,:,1), ene_virial, npt)
      virial(1,1,1) = virial(1,1,1) + ene_virial(1)
      virial(2,2,1) = virial(2,2,1) + ene_virial(2)
      virial(3,3,1) = virial(3,3,1) + ene_virial(3)
    end if
    call timer(TimerPmeReal, TimerOff)
    call timer(TimerNonBond, TimerOff)

  end subroutine compute_energy_nonbond_pme_wait
#endif

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_energy_genesis
  !> @brief        output energy in GENESIS style
  !! @authors      TM, CK
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
    real(dp)                 :: values(999)
    real(dp)                 :: ene_restraint


    write(title,'(A16)') 'STEP'
    write(frmt,'(A2,I2,A)') '(A',clength,')'
    write(frmt_cont,'(A2,I2,A3)') '(A',clength,',$)'
    write(frmt_res,'(A2,I2,A6)') '(A',clength-3,',I3.3)'
    write(rfrmt,'(A2,I2,A1,I1,A1)') '(F',clength,'.',flength,')'
    write(rfrmt_cont,'(A2,I2,A1,I1,A3)') '(F',clength,'.',flength,',$)'

    ifm = 1

    if (enefunc%num_bond_all > 0) then
      write(category(ifm),frmt) 'BOND'
      values(ifm) = energy%bond
      ifm = ifm+1
    endif

    if (enefunc%num_angl_all > 0) then
      write(category(ifm),frmt) 'ANGLE'
      values(ifm) = energy%angle
      ifm = ifm+1

      if (enefunc%forcefield == ForcefieldCHARMM) then
        write(category(ifm),frmt) 'UREY-BRADLEY'
        values(ifm) = energy%urey_bradley
        ifm = ifm+1
      endif
    endif

    if (enefunc%num_dihe_all > 0 .or. enefunc%num_rb_dihe_all > 0) then
      write(category(ifm),frmt) 'DIHEDRAL'
      values(ifm) = energy%dihedral
      ifm = ifm+1
    endif

    if (enefunc%num_impr_all > 0 ) then
      write(category(ifm),frmt) 'IMPROPER'
      values(ifm) = energy%improper
      ifm = ifm+1
    endif

    if (enefunc%forcefield == ForcefieldCHARMM) then

      if (enefunc%num_cmap_all > 0 ) then

        write(category(ifm),frmt) 'CMAP'
        values(ifm) = energy%cmap
        ifm = ifm+1
      endif
    endif

    if (enefunc%forcefield == ForcefieldAAGO) then

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
      endif

      write(category(ifm),frmt) 'ELECT'
      values(ifm) = energy%electrostatic
      ifm = ifm+1

    end if

    if (enefunc%num_restraintfuncs > 0) then
      if (enefunc%restraint_rmsd) then
        write(category(ifm),frmt) 'RMSD'
        values(ifm) = energy%rmsd
        ifm = ifm+1
      endif

      if (enefunc%restraint_emfit) then
        write(category(ifm),frmt) 'EMCORR'
        values(ifm) = energy%emcorr
        ifm = ifm+1
      end if

      ene_restraint =   energy%restraint_distance &
                      + energy%restraint_position &
                      + energy%restraint_rmsd     &
                      + energy%restraint_emfit
      write(category(ifm),frmt) 'RESTRAINT_TOTAL'
      values(ifm) = ene_restraint
      ifm = ifm+1

    endif

    if (enefunc%gamd_use) then
      if (enefunc%gamd%boost_pot) then
        write(category(ifm),frmt) 'POTENTIAL_GAMD'
        values(ifm) = energy%total_gamd
        ifm = ifm+1
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
      end if
    end if

    if (etitle) then

      write(MsgOut,'(A,$)') title

      do i = 1, ifm-1

        if (i == ifm-1) then
          write(MsgOut,frmt) category(i)
        else
          write(MsgOut,frmt_cont) category(i)
        endif
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
      endif
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
    real(dp)                 :: time, totener, totke, energy_, temperature
    real(dp)                 :: grms, hfctote, hfcke, ehfcor, virke
    real(dp)                 :: hbonds, asp, user
    real(dp)                 :: imnbvdw, imelec, imhbnd, rxnfield, extelec
    real(dp)                 :: ewksum, ewself, ewexcl, ewqcor, ewutil
    real(dp)                 :: vire, viri, presse, pressi, volume
    real(dp)                 :: cdihe, cintcr, noe

    real(dp),        pointer :: bonds, angles, urey_b, dihedrals, impropers
    real(dp),        pointer :: vdwaals, elec, cmaps, disp_corr
    real(dp),        pointer :: posicon, restdist


    time        = 0.0_dp
    hbonds      = 0.0_dp
    totener     = 0.0_dp
    totke       = 0.0_dp
    energy_     = 0.0_dp
    temperature = 0.0_dp
    asp         = 0.0_dp
    user        = 0.0_dp
    imnbvdw     = 0.0_dp
    imelec      = 0.0_dp
    imhbnd      = 0.0_dp
    rxnfield    = 0.0_dp
    extelec     = 0.0_dp
    ewksum      = 0.0_dp
    ewself      = 0.0_dp
    ewexcl      = 0.0_dp
    ewqcor      = 0.0_dp
    ewutil      = 0.0_dp
    vire        = 0.0_dp
    viri        = 0.0_dp
    presse      = 0.0_dp
    pressi      = 0.0_dp
    volume      = 0.0_dp
    grms        = 0.0_dp
    hfctote     = 0.0_dp
    hfcke       = 0.0_dp
    ehfcor      = 0.0_dp
    virke       = 0.0_dp
    volume      = 0.0_dp
    noe         = 0.0_dp
    cdihe       = 0.0_dp
    cintcr      = 0.0_dp

    ! write title if necessary
    !
    if (etitle) then
      write(MsgOut,'(A)') 'Output_Energy> CHARMM_Style is used'
      write(MsgOut,'(A)') ' '
      write(MsgOut,'(A79)') 'DYNA DYN: Step         Time      TOTEner        TOTKe       ENERgy  TEMPerature'
      write(MsgOut,'(A79)') 'DYNA PROP:             GRMS      HFCTote        HFCKe       EHFCor        VIRKe'
      write(MsgOut,'(A79)') 'DYNA INTERN:          BONDs       ANGLes       UREY-b    DIHEdrals    IMPRopers'

      if (enefunc%num_cmap_all > 0) then
        write(MsgOut,'(A79)') 'DYNA CROSS:           CMAPs                                                    '
      end if
      if (enefunc%dispersion_corr /= Disp_Corr_NONE) then
        write(MsgOut,'(A79)') 'DYNA  DISP:       Disp-Corr                                                    '
      endif

      write(MsgOut,'(A79)') 'DYNA EXTERN:        VDWaals         ELEC       HBONds          ASP         USER'

      write(MsgOut,'(A79)') 'DYNA IMAGES:        IMNBvdw       IMELec       IMHBnd       RXNField    EXTElec'
      write(MsgOut,'(A79)') 'DYNA EWALD:          EWKSum       EWSElf       EWEXcl       EWQCor       EWUTil'

      if (enefunc%restraint) then
        write(MsgOut,'(A79)') 'DYNA CONSTR:       HARMonic    CDIHedral          CIC     RESDistance       NOE'
      end if

      write(MsgOut,'(A79)') 'DYNA PRESS:            VIRE         VIRI       PRESSE       PRESSI       VOLUme'
      write(MsgOut,'(A79)') ' ----------       ---------    ---------    ---------    ---------    ---------'
      etitle = .false.

    end if


    bonds      => energy%bond
    angles     => energy%angle
    urey_b     => energy%urey_bradley
    dihedrals  => energy%dihedral
    impropers  => energy%improper
    elec       => energy%electrostatic
    vdwaals    => energy%van_der_waals
    cmaps      => energy%cmap
    restdist   => energy%restraint_distance
    posicon    => energy%restraint_position 
    disp_corr  => energy%disp_corr_energy

    ! write energy in CHARMM-style
    !
    write(MsgOut,'(A5,I9,5F13.5)') 'DYNA>', step, time, totener, totke, energy_, temperature
    write(MsgOut,'(A14,5F13.5)')   'DYNA PROP>    ', grms, hfctote, hfcke, ehfcor, virke
    write(MsgOut,'(A14,5F13.5)')   'DYNA INTERN>  ', bonds, angles, urey_b, dihedrals, impropers

    if (enefunc%num_cmap_all > 0) then
      write(MsgOut,'(A14, F13.5)')   'DYNA CROSS>   ', cmaps
    end if
    if (enefunc%dispersion_corr /= Disp_Corr_NONE) then
      write(MsgOut,'(A14,F13.5)')    'DYNA DISP>    ', disp_corr
    end if

    write(MsgOut,'(A14,5F13.5)')   'DYNA EXTERN>  ', vdwaals, elec, hbonds, asp, user

    write(MsgOut,'(A14,5F13.5)')   'DYNA IMAGES>  ', imnbvdw, imelec, imhbnd, rxnfield, extelec
    write(MsgOut,'(A14,5F13.5)')   'DYNA EWALD>   ', ewksum, ewself, ewexcl, ewqcor, ewutil

    if (enefunc%restraint) then
      write(MsgOut,'(A14,5F13.5)')   'DYNA CONSTR>  ', posicon, cdihe, cintcr, restdist, noe
    end if

    write(MsgOut,'(A14,5F13.5)')   'DYNA PRESS>   ', vire, viri, presse, pressi, volume
    write(MsgOut,'(A79)') ' ----------       ---------    ---------    ---------    ---------    ---------'
    write(MsgOut,'(A)') ' '

    return

  end subroutine output_energy_charmm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_energy_namd
  !> @brief        output energy in NAMD style
  !! @authors      YS, CK
  !! @param[in]    step    : step 
  !! @param[in]    energy : energy information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_energy_namd(step, energy)

    ! formal arguments
    integer,                 intent(in)    :: step
    type(s_energy),  target, intent(in)    :: energy

    ! local variables
    real(dp)                 :: ebond, eangle
    real(dp)                 :: misc, kinetic
    real(dp)                 :: total, temp, total2, total3, tempavg
    real(dp)                 :: pressure, gpressure, volume
    real(dp)                 :: pressavg, gpressavg

    real(dp),        pointer :: edihed, eimprp
    real(dp),        pointer :: eelect, evdw
    real(dp),        pointer :: eboundary


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


    ebond     =  energy%bond  + energy%restraint_distance
    eangle    =  energy%angle + energy%urey_bradley
    edihed    => energy%dihedral
    eimprp    => energy%improper
    eelect    => energy%electrostatic
    evdw      => energy%van_der_waals
    eboundary => energy%restraint_position

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
  !! @param[in]    energy  : energy information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_energy_gromacs(step, energy)

    ! formal arguments
    integer,                 intent(in)    :: step
    type(s_energy),  target, intent(in)    :: energy

    ! local variables
    real(dp)                 :: time, bonds, angles, urey_b
    real(dp)                 :: dihedrals, impropers, vdwaals, elec
    real(dp)                 :: restdist, posicon
    real(dp)                 :: energy_, totke, totener
    real(dp)                 :: temperature, pressi, presse
    real(dp)                 :: disp_corr
    integer                  :: i


    time        = 0.0_dp
    bonds       = energy%bond          * CAL2JOU
    angles      = energy%angle         * CAL2JOU
    urey_b      = energy%urey_bradley  * CAL2JOU
    dihedrals   = energy%dihedral      * CAL2JOU
    impropers   = energy%improper      * CAL2JOU
    vdwaals     = energy%van_der_waals * CAL2JOU
    elec        = energy%electrostatic * CAL2JOU
    disp_corr   = energy%disp_corr_energy * CAL2JOU
    restdist    = energy%restraint_distance
    posicon     = energy%restraint_position+energy%restraint_rmsd
    energy_     = 0.0_dp
    totke       = 0.0_dp
    totener     = 0.0_dp
    temperature = 0.0_dp
    pressi      = 0.0_dp
    presse      = 0.0_dp


    write(MsgOut,'(3A15)') &
         'Step', 'Time', 'Lambda'
    write(MsgOut,'(I15,2F15.5)') &
         step, time, 0.0_dp
    write(MsgOut,'(A)') &
         ' '
    write(MsgOut,'(A)') &
         '   Energies (kJ/mol)'

    write(MsgOut,'(5A15)') &
         'Bond', 'Angle', 'Urey-bradley', 'Dihedral', 'Improper Dih.'
    write(MsgOut,'(5ES15.5E2)') &
         bonds,angles,urey_b,dihedrals,impropers

    write(MsgOut,'(4A15)') &
       'LJ (1-4,SR', ' Coulomb(1-4,SR', 'Disper. corr.', 'Position Rest.'
    write(MsgOut,'(4ES15.5E2)') &
         vdwaals,elec,disp_corr,posicon
    write(MsgOut,'(2A15)') &
        'Potential', 'Kinetic En.'
    write(MsgOut,'(2ES15.5E2)') &
        energy_,totke


    write(MsgOut,'(5A15)') &
         'Total Energy', 'Temperature', 'Pressure(int.)', 'Pressure(ext.)'
    write(MsgOut,'(5ES15.5E2)') &
         totener,temperature,pressi,presse
    
    write(MsgOut,'(A)')  ' '

    return

  end subroutine output_energy_gromacs

  !======1=========2=========3=========4=========5=========6=========7=========8
  ! 
  !  Subroutine    reduce_ene
  !> @brief        reduce energy and virial
  !! @authors      JJ
  !! @param[inout] energy : energy information
  !! @param[inout] virial : virial term of 
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine reduce_ene(energy, virial)

    ! formal arguments
    type(s_energy),          intent(inout) :: energy
    real(dp),                intent(inout) :: virial(3,3)
    
#ifdef HAVE_MPI_GENESIS 

    ! local variables      
    real(dp)                 :: before_reduce(19), after_reduce(19)
    integer                  :: i, j, n


    ! Allreduce virial and energy components
    !
    n = 0
    do i = 1, 3
      do j = 1, 3
        n = n + 1
        before_reduce(n) = virial(i,j)
      end do
    end do

    before_reduce(10) = energy%bond
    before_reduce(11) = energy%angle
    before_reduce(12) = energy%urey_bradley
    before_reduce(13) = energy%dihedral
    before_reduce(14) = energy%improper
    before_reduce(15) = energy%cmap
    before_reduce(16) = energy%electrostatic
    before_reduce(17) = energy%van_der_waals
    before_reduce(18) = energy%restraint_position
    before_reduce(19) = energy%total

    call mpi_allreduce(before_reduce, after_reduce, 19, mpi_real8,  &
                       mpi_sum, mpi_comm_country, ierror)

    n = 0
    do i = 1, 3
      do j = 1, 3
        n = n + 1
        virial(i,j) = after_reduce(n)
      end do
    end do

    energy%bond               = after_reduce(10)
    energy%angle              = after_reduce(11)
    energy%urey_bradley       = after_reduce(12)
    energy%dihedral           = after_reduce(13)
    energy%improper           = after_reduce(14)
    energy%cmap               = after_reduce(15)
    energy%electrostatic      = after_reduce(16)
    energy%van_der_waals      = after_reduce(17)
    energy%restraint_position = after_reduce(18)
    energy%total              = after_reduce(19)

#endif

    return

  end subroutine reduce_ene

  subroutine reduce_ene_go(energy)

    ! formal arguments
    type(s_energy),          intent(inout) :: energy

#ifdef HAVE_MPI_GENESIS

    ! local variables
    real(dp)                 :: before_reduce(8), after_reduce(8)
    integer                  :: i, j, n


    ! Allreduce energy components
    !
    before_reduce(1) = energy%bond
    before_reduce(2) = energy%angle
    before_reduce(3) = energy%dihedral
    before_reduce(4) = energy%improper
    before_reduce(5) = energy%contact
    before_reduce(6) = energy%noncontact
    before_reduce(7) = energy%restraint_position
    before_reduce(8) = energy%total

    call mpi_reduce(before_reduce, after_reduce, 8, mpi_real8,  &
                    mpi_sum, 0, mpi_comm_country, ierror)

    energy%bond               = after_reduce(1)
    energy%angle              = after_reduce(2)
    energy%dihedral           = after_reduce(3)
    energy%improper           = after_reduce(4)
    energy%contact            = after_reduce(5)
    energy%noncontact         = after_reduce(6)
    energy%restraint_position = after_reduce(7)
    energy%total              = after_reduce(8)

#endif

    return

  end subroutine reduce_ene_go

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
    integer                  :: i, j, k, ifunc
    integer                  :: dimno_i, dimno_j
    integer                  :: atom_i, atom_j
    real(dp)                 :: etmp, stmp, dtmp
    real(dp)                 :: d(1:3)
    real(dp),    allocatable :: collection(:)

    if (enefunc%rpath_pos_func > 0) then

      allocate(collection(enefunc%stats_dimension))

      collection(:) = 0.0_dp

      call mpi_reduce(enefunc%stats_delta, collection, enefunc%stats_dimension,&
                      mpi_real8, mpi_sum, 0, mpi_comm_country, ierror)

    end if

    if (.not. replica_main_rank) then
      if (allocated(collection)) deallocate(collection)
      return
    endif

    do i = 1, enefunc%stats_dimension

      if (enefunc%rpath_pos_func > 0) then

        ifunc = enefunc%rpath_pos_func
        dtmp  = collection(i)
        enefunc%stats_force(i) = enefunc%stats_force(i) + &
        2.0_dp * real(enefunc%restraint_const(1,ifunc),dp) * dtmp

      else

        ifunc = enefunc%rpath_rest_function(i)
        dtmp  = enefunc%stats_delta(i)
        enefunc%stats_force(i) = enefunc%stats_force(i) + &
        2.0_dp * real(enefunc%restraint_const(1,ifunc),dp) * dtmp

      end if


    end do

    ifunc = enefunc%rpath_rest_function(1)

    if (enefunc%restraint_kind(ifunc) == RestraintsFuncPOSI .or. &
        enefunc%restraint_kind(ifunc) == RestraintsFuncPC .or. &
        enefunc%restraint_kind(ifunc) == RestraintsFuncPCCOM) then
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
    end if

    if (enefunc%rpath_pos_func > 0) then
      deallocate(collection)
    end if

    return

  end subroutine compute_stats

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_fep
  !> @brief        compute potential energy for FEP calculations
  !! @authors      NK
  !! @param[in]    domain        : domain information
  !! @param[in]    enefunc       : potential energy functions information
  !! @param[in]    pairlist      : pair list information
  !! @param[in]    boundary      : boundary information
  !! @param[in]    coord         : coordinates of target systems
  !! @param[in]    npt           : flag for NPT or not
  !! @param[in]    reduce        : flag for reduce energy and virial
  !! @param[in]    nonb_ene      : flag for calculate nonbonded energy
  !! @param[in]    nonb_limiter  : flag for nonbond limiter
  !! @param[inout] energy        : energy information
  !! @param[input] coord_pbc     : !TODO
  !! @param[inout] force         : forces of target systems
  !! @param[inout] force_omp     : temprary forces of target systems
  !! @param[inout] force_pbc     : !TODO
  !! @param[inout] virial_cell   : !TODO
  !! @param[inout] virial        : virial term of target systems
  !! @param[inout] virial_ext    : extern virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_fep(domain, enefunc, pairlist, boundary, coord,  &
                            npt, reduce, nonb_ene, merge_force,          &
                            nonb_limiter, energy,  &
                            coord_pbc, force, force_long, force_omp,     &
                            force_pbc, virial_cell, virial, virial_long, &
                            virial_ext)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(dp),                intent(in)    :: coord(:,:,:)
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: reduce
    logical,                 intent(in)    :: nonb_ene
    logical,                 intent(in)    :: merge_force
    logical,                 intent(in)    :: nonb_limiter
    type(s_energy),          intent(inout) :: energy
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(dp),                intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: force_long(:,:,:)
    real(wp),                intent(inout) :: force_omp(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3)
    real(dp),                intent(inout) :: virial_long(3,3)
    real(dp),                intent(inout) :: virial_ext(3,3)

    ! local variables
    real(dp)                 :: volume

    call timer(TimerEnergy, TimerOn)

    if (enefunc%gamd_use) &
      call error_msg('Compute_Energy> GaMD is not available in FEP')

    select case (enefunc%forcefield)

    case (ForcefieldCHARMM)

      call compute_energy_charmm_fep( &
                              domain, enefunc, pairlist, boundary, coord, &
                              npt, reduce, nonb_ene, merge_force,         &
                              nonb_limiter,                               &
                              energy, coord_pbc, force, force_long,       &
                              force_omp, force_pbc, virial_cell, virial,  &
                              virial_long, virial_ext)

    case (ForcefieldAMBER)
      call compute_energy_amber_fep( &
                              domain, enefunc, pairlist, boundary, coord, &
                              npt, reduce, nonb_ene, merge_force,         &
                              nonb_limiter,                               &
                              energy, coord_pbc, force, force_long,       &
                              force_omp, force_pbc, virial_cell, virial,  &
                              virial_long, virial_ext)

    case (ForcefieldGROAMBER)

      call compute_energy_gro_amber_fep( &
                              domain, enefunc, pairlist, boundary, coord, &
                              npt, reduce, nonb_ene, merge_force,         &
                              nonb_limiter,                               &
                              energy, coord_pbc, force, force_long,       &
                              force_omp, force_pbc, virial_cell, virial,  &
                              virial_long, virial_ext)

    case (ForcefieldGROMARTINI)

      call error_msg('Compute_Energy> gro_martini is not available in FEP')

    case (ForcefieldAAGO)

      call error_msg('Compute_Energy> AAGO is not available in FEP')

    end select

    ! Dispersion correction
    if (enefunc%dispersion_corr /= Disp_corr_NONE) then
      volume =  boundary%box_size_x_ref * &
                boundary%box_size_y_ref * &
                boundary%box_size_z_ref

      energy%disp_corr_energy = enefunc%dispersion_energy_preserve &
        + enefunc%lambljA*enefunc%dispersion_energy_vanish &
        + enefunc%lambljB*enefunc%dispersion_energy_appear
      energy%disp_corr_energy = energy%disp_corr_energy / volume

      if (enefunc%dispersion_corr == Disp_corr_EPress) then
        energy%disp_corr_virial = enefunc%dispersion_virial_preserve &
          + enefunc%lambljA*enefunc%dispersion_virial_vanish &
          + enefunc%lambljB*enefunc%dispersion_virial_appear
        energy%disp_corr_virial= energy%disp_corr_virial / volume
        if (replica_main_rank .or. main_rank) then
          virial(1,1) = virial(1,1) + real(energy%disp_corr_virial,dp)
          virial(2,2) = virial(2,2) + real(energy%disp_corr_virial,dp)
          virial(3,3) = virial(3,3) + real(energy%disp_corr_virial,dp)
        endif
      end if

    end if

    if (enefunc%rpath_sum_mf_flag) then
      call compute_stats(enefunc)
    end if

    call timer(TimerEnergy, TimerOff)

    return

  end subroutine compute_energy_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_energy_fep
  !> @brief        output energy for FEP calculations
  !! @authors      NK
  !! @param[in]    step    : step
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    energy  : energy information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_energy_fep(step, enefunc, energy)

    ! formal arguments
    integer,                 intent(in)    :: step
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_energy),          intent(in)    :: energy

    if (.not. main_rank) return

    select case (enefunc%output_style)

    case (OutputStyleGENESIS)

      call output_energy_genesis_fep(step, enefunc, energy)

    case (OutputStyleCHARMM)

      call output_energy_charmm(step, enefunc, energy)

    case (OutputStyleNAMD)

      call output_energy_namd(step, energy)

    case (OutputStyleGROMACS)

      call output_energy_gromacs(step, energy)

    end select

    return

  end subroutine output_energy_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_charmm_fep
  !> @brief        compute potential energy with charmm force field for FEP
  !! @authors      NK, HO
  !! @param[in]    domain        : domain information
  !! @param[in]    enefunc       : potential energy functions information
  !! @param[in]    pairlist      : pair list information
  !! @param[in]    boundary      : boundary information
  !! @param[in]    coord         : coordinates of target systems
  !! @param[in]    npt           : flag for NPT or not
  !! @param[in]    reduce        : flag for reduce energy and virial
  !! @param[in]    nonb_ene      : flag for calculate nonbonded energy
  !! @param[in]    nonb_limiter  : flag for nonbond limiter
  !! @param[inout] energy        : energy information
  !! @param[input] coord_pbc     : !TODO
  !! @param[inout] force         : forces of target systems
  !! @param[inout] force_omp     : temprary forces of target systems
  !! @param[inout] force_pbc     : !TODO
  !! @param[inout] virial_cell   : !TODO
  !! @param[inout] virial        : virial term of target systems
  !! @param[inout] virial_ext    : extern virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_charmm_fep(domain, enefunc, pairlist, boundary,    &
                                   coord, npt, reduce, nonb_ene, merge_force,  &
                                   nonb_limiter,                               &
                                   energy, coord_pbc, force, force_long,       &
                                   force_omp, force_pbc, virial_cell, virial,  &
                                   virial_long, virial_ext)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(dp),                intent(in)    :: coord(:,:,:)
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: reduce
    logical,                 intent(in)    :: nonb_ene
    logical,                 intent(in)    :: merge_force
    logical,                 intent(in)    :: nonb_limiter
    type(s_energy),          intent(inout) :: energy
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(dp),                intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: force_long(:,:,:)
    real(wp),                intent(inout) :: force_omp(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3)
    real(dp),                intent(inout) :: virial_long(3,3)
    real(dp),                intent(inout) :: virial_ext(3,3)

    ! local variable
    real(dp)                 :: virial_omp(3,3,nthread)
    real(dp)                 :: virial_ext_omp(3,3,nthread)
    real(dp)                 :: elec_omp    (nthread)
    real(dp)                 :: evdw_omp    (nthread)
    real(dp)                 :: ebond_omp   (nthread)
    real(dp)                 :: eangle_omp  (nthread)
    real(dp)                 :: eurey_omp   (nthread)
    real(dp)                 :: edihed_omp  (nthread)
    real(dp)                 :: eimprop_omp (nthread)
    real(dp)                 :: ecmap_omp   (nthread)
    real(dp)                 :: eposi_omp   (nthread)
    real(dp)                 :: trans(1:3)
    real(dp)                 :: force_tmp(1:3)
    integer                  :: ncell, natom, id, i, j, k, ix, ic, jc
    integer                  :: omp_get_thread_num
    integer                  :: dimno_i, dimno_j
    integer                  :: atom_i, atom_j

    ! number of cells and atoms
    !
    ncell = domain%num_cell_local + domain%num_cell_boundary
    natom = domain%max_num_atom

    call timer(TimerTest1, TimerOn)
    ! initialization of energy and forces
    !
    call init_energy(energy)

    force      (1:3,1:natom,1:ncell) = 0.0_dp
    force_long (1:3,1:natom,1:ncell) = 0.0_dp
    virial     (1:3,1:3)             = 0.0_dp
    virial_long(1:3,1:3)             = 0.0_dp
    virial_ext (1:3,1:3)             = 0.0_dp

#ifndef USE_GPU
    !$omp parallel do
    do id = 1, nthread
      force_omp(1:3,1:natom,1:ncell,id) = 0.0_wp
      force_pbc(1:3,1:natom,1:ncell,id) = 0.0_wp
    end do
    !$omp end parallel do
#else
    !$omp parallel do
    do id = 1, nthread
      force_omp(1:3,1:natom,1:ncell,id) = 0.0_wp
    end do
    !$omp end parallel do
    !$omp parallel do
    do i = 1, ncell
      force_pbc(1:3,1:natom,i,1) = 0.0_wp
    end do
    !$omp end parallel do
#endif
    virial_omp    (1:3,1:3,1:nthread) = 0.0_dp
    virial_ext_omp(1:3,1:3,1:nthread) = 0.0_dp
    virial_cell   (1:3,1:maxcell) = 0.0_dp
    ebond_omp     (1:nthread) = 0.0_dp
    eangle_omp    (1:nthread) = 0.0_dp
    eurey_omp     (1:nthread) = 0.0_dp
    edihed_omp    (1:nthread) = 0.0_dp
    eimprop_omp   (1:nthread) = 0.0_dp
    ecmap_omp     (1:nthread) = 0.0_dp
    elec_omp      (1:nthread) = 0.0_dp
    evdw_omp      (1:nthread) = 0.0_dp
    eposi_omp     (1:nthread) = 0.0_dp

    call timer(TimerTest1, TimerOff)

    ! setup for emfit
    !
    if (enefunc%do_emfit) then
      call emfit_pre(domain, boundary, coord)
    end if

    ! nonbonded energy
    !
    select case (boundary%type)

    case (BoundaryTypePBC)

      if (enefunc%pme_use) then

        call compute_energy_nonbond_pme_fep(domain, enefunc, pairlist, &
                              boundary, npt, nonb_ene, nonb_limiter, &
                              coord_pbc,                &
                              force_long,               &
                              force_omp, force_pbc,     &
                              virial_cell, virial_omp,  &
                              elec_omp, evdw_omp)

        if (.not.merge_force) then
          do id = 1, nthread
            do k = 1, 3
              virial_long(k,k) = virial_long(k,k) + virial_omp(k,k,id)
            end do
          end do
          virial_omp(1:3,1:3,1:nthread) = 0.0_dp
        end if

      else

        call compute_energy_nonbond_cutoff_fep(domain, enefunc, pairlist, &
                              nonb_ene, force_pbc, virial_omp, &
                              elec_omp, evdw_omp)

      end if

    case default

      call error_msg('Compute_Energy_Charmm> Unknown boundary condition')

    end select

    if (real_calc) then

      ! bond energy
      !
      call compute_energy_bond_fep(domain, enefunc, coord,  &
        force_omp, ebond_omp)

      ! angle energy
      !
      call compute_energy_angle_fep(domain, enefunc, coord,  &
        force_omp, eangle_omp, eurey_omp)

      ! dihedral energy
      !
      if (enefunc%local_restraint) then
        call compute_energy_dihed_localres_fep(domain, enefunc, coord, &
          force_omp, edihed_omp)
      else
        call compute_energy_dihed_fep(domain, enefunc, coord, &
          force_omp, edihed_omp)
      end if

      ! improper energy
      !
      call compute_energy_improp_fep(domain, enefunc, coord, &
        force_omp, eimprop_omp)

      ! cmap energy
      !
      call compute_energy_cmap_fep(domain, enefunc, coord, &
        force_omp, ecmap_omp)

      ! 1-4 interaction
      !
      call timer(TimerTest2, TimerOn)
      if (enefunc%pme_use) then

        call compute_energy_nonbond14_table_linear_fep(domain, enefunc, &
                              force_omp, elec_omp, evdw_omp)

        call pme_bond_corr_linear_fep(domain, enefunc, force_omp, elec_omp)

      end if
      call timer(TimerTest2, TimerOff)

      call timer(TimerTest3, TimerOn)

      ! virial for bonding
      !
      if (nonb_ene .or. npt) then

        !$omp parallel default(shared) private(id, i, ix, k)
#ifdef OMP
        id = omp_get_thread_num()
#else
        id = 0
#endif
        do i = 1, ncell
          do ix = 1, domain%num_atom(i)
            do k = 1, 3
              virial_omp(k,k,id+1) = virial_omp(k,k,id+1) + &
                                     coord(k,ix,i)*force_omp(k,ix,i,id+1)
            end do
          end do
        end do
        !$omp end parallel

      end if
      call timer(TimerTest3, TimerOff)

      ! restraint energy
      !
      if (enefunc%restraint) &
        call compute_energy_restraints_fep(.true., .true., domain, boundary, &
                              enefunc, coord, &
                              force_omp, virial_omp, virial_ext_omp,           &
                              eposi_omp, energy%restraint_rmsd,                &
                              energy%rmsd, energy%restraint_distance,          &
                              energy%restraint_emfit, energy%emcorr)

    end if

    ! finish GPU
    !
#ifdef USE_GPU
    if (enefunc%pme_use) then
      call compute_energy_nonbond_pme_wait( &
           domain, enefunc, pairlist, boundary, &
           npt, nonb_ene, coord_pbc, force_omp, &
           force_pbc, virial_cell, virial_omp,  &
           elec_omp, evdw_omp)
    end if
#endif

    call timer(TimerTest4, TimerOn)

    ! virial with periodic boundary condition
    !
    if (enefunc%pme_use .and. (nonb_ene .or. npt)) then

      !$omp parallel default(shared) private(id, i, ix, ic, jc, trans, k)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
#ifndef USE_GPU
      do i = 1, ncell
        do ix = 1, domain%num_atom(i)
          do k = 1, 3
            virial_omp(k,k,id+1) = virial_omp(k,k,id+1) + &
                domain%translated(k,ix,i)*force_pbc(k,ix,i,id+1)
          end do
        end do
      end do
      do i = id+1, maxcell, nthread
        ic = domain%cell_pairlist1(1,i)
        jc = domain%cell_pairlist1(2,i)
        if (domain%virial_check(jc,ic) == 1) then
          trans(1:3) = domain%cell_move(1:3,jc,ic) * domain%system_size(1:3)
          do k = 1, 3
            virial_omp(k,k,id+1) = virial_omp(k,k,id+1)   &
                                   - trans(k)*virial_cell(k,i)
          end do
        end if
      end do
#else
      !$omp do schedule(dynamic,1)
      do i = 1, ncell
        do ix = 1, domain%num_atom(i)
          do k = 1, 3
            virial_omp(k,k,id+1) = virial_omp(k,k,id+1) + &
                                   domain%translated(k,ix,i)*force_pbc(k,ix,i,1)
          end do
        end do
      end do
#endif
      !$omp end parallel

    end if
    call timer(TimerTest4, TimerOff)

    call timer(TimerTest5, TimerOn)
    ! gather values
    !
    !$omp parallel do default(shared) private(id, i, ix, force_tmp) &
    !$omp schedule(dynamic,1)
    do i = 1, ncell
      do ix = 1, domain%num_atom(i)
        force_tmp(1:3) = force_omp(1:3,ix,i,1)
        force_tmp(1:3) = force_tmp(1:3) + force_pbc(1:3,ix,i,1)
        if (merge_force) force_tmp(1:3) = force_tmp(1:3) + force_long(1:3,ix,i)
#ifndef USE_GPU
        do id = 2, nthread
          force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,id)
          force_tmp(1:3) = force_tmp(1:3) + force_pbc(1:3,ix,i,id)
        end do
#else
        do id = 2, nthread
          force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,id)
        end do
#endif
        force(1:3,ix,i) = force_tmp(1:3)
      end do
    end do
    !$omp end parallel do

    do id = 1, nthread

      virial    (1:3,1:3) = virial    (1:3,1:3) + virial_omp    (1:3,1:3,id)
      virial_ext(1:3,1:3) = virial_ext(1:3,1:3) + virial_ext_omp(1:3,1:3,id)

      energy%bond               = energy%bond               + ebond_omp(id)
      energy%angle              = energy%angle              + eangle_omp(id)
      energy%urey_bradley       = energy%urey_bradley       + eurey_omp(id)
      energy%dihedral           = energy%dihedral           + edihed_omp(id)
      energy%improper           = energy%improper           + eimprop_omp(id)
      energy%cmap               = energy%cmap               + ecmap_omp(id)
      energy%electrostatic      = energy%electrostatic      + elec_omp(id)
      energy%van_der_waals      = energy%van_der_waals      + evdw_omp(id)
      energy%restraint_position = energy%restraint_position + eposi_omp(id)
    end do

    call timer(TimerTest5, TimerOff)

    ! total energy
    !
    energy%total = energy%bond          &
                 + energy%angle         &
                 + energy%urey_bradley  &
                 + energy%dihedral      &
                 + energy%cmap          &
                 + energy%improper      &
                 + energy%electrostatic &
                 + energy%van_der_waals &
                 + energy%restraint_position

    if (reduce) &
      call reduce_ene(energy, virial)

    return

  end subroutine compute_energy_charmm_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_amber_fep
  !> @brief        compute potential energy with AMBER99 force field
  !! @authors      NK, HO
  !! @param[in]    domain        : domain information
  !! @param[in]    enefunc       : potential energy functions information
  !! @param[in]    pairlist      : pair list information
  !! @param[in]    boundary      : boundary information
  !! @param[in]    coord         : coordinates of target systems
  !! @param[in]    npt           : flag for NPT or not
  !! @param[in]    reduce        : flag for reduce energy and virial
  !! @param[in]    nonb_ene      : flag for calculate nonbonded energy
  !! @param[in]    nonb_limiter  : flag for nonbond limiter
  !! @param[inout] energy        : energy information
  !! @param[input] coord_pbc     : !TODO
  !! @param[inout] force         : forces of target systems
  !! @param[inout] force_omp     : temprary forces of target systems
  !! @param[inout] force_pbc     : !TODO
  !! @param[inout] virial_cell   : !TODO
  !! @param[inout] virial        : virial term of target systems
  !! @param[inout] virial_ext    : extern virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_amber_fep(domain, enefunc, pairlist, boundary, coord,  &
                                  npt, reduce, nonb_ene, merge_force,          &
                                  nonb_limiter,                                &
                                  energy, coord_pbc, force, force_long,        &
                                  force_omp, force_pbc, virial_cell, virial,   &
                                  virial_long, virial_ext)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(dp),                intent(in)    :: coord(:,:,:)
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: reduce
    logical,                 intent(in)    :: nonb_ene
    logical,                 intent(in)    :: merge_force
    logical,                 intent(in)    :: nonb_limiter
    type(s_energy),          intent(inout) :: energy
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(dp),                intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: force_long(:,:,:)
    real(wp),                intent(inout) :: force_omp(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3)
    real(dp),                intent(inout) :: virial_long(3,3)
    real(dp),                intent(inout) :: virial_ext(3,3)

    ! local variable
    real(dp)                 :: virial_omp(3,3,nthread)
    real(dp)                 :: virial_ext_omp(3,3,nthread)
    real(dp)                 :: elec_omp   (nthread)
    real(dp)                 :: evdw_omp   (nthread)
    real(dp)                 :: ebond_omp  (nthread)
    real(dp)                 :: eangle_omp (nthread)
    real(dp)                 :: eurey_omp  (nthread)
    real(dp)                 :: edihed_omp (nthread)
    real(dp)                 :: eimprop_omp(nthread)
    real(dp)                 :: ecmap_omp  (nthread)
    real(dp)                 :: eposi_omp  (nthread)
    real(dp)                 :: trans(1:3)
    real(dp)                 :: force_tmp(1:3)
    integer                  :: ncell, natom, id, i, ix, ic, jc, k
    integer                  :: omp_get_thread_num

    ! number of cells and atoms
    !
    ncell = domain%num_cell_local + domain%num_cell_boundary
    natom = domain%max_num_atom

    ! initialization of energy and forces
    !
    call init_energy(energy)

    force      (1:3,1:natom,1:ncell) = 0.0_dp
    force_long (1:3,1:natom,1:ncell) = 0.0_dp
    virial     (1:3,1:3)             = 0.0_dp
    virial_long(1:3,1:3)             = 0.0_dp
    virial_ext (1:3,1:3)             = 0.0_dp

#ifndef USE_GPU
    force_omp     (1:3,1:natom,1:ncell,1:nthread) = 0.0_wp
    force_pbc     (1:3,1:natom,1:ncell,1:nthread) = 0.0_wp
#else
    force_omp     (1:3,1:natom,1:ncell,1:nthread) = 0.0_wp
    force_pbc     (1:3,1:natom,1:ncell,1) = 0.0_wp
#endif
    virial_omp    (1:3,1:3,1:nthread) = 0.0_dp
    virial_ext_omp(1:3,1:3,1:nthread) = 0.0_dp
    virial_cell   (1:3,1:maxcell) = 0.0_dp
    ebond_omp     (1:nthread) = 0.0_dp
    eangle_omp    (1:nthread) = 0.0_dp
    eurey_omp     (1:nthread) = 0.0_dp
    edihed_omp    (1:nthread) = 0.0_dp
    eimprop_omp   (1:nthread) = 0.0_dp
    ecmap_omp     (1:nthread) = 0.0_dp
    elec_omp      (1:nthread) = 0.0_dp
    evdw_omp      (1:nthread) = 0.0_dp
    eposi_omp     (1:nthread) = 0.0_dp

    ! setup for emfit
    !
    if (enefunc%do_emfit) then
      call emfit_pre(domain, boundary, coord)
    end if

    ! nonbonded energy
    !
    select case (boundary%type)

    case (BoundaryTypePBC)

      if (enefunc%pme_use) then

        call compute_energy_nonbond_pme_fep(domain, enefunc, pairlist, boundary, &
                              npt, nonb_ene, nonb_limiter, coord_pbc, &
                              force_long,               &
                              force_omp, force_pbc,     &
                              virial_cell, virial_omp,  &
                              elec_omp, evdw_omp)

        if (.not.merge_force) then
          do id = 1, nthread
            do k = 1, 3
              virial_long(k,k) = virial_long(k,k) + virial_omp(k,k,id)
            end do
          end do
          virial_omp(1:3,1:3,1:nthread) = 0.0_dp
        end if

      else

        call compute_energy_nonbond_cutoff_fep(domain, enefunc, pairlist, &
                              nonb_ene, force_pbc, virial_omp, &
                              elec_omp, evdw_omp)

      end if

    case default

      call error_msg('Compute_Energy_Amber> Unknown boundary condition')

    end select

    if (real_calc) then

      ! bond energy
      !
      call compute_energy_bond_fep(domain, enefunc, coord,  &
        force_omp, ebond_omp)

      ! angle energy
      !
      call compute_energy_angle_fep(domain, enefunc, coord,  &
        force_omp, eangle_omp, eurey_omp)

      ! dihedral energy
      !
      if (enefunc%local_restraint) then
        call compute_energy_dihed_localres_fep(domain, enefunc, coord, &
          force_omp, edihed_omp)
      else
        call compute_energy_dihed_fep(domain, enefunc, coord, &
          force_omp, edihed_omp)
      end if

      ! improper energy
      !
      call compute_energy_improp_cos_fep(domain, enefunc, coord, &
        force_omp, eimprop_omp)

      ! 1-4 interaction with linear table
      !
      if (enefunc%pme_use) then

        call compute_energy_nonbond14_table_linear_fep(domain, enefunc, &
                              force_omp, elec_omp, evdw_omp)

        call pme_bond_corr_linear_fep(domain, enefunc, force_omp, elec_omp)

      end if

      ! virial for bonding
      !
      if (nonb_ene .or. npt) then

        !$omp parallel default(shared) private(id, i, ix, k)
#ifdef OMP
        id = omp_get_thread_num()
#else
        id = 0
#endif
        do i = 1, ncell
          do ix = 1, domain%num_atom(i)
            do k = 1, 3
              virial_omp(k,k,id+1) = virial_omp(k,k,id+1) + &
                                     coord(k,ix,i)*force_omp(k,ix,i,id+1)
            end do
          end do
        end do
        !$omp end parallel

      end if

      ! restraint energy
      !
      if (enefunc%restraint) &
        call compute_energy_restraints_fep(.true., .true., domain, boundary, &
                              enefunc, coord, &
                              force_omp, virial_omp, virial_ext_omp,           &
                              eposi_omp, energy%restraint_rmsd,                &
                              energy%rmsd, energy%restraint_distance,          &
                              energy%restraint_emfit, energy%emcorr)

    end if

    ! finish GPU
    !
#ifdef USE_GPU
    if (enefunc%pme_use) then
      call compute_energy_nonbond_pme_wait( &
           domain, enefunc, pairlist, boundary, &
           npt, nonb_ene, coord_pbc, force_omp, &
           force_pbc, virial_cell, virial_omp,  &
           elec_omp, evdw_omp)
    end if
#endif

    ! virial with periodic boundary condition
    !
    if (enefunc%pme_use .and. (nonb_ene .or. npt)) then

      !$omp parallel default(shared) private(id, i, ix, ic, jc, trans, k)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
#ifndef USE_GPU
      do i = 1, ncell
        do ix = 1, domain%num_atom(i)
          do k = 1, 3
            virial_omp(k,k,id+1) = virial_omp(k,k,id+1) + &
                domain%translated(k,ix,i)*force_pbc(k,ix,i,id+1)
          end do
        end do
      end do
      do i = id+1, maxcell, nthread
        ic = domain%cell_pairlist1(1,i)
        jc = domain%cell_pairlist1(2,i)
        if (domain%virial_check(jc,ic) == 1) then
          trans(1:3) = domain%cell_move(1:3,jc,ic) * domain%system_size(1:3)
          do k = 1, 3
            virial_omp(k,k,id+1) = virial_omp(k,k,id+1) &
                                 - trans(k)*virial_cell(k,i)
          end do
        end if
      end do
#else
      !$omp do schedule(dynamic,1)
      do i = 1, ncell
        do ix = 1, domain%num_atom(i)
          do k = 1, 3
            virial_omp(k,k,id+1) = virial_omp(k,k,id+1) + &
                domain%translated(k,ix,i)*force_pbc(k,ix,i,1)
          end do
        end do
      end do
#endif
      !$omp end parallel

    end if

    ! gather values
    !
    !$omp parallel do default(shared) private(id, i, ix, force_tmp) &
    !$omp schedule(dynamic,1)
    do i = 1, ncell
      do ix = 1, domain%num_atom(i)
        force_tmp(1:3) = force_omp(1:3,ix,i,1)
        force_tmp(1:3) = force_tmp(1:3) + force_pbc(1:3,ix,i,1)
        if (merge_force) force_tmp(1:3) = force_tmp(1:3) + force_long(1:3,ix,i)
#ifndef USE_GPU
        do id = 2, nthread
          force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,id)
          force_tmp(1:3) = force_tmp(1:3) + force_pbc(1:3,ix,i,id)
        end do
#else
        do id = 2, nthread
          force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,id)
        end do
#endif
        force(1:3,ix,i) = force_tmp(1:3)
      end do
    end do
    !$omp end parallel do

    do id = 1, nthread

      virial    (1:3,1:3) = virial    (1:3,1:3) + virial_omp    (1:3,1:3,id)
      virial_ext(1:3,1:3) = virial_ext(1:3,1:3) + virial_ext_omp(1:3,1:3,id)

      energy%bond               = energy%bond               + ebond_omp(id)
      energy%angle              = energy%angle              + eangle_omp(id)
      energy%urey_bradley       = energy%urey_bradley       + eurey_omp(id)
      energy%dihedral           = energy%dihedral           + edihed_omp(id)
      energy%improper           = energy%improper           + eimprop_omp(id)
      energy%electrostatic      = energy%electrostatic      + elec_omp(id)
      energy%van_der_waals      = energy%van_der_waals      + evdw_omp(id)
      energy%restraint_position = energy%restraint_position + eposi_omp(id)
    end do

    ! total energy
    !
    energy%total = energy%bond          &
                 + energy%angle         &
                 + energy%urey_bradley  &
                 + energy%dihedral      &
                 + energy%improper      &
                 + energy%electrostatic &
                 + energy%van_der_waals &
                 + energy%restraint_position

    if (reduce) &
      call reduce_ene(energy, virial)

    return

  end subroutine compute_energy_amber_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_gro_amber_fep
  !> @brief        compute potential energy with GROMACS-AMBER force field
  !                for FEP
  !! @authors      HO
  !! @param[in]    domain      : domain information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[in]    pairlist    : pair list information
  !! @param[in]    boundary    : boundary information
  !! @param[in]    coord       : coordinates of target systems
  !! @param[in]    npt         : flag for NPT or not
  !! @param[in]    reduce      : flag for reduce energy and virial
  !! @param[in]    nonb_ene    : flag for calculate nonbonded energy
  !! @param[in]    nonb_limiter  : flag for nonbond limiter
  !! @param[inout] energy      : energy information
  !! @param[input] coord_pbc   : !TODO
  !! @param[inout] force       : forces of target systems
  !! @param[inout] force_omp   : temprary forces of target systems
  !! @param[inout] force_pbc   : !TODO
  !! @param[inout] virial_cell : !TODO
  !! @param[inout] virial      : virial term of target systems
  !! @param[inout] virial_ext  : extern virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_gro_amber_fep(domain, enefunc, pairlist, boundary, &
                                   coord, npt, reduce, nonb_ene, merge_force,  &
                                   nonb_limiter,                               &
                                   energy, coord_pbc, force, force_long,       &
                                   force_omp, force_pbc, virial_cell, virial,  &
                                   virial_long, virial_ext)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(dp),                intent(in)    :: coord(:,:,:)
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: reduce
    logical,                 intent(in)    :: nonb_ene
    logical,                 intent(in)    :: merge_force
    logical,                 intent(in)    :: nonb_limiter
    type(s_energy),          intent(inout) :: energy
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(dp),                intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: force_long(:,:,:)
    real(wp),                intent(inout) :: force_omp(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3)
    real(dp),                intent(inout) :: virial_long(3,3)
    real(dp),                intent(inout) :: virial_ext(3,3)

    ! local variable
    real(dp)                 :: virial_omp(3,3,nthread)
    real(dp)                 :: virial_ext_omp(3,3,nthread)
    real(dp)                 :: elec_omp    (nthread)
    real(dp)                 :: evdw_omp    (nthread)
    real(dp)                 :: ebond_omp   (nthread)
    real(dp)                 :: eangle_omp  (nthread)
    real(dp)                 :: eurey_omp   (nthread)
    real(dp)                 :: edihed_omp  (nthread)
    real(dp)                 :: erbdihed_omp(nthread)
    real(dp)                 :: eposi_omp   (nthread)
    real(dp)                 :: trans(1:3)
    real(dp)                 :: force_tmp(1:3)
    integer                  :: ncell, natom, id, i, ix, ic, jc, k
    integer                  :: omp_get_thread_num

    ! number of cells and atoms
    !
    ncell = domain%num_cell_local + domain%num_cell_boundary
    natom = domain%max_num_atom

    ! initialization of energy and forces
    !
    call init_energy(energy)

    force      (1:3,1:natom,1:ncell)   = 0.0_dp
    force_long (1:3,1:natom,1:ncell)   = 0.0_dp
    virial     (1:3,1:3)               = 0.0_dp
    virial_long(1:3,1:3)               = 0.0_dp
    virial_ext (1:3,1:3)               = 0.0_dp

#ifndef USE_GPU
    force_omp     (1:3,1:natom,1:ncell,1:nthread) = 0.0_wp
    force_pbc     (1:3,1:natom,1:ncell,1:nthread) = 0.0_wp
#else
    force_omp     (1:3,1:natom,1:ncell,1:nthread) = 0.0_wp
    force_pbc     (1:3,1:natom,1:ncell,1) = 0.0_wp
#endif
    virial_omp    (1:3,1:3,1:nthread) = 0.0_dp
    virial_ext_omp(1:3,1:3,1:nthread) = 0.0_dp
    virial_cell   (1:3,1:maxcell)     = 0.0_dp
    ebond_omp     (1:nthread)         = 0.0_dp
    eangle_omp    (1:nthread)         = 0.0_dp
    eurey_omp     (1:nthread)         = 0.0_dp
    edihed_omp    (1:nthread)         = 0.0_dp
    erbdihed_omp  (1:nthread)         = 0.0_dp
    elec_omp      (1:nthread)         = 0.0_dp
    evdw_omp      (1:nthread)         = 0.0_dp
    eposi_omp     (1:nthread)         = 0.0_dp

    ! setup for emfit
    !
    if (enefunc%do_emfit) then
      call emfit_pre(domain, boundary, coord)
    end if

    ! nonbonded energy
    !
    select case (boundary%type)

    case (BoundaryTypePBC)

      if (enefunc%pme_use) then

        call compute_energy_nonbond_pme_fep(domain, enefunc, pairlist, boundary, &
                              npt, nonb_ene, nonb_limiter, coord_pbc,        &
                              force_long,               &
                              force_omp, force_pbc,     &
                              virial_cell, virial_omp,  &
                              elec_omp, evdw_omp)
        if (.not.merge_force) then
          do id = 1, nthread
            do k = 1, 3
              virial_long(k,k) = virial_long(k,k) + virial_omp(k,k,id)
            end do
          end do
          virial_omp(1:3,1:3,1:nthread) = 0.0_dp
        end if

      else

        call compute_energy_nonbond_cutoff_fep(domain, enefunc, pairlist, &
                              nonb_ene, force_pbc, virial_omp, &
                              elec_omp, evdw_omp)
      end if

    case default

      call error_msg('Compute_Energy_Gro_Amber> Unknown boundary condition')

    end select

    if (real_calc) then
    
      ! bond energy
      !
      call compute_energy_bond_fep(domain, enefunc, coord, &
                              force_omp, ebond_omp)

      ! angle energy
      !
      call compute_energy_angle_fep(domain, enefunc, coord, &
                              force_omp, eangle_omp, eurey_omp)

      ! dihedral energy
      !
      call compute_energy_dihed_fep(domain, enefunc, coord, &
                              force_omp, edihed_omp)
  
      ! Ryckaert-Bellemans dihedral energy
      !
      call compute_energy_rb_dihed_fep(domain, enefunc, coord, &
                              force_omp, erbdihed_omp)

      if (enefunc%pme_use) then
 
        call compute_energy_nonbond14_table_linear_fep(domain, enefunc, force_omp, &
                              elec_omp, evdw_omp)

        call pme_bond_corr_linear_fep(domain, enefunc, force_omp, elec_omp)
 
      end if

      ! virial for bonding
      !
      if (nonb_ene .or. npt) then

        !$omp parallel default(shared) private(id, i, ix, k)
#ifdef OMP
        id = omp_get_thread_num()
#else
        id = 0
#endif
        do i = 1, ncell
          do ix = 1, domain%num_atom(i)
            do k = 1, 3
              virial_omp(k,k,id+1) = virial_omp(k,k,id+1) + &
                                     coord(k,ix,i)*force_omp(k,ix,i,id+1)
            end do
          end do
        end do
        !$omp end parallel

      end if

      ! restraint energy
      !
      if (enefunc%restraint) &
        call compute_energy_restraints_fep(.true., .true., domain, boundary, enefunc, coord, &
                              force_omp, virial_omp, virial_ext_omp,           &
                              eposi_omp, energy%restraint_rmsd,                &
                              energy%rmsd, energy%restraint_distance,          &
                              energy%restraint_emfit, energy%emcorr)
    end if

    ! finish GPU
    !
#ifdef USE_GPU
    if (enefunc%pme_use) then
      call compute_energy_nonbond_pme_wait( &
           domain, enefunc, pairlist, boundary, &
           npt, nonb_ene, coord_pbc, force_omp, &
           force_pbc, virial_cell, virial_omp,  &
           elec_omp, evdw_omp)
    end if
#endif

    ! virial with periodic boundary condition
    !
    if (enefunc%pme_use .and. (nonb_ene .or. npt)) then

      !$omp parallel default(shared) private(id, i, ix, ic, jc, trans, k)
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
#ifndef USE_GPU
      do i = 1, ncell
        do ix = 1, domain%num_atom(i)
          do k = 1, 3
            virial_omp(k,k,id+1) = virial_omp(k,k,id+1) + &
                domain%translated(k,ix,i)*force_pbc(k,ix,i,id+1)
          end do
        end do
      end do
      do i = id+1, maxcell, nthread
        ic = domain%cell_pairlist1(1,i)
        jc = domain%cell_pairlist1(2,i)
        if (domain%virial_check(jc,ic) == 1) then
          trans(1:3) = domain%cell_move(1:3,jc,ic) * domain%system_size(1:3)
          do k = 1, 3
            virial_omp(k,k,id+1) = virial_omp(k,k,id+1) - &
                                   trans(k)*virial_cell(k,i)
          end do
        end if
      end do
#else
      !$omp do schedule(dynamic,1)
      do i = 1, ncell
        do ix = 1, domain%num_atom(i)
          do k = 1, 3
            virial_omp(k,k,id+1) = virial_omp(k,k,id+1) + &
                                   domain%translated(k,ix,i)*force_pbc(k,ix,i,1)
          end do
        end do
      end do
#endif
      !$omp end parallel

    end if

    ! gather values
    !
    !$omp parallel do default(shared) private(id, i, ix, force_tmp) &
    !$omp schedule(dynamic,1)
    do i = 1, ncell
      do ix = 1, domain%num_atom(i)
        force_tmp(1:3) = force_omp(1:3,ix,i,1)
        force_tmp(1:3) = force_tmp(1:3) + force_pbc(1:3,ix,i,1)
        if (merge_force) force_tmp(1:3) = force_tmp(1:3) + force_long(1:3,ix,i)
#ifndef USE_GPU
        do id = 2, nthread
          force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,id)
          force_tmp(1:3) = force_tmp(1:3) + force_pbc(1:3,ix,i,id)
        end do
#else
        do id = 2, nthread
          force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,id)
        end do
#endif
        force(1:3,ix,i) = force_tmp(1:3)
      end do
    end do
    !$omp end parallel do

    do id = 1, nthread

      virial    (1:3,1:3) = virial    (1:3,1:3) + virial_omp    (1:3,1:3,id)
      virial_ext(1:3,1:3) = virial_ext(1:3,1:3) + virial_ext_omp(1:3,1:3,id)

      energy%bond               = energy%bond               + ebond_omp(id)
      energy%angle              = energy%angle              + eangle_omp(id)
      energy%urey_bradley       = energy%urey_bradley       + eurey_omp(id)
      energy%dihedral           = energy%dihedral           + edihed_omp(id)
      energy%dihedral           = energy%dihedral           + erbdihed_omp(id)
      energy%electrostatic      = energy%electrostatic      + elec_omp(id)
      energy%van_der_waals      = energy%van_der_waals      + evdw_omp(id)
      energy%restraint_position = energy%restraint_position + eposi_omp(id)

    end do

    ! total energy
    !
    energy%total = energy%bond          &
                 + energy%angle         &
                 + energy%urey_bradley  &
                 + energy%dihedral      &
                 + energy%improper      &
                 + energy%electrostatic &
                 + energy%van_der_waals &
                 + energy%restraint_position

    if (reduce) &
      call reduce_ene(energy, virial)

    return

  end subroutine compute_energy_gro_amber_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_energy_genesis_fep
  !> @brief        output energy in GENESIS style for FEP calculation
  !! @authors      NK
  !! @param[in]    step    : step
  !! @param[in]    enefunc : information of potential functions
  !! @param[in]    energy  : energy information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_energy_genesis_fep(step, enefunc, energy)

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
    real(dp)                 :: values(999)
    real(dp)                 :: ene_restraint


    write(title,'(A16)') 'STEP'
    write(frmt,'(A2,I2,A)') '(A',clength,')'
    write(frmt_cont,'(A2,I2,A3)') '(A',clength,',$)'
    write(frmt_res,'(A2,I2,A6)') '(A',clength-3,',I3.3)'
    write(rfrmt,'(A2,I2,A1,I1,A1)') '(F',clength,'.',flength,')'
    write(rfrmt_cont,'(A2,I2,A1,I1,A3)') '(F',clength,'.',flength,',$)'

    ifm = 1

    if (enefunc%num_bond_all > 0) then
      write(category(ifm),frmt) 'BOND'
      values(ifm) = energy%bond
      ifm = ifm+1
    endif

    if (enefunc%num_angl_all > 0) then
      write(category(ifm),frmt) 'ANGLE'
      values(ifm) = energy%angle
      ifm = ifm+1

      if (enefunc%forcefield == ForcefieldCHARMM) then
        write(category(ifm),frmt) 'UREY-BRADLEY'
        values(ifm) = energy%urey_bradley
        ifm = ifm+1
      endif
    endif

    if (enefunc%num_dihe_all > 0 .or. enefunc%num_rb_dihe_all > 0) then
      write(category(ifm),frmt) 'DIHEDRAL'
      values(ifm) = energy%dihedral
      ifm = ifm+1
    endif

    if (enefunc%num_impr_all > 0 ) then
      write(category(ifm),frmt) 'IMPROPER'
      values(ifm) = energy%improper
      ifm = ifm+1
    endif

    if (enefunc%forcefield == ForcefieldCHARMM) then

      if (enefunc%num_cmap_all > 0 ) then

        write(category(ifm),frmt) 'CMAP'
        values(ifm) = energy%cmap
        ifm = ifm+1
      endif
    endif

    if (enefunc%forcefield == ForcefieldAAGO) then

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
      endif

      write(category(ifm),frmt) 'ELECT'
      values(ifm) = energy%electrostatic
      ifm = ifm+1

    end if

    if (enefunc%num_restraintfuncs > 0) then
      if (enefunc%restraint_rmsd) then
        write(category(ifm),frmt) 'RMSD'
        values(ifm) = energy%rmsd
        ifm = ifm+1
      endif

      if (enefunc%restraint_emfit) then
        write(category(ifm),frmt) 'EMCORR'
        values(ifm) = energy%emcorr
        ifm = ifm+1
      end if

      ene_restraint =   energy%restraint_distance &
                      + energy%restraint_position &
                      + energy%restraint_rmsd     &
                      + energy%restraint_emfit
      write(category(ifm),frmt) 'RESTRAINT_TOTAL'
      values(ifm) = ene_restraint
      ifm = ifm+1

    endif

    ! FEP
    write(MsgOut,'(A,F5.2)') "lambljA   = ", enefunc%lambljA
    write(MsgOut,'(A,F5.2)') "lambljB   = ", enefunc%lambljB
    write(MsgOut,'(A,F5.2)') "lambelA   = ", enefunc%lambelA
    write(MsgOut,'(A,F5.2)') "lambelB   = ", enefunc%lambelB
    write(MsgOut,'(A,F5.2)') "lambbondA = ", enefunc%lambbondA
    write(MsgOut,'(A,F5.2)') "lambbondB = ", enefunc%lambbondB

    if (etitle) then

      write(MsgOut,'(A,$)') title

      do i = 1, ifm-1

        if (i == ifm-1) then
          write(MsgOut,frmt) category(i)
        else
          write(MsgOut,frmt_cont) category(i)
        endif
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
      endif
    end do

    write(MsgOut,'(A)') ''

    return

  end subroutine output_energy_genesis_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_short_fep
  !> @brief        compute potential energy with short range interaction for FEP
  !! @authors      HO
  !! @param[in]    domain      : domain information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[in]    pairlist    : pair list information
  !! @param[in]    boundary    : boundary information
  !! @param[in]    coord       : coordinates of target systems
  !! @param[in]    npt         : flag for npt
  !! @param[in]    nonb_ene    : flag for calculate nonbonded energy
  !! @param[inout] energy      : energy information
  !! @param[input] coord_pbc   : !TODO
  !! @param[inout] force       : forces of target systems
  !! @param[inout] force_omp   : temprary forces of target systems
  !! @param[inout] force_pbc   : !TODO
  !! @param[inout] virial_cell : !TODO
  !! @param[inout] virial      : virial term of target systems
  !! @param[inout] virial_ext  : extern virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_short_fep(domain, enefunc, pairlist, boundary, coord, &
                                  npt, nonb_ene, energy, coord_pbc, &
                                  force, force_omp, force_pbc, &
                                  virial_cell, virial, virial_ext)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(dp),                intent(in)    :: coord(:,:,:)
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: nonb_ene
    type(s_energy),          intent(inout) :: energy
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(dp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: force_omp(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3)
    real(dp),                intent(inout) :: virial_ext(3,3)


    call timer(TimerEnergy, TimerOn)


    select case (enefunc%forcefield)

    case (ForcefieldCHARMM)

      call compute_energy_charmm_short_fep( &
                              domain, enefunc, pairlist, boundary, coord, &
                              npt, nonb_ene, energy, coord_pbc,           &
                              force, force_omp, force_pbc,                &
                              virial_cell, virial, virial_ext)

    case (ForcefieldAMBER)

      call compute_energy_amber_short_fep( &
                              domain, enefunc, pairlist, boundary, coord, &
                              npt, nonb_ene, energy, coord_pbc,           &
                              force, force_omp, force_pbc,                &
                              virial_cell, virial, virial_ext)

    case (ForcefieldGROAMBER)

      call compute_energy_gro_amber_short_fep( &
                              domain, enefunc, pairlist, boundary, coord, &
                              npt, nonb_ene, energy, coord_pbc,           &
                              force, force_omp, force_pbc,                &
                              virial_cell, virial, virial_ext)

    case (ForcefieldGROMARTINI)

      call error_msg('Compute_Energy_Short> gro_martini is not available in FEP')

    case (ForcefieldAAGO)

      call error_msg('Compute_Energy_Short> AAGO is not available in FEP')

    end select


    call timer(TimerEnergy, TimerOff)

    return

  end subroutine compute_energy_short_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_charmm_short_fep
  !> @brief        compute potential energy with charmm force field for FEP
  !! @authors      HO
  !! @param[in]    domain      : domain information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[in]    pairlist    : pair list information
  !! @param[in]    boundary    : boundary information
  !! @param[in]    coord       : coordinates of target systems
  !! @param[in]    npt         : flag for npt
  !! @param[in]    nonb_ene    : flag for calculate nonbonded energy
  !! @param[inout] energy      : energy information
  !! @param[input] coord_pbc   : !TODO
  !! @param[inout] force       : forces of target systems
  !! @param[inout] force_omp   : temprary forces of target systems
  !! @param[inout] force_pbc   : !TODO
  !! @param[inout] virial_cell : !TODO
  !! @param[inout] virial      : virial term of target systems
  !! @param[inout] virial_ext  : extern virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_charmm_short_fep(domain, enefunc, pairlist, boundary, &
                                   coord, npt, nonb_ene, energy, coord_pbc, &
                                   force, force_omp, force_pbc, &
                                   virial_cell, virial, virial_ext)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(dp),                intent(in)    :: coord(:,:,:)
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: nonb_ene
    type(s_energy),          intent(inout) :: energy
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(dp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: force_omp(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3)
    real(dp),                intent(inout) :: virial_ext(3,3)

    ! local variable
    real(dp)                 :: virial_omp(3,3,nthread)
    real(dp)                 :: virial_ext_omp(3,3,nthread)
    real(dp)                 :: elec_omp   (nthread)
    real(dp)                 :: evdw_omp   (nthread)
    real(dp)                 :: ebond_omp  (nthread)
    real(dp)                 :: eangle_omp (nthread)
    real(dp)                 :: eurey_omp  (nthread)
    real(dp)                 :: edihed_omp (nthread)
    real(dp)                 :: eimprop_omp(nthread)
    real(dp)                 :: ecmap_omp  (nthread)
    real(dp)                 :: eposi_omp  (nthread)
    real(dp)                 :: trans(1:3)
    real(dp)                 :: force_tmp(1:3)
    integer                  :: ncell, natom, id, i, k, ix, ic, jc
    integer                  :: omp_get_thread_num

    ! number of cells and atoms
    !
    ncell = domain%num_cell_local + domain%num_cell_boundary
    natom = domain%max_num_atom

    ! initialization of energy and forces
    !
    call init_energy(energy)

    force     (1:3,1:natom,1:ncell) = 0.0_dp
    virial    (1:3,1:3)             = 0.0_dp
    virial_ext(1:3,1:3)             = 0.0_dp

#ifndef USE_GPU
    force_omp     (1:3,1:natom,1:ncell,1:nthread) = 0.0_wp
    force_pbc     (1:3,1:natom,1:ncell,1:nthread) = 0.0_wp
#else
    force_omp     (1:3,1:natom,1:ncell,1:nthread) = 0.0_wp
    force_pbc     (1:3,1:natom,1:ncell,1) = 0.0_wp
#endif
    virial_omp    (1:3,1:3,1:nthread) = 0.0_dp
    virial_ext_omp(1:3,1:3,1:nthread) = 0.0_dp
    virial_cell   (1:3,1:maxcell) = 0.0_dp
    ebond_omp     (1:nthread) = 0.0_dp
    eangle_omp    (1:nthread) = 0.0_dp
    eurey_omp     (1:nthread) = 0.0_dp
    edihed_omp    (1:nthread) = 0.0_dp
    eimprop_omp   (1:nthread) = 0.0_dp
    ecmap_omp     (1:nthread) = 0.0_dp
    elec_omp      (1:nthread) = 0.0_dp
    evdw_omp      (1:nthread) = 0.0_dp
    eposi_omp     (1:nthread) = 0.0_dp

    ! nonbonded energy
    !
    select case (boundary%type)

    case (BoundaryTypePBC)

      if (enefunc%pme_use) then

        call compute_energy_nonbond_pme_short_fep(domain, enefunc, pairlist, &
                              npt, nonb_ene, coord_pbc,  &
                              force_omp, force_pbc, &
                              virial_cell, virial_omp, elec_omp, evdw_omp)

      else

        call error_msg( &
        'Compute_Energy_Charmm_Short>  Multi-step is available only with PME')

      end if

    case default

      call error_msg('Compute_Energy_Charmm_Short> Unknown boundary condition')

    end select
 
    ! bond energy
    !
    call compute_energy_bond_fep(domain, enefunc, coord,  &
                            force_omp, ebond_omp)

    ! angle energy
    !
    call compute_energy_angle_fep(domain, enefunc, coord,  &
      force_omp, eangle_omp, eurey_omp)

    ! dihedral energy
    !
    if (enefunc%local_restraint) then
      call compute_energy_dihed_localres_fep(domain, enefunc, coord, &
        force_omp, edihed_omp)
    else
      call compute_energy_dihed_fep(domain, enefunc, coord, &
        force_omp, edihed_omp)
    end if

    ! improper energy
    !
    call compute_energy_improp_fep(domain, enefunc, coord, &
                            force_omp, eimprop_omp)

    ! cmap energy
    !
    call compute_energy_cmap_fep(domain, enefunc, coord, &
                            force_omp, ecmap_omp)

    ! 1-4 interaction with linear table
    !
    call compute_energy_nonbond14_table_linear_fep(domain, enefunc, &
                          force_omp, elec_omp, evdw_omp)

    call pme_bond_corr_linear_fep(domain, enefunc, force_omp, elec_omp)

    ! restraint energy
    !
    if (enefunc%restraint) &
      call compute_energy_restraints_fep(.true., .true., domain, boundary, &
                              enefunc, coord, &
                              force_omp, virial_omp, virial_ext_omp,           &
                              eposi_omp, energy%restraint_rmsd,                &
                              energy%rmsd, energy%restraint_distance,          &
                              energy%restraint_emfit, energy%emcorr)

    ! finish GPU
    !
#ifdef USE_GPU
    if (enefunc%pme_use) then
      call compute_energy_nonbond_pme_wait( &
           domain, enefunc, pairlist, boundary, &
           npt, nonb_ene, coord_pbc, force_omp, &
           force_pbc, virial_cell, virial_omp,  &
           elec_omp, evdw_omp)
    end if
#endif

    ! gather values
    !
    !$omp parallel do default(shared) private(id, i, ix, force_tmp) &
    !$omp schedule(dynamic,1)
    do i = 1, ncell
      do ix = 1, domain%num_atom(i)
        force_tmp(1:3) = force_omp(1:3,ix,i,1)
        force_tmp(1:3) = force_tmp(1:3) + force_pbc(1:3,ix,i,1)
#ifndef USE_GPU
        do id = 2, nthread
          force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,id)
          force_tmp(1:3) = force_tmp(1:3) + force_pbc(1:3,ix,i,id)
        end do
#else
        do id = 2, nthread
          force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,id)
        end do
#endif
        force(1:3,ix,i) = force_tmp(1:3)
      end do
    end do
    !$omp end parallel do

    return

  end subroutine compute_energy_charmm_short_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_amber_short_fep
  !> @brief        compute potential energy with AMBER force field for FEP
  !! @authors      HO
  !! @param[in]    domain      : domain information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[in]    pairlist    : pair list information
  !! @param[in]    boundary    : boundary information
  !! @param[in]    coord       : coordinates of target systems
  !! @param[in]    npt         : flag for npt
  !! @param[in]    nonb_ene    : flag for calculate nonbonded energy
  !! @param[inout] energy      : energy information
  !! @param[input] coord_pbc   : !TODO
  !! @param[inout] force       : forces of target systems
  !! @param[inout] force_omp   : temprary forces of target systems
  !! @param[inout] force_pbc   : !TODO
  !! @param[inout] virial_cell : !TODO
  !! @param[inout] virial      : virial term of target systems
  !! @param[inout] virial_ext  : extern virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_amber_short_fep(domain, enefunc, pairlist, boundary, &
                              coord, npt, nonb_ene, energy, coord_pbc, &
                              force, force_omp, force_pbc, &
                              virial_cell, virial, virial_ext)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(dp),                intent(in)    :: coord(:,:,:)
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: nonb_ene
    type(s_energy),          intent(inout) :: energy
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(dp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: force_omp(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3)
    real(dp),                intent(inout) :: virial_ext(3,3)

    ! local variable
    real(dp)                 :: virial_omp(3,3,nthread)
    real(dp)                 :: virial_ext_omp(3,3,nthread)
    real(dp)                 :: elec_omp   (nthread)
    real(dp)                 :: evdw_omp   (nthread)
    real(dp)                 :: ebond_omp  (nthread)
    real(dp)                 :: eangle_omp (nthread)
    real(dp)                 :: eurey_omp  (nthread)
    real(dp)                 :: edihed_omp (nthread)
    real(dp)                 :: eimprop_omp(nthread)
    real(dp)                 :: ecmap_omp  (nthread)
    real(dp)                 :: eposi_omp  (nthread)
    real(dp)                 :: trans(1:3)
    real(dp)                 :: force_tmp(1:3)
    integer                  :: ncell, natom, id, i, ix, ic, jc, k
    integer                  :: omp_get_thread_num

    ! number of cells and atoms
    !
    ncell = domain%num_cell_local + domain%num_cell_boundary
    natom = domain%max_num_atom

    ! initialization of energy and forces
    !
    call init_energy(energy)

    force     (1:3,1:natom,1:ncell) = 0.0_dp
    virial    (1:3,1:3)             = 0.0_dp
    virial_ext(1:3,1:3)             = 0.0_dp

#ifndef USE_GPU
    force_omp     (1:3,1:natom,1:ncell,1:nthread) = 0.0_wp
    force_pbc     (1:3,1:natom,1:ncell,1:nthread) = 0.0_wp
#else
    force_omp     (1:3,1:natom,1:ncell,1:nthread) = 0.0_wp
    force_pbc     (1:3,1:natom,1:ncell,1) = 0.0_wp
#endif
    virial_omp    (1:3,1:3,1:nthread) = 0.0_dp
    virial_ext_omp(1:3,1:3,1:nthread) = 0.0_dp
    virial_cell   (1:3,1:maxcell) = 0.0_dp
    ebond_omp     (1:nthread) = 0.0_dp
    eangle_omp    (1:nthread) = 0.0_dp
    eurey_omp     (1:nthread) = 0.0_dp
    edihed_omp    (1:nthread) = 0.0_dp
    eimprop_omp   (1:nthread) = 0.0_dp
    ecmap_omp     (1:nthread) = 0.0_dp
    elec_omp      (1:nthread) = 0.0_dp
    evdw_omp      (1:nthread) = 0.0_dp
    eposi_omp     (1:nthread) = 0.0_dp

    ! nonbonded energy
    !
    select case (boundary%type)

    case (BoundaryTypePBC)

      if (enefunc%pme_use) then

        call compute_energy_nonbond_pme_short_fep(domain, enefunc, pairlist, &
                              npt, nonb_ene, coord_pbc,  &
                              force_omp, force_pbc, &
                              virial_cell, virial_omp, elec_omp, evdw_omp)

      else

        call error_msg( &
        'Compute_Energy_Amber_Short>  Multi-step is available only with PME')

      end if

    case default

      call error_msg('Compute_Energy_Amber_Short> Unknown boundary condition')

    end select

    ! bond energy
    !
    call compute_energy_bond_fep(domain, enefunc, coord,  &
      force_omp, ebond_omp)

    ! angle energy
    !
    call compute_energy_angle_fep(domain, enefunc, coord,  &
      force_omp, eangle_omp, eurey_omp)

    ! dihedral energy
    !
    if (enefunc%local_restraint) then
      call compute_energy_dihed_localres_fep(domain, enefunc, coord, &
        force_omp, edihed_omp)
    else
      call compute_energy_dihed_fep(domain, enefunc, coord, &
        force_omp, edihed_omp)
    end if

    ! improper energy
    !
    call compute_energy_improp_cos_fep(domain, enefunc, coord, &
      force_omp, eimprop_omp)

    ! 1-4 interaction with linear table
    !
    call compute_energy_nonbond14_table_linear_fep(domain, enefunc, &
                          force_omp, elec_omp, evdw_omp)

    call pme_bond_corr_linear_fep(domain, enefunc, force_omp, elec_omp)

    ! restraint energy
    !
    if (enefunc%restraint) &
      call compute_energy_restraints_fep(.true., .true., domain, boundary, &
                              enefunc, coord, &
                              force_omp, virial_omp, virial_ext_omp,           &
                              eposi_omp, energy%restraint_rmsd,                &
                              energy%rmsd, energy%restraint_distance,          &
                              energy%restraint_emfit, energy%emcorr)

    ! finish GPU
    !
#ifdef USE_GPU
    if (enefunc%pme_use) then
      call compute_energy_nonbond_pme_wait( &
           domain, enefunc, pairlist, boundary, &
           npt, nonb_ene, coord_pbc, force_omp, &
           force_pbc, virial_cell, virial_omp,  &
           elec_omp, evdw_omp)
    end if
#endif

    ! gather values
    !
    !$omp parallel do default(shared) private(id, i, ix, force_tmp) &
    !$omp schedule(dynamic,1)
    do i = 1, ncell
      do ix = 1, domain%num_atom(i)
        force_tmp(1:3) = force_omp(1:3,ix,i,1)
        force_tmp(1:3) = force_tmp(1:3) + force_pbc(1:3,ix,i,1)
#ifndef USE_GPU
        do id = 2, nthread
          force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,id)
          force_tmp(1:3) = force_tmp(1:3) + force_pbc(1:3,ix,i,id)
        end do
#else
        do id = 2, nthread
          force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,id)
        end do
#endif
        force(1:3,ix,i) = force_tmp(1:3)
      end do
    end do
    !$omp end parallel do

    return

  end subroutine compute_energy_amber_short_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_gro_amber_short_fep
  !> @brief        compute potential energy with GROMACS-AMBER force field
  !                for FEP
  !! @authors      HO
  !! @param[in]    domain      : domain information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[in]    pairlist    : pair list information
  !! @param[in]    boundary    : boundary information
  !! @param[in]    coord       : coordinates of target systems
  !! @param[in]    nonb_ene    : flag for calculate nonbonded energy
  !! @param[inout] energy      : energy information
  !! @param[input] coord_pbc   : !TODO
  !! @param[inout] force       : forces of target systems
  !! @param[inout] force_omp   : temprary forces of target systems
  !! @param[inout] force_pbc   : !TODO
  !! @param[inout] virial_cell : !TODO
  !! @param[inout] virial      : virial term of target systems
  !! @param[inout] virial_ext  : extern virial term of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_gro_amber_short_fep(domain, enefunc, pairlist,         &
                            boundary, coord, npt, nonb_ene, energy, coord_pbc, &
                            force, force_omp, force_pbc,                       &
                            virial_cell, virial, virial_ext)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    real(dp),                intent(in)    :: coord(:,:,:)
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: nonb_ene
    type(s_energy),          intent(inout) :: energy
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(dp),                intent(inout) :: force(:,:,:)
    real(wp),                intent(inout) :: force_omp(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3)
    real(dp),                intent(inout) :: virial_ext(3,3)

    ! local variable
    real(dp)                 :: virial_omp(3,3,nthread)
    real(dp)                 :: virial_ext_omp(3,3,nthread)
    real(dp)                 :: elec_omp    (nthread)
    real(dp)                 :: evdw_omp    (nthread)
    real(dp)                 :: ebond_omp   (nthread)
    real(dp)                 :: eangle_omp  (nthread)
    real(dp)                 :: eurey_omp   (nthread)
    real(dp)                 :: edihed_omp  (nthread)
    real(dp)                 :: erbdihed_omp(nthread)
    real(dp)                 :: eposi_omp   (nthread)
    real(dp)                 :: trans(1:3)
    real(dp)                 :: force_tmp(1:3)
    integer                  :: ncell, natom, id, i, ix, ic, jc, k
    integer                  :: omp_get_thread_num

    ! number of cells and atoms
    !
    ncell = domain%num_cell_local + domain%num_cell_boundary
    natom = domain%max_num_atom

    ! initialization of energy and forces
    !
    call init_energy(energy)

    force     (1:3,1:natom,1:ncell) = 0.0_dp
    virial    (1:3,1:3)             = 0.0_dp
    virial_ext(1:3,1:3)             = 0.0_dp

#ifndef USE_GPU
    force_omp     (1:3,1:natom,1:ncell,1:nthread) = 0.0_wp
    force_pbc     (1:3,1:natom,1:ncell,1:nthread) = 0.0_wp
#else
    force_omp     (1:3,1:natom,1:ncell,1:nthread) = 0.0_wp
    force_pbc     (1:3,1:natom,1:ncell,1) = 0.0_wp
#endif
    virial_omp    (1:3,1:3,1:nthread) = 0.0_dp
    virial_ext_omp(1:3,1:3,1:nthread) = 0.0_dp
    virial_cell   (1:3,1:maxcell) = 0.0_dp
    ebond_omp     (1:nthread) = 0.0_dp
    eangle_omp    (1:nthread) = 0.0_dp
    eurey_omp     (1:nthread) = 0.0_dp
    edihed_omp    (1:nthread) = 0.0_dp
    erbdihed_omp  (1:nthread) = 0.0_dp
    elec_omp      (1:nthread) = 0.0_dp
    evdw_omp      (1:nthread) = 0.0_dp
    eposi_omp     (1:nthread) = 0.0_dp

    ! nonbonded energy
    !
    select case (boundary%type)

    case (BoundaryTypePBC)

      if (enefunc%pme_use) then

        call compute_energy_nonbond_pme_short_fep(domain, enefunc, pairlist, &
                              npt, nonb_ene, coord_pbc,  &
                              force_omp, force_pbc, &
                              virial_cell, virial_omp, elec_omp, evdw_omp)

      else

        call error_msg( &
        'Compute_Energy_Gro_Amber_Short>  Multi-step is available only with PME')

      end if

    case default

      call error_msg('Compute_Energy_Gro_Amber_Short> Unknown boundary condition')

    end select

    ! bond energy
    !
    call compute_energy_bond_fep(domain, enefunc, coord, &
                             force_omp, ebond_omp)

    ! angle energy
    !
    call compute_energy_angle_fep(domain, enefunc, coord, &
                              force_omp, eangle_omp, eurey_omp)

    ! dihedral energy
    !
    call compute_energy_dihed_fep(domain, enefunc, coord, &
                              force_omp, edihed_omp)
  
    ! Ryckaert-Bellemans dihedral energy
    !
    call compute_energy_rb_dihed_fep(domain, enefunc, coord, &
                                 force_omp, erbdihed_omp)

    ! 1-4 interaction with linear table
    !
    call compute_energy_nonbond14_table_linear_fep(domain, enefunc, force_omp, &
                                               elec_omp, evdw_omp)
    call pme_bond_corr_linear_fep(domain, enefunc, force_omp, elec_omp)

    ! restraint energy
    !
    if (enefunc%restraint) &
      call compute_energy_restraints_fep(.true., .true., domain, boundary, enefunc, coord, &
                              force_omp, virial_omp, virial_ext_omp,         &
                              eposi_omp, energy%restraint_rmsd,              &
                              energy%rmsd, energy%restraint_distance,        &
                              energy%restraint_emfit, energy%emcorr)

    ! finish GPU
    !
#ifdef USE_GPU
    if (enefunc%pme_use) then
      call compute_energy_nonbond_pme_wait( &
           domain, enefunc, pairlist, boundary, &
           npt, nonb_ene, coord_pbc, force_omp, &
           force_pbc, virial_cell, virial_omp,  &
           elec_omp, evdw_omp)
    end if
#endif

    ! gather values
    !
    !$omp parallel do default(shared) private(id, i, ix, force_tmp) &
    !$omp schedule(dynamic,1)
    do i = 1, ncell
      do ix = 1, domain%num_atom(i)
        force_tmp(1:3) = force_omp(1:3,ix,i,1)
        force_tmp(1:3) = force_tmp(1:3) + force_pbc(1:3,ix,i,1)
#ifndef USE_GPU
        do id = 2, nthread
          force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,id)
          force_tmp(1:3) = force_tmp(1:3) + force_pbc(1:3,ix,i,id)
        end do
#else
        do id = 2, nthread
          force_tmp(1:3) = force_tmp(1:3) + force_omp(1:3,ix,i,id)
        end do
#endif
        force(1:3,ix,i) = force_tmp(1:3)
      end do
    end do
    !$omp end parallel do

    return

  end subroutine compute_energy_gro_amber_short_fep

end module sp_energy_mod
