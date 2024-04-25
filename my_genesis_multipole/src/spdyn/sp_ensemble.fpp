!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_ensemble_mod
!> @brief   setup parameters for temperature and pressure control
!! @authors Takaharu Mori (TM), Jaewoon Jung (JJ), Norio Takase (NT)
!           Tadashi Ando (TA)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_ensemble_mod

  use sp_dynamics_str_mod
  use sp_ensemble_str_mod
  use fileio_control_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_ens_info

    integer                   :: ensemble        = EnsembleNVE
    real(wp)                  :: temperature     = 298.15_wp
    real(wp)                  :: pressure        = 1.0_wp
    real(wp)                  :: gamma           = 0.0_wp
    integer                   :: tpcontrol       = TpcontrolNO

    ! parameters for Berendsen and Bussi
    real(wp)                  :: tau_t           = 5.0_wp
    real(wp)                  :: tau_p           = 5.0_wp
    real(wp)                  :: compressibility = 4.63e-5_wp

    ! parameters for Langevin
    real(wp)                  :: gamma_t         = 1.0_wp
    real(wp)                  :: gamma_p         = 0.1_wp

    ! parameters for Nose-Hoover chain
    !
    integer                   :: nhchain         = 5
    integer                   :: nhmultistep     = 4

    ! parameters for NPT
    integer                   :: isotropy        = IsotropyISO

  end type s_ens_info

  ! subroutines
  public  :: show_ctrl_ensemble
  public  :: read_ctrl_ensemble
  public  :: setup_ensemble

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_ensemble
  !> @brief        show ENSEMBLE section usage
  !! @authors      NT
  !! @param[in]    show_all : show all usage or not
  !! @param[in]    run_mode : run mode string : "md", "min"
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_ensemble(show_all, run_mode)

    ! formal arguments
    logical,                 intent(in)    :: show_all
    character(*),            intent(in)    :: run_mode


    if (show_all) then

      select case (run_mode)

      case ('md', 'remd')

        write(MsgOut,'(A)') '[ENSEMBLE]'
        write(MsgOut,'(A)') 'ensemble      = NVE       # [NVE,NVT,NPT,NPAT,NPgT]'
        write(MsgOut,'(A)') 'tpcontrol     = NO        # [NO,BERENDSEN,BUSSI,LANGEVIN]'
        write(MsgOut,'(A)') '# temperature   = 298.15    # initial and target temperature (K)'
        write(MsgOut,'(A)') '# pressure      = 1.0       # target pressure (atm)'
        write(MsgOut,'(A)') '# gamma         = 0.0       # target surface tension (dyn/cm)'
        write(MsgOut,'(A)') '# tau_t         = 5.0       # temperature coupling time (ps) in [BERENDSEN,NOSE-HOOVER,BUSSI]'
        write(MsgOut,'(A)') '# tau_p         = 5.0       # pressure coupling time (ps)    in [BERENDSEN,BUSSI]'
        write(MsgOut,'(A)') '# compressibility = 4.63e-5 # compressibility (atm-1) in [BERENDSEN]'
        write(MsgOut,'(A)') '# gamma_t       = 1.0       # thermostat friction (ps-1) in [LANGEVIN]'
        write(MsgOut,'(A)') '# gamma_p       = 0.1       # barostat friction (ps-1)   in [LANGEVIN]'
        write(MsgOut,'(A)') '# isotropy      = ISO       # [ISO,SEMI-ISO,ANISO,XY-FIXED]'
        write(MsgOut,'(A)') ' '

      end select

    else

      select case (run_mode)

      case ('md', 'remd')

        write(MsgOut,'(A)') '[ENSEMBLE]'
        write(MsgOut,'(A)') 'ensemble      = NVE       # [NVE,NVT,NPT,NPAT,NPgT]'
        write(MsgOut,'(A)') 'tpcontrol     = NO        # [NO,BERENDSEN,BUSSI,LANGEVIN]'
        write(MsgOut,'(A)') ' '

      end select

    end if

    return

  end subroutine show_ctrl_ensemble
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_ensemble
  !> @brief        read ENSEMBLE section in the control file
  !! @authors      TM, JJ
  !! @param[in]    handle   : unit number of control file
  !! @param[out]   ens_info : ENSEMBLE section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_ensemble(handle, ens_info)

    ! parameters
    character(*),            parameter     :: Section = 'Ensemble'

    ! formal arguments
    integer,                 intent(in)    :: handle
    type(s_ens_info),        intent(inout) :: ens_info



    ! read parameters from control file
    !
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_type(handle, Section, 'ensemble',          &
                            ens_info%ensemble, EnsembleTypes)
    call read_ctrlfile_type(handle, Section, 'tpcontrol',         &
                            ens_info%tpcontrol, TpcontrolTypes)
    call read_ctrlfile_real(handle, Section, 'temperature',       &
                            ens_info%temperature)
    call read_ctrlfile_real(handle, Section, 'pressure',          &
                            ens_info%pressure)
    call read_ctrlfile_real(handle, Section, 'gamma',             &
                            ens_info%gamma)
    call read_ctrlfile_real(handle, Section, 'tau_t',             &
                            ens_info%tau_t)
    call read_ctrlfile_real(handle, Section, 'tau_p',             &
                            ens_info%tau_p)
    call read_ctrlfile_real(handle, Section, 'compressibility',   &
                            ens_info%compressibility)
    call read_ctrlfile_real(handle, Section, 'gamma_t',           &
                            ens_info%gamma_t)
    call read_ctrlfile_real(handle, Section, 'gamma_p',           &
                            ens_info%gamma_p)
    call read_ctrlfile_type(handle, Section, 'isotropy',          &
                            ens_info%isotropy, IsotropyTypes)
    call read_ctrlfile_integer(handle, Section, 'nhchain',        &
                            ens_info%nhchain)
    call read_ctrlfile_integer(handle, Section, 'nhmultistep',    &
                            ens_info%nhmultistep)

    call end_ctrlfile_section(handle)

    ! error check
    !
    if (ens_info%ensemble == EnsembleNVE .and.  &
        ens_info%tpcontrol /= TpcontrolNo) &
      call error_msg('Read_Ctrl_Ensemble> Error in specifying "tpcontrol"')

    if (ens_info%ensemble == EnsembleNPAT .and.   &
        ens_info%isotropy /= IsotropyXY_Fixed)    &
      call error_msg('Read_Ctrl_Ensemble> NPAT should be XY-FIXED')

    if (ens_info%ensemble == EnsembleNPT .and.   &
        ens_info%isotropy == IsotropyXY_Fixed) &
          call error_msg('Read_Ctrl_Ensemble> Cannot specify XY-FIXED with NPT')

    if (ens_info%ensemble == EnsembleNPgT .and. &
        ens_info%isotropy /= IsotropySEMI_ISO) &
      call error_msg('Read_Ctrl_Ensemble> NPgT is available with SEMI-ISO')

    if (ens_info%ensemble == EnsembleNPAT .and.  &
        ens_info%tpcontrol /= TpcontrolLangevin) &
      call error_msg('Read_Ctrl_Ensemble> NPAT is available with Langevin')

    if (ens_info%ensemble == EnsembleNPgT .and.  &
        ens_info%tpcontrol /= TpcontrolLangevin) &
      call error_msg('Read_Ctrl_Ensemble> NPgT is available with Langevin')

    ! write parameters to MsgOut
    !
    if (main_rank) then

      write(MsgOut,'(A)') 'Read_Ctrl_Ensemble> Parameters for Ensemble'

      select case (ens_info%ensemble)

      case (EnsembleNVE)
        write(MsgOut,'(A20,A10)')   '  ensemble        = ', &
                                    trim(EnsembleTypes(ens_info%ensemble))
        write(MsgOut,'(A20,F10.3)') '  initial_temp    = ', &
                                    ens_info%temperature
        write(MsgOut,'(A20,A10)')   '  tpcontrol       = ', &
                                    trim(TpcontrolTypes(ens_info%tpcontrol))

      case (EnsembleNVT)
        write(MsgOut,'(A20,A10)')   '  ensemble        = ', &
                                    trim(EnsembleTypes(ens_info%ensemble))
        write(MsgOut,'(A20,F10.3)') '  temperature     = ', &
                                    ens_info%temperature
        write(MsgOut,'(A20,A10)')   '  tpcontrol       = ', &
                                    trim(TpcontrolTypes(ens_info%tpcontrol))

        select case (ens_info%tpcontrol)

        case (TpcontrolGauss)
          ! do nothing

        case (TpcontrolEvans)
          ! do nothing

        case (TpcontrolBerendsen)
          write(MsgOut,'(A20,F10.3)') '  tau_t           = ', ens_info%tau_t

        case (TpcontrolLangevin)
          write(MsgOut,'(A20,F10.3)') '  gamma_t         = ', ens_info%gamma_t

        case (TpcontrolNoseHoover)
          write(MsgOut,'(A20,F10.3)') '  tau_t           = ', ens_info%tau_t

        case (TpcontrolAndersen)

        case (TpcontrolBussi)
          write(MsgOut,'(A20,F10.3)') '  tau_t           = ', ens_info%tau_t

        case default
          call error_msg('Read_Ctrl_Ensemble> Error in specifying "tpcontrol"')

        end select

      case (EnsembleNPT)

        select case (ens_info%tpcontrol)

        case (TpcontrolBerendsen)
          write(MsgOut,'(A20,A10)')   '  ensemble        = ', &
                                      trim(EnsembleTypes(ens_info%ensemble))
          write(MsgOut,'(A20,F10.3)') '  temperature     = ', &
                                      ens_info%temperature
          write(MsgOut,'(A20,F10.3)') '  pressure        = ', &
                                      ens_info%pressure
          write(MsgOut,'(A20,A20)')   '  tpcontrol       = ', &
                                      trim(TpcontrolTypes(ens_info%tpcontrol))
          write(MsgOut,'(A20,F10.3)') '  tau_t           = ', &
                                      ens_info%tau_t
          write(MsgOut,'(A20,F10.3)') '  tau_p           = ', &
                                      ens_info%tau_p
          write(MsgOut,'(A20,E10.3)') '  compressibility = ', &
                                      ens_info%compressibility
          write(MsgOut,'(A20,A10)')   '  isotropy        = ', &
                                      trim(IsotropyTypes(ens_info%isotropy))

        case (TpcontrolLangevin)
          write(MsgOut,'(A20,A10)')   '  ensemble        = ', &
                                      trim(EnsembleTypes(ens_info%ensemble))
          write(MsgOut,'(A20,F10.3)') '  temperature     = ', &
                                      ens_info%temperature
          write(MsgOut,'(A20,F10.3)') '  pressure        = ', &
                                      ens_info%pressure
          write(MsgOut,'(A20,A10)')   '  tpcontrol       = ', &
                                      trim(TpcontrolTypes(ens_info%tpcontrol))
          write(MsgOut,'(A20,F10.3)') '  gamma_t         = ', &
                                      ens_info%gamma_t
          write(MsgOut,'(A20,F10.3)') '  gamma_p         = ', &
                                      ens_info%gamma_p
          write(MsgOut,'(A20,A10)')   '  isotropy        = ', &
                                      trim(IsotropyTypes(ens_info%isotropy))

        case (TpcontrolMTK)
          write(MsgOut,'(A20,A10)')   '  ensemble        = ', &
                                      trim(EnsembleTypes(ens_info%ensemble))
          write(MsgOut,'(A20,F10.3)') '  temperature     = ', &
                                      ens_info%temperature
          write(MsgOut,'(A20,F10.3)') '  pressure        = ', &
                                      ens_info%pressure
          write(MsgOut,'(A20,A10)')   '  tpcontrol       = ', &
                                      trim(TpcontrolTypes(ens_info%tpcontrol))

        case (TpcontrolBussi)
          write(MsgOut,'(A20,A10)')   '  ensemble        = ', &
                                      trim(EnsembleTypes(ens_info%ensemble))
          write(MsgOut,'(A20,F10.3)') '  temperature     = ', &
                                      ens_info%temperature
          write(MsgOut,'(A20,F10.3)') '  pressure        = ', &
                                      ens_info%pressure
          write(MsgOut,'(A20,A20)')   '  tpcontrol       = ', &
                                      trim(TpcontrolTypes(ens_info%tpcontrol))
          write(MsgOut,'(A20,F10.3)') '  tau_t           = ', &
                                      ens_info%tau_t
          write(MsgOut,'(A20,F10.3)') '  tau_p           = ', &
                                      ens_info%tau_p
          write(MsgOut,'(A20,A10)')   '  isotropy        = ', &
                                      trim(IsotropyTypes(ens_info%isotropy))

        case default
          call error_msg('Read_Ctrl_Ensemble> Error in specifying "tpcontrol"')

        end select

      case (EnsembleNPAT)

        select case (ens_info%tpcontrol)

        case (TpcontrolBerendsen)
          write(MsgOut,'(A20,A10)')   '  ensemble        = ', &
                                      trim(EnsembleTypes(ens_info%ensemble))
          write(MsgOut,'(A20,F10.3)') '  temperature     = ', &
                                      ens_info%temperature
          write(MsgOut,'(A20,F10.3)') '  pressure        = ', &
                                      ens_info%pressure
          write(MsgOut,'(A20,A10)')   '  tpcontrol       = ', &
                                      trim(TpcontrolTypes(ens_info%tpcontrol))
          write(MsgOut,'(A20,F10.3)') '  tau_t           = ', &
                                      ens_info%tau_t
          write(MsgOut,'(A20,F10.3)') '  tau_p           = ', &
                                      ens_info%tau_p
          write(MsgOut,'(A20,E10.3)') '  compressibility = ', &
                                      ens_info%compressibility
          write(MsgOut,'(A20,A10)')   '  isotropy        = ', &
                                      trim(IsotropyTypes(ens_info%isotropy))

        case (TpcontrolLangevin)
          write(MsgOut,'(A20,A10)')   '  ensemble        = ', &
                                      trim(EnsembleTypes(ens_info%ensemble))
          write(MsgOut,'(A20,F10.3)') '  temperature     = ', &
                                      ens_info%temperature
          write(MsgOut,'(A20,F10.3)') '  pressure        = ', &
                                      ens_info%pressure
          write(MsgOut,'(A20,A10)')   '  tpcontrol       = ', &
                                      trim(TpcontrolTypes(ens_info%tpcontrol))
          write(MsgOut,'(A20,F10.3)') '  gamma_t         = ', &
                                      ens_info%gamma_t
          write(MsgOut,'(A20,F10.3)') '  gamma_p         = ', &
                                      ens_info%gamma_p
          write(MsgOut,'(A20,A10)')   '  isotropy        = ', &
                                      trim(IsotropyTypes(ens_info%isotropy))

        case default
          call error_msg('Read_Ctrl_Ensemble> Error in specifying "tpcontrol"')

        end select

      case (EnsembleNPgT)

        select case (ens_info%tpcontrol)

        case (TpcontrolBerendsen)
          write(MsgOut,'(A20,A10)')   '  ensemble        = ', &
                                      trim(EnsembleTypes(ens_info%ensemble))
          write(MsgOut,'(A20,F10.3)') '  temperature     = ', &
                                      ens_info%temperature
          write(MsgOut,'(A20,F10.3)') '  pressure        = ', &
                                      ens_info%pressure
          write(MsgOut,'(A20,F10.3)') '  gamma           = ', &
                                      ens_info%gamma
          write(MsgOut,'(A20,A10)')   '  tpcontrol       = ', &
                                      trim(TpcontrolTypes(ens_info%tpcontrol))
          write(MsgOut,'(A20,F10.3)') '  tau_t           = ', &
                                      ens_info%tau_t
          write(MsgOut,'(A20,F10.3)') '  tau_p           = ', &
                                      ens_info%tau_p
          write(MsgOut,'(A20,E10.3)') '  compressibility = ', &
                                      ens_info%compressibility
          write(MsgOut,'(A20,A10)')   '  isotropy        = ', &
                                      trim(IsotropyTypes(ens_info%isotropy))

        case (TpcontrolLangevin)
          write(MsgOut,'(A20,A10)')   '  ensemble        = ', &
                                      trim(EnsembleTypes(ens_info%ensemble))
          write(MsgOut,'(A20,F10.3)') '  temperature     = ', &
                                      ens_info%temperature
          write(MsgOut,'(A20,F10.3)') '  pressure        = ', &
                                      ens_info%pressure
          write(MsgOut,'(A20,F10.3)') '  gamma           = ', &
                                      ens_info%gamma
          write(MsgOut,'(A20,A10)')   '  tpcontrol       = ', &
                                      trim(TpcontrolTypes(ens_info%tpcontrol))
          write(MsgOut,'(A20,F10.3)') '  gamma_t         = ', &
                                      ens_info%gamma_t
          write(MsgOut,'(A20,F10.3)') '  gamma_p         = ', &
                                      ens_info%gamma_p
          write(MsgOut,'(A20,A10)')   '  isotropy        = ', &
                                      trim(IsotropyTypes(ens_info%isotropy))

        case default
          call error_msg('Read_Ctrl_Ensemble> Error in specifying "tpcontrol"')

        end select

      end select

      write(MsgOut,'(A)') ' '

    end if
    
    return

  end subroutine read_ctrl_ensemble

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_ensemble
  !> @brief        setup structure of ensemble
  !! @authors      TM, JJ
  !! @param[in]    ens_info : ENSEMBLE section control parameters information
  !! @param[in]    dynamics : dynamics information
  !! @param[out]   ensemble : ensemble information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_ensemble(ens_info, dynamics, ensemble)

    ! formal arguments
    type(s_ens_info),        intent(in)    :: ens_info
    type(s_dynamics),        intent(in)    :: dynamics
    type(s_ensemble),        intent(inout) :: ensemble


    ! initialize
    !
    call init_ensemble(ensemble)


    ! setup variables
    !
    ensemble%ensemble        = ens_info%ensemble
    ensemble%temperature     = ens_info%temperature
    ensemble%pressure        = ens_info%pressure
    ensemble%gamma           = ens_info%gamma
    ensemble%tpcontrol       = ens_info%tpcontrol
    ensemble%tau_t           = ens_info%tau_t
    ensemble%tau_p           = ens_info%tau_p
    ensemble%compressibility = ens_info%compressibility
    ensemble%gamma_t         = ens_info%gamma_t
    ensemble%gamma_p         = ens_info%gamma_p
    ensemble%isotropy        = ens_info%isotropy
    ensemble%nhchain         = ens_info%nhchain
    ensemble%nhmultistep     = ens_info%nhmultistep

    if (ensemble%ensemble == EnsembleNPT  .or. &
        ensemble%ensemble == EnsembleNPAT .or. &
        ensemble%ensemble == EnsembleNPgT) then
      ensemble%use_barostat = .true.
    else
      ensemble%use_barostat = .false.
    end if
    !
    ! error check for tpcontrol & integrator
    !
    if (dynamics%integrator == IntegratorLEAP .and. &
        ensemble%tpcontrol           == TpcontrolBUSSI) then
      call error_msg('Setup_Ensemble> BUSSI is not allowed in Leap-frog')
    end if

    if (dynamics%integrator == IntegratorVRES) then
      if(((ensemble%tpcontrol  /= TpcontrolLangevin).and. &
          (ensemble%tpcontrol  /= TpcontrolBussi)) .and.  &
        ensemble%ensemble /= EnsembleNVE)                 &
      call error_msg('Setup_Ensemble> RESPA is allowed only in Langevin')
    end if
    if (dynamics%integrator == IntegratorVVER .and. &
        ensemble%use_barostat .and. &
        ensemble%tpcontrol   == TpcontrolBerendsen) then
      call error_msg('Setup_Ensemble> Berendsen+VVER is allowed only in NVT')
    end if

    if (dynamics%integrator == IntegratorVVER .and. &
        ensemble%tpcontrol  == TpcontrolBerendsen .and.  &
        ensemble%ensemble /= EnsembleNVT) then
      call error_msg('Setup_Ensemble> Berendsen works only in NVT in VVER')
    end if

    if (dynamics%integrator == IntegratorVRES .and. &
        (ensemble%ensemble == EnsembleNPAT .or.  &
        ensemble%ensemble == EnsembleNPgT)) then
      call error_msg('Setup_Ensemble> RESPA+NPAT/NPgT is not allowed')
    end if

    return

  end subroutine setup_ensemble

end module sp_ensemble_mod
