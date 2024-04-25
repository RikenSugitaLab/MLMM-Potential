!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_alchemy_mod
!> @brief   Setup FEP calculation
!! @authors Hiraku Oshima (HO), Nobuhiko Kato (NK)
!
!  (c) Copyright 2019 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_alchemy_mod

  use sp_energy_mod
  use sp_energy_str_mod
  use sp_output_str_mod
  use sp_output_mod
  use fileio_control_mod
  use string_mod
  use random_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use sp_enefunc_str_mod
  use sp_alchemy_str_mod
  use sp_fep_utils_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! structures
  type, public :: s_alch_info
    integer                :: fepout_period  = 100
    integer                :: equilsteps     = 0
    character(MaxLine)     :: singleA        = ''
    character(MaxLine)     :: singleB        = ''
    character(MaxLine)     :: dualA          = ''
    character(MaxLine)     :: dualB          = ''
    real(wp)               :: sc_alpha       = 5.0_wp
    real(wp)               :: sc_beta        = 0.5_wp
    character(MaxLine)     :: lambljA        = ''
    character(MaxLine)     :: lambljB        = ''
    character(MaxLine)     :: lambelA        = ''
    character(MaxLine)     :: lambelB        = ''
    character(MaxLine)     :: lambbondA      = ''
    character(MaxLine)     :: lambbondB      = ''
    character(MaxLine)     :: lambrest       = ''
    integer                :: fep_direction  = FEP_Bothsides
    integer                :: fep_topology   = FEP_Hybrid
    integer                :: fep_md_type    = FEP_Serial
    integer                :: ref_lambid     = 0
  end type s_alch_info

  ! subroutines
  public  :: show_ctrl_alchemy
  public  :: read_ctrl_alchemy
  public  :: setup_alchemy_remd
  public  :: setup_alchemy_min
  public  :: setup_alchemy_md

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_alchemy
  !> @brief        show ALCHEMY section usage
  !! @authors      HO
  !! @param[in]    show_all : show all usage or not
  !! @param[in]    run_mode : run mode string : "md", "min"
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine show_ctrl_alchemy(show_all, run_mode)

    ! formal arguments
    logical,                 intent(in)    :: show_all
    character(*),            intent(in)    :: run_mode

    select case (run_mode)

    case ('min')

      write(MsgOut,'(A)') '# [ALCHEMY]'
      write(MsgOut,'(A)') '# fep_topology   = Hybrid # Fep topology [Hybrid, Dual]'
      write(MsgOut,'(A)') '# singleA        = 1      # Group number for singleA atoms'
      write(MsgOut,'(A)') '# singleB        = 2      # Group number for singleB atoms'
      write(MsgOut,'(A)') '# dualA          = 3      # Group number for dualA atoms'
      write(MsgOut,'(A)') '# dualB          = 4      # Group number for dualB atoms'
      write(MsgOut,'(A)') '# sc_alpha       = 5.0    # LJ soft-core parameter'
      write(MsgOut,'(A)') '# sc_beta        = 0.5    # Electrostatic soft-core parameter'
      write(MsgOut,'(A)') '# lambljA        = 1.0    # lambda for LJ in region A'
      write(MsgOut,'(A)') '# lambljB        = 0.0    # lambda for LJ in region B'
      write(MsgOut,'(A)') '# lambelA        = 1.0    # lambda for electrostatic in region A'
      write(MsgOut,'(A)') '# lambelB        = 0.0    # lambda for electrostatic in region B'
      write(MsgOut,'(A)') '# lambbondA      = 1.0    # lambda for internal bonds in singleA'
      write(MsgOut,'(A)') '# lambbondB      = 0.0    # lambda for internal bonds in singleB'
      write(MsgOut,'(A)') '# lambrest       = 1.0    # lambda for restraint energy'
      write(MsgOut,'(A)') ''

    case ('md', 'remd')

      write(MsgOut,'(A)') '# [ALCHEMY]'
      write(MsgOut,'(A)') '# fep_direction  = Bothsides              # Fep direction [Bothsides, Forward, Reverse]'
      write(MsgOut,'(A)') '# fep_topology   = Hybrid                 # Fep topology [Hybrid, Dual]'
      write(MsgOut,'(A)') '# fep_md_type    = Serial                 # FEP-MD or FEP-REMD [Serial, Single, Parallel]'
      write(MsgOut,'(A)') '# singleA        = 1                      # Group number for singleA atoms'
      write(MsgOut,'(A)') '# singleB        = 2                      # Group number for singleB atoms'
      write(MsgOut,'(A)') '# dualA          = 3                      # Group number for dualA atoms'
      write(MsgOut,'(A)') '# dualB          = 4                      # Group number for dualB atoms'
      write(MsgOut,'(A)') '# fepout_period  = 100                    # output period for energy difference'
      write(MsgOut,'(A)') '# equilsteps     = 0                      # number of equilibration steps'
      write(MsgOut,'(A)') '# sc_alpha       = 5.0                    # LJ soft-core parameter'
      write(MsgOut,'(A)') '# sc_beta        = 0.5                    # Electrostatic soft-core parameter'
      write(MsgOut,'(A)') '# lambljA        = 1.0 0.75 0.5 0.25 0.0  # lambda for LJ in region A'
      write(MsgOut,'(A)') '# lambljB        = 0.0 0.25 0.5 0.75 1.0  # lambda for LJ in region B'
      write(MsgOut,'(A)') '# lambelA        = 1.0 0.75 0.5 0.25 0.0  # lambda for electrostatic in region A'
      write(MsgOut,'(A)') '# lambelB        = 0.0 0.25 0.5 0.75 1.0  # lambda for electrostatic in region B'
      write(MsgOut,'(A)') '# lambbondA      = 1.0 0.75 0.5 0.25 0.0  # lambda for internal bonds in singleA'
      write(MsgOut,'(A)') '# lambbondB      = 0.0 0.25 0.5 0.75 1.0  # lambda for internal bonds in singleB'
      write(MsgOut,'(A)') '# lambrest       = 1.0 1.0  1.0 1.0  1.0  # lambda for restraint energy'
      write(MsgOut,'(A)') '# ref_lambid     = 0                      # Reference window id for single FEP-MD'
      write(MsgOut,'(A)') ''

    end select

    return

  end subroutine show_ctrl_alchemy

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_alchemy
  !> @brief        read ALCHEMY section in the control file for FEP calculation
  !! @authors      HO, NK
  !! @param[in]    handle   : unit number of control file
  !! @param[out]   alch_info : ALCHEMY section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_alchemy(handle, alch_info)

    ! parameters
    character(*),            parameter     :: Section = 'ALCHEMY'

    ! formal arguments
    integer,                 intent(in)    :: handle
    type(s_alch_info),       intent(inout) :: alch_info

    ! local variables
    integer                  :: i, trep
    character(30)            :: partmp


    ! read parameters from control file
    !
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_integer(handle, Section, 'fepout_period', &
                               alch_info%fepout_period)
    call read_ctrlfile_string (handle, Section, 'singleA',  &
                               alch_info%singleA)
    call read_ctrlfile_string (handle, Section, 'singleB',  &
                               alch_info%singleB)
    call read_ctrlfile_string (handle, Section, 'dualA',  &
                               alch_info%dualA)
    call read_ctrlfile_string (handle, Section, 'dualB',  &
                               alch_info%dualB)
    call read_ctrlfile_integer(handle, Section, 'equilsteps',     &
                               alch_info%equilsteps)
    call read_ctrlfile_real   (handle, Section, 'sc_alpha',      &
                               alch_info%sc_alpha)
    call read_ctrlfile_real   (handle, Section, 'sc_beta',      &
                               alch_info%sc_beta)
    call read_ctrlfile_string (handle, Section, 'lambljA',  &
                               alch_info%lambljA)
    call read_ctrlfile_string (handle, Section, 'lambljB',  &
                               alch_info%lambljB)
    call read_ctrlfile_string (handle, Section, 'lambelA',  &
                               alch_info%lambelA)
    call read_ctrlfile_string (handle, Section, 'lambelB',  &
                               alch_info%lambelB)
    call read_ctrlfile_string (handle, Section, 'lambbondA',  &
                               alch_info%lambbondA)
    call read_ctrlfile_string (handle, Section, 'lambbondB',  &
                               alch_info%lambbondB)
    call read_ctrlfile_string (handle, Section, 'lambrest',  &
                               alch_info%lambrest)
    call read_ctrlfile_type   (handle, Section, 'fep_direction',    &
                               alch_info%fep_direction, FEPDirectionTypes)
    call read_ctrlfile_type   (handle, Section, 'fep_topology',    &
                               alch_info%fep_topology, FEPTopologyTypes)
    call read_ctrlfile_type   (handle, Section, 'fep_md_type',    &
                               alch_info%fep_md_type, FEPMDTypes)
    call read_ctrlfile_integer(handle, Section, 'ref_lambid',     &
                               alch_info%ref_lambid)

    call end_ctrlfile_section(handle)


    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Read_Ctrl_ALCHEMY> Alchemy information'
      write(MsgOut,'(A,A)')    '  fep_direction  = ', &
                                  FEPDirectionTypes(alch_info%fep_direction)
      write(MsgOut,'(A,A)')    '  fep_topology  = ', &
                                  FEPTopologyTypes(alch_info%fep_topology)
      write(MsgOut,'(A,A)')    '  fep_md_type  = ', &
                                  FEPMDTypes(alch_info%fep_md_type)
      write(MsgOut,'(A,A)')    '  singleA group = ', &
                                  trim(alch_info%singleA)
      write(MsgOut,'(A,A)')    '  singleB group = ', &
                                  trim(alch_info%singleB)
      write(MsgOut,'(A,A)')    '  dualA group = ', &
                                  trim(alch_info%dualA)
      write(MsgOut,'(A,A)')    '  dualB group = ', &
                                  trim(alch_info%dualB)
      write(MsgOut,'(A,I10)')  '  fepout_period = ', &
                                  alch_info%fepout_period
      write(MsgOut,'(A,I10)')  '  equilsteps     = ', &
                                  alch_info%equilsteps
      write(MsgOut,'(A,F5.2)') '  sc_alpha       = ', &
                                  alch_info%sc_alpha
      write(MsgOut,'(A,F5.2)') '  sc_beta        = ', &
                                  alch_info%sc_beta
      write(MsgOut,'(A,A)')    '  lambljA  = ', &
                                  trim(alch_info%lambljA)
      write(MsgOut,'(A,A)')    '  lambljB  = ', &
                                  trim(alch_info%lambljB)
      write(MsgOut,'(A,A)')    '  lambelA  = ', &
                                  trim(alch_info%lambelA)
      write(MsgOut,'(A,A)')    '  lambelB  = ', &
                                  trim(alch_info%lambelB)
      write(MsgOut,'(A,A)')    '  lambbondA  = ', &
                                  trim(alch_info%lambbondA)
      write(MsgOut,'(A,A)')    '  lambbondB  = ', &
                                  trim(alch_info%lambbondB)
      write(MsgOut,'(A,A)')    '  lambrest = ', &
                                  trim(alch_info%lambrest)
      write(MsgOut,'(A)') ''
    end if

    return

  end subroutine read_ctrl_alchemy

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_alchemy_remd
  !> @brief        setup ALCHEMY for FEP calculations
  !! @authors      HO
  !! @param[in]    alch_info  : ALCHEMY section control parameters information
  !! @param[inout] alchemy    : ALCHEMY information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_alchemy_remd(alch_info, enefunc, alchemy)

    ! formal arguments
    type(s_alch_info),       intent(in)    :: alch_info
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_alchemy),         intent(inout) :: alchemy

    ! local variables
    integer                  :: i, j, k, m
    real(wp), allocatable    :: lamb_temp(:)

    ! initialize
    !
    alchemy%sc_alpha         = alch_info%sc_alpha
    alchemy%sc_beta          = alch_info%sc_beta
    alchemy%num_fep_neighbor = 2
    alchemy%fep_direction    = FEP_Bothsides
    alchemy%fep_topology     = alch_info%fep_topology
    alchemy%equilsteps       = alch_info%equilsteps
    if (main_rank) then
      if (alch_info%fep_md_type /= FEP_Parallel) then
        write(MsgOut,'(A)') 'Setup_Alchemy_Remd> fep_md_type was changed to Parallel.'
        write(MsgOut,'(A)') 'If you need serial FEP MD calculations, please remove [REMD] section.'
        write(MsgOut,'(A)') ''
      end if
    end if
    alchemy%fep_md_type       = FEP_Parallel

    ! setup lambda of each window for FEP
    ! split alch_info%lamb?_? (char) into lamb?_? (real)
    !
    alchemy%num_fep_windows  = split_num(alch_info%lambljA)

    call alloc_alchemy(alchemy, alchemy%num_fep_windows)

    allocate(lamb_temp(alchemy%num_fep_windows))

    call split(alchemy%num_fep_windows, alchemy%num_fep_windows, alch_info%lambljA, lamb_temp)
    do j = 1, alchemy%num_fep_windows
      alchemy%lambljA(j) = lamb_temp(j)
    end do

    if (split_num(alch_info%lambljB) /= (alchemy%num_fep_windows)) &
        call error_msg('Setup_Alchemy_Remd> &
        "num_fep_windows" and number of data in "lambljB" are inconsistent')
    call split(alchemy%num_fep_windows, alchemy%num_fep_windows, alch_info%lambljB, lamb_temp)
    do j = 1, alchemy%num_fep_windows
      alchemy%lambljB(j) = lamb_temp(j)
    end do

    if (split_num(alch_info%lambelA) /= (alchemy%num_fep_windows)) &
        call error_msg('Setup_Alchemy_Remd> &
        "num_fep_windows" and number of data in "lambelA" are inconsistent')
    call split(alchemy%num_fep_windows, alchemy%num_fep_windows, alch_info%lambelA, lamb_temp)
    do j = 1, alchemy%num_fep_windows
      alchemy%lambelA(j) = lamb_temp(j)
    end do

    if (split_num(alch_info%lambelB) /= (alchemy%num_fep_windows)) &
        call error_msg('Setup_Alchemy_Remd> &
        "num_fep_windows" and number of data in "lambelB" are inconsistent')
    call split(alchemy%num_fep_windows, alchemy%num_fep_windows, alch_info%lambelB, lamb_temp)
    do j = 1, alchemy%num_fep_windows
      alchemy%lambelB(j) = lamb_temp(j)
    end do

    if (split_num(alch_info%lambbondA) /= (alchemy%num_fep_windows)) &
        call error_msg('Setup_Alchemy_Remd> &
        "num_fep_windows" and number of data in "lambbondA" are inconsistent')
    call split(alchemy%num_fep_windows, alchemy%num_fep_windows, alch_info%lambbondA, lamb_temp)
    do j = 1, alchemy%num_fep_windows
      alchemy%lambbondA(j) = lamb_temp(j)
    end do

    if (split_num(alch_info%lambbondB) /= (alchemy%num_fep_windows)) &
        call error_msg('Setup_Alchemy_Remd> &
        "num_fep_windows" and number of data in "lambbondB" are inconsistent')
    call split(alchemy%num_fep_windows, alchemy%num_fep_windows, alch_info%lambbondB, lamb_temp)
    do j = 1, alchemy%num_fep_windows
      alchemy%lambbondB(j) = lamb_temp(j)
    end do

    if (alch_info%lambrest == '') then
      if (main_rank) then
        write(MsgOut,'(A)') 'Setup_Alchemy_Remd> All values of lambrest &
          are set to 1.0 because lambrest is not specified.'
        write(MsgOut,'(A)') ''
      end if
      do j = 1, alchemy%num_fep_windows
        alchemy%lambrest(j) = 1.0_wp
      end do
    else
      if (split_num(alch_info%lambrest) /= (alchemy%num_fep_windows)) &
          call error_msg('Setup_Alchemy_Remd> &
          "num_fep_windows" and number of data in "lambrest" are inconsistent')
      call split(alchemy%num_fep_windows, alchemy%num_fep_windows, alch_info%lambrest, lamb_temp)
      do j = 1, alchemy%num_fep_windows
        alchemy%lambrest(j) = lamb_temp(j)
      end do
    end if

    deallocate(lamb_temp)

    enefunc%sc_alpha           = alchemy%sc_alpha
    enefunc%sc_beta            = alchemy%sc_beta
    enefunc%num_fep_neighbor   = alchemy%num_fep_neighbor
    enefunc%fep_direction      = alchemy%fep_direction
    enefunc%lambljA            = alchemy%lambljA(1)
    enefunc%lambljB            = alchemy%lambljB(1)
    enefunc%lambelA            = alchemy%lambelA(1)
    enefunc%lambelB            = alchemy%lambelB(1)
    enefunc%lambbondA          = alchemy%lambbondA(1)
    enefunc%lambbondB          = alchemy%lambbondB(1)
    enefunc%lambrest           = alchemy%lambrest(1)

    ! Initialize table of lambda and softcore in FEP
    call set_lambda_table_fep(enefunc)

    ! write the summary of setup
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Setup_Alchemy_Remd> Alchemy information'
      write(MsgOut,'(A,I10)')   '  num_fep_windows = ', alchemy%num_fep_windows
      write(MsgOut,'(A,F10.4)') '  sc_alpha        = ', alchemy%sc_alpha
      write(MsgOut,'(A,F10.4)') '  sc_beta         = ', alchemy%sc_beta
      write(MsgOut,'(A,I10)')   '  equilsteps      = ', alchemy%equilsteps
      write(MsgOut,'(A,A)')     '  fep_direction   = ', FEPDirectionTypes(alchemy%fep_direction)
      write(MsgOut,'(A,A)')     '  fep_md_type     = ', FEPMDTypes(alchemy%fep_md_type)
      write(MsgOut,'(A)') ''
    end if

    return

  end subroutine setup_alchemy_remd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_alchemy_min
  !> @brief        setup ALCHEMY for FEP calculations
  !! @authors      HO
  !! @param[in]    alch_info  : ALCHEMY section control parameters information
  !! @param[inout] alchemy    : ALCHEMY information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_alchemy_min(alch_info, enefunc, alchemy)

    ! formal arguments
    type(s_alch_info),       intent(in)    :: alch_info
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_alchemy),         intent(inout) :: alchemy

    ! local variables
    integer                  :: i, j, k, m
    real(wp), allocatable    :: lamb_temp(:)

    if (main_rank) then
    end if

    ! initialize
    !
    alchemy%equilsteps       = alch_info%equilsteps
    alchemy%sc_alpha         = alch_info%sc_alpha
    alchemy%sc_beta          = alch_info%sc_beta
    alchemy%num_fep_neighbor = 0
    alchemy%fep_direction    = FEP_Nodirection
    alchemy%fep_topology     = alch_info%fep_topology

    ! setup lambda of each window for FEP
    ! split alch_info%lamb?_? (char) into lamb?_? (real)
    !
    alchemy%num_fep_windows  = split_num(alch_info%lambljA)
    if (alchemy%num_fep_windows < 1) &
      call error_msg('Setup_Alchemy_Min> number of data in &
                      "lambljA" should be > 0')

    if (alchemy%num_fep_windows > 1) then
        call alloc_alchemy(alchemy, alchemy%num_fep_windows)

        allocate(lamb_temp(alchemy%num_fep_windows))

        call split(alchemy%num_fep_windows, alchemy%num_fep_windows, alch_info%lambljA, lamb_temp)
        alchemy%lambljA(1) = lamb_temp(1)

        if (split_num(alch_info%lambljB) < 1) &
          call error_msg('Setup_Alchemy_Min> number of data in &
                          "lambljB" should be > 0')
        call split(alchemy%num_fep_windows, alchemy%num_fep_windows, alch_info%lambljB, lamb_temp)
        alchemy%lambljB(1) = lamb_temp(1)

        if (split_num(alch_info%lambelA) < 1) &
          call error_msg('Setup_Alchemy_Min> number of data in &
                          "lambelA" should be > 0')
        call split(alchemy%num_fep_windows, alchemy%num_fep_windows, alch_info%lambelA, lamb_temp)
        alchemy%lambelA(1) = lamb_temp(1)

        if (split_num(alch_info%lambelB) < 1) &
          call error_msg('Setup_Alchemy_Min> number of data in &
                          "lambelB" should be > 0')
        call split(alchemy%num_fep_windows, alchemy%num_fep_windows, alch_info%lambelB, lamb_temp)
        alchemy%lambelB(1) = lamb_temp(1)

        if (split_num(alch_info%lambbondA) < 1) &
          call error_msg('Setup_Alchemy_Min> number of data in &
                          "lambbondA" should be > 0')
        call split(alchemy%num_fep_windows, alchemy%num_fep_windows, alch_info%lambbondA, lamb_temp)
        alchemy%lambbondA(1) = lamb_temp(1)

        if (split_num(alch_info%lambbondB) < 1) &
          call error_msg('Setup_Alchemy_Min> number of data in &
                          "lambbondB" should be > 0')
        call split(alchemy%num_fep_windows, alchemy%num_fep_windows, alch_info%lambbondB, lamb_temp)
        alchemy%lambbondB(1) = lamb_temp(1)

        if (alch_info%lambrest == '') then
          if (main_rank) then
            write(MsgOut,'(A)') 'Setup_Alchemy_Min> lambrest &
              is set to 1.0 because lambrest is not specified.'
            write(MsgOut,'(A)') ''
          end if
          alchemy%lambrest(1) = 1.0_wp
        else
          if (split_num(alch_info%lambrest) < 1) &
            call error_msg('Setup_Alchemy_Min> number of data in &
                            "lambrest" should be > 0')
          call split(alchemy%num_fep_windows, alchemy%num_fep_windows, alch_info%lambrest, lamb_temp)
          alchemy%lambrest(1) = lamb_temp(1)
        end if

        deallocate(lamb_temp)
    else
        alchemy%num_fep_windows = 1
        call alloc_alchemy(alchemy, alchemy%num_fep_windows)
        read(alch_info%lambljA,*) alchemy%lambljA(1)
        read(alch_info%lambljB,*) alchemy%lambljB(1)
        read(alch_info%lambelA,*) alchemy%lambelA(1)
        read(alch_info%lambelB,*) alchemy%lambelB(1)
        read(alch_info%lambbondA,*) alchemy%lambbondA(1)
        read(alch_info%lambbondB,*) alchemy%lambbondB(1)
        if (alch_info%lambrest == '') then
          if (main_rank) then
            write(MsgOut,'(A)') 'Setup_Alchemy_Min> lambrest &
              is set to 1.0 because lambrest is not specified.'
            write(MsgOut,'(A)') ''
          end if
          alchemy%lambrest(1) = 1.0_wp
        else
          read(alch_info%lambrest,*) alchemy%lambrest(1)
        end if
    end if

    alchemy%num_fep_windows = 1

    enefunc%sc_alpha         = alchemy%sc_alpha
    enefunc%sc_beta          = alchemy%sc_beta
    enefunc%num_fep_neighbor = alchemy%num_fep_neighbor
    enefunc%fep_direction    = alchemy%fep_direction
    enefunc%lambljA          = alchemy%lambljA(1)
    enefunc%lambljB          = alchemy%lambljB(1)
    enefunc%lambelA          = alchemy%lambelA(1)
    enefunc%lambelB          = alchemy%lambelB(1)
    enefunc%lambbondA        = alchemy%lambbondA(1)
    enefunc%lambbondB        = alchemy%lambbondB(1)
    enefunc%lambrest         = alchemy%lambrest(1)

    ! Initialize table of lambda and softcore in FEP
    call set_lambda_table_fep(enefunc)

    ! write the summary of setup
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Setup_Alchemy_Min> Alchemy information'
      write(MsgOut,'(A,F10.5)') '  sc_alpha       = ', alchemy%sc_alpha
      write(MsgOut,'(A,F10.5)') '  sc_beta        = ', alchemy%sc_beta
      write(MsgOut,'(A,F10.5)') '  lambljA  = ', alchemy%lambljA(1)
      write(MsgOut,'(A,F10.5)') '  lambljB  = ', alchemy%lambljB(1)
      write(MsgOut,'(A,F10.5)') '  lambelA  = ', alchemy%lambelA(1)
      write(MsgOut,'(A,F10.5)') '  lambelB  = ', alchemy%lambelB(1)
      write(MsgOut,'(A,F10.5)') '  lambbondA  = ', alchemy%lambbondA(1)
      write(MsgOut,'(A,F10.5)') '  lambbondB  = ', alchemy%lambbondB(1)
      write(MsgOut,'(A,F10.5)') '  lambrest = ', alchemy%lambrest(1)
      write(MsgOut,'(A,A)')     '  fep_direction  = ', FEPDirectionTypes(alchemy%fep_direction)
      write(MsgOut,'(A)') ''
    end if

    return

  end subroutine setup_alchemy_min

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_alchemy_md
  !> @brief        setup ALCHEMY for FEP calculations
  !! @authors      HO
  !! @param[in]    alch_info  : ALCHEMY section control parameters information
  !! @param[inout] alchemy    : ALCHEMY information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_alchemy_md(alch_info, enefunc, alchemy)

    ! formal arguments
    type(s_alch_info),       intent(in)    :: alch_info
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_alchemy),         intent(inout) :: alchemy

    ! local variables
    integer                  :: i, j, k, m
    real(wp), allocatable    :: lamb_temp(:)

    ! initialize
    !
    alchemy%sc_alpha      = alch_info%sc_alpha
    alchemy%sc_beta       = alch_info%sc_beta
    alchemy%fep_direction = alch_info%fep_direction
    alchemy%fep_topology  = alch_info%fep_topology
    alchemy%equilsteps    = alch_info%equilsteps
    if (alch_info%fep_md_type == FEP_Parallel) then
      if (main_rank) then
        write(MsgOut,'(A)') 'Setup_Alchemy_Md> fep_md_type was changed to Serial.'
        write(MsgOut,'(A)') 'If you need FEP REMD, please add [REMD] section.'
        write(MsgOut,'(A)') ''
      end if
      alchemy%fep_md_type = FEP_Serial
    else if (alch_info%fep_md_type == FEP_Serial) then
      alchemy%fep_md_type = FEP_Serial
    else if (alch_info%fep_md_type == FEP_Single) then
      alchemy%fep_md_type = FEP_Single
      if (alch_info%ref_lambid == 0) then
        if (main_rank) then
          call error_msg('Setup_Alchemy_Md> &
            &ref_lambid should be set to window id (> 0) if fep_md_type is Single.')
        end if
      end if
      alchemy%ref_lambid = alch_info%ref_lambid
    end if
    if (alch_info%fep_direction == FEP_Bothsides) then
      alchemy%num_fep_neighbor = 2
    else if (alch_info%fep_direction == FEP_Forward) then
      alchemy%num_fep_neighbor = 1
    else if (alch_info%fep_direction == FEP_Reverse) then
      alchemy%num_fep_neighbor = 1
    else
      if (main_rank) then
        call error_msg('Setup_Alchemy_Md> &
          fep_direction should be Bothsides, Forward, or Reverse.')
      end if
    end if

    ! setup lambda of each window for FEP
    ! split alch_info%lamb?_? (char) into lamb?_? (real)
    !
    alchemy%num_fep_windows  = split_num(alch_info%lambljA)

    call alloc_alchemy(alchemy, alchemy%num_fep_windows)

    allocate(lamb_temp(alchemy%num_fep_windows))

    call split(alchemy%num_fep_windows, alchemy%num_fep_windows, alch_info%lambljA, lamb_temp)
    do j = 1, alchemy%num_fep_windows
      alchemy%lambljA(j) = lamb_temp(j)
    end do

    if (split_num(alch_info%lambljB) /= (alchemy%num_fep_windows)) &
        call error_msg('Setup_Alchemy_Md> &
        "num_fep_windows" and number of data in "lambljB" are inconsistent')
    call split(alchemy%num_fep_windows, alchemy%num_fep_windows, alch_info%lambljB, lamb_temp)
    do j = 1, alchemy%num_fep_windows
      alchemy%lambljB(j) = lamb_temp(j)
    end do

    if (split_num(alch_info%lambelA) /= (alchemy%num_fep_windows)) &
        call error_msg('Setup_Alchemy_Md> &
        "num_fep_windows" and number of data in "lambelA" are inconsistent')
    call split(alchemy%num_fep_windows, alchemy%num_fep_windows, alch_info%lambelA, lamb_temp)
    do j = 1, alchemy%num_fep_windows
      alchemy%lambelA(j) = lamb_temp(j)
    end do

    if (split_num(alch_info%lambelB) /= (alchemy%num_fep_windows)) &
        call error_msg('Setup_Alchemy_Md> &
        "num_fep_windows" and number of data in "lambelB" are inconsistent')
    call split(alchemy%num_fep_windows, alchemy%num_fep_windows, alch_info%lambelB, lamb_temp)
    do j = 1, alchemy%num_fep_windows
      alchemy%lambelB(j) = lamb_temp(j)
    end do

    if (split_num(alch_info%lambbondA) /= (alchemy%num_fep_windows)) &
        call error_msg('Setup_Alchemy_Md> &
        "num_fep_windows" and number of data in "lambbondA" are inconsistent')
    call split(alchemy%num_fep_windows, alchemy%num_fep_windows, alch_info%lambbondA, lamb_temp)
    do j = 1, alchemy%num_fep_windows
      alchemy%lambbondA(j) = lamb_temp(j)
    end do

    if (split_num(alch_info%lambbondB) /= (alchemy%num_fep_windows)) &
        call error_msg('Setup_Alchemy_Md> &
        "num_fep_windows" and number of data in "lambbondB" are inconsistent')
    call split(alchemy%num_fep_windows, alchemy%num_fep_windows, alch_info%lambbondB, lamb_temp)
    do j = 1, alchemy%num_fep_windows
      alchemy%lambbondB(j) = lamb_temp(j)
    end do

    if (alch_info%lambrest == '') then
      if (main_rank) then
        write(MsgOut,'(A)') 'Setup_Alchemy_Md> All values of lambrest &
          are set to 1.0 because lambrest is not specified.'
        write(MsgOut,'(A)') ''
      end if
      do j = 1, alchemy%num_fep_windows
        alchemy%lambrest(j) = 1.0_wp
      end do
    else
      if (split_num(alch_info%lambrest) /= (alchemy%num_fep_windows)) &
          call error_msg('Setup_Alchemy_Md> &
          "num_fep_windows" and number of data in "lambrest" are inconsistent')
      call split(alchemy%num_fep_windows, alchemy%num_fep_windows, alch_info%lambrest, lamb_temp)
      do j = 1, alchemy%num_fep_windows
        alchemy%lambrest(j) = lamb_temp(j)
      end do
    end if

    deallocate(lamb_temp)

    enefunc%sc_alpha           = alchemy%sc_alpha
    enefunc%sc_beta            = alchemy%sc_beta
    enefunc%num_fep_neighbor   = alchemy%num_fep_neighbor
    enefunc%fep_direction      = alchemy%fep_direction
    if (alchemy%fep_md_type == FEP_Single) then
      enefunc%lambljA            = alchemy%lambljA(alchemy%ref_lambid)
      enefunc%lambljB            = alchemy%lambljB(alchemy%ref_lambid)
      enefunc%lambelA            = alchemy%lambelA(alchemy%ref_lambid)
      enefunc%lambelB            = alchemy%lambelB(alchemy%ref_lambid)
      enefunc%lambbondA          = alchemy%lambbondA(alchemy%ref_lambid)
      enefunc%lambbondB          = alchemy%lambbondB(alchemy%ref_lambid)
      enefunc%lambrest           = alchemy%lambrest(alchemy%ref_lambid)
    else
      enefunc%lambljA            = alchemy%lambljA(1)
      enefunc%lambljB            = alchemy%lambljB(1)
      enefunc%lambelA            = alchemy%lambelA(1)
      enefunc%lambelB            = alchemy%lambelB(1)
      enefunc%lambbondA          = alchemy%lambbondA(1)
      enefunc%lambbondB          = alchemy%lambbondB(1)
      enefunc%lambrest           = alchemy%lambrest(1)
    end if

    ! Initialize table of lambda and softcore in FEP
    call set_lambda_table_fep(enefunc)

    ! write the summary of setup
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Setup_Alchemy_Md> Alchemy information'
      write(MsgOut,'(A,I10)')   '  num_fep_windows = ', alchemy%num_fep_windows
      write(MsgOut,'(A,F10.4)') '  sc_alpha        = ', alchemy%sc_alpha
      write(MsgOut,'(A,F10.4)') '  sc_beta         = ', alchemy%sc_beta
      write(MsgOut,'(A,I10)')   '  equilsteps      = ', alchemy%equilsteps
      write(MsgOut,'(A,A)')     '  fep_direction   = ', FEPDirectionTypes(alchemy%fep_direction)
      write(MsgOut,'(A,A)')     '  fep_md_type     = ', FEPMDTypes(alchemy%fep_md_type)
      if (alchemy%fep_md_type == FEP_Single) then
        write(MsgOut,'(A,I10)') '  ref_lambid      = ', alchemy%ref_lambid
      end if
      write(MsgOut,'(A)') ''
      write(MsgOut,'(A)') '  FEP Windows'
      do j = 1, alchemy%num_fep_windows
          write(MsgOut,'(A,I4)') '    Window index = ', j

          write(MsgOut,'(6X,A,F10.5,5X)') &
          ' lambljA = ', alchemy%lambljA(j)

          write(MsgOut,'(6X,A,F10.5,5X)') &
          ' lambljB = ', alchemy%lambljB(j)

          write(MsgOut,'(6X,A,F10.5,5X)') &
          ' lambelA = ', alchemy%lambelA(j)

          write(MsgOut,'(6X,A,F10.5,5X)') &
          ' lambelB = ', alchemy%lambelB(j)

          write(MsgOut,'(6X,A,F10.5,5X)') &
          ' lambbondA = ', alchemy%lambbondA(j)

          write(MsgOut,'(6X,A,F10.5,5X)') &
          ' lambbondB = ', alchemy%lambbondB(j)

          write(MsgOut,'(6X,A,F10.5,5X)') &
          ' lambrest = ', alchemy%lambrest(j)

          write(MsgOut,'(A)') ''
      end do

      if (alchemy%fep_md_type == FEP_Serial) then
        write(MsgOut,'(A)') '  Serial FEP MD simulations will be performed by changing FEP windows'
      else
        write(MsgOut,'(A,I4)') '  Single FEP MD simulation will be performed only at FEP window ', alchemy%ref_lambid
      end if
      write(MsgOut,'(A)') ''
    end if

    return

  end subroutine setup_alchemy_md

end module sp_alchemy_mod
