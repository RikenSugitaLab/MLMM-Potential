!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   timers_mod
!> @brief   process timer (msec)
!! @authors Yuji Sugita (YS), Takaharu Mori (TM), Naoyuki Miyashita (NM), 
!!          Norio Takase (NT), Kiyoshi Yagi (KY)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module timers_mod

  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! parameters
  integer, public, parameter :: TimerOff        = 0
  integer, public, parameter :: TimerOn         = 1

  integer, public, parameter :: TimerTotal      = 1
  integer, public, parameter :: TimerEnergy     = 2
  integer, public, parameter :: TimerBond       = 3
  integer, public, parameter :: TimerAngle      = 4
  integer, public, parameter :: TimerDihedral   = 5
  integer, public, parameter :: TimerRestraint  = 6
  integer, public, parameter :: TimerNonBond    = 7
  integer, public, parameter :: TimerPmeReal    = 8
  integer, public, parameter :: TimerPmeRecip   = 9
  integer, public, parameter :: TimerPairList   = 10
  integer, public, parameter :: TimerDynamics   = 11
  integer, public, parameter :: TimerConstraint = 12
  integer, public, parameter :: TimerUpdate     = 13
  integer, public, parameter :: TimerComm1      = 14
  integer, public, parameter :: TimerComm2      = 15
  integer, public, parameter :: TimerComm3      = 16
  integer, public, parameter :: TimerIntegrator = 17

  integer, public, parameter :: TimerTest1      = 18
  integer, public, parameter :: TimerTest2      = 19
  integer, public, parameter :: TimerTest3      = 20
  integer, public, parameter :: TimerTest4      = 21
  integer, public, parameter :: TimerTest5      = 22
  integer, public, parameter :: TimerTest6      = 23
  integer, public, parameter :: TimerTest7      = 24
  integer, public, parameter :: TimerTest8      = 25
  integer, public, parameter :: TimerTest9      = 26
  integer, public, parameter :: TimerTest10     = 27

  integer, public, parameter :: TimerRespa      = 28
  integer, public, parameter :: TimerEmfit      = 29
  integer, public, parameter :: TimerEnergyEmfit= 30
  integer, public, parameter :: TimerQMMM       = 31
  integer, public, parameter :: TimerSolvation  = 32
  integer, public, parameter :: TimerGB         = 33
  integer, public, parameter :: TimerSA         = 34

  integer,         parameter :: InvalidID       = -1
  integer,         parameter :: NumTimers = 34  !< total number of timers
  integer,         parameter :: MaxProc   = 100 !< maximum number of processes

  ! variables
  real(dp),        public    :: total_time(NumTimers)
  real(dp),        save      :: proc_time_st(NumTimers)

  ! functions
  public  :: timer
  public  :: output_time
  private :: get_unix_time
  private :: days

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    timer
  !> @brief        integrate process time
  !! @authors      TM
  !! @param[in]    mode : timer index
  !! @param[in]    set  : 0 = turn off timer; 1 = turn on timer
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine timer(mode, set)

    ! formal arguments
    integer,                 intent(in) :: mode
    integer,                 intent(in) :: set


    select case (set)

    case (TimerOn)

      proc_time_st(mode) = get_unix_time()

    case (TimerOff)

      total_time(mode) = total_time(mode) &
                       + (get_unix_time()-proc_time_st(mode))*1.0d-6

    end select

    return

  end subroutine timer

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_time
  !> @brief        output timer profile
  !! @authors      TM
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_time

    ! local variables
    real(dp)                      :: tottime_ebond, tottime_enbond
    real(dp)                      :: tottime_dynamc, avetime_setup
    integer                       :: i, j, alloc_stat, dealloc_stat
    integer                       :: num_realnodes, num_recipnodes

    real(dp)        , allocatable :: tottime(:,:), avetime(:), sumtime(:)
    real(dp)        , allocatable :: maxtime(:), mintime(:)
    logical,          allocatable :: node_real(:), node_recip(:)


    ! allocate local array
    !
    allocate(node_real (nproc_world),            &
             node_recip(nproc_world),            &
             tottime   (NumTimers, nproc_world), &
             sumtime   (NumTimers),              &
             avetime   (NumTimers),              &
             mintime   (NumTimers),              &
             maxtime   (NumTimers),              &
             stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

#ifdef HAVE_MPI_GENESIS
    if (change_world)then !added by YKL
      call mpi_allgather(total_time, NumTimers, mpi_double_precision, &
                         tottime,    NumTimers, mpi_double_precision, &
                         mpi_comm_world2, ierror)
      call mpi_allgather(real_calc, 1, mpi_logical,                   &
                         node_real, 1, mpi_logical,                   &
                         mpi_comm_world2, ierror)
      call mpi_allgather(reciprocal_calc, 1, mpi_logical,             &
                         node_recip,      1, mpi_logical,             &
                         mpi_comm_world2, ierror)
    else
      call mpi_allgather(total_time, NumTimers, mpi_double_precision, &
                         tottime,    NumTimers, mpi_double_precision, &
                         mpi_comm_world, ierror)
      call mpi_allgather(real_calc, 1, mpi_logical,                   &
                         node_real, 1, mpi_logical,                   &
                         mpi_comm_world, ierror)
      call mpi_allgather(reciprocal_calc, 1, mpi_logical,             &
                         node_recip,      1, mpi_logical,             &
                         mpi_comm_world, ierror)
    end if
#else
    tottime(1:NumTimers,1) = total_time(1:NumTimers)
    node_real(1)  = real_calc
    node_recip(1) = reciprocal_calc
#endif

    ! count the number of real and resiprocal part nodes
    !
    if (main_rank) then
      num_realnodes  = 0
      num_recipnodes = 0
      do i = 1, nproc_world
        if (node_real(i))  num_realnodes  = num_realnodes  + 1
        if (node_recip(i)) num_recipnodes = num_recipnodes + 1
      end do

    ! calculate total time over all nodes
    !
      sumtime(1:NumTimers) = 0.0_dp
      do i = 1, NumTimers
        maxtime(i) = tottime(i,1)
        mintime(i) = tottime(i,1)
        do j = 1, nproc_world
          sumtime(i) = sumtime(i) + tottime(i,j)
          if (maxtime(i) < tottime(i,j)) then
            maxtime(i) = tottime(i,j)
          endif
          if (mintime(i) > tottime(i,j)) then
            mintime(i) = tottime(i,j)
          endif
        end do
      end do

      ! average over all ranks
      !
      avetime(TimerTotal)      = sumtime(TimerTotal)      / nproc_world
      avetime(TimerEnergy)     = sumtime(TimerEnergy)     / nproc_world
      avetime(TimerBond)       = sumtime(TimerBond)       / num_realnodes
      avetime(TimerAngle)      = sumtime(TimerAngle)      / num_realnodes
      avetime(TimerDihedral)   = sumtime(TimerDihedral)   / num_realnodes
      avetime(TimerRestraint)  = sumtime(TimerRestraint)  / num_realnodes
      avetime(TimerNonBond)    = sumtime(TimerNonBond)    / nproc_world
      avetime(TimerPmeReal)    = sumtime(TimerPmeReal)    / num_realnodes
      if (num_recipnodes == 0) then
        avetime(TimerPmeRecip) = 0.0_dp
      else
        avetime(TimerPmeRecip) = sumtime(TimerPmeRecip)   / num_recipnodes
      end if
      avetime(TimerPairList)   = sumtime(TimerPairList)   / num_realnodes
      avetime(TimerDynamics)   = sumtime(TimerDynamics)   / nproc_world
      avetime(TimerConstraint) = sumtime(TimerConstraint) / nproc_world
      avetime(TimerUpdate)     = sumtime(TimerUpdate)     / nproc_world
      avetime(TimerComm1)      = sumtime(TimerComm1)      / nproc_world
      avetime(TimerComm2)      = sumtime(TimerComm2)      / nproc_world
      avetime(TimerComm3)      = sumtime(TimerComm3)      / nproc_world
      avetime(TimerIntegrator) = sumtime(TimerIntegrator) / nproc_world
      avetime(TimerSolvation)  = sumtime(TimerSolvation)  / nproc_world
      avetime(TimerGB)         = sumtime(TimerGB)         / nproc_world
      avetime(TimerSA)         = sumtime(TimerSA)         / nproc_world
    
      avetime(TimerQMMM ) = sumtime(TimerQMMM ) / nproc_world
    
      avetime(TimerTest1) = sumtime(TimerTest1) / nproc_world
      avetime(TimerTest2) = sumtime(TimerTest2) / nproc_world
      avetime(TimerTest3) = sumtime(TimerTest3) / nproc_world
      avetime(TimerTest4) = sumtime(TimerTest4) / nproc_world
      avetime(TimerTest5) = sumtime(TimerTest5) / nproc_world
      avetime(TimerTest6) = sumtime(TimerTest6) / nproc_world
      avetime(TimerTest7) = sumtime(TimerTest7) / nproc_world
      avetime(TimerTest8) = sumtime(TimerTest8) / nproc_world
      avetime(TimerTest9) = sumtime(TimerTest9) / nproc_world
      avetime(TimerTest10) = sumtime(TimerTest10) / nproc_world

      ! write detailed timer profile of the main_rank
      !
!    if (main_rank) then
      avetime_setup   = avetime(TimerTotal) - avetime(TimerDynamics)

      write(MsgOut,'(a)') 'Output_Time> Averaged timer profile (Min, Max)'
      write(MsgOut,'(a,f12.3)') '  total time      =', avetime(TimerTotal)
      write(MsgOut,'(a,f12.3)') '    setup         =', avetime_setup
      write(MsgOut,'(a,f12.3)') '    dynamics      =', avetime(TimerDynamics)
      write(MsgOut,'(a,f12.3)') '      energy      =', avetime(TimerEnergy)
      write(MsgOut,'(a,f12.3)') '      integrator  =', avetime(TimerIntegrator)
      write(MsgOut,'(a,f12.3,a,f12.3,a,f12.3,a)')                              &
                                '      pairlist    =',avetime(TimerPairList),  &
                                ' (',mintime(TimerPairList),',',               &
                                 maxtime(TimerPairList),')'
      write(MsgOut,'(a)')       '  energy           '
      write(MsgOut,'(a,f12.3,a,f12.3,a,f12.3,a)')                              &
                                '    bond          =',avetime(TimerBond),      &
                                ' (',mintime(TimerBond),',',                   &
                                 maxtime(TimerBond),')'
      write(MsgOut,'(a,f12.3,a,f12.3,a,f12.3,a)')                              &
                                '    angle         =',avetime(TimerAngle),     &
                                ' (',mintime(TimerAngle),',',                  &
                                 maxtime(TimerAngle),')'
      write(MsgOut,'(a,f12.3,a,f12.3,a,f12.3,a)')                              &
                                '    dihedral      =',avetime(TimerDihedral),  &
                                ' (',mintime(TimerDihedral),',',               &
                                 maxtime(TimerDihedral),')'
      write(MsgOut,'(a,f12.3,a,f12.3,a,f12.3,a)')                              &
                                '    nonbond       =',avetime(TimerNonBond),   &
                                ' (',mintime(TimerNonBond),',',                &
                                 maxtime(TimerNonBond),')'
      write(MsgOut,'(a,f12.3,a,f12.3,a,f12.3,a)')                              &
                                '      pme real    =',avetime(TimerPmeReal),   &
                                ' (',mintime(TimerPmeReal),',',                &
                                 maxtime(TimerPmeReal),')'
      write(MsgOut,'(a,f12.3,a,f12.3,a,f12.3,a)')                              &
                                '      pme recip   =',avetime(TimerPmeRecip),  &
                                ' (',mintime(TimerPmeRecip),',',               &
                                 maxtime(TimerPmeRecip),')'
      write(MsgOut,'(a,f12.3,a,f12.3,a,f12.3,a)')                              &
                                '    solvation     =',avetime(TimerSolvation), &
                                ' (',mintime(TimerSolvation),',',              &
                                 maxtime(TimerSolvation),')'
      write(MsgOut,'(a,f12.3,a,f12.3,a,f12.3,a)')                              &
                                '      polar       =',avetime(TimerGB),        &
                                ' (',mintime(TimerGB),',',                     &
                                 maxtime(TimerGB),')'
      write(MsgOut,'(a,f12.3,a,f12.3,a,f12.3,a)')                              &
                                '      non-polar   =',avetime(TimerSA),        &
                                ' (',mintime(TimerSA),',',                     &
                                 maxtime(TimerSA),')'
!     write(MsgOut,'(a,f12.3,a,f12.3,a,f12.3,a)')                              &
!                               '      cpu_real    =',avetime(TimerTest10),    &
!                               ' (',mintime(TimerTest10),',',                 &
!                                maxtime(TimerTest10),')'
      write(MsgOut,'(a,f12.3,a,f12.3,a,f12.3,a)')                              &
                                '    restraint     =',avetime(TimerRestraint), &
                                ' (',mintime(TimerRestraint),',',              &
                                 maxtime(TimerRestraint),')'
      write(MsgOut,'(a,f12.3,a,f12.3,a,f12.3,a)')                              &
                                '    qmmm          =',avetime(TimerQMMM),      &
                                ' (',mintime(TimerQMMM),',',                   &
                                 maxtime(TimerQMMM),')'
      write(MsgOut,'(a)')       '  integrator       '
      write(MsgOut,'(a,f12.3,a,f12.3,a,f12.3,a)')                              &
                                '    constraint    =',                         &
                                avetime(TimerConstraint),                      &
                                ' (',mintime(TimerConstraint),',',             &
                                 maxtime(TimerConstraint),')'
      write(MsgOut,'(a,f12.3,a,f12.3,a,f12.3,a)')                              &
                                '    update        =',                         &
                                avetime(TimerUpdate),                          &
                                ' (',mintime(TimerUpdate),',',                 &
                                 maxtime(TimerUpdate),')'
      write(MsgOut,'(a,f12.3,a,f12.3,a,f12.3,a)')                              &
                                '    comm_coord    =',                         &
                                avetime(TimerComm1),                           &
                                ' (',mintime(TimerComm1),',',                  &
                                 maxtime(TimerComm1),')'
      write(MsgOut,'(a,f12.3,a,f12.3,a,f12.3,a)')                              &
                                '    comm_force    =',                         &
                                avetime(TimerComm2),                           &
                                ' (',mintime(TimerComm2),',',                  &
                                 maxtime(TimerComm2),')'
      write(MsgOut,'(a,f12.3,a,f12.3,a,f12.3,a)')                              &
                                '    comm_migrate  =',                         &
                                avetime(TimerComm3),                           &
                                ' (',mintime(TimerComm3),',',                  &
                                 maxtime(TimerComm3),')'
      write(MsgOut,'(a)') ''

    end if


    ! deallocate local array
    !
    deallocate(node_real, node_recip, tottime, sumtime, avetime, &
               mintime, maxtime, stat=dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine output_time

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      get_unix_time
  !> @brief        get the local unix time (not GMT)
  !! @authors      JJ, ST(FJ), NT
  !! @return       milliseconds
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_unix_time()

    ! return value
    real(dp)               :: get_unix_time

    ! parameters
    real(dp)               :: usec

#ifdef KCOMP
    call gettod(usec)
#elif __INTEL_COMPILER
    call clockx(usec)
#else
    integer,   save :: init_flg = 0
    integer*8 :: CountPerSec, CountMax
    real(dp), save :: ratio_usec
    integer*8 time_c

    if( init_flg == 0 ) then
      init_flg = 1
      call system_clock(time_c, CountPerSec, CountMax)
      ratio_usec = 1.0d+6 / dble(CountPerSec)
    else
      call system_clock(time_c)
    endif

    usec = dble(time_c)*ratio_usec
#endif

     get_unix_time = usec

    return

  end function get_unix_time

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      days
  !> @brief        get number of days from AD.1
  !! @authors      NT
  !! @param[in]    year : A.D. year 
  !! @return       number of days
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function days(year)

    ! return values
    integer(8)               :: days

    ! formal arguments
    integer(8),              intent(in)    :: year


    days = year * 365
    days = days + year / 4
    days = days - year / 100
    days = days + year / 400

    return

  end function days

end module timers_mod
