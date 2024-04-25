!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_pairlist_mod
!> @brief   set pairlist for nonbonded interactions
!! @authors Jaewoon Jung (JJ), Takaharu Mori (TM), Takashi Imai (TI), 
!!          Motoshi Kamiya (MK), Yuji Sugita (YS), Kiyoshi Yagi (KY)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_pairlist_mod

  use at_boundary_str_mod
  use at_pairlist_str_mod
  use at_enefunc_str_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  real(wp),         parameter :: TableWaterMargin = 2.0_wp
  real(wp),         parameter :: FactNumNb15   = 1.5_wp
  real(wp),         parameter :: FactThreshold = 0.95_wp

  ! subroutines
  public  :: setup_pairlist
  public  :: update_pairlist
  private :: update_pairlist_nobc
  private :: update_pairlist_nobc_gbsa
  private :: update_ecqm15_nonb
!  private :: update_pairlist_pbc
!  private :: update_pairlist_pbc_cell
  private :: update_pairlist_solute_solute
  private :: update_pairlist_solute_water
  private :: update_pairlist_water_water
  private :: check_pairlist_memory_size

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_pairlist
  !> @brief        initialize/allocate/setup pairlist for nonbonded interactions
  !! @authors      YS, TI, JJ, TM
  !! @param[in]    enefunc  : information of potential energy function
  !! @param[in]    boundary : information of boundary condition
  !! @param[in]    coord    : atomic coordinate
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_pairlist(enefunc, boundary, coord, pairlist)

    ! formal arguments
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_boundary),        intent(in)    :: boundary
    real(wp),                intent(in)    :: coord(:,:)
    type(s_pairlist),        intent(inout) :: pairlist
    
    ! local variables
    integer                  :: natom, nthread
#ifdef OMP
    integer                  :: omp_get_max_threads
#endif


    ! do not make pairlist, if my_node does not compute 
    !   non-bonded energy in the real part
    !
    if (.not. real_calc) &
      return

    ! initialize pairlist
    !
    call init_pairlist(pairlist)

    pairlist%pairlistdist = enefunc%pairlistdist

    ! allocate pairlist
    !
#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif
    call alloc_pairlist(pairlist, PairListNthreads, nthread)
    if (enefunc%gbsa_use) then
      call alloc_pairlist(pairlist, PairListNthreadsGbsa, nthread)
    end if

    natom = size(coord(1,:))

    select case (boundary%type)

    case (BoundaryTypeNOBC)

      call alloc_pairlist(pairlist, PairListAtomNobc, natom)
      if (enefunc%gbsa_use) then
        call alloc_pairlist(pairlist, PairListAtomNobcGbsa, natom)
      end if

    case (BoundaryTypePBC)

      call alloc_pairlist(pairlist, PairListPbcSolute, enefunc%table%num_solute)
      call alloc_pairlist(pairlist, PairListPbcWater,  enefunc%table%num_water)

    end select

    ! make pairlist
    !
    pairlist%allocation = .true.
    call update_pairlist(enefunc, boundary, coord, pairlist)
    pairlist%allocation = .false.

    return

  end subroutine setup_pairlist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist
  !> @brief        update pairlist for nonbonded interactions
  !! @authors      YS, TM, TI, JJ
  !! @param[in]    enefunc  : information of potential energy function
  !! @param[in]    boundary : information of boundary condition
  !! @param[in]    coord    : atomic coordinate
  !! @param[inout] pairlist : information of pairlist
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist(enefunc, boundary, coord, pairlist)

    ! formal arguments
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_boundary),        intent(in)    :: boundary
    real(wp),                intent(in)    :: coord(:,:)
    type(s_pairlist),        intent(inout) :: pairlist


    call timer(TimerPairList, TimerOn)

    if (pairlist%allocation) then
      pairlist%allocate_nobc   = .true.
      pairlist%allocate_pbc    = .true.
      pairlist%allocate_solsol = .true.
      pairlist%allocate_solwat = .true.
      pairlist%allocate_watwat = .true.
    end if

    select case (boundary%type)

    case (BoundaryTypeNOBC)

      if (enefunc%gbsa_use) then
        call update_pairlist_nobc_gbsa(enefunc, coord, pairlist)
      else
        call update_pairlist_nobc     (enefunc, coord, pairlist)
      end if

      if(enefunc%qmmm%do_qmmm .and. enefunc%qmmm%num_qmmmbonds > 0)  &
           call update_ecqm15_nonb(enefunc, coord, pairlist)

    case (BoundaryTypePBC)

!      if (enefunc%table%table) then
        call update_pairlist_solute_solute(enefunc, boundary, coord, pairlist)
        call update_pairlist_solute_water (enefunc, boundary, coord, pairlist)
        call update_pairlist_water_water  (enefunc, boundary, coord, pairlist)
!      else
!        if (boundary%use_cell_linked_list) then
!          call update_pairlist_pbc_cell(enefunc, boundary, coord, pairlist)
!        else
!          call update_pairlist_pbc(enefunc, boundary, coord, pairlist)
!        end if
!      end if

    end select

    call timer(TimerPairList, TimerOff)

    return

  end subroutine update_pairlist

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_nobc
  !> @brief        update pairlist when no boundary condition is applied
  !! @authors      YS, JJ, TM
  !! @param[in]    enefunc  : information of potential energy function
  !! @param[in]    coord    : atomic coordinate
  !! @param[inout] pairlist : information of pairlist
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_nobc(enefunc, coord, pairlist)
  
    ! formal arguments
    type(s_enefunc),          intent(in)    :: enefunc
    real(wp),                 intent(in)    :: coord(:,:)
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                 :: pairdist2, dij(1:3), rij2
    integer                  :: i, j, k, natom, n, nloops
    integer                  :: num_excl, ini_excl, fin_excl
    integer                  :: num_nb14, ini_nb14, fin_nb14
    integer                  :: num_nb15_max
    integer                  :: id, my_id, nthread
#ifdef OMP
    integer                  :: omp_get_thread_num, omp_get_max_threads
#endif
    logical                  :: nb15_calc
    logical                  :: proceed, do_allocate

    integer,     pointer     :: num_nb15_pre(:), num_nb15(:)


    num_nb15_pre => pairlist%num_nb15_pre
    num_nb15     => pairlist%num_nb15
    pairdist2    =  pairlist%pairlistdist * pairlist%pairlistdist
    natom        =  size(coord(1,:))

    ! initialize
    !
#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif
    pairlist%num_nb15_calc(1:natom,1:nthread) = 0

    ! if pairlist%allocation = .true. , allocation + make pairlist
    ! if pairlist%allocation = .false., make pairlist only
    !
    if (pairlist%allocate_nobc) then
      nloops = 2
      do_allocate = .true.
    else
      nloops = 1
      do_allocate = .false.
    end if
    
    do n = 1, nloops

      num_excl = 0
      num_nb14 = 0
      num_nb15(1:nthread) = 0

      if (.not. do_allocate) &
        num_nb15_pre(1:nthread) = 0

      !$omp parallel                                                       &
      !$omp private(id, my_id, i, ini_excl, fin_excl, ini_nb14, fin_nb14,  &
      !$omp         proceed, j, nb15_calc, k, dij, rij2)                   &
      !$omp firstprivate(num_excl, num_nb14, do_allocate)
      !
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      my_id = my_city_rank * nthread + id
      id = id + 1

      do i = 1, natom-1

        ini_excl = num_excl + 1
        fin_excl = num_excl + enefunc%num_nonb_excl(i)
        num_excl = fin_excl

        proceed = .true.
        if (mod(i-1,nproc_city*nthread) /= my_id) proceed = .false.

        if (proceed) then
          do j = i+1, natom

            ! compute distance
            !
            dij(1:3) = coord(1:3,i) - coord(1:3,j)
            rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            ! store interaction table
            !
            if (rij2 < pairdist2) then

              ! 1-2 and 1-3 interactions are removed
              !
              nb15_calc = .true.
              do k = ini_excl, fin_excl
                if (j == enefunc%nonb_excl_list(k)) then
                  nb15_calc = .false.
                  exit
                end if
              end do

              ! 1-4 interactions are removed
              !
              if (enefunc%forcefield/=ForcefieldKBGO) then
                do k = 1, enefunc%num_nb14_calc(i)
                  if (j == enefunc%nb14_calc_list(k,i)) then
                    nb15_calc = .false.
                    exit
                  end if
                end do
              endif
              if (.not. nb15_calc) cycle

              num_nb15(id) = num_nb15(id) + 1

              if (.not. do_allocate) &
                pairlist%nb15_calc_list(num_nb15(id),id) = j
            end if

          end do
        end if

        if (.not. do_allocate) then
          pairlist%num_nb15_calc(i,id) = num_nb15(id) - num_nb15_pre(id)
          num_nb15_pre(id)             = num_nb15(id)
        end if

      end do
      !$omp end parallel

      ! allocate memory of pairlist
      !
      if (do_allocate) then
        num_nb15_max = 0
        do i = 1, nthread
          num_nb15_max = max(num_nb15_max, num_nb15(i))
        end do

        num_nb15_max = int(real(num_nb15_max,wp)*FactNumNb15)
        call alloc_pairlist(pairlist, PairListIntNobc, num_nb15_max)
        pairlist%num_nb15_max = num_nb15_max

        do_allocate = .false.
      end if

    end do

    ! check memory size
    !
    call check_pairlist_memory_size(num_nb15, pairlist%num_nb15_max, &
                                    pairlist%allocate_nobc)

    return

  end subroutine update_pairlist_nobc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_nobc_gbsa
  !> @brief        update pairlist under NOBC for GBSA
  !! @authors      TM
  !! @param[in]    enefunc  : information of potential energy function
  !! @param[in]    coord    : atomic coordinate
  !! @param[inout] pairlist : information of pairlist
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_nobc_gbsa(enefunc, coord, pairlist)

    ! formal arguments
    type(s_enefunc),          intent(in)    :: enefunc
    real(wp),                 intent(in)    :: coord(:,:)
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                 :: pairdist2, dij(1:3), rij2
    integer                  :: i, j, k, natom, n, nloops
    integer                  :: num_excl, ini_excl, fin_excl
    integer                  :: num_nb14, ini_nb14, fin_nb14
    integer                  :: num_nb15_max
    integer                  :: num_all_max
    integer                  :: id, my_id, nthread
#ifdef OMP
    integer                  :: omp_get_thread_num, omp_get_max_threads
#endif
    logical                  :: nb15_calc
    logical                  :: proceed, do_allocate

    integer,     pointer     :: num_nb15_pre(:), num_nb15(:)
    integer,     pointer     :: num_all_pre(:), num_all(:)

    num_nb15_pre => pairlist%num_nb15_pre
    num_nb15     => pairlist%num_nb15
    num_all_pre  => pairlist%num_all_pre
    num_all      => pairlist%num_all
    pairdist2    =  pairlist%pairlistdist * pairlist%pairlistdist
    natom        =  size(coord(1,:))

    ! initialize
    !
#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif
    pairlist%num_nb15_calc(1:natom,1:nthread) = 0
    pairlist%num_all_calc (1:natom,1:nthread) = 0

    ! if pairlist%allocation = .true. , allocation + make pairlist
    ! if pairlist%allocation = .false., make pairlist only
    !
    if (pairlist%allocate_nobc) then
      nloops = 2
      do_allocate = .true.
    else
      nloops = 1
      do_allocate = .false.
    end if

    do n = 1, nloops

      num_excl = 0
      num_nb14 = 0
      num_nb15(1:nthread) = 0
      num_all (1:nthread) = 0

      if (.not. do_allocate) then
        num_nb15_pre(1:nthread) = 0
        num_all_pre (1:nthread) = 0
      end if

      !$omp parallel                                                       &
      !$omp private(id, my_id, i, ini_excl, fin_excl, ini_nb14, fin_nb14,  &
      !$omp         proceed, j, nb15_calc, k, dij, rij2)                   &
      !$omp firstprivate(num_excl, num_nb14, do_allocate)
      !
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      my_id = my_city_rank * nthread + id
      id = id + 1

      do i = 1, natom-1

        ini_excl = num_excl + 1
        fin_excl = num_excl + enefunc%num_nonb_excl(i)
        num_excl = fin_excl

        proceed = .true.
        if (mod(i-1,nproc_city*nthread) /= my_id) proceed = .false.

        if (proceed) then
          do j = i+1, natom

            ! compute distance
            !
            dij(1:3) = coord(1:3,i) - coord(1:3,j)
            rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            ! store interaction table
            !
            if (rij2 < pairdist2) then

              ! for GB
              !
              num_all(id) = num_all(id) + 1

              if (.not. do_allocate) &
                pairlist%all_calc_list(num_all(id),id) = j

              ! 1-2 and 1-3 interactions are removed
              !
              nb15_calc = .true.
              do k = ini_excl, fin_excl
                if (j == enefunc%nonb_excl_list(k)) then
                  nb15_calc = .false.
                  exit
                end if
              end do

              ! 1-4 interactions are removed
              !
              do k = 1, enefunc%num_nb14_calc(i)
                if (j == enefunc%nb14_calc_list(k,i)) then
                  nb15_calc = .false.
                  exit
                end if
              end do
              if (.not. nb15_calc) cycle

              num_nb15(id) = num_nb15(id) + 1

              if (.not. do_allocate) &
                pairlist%nb15_calc_list(num_nb15(id),id) = j
            end if

          end do
        end if

        if (.not. do_allocate) then
          pairlist%num_nb15_calc(i,id) = num_nb15(id) - num_nb15_pre(id)
          num_nb15_pre(id)             = num_nb15(id)

          pairlist%num_all_calc(i,id) = num_all(id) - num_all_pre(id)
          num_all_pre(id)             = num_all(id)
        end if

      end do
      !$omp end parallel

      ! allocate memory of pairlist
      !
      if (do_allocate) then

        num_nb15_max = 0
        do i = 1, nthread
          num_nb15_max = max(num_nb15_max, num_nb15(i))
        end do
        num_nb15_max = int(real(num_nb15_max,wp)*FactNumNb15)
        call alloc_pairlist(pairlist, PairListIntNobc, num_nb15_max)
        pairlist%num_nb15_max = num_nb15_max

        num_all_max = 0
        do i = 1, nthread
          num_all_max = max(num_all_max, num_all(i))
        end do
        num_all_max = int(real(num_all_max,wp)*FactNumNb15)
        call alloc_pairlist(pairlist, PairListIntNobcGbsa, num_all_max)
        pairlist%num_all_max = num_all_max

        do_allocate = .false.
      end if

    end do

    ! check memory size
    !
    call check_pairlist_memory_size(num_nb15, pairlist%num_nb15_max, &
                                    pairlist%allocate_nobc)

    call check_pairlist_memory_size(num_all, pairlist%num_all_max, &
                                    pairlist%allocate_nobc)

    return

  end subroutine update_pairlist_nobc_gbsa

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_ecqm15_nonb
  !> @brief        create ecqm_nb15_list
  !! @authors      KY
  !! @param[in] enefunc  : potential energy functions information
  !! @param[in] coord    : atomic coordinate
  !! @param[in] pairlist : information of pairlist
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_ecqm15_nonb(enefunc, coord, pairlist)

    ! formal arguments
    type(s_enefunc),  target, intent(in)   :: enefunc
    real(wp),                 intent(in)   :: coord(:,:)
    type(s_pairlist), target, intent(inout):: pairlist

    ! local variables
    real(wp)                 :: pairdist2, dij(1:3), rij2
    integer                  :: i, j, k, ii, jj
    logical                  :: nb15_calc

    integer,     pointer     :: nonb_excl_list(:,:), nb14_calc_list(:,:)
    integer                  :: qm_natoms, ec_natoms, iqm, iec
    integer,     pointer     :: qmatom_id(:), ecatom_id(:)
    integer, allocatable     :: num_nb15(:), nb15_list(:,:)
    integer                  :: kalloc_stat, kdealloc_stat
    integer                  :: n_init, n_final

    integer                  :: id, nthread
#ifdef OMP
    integer                  :: omp_get_thread_num, omp_get_max_threads
#endif

    !dbg write(MsgOut,'("QM-EC_15 pair")')

    qm_natoms =  enefunc%qmmm%qm_natoms
    qmatom_id => enefunc%qmmm%qmatom_id
    ec_natoms =  enefunc%qmmm%ec_natoms
    ecatom_id => enefunc%qmmm%ecatom_id

    nonb_excl_list => enefunc%table%nonb_excl_list
    nb14_calc_list => enefunc%nb14_calc_list
    pairdist2      =  pairlist%pairlistdist * pairlist%pairlistdist

    if (.not. allocated(pairlist%ecqm_nb15_list)) then
      call alloc_pairlist(pairlist, PairListEcqm, qm_natoms*ec_natoms)
    end if

#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif

    !dbg write(MsgOut,'("qm_natom, ec_natom",2i4)') qm_natoms, ec_natoms

    allocate(num_nb15(nthread), stat=kalloc_stat)
    if (kalloc_stat /=0) call error_msg_alloc
    allocate(nb15_list(2, qm_natoms*ec_natoms), stat=kalloc_stat)
    if (kalloc_stat /=0) call error_msg_alloc

    !$omp parallel                                             &
    !$omp private(id, iqm, iec, i, j, k, dij, rij2, nb15_calc, &
    !$omp         nb15_list, n_init, n_final)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    id = id + 1

    num_nb15(id) = 0
    do iqm = 1, qm_natoms
      if (mod(iqm-1,nthread) /= id-1) cycle
      !dbg write(MsgOut,'("id, nthread, mod(iqm-1,nthread)",3i4)') id, nthread, mod(iqm-1,nthread)

      do iec = 1, ec_natoms
        if (qmatom_id(iqm) < ecatom_id(iec)) then
          i = qmatom_id(iqm)
          j = ecatom_id(iec)
        else
          i = ecatom_id(iec)
          j = qmatom_id(iqm)
        end if

        ! compute distance
        !
        dij(1:3) = coord(1:3,i) - coord(1:3,j)
        rij2     = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

        ! store interaction table
        !
        if (rij2 < pairdist2) then

          ! 1-2 and 1-3 interactions are removed
          !
          nb15_calc = .true.
          do k = 1, enefunc%num_nonb_excl(i)
            if (j .eq. nonb_excl_list(k,i)) then
              nb15_calc = .false.
              exit
            end if
          end do
          if (.not. nb15_calc) cycle

          ! 1-4 interactions are removed
          !
          if (enefunc%forcefield/=ForcefieldKBGO) then
            do k = 1, enefunc%num_nb14_calc(i)
              if (j == nb14_calc_list(k,i)) then
                nb15_calc = .false.
                exit
              end if
            end do
          endif
          if (.not. nb15_calc) cycle

          ! NOTE: i < j
          num_nb15(id) = num_nb15(id) + 1
          nb15_list(1,num_nb15(id)) = i
          nb15_list(2,num_nb15(id)) = j
        end if

      end do
    end do

    !$omp barrier

    !dbg write(MsgOut,'("id,num_nb15",5i4)') id, num_nb15

    n_init  = 1
    do i = 1, id-1
       n_init = n_init + num_nb15(i)
    end do
    n_final = n_init + num_nb15(id) - 1
    !dbg write(MsgOut,'("id,n_init,n_final",3i4)') id, n_init, n_final

    if (id == nthread)  pairlist%ecqm_num_nb15 = n_final 
    if (num_nb15(id) > 0) then
      pairlist%ecqm_nb15_list(:,n_init:n_final) = nb15_list(:,1:num_nb15(id))
    end if

    !$omp end parallel

    deallocate(nb15_list, num_nb15, stat=kdealloc_stat)
    if (kdealloc_stat /=0) call error_msg_dealloc

    !dbg write(MsgOut,'("ecqm_num_nb15",i4)') pairlist%ecqm_num_nb15
    !dbg do i = 1, pairlist%ecqm_num_nb15
    !dbg   write(MsgOut,'(2i4)') pairlist%ecqm_nb15_list(:,i)
    !dbg end do

  end subroutine update_ecqm15_nonb

! !======1=========2=========3=========4=========5=========6=========7=========8
! !
! !  Subroutine    update_pairlist_pbc
! !> @brief        update pairlist under PBC without cell-linked list 
! !! @authors      TM, JJ
! !! @param[in]    enefunc  : information of potential energy function
! !! @param[in]    boundary : information of boundary
! !! @param[in]    coord    : atomic coordinate
! !! @param[inout] pairlist : information of pairlist
! !! @note         pairlist%trans_vector is not considered
! !
! !======1=========2=========3=========4=========5=========6=========7=========8

! subroutine update_pairlist_pbc(enefunc, boundary, coord, pairlist)
! 
!   ! formal arguments
!   type(s_enefunc),          intent(in)    :: enefunc
!   type(s_boundary), target, intent(in)    :: boundary
!   real(wp),                 intent(in)    :: coord(:,:)
!   type(s_pairlist), target, intent(inout) :: pairlist

!   ! local variables
!   real(wp)                  :: pairdist2, dij(1:3), rij2
!   real(wp)                  :: dij_pbc(1:3)
!   real(wp)                  :: box_size(1:3)
!   integer                   :: i, j, k, n, natom, nloops
!   integer                   :: num_excl, ini_excl, fin_excl
!   integer                   :: num_nb14, ini_nb14, fin_nb14
!   integer                   :: num_nb15_max
!   integer                   :: id, my_id, nthread
!#ifdef OMP
!   integer                   :: omp_get_thread_num, omp_get_max_threads
!#endif
!   logical                   :: nb15_calc
!   logical                   :: proceed, do_allocate

!   integer,      pointer     :: num_nb15_pre(:), num_nb15(:)


!   num_nb15_pre => pairlist%num_nb15_pre
!   num_nb15     => pairlist%num_nb15
!   box_size(1)  =  boundary%box_size_x
!   box_size(2)  =  boundary%box_size_y
!   box_size(3)  =  boundary%box_size_z
!   pairdist2    =  pairlist%pairlistdist * pairlist%pairlistdist
!   natom        =  size(coord(1,:))

!   ! initialize
!   !
!#ifdef OMP
!   nthread = omp_get_max_threads()
!#else
!   nthread = 1
!#endif
!   pairlist%num_nb15_calc(1:natom,1:nthread) = 0


!   ! if pairlist%allocation = .true. , allocation + make pairlist
!   ! if pairlist%allocation = .false., make pairlist only
!   !
!   if (pairlist%allocate_pbc) then
!     nloops = 2
!     do_allocate = .true.
!   else
!     nloops = 1
!     do_allocate = .false.
!   end if

!   do n = 1, nloops

!     num_excl = 0
!     num_nb14 = 0
!     num_nb15(1:nthread) = 0

!     if (.not. do_allocate) &
!       num_nb15_pre(1:nthread) = 0

!     !$omp parallel                                                       &
!     !$omp private(id, my_id, i, ini_excl, fin_excl, ini_nb14, fin_nb14,  &
!     !$omp         proceed, j, nb15_calc, k, dij, dij_pbc, rij2) &
!     !$omp firstprivate(num_excl, num_nb14, do_allocate)
!     !
!#ifdef OMP
!     id = omp_get_thread_num()
!#else
!     id = 0
!#endif
!     my_id = my_city_rank * nthread + id
!     id    = id + 1

!     do i = 1, natom-1

!       ini_excl = num_excl + 1
!       fin_excl = num_excl + enefunc%num_nonb_excl(i)
!       num_excl = fin_excl

!       ini_nb14 = num_nb14 + 1
!       fin_nb14 = num_nb14 + enefunc%num_nb14_calc(i)
!       num_nb14 = fin_nb14

!       proceed = .true.
!       if (mod(i-1,nproc_city*nthread) /= my_id) proceed = .false.

!       if (proceed) then

!         do j = i+1, natom

!           ! compute distance
!           !
!           dij(1:3) = coord(1:3,i) - coord(1:3,j)

!           ! consider the periodic boundary
!           !
!           dij_pbc(1:3) = box_size(1:3)*anint(dij(1:3)/box_size(1:3))
!           dij(1:3)  = dij(1:3) - dij_pbc(1:3)
!           rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

!           ! store interaction table
!           !
!           if (rij2 < pairdist2) then

!             ! 1-2 and 1-3 interactions are removed
!             !
!             nb15_calc = .true.
!             do k = ini_excl, fin_excl
!               if (j == enefunc%nonb_excl_list(k)) then
!                 nb15_calc = .false.
!                 exit
!               end if
!             end do

!             ! 1-4 interactions are removed
!             !
!             do k = ini_nb14, fin_nb14
!               if (j == enefunc%nb14_calc_list(k)) then
!                 nb15_calc = .false.
!                 exit
!               end if
!             end do
!             if (.not. nb15_calc) cycle

!             num_nb15(id) = num_nb15(id) + 1

!             if (.not. do_allocate) &
!               pairlist%nb15_calc_list(num_nb15(id),id) = j
!           end if

!         end do
!       end if

!       if (.not. do_allocate) then
!         pairlist%num_nb15_calc(i,id) = num_nb15(id) - num_nb15_pre(id)
!         num_nb15_pre(id)             = num_nb15(id)
!       end if

!     end do
!     !$omp end parallel

!     ! allocate memory of pairlist
!     !
!     if (do_allocate) then
!       num_nb15_max = 0
!       do i = 1, nthread
!         num_nb15_max = max(num_nb15_max, num_nb15(i))
!       end do

!       num_nb15_max = int(real(num_nb15_max,wp)*FactNumNb15)
!       call alloc_pairlist(pairlist, PairListIntPbc, num_nb15_max)
!       pairlist%num_nb15_max = num_nb15_max

!       do_allocate = .false.
!     end if

!   end do

!   ! check memory size
!   !
!   call check_pairlist_memory_size(num_nb15, pairlist%num_nb15_max, &
!                                   pairlist%allocate_pbc)

!   return

! end subroutine update_pairlist_pbc

! !======1=========2=========3=========4=========5=========6=========7=========8
! !
! !  Subroutine    update_pairlist_pbc_cell
! !> @brief        update pairlist when periodic boundary condition is applied
! !! @authors      TI, JJ, TM
! !! @param[in]    enefunc  : information of potential energy function
! !! @param[in]    boundary : information of boundary condition
! !! @param[in]    coord    : atomic coordinate
! !! @param[inout] pairlist : information of pairlist
! !! @note         pairlist%trans_vector is not used
! !
! !======1=========2=========3=========4=========5=========6=========7=========8

! subroutine update_pairlist_pbc_cell(enefunc, boundary, coord, pairlist)

!   ! formal arguments
!   type(s_enefunc),          intent(in)    :: enefunc
!   type(s_boundary), target, intent(in)    :: boundary
!   real(wp),                 intent(in)    :: coord(:,:)
!   type(s_pairlist), target, intent(inout) :: pairlist

!   ! local variables
!   real(wp)                  :: dij(1:3), rij2, pairdist2
!   real(wp)                  :: origin(1:3), dij_pbc(1:3)
!   real(wp)                  :: r_shift(1:3)
!   real(wp)                  :: csize_inv(1:3)
!   real(wp)                  :: box_size(1:3), half_box_size(1:3)
!   integer                   :: i, j, k, n
!   integer                   :: natom, nloops
!   integer                   :: num_excl, ini_excl, fin_excl
!   integer                   :: num_nb14, ini_nb14, fin_nb14
!   integer                   :: ic(1:3), icel, inbc
!   integer                   :: num_nb15_max
!   integer                   :: id, my_id, nthread
!#ifdef OMP
!   integer                   :: omp_get_max_threads, omp_get_thread_num
!#endif
!   logical                   :: nb15_calc
!   logical                   :: proceed, do_allocate
!   
!   real(wp),         pointer :: bsize_x, bsize_y, bsize_z
!   integer,          pointer :: list(:), ic_atm(:), head(:)
!   integer,          pointer :: ncel_x, ncel_y, ncel_z, ncel
!   integer,          pointer :: num_nb15_pre(:), num_nb15(:)


!   ! use pointers
!   !
!   natom = size(coord(1,:))

!   num_nb15_pre => pairlist%num_nb15_pre
!   num_nb15     => pairlist%num_nb15
!   list         => pairlist%cell_linked_list
!   ic_atm       => pairlist%atom_cell_index
!   head         => boundary%cell_head_atom
!   ncel_x       => boundary%num_cells_x
!   ncel_y       => boundary%num_cells_y
!   ncel_z       => boundary%num_cells_z
!   ncel         => boundary%num_cells

!   origin(1)    = boundary%origin_x
!   origin(2)    = boundary%origin_y
!   origin(3)    = boundary%origin_z
!   box_size(1)  = boundary%box_size_x
!   box_size(2)  = boundary%box_size_y
!   box_size(3)  = boundary%box_size_z
!   csize_inv(1) = 1.0_wp / boundary%cell_size_x
!   csize_inv(2) = 1.0_wp / boundary%cell_size_y
!   csize_inv(3) = 1.0_wp / boundary%cell_size_z
!   half_box_size(1:3) = box_size(1:3) * 0.5_wp

!   pairdist2    = pairlist%pairlistdist * pairlist%pairlistdist

!   ! initialize
!   !
!#ifdef OMP
!   nthread = omp_get_max_threads()
!#else
!   nthread = 1
!#endif
!   pairlist%num_nb15_calc(1:natom,1:nthread) = 0


!   ! make cell linked list
!   !
!   head(1:ncel) = 0
!   do i = 1, natom

!     ! coordinate shifted to the first quadrant and set into the boundary box
!     !
!     r_shift(1:3) = coord(1:3,i) - origin(1:3)
!     r_shift(1:3) = r_shift(1:3) + half_box_size(1:3) - &
!                    box_size(1:3)*anint(r_shift(1:3)/box_size(1:3))

!     ! assign which cell
!     !
!     ic(1:3) = int(r_shift(1:3)*csize_inv(1:3))
!     
!     ! upper limit of ic[x,y,z] should be ncel_[x,y,z]
!     !
!     if (ic(1) >= ncel_x) ic(1) = ncel_x-1
!     if (ic(2) >= ncel_y) ic(2) = ncel_y-1
!     if (ic(3) >= ncel_z) ic(3) = ncel_z-1
!     icel = 1 + ic(1) + ic(2)*ncel_x + ic(3)*ncel_x*ncel_y

!     ! store which cell atom i is placed
!     !
!     ic_atm(i)  = icel
!     list(i)    = head(icel)
!     head(icel) = i

!   end do


!   ! if pairlist%allocation = .true. , allocation + make pairlist
!   ! if pairlist%allocation = .false., make pairlist only
!   !
!   if (pairlist%allocate_pbc) then
!     nloops = 2
!     do_allocate = .true.
!   else
!     nloops = 1
!     do_allocate = .false.
!   end if

!   do n = 1, nloops

!     num_excl = 0
!     num_nb14 = 0
!     num_nb15(1:nthread) = 0

!     if (.not. do_allocate) &
!       num_nb15_pre(1:nthread) = 0

!     !$omp parallel                                                           &
!     !$omp private(id, my_id, i, ini_excl, fin_excl, ini_nb14, fin_nb14,      &
!     !$omp         proceed, inbc, j, nb15_calc, k, dij, dij_pbc,rij2)         &
!     !$omp firstprivate(num_excl, num_nb14, do_allocate)
!     !
!#ifdef OMP
!     id = omp_get_thread_num()
!#else
!     id = 0
!#endif
!     my_id = my_city_rank * nthread + id
!     id    = id + 1

!     do i = 1, natom - 1

!       ini_excl = num_excl + 1
!       fin_excl = num_excl + enefunc%num_nonb_excl(i)
!       num_excl = fin_excl

!       ini_nb14 = num_nb14 + 1
!       fin_nb14 = num_nb14 + enefunc%num_nb14_calc(i)
!       num_nb14 = fin_nb14

!       proceed = .true.
!       if (mod(i-1,nproc_city*nthread) /= my_id) proceed = .false.

!       if (proceed) then

!         ! serch only in the self and neighboring cells
!         !
!         do inbc = 1, 27

!           j = head(boundary%neighbor_cells(inbc,ic_atm(i)))

!           do while (j > 0)

!             if (j <= i) then
!               j = list(j)
!               cycle
!             end if

!             ! compute distance
!             !
!             dij(1:3) = coord(1:3,i) - coord(1:3,j)

!             ! consider the periodic boundary
!             !
!             dij_pbc(1:3) = box_size(1:3)*anint(dij(1:3)/box_size(1:3))
!             dij(1:3)     = dij(1:3) - dij_pbc(1:3)
!             rij2         = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

!             ! cout the number of interactions
!             !
!             if (rij2 < pairdist2) then

!               ! 1-2 and 1-3 interactions are removed
!               !
!               nb15_calc = .true.
!               do k = ini_excl, fin_excl
!                 if (j == enefunc%nonb_excl_list(k)) then
!                   nb15_calc = .false.
!                   exit
!                 end if
!               end do

!               ! 1-4 interactions are removed
!               !
!               do k = ini_nb14, fin_nb14
!                 if (j == enefunc%nb14_calc_list(k)) then
!                   nb15_calc = .false.
!                   exit
!                 end if
!               end do
!               if (.not. nb15_calc) then
!                 j = list(j)
!                 cycle
!               end if

!               num_nb15(id) = num_nb15(id) + 1

!               if (.not. do_allocate) &
!                 pairlist%nb15_calc_list(num_nb15(id),id) = j
!             end if

!             j = list(j)

!           end do
!         end do
!       end if

!       if (.not. do_allocate) then
!         pairlist%num_nb15_calc(i,id) = num_nb15(id) - num_nb15_pre(id)
!         num_nb15_pre(id)             = num_nb15(id)
!       end if

!     end do
!     !$omp end parallel

!     ! allocate memory of pairlist
!     !
!     if (do_allocate) then
!       num_nb15_max = 0
!       do i = 1, nthread
!         num_nb15_max = max(num_nb15_max, num_nb15(i))
!       end do

!       num_nb15_max = int(real(num_nb15_max,wp)*FactNumNb15)
!       call alloc_pairlist(pairlist, PairListIntPbc, num_nb15_max)
!       pairlist%num_nb15_max = num_nb15_max

!       do_allocate = .false.
!     end if

!   end do

!   ! check memory size
!   !
!   call check_pairlist_memory_size(num_nb15, pairlist%num_nb15_max, &
!                                   pairlist%allocate_pbc)

!   return

! end subroutine update_pairlist_pbc_cell

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_solute_solute
  !> @brief        update pairlist between solutes
  !! @authors      JJ, TM, MK
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : boundary conditions information
  !! @param[in]    coord    : coordinates
  !! @param[inout] pairlist : pairlist information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_solute_solute(enefunc, boundary, coord, pairlist)

    ! formal arguments
    type(s_enefunc),          intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    real(wp),                 intent(in)    :: coord(:,:)
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                  :: dij(1:3), rij2, pairdist2
    real(wp)                  :: dij_pbc(1:3), origin(1:3)
    real(wp)                  :: r_shift(1:3)
    real(wp)                  :: csize_inv(1:3)
    real(wp)                  :: box_size(1:3), half_box_size(1:3)
    integer                   :: ii, jj, i, j, k, ik, jk, n, nloops
    integer                   :: natom
    integer                   :: ic(1:3), icel, inbc
    integer                   :: num_nb15_max
    integer                   :: id, my_id, nthread
    integer                   :: ncell(1:3)
#ifdef OMP
    integer                   :: omp_get_max_threads, omp_get_thread_num
#endif
    logical                   :: nb15_calc
    logical                   :: do_allocate

    integer,          pointer :: num_nb15_pre(:), num_nb15(:), num_nb15_calc(:)
    integer,          pointer :: list(:), ic_atm(:), head(:)
    integer,          pointer :: ncel


    natom   = enefunc%table%num_solute

    num_nb15_pre => pairlist%num_nb15_pre
    num_nb15     => pairlist%num_nb15
    num_nb15_calc => pairlist%table%num_nb15_calc
    list         => pairlist%table%cell_linked_list
    ic_atm       => pairlist%table%atom_cell_index
    head         => boundary%cell_head_atom
    ncel         => boundary%num_cells

    box_size(1)  = boundary%box_size_x
    box_size(2)  = boundary%box_size_y
    box_size(3)  = boundary%box_size_z
    origin(1)    = boundary%origin_x
    origin(2)    = boundary%origin_y
    origin(3)    = boundary%origin_z
    ncell(1)     = boundary%num_cells_x
    ncell(2)     = boundary%num_cells_y
    ncell(3)     = boundary%num_cells_z
    csize_inv(1) = 1.0_wp / boundary%cell_size_x
    csize_inv(2) = 1.0_wp / boundary%cell_size_y
    csize_inv(3) = 1.0_wp / boundary%cell_size_z
    half_box_size(1:3) = box_size(1:3)*0.5_wp

    pairdist2 = pairlist%pairlistdist * pairlist%pairlistdist

    ! initialize
    !
#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif


    ! make cell linked list for solute
    !
    head(1:ncel) = 0
    do ii = 1, natom

      ! coordinate shifted to the first quadrant and set into the boundary box
      !
      i = enefunc%table%solute_list(ii)
      r_shift(1:3) = coord(1:3,i) - origin(1:3)
      r_shift(1:3) = r_shift(1:3) + half_box_size(1:3) -  &
                     box_size(1:3)*anint(r_shift(1:3)/box_size(1:3))

      ! assign which cell
      !
      ic(1:3) = int(r_shift(1:3)*csize_inv(1:3))

      ! upper limit of ic[x,y,z] should be ncel_[x,y,z]
      !
      do k = 1, 3
        if (ic(k) >= ncell(k)) ic(k) = ncell(k)-1
      end do
      icel = 1 + ic(1) + ic(2)*ncell(1) + ic(3)*ncell(1)*ncell(2)

      ! store which cell atom i is placed
      !
      ic_atm(ii) = icel
      list(ii)   = head(icel)
      head(icel) = ii

    end do

    num_nb15_calc(1:natom) = 0

    ! if pairlist%allocation = .true. , allocation + make pairlist
    ! if pairlist%allocation = .false., make pairlist only
    !
    if (pairlist%allocate_solsol) then
      nloops = 2
      do_allocate = .true.
    else
      nloops = 1
      do_allocate = .false.
    end if

    do n = 1, nloops

      !$omp parallel                                      &
      !$omp private(id, my_id, ii, jj, i, j, k, inbc,     &
      !$omp         nb15_calc, dij, dij_pbc, rij2, icel ) 
      !
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      my_id = my_city_rank * nthread + id
      id    = id + 1

      do icel = 1, ncel

        ii = head(icel)
        if ( ii <= 0 ) cycle

        do while ( ii > 0 )

          i = enefunc%table%solute_list(ii)

          if (mod(ii-1,nproc_city*nthread) == my_id) then

            num_nb15_calc(ii) = 0

            ! serch only in the self and neighboring cells 
            !
            do inbc = 1, boundary%num_neighbor_cells
              jj = head(boundary%neighbor_cells(inbc,icel))

              do while (jj > 0)

                if (jj <= ii) then
                  jj = list(jj)
                  cycle
                end if

                ! 1-2 and 1-3 interactions are removed
                !
                nb15_calc = .true.
                do k = 1, enefunc%table%num_nonb_excl(ii)
                  if (jj == enefunc%table%nonb_excl_list(k,ii)) then
                    nb15_calc = .false.
                    exit
                  end if
                end do

                ! 1-4 interactions are removed
                !
                do k = 1, enefunc%num_nb14_calc(ii)
                  if (jj == enefunc%nb14_calc_list(k,ii)) then
                    nb15_calc = .false.
                    exit
                  end if
                end do

                if (.not.nb15_calc) then
                  jj = list(jj)
                  cycle
                end if

                j = enefunc%table%solute_list(jj)

                ! compute distance
                !
                dij(1:3) = coord(1:3,i) - coord(1:3,j)

                ! consider the periodic boundary
                !
                dij_pbc(1:3) = box_size(1:3)*anint(dij(1:3)/box_size(1:3))
                dij(1:3)     = dij(1:3) - dij_pbc(1:3)
                rij2         = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

                ! cout the number of interactions
                !
                if (rij2 < pairdist2) then
                  num_nb15_calc(ii) = num_nb15_calc(ii) + 1

                  if (.not. do_allocate) &
                    pairlist%table%nb15_calc_list(num_nb15_calc(ii),ii) = j
                end if

                jj = list(jj)

              end do
            end do
          end if

          ii = list(ii)

        end do
      end do
      !$omp end parallel

      ! allocate memory of pairlist
      !
      if (do_allocate) then
        num_nb15_max = max(1,maxval(num_nb15_calc(1:natom)))
        num_nb15_max = int(real(num_nb15_max,wp)*FactNumNb15)
        call alloc_pairlist2(pairlist, PairListPbcSoluteSolute,           &
                             num_nb15_max, natom)
        pairlist%num_nb15_max = num_nb15_max

        do_allocate = .false.
      end if

    end do

    ! check memory size
    !
    call check_pairlist_memory_size(num_nb15_calc, pairlist%num_nb15_max, &
                                    pairlist%allocate_solsol)

    return

  end subroutine update_pairlist_solute_solute

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_solute_water
  !> @brief        update pairlist between solute and water
  !! @authors      JJ, TM, MK
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : boundary conditions information
  !! @param[in]    coord    : coordinates
  !! @param[inout] pairlist : pairlist information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_solute_water(enefunc, boundary, coord, pairlist)

    ! formal arguments
    type(s_enefunc),          intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    real(wp),                 intent(in)    :: coord(:,:)
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                  :: dij(1:3), rij2, pairdist2
    real(wp)                  :: dij_pbc(1:3), origin(1:3)
    real(wp)                  :: r_shift(1:3)
    real(wp)                  :: csize_inv(1:3)
    real(wp)                  :: box_size(1:3), half_box_size(1:3)
    integer                   :: ii, jj, i, j, k, ik, jk, n, nloops
    integer                   :: natom, nwater
    integer                   :: num_excl, ini_excl, fin_excl
    integer                   :: num_nb14, ini_nb14, fin_nb14
    integer                   :: ic(1:3), icel, inbc
    integer                   :: num_nb15_max
    integer                   :: id, my_id, nthread
    integer                   :: ncell(1:3)
#ifdef OMP
    integer                   :: omp_get_max_threads, omp_get_thread_num
#endif
    logical                   :: nb15_calc, do_allocate

    integer,          pointer :: num_nb15_pre(:), num_nb15_calcw(:)
    integer,          pointer :: ic_atm(:), head(:), list(:)
    integer,          pointer :: headw(:), listw(:)
    integer,          pointer :: ncel

    natom   = enefunc%table%num_solute
    nwater  = enefunc%table%num_water

    num_nb15_pre => pairlist%num_nb15_pre
    num_nb15_calcw => pairlist%table%num_nb15_calcw
    ic_atm       => pairlist%table%atom_cell_index
    list         => pairlist%table%cell_linked_list
    listw        => pairlist%table%cell_linked_listw
    head         => boundary%cell_head_atom
    headw        => boundary%cell_head_atomw
    ncel         => boundary%num_cells

    box_size(1)  = boundary%box_size_x
    box_size(2)  = boundary%box_size_y
    box_size(3)  = boundary%box_size_z
    origin(1)    = boundary%origin_x
    origin(2)    = boundary%origin_y
    origin(3)    = boundary%origin_z
    ncell(1)     = boundary%num_cells_x
    ncell(2)     = boundary%num_cells_y
    ncell(3)     = boundary%num_cells_z
    csize_inv(1) = 1.0_wp / boundary%cell_size_x
    csize_inv(2) = 1.0_wp / boundary%cell_size_y
    csize_inv(3) = 1.0_wp / boundary%cell_size_z
    half_box_size(1:3) = box_size(1:3)*0.5_wp

    pairdist2 = pairlist%pairlistdist * pairlist%pairlistdist

    ! initialize
    !
#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif
    num_nb15_calcw(1:natom) = 0


    ! make cell linked list of water oxygen
    !
    headw(1:ncel) = 0
    ik = 0
    do ii = 1, nwater

      ! water oxygen index
      !
      do k = 1, 3

        ik = ik + 1

        ! coordinate shifted to the first quadrant and set into the boundary box
        !
        i = enefunc%table%water_list(k,ii)
        r_shift(1:3) = coord(1:3,i) - origin(1:3)
        r_shift(1:3) = r_shift(1:3) + half_box_size(1:3) - &
                       box_size(1:3)*anint(r_shift(1:3)/box_size(1:3))

        ! assign which cell
        !
        ic(1:3) = int(r_shift(1:3)*csize_inv(1:3))

        ! upper limit of ic[x,y,z] should be ncel_[x,y,z]
        !
        if (ic(1) >= ncell(1)) ic(1) = ncell(1)-1
        if (ic(2) >= ncell(2)) ic(2) = ncell(2)-1
        if (ic(3) >= ncell(3)) ic(3) = ncell(3)-1
        icel = 1 + ic(1) + ic(2)*ncell(1) + ic(3)*ncell(1)*ncell(2)

        ! store which cell atom i is placed
        !
        listw(ik)   = headw(icel)
        headw(icel)  = ik
      end do

    end do

    num_nb15_calcw(1:natom) = 0

    ! if pairlist%allocation = .true. , allocation + make pairlist
    ! if pairlist%allocation = .false., make pairlist only
    !
    if (pairlist%allocate_solwat) then
      nloops = 2
      do_allocate = .true.
    else
      nloops = 1
      do_allocate = .false.
    end if

    do n = 1, nloops

      if (.not. do_allocate) &
        num_nb15_calcw(1:natom) = 0

      !$omp parallel                                                           &
      !$omp private(id, my_id, ii, jj, i, j, k, jk, inbc, ini_excl, fin_excl,  &
      !$omp         ini_nb14, fin_nb14, nb15_calc, dij, dij_pbc, rij2,         &
      !$omp         icel)&
      !$omp firstprivate(do_allocate)
      !
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      my_id = my_city_rank * nthread + id
      id    = id + 1

      do icel = 1, ncel

        ii = head(icel)
        if ( ii <= 0 ) cycle

        do while ( ii > 0 )

          i = enefunc%table%solute_list(ii)

          if (mod(ii-1,nproc_city*nthread) == my_id) then

            ! serch only in the self and neighboring cells 
            !
            do inbc = 1, boundary%num_neighbor_cells
              jk = headw(boundary%neighbor_cells(inbc,icel))

              do while (jk > 0)

                ! water list
                !
                jj = (jk+2) / 3
                k = mod(jk+2,3) + 1
                j = enefunc%table%water_list(k,jj)

                ! compute distance
                !
                dij(1:3) = coord(1:3,i) - coord(1:3,j)

                ! consider the periodic boundary
                !
                dij_pbc(1:3) = box_size(1:3)*anint(dij(1:3)/box_size(1:3))
                dij(1:3)     = dij(1:3) - dij_pbc(1:3)
                rij2         = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

                ! cout the number of interactions
                !
                if (rij2 < pairdist2) then
                  num_nb15_calcw(ii) = num_nb15_calcw(ii) + 1
                  if (.not. do_allocate) &
                    pairlist%table%nb15_calc_listw(num_nb15_calcw(ii),ii) = j
                end if

                jk = listw(jk)

              end do
            end do
          end if

          ii = list(ii)

        end do
      end do
      !$omp end parallel

      ! allocate memory of pairlist
      !
      if (do_allocate) then
        num_nb15_max = max(1,maxval(num_nb15_calcw(1:natom)))
        num_nb15_max = int(real(num_nb15_max,wp)*FactNumNb15)
        call alloc_pairlist2(pairlist, PairListPbcSoluteWater,   &
                             num_nb15_max, max(1,natom))
        pairlist%table%num_nb15_max = num_nb15_max

        do_allocate = .false.
      end if

    end do

    ! check memory size
    !
    call check_pairlist_memory_size(num_nb15_calcw,              &
                                    pairlist%table%num_nb15_max, &
                                    pairlist%allocate_solwat)

    return

  end subroutine update_pairlist_solute_water

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_water_water
  !> @brief        update pairlist between water and water
  !! @authors      JJ, TM, MK
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : boundary conditions information
  !! @param[in]    coord    : coordinates
  !! @param[inout] pairlist : pairlist information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_water_water(enefunc, boundary, coord, pairlist)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_boundary), target, intent(in)    :: boundary
    real(wp),                 intent(in)    :: coord(:,:)
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                  :: dij(1:3), rij2, pairdist2
    real(wp)                  :: dij_pbc(1:3), origin(1:3)
    real(wp)                  :: r_shift(1:3)
    real(wp)                  :: csize_inv(1:3)
    real(wp)                  :: box_size(1:3), half_box_size(1:3)
    integer                   :: ii, jj, i, j, k, n, nloops
    integer                   :: nwater
    integer                   :: ic(1:3), ncell(1:3), icel, inbc
    integer                   :: num_nb15_max
    integer                   :: id, my_id, nthread
#ifdef OMP
    integer                   :: omp_get_max_threads, omp_get_thread_num
#endif
    logical                   :: do_allocate

    integer,          pointer :: num_nb15_calc_water(:)
    integer,          pointer :: listw(:), ic_atmw(:), headw(:)
    integer,          pointer :: ncel


    nwater       =  enefunc%table%num_water
    num_nb15_calc_water => pairlist%table%num_nb15_calc_water
    ic_atmw      => pairlist%table%atom_cell_index_water
    listw        => pairlist%table%cell_linked_list_water
    headw        => boundary%cell_head_atomw
    ncel         => boundary%num_cells

    box_size(1)  =  boundary%box_size_x
    box_size(2)  =  boundary%box_size_y
    box_size(3)  =  boundary%box_size_z
    csize_inv(1) =  1.0_wp / boundary%cell_size_x
    csize_inv(2) =  1.0_wp / boundary%cell_size_y
    csize_inv(3) =  1.0_wp / boundary%cell_size_z
    ncell(1)     =  boundary%num_cells_x
    ncell(2)     =  boundary%num_cells_y
    ncell(3)     =  boundary%num_cells_z
    origin(1)    =  boundary%origin_x
    origin(2)    =  boundary%origin_y
    origin(3)    =  boundary%origin_z
    half_box_size(1:3) = box_size(1:3)*0.5_wp

    pairdist2 = (pairlist%pairlistdist + TableWaterMargin) &
              * (pairlist%pairlistdist + TableWaterMargin)

    ! initialize
    !
#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif


    ! make cell linked list for water
    !
    headw(1:ncel) = 0
    do ii = 1, nwater

      ! coordinate shifted to the first quadrant and set into the boundary box
      !
      i = enefunc%table%water_list(1,ii)
      r_shift(1:3) = coord(1:3,i) - origin(1:3)
      r_shift(1:3) = r_shift(1:3) + half_box_size(1:3) - &
                     box_size(1:3)*anint(r_shift(1:3)/box_size(1:3))

      ! assign which cell
      !
      ic(1:3) = int(r_shift(1:3)*csize_inv(1:3))

      ! upper limit of ic[x,y,z] should be ncel_[x,y,z]
      !
      do k = 1, 3
        if (ic(k) >= ncell(k)) ic(k) = ncell(k)-1
      end do
      icel = 1 + ic(1) + ic(2)*ncell(1) + ic(3)*ncell(1)*ncell(2)

      ! store which cell atom i is placed
      !
      ic_atmw(ii) = icel
      listw(ii)   = headw(icel)
      headw(icel)  = ii

    end do

    num_nb15_calc_water(1:nwater) = 0


    ! if pairlist%allocation = .true. , allocation + make pairlist
    ! if pairlist%allocation = .false., make pairlist only
    !
    if (pairlist%allocate_watwat) then
      nloops = 2
      do_allocate = .true.
    else
      nloops = 1
      do_allocate = .false.
    end if

    do n = 1, nloops

      if (.not. do_allocate) &
        num_nb15_calc_water(1:nwater) = 0

      !$omp parallel                               &
      !$omp private(id, my_id, ii, jj, i, j, inbc, &
      !$omp         dij, dij_pbc, rij2, icel)      &
      !$omp firstprivate(do_allocate)
      !
#ifdef OMP
      id = omp_get_thread_num()
#else
      id = 0
#endif
      my_id = my_city_rank * nthread + id
      id    = id + 1

      do icel = 1, ncel

        ii = headw(icel)
        if ( ii <= 0 ) cycle

        do while ( ii > 0 )

          if (mod(ii-1,nproc_city*nthread) == my_id) then

            i = enefunc%table%water_list(1,ii)

            ! serch only in the self and neighboring cells
            !
            do inbc = 1, boundary%num_neighbor_cells

              jj = headw(boundary%neighbor_cells(inbc,icel))

              do while (jj > 0)

                if (jj <= ii) then
                  jj = listw(jj)
                  cycle
                end if

                j = enefunc%table%water_list(1,jj)

                ! compute distance
                !
                dij(1:3) = coord(1:3,i) - coord(1:3,j)

                ! consider the periodic boundary
                !
                dij_pbc(1:3) = box_size(1:3)*anint(dij(1:3)/box_size(1:3))
                dij(1:3)     = dij(1:3) - dij_pbc(1:3)
                rij2         = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

                ! cout the number of interactions
                !
                if (rij2 < pairdist2) then
                  num_nb15_calc_water(ii) = num_nb15_calc_water(ii) + 1

                  if (.not. do_allocate) &
                    pairlist%table%nb15_calc_list_water(num_nb15_calc_water(ii),ii) = jj
                end if

                jj = listw(jj)
              end do
            end do
          end if

          ii = listw(ii)

        end do
      end do
      !$omp end parallel

      ! allocate memory of pairlist
      !
      if (do_allocate) then
        num_nb15_max = max(1,maxval(num_nb15_calc_water(1:nwater)))
        num_nb15_max = int(real(num_nb15_max,wp)*FactNumNb15)
        call alloc_pairlist2(pairlist, PairListPbcWaterWater,          &
                             num_nb15_max, max(1,nwater))
        pairlist%table%water_nb15_max = num_nb15_max

        do_allocate = .false.
      end if

    end do

    ! check memory size
    !
    call check_pairlist_memory_size(num_nb15_calc_water,           &
                                    pairlist%table%water_nb15_max, &
                                    pairlist%allocate_watwat)

    return

  end subroutine update_pairlist_water_water

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    check_pairlist_memory_size
  !> @brief        check memory size of pairlist
  !! @authors      TM
  !! @param[in]    num_nb15       : number of 1-5 nonbonded pairs
  !! @param[in]    num_nb15_limit : already allocated size
  !! @param[out]   alloc_list     : allocate list or not
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_pairlist_memory_size(num_nb15_calc, num_nb15_limit,         &
                                        alloc_list)

    ! formal arguments
    integer,          intent(in)    :: num_nb15_calc(:)
    integer,          intent(in)    :: num_nb15_limit
    logical,          intent(out)   :: alloc_list

    ! local variables
    integer                  :: threshold
    integer                  :: num_nb15_max

    num_nb15_max = 0
    num_nb15_max = max(num_nb15_max, maxval(num_nb15_calc(:)))

    ! check memory size
    !   if the num_nb15_max is larger than FactThreshold*(allocated size),
    !   allocate memory of pairlist at the next update
    !
    threshold = int(real(num_nb15_limit,wp)*FactThreshold)

    if (num_nb15_max > threshold) then
      alloc_list = .true.
    else
      alloc_list = .false.
    end if

    return

  end subroutine check_pairlist_memory_size

end module at_pairlist_mod
