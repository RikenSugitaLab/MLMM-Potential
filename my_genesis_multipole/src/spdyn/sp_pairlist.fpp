!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_pairlist_mod
!> @brief   set pairlist for nonbonded interactions
!! @authors Jaewoon Jung (JJ), Kiyotaka Sakamoto (KS)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_pairlist_mod

  use sp_pairlist_str_mod
  use sp_boundary_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use molecules_str_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: setup_pairlist
  public  :: update_pairlist_pbc
  public  :: update_pairlist_nobc
  public  :: update_pairlist_pbc_check
#ifdef USE_GPU
  public  :: update_pairlist_pbc_univ
#endif
  !FEP
  public  :: setup_pairlist_fep
  public  :: update_pairlist_pbc_fep
  public  :: update_pairlist_pbc_check_fep

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_pairlist
  !> @brief        initialize/allocate/setup pairlist for nonbonded interactions
  !! @authors      JJ
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    domain   : domain information
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_pairlist(enefunc, domain, pairlist)

    ! formal arguments
    type(s_enefunc),         intent(in)    :: enefunc 
    type(s_domain),          intent(inout) :: domain 
    type(s_pairlist),        intent(inout) :: pairlist
    

    call init_pairlist(pairlist)

    pairlist%pairlistdist = enefunc%pairlistdist

    if (enefunc%forcefield == ForcefieldAAGO) then

      call alloc_pairlist(pairlist, PairListNOBC, &
        domain%num_cell_local+domain%num_cell_boundary)

      call timer(TimerPairList, TimerOn)
      call update_pairlist_nobc(enefunc, domain, pairlist)
      call timer(TimerPairList, TimerOff)

    else

#ifdef KCOMP
      if (enefunc%forcefield /= ForcefieldGROMARTINI) then
        call alloc_pairlist(pairlist, PairListNoTable, domain%num_cell_local)
      else
        call alloc_pairlist(pairlist, PairListNoTable_CG, domain%num_cell_local)
      endif
#else
      call alloc_pairlist(pairlist, PairListNoTable, domain%num_cell_local)
#endif

      call timer(TimerPairList, TimerOn)
#ifndef USE_GPU
      if (enefunc%pairlist_check) then
        call update_pairlist_pbc_check(enefunc, domain, pairlist)
      else
        call update_pairlist_pbc(enefunc, domain, pairlist)
      endif
#else
      call update_pairlist_pbc_univ(enefunc, domain, pairlist)
#endif
      call timer(TimerPairList, TimerOff)

    end if

    return

  end subroutine setup_pairlist


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_pbc
  !> @brief        update pairlist in each domain with periodic boundary 
  !!               condition 
  !! @authors      JJ
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    domain   : domain information
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_pbc(enefunc, domain, pairlist)
  
    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_domain),   target, intent(in)    :: domain  
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                  :: pairdist2
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: rtmp(1:3)

    integer                   :: i, j, ij, k, ix, iy
    integer                   :: num_nb15_pre, num_nb15
    integer                   :: num_nb15_total, num_nb15_total1
    integer                   :: num_excl, ini_excl, fin_excl
    integer                   :: num_nb14, ini_nb14, fin_nb14
    logical                   :: nb15_calc
    integer                   :: id, omp_get_thread_num, num_nb15_cell

    real(dp),         pointer :: coord(:,:,:)
    real(wp),         pointer :: cell_move(:,:,:), system_size(:)
    real(wp),         pointer :: trans1(:,:,:), trans2(:,:,:)
    integer,          pointer :: num_nonb_excl1(:,:), num_nb14_calc1(:,:)
    integer,          pointer :: num_nonb_excl(:,:), num_nb14_calc(:,:)
    integer,          pointer :: nonb_excl_list1(:,:), nb14_calc_list1(:,:)
    integer,          pointer :: nonb_excl_list(:,:), nb14_calc_list(:,:)
    integer,          pointer :: natom(:), cell_pairlist1(:,:)
    integer,          pointer :: ncell, nboundary
    integer,          pointer :: nb15_cell(:), nb15_list(:,:)
    integer,          pointer :: nb15_calc_list1(:,:), nb15_calc_list(:,:)
    integer,          pointer :: num_nb15_calc1(:,:), num_nb15_calc(:,:)

    num_nonb_excl1  => enefunc%num_nonb_excl1
    num_nb14_calc1  => enefunc%num_nb14_calc1
    nonb_excl_list1 => enefunc%nonb_excl_list1
    nb14_calc_list1 => enefunc%nb14_calc_list1
    num_nonb_excl   => enefunc%num_nonb_excl
    num_nb14_calc   => enefunc%num_nb14_calc
    nonb_excl_list  => enefunc%nonb_excl_list
    nb14_calc_list  => enefunc%nb14_calc_list

    ncell           => domain%num_cell_local
    nboundary       => domain%num_cell_boundary
    natom           => domain%num_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    trans2          => domain%translated
    cell_pairlist1  => domain%cell_pairlist1
    cell_move       => domain%cell_move
    system_size     => domain%system_size

    num_nb15_calc1  => pairlist%num_nb15_calc1
    num_nb15_calc   => pairlist%num_nb15_calc
    nb15_list       => pairlist%nb15_list
    nb15_cell       => pairlist%nb15_cell
    nb15_calc_list1 => pairlist%nb15_calc_list1
    nb15_calc_list  => pairlist%nb15_calc_list

    pairdist2       =  pairlist%pairlistdist * pairlist%pairlistdist
    num_nb15_total  =  0

    !$omp parallel                                                             &
    !$omp private(id, i, ix, ini_excl, num_excl, ini_nb14, num_nb14, num_nb15, &
    !$omp         num_nb15_pre, fin_excl, fin_nb14, iy, k, nb15_calc, ij, j,   &
    !$omp         rtmp, dij, rij2, num_nb15_cell)                              &
    !$omp reduction(+:num_nb15_total)

#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell+nboundary, nthread
      do ix = 1, natom(i)
        trans2(1:3,ix,i) = coord(1:3,ix,i) + trans1(1:3,ix,i)
      end do
    end do

    !$omp barrier

    ! make a pairlist in the same cell
    !
    do i = id+1, ncell, nthread

      ini_excl = 0
      num_excl = 0
      ini_nb14 = 0
      num_nb14 = 0
      num_nb15 = 0
      num_nb15_pre = 0

#ifdef DEBUG
      if (natom(i) > MaxAtom) &
        call error_msg( &
             'Debug: Update_Pairlist_Pbc> natom(cell) is exceed MaxAtom')
#endif

      do ix = 1, natom(i) - 1

        ini_excl = num_excl + 1
        fin_excl = num_excl + num_nonb_excl1(ix,i)
        num_excl = fin_excl

        ini_nb14 = num_nb14 + 1
        fin_nb14 = num_nb14 + num_nb14_calc1(ix,i)
        num_nb14 = fin_nb14

        do iy = ix + 1, natom(i)

          nb15_calc = .true.
          do k = ini_excl, fin_excl
            if (iy == nonb_excl_list1(k,i)) then
              nb15_calc = .false.
              exit
            end if
          end do

          do k = ini_nb14, fin_nb14
            if (iy == nb14_calc_list1(k,i)) then
              nb15_calc = .false.
              exit
            end if
          end do

          if (nb15_calc) then

            ! store interaction table
            !
            num_nb15 = num_nb15 + 1
#ifdef DEBUG
            if (num_nb15 > MaxNb15_chk) &
              call error_msg( &
                   'Debug: Update_Pairlist_Pbc> num_nb15 is exceed MaxNb15')
#endif
            nb15_calc_list1(num_nb15,i) = iy
          end if
        end do

        num_nb15_calc1(ix,i) = num_nb15 - num_nb15_pre
        num_nb15_total = num_nb15_total + num_nb15 - num_nb15_pre
        num_nb15_pre = num_nb15

      end do
    end do 

    ! Make a pairlist between different cells
    !
    do ij = id+1, maxcell_near, nthread

      i = cell_pairlist1(1,ij)
      j = cell_pairlist1(2,ij)

#ifdef DEBUG
      if (natom(i) > MaxAtom .or. natom(j) > MaxAtom) &
        call error_msg( &
             'Debug: Update_Pairlist_Pbc> natom(cell) is exceed MaxAtom')
#endif

      ini_excl = 0
      num_excl = 0
      ini_nb14 = 0
      num_nb14 = 0
      num_nb15 = 0
      num_nb15_pre = 0
      num_nb15_cell = 0

      do ix = 1, natom(i)

        ini_excl = num_excl + 1
        fin_excl = num_excl + num_nonb_excl(ix,ij)
        num_excl = fin_excl

        ini_nb14 = num_nb14 + 1
        fin_nb14 = num_nb14 + num_nb14_calc(ix,ij)
        num_nb14 = fin_nb14

        rtmp(1:3) = trans2(1:3,ix,i) + cell_move(1:3,j,i)*system_size(1:3)

        do iy = 1, natom(j)

          nb15_calc = .true.

          do k = ini_excl, fin_excl
            if (iy == nonb_excl_list(k,ij)) then
              nb15_calc = .false.
              exit
            end if
          end do

          do k = ini_nb14, fin_nb14
            if (iy == nb14_calc_list(k,ij)) then
              nb15_calc = .false.
              exit
            end if
          end do

          if (nb15_calc) then

            dij(1:3) = rtmp(1:3) - trans2(1:3,iy,j) 

            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            ! store interaction table
            !
            if (rij2 < pairdist2) then

              num_nb15 = num_nb15 + 1

#ifdef DEBUG
              if (num_nb15 > MaxNb15_chk) &
                call error_msg( &
                     'Debug: Update_Pairlist_Pbc> num_nb15 is exceed MaxNb15')
#endif
              nb15_calc_list(num_nb15,ij) = iy
            end if
          end if
        end do

        num_nb15_calc(ix,ij) = num_nb15 - num_nb15_pre

        if (num_nb15 /= num_nb15_pre) then
          num_nb15_cell = num_nb15_cell + 1
          nb15_list(num_nb15_cell,ij) = ix
        end if

        num_nb15_total = num_nb15_total + num_nb15 - num_nb15_pre
        num_nb15_pre = num_nb15

      end do

      nb15_cell(ij) = num_nb15_cell

    end do

    do ij = id+maxcell_near+1, maxcell, nthread

      i = cell_pairlist1(1,ij)
      j = cell_pairlist1(2,ij)

#ifdef DEBUG
      if (natom(i) > MaxAtom .or. natom(j) > MaxAtom) &
        call error_msg( &
             'Debug: Update_Pairlist_Pbc> natom(cell) is exceed MaxAtom')
#endif

      num_nb15 = 0
      num_nb15_pre = 0
      num_nb15_cell = 0

      do ix = 1, natom(i)

        rtmp(1:3) = trans2(1:3,ix,i) + cell_move(1:3,j,i)*system_size(1:3)

        do iy = 1, natom(j)

          nb15_calc = .true.

          if (nb15_calc) then

            dij(1:3) = rtmp(1:3) - trans2(1:3,iy,j)

            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            ! store interaction table
            !
            if (rij2 < pairdist2) then 

              num_nb15 = num_nb15 + 1

#ifdef DEBUG
              if (num_nb15 > MaxNb15_chk) &
                call error_msg( &
                     'Debug: Update_Pairlist_Pbc> num_nb15 is exceed MaxNb15')
#endif
              nb15_calc_list(num_nb15,ij) = iy
            end if
          end if
        end do

        num_nb15_calc(ix,ij) = num_nb15 - num_nb15_pre

        if (num_nb15 /= num_nb15_pre) then
          num_nb15_cell = num_nb15_cell + 1
          nb15_list(num_nb15_cell,ij) = ix
        end if

        num_nb15_total = num_nb15_total + num_nb15 - num_nb15_pre
        num_nb15_pre = num_nb15

      end do

      nb15_cell(ij) = num_nb15_cell

    end do

    !$omp end parallel

#ifdef HAVE_MPI_GENESIS
    call mpi_reduce(num_nb15_total, num_nb15_total1, 1, mpi_integer, mpi_sum, &
         0, mpi_comm_country, ierror)
#endif

    return

  end subroutine update_pairlist_pbc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_pbc_check
  !> @brief        update pairlist in each domain with periodic boundary 
  !!               condition 
  !! @authors      JJ,CK
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    domain   : domain information
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_pbc_check(enefunc, domain, pairlist)
  
    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_domain),   target, intent(in)    :: domain  
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                  :: pairdist2
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: rtmp(1:3)

    integer                   :: i, j, ij, k, ix, iy
    integer                   :: num_nb15_pre, num_nb15
    integer                   :: num_nb15_total, num_nb15_total1
    integer                   :: num_excl, ini_excl, fin_excl
    integer                   :: num_nb14, ini_nb14, fin_nb14
    integer                   :: small_contact
    logical                   :: nb15_calc
    integer                   :: id, omp_get_thread_num, num_nb15_cell

    real(dp),         pointer :: coord(:,:,:)
    real(wp),         pointer :: cell_move(:,:,:), system_size(:)
    real(wp),         pointer :: trans1(:,:,:), trans2(:,:,:)
    real(wp),         pointer :: err_minimum_contact
    integer,          pointer :: num_nonb_excl1(:,:), num_nb14_calc1(:,:)
    integer,          pointer :: num_nonb_excl(:,:), num_nb14_calc(:,:)
    integer,          pointer :: nonb_excl_list1(:,:), nb14_calc_list1(:,:)
    integer,          pointer :: nonb_excl_list(:,:), nb14_calc_list(:,:)
    integer,          pointer :: natom(:), cell_pairlist1(:,:)
    integer,          pointer :: ncell, nboundary
    integer,          pointer :: nb15_cell(:), nb15_list(:,:)
    integer,          pointer :: nb15_calc_list1(:,:), nb15_calc_list(:,:)
    integer,          pointer :: num_nb15_calc1(:,:), num_nb15_calc(:,:)
    integer,          pointer :: table_order
    logical,          pointer :: nonb_limiter, table

    num_nonb_excl1  => enefunc%num_nonb_excl1
    num_nb14_calc1  => enefunc%num_nb14_calc1
    nonb_excl_list1 => enefunc%nonb_excl_list1
    nb14_calc_list1 => enefunc%nb14_calc_list1
    num_nonb_excl   => enefunc%num_nonb_excl
    num_nb14_calc   => enefunc%num_nb14_calc
    nonb_excl_list  => enefunc%nonb_excl_list
    nb14_calc_list  => enefunc%nb14_calc_list

    table               => enefunc%table%table
    table_order         => enefunc%table%table_order
    nonb_limiter        => enefunc%nonb_limiter
    err_minimum_contact => enefunc%err_minimum_contact

    ncell           => domain%num_cell_local
    nboundary       => domain%num_cell_boundary
    natom           => domain%num_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    trans2          => domain%translated
    cell_pairlist1  => domain%cell_pairlist1
    cell_move       => domain%cell_move
    system_size     => domain%system_size

    num_nb15_calc1  => pairlist%num_nb15_calc1
    num_nb15_calc   => pairlist%num_nb15_calc
    nb15_list       => pairlist%nb15_list
    nb15_cell       => pairlist%nb15_cell
    nb15_calc_list1 => pairlist%nb15_calc_list1
    nb15_calc_list  => pairlist%nb15_calc_list

    pairdist2       =  pairlist%pairlistdist * pairlist%pairlistdist
    num_nb15_total  =  0
    small_contact   =  0

    !$omp parallel                                                             &
    !$omp private(id, i, ix, ini_excl, num_excl, ini_nb14, num_nb14, num_nb15, &
    !$omp         num_nb15_pre, fin_excl, fin_nb14, iy, k, nb15_calc, ij, j,   &
    !$omp         rtmp, dij, rij2, num_nb15_cell)                              &
    !$omp reduction(+:num_nb15_total) reduction(+:small_contact)

#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell+nboundary, nthread
      do ix = 1, natom(i)
        trans2(1:3,ix,i) = coord(1:3,ix,i) + trans1(1:3,ix,i)
      end do
    end do

    !$omp barrier

    ! make a pairlist in the same cell
    !
    do i = id+1, ncell, nthread

      ini_excl = 0
      num_excl = 0
      ini_nb14 = 0
      num_nb14 = 0
      num_nb15 = 0
      num_nb15_pre = 0

      if (natom(i) > MaxAtom) &
        call error_msg( &
          'Debug: Update_Pairlist_Pbc_Check> natom(cell) is exceed MaxAtom')

      do ix = 1, natom(i) - 1

        ini_excl = num_excl + 1
        fin_excl = num_excl + num_nonb_excl1(ix,i)
        num_excl = fin_excl

        ini_nb14 = num_nb14 + 1
        fin_nb14 = num_nb14 + num_nb14_calc1(ix,i)
        num_nb14 = fin_nb14

        do iy = ix + 1, natom(i)

          nb15_calc = .true.
          do k = ini_excl, fin_excl
            if (iy == nonb_excl_list1(k,i)) then
              nb15_calc = .false.
              exit
            end if
          end do

          do k = ini_nb14, fin_nb14
            if (iy == nb14_calc_list1(k,i)) then
              nb15_calc = .false.
              exit
            end if
          end do

          if (nb15_calc) then

            ! store interaction table
            !
            num_nb15 = num_nb15 + 1
            if (table .and. table_order == 1 )then
              dij(1:3) = trans2(1:3,ix,i) - trans2(1:3,iy,i)
              rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
              if (rij2 < err_minimum_contact) then
                if (.not. nonb_limiter) &
                  call error_msg( &
                  'Debug: Update_Pairlist_Pbc_Check> too small contact')
                small_contact = small_contact + 1 
              endif
            endif
            if (num_nb15 > MaxNb15_chk) &
              call error_msg( &
                 'Debug: Update_Pairlist_Pbc_Check> num_nb15 is exceed MaxNb15')

            nb15_calc_list1(num_nb15,i) = iy
          end if
        end do

        num_nb15_calc1(ix,i) = num_nb15 - num_nb15_pre
        num_nb15_total = num_nb15_total + num_nb15 - num_nb15_pre
        num_nb15_pre = num_nb15

      end do
    end do 

    ! Make a pairlist between different cells
    !
    do ij = id+1, maxcell_near, nthread

      i = cell_pairlist1(1,ij)
      j = cell_pairlist1(2,ij)

      if (natom(i) > MaxAtom .or. natom(j) > MaxAtom) &
        call error_msg( &
             'Debug: Update_Pairlist_Pbc_Check> natom(cell) is exceed MaxAtom')

      ini_excl = 0
      num_excl = 0
      ini_nb14 = 0
      num_nb14 = 0
      num_nb15 = 0
      num_nb15_pre = 0
      num_nb15_cell = 0

      do ix = 1, natom(i)

        ini_excl = num_excl + 1
        fin_excl = num_excl + num_nonb_excl(ix,ij)
        num_excl = fin_excl

        ini_nb14 = num_nb14 + 1
        fin_nb14 = num_nb14 + num_nb14_calc(ix,ij)
        num_nb14 = fin_nb14

        rtmp(1:3) = trans2(1:3,ix,i) + cell_move(1:3,j,i)*system_size(1:3)

        do iy = 1, natom(j)

          nb15_calc = .true.

          do k = ini_excl, fin_excl
            if (iy == nonb_excl_list(k,ij)) then
              nb15_calc = .false.
              exit
            end if
          end do

          do k = ini_nb14, fin_nb14
            if (iy == nb14_calc_list(k,ij)) then
              nb15_calc = .false.
              exit
            end if
          end do

          if (nb15_calc) then

            dij(1:3) = rtmp(1:3) - trans2(1:3,iy,j) 

            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            ! store interaction table
            !
            if (rij2 < pairdist2) then

              num_nb15 = num_nb15 + 1

              if (table .and. table_order == 1 )then
                if (rij2 < err_minimum_contact) then
                  if (.not. nonb_limiter) &
                    call error_msg( &
                    'Debug: Update_Pairlist_Pbc_Check> too small contact')
                  small_contact = small_contact + 1 
                endif
              endif
              if (num_nb15 > MaxNb15_chk) &
                call error_msg( &
                 'Debug: Update_Pairlist_Pbc_Check> num_nb15 is exceed MaxNb15')

              nb15_calc_list(num_nb15,ij) = iy
            end if
          end if
        end do

        num_nb15_calc(ix,ij) = num_nb15 - num_nb15_pre

        if (num_nb15 /= num_nb15_pre) then
          num_nb15_cell = num_nb15_cell + 1
          nb15_list(num_nb15_cell,ij) = ix
        end if

        num_nb15_total = num_nb15_total + num_nb15 - num_nb15_pre
        num_nb15_pre = num_nb15

      end do

      nb15_cell(ij) = num_nb15_cell

    end do

    do ij = id+maxcell_near+1, maxcell, nthread

      i = cell_pairlist1(1,ij)
      j = cell_pairlist1(2,ij)

      if (natom(i) > MaxAtom .or. natom(j) > MaxAtom) &
        call error_msg( &
           'Debug: Update_Pairlist_Pbc_Check> natom(cell) is exceed MaxAtom')

      num_nb15 = 0
      num_nb15_pre = 0
      num_nb15_cell = 0

      do ix = 1, natom(i)

        rtmp(1:3) = trans2(1:3,ix,i) + cell_move(1:3,j,i)*system_size(1:3)

        do iy = 1, natom(j)

          nb15_calc = .true.

          if (nb15_calc) then

            dij(1:3) = rtmp(1:3) - trans2(1:3,iy,j)

            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            ! store interaction table
            !
            if (rij2 < pairdist2) then 

              num_nb15 = num_nb15 + 1

              if (num_nb15 > MaxNb15_chk) &
                call error_msg( &
                 'Debug: Update_Pairlist_Pbc_Check> num_nb15 is exceed MaxNb15')

              nb15_calc_list(num_nb15,ij) = iy
            end if
          end if
        end do

        num_nb15_calc(ix,ij) = num_nb15 - num_nb15_pre

        if (num_nb15 /= num_nb15_pre) then
          num_nb15_cell = num_nb15_cell + 1
          nb15_list(num_nb15_cell,ij) = ix
        end if

        num_nb15_total = num_nb15_total + num_nb15 - num_nb15_pre
        num_nb15_pre = num_nb15

      end do

      nb15_cell(ij) = num_nb15_cell

    end do

    !$omp end parallel

    if (small_contact > 0) then
      write(MsgOut, *) "Warning: small contacts exist in inner loop"
    endif

#ifdef HAVE_MPI_GENESIS
    call mpi_reduce(num_nb15_total, num_nb15_total1, 1, mpi_integer, mpi_sum, &
         0, mpi_comm_country, ierror)
#endif

    return

  end subroutine update_pairlist_pbc_check

#ifdef USE_GPU
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_pbc_univ
  !> @brief        update pairlist in each domain with periodic boundary
  !!               condition for GPU
  !! @authors      JJ, Naruse, KS
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    domain   : domain information
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_pbc_univ(enefunc, domain, pairlist)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_domain),   target, intent(inout) :: domain
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                  :: pairdist2
    real(wp)                  :: cutoffdist2
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: rtmp(1:3)

    integer                   :: i, j, ij, k, ix, iy
    integer                   :: num_nb15_pre, num_nb15
    integer                   :: num_nb15_total, num_nb15_total1
    integer                   :: num_excl, ini_excl, fin_excl
    integer                   :: num_nb14, ini_nb14, fin_nb14
    logical                   :: nb15_calc
    integer                   :: id, omp_get_thread_num

    real(dp),         pointer :: coord(:,:,:)
    real(wp),         pointer :: cell_move(:,:,:), system_size(:)
    real(wp),         pointer :: trans1(:,:,:), trans2(:,:,:)
    integer,          pointer :: num_nonb_excl1(:,:), num_nb14_calc1(:,:)
    integer,          pointer :: num_nonb_excl(:,:), num_nb14_calc(:,:)
    integer,          pointer :: nonb_excl_list1(:,:), nb14_calc_list1(:,:)
    integer,          pointer :: nonb_excl_list(:,:), nb14_calc_list(:,:)
    integer,          pointer :: natom(:), cell_pairlist1(:,:)
    integer,          pointer :: ncell, nboundary
    real(wp)                  :: dij1, dij2, dij3
    real(wp)                  :: rtmp1, rtmp2, rtmp3

    ! for GPU
    integer                   :: univ_ij
    integer                   :: num_all, num_target, num_calc
    integer                   :: max_atoms, max_load
    integer,      allocatable :: index(:), count(:)
    integer                   :: id_load, iix, iiy
    integer                   :: ix_natom, iy_natom
    integer                   :: idx, base_idx
    integer                   :: load_ij, load_ji
    integer                   :: all_iter_32x1, all_iter_16x2, all_iter_8x4
    integer                   :: act_iter
    real(8)                   :: rate_iter
    integer(1)                :: mark
    integer                   :: max_natom, univ_mask2_size
    integer(1),       pointer :: univ_ix_list(:,:), univ_iy_list(:,:)
    integer,          pointer :: univ_ix_natom(:), univ_iy_natom(:)
    integer,          pointer :: univ_cell_pairlist1(:,:)
    integer,          pointer :: univ_cell_pairlist2(:,:)

    integer                   :: ncell_local, ncell_bound
    integer                   :: ret_shape(1), ncell_max
    integer                   :: num
    integer                   :: univ_ix_num(MaxAtom), univ_iy_num(MaxAtom)

    integer                   :: pack_univ_mask2_size
    integer                   :: return_size

    num_nonb_excl1      => enefunc%num_nonb_excl1
    num_nb14_calc1      => enefunc%num_nb14_calc1
    nonb_excl_list1     => enefunc%nonb_excl_list1
    nb14_calc_list1     => enefunc%nb14_calc_list1
    num_nonb_excl       => enefunc%num_nonb_excl
    num_nb14_calc       => enefunc%num_nb14_calc
    nonb_excl_list      => enefunc%nonb_excl_list
    nb14_calc_list      => enefunc%nb14_calc_list

    ncell               => domain%num_cell_local
    nboundary           => domain%num_cell_boundary
    natom               => domain%num_atom
    coord               => domain%coord
    trans1              => domain%trans_vec
    trans2              => domain%translated
    cell_move           => domain%cell_move
    system_size         => domain%system_size
    cell_pairlist1      => domain%cell_pairlist1

    pairdist2           =  pairlist%pairlistdist * pairlist%pairlistdist
    cutoffdist2         =  enefunc%cutoffdist * enefunc%cutoffdist

    univ_ix_list        => pairlist%univ_ix_list
    univ_iy_list        => pairlist%univ_iy_list
    univ_ix_natom       => pairlist%univ_ix_natom
    univ_iy_natom       => pairlist%univ_iy_natom
    univ_cell_pairlist1 => domain%univ_cell_pairlist1
    univ_cell_pairlist2 => domain%univ_cell_pairlist2
    ncell_local         =  domain%num_cell_local
    ncell_bound         =  domain%num_cell_boundary

    ! allocation of mask
    !
    max_natom = 0
    do i = 1, ncell+nboundary
      if ( max_natom < natom(i) ) then
        max_natom = natom(i)
      endif
    enddo
#ifdef DEBUG
    if ( max_natom > MaxAtom ) then
      call error_msg('Debug: Update_Pairlist_Pbc> natom(cell) is exceed MaxAtom')
    endif
#endif
    max_load = 16 * ((max_natom+7)/8 * (max_natom+3)/4) + max_natom
    allocate( index(0:max_load), count(0:max_load) )

    univ_mask2_size = max_natom * max_natom
    if ( pairlist%univ_mask2_size < univ_mask2_size ) then
      if ( allocated(pairlist%univ_mask2) ) then
        deallocate( pairlist%univ_mask2)
        call unset_pinned_memory(pairlist%pack_univ_mask2)
        deallocate( pairlist%pack_univ_mask2)
      endif
      max_natom = max_natom + max_natom/10  ! to reduce number of memory allocation
      univ_mask2_size = max_natom * max_natom
      allocate(pairlist%univ_mask2(univ_mask2_size,univ_ncell_near))
      pairlist%univ_mask2_size = univ_mask2_size

      pack_univ_mask2_size = (univ_mask2_size * univ_ncell_near + 7)/8
      allocate(pairlist%pack_univ_mask2(pack_univ_mask2_size))
      call set_pinned_memory(pairlist%pack_univ_mask2, pack_univ_mask2_size)
      pairlist%pack_univ_mask2_size = pack_univ_mask2_size
    endif
    univ_natom_max = max_natom

    !
    ! ******** Initialization of GPU ********
    !
    call gpu_init()

    !
    ! ******** make pairlist on GPU ********
    !
    ret_shape = shape(natom)
    ncell_max = ret_shape(1)

    call gpu_launch_build_pairlist( &
         trans2, coord, trans1, cell_move, natom, univ_cell_pairlist1, &
         univ_ix_list, univ_iy_list, univ_ix_natom, univ_iy_natom,     &
         MaxAtom, ncell_local, ncell_bound, ncell_max, univ_maxcell,   &
         univ_maxcell1, pairdist2, cutoffdist2,                        &
         system_size(1), system_size(2), system_size(3) )

    ! Initialization of mask on CPU
    !
    !$omp parallel do default(shared) private(univ_ij, i, j, ix, iy, idx)
    do univ_ij =1, univ_ncell_near
       i = univ_cell_pairlist1(1,univ_ij)
       j = univ_cell_pairlist1(2,univ_ij)
       do ix = 1, natom(i)
          do iy = 1, natom(j)
             idx = 1 + (iy-1) + (ix-1)*max_natom
             pairlist%univ_mask2(idx, univ_ij) = 1
          enddo
       enddo
    enddo
    !$omp end parallel do

    ! Mask calculation on CPU by exclusion
    !
    !$omp parallel default(shared)                    &
    !$omp private(id, idx, i, ix, iy, k, ij, univ_ij, &
    !$omp         ini_excl, fin_excl, num_excl,       &
    !$omp         ini_nb14, fin_nb14, num_nb14)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell, nthread

      ini_excl = 0
      num_excl = 0
      ini_nb14 = 0
      num_nb14 = 0

      do ix = 1, natom(i)-1

        ini_excl = num_excl + 1
        fin_excl = num_excl + num_nonb_excl1(ix,i)
        num_excl = fin_excl

        ini_nb14 = num_nb14 + 1
        fin_nb14 = num_nb14 + num_nb14_calc1(ix,i)
        num_nb14 = fin_nb14

        do k = ini_excl, fin_excl
          iy = nonb_excl_list1(k,i)
          idx = iy + (ix-1)*max_natom
          pairlist%univ_mask2(idx,i) = 0
          idx = ix + (iy-1)*max_natom
          pairlist%univ_mask2(idx,i) = 0
        end do
        do k = ini_nb14, fin_nb14
          iy = nb14_calc_list1(k,i)
          idx = iy + (ix-1)*max_natom
          pairlist%univ_mask2(idx,i) = 0
          idx = ix + (iy-1)*max_natom
          pairlist%univ_mask2(idx,i) = 0
        end do
        idx = ix + (ix-1)*max_natom
        pairlist%univ_mask2(idx,i) = 0

      end do

      if (natom(i) >= 1) then
        idx = natom(i) + (natom(i)-1)*max_natom
        pairlist%univ_mask2(idx,i) = 0
      end if

    end do

    do ij = id+1, maxcell_near, nthread

      i = cell_pairlist1(1,ij)
      j = cell_pairlist1(2,ij)
      univ_ij = ij + ncell

      ini_excl = 0
      num_excl = 0
      ini_nb14 = 0
      num_nb14 = 0

      do ix = 1, natom(i)

        ini_excl = num_excl + 1
        fin_excl = num_excl + num_nonb_excl(ix,ij)
        num_excl = fin_excl

        ini_nb14 = num_nb14 + 1
        fin_nb14 = num_nb14 + num_nb14_calc(ix,ij)
        num_nb14 = fin_nb14

        do k = ini_excl, fin_excl
          iy = nonb_excl_list(k,ij)
          idx = iy + (ix-1)*max_natom
          pairlist%univ_mask2(idx,univ_ij) = 0
        end do
        do k = ini_nb14, fin_nb14
          iy = nb14_calc_list(k,ij)
          idx = iy + (ix-1)*max_natom
          pairlist%univ_mask2(idx,univ_ij) = 0
        end do

      end do

    end do
    !$omp end parallel

    call gpu_alocate_packdata(pairlist%pack_univ_mask2, &
                              pairlist%univ_mask2_size, univ_ncell_near)
    call send_mask_data(pairlist%univ_mask2, pairlist%univ_mask2_size,  &
                        pairlist%pack_univ_mask2_size, univ_ncell_near, &
                        pairlist%pack_univ_mask2)

    call gpu_wait_build_pairlist(univ_ix_list, univ_iy_list, univ_ix_natom, &
                                 univ_iy_natom, pairlist%univ_mask2,        &
                                 pairlist%univ_mask2_size, univ_ncell_near )

    ! Sort pairlist
    !
    !$omp parallel default(shared) &
    !$omp private(univ_ij, ix_natom, iy_natom, base_idx, idx, id_load)
    !$omp do schedule(static,1)
    do univ_ij =1, univ_maxcell
       ix_natom = univ_ix_natom(univ_ij)
       iy_natom = univ_iy_natom(univ_ij)

       id_load = 16*(((ix_natom+7)/8) * ((iy_natom+3)/4)) + (ix_natom+7)/8;

       if ( univ_ij <= ncell ) then
          ! self cell-pair
          id_load = max_load - max_natom + (ix_natom + iy_natom)/2
       else if ( univ_ij <= univ_ncell_near) then
          ! near
          if ( id_load > max_load - max_natom ) then
             id_load = max_load - max_natom
          endif
       else
          ! far
          if ( id_load > max_load - max_natom - 1 ) then
             id_load = max_load - max_natom - 1
          endif
       endif
       pairlist%univ_ij_load(univ_ij) = id_load
    enddo
    !$omp end parallel

    count(:) = 0
    do univ_ij = 1, univ_maxcell
       id_load = pairlist%univ_ij_load(univ_ij)
       count(id_load) = count(id_load) + 1
    enddo

    pairlist%univ_ncell_nonzero = univ_maxcell - count(0)

    index(max_load) = 0
    do id_load = max_load-1, 0, -1
       index(id_load) = index(id_load+1) + count(id_load+1)
    enddo

    do univ_ij = 1, univ_maxcell
       id_load = pairlist%univ_ij_load(univ_ij)
       index(id_load) = index(id_load) + 1
       pairlist%univ_ij_sort_list( index(id_load) ) = univ_ij
    enddo

    deallocate( index, count )

    return

  end subroutine update_pairlist_pbc_univ

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    send_mask_data
  !> @brief        Compress and Send mask data
  !! @authors      KS
  !! @param[in]    univ_mask2           : mask data
  !! @param[in]    univ_mask2_size      : size of mask2 data
  !! @param[in]    pack_univ_mask2_size : size of compressed mask2 data
  !! @param[in]    univ_ncell_near      : univ_ncell_near
  !! @param[out]   pack_univ_mask2      : compressed mask2 data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine send_mask_data(univ_mask2, univ_mask2_size, &
                            pack_univ_mask2_size, univ_ncell_near, &
                            pack_univ_mask2)

    integer(kind=1), intent(in)    :: univ_mask2(*)
    integer,         intent(in)    :: univ_mask2_size
    integer,         intent(in)    :: pack_univ_mask2_size
    integer,         intent(in)    :: univ_ncell_near
    integer(kind=1), intent(out)   :: pack_univ_mask2(*)

    integer                        :: all_size
    integer                        :: divided_index


    all_size = univ_mask2_size*univ_ncell_near
    divided_index = int(all_size/2/8) * 8

    call pack_mask_data(univ_mask2, 1, divided_index, pack_univ_mask2)
    call gpu_copy_mask2(pack_univ_mask2, univ_mask2_size, univ_ncell_near, &
                        0, divided_index/8 - 1, 1)

    call pack_mask_data(univ_mask2, divided_index+1, all_size, pack_univ_mask2)
    call gpu_copy_mask2(pack_univ_mask2, univ_mask2_size, univ_ncell_near, &
                        divided_index/8, pack_univ_mask2_size - 1, 2)

    return

  end subroutine send_mask_data


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pack_mask_data
  !> @brief        Compress mask2 data
  !! @authors      KS
  !! @param[in]    mask2       : mask2 data
  !! @param[in]    start_index : Start index
  !! @param[in]    end_index   : End index
  !! @param[inout] pack_mask   : Compressed mask2 data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  subroutine pack_mask_data(mask2, start_index, end_index, pack_mask)

    ! formal arguments
    integer(kind=1), intent(in)    :: mask2(*)
    integer,         intent(in)    :: start_index
    integer,         intent(in)    :: end_index
    integer(kind=1), intent(inout) :: pack_mask(*)

    integer                        :: i,j,k
    integer                        :: masksize
    integer                        :: packblock_size
    integer                        :: mod_size
    integer                        :: packblock_offset

    integer(kind=1)                :: packed
    integer                        :: offset_index
    integer                        :: base_index
    integer(kind=1)                :: setval


    ! Size of data to pack
    masksize = end_index - start_index + 1
    ! Size of data after packing
    packblock_size = masksize/8
    packblock_offset = start_index/8

    ! Loop of packblock
    !$omp parallel do default(shared) private(packed,setval,offset_index)
    do k = 0, packblock_size - 1
      setval = 1
      packed = 0
      do j = 0, 7
        offset_index = k * 8 + j + start_index
        if (mask2(offset_index) == 1) packed = IOR(packed, setval)
        setval = ISHFT(setval, 1)
      end do
      pack_mask(packblock_offset + k + 1) = packed
    end do
    !$omp end parallel do

    ! Check if there is a remainder
    mod_size = mod(masksize, 8)
    if( mod_size > 0 ) then
      setval = 1
      packed = 0
      do j = 0, mod_size - 1
        offset_index = packblock_size * 8 + j + start_index
        if( mask2(offset_index) == 1 ) packed = IOR(packed, setval)
        setval = ISHFT(setval, 1)
      end do
      pack_mask(packblock_offset + packblock_size + 1) = packed
    end if

    return

  end subroutine pack_mask_data

#endif


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_nobc
  !> @brief        update pairlist in each domain without periodic boundary
  !!               condition
  !! @authors      JJ
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    domain   : domain information
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_nobc(enefunc, domain, pairlist)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_domain),   target, intent(in)    :: domain
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                  :: pairdist2
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: rtmp(1:3)

    integer                   :: i, j, ij, k, ix, iy
    integer                   :: num_nb15
    integer                   :: num_excl, ini_excl, fin_excl
    integer                   :: num_nb14, ini_nb14, fin_nb14
    logical                   :: nb15_calc
    integer                   :: id, omp_get_thread_num
    integer                   :: ncell, nboundary

    real(dp),         pointer :: coord(:,:,:)
    integer,          pointer :: natom(:), cell_pairlist1(:,:)
    integer,          pointer :: nb15_calc_list_nobc(:,:,:,:)
    integer,          pointer :: num_nb15_nobc(:,:)
    integer(1),       pointer :: exclusion_mask(:,:,:), exclusion_mask1(:,:,:)

    ncell           =  domain%num_cell_local
    nboundary       =  domain%num_cell_boundary
    natom           => domain%num_atom
    coord           => domain%coord
    cell_pairlist1  => domain%cell_pairlist1
    exclusion_mask  => enefunc%exclusion_mask
    exclusion_mask1 => enefunc%exclusion_mask1

    num_nb15_nobc   => pairlist%num_nb15_nobc
    nb15_calc_list_nobc => pairlist%nb15_calc_list_nobc

    pairdist2       =  pairlist%pairlistdist * pairlist%pairlistdist

    num_nb15_nobc(1:Maxatom,1:(ncell+nboundary)) = 0
    !$omp parallel                                                       &
    !$omp private(id, i, ix, iy, k, ij, j, num_excl, num_nb14, num_nb15, &
    !$omp         nb15_calc, rtmp, dij, rij2)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    ! make a pairlist in the same cell
    !
    do i = id+1, ncell, nthread

      do ix = 1, natom(i) - 1

        num_nb15 = 0

        do iy = ix + 1, natom(i)

          nb15_calc = .true.

          if (exclusion_mask1(iy,ix,i) /= 1) then

            ! store interaction table
            !
            num_nb15 = num_nb15 + 1
            nb15_calc_list_nobc(1,num_nb15,ix,i) = i
            nb15_calc_list_nobc(2,num_nb15,ix,i) = iy
          end if
        end do
        num_nb15_nobc(ix,i) = num_nb15

      end do
    end do

    ! Make a pairlist between different cells
    !
    do ij = 1, maxcell_near

      i = cell_pairlist1(1,ij)
      j = cell_pairlist1(2,ij)

      if (mod(i-1,nthread) == id) then

        do ix = 1, natom(i)

          rtmp(1:3) = coord(1:3,ix,i)
          num_nb15 = num_nb15_nobc(ix,i)

          do iy = 1, natom(j)

            if (exclusion_mask(iy,ix,ij) /= 1) then

              dij(1:3) = rtmp(1:3) - coord(1:3,iy,j)
              rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

              ! store interaction table
              !
              if (rij2 < pairdist2) then

                num_nb15 = num_nb15 + 1
                nb15_calc_list_nobc(1,num_nb15,ix,i) = j
                nb15_calc_list_nobc(2,num_nb15,ix,i) = iy
              end if
            end if
          end do

          num_nb15_nobc(ix,i) = num_nb15

        end do

      end if

    end do

    do ij = maxcell_near+1, maxcell

      i = cell_pairlist1(1,ij)
      j = cell_pairlist1(2,ij)

      if (mod(i-1,nthread) == id) then

        do ix = 1, natom(i)

          rtmp(1:3) = coord(1:3,ix,i)
          num_nb15  = num_nb15_nobc(ix,i)

          do iy = 1, natom(j)

            dij(1:3) = rtmp(1:3) - coord(1:3,iy,j)
            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            ! store interaction table
            !
            if (rij2 < pairdist2) then

              num_nb15 = num_nb15 + 1
              nb15_calc_list_nobc(1,num_nb15,ix,i) = j
              nb15_calc_list_nobc(2,num_nb15,ix,i) = iy
            end if
          end do

          num_nb15_nobc(ix,i) = num_nb15

        end do

      end if

    end do

    !$omp end parallel
    return

  end subroutine update_pairlist_nobc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_pairlist_fep
  !> @brief        initialize/allocate/setup pairlist for nonbonded interactions
  !>                for FEP calculations
  !! @authors      NK, HO
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    domain   : domain information
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_pairlist_fep(enefunc, domain, pairlist)

    ! formal arguments
    type(s_enefunc),         intent(in)    :: enefunc 
    type(s_domain),          intent(inout) :: domain 
    type(s_pairlist),        intent(inout) :: pairlist
    

    call init_pairlist(pairlist)

    pairlist%pairlistdist = enefunc%pairlistdist

    if (enefunc%forcefield == ForcefieldAAGO) then

      call alloc_pairlist(pairlist, PairListNOBC, &
        domain%num_cell_local+domain%num_cell_boundary)

      call timer(TimerPairList, TimerOn)
      call update_pairlist_nobc(enefunc, domain, pairlist)
      call timer(TimerPairList, TimerOff)

    else

#ifdef KCOMP
      if (enefunc%forcefield /= ForcefieldGROMARTINI) then
        call alloc_pairlist(pairlist, PairListNoTable_FEP, domain%num_cell_local)
      else
        call alloc_pairlist(pairlist, PairListNoTable_CG, domain%num_cell_local)
      endif
#else
      call alloc_pairlist(pairlist, PairListNoTable_FEP, domain%num_cell_local)
#endif

      call timer(TimerPairList, TimerOn)
#ifndef USE_GPU
      if (enefunc%pairlist_check) then
        call update_pairlist_pbc_check_fep(enefunc, domain, pairlist)
      else
        call update_pairlist_pbc_fep(enefunc, domain, pairlist)
      endif
#else
      call update_pairlist_pbc_univ(enefunc, domain, pairlist)
#endif
      call timer(TimerPairList, TimerOff)

    end if

    return

  end subroutine setup_pairlist_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_pbc_fep
  !> @brief        update pairlist in each domain for FEP calculations
  !!               with periodic boundary condition
  !! @authors      NK, HO
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    domain   : domain information
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_pbc_fep(enefunc, domain, pairlist)
  
    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_domain),   target, intent(in)    :: domain  
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                  :: pairdist2
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: rtmp(1:3)

    integer                   :: i, j, ij, k, ix, iy
    integer                   :: num_nb15_pre, num_nb15
    integer                   :: num_nb15_total, num_nb15_total1
    integer                   :: num_excl, ini_excl, fin_excl
    integer                   :: num_nb14, ini_nb14, fin_nb14
    logical                   :: nb15_calc
    integer                   :: id, omp_get_thread_num, num_nb15_cell

    real(dp),         pointer :: coord(:,:,:)
    real(wp),         pointer :: cell_move(:,:,:), system_size(:)
    real(wp),         pointer :: trans1(:,:,:), trans2(:,:,:)
    integer,          pointer :: num_nonb_excl1(:,:), num_nb14_calc1(:,:)
    integer,          pointer :: num_nonb_excl(:,:), num_nb14_calc(:,:)
    integer,          pointer :: nonb_excl_list1(:,:), nb14_calc_list1(:,:)
    integer,          pointer :: nonb_excl_list(:,:), nb14_calc_list(:,:)
    integer,          pointer :: natom(:), cell_pairlist1(:,:)
    integer,          pointer :: ncell, nboundary
    integer,          pointer :: nb15_cell(:), nb15_list(:,:)
    integer,          pointer :: nb15_calc_list1(:,:), nb15_calc_list(:,:)
    integer,          pointer :: num_nb15_calc1(:,:), num_nb15_calc(:,:)

    ! FEP
    integer                   :: num_nb15_pre_fep, num_nb15_fep
    integer                   :: num_nb15_cell_fep
    integer                   :: num_nb15_total_fep
    integer                   :: fg1, fg2
    integer,          pointer :: fepgrp(:,:)
    integer,          pointer :: nb15_cell_fep(:), nb15_list_fep(:,:)
    integer,          pointer :: nb15_calc_list1_fep(:,:), nb15_calc_list_fep(:,:)
    integer,          pointer :: num_nb15_calc1_fep(:,:), num_nb15_calc_fep(:,:)

    num_nonb_excl1  => enefunc%num_nonb_excl1
    num_nb14_calc1  => enefunc%num_nb14_calc1
    nonb_excl_list1 => enefunc%nonb_excl_list1
    nb14_calc_list1 => enefunc%nb14_calc_list1
    num_nonb_excl   => enefunc%num_nonb_excl
    num_nb14_calc   => enefunc%num_nb14_calc
    nonb_excl_list  => enefunc%nonb_excl_list
    nb14_calc_list  => enefunc%nb14_calc_list

    ncell           => domain%num_cell_local
    nboundary       => domain%num_cell_boundary
    natom           => domain%num_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    trans2          => domain%translated
    cell_pairlist1  => domain%cell_pairlist1
    cell_move       => domain%cell_move
    system_size     => domain%system_size

    num_nb15_calc1  => pairlist%num_nb15_calc1
    num_nb15_calc   => pairlist%num_nb15_calc
    nb15_list       => pairlist%nb15_list
    nb15_cell       => pairlist%nb15_cell
    nb15_calc_list1 => pairlist%nb15_calc_list1
    nb15_calc_list  => pairlist%nb15_calc_list

    pairdist2       =  pairlist%pairlistdist * pairlist%pairlistdist
    num_nb15_total  =  0

    ! FEP
    num_nb15_calc1_fep  => pairlist%num_nb15_calc1_fep
    num_nb15_calc_fep   => pairlist%num_nb15_calc_fep
    nb15_list_fep       => pairlist%nb15_list_fep
    nb15_cell_fep       => pairlist%nb15_cell_fep
    nb15_calc_list1_fep => pairlist%nb15_calc_list1_fep
    nb15_calc_list_fep  => pairlist%nb15_calc_list_fep
    fepgrp              => domain%fepgrp
    num_nb15_total_fep  =  0

    !$omp parallel                                                             &
    !$omp private(id, i, ix, ini_excl, num_excl, ini_nb14, num_nb14, num_nb15, &
    !$omp         num_nb15_pre, fin_excl, fin_nb14, iy, k, nb15_calc, ij, j,   &
    !$omp         rtmp, dij, rij2, num_nb15_cell,                              &
    !$omp         num_nb15_fep, num_nb15_pre_fep,                              &
    !$omp         num_nb15_cell_fep, fg1, fg2)                                 &
    !$omp reduction(+:num_nb15_total) reduction(+:num_nb15_total_fep)

#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell+nboundary, nthread
      do ix = 1, natom(i)
        trans2(1:3,ix,i) = coord(1:3,ix,i) + trans1(1:3,ix,i)
      end do
    end do

    !$omp barrier

    ! make a pairlist in the same cell
    !
    do i = id+1, ncell, nthread

      ini_excl = 0
      num_excl = 0
      ini_nb14 = 0
      num_nb14 = 0
      num_nb15 = 0
      num_nb15_pre = 0
      ! FEP: for perturbed interactions
      num_nb15_fep = 0
      num_nb15_pre_fep = 0
#ifdef DEBUG
      if (natom(i) > MaxAtom) &
        call error_msg( &
             'Debug: Update_Pairlist_Pbc> natom(cell) is exceed MaxAtom')
#endif

      do ix = 1, natom(i) - 1

        ini_excl = num_excl + 1
        fin_excl = num_excl + num_nonb_excl1(ix,i)
        num_excl = fin_excl

        ini_nb14 = num_nb14 + 1
        fin_nb14 = num_nb14 + num_nb14_calc1(ix,i)
        num_nb14 = fin_nb14
        do iy = ix + 1, natom(i)

          ! FEP: skip partA-partB interactions in FEP
          fg1 = fepgrp(ix,i)
          fg2 = fepgrp(iy,i)
          if (enefunc%fepgrp_nonb(fg1,fg2) == 0) cycle

          nb15_calc = .true.
          do k = ini_excl, fin_excl
            if (iy == nonb_excl_list1(k,i)) then
              nb15_calc = .false.
              exit
            end if
          end do

          do k = ini_nb14, fin_nb14
            if (iy == nb14_calc_list1(k,i)) then
              nb15_calc = .false.
              exit
            end if
          end do

          if (nb15_calc) then

            ! store interaction table
            !
            if (enefunc%fepgrp_nonb(fg1,fg2) == 5) then
              num_nb15 = num_nb15 + 1
#ifdef DEBUG
              if (num_nb15 > MaxNb15) &
                call error_msg('Debug: Update_Pairlist_Pbc> num_nb15 is exceed MaxNb15')
#endif
              nb15_calc_list1(num_nb15,i) = iy
            else
              ! FEP: for perturbed interactions
              num_nb15_fep = num_nb15_fep + 1
#ifdef DEBUG
              if (num_nb15 > MaxNb15) &
                call error_msg('Debug: Update_Pairlist_Pbc> num_nb15 is exceed MaxNb15')
#endif
              nb15_calc_list1_fep(num_nb15_fep,i) = iy
            end if
          end if
        end do

        num_nb15_calc1(ix,i) = num_nb15 - num_nb15_pre
        num_nb15_total = num_nb15_total + num_nb15 - num_nb15_pre
        num_nb15_pre = num_nb15
        ! FEP: for perturbed interactions
        num_nb15_calc1_fep(ix,i) = num_nb15_fep - num_nb15_pre_fep
        num_nb15_total_fep = num_nb15_total_fep + num_nb15_fep - num_nb15_pre_fep
        num_nb15_pre_fep = num_nb15_fep
      end do
    end do 

    ! Make a pairlist between different cells
    !
    do ij = id+1, maxcell_near, nthread

      i = cell_pairlist1(1,ij)
      j = cell_pairlist1(2,ij)

#ifdef DEBUG
      if (natom(i) > MaxAtom .or. natom(j) > MaxAtom) &
        call error_msg( &
             'Debug: Update_Pairlist_Pbc> natom(cell) is exceed MaxAtom')
#endif

      ini_excl = 0
      num_excl = 0
      ini_nb14 = 0
      num_nb14 = 0
      num_nb15 = 0
      num_nb15_pre = 0
      num_nb15_cell = 0
      ! FEP: for perturbed interactions
      num_nb15_fep = 0
      num_nb15_pre_fep = 0
      num_nb15_cell_fep = 0
      do ix = 1, natom(i)

        ini_excl = num_excl + 1
        fin_excl = num_excl + num_nonb_excl(ix,ij)
        num_excl = fin_excl

        ini_nb14 = num_nb14 + 1
        fin_nb14 = num_nb14 + num_nb14_calc(ix,ij)
        num_nb14 = fin_nb14

        rtmp(1:3) = trans2(1:3,ix,i) + cell_move(1:3,j,i)*system_size(1:3)
        do iy = 1, natom(j)

          ! FEP: skip partA-partB interactions in FEP
          fg1 = fepgrp(ix,i)
          fg2 = fepgrp(iy,j)
          if (enefunc%fepgrp_nonb(fg1,fg2) == 0) cycle

          nb15_calc = .true.

          do k = ini_excl, fin_excl
            if (iy == nonb_excl_list(k,ij)) then
              nb15_calc = .false.
              exit
            end if
          end do

          do k = ini_nb14, fin_nb14
            if (iy == nb14_calc_list(k,ij)) then
              nb15_calc = .false.
              exit
            end if
          end do

          if (nb15_calc) then

            dij(1:3) = rtmp(1:3) - trans2(1:3,iy,j) 

            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            ! store interaction table
            !
            if (rij2 < pairdist2) then
              if (enefunc%fepgrp_nonb(fg1,fg2) == 5) then
                num_nb15 = num_nb15 + 1
#ifdef DEBUG
                if (num_nb15 > MaxNb15) &
                  call error_msg( &
                       'Debug: Update_Pairlist_Pbc> num_nb15 is exceed MaxNb15')
#endif
                nb15_calc_list(num_nb15,ij) = iy
              else
                ! FEP: for perturbed interactions
                num_nb15_fep = num_nb15_fep + 1
#ifdef DEBUG
                if (num_nb15 > MaxNb15) &
                  call error_msg( &
                       'Debug: Update_Pairlist_Pbc> num_nb15 is exceed MaxNb15')
#endif
                nb15_calc_list_fep(num_nb15_fep,ij) = iy
              end if
            end if
          end if
        end do

        num_nb15_calc(ix,ij) = num_nb15 - num_nb15_pre
        if (num_nb15 /= num_nb15_pre) then
          num_nb15_cell = num_nb15_cell + 1
          nb15_list(num_nb15_cell,ij) = ix
        end if
        num_nb15_total = num_nb15_total + num_nb15 - num_nb15_pre
        num_nb15_pre = num_nb15

        ! FEP: for perturbed interactions
        num_nb15_calc_fep(ix,ij) = num_nb15_fep - num_nb15_pre_fep
        if (num_nb15_fep /= num_nb15_pre_fep) then
          num_nb15_cell_fep = num_nb15_cell_fep + 1
          nb15_list_fep(num_nb15_cell_fep,ij) = ix
        end if
        num_nb15_total_fep = num_nb15_total_fep + num_nb15_fep - num_nb15_pre_fep
        num_nb15_pre_fep = num_nb15_fep
      end do

      nb15_cell(ij) = num_nb15_cell
      ! FEP: for perturbed interactions
      nb15_cell_fep(ij) = num_nb15_cell_fep
    end do

    do ij = id+maxcell_near+1, maxcell, nthread

      i = cell_pairlist1(1,ij)
      j = cell_pairlist1(2,ij)

#ifdef DEBUG
      if (natom(i) > MaxAtom .or. natom(j) > MaxAtom) &
        call error_msg( &
             'Debug: Update_Pairlist_Pbc> natom(cell) is exceed MaxAtom')
#endif

      num_nb15 = 0
      num_nb15_pre = 0
      num_nb15_cell = 0
      ! FEP: for perturbed interactions
      num_nb15_fep = 0
      num_nb15_pre_fep = 0
      num_nb15_cell_fep = 0
      do ix = 1, natom(i)

        rtmp(1:3) = trans2(1:3,ix,i) + cell_move(1:3,j,i)*system_size(1:3)
        do iy = 1, natom(j)

          ! FEP: skip partA-partB interactions in FEP
          fg1 = fepgrp(ix,i)
          fg2 = fepgrp(iy,j)
          if (enefunc%fepgrp_nonb(fg1,fg2) == 0) cycle

          dij(1:3) = rtmp(1:3) - trans2(1:3,iy,j)
          rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          ! store interaction table
          !
          if (rij2 < pairdist2) then 
            if (enefunc%fepgrp_nonb(fg1,fg2) == 5) then
              num_nb15 = num_nb15 + 1
#ifdef DEBUG
              if (num_nb15 > MaxNb15) &
                call error_msg( &
                     'Debug: Update_Pairlist_Pbc> num_nb15 is exceed MaxNb15')
#endif
              nb15_calc_list(num_nb15,ij) = iy
            else
              ! FEP: for perturbed interactions
              num_nb15_fep = num_nb15_fep + 1
#ifdef DEBUG
              if (num_nb15 > MaxNb15) &
                call error_msg( &
                     'Debug: Update_Pairlist_Pbc> num_nb15 is exceed MaxNb15')
#endif
              nb15_calc_list_fep(num_nb15_fep,ij) = iy
            end if
          end if
        end do

        num_nb15_calc(ix,ij) = num_nb15 - num_nb15_pre
        if (num_nb15 /= num_nb15_pre) then
          num_nb15_cell = num_nb15_cell + 1
          nb15_list(num_nb15_cell,ij) = ix
        end if
        num_nb15_total = num_nb15_total + num_nb15 - num_nb15_pre
        num_nb15_pre = num_nb15

        ! FEP: for perturbed interactions
        num_nb15_calc_fep(ix,ij) = num_nb15_fep - num_nb15_pre_fep
        if (num_nb15_fep /= num_nb15_pre_fep) then
          num_nb15_cell_fep = num_nb15_cell_fep + 1
          nb15_list_fep(num_nb15_cell_fep,ij) = ix
        end if
        num_nb15_total_fep = num_nb15_total_fep + num_nb15_fep - num_nb15_pre_fep
        num_nb15_pre_fep = num_nb15_fep
      end do

      nb15_cell(ij) = num_nb15_cell
      ! FEP: for perturbed interactions
      nb15_cell_fep(ij) = num_nb15_cell_fep
    end do

    !$omp end parallel

#ifdef HAVE_MPI_GENESIS
    call mpi_reduce(num_nb15_total, num_nb15_total1, 1, mpi_integer, mpi_sum, &
         0, mpi_comm_country, ierror)
#endif

    return

  end subroutine update_pairlist_pbc_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_pairlist_pbc_check_fep
  !> @brief        update pairlist in each domain for FEP calculations
  !!               with periodic boundary condition
  !! @authors      NK, HO
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    domain   : domain information
  !! @param[inout] pairlist : pair-list information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_pairlist_pbc_check_fep(enefunc, domain, pairlist)
  
    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_domain),   target, intent(in)    :: domain  
    type(s_pairlist), target, intent(inout) :: pairlist

    ! local variables
    real(wp)                  :: pairdist2
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: rtmp(1:3)

    integer                   :: i, j, ij, k, ix, iy
    integer                   :: num_nb15_pre, num_nb15
    integer                   :: num_nb15_total, num_nb15_total1
    integer                   :: num_excl, ini_excl, fin_excl
    integer                   :: num_nb14, ini_nb14, fin_nb14
    integer                   :: small_contact
    logical                   :: nb15_calc
    integer                   :: id, omp_get_thread_num, num_nb15_cell

    real(dp),         pointer :: coord(:,:,:)
    real(wp),         pointer :: cell_move(:,:,:), system_size(:)
    real(wp),         pointer :: trans1(:,:,:), trans2(:,:,:)
    real(wp),         pointer :: err_minimum_contact
    integer,          pointer :: num_nonb_excl1(:,:), num_nb14_calc1(:,:)
    integer,          pointer :: num_nonb_excl(:,:), num_nb14_calc(:,:)
    integer,          pointer :: nonb_excl_list1(:,:), nb14_calc_list1(:,:)
    integer,          pointer :: nonb_excl_list(:,:), nb14_calc_list(:,:)
    integer,          pointer :: natom(:), cell_pairlist1(:,:)
    integer,          pointer :: ncell, nboundary
    integer,          pointer :: nb15_cell(:), nb15_list(:,:)
    integer,          pointer :: nb15_calc_list1(:,:), nb15_calc_list(:,:)
    integer,          pointer :: num_nb15_calc1(:,:), num_nb15_calc(:,:)
    integer,          pointer :: table_order
    logical,          pointer :: nonb_limiter, table

    ! FEP
    integer                   :: num_nb15_pre_fep, num_nb15_fep
    integer                   :: num_nb15_cell_fep
    integer                   :: num_nb15_total_fep
    integer                   :: fg1, fg2
    integer,          pointer :: fepgrp(:,:)
    integer,          pointer :: nb15_cell_fep(:), nb15_list_fep(:,:)
    integer,          pointer :: nb15_calc_list1_fep(:,:), nb15_calc_list_fep(:,:)
    integer,          pointer :: num_nb15_calc1_fep(:,:), num_nb15_calc_fep(:,:)

    num_nonb_excl1  => enefunc%num_nonb_excl1
    num_nb14_calc1  => enefunc%num_nb14_calc1
    nonb_excl_list1 => enefunc%nonb_excl_list1
    nb14_calc_list1 => enefunc%nb14_calc_list1
    num_nonb_excl   => enefunc%num_nonb_excl
    num_nb14_calc   => enefunc%num_nb14_calc
    nonb_excl_list  => enefunc%nonb_excl_list
    nb14_calc_list  => enefunc%nb14_calc_list

    table               => enefunc%table%table
    table_order         => enefunc%table%table_order
    nonb_limiter        => enefunc%nonb_limiter
    err_minimum_contact => enefunc%err_minimum_contact

    ncell           => domain%num_cell_local
    nboundary       => domain%num_cell_boundary
    natom           => domain%num_atom
    coord           => domain%coord
    trans1          => domain%trans_vec
    trans2          => domain%translated
    cell_pairlist1  => domain%cell_pairlist1
    cell_move       => domain%cell_move
    system_size     => domain%system_size

    num_nb15_calc1  => pairlist%num_nb15_calc1
    num_nb15_calc   => pairlist%num_nb15_calc
    nb15_list       => pairlist%nb15_list
    nb15_cell       => pairlist%nb15_cell
    nb15_calc_list1 => pairlist%nb15_calc_list1
    nb15_calc_list  => pairlist%nb15_calc_list

    pairdist2       =  pairlist%pairlistdist * pairlist%pairlistdist
    num_nb15_total  =  0
    small_contact   =  0

    ! FEP
    num_nb15_calc1_fep  => pairlist%num_nb15_calc1_fep
    num_nb15_calc_fep   => pairlist%num_nb15_calc_fep
    nb15_list_fep       => pairlist%nb15_list_fep
    nb15_cell_fep       => pairlist%nb15_cell_fep
    nb15_calc_list1_fep => pairlist%nb15_calc_list1_fep
    nb15_calc_list_fep  => pairlist%nb15_calc_list_fep
    fepgrp              => domain%fepgrp
    num_nb15_total_fep  =  0

    !$omp parallel                                                             &
    !$omp private(id, i, ix, ini_excl, num_excl, ini_nb14, num_nb14, num_nb15, &
    !$omp         num_nb15_pre, fin_excl, fin_nb14, iy, k, nb15_calc, ij, j,   &
    !$omp         rtmp, dij, rij2, num_nb15_cell,                              &
    !$omp         num_nb15_fep, num_nb15_pre_fep,                              &
    !$omp         num_nb15_cell_fep, fg1, fg2)                                 &
    !$omp reduction(+:num_nb15_total) reduction(+:small_contact)               &
    !$omp reduction(+:num_nb15_total_fep)

#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell+nboundary, nthread
      do ix = 1, natom(i)
        trans2(1:3,ix,i) = coord(1:3,ix,i) + trans1(1:3,ix,i)
      end do
    end do

    !$omp barrier

    ! make a pairlist in the same cell
    !
    do i = id+1, ncell, nthread

      ini_excl = 0
      num_excl = 0
      ini_nb14 = 0
      num_nb14 = 0
      num_nb15 = 0
      num_nb15_pre = 0
      ! FEP: for perturbed interactions
      num_nb15_fep = 0
      num_nb15_pre_fep = 0

      if (natom(i) > MaxAtom) &
        call error_msg( &
          'Debug: Update_Pairlist_Pbc_Check_Fep> natom(cell) is exceed MaxAtom')

      do ix = 1, natom(i) - 1

        ini_excl = num_excl + 1
        fin_excl = num_excl + num_nonb_excl1(ix,i)
        num_excl = fin_excl

        ini_nb14 = num_nb14 + 1
        fin_nb14 = num_nb14 + num_nb14_calc1(ix,i)
        num_nb14 = fin_nb14

        do iy = ix + 1, natom(i)

          ! FEP: skip partA-partB interactions in FEP
          fg1 = fepgrp(ix,i)
          fg2 = fepgrp(iy,i)
          if (enefunc%fepgrp_nonb(fg1,fg2) == 0) cycle

          nb15_calc = .true.
          do k = ini_excl, fin_excl
            if (iy == nonb_excl_list1(k,i)) then
              nb15_calc = .false.
              exit
            end if
          end do

          do k = ini_nb14, fin_nb14
            if (iy == nb14_calc_list1(k,i)) then
              nb15_calc = .false.
              exit
            end if
          end do

          if (nb15_calc) then

            ! store interaction table
            !
            if (enefunc%fepgrp_nonb(fg1,fg2) == 5) then
              num_nb15 = num_nb15 + 1
              if (table .and. table_order == 1 )then
                dij(1:3) = trans2(1:3,ix,i) - trans2(1:3,iy,i)
                rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
                if (rij2 < err_minimum_contact) then
                  if (.not. nonb_limiter) &
                    call error_msg( &
                    'Debug: Update_Pairlist_Pbc_Check_Fep> too small contact')
                  small_contact = small_contact + 1 
                endif
              endif
              if (num_nb15 > MaxNb15_chk) &
                call error_msg( &
                   'Debug: Update_Pairlist_Pbc_Check_Fep> num_nb15 is exceed MaxNb15')

              nb15_calc_list1(num_nb15,i) = iy
            else
              ! FEP: for perturbed interactions
              num_nb15_fep = num_nb15_fep + 1
              if (table .and. table_order == 1 )then
                dij(1:3) = trans2(1:3,ix,i) - trans2(1:3,iy,i)
                rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
                if (rij2 < err_minimum_contact) then
                  if (.not. nonb_limiter) &
                    call error_msg( &
                    'Debug: Update_Pairlist_Pbc_Check_Fep> too small contact')
                  small_contact = small_contact + 1 
                endif
              endif
              if (num_nb15_fep > MaxNb15_chk) &
                call error_msg( &
                   'Debug: Update_Pairlist_Pbc_Check_Fep> num_nb15 is exceed MaxNb15')

              nb15_calc_list1_fep(num_nb15_fep,i) = iy
            end if

          end if
        end do

        num_nb15_calc1(ix,i) = num_nb15 - num_nb15_pre
        num_nb15_total = num_nb15_total + num_nb15 - num_nb15_pre
        num_nb15_pre = num_nb15

        ! FEP: for perturbed interactions
        num_nb15_calc1_fep(ix,i) = num_nb15_fep - num_nb15_pre_fep
        num_nb15_total_fep = num_nb15_total_fep + num_nb15_fep - num_nb15_pre_fep
        num_nb15_pre_fep = num_nb15_fep

      end do
    end do 

    ! Make a pairlist between different cells
    !
    do ij = id+1, maxcell_near, nthread

      i = cell_pairlist1(1,ij)
      j = cell_pairlist1(2,ij)

      if (natom(i) > MaxAtom .or. natom(j) > MaxAtom) &
        call error_msg( &
             'Debug: Update_Pairlist_Pbc_Check_Fep> natom(cell) is exceed MaxAtom')

      ini_excl = 0
      num_excl = 0
      ini_nb14 = 0
      num_nb14 = 0
      num_nb15 = 0
      num_nb15_pre = 0
      num_nb15_cell = 0
      ! FEP: for perturbed interactions
      num_nb15_fep = 0
      num_nb15_pre_fep = 0
      num_nb15_cell_fep = 0

      do ix = 1, natom(i)

        ini_excl = num_excl + 1
        fin_excl = num_excl + num_nonb_excl(ix,ij)
        num_excl = fin_excl

        ini_nb14 = num_nb14 + 1
        fin_nb14 = num_nb14 + num_nb14_calc(ix,ij)
        num_nb14 = fin_nb14

        rtmp(1:3) = trans2(1:3,ix,i) + cell_move(1:3,j,i)*system_size(1:3)

        do iy = 1, natom(j)

          ! FEP: skip partA-partB interactions in FEP
          fg1 = fepgrp(ix,i)
          fg2 = fepgrp(iy,j)
          if (enefunc%fepgrp_nonb(fg1,fg2) == 0) cycle

          nb15_calc = .true.

          do k = ini_excl, fin_excl
            if (iy == nonb_excl_list(k,ij)) then
              nb15_calc = .false.
              exit
            end if
          end do

          do k = ini_nb14, fin_nb14
            if (iy == nb14_calc_list(k,ij)) then
              nb15_calc = .false.
              exit
            end if
          end do

          if (nb15_calc) then

            dij(1:3) = rtmp(1:3) - trans2(1:3,iy,j) 

            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            ! store interaction table
            !
            if (rij2 < pairdist2) then

              if (enefunc%fepgrp_nonb(fg1,fg2) == 5) then
                num_nb15 = num_nb15 + 1
                if (table .and. table_order == 1 )then
                  if (rij2 < err_minimum_contact) then
                    if (.not. nonb_limiter) &
                      call error_msg( &
                      'Debug: Update_Pairlist_Pbc_Check_Fep> too small contact')
                    small_contact = small_contact + 1 
                  endif
                endif
                if (num_nb15 > MaxNb15_chk) &
                  call error_msg( &
                   'Debug: Update_Pairlist_Pbc_Check_Fep> num_nb15 is exceed MaxNb15')

                nb15_calc_list(num_nb15,ij) = iy
              else
                ! FEP: for perturbed interactions
                num_nb15_fep = num_nb15_fep + 1
                if (table .and. table_order == 1 )then
                  if (rij2 < err_minimum_contact) then
                    if (.not. nonb_limiter) &
                      call error_msg( &
                      'Debug: Update_Pairlist_Pbc_Check_Fep> too small contact')
                    small_contact = small_contact + 1 
                  endif
                endif
                if (num_nb15_fep > MaxNb15_chk) &
                  call error_msg( &
                   'Debug: Update_Pairlist_Pbc_Check_Fep> num_nb15 is exceed MaxNb15')

                nb15_calc_list_fep(num_nb15_fep,ij) = iy
              end if
            end if
          end if
        end do

        num_nb15_calc(ix,ij) = num_nb15 - num_nb15_pre
        if (num_nb15 /= num_nb15_pre) then
          num_nb15_cell = num_nb15_cell + 1
          nb15_list(num_nb15_cell,ij) = ix
        end if
        num_nb15_total = num_nb15_total + num_nb15 - num_nb15_pre
        num_nb15_pre = num_nb15

        ! FEP: for perturbed interactions
        num_nb15_calc_fep(ix,ij) = num_nb15_fep - num_nb15_pre_fep
        if (num_nb15_fep /= num_nb15_pre_fep) then
          num_nb15_cell_fep = num_nb15_cell_fep + 1
          nb15_list_fep(num_nb15_cell_fep,ij) = ix
        end if
        num_nb15_total_fep = num_nb15_total_fep + num_nb15_fep - num_nb15_pre_fep
        num_nb15_pre_fep = num_nb15_fep
      end do

      nb15_cell(ij) = num_nb15_cell

      ! FEP: for perturbed interactions
      nb15_cell_fep(ij) = num_nb15_cell_fep

    end do

    do ij = id+maxcell_near+1, maxcell, nthread

      i = cell_pairlist1(1,ij)
      j = cell_pairlist1(2,ij)

      if (natom(i) > MaxAtom .or. natom(j) > MaxAtom) &
        call error_msg( &
           'Debug: Update_Pairlist_Pbc_Check_Fep> natom(cell) is exceed MaxAtom')

      num_nb15 = 0
      num_nb15_pre = 0
      num_nb15_cell = 0
      ! FEP: for perturbed interactions
      num_nb15_fep = 0
      num_nb15_pre_fep = 0
      num_nb15_cell_fep = 0

      do ix = 1, natom(i)

        rtmp(1:3) = trans2(1:3,ix,i) + cell_move(1:3,j,i)*system_size(1:3)

        do iy = 1, natom(j)

          ! FEP: skip partA-partB interactions in FEP
          fg1 = fepgrp(ix,i)
          fg2 = fepgrp(iy,j)
          if (enefunc%fepgrp_nonb(fg1,fg2) == 0) cycle

          nb15_calc = .true.

          if (nb15_calc) then

            dij(1:3) = rtmp(1:3) - trans2(1:3,iy,j)

            rij2 = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

            ! store interaction table
            !
            if (rij2 < pairdist2) then 

              if (enefunc%fepgrp_nonb(fg1,fg2) == 5) then
                num_nb15 = num_nb15 + 1

                if (num_nb15 > MaxNb15_chk) &
                  call error_msg( &
                   'Debug: Update_Pairlist_Pbc_Check_Fep> num_nb15 is exceed MaxNb15')

                nb15_calc_list(num_nb15,ij) = iy
              else
                ! FEP: for perturbed interactions
                num_nb15_fep = num_nb15_fep + 1

                if (num_nb15_fep > MaxNb15_chk) &
                  call error_msg( &
                   'Debug: Update_Pairlist_Pbc_Check_Fep> num_nb15 is exceed MaxNb15')

                nb15_calc_list_fep(num_nb15_fep,ij) = iy
              end if
            end if
          end if
        end do

        num_nb15_calc(ix,ij) = num_nb15 - num_nb15_pre
        if (num_nb15 /= num_nb15_pre) then
          num_nb15_cell = num_nb15_cell + 1
          nb15_list(num_nb15_cell,ij) = ix
        end if
        num_nb15_total = num_nb15_total + num_nb15 - num_nb15_pre
        num_nb15_pre = num_nb15

        ! FEP: for perturbed interactions
        num_nb15_calc_fep(ix,ij) = num_nb15_fep - num_nb15_pre_fep
        if (num_nb15_fep /= num_nb15_pre_fep) then
          num_nb15_cell_fep = num_nb15_cell_fep + 1
          nb15_list_fep(num_nb15_cell_fep,ij) = ix
        end if
        num_nb15_total_fep = num_nb15_total_fep + num_nb15_fep - num_nb15_pre_fep
        num_nb15_pre_fep = num_nb15_fep

      end do

      nb15_cell(ij) = num_nb15_cell
      ! FEP: for perturbed interactions
      nb15_cell_fep(ij) = num_nb15_cell_fep

    end do

    !$omp end parallel

    if (small_contact > 0) then
      write(MsgOut, *) "Warning: small contacts exist in inner loop"
    endif

#ifdef HAVE_MPI_GENESIS
    call mpi_reduce(num_nb15_total, num_nb15_total1, 1, mpi_integer, mpi_sum, &
         0, mpi_comm_country, ierror)
#endif

    return

  end subroutine update_pairlist_pbc_check_fep

end module sp_pairlist_mod
