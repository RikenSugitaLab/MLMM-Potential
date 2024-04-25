!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_boundary_str_mod
!> @brief   structure of boundary
!! @authors Jaewoon Jung (JJ), Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_boundary_str_mod

  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_boundary

    integer                       :: type
    real(wp)                      :: origin_x
    real(wp)                      :: origin_y
    real(wp)                      :: origin_z

    integer                       :: num_cells_x
    integer                       :: num_cells_y
    integer                       :: num_cells_z
    real(dp)                      :: box_size_x
    real(dp)                      :: box_size_y
    real(dp)                      :: box_size_z
    real(dp)                      :: box_size_x_ref
    real(dp)                      :: box_size_y_ref
    real(dp)                      :: box_size_z_ref
    real(dp)                      :: cell_size_x
    real(dp)                      :: cell_size_y
    real(dp)                      :: cell_size_z
    real(dp)                      :: box_size_x_max
    real(dp)                      :: box_size_y_max
    real(dp)                      :: box_size_z_max
    real(dp)                      :: box_size_x_min
    real(dp)                      :: box_size_y_min
    real(dp)                      :: box_size_z_min

    integer                       :: num_domain(3)

    !  BoundaryCells
    integer,          allocatable :: neighbor_cells(:,:)

  end type s_boundary

  ! parameters for allocatable variables
  integer,      public, parameter :: BoundaryCells        = 1

  ! parameters
  integer,      public, parameter :: BoundaryTypeNOBC     = 1
  integer,      public, parameter :: BoundaryTypePBC      = 2

  character(*), public, parameter :: BoundaryTypeTypes(2) = (/'NOBC',&
                                                              'PBC '/)

  ! subroutines
  public :: init_boundary
  public :: alloc_boundary
  public :: dealloc_boundary
  public :: dealloc_boundary_all

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_boundary
  !> @brief        initialize boundary conditions information
  !! @authors      JJ
  !! @param[out]   boundary : boundary conditions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_boundary(boundary)

    ! formal arguments
    type(s_boundary),        intent(inout) :: boundary


    boundary%type            = 0
    boundary%origin_x        = 0.0_wp
    boundary%origin_y        = 0.0_wp
    boundary%origin_z        = 0.0_wp
    boundary%num_cells_x     = 0
    boundary%num_cells_y     = 0
    boundary%num_cells_z     = 0
    boundary%box_size_x      = 0.0_wp
    boundary%box_size_y      = 0.0_wp
    boundary%box_size_z      = 0.0_wp
    boundary%box_size_x_ref  = 0.0_wp
    boundary%box_size_y_ref  = 0.0_wp
    boundary%box_size_z_ref  = 0.0_wp
    boundary%cell_size_x     = 0.0_wp
    boundary%cell_size_y     = 0.0_wp
    boundary%cell_size_z     = 0.0_wp

    return

  end subroutine init_boundary

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_boundary
  !> @brief        allocate boundary with domain decomposition (half size cell)
  !! @authors      JJ
  !! @param[inout] boundary : boundary conditions information
  !! @param[in]    variable : selected variable
  !! @param[in]    var_size : size of the selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_boundary(boundary, variable, var_size)

    ! formal arguments
    type(s_boundary),        intent(inout) :: boundary
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size

    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat


    alloc_stat   = 0
    dealloc_stat = 0

    ! allocate selected variables
    !
    select case (variable)

    case(BoundaryCells)

      if (allocated(boundary%neighbor_cells)) then
        if (size(boundary%neighbor_cells(1,:)) /= var_size) &
          deallocate(boundary%neighbor_cells, stat = dealloc_stat)
      end if

      if (.not. allocated(boundary%neighbor_cells)) &
        allocate(boundary%neighbor_cells(125,var_size), stat = alloc_stat)

      boundary%neighbor_cells(1:125,1:var_size) = 0

    case default

      call error_msg('Alloc_Boundary> bad variable')

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_boundary

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_boundary
  !> @brief        deallocate boundary conditions information
  !! @authors      TM
  !! @param[inout] boundary : boundary conditions information
  !! @param[in]    variable : selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_boundary(boundary, variable)

    ! formal arguments
    type(s_boundary),        intent(inout) :: boundary
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat


    alloc_stat   = 0
    dealloc_stat = 0

    ! deallocate selected variables
    !
    select case (variable)

    case(BoundaryCells)

      if (allocated(boundary%neighbor_cells)) then
        deallocate (boundary%neighbor_cells, &
                    stat = dealloc_stat)
      end if

    case default

      call error_msg('Dealloc_Boundary> bad variable')

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_boundary

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_boundary_all
  !> @brief        deallocate all boundary conditions information
  !! @authors      TM
  !! @param[inout] boundary : boundary conditions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_boundary_all(boundary)

    ! formal arguments
    type(s_boundary),        intent(inout) :: boundary


    call dealloc_boundary(boundary, BoundaryCells)

    return

  end subroutine dealloc_boundary_all

end module sp_boundary_str_mod
