!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_boundary_str_mod
!> @brief   structure of boundary
!! @authors Takaharu Mori (TM), Norio Takase (NT), Kiyoshi Yagi (KY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_boundary_str_mod

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

    logical                       :: use_cell_linked_list
    logical                       :: wrap_all
    integer                       :: num_cells_x
    integer                       :: num_cells_y
    integer                       :: num_cells_z
    real(wp)                      :: box_size_x
    real(wp)                      :: box_size_y
    real(wp)                      :: box_size_z
    real(wp)                      :: box_size_x_ref
    real(wp)                      :: box_size_y_ref
    real(wp)                      :: box_size_z_ref
    real(wp)                      :: cell_size_x
    real(wp)                      :: cell_size_y
    real(wp)                      :: cell_size_z
    real(wp)                      :: pairlist_grid

    integer                       :: num_cells
    integer, allocatable          :: neighbor_cells(:,:)
    integer, allocatable          :: cell_head_atom(:)
    integer, allocatable          :: cell_head_atomw(:)

    integer                       :: num_neighbor_cells
    integer, allocatable          :: neighbor_cell_common_x(:)
    integer, allocatable          :: neighbor_cell_common_y(:)
    integer, allocatable          :: neighbor_cell_common_z(:)

    ! spherical boundary condition
    logical                       :: sph_pot
    integer                       :: nfunctions
    real(wp), allocatable         :: radius(:)
    real(wp), allocatable         :: center(:,:)
    real(wp), allocatable         :: const(:)
    integer, allocatable          :: exponent(:)
    real(wp)                      :: fix_layer
    integer                       :: num_fixatm
    logical, allocatable          :: fixatm(:)

    ! currently not in use
    integer                       :: num_nospot
    integer, allocatable          :: nospotlist(:)
    integer, allocatable          :: atomlist(:)

  end type s_boundary

  ! parameters for allocatable variables
  integer,      public, parameter :: BoundaryCells        = 1
  integer,      public, parameter :: BoundarySphericalPot = 2

  ! parameters
  integer,      public, parameter :: BoundaryTypeNOBC     = 1
  integer,      public, parameter :: BoundaryTypePBC      = 2

  character(*), public, parameter :: BoundaryTypeTypes(2) = (/'NOBC', &
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
  !! @authors      NT
  !! @param[out]   boundary : boundary conditions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_boundary(boundary)

    ! formal arguments
    type(s_boundary),        intent(inout) :: boundary


    boundary%type                 = 0
    boundary%origin_x             = 0.0_wp
    boundary%origin_y             = 0.0_wp
    boundary%origin_z             = 0.0_wp
    boundary%use_cell_linked_list = .false.
    boundary%wrap_all             = .false.
    boundary%num_cells_x          = 0
    boundary%num_cells_y          = 0
    boundary%num_cells_z          = 0
    boundary%box_size_x           = 0.0_wp
    boundary%box_size_y           = 0.0_wp
    boundary%box_size_z           = 0.0_wp
    boundary%box_size_x_ref       = 0.0_wp
    boundary%box_size_y_ref       = 0.0_wp
    boundary%box_size_z_ref       = 0.0_wp
    boundary%cell_size_x          = 0.0_wp
    boundary%cell_size_y          = 0.0_wp
    boundary%cell_size_z          = 0.0_wp
    boundary%num_cells            = 0
    boundary%pairlist_grid        = 0.0_wp
    boundary%sph_pot              = .false.
    boundary%nfunctions           = 0
    boundary%num_nospot           = 0
    boundary%fix_layer            = 0.0_wp

    return

  end subroutine init_boundary

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_boundary
  !> @brief        allocate boundary conditions information
  !! @authors      TM, KY
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

      if (allocated(boundary%neighbor_cell_common_x)) then

        if (size(boundary%neighbor_cell_common_x) <        &
            boundary%num_neighbor_cells ) then
          deallocate(boundary%neighbor_cells,              &
                     boundary%cell_head_atom,              &
                     boundary%cell_head_atomw,             &
                     boundary%neighbor_cell_common_x,      &
                     boundary%neighbor_cell_common_y,      &
                     boundary%neighbor_cell_common_z,      &
                     stat = dealloc_stat)
        endif
      endif

      if (allocated(boundary%neighbor_cells)) then

        if (size(boundary%neighbor_cells(1,:)) < var_size) then
          deallocate(boundary%neighbor_cells,  &
                     boundary%cell_head_atom,  &
                     boundary%cell_head_atomw, &
                     stat = dealloc_stat)
           if (dealloc_stat /= 0) call error_msg_dealloc
        endif
      end if

      if (.not. allocated(boundary%neighbor_cells)) then

        allocate(boundary%neighbor_cells(boundary%num_neighbor_cells,var_size),&
                 boundary%cell_head_atom(var_size),                            &
                 boundary%cell_head_atomw(var_size),                           &
                 stat = alloc_stat)
        if (alloc_stat /= 0)   call error_msg_alloc
      endif

      if (.not. allocated(boundary%neighbor_cell_common_x)) then

        allocate(boundary%neighbor_cell_common_x(boundary%num_neighbor_cells), &
                 boundary%neighbor_cell_common_y(boundary%num_neighbor_cells), &
                 boundary%neighbor_cell_common_z(boundary%num_neighbor_cells), &
                 stat = alloc_stat)
      endif

      boundary%neighbor_cells(1:boundary%num_neighbor_cells,1:var_size) = 0
      boundary%neighbor_cell_common_x(1:boundary%num_neighbor_cells)    = 0
      boundary%neighbor_cell_common_y(1:boundary%num_neighbor_cells)    = 0
      boundary%neighbor_cell_common_z(1:boundary%num_neighbor_cells)    = 0
      boundary%cell_head_atom (1:var_size) = 0
      boundary%cell_head_atomw(1:var_size) = 0

    case (BoundarySphericalPot)

      if (allocated(boundary%radius)) then
        if (size(boundary%radius) == var_size) return
        
        deallocate(boundary%radius, &
                   boundary%center, &
                   boundary%const,  &
                   boundary%exponent, &
                   stat = dealloc_stat)
         if (dealloc_stat /= 0) call error_msg_dealloc
      end if

      allocate(boundary%radius(var_size),    &
               boundary%center(3, var_size), &
               boundary%const(var_size),     &
               boundary%exponent(var_size),  &
               stat = alloc_stat)

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
  !! @authors      TM, KY
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
        deallocate (boundary%neighbor_cells,  &
                    boundary%neighbor_cell_common_x, &
                    boundary%neighbor_cell_common_y, &
                    boundary%neighbor_cell_common_z, &
                    boundary%cell_head_atom,  &
                    boundary%cell_head_atomw, &
                    stat = dealloc_stat)
      end if

    case (BoundarySphericalPot)

      if (allocated(boundary%radius)) then
        deallocate(boundary%radius, &
                   boundary%center, &
                   boundary%const,  &
                   boundary%exponent, &
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
  !! @authors      TM, KY
  !! @param[inout] boundary : boundary conditions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_boundary_all(boundary)

    ! formal arguments
    type(s_boundary),        intent(inout) :: boundary


    call dealloc_boundary(boundary, BoundaryCells)
    call dealloc_boundary(boundary, BoundarySphericalPot)

    return

  end subroutine dealloc_boundary_all

end module at_boundary_str_mod
