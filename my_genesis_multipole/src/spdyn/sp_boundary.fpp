!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_boundary_mod
!> @brief   utilities for boundary conditions
!! @authors Jaewoon Jung (JJ), Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_boundary_mod

  use sp_ensemble_str_mod
  use sp_boundary_str_mod
  use molecules_str_mod
  use fileio_rst_mod
  use fileio_control_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_pbc_info
    real(wp)            :: box_size_x    = 0.0_wp
    real(wp)            :: box_size_y    = 0.0_wp
    real(wp)            :: box_size_z    = 0.0_wp
    integer             :: num_cells_x   = 0
    integer             :: num_cells_y   = 0
    integer             :: num_cells_z   = 0
  end type s_pbc_info

  type, public :: s_boundary_info
    integer             :: type          = BoundaryTypePBC
    type(s_pbc_info)    :: pbc_info
    real(wp)            :: origin_x      = 0.0_wp
    real(wp)            :: origin_y      = 0.0_wp
    real(wp)            :: origin_z      = 0.0_wp
    integer             :: domain_x      = 0
    integer             :: domain_y      = 0
    integer             :: domain_z      = 0
    real(wp)            :: box_size_x_max = 100000000.0_wp
    real(wp)            :: box_size_y_max = 100000000.0_wp
    real(wp)            :: box_size_z_max = 100000000.0_wp
    real(wp)            :: box_size_x_min = 0.0_wp
    real(wp)            :: box_size_y_min = 0.0_wp
    real(wp)            :: box_size_z_min = 0.0_wp
  end type s_boundary_info

  ! subroutines
  public  :: show_ctrl_boundary
  public  :: read_ctrl_boundary
  public  :: setup_boundary
  public  :: setup_boundary_pio
  public  :: setup_processor_number
  public  :: setup_boundary_cell

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_boundary
  !> @brief        show BOUNDARY section usage
  !! @authors      NT
  !! @param[in]    show_all : show all usage or not
  !! @param[in]    run_mode : run mode string : "md", "min"
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_boundary(show_all, run_mode)

    ! formal arguments
    logical,                 intent(in)    :: show_all
    character(*),            intent(in)    :: run_mode


    if (show_all) then

      select case (run_mode)

      case ('md', 'min', 'remd')

        write(MsgOut,'(A)') '[BOUNDARY]'
        write(MsgOut,'(A)') 'type          = PBC       # [PBC, NOBC]'
        write(MsgOut,'(A)') '# box_size_x    = 0.0       # box size (x) in [PBC, NOBC]'
        write(MsgOut,'(A)') '# box_size_y    = 0.0       # box size (y) in [PBC, NOBC]'
        write(MsgOut,'(A)') '# box_size_z    = 0.0       # box size (z) in [PBC, NOBC]'
        write(MsgOut,'(A)') 'domain_x      = 0         # domain size (x)'
        write(MsgOut,'(A)') 'domain_y      = 0         # domain size (y)'
        write(MsgOut,'(A)') 'domain_z      = 0         # domain size (z)'
        write(MsgOut,'(A)') ' '

      end select

    else

      select case (run_mode)

      case ('md', 'min', 'remd')

        write(MsgOut,'(A)') '[BOUNDARY]'
        write(MsgOut,'(A)') 'type          = PBC       # [PBC, NOBC]'
        write(MsgOut,'(A)') ' '

      end select

    end if


    return

  end subroutine show_ctrl_boundary
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_boundary
  !> @brief        read BOUNDARY section in the control file
  !! @authors      JJ
  !! @param[in]    handle     : unit number
  !! @param[out]   bound_info : BOUNDARY section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_boundary(handle, bound_info)

    ! parameters
    character(*),            parameter     :: Section = 'Boundary'

    ! formal arguments
    integer,                 intent(in)    :: handle
    type(s_boundary_info),   intent(inout) :: bound_info


    ! read parameters from control file
    ! 
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_type(handle, Section, 'type', &
                            bound_info%type, BoundaryTypeTypes)


    select case (bound_info%type)

    case (BoundaryTypePBC)
      call read_ctrlfile_real (handle, Section, 'box_size_x',  &
                               bound_info%pbc_info%box_size_x)
      call read_ctrlfile_real (handle, Section, 'box_size_y',  &
                               bound_info%pbc_info%box_size_y)
      call read_ctrlfile_real (handle, Section, 'box_size_z',  &
                               bound_info%pbc_info%box_size_z)

    case (BoundaryTypeNOBC)
      call read_ctrlfile_real (handle, Section, 'box_size_x_max',  &
                               bound_info%box_size_x_max)
      call read_ctrlfile_real (handle, Section, 'box_size_y_max',  &
                               bound_info%box_size_y_max)
      call read_ctrlfile_real (handle, Section, 'box_size_z_max',  &
                               bound_info%box_size_z_max)
      call read_ctrlfile_real (handle, Section, 'box_size_x_min',  &
                               bound_info%box_size_x_min)
      call read_ctrlfile_real (handle, Section, 'box_size_y_min',  &
                               bound_info%box_size_y_min)
      call read_ctrlfile_real (handle, Section, 'box_size_z_min',  &
                               bound_info%box_size_z_min)
    end select

    call read_ctrlfile_integer(handle, Section, 'domain_x', &
                               bound_info%domain_x)
    call read_ctrlfile_integer(handle, Section, 'domain_y', &
                               bound_info%domain_y)
    call read_ctrlfile_integer(handle, Section, 'domain_z', &
                               bound_info%domain_z)

    call end_ctrlfile_section(handle)

    ! write parameters to MsgOut
    !
    if (main_rank) then

      write(MsgOut,'(A)') 'Read_Ctrl_Boundary> Parameters of Boundary Condition'
      write(MsgOut,'(A20,A10)') &
            '  type            = ', BoundaryTypeTypes(bound_info%type)

      select case (bound_info%type)

      case (BoundaryTypePBC)

        write(MsgOut,'(A20,3F10.3)')                                &
            '  box_size(x,y,z) = ', bound_info%pbc_info%box_size_x, &
                                    bound_info%pbc_info%box_size_y, &
                                    bound_info%pbc_info%box_size_z

      case (BoundaryTypeNOBC)

        write(MsgOut,'(A24,3F10.3)')                                    &
            '  box_size_max(x,y,z) = ', bound_info%box_size_x_max, &
                                        bound_info%box_size_y_max, &
                                        bound_info%box_size_z_max
        write(MsgOut,'(A20,3F10.3)')                                    &
            '  box_size_max(x,y,z) = ', bound_info%box_size_x_min, &
                                        bound_info%box_size_y_min, &
                                        bound_info%box_size_z_min
      end select

      if (bound_info%domain_x /= 0 .and. &
          bound_info%domain_y /= 0 .and. &
          bound_info%domain_z /= 0) &
        write(MsgOut,'(A20,3I10)')                         &
              '  domain (x,y,z)  = ', bound_info%domain_x, &
                                      bound_info%domain_y, &
                                      bound_info%domain_z

      write(MsgOut,'(A)') ' '

    end if

    return

  end subroutine read_ctrl_boundary

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_boundary
  !> @brief        set essential variables for boundary condition
  !! @authors      JJ
  !! @param[in]    bound_info  : BOUNDARY section control parameters information
  !! @param[in]    table        : flag for use table or not
  !! @param[in]    pairlistdist : pair-list distance
  !! @pmara[in]    water_model  : water model
  !! @param[in]    ensemble     : type of ensemble 
  !! @param[in]    rigid_bond   : flag for rigid-bond
  !! @param[in]    rst          : restart file information
  !! @param[in]    dsize_cg     : flag for reset domain size for CG-model
  !! @param[in]    dmin_size_cg : minimum domain size for CG-model
  !! @param[inout] boundary     : boundary information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_boundary(bound_info, table, pairlistdist, water_model, &
                            ensemble, rigid_bond, dsize_cg, dmin_size_cg, &
                            molecule, rst, boundary)

    ! formal arguments
    type(s_boundary_info),   intent(in)    :: bound_info
    logical,                 intent(in)    :: table
    real(wp),                intent(in)    :: pairlistdist
    character(*),            intent(in)    :: water_model
    integer,                 intent(in)    :: ensemble
    logical,                 intent(in)    :: rigid_bond
    logical,                 intent(in)    :: dsize_cg
    real(wp),                intent(in)    :: dmin_size_cg
    type(s_rst),             intent(in)    :: rst
    type(s_molecule),        intent(in)    :: molecule
    type(s_boundary),        intent(inout) :: boundary

    ! local variables
    integer                  :: i
    real(dp)                 :: coord_min(3), coord_max(3), box_size(3)
    real(dp)                 :: center(3)

    call init_boundary(boundary)

    if (.not.((bound_info%domain_x == 2 .and. &
               bound_info%domain_y == 1 .and. &
               bound_info%domain_z == 1) .or. &
              (bound_info%domain_x == 2 .and. &
               bound_info%domain_y == 2 .and. &
               bound_info%domain_z == 1) .or. &
              (bound_info%domain_x == 1 .and. &
               bound_info%domain_y == 1 .and. &
               bound_info%domain_z == 1))) then

      if (bound_info%domain_x == 1 .or. &
          bound_info%domain_y == 1 .or. &
          bound_info%domain_z == 1 ) then
        call error_msg('Setup_Boundary> other than (2,1,1)/(2,2,1)/(1,1,1), &
                       &domain[x,y,z] should be larger than 1')
      end if

    end if

    select case (bound_info%type)

    case (BoundaryTypeNOBC)

      boundary%type           = bound_info%type
      boundary%origin_x       = bound_info%origin_x 
      boundary%origin_y       = bound_info%origin_y 
      boundary%origin_z       = bound_info%origin_z
      boundary%box_size_x_max = bound_info%box_size_x_max
      boundary%box_size_y_max = bound_info%box_size_y_max
      boundary%box_size_z_max = bound_info%box_size_z_max
      boundary%box_size_x_min = bound_info%box_size_x_min
      boundary%box_size_y_min = bound_info%box_size_y_min
      boundary%box_size_z_min = bound_info%box_size_z_min


      ! Decidie the system size
      !
      coord_min(1:3)          =  1000000000000.0_dp 
      coord_max(1:3)          = -1000000000000.0_dp 
      if (rst%rstfile_type /= RstfileTypeUndef) then
        do i = 1, rst%num_atoms
          coord_min(1:3) = min(coord_min(1:3), rst%coord(1:3,i))
          coord_max(1:3) = max(coord_max(1:3), rst%coord(1:3,i))
        end do
      else
        do i = 1, molecule%num_atoms
          coord_min(1:3) = min(coord_min(1:3), molecule%atom_coord(1:3,i))
          coord_max(1:3) = max(coord_max(1:3), molecule%atom_coord(1:3,i))
        end do
      end if
      box_size(1:3) = max(-coord_min(1:3)+0.1_wp,coord_max(1:3)+0.1_wp)
      boundary%box_size_x     = box_size(1)*2.0_wp
      boundary%box_size_y     = box_size(2)*2.0_wp
      boundary%box_size_z     = box_size(3)*2.0_wp
      if (boundary%box_size_x > boundary%box_size_x_max) &
        call error_msg('calculated box_size_x is greater than box_size_x_max')
      if (boundary%box_size_y > boundary%box_size_y_max) &
        call error_msg('calculated box_size_y is greater than box_size_y_max')
      if (boundary%box_size_z > boundary%box_size_z_max) &
        call error_msg('calculated box_size_z is greater than box_size_z_max')
      if (boundary%box_size_x < boundary%box_size_x_min) &
        write(Msgout, '(A)') 'WARNING : calculated box size_x is less than box_size_x_min'
      if (boundary%box_size_y < boundary%box_size_y_min) &
        write(Msgout, '(A)') 'WARNING : calculated box size_y is less than box_size_y_min'
      if (boundary%box_size_z < boundary%box_size_z_min) &
        write(Msgout, '(A)') 'WARNING : calculated box size_z is less than box_size_z_min'

      boundary%box_size_x_ref = boundary%box_size_x
      boundary%box_size_y_ref = boundary%box_size_y
      boundary%box_size_z_ref = boundary%box_size_z

    case (BoundaryTypePBC)

      boundary%type           = bound_info%type
      boundary%origin_x       = bound_info%origin_x
      boundary%origin_y       = bound_info%origin_y
      boundary%origin_z       = bound_info%origin_z
      boundary%box_size_x     = bound_info%pbc_info%box_size_x
      boundary%box_size_y     = bound_info%pbc_info%box_size_y
      boundary%box_size_z     = bound_info%pbc_info%box_size_z
      boundary%box_size_x_ref = boundary%box_size_x
      boundary%box_size_y_ref = boundary%box_size_y
      boundary%box_size_z_ref = boundary%box_size_z
    
      if (rst%rstfile_type /= RstfileTypeUndef) then

        boundary%box_size_x     = rst%box_size_x
        boundary%box_size_y     = rst%box_size_y
        boundary%box_size_z     = rst%box_size_z
        boundary%box_size_x_ref = boundary%box_size_x
        boundary%box_size_y_ref = boundary%box_size_y
        boundary%box_size_z_ref = boundary%box_size_z

      end if

    end select

    call setup_processor_number(bound_info, &
                                table, pairlistdist, water_model, &
                                ensemble, rigid_bond, boundary)

    call setup_boundary_cell   (table, pairlistdist, water_model, &
                                ensemble, rigid_bond, dsize_cg, dmin_size_cg,  &
                                boundary)

    return

  end subroutine setup_boundary

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_boundary_pio
  !> @brief        set essential variables for boundary condition
  !! @authors      NT
  !! @param[in]    bound_info  : BOUNDARY section control parameters information
  !! @param[in]    table        : flag for use table or not
  !! @param[in]    pairlistdist : pair-list distance
  !! @pmara[in]    water_model  : water model
  !! @param[in]    ensemble     : type of ensemble 
  !! @param[in]    rigid_bond   : flag for rigid-bond
  !! @param[in]    dsize_cg     : flag for reset domain size for CG-model
  !! @param[in]    dmin_size_cg : minimum domain size for CG-model
  !! @param[inout] boundary     : boundary information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_boundary_pio(bound_info, table, pairlistdist, water_model, &
                                ensemble, rigid_bond, dsize_cg, dmin_size_cg, &
                                boundary)

    ! formal arguments
    type(s_boundary_info),   intent(in)    :: bound_info
    logical,                 intent(in)    :: table
    real(wp),                intent(in)    :: pairlistdist
    character(*),            intent(in)    :: water_model
    integer,                 intent(in)    :: ensemble
    logical,                 intent(in)    :: rigid_bond
    logical,                 intent(in)    :: dsize_cg
    real(wp),                intent(in)    :: dmin_size_cg
    type(s_boundary),        intent(inout) :: boundary


    if (.not.((bound_info%domain_x == 2 .and. &
               bound_info%domain_y == 1 .and. &
               bound_info%domain_z == 1) .or. &
              (bound_info%domain_x == 2 .and. &
               bound_info%domain_y == 2 .and. &
               bound_info%domain_z == 1))) then

      if (bound_info%domain_x == 1 .or. &
          bound_info%domain_y == 1 .or. &
          bound_info%domain_z == 1 ) then
        call error_msg('Setup_Boundary_Pio> other than (2,1,1)/(2,2,1)/(1,1,1), &
                       &domain[x,y,z] should be larger than 1')
      end if

    end if

    boundary%type = bound_info%type

    call setup_processor_number(bound_info, table, pairlistdist, water_model, &
                                ensemble, rigid_bond, boundary)

    call setup_boundary_cell(table, pairlistdist, water_model,              &
                             ensemble, rigid_bond, dsize_cg, dmin_size_cg,  &
                             boundary)

    return

  end subroutine setup_boundary_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_processor_number
  !> @brief        define the processor number in each dimension
  !! @authors      JJ
  !! @param[in]    bound_info : BOUNDARY section control parameters information
  !! @param[in]    table        : flag for use table or not
  !! @param[in]    pairlistdist : pair-list distance
  !! @pmara[in]    water_model  : water model
  !! @param[in]    ensemble     : type of ensemble 
  !! @param[in]    rigid_bond   : flag for rigid-bond
  !! @param[inout] boundary     : boundary information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_processor_number(bound_info, table, pairlistdist, &
                                    water_model, ensemble, rigid_bond, boundary)

    ! formal arguments
    type(s_boundary_info),   intent(in)    :: bound_info
    logical,                 intent(in)    :: table
    real(wp),                intent(in)    :: pairlistdist
    character(*),            intent(in)    :: water_model
    integer,                 intent(in)    :: ensemble
    logical,                 intent(in)    :: rigid_bond
    type(s_boundary),        intent(inout) :: boundary

    ! local variable
    real(wp)                 :: bsize_x, bsize_y, bsize_z
    real(wp)                 :: size_x, size_y, size_z, maxsize(0:100)
    integer                  :: total_proc
    integer                  :: nx, ny, nz, nx1, ny1,nz1, i, j, k, itype
    integer                  :: nc(3,100), cell_size(3,100)
    logical                  :: extend, extend1
    logical                  :: defined_proc

    extend = .false.
    if ( table  .and. &
        (water_model(1:4) == 'TIP3' .or. &
         water_model(1:3) == 'WAT'  .or. &
         water_model(1:3) == 'SOL') )    & 
       extend = .true.

    if (rigid_bond) &
      extend = .true.

    extend1 = .false.
    if (ensemble == EnsembleNPT  .or. &
        ensemble == EnsembleNPAT .or. &
        ensemble == EnsembleNPgT ) &
      extend1 = .true.

    ! check processor number based on the num of domains
    !
    if (bound_info%domain_x /= 0 .and. &
        bound_info%domain_y /= 0 .and. &
        bound_info%domain_z /= 0) then

      total_proc = bound_info%domain_x * &
                   bound_info%domain_y * &
                   bound_info%domain_z

      if (total_proc == nproc_country) then

        boundary%num_domain(1) = bound_info%domain_x
        boundary%num_domain(2) = bound_info%domain_y
        boundary%num_domain(3) = bound_info%domain_z
        return

      else
        call error_msg('Setup_Processor_Number> # of process is not '// &
                       'domain_x * domain_y * domain_z ')
      end if

    end if

    bsize_x = boundary%box_size_x
    bsize_y = boundary%box_size_y
    bsize_z = boundary%box_size_z

    if (extend) then
      if (extend1) then
        nx  = int(bsize_x / (pairlistdist+2.6_wp))
        ny  = int(bsize_y / (pairlistdist+2.6_wp))
        nz  = int(bsize_z / (pairlistdist+2.6_wp))
        nx1 = int(bsize_x / ((pairlistdist+2.6_wp)/2.0_wp))
        ny1 = int(bsize_y / ((pairlistdist+2.6_wp)/2.0_wp))
        nz1 = int(bsize_z / ((pairlistdist+2.6_wp)/2.0_wp))
      else
        nx  = int(bsize_x / (pairlistdist+2.0_wp))
        ny  = int(bsize_y / (pairlistdist+2.0_wp))
        nz  = int(bsize_z / (pairlistdist+2.0_wp))
        nx1 = int(bsize_x / ((pairlistdist+2.0_wp)/2.0_wp))
        ny1 = int(bsize_y / ((pairlistdist+2.0_wp)/2.0_wp))
        nz1 = int(bsize_z / ((pairlistdist+2.0_wp)/2.0_wp))
      end if
    else
      if (extend1) then
        nx  = int(bsize_x / (pairlistdist+0.6_wp))
        ny  = int(bsize_y / (pairlistdist+0.6_wp))
        nz  = int(bsize_z / (pairlistdist+0.6_wp))
        nx1 = int(bsize_x / ((pairlistdist+0.6_wp)/2.0_wp))
        ny1 = int(bsize_y / ((pairlistdist+0.6_wp)/2.0_wp))
        nz1 = int(bsize_z / ((pairlistdist+0.6_wp)/2.0_wp))
      else
        nx  = int(bsize_x / pairlistdist)
        ny  = int(bsize_y / pairlistdist)
        nz  = int(bsize_z / pairlistdist)
        nx1 = int(2.0_wp * bsize_x / pairlistdist)
        ny1 = int(2.0_wp * bsize_y / pairlistdist)
        nz1 = int(2.0_wp * bsize_z / pairlistdist)
      end if
    end if

    if (nx*ny*nz >= nproc_city) then

      itype = 0
      if (mod(nproc_city,8) == 0) then
        do k = 2, nz
          do j = 2, ny
            do i = 2, nx
              if (i*j*k == nproc_city) then
                itype = itype + 1
                nc(1,itype) = i
                nc(2,itype) = j
                nc(3,itype) = k
                cell_size(1,itype) = nx/i * i
                cell_size(2,itype) = ny/j * j
                cell_size(3,itype) = nz/k * k
              end if
            end do
          end do
        end do
      else if (nproc_city == 1) then
        itype = itype + 1
        nc(1,itype) = 1
        nc(2,itype) = 1
        nc(3,itype) = 1
        cell_size(1,itype) = nx
        cell_size(2,itype) = ny
        cell_size(3,itype) = nz
      else if (nproc_city == 2) then
        itype = itype + 1
        nc(1,itype) = 2
        nc(2,itype) = 1
        nc(3,itype) = 1
        cell_size(1,itype) = nx/2 * 2
        cell_size(2,itype) = ny
        cell_size(3,itype) = nz
      else if (nproc_city == 4) then
        itype = itype + 1
        nc(1,itype) = 2
        nc(2,itype) = 2
        nc(3,itype) = 1
        cell_size(1,itype) = nx/2 * 2
        cell_size(2,itype) = ny/2 * 2
        cell_size(3,itype) = nz
      else
        defined_proc = .false.
        do k = 2, nz
          do j = 2, ny
            do i = 2, nx
              if (i*j*k == nproc_city) then
                itype = itype + 1
                nc(1,itype) = i
                nc(2,itype) = j
                nc(3,itype) = k
                cell_size(1,itype) = nx/i * i
                cell_size(2,itype) = ny/j * j
                cell_size(3,itype) = nz/k * k
              end if
            end do
          end do
        end do
        if (.not. defined_proc) &
          call error_msg('Setup_Processor_Number> MPI Process number can not'//&
                'be defined, please set them manualy')
      end if
    else
      call error_msg('Setup_Processor_Number> Cannot define domains '//         &
                     'and cells. '//                                          &
                     'Smaller MPI processors, or shorter pairlistdist, or '//  &
                     'larger boxsize should be used (see "Chapter: Trouble shooting" in the user manual).')
    end if

    if (itype == 0) then
      call error_msg('Setup_Processor_Number> Cannot define domains '//         &
                     'and cells. '//                                           &
                     'Smaller or adjusted MPI processors, or shorter'//        &
                     ' pairlistdist, or larger boxsize should be used (see "Chapter: Trouble shooting" in the user manual).')
    end if

    k = 0
    maxsize(0) = 100000000000.0_wp
    do i = 1, itype
      size_x = bsize_x/cell_size(1,i)
      size_y = bsize_y/cell_size(2,i)
      size_z = bsize_z/cell_size(3,i)
      maxsize(i) = size_x*size_y*size_z
      if (maxsize(i) < maxsize(k)) &
        k = i
    end do

    boundary%num_domain(1) = nc(1,k)
    boundary%num_domain(2) = nc(2,k)
    boundary%num_domain(3) = nc(3,k)

    return

  end subroutine setup_processor_number

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_boundary_cell
  !> @brief        setup boundary cell information
  !! @authors      NT
  !! @param[in]    table        : flag for use table or not
  !! @param[in]    pairlistdist : pair-list distance
  !! @pmara[in]    water_model  : water model
  !! @param[in]    ensemble     : type of ensemble 
  !! @param[in]    rigid_bond   : flag for rigid-bond
  !! @param[inout] boundary     : boundary information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_boundary_cell(table, pairlistdist, &
                                 water_model, ensemble, rigid_bond,  &
                                 dsize_cg, dmin_size_cg,  boundary)

    ! formal arguments
    logical,                 intent(in)    :: table
    real(wp),                intent(in)    :: pairlistdist
    character(*),            intent(in)    :: water_model
    integer,                 intent(in)    :: ensemble
    logical,                 intent(in)    :: rigid_bond
    logical,                 intent(in)    :: dsize_cg
    real(wp),                intent(in)    :: dmin_size_cg
    type(s_boundary),        intent(inout) :: boundary

    ! local variables
    real(wp)                 :: csize_x, csize_y, csize_z, cutoff
    real(wp)                 :: bsize_x, bsize_y, bsize_z
    integer                  :: i, j, k, ncell
    integer                  :: inb,jnb,knb,inbs,jnbs,knbs,lc,lcnb,ln
    integer                  :: ncell_x, ncell_y, ncell_z
    logical                  :: extend, extend1


    extend  = .false.
    if ( table  .and. &
        (water_model(1:4) == 'TIP3' .or. &
         water_model(1:3) == 'WAT'  .or. &
         water_model(1:3) == 'SOL') )    &
      extend = .true.

    if (rigid_bond) &
      extend = .true.

    extend1 = .false.
    if (ensemble == EnsembleNPT  .or. &
        ensemble == EnsembleNPAT .or. &
        ensemble == EnsembleNPgT) &
      extend1 = .true.

    if (extend) then
      if (extend1) then
        cutoff = pairlistdist + 2.6_wp
      else
        cutoff = pairlistdist + 2.0_wp
      end if
    else
      if (extend1) then
        cutoff = pairlistdist + 0.6_wp
      else
        cutoff = pairlistdist
      end if
    end if

    if (dsize_cg) then
      cutoff=max(dmin_size_cg,pairlistdist)
    endif

    bsize_x = boundary%box_size_x
    bsize_y = boundary%box_size_y
    bsize_z = boundary%box_size_z

    if ( (boundary%num_cells_x /= 0) .and. &
         (boundary%num_cells_y /= 0) .and. &
         (boundary%num_cells_z /= 0)) then
      ncell_x = boundary%num_cells_x
      ncell_y = boundary%num_cells_y
      ncell_z = boundary%num_cells_z
    else
      ncell_x = int(bsize_x/(cutoff/2.0_wp))
      ncell_y = int(bsize_y/(cutoff/2.0_wp))
      ncell_z = int(bsize_z/(cutoff/2.0_wp))
    end if

#ifdef DEBUG
    if (main_rank) then
      write(MsgOut,'(a,f15.8)')  'Debuging > cutoff', cutoff
      write(MsgOut,'(a,3f15.8)') 'Debuging > bsize_[x,y,z]', &
                                 bsize_x, bsize_y, bsize_z
      write(MsgOut,'(a,3i8)')    'Debuging > ncell_[x,y,z]', &
                                 ncell_x, ncell_y, ncell_z
    end if
#endif

    k = mod(ncell_x, boundary%num_domain(1))
    if (k /= 0) &
      ncell_x = ncell_x - k

    k = mod(ncell_y, boundary%num_domain(2))
    if (k /= 0) &
      ncell_y = ncell_y - k

    k = mod(ncell_z, boundary%num_domain(3))
    if (k /= 0) &
      ncell_z = ncell_z - k

    if (ncell_x < 5 .or. ncell_y < 5 .or. ncell_z < 5) &
      call error_msg('Setup_Boundary_Cell> too small boxsize/pairlistdist. '//&
                   'shorter pairlistdist or larger boxsize or less MPI processors'//&
                   ' should be used.')

    csize_x = bsize_x/real(ncell_x, wp)
    csize_y = bsize_y/real(ncell_y, wp)
    csize_z = bsize_z/real(ncell_z, wp)
    ncell   = ncell_x*ncell_y*ncell_z

    boundary%num_cells_x = ncell_x
    boundary%num_cells_y = ncell_y
    boundary%num_cells_z = ncell_z
    boundary%cell_size_x = csize_x
    boundary%cell_size_y = csize_y
    boundary%cell_size_z = csize_z

    ! prepare cell neighbor list
    !
    call alloc_boundary(boundary, BoundaryCells, ncell)

    do k = 0, ncell_z-1
    do j = 0, ncell_y-1
    do i = 0, ncell_x-1

      lc = 1 + i + j*ncell_x + k*ncell_x*ncell_y
      ln = 0

      do knb = k-1, k+1

        if (knb == -1) then
          knbs = ncell_z - 1
        else if (knb == ncell_z) then
          knbs = 0
        else
          knbs = knb
        end if

        do jnb = j-1, j+1

          if (jnb == -1) then
            jnbs = ncell_y - 1
          else if (jnb == ncell_y) then
           jnbs = 0
          else
            jnbs = jnb
          end if

          do inb = i-1, i+1
           if (inb == -1) then
              inbs = ncell_x - 1
            else if (inb == ncell_x) then
              inbs = 0
            else
              inbs = inb
            end if

            lcnb = 1 + inbs + jnbs*ncell_x + knbs*ncell_x*ncell_y
            ln   = ln + 1

            boundary%neighbor_cells(ln,lc) = lcnb
          end do

        end do
      end do

      do knb = k-2, k+2

        if (knb == -2) then
          knbs = ncell_z - 2
        else if (knb == -1) then
          knbs = ncell_z - 1
        else if (knb == ncell_z) then
          knbs = 0
        else if (knb == (ncell_z+1)) then
          knbs = 1
        else
          knbs = knb
        end if

        do jnb = j-2, j+2

          if (jnb == -2) then
            jnbs = ncell_y - 2
          else if (jnb == -1) then
            jnbs = ncell_y - 1
          else if (jnb == ncell_y) then
            jnbs = 0
          else if (jnb == (ncell_y+1)) then
            jnbs = 1
          else
            jnbs = jnb
          end if

          do inb = i-2, i+2
            if (inb == -2) then
              inbs = ncell_x - 2
            else if (inb == -1) then
              inbs = ncell_x - 1
            else if (inb == ncell_x) then
              inbs = 0
            else if (inb == (ncell_x+1)) then
              inbs = 1
            else
              inbs = inb
            end if

            if (abs(inb-i) >= 2 .or. abs(jnb-j) >= 2 .or. abs(knb-k) >= 2) then
              lcnb = 1 + inbs + jnbs*ncell_x + knbs*ncell_x*ncell_y
              ln   = ln + 1
              boundary%neighbor_cells(ln,lc) = lcnb
            end if

            boundary%neighbor_cells(ln,lc) = lcnb
          end do

        end do
      end do
    end do
    end do
    end do

    if (main_rank) then
      write(MsgOut,'(A)') &
           'Setup_Boundary_Cell> Set Variables for Boundary Condition'

      write(MsgOut,'(A20,3I10)')                 &
           '  domains (x,y,z) = ', boundary%num_domain(1), &
                                   boundary%num_domain(2), &
                                   boundary%num_domain(3)

      write(MsgOut,'(A20,3I10)')            &
           '  ncells (x,y,z)  = ', ncell_x, &
                                   ncell_y, &
                                   ncell_z
      write(MsgOut,'(A)') ' '
    end if

    if (ncell_x <= boundary%num_domain(1) .or. &
        ncell_y <= boundary%num_domain(2) .or. &
        ncell_z <= boundary%num_domain(3))     &
      call error_msg( &
          'Setup_Boundary_Cell> ncell_[x,y,z] should be greater than or equal to '//&
          '2*domain_[x,y,z]. Please reduce MPI and increase OpenMP to use the '//   &
          'same number of processors.')
    return

  end subroutine setup_boundary_cell

end module sp_boundary_mod
