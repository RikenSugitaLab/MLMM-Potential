! The main program which runs our driver test case potentials
!
! Copyright (C) 2019, Bernat Font Garcia

! Permission is hereby granted, free of charge, to any person obtaining
! a copy of this software and associated documentation files (the
! "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish,
! distribute, sublicense, and/or sell copies of the Software, and to
! permit persons to whom the Software is furnished to do so, subject to
! the following conditions:
!
! The above copyright notice and this permission notice shall be included
! in all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
! IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
! CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
! TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
! SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!
! ./driver.x -h localhost -p 31415

module at_driver_mod
    use at_enefunc_str_mod
    use molecules_str_mod
    use constants_mod
    use at_energy_str_mod
    use messages_mod
    use at_fsockets_mod, only : open_socket, writebuffer, readbuffer
    use, intrinsic :: iso_c_binding

    implicit none

    ! socket variables
    integer, parameter :: msglen=12   ! length of the headers of the driver/wrapper communication protocol

    contains

    subroutine ml_model( molecule, coord, qmmm, energy)

      type(s_molecule), target, INTENT(IN)    :: molecule
      real(wp),                 INTENT(IN)    :: coord(:,:)
      type(s_energy),   target, INTENT(INOUT) :: energy
      type(s_qmmm),     target, INTENT(INOUT) :: qmmm

      ! socket communication buffers
      character(len=12) :: header  !12
      logical :: isinit=.false.   ! The driver has been initialised by the server
      logical :: hasdata=.true.   ! The driver has finished computing and can send data to server
      real(wp), ALLOCATABLE :: smg(:), rmg(:), coord_qm(:,:)
      real(wp), ALLOCATABLE :: qm_force_(:,:), qm_charge_(:)
      integer,  ALLOCATABLE :: Z(:) 
      integer :: rmg_size, smg_size, dof, i, n_atoms

      integer, pointer      :: socket
      integer, pointer      :: max_order

      !qmmm%qm_debug = .true.
      socket => qmmm%socket
      max_order=> qmmm%max_order

      n_atoms = qmmm%qm_natoms + qmmm%num_qmmmbonds
      dof = n_atoms*3
      smg_size = dof * 2
      smg_size = smg_size + n_atoms 

      if (max_order >= 1) smg_size = smg_size + n_atoms*9
      if (max_order >= 2) smg_size = smg_size + n_atoms*27
      if (max_order >= 3) smg_size = smg_size + n_atoms*81

      rmg_size = n_atoms + 1 + dof
      rmg_size = rmg_size + n_atoms
      rmg_size = rmg_size + dof
      if (max_order > 1) rmg_size = rmg_size + n_atoms*9
      if (max_order > 2) rmg_size = rmg_size + n_atoms*27

      ALLOCATE(smg( smg_size ), rmg( rmg_size ), coord_qm(3, n_atoms), Z(n_atoms))
      ALLOCATE(qm_force_(3,n_atoms), qm_charge_(n_atoms))

      do i = 1, n_atoms
        if (i<=qmmm%qm_natoms) then
          coord_qm(:,i) = coord(:,qmmm%qmatom_id(i))
          Z(i) = qmmm%qm_atomic_no(i)
        else
          coord_qm(:,i) = qmmm%linkatom_coord(:,i-qmmm%qm_natoms) 
          Z(i) = 0
        endif
      end do

      smg(1:dof) = reshape(coord_qm, [dof])
      smg(dof+1:dof+n_atoms) = qmmm%T0
      smg(dof+1+n_atoms:dof+n_atoms+dof) = reshape(qmmm%T1, [dof])

      if (max_order>=1) then
        smg(dof+1+n_atoms+dof:dof+n_atoms+dof+dof*3) = reshape(qmmm%T2, [dof*3])
      end if

      if (max_order>=2) then
        smg(dof+1+n_atoms+dof+dof*3:dof+n_atoms+dof+dof*3+dof*9) = reshape(qmmm%T3, [dof*9])
      end if

      if (max_order>=3) then
        smg(dof+1+n_atoms+dof*13:dof+n_atoms+dof*13+dof*27) = reshape(qmmm%T4, [dof*27])
      end if
      
      hasdata = .true.
      do while (.true.) ! loops forever (or until the wrapper ends!)
          ! reads  qmmm%qm_atomic_no(i)from the socket one message header
          call readbuffer(socket, header, msglen)
          if (qmmm%qm_debug) write(MsgOut,'(''QMMM_debug> message from server:'', a)') trim(header)
          !write(*,*) "has_data", hasdata
          if (qmmm%qm_debug) write(MsgOut,'(''QMMM_debug> has_data:'', l2)') hasdata
          if (trim(header) == "STATUS") then ! the wrapper is inquiring on what we are doing
              if (.not. isinit) then
                  call writebuffer(socket, "NEEDINIT    ", msglen)  ! signals that we need initialization
                  call writebuffer(socket, "DATAREADY   ", msglen)
                  call writebuffer(socket, Z, size(Z))  ! we are idling and eager to compute something
                  call writebuffer(socket, qmmm%qm_classic_charge, size(Z))  ! we are idling and eager to compute something
                 ! write(*,*) "@ message to server: NEEDINIT"
                  if (qmmm%qm_debug) write(MsgOut,'("QMMM_debug> message to server: NEEDINIT")') 
              elseif (hasdata) then
                  call writebuffer(socket, "HAVEDATA    ", msglen)  ! signals that we are done computing and can data
                  if (qmmm%qm_debug) write(MsgOut,'("QMMM_debug> message to server: HAVEDATA")') 
                  !write(*,*) "@ message to server: HAVEDATA"
              else
                  call writebuffer(socket, "READY       ", msglen)  ! we are idling and eager to compute something
                  if (qmmm%qm_debug) write(MsgOut,'("QMMM_debug> message to server: READY")') 
                  !write(*,*) "@ message to server: READY"
              endif

          elseif (trim(header) == "INIT") then     ! the driver is kindly sending a string for initialization
              if (qmmm%qm_debug) write(MsgOut,'("QMMM_debug> Initializing system from server")') 
              !write(*,*) " Initializing system from server"
              isinit=.true. ! we actually do nothing with this string, thanks anyway. could be used to pass some information (e.g. the input parameters, or the index of the replica, from the driver

          elseif (trim(header) == "SENDDATA") then  ! Server wants to send data to the driver
              if (.not. isinit) then
                  if (qmmm%qm_debug) write(MsgOut,'("QMMM_debug> Driver not inilialized")') 
                  !write(*,*) "Driver not iniliasied."
              elseif (hasdata) then
                  if (qmmm%qm_debug) write(MsgOut,'("QMMM_debug> Driver has data to send back to server")') 
                  !write(*,*) "Driver has data to send back to server"
              else ! Driver is ready to receive data
                  if (.not. allocated(rmg)) allocate(rmg(rmg_size))
                  call readbuffer(socket, rmg, size(rmg))

                  energy%qm_ene = rmg(1)
                  qm_force_ = reshape(rmg(2:dof+1), shape(qmmm%T1))
                  qm_charge_ = reshape(rmg(dof+2:dof+1+n_atoms), &
                                            shape(qmmm%T0))

                  qmmm%qm_force = qm_force_(:,1:qmmm%qm_natoms)
                  qmmm%qm_charge = qm_charge_(1:qmmm%qm_natoms)
                  if (qmmm%num_qmmmbonds > 0) then
                    qmmm%linkatom_force = qm_force_(:,qmmm%qm_natoms+1:n_atoms)
                    qmmm%linkatom_charge = qm_charge_(qmmm%qm_natoms+1:n_atoms)
                  endif

                  qmmm%T0_grad = reshape(rmg(dof+2+n_atoms:dof+1+2*n_atoms), &
                                         shape(qmmm%T0_grad))
                  qmmm%T1_grad = reshape(rmg(dof+2+2*n_atoms:dof+1+2*n_atoms+dof), &
                                         shape(qmmm%T1_grad))

                  if (max_order > 1) then
                    qmmm%T2_grad = reshape(rmg(dof+2+2*n_atoms+dof:dof+1+2*n_atoms+dof+3*dof), &
                                           shape(qmmm%T2_grad))
                  end if

                  if (max_order > 2) then
                    qmmm%T3_grad = reshape(rmg(dof+2+2*n_atoms+4*dof:rmg_size), &
                                           shape(qmmm%T3_grad))
                  end if

                  hasdata = .true.
                  exit
              end if

          elseif (trim(header) == "GETDATA") then  ! Server signaling driver to send data
              if (.not. isinit) then
                 ! write(*,*) "Driver not iniliasied."
                  if (qmmm%qm_debug) write(MsgOut,'("QMMM_debug> Driver not inilialized")') 
              elseif (.not. hasdata) then
                  if (qmmm%qm_debug) write(MsgOut,'("QMMM_debug> Driver does not have data to send")') 
                  !write(*,*) "Driver does not have data to send"
              else
                  call writebuffer(socket, "DATAREADY   ", msglen)
                  if (qmmm%qm_debug) write(MsgOut,'("QMMM_debug> message to server: DATAREADY")') 
                  !write(*,*) "@ message to server: DATAREADY"

                  call writebuffer(socket, smg, size(smg)) ! writing data

                  hasdata = .false.
              end if
          else
              write(MsgOut, '('' unexpected header: '', a)') header
              stop "ended"
          endif
      enddo

    end subroutine ml_model


    subroutine helpmessage ! Help banner
        write(*,*) " syntax: driver.x -h hostname -p port "
    end subroutine helpmessage

end module
