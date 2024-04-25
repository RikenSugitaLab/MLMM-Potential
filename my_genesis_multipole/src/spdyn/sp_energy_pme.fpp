!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_energy_pme
!> @brief   Smooth particle mesh ewald method
!! @authors Jaewoon Jung (JJ)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_energy_pme_mod

  use sp_boundary_str_mod
  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use math_libs_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use fft3d_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  real(wp),         save :: el_fact      ! e^2/(4*pi*eps) A*kcal/mol
  real(wp),         save :: alpha        ! alpha
  real(wp),         save :: alpha2m      ! -alpha^2
  real(wp),         save :: alpha2sp     ! 2*alpha/sqrt(pi)
  real(wp),         save :: vol_fact2    ! (2*pi/V)*el_fact
  real(wp),         save :: vol_fact4    ! (4*pi/V)*el_fact
  real(dp), public, save :: u_self       ! Ewald self energy
  real(wp),         save :: box(3)       ! box size
  real(wp),         save :: box_inv(3)   ! Inverse of box size
  real(wp),         save :: bs_fact      ! B-spline factor (1/(n-1)...2)
  real(wp),         save :: bs_fact3     ! bs_fact^3
  real(wp),         save :: bs_fact3d    ! bs_fact^3*(n-1)
  real(wp),         save :: r_scale(3)   ! coordinate-scaling factor (I/L)
  integer,          save :: n_bs         ! Order of B-spline
  integer,          save :: ngrid(3)     ! Number of grid
  integer,          save :: nx           ! process number in x dimension
  integer,          save :: ny           ! process number in y dimension
  integer,          save :: nz           ! process number in z dimension
  integer,          save :: x_start
  integer,          save :: x_end
  integer,          save :: x_local  
  integer,          save :: x_start1
  integer,          save :: x_end1
  integer,          save :: x_local1 
  integer,          save :: y_start
  integer,          save :: y_end
  integer,          save :: y_local  
  integer,          save :: z_start
  integer,          save :: z_end
  integer,          save :: z_local  
  integer,          save :: nlocalx
  integer,          save :: nlocalx1
  integer,          save :: nlocaly
  integer,          save :: nlocalz
  integer,          save :: maxproc
  integer,          save :: ngridmax
  integer,          save :: fft_scheme
  !FEP
  real(dp), public, save :: u_self_preserve ! Ewald self energy
  real(dp), public, save :: u_self_appear   ! Ewald self energy
  real(dp), public, save :: u_self_vanish   ! Ewald self energy

  real(wp),         save, allocatable :: b2(:,:)        ! b^2(hx,hy,hz)
  real(wp),         save, allocatable :: gx(:)
  real(wp),         save, allocatable :: gy(:)
  real(wp),         save, allocatable :: gz(:,:)
  real(wp),         save, allocatable :: vir_fact(:,:,:)! -2*(1+G^2/4a)/G^2
  real(wp),         save, allocatable :: theta(:,:,:,:) ! F^-1[Theta](hz,hy,hx,procz)
  real(wp),         save, allocatable :: f(:,:,:)
  real(wp),         save, allocatable :: v(:,:,:)
  integer,          save, allocatable :: vi(:,:,:)
  real(wp),         save, allocatable :: bsc(:,:,:,:)
  real(wp),         save, allocatable :: bscd(:,:,:,:)
  real(wp),         save, allocatable :: qdf(:,:,:,:)
  real(wp),         save, allocatable :: qdf_real(:)
  real(wp),         save, allocatable :: qdf_work(:,:)
  complex(wp),      save, allocatable :: ftqdf(:)
  complex(wp),      save, allocatable :: ftqdf_work(:,:)

  !$omp threadprivate(f)
  !$omp threadprivate(bsc)
  !$omp threadprivate(bscd)

  complex(wp),      save, allocatable :: ftqdf2(:)
  complex(wp),      save, allocatable :: ftqdf3(:)
  complex(wp),      save, allocatable :: ftqdf4(:)
  complex(wp),      save, allocatable :: ftqdf5(:)
  complex(wp),      save, allocatable :: ftqdf_work2(:,:,:)
  complex(wp),      save, allocatable :: ftqdf_work3(:,:,:)
  complex(wp),      save, allocatable :: c_work(:,:,:), c_work3(:,:,:)
  complex(wp),      save, allocatable :: c_work1(:), c_work2(:)
  real(wp),         save, allocatable :: theta_local(:)
  real(wp),         save, allocatable :: r_work(:), r_work1(:)
  real(wp),         save, allocatable :: r_work2(:), r_work3(:,:,:)

  ! parameters
  integer, parameter      :: NumIndex        = 45
  integer, parameter      :: Index(NumIndex) = (/ &
                              1,   2,   3,   4,   5,   6,   8,   9,  10,  12, &
                             15,  16,  18,  20,  24,  25,  27,  30,  32,  36, &
                             40,  45,  48,  50,  54,  60,  64,  72,  75,  80, &
                             81,  90,  96, 100, 120, 125, 128, 135, 144, 150, &
                            160, 162, 180, 192, 200/)

  ! subroutines
  public  :: setup_pme
  public  :: dealloc_pme
  public  :: pme_pre
  public  :: pme_recip
  !FEP
  public  :: pme_pre_fep
  public  :: pme_recip_fep

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_pme
  !> @brief        Setup for PME with domain decomposition
  !! @authors      JJ
  !! @param[in]    domain   : domain information
  !! @param[in]    boundary : boundary information
  !! @param[inout] enefunc  : energy function
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_pme(domain, boundary, enefunc)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_boundary),        intent(in)    :: boundary
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    real(wp)                 :: bcos, bsin, fact, grid_space
    integer                  :: i, j, k, js
    integer                  :: localmax, ncell
    integer                  :: iproc(3), index_x, index_y, index_z
    integer                  :: remainder, quotient, expo
    integer                  :: alloc_stat

    real(wp),    allocatable :: bs(:)

    if (boundary%type /= BoundaryTypePBC) &
      call error_msg('Setup_PME> Error! BoundaryType is not PBC.')

    ncell    =  domain%num_cell_local + domain%num_cell_boundary

    ! FFT scheme
    !
    fft_scheme = enefunc%fft_scheme

    ! Setting the number of grids (even number)
    !
    ngrid(1) = enefunc%pme_ngrid_x
    ngrid(2) = enefunc%pme_ngrid_y
    ngrid(3) = enefunc%pme_ngrid_z
    ngridmax = max(ngrid(1),ngrid(2),ngrid(3))

    if (ngrid(1) == 0 .and. ngrid(2) == 0 .and. ngrid(3) ==0) then

      grid_space = enefunc%pme_max_spacing
      ngrid(1) = int(boundary%box_size_x/grid_space)
      ngrid(2) = int(boundary%box_size_y/grid_space)
      ngrid(3) = int(boundary%box_size_z/grid_space)
      ngridmax = max(ngrid(1),ngrid(2),ngrid(3))

    end if

    ! grid partition
    !
    if (reciprocal_calc) then

      nx = boundary%num_domain(1)
      ny = boundary%num_domain(2)
      nz = boundary%num_domain(3)

      ! Check grid point in x direction
      !
      if ((fft_scheme == FFT_1dallgather) .or. &
          (fft_scheme == FFT_1dalltoall)) then

        remainder = mod(ngrid(1),2*nx)
        if (remainder /= 0) ngrid(1) = ngrid(1) + 2*nx - remainder

        quotient = ngrid(1)/(2*nx)
        if (quotient <= Index(NumIndex)) then
          do i = 1, NumIndex
            if (quotient <= Index(i)) exit
          end do
          quotient = Index(i)
          ngrid(1) = (2*nx) * quotient
        else
          expo = int(log(real(quotient,wp))/log(real(2,wp)))
          if (2**expo >= quotient) then
            quotient = 2**expo
          else
            quotient = 2**(expo+1)
          end if
          ngrid(1) = (2*nx) * quotient
        end if

      else if (fft_scheme == FFT_2dalltoall) then

        remainder = mod(ngrid(1),2*nx*ny)
        if (remainder /= 0) ngrid(1) = ngrid(1) + 2*nx*ny - remainder

        quotient = ngrid(1)/(2*nx*ny)
        if (quotient <= Index(NumIndex)) then
          do i = 1, NumIndex
            if (quotient <= Index(i)) exit
          end do
          quotient = Index(i)
          ngrid(1) = (2*nx*ny) * quotient
        else
          expo = int(log(real(quotient,wp))/log(real(2,wp)))
          if (2**expo >= quotient) then
            quotient = 2**expo
          else
            quotient = 2**(expo+1)
          end if
          ngrid(1) = (2*nx*ny) * quotient
        end if

      end if

      ! check grid points in y direction
      !
      if ((fft_scheme == FFT_1dalltoall) .or. &
          (fft_scheme == FFT_2dalltoall)) then

        if (mod(nz,2) == 0) then

          remainder = mod(ngrid(2),ny*nz)
          if (remainder /= 0) ngrid(2) = ngrid(2) + nz*ny - remainder

          quotient = ngrid(2)/(nz*ny)
          if (quotient <= Index(NumIndex)) then
            do i = 1, NumIndex
              if (quotient <= Index(i)) exit
            end do
            quotient = Index(i)
            ngrid(2) = (nz*ny) * quotient
          else
            expo = int(log(real(quotient,wp))/log(real(2,wp)))
            if (2**expo >= quotient) then
              quotient = 2**expo
            else
              quotient = 2**(expo+1)
            end if
            ngrid(2) = (nz*ny) * quotient
          end if

        else

          remainder = mod(ngrid(2),ny*nz*2)
          if (remainder /= 0) ngrid(2) = ngrid(2) + 2*nz*ny - remainder

          quotient = ngrid(2)/(2*nz*ny)
          if (quotient <= Index(NumIndex)) then
            do i = 1, NumIndex
              if (quotient <= Index(i)) exit
            end do
            quotient = Index(i)
            ngrid(2) = (2*nz*ny) * quotient
          else
            expo = int(log(real(quotient,wp))/log(real(2,wp)))
            if (2**expo >= quotient) then
              quotient = 2**expo
            else
              quotient = 2**(expo+1)
            end if
            ngrid(2) = (2*nz*ny) * quotient
          end if

        end if

      else if (fft_scheme == FFT_1dallgather) then

        remainder = mod(ngrid(2),2*ny)
        if (remainder /= 0) ngrid(2) = ngrid(2) + 2*ny - remainder

        quotient = ngrid(2)/(2*ny)
        if (quotient <= Index(NumIndex)) then
          do i = 1, NumIndex
            if (quotient <= Index(i)) exit
          end do
          quotient = Index(i)
          ngrid(2) = (2*ny) * quotient

        else
          expo = int(log(real(quotient,wp))/log(real(2,wp)))
          if (2**expo >= quotient) then
            quotient = 2**expo
          else
            quotient = 2**(expo+1)
          end if
          ngrid(2) = (2*ny) * quotient
        end if

      end if

      ! Check grid point in z direction
      !
      if ((fft_scheme == FFT_1dalltoall) .or. &
          (fft_scheme == FFT_2dalltoall)) then

        remainder = mod(ngrid(3),nz*nx)
        if (remainder /= 0) ngrid(3) = ngrid(3) + nz*nx - remainder

        quotient = ngrid(3)/(nz*nx)
        if (quotient <= Index(NumIndex)) then
          do i = 1, NumIndex
            if (quotient <= Index(i)) exit
          end do
          quotient = Index(i)
          ngrid(3) = (nz*nx) * quotient
        else
          expo = int(log(real(quotient,wp))/log(real(2,wp)))
          if (2**expo >= quotient) then
            quotient = 2**expo
          else
            quotient = 2**(expo+1)
          end if
          ngrid(3) = (nz*nx) * quotient
        end if

        if (fft_scheme == FFT_1dalltoall) then

          remainder = mod(ngrid(3),nz*ny)
          if (remainder /= 0) ngrid(3) = ngrid(3) + nz*ny - remainder

          quotient = ngrid(3)/(nz*ny)
          if (quotient <= Index(NumIndex)) then
            do i = 1, NumIndex
              if (quotient <= Index(i)) exit
            end do
            quotient = Index(i)
            ngrid(3) = (nz*ny) * quotient
          else
            expo = int(log(real(quotient,wp))/log(real(2,wp)))
            if (2**expo >= quotient) then
              quotient = 2**expo
            else
              quotient = 2**(expo+1)
            end if
            ngrid(3) = (nz*ny) * quotient
          end if
        end if

      else if (fft_scheme == FFT_1dallgather) then

        remainder = mod(ngrid(3),nz)
        if (remainder /= 0) ngrid(3) = ngrid(3) + nz - remainder

        quotient = ngrid(3)/nz
        if (quotient <= Index(NumIndex)) then
          do i = 1, NumIndex
            if (quotient <= Index(i)) exit
          end do
          quotient = Index(i)
          ngrid(3) = nz * quotient
        else
          expo = int(log(real(quotient,wp))/log(real(2,wp)))
          if (2**expo >= quotient) then
            quotient = 2**expo
          else
            quotient = 2**(expo+1)
          end if
          ngrid(3) = nz * quotient
        end if

      end if

      if ((enefunc%pme_ngrid_x /= ngrid(1)) .or. &
          (enefunc%pme_ngrid_y /= ngrid(2)) .or. &
          (enefunc%pme_ngrid_z /= ngrid(3))) then

        if (main_rank) then
          write(MsgOut,'(A)') &
            ''
          if (enefunc%pme_ngrid_x == 0 .and. enefunc%pme_ngrid_y == 0 .and. &
              enefunc%pme_ngrid_z == 0) then
            write(MsgOut,'(A)') &
              'Setup_PME> Proper PME grid number was generated automatically'
          else
            write(MsgOut,'(A)') &
              '  WARNING: PME grid number is different from the input'
          end if
          write(MsgOut,'(A20,3I10)') &
            '  pme_ngrid(x,y,z)= ', ngrid(1), ngrid(2), ngrid(3)
          write(MsgOut,'(A)') &
            ''
        end if
      end if

      ngridmax = max(ngrid(1),ngrid(2),ngrid(3))

      ! index for each processor
      !

      iproc(1) = mod(my_city_rank, boundary%num_domain(1))
      iproc(2) = mod(my_city_rank/boundary%num_domain(1),boundary%num_domain(2))
      iproc(3) = my_city_rank/(boundary%num_domain(1)*boundary%num_domain(2))

      nlocalx  = int((ngrid(1)+nx-1)/nx)
      nlocaly  = int((ngrid(2)+ny-1)/ny)
      nlocalz  = int((ngrid(3)+nz-1)/nz)

      x_start  = (iproc(1))*nlocalx + 1
      y_start  = (iproc(2))*nlocaly + 1
      z_start  = (iproc(3))*nlocalz + 1

      x_end    = (iproc(1)+1)*nlocalx 
      y_end    = (iproc(2)+1)*nlocaly 
      z_end    = (iproc(3)+1)*nlocalz 

      x_local  = nlocalx
      y_local  = nlocaly
      z_local  = nlocalz

      localmax = max(nlocalx*nlocaly,nlocalx*nlocalz,nlocaly*nlocalz)
      index_x  = iproc(2)*nlocaly*nlocalz + iproc(3)
      index_y  = iproc(3)*nlocalx*nlocalz + iproc(1)
      index_z  = iproc(1)*nlocalx*nlocaly + iproc(2)

      ! new communicator according to the grid index
      !
#ifdef HAVE_MPI_GENESIS
      call mpi_comm_split(mpi_comm_city, index_x,my_city_rank,grid_commx,ierror)
      call mpi_comm_size (grid_commx, nprocx, ierror)
      call mpi_comm_rank (grid_commx, my_x_rank, ierror)
      call mpi_comm_split(mpi_comm_city, index_y,my_city_rank,grid_commy,ierror)
      call mpi_comm_size (grid_commy, nprocy, ierror)
      call mpi_comm_rank (grid_commy, my_y_rank, ierror)
      call mpi_comm_split(mpi_comm_city, index_z,my_city_rank,grid_commz,ierror)
      call mpi_comm_size (grid_commz, nprocz, ierror)
      call mpi_comm_rank (grid_commz, my_z_rank, ierror)
      call mpi_comm_split(mpi_comm_city, my_city_rank/(nx*ny), my_city_rank, &
                          grid_commxy, ierror)
      call mpi_comm_size(grid_commxy, nprocxy, ierror)
      call mpi_comm_rank(grid_commxy, my_xy_rank, ierror)
#endif
      maxproc = max(nprocx,nprocy,nprocz,1)

      nlocalx1 = nlocalx/2 + 1
      if (iproc(1) == 0) then
        x_local1 = nlocalx1
        x_start1 = 1
        x_end1   = x_local1
      else
        x_local1 = nlocalx1 - 1
        x_start1 = iproc(1)*nlocalx/2 + 2
        x_end1   = x_start1 + x_local1 - 1
      end if
    end if

    ! Setting common parameters
    !
    el_fact  = ELECOEF/enefunc%dielec_const
    n_bs     = enefunc%pme_nspline
    alpha    = enefunc%pme_alpha
    alpha2m  = -alpha**2
    alpha2sp = 2.0_wp*alpha/sqrt(PI)

    if (reciprocal_calc) then

      j = 1
      do i = 1, n_bs - 2
        j = j*(n_bs - i)
      end do

      bs_fact   = 1.0_wp/real(j,wp)
      bs_fact3  = bs_fact**3
      bs_fact3d = bs_fact3 * real(n_bs - 1,wp)

      ! Preparing b2=b(h)^2, h shifted
      !
      allocate(bs(n_bs), b2(ngridmax,3))

      call b_spline_coef(n_bs, 1.0_wp, bs)

      do i = 1, n_bs - 1
        bs(i) = bs(i)*bs_fact
      end do

      do k = 1, 3

        fact = 2.0_wp * PI/real(ngrid(k),wp)

        do j = 1, ngrid(k)/2 + 1

          js = j - 1
          bcos = 0.0_wp
          bsin = 0.0_wp

          do i = 0, n_bs - 2
            bcos = bcos + bs(i+1) * cos(real(js*i,wp)*fact)
            bsin = bsin + bs(i+1) * sin(real(js*i,wp)*fact)
          end do

          b2(j, k) = 1.0_wp/(bcos**2 + bsin**2)

        end do

        do j = ngrid(k)/2 + 2, ngrid(k)

          js = j - 1 - ngrid(k)
          bcos = 0.0_wp
          bsin = 0.0_wp

          do i = 0, n_bs - 2
            bcos = bcos + bs(i+1) * cos(real(js*i,wp)*fact)
            bsin = bsin + bs(i+1) * sin(real(js*i,wp)*fact)
          end do

          b2(j, k) = 1.0_wp/(bcos**2 + bsin**2)

        end do

      end do

      if (mod(n_bs,2) == 1) then
        do k = 1, 3
          b2(ngrid(k)/2 + 1,k) = 0.0_wp
        end do
      end if

      ! Prepareing theta = F^-1[theta](h), h shifted
      !
      if (fft_scheme == FFT_1dallgather) then
        allocate(gx(x_local1), gy(nlocaly), gz(nlocalz,nprocz), &
                 vir_fact(nlocalz, nlocaly, x_local1),          &
                 v(3,MaxAtom,ncell),                            &
                 vi(3,MaxAtom,ncell),                           &
                 qdf(nlocalx,nlocaly,nlocalz,nthread),          &
                 qdf_real(nlocalx*nlocaly*nlocalz),             &
                 theta(nlocalz,nlocaly,x_local1,nprocz),        &
                 qdf_work(nlocalx*nlocaly*nlocalz,maxproc),     &
                 ftqdf(nlocalx1*nlocaly*nlocalz),               &
                 ftqdf_work(nlocalx1*nlocaly*nlocalz,maxproc),  &
                 stat = alloc_stat)
      else if (fft_scheme == FFT_1dalltoall) then
        allocate(gx(x_local1), gy(nlocaly), gz(ngrid(3),1),     &
                 vir_fact(ngrid(3), nlocaly, x_local1),         &
                 vi(3,MaxAtom,ncell),                           &
                 qdf(nlocalx,nlocaly,nlocalz,nthread),          &
                 qdf_real(nlocalx*nlocaly*nlocalz),             &
                 theta(ngrid(3),nlocaly,x_local1,1),            &
                 qdf_work(nlocalx*nlocaly*nlocalz,maxproc),     &
                 ftqdf(nlocalx1*nlocaly*nlocalz),               &
                 ftqdf_work(nlocalx1*nlocaly*nlocalz*maxproc,1),&
                 stat = alloc_stat)
      else if (fft_scheme == FFT_2dalltoall) then
        allocate(gx(x_local1), gy(nlocaly), gz(nlocalz,nprocz), &
                 vir_fact(nlocalz, nlocaly, x_local1),          &
                 vi(3,MaxAtom,ncell),                           &
                 theta_local(nlocalx*nlocaly*nlocalz),          &
                 c_work(nlocalz,nlocalx,nlocaly),               &
                 c_work1(nlocalx*nlocaly*nlocalz),              &
                 c_work2(nlocalx*nlocaly*nlocalz),              &
                 c_work3(nlocaly,nlocalx,nlocalz),              &
                 r_work (nlocalx*nlocaly*nlocalz),              &
                 r_work1(nlocalx*nlocaly*nlocalz),              &
                 r_work2(nlocalx*nlocaly*nlocalz),              &
                 r_work3(nlocaly,nlocalx,nlocalz),              &
                 theta(1,1,1,1),                                &
                 qdf(nlocalx,nlocaly,nlocalz,nthread),          &
                 qdf_real(nlocalx*nlocaly*nlocalz),             &
                 qdf_work(nlocalx*nlocaly*nlocalz,1),           &
                 ftqdf_work(nlocalx1*nlocaly*nlocalz,1),        &
                 stat = alloc_stat)
      end if
      if (alloc_stat /= 0) call error_msg_alloc

      if (fft_scheme == FFT_1dalltoall) then
        allocate(ftqdf2(nlocalx1*nlocaly*nlocalz),              &
                 ftqdf3(nlocalx1*nlocaly*nlocalz),              &
                 ftqdf4(nlocalx1*nlocaly*nlocalz),              &
                 ftqdf5(nlocalx1*nlocaly*nlocalz),              &
                 ftqdf_work2(ngrid(3),nlocalx1,nlocaly),        &
                 ftqdf_work3(ngrid(3),nlocalx1,nlocaly),        &
                 stat = alloc_stat)
        if (alloc_stat /= 0) call error_msg_alloc
      end if

      !$omp parallel private(alloc_stat)

      allocate(f(3,MaxAtom,ncell),                                   &
               bsc(n_bs,3,MaxAtom,ncell),bscd(n_bs,3,MaxAtom,ncell), &
               stat = alloc_stat)
      if (alloc_stat /= 0)   call error_msg_alloc

      !$omp end parallel

      deallocate(bs)

    end if

    return

  end subroutine setup_pme

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_pme
  !> @brief        deallocate all arrays used in PME
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_pme

    ! local variables
    integer :: dealloc_stat


    if (.not. allocated(f)) &
      return 

    !$omp parallel private(dealloc_stat)

    deallocate(f, bsc, bscd, stat = dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc

    !$omp end parallel

    deallocate(gx, gy, gz, vir_fact, vi, qdf, qdf_real, &
               theta, qdf_work, ftqdf_work,      &
               stat = dealloc_stat)

    if (dealloc_stat /= 0) call error_msg_dealloc
    if (fft_scheme == FFT_2dalltoall) then
      deallocate(theta_local, c_work, c_work1, c_work2,         &
                 c_work3, r_work, r_work1, r_work2, r_work3,    &
                 stat = dealloc_stat)
      if (dealloc_stat /= 0) call error_msg_dealloc
    else
      deallocate(ftqdf, stat = dealloc_stat)
      if (dealloc_stat /= 0) call error_msg_dealloc
    end if

    if (fft_scheme == FFT_1dalltoall) then
      deallocate(ftqdf2, ftqdf3, ftqdf4, ftqdf5, &
                 ftqdf_work2, ftqdf_work3,       &
                 stat = dealloc_stat)
      if (dealloc_stat /= 0) call error_msg_dealloc
    end if

    deallocate(b2, stat = dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_pme

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_pre
  !> @brief        Prepare functions for PME calculation
  !! @authors      JJ
  !! @param[in]    domain   : domain information
  !! @param[in]    boundary : boundary information
  !! @note         Extracted from setup_pme for NPT calculation
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_pre(domain, boundary)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_boundary),        intent(in)    :: boundary

    ! local variables
    integer                  :: i, j, k, is, js, ks, ix, iy, iz, iproc
    integer                  :: ic, kk, iprocx, iprocy, ixs, izs, k_f, k_t
    integer                  :: nix, nix1, niy, niz, iorg
    integer                  :: ix_start, ix_end, iz_start
    integer                  :: iz_end, iy_start, iy_end
    real(wp)                 :: fact, gfact(3), g2


    box(1) = boundary%box_size_x
    box(2) = boundary%box_size_y
    box(3) = boundary%box_size_z

    do k = 1, 3
      box_inv(k) = 1.0_wp / box(k)
      r_scale(k) = real(ngrid(k),wp) * box_inv(k)
    end do

    vol_fact2 = 2.0_wp * PI * el_fact * box_inv(1)*box_inv(2)*box_inv(3)
    vol_fact4 = 2.0_wp * vol_fact2

    if (fft_scheme == FFT_2dalltoall) then
      theta_local(1:nlocalx*nlocaly*nlocalz) = 0.0_wp
    end if

    ! Prepareing theta = F^-1[theta](h), h shifted
    ! Gx, Gy, Gz, and vir_fact=-2(1+G^2/4a)/G^2 are also prepared 
    ! for virial calculation
    !
    do k = 1, 3
      gfact(k) = 2.0_wp * PI * box_inv(k)
    end do

    fact = 0.25_wp / alpha2m

    if (fft_scheme == FFT_2dalltoall) then

      do i = x_start1, x_end1

        ix = i - x_start1 + 1
        is = i - 1
        gx(ix) = gfact(1) * real(is,wp)

        do j = y_start, y_end

          iy = j - y_start + 1
          if (j <= ngrid(2)/2+1) then
            js = j - 1
          else
            js = j - 1 - ngrid(2)
          end if

          gy(iy) = gfact(2) * real(js,wp)

          do iproc = 1, nprocz
            do iz = 1, nlocalz
              k = (iproc-1)*nlocalz + iz
              if (k <= ngrid(3)/2+1) then
                ks = k - 1
              else
                ks = k - 1 - ngrid(3)
              end if
              gz(iz,iproc) = gfact(3) * real(ks,wp)
              g2 = gx(ix)**2 + gy(iy)**2 + gz(iz,iproc)**2
              if (g2 > EPS) then
                if (iproc == (my_z_rank+1)) &
                  vir_fact(iz,iy,ix) = -2.0_wp * (1.0_wp - g2 * fact)/g2
                if (iproc == (my_z_rank+1)) then
                  kk = ix + (iy-1)*nlocalx1 + (iz-1)*nlocalx1*nlocaly
                  theta_local(kk) = b2(i,1) * b2(j,2) * b2(k,3)             & 
                                    * exp(g2 * fact)/g2
                endif
              end if

            end do

          end do

        end do

      end do

      if (my_x_rank == 0 .and. my_y_rank == 0 .and. my_z_rank == 0) &
        theta_local(1) = 0.0d0

      nix = (nlocalx1-1)/nprocy
      nix1 = nix+1
      niy = nlocaly/nprocz
      niz = nlocalz/nprocx

      call mpi_alltoall(theta_local, nlocalx1*nlocaly*niz, mpi_wp_real, &
                        r_work1,     nlocalx1*nlocaly*niz, mpi_wp_real, &
                        grid_commx, ierror)

      do iz = 1, nlocalz
        do iy = 1, nlocaly
          do ix = 1, nlocalx1
            k_f = ix + (iy-1)*nlocalx1 + (iz-1)*nlocalx1*nlocaly
            k_t = iy + (ix-1)*nlocaly  + (iz-1)*nlocalx1*nlocaly
            r_work2(k_t) = r_work1(k_f)
          enddo
        enddo
      enddo

      r_work1 = 0.0_wp
      do iprocx = 0, nprocx-1
        iz_start = iprocx*niz+1
        iz_end = iz_start+niz-1
        do iprocy = 0, nprocy-1
          ix_start = iprocy*nix+1
          ix_end   = ix_start+nix-1
          if(iprocy == nprocy-1) ix_end = ix_start+nix1-1
          iorg = nlocaly*nix1*niz*(iprocx+nprocx*iprocy)
          do iz = iz_start, iz_end
            izs = iz - iz_start
            do ix = ix_start, ix_end
              ixs = ix - ix_start
              k = iorg+nlocaly*ixs+nlocaly*nix1*izs
              do iy = 1, nlocaly
                kk = iy + (ix-1)*nlocaly  + (iz-1)*nlocalx1*nlocaly
                r_work1(k+iy) = r_work2(kk)
              enddo
            enddo
          enddo
        enddo
      enddo

      call mpi_alltoall(r_work1, nlocaly*nix1*niz, mpi_wp_real, &
                        r_work2, nlocaly*nix1*niz, mpi_wp_real, &
                        grid_commxy, ierror)

      do iprocx = 0, nprocx-1
        iz_start = iprocx*niz+1
        iz_end = iz_start+niz-1
        do iprocy = 0, nprocy-1
          ix_start = iprocy*nix1+1
          ix_end   = ix_start+nix1-1
          iorg = nlocaly*nix1*niz*(iprocx+nprocx*iprocy)
          do iz = iz_start, iz_end
            izs = iz - iz_start
            do ix = ix_start, ix_end
              ixs = ix - ix_start
              k = iorg+nlocaly*ixs+nlocaly*nix1*izs
              do iy = 1, nlocaly
                r_work3(iy,ix,iz) = r_work2(k+iy)
              enddo
            enddo
          enddo
        enddo
      enddo

      do iz = 1, nlocalz
        do ix = 1, nix1
          do iproc = 1, nprocy
            k = (iproc-1)*nlocaly
            ixs = (iproc-1)*nix1 + ix
            do iy = 1, nlocaly
              kk = iz + nlocalz*(ixs-1) + nlocalz*nix1*nprocy*(iy-1)
              r_work2(kk) = r_work3(iy,ixs,iz)
            end do
          end do
        end do
      end do

      call mpi_alltoall(r_work2, nix1*nprocy*niy*nlocalz, mpi_wp_real, &
                        r_work1, nix1*nprocy*niy*nlocalz, mpi_wp_real, &
                        grid_commz, ierror)

      do iy = 1,nlocaly
        do ix = 1,nix1*nprocy
          do iz = 1,nlocalz
             kk = iz + nlocalz*(ix-1) + nlocalz*nix1*nprocy*(iy-1)
             theta_local(kk) = r_work1(kk)
           enddo
        enddo
      enddo
      if (my_x_rank == 0 .and. my_y_rank == 0) &
        theta(1,1,1,1) = 0.0_wp

    else if (fft_scheme == FFT_1dalltoall) then

      niy = nlocaly / nprocz
      iy_start = y_start + niy*my_z_rank
      iy_end   = iy_start + niy - 1

      do i = x_start1, x_end1

        ix = i - x_start1 + 1
        is = i - 1
        gx(ix) = gfact(1) * real(is,wp)

        do j = iy_start, iy_end

          iy = j - iy_start + 1

          if (j <= ngrid(2)/2+1) then
            js = j - 1
          else
            js = j - 1 - ngrid(2)
          end if
          gy(iy) = gfact(2) * real(js,wp)

          do k = 1, ngrid(3)
            if (k <= ngrid(3)/2+1) then
              ks = k - 1
            else
              ks = k - 1 - ngrid(3)
            end if
            gz(k,1) = gfact(3) * real(ks,wp)

            g2 = gx(ix)**2 + gy(iy)**2 + gz(k,1)**2
            if (g2 > EPS) then
              vir_fact(k,iy,ix) = -2.0_wp * (1.0_wp - g2 * fact) / g2
              theta(k,iy,ix,1) = b2(i,1)*b2(j,2)*b2(k,3)*exp(g2*fact) / g2
            end if

          end do

        end do

      end do
      if (my_x_rank == 0 .and. my_y_rank == 0 .and. my_z_rank == 0) &
        theta(1,1,1,1) = 0.0_wp

    else

      do i = x_start1, x_end1

        ix = i - x_start1 + 1
        is = i - 1
        gx(ix) = gfact(1) * real(is,wp)

        do j = y_start, y_end

          iy = j - y_start + 1

          if (j <= ngrid(2)/2+1) then
            js = j - 1
          else
            js = j - 1 - ngrid(2)
          end if

          gy(iy) = gfact(2) * real(js,wp)

          do iproc = 1, nprocz
            do iz = 1, nlocalz
              k = (iproc-1)*nlocalz + iz
              if (k <= ngrid(3)/2+1) then
                ks = k - 1
              else
                ks = k - 1 - ngrid(3)
              end if
              gz(iz,iproc) = gfact(3) * real(ks,wp)
              g2 = gx(ix)**2 + gy(iy)**2 + gz(iz,iproc)**2
              if (g2 > EPS) then
                if (iproc == (my_z_rank+1)) &
                  vir_fact(iz,iy,ix) = -2.0_wp * (1.0_wp - g2 * fact)/g2
                theta(iz,iy,ix,iproc) = &
                     b2(i,1) * b2(j,2) * b2(k,3) * exp(g2 * fact)/g2
              end if

            end do

          end do

        end do

      end do
      if (my_x_rank == 0 .and. my_y_rank == 0) &
        theta(1,1,1,1) = 0.0_wp

    end if

    if (my_city_rank == 0) &
      vir_fact(1,1,1) = 0.0_wp

    ! Calculating self energy
    !
    u_self = 0.0_dp

    do i = 1, domain%num_cell_local
      do ix = 1, domain%num_atom(i)
        u_self = u_self + domain%charge(ix,i)**2
      end do 
    end do

    u_self = - u_self * el_fact * alpha/sqrt(PI)

    return

  end subroutine pme_pre

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_recip
  !> @brief        Calculate PME reciprocal part with domain decomposition
  !! @authors      JJ
  !! @param[in]    domain : domain information
  !! @param[inout] force  : forces of target systems
  !! @param[inout] virial : virial term of target systems
  !! @param[inout] eelec  : electrostatic energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_recip(domain, force, virial, eelec)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    real(dp),                intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: virial(3,3,nthread)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: elec_temp
    real(wp)                 :: vr(3), dv(3), bsc_tmp(3,n_bs), bscd_tmp(3,n_bs)
    real(wp)                 :: vxx, vyy, vzz, temp
    real(wp)                 :: tqq
    real(wp)                 :: f_1, f_2, f_3
    integer                  :: i, k, k1, icel, iproc, id
    integer                  :: ix, iv, iy, iz, ixs, iys, izs
    integer                  :: ixyz, nix, nix1, niy, niz
    integer                  :: ixx, iyy, izz, ii(3)
    integer                  :: vx_tmp(n_bs), vy_tmp(n_bs), vz_tmp(n_bs)
    integer                  :: ncell, omp_get_thread_num
    integer                  :: kk, iorg, k_f, k_t
    integer                  :: iprocx, iprocy
    integer                  :: ix_start, ix_end, iz_start, iz_end

    complex(wp), allocatable :: work1(:), work2(:)
    real(dp),        pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    integer,         pointer :: natom(:)


    coord  => domain%coord
    charge => domain%charge
    natom  => domain%num_atom
 
    ncell  = domain%num_cell_local + domain%num_cell_boundary

    ! Initializing the energy and force
    !
    qdf_real(1:nlocalx*nlocaly*nlocalz) = 0.0_wp

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, iv, ix, iy, iz, icel, k, ii, dv, izz, izs, iyy, iys,  &
    !$omp         iz_start, iz_end, ix_start, ix_end, iorg, k_f, k_t,          &
    !$omp         ixx, ixs, vxx, vyy, vzz, iproc, tqq, k1, temp, vr, &
    !$omp         work1, work2, elec_temp, f_1, f_2, f_3, ixyz, vx_tmp, vy_tmp,&
    !$omp         vz_tmp, nix, nix1, niy, niz, bsc_tmp, bscd_tmp, kk, iprocx,  &
    !$omp         iprocy)
    !

#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    !$omp barrier

    ! initialization
    !
    allocate(work1(ngridmax), work2(2*ngridmax))

    ! Initialization of Q-fanction
    !
!#ifndef USE_GPU
    do iz = 1, nlocalz
      do iy = 1, nlocaly
        do ix = 1, nlocalx
          qdf(ix,iy,iz,id+1) = 0.0_wp
        end do
      end do
    end do

    ! Calculateing Q-fanction
    !
    do icel = id+1, ncell, nthread
      do i = 1, natom(icel)
        vr(1:3) = coord(1:3,i,icel) * r_scale(1:3)
        vr(1:3) = vr(1:3) + real(ngrid(1:3)/2,wp) &
                - real(ngrid(1:3),wp)*anint(vr(1:3)/real(ngrid(1:3),wp))

        do k = 1, 3
          ii(k) = int(vr(k))
          if (ii(k) >= ngrid(k)) then
            vr(k) = real(ngrid(k),wp) - 0.00001_wp
            ii(k) = ngrid(k)-1
          endif
          vi(k,i,icel) = ii(k)
          dv(k) = vr(k) - real(ii(k),wp)
          call b_spline_dev_coef(n_bs, dv(k), bsc(1,k,i,icel), bscd(1,k,i,icel))
        end do

        nix = 0 
        niy = 0 
        niz = 0 

        do ixyz = 1, n_bs
          ixs = ii(1) - ixyz + 2
          iys = ii(2) - ixyz + 2
          izs = ii(3) - ixyz + 2
          if (ixs <= 0) ixs = ixs + ngrid(1)
          if (iys <= 0) iys = iys + ngrid(2)
          if (izs <= 0) izs = izs + ngrid(3)

          if (ixs >= x_start .and. ixs <= x_end) then
            nix = nix + 1
            vx_tmp(nix) = ixs - x_start + 1
            bsc_tmp(1,nix) = bsc(ixyz,1,i,icel)
          end if

          if (iys >= y_start .and. iys <= y_end) then
            niy = niy + 1
            vy_tmp(niy) = iys - y_start + 1
            bsc_tmp(2,niy) = bsc(ixyz,2,i,icel)
          end if

          if (izs >= z_start .and. izs <= z_end) then
            niz = niz + 1
            vz_tmp(niz) = izs - z_start + 1
            bsc_tmp(3,niz) = bsc(ixyz,3,i,icel)
          end if
        end do

        do iz = 1, niz
          izs = vz_tmp(iz)
          do iy = 1, niy
            iys = vy_tmp(iy)
            do ix = 1, nix
              ixs = vx_tmp(ix)
              qdf(ixs,iys,izs,id+1) = qdf(ixs,iys,izs,id+1)           &
                                    + bsc_tmp(1,ix)*bsc_tmp(2,iy)     &
                                     * bsc_tmp(3,iz)*charge(i,icel)   &
                                     * bs_fact3
            end do
          end do
        end do
      end do
    end do

    !$omp barrier
    do iproc = 1, nthread
      !$omp do
      do iz = 1, nlocalz
        do iy = 1, nlocaly
          k = (iy-1)*nlocalx + (iz-1)*nlocalx*nlocaly
          do iv = 1, nlocalx, 16
            do ix = iv, min0(iv+15,nlocalx)
              qdf_real(k+ix) = qdf_real(k+ix) + qdf(ix,iy,iz,iproc)
            end do
          end do
        end do
      end do
    end do

    !$omp barrier

!#else
!    !$omp master
!    call gpu_pme_recip_build_qdf( qdf_real, coord, charge, natom, &
!         MaxAtom, ncell, nlocalx, nlocaly, nlocalz, ngrid, &
!         x_start, x_end, y_start, y_end, z_start, z_end, &
!         r_scale, &
!         bs_fact3, &
!         vi, bsc, bscd )
!    !$omp end master
!#endif

    ! FFT (forward)
    !
    if ( (fft_scheme == FFT_1dalltoall) .or.   &
         (fft_scheme == FFT_2dalltoall) ) then

      !$omp barrier
      !$omp master
#ifdef HAVE_MPI_GENESIS
      call mpi_alltoall(qdf_real, nlocalx*nlocaly*nlocalz/nprocx, mpi_wp_real, &
                        qdf_work, nlocalx*nlocaly*nlocalz/nprocx, mpi_wp_real, &
                        grid_commx, ierror)
#else
      qdf_work(1:nlocalx*nlocaly*nlocalz,1) = qdf_real(1:nlocalx*nlocaly*nlocalz)
#endif
      !$omp end master
      !$omp barrier

      if (fft_scheme == FFT_1dalltoall) then
        call fft3d_1d_alltoall(qdf_work, ftqdf, ftqdf2, ftqdf3, ftqdf4, ftqdf5,    &
                         ftqdf_work, ftqdf_work, ftqdf_work, ftqdf_work,           &
                         ftqdf_work2, ftqdf_work3, work1, work2, nlocalx, nlocaly, &
                         nlocalz, nlocalx1, x_local1, ngrid(1), ngrid(2), ngrid(3),&
                         nprocx, nprocy, nprocz, id, nthread, my_x_rank, my_y_rank,&
                         my_z_rank, x_start1, x_end1, y_start, y_end, z_start,     &
                         z_end, grid_commx, grid_commy, grid_commz)

      else
        call fft3d_2d_alltoall(qdf_work, ftqdf_work, work1, work2, nlocalx,        &
                         nlocaly, nlocalz, nlocalx1, x_local1, ngrid(1), ngrid(2), &
                         ngrid(3), nprocx, nprocy, nprocz, id, nthread, my_x_rank, &
                         my_y_rank, my_z_rank, x_start1, x_end1, y_start, y_end,   &
                         z_start, z_end, grid_commx, grid_commy, grid_commz,       &
                         grid_commxy, c_work, c_work1, c_work2, c_work3)
      end if

    else if (fft_scheme == FFT_1dallgather) then

      !$omp barrier
      !$omp master
#ifdef HAVE_MPI_GENESIS
      call mpi_allgather(qdf_real, nlocalx*nlocaly*nlocalz, mpi_wp_real,         &
                         qdf_work, nlocalx*nlocaly*nlocalz, mpi_wp_real,         &
                         grid_commx, ierror)
#else
      qdf_work(1:nlocalx*nlocaly*nlocalz,1) = qdf_real(1:nlocalx*nlocaly*nlocalz)
#endif
      !$omp end master
      !$omp barrier

      call fft3d(qdf_work, ftqdf, ftqdf, ftqdf, ftqdf_work, ftqdf_work,          &
                 ftqdf_work, work1, work2, nlocalx, nlocaly,                     &
                 nlocalz, nlocalx1, x_local1, ngrid(1), ngrid(2), ngrid(3),      &
                 nprocx, nprocy, nprocz, id, nthread, my_x_rank, x_start1,       &
                 x_end1, y_start, y_end, z_start, z_end, grid_commx, grid_commy, &
                 grid_commz)
    end if

    ! Energy calculation
    !
    iproc = my_z_rank + 1

    if (fft_scheme == FFT_2dalltoall) then

      !$omp barrier
      nix = (nlocalx1-1)/nprocy
      nix1 = nix+1
      niy = nlocaly/nprocz
      niz = nlocalz/nprocx

      !$omp do
      do iy = 1, nlocaly
        do ix = 1, nix1*nprocy
          iorg = nlocalz*(ix-1) + nlocalz*nix1*nprocy*(iy-1)
          elec_temp = 0.0_wp
          do iz = 1, nlocalz
            k = iz + iorg
            r_work(k) = real(c_work(iz,ix,iy),wp)**2 + imag(c_work(iz,ix,iy))**2
            c_work(iz,ix,iy) = c_work(iz,ix,iy) * cmplx(theta_local(k),0.0_wp,wp)
            tqq = r_work(k)
            tqq = tqq * theta_local(k) * vol_fact4
            elec_temp = elec_temp + tqq
          end do
          eelec(id+1) = eelec(id+1) + elec_temp
        end do
      end do

      if(my_x_rank == 0 .and. my_y_rank == 0) then

        !$omp do
        do iy = 1, nlocaly
          do iprocy = 0, nprocy-1
            ix = iprocy*nix1+1
            iorg = nlocalz*(ix-1) + nlocalz*nix1*nprocy*(iy-1)
            elec_temp = 0.0_wp
            do iz = 1, nlocalz
              k = iz + iorg
              tqq = r_work(k)
              tqq = tqq * theta_local(k) * vol_fact2
              elec_temp = elec_temp - tqq
            end do
            eelec(id+1) = eelec(id+1) + elec_temp
          end do
        end do
      end if

      if(my_x_rank == (nprocx-1) .and. my_y_rank == (nprocy-1)) then

        !$omp do
        do iy = 1, nlocaly
          do iprocy = 0, nprocy-1
            ix = (iprocy+1)*nix1-1
            iorg = nlocalz*(ix-1) + nlocalz*nix1*nprocy*(iy-1)
            elec_temp = 0.0_wp
            do iz = 1, nlocalz
              k = iz + iorg
              tqq = r_work(k)
              tqq = tqq * theta_local(k) * vol_fact2
              elec_temp = elec_temp - tqq
            end do
            eelec(id+1) = eelec(id+1) + elec_temp
          end do
        end do
      end if

    else if (fft_scheme == FFT_1dalltoall) then

      niy = nlocaly / nprocz
      do ix = 1, x_local1
        do iy = id+1, niy, nthread
          elec_temp = 0.0_wp
          vxx = 0.0_wp
          vyy = 0.0_wp
          vzz = 0.0_wp
          do iz = 1, ngrid(3)
            tqq = real(ftqdf_work2(iz,ix,iy),wp)**2 + &
                  imag(ftqdf_work2(iz,ix,iy))**2
            tqq = tqq * theta(iz,iy,ix,1) * vol_fact4
            elec_temp = elec_temp + tqq

            ! virial
            ! 
            vxx = vxx + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gx(ix) * gx(ix))
            vyy = vyy + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gy(iy) * gy(iy))
            vzz = vzz +                                                        &
                  tqq * (1.0_wp + vir_fact(iz,iy,ix) * gz(iz,1) * gz(iz,1))

          end do
          eelec(id+1) = eelec(id+1) + elec_temp
          virial(1,1,id+1) = virial(1,1,id+1) + vxx
          virial(2,2,id+1) = virial(2,2,id+1) + vyy
          virial(3,3,id+1) = virial(3,3,id+1) + vzz
        end do
      end do

      if (my_x_rank == 0) then
        do iy = id+1, niy, nthread
          elec_temp = 0.0_wp
          vxx = 0.0_wp
          vyy = 0.0_wp
          vzz = 0.0_wp
          do iz = 1, ngrid(3)
            tqq = real(ftqdf_work2(iz,1,iy),wp)**2 + &
                  imag(ftqdf_work2(iz,1,iy))**2
            tqq = tqq * theta(iz,iy,1,1) * vol_fact2
            elec_temp = elec_temp - tqq
            vxx = vxx - tqq * (1.0_wp + vir_fact(iz,iy,1) * gx(1) * gx(1))
            vyy = vyy - tqq * (1.0_wp + vir_fact(iz,iy,1) * gy(iy) * gy(iy))
            vzz = vzz - tqq * (1.0_wp + vir_fact(iz,iy,1) * gz(iz,1) * gz(iz,1))
          end do
          eelec(id+1) = eelec(id+1) + elec_temp
          virial(1,1,id+1) = virial(1,1,id+1) + vxx
          virial(2,2,id+1) = virial(2,2,id+1) + vyy
          virial(3,3,id+1) = virial(3,3,id+1) + vzz
        end do
      end if

      if (my_x_rank == (nprocx-1)) then
        ix = x_local1
        do iy = id+1, niy, nthread
          elec_temp = 0.0_wp
          vxx = 0.0_wp
          vyy = 0.0_wp
          vzz = 0.0_wp
          do iz = 1, ngrid(3)
            tqq = real(ftqdf_work2(iz,ix,iy),wp)**2 + &
                  imag(ftqdf_work2(iz,ix,iy))**2
            tqq = tqq * theta(iz,iy,ix,1) * vol_fact2
            elec_temp = elec_temp - tqq
            vxx = vxx - tqq * (1.0_wp + vir_fact(iz,iy,ix) * gx(ix) * gx(ix))
            vyy = vyy - tqq * (1.0_wp + vir_fact(iz,iy,ix) * gy(iy) * gy(iy))
            vzz = vzz                                                          &
                  - tqq * (1.0_wp + vir_fact(iz,iy,ix) * gz(iz,1) * gz(iz,1))
          end do
          eelec(id+1) = eelec(id+1) + elec_temp
          virial(1,1,id+1) = virial(1,1,id+1) + vxx
          virial(2,2,id+1) = virial(2,2,id+1) + vyy
          virial(3,3,id+1) = virial(3,3,id+1) + vzz
        end do
      end if


      ! F^-1[Th]*F^-1[Q] (=X)
      !
      !$omp barrier
      do ix = 1, x_local1
        do iy = id+1, niy, nthread
          k = (iy-1)*ngrid(3) + (ix-1)*ngrid(3)*niy
          do iz = 1, ngrid(3)
            ftqdf_work(k+iz,1) = ftqdf_work2(iz,ix,iy)   &
                               * cmplx(theta(iz,iy,ix,1),0.0_wp,wp)
          end do
        end do
      end do

    else if (fft_scheme == FFT_1dallgather) then

      do ix = 1, x_local1
        do iy = id+1, nlocaly, nthread
          elec_temp = 0.0_wp
          vxx = 0.0_wp
          vyy = 0.0_wp
          vzz = 0.0_wp
          do iz = 1, nlocalz
            k = iz + (iy-1)*nlocalz + (ix-1)*nlocalz*nlocaly
            tqq = real(ftqdf_work(k,iproc),wp)**2 + &
                  imag(ftqdf_work(k,iproc))**2
            tqq = tqq * theta(iz,iy,ix,iproc) * vol_fact4
            elec_temp = elec_temp + tqq

            ! virial
            ! 
            vxx = vxx + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gx(ix) * gx(ix))
            vyy = vyy + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gy(iy) * gy(iy))
            vzz = vzz + tqq                                                     &
                 * (1.0_wp + vir_fact(iz,iy,ix) * gz(iz,iproc) * gz(iz,iproc))
          end do
          eelec(id+1) = eelec(id+1) + elec_temp
          virial(1,1,id+1) = virial(1,1,id+1) + vxx
          virial(2,2,id+1) = virial(2,2,id+1) + vyy
          virial(3,3,id+1) = virial(3,3,id+1) + vzz
        end do
      end do

      if (my_x_rank == 0) then
        do iy = id+1, nlocaly, nthread
          elec_temp = 0.0_wp
          vxx = 0.0_wp
          vyy = 0.0_wp
          vzz = 0.0_wp
          do iz = 1, nlocalz
            k = iz + (iy-1)*nlocalz
            tqq = real(ftqdf_work(k,iproc),wp)**2 + imag(ftqdf_work(k,iproc))**2
            tqq = tqq * theta(iz,iy,1,iproc) * vol_fact2
            elec_temp = elec_temp - tqq
            vxx = vxx - tqq * (1.0_wp + vir_fact(iz,iy,1) * gx(1) * gx(1))
            vyy = vyy - tqq * (1.0_wp + vir_fact(iz,iy,1) * gy(iy) * gy(iy))
            vzz = vzz - tqq                                                    &
                 * (1.0_wp + vir_fact(iz,iy,1) * gz(iz,iproc) * gz(iz,iproc))
          end do
          eelec(id+1) = eelec(id+1) + elec_temp
          virial(1,1,id+1) = virial(1,1,id+1) + vxx
          virial(2,2,id+1) = virial(2,2,id+1) + vyy
          virial(3,3,id+1) = virial(3,3,id+1) + vzz
        end do
      end if

      if (my_x_rank == (nprocx-1)) then
        ix = x_local1
        do iy = id+1, nlocaly, nthread
          elec_temp = 0.0_wp
          vxx = 0.0_wp
          vyy = 0.0_wp
          vzz = 0.0_wp
          do iz = 1, nlocalz
            k = iz + (iy-1)*nlocalz + (x_local1-1)*nlocalz*nlocaly
            tqq = real(ftqdf_work(k,iproc),wp)**2 + &
                  imag(ftqdf_work(k,iproc))**2
            tqq = tqq * theta(iz,iy,ix,iproc) * vol_fact2
            elec_temp = elec_temp - tqq
            vxx = vxx - tqq * (1.0_wp + vir_fact(iz,iy,ix) * gx(ix) * gx(ix))
            vyy = vyy - tqq * (1.0_wp + vir_fact(iz,iy,ix) * gy(iy) * gy(iy))
            vzz = vzz - tqq                                                    &
                 * (1.0_wp + vir_fact(iz,iy,ix) * gz(iz,iproc) * gz(iz,iproc))
          end do
          eelec(id+1) = eelec(id+1) + elec_temp
          virial(1,1,id+1) = virial(1,1,id+1) + vxx
          virial(2,2,id+1) = virial(2,2,id+1) + vyy
          virial(3,3,id+1) = virial(3,3,id+1) + vzz
        end do
      end if

      ! F^-1[Th]*F^-1[Q] (=X)
      !
      !$omp barrier
      do ix = 1, x_local1
        do iy = id+1, niy, nthread
          k = (iy-1)*nlocalz + (ix-1)*nlocalz*nlocaly
          do iproc = 1, nprocz
            ftqdf_work(k+1:k+nlocalz,iproc) = ftqdf_work(k+1:k+nlocalz,iproc)  &
                                          * cmplx(theta(1:nlocalz,iy,ix,iproc),0.0_wp,wp)
          end do
        end do
      end do
      !$omp barrier

    end if

    ! FFT(backward)
    ! 
    if (fft_scheme == FFT_1dalltoall) then
      call bfft3d_1d_alltoall(qdf_work, qdf_real, ftqdf, ftqdf2, ftqdf3,   &
                  ftqdf_work3, ftqdf_work2, ftqdf_work, work1, work2,      &
                  nlocalx, nlocaly, nlocalz, nlocalx1, x_local1, ngrid(1), &
                  ngrid(2), ngrid(3), nprocx, nprocy, nprocz, id, nthread, &
                  my_x_rank, my_y_rank, my_z_rank, x_start1, x_end1,       &
                  y_start, y_end, z_start, z_end, grid_commx, grid_commy,  &
                  grid_commz, niy)
    else if (fft_scheme == FFT_2dalltoall) then
      call bfft3d_2d_alltoall(qdf_work, qdf_real, ftqdf_work, work1,       &
                  work2, nlocalx, nlocaly, nlocalz, nlocalx1, x_local1,    &
                  ngrid(1), ngrid(2), ngrid(3), nprocx, nprocy, nprocz,    &
                  id, nthread, my_x_rank, my_y_rank, my_z_rank, x_start1,  &
                  x_end1, y_start, y_end, z_start, z_end, grid_commx,      &
                  grid_commy, grid_commz, grid_commxy, c_work, c_work1,    &
                  c_work2,c_work3)
    else if (fft_scheme == FFT_1dallgather) then
      call bfft3d(qdf_work, qdf_real, ftqdf, ftqdf, ftqdf, ftqdf_work,     &
                  ftqdf_work, ftqdf_work, work1, work2, nlocalx, nlocaly,  &
                  nlocalz, nlocalx1, x_local1, ngrid(1), ngrid(2),         &
                  ngrid(3), nprocx, nprocy, nprocz, id, nthread,           &
                  my_x_rank, x_start1, x_end1, y_start, y_end, z_start,    &
                  z_end, grid_commx, grid_commy,  grid_commz)
    end if

    !$omp barrier

!#ifndef USE_GPU
    ! Gradient on CPU 
    !
    ! X is saved on qdf
    !
    do iz = id+1, nlocalz, nthread
      do iy = 1, nlocaly
        k = (iy-1)*nlocalx + (iz-1)*nlocalx*nlocaly
        do ix = 1, nlocalx
          qdf(ix,iy,iz,1) = qdf_real(k+ix)
        end do
      end do
    end do

    !$omp barrier

    do icel = id+1, ncell, nthread
      do i = 1, natom(icel)

        f_1 = 0.0_wp
        f_2 = 0.0_wp
        f_3 = 0.0_wp

        do k = 1, 3
          ii(k) = vi(k,i,icel)
        end do

        nix = 0
        niy = 0
        niz = 0

        do ixyz = 1, n_bs

          ixs = ii(1) - (ixyz-1) + 1
          iys = ii(2) - (ixyz-1) + 1
          izs = ii(3) - (ixyz-1) + 1
          if (ixs <= 0) ixs = ixs + ngrid(1)
          if (iys <= 0) iys = iys + ngrid(2)
          if (izs <= 0) izs = izs + ngrid(3)

          if (ixs >= x_start .and. ixs <= x_end) then
            nix = nix + 1
            vx_tmp(nix) = ixs - x_start + 1
            bsc_tmp(1,nix)  = bsc(ixyz,1,i,icel)
            bscd_tmp(1,nix) = bscd(ixyz,1,i,icel)
          end if

          if (iys >= y_start .and. iys <= y_end) then
            niy = niy + 1
            vy_tmp(niy) = iys - y_start + 1
            bsc_tmp(2,niy)  = bsc(ixyz,2,i,icel)
            bscd_tmp(2,niy) = bscd(ixyz,2,i,icel)
          end if

          if (izs >= z_start .and. izs <= z_end) then
            niz = niz + 1
            vz_tmp(niz) = izs - z_start + 1
            bsc_tmp(3,niz)  = bsc(ixyz,3,i,icel)
            bscd_tmp(3,niz) = bscd(ixyz,3,i,icel)
          end if
        end do

        if (nix*niy*niz == 0) cycle

        do iz = 1, niz
          izs = vz_tmp(iz)
          do iy = 1, niy
            iys = vy_tmp(iy)
            do ix = 1, nix
              ixs = vx_tmp(ix)

              f_1 = f_1                                  &
                  + bscd_tmp(1,ix)*bsc_tmp(2,iy)         &
                  * bsc_tmp(3,iz)*qdf(ixs,iys,izs,1)
              f_2 = f_2                                  &
                  + bsc_tmp(1,ix)*bscd_tmp(2,iy)         &
                  * bsc_tmp(3,iz)*qdf(ixs,iys,izs,1)
              f_3 = f_3                                  &
                  + bsc_tmp(1,ix)*bsc_tmp(2,iy)          &
                  * bscd_tmp(3,iz)*qdf(ixs,iys,izs,1)
            end do
          end do
        end do

        force(1,i,icel) = force(1,i,icel)                          &
             - f_1 * charge(i,icel) * r_scale(1) * vol_fact4 * bs_fact3d
        force(2,i,icel) = force(2,i,icel)                          &
             - f_2 * charge(i,icel) * r_scale(2) * vol_fact4 * bs_fact3d
        force(3,i,icel) = force(3,i,icel)                          &
             - f_3 * charge(i,icel) * r_scale(3) * vol_fact4 * bs_fact3d
 
      end do
    end do

!#else
!    !
!    ! Gradient on GPU 
!    !
!    !$omp master
!    call gpu_pme_recip_interpolate_force( &
!         force, qdf_real, charge, natom, &
!         MaxAtom, ncell, nlocalx, nlocaly, nlocalz, ngrid, &
!         x_start, x_end, y_start, y_end, z_start, z_end, &
!         r_scale, bs_fact3d, vol_fact4, &
!         vi, bsc, bscd )
!    !$omp end master
!
!    !$omp barrier
!
!#endif

    !$omp barrier

    deallocate(work1, work2)

    !$omp end parallel
 
    return

  end subroutine pme_recip

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_pre_fep
  !> @brief        Prepare functions for PME calculation for FEP calculation
  !! @authors      NK
  !! @param[in]    domain   : domain information
  !! @param[in]    boundary : boundary information
  !! @note         Extracted from setup_pme for NPT calculation
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_pre_fep(domain, boundary)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_boundary),        intent(in)    :: boundary
    ! local variables
    integer                  :: i, j, k, is, js, ks, ix, iy, iz, iproc
    integer                  :: ic, kk, iprocx, iprocy, ixs, izs, k_f, k_t
    integer                  :: nix, nix1, niy, niz, iorg
    integer                  :: ix_start, ix_end, iz_start
    integer                  :: iz_end, iy_start, iy_end
    real(wp)                 :: fact, gfact(3), g2


    box(1) = boundary%box_size_x
    box(2) = boundary%box_size_y
    box(3) = boundary%box_size_z

    do k = 1, 3
      box_inv(k) = 1.0_wp / box(k)
      r_scale(k) = real(ngrid(k),wp) * box_inv(k)
    end do

    vol_fact2 = 2.0_wp * PI * el_fact * box_inv(1)*box_inv(2)*box_inv(3)
    vol_fact4 = 2.0_wp * vol_fact2

    if (fft_scheme == FFT_2dalltoall) then
      theta_local(1:nlocalx*nlocaly*nlocalz) = 0.0_wp
    end if

    ! Prepareing theta = F^-1[theta](h), h shifted
    ! Gx, Gy, Gz, and vir_fact=-2(1+G^2/4a)/G^2 are also prepared 
    ! for virial calculation
    !
    do k = 1, 3
      gfact(k) = 2.0_wp * PI * box_inv(k)
    end do

    fact = 0.25_wp / alpha2m

    if (fft_scheme == FFT_2dalltoall) then

      do i = x_start1, x_end1

        ix = i - x_start1 + 1
        is = i - 1
        gx(ix) = gfact(1) * real(is,wp)

        do j = y_start, y_end

          iy = j - y_start + 1
          if (j <= ngrid(2)/2+1) then
            js = j - 1
          else
            js = j - 1 - ngrid(2)
          end if

          gy(iy) = gfact(2) * real(js,wp)

          do iproc = 1, nprocz
            do iz = 1, nlocalz
              k = (iproc-1)*nlocalz + iz
              if (k <= ngrid(3)/2+1) then
                ks = k - 1
              else
                ks = k - 1 - ngrid(3)
              end if
              gz(iz,iproc) = gfact(3) * real(ks,wp)
              g2 = gx(ix)**2 + gy(iy)**2 + gz(iz,iproc)**2
              if (g2 > EPS) then
                if (iproc == (my_z_rank+1)) &
                  vir_fact(iz,iy,ix) = -2.0_wp * (1.0_wp - g2 * fact)/g2
                if (iproc == (my_z_rank+1)) then
                  kk = ix + (iy-1)*nlocalx1 + (iz-1)*nlocalx1*nlocaly
                  theta_local(kk) = b2(i,1) * b2(j,2) * b2(k,3)             & 
                                    * exp(g2 * fact)/g2
                endif
              end if

            end do

          end do

        end do

      end do

      if (my_x_rank == 0 .and. my_y_rank == 0 .and. my_z_rank == 0) &
        theta_local(1) = 0.0d0

      nix = (nlocalx1-1)/nprocy
      nix1 = nix+1
      niy = nlocaly/nprocz
      niz = nlocalz/nprocx

      call mpi_alltoall(theta_local, nlocalx1*nlocaly*niz, mpi_wp_real, &
                        r_work1,     nlocalx1*nlocaly*niz, mpi_wp_real, &
                        grid_commx, ierror)

      do iz = 1, nlocalz
        do iy = 1, nlocaly
          do ix = 1, nlocalx1
            k_f = ix + (iy-1)*nlocalx1 + (iz-1)*nlocalx1*nlocaly
            k_t = iy + (ix-1)*nlocaly  + (iz-1)*nlocalx1*nlocaly
            r_work2(k_t) = r_work1(k_f)
          enddo
        enddo
      enddo

      r_work1 = 0.0_wp
      do iprocx = 0, nprocx-1
        iz_start = iprocx*niz+1
        iz_end = iz_start+niz-1
        do iprocy = 0, nprocy-1
          ix_start = iprocy*nix+1
          ix_end   = ix_start+nix-1
          if(iprocy == nprocy-1) ix_end = ix_start+nix1-1
          iorg = nlocaly*nix1*niz*(iprocx+nprocx*iprocy)
          do iz = iz_start, iz_end
            izs = iz - iz_start
            do ix = ix_start, ix_end
              ixs = ix - ix_start
              k = iorg+nlocaly*ixs+nlocaly*nix1*izs
              do iy = 1, nlocaly
                kk = iy + (ix-1)*nlocaly  + (iz-1)*nlocalx1*nlocaly
                r_work1(k+iy) = r_work2(kk)
              enddo
            enddo
          enddo
        enddo
      enddo

      call mpi_alltoall(r_work1, nlocaly*nix1*niz, mpi_wp_real, &
                        r_work2, nlocaly*nix1*niz, mpi_wp_real, &
                        grid_commxy, ierror)

      do iprocx = 0, nprocx-1
        iz_start = iprocx*niz+1
        iz_end = iz_start+niz-1
        do iprocy = 0, nprocy-1
          ix_start = iprocy*nix1+1
          ix_end   = ix_start+nix1-1
          iorg = nlocaly*nix1*niz*(iprocx+nprocx*iprocy)
          do iz = iz_start, iz_end
            izs = iz - iz_start
            do ix = ix_start, ix_end
              ixs = ix - ix_start
              k = iorg+nlocaly*ixs+nlocaly*nix1*izs
              do iy = 1, nlocaly
                r_work3(iy,ix,iz) = r_work2(k+iy)
              enddo
            enddo
          enddo
        enddo
      enddo

      do iz = 1, nlocalz
        do ix = 1, nix1
          do iproc = 1, nprocy
            k = (iproc-1)*nlocaly
            ixs = (iproc-1)*nix1 + ix
            do iy = 1, nlocaly
              kk = iz + nlocalz*(ixs-1) + nlocalz*nix1*nprocy*(iy-1)
              r_work2(kk) = r_work3(iy,ixs,iz)
            end do
          end do
        end do
      end do

      call mpi_alltoall(r_work2, nix1*nprocy*niy*nlocalz, mpi_wp_real, &
                        r_work1, nix1*nprocy*niy*nlocalz, mpi_wp_real, &
                        grid_commz, ierror)

      do iy = 1,nlocaly
        do ix = 1,nix1*nprocy
          do iz = 1,nlocalz
             kk = iz + nlocalz*(ix-1) + nlocalz*nix1*nprocy*(iy-1)
             theta_local(kk) = r_work1(kk)
           enddo
        enddo
      enddo
      if (my_x_rank == 0 .and. my_y_rank == 0) &
        theta(1,1,1,1) = 0.0_wp

    else if (fft_scheme == FFT_1dalltoall) then

      niy = nlocaly / nprocz
      iy_start = y_start + niy*my_z_rank
      iy_end   = iy_start + niy - 1

      do i = x_start1, x_end1

        ix = i - x_start1 + 1
        is = i - 1
        gx(ix) = gfact(1) * real(is,wp)

        do j = iy_start, iy_end

          iy = j - iy_start + 1

          if (j <= ngrid(2)/2+1) then
            js = j - 1
          else
            js = j - 1 - ngrid(2)
          end if
          gy(iy) = gfact(2) * real(js,wp)

          do k = 1, ngrid(3)
            if (k <= ngrid(3)/2+1) then
              ks = k - 1
            else
              ks = k - 1 - ngrid(3)
            end if
            gz(k,1) = gfact(3) * real(ks,wp)

            g2 = gx(ix)**2 + gy(iy)**2 + gz(k,1)**2
            if (g2 > EPS) then
              vir_fact(k,iy,ix) = -2.0_wp * (1.0_wp - g2 * fact) / g2
              theta(k,iy,ix,1) = b2(i,1)*b2(j,2)*b2(k,3)*exp(g2*fact) / g2
            end if

          end do

        end do

      end do
      if (my_x_rank == 0 .and. my_y_rank == 0 .and. my_z_rank == 0) &
        theta(1,1,1,1) = 0.0_wp

    else

      do i = x_start1, x_end1

        ix = i - x_start1 + 1
        is = i - 1
        gx(ix) = gfact(1) * real(is,wp)

        do j = y_start, y_end

          iy = j - y_start + 1

          if (j <= ngrid(2)/2+1) then
            js = j - 1
          else
            js = j - 1 - ngrid(2)
          end if

          gy(iy) = gfact(2) * real(js,wp)

          do iproc = 1, nprocz
            do iz = 1, nlocalz
              k = (iproc-1)*nlocalz + iz
              if (k <= ngrid(3)/2+1) then
                ks = k - 1
              else
                ks = k - 1 - ngrid(3)
              end if
              gz(iz,iproc) = gfact(3) * real(ks,wp)
              g2 = gx(ix)**2 + gy(iy)**2 + gz(iz,iproc)**2
              if (g2 > EPS) then
                if (iproc == (my_z_rank+1)) &
                  vir_fact(iz,iy,ix) = -2.0_wp * (1.0_wp - g2 * fact)/g2
                theta(iz,iy,ix,iproc) = &
                     b2(i,1) * b2(j,2) * b2(k,3) * exp(g2 * fact)/g2
              end if

            end do

          end do

        end do

      end do
      if (my_x_rank == 0 .and. my_y_rank == 0) &
        theta(1,1,1,1) = 0.0_wp

    end if

    if (my_city_rank == 0) &
      vir_fact(1,1,1) = 0.0_wp

    ! Calculating self energy
    !
    u_self_preserve = 0.0_dp
    u_self_appear   = 0.0_dp
    u_self_vanish   = 0.0_dp

    do i = 1, domain%num_cell_local
      do ix = 1, domain%num_atom_preserve(i)
        u_self_preserve = u_self_preserve &
                          + domain%charge(domain%pmelist_preserve(ix,i),i)**2
      end do
      do ix = 1, domain%num_atom_appear_gr(i)
        u_self_appear = u_self_appear &
                        + domain%charge(domain%pmelist_appear_gr(ix,i),i)**2
      end do
      do ix = 1, domain%num_atom_vanish_gr(i)
        u_self_vanish = u_self_vanish &
                        + domain%charge(domain%pmelist_vanish_gr(ix,i),i)**2
      end do
    end do

    u_self_preserve = - u_self_preserve * el_fact * alpha/sqrt(PI)
    u_self_appear   = - u_self_appear * el_fact * alpha/sqrt(PI)
    u_self_vanish   = - u_self_vanish * el_fact * alpha/sqrt(PI)


    return

  end subroutine pme_pre_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pme_recip_fep
  !> @brief        Calculate PME reciprocal part with domain decomposition
  !                for FEP calculation
  !! @authors      NK
  !! @param[in]    domain : domain information
  !! @param[in]    enefunc: enefunc information
  !! @param[in]    flg_fep: flag for FEP calculations
  !! @param[inout] force  : forces of target systems
  !! @param[inout] virial : virial term of target systems
  !! @param[inout] eelec  : electrostatic energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pme_recip_fep(domain, enefunc, flg_fep, force, virial, eelec)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    integer,                 intent(in)    :: flg_fep
    real(dp),                intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: virial(3,3,nthread)
    real(dp),                intent(inout) :: eelec(nthread)

    ! local variables
    real(wp)                 :: elec_temp
    real(wp)                 :: vr(3), dv(3), bsc_tmp(3,n_bs), bscd_tmp(3,n_bs)
    real(wp)                 :: vxx, vyy, vzz, temp
    real(wp)                 :: tqq
    real(wp)                 :: f_1, f_2, f_3
    integer                  :: i, k, k1, icel, iproc, id
    integer                  :: ix, iv, iy, iz, ixs, iys, izs, ip
    integer                  :: ixyz, nix, nix1, niy, niz
    integer                  :: ixx, iyy, izz, ii(3)
    integer                  :: vx_tmp(n_bs), vy_tmp(n_bs), vz_tmp(n_bs)
    integer                  :: ncell, omp_get_thread_num
    integer                  :: kk, iorg, k_f, k_t
    integer                  :: iprocx, iprocy
    integer                  :: ix_start, ix_end, iz_start, iz_end

    complex(wp), allocatable :: work1(:), work2(:)
    real(dp),        pointer :: coord(:,:,:)
    real(wp),        pointer :: charge(:,:)
    integer,         pointer :: natom(:)
    integer,         pointer :: pmelist(:,:)

    real(dp)                 :: eelec_temp(nthread)
    real(wp)                 :: lambel


    coord  => domain%coord
    charge => domain%charge
    ncell  = domain%num_cell_local + domain%num_cell_boundary

    select case (flg_fep)

    case(FEP_PRESERVE)
      natom   => domain%num_atom_preserve
      pmelist => domain%pmelist_preserve
      lambel  =  1.0_wp - enefunc%lambelB - enefunc%lambelA

    case(FEP_APPEAR)
      natom   => domain%num_atom_appear_gr
      pmelist => domain%pmelist_appear_gr
      lambel  =  enefunc%lambelB

    case(FEP_VANISH)
      natom   => domain%num_atom_vanish_gr
      pmelist => domain%pmelist_vanish_gr
      lambel  =  enefunc%lambelA

    end select

    ! Initializing the energy and force
    !
    qdf_real(1:nlocalx*nlocaly*nlocalz) = 0.0_wp
    eelec_temp(1:nthread) = 0.0_dp

    !$omp parallel default(shared)                                             &
    !$omp private(id, i, iv, ix, iy, iz, icel, k, ii, dv, izz, izs, iyy, iys,  &
    !$omp         iz_start, iz_end, ix_start, ix_end, iorg, k_f, k_t,          &
    !$omp         ixx, ixs, vxx, vyy, vzz, iproc, tqq, k1, temp, vr, &
    !$omp         work1, work2, elec_temp, f_1, f_2, f_3, ixyz, vx_tmp, vy_tmp,&
    !$omp         vz_tmp, nix, nix1, niy, niz, bsc_tmp, bscd_tmp, kk, iprocx,  &
    !$omp         iprocy, ip)
    !

#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    !$omp barrier

    ! initialization
    !
    allocate(work1(ngridmax), work2(2*ngridmax))

    ! Initialization of Q-fanction
    !
!#ifndef USE_GPU
    do iz = 1, nlocalz
      do iy = 1, nlocaly
        do ix = 1, nlocalx
          qdf(ix,iy,iz,id+1) = 0.0_wp
        end do
      end do
    end do


    ! Calculateing Q-fanction
    !
    do icel = id+1, ncell, nthread
      do i = 1, natom(icel)
        ip = pmelist(i,icel)

        vr(1:3) = coord(1:3,ip,icel) * r_scale(1:3)
        vr(1:3) = vr(1:3) + real(ngrid(1:3)/2,wp) &
                - real(ngrid(1:3),wp)*anint(vr(1:3)/real(ngrid(1:3),wp))

        do k = 1, 3
          ii(k) = int(vr(k))
          if (ii(k) >= ngrid(k)) then
            vr(k) = real(ngrid(k),wp) - 0.00001_wp
            ii(k) = ngrid(k)-1
          endif
          vi(k,ip,icel) = ii(k)
          dv(k) = vr(k) - real(ii(k),wp)
          call b_spline_dev_coef(n_bs, dv(k), bsc(1,k,ip,icel), &
                                   bscd(1,k,ip,icel))
        end do

        nix = 0 
        niy = 0 
        niz = 0 

        do ixyz = 1, n_bs
          ixs = ii(1) - ixyz + 2
          iys = ii(2) - ixyz + 2
          izs = ii(3) - ixyz + 2
          if (ixs <= 0) ixs = ixs + ngrid(1)
          if (iys <= 0) iys = iys + ngrid(2)
          if (izs <= 0) izs = izs + ngrid(3)

          if (ixs >= x_start .and. ixs <= x_end) then
            nix = nix + 1
            vx_tmp(nix) = ixs - x_start + 1
            bsc_tmp(1,nix) = bsc(ixyz,1,ip,icel)
          end if

          if (iys >= y_start .and. iys <= y_end) then
            niy = niy + 1
            vy_tmp(niy) = iys - y_start + 1
            bsc_tmp(2,niy) = bsc(ixyz,2,ip,icel)
          end if

          if (izs >= z_start .and. izs <= z_end) then
            niz = niz + 1
            vz_tmp(niz) = izs - z_start + 1
            bsc_tmp(3,niz) = bsc(ixyz,3,ip,icel)
          end if
        end do

        do iz = 1, niz
          izs = vz_tmp(iz)
          do iy = 1, niy
            iys = vy_tmp(iy)
            do ix = 1, nix
              ixs = vx_tmp(ix)
              qdf(ixs,iys,izs,id+1) = qdf(ixs,iys,izs,id+1)           &
                                    + bsc_tmp(1,ix)*bsc_tmp(2,iy)     &
                                     * bsc_tmp(3,iz)*charge(ip,icel)  &
                                     * bs_fact3
            end do
          end do
        end do
      end do
    end do

    !$omp barrier
    do iproc = 1, nthread
      !$omp do
      do iz = 1, nlocalz
        do iy = 1, nlocaly
          k = (iy-1)*nlocalx + (iz-1)*nlocalx*nlocaly
          do iv = 1, nlocalx, 16
            do ix = iv, min0(iv+15,nlocalx)
              qdf_real(k+ix) = qdf_real(k+ix) + qdf(ix,iy,iz,iproc)
            end do
          end do
        end do
      end do
    end do

    !$omp barrier

!#else
!    !$omp master
!    call gpu_pme_recip_build_qdf( qdf_real, coord, charge, natom, &
!         MaxAtom, ncell, nlocalx, nlocaly, nlocalz, ngrid, &
!         x_start, x_end, y_start, y_end, z_start, z_end, &
!         r_scale, &
!         bs_fact3, &
!         vi, bsc, bscd )
!    !$omp end master
!#endif

    ! FFT (forward)
    !
    if ( (fft_scheme == FFT_1dalltoall) .or.   &
         (fft_scheme == FFT_2dalltoall) ) then

      !$omp barrier
      !$omp master
#ifdef HAVE_MPI_GENESIS
      call mpi_alltoall(qdf_real, nlocalx*nlocaly*nlocalz/nprocx, mpi_wp_real, &
                        qdf_work, nlocalx*nlocaly*nlocalz/nprocx, mpi_wp_real, &
                        grid_commx, ierror)
#else
      qdf_work(1:nlocalx*nlocaly*nlocalz,1) = qdf_real(1:nlocalx*nlocaly*nlocalz)
#endif
      !$omp end master
      !$omp barrier

      if (fft_scheme == FFT_1dalltoall) then
        call fft3d_1d_alltoall(qdf_work, ftqdf, ftqdf2, ftqdf3, ftqdf4, ftqdf5,    &
                         ftqdf_work, ftqdf_work, ftqdf_work, ftqdf_work,           &
                         ftqdf_work2, ftqdf_work3, work1, work2, nlocalx, nlocaly, &
                         nlocalz, nlocalx1, x_local1, ngrid(1), ngrid(2), ngrid(3),&
                         nprocx, nprocy, nprocz, id, nthread, my_x_rank, my_y_rank,&
                         my_z_rank, x_start1, x_end1, y_start, y_end, z_start,     &
                         z_end, grid_commx, grid_commy, grid_commz)

      else
        call fft3d_2d_alltoall(qdf_work, ftqdf_work, work1, work2, nlocalx,        &
                         nlocaly, nlocalz, nlocalx1, x_local1, ngrid(1), ngrid(2), &
                         ngrid(3), nprocx, nprocy, nprocz, id, nthread, my_x_rank, &
                         my_y_rank, my_z_rank, x_start1, x_end1, y_start, y_end,   &
                         z_start, z_end, grid_commx, grid_commy, grid_commz,       &
                         grid_commxy, c_work, c_work1, c_work2, c_work3)
      end if

    else if (fft_scheme == FFT_1dallgather) then

      !$omp barrier
      !$omp master
#ifdef HAVE_MPI_GENESIS
      call mpi_allgather(qdf_real, nlocalx*nlocaly*nlocalz, mpi_wp_real,         &
                         qdf_work, nlocalx*nlocaly*nlocalz, mpi_wp_real,         &
                         grid_commx, ierror)
#else
      qdf_work(1:nlocalx*nlocaly*nlocalz,1) = qdf_real(1:nlocalx*nlocaly*nlocalz)
#endif
      !$omp end master
      !$omp barrier

      call fft3d(qdf_work, ftqdf, ftqdf, ftqdf, ftqdf_work, ftqdf_work,          &
                 ftqdf_work, work1, work2, nlocalx, nlocaly,                     &
                 nlocalz, nlocalx1, x_local1, ngrid(1), ngrid(2), ngrid(3),      &
                 nprocx, nprocy, nprocz, id, nthread, my_x_rank, x_start1,       &
                 x_end1, y_start, y_end, z_start, z_end, grid_commx, grid_commy, &
                 grid_commz)
    end if

    ! Energy calculation
    !
    iproc = my_z_rank + 1

    if (fft_scheme == FFT_2dalltoall) then

      !$omp barrier
      nix = (nlocalx1-1)/nprocy
      nix1 = nix+1
      niy = nlocaly/nprocz
      niz = nlocalz/nprocx

      !$omp do
      do iy = 1, nlocaly
        do ix = 1, nix1*nprocy
          iorg = nlocalz*(ix-1) + nlocalz*nix1*nprocy*(iy-1)
          elec_temp = 0.0_wp
          do iz = 1, nlocalz
            k = iz + iorg
            r_work(k) = real(c_work(iz,ix,iy),wp)**2 + imag(c_work(iz,ix,iy))**2
            c_work(iz,ix,iy) = c_work(iz,ix,iy) * cmplx(theta_local(k),0.0_wp,wp)
            tqq = r_work(k)
            tqq = tqq * theta_local(k) * vol_fact4
            elec_temp = elec_temp + tqq
          end do
          eelec_temp(id+1) = eelec_temp(id+1) + elec_temp
        end do
      end do

      if(my_x_rank == 0 .and. my_y_rank == 0) then

        !$omp do
        do iy = 1, nlocaly
          do iprocy = 0, nprocy-1
            ix = iprocy*nix1+1
            iorg = nlocalz*(ix-1) + nlocalz*nix1*nprocy*(iy-1)
            elec_temp = 0.0_wp
            do iz = 1, nlocalz
              k = iz + iorg
              tqq = r_work(k)
              tqq = tqq * theta_local(k) * vol_fact2
              elec_temp = elec_temp - tqq
            end do
            eelec_temp(id+1) = eelec_temp(id+1) + elec_temp
          end do
        end do
      end if

      if(my_x_rank == (nprocx-1) .and. my_y_rank == (nprocy-1)) then

        !$omp do
        do iy = 1, nlocaly
          do iprocy = 0, nprocy-1
            ix = (iprocy+1)*nix1-1
            iorg = nlocalz*(ix-1) + nlocalz*nix1*nprocy*(iy-1)
            elec_temp = 0.0_wp
            do iz = 1, nlocalz
              k = iz + iorg
              tqq = r_work(k)
              tqq = tqq * theta_local(k) * vol_fact2
              elec_temp = elec_temp - tqq
            end do
            eelec_temp(id+1) = eelec_temp(id+1) + elec_temp
          end do
        end do
      end if

    else if (fft_scheme == FFT_1dalltoall) then

      niy = nlocaly / nprocz
      do ix = 1, x_local1
        do iy = id+1, niy, nthread
          elec_temp = 0.0_wp
          vxx = 0.0_wp
          vyy = 0.0_wp
          vzz = 0.0_wp
          do iz = 1, ngrid(3)
            tqq = real(ftqdf_work2(iz,ix,iy),wp)**2 + &
                  imag(ftqdf_work2(iz,ix,iy))**2
            tqq = tqq * theta(iz,iy,ix,1) * vol_fact4
            elec_temp = elec_temp + tqq
            tqq = tqq * lambel

            ! virial
            ! 
            vxx = vxx + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gx(ix) * gx(ix))
            vyy = vyy + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gy(iy) * gy(iy))
            vzz = vzz +                                                        &
                  tqq * (1.0_wp + vir_fact(iz,iy,ix) * gz(iz,1) * gz(iz,1))

          end do
          eelec_temp(id+1) = eelec_temp(id+1) + elec_temp
          virial(1,1,id+1) = virial(1,1,id+1) + vxx
          virial(2,2,id+1) = virial(2,2,id+1) + vyy
          virial(3,3,id+1) = virial(3,3,id+1) + vzz
        end do
      end do

      if (my_x_rank == 0) then
        do iy = id+1, niy, nthread
          elec_temp = 0.0_wp
          vxx = 0.0_wp
          vyy = 0.0_wp
          vzz = 0.0_wp
          do iz = 1, ngrid(3)
            tqq = real(ftqdf_work2(iz,1,iy),wp)**2 + &
                  imag(ftqdf_work2(iz,1,iy))**2
            tqq = tqq * theta(iz,iy,1,1) * vol_fact2
            elec_temp = elec_temp - tqq
            tqq = tqq * lambel
            vxx = vxx - tqq * (1.0_wp + vir_fact(iz,iy,1) * gx(1) * gx(1))
            vyy = vyy - tqq * (1.0_wp + vir_fact(iz,iy,1) * gy(iy) * gy(iy))
            vzz = vzz - tqq * (1.0_wp + vir_fact(iz,iy,1) * gz(iz,1) * gz(iz,1))
          end do
          eelec_temp(id+1) = eelec_temp(id+1) + elec_temp
          virial(1,1,id+1) = virial(1,1,id+1) + vxx
          virial(2,2,id+1) = virial(2,2,id+1) + vyy
          virial(3,3,id+1) = virial(3,3,id+1) + vzz
        end do
      end if

      if (my_x_rank == (nprocx-1)) then
        ix = x_local1
        do iy = id+1, niy, nthread
          elec_temp = 0.0_wp
          vxx = 0.0_wp
          vyy = 0.0_wp
          vzz = 0.0_wp
          do iz = 1, ngrid(3)
            tqq = real(ftqdf_work2(iz,ix,iy),wp)**2 + &
                  imag(ftqdf_work2(iz,ix,iy))**2
            tqq = tqq * theta(iz,iy,ix,1) * vol_fact2
            elec_temp = elec_temp - tqq
            tqq = tqq * lambel
            vxx = vxx - tqq * (1.0_wp + vir_fact(iz,iy,ix) * gx(ix) * gx(ix))
            vyy = vyy - tqq * (1.0_wp + vir_fact(iz,iy,ix) * gy(iy) * gy(iy))
            vzz = vzz                                                          &
                  - tqq * (1.0_wp + vir_fact(iz,iy,ix) * gz(iz,1) * gz(iz,1))
          end do
          eelec_temp(id+1) = eelec_temp(id+1) + elec_temp
          virial(1,1,id+1) = virial(1,1,id+1) + vxx
          virial(2,2,id+1) = virial(2,2,id+1) + vyy
          virial(3,3,id+1) = virial(3,3,id+1) + vzz
        end do
      end if


      ! F^-1[Th]*F^-1[Q] (=X)
      !
      !$omp barrier
      do ix = 1, x_local1
        do iy = id+1, niy, nthread
          k = (iy-1)*ngrid(3) + (ix-1)*ngrid(3)*niy
          do iz = 1, ngrid(3)
            ftqdf_work(k+iz,1) = ftqdf_work2(iz,ix,iy)   &
                               * cmplx(theta(iz,iy,ix,1),0.0_wp,wp)
          end do
        end do
      end do

    else if (fft_scheme == FFT_1dallgather) then

      do ix = 1, x_local1
        do iy = id+1, nlocaly, nthread
          elec_temp = 0.0_wp
          vxx = 0.0_wp
          vyy = 0.0_wp
          vzz = 0.0_wp
          do iz = 1, nlocalz
            k = iz + (iy-1)*nlocalz + (ix-1)*nlocalz*nlocaly
            tqq = real(ftqdf_work(k,iproc),wp)**2 + &
                  imag(ftqdf_work(k,iproc))**2
            tqq = tqq * theta(iz,iy,ix,iproc) * vol_fact4
            elec_temp = elec_temp + tqq
            tqq = tqq * lambel

            ! virial
            ! 
            vxx = vxx + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gx(ix) * gx(ix))
            vyy = vyy + tqq * (1.0_wp + vir_fact(iz,iy,ix) * gy(iy) * gy(iy))
            vzz = vzz + tqq                                                     &
                 * (1.0_wp + vir_fact(iz,iy,ix) * gz(iz,iproc) * gz(iz,iproc))
          end do
          eelec_temp(id+1) = eelec_temp(id+1) + elec_temp
          virial(1,1,id+1) = virial(1,1,id+1) + vxx
          virial(2,2,id+1) = virial(2,2,id+1) + vyy
          virial(3,3,id+1) = virial(3,3,id+1) + vzz
        end do
      end do

      if (my_x_rank == 0) then
        do iy = id+1, nlocaly, nthread
          elec_temp = 0.0_wp
          vxx = 0.0_wp
          vyy = 0.0_wp
          vzz = 0.0_wp
          do iz = 1, nlocalz
            k = iz + (iy-1)*nlocalz
            tqq = real(ftqdf_work(k,iproc),wp)**2 + imag(ftqdf_work(k,iproc))**2
            tqq = tqq * theta(iz,iy,1,iproc) * vol_fact2
            elec_temp = elec_temp - tqq
            tqq = tqq * lambel
            vxx = vxx - tqq * (1.0_wp + vir_fact(iz,iy,1) * gx(1) * gx(1))
            vyy = vyy - tqq * (1.0_wp + vir_fact(iz,iy,1) * gy(iy) * gy(iy))
            vzz = vzz - tqq                                                    &
                 * (1.0_wp + vir_fact(iz,iy,1) * gz(iz,iproc) * gz(iz,iproc))
          end do
          eelec_temp(id+1) = eelec_temp(id+1) + elec_temp
          virial(1,1,id+1) = virial(1,1,id+1) + vxx
          virial(2,2,id+1) = virial(2,2,id+1) + vyy
          virial(3,3,id+1) = virial(3,3,id+1) + vzz
        end do
      end if

      if (my_x_rank == (nprocx-1)) then
        ix = x_local1
        do iy = id+1, nlocaly, nthread
          elec_temp = 0.0_wp
          vxx = 0.0_wp
          vyy = 0.0_wp
          vzz = 0.0_wp
          do iz = 1, nlocalz
            k = iz + (iy-1)*nlocalz + (x_local1-1)*nlocalz*nlocaly
            tqq = real(ftqdf_work(k,iproc),wp)**2 + &
                  imag(ftqdf_work(k,iproc))**2
            tqq = tqq * theta(iz,iy,ix,iproc) * vol_fact2
            elec_temp = elec_temp - tqq
            tqq = tqq * lambel
            vxx = vxx - tqq * (1.0_wp + vir_fact(iz,iy,ix) * gx(ix) * gx(ix))
            vyy = vyy - tqq * (1.0_wp + vir_fact(iz,iy,ix) * gy(iy) * gy(iy))
            vzz = vzz - tqq                                                    &
                 * (1.0_wp + vir_fact(iz,iy,ix) * gz(iz,iproc) * gz(iz,iproc))
          end do
          eelec_temp(id+1) = eelec_temp(id+1) + elec_temp
          virial(1,1,id+1) = virial(1,1,id+1) + vxx
          virial(2,2,id+1) = virial(2,2,id+1) + vyy
          virial(3,3,id+1) = virial(3,3,id+1) + vzz
        end do
      end if

      ! F^-1[Th]*F^-1[Q] (=X)
      !
      !$omp barrier
      do ix = 1, x_local1
        do iy = id+1, niy, nthread
          k = (iy-1)*nlocalz + (ix-1)*nlocalz*nlocaly
          do iproc = 1, nprocz
            ftqdf_work(k+1:k+nlocalz,iproc) = ftqdf_work(k+1:k+nlocalz,iproc)  &
                                          * cmplx(theta(1:nlocalz,iy,ix,iproc),0.0_wp,wp)
          end do
        end do
      end do
      !$omp barrier

    end if

    ! FFT(backward)
    ! 
    if (fft_scheme == FFT_1dalltoall) then
      call bfft3d_1d_alltoall(qdf_work, qdf_real, ftqdf, ftqdf2, ftqdf3,   &
                  ftqdf_work3, ftqdf_work2, ftqdf_work, work1, work2,      &
                  nlocalx, nlocaly, nlocalz, nlocalx1, x_local1, ngrid(1), &
                  ngrid(2), ngrid(3), nprocx, nprocy, nprocz, id, nthread, &
                  my_x_rank, my_y_rank, my_z_rank, x_start1, x_end1,       &
                  y_start, y_end, z_start, z_end, grid_commx, grid_commy,  &
                  grid_commz, niy)
    else if (fft_scheme == FFT_2dalltoall) then
      call bfft3d_2d_alltoall(qdf_work, qdf_real, ftqdf_work, work1,       &
                  work2, nlocalx, nlocaly, nlocalz, nlocalx1, x_local1,    &
                  ngrid(1), ngrid(2), ngrid(3), nprocx, nprocy, nprocz,    &
                  id, nthread, my_x_rank, my_y_rank, my_z_rank, x_start1,  &
                  x_end1, y_start, y_end, z_start, z_end, grid_commx,      &
                  grid_commy, grid_commz, grid_commxy, c_work, c_work1,    &
                  c_work2,c_work3)
    else if (fft_scheme == FFT_1dallgather) then
      call bfft3d(qdf_work, qdf_real, ftqdf, ftqdf, ftqdf, ftqdf_work,     &
                  ftqdf_work, ftqdf_work, work1, work2, nlocalx, nlocaly,  &
                  nlocalz, nlocalx1, x_local1, ngrid(1), ngrid(2),         &
                  ngrid(3), nprocx, nprocy, nprocz, id, nthread,           &
                  my_x_rank, x_start1, x_end1, y_start, y_end, z_start,    &
                  z_end, grid_commx, grid_commy,  grid_commz)
    end if

    !$omp barrier

!#ifndef USE_GPU
    ! Gradient on CPU 
    !
    ! X is saved on qdf
    !
    do iz = id+1, nlocalz, nthread
      do iy = 1, nlocaly
        k = (iy-1)*nlocalx + (iz-1)*nlocalx*nlocaly
        do ix = 1, nlocalx
          qdf(ix,iy,iz,1) = qdf_real(k+ix)
        end do
      end do
    end do

    !$omp barrier

    do icel = id+1, ncell, nthread
      do i = 1, natom(icel)

        ip = pmelist(i,icel)

        f_1 = 0.0_wp
        f_2 = 0.0_wp
        f_3 = 0.0_wp

        do k = 1, 3
          ii(k) = vi(k,ip,icel)
        end do

        nix = 0
        niy = 0
        niz = 0

        do ixyz = 1, n_bs

          ixs = ii(1) - (ixyz-1) + 1
          iys = ii(2) - (ixyz-1) + 1
          izs = ii(3) - (ixyz-1) + 1
          if (ixs <= 0) ixs = ixs + ngrid(1)
          if (iys <= 0) iys = iys + ngrid(2)
          if (izs <= 0) izs = izs + ngrid(3)

          if (ixs >= x_start .and. ixs <= x_end) then
            nix = nix + 1
            vx_tmp(nix) = ixs - x_start + 1
            bsc_tmp(1,nix)  = bsc(ixyz,1,ip,icel)
            bscd_tmp(1,nix) = bscd(ixyz,1,ip,icel)
          end if

          if (iys >= y_start .and. iys <= y_end) then
            niy = niy + 1
            vy_tmp(niy) = iys - y_start + 1
            bsc_tmp(2,niy)  = bsc(ixyz,2,ip,icel)
            bscd_tmp(2,niy) = bscd(ixyz,2,ip,icel)
          end if

          if (izs >= z_start .and. izs <= z_end) then
            niz = niz + 1
            vz_tmp(niz) = izs - z_start + 1
            bsc_tmp(3,niz)  = bsc(ixyz,3,ip,icel)
            bscd_tmp(3,niz) = bscd(ixyz,3,ip,icel)
          end if
        end do

        if (nix*niy*niz == 0) cycle

        do iz = 1, niz
          izs = vz_tmp(iz)
          do iy = 1, niy
            iys = vy_tmp(iy)
            do ix = 1, nix
              ixs = vx_tmp(ix)

              f_1 = f_1                                  &
                  + bscd_tmp(1,ix)*bsc_tmp(2,iy)         &
                  * bsc_tmp(3,iz)*qdf(ixs,iys,izs,1)
              f_2 = f_2                                  &
                  + bsc_tmp(1,ix)*bscd_tmp(2,iy)         &
                  * bsc_tmp(3,iz)*qdf(ixs,iys,izs,1)
              f_3 = f_3                                  &
                  + bsc_tmp(1,ix)*bsc_tmp(2,iy)          &
                  * bscd_tmp(3,iz)*qdf(ixs,iys,izs,1)
            end do
          end do
        end do

        force(1,ip,icel) = force(1,ip,icel) &
             - lambel * f_1 * charge(ip,icel) &
             * r_scale(1) * vol_fact4 * bs_fact3d
        force(2,ip,icel) = force(2,ip,icel) &
             - lambel * f_2 * charge(ip,icel) &
             * r_scale(2) * vol_fact4 * bs_fact3d
        force(3,ip,icel) = force(3,ip,icel) &
             - lambel * f_3 * charge(ip,icel) &
             * r_scale(3) * vol_fact4 * bs_fact3d

      end do
    end do

    eelec(id+1) = eelec(id+1) + lambel*eelec_temp(id+1)

!#else
!    !
!    ! Gradient on GPU 
!    !
!    !$omp master
!    call gpu_pme_recip_interpolate_force( &
!         force, qdf_real, charge, natom, &
!         MaxAtom, ncell, nlocalx, nlocaly, nlocalz, ngrid, &
!         x_start, x_end, y_start, y_end, z_start, z_end, &
!         r_scale, bs_fact3d, vol_fact4, &
!         vi, bsc, bscd )
!    !$omp end master
!
!    !$omp barrier
!
!#endif

    !$omp barrier

    deallocate(work1, work2)

    !$omp end parallel
 
    return

  end subroutine pme_recip_fep

end module sp_energy_pme_mod
