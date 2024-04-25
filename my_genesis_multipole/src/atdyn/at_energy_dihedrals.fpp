!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_energy_dihedrals_mod
!> @brief   calculate dihedral energy
!! @authors Chigusa Kobayashi(CK), Takao Yoda (TY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_energy_dihedrals_mod

  use messages_mod
  use at_enefunc_str_mod
  use dihedral_libs_mod
  use timers_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! subroutine
  public :: compute_energy_dihed
  public :: compute_energy_rb_dihed
  public :: compute_energy_improp
  public :: compute_energy_improp_cos
  public :: compute_energy_cmap

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_dihed
  !> @brief        calculate dihedral energy
  !! @authors      CK
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[inout] edihe   : dihedral energy of target systems
  !! @note         Blondel and Karplus, J. Comput. Chem., 17, 1132-1141 (1996)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_dihed(enefunc, coord, force, virial, edihe)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: edihe

    ! local variables
    integer                  :: i, j, k, id, krot
    integer                  :: i1, i2, i3, i4
    integer                  :: istart, iend
    integer                  :: aindex(1:4)
    real(wp)                 :: cospha, sinpha, grad_coef, tmp
    real(wp)                 :: cos_dih, sin_dih
    real(wp)                 :: cosnt, sinnt, vtmp
    real(wp)                 :: grad(1:9), v(1:3,1:3)

    real(wp),        pointer :: fc(:), phase(:)
    real(wp),        pointer :: work(:,:)
    integer,         pointer :: list(:,:)
    integer,         pointer :: nperiod(:)


    call timer(TimerDihedral, TimerOn)

    ! use pointers
    !
    istart  = enefunc%istart_dihedral
    iend    = enefunc%iend_dihedral
    list    => enefunc%dihe_list
    fc      => enefunc%dihe_force_const
    nperiod => enefunc%dihe_periodicity
    phase   => enefunc%dihe_phase
    work    => enefunc%work
    edihe = 0.0_wp

    ! calculation of dihedral energy and gradient
    !
    !$omp parallel do default(none)                                            &
    !$omp private(i, j, k, id, aindex, cos_dih, sin_dih, grad, v, grad_coef,   &
    !$omp         cospha, sinpha, cosnt, sinnt, krot, tmp, vtmp)               &
    !$omp shared(istart, iend, list, fc, nperiod, phase, work, coord)          &
    !$omp reduction(+:edihe) reduction(+:virial)
    !
    do i = istart, iend

      aindex(1:4) = list(1:4,i)
      call calculate_dihedral(aindex, coord, cos_dih, sin_dih, grad, v)

      cosnt = 1.0_wp
      sinnt = 0.0_wp
      krot = 0
      do while (krot < nperiod(i))
        tmp   = cosnt*cos_dih - sinnt*sin_dih
        sinnt = sinnt*cos_dih + cosnt*sin_dih
        cosnt = tmp
        krot = krot+1
      end do

      cospha = cos(phase(i))
      sinpha = sin(phase(i))

      edihe = edihe + fc(i) * (1.0_wp + cospha*cosnt + sinnt*sinpha)

      grad_coef = fc(i) * real(nperiod(i),wp) * (cospha*sinnt - cosnt*sinpha)
      work(1:9,i) = grad_coef*grad(1:9)

      ! virial
      !
      do j = 1, 3
        do k = j+1, 3
          vtmp = grad_coef*v(k,j)
          virial(k,j) = virial(k,j) + vtmp
          virial(j,k) = virial(j,k) + vtmp
        end do
        vtmp = grad_coef*v(j,j)
        virial(j,j) = virial(j,j) + vtmp
      end do

    end do
    !$omp end parallel do

    ! store force
    !
    do i = istart, iend
      force(1:3,list(1,i)) = force(1:3,list(1,i)) - work(1:3,i)
      force(1:3,list(2,i)) = force(1:3,list(2,i)) + work(1:3,i) - work(4:6,i)
      force(1:3,list(3,i)) = force(1:3,list(3,i)) + work(4:6,i) + work(7:9,i)
      force(1:3,list(4,i)) = force(1:3,list(4,i)) - work(7:9,i)
    end do

    call timer(TimerDihedral, TimerOff)

    return

  end subroutine compute_energy_dihed

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_rb_dihed
  !> @brief        calculate Ryckaert-Bellemans dihedral energy
  !! @authors      CK
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[inout] edihe   : dihedral energy of target systems
  !! @note         Blondel and Karplus, J. Comput. Chem., 17, 1132-1141 (1996)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_rb_dihed(enefunc, coord, force, virial, edihe)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: edihe

    ! local variables
    integer                  :: i, j, k, istart, iend, icn
    integer                  :: aindex(1:4)
    real(wp)                 :: cospha, sinpha, grad_coef, coef
    real(wp)                 :: cos_dih, sin_dih, vtmp
    real(wp)                 :: grad(1:9), v(1:3,1:3)

    real(wp),        pointer :: fc(:,:)
    real(wp),        pointer :: work(:,:)
    integer,         pointer :: list(:,:)


    call timer(TimerDihedral, TimerOn)

    list    => enefunc%rb_dihe_list
    fc      => enefunc%rb_dihe_c
    work    => enefunc%work

    istart  = enefunc%istart_rb_dihed
    iend    = enefunc%iend_rb_dihed

    ! calculation of dihedral energy and gradient
    !
    !$omp parallel do default(none)                                            &
    !$omp private(i, j, k,  aindex, cos_dih, sin_dih, grad, v, grad_coef,   &
    !$omp         icn,  vtmp, coef)               &
    !$omp shared(istart, iend, list, fc, work, coord)          &
    !$omp reduction(+:edihe) reduction(+:virial)
    !
    do i = istart, iend

      aindex(1:4) = list(1:4,i)
      call calculate_dihedral(aindex, coord, cos_dih, sin_dih, grad, v)

!     psi = phi - pi
      cos_dih = -cos_dih
      sin_dih = -sin_dih

      coef = 0.0_wp
      do icn = 1, 6
        coef = coef + real(icn-1,wp) * fc(icn,i) * cos_dih**(icn-2)
        edihe = edihe + fc(icn,i) * cos_dih**(icn-1)
      end do

      grad_coef = sin_dih * coef
      work(1:9,i) = grad_coef*grad(1:9)

      ! virial 
      !
      do j = 1, 3
        do k = j + 1, 3
          vtmp = grad_coef*v(k,j)
          virial(k,j) = virial(k,j) - vtmp
          virial(j,k) = virial(j,k) - vtmp
        end do
        vtmp = grad_coef*v(j,j)
        virial(j,j) = virial(j,j) - vtmp
      end do

    end do
    !$omp end parallel do

    ! store force
    !
    do i = istart, iend
      force(1:3,list(1,i)) = force(1:3,list(1,i)) - work(1:3,i)
      force(1:3,list(2,i)) = force(1:3,list(2,i)) + work(1:3,i) - work(4:6,i)
      force(1:3,list(3,i)) = force(1:3,list(3,i)) + work(4:6,i) + work(7:9,i)
      force(1:3,list(4,i)) = force(1:3,list(4,i)) - work(7:9,i)
    end do

    call timer(TimerDihedral, TimerOff)

    return

  end subroutine compute_energy_rb_dihed

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_improp
  !> @brief        calculate improper dihedral energy
  !! @authors      CK
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[inout] eimpr   : improper dihedral energy of target systems
  !! @note         Blondel and Karplus, J. Comput. Chem., 17, 1132-1141 (1996)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_improp(enefunc, coord, force, virial, eimpr)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: eimpr

    ! local variables
    integer                  :: i, j, k
    integer                  :: i1, i2, i3, i4
    integer                  :: istart, iend
    integer                  :: aindex(1:4)
    real(wp)                 :: cospha, sinpha, grad_coef, diffphi, vtmp
    real(wp)                 :: cos_dih, sin_dih
    real(wp)                 :: cosdif, sindif
    real(wp)                 :: grad(1:9), v(1:3,1:3)

    real(wp),        pointer :: fc(:), phase(:)
    real(wp),        pointer :: work(:,:)
    integer,         pointer :: list(:,:)
    integer,         pointer :: nperiod(:)


    call timer(TimerDihedral, TimerOn)

    istart = enefunc%istart_improper
    iend   = enefunc%iend_improper

    list   => enefunc%impr_list
    fc     => enefunc%impr_force_const
    phase  => enefunc%impr_phase
    work   => enefunc%work

    ! calculation of dihedral energy and gradient
    !
    !$omp parallel do default(none)                                            &
    !$omp private(i, j, k,  aindex, cos_dih, sin_dih, grad, v, grad_coef,   &
    !$omp         cospha, sinpha, diffphi, cosdif, sindif, vtmp)               &
    !$omp shared(istart, iend, list, fc, nperiod, phase, work, coord)          &
    !$omp reduction(+:eimpr) reduction(+:virial)
    !
    do i = istart, iend

      aindex(1:4) = list(1:4,i)
      call calculate_dihedral(aindex, coord, cos_dih, sin_dih, grad, v)

      cospha = cos(phase(i))
      sinpha = sin(phase(i))

      cosdif = cos_dih*cospha + sin_dih*sinpha
      sindif = cos_dih*sinpha - sin_dih*cospha

      if (cosdif > 1.0E-1_wp) then
        diffphi = asin(sindif)
      else
        diffphi = sign(1.0_wp,sindif)*acos(cosdif)
      endif

      eimpr     = eimpr + fc(i)*diffphi*diffphi

      grad_coef = 2.0_wp*fc(i)*diffphi
      work(1:9,i) = grad_coef*grad(1:9)

      ! virial
      !
      do j = 1, 3
        do k = j+1, 3
          vtmp = grad_coef*v(k,j)
          virial(k,j) = virial(k,j) + vtmp
          virial(j,k) = virial(j,k) + vtmp
        end do
        vtmp = grad_coef*v(j,j)
        virial(j,j) = virial(j,j) + vtmp
      end do
    end do
    !$omp end parallel do

    ! store force
    !
    do i = istart, iend
      i1 = list(1,i)
      i2 = list(2,i)
      i3 = list(3,i)
      i4 = list(4,i)
      force(1,i1) = force(1,i1) - work(1,i)
      force(2,i1) = force(2,i1) - work(2,i)
      force(3,i1) = force(3,i1) - work(3,i)
      force(1,i2) = force(1,i2) + work(1,i) - work(4,i)
      force(2,i2) = force(2,i2) + work(2,i) - work(5,i)
      force(3,i2) = force(3,i2) + work(3,i) - work(6,i)
      force(1,i3) = force(1,i3) + work(4,i) + work(7,i)
      force(2,i3) = force(2,i3) + work(5,i) + work(8,i)
      force(3,i3) = force(3,i3) + work(6,i) + work(9,i)
      force(1,i4) = force(1,i4) - work(7,i)
      force(2,i4) = force(2,i4) - work(8,i)
      force(3,i4) = force(3,i4) - work(9,i)
    end do

    call timer(TimerDihedral, TimerOff)

    return

  end subroutine compute_energy_improp

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_improp_cos
  !> @brief        calculate improper dihedral energy
  !! @authors      CK
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[inout] eimpr   : improper dihedral energy of target systems
  !! @note         Blondel and Karplus, J. Comput. Chem., 17, 1132-1141 (1996)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_improp_cos(enefunc, coord, force, virial, eimpr)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: eimpr

    ! local variables
    integer                  :: i, j, k,  krot
    integer                  :: i1, i2, i3, i4
    integer                  :: istart, iend
    integer                  :: aindex(1:4)
    real(wp)                 :: cospha, sinpha, grad_coef, diffphi, tmp
    real(wp)                 :: cos_dih, sin_dih
    real(wp)                 :: cosdif, sindif
    real(wp)                 :: cosnt, sinnt, vtmp
    real(wp)                 :: grad(1:9), v(1:3,1:3)

    real(wp),        pointer :: fc(:), phase(:)
    real(wp),        pointer :: work(:,:)
    integer,         pointer :: list(:,:)
    integer,         pointer :: nperiod(:)


    call timer(TimerDihedral, TimerOn)

    istart = enefunc%istart_improper
    iend   = enefunc%iend_improper

    list    => enefunc%impr_list
    fc      => enefunc%impr_force_const
    nperiod => enefunc%impr_periodicity
    phase   => enefunc%impr_phase
    work    => enefunc%work

    ! calculation of dihedral energy and gradient
    !
    !$omp parallel do default(none)                                            &
    !$omp private(i, j, k, aindex, cos_dih, sin_dih, grad, v, grad_coef,   &
    !$omp         cospha, sinpha, cosnt, sinnt, krot, tmp, vtmp)               &
    !$omp shared(istart, iend, list, fc, nperiod, phase, work, coord)          &
    !$omp reduction(+:eimpr) reduction(+:virial)
    !
    do i = istart, iend

      aindex(1:4) = list(1:4,i)
      call calculate_dihedral(aindex, coord, cos_dih, sin_dih, grad, v)

      cospha = cos(phase(i))
      sinpha = sin(phase(i))

      cosnt = 1.0_wp
      sinnt = 0.0_wp
      krot = 0
      do while (krot < nperiod(i))
        tmp   = cosnt*cos_dih - sinnt*sin_dih
        sinnt = sinnt*cos_dih + cosnt*sin_dih
        cosnt = tmp
        krot = krot+1
      end do

      cospha = cos(phase(i))
      sinpha = sin(phase(i))

      eimpr = eimpr + fc(i) * (1.0_wp + cospha*cosnt + sinnt*sinpha)

      grad_coef = fc(i) * real(nperiod(i),wp) * (cospha*sinnt - cosnt*sinpha)
      work(1:9,i) = grad_coef*grad(1:9)

      ! virial
      !
      do j = 1, 3
        do k = j+1, 3
          vtmp = grad_coef*v(k,j)
          virial(k,j) = virial(k,j) + vtmp
          virial(j,k) = virial(j,k) + vtmp
        end do
        vtmp = grad_coef*v(j,j)
        virial(j,j) = virial(j,j) + vtmp
      end do
    end do
    !$omp end parallel do

    ! store force
    !
    do i = istart, iend
      i1 = list(1,i)
      i2 = list(2,i)
      i3 = list(3,i)
      i4 = list(4,i)
      force(1,i1) = force(1,i1) - work(1,i)
      force(2,i1) = force(2,i1) - work(2,i)
      force(3,i1) = force(3,i1) - work(3,i)
      force(1,i2) = force(1,i2) + work(1,i) - work(4,i)
      force(2,i2) = force(2,i2) + work(2,i) - work(5,i)
      force(3,i2) = force(3,i2) + work(3,i) - work(6,i)
      force(1,i3) = force(1,i3) + work(4,i) + work(7,i)
      force(2,i3) = force(2,i3) + work(5,i) + work(8,i)
      force(3,i3) = force(3,i3) + work(6,i) + work(9,i)
      force(1,i4) = force(1,i4) - work(7,i)
      force(2,i4) = force(2,i4) - work(8,i)
      force(3,i4) = force(3,i4) - work(9,i)
    end do

    call timer(TimerDihedral, TimerOff)

    return

  end subroutine compute_energy_improp_cos

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_cmap
  !> @brief        calculate cmap energy
  !! @authors      CK, TY
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[inout] ecmap   : cmap energy of target systems
  !! @note         A.D.MacKerell et al., J.Comput.Chem., 25, 1400-1415 (2004).
  !! @note         Blondel and Karplus, J. Comput. Chem., 17, 1132-1141 (1996)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_cmap(enefunc, coord, force, virial, ecmap)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: ecmap

    ! local variables
    real(wp)                 :: cos_dih, sin_dih, dihed
    real(wp)                 :: delta, inv_delta, ctmp, vtmp
    integer                  :: aindex(1:8)
    integer                  :: ncmap, istart, iend
    integer                  :: icmap
    integer                  :: igrid(1:2), ngrid0, i, j, m, igr, ip, ictype
    real(wp)                 :: work(1:9,1:2), grad_coef(1:2), gradtot(1:2)
    real(wp)                 :: dgrid(1:2)
    real(wp)                 :: gradphi(1:9), vphi(1:3,1:3)
    real(wp)                 :: gradpsi(1:9), vpsi(1:3,1:3)
    real(wp)                 :: grid_power(1:2,1:4), dgrid_power(1:2,1:4)

    real(wp),        pointer :: coef(:,:,:,:,:)
    real(wp),        pointer :: f(:,:,:)
    integer,         pointer :: resol(:), list(:,:), ctype(:)


    call timer(TimerDihedral, TimerOn)

    ! use pointers
    !
    istart = enefunc%istart_cmap
    iend   = enefunc%iend_cmap
    ncmap  =  enefunc%num_cmaps
    list   => enefunc%cmap_list
    coef   => enefunc%cmap_coef
    resol  => enefunc%cmap_resolution
    ctype  => enefunc%cmap_type
    f      => enefunc%cmap_force

    ! calculate cmap energy and gradient
    !
    !$omp parallel do default(none)                                       &
    !$omp private(icmap, i, j, aindex, ngrid0, delta, inv_delta, ip,      &
    !$omp         cos_dih, sin_dih, dihed, igr, dgrid, igrid, vtmp,       &
    !$omp         grid_power, dgrid_power, work, m,  ctmp, grad_coef,     &
    !$omp         gradpsi, gradphi, vpsi, vphi, gradtot, ictype)          &
    !$omp shared(istart, iend, list, coord, coef, resol,  f, ctype)       &
    !$omp reduction(+:ecmap) reduction(+:virial)
    do icmap = istart, iend

      ictype    = ctype(icmap)
      ngrid0    = resol(ictype)
      delta     = 360.0_wp/real(ngrid0,wp)
      inv_delta = 1.0_wp/delta

      ! calc dih1 and gradient
      !
      aindex(1:8) = list(1:8,icmap)

      ! dih1 is backbone phi
      call calculate_dihedral(aindex(1:4), coord, cos_dih, sin_dih,          &
                              gradphi, vphi)

      gradphi(1:9) = gradphi(1:9)/RAD

      if (abs(cos_dih) > 1.0E-1_wp) then
        dihed = asin(sin_dih)/RAD
        if (cos_dih < 0.0_wp) then
          if (dihed > 0.0_wp) then
            dihed=180.0_wp-dihed
          else
            dihed=-180.0_wp-dihed
          endif
        endif
      else
        dihed = sign(1.0_wp,sin_dih)*acos(cos_dih)/RAD
      endif

      if (dihed < -180.0_wp) then
        dihed = dihed + 360.0_wp
      else if (dihed > 180.0_wp) then
        dihed = dihed - 360.0_wp
      endif

      igr      = int((dihed+180.0_wp)*inv_delta)
      dgrid(1) = (dihed - (delta*real(igr,wp) - 180.0_wp))*inv_delta
      igrid(1) = igr + 1

      !dbg
      if (igrid(1) < 1 .or. igrid(1) > 24) then
        write(6,'("debug-charmm-6-2")')
        write(6,'(2e20.4)') dihed, inv_delta
        write(6,'(i8)') igrid(1)
      end if
      !dbg

      ! dih2 is backbone psi
      call calculate_dihedral(aindex(5:8), coord, cos_dih, sin_dih,          &
                              gradpsi, vpsi)

      gradpsi(1:9) = gradpsi(1:9)/RAD

      if (abs(cos_dih) > 1.0E-1_wp) then
        dihed = asin(sin_dih)/RAD
        if (cos_dih < 0.0_wp) then
          if (dihed > 0.0_wp) then
            dihed=180.0_wp-dihed
          else
            dihed=-180.0_wp-dihed
          endif
        endif
      else
        dihed = sign(1.0_wp,sin_dih)*acos(cos_dih)/RAD
      endif

      if (dihed < -180.0_wp) then
        dihed = dihed + 360.0_wp
      else if (dihed > 180.0_wp) then
        dihed = dihed - 360.0_wp
      endif

      igr      = int((dihed+180.0_wp)*inv_delta)
      dgrid(2) = (dihed - (delta*real(igr,wp) - 180.0_wp))*inv_delta
      igrid(2) = igr + 1

      !dbg
      if (igrid(2) < 1 .or. igrid(2) > 24) then
        write(6,'("debug-charmm-6-3")')
        write(6,'(2e20.4)') dihed, inv_delta
        write(6,'(i8)') igrid(2)
      end if
      !dbg

      grid_power(1:2,1) = 1.0_wp

      ip = 1
      do while(ip < 4)
        ip = ip + 1
        grid_power(1:2,ip) = grid_power(1:2,ip-1)*dgrid(1:2)
      end do

      dgrid_power(1:2,1) = 0.0_wp
      dgrid_power(1:2,2) = inv_delta

      ip = 2
      do while(ip < 4)
        ip = ip + 1
        dgrid_power(1:2,ip) = grid_power(1:2,ip-1)*real(ip-1,wp)*inv_delta
      end do

      ! calculate Ecmap and gradient
      !
      work(1:9,1:2) = 0.0_wp
      gradtot(1:2) = 0.0_wp

      do j = 1, 4
        do i = 1, 4
          ctmp = coef(i,j,igrid(2),igrid(1),ictype)

          ! cmap energy
          !
          ecmap = ecmap + grid_power(2,i)*grid_power(1,j)*ctmp

          ! for gradient
          !
          grad_coef(1) = -dgrid_power(1,j)*grid_power(2,i)*ctmp
          grad_coef(2) = -dgrid_power(2,i)*grid_power(1,j)*ctmp

          gradtot(1:2) = gradtot(1:2) + grad_coef(1:2)
        end do
      end do
      work(1:9,1) = gradtot(1)*gradphi(1:9)
      work(1:9,2) = gradtot(2)*gradpsi(1:9)

      f(1:3,1,icmap) = -work(1:3,1)
      f(1:3,2,icmap) =  work(1:3,1) - work(4:6,1)
      f(1:3,3,icmap) =  work(4:6,1) + work(7:9,1)
      f(1:3,4,icmap) = -work(7:9,1)

      f(1:3,5,icmap) = -work(1:3,2)
      f(1:3,6,icmap) =  work(1:3,2) - work(4:6,2)
      f(1:3,7,icmap) =  work(4:6,2) + work(7:9,2)
      f(1:3,8,icmap) = -work(7:9,2)

      do i  = 1, 3
        do j  = i+1, 3
          vtmp =  (gradtot(1)*vphi(j,i) + gradtot(2)*vpsi(j,i))/RAD
          virial(j,i) =  virial(j,i)  + vtmp
          virial(i,j) =  virial(i,j)  + vtmp
        end do
        vtmp =  (gradtot(1)*vphi(i,i) + gradtot(2)*vpsi(i,i))/RAD
        virial(i,i) =  virial(i,i)  + vtmp
      end do
    end do
    !$omp end parallel do

    ! store force
    !
    do icmap = istart, iend
      aindex(1:8) = list(1:8, icmap)

      force(1:3, aindex(1)) = force(1:3, aindex(1)) + f(1:3,1,icmap)
      force(1:3, aindex(2)) = force(1:3, aindex(2)) + f(1:3,2,icmap)
      force(1:3, aindex(3)) = force(1:3, aindex(3)) + f(1:3,3,icmap)
      force(1:3, aindex(4)) = force(1:3, aindex(4)) + f(1:3,4,icmap)
      force(1:3, aindex(5)) = force(1:3, aindex(5)) + f(1:3,5,icmap)
      force(1:3, aindex(6)) = force(1:3, aindex(6)) + f(1:3,6,icmap)
      force(1:3, aindex(7)) = force(1:3, aindex(7)) + f(1:3,7,icmap)
      force(1:3, aindex(8)) = force(1:3, aindex(8)) + f(1:3,8,icmap)
    end do

    call timer(TimerDihedral, TimerOff)

    return

  end subroutine compute_energy_cmap

end module at_energy_dihedrals_mod
