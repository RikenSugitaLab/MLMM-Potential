!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_fep_utils_mod
!> @brief   Utilities for FEP calculations
!! @authors Nobuhiko Kato (NK), Hiraku Oshima (HO)
!
!  (c) Copyright 2020 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_fep_utils_mod

  use constants_mod
  use sp_domain_str_mod
  use sp_enefunc_str_mod

  implicit none
  private

  ! subroutines
  public  :: sync_single_fep
  public  :: add_single_fep
  public  :: check_fep
  public  :: set_lambda_table_fep

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    sync_single_fep
  !> @brief        synchronize arrays (e.g. coord, vel, and force) of singleA
  !                and singleB in FEP
  !! @authors      HO
  !! @param[in]    domain  : domain information
  !! @param[inout] array   : given array (coord, vel, and force)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine sync_single_fep(domain, array)

    ! formal arguments
    type(s_domain), target, intent(in) :: domain
    real(dp),            intent(inout) :: array(:,:,:)

    ! local variables
    integer                  :: i, ncell
    integer,         pointer :: natom(:)
    integer,         pointer :: natom_singleB(:)
    integer,         pointer :: id_singleA(:,:,:)
    integer,         pointer :: id_singleB(:,:,:)
    integer                  :: iA, iB, iatomA, iatomB, icel

    natom         => domain%num_atom
    ncell         =  domain%num_cell_local + domain%num_cell_boundary
    natom_singleB => domain%num_atom_singleB
    id_singleA    => domain%id_singleA
    id_singleB    => domain%id_singleB

    do icel = 1, ncell
      do i = 1, natom_singleB(icel)
        iA = id_singleA(i, icel, 3)
        iB = id_singleB(i, icel, 3)
        iatomA = id_singleA(iA, icel, 1)
        iatomB = id_singleB(iB, icel, 1)
        array(1:3,iatomB,icel) = array(1:3,iatomA,icel)
      end do
    end do

    return

  end subroutine sync_single_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    add_single_fep
  !> @brief        sum two arrays (e.g. force) of singleA and singleB in FEP
  !! @authors      HO
  !! @param[in]    domain  : domain information
  !! @param[inout] array   : given array (force)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine add_single_fep(domain, array)

    ! formal arguments
    type(s_domain), target, intent(in) :: domain
    real(dp),            intent(inout) :: array(:,:,:)

    ! local variables
    integer                  :: i, ncell
    integer,         pointer :: natom(:)
    integer,         pointer :: natom_singleB(:)
    integer,         pointer :: id_singleA(:,:,:)
    integer,         pointer :: id_singleB(:,:,:)
    integer                  :: iA, iB, iatomA, iatomB, icel
    real(dp)                 :: arrayA_temp(1:3), arrayB_temp(1:3)

    natom         => domain%num_atom
    ncell         =  domain%num_cell_local + domain%num_cell_boundary
    natom_singleB => domain%num_atom_singleB
    id_singleA    => domain%id_singleA
    id_singleB    => domain%id_singleB

    do icel = 1, ncell
      do i = 1, natom_singleB(icel)
        iA = id_singleA(i, icel, 3)
        iB = id_singleB(i, icel, 3)
        iatomA = id_singleA(iA, icel, 1)
        iatomB = id_singleB(iB, icel, 1)
        arrayA_temp(:) = array(:,iatomA,icel)
        arrayB_temp(:) = array(:,iatomB,icel)
        array(:,iatomA,icel) = array(:,iatomA,icel) + arrayB_temp(:)
        array(:,iatomB,icel) = array(:,iatomB,icel) + arrayA_temp(:)
      end do
    end do

    return

  end subroutine add_single_fep

  !======1=========2=========3=========4=========5=========6=========7=========8

!  subroutine add_single_mass_weighted_fep(domain, array)
!
!    ! formal arguments
!    type(s_domain), target, intent(in) :: domain
!    real(dp),            intent(inout) :: array(:,:,:)
!
!    ! local variables
!    integer                  :: i, ncell
!    integer,         pointer :: natom(:)
!    integer,         pointer :: natom_singleB(:)
!    integer,         pointer :: id_singleA(:,:,:)
!    integer,         pointer :: id_singleB(:,:,:)
!    integer                  :: iA, iB, iatomA, iatomB, icel
!    real(dp)                 :: arrayA_temp(1:3), arrayB_temp(1:3)
!    real(dp),        pointer :: mass(:,:)
!    real(dp)                 :: massA, massB, mass_merged
!    real(dp)                 :: lambA, lambB
!
!    natom         => domain%num_atom
!    ncell         =  domain%num_cell_local + domain%num_cell_boundary
!    natom_singleB => domain%num_atom_singleB
!    id_singleA    => domain%id_singleA
!    id_singleB    => domain%id_singleB
!    mass          => domain%mass
!    lambA         =  real(domain%lambbondA,dp)
!    lambB         =  real(domain%lambbondB,dp)
!
!    do icel = 1, ncell
!      do i = 1, natom_singleB(icel)
!        iA = id_singleA(i, icel, 3)
!        iB = id_singleB(i, icel, 3)
!        iatomA = id_singleA(iA, icel, 1)
!        iatomB = id_singleB(iB, icel, 1)
!        arrayA_temp(:) = array(:,iatomA,icel)
!        arrayB_temp(:) = array(:,iatomB,icel)
!        array(:,iatomA,icel) = array(:,iatomA,icel) + arrayB_temp(:)
!        array(:,iatomB,icel) = array(:,iatomB,icel) + arrayA_temp(:)
!
!        massA = mass(iatomA,icel)
!        massB = mass(iatomB,icel)
!        mass_merged = lambA * massA + lambB * massB
!        array(:,iatomA,icel) = array(:,iatomA,icel) * massA / mass_merged
!        array(:,iatomB,icel) = array(:,iatomB,icel) * massB / mass_merged
!      end do
!    end do
!
!    return
!
!  end subroutine add_single_mass_weighted_fep

  !======1=========2=========3=========4=========5=========6=========7=========8

!  subroutine scale_mass_single_fep(domain)
!
!    ! formal arguments
!    type(s_domain), target, intent(in) :: domain
!
!    ! local variables
!    integer                  :: i, ncell
!    integer,         pointer :: natom(:)
!    integer,         pointer :: natom_singleB(:)
!    integer,         pointer :: id_singleA(:,:,:)
!    integer,         pointer :: id_singleB(:,:,:)
!    integer                  :: iA, iB, iatomA, iatomB, icel
!    real(dp),        pointer :: mass(:,:), mass_fep(:,:)
!    real(dp)                 :: massA, massB, mass_merged
!    real(dp)                 :: lambA, lambB
!
!    natom         => domain%num_atom
!    ncell         =  domain%num_cell_local + domain%num_cell_boundary
!    natom_singleB => domain%num_atom_singleB
!    id_singleA    => domain%id_singleA
!    id_singleB    => domain%id_singleB
!    mass          => domain%mass
!    mass_fep      => domain%mass_fep
!    lambA         =  real(domain%lambbondA,dp)
!    lambB         =  real(domain%lambbondB,dp)
!
!    mass_fep(:,:) = mass(:,:)
!    do icel = 1, ncell
!      do i = 1, natom_singleB(icel)
!        iA = id_singleA(i, icel, 3)
!        iB = id_singleB(i, icel, 3)
!        iatomA = id_singleA(iA, icel, 1)
!        iatomB = id_singleB(iB, icel, 1)
!        massA = mass(iatomA,icel)
!        massB = mass(iatomB,icel)
!        mass_merged = lambA * massA + lambB * massB
!        mass_fep(iatomA,icel) = mass_merged
!        mass_fep(iatomB,icel) = mass_merged
!      end do
!    end do
!
!    return
!
!  end subroutine scale_mass_single_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    check_fep
  !> @brief        check coord, vel, and force of singleA and singleB
  !! @authors      HO
  !! @param[in]    domain  : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_fep(domain)

    ! formal arguments
    type(s_domain),  target, intent(inout) :: domain

    ! local variables
    integer                  :: i, ncell
    real(dp),        pointer :: coord(:,:,:)
    real(dp),        pointer :: vel(:,:,:)
    real(dp),        pointer :: force(:,:,:)
    integer,         pointer :: natom(:)
    integer,         pointer :: natom_singleB(:)
    integer,         pointer :: id_singleA(:,:,:)
    integer,         pointer :: id_singleB(:,:,:)
    integer                  :: iA, iB, iatomA, iatomB, icel, idA, idB

    natom         => domain%num_atom
    coord         => domain%coord
    vel           => domain%velocity
    force         => domain%force
    ncell         =  domain%num_cell_local + domain%num_cell_boundary
    natom_singleB => domain%num_atom_singleB
    id_singleA    => domain%id_singleA
    id_singleB    => domain%id_singleB

    do icel = 1, ncell
      do i = 1, natom_singleB(icel)
        iA = id_singleA(i, icel, 3)
        iB = id_singleB(i, icel, 3)
        iatomA = id_singleA(iA, icel, 1)
        iatomB = id_singleB(iB, icel, 1)
        idA = domain%id_l2g(iatomA,icel)
        idB = domain%id_l2g(iatomB,icel)

        write(*,*) "[debug001]_coordX", idA, idB, &
          coord(1,iatomA,icel)-coord(1,iatomB,icel)
        write(*,*) "[debug001]_coordY", idA, idB, &
          coord(2,iatomA,icel)-coord(2,iatomB,icel)
        write(*,*) "[debug001]_coordZ", idA, idB, &
          coord(3,iatomA,icel)-coord(3,iatomB,icel)

        write(*,*) "[debug001]_velX", idA, idB, &
          vel(1,iatomA,icel)-vel(1,iatomB,icel)
        write(*,*) "[debug001]_velY", idA, idB, &
          vel(2,iatomA,icel)-vel(2,iatomB,icel)
        write(*,*) "[debug001]_velZ", idA, idB, &
          vel(3,iatomA,icel)-vel(3,iatomB,icel)

        write(*,*) "[debug001]_forceX", idA, idB, &
          force(1,iatomA,icel)-force(1,iatomB,icel)
        write(*,*) "[debug001]_forceY", idA, idB, &
          force(2,iatomA,icel)-force(2,iatomB,icel)
        write(*,*) "[debug001]_forceZ", idA, idB, &
          force(3,iatomA,icel)-force(3,iatomB,icel)

        if ((coord(1,iatomA,icel)-coord(1,iatomB,icel))/=0.0_dp) stop

      end do
    end do

    return

  end subroutine check_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    set_lambda_table_fep
  !> @brief        make lambda table for calculations of bonded and nonbonded
  !                energies in FEP
  !! @authors      HO
  !! @param[in]    enefunc  : enefunc information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine set_lambda_table_fep(enefunc)

    ! formal arguments
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: fg1, fg2
    integer                  :: i1, i2, i3, i4
    integer                  :: i5, i6, i7, i8
    integer                  :: idx
    integer,         pointer :: fepgrp_nonb(:,:)
    real(wp),        pointer :: table_nonb_lambda(:,:,:)
    real(wp),        pointer :: table_bond_lambda(:,:)
    real(wp),        pointer :: table_angl_lambda(:,:,:)
    real(wp),        pointer :: table_dihe_lambda(:,:,:,:)
    real(wp),        pointer :: table_cmap_lambda(:)

    fepgrp_nonb       => enefunc%fepgrp_nonb
    table_nonb_lambda => enefunc%table_nonb_lambda
    table_bond_lambda => enefunc%table_bond_lambda
    table_angl_lambda => enefunc%table_angl_lambda
    table_dihe_lambda => enefunc%table_dihe_lambda
    table_cmap_lambda => enefunc%table_cmap_lambda

    ! Make lambda table for calculation of nonbonded energy
    !
    fepgrp_nonb(:,:) = 0
    table_nonb_lambda(:,:,:) = 0.0_wp
    do fg1 = 1, 5
      do fg2 = 1, 5

        ! Make table of lambda and softcore in FEP
        !
        if ((fg1==5).and.(fg2==5)) then
          ! preserved-preserved
          fepgrp_nonb(fg1,fg2) = 5
          table_nonb_lambda(1,fg1,fg2) = 1.0_wp
          table_nonb_lambda(2,fg1,fg2) = 1.0_wp
        else if (((fg1==1).and.(fg2==1)) .or. &
          ((fg1==1).and.(fg2==5)) .or. &
          ((fg1==5).and.(fg2==1))) then
          ! singleA-singleA and singleA-preserved
          fepgrp_nonb(fg1,fg2) = 1
          table_nonb_lambda(1,fg1,fg2) = enefunc%lambljA
          table_nonb_lambda(2,fg1,fg2) = enefunc%lambelA
        else if (((fg1==2).and.(fg2==2)) .or. &
          ((fg1==2).and.(fg2==5)) .or. &
          ((fg1==5).and.(fg2==2))) then
          ! singleB-singleB and singleB-preserved
          fepgrp_nonb(fg1,fg2) = 2
          table_nonb_lambda(1,fg1,fg2) = enefunc%lambljB
          table_nonb_lambda(2,fg1,fg2) = enefunc%lambelB
        else if (((fg1==3).and.(fg2==3)) .or. &
          ((fg1==3).and.(fg2==5)) .or. &
          ((fg1==5).and.(fg2==3)) .or. &
          ((fg1==3).and.(fg2==1)) .or. &
          ((fg1==1).and.(fg2==3))) then
          ! dualA-dualA and dualA-other
          fepgrp_nonb(fg1,fg2) = 3
          table_nonb_lambda(1,fg1,fg2) = enefunc%lambljA
          table_nonb_lambda(2,fg1,fg2) = enefunc%lambelA
          table_nonb_lambda(3,fg1,fg2) = enefunc%sc_alpha*(1.0_wp-enefunc%lambljA)
          table_nonb_lambda(4,fg1,fg2) = enefunc%sc_beta*(1.0_wp-enefunc%lambelA)
        else if (((fg1==4).and.(fg2==4)) .or. &
          ((fg1==4).and.(fg2==5)) .or. &
          ((fg1==5).and.(fg2==4)) .or. &
          ((fg1==4).and.(fg2==2)) .or. &
          ((fg1==2).and.(fg2==4))) then
          ! dualB-dualB and dualB-other
          fepgrp_nonb(fg1,fg2) = 4
          table_nonb_lambda(1,fg1,fg2) = enefunc%lambljB
          table_nonb_lambda(2,fg1,fg2) = enefunc%lambelB
          table_nonb_lambda(3,fg1,fg2) = enefunc%sc_alpha*(1.0_wp-enefunc%lambljB)
          table_nonb_lambda(4,fg1,fg2) = enefunc%sc_beta*(1.0_wp-enefunc%lambelB)
        end if

        ! Interactons between partA and partB are set to 0,
        ! while others are set to 1.
        !
        if (((fg1==1).and.(fg2==2)) .or. ((fg1==2).and.(fg2==1)) .or. &
          ((fg1==1).and.(fg2==4)) .or. ((fg1==4).and.(fg2==1)) .or. &
          ((fg1==3).and.(fg2==2)) .or. ((fg1==2).and.(fg2==3)) .or. &
          ((fg1==3).and.(fg2==4)) .or. ((fg1==4).and.(fg2==3))) then
          fepgrp_nonb(fg1,fg2) = 0
          table_nonb_lambda(5,fg1,fg2) = 0.0_wp
        else
          table_nonb_lambda(5,fg1,fg2) = 1.0_wp
        end if
      end do
    end do
#ifdef USE_GPU
    call gpu_upload_table_nonb_lambda( table_nonb_lambda )
#endif /* USE_GPU */


    ! Make lambda table for calculation of bonded energy
    !
    table_bond_lambda(:,:) = 0.0_wp
    table_angl_lambda(:,:,:) = 0.0_wp
    table_dihe_lambda(:,:,:,:) = 0.0_wp
    table_cmap_lambda(:) = 0.0_wp

    ! Bond table
    !
    do i2 = 1, 5
      do i1 = 1, 5
        ! preserve: all 5
        if ((i1 == 5) .and. (i2 == 5)) then
          table_bond_lambda(i1,i2) = 1.0_wp
        end if
        ! singleA: all 1
        if ((i1 == 1) .and. (i2 == 1)) then
          table_bond_lambda(i1,i2) = enefunc%lambbondA
        end if
        ! singleB: all 2
        if ((i1 == 2) .and. (i2 == 2)) then
          table_bond_lambda(i1,i2) = enefunc%lambbondB
        end if
        ! dualA: all 3, 3 and 5, 3 and 1
        if (((i1 == 3) .or. (i1 == 5) .or. (i1 == 1)) .and. &
          ((i2 == 3) .or. (i2 == 5) .or. (i2 == 1))) then
          if ((i1 == 3) .or. (i2 == 3)) then
            table_bond_lambda(i1,i2) = 1.0_wp
          end if
        end if
        ! dualB: all 4, 4 and 5, 4 and 2
        if (((i1 == 4) .or. (i1 == 5) .or. (i1 == 2)) .and. &
          ((i2 == 4) .or. (i2 == 5) .or. (i2 == 2))) then
          if ((i1 == 4) .or. (i2 == 4)) then
            table_bond_lambda(i1,i2) = 1.0_wp
          end if
        end if
      end do
    end do

    ! Angle table
    !
    do i3 = 1, 5
      do i2 = 1, 5
        do i1 = 1, 5
          ! preserve: all 5
          if ((i1 == 5) .and. (i2 == 5) .and. (i3 == 5)) then
            table_angl_lambda(i1,i2,i3) = 1.0_wp
          end if
          ! singleA: all 1
          if ((i1 == 1) .and. (i2 == 1) .and. (i3 == 1)) then
            table_angl_lambda(i1,i2,i3) = enefunc%lambbondA
          end if
          ! singleB: all 2
          if ((i1 == 2) .and. (i2 == 2) .and. (i3 == 2)) then
            table_angl_lambda(i1,i2,i3) = enefunc%lambbondB
          end if
          ! dualA: all 3, 3 and 5, 3 and 1
          if (((i1 == 3) .or. (i1 == 5) .or. (i1 == 1)) .and. &
            ((i2 == 3) .or. (i2 == 5) .or. (i2 == 1)) .and. &
            ((i3 == 3) .or. (i3 == 5) .or. (i3 == 1))) then
            if ((i1 == 3) .or. (i2 == 3) .or. (i3 == 3)) then
              table_angl_lambda(i1,i2,i3) = 1.0_wp
            end if
          end if
          ! dualB: all 4, 4 and 5, 4 and 2
          if (((i1 == 4) .or. (i1 == 5) .or. (i1 == 2)) .and. &
            ((i2 == 4) .or. (i2 == 5) .or. (i2 == 2)) .and. &
            ((i3 == 4) .or. (i3 == 5) .or. (i3 == 2))) then
            if ((i1 == 4) .or. (i2 == 4) .or. (i3 == 4)) then
              table_angl_lambda(i1,i2,i3) = 1.0_wp
            end if
          end if
        end do
      end do
    end do

    ! Dihedral table
    !
    do i4 = 1, 5
      do i3 = 1, 5
        do i2 = 1, 5
          do i1 = 1, 5
            ! preserve: all 5
            if ((i1 == 5) .and. (i2 == 5) .and. (i3 == 5) .and. (i4 == 5)) then
              table_dihe_lambda(i1,i2,i3,i4) = 1.0_wp
            end if
            ! singleA: all 1
            if ((i1 == 1) .and. (i2 == 1) .and. (i3 == 1) .and. (i4 == 1)) then
              table_dihe_lambda(i1,i2,i3,i4) = enefunc%lambbondA
            end if
            ! singleB: all 2
            if ((i1 == 2) .and. (i2 == 2) .and. (i3 == 2) .and. (i4 == 2)) then
              table_dihe_lambda(i1,i2,i3,i4) = enefunc%lambbondB
            end if
            ! dualA: all 3, 3 and 5, 3 and 1
            if (((i1 == 3) .or. (i1 == 5) .or. (i1 == 1)) .and. &
              ((i2 == 3) .or. (i2 == 5) .or. (i2 == 1)) .and. &
              ((i3 == 3) .or. (i3 == 5) .or. (i3 == 1)) .and. &
              ((i4 == 3) .or. (i4 == 5) .or. (i4 == 1))) then
              if ((i1 == 3) .or. (i2 == 3) .or. (i3 == 3) .or. (i4 == 3)) then
                table_dihe_lambda(i1,i2,i3,i4) = 1.0_wp
              end if
            end if
            ! dualB: all 4, 4 and 5, 4 and 2
            if (((i1 == 4) .or. (i1 == 5) .or. (i1 == 2)) .and. &
              ((i2 == 4) .or. (i2 == 5) .or. (i2 == 2)) .and. &
              ((i3 == 4) .or. (i3 == 5) .or. (i3 == 2)) .and. &
              ((i4 == 4) .or. (i4 == 5) .or. (i4 == 2))) then
              if ((i1 == 4) .or. (i2 == 4) .or. (i3 == 4) .or. (i4 == 4)) then
                table_dihe_lambda(i1,i2,i3,i4) = 1.0_wp
              end if
            end if
          end do
        end do
      end do
    end do

    ! CMAP table
    !
    do i8 = 1, 5
      do i7 = 1, 5
        do i6 = 1, 5
          do i5 = 1, 5
            do i4 = 1, 5
              do i3 = 1, 5
                do i2 = 1, 5
                  do i1 = 1, 5
                    idx = i1 + 5*(i2-1 + 5*(i3-1 + 5*(i4-1 + 5*(i5-1 + 5*(i6-1 + 5*(i7-1 + 5*(i8-1)))))))
                    ! preserve: all 5
                    if ((i1 == 5) .and. (i2 == 5) .and. (i3 == 5) .and. (i4 == 5) .and. &
                        (i5 == 5) .and. (i6 == 5) .and. (i7 == 5) .and. (i8 == 5)) then
                      table_cmap_lambda(idx) = 1.0_wp
                    end if
                    ! singleA: all 1
                    if ((i1 == 1) .and. (i2 == 1) .and. (i3 == 1) .and. (i4 == 1) .and. &
                        (i5 == 1) .and. (i6 == 1) .and. (i7 == 1) .and. (i8 == 1)) then
                      table_cmap_lambda(idx) = enefunc%lambbondA
                    end if
                    ! singleB: all 2
                    if ((i1 == 2) .and. (i2 == 2) .and. (i3 == 2) .and. (i4 == 2) .and. &
                        (i5 == 2) .and. (i6 == 2) .and. (i7 == 2) .and. (i8 == 2)) then
                      table_cmap_lambda(idx) = enefunc%lambbondB
                    end if
                    ! dualA: all 3, 3 and 5, 3 and 1
                    if (((i1 == 3) .or. (i1 == 5) .or. (i1 == 1)) .and. &
                        ((i2 == 3) .or. (i2 == 5) .or. (i2 == 1)) .and. &
                        ((i3 == 3) .or. (i3 == 5) .or. (i3 == 1)) .and. &
                        ((i4 == 3) .or. (i4 == 5) .or. (i4 == 1)) .and. &
                        ((i5 == 3) .or. (i5 == 5) .or. (i5 == 1)) .and. &
                        ((i6 == 3) .or. (i6 == 5) .or. (i6 == 1)) .and. &
                        ((i7 == 3) .or. (i7 == 5) .or. (i7 == 1)) .and. &
                        ((i8 == 3) .or. (i8 == 5) .or. (i8 == 1))) then
                      if ((i1 == 3) .or. (i2 == 3) .or. (i3 == 3) .or. (i4 == 3) .or. &
                          (i5 == 3) .or. (i6 == 3) .or. (i7 == 3) .or. (i8 == 3)) then
                        table_cmap_lambda(idx) = 1.0_wp
                      end if
                    end if
                    ! dualB: all 4, 4 and 5, 4 and 2
                    if (((i1 == 4) .or. (i1 == 5) .or. (i1 == 2)) .and. &
                        ((i2 == 4) .or. (i2 == 5) .or. (i2 == 2)) .and. &
                        ((i3 == 4) .or. (i3 == 5) .or. (i3 == 2)) .and. &
                        ((i4 == 4) .or. (i4 == 5) .or. (i4 == 2)) .and. &
                        ((i5 == 4) .or. (i5 == 5) .or. (i5 == 2)) .and. &
                        ((i6 == 4) .or. (i6 == 5) .or. (i6 == 2)) .and. &
                        ((i7 == 4) .or. (i7 == 5) .or. (i7 == 2)) .and. &
                        ((i8 == 4) .or. (i8 == 5) .or. (i8 == 2))) then
                      if ((i1 == 4) .or. (i2 == 4) .or. (i3 == 4) .or. (i4 == 4) .or. &
                          (i5 == 4) .or. (i6 == 4) .or. (i7 == 4) .or. (i8 == 4)) then
                        table_cmap_lambda(idx) = 1.0_wp
                      end if
                    end if
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do

    return

  end subroutine set_lambda_table_fep

end module sp_fep_utils_mod
