!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_enefunc_gamd_mod
!> @brief   define potential energy functions for GaMD
!           Refs: Miao Y. et al., J. Chem. Theory Comput. 11, 3584 (2015)
!                 Pang Y. et al., J. Chem. Theory Comput. 13, 9 (2017)
!! @authors Hiraku Oshima (HO)
!  
!  (c) Copyright 2019 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_enefunc_gamd_mod
  use at_enefunc_str_mod
  use molecules_str_mod
  use fileio_control_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! public subroutines
  public  :: init_enefunc_gamd
  public  :: alloc_enefunc_gamd
  public  :: setup_enefunc_gamd

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_enefunc_gamd
  !> @brief        initialize parameters of GaMD
  !! @authors      HO
  !! @param[inout] gamd    : GaMD information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_enefunc_gamd(gamd)

    ! formal arguments
    type(s_enefunc_gamd), intent(inout) :: gamd

    gamd%gamd_stat               = .false.
    gamd%gamd_boost              = .true.
    gamd%boost_pot               = .false.
    gamd%boost_dih               = .false.
    gamd%boost_dual              = .false.
    gamd%boost_qm                = .false. !added by ykl
    gamd%boost_qm_only           = .false.
    gamd%thresh_lower            = .true.
    gamd%thresh_upper            = .true.
    gamd%ene_pot_max             = -1.0d99
    gamd%ene_pot_min             =  1.0d99
    gamd%ene_pot_ave             =  0.0_wp
    gamd%ene_pot_ave2            =  0.0_wp
    gamd%ene_pot_dev             =  0.0_wp
    gamd%ene_dih_max             = -1.0d99
    gamd%ene_dih_min             =  1.0d99
    gamd%ene_dih_ave             =  0.0_wp
    gamd%ene_dih_ave2            =  0.0_wp
    gamd%ene_dih_dev             =  0.0_wp
    gamd%ene_qm_max              = -1.0d99
    gamd%ene_qm_min              =  1.0d99
    gamd%ene_qm_ave              =  0.0_wp
    gamd%ene_qm_ave2             =  0.0_wp
    gamd%ene_qm_dev              =  0.0_wp
    gamd%count_pot               =  1
    gamd%count_dih               =  1
    gamd%count_qm                =  1
    gamd%ene_pot_th              = 0.0_wp
    gamd%ene_dih_th              = 0.0_wp
    gamd%ene_qm_th               = 0.0_wp

    return

  end subroutine init_enefunc_gamd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_enefunc_gamd
  !> @brief        allocate additional array for GaMD
  !! @authors      HO
  !! @param[inout] gamd    : GaMD information
  !! @param[inout] domain  : domain information
  !! @param[inout] enefunc : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_enefunc_gamd(gamd, molecule, enefunc)

    ! formal arguments
    type(s_enefunc_gamd),     intent(inout) :: gamd
    type(s_molecule), target, intent(inout) :: molecule
    type(s_enefunc),          intent(inout) :: enefunc

    ! local varialbles
    integer :: natoms

    natoms = molecule%num_atoms

    if (gamd%boost_pot) then

      call alloc_enefunc(enefunc, EneFuncGamdRest, natoms)

    else if (gamd%boost_dih) then

      call alloc_enefunc(enefunc, EneFuncGamdDih, natoms)

    else if (gamd%boost_dual) then

      call alloc_enefunc(enefunc, EneFuncGamdRest, natoms)
      call alloc_enefunc(enefunc, EneFuncGamdDih, natoms)

    end if
    
    if (gamd%boost_qm .or. gamd%boost_qm_only) then !added by ykl
      call alloc_enefunc(enefunc, EneFuncGamdQm, natoms)
    end if

  end subroutine alloc_enefunc_gamd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_gamd
  !> @brief        define parameters of GaMD boost potential
  !! @authors      HO, ykl
  !! @param[inout] gamd : GaMD information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_gamd(gamd)

    ! formal arguments
    type(s_enefunc_gamd), intent(inout) :: gamd

    if (gamd%boost_pot) then

      if (gamd%thresh_lower) then
        call setup_enefunc_gamd_lower( &
          gamd%ene_pot_th, gamd%k_pot, gamd%k0_pot, &
          gamd%ene_pot_max, gamd%ene_pot_min, gamd%ene_pot_ave, &
          gamd%ene_pot_dev, gamd%sigma0_pot)
        
        if (gamd%boost_qm) then !added by ykl
          call setup_enefunc_gamd_lower( &
            gamd%ene_qm_th, gamd%k_qm, gamd%k0_qm, &
            gamd%ene_qm_max, gamd%ene_qm_min, gamd%ene_qm_ave, &
            gamd%ene_qm_dev, gamd%sigma0_qm)
        end if
      else if (gamd%thresh_upper) then
        call setup_enefunc_gamd_upper( &
          gamd%ene_pot_th, gamd%k_pot, gamd%k0_pot, &
          gamd%ene_pot_max, gamd%ene_pot_min, gamd%ene_pot_ave, &
          gamd%ene_pot_dev, gamd%sigma0_pot)

        if (gamd%boost_qm) then
          call setup_enefunc_gamd_upper( &
            gamd%ene_qm_th, gamd%k_qm, gamd%k0_qm, &
            gamd%ene_qm_max, gamd%ene_qm_min, gamd%ene_qm_ave, &
            gamd%ene_qm_dev, gamd%sigma0_qm)
        end if
      end if

    else if (gamd%boost_dih) then

      if (gamd%thresh_lower) then
        call setup_enefunc_gamd_lower( &
          gamd%ene_dih_th, gamd%k_dih, gamd%k0_dih, &
          gamd%ene_dih_max, gamd%ene_dih_min, gamd%ene_dih_ave, &
          gamd%ene_dih_dev, gamd%sigma0_dih)
      else if (gamd%thresh_upper) then
        call setup_enefunc_gamd_upper( &
          gamd%ene_dih_th, gamd%k_dih, gamd%k0_dih, &
          gamd%ene_dih_max, gamd%ene_dih_min, gamd%ene_dih_ave, &
          gamd%ene_dih_dev, gamd%sigma0_dih)
      end if

    else if (gamd%boost_dual) then

      if (gamd%thresh_lower) then
        call setup_enefunc_gamd_lower( &
          gamd%ene_pot_th, gamd%k_pot, gamd%k0_pot, &
          gamd%ene_pot_max, gamd%ene_pot_min, gamd%ene_pot_ave, &
          gamd%ene_pot_dev, gamd%sigma0_pot)

        call setup_enefunc_gamd_lower( &
          gamd%ene_dih_th, gamd%k_dih, gamd%k0_dih, &
          gamd%ene_dih_max, gamd%ene_dih_min, gamd%ene_dih_ave, &
          gamd%ene_dih_dev, gamd%sigma0_dih)

        if (gamd%boost_qm) then
          call setup_enefunc_gamd_lower( &
            gamd%ene_qm_th, gamd%k_qm, gamd%k0_qm, &
            gamd%ene_qm_max, gamd%ene_qm_min, gamd%ene_qm_ave, &
            gamd%ene_qm_dev, gamd%sigma0_qm)
        end if
      else if (gamd%thresh_upper) then
        call setup_enefunc_gamd_upper( &
          gamd%ene_pot_th, gamd%k_pot, gamd%k0_pot, &
          gamd%ene_pot_max, gamd%ene_pot_min, gamd%ene_pot_ave, &
          gamd%ene_pot_dev, gamd%sigma0_pot)

        call setup_enefunc_gamd_upper( &
          gamd%ene_dih_th, gamd%k_dih, gamd%k0_dih, &
          gamd%ene_dih_max, gamd%ene_dih_min, gamd%ene_dih_ave, &
          gamd%ene_dih_dev, gamd%sigma0_dih)

        if (gamd%boost_qm) then
          call setup_enefunc_gamd_upper( &
            gamd%ene_qm_th, gamd%k_qm, gamd%k0_qm, &
            gamd%ene_qm_max, gamd%ene_qm_min, gamd%ene_qm_ave, &
            gamd%ene_qm_dev, gamd%sigma0_qm)
        end if
      end if

    else if (gamd%boost_qm_only) then
      if (gamd%thresh_lower) then
          call setup_enefunc_gamd_lower( &
            gamd%ene_qm_th, gamd%k_qm, gamd%k0_qm, &
            gamd%ene_qm_max, gamd%ene_qm_min, gamd%ene_qm_ave, &
            gamd%ene_qm_dev, gamd%sigma0_qm)
      else if (gamd%thresh_upper) then
          call setup_enefunc_gamd_upper( &
            gamd%ene_qm_th, gamd%k_qm, gamd%k0_qm, &
            gamd%ene_qm_max, gamd%ene_qm_min, gamd%ene_qm_ave, &
            gamd%ene_qm_dev, gamd%sigma0_qm)
      end if
    end if

  end subroutine setup_enefunc_gamd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_gamd_lower
  !> @brief        define parameters of lower threshold in GaMD
  !! @authors      HO
  !! @param[inout] E_th   : threshold of GaMD
  !! @param[inout] k      : force constant
  !! @param[inout] k0     : normalized force constant
  !! @param[in]    E_max  : maximum of energy
  !! @param[in]    E_min  : minimum of energy
  !! @param[in]    E_ave  : average of energy
  !! @param[in]    E_dev  : standard deviation of energy
  !! @param[in]    sigma0 : control parameter of GaMD
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_gamd_lower(E_th, k, k0, E_max, E_min, E_ave, &
      E_dev, sigma0)

    ! formal arguments
    real(wp),                 intent(inout) :: E_th
    real(wp),                 intent(inout) :: k
    real(wp),                 intent(inout) :: k0
    real(wp),                    intent(in) :: E_max
    real(wp),                    intent(in) :: E_min
    real(wp),                    intent(in) :: E_ave
    real(wp),                    intent(in) :: E_dev
    real(wp),                    intent(in) :: sigma0

    E_th = E_max
    k0   = (sigma0/E_dev) * (E_max-E_min) / (E_ave-E_min)
    k0   = min(1.0, k0)
    k = k0 / (E_max-E_min)

  end subroutine setup_enefunc_gamd_lower

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_gamd_upper
  !> @brief        define parameters of upper threshold in GaMD
  !! @authors      HO
  !! @param[inout] E_th   : threshold of GaMD
  !! @param[inout] k      : force constant
  !! @param[inout] k0     : normalized force constant
  !! @param[in]    E_max  : maximum of energy
  !! @param[in]    E_min  : minimum of energy
  !! @param[in]    E_ave  : average of energy
  !! @param[in]    E_dev  : standard deviation of energy
  !! @param[in]    sigma0 : control parameter of GaMD
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_gamd_upper(E_th, k, k0, E_max, E_min, E_ave, &
      E_dev, sigma0)

    ! formal arguments
    real(wp),                 intent(inout) :: E_th
    real(wp),                 intent(inout) :: k
    real(wp),                 intent(inout) :: k0
    real(wp),                    intent(in) :: E_max
    real(wp),                    intent(in) :: E_min
    real(wp),                    intent(in) :: E_ave
    real(wp),                    intent(in) :: E_dev
    real(wp),                    intent(in) :: sigma0

    k0   = (1.0_wp - sigma0/E_dev) * (E_max-E_min) / (E_ave-E_min)
    if ((k0 > 0.0_wp) .and. (k0 <= 1.0_wp)) then
      E_th = E_min + (E_max-E_min)/k0
    else
      E_th = E_max
      k0   = (sigma0/E_dev) * (E_max-E_min) / (E_ave-E_min)
      k0   = min(1.0, k0)
    end if
    k = k0 / (E_max-E_min)

  end subroutine setup_enefunc_gamd_upper

end module at_enefunc_gamd_mod
