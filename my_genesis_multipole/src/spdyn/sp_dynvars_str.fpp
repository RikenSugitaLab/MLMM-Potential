!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_dynvars_str_mod
!> @brief   structure of dynvars
!! @authors Yuji Sugita (YS), Takaharu Mori (TM), Chigusa Kobayashi (CK)
!! @note    dynvars include coord, force, energy, boxsize, pressure.
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_dynvars_str_mod

  use sp_energy_str_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_dynvars
    type(s_energy)                :: energy

    real(dp)                      :: virial(3,3)
    real(dp)                      :: virial_long(3,3)
    real(dp)                      :: virial_const(3,3)
    real(dp)                      :: virial_extern(3,3)

    integer                       :: step
    real(dp)                      :: time
    real(dp)                      :: total_pene
    real(dp)                      :: total_kene
    real(dp)                      :: total_energy
    real(dp)                      :: temperature
    real(dp)                      :: rms_gradient
    real(dp)                      :: hfc_kene
    real(dp)                      :: virial_kene

    real(dp)                      :: volume
    real(dp)                      :: internal_virial
    real(dp)                      :: internal_pressure
    real(dp)                      :: external_virial
    real(dp)                      :: external_pressure

    real(dp)                      :: thermostat_momentum
    real(dp)                      :: thermostat_momentum_ref
    real(dp)                      :: barostat_momentum(3)
    real(dp)                      :: barostat_momentum_ref(3)

    real(dp)                      :: nh_mass(1:10)
    real(dp)                      :: nh_velocity(1:10)
    real(dp)                      :: nh_force(1:10)
    real(dp)                      :: nh_coef(1:10)
    real(dp)                      :: nh_baro_velocity(1:10)
    real(dp)                      :: nh_baro_force(1:10)
    real(dp)                      :: nh_baro_coef(1:10)

    real(dp)                      :: box_size_x
    real(dp)                      :: box_size_y
    real(dp)                      :: box_size_z
    real(dp)                      :: kinetic_ref(3)
    real(dp)                      :: pressure_xx
    real(dp)                      :: pressure_yy
    real(dp)                      :: pressure_zz
  end type s_dynvars

  ! subroutines
  public :: init_dynvars

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_dynvars
  !> @brief        initialize energy terms in dynvars
  !! @authors      YS, TM, CK
  !! @param[out]   dynvars : dynamic variable information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_dynvars(dynvars)

    ! formal arguments
    type(s_dynvars),         intent(inout) :: dynvars


    call init_energy (dynvars%energy)

    dynvars%virial(1:3,1:3)            = 0.0_dp
    dynvars%virial_const(1:3,1:3)      = 0.0_dp
    dynvars%virial_extern(1:3,1:3)     = 0.0_dp

    dynvars%step                       = 0
    dynvars%time                       = 0.0_dp
    dynvars%total_pene                 = 0.0_dp
    dynvars%total_kene                 = 0.0_dp
    dynvars%total_energy               = 0.0_dp
    dynvars%temperature                = 0.0_dp
    dynvars%rms_gradient               = 0.0_dp
    dynvars%hfc_kene                   = 0.0_dp
    dynvars%virial_kene                = 0.0_dp

    dynvars%volume                     = 0.0_dp
    dynvars%internal_virial            = 0.0_dp
    dynvars%internal_pressure          = 0.0_dp
    dynvars%external_virial            = 0.0_dp
    dynvars%external_pressure          = 0.0_dp

    dynvars%thermostat_momentum        = 0.0_dp
    dynvars%thermostat_momentum_ref    = 0.0_dp
    dynvars%barostat_momentum(1:3)     = 0.0_dp
    dynvars%barostat_momentum_ref(1:3) = 0.0_dp

    dynvars%box_size_x                 = 0.0_dp
    dynvars%box_size_y                 = 0.0_dp
    dynvars%box_size_z                 = 0.0_dp
    dynvars%pressure_xx                = 0.0_dp
    dynvars%pressure_yy                = 0.0_dp
    dynvars%pressure_zz                = 0.0_dp

    return

  end subroutine init_dynvars

end module sp_dynvars_str_mod
