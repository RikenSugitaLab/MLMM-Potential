!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   table_libs_mod
!> @brief   libralies for tables
!! @authors Chigusa Kobayashi (CK), Jaewoon Jung (JJ)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module table_libs_mod
  use math_libs_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  public :: table_pme_linear_noswitch
  public :: table_pme_linear_switch
  public :: table_pme_linear_fswitch
  public :: table_pme_linear_groshift
  public :: table_pme_linear_groswitch
  public :: table_water_pme_linear_noswitch
  public :: table_water_pme_linear_switch
  public :: table_water_pme_linear_fswitch
  public :: table_water_pme_linear_groshift
  public :: table_water_pme_linear_groswitch
  public :: table_cutoff_cubic_noswitch
  public :: table_cutoff_cubic_switch
  public :: table_cutoff_cubic_fswitch
  public :: table_cutoff_cubic_groshift
  public :: table_cutoff_cubic_groswitch
  public :: table_cutoff_cubic_grodoubleshift

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    table_pme_linear_noswitch
  !> @brief        table for pme and linear  without switch
  !! @authors      CK, JJ
  !! @param[in]    switch_int  : integer switchdist
  !! @param[in]    cutoff_int  : integer cutoff distance
  !! @param[in]    density     : table density
  !! @param[in]    cutoff2     : square cutoff distance
  !! @param[in]    switchdist2 : square switch distance
  !! @param[in]    el_fact     : inverse of dielectoric constant
  !! @param[in]    alpha       : ewald factor
  !! @param[in]    alpha2sp    : ewald factor
  !! @param[in]    alpha2m     : ewald factor
  !! @param[out]   table_ene   : energy table
  !! @param[out]   table_grad  : gradient table
  !! @param[out]   table_ecor  : energy table (bond correction)
  !! @param[out]   table_decor : gradient table (bond correction)
  !! @note         J.Jung et al., J.Comput.Chem., 34, 2412-2420 (2013)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine table_pme_linear_noswitch(switch_int, cutoff_int, density,        &
                                     cutoff2, switchdist2, el_fact, alpha,     &
                                     alpha2sp, alpha2m,                        &
                                     table_ene, table_grad, table_ecor,        &
                                     table_decor)
    ! formal arguments
    integer,                 intent(in)    :: switch_int
    integer,                 intent(in)    :: cutoff_int
    real(wp),                intent(in)    :: density
    real(wp),                intent(in)    :: cutoff2
    real(wp),                intent(in)    :: switchdist2
    real(wp),                intent(in)    :: el_fact
    real(wp),                intent(in)    :: alpha
    real(wp),                intent(in)    :: alpha2m
    real(wp),                intent(in)    :: alpha2sp
    real(wp),                intent(inout) :: table_ene(:)
    real(wp),                intent(inout) :: table_grad(:)
    real(wp),                intent(inout) :: table_ecor(:)
    real(wp),                intent(inout) :: table_decor(:)

    ! local variables
    integer                  :: i
    real(wp)                 :: rij, rij2
    real(wp)                 :: inv_r2, inv_rij, inv_r6, inv_r12

    do i = int(density)+1, cutoff_int

      rij2     = cutoff2*density/real(i,wp)
      rij      = sqrt(rij2)

      inv_r2   = 1.0_wp / rij2
      inv_rij  = 1.0_wp / rij
      inv_r6   = inv_r2 * inv_r2 * inv_r2
      inv_r12  = inv_r6 * inv_r6

      table_ene(3*i-2)  = inv_r12
      table_ene(3*i-1)  = inv_r6
      table_grad(3*i-2) = -12.0_wp*inv_r12*inv_r2
      table_grad(3*i-1) = -6.0_wp*inv_r6*inv_r2
      table_ene(3*i)    = el_fact*erfc04(alpha*rij)*inv_rij
      table_grad(3*i)   = -(table_ene(3*i)+ &
                            alpha2sp*el_fact*exp(alpha2m*rij2))*inv_r2
      table_ecor(i)     = -erf_n(alpha*rij)*inv_rij*el_fact
      table_decor(i)    = -(table_ecor(i)+ &
                            alpha2sp*el_fact*exp(alpha2m*rij2))/rij2
    end do

    return

  end subroutine table_pme_linear_noswitch

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    table_pme_linear_switch
  !> @brief        table for pme and linear 
  !! @authors      CK, JJ
  !! @param[in]    switch_int  : integer switchdist
  !! @param[in]    cutoff_int  : integer cutoff distance
  !! @param[in]    density     : table density
  !! @param[in]    cutoff2     : square cutoff distance
  !! @param[in]    switchdist2 : square switch distance
  !! @param[in]    el_fact     : inverse of dielectoric constant
  !! @param[in]    alpha       : ewald factor
  !! @param[in]    alpha2sp    : ewald factor
  !! @param[in]    alpha2m     : ewald factor
  !! @param[out]   table_ene   : energy table
  !! @param[out]   table_grad  : gradient table
  !! @param[out]   table_ecor  : energy table (bond correction)
  !! @param[out]   table_decor : gradient table (bond correction)
  !! @note         J.Jung et al., J.Comput.Chem., 34, 2412-2420 (2013)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine table_pme_linear_switch(switch_int, cutoff_int, density,          &
                                     cutoff2, switchdist2, el_fact, alpha,     &
                                     alpha2sp, alpha2m,                        &
                                     table_ene, table_grad, table_ecor,        &
                                     table_decor)
    ! formal arguments
    integer,                 intent(in)    :: switch_int
    integer,                 intent(in)    :: cutoff_int
    real(wp),                intent(in)    :: density
    real(wp),                intent(in)    :: cutoff2
    real(wp),                intent(in)    :: switchdist2
    real(wp),                intent(in)    :: el_fact
    real(wp),                intent(in)    :: alpha
    real(wp),                intent(in)    :: alpha2m
    real(wp),                intent(in)    :: alpha2sp
    real(wp),                intent(inout) :: table_ene(:)
    real(wp),                intent(inout) :: table_grad(:)
    real(wp),                intent(inout) :: table_ecor(:)
    real(wp),                intent(inout) :: table_decor(:)

    ! local variables
    integer                  :: i
    real(wp)                 :: rij, rij2
    real(wp)                 :: inv_r2, inv_rij, inv_r6, inv_r12
    real(wp)                 :: lj1, lj2, lj4, esw, dsw
    real(wp)                 :: switch, dswitch

    lj1 = 1.0_wp/(cutoff2-switchdist2)
    lj1 = lj1*lj1*lj1

    do i = switch_int+1, cutoff_int
      rij2    = cutoff2*density/real(i,wp)
      rij     = sqrt(rij2)
      inv_r2  = 1.0_wp / rij2
      inv_rij = 1.0_wp / rij
      inv_r6  = inv_r2 * inv_r2 * inv_r2
      inv_r12 = inv_r6 * inv_r6

      table_ene(3*i-2)  = inv_r12
      table_ene(3*i-1)  = inv_r6
      table_grad(3*i-2) = -12.0_wp*inv_r12*inv_r2
      table_grad(3*i-1) = -6.0_wp*inv_r6*inv_r2
      table_ene(3*i)    = el_fact*erfc04(alpha*rij)*inv_rij
      table_grad(3*i)   = -(table_ene(3*i)+ &
                            alpha2sp*el_fact*exp(alpha2m*rij2))*inv_r2
      table_ecor(i)     = -erf_n(alpha*rij)*inv_rij*el_fact
      table_decor(i)    = -(table_ecor(i)+ &
                            alpha2sp*el_fact*exp(alpha2m*rij2))/rij2
    end do

    do i = int(density), switch_int
      rij2     = cutoff2*density/real(i,wp)
      rij      = sqrt(rij2)

      lj2      = cutoff2 - rij2
      lj4      = lj2*(cutoff2 + 2.0_wp*rij2 - 3.0_wp*switchdist2)
      switch   = lj2*lj4*lj1
      dswitch  = 4.0_wp*lj1*rij  *(lj2*lj2 - lj4)

      inv_r2   = 1.0_wp / rij2
      inv_rij  = 1.0_wp / rij
      inv_r6   = inv_r2 * inv_r2 * inv_r2
      inv_r12  = inv_r6 * inv_r6

      table_ene(3*i-2)  = switch*inv_r12
      table_ene(3*i-1)  = switch*inv_r6
      table_grad(3*i-2) = dswitch*inv_r12*inv_rij              &
                         - switch*12.0_wp*inv_r12*inv_r2
      table_grad(3*i-1) = dswitch*inv_r6*inv_rij               &
                         - switch*6.0_wp*inv_r6*inv_r2
      table_ene(3*i)    = el_fact*erfc04(alpha*rij)*inv_rij
      table_grad(3*i)   = -(table_ene(3*i)+ &
                            alpha2sp*el_fact*exp(alpha2m*rij2))*inv_r2
    end do

    return

  end subroutine table_pme_linear_switch

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    table_pme_linear_fswitch
  !> @brief        table for pme and linear and force switch
  !! @authors      CK, JJ
  !! @param[in]    switch_int  : integer switchdist
  !! @param[in]    cutoff_int  : integer cutoff distance
  !! @param[in]    density     : table density
  !! @param[in]    cutoff2     : square cutoff distance
  !! @param[in]    switchdist2 : square switch distance
  !! @param[in]    el_fact     : inverse of dielectoric constant
  !! @param[in]    alpha       : ewald factor
  !! @param[in]    alpha2sp    : ewald factor
  !! @param[in]    alpha2m     : ewald factor
  !! @param[out]   table_ene   : energy table
  !! @param[out]   table_grad  : gradient table
  !! @param[out]   table_ecor  : energy table (bond correction)
  !! @param[out]   table_decor : gradient table (bond correction)
  !! @note         J.Jung et al., J.Comput.Chem., 34, 2412-2420 (2013)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine table_pme_linear_fswitch(switch_int, cutoff_int, density,      &
                                     cutoff2, switchdist2, el_fact, alpha,     &
                                     alpha2sp, alpha2m,                        &
                                     table_ene, table_grad, table_ecor,        &
                                     table_decor)
    ! formal arguments
    integer,                 intent(in)    :: switch_int
    integer,                 intent(in)    :: cutoff_int
    real(wp),                intent(in)    :: density
    real(wp),                intent(in)    :: cutoff2
    real(wp),                intent(in)    :: switchdist2
    real(wp),                intent(in)    :: el_fact
    real(wp),                intent(in)    :: alpha
    real(wp),                intent(in)    :: alpha2m
    real(wp),                intent(in)    :: alpha2sp
    real(wp),                intent(inout) :: table_ene(:)
    real(wp),                intent(inout) :: table_grad(:)
    real(wp),                intent(inout) :: table_ecor(:)
    real(wp),                intent(inout) :: table_decor(:)

    ! local variables
    integer                  :: i
    real(wp)                 :: rij, rij2
    real(wp)                 :: inv_r2, inv_rij, inv_r3, inv_r6, inv_r12
    real(wp)                 :: lj1, lj2, lj4, lj6, lj12
    real(wp)                 :: switch, dswitch

    lj1 = 1.0_wp / (cutoff2*switchdist2)
    lj1 = lj1*lj1*lj1
    lj2 = sqrt(lj1)

    do i = switch_int+1, cutoff_int
      rij2    = cutoff2*density/real(i,wp)
      rij     = sqrt(rij2)
      inv_r2  = 1.0_wp / rij2
      inv_rij = 1.0_wp / rij
      inv_r6  = inv_r2 * inv_r2 * inv_r2
      inv_r12 = inv_r6 * inv_r6

      table_ene(3*i-2)  = inv_r12 - lj1
      table_ene(3*i-1)  = inv_r6 - lj2
      table_grad(3*i-2) = -12.0_wp*inv_r12*inv_r2
      table_grad(3*i-1) = -6.0_wp*inv_r6*inv_r2
      table_ene(3*i)    = el_fact*erfc04(alpha*rij)*inv_rij
      table_grad(3*i)   = -(table_ene(3*i)+ &
                            alpha2sp*el_fact*exp(alpha2m*rij2))*inv_r2
      table_ecor(i)     = -erf_n(alpha*rij)*inv_rij*el_fact
      table_decor(i)    = -(table_ecor(i)+ &
                            alpha2sp*el_fact*exp(alpha2m*rij2))/rij2
    end do

    lj1  = cutoff2*cutoff2*cutoff2
    lj2  = cutoff2*sqrt(cutoff2)
    lj12 = lj1 / (lj1-switchdist2*switchdist2*switchdist2)
    lj6  = lj2 / (lj2-switchdist2*sqrt(switchdist2))
    lj1  = 1.0_wp / lj1
    lj2  = 1.0_wp / lj2

    do i = int(density), switch_int

      rij2    = cutoff2*density/real(i,wp)
      rij     = sqrt(rij2)
      inv_r2  = 1.0_wp / rij2
      inv_rij = 1.0_wp / rij
      inv_r6  = inv_r2 * inv_r2 * inv_r2
      inv_r3  = sqrt(inv_r6)

      table_ene(3*i-2)  = lj12*(inv_r6-lj1)**2
      table_ene(3*i-1)  = lj6*(inv_r3-lj2)**2
      table_grad(3*i-2) = -12.0_wp*lj12*(inv_r6-lj1)*inv_r6*inv_r2
      table_grad(3*i-1) = -6.0_wp*lj6*(inv_r3-lj2)*inv_r3*inv_r2
      table_ene(3*i)    = el_fact*erfc04(alpha*rij)*inv_rij
      table_grad(3*i)   = -(table_ene(3*i)+ &
                            alpha2sp*el_fact*exp(alpha2m*rij2))*inv_r2
    end do

    return

  end subroutine table_pme_linear_fswitch

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    table_pme_linear_groshift
  !> @brief        table for pme and linear and gromacs shift
  !! @authors      CK, JJ
  !! @param[in]    switch_int  : integer switchdist
  !! @param[in]    cutoff_int  : integer cutoff distance
  !! @param[in]    density     : table density
  !! @param[in]    cutoff      : cutoff distance
  !! @param[in]    switchdist  : switch distance
  !! @param[in]    el_fact     : inverse of dielectoric constant
  !! @param[in]    alpha       : ewald factor
  !! @param[in]    alpha2sp    : ewald factor
  !! @param[in]    alpha2m     : ewald factor
  !! @param[out]   table_ene   : energy table
  !! @param[out]   table_grad  : gradient table
  !! @param[out]   table_ecor  : energy table (bond correction)
  !! @param[out]   table_decor : gradient table (bond correction)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine table_pme_linear_groshift(switch_int, cutoff_int, density,        &
                                       cutoff, switchdist, el_fact, alpha,     &
                                       alpha2sp, alpha2m,                      &
                                       table_ene, table_grad, table_ecor,      &
                                       table_decor)
    ! formal arguments
    integer,                 intent(in)    :: switch_int
    integer,                 intent(in)    :: cutoff_int
    real(wp),                intent(in)    :: density
    real(wp),                intent(in)    :: cutoff
    real(wp),                intent(in)    :: switchdist
    real(wp),                intent(in)    :: el_fact
    real(wp),                intent(in)    :: alpha
    real(wp),                intent(in)    :: alpha2m
    real(wp),                intent(in)    :: alpha2sp
    real(wp),                intent(inout) :: table_ene(:)
    real(wp),                intent(inout) :: table_grad(:)
    real(wp),                intent(inout) :: table_ecor(:)
    real(wp),                intent(inout) :: table_decor(:)

    ! local variables
    integer                  :: i
    real(wp)                 :: rij, rij2
    real(wp)                 :: inv_r2, inv_rij, inv_r6, inv_r12
    real(wp)                 :: lj1, lj2, lj4, cutoff2, switchdist2
    real(wp)                 :: switch, dswitch
    real(wp)                 :: Coef_A12, Coef_B12, Coef_C12
    real(wp)                 :: Coef_A06, Coef_B06, Coef_C06
    real(wp)                 :: term_lj12_OO, term_lj12_HO, term_lj12_HH
    real(wp)                 :: term_lj6_OO, term_lj6_HO, term_lj6_HH
    real(wp)                 :: F12, F06, P12, P06

    cutoff2 = cutoff*cutoff
    switchdist2 = switchdist*switchdist

    lj1 = 1.0_wp/(cutoff2-switchdist2)

    Coef_A12  = 16.0_wp*cutoff - 13.0_wp*switchdist
    Coef_A12  = - Coef_A12 / cutoff**14 / (cutoff-switchdist)**2
    Coef_B12  = 15.0_wp*cutoff - 13.0_wp*switchdist
    Coef_B12  = Coef_B12 / cutoff**14 / (cutoff-switchdist)**3
    Coef_A06  = 10.0_wp*cutoff - 7.0_wp*switchdist
    Coef_A06  = - Coef_A06 / cutoff**8 / (cutoff-switchdist)**2
    Coef_B06  = 9.0_wp*cutoff - 7.0_wp*switchdist
    Coef_B06  = Coef_B06 / cutoff**8 / (cutoff-switchdist)**3
    Coef_C12  = 1.0_wp/cutoff**12
    Coef_C12  = Coef_C12 - 12.0_wp*Coef_A12/3.0_wp*(cutoff-switchdist)**3
    Coef_C12  = Coef_C12 - 12.0_wp*Coef_B12/4.0_wp*(cutoff-switchdist)**4
    Coef_C06  = 1.0_wp/cutoff**6
    Coef_C06  = Coef_C06 - 6.0_wp*Coef_A06/3.0_wp*(cutoff-switchdist)**3
    Coef_C06  = Coef_C06 - 6.0_wp*Coef_B06/4.0_wp*(cutoff-switchdist)**4

    if (switchdist .le. EPS) then

      do i = int(density), cutoff_int

        rij2      = cutoff2*density/real(i,wp)
        rij       = sqrt(rij2)
        F12       = Coef_A12*(rij-switchdist)**2 + Coef_B12*(rij-switchdist)**3
        F06       = Coef_A06*(rij-switchdist)**2 + Coef_B06*(rij-switchdist)**3
        P12       = - 12.0_wp*Coef_A12/3.0_wp*(rij-switchdist)**3
        P12       = P12 - 12.0_wp*Coef_B12/4.0_wp*(rij-switchdist)**4 - Coef_C12
        P06       = - 6.0_wp*Coef_A06/3.0_wp*(rij-switchdist)**3
        P06       = P06 - 6.0_wp*Coef_B06/4.0_wp*(rij-switchdist)**4 - Coef_C06
        inv_r2    = 1.0_wp/rij2
        inv_rij   = 1.0_wp/rij
        inv_r6    = inv_r2 * inv_r2 * inv_r2
        inv_r12   = inv_r6 * inv_r6
        table_ene(3*i-2)  = inv_r12 + P12
        table_ene(3*i-1)  = inv_r6 + P06
        table_grad(3*i-2) = -12.0_wp*inv_r12*inv_r2 - 12.0_wp*F12*inv_rij
        table_grad(3*i-1) = -6.0_wp*inv_r6*inv_r2 - 6.0_wp*F06*inv_rij
        table_ene(3*i)    = el_fact*erfc04(alpha*rij)*inv_rij
        table_grad(3*i)   = -(table_ene(3*i)+ &
                            alpha2sp*el_fact*exp(alpha2m*rij2))*inv_r2
        table_ecor(i)     = -erf_n(alpha*rij)*inv_rij*el_fact
        table_decor(i)    = -(table_ecor(i)+ &
                            alpha2sp*el_fact*exp(alpha2m*rij2))/rij2
      end do

    else

      do i = int(density)+1, switch_int
  
        rij2      = cutoff2*density/real(i,wp)
        rij       = sqrt(rij2)
        F12       = Coef_A12*(rij-switchdist)**2 + Coef_B12*(rij-switchdist)**3
        F06       = Coef_A06*(rij-switchdist)**2 + Coef_B06*(rij-switchdist)**3
        P12       = - 12.0_wp*Coef_A12/3.0_wp*(rij-switchdist)**3
        P12       = P12 - 12.0_wp*Coef_B12/4.0_wp*(rij-switchdist)**4 - Coef_C12
        P06       = - 6.0_wp*Coef_A06/3.0_wp*(rij-switchdist)**3
        P06       = P06 - 6.0_wp*Coef_B06/4.0_wp*(rij-switchdist)**4 - Coef_C06
        inv_r2    = 1.0_wp/rij2
        inv_rij   = 1.0_wp/rij
        inv_r6    = inv_r2 * inv_r2 * inv_r2
        inv_r12   = inv_r6 * inv_r6
        table_ene(3*i-2)  = inv_r12 + P12
        table_ene(3*i-1)  = inv_r6 + P06
        table_grad(3*i-2) = -12.0_wp*inv_r12*inv_r2 - 12.0_wp*F12*inv_rij
        table_grad(3*i-1) = -6.0_wp*inv_r6*inv_r2 - 6.0_wp*F06*inv_rij
        table_ene(3*i)    = el_fact*erfc04(alpha*rij)*inv_rij
        table_grad(3*i)   = -(table_ene(3*i)+ &
                            alpha2sp*el_fact*exp(alpha2m*rij2))*inv_r2
        table_ecor(i)     = -erf_n(alpha*rij)*inv_rij*el_fact
        table_decor(i)    = -(table_ecor(i)+ &
                            alpha2sp*el_fact*exp(alpha2m*rij2))/rij2
      end do

      do i = switch_int+1, cutoff_int

        rij2     = cutoff2*density/real(i,wp)
        rij      = sqrt(rij2)

        inv_r2   = 1.0_wp / rij2
        inv_rij  = 1.0_wp / rij
        inv_r6   = inv_r2 * inv_r2 * inv_r2
        inv_r12  = inv_r6 * inv_r6

        table_ene(3*i-2)  = inv_r12 - Coef_C12
        table_ene(3*i-1)  = inv_r6 - Coef_C06
        table_grad(3*i-2) = -12.0_wp*inv_r12*inv_r2
        table_grad(3*i-1) = -6.0_wp*inv_r6*inv_r2
        table_ene(3*i)    = el_fact*erfc04(alpha*rij)*inv_rij
        table_grad(3*i)   = -(table_ene(3*i)+ &
                            alpha2sp*el_fact*exp(alpha2m*rij2))*inv_r2
        table_ecor(i)     = -erf_n(alpha*rij)*inv_rij*el_fact
        table_decor(i)    = -(table_ecor(i)+ &
                            alpha2sp*el_fact*exp(alpha2m*rij2))/rij2
      end do

    end if

    return

  end subroutine table_pme_linear_groshift

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    table_pme_linear_groswitch
  !> @brief        table for pme and gromacs switch with linear interpolation
  !! @authors      CK, JJ
  !! @param[in]    switch_int  : integer switchdist
  !! @param[in]    cutoff_int  : integer cutoff distance
  !! @param[in]    density     : table density
  !! @param[in]    cutoff      : cutoff distance
  !! @param[in]    switchdist  : switch distance
  !! @param[in]    el_fact     : inverse of dielectoric constant
  !! @param[in]    alpha       : ewald factor
  !! @param[in]    alpha2sp    : ewald factor
  !! @param[in]    alpha2m     : ewald factor
  !! @param[out]   table_ene   : energy table
  !! @param[out]   table_grad  : gradient table
  !! @param[out]   table_ecor  : energy table (bond correction)
  !! @param[out]   table_decor : gradient table (bond correction)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine table_pme_linear_groswitch(switch_int, cutoff_int, density,       &
                                        cutoff, switchdist, el_fact, alpha,    &
                                        alpha2sp, alpha2m,                     &
                                        table_ene, table_grad, table_ecor,     &
                                        table_decor, eswitch, vswitch)
    ! formal arguments
    integer,                 intent(in)    :: switch_int
    integer,                 intent(in)    :: cutoff_int
    real(wp),                intent(in)    :: density
    real(wp),                intent(in)    :: cutoff
    real(wp),                intent(in)    :: switchdist
    real(wp),                intent(in)    :: el_fact
    real(wp),                intent(in)    :: alpha
    real(wp),                intent(in)    :: alpha2m
    real(wp),                intent(in)    :: alpha2sp
    real(wp),                intent(inout) :: table_ene(:)
    real(wp),                intent(inout) :: table_grad(:)
    real(wp),                intent(inout) :: table_ecor(:)
    real(wp),                intent(inout) :: table_decor(:)
    real(wp),                intent(inout) :: eswitch, vswitch

    ! local variables
    integer                  :: i, i_min, i_max
    real(wp)                 :: rij, rij2, rij1, rij_min, rij_max, delta_rij
    real(wp)                 :: inv_r2, inv_rij, inv_r6, inv_r12, delta_f
    real(wp)                 :: cutoff2, switchdist2
    real(wp)                 :: switch, dswitch
    real(wp)                 :: coef_A, coef_B, coef_C
    real(wp)                 :: coef_A1, coef_B1, coef_C1

    cutoff2 = cutoff*cutoff
    switchdist2 = switchdist*switchdist

    coef_A    = -10.0_wp / (cutoff-switchdist)**3
    coef_B    =  15.0_wp / (cutoff-switchdist)**4
    coef_C    = - 6.0_wp / (cutoff-switchdist)**5
    coef_A1   = -30.0_wp / (cutoff-switchdist)**3
    coef_B1   =  60.0_wp / (cutoff-switchdist)**4
    coef_C1   = -30.0_wp / (cutoff-switchdist)**5

    do i = int(density), switch_int

      rij2      = cutoff2*density/real(i,wp)
      rij       = sqrt(rij2)
      inv_rij   = 1.0_wp / rij
      rij1      = rij - switchdist
      switch    = 1.0_wp + coef_A*rij1**3 + coef_B*rij1**4 + coef_C*rij1**5
      dswitch   = coef_A1*rij1**2 + coef_B*rij1**3 + coef_C*rij1**4
      inv_r2    = 1.0_wp/rij2
      inv_rij   = 1.0_wp/rij
      inv_r6    = inv_r2 * inv_r2 * inv_r2
      inv_r12   = inv_r6 * inv_r6
      table_ene(3*i-2)  = inv_r12 * switch
      table_ene(3*i-1)  = inv_r6 * switch
      table_grad(3*i-2) = -12.0_wp*inv_r12*inv_r2*switch  &
                          + dswitch*inv_r12*inv_rij
      table_grad(3*i-1) = -6.0_wp*inv_r6*inv_r2*switch    &
                          + dswitch*inv_r6*inv_rij
      table_ene(3*i)    = el_fact*erfc04(alpha*rij)*inv_rij
      table_grad(3*i)   = -(table_ene(3*i)+ &
                            alpha2sp*el_fact*exp(alpha2m*rij2))*inv_r2
    end do

    do i = switch_int+1, cutoff_int

      rij2     = cutoff2*density/real(i,wp)
      rij      = sqrt(rij2)

      inv_r2   = 1.0_wp / rij2
      inv_rij  = 1.0_wp / rij
      inv_r6   = inv_r2 * inv_r2 * inv_r2
      inv_r12  = inv_r6 * inv_r6

      table_ene(3*i-2)  = inv_r12
      table_ene(3*i-1)  = inv_r6
      table_grad(3*i-2) = -12.0_wp*inv_r12*inv_r2
      table_grad(3*i-1) = -6.0_wp*inv_r6*inv_r2
      table_ene(3*i)    = el_fact*erfc04(alpha*rij)*inv_rij
      table_grad(3*i)   = -(table_ene(3*i)+ &
                            alpha2sp*el_fact*exp(alpha2m*rij2))*inv_r2
      table_ecor(i)     = -erf_n(alpha*rij)*inv_rij*el_fact
      table_decor(i)    = -(table_ecor(i)+ &
                            alpha2sp*el_fact*exp(alpha2m*rij2))/rij2
    end do

    ! dispersion correction
    !
    eswitch = 0.0_wp
    vswitch = 0.0_wp
    i_min   = int(switchdist * 1000.0_wp)
    i_max   = int(cutoff * 1000.0_wp)
    do i = i_min, i_max-1
      rij       = (real(i,wp)+0.5) / 1000.0_wp
      rij2      = rij*rij
      inv_rij   = 1.0_wp / rij
      rij1      = rij - switchdist
      switch    = 1.0_wp + coef_A*rij1**3 + coef_B*rij1**4 + coef_C*rij1**5
      dswitch   = coef_A1*rij1**2 + coef_B*rij1**3 + coef_C*rij1**4
      inv_r2    = 1.0_wp/rij2
      inv_r6    = inv_r2 * inv_r2 * inv_r2
      inv_r12   = inv_r6 * inv_r6
      rij_min   = real(i,wp) / 1000.0_wp
      rij_max   = real(i+1,wp) / 1000.0_wp
      delta_rij = rij_max - rij_min
      eswitch   = eswitch-rij2*inv_r6*(1.0_wp-switch)*delta_rij
!     delta_f   = 6.0_wp*inv_r6 - 6.0_wp*inv_r6*switch + dswitch*inv_r6*rij
      delta_f   = 6.0_wp*inv_r6 - 6.0_wp*inv_r6*switch 
      vswitch   = vswitch + 0.5_wp*delta_f*delta_rij*rij*rij
    end do
    vswitch   = vswitch / 3.0_wp

    return

  end subroutine table_pme_linear_groswitch


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    table_water_pme_linear_noswitch
  !> @brief        table for pme and linear without switch
  !! @authors      CK, JJ
  !! @param[in]    switch_int  : integer switchdist
  !! @param[in]    cutoff_int  : integer cutoff distance
  !! @param[in]    density     : table density
  !! @param[in]    cutoff2     : square cutoff distance
  !! @param[in]    switchdist2 : square switch distance
  !! @param[in]    el_fact     : inverse of dielectoric constant
  !! @param[in]    alpha       : ewald factor
  !! @param[in]    alpha2sp    : ewald factor
  !! @param[in]    alpha2m     : ewald factor
  !! @param[in]    lj6_OO      : LJ6 oxygen-oxygen
  !! @param[in]    lj6_OH      : LJ6 oxygen-hydrogen
  !! @param[in]    lj6_HH      : LJ6 hydrogen-hydrogen
  !! @param[in]    lj12_OO     : LJ12 oxygen-oxygen
  !! @param[in]    lj12_OH     : LJ12 oxygen-hydrogen
  !! @param[in]    lj12_HH     : LJ12 hydrogen-hydrogen
  !! @param[in]    cc_OO       : charge oxygen-oxygen
  !! @param[in]    cc_OH       : charge oxygen-hydrogen
  !! @param[in]    cc_HH       : charge hydrogen-hydrogen
  !! @param[out]   table_ene_WW: energy table
  !! @param[out]   table_de_WW  : gradient table
  !! @note         J.Jung et al., J.Comput.Chem., 34, 2412-2420 (2013)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine table_water_pme_linear_noswitch(switch_int, cutoff_int, density,  &
                                     cutoff2, switchdist2, el_fact, alpha,     &
                                     lj6_OO,  lj6_OH,  lj6_HH,                 &
                                     lj12_OO, lj12_OH, lj12_HH,                &
                                     cc_OO, cc_OH, cc_HH,                      &
                                     alpha2sp, alpha2m,                        &
                                     table_ene_WW, table_de_WW)
    ! formal arguments
    integer,                 intent(in)    :: switch_int
    integer,                 intent(in)    :: cutoff_int
    real(wp),                intent(in)    :: density
    real(wp),                intent(in)    :: cutoff2
    real(wp),                intent(in)    :: switchdist2
    real(wp),                intent(in)    :: el_fact
    real(wp),                intent(in)    :: alpha
    real(wp),                intent(in)    :: alpha2m
    real(wp),                intent(in)    :: alpha2sp
    real(wp),                intent(in)    :: lj6_OO
    real(wp),                intent(in)    :: lj6_OH
    real(wp),                intent(in)    :: lj6_HH
    real(wp),                intent(in)    :: lj12_OO
    real(wp),                intent(in)    :: lj12_OH
    real(wp),                intent(in)    :: lj12_HH
    real(wp),                intent(in)    :: cc_OO
    real(wp),                intent(in)    :: cc_OH
    real(wp),                intent(in)    :: cc_HH
    real(wp),                intent(inout) :: table_ene_WW(:,:)
    real(wp),                intent(inout) :: table_de_WW(:,:)

    ! local variables
    integer                  :: i
    real(wp)                 :: rij, rij2
    real(wp)                 :: inv_r2, inv_rij, inv_r6, inv_r12
    real(wp)                 :: term_lj12_OO, term_lj12_OH, term_lj12_HH
    real(wp)                 :: term_lj6_OO, term_lj6_OH, term_lj6_HH
    real(wp)                 :: term_lj_OO, term_lj_OH, term_lj_HH

    do i = int(density)+1, cutoff_int
      rij2     = cutoff2*density/real(i,wp)
      rij      = sqrt(rij2)
      inv_r2   = 1.0_wp/rij2
      inv_rij  = 1.0_wp / rij
      inv_r6   = inv_r2 * inv_r2 * inv_r2
      inv_r12  = inv_r6 * inv_r6

      term_lj12_OO  = lj12_OO * inv_r12
      term_lj12_OH  = lj12_OH * inv_r12
      term_lj12_HH  = lj12_HH * inv_r12
      term_lj6_OO   = lj6_OO * inv_r6
      term_lj6_OH   = lj6_OH * inv_r6
      term_lj6_HH   = lj6_HH * inv_r6

      term_lj_OO = term_lj12_OO - term_lj6_OO
      term_lj_OH = term_lj12_OH - term_lj6_OH
      term_lj_HH = term_lj12_HH - term_lj6_HH

      table_ene_WW(6*i-5,1) = term_lj_OO
      table_ene_WW(6*i-5,2) = term_lj_OH
      table_ene_WW(6*i-5,3) = term_lj_HH

      table_de_WW(i,1)   = -(12.0_wp*term_lj12_OO-6.0_wp*term_lj6_OO)*inv_r2
      table_de_WW(i,2)   = -(12.0_wp*term_lj12_OH-6.0_wp*term_lj6_OH)*inv_r2
      table_de_WW(i,3)   = -(12.0_wp*term_lj12_HH-6.0_wp*term_lj6_HH)*inv_r2

      table_ene_WW(6*i-4,1) = cc_OO*erfc04(alpha*rij)*inv_rij
      table_ene_WW(6*i-4,2) = cc_OH*erfc04(alpha*rij)*inv_rij
      table_ene_WW(6*i-4,3) = cc_HH*erfc04(alpha*rij)*inv_rij

      table_de_WW(i,1)   = table_de_WW(i,1) &
              -(table_ene_WW(6*i-4,1)+alpha2sp*cc_OO*exp(alpha2m*rij2))*inv_r2
      table_de_WW(i,2)   = table_de_WW(i,2) &
              -(table_ene_WW(6*i-4,2)+alpha2sp*cc_OH*exp(alpha2m*rij2))*inv_r2
      table_de_WW(i,3)   = table_de_WW(i,3) &
              -(table_ene_WW(6*i-4,3)+alpha2sp*cc_HH*exp(alpha2m*rij2))*inv_r2

      table_ene_WW(6*i-3,1:3) = table_de_WW(i,1:3)

    end do

    do i = int(density), cutoff_int
      table_ene_WW(6*i-2,1:3) = table_ene_WW(6*(i+1)-5,1:3) - &
                                table_ene_WW(6*i-5,1:3)
      table_ene_WW(6*i-1,1:3) = table_ene_WW(6*(i+1)-4,1:3) - &
                                table_ene_WW(6*i-4,1:3)
      table_ene_WW(6*i,1:3)   = table_ene_WW(6*(i+1)-3,1:3) - &
                                table_ene_WW(6*i-3,1:3)
    end do

    return

  end subroutine table_water_pme_linear_noswitch

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    table_water_pme_linear_switch
  !> @brief        table for pme and linear and switch
  !! @authors      CK, JJ
  !! @param[in]    switch_int  : integer switchdist
  !! @param[in]    cutoff_int  : integer cutoff distance
  !! @param[in]    density     : table density
  !! @param[in]    cutoff2     : square cutoff distance
  !! @param[in]    switchdist2 : square switch distance
  !! @param[in]    el_fact     : inverse of dielectoric constant
  !! @param[in]    alpha       : ewald factor
  !! @param[in]    alpha2sp    : ewald factor
  !! @param[in]    alpha2m     : ewald factor
  !! @param[in]    lj6_OO      : LJ6 oxygen-oxygen
  !! @param[in]    lj6_OH      : LJ6 oxygen-hydrogen
  !! @param[in]    lj6_HH      : LJ6 hydrogen-hydrogen
  !! @param[in]    lj12_OO     : LJ12 oxygen-oxygen
  !! @param[in]    lj12_OH     : LJ12 oxygen-hydrogen
  !! @param[in]    lj12_HH     : LJ12 hydrogen-hydrogen
  !! @param[in]    cc_OO       : charge oxygen-oxygen
  !! @param[in]    cc_OH       : charge oxygen-hydrogen
  !! @param[in]    cc_HH       : charge hydrogen-hydrogen
  !! @param[out]   table_ene_WW: energy table
  !! @param[out]   table_de_WW  : gradient table
  !! @note         J.Jung et al., J.Comput.Chem., 34, 2412-2420 (2013)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine table_water_pme_linear_switch(switch_int, cutoff_int, density,    &
                                     cutoff2, switchdist2, el_fact, alpha,     &
                                     lj6_OO,  lj6_OH,  lj6_HH,                 &
                                     lj12_OO, lj12_OH, lj12_HH,                &
                                     cc_OO, cc_OH, cc_HH,                      &
                                     alpha2sp, alpha2m,                        &
                                     table_ene_WW, table_de_WW)
    ! formal arguments
    integer,                 intent(in)    :: switch_int
    integer,                 intent(in)    :: cutoff_int
    real(wp),                intent(in)    :: density
    real(wp),                intent(in)    :: cutoff2
    real(wp),                intent(in)    :: switchdist2
    real(wp),                intent(in)    :: el_fact
    real(wp),                intent(in)    :: alpha
    real(wp),                intent(in)    :: alpha2m
    real(wp),                intent(in)    :: alpha2sp
    real(wp),                intent(in)    :: lj6_OO
    real(wp),                intent(in)    :: lj6_OH
    real(wp),                intent(in)    :: lj6_HH
    real(wp),                intent(in)    :: lj12_OO
    real(wp),                intent(in)    :: lj12_OH
    real(wp),                intent(in)    :: lj12_HH
    real(wp),                intent(in)    :: cc_OO
    real(wp),                intent(in)    :: cc_OH
    real(wp),                intent(in)    :: cc_HH
    real(wp),                intent(inout) :: table_ene_WW(:,:)
    real(wp),                intent(inout) :: table_de_WW(:,:)

    ! local variables
    integer                  :: i
    real(wp)                 :: rij, rij2
    real(wp)                 :: inv_r2, inv_rij, inv_r6, inv_r12
    real(wp)                 :: lj1, lj2, lj4
    real(wp)                 :: term_lj12_OO, term_lj12_OH, term_lj12_HH
    real(wp)                 :: term_lj6_OO, term_lj6_OH, term_lj6_HH
    real(wp)                 :: term_lj_OO, term_lj_OH, term_lj_HH
    real(wp)                 :: switch, dswitch

    lj1 = 1.0_wp/(cutoff2-switchdist2)
    lj1 = lj1*lj1*lj1

    do i = switch_int+1, cutoff_int

      rij2     = cutoff2*density/real(i,wp)
      rij    = sqrt(rij2)
      switch  = 1.0_wp
      dswitch = 0.0_wp

      inv_r2  = 1.0_wp / rij2
      inv_rij = 1.0_wp / rij
      inv_r6  = inv_r2 * inv_r2 * inv_r2
      inv_r12 = inv_r6 * inv_r6

      term_lj12_OO = lj12_OO * inv_r12
      term_lj12_OH = lj12_OH * inv_r12
      term_lj12_HH = lj12_HH * inv_r12
      term_lj6_OO  = lj6_OO * inv_r6
      term_lj6_OH  = lj6_OH * inv_r6
      term_lj6_HH  = lj6_HH * inv_r6

      term_lj_OO = term_lj12_OO - term_lj6_OO
      term_lj_OH = term_lj12_OH - term_lj6_OH
      term_lj_HH = term_lj12_HH - term_lj6_HH

      table_ene_WW(6*i-5,1) = switch*term_lj_OO
      table_ene_WW(6*i-5,2) = switch*term_lj_OH
      table_ene_WW(6*i-5,3) = switch*term_lj_HH

      table_de_WW(i,1)   = (  dswitch * term_lj_OO * rij  &
              - switch * (12.0_wp*term_lj12_OO-6.0_wp*term_lj6_OO))*inv_r2
      table_de_WW(i,2)   = (  dswitch * term_lj_OH * rij  &
              - switch * (12.0_wp*term_lj12_OH-6.0_wp*term_lj6_OH))*inv_r2
      table_de_WW(i,3)   = (  dswitch * term_lj_HH * rij  &
              - switch * (12.0_wp*term_lj12_HH-6.0_wp*term_lj6_HH))*inv_r2

      table_ene_WW(6*i-4,1) = cc_OO*erfc04(alpha*rij)*inv_rij
      table_ene_WW(6*i-4,2) = cc_OH*erfc04(alpha*rij)*inv_rij
      table_ene_WW(6*i-4,3) = cc_HH*erfc04(alpha*rij)*inv_rij

      table_de_WW(i,1)   = table_de_WW(i,1)       &
                         - (table_ene_WW(6*i-4,1) &
                         + alpha2sp*cc_OO*exp(alpha2m*rij2))*inv_r2
      table_de_WW(i,2)   = table_de_WW(i,2)       &
                         - (table_ene_WW(6*i-4,2) &
                         + alpha2sp*cc_OH*exp(alpha2m*rij2))*inv_r2
      table_de_WW(i,3)   = table_de_WW(i,3)       &
                         - (table_ene_WW(6*i-4,3) &
                         +alpha2sp*cc_HH*exp(alpha2m*rij2))*inv_r2

      table_ene_WW(6*i-3,1:3) = table_de_WW(i,1:3)

    end do

    do i = int(density), switch_int
      rij2     = cutoff2*density/real(i,wp)
      rij      = sqrt(rij2)
      lj2      = cutoff2 - rij2
      lj4      = lj2*(cutoff2 + 2.0_wp*rij2 - 3.0_wp*switchdist2)
      switch   = lj2*lj4*lj1
      dswitch  = 4.0_wp*lj1*rij*(lj2*lj2 - lj4)
      inv_r2   = 1.0_wp/rij2
      inv_rij  = 1.0_wp / rij
      inv_r6   = inv_r2 * inv_r2 * inv_r2
      inv_r12  = inv_r6 * inv_r6

      term_lj12_OO  = lj12_OO * inv_r12
      term_lj12_OH  = lj12_OH * inv_r12
      term_lj12_HH  = lj12_HH * inv_r12
      term_lj6_OO   = lj6_OO * inv_r6
      term_lj6_OH   = lj6_OH * inv_r6
      term_lj6_HH   = lj6_HH * inv_r6

      term_lj_OO = term_lj12_OO - term_lj6_OO
      term_lj_OH = term_lj12_OH - term_lj6_OH
      term_lj_HH = term_lj12_HH - term_lj6_HH

      table_ene_WW(6*i-5,1) = switch*term_lj_OO
      table_ene_WW(6*i-5,2) = switch*term_lj_OH
      table_ene_WW(6*i-5,3) = switch*term_lj_HH

      table_de_WW(i,1)   = (  dswitch * term_lj_OO * rij     &
              - switch*(12.0_wp*term_lj12_OO-6.0_wp*term_lj6_OO))*inv_r2
      table_de_WW(i,2)   = (  dswitch * term_lj_OH * rij     &
              - switch * (12.0_wp*term_lj12_OH-6.0_wp*term_lj6_OH))*inv_r2
      table_de_WW(i,3)   = (  dswitch * term_lj_HH * rij     &
              - switch * (12.0_wp*term_lj12_HH-6.0_wp*term_lj6_HH))*inv_r2

      table_ene_WW(6*i-4,1) = cc_OO*erfc04(alpha*rij)*inv_rij
      table_ene_WW(6*i-4,2) = cc_OH*erfc04(alpha*rij)*inv_rij
      table_ene_WW(6*i-4,3) = cc_HH*erfc04(alpha*rij)*inv_rij

      table_de_WW(i,1)   = table_de_WW(i,1) &
              -(table_ene_WW(6*i-4,1)+alpha2sp*cc_OO*exp(alpha2m*rij2))*inv_r2
      table_de_WW(i,2)   = table_de_WW(i,2) &
              -(table_ene_WW(6*i-4,2)+alpha2sp*cc_OH*exp(alpha2m*rij2))*inv_r2
      table_de_WW(i,3)   = table_de_WW(i,3) &
              -(table_ene_WW(6*i-4,3)+alpha2sp*cc_HH*exp(alpha2m*rij2))*inv_r2

      table_ene_WW(6*i-3,1:3) = table_de_WW(i,1:3)

    end do

    do i = int(density), cutoff_int
      table_ene_WW(6*i-2,1:3) = table_ene_WW(6*(i+1)-5,1:3) - &
                                table_ene_WW(6*i-5,1:3)
      table_ene_WW(6*i-1,1:3) = table_ene_WW(6*(i+1)-4,1:3) - &
                                table_ene_WW(6*i-4,1:3)
      table_ene_WW(6*i,1:3)   = table_ene_WW(6*(i+1)-3,1:3) - &
                                table_ene_WW(6*i-3,1:3)
    end do

    return

  end subroutine table_water_pme_linear_switch

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    table_water_pme_linear_fswitch
  !> @brief        table for pme and linear and force switch
  !! @authors      CK, JJ
  !! @param[in]    switch_int  : integer switchdist
  !! @param[in]    cutoff_int  : integer cutoff distance
  !! @param[in]    density     : table density
  !! @param[in]    cutoff2     : square cutoff distance
  !! @param[in]    switchdist2 : square switch distance
  !! @param[in]    el_fact     : inverse of dielectoric constant
  !! @param[in]    alpha       : ewald factor
  !! @param[in]    alpha2sp    : ewald factor
  !! @param[in]    alpha2m     : ewald factor
  !! @param[in]    lj6_OO      : LJ6 oxygen-oxygen
  !! @param[in]    lj6_OH      : LJ6 oxygen-hydrogen
  !! @param[in]    lj6_HH      : LJ6 hydrogen-hydrogen
  !! @param[in]    lj12_OO     : LJ12 oxygen-oxygen
  !! @param[in]    lj12_OH     : LJ12 oxygen-hydrogen
  !! @param[in]    lj12_HH     : LJ12 hydrogen-hydrogen
  !! @param[in]    cc_OO       : charge oxygen-oxygen
  !! @param[in]    cc_OH       : charge oxygen-hydrogen
  !! @param[in]    cc_HH       : charge hydrogen-hydrogen
  !! @param[out]   table_ene_WW: energy table
  !! @param[out]   table_de_WW  : gradient table
  !! @note         J.Jung et al., J.Comput.Chem., 34, 2412-2420 (2013)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8


  subroutine table_water_pme_linear_fswitch(switch_int, cutoff_int, density,   &
                                     cutoff2, switchdist2, el_fact, alpha,     &
                                     lj6_OO,  lj6_OH,  lj6_HH,                 &
                                     lj12_OO, lj12_OH, lj12_HH,                &
                                     cc_OO, cc_OH, cc_HH,                      &
                                     alpha2sp, alpha2m,                        &
                                     table_ene_WW, table_de_WW)
    ! formal arguments
    integer,                 intent(in)    :: switch_int
    integer,                 intent(in)    :: cutoff_int
    real(wp),                intent(in)    :: density
    real(wp),                intent(in)    :: cutoff2
    real(wp),                intent(in)    :: switchdist2
    real(wp),                intent(in)    :: el_fact
    real(wp),                intent(in)    :: alpha
    real(wp),                intent(in)    :: alpha2m
    real(wp),                intent(in)    :: alpha2sp
    real(wp),                intent(in)    :: lj6_OO
    real(wp),                intent(in)    :: lj6_OH
    real(wp),                intent(in)    :: lj6_HH
    real(wp),                intent(in)    :: lj12_OO
    real(wp),                intent(in)    :: lj12_OH
    real(wp),                intent(in)    :: lj12_HH
    real(wp),                intent(in)    :: cc_OO
    real(wp),                intent(in)    :: cc_OH
    real(wp),                intent(in)    :: cc_HH
    real(wp),                intent(inout) :: table_ene_WW(:,:)
    real(wp),                intent(inout) :: table_de_WW(:,:)

    ! local variables
    integer                  :: i
    real(wp)                 :: rij, rij2
    real(wp)                 :: inv_r2, inv_rij, inv_r3, inv_r6, inv_r12
    real(wp)                 :: lj1, lj2, lj4, lj6, lj12
    real(wp)                 :: term_lj12_OO, term_lj12_OH, term_lj12_HH
    real(wp)                 :: term_lj6_OO, term_lj6_OH, term_lj6_HH
    real(wp)                 :: term_lj_OO, term_lj_OH, term_lj_HH
    real(wp)                 :: term_lj12, term_lj6
    real(wp)                 :: switch, dswitch

    lj1 = 1.0_wp / (cutoff2*switchdist2)
    lj1 = lj1*lj1*lj1
    lj2 = sqrt(lj1)

    do i = switch_int+1, cutoff_int

      rij2     = cutoff2*density/real(i,wp)
      rij    = sqrt(rij2)

      inv_r2  = 1.0_wp / rij2
      inv_rij = 1.0_wp / rij
      inv_r6  = inv_r2 * inv_r2 * inv_r2
      inv_r12 = inv_r6 * inv_r6

      term_lj12_OO = lj12_OO * (inv_r12-lj1)
      term_lj12_OH = lj12_OH * (inv_r12-lj1)
      term_lj12_HH = lj12_HH * (inv_r12-lj1)
      term_lj6_OO  = lj6_OO * (inv_r6-lj2)
      term_lj6_OH  = lj6_OH * (inv_r6-lj2)
      term_lj6_HH  = lj6_HH * (inv_r6-lj2)

      term_lj_OO = term_lj12_OO - term_lj6_OO
      term_lj_OH = term_lj12_OH - term_lj6_OH
      term_lj_HH = term_lj12_HH - term_lj6_HH

      table_ene_WW(6*i-5,1) = term_lj_OO
      table_ene_WW(6*i-5,2) = term_lj_OH
      table_ene_WW(6*i-5,3) = term_lj_HH

      term_lj12_OO = lj12_OO * inv_r12
      term_lj12_OH = lj12_OH * inv_r12
      term_lj12_HH = lj12_HH * inv_r12
      term_lj6_OO  = lj6_OO * inv_r6
      term_lj6_OH  = lj6_OH * inv_r6
      term_lj6_HH  = lj6_HH * inv_r6

      table_de_WW(i,1)   = -(12.0_wp*term_lj12_OO-6.0_wp*term_lj6_OO)*inv_r2
      table_de_WW(i,2)   = -(12.0_wp*term_lj12_OH-6.0_wp*term_lj6_OH)*inv_r2
      table_de_WW(i,3)   = -(12.0_wp*term_lj12_HH-6.0_wp*term_lj6_HH)*inv_r2

      table_ene_WW(6*i-4,1) = cc_OO*erfc04(alpha*rij)*inv_rij
      table_ene_WW(6*i-4,2) = cc_OH*erfc04(alpha*rij)*inv_rij
      table_ene_WW(6*i-4,3) = cc_HH*erfc04(alpha*rij)*inv_rij

      table_de_WW(i,1)   = table_de_WW(i,1)                      &
                         - (table_ene_WW(6*i-4,1)                &
                         + alpha2sp*cc_OO*exp(alpha2m*rij2))*inv_r2
      table_de_WW(i,2)   = table_de_WW(i,2)                      &
                         - (table_ene_WW(6*i-4,2)                &
                         + alpha2sp*cc_OH*exp(alpha2m*rij2))*inv_r2
      table_de_WW(i,3)   = table_de_WW(i,3)                      &
                         - (table_ene_WW(6*i-4,3)                &
                         +alpha2sp*cc_HH*exp(alpha2m*rij2))*inv_r2

      table_ene_WW(6*i-3,1:3) = table_de_WW(i,1:3)

    end do

    lj1  = cutoff2*cutoff2*cutoff2
    lj2  = cutoff2*sqrt(cutoff2)
    lj12 = lj1 / (lj1-switchdist2*switchdist2*switchdist2)
    lj6  = lj2 / (lj2-switchdist2*sqrt(switchdist2))
    lj1  = 1.0_wp / lj1
    lj2  = 1.0_wp / lj2

    do i = int(density), switch_int
      rij2      = cutoff2*density/real(i,wp)
      rij       = sqrt(rij2)
      inv_r2    = 1.0_wp / rij2
      inv_rij   = 1.0_wp / rij
      inv_r6    = inv_r2 * inv_r2 * inv_r2
      inv_r3    = sqrt(inv_r6)
      term_lj12 = lj12*(inv_r6-lj1)**2
      term_lj6  = lj6*(inv_r3-lj2)**2

      term_lj12_OO = lj12_OO * term_lj12
      term_lj12_OH = lj12_OH * term_lj12
      term_lj12_HH = lj12_HH * term_lj12
      term_lj6_OO   = lj6_OO * term_lj6
      term_lj6_OH   = lj6_OH * term_lj6
      term_lj6_HH   = lj6_HH * term_lj6

      term_lj_OO = term_lj12_OO - term_lj6_OO
      term_lj_OH = term_lj12_OH - term_lj6_OH
      term_lj_HH = term_lj12_HH - term_lj6_HH

      table_ene_WW(6*i-5,1) = term_lj_OO
      table_ene_WW(6*i-5,2) = term_lj_OH
      table_ene_WW(6*i-5,3) = term_lj_HH

      term_lj12 = -12.0_wp*lj12*(inv_r6-lj1)*inv_r6*inv_r2
      term_lj6  = -6.0_wp*lj6*(inv_r3-lj2)*inv_r3*inv_r2

      term_lj12_OO = lj12_OO * term_lj12
      term_lj12_OH = lj12_OH * term_lj12
      term_lj12_HH = lj12_HH * term_lj12
      term_lj6_OO  = lj6_OO * term_lj6
      term_lj6_OH  = lj6_OH * term_lj6
      term_lj6_HH  = lj6_HH * term_lj6

      term_lj_OO = term_lj12_OO - term_lj6_OO
      term_lj_OH = term_lj12_OH - term_lj6_OH
      term_lj_HH = term_lj12_HH - term_lj6_HH

      table_de_WW(i,1)   = term_lj_OO
      table_de_WW(i,2)   = term_lj_OH
      table_de_WW(i,3)   = term_lj_HH

      table_ene_WW(6*i-4,1) = cc_OO*erfc04(alpha*rij)*inv_rij
      table_ene_WW(6*i-4,2) = cc_OH*erfc04(alpha*rij)*inv_rij
      table_ene_WW(6*i-4,3) = cc_HH*erfc04(alpha*rij)*inv_rij

      table_de_WW(i,1) = table_de_WW(i,1) &
          -(table_ene_WW(6*i-4,1)+alpha2sp*cc_OO*exp(alpha2m*rij2))*inv_r2
      table_de_WW(i,2) = table_de_WW(i,2) &
          -(table_ene_WW(6*i-4,2)+alpha2sp*cc_OH*exp(alpha2m*rij2))*inv_r2
      table_de_WW(i,3) = table_de_WW(i,3) &
          -(table_ene_WW(6*i-4,3)+alpha2sp*cc_HH*exp(alpha2m*rij2))*inv_r2

      table_ene_WW(6*i-3,1:3) = table_de_WW(i,1:3)

    end do

    do i = int(density), cutoff_int
      table_ene_WW(6*i-2,1:3) = table_ene_WW(6*(i+1)-5,1:3) - &
                                table_ene_WW(6*i-5,1:3)
      table_ene_WW(6*i-1,1:3) = table_ene_WW(6*(i+1)-4,1:3) - &
                                table_ene_WW(6*i-4,1:3)
      table_ene_WW(6*i,1:3)   = table_ene_WW(6*(i+1)-3,1:3) - &
                                table_ene_WW(6*i-3,1:3)
    end do

    return

  end subroutine table_water_pme_linear_fswitch

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    table_water_pme_linear_groshift 
  !> @brief        table for pme and linear and gromacs shift water
  !! @authors      CK, JJ
  !! @param[in]    switch_int  : integer switchdist
  !! @param[in]    cutoff_int  : integer cutoff distance
  !! @param[in]    density     : table density
  !! @param[in]    cutoff      : cutoff distance
  !! @param[in]    switchdist  : switch distance
  !! @param[in]    el_fact     : inverse of dielectoric constant
  !! @param[in]    alpha       : ewald factor
  !! @param[in]    alpha2sp    : ewald factor
  !! @param[in]    alpha2m     : ewald factor
  !! @param[in]    lj6_OO      : LJ6 oxygen-oxygen
  !! @param[in]    lj6_OH      : LJ6 oxygen-hydrogen
  !! @param[in]    lj6_HH      : LJ6 hydrogen-hydrogen
  !! @param[in]    lj12_OO     : LJ12 oxygen-oxygen
  !! @param[in]    lj12_OH     : LJ12 oxygen-hydrogen
  !! @param[in]    lj12_HH     : LJ12 hydrogen-hydrogen
  !! @param[in]    cc_OO       : charge oxygen-oxygen
  !! @param[in]    cc_OH       : charge oxygen-hydrogen
  !! @param[in]    cc_HH       : charge hydrogen-hydrogen
  !! @param[out]   table_ene_WW: energy table
  !! @param[out]   table_de_WW  : gradient table
  !! @note         J.Jung et al., J.Comput.Chem., 34, 2412-2420 (2013)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine table_water_pme_linear_groshift(switch_int, cutoff_int, density, &
                                     cutoff, switchdist, el_fact, alpha,     &
                                     lj6_OO,  lj6_OH,  lj6_HH,                 &
                                     lj12_OO, lj12_OH, lj12_HH,                &
                                     cc_OO, cc_OH, cc_HH,                      &
                                     alpha2sp, alpha2m,                        &
                                     table_ene_WW, table_de_WW)
    ! formal arguments
    integer,                 intent(in)    :: switch_int
    integer,                 intent(in)    :: cutoff_int
    real(wp),                intent(in)    :: density
    real(wp),                intent(in)    :: cutoff
    real(wp),                intent(in)    :: switchdist
    real(wp),                intent(in)    :: el_fact
    real(wp),                intent(in)    :: alpha
    real(wp),                intent(in)    :: alpha2m
    real(wp),                intent(in)    :: alpha2sp
    real(wp),                intent(in)    :: lj6_OO
    real(wp),                intent(in)    :: lj6_OH
    real(wp),                intent(in)    :: lj6_HH
    real(wp),                intent(in)    :: lj12_OO
    real(wp),                intent(in)    :: lj12_OH
    real(wp),                intent(in)    :: lj12_HH
    real(wp),                intent(in)    :: cc_OO
    real(wp),                intent(in)    :: cc_OH
    real(wp),                intent(in)    :: cc_HH
    real(wp),                intent(inout) :: table_ene_WW(:,:)
    real(wp),                intent(inout) :: table_de_WW(:,:)

    ! local variables
    integer                  :: i
    real(wp)                 :: rij, rij2
    real(wp)                 :: inv_r2, inv_rij, inv_r6, inv_r12
    real(wp)                 :: lj1, lj2, lj4, cutoff2, switchdist2
    real(wp)                 :: switch, dswitch
    real(wp)                 :: Coef_A12, Coef_B12, Coef_C12
    real(wp)                 :: Coef_A06, Coef_B06, Coef_C06
    real(wp)                 :: term_lj12_OO, term_lj12_OH, term_lj12_HH
    real(wp)                 :: term_lj6_OO, term_lj6_OH, term_lj6_HH
    real(wp)                 :: term_lj_OO, term_lj_OH, term_lj_HH
    real(wp)                 :: F12, F06, P12, P06

    cutoff2 = cutoff*cutoff
    switchdist2 = switchdist*switchdist

    lj1 = 1.0_wp/(cutoff2-switchdist2)
    lj1 = lj1*lj1*lj1

    Coef_A12  = 16.0_wp*cutoff - 13.0_wp*switchdist
    Coef_A12  = - Coef_A12 / cutoff**14 / (cutoff-switchdist)**2
    Coef_B12  = 15.0_wp*cutoff - 13.0_wp*switchdist
    Coef_B12  = Coef_B12 / cutoff**14 / (cutoff-switchdist)**3
    Coef_A06  = 10.0_wp*cutoff - 7.0_wp*switchdist
    Coef_A06  = - Coef_A06 / cutoff**8 / (cutoff-switchdist)**2
    Coef_B06  = 9.0_wp*cutoff - 7.0_wp*switchdist
    Coef_B06  = Coef_B06 / cutoff**8 / (cutoff-switchdist)**3
    Coef_C12  = 1.0_wp/cutoff**12
    Coef_C12  = Coef_C12 - 12.0_wp*Coef_A12/3.0_wp*(cutoff-switchdist)**3
    Coef_C12  = Coef_C12 - 12.0_wp*Coef_B12/4.0_wp*(cutoff-switchdist)**4
    Coef_C06  = 1.0_wp/cutoff**6
    Coef_C06  = Coef_C06 - 6.0_wp*Coef_A06/3.0_wp*(cutoff-switchdist)**3
    Coef_C06  = Coef_C06 - 6.0_wp*Coef_B06/4.0_wp*(cutoff-switchdist)**4

    if (switchdist .le. EPS) then

      do i = int(density), cutoff_int
        rij2     = cutoff2*density/real(i,wp)
        rij      = sqrt(rij2)
        F12       = Coef_A12*(rij-switchdist)**2 + Coef_B12*(rij-switchdist)**3
        F06       = Coef_A06*(rij-switchdist)**2 + Coef_B06*(rij-switchdist)**3
        P12       = - 12.0_wp*Coef_A12/3.0_wp*(rij-switchdist)**3
        P12       = P12 - 12.0_wp*Coef_B12/4.0_wp*(rij-switchdist)**4 - Coef_C12
        P06       = - 6.0_wp*Coef_A06/3.0_wp*(rij-switchdist)**3
        P06       = P06 - 6.0_wp*Coef_B06/4.0_wp*(rij-switchdist)**4 - Coef_C06
        inv_r2    = 1.0_wp/rij2
        inv_rij   = 1.0_wp/rij
        inv_r6    = inv_r2 * inv_r2 * inv_r2
        inv_r12   = inv_r6 * inv_r6

        term_lj12_OO  = lj12_OO * (inv_r12+P12)
        term_lj12_OH  = lj12_OH * (inv_r12+P12)
        term_lj12_HH  = lj12_HH * (inv_r12+P12)
        term_lj6_OO   = lj6_OO * (inv_r6+P06)
        term_lj6_OH   = lj6_OH * (inv_r6+P06)
        term_lj6_HH   = lj6_HH * (inv_r6+P06)

        table_ene_WW(6*i-5,1) = term_lj12_OO - term_lj6_OO
        table_ene_WW(6*i-5,2) = term_lj12_OH - term_lj6_OH
        table_ene_WW(6*i-5,3) = term_lj12_HH - term_lj6_HH

        term_lj12_OO  = -12.0_wp * lj12_OO * (inv_r12*inv_r2+F12*inv_rij)
        term_lj12_OH  = -12.0_wp * lj12_OH * (inv_r12*inv_r2+F12*inv_rij)
        term_lj12_HH  = -12.0_wp * lj12_HH * (inv_r12*inv_r2+F12*inv_rij)
        term_lj6_OO   = -6.0_wp  * lj6_OO  * (inv_r6*inv_r2+F06*inv_rij)
        term_lj6_OH   = -6.0_wp  * lj6_OH  * (inv_r6*inv_r2+F06*inv_rij)
        term_lj6_HH   = -6.0_wp  * lj6_HH  * (inv_r6*inv_r2+F06*inv_rij)

        table_de_WW(i,1)   = term_lj12_OO - term_lj6_OO
        table_de_WW(i,2)   = term_lj12_OH - term_lj6_OH
        table_de_WW(i,3)   = term_lj12_HH - term_lj6_HH
  
        table_ene_WW(6*i-4,1) = cc_OO*erfc04(alpha*rij)*inv_rij
        table_ene_WW(6*i-4,2) = cc_OH*erfc04(alpha*rij)*inv_rij
        table_ene_WW(6*i-4,3) = cc_HH*erfc04(alpha*rij)*inv_rij
  
        table_de_WW(i,1)   = table_de_WW(i,1) &
                -(table_ene_WW(6*i-4,1)+alpha2sp*cc_OO*exp(alpha2m*rij2))*inv_r2
        table_de_WW(i,2)   = table_de_WW(i,2) &
                -(table_ene_WW(6*i-4,2)+alpha2sp*cc_OH*exp(alpha2m*rij2))*inv_r2
        table_de_WW(i,3)   = table_de_WW(i,3) &
                -(table_ene_WW(6*i-4,3)+alpha2sp*cc_HH*exp(alpha2m*rij2))*inv_r2
  
        table_ene_WW(6*i-3,1:3) = table_de_WW(i,1:3)

      end do

    else

      do i = int(density), switch_int
        rij2     = cutoff2*density/real(i,wp)
        rij      = sqrt(rij2)
        F12       = Coef_A12*(rij-switchdist)**2 + Coef_B12*(rij-switchdist)**3
        F06       = Coef_A06*(rij-switchdist)**2 + Coef_B06*(rij-switchdist)**3
        P12       = - 12.0_wp*Coef_A12/3.0_wp*(rij-switchdist)**3
        P12       = P12 - 12.0_wp*Coef_B12/4.0_wp*(rij-switchdist)**4 - Coef_C12
        P06       = - 6.0_wp*Coef_A06/3.0_wp*(rij-switchdist)**3
        P06       = P06 - 6.0_wp*Coef_B06/4.0_wp*(rij-switchdist)**4 - Coef_C06
        inv_r2    = 1.0_wp/rij2
        inv_rij   = 1.0_wp/rij
        inv_r6    = inv_r2 * inv_r2 * inv_r2
        inv_r12   = inv_r6 * inv_r6

        term_lj12_OO  = lj12_OO * (inv_r12+P12)
        term_lj12_OH  = lj12_OH * (inv_r12+P12)
        term_lj12_HH  = lj12_HH * (inv_r12+P12)
        term_lj6_OO   = lj6_OO * (inv_r6+P06)
        term_lj6_OH   = lj6_OH * (inv_r6+P06)
        term_lj6_HH   = lj6_HH * (inv_r6+P06)

        table_ene_WW(6*i-5,1) = term_lj12_OO - term_lj6_OO
        table_ene_WW(6*i-5,2) = term_lj12_OH - term_lj6_OH
        table_ene_WW(6*i-5,3) = term_lj12_HH - term_lj6_HH

        term_lj12_OO  = -12.0_wp * lj12_OO * (inv_r12*inv_r2+F12*inv_rij)
        term_lj12_OH  = -12.0_wp * lj12_OH * (inv_r12*inv_r2+F12*inv_rij)
        term_lj12_HH  = -12.0_wp * lj12_HH * (inv_r12*inv_r2+F12*inv_rij)
        term_lj6_OO   = -6.0_wp  * lj6_OO  * (inv_r6*inv_r2+F06*inv_rij)
        term_lj6_OH   = -6.0_wp  * lj6_OH  * (inv_r6*inv_r2+F06*inv_rij)
        term_lj6_HH   = -6.0_wp  * lj6_HH  * (inv_r6*inv_r2+F06*inv_rij)

        table_de_WW(i,1)   = term_lj12_OO - term_lj6_OO
        table_de_WW(i,2)   = term_lj12_OH - term_lj6_OH
        table_de_WW(i,3)   = term_lj12_HH - term_lj6_HH

        table_ene_WW(6*i-4,1) = cc_OO*erfc04(alpha*rij)*inv_rij
        table_ene_WW(6*i-4,2) = cc_OH*erfc04(alpha*rij)*inv_rij
        table_ene_WW(6*i-4,3) = cc_HH*erfc04(alpha*rij)*inv_rij

        table_de_WW(i,1)   = table_de_WW(i,1) &
                -(table_ene_WW(6*i-4,1)+alpha2sp*cc_OO*exp(alpha2m*rij2))*inv_r2
        table_de_WW(i,2)   = table_de_WW(i,2) &
                -(table_ene_WW(6*i-4,2)+alpha2sp*cc_OH*exp(alpha2m*rij2))*inv_r2
        table_de_WW(i,3)   = table_de_WW(i,3) &
                -(table_ene_WW(6*i-4,3)+alpha2sp*cc_HH*exp(alpha2m*rij2))*inv_r2

        table_ene_WW(6*i-3,1:3) = table_de_WW(i,1:3)

      end do

      do i = switch_int+1, cutoff_int

        rij2     = cutoff2*density/real(i,wp)
        rij    = sqrt(rij2)

        inv_r2  = 1.0_wp / rij2
        inv_rij = 1.0_wp / rij
        inv_r6  = inv_r2 * inv_r2 * inv_r2
        inv_r12 = inv_r6 * inv_r6
  
        term_lj12_OO = lj12_OO * (inv_r12-Coef_C12)
        term_lj12_OH = lj12_OH * (inv_r12-Coef_C12)
        term_lj12_HH = lj12_HH * (inv_r12-Coef_C12)
        term_lj6_OO  = lj6_OO * (inv_r6-Coef_C06)
        term_lj6_OH  = lj6_OH * (inv_r6-Coef_C06)
        term_lj6_HH  = lj6_HH * (inv_r6-Coef_C06)

        term_lj_OO = term_lj12_OO - term_lj6_OO
        term_lj_OH = term_lj12_OH - term_lj6_OH
        term_lj_HH = term_lj12_HH - term_lj6_HH
  
        table_ene_WW(6*i-5,1) = term_lj_OO
        table_ene_WW(6*i-5,2) = term_lj_OH
        table_ene_WW(6*i-5,3) = term_lj_HH
  
        table_de_WW(i,1)   = -(12.0_wp*term_lj12_OO-6.0_wp*term_lj6_OO)*inv_r2
        table_de_WW(i,2)   = -(12.0_wp*term_lj12_OH-6.0_wp*term_lj6_OH)*inv_r2
        table_de_WW(i,3)   = -(12.0_wp*term_lj12_HH-6.0_wp*term_lj6_HH)*inv_r2
  
        table_ene_WW(6*i-4,1) = cc_OO*erfc04(alpha*rij)*inv_rij
        table_ene_WW(6*i-4,2) = cc_OH*erfc04(alpha*rij)*inv_rij
        table_ene_WW(6*i-4,3) = cc_HH*erfc04(alpha*rij)*inv_rij
  
        table_de_WW(i,1)   = table_de_WW(i,1)       &
                           - (table_ene_WW(6*i-4,1) &
                           + alpha2sp*cc_OO*exp(alpha2m*rij2))*inv_r2
        table_de_WW(i,2)   = table_de_WW(i,2)       &
                           - (table_ene_WW(6*i-4,2) &
                           + alpha2sp*cc_OH*exp(alpha2m*rij2))*inv_r2
        table_de_WW(i,3)   = table_de_WW(i,3)       &
                           - (table_ene_WW(6*i-4,3) &
                           +alpha2sp*cc_HH*exp(alpha2m*rij2))*inv_r2
  
        table_ene_WW(6*i-3,1:3) = table_de_WW(i,1:3)
  
      end do

    end if


    do i = int(density), cutoff_int
      table_ene_WW(6*i-2,1:3) = table_ene_WW(6*(i+1)-5,1:3) - &
                                table_ene_WW(6*i-5,1:3)
      table_ene_WW(6*i-1,1:3) = table_ene_WW(6*(i+1)-4,1:3) - &
                                table_ene_WW(6*i-4,1:3)
      table_ene_WW(6*i,1:3)   = table_ene_WW(6*(i+1)-3,1:3) - &
                                table_ene_WW(6*i-3,1:3)
    end do

    return

  end subroutine table_water_pme_linear_groshift

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    table_water_pme_linear_groswitch
  !> @brief        table for pme and linear and switch
  !! @authors      CK, JJ
  !! @param[in]    switch_int  : integer switchdist
  !! @param[in]    cutoff_int  : integer cutoff distance
  !! @param[in]    density     : table density
  !! @param[in]    cutoff2     : square cutoff distance
  !! @param[in]    switchdist2 : square switch distance
  !! @param[in]    el_fact     : inverse of dielectoric constant
  !! @param[in]    alpha       : ewald factor
  !! @param[in]    alpha2sp    : ewald factor
  !! @param[in]    alpha2m     : ewald factor
  !! @param[in]    lj6_OO      : LJ6 oxygen-oxygen
  !! @param[in]    lj6_OH      : LJ6 oxygen-hydrogen
  !! @param[in]    lj6_HH      : LJ6 hydrogen-hydrogen
  !! @param[in]    lj12_OO     : LJ12 oxygen-oxygen
  !! @param[in]    lj12_OH     : LJ12 oxygen-hydrogen
  !! @param[in]    lj12_HH     : LJ12 hydrogen-hydrogen
  !! @param[in]    cc_OO       : charge oxygen-oxygen
  !! @param[in]    cc_OH       : charge oxygen-hydrogen
  !! @param[in]    cc_HH       : charge hydrogen-hydrogen
  !! @param[out]   table_ene_WW: energy table
  !! @param[out]   table_de_WW  : gradient table
  !! @note         J.Jung et al., J.Comput.Chem., 34, 2412-2420 (2013)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine table_water_pme_linear_groswitch(switch_int, cutoff_int, density, &
                                     cutoff, switchdist, el_fact, alpha,       &
                                     lj6_OO,  lj6_OH,  lj6_HH,                 &
                                     lj12_OO, lj12_OH, lj12_HH,                &
                                     cc_OO, cc_OH, cc_HH,                      &
                                     alpha2sp, alpha2m,                        &
                                     table_ene_WW, table_de_WW)
    ! formal arguments
    integer,                 intent(in)    :: switch_int
    integer,                 intent(in)    :: cutoff_int
    real(wp),                intent(in)    :: density
    real(wp),                intent(in)    :: cutoff
    real(wp),                intent(in)    :: switchdist
    real(wp),                intent(in)    :: el_fact
    real(wp),                intent(in)    :: alpha
    real(wp),                intent(in)    :: alpha2m
    real(wp),                intent(in)    :: alpha2sp
    real(wp),                intent(in)    :: lj6_OO
    real(wp),                intent(in)    :: lj6_OH
    real(wp),                intent(in)    :: lj6_HH
    real(wp),                intent(in)    :: lj12_OO
    real(wp),                intent(in)    :: lj12_OH
    real(wp),                intent(in)    :: lj12_HH
    real(wp),                intent(in)    :: cc_OO
    real(wp),                intent(in)    :: cc_OH
    real(wp),                intent(in)    :: cc_HH
    real(wp),                intent(inout) :: table_ene_WW(:,:)
    real(wp),                intent(inout) :: table_de_WW(:,:)

    ! local variables
    integer                  :: i
    real(wp)                 :: rij, rij2, rij1
    real(wp)                 :: inv_r2, inv_rij, inv_r6, inv_r12
    real(wp)                 :: lj1, lj2, lj4
    real(wp)                 :: term_lj12_OO, term_lj12_OH, term_lj12_HH
    real(wp)                 :: term_lj6_OO, term_lj6_OH, term_lj6_HH
    real(wp)                 :: term_lj_OO, term_lj_OH, term_lj_HH
    real(wp)                 :: cutoff2
    real(wp)                 :: switch, dswitch
    real(wp)                 :: coef_A, coef_B, coef_C
    real(wp)                 :: coef_A1, coef_B1, coef_C1

    cutoff2 = cutoff*cutoff

    do i = switch_int+1, cutoff_int

      rij2     = cutoff2*density/real(i,wp)
      rij    = sqrt(rij2)
      switch  = 1.0_wp
      dswitch = 0.0_wp

      inv_r2  = 1.0_wp / rij2
      inv_rij = 1.0_wp / rij
      inv_r6  = inv_r2 * inv_r2 * inv_r2
      inv_r12 = inv_r6 * inv_r6

      term_lj12_OO = lj12_OO * inv_r12
      term_lj12_OH = lj12_OH * inv_r12
      term_lj12_HH = lj12_HH * inv_r12
      term_lj6_OO  = lj6_OO * inv_r6
      term_lj6_OH  = lj6_OH * inv_r6
      term_lj6_HH  = lj6_HH * inv_r6

      term_lj_OO = term_lj12_OO - term_lj6_OO
      term_lj_OH = term_lj12_OH - term_lj6_OH
      term_lj_HH = term_lj12_HH - term_lj6_HH

      table_ene_WW(6*i-5,1) = term_lj_OO
      table_ene_WW(6*i-5,2) = term_lj_OH
      table_ene_WW(6*i-5,3) = term_lj_HH

      table_de_WW(i,1)   = - (12.0_wp*term_lj12_OO-6.0_wp*term_lj6_OO)*inv_r2
      table_de_WW(i,2)   = - (12.0_wp*term_lj12_OH-6.0_wp*term_lj6_OH)*inv_r2
      table_de_WW(i,3)   = - (12.0_wp*term_lj12_HH-6.0_wp*term_lj6_HH)*inv_r2

      table_ene_WW(6*i-4,1) = cc_OO*erfc04(alpha*rij)*inv_rij
      table_ene_WW(6*i-4,2) = cc_OH*erfc04(alpha*rij)*inv_rij
      table_ene_WW(6*i-4,3) = cc_HH*erfc04(alpha*rij)*inv_rij

      table_de_WW(i,1)   = table_de_WW(i,1)       &
                         - (table_ene_WW(6*i-4,1) &
                         + alpha2sp*cc_OO*exp(alpha2m*rij2))*inv_r2
      table_de_WW(i,2)   = table_de_WW(i,2)       &
                         - (table_ene_WW(6*i-4,2) &
                         + alpha2sp*cc_OH*exp(alpha2m*rij2))*inv_r2
      table_de_WW(i,3)   = table_de_WW(i,3)       &
                         - (table_ene_WW(6*i-4,3) &
                         +alpha2sp*cc_HH*exp(alpha2m*rij2))*inv_r2

      table_ene_WW(6*i-3,1:3) = table_de_WW(i,1:3)

    end do

    coef_A    = -10.0_wp / (cutoff-switchdist)**3
    coef_B    =  15.0_wp / (cutoff-switchdist)**4
    coef_C    = - 6.0_wp / (cutoff-switchdist)**5
    coef_A1   = -30.0_wp / (cutoff-switchdist)**3
    coef_B1   =  60.0_wp / (cutoff-switchdist)**4
    coef_C1   = -30.0_wp / (cutoff-switchdist)**5

    do i = int(density), switch_int
      rij2     = cutoff2*density/real(i,wp)
      rij      = sqrt(rij2)
      rij1      = rij - switchdist
      switch   = 1.0_wp + coef_A*rij1**3 + coef_B*rij1**4 + coef_C*rij1**5
      dswitch  = coef_A1*rij1**2 + coef_B*rij1**3 + coef_C*rij1**4
      inv_r2   = 1.0_wp/rij2
      inv_rij  = 1.0_wp / rij
      inv_r6   = inv_r2 * inv_r2 * inv_r2
      inv_r12  = inv_r6 * inv_r6

      term_lj12_OO  = lj12_OO * inv_r12
      term_lj12_OH  = lj12_OH * inv_r12
      term_lj12_HH  = lj12_HH * inv_r12
      term_lj6_OO   = lj6_OO * inv_r6
      term_lj6_OH   = lj6_OH * inv_r6
      term_lj6_HH   = lj6_HH * inv_r6

      term_lj_OO = term_lj12_OO - term_lj6_OO
      term_lj_OH = term_lj12_OH - term_lj6_OH
      term_lj_HH = term_lj12_HH - term_lj6_HH

      table_ene_WW(6*i-5,1) = switch*term_lj_OO
      table_ene_WW(6*i-5,2) = switch*term_lj_OH
      table_ene_WW(6*i-5,3) = switch*term_lj_HH

      table_de_WW(i,1)   = (  dswitch * term_lj_OO * rij     &
              - switch*(12.0_wp*term_lj12_OO-6.0_wp*term_lj6_OO))*inv_r2
      table_de_WW(i,2)   = (  dswitch * term_lj_OH * rij     &
              - switch * (12.0_wp*term_lj12_OH-6.0_wp*term_lj6_OH))*inv_r2
      table_de_WW(i,3)   = (  dswitch * term_lj_HH * rij     &
              - switch * (12.0_wp*term_lj12_HH-6.0_wp*term_lj6_HH))*inv_r2

      table_ene_WW(6*i-4,1) = cc_OO*erfc04(alpha*rij)*inv_rij
      table_ene_WW(6*i-4,2) = cc_OH*erfc04(alpha*rij)*inv_rij
      table_ene_WW(6*i-4,3) = cc_HH*erfc04(alpha*rij)*inv_rij

      table_de_WW(i,1)   = table_de_WW(i,1) &
              -(table_ene_WW(6*i-4,1)+alpha2sp*cc_OO*exp(alpha2m*rij2))*inv_r2
      table_de_WW(i,2)   = table_de_WW(i,2) &
              -(table_ene_WW(6*i-4,2)+alpha2sp*cc_OH*exp(alpha2m*rij2))*inv_r2
      table_de_WW(i,3)   = table_de_WW(i,3) &
              -(table_ene_WW(6*i-4,3)+alpha2sp*cc_HH*exp(alpha2m*rij2))*inv_r2

      table_ene_WW(6*i-3,1:3) = table_de_WW(i,1:3)

    end do

    do i = int(density), cutoff_int
      table_ene_WW(6*i-2,1:3) = table_ene_WW(6*(i+1)-5,1:3) - &
                                table_ene_WW(6*i-5,1:3)
      table_ene_WW(6*i-1,1:3) = table_ene_WW(6*(i+1)-4,1:3) - &
                                table_ene_WW(6*i-4,1:3)
      table_ene_WW(6*i,1:3)   = table_ene_WW(6*(i+1)-3,1:3) - &
                                table_ene_WW(6*i-3,1:3)
    end do

    return

  end subroutine table_water_pme_linear_groswitch

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    table_cutoff_cubic_noswitch
  !> @brief        table for pme and cubic  without switch
  !! @authors      CK, JJ
  !! @param[in]    switch_int  : integer switchdist
  !! @param[in]    cutoff_int  : integer cutoff distance
  !! @param[in]    cutoff_int2 : integer cutoff distance (buffering region)
  !! @param[in]    density     : table density
  !! @param[in]    cutoff2     : square cutoff distance
  !! @param[in]    switchdist2 : square switch distance
  !! @param[out]   table_ene   : energy table
  !! @param[out]   table_grad  : gradient table
  !! @note         J.Jung et al., J.Comput.Chem., 34, 2412-2420 (2013)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine table_cutoff_cubic_noswitch(switch_int, cutoff_int, cutoff_int2, &
                                  density, el_fact,                           &
                                  cutoff2, switchdist2,  table_ene, table_grad) 
    ! formal arguments
    integer,                 intent(in)    :: switch_int
    integer,                 intent(in)    :: cutoff_int
    integer,                 intent(in)    :: cutoff_int2
    real(wp),                intent(in)    :: density
    real(wp),                intent(in)    :: el_fact
    real(wp),                intent(in)    :: cutoff2
    real(wp),                intent(in)    :: switchdist2
    real(wp),                intent(inout) :: table_ene(:)
    real(wp),                intent(inout) :: table_grad(:)

    ! local variables
    integer                  :: i, icoun
    real(wp)                 :: rij, rij2
    real(wp)                 :: inv_r2, inv_rij, inv_r6, inv_r12
    real(wp)                 :: lj2
    real(wp)                 :: cutoff2_inv

    real(wp),    allocatable :: elec_array(:), delec_array(:)
    real(wp),    allocatable :: lj12_array(:), dlj12_array(:)
    real(wp),    allocatable :: lj6_array(:), dlj6_array(:)

    cutoff2_inv = 1.0_wp/cutoff2

    allocate(elec_array(cutoff_int2), delec_array(cutoff_int2))
    allocate(lj12_array(cutoff_int2), dlj12_array(cutoff_int2))
    allocate(lj6_array (cutoff_int2), dlj6_array (cutoff_int2))

    do i = 2, switch_int

      rij2     = real(i,wp)/density
      rij      = sqrt(rij2)
      inv_r2   = 1.0_wp / rij2
      inv_rij  = 1.0_wp / rij
      inv_r6   = inv_r2 * inv_r2 * inv_r2
      inv_r12  = inv_r6 * inv_r6

      table_ene(6*i-5)  = inv_r12
      table_ene(6*i-3)  = inv_r6
      table_grad(6*i-5) = -12.0_wp*inv_r12*inv_r2
      table_grad(6*i-3) = -6.0_wp*inv_r6*inv_r2
      table_ene(6*i-1)  = el_fact * inv_rij
      table_grad(6*i-1) = -el_fact * inv_r2 * inv_rij

      lj12_array(i-1)  = table_ene(6*i-5)  - table_ene(6*i-11)
      lj6_array(i-1)   = table_ene(6*i-3)  - table_ene(6*i-9)
      dlj12_array(i-1) = table_grad(6*i-5) - table_grad(6*i-11)
      dlj6_array(i-1)  = table_grad(6*i-3) - table_grad(6*i-9)
      elec_array(i-1)  = table_ene(6*i-1)  - table_ene(6*i-7)
      delec_array(i-1) = table_grad(6*i-1) - table_grad(6*i-7)
    end do

    icoun = max(2, switch_int+1)

    do i = icoun, cutoff_int
!    do i = switch_int+1, cutoff_int

      rij2     = real(i,wp)/density
      rij      = sqrt(rij2)
      lj2      = cutoff2 - rij2
      inv_r2   = 1.0_wp / rij2
      inv_rij  = 1.0_wp / rij
      inv_r6   = inv_r2 * inv_r2 * inv_r2
      inv_r12  = inv_r6 * inv_r6

      table_ene(6*i-5)  = inv_r12
      table_ene(6*i-3)  = inv_r6
      table_grad(6*i-5) = -12.0_wp*inv_r12*inv_r2
      table_grad(6*i-3) = -6.0_wp*inv_r6*inv_r2
      table_ene(6*i-1)  = el_fact * inv_rij
      table_grad(6*i-1) = -el_fact * inv_r2 * inv_rij

      lj12_array(i-1)  = table_ene(6*i-5)  - table_ene(6*i-11)
      lj6_array(i-1)   = table_ene(6*i-3)  - table_ene(6*i-9)
      dlj12_array(i-1) = table_grad(6*i-5) - table_grad(6*i-11)
      dlj6_array(i-1)  = table_grad(6*i-3) - table_grad(6*i-9)
      elec_array(i-1)  = table_ene(6*i-1)  - table_ene(6*i-7)
      delec_array(i-1) = table_grad(6*i-1) - table_grad(6*i-7)
    end do

    table_ene (6*cutoff_int-5:6*cutoff_int) = 0.0_wp
    table_grad(6*cutoff_int-5:6*cutoff_int) = 0.0_wp

    do i = 2, cutoff_int-1
      table_ene(6*i-4)  = (lj12_array(i-1) + lj12_array(i))   / 2.0_wp
      table_grad(6*i-4) = (dlj12_array(i-1) + dlj12_array(i)) / 2.0_wp
      table_ene(6*i-2)  = (lj6_array(i-1) + lj6_array(i))     / 2.0_wp
      table_grad(6*i-2) = (dlj6_array(i-1) + dlj6_array(i))   / 2.0_wp
      table_ene(6*i)    = (elec_array(i-1) + elec_array(i))   / 2.0_wp
      table_grad(6*i)   = (delec_array(i-1) + delec_array(i)) / 2.0_wp
    end do

    table_ene(2)               = lj12_array(1)
    table_ene(4)               = lj6_array(1)
    table_ene(6)               = elec_array(1)
    table_ene(6*cutoff_int-4)  = lj12_array(cutoff_int-1)
    table_ene(6*cutoff_int-2)  = lj6_array(cutoff_int-1)
    table_ene(6*cutoff_int)    = elec_array(cutoff_int-1)
    table_grad(2)              = dlj12_array(1)
    table_grad(4)              = dlj6_array(1)
    table_grad(6)              = delec_array(1)
    table_grad(6*cutoff_int-4) = dlj12_array(cutoff_int-1)
    table_grad(6*cutoff_int-2) = dlj6_array(cutoff_int-1)
    table_grad(6*cutoff_int)   = delec_array(cutoff_int-1)

    deallocate(elec_array, delec_array)
    deallocate(lj12_array, dlj12_array)
    deallocate(lj6_array, dlj6_array)

    return

  end subroutine table_cutoff_cubic_noswitch

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    table_cutoff_cubic_switch
  !> @brief        table for pme and cubic 
  !! @authors      CK, JJ
  !! @param[in]    switch_int  : integer switchdist
  !! @param[in]    cutoff_int  : integer cutoff distance
  !! @param[in]    cutoff_int2 : integer cutoff distance (buffering region)
  !! @param[in]    density     : table density
  !! @param[in]    cutoff2     : square cutoff distance
  !! @param[in]    switchdist2 : square switch distance
  !! @param[out]   table_ene   : energy table
  !! @param[out]   table_grad  : gradient table
  !! @note         J.Jung et al., J.Comput.Chem., 34, 2412-2420 (2013)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine table_cutoff_cubic_switch(switch_int, cutoff_int, cutoff_int2,   &
                                  density, el_fact,                           &
                                  cutoff2, switchdist2,  table_ene, table_grad) 
    ! formal arguments
    integer,                 intent(in)    :: switch_int
    integer,                 intent(in)    :: cutoff_int
    integer,                 intent(in)    :: cutoff_int2
    real(wp),                intent(in)    :: density
    real(wp),                intent(in)    :: el_fact
    real(wp),                intent(in)    :: cutoff2
    real(wp),                intent(in)    :: switchdist2
    real(wp),                intent(inout) :: table_ene(:)
    real(wp),                intent(inout) :: table_grad(:)

    ! local variables
    integer                  :: i, icoun
    real(wp)                 :: rij, rij2
    real(wp)                 :: inv_r2, inv_rij, inv_r6, inv_r12
    real(wp)                 :: lj1, lj2, lj4 
    real(wp)                 :: switch, dswitch, cutoff2_inv
    real(wp)                 :: shift,  dshift

    real(wp),    allocatable :: elec_array(:), delec_array(:)
    real(wp),    allocatable :: lj12_array(:), dlj12_array(:)
    real(wp),    allocatable :: lj6_array(:), dlj6_array(:)

    lj1 = 1.0_wp/(cutoff2-switchdist2)
    lj1 = lj1*lj1*lj1

    cutoff2_inv = 1.0_wp/cutoff2

    allocate(elec_array(cutoff_int2), delec_array(cutoff_int2))
    allocate(lj12_array(cutoff_int2), dlj12_array(cutoff_int2))
    allocate(lj6_array (cutoff_int2), dlj6_array (cutoff_int2))

    do i = 2, switch_int

      rij2     = real(i,wp)/density
      rij      = sqrt(rij2)
      switch   = 1.0_wp
      dswitch  = 0.0_wp
      shift    = (1.0_wp - (rij2*cutoff2_inv))
      dshift   = -4.0_wp*shift*rij*cutoff2_inv
      shift    = shift * shift
      inv_r2   = 1.0_wp / rij2
      inv_rij  = 1.0_wp / rij
      inv_r6   = inv_r2 * inv_r2 * inv_r2
      inv_r12  = inv_r6 * inv_r6

      table_ene(6*i-5)  = switch*inv_r12
      table_ene(6*i-3)  = switch*inv_r6
      table_grad(6*i-5) = dswitch*inv_r12*inv_rij          &
                         - switch*12.0_wp*inv_r12*inv_r2
      table_grad(6*i-3) = dswitch*inv_r6*inv_rij           &
                         - switch*6.0_wp*inv_r6*inv_r2
      table_ene(6*i-1)  = el_fact * shift * inv_rij
      table_grad(6*i-1) = el_fact * (dshift*inv_rij - shift*inv_r2) * inv_rij

      lj12_array(i-1)  = table_ene(6*i-5)  - table_ene(6*i-11)
      lj6_array(i-1)   = table_ene(6*i-3)  - table_ene(6*i-9)
      dlj12_array(i-1) = table_grad(6*i-5) - table_grad(6*i-11)
      dlj6_array(i-1)  = table_grad(6*i-3) - table_grad(6*i-9)
      elec_array(i-1)  = table_ene(6*i-1)  - table_ene(6*i-7)
      delec_array(i-1) = table_grad(6*i-1) - table_grad(6*i-7)
    end do

    icoun = max(2, switch_int+1)

    do i = icoun, cutoff_int
!    do i = switch_int+1, cutoff_int

      rij2     = real(i,wp)/density
      rij      = sqrt(rij2)
      lj2      = cutoff2 - rij2
      lj4      = lj2*(cutoff2 + 2.0_wp*rij2 - 3.0_wp*switchdist2)
      switch   = lj2*lj4*lj1
      dswitch  = 4.0_wp*lj1*rij*(lj2*lj2 - lj4)
      shift    = (1.0_wp - (rij2*cutoff2_inv))
      dshift   = -4.0_wp*shift*rij*cutoff2_inv
      shift    = shift * shift
      inv_r2   = 1.0_wp / rij2
      inv_rij  = 1.0_wp / rij
      inv_r6   = inv_r2 * inv_r2 * inv_r2
      inv_r12  = inv_r6 * inv_r6

      table_ene(6*i-5)  = switch*inv_r12
      table_ene(6*i-3)  = switch*inv_r6
      table_grad(6*i-5) = dswitch*inv_r12*inv_rij        &
                         - switch*12.0_wp*inv_r12*inv_r2
      table_grad(6*i-3) = dswitch*inv_r6*inv_rij         &
                         - switch*6.0_wp*inv_r6*inv_r2
      table_ene(6*i-1)  = el_fact * shift * inv_rij
      table_grad(6*i-1) = el_fact * (dshift*inv_rij - shift*inv_r2) * inv_rij

      lj12_array(i-1)  = table_ene(6*i-5)  - table_ene(6*i-11)
      lj6_array(i-1)   = table_ene(6*i-3)  - table_ene(6*i-9)
      dlj12_array(i-1) = table_grad(6*i-5) - table_grad(6*i-11)
      dlj6_array(i-1)  = table_grad(6*i-3) - table_grad(6*i-9)
      elec_array(i-1)  = table_ene(6*i-1)  - table_ene(6*i-7)
      delec_array(i-1) = table_grad(6*i-1) - table_grad(6*i-7)
    end do

    table_ene (6*cutoff_int-5:6*cutoff_int) = 0.0_wp
    table_grad(6*cutoff_int-5:6*cutoff_int) = 0.0_wp

    do i = 2, cutoff_int-1
      table_ene(6*i-4)  = (lj12_array(i-1) + lj12_array(i))   / 2.0_wp
      table_grad(6*i-4) = (dlj12_array(i-1) + dlj12_array(i)) / 2.0_wp
      table_ene(6*i-2)  = (lj6_array(i-1) + lj6_array(i))     / 2.0_wp
      table_grad(6*i-2) = (dlj6_array(i-1) + dlj6_array(i))   / 2.0_wp
      table_ene(6*i)    = (elec_array(i-1) + elec_array(i))   / 2.0_wp
      table_grad(6*i)   = (delec_array(i-1) + delec_array(i)) / 2.0_wp
    end do

    table_ene(2)               = lj12_array(1)
    table_ene(4)               = lj6_array(1)
    table_ene(6)               = elec_array(1)
    table_ene(6*cutoff_int-4)  = lj12_array(cutoff_int-1)
    table_ene(6*cutoff_int-2)  = lj6_array(cutoff_int-1)
    table_ene(6*cutoff_int)    = elec_array(cutoff_int-1)
    table_grad(2)              = dlj12_array(1)
    table_grad(4)              = dlj6_array(1)
    table_grad(6)              = delec_array(1)
    table_grad(6*cutoff_int-4) = dlj12_array(cutoff_int-1)
    table_grad(6*cutoff_int-2) = dlj6_array(cutoff_int-1)
    table_grad(6*cutoff_int)   = delec_array(cutoff_int-1)

    deallocate(elec_array, delec_array)
    deallocate(lj12_array, dlj12_array)
    deallocate(lj6_array, dlj6_array)

    return

  end subroutine table_cutoff_cubic_switch

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    table_cutoff_cubic_fswitch
  !> @brief        table for pme and cubic force switch
  !! @authors      CK, JJ
  !! @param[in]    switch_int  : integer switchdist
  !! @param[in]    cutoff_int  : integer cutoff distance
  !! @param[in]    cutoff_int2 : integer cutoff distance (buffering region)
  !! @param[in]    density     : table density
  !! @param[in]    cutoff2     : square cutoff distance
  !! @param[in]    switchdist2 : square switch distance
  !! @param[out]   table_ene   : energy table
  !! @param[out]   table_grad  : gradient table
  !! @note         J.Jung et al., J.Comput.Chem., 34, 2412-2420 (2013)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine table_cutoff_cubic_fswitch(switch_int, cutoff_int, cutoff_int2,   &
                                  density, el_fact,                           &
                                  cutoff2, switchdist2,  table_ene, table_grad) 
    ! formal arguments
    integer,                 intent(in)    :: switch_int
    integer,                 intent(in)    :: cutoff_int
    integer,                 intent(in)    :: cutoff_int2
    real(wp),                intent(in)    :: density
    real(wp),                intent(in)    :: el_fact
    real(wp),                intent(in)    :: cutoff2
    real(wp),                intent(in)    :: switchdist2
    real(wp),                intent(inout) :: table_ene(:)
    real(wp),                intent(inout) :: table_grad(:)

    ! local variables
    integer                  :: i, icoun
    real(wp)                 :: rij, rij2
    real(wp)                 :: inv_r2, inv_rij, inv_r3, inv_r6, inv_r12
    real(wp)                 :: lj1, lj2, lj4, lj6, lj12
    real(wp)                 :: switch, dswitch, cutoff2_inv
    real(wp)                 :: shift,  dshift

    real(wp),    allocatable :: elec_array(:), delec_array(:)
    real(wp),    allocatable :: lj12_array(:), dlj12_array(:)
    real(wp),    allocatable :: lj6_array(:), dlj6_array(:)

    lj1 = 1.0_wp / (cutoff2*switchdist2)
    lj1 = lj1*lj1*lj1
    lj2 = sqrt(lj1)

    cutoff2_inv = 1.0_wp/cutoff2

    allocate(elec_array(cutoff_int2), delec_array(cutoff_int2))
    allocate(lj12_array(cutoff_int2), dlj12_array(cutoff_int2))
    allocate(lj6_array (cutoff_int2), dlj6_array (cutoff_int2))

    do i = 2, switch_int
      rij2    = real(i,wp)/density
      rij     = sqrt(rij2)
      shift   = (1.0_wp - (rij2*cutoff2_inv))
      dshift  = -4.0_wp*shift*rij*cutoff2_inv
      shift   = shift * shift
      inv_r2  = 1.0_wp / rij2
      inv_rij = 1.0_wp / rij
      inv_r6  = inv_r2 * inv_r2 * inv_r2
      inv_r12 = inv_r6 * inv_r6

      table_ene(6*i-5)  = inv_r12 - lj1
      table_ene(6*i-3)  = inv_r6 - lj2
      table_grad(6*i-5) = -12.0_wp*inv_r12*inv_r2
      table_grad(6*i-3) = -6.0_wp*inv_r6*inv_r2
      table_ene(6*i-1)  = el_fact * shift * inv_rij
      table_grad(6*i-1) = el_fact * (dshift*inv_rij - shift*inv_r2) * inv_rij

      lj12_array(i-1)  = table_ene(6*i-5)  - table_ene(6*i-11)
      lj6_array(i-1)   = table_ene(6*i-3)  - table_ene(6*i-9)
      dlj12_array(i-1) = table_grad(6*i-5) - table_grad(6*i-11)
      dlj6_array(i-1)  = table_grad(6*i-3) - table_grad(6*i-9)
      elec_array(i-1)  = table_ene(6*i-1)  - table_ene(6*i-7)
      delec_array(i-1) = table_grad(6*i-1) - table_grad(6*i-7)
    end do

    lj1  = cutoff2*cutoff2*cutoff2
    lj2  = cutoff2*sqrt(cutoff2)
    lj12 = lj1 / (lj1-switchdist2*switchdist2*switchdist2)
    lj6  = lj2 / (lj2-switchdist2*sqrt(switchdist2))
    lj1  = 1.0_wp / lj1
    lj2  = 1.0_wp / lj2

    icoun = max(2, switch_int+1)

    do i = icoun, cutoff_int
!    do i = switch_int+1, cutoff_int
      rij2     = real(i,wp)/density
      rij      = sqrt(rij2)
      shift    = (1.0_wp - (rij2*cutoff2_inv))
      dshift   = -4.0_wp*shift*rij*cutoff2_inv
      shift    = shift * shift
      inv_r2   = 1.0_wp / rij2
      inv_rij  = 1.0_wp / rij
      inv_r6   = inv_r2 * inv_r2 * inv_r2
      inv_r3   = sqrt(inv_r6)

      table_ene(6*i-5)  = lj12*(inv_r6-lj1)**2
      table_ene(6*i-3)  = lj6*(inv_r3-lj2)**2
      table_grad(6*i-5) = -12.0_wp*lj12*(inv_r6-lj1)*inv_r6*inv_r2
      table_grad(6*i-3) = -6.0_wp*lj6*(inv_r3-lj2)*inv_r3*inv_r2
      table_ene(6*i-1)  = el_fact * shift * inv_rij
      table_grad(6*i-1) = el_fact * (dshift*inv_rij - shift*inv_r2) * inv_rij

      lj12_array(i-1)  = table_ene(6*i-5)  - table_ene(6*i-11)
      lj6_array(i-1)   = table_ene(6*i-3)  - table_ene(6*i-9)
      dlj12_array(i-1) = table_grad(6*i-5) - table_grad(6*i-11)
      dlj6_array(i-1)  = table_grad(6*i-3) - table_grad(6*i-9)
      elec_array(i-1)  = table_ene(6*i-1)  - table_ene(6*i-7)
      delec_array(i-1) = table_grad(6*i-1) - table_grad(6*i-7)
    end do

    table_ene (6*cutoff_int-5:6*cutoff_int) = 0.0_wp
    table_grad(6*cutoff_int-5:6*cutoff_int) = 0.0_wp

    do i = 2, cutoff_int-1
      table_ene(6*i-4)  = (lj12_array(i-1) + lj12_array(i))   / 2.0_wp
      table_grad(6*i-4) = (dlj12_array(i-1) + dlj12_array(i)) / 2.0_wp
      table_ene(6*i-2)  = (lj6_array(i-1) + lj6_array(i))     / 2.0_wp
      table_grad(6*i-2) = (dlj6_array(i-1) + dlj6_array(i))   / 2.0_wp
      table_ene(6*i)    = (elec_array(i-1) + elec_array(i))   / 2.0_wp
      table_grad(6*i)   = (delec_array(i-1) + delec_array(i)) / 2.0_wp
    end do

    table_ene(2)               = lj12_array(1)
    table_ene(4)               = lj6_array(1)
    table_ene(6)               = elec_array(1)
    table_ene(6*cutoff_int-4)  = lj12_array(cutoff_int-1)
    table_ene(6*cutoff_int-2)  = lj6_array(cutoff_int-1)
    table_ene(6*cutoff_int)    = elec_array(cutoff_int-1)
    table_grad(2)              = dlj12_array(1)
    table_grad(4)              = dlj6_array(1)
    table_grad(6)              = delec_array(1)
    table_grad(6*cutoff_int-4) = dlj12_array(cutoff_int-1)
    table_grad(6*cutoff_int-2) = dlj6_array(cutoff_int-1)
    table_grad(6*cutoff_int)   = delec_array(cutoff_int-1)

    deallocate(elec_array, delec_array)
    deallocate(lj12_array, dlj12_array)
    deallocate(lj6_array, dlj6_array)

    return

  end subroutine table_cutoff_cubic_fswitch

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    table_cutoff_cubic_groshift
  !> @brief        table for pme and cubic force switch
  !! @authors      CK, JJ
  !! @param[in]    switch_int  : integer switchdist
  !! @param[in]    cutoff_int  : integer cutoff distance
  !! @param[in]    cutoff_int2 : integer cutoff distance (buffering region)
  !! @param[in]    density     : table density
  !! @param[in]    el_fact     : inverse of dielectoric constant
  !! @param[in]    cutoff      : cutoff distance
  !! @param[in]    switchdist  : switch distance
  !! @param[out]   table_ene   : energy table
  !! @param[out]   table_grad  : gradient table
  !! @note         J.Jung et al., J.Comput.Chem., 34, 2412-2420 (2013)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine table_cutoff_cubic_groshift(switch_int, cutoff_int, cutoff_int2, &
                                  density, el_fact,                           &
                                  cutoff, switchdist,  table_ene, table_grad) 
    ! formal arguments
    integer,                 intent(in)    :: switch_int
    integer,                 intent(in)    :: cutoff_int
    integer,                 intent(in)    :: cutoff_int2
    real(wp),                intent(in)    :: density
    real(wp),                intent(in)    :: el_fact
    real(wp),                intent(in)    :: cutoff
    real(wp),                intent(in)    :: switchdist
    real(wp),                intent(inout) :: table_ene(:)
    real(wp),                intent(inout) :: table_grad(:)

    ! local variables
    integer                  :: i, icoun
    real(wp)                 :: rij, rij2
    real(wp)                 :: inv_r2, inv_rij, inv_r6, inv_r12
    real(wp)                 :: lj1, lj2, lj4
    real(wp)                 :: switch, dswitch, cutoff2, switchdist2
    real(wp)                 :: Coef_A12, Coef_B12, Coef_C12
    real(wp)                 :: Coef_A06, Coef_B06, Coef_C06
    real(wp)                 :: term_lj12_OO, term_lj12_HO, term_lj12_HH
    real(wp)                 :: term_lj6_OO, term_lj6_HO, term_lj6_HH
    real(wp)                 :: F12, F06, P12, P06

    real(wp),    allocatable :: elec_array(:), delec_array(:)
    real(wp),    allocatable :: lj12_array(:), dlj12_array(:)
    real(wp),    allocatable :: lj6_array(:), dlj6_array(:)

    cutoff2 = cutoff*cutoff
    switchdist2 = switchdist*switchdist

    allocate(elec_array(cutoff_int2), delec_array(cutoff_int2))
    allocate(lj12_array(cutoff_int2), dlj12_array(cutoff_int2))
    allocate(lj6_array (cutoff_int2), dlj6_array (cutoff_int2))

    Coef_A12  = 16.0_wp*cutoff - 13.0_wp*switchdist
    Coef_A12  = - Coef_A12 / cutoff**14 / (cutoff-switchdist)**2
    Coef_B12  = 15.0_wp*cutoff - 13.0_wp*switchdist
    Coef_B12  = Coef_B12 / cutoff**14 / (cutoff-switchdist)**3
    Coef_A06  = 10.0_wp*cutoff - 7.0_wp*switchdist
    Coef_A06  = - Coef_A06 / cutoff**8 / (cutoff-switchdist)**2
    Coef_B06  = 9.0_wp*cutoff - 7.0_wp*switchdist
    Coef_B06  = Coef_B06 / cutoff**8 / (cutoff-switchdist)**3
    Coef_C12  = 1.0_wp/cutoff**12
    Coef_C12  = Coef_C12 - 12.0_wp*Coef_A12/3.0_wp*(cutoff-switchdist)**3
    Coef_C12  = Coef_C12 - 12.0_wp*Coef_B12/4.0_wp*(cutoff-switchdist)**4
    Coef_C06  = 1.0_wp/cutoff**6
    Coef_C06  = Coef_C06 - 6.0_wp*Coef_A06/3.0_wp*(cutoff-switchdist)**3
    Coef_C06  = Coef_C06 - 6.0_wp*Coef_B06/4.0_wp*(cutoff-switchdist)**4

    do i = 2, switch_int

      rij2     = real(i,wp)/density
      rij      = sqrt(rij2)
      inv_r2   = 1.0_wp / rij2
      inv_rij  = 1.0_wp / rij
      inv_r6   = inv_r2 * inv_r2 * inv_r2
      inv_r12  = inv_r6 * inv_r6

      table_ene(6*i-5)  = inv_r12 - Coef_C12
      table_ene(6*i-3)  = inv_r6 - Coef_C06
      table_grad(6*i-5) = -12.0_wp*inv_r12*inv_r2
      table_grad(6*i-3) = -6.0_wp*inv_r6*inv_r2
      table_ene(6*i-1)  = el_fact * inv_rij
      table_grad(6*i-1) = -el_fact * inv_r2 * inv_rij

      lj12_array(i-1)  = table_ene(6*i-5)  - table_ene(6*i-11)
      lj6_array(i-1)   = table_ene(6*i-3)  - table_ene(6*i-9)
      dlj12_array(i-1) = table_grad(6*i-5) - table_grad(6*i-11)
      dlj6_array(i-1)  = table_grad(6*i-3) - table_grad(6*i-9)
      elec_array(i-1)  = table_ene(6*i-1)  - table_ene(6*i-7)
      delec_array(i-1) = table_grad(6*i-1) - table_grad(6*i-7)
    end do

    icoun = max(2, switch_int+1)

    do i = icoun, cutoff_int
!    do i = switch_int+1, cutoff_int
      rij2     = real(i,wp)/density
      rij      = sqrt(rij2)
      inv_r2   = 1.0_wp / rij2
      inv_rij  = 1.0_wp / rij
      inv_r6   = inv_r2 * inv_r2 * inv_r2
      inv_r12 = inv_r6 * inv_r6

      F12       = Coef_A12*(rij-switchdist)**2 + Coef_B12*(rij-switchdist)**3
      F06       = Coef_A06*(rij-switchdist)**2 + Coef_B06*(rij-switchdist)**3
      P12       = - 12.0_wp*Coef_A12/3.0_wp*(rij-switchdist)**3
      P12       = P12 - 12.0_wp*Coef_B12/4.0_wp*(rij-switchdist)**4 - Coef_C12
      P06       = - 6.0_wp*Coef_A06/3.0_wp*(rij-switchdist)**3
      P06       = P06 - 6.0_wp*Coef_B06/4.0_wp*(rij-switchdist)**4 - Coef_C06

      table_ene(6*i-5)  = inv_r12 + P12
      table_ene(6*i-3)  = inv_r6  + P06
      table_grad(6*i-5) = -12.0_wp*inv_r12*inv_r2 - 12.0_wp*F12*inv_rij
      table_grad(6*i-3) = -6.0_wp*inv_r6*inv_r2- 6.0_wp*F06*inv_rij
      table_ene(6*i-1)  = el_fact*inv_rij
      table_grad(6*i-1) = -el_fact*inv_rij*inv_r2

      lj12_array(i-1)  = table_ene(6*i-5)  - table_ene(6*i-11)
      lj6_array(i-1)   = table_ene(6*i-3)  - table_ene(6*i-9)
      dlj12_array(i-1) = table_grad(6*i-5) - table_grad(6*i-11)
      dlj6_array(i-1)  = table_grad(6*i-3) - table_grad(6*i-9)
      elec_array(i-1)  = table_ene(6*i-1)  - table_ene(6*i-7)
      delec_array(i-1) = table_grad(6*i-1) - table_grad(6*i-7)
    end do

    table_ene (6*cutoff_int-5:6*cutoff_int) = 0.0_wp
    table_grad(6*cutoff_int-5:6*cutoff_int) = 0.0_wp

    do i = 2, cutoff_int-1
      table_ene(6*i-4)  = (lj12_array(i-1) + lj12_array(i))   / 2.0_wp
      table_grad(6*i-4) = (dlj12_array(i-1) + dlj12_array(i)) / 2.0_wp
      table_ene(6*i-2)  = (lj6_array(i-1) + lj6_array(i))     / 2.0_wp
      table_grad(6*i-2) = (dlj6_array(i-1) + dlj6_array(i))   / 2.0_wp
      table_ene(6*i)    = (elec_array(i-1) + elec_array(i))   / 2.0_wp
      table_grad(6*i)   = (delec_array(i-1) + delec_array(i)) / 2.0_wp
    end do

    table_ene(2)               = lj12_array(1)
    table_ene(4)               = lj6_array(1)
    table_ene(6)               = elec_array(1)
    table_ene(6*cutoff_int-4)  = lj12_array(cutoff_int-1)
    table_ene(6*cutoff_int-2)  = lj6_array(cutoff_int-1)
    table_ene(6*cutoff_int)    = elec_array(cutoff_int-1)
    table_grad(2)              = dlj12_array(1)
    table_grad(4)              = dlj6_array(1)
    table_grad(6)              = delec_array(1)
    table_grad(6*cutoff_int-4) = dlj12_array(cutoff_int-1)
    table_grad(6*cutoff_int-2) = dlj6_array(cutoff_int-1)
    table_grad(6*cutoff_int)   = delec_array(cutoff_int-1)

    deallocate(elec_array, delec_array)
    deallocate(lj12_array, dlj12_array)
    deallocate(lj6_array, dlj6_array)

    return

  end subroutine table_cutoff_cubic_groshift

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    table_cutoff_cubic_groswitch
  !> @brief        table for cutoff and gromacs swith function for vdW
  !! @authors      CK, JJ
  !! @param[in]    switch_int  : integer switchdist
  !! @param[in]    cutoff_int  : integer cutoff distance
  !! @param[in]    cutoff_int2 : integer cutoff distance (buffering region)
  !! @param[in]    density     : table density
  !! @param[in]    el_fact     : inverse of dielectoric constant
  !! @param[in]    cutoff      : cutoff distance
  !! @param[in]    switchdist  : switch distance
  !! @param[out]   table_ene   : energy table
  !! @param[out]   table_grad  : gradient table
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine table_cutoff_cubic_groswitch(switch_int, cutoff_int, cutoff_int2, &
                                  density, el_fact,                           &
                                  cutoff, switchdist,  table_ene, table_grad) 
    ! formal arguments
    integer,                 intent(in)    :: switch_int
    integer,                 intent(in)    :: cutoff_int
    integer,                 intent(in)    :: cutoff_int2
    real(wp),                intent(in)    :: density
    real(wp),                intent(in)    :: el_fact
    real(wp),                intent(in)    :: cutoff
    real(wp),                intent(in)    :: switchdist
    real(wp),                intent(inout) :: table_ene(:)
    real(wp),                intent(inout) :: table_grad(:)

    ! local variables
    integer                  :: i, icoun
    real(wp)                 :: rij, rij1, rij2
    real(wp)                 :: inv_r2, inv_rij, inv_r6, inv_r12
    real(wp)                 :: switch, dswitch, cutoff2, switchdist2
    real(wp)                 :: coef_A, coef_B, coef_C
    real(wp)                 :: coef_A1, coef_B1, coef_C1

    real(wp),    allocatable :: elec_array(:), delec_array(:)
    real(wp),    allocatable :: lj12_array(:), dlj12_array(:)
    real(wp),    allocatable :: lj6_array(:), dlj6_array(:)

    cutoff2 = cutoff*cutoff
    switchdist2 = switchdist*switchdist

    allocate(elec_array(cutoff_int2), delec_array(cutoff_int2))
    allocate(lj12_array(cutoff_int2), dlj12_array(cutoff_int2))
    allocate(lj6_array (cutoff_int2), dlj6_array (cutoff_int2))

    do i = 2, switch_int
      rij2    = real(i,wp)/density
      rij     = sqrt(rij2)
      inv_r2  = 1.0_wp / rij2
      inv_rij = 1.0_wp / rij
      inv_r6  = inv_r2 * inv_r2 * inv_r2
      inv_r12 = inv_r6 * inv_r6

      table_ene(6*i-5)  = inv_r12 
      table_ene(6*i-3)  = inv_r6
      table_grad(6*i-5) = -12.0_wp*inv_r12*inv_r2
      table_grad(6*i-3) = -6.0_wp*inv_r6*inv_r2
      table_ene(6*i-1)  = el_fact*inv_rij
      table_grad(6*i-1) = -el_fact*inv_rij*inv_r2

      lj12_array(i-1)  = table_ene(6*i-5)  - table_ene(6*i-11)
      lj6_array(i-1)   = table_ene(6*i-3)  - table_ene(6*i-9)
      dlj12_array(i-1) = table_grad(6*i-5) - table_grad(6*i-11)
      dlj6_array(i-1)  = table_grad(6*i-3) - table_grad(6*i-9)
      elec_array(i-1)  = table_ene(6*i-1)  - table_ene(6*i-7)
      delec_array(i-1) = table_grad(6*i-1) - table_grad(6*i-7)
    end do

    coef_A    = -10.0_wp / (cutoff-switchdist)**3
    coef_B    =  15.0_wp / (cutoff-switchdist)**4
    coef_C    = - 6.0_wp / (cutoff-switchdist)**5
    coef_A1   = -30.0_wp / (cutoff-switchdist)**3
    coef_B1   =  60.0_wp / (cutoff-switchdist)**4
    coef_C1   = -30.0_wp / (cutoff-switchdist)**5

    icoun = max(2, switch_int+1)

    do i = icoun, cutoff_int
!    do i = switch_int+1, cutoff_int
      rij2     = real(i,wp)/density
      rij      = sqrt(rij2)
      rij1     = rij - switchdist
      switch    = 1.0_wp + coef_A*rij1**3 + coef_B*rij1**4 + coef_C*rij1**5
      dswitch   = coef_A1*rij1**2 + coef_B*rij1**3 + coef_C*rij1**4
      inv_r2   = 1.0_wp / rij2
      inv_rij  = 1.0_wp / rij
      inv_r6   = inv_r2 * inv_r2 * inv_r2
      inv_r12 = inv_r6 * inv_r6

      table_ene(6*i-5)  = inv_r12 * switch
      table_ene(6*i-3)  = inv_r6  * switch
      table_grad(6*i-5) = -12.0_wp*inv_r12*inv_r2*switch  &
                         + dswitch*inv_r12*inv_rij
      table_grad(6*i-3) = -6.0_wp*inv_r6*inv_r2*switch    &
                         + dswitch*inv_r6*inv_rij
      table_ene(6*i-1)  = el_fact*inv_rij
      table_grad(6*i-1) = -el_fact*inv_rij*inv_r2

      lj12_array(i-1)  = table_ene(6*i-5)  - table_ene(6*i-11)
      lj6_array(i-1)   = table_ene(6*i-3)  - table_ene(6*i-9)
      dlj12_array(i-1) = table_grad(6*i-5) - table_grad(6*i-11)
      dlj6_array(i-1)  = table_grad(6*i-3) - table_grad(6*i-9)
      elec_array(i-1)  = table_ene(6*i-1)  - table_ene(6*i-7)
      delec_array(i-1) = table_grad(6*i-1) - table_grad(6*i-7)
    end do

    table_ene (6*cutoff_int-5:6*cutoff_int) = 0.0_wp
    table_grad(6*cutoff_int-5:6*cutoff_int) = 0.0_wp

    do i = 2, cutoff_int-1
      table_ene(6*i-4)  = (lj12_array(i-1) + lj12_array(i))   / 2.0_wp
      table_grad(6*i-4) = (dlj12_array(i-1) + dlj12_array(i)) / 2.0_wp
      table_ene(6*i-2)  = (lj6_array(i-1) + lj6_array(i))     / 2.0_wp
      table_grad(6*i-2) = (dlj6_array(i-1) + dlj6_array(i))   / 2.0_wp
      table_ene(6*i)    = (elec_array(i-1) + elec_array(i))   / 2.0_wp
      table_grad(6*i)   = (delec_array(i-1) + delec_array(i)) / 2.0_wp
    end do

    table_ene(2)               = lj12_array(1)
    table_ene(4)               = lj6_array(1)
    table_ene(6)               = elec_array(1)
    table_ene(6*cutoff_int-4)  = lj12_array(cutoff_int-1)
    table_ene(6*cutoff_int-2)  = lj6_array(cutoff_int-1)
    table_ene(6*cutoff_int)    = elec_array(cutoff_int-1)
    table_grad(2)              = dlj12_array(1)
    table_grad(4)              = dlj6_array(1)
    table_grad(6)              = delec_array(1)
    table_grad(6*cutoff_int-4) = dlj12_array(cutoff_int-1)
    table_grad(6*cutoff_int-2) = dlj6_array(cutoff_int-1)
    table_grad(6*cutoff_int)   = delec_array(cutoff_int-1)

    deallocate(elec_array, delec_array)
    deallocate(lj12_array, dlj12_array)
    deallocate(lj6_array, dlj6_array)

    return

  end subroutine table_cutoff_cubic_groswitch

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    table_cutoff_cubic_grodoubleshift
  !> @brief        table for pme and cubic force switch
  !! @authors      CK, JJ
  !! @param[in]    switch_int  : integer switchdist
  !! @param[in]    cutoff_int  : integer cutoff distance
  !! @param[in]    cutoff_int2 : integer cutoff distance (buffering region)
  !! @param[in]    density     : table density
  !! @param[in]    el_fact     : inverse of dielectoric constant
  !! @param[in]    cutoff      : cutoff distance
  !! @param[in]    switchdist  : switch distance
  !! @param[out]   table_ene   : energy table
  !! @param[out]   table_grad  : gradient table
  !! @note         J.Jung et al., J.Comput.Chem., 34, 2412-2420 (2013)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine table_cutoff_cubic_grodoubleshift(switch_int,                    &
                                  cutoff_int, cutoff_int2,                    &
                                  density, el_fact,                           &
                                  cutoff, switchdist,  table_ene, table_grad) 
    ! formal arguments
    integer,                 intent(in)    :: switch_int
    integer,                 intent(in)    :: cutoff_int
    integer,                 intent(in)    :: cutoff_int2
    real(wp),                intent(in)    :: density
    real(wp),                intent(in)    :: el_fact
    real(wp),                intent(in)    :: cutoff
    real(wp),                intent(in)    :: switchdist
    real(wp),                intent(inout) :: table_ene(:)
    real(wp),                intent(inout) :: table_grad(:)

    ! local variables
    integer                  :: i
    integer                  :: icoun
    real(wp)                 :: rij, rij2
    real(wp)                 :: inv_r2, inv_rij, inv_r6, inv_r12
    real(wp)                 :: lj1, lj2, lj4
    real(wp)                 :: switch, dswitch, cutoff2, switchdist2
    real(wp)                 :: Coef_A12, Coef_B12, Coef_C12
    real(wp)                 :: Coef_A06, Coef_B06, Coef_C06
    real(wp)                 :: Coef_A01, Coef_B01, Coef_C01
    real(wp)                 :: term_lj12_OO, term_lj12_HO, term_lj12_HH
    real(wp)                 :: term_lj6_OO, term_lj6_HO, term_lj6_HH
    real(wp)                 :: F12, F06, F01, P12, P06, P01

    real(wp),    allocatable :: elec_array(:), delec_array(:)
    real(wp),    allocatable :: lj12_array(:), dlj12_array(:)
    real(wp),    allocatable :: lj6_array(:), dlj6_array(:)

    cutoff2 = cutoff*cutoff
    switchdist2 = switchdist*switchdist

    allocate(elec_array(cutoff_int2), delec_array(cutoff_int2))
    allocate(lj12_array(cutoff_int2), dlj12_array(cutoff_int2))
    allocate(lj6_array (cutoff_int2), dlj6_array (cutoff_int2))

    Coef_A12  = 16.0_wp*cutoff - 13.0_wp*switchdist
    Coef_A12  = - Coef_A12 / cutoff**14 / (cutoff-switchdist)**2
    Coef_B12  = 15.0_wp*cutoff - 13.0_wp*switchdist
    Coef_B12  = Coef_B12 / cutoff**14 / (cutoff-switchdist)**3
    Coef_A06  = 10.0_wp*cutoff - 7.0_wp*switchdist
    Coef_A06  = - Coef_A06 / cutoff**8 / (cutoff-switchdist)**2
    Coef_B06  = 9.0_wp*cutoff - 7.0_wp*switchdist
    Coef_B06  = Coef_B06 / cutoff**8 / (cutoff-switchdist)**3
    Coef_C12  = 1.0_wp/cutoff**12
    Coef_C12  = Coef_C12 - 12.0_wp*Coef_A12/3.0_wp*(cutoff-switchdist)**3
    Coef_C12  = Coef_C12 - 12.0_wp*Coef_B12/4.0_wp*(cutoff-switchdist)**4
    Coef_C06  = 1.0_wp/cutoff**6
    Coef_C06  = Coef_C06 - 6.0_wp*Coef_A06/3.0_wp*(cutoff-switchdist)**3
    Coef_C06  = Coef_C06 - 6.0_wp*Coef_B06/4.0_wp*(cutoff-switchdist)**4

!  switch=0 in coulomb
    Coef_A01  = 5.0_wp*cutoff
    Coef_A01  = - Coef_A01 / cutoff**3 / (cutoff)**2
    Coef_B01  = 4.0_wp*cutoff
    Coef_B01  = Coef_B01 / cutoff**3 / (cutoff)**3
    Coef_C01  = 1.0_wp/cutoff
    Coef_C01  = Coef_C01 - Coef_A01/3.0_wp*(cutoff)**3
    Coef_C01  = Coef_C01 - Coef_B01/4.0_wp*(cutoff)**4

    do i = 2, switch_int
      rij2    = real(i,wp)/density
      rij     = sqrt(rij2)
      inv_r2  = 1.0_wp / rij2
      inv_rij = 1.0_wp / rij
      inv_r6  = inv_r2 * inv_r2 * inv_r2
      inv_r12 = inv_r6 * inv_r6

!  switch=0 in coulomb
      F01       = Coef_A01*(rij)**2 + Coef_B01*(rij)**3
      P01       = - Coef_A01/3.0_wp*(rij)**3
      P01       = P01 - Coef_B01/4.0_wp*(rij)**4 - Coef_C01

      table_ene(6*i-5)  = inv_r12 - Coef_C12
      table_ene(6*i-3)  = inv_r6 - Coef_C06
      table_grad(6*i-5) = -12.0_wp*inv_r12*inv_r2
      table_grad(6*i-3) = -6.0_wp*inv_r6*inv_r2
      table_ene(6*i-1)  = el_fact*(inv_rij + P01)
      table_grad(6*i-1) = -el_fact*(inv_rij*inv_r2 + F01*inv_rij)

      lj12_array(i-1)  = table_ene(6*i-5)  - table_ene(6*i-11)
      lj6_array(i-1)   = table_ene(6*i-3)  - table_ene(6*i-9)
      dlj12_array(i-1) = table_grad(6*i-5) - table_grad(6*i-11)
      dlj6_array(i-1)  = table_grad(6*i-3) - table_grad(6*i-9)
      elec_array(i-1)  = table_ene(6*i-1)  - table_ene(6*i-7)
      delec_array(i-1) = table_grad(6*i-1) - table_grad(6*i-7)
    end do

    icoun = max(2, switch_int+1)

    do i = icoun, cutoff_int
!    do i = switch_int+1, cutoff_int
      rij2     = real(i,wp)/density
      rij      = sqrt(rij2)
      inv_r2   = 1.0_wp / rij2
      inv_rij  = 1.0_wp / rij
      inv_r6   = inv_r2 * inv_r2 * inv_r2
      inv_r12 = inv_r6 * inv_r6

      F12       = Coef_A12*(rij-switchdist)**2 + Coef_B12*(rij-switchdist)**3
      F06       = Coef_A06*(rij-switchdist)**2 + Coef_B06*(rij-switchdist)**3
      P12       = - 12.0_wp*Coef_A12/3.0_wp*(rij-switchdist)**3
      P12       = P12 - 12.0_wp*Coef_B12/4.0_wp*(rij-switchdist)**4 - Coef_C12
      P06       = - 6.0_wp*Coef_A06/3.0_wp*(rij-switchdist)**3
      P06       = P06 - 6.0_wp*Coef_B06/4.0_wp*(rij-switchdist)**4 - Coef_C06

!  switch=0 in coulomb
      F01       = Coef_A01*(rij)**2 + Coef_B01*(rij)**3
      P01       = - Coef_A01/3.0_wp*(rij)**3
      P01       = P01 - Coef_B01/4.0_wp*(rij)**4 - Coef_C01

      table_ene(6*i-5)  = inv_r12 + P12
      table_ene(6*i-3)  = inv_r6  + P06
      table_grad(6*i-5) = -12.0_wp*inv_r12*inv_r2 - 12.0_wp*F12*inv_rij
      table_grad(6*i-3) = -6.0_wp*inv_r6*inv_r2- 6.0_wp*F06*inv_rij
      table_ene(6*i-1)  = el_fact*(inv_rij + P01)
      table_grad(6*i-1) = -el_fact*(inv_rij*inv_r2 + F01*inv_rij)

      lj12_array(i-1)  = table_ene(6*i-5)  - table_ene(6*i-11)
      lj6_array(i-1)   = table_ene(6*i-3)  - table_ene(6*i-9)
      dlj12_array(i-1) = table_grad(6*i-5) - table_grad(6*i-11)
      dlj6_array(i-1)  = table_grad(6*i-3) - table_grad(6*i-9)
      elec_array(i-1)  = table_ene(6*i-1)  - table_ene(6*i-7)
      delec_array(i-1) = table_grad(6*i-1) - table_grad(6*i-7)
    end do

    table_ene (6*cutoff_int-5:6*cutoff_int) = 0.0_wp
    table_grad(6*cutoff_int-5:6*cutoff_int) = 0.0_wp

    do i = 2, cutoff_int-1
      table_ene(6*i-4)  = (lj12_array(i-1) + lj12_array(i))   / 2.0_wp
      table_grad(6*i-4) = (dlj12_array(i-1) + dlj12_array(i)) / 2.0_wp
      table_ene(6*i-2)  = (lj6_array(i-1) + lj6_array(i))     / 2.0_wp
      table_grad(6*i-2) = (dlj6_array(i-1) + dlj6_array(i))   / 2.0_wp
      table_ene(6*i)    = (elec_array(i-1) + elec_array(i))   / 2.0_wp
      table_grad(6*i)   = (delec_array(i-1) + delec_array(i)) / 2.0_wp
    end do

    table_ene(2)               = lj12_array(1)
    table_ene(4)               = lj6_array(1)
    table_ene(6)               = elec_array(1)
    table_ene(6*cutoff_int-4)  = lj12_array(cutoff_int-1)
    table_ene(6*cutoff_int-2)  = lj6_array(cutoff_int-1)
    table_ene(6*cutoff_int)    = elec_array(cutoff_int-1)
    table_grad(2)              = dlj12_array(1)
    table_grad(4)              = dlj6_array(1)
    table_grad(6)              = delec_array(1)
    table_grad(6*cutoff_int-4) = dlj12_array(cutoff_int-1)
    table_grad(6*cutoff_int-2) = dlj6_array(cutoff_int-1)
    table_grad(6*cutoff_int)   = delec_array(cutoff_int-1)

    deallocate(elec_array, delec_array)
    deallocate(lj12_array, dlj12_array)
    deallocate(lj6_array, dlj6_array)

    return

  end subroutine table_cutoff_cubic_grodoubleshift

end module table_libs_mod
