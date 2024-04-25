!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_enefunc_table_mod
!> @brief   define lookup table for energy calculation
!! @authors Jaewoon Jung (JJ), Chigusa Kobayashi (CK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_enefunc_table_mod

  use sp_boundary_mod
  use sp_energy_pme_mod
  use sp_energy_mod
  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use molecules_str_mod
  use table_libs_mod
  use mpi_parallel_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! variables
  real(wp),         save      :: el_fact    ! e^2/(4*pi*eps) A*kcal/mol
  real(wp),         save      :: alpha      ! alpha
  real(wp),         save      :: alpha2m    ! -alpha^2
  real(wp),         save      :: alpha2sp   ! 2*alpha/sqrt(pi)

  public  :: setup_enefunc_table
  private :: setup_table_general_pme_linear
  private :: setup_table_general_cutoff_cubic
  private :: setup_table_water_pme_linear

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_table
  !> @brief        define lookup table of potential energy function
  !! @authors      JJ
  !! @param[in]    ene_info : ENERGY section control parameters information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_table(ene_info, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    real(wp)                 :: cutoff, cutoff2
    real(wp)                 :: cutoff2_water, cutoff2_water1, density
    real(wp)                 :: mind
    integer                  :: cutoff_int, cutoff_int2, cutoff_int1


    if (.not. ene_info%table) &
      return

    ! set cutoff distance etc.
    !
    enefunc%table%table       = ene_info%table
    enefunc%table%table_order = ene_info%table_order
    enefunc%table%water_model = ene_info%water_model
    enefunc%table%density     = ene_info%table_density

    ! setting common parameters
    !
    el_fact  = ELECOEF/ene_info%dielec_const
    alpha    = ene_info%pme_alpha
    alpha2m  = -alpha**2
    alpha2sp = 2.0_wp*alpha/sqrt(PI)

    ! set parameters for table
    !
    cutoff   = ene_info%cutoffdist
    density  = ene_info%table_density

    cutoff2        = cutoff * cutoff
    cutoff2_water  = (cutoff+4.5_wp) * (cutoff+4.5_wp)
    cutoff2_water1 = (cutoff+1.5_wp) * (cutoff+1.5_wp)
    cutoff_int     = int(cutoff2*density)
    cutoff_int2    = int(cutoff2_water*density)
    cutoff_int1    = int(cutoff2_water1*density)

    if (enefunc%table%table_order == 1) then

      cutoff_int = int(2.1_wp*cutoff2_water1*density)
      cutoff_int2 = cutoff_int+1

    end if

    enefunc%table%cutoff_int = cutoff_int2

    if (enefunc%contact_check .or. enefunc%nonb_limiter) then
      if (enefunc%table%table_order /= 1) then
        if (main_rank) &
          write(MsgOut,'(A)') 'Setup_Enefunc_Table> Warning:'//&
          'concact_check is only in table_order =1'
      else
        mind=cutoff2*density/real(cutoff_int2,wp)+0.001_wp
        enefunc%minimum_contact=max(mind, enefunc%minimum_contact)
      endif
    endif

    if (enefunc%table%water_table) then
      call alloc_enefunc(enefunc, EneFuncTblWatDomain, cutoff_int2, 1)

    else
      call alloc_enefunc(enefunc, EneFuncTableDomain, cutoff_int2, 1)

    end if

    if (enefunc%pme_use) then

        call setup_table_general_pme_linear(ene_info, cutoff2, cutoff_int,  &
                                            cutoff_int2, enefunc)

        if (enefunc%table%water_table) then

          call setup_table_water_pme_linear(ene_info, cutoff2, cutoff_int,  &
                                            cutoff_int2, enefunc)

        end if

    else

        call setup_table_general_cutoff_cubic(ene_info, cutoff2, cutoff_int,&
                                              cutoff_int2, enefunc)

        if (enefunc%table%water_table) then
          call error_msg('Setup_Enefunc_Table> watermodel is not'//         &
                         ' allowed in cutoff')
        end if


    end if

    return

  end subroutine setup_enefunc_table

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_table_general_pme_linear
  !> @brief        define lookup table of potential energy function
  !!               (linear)
  !! @authors      JJ, CK
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[in]    cutoff2     : square cutoff distance
  !! @param[in]    cutoff_int  : integer cutoff distance
  !! @param[in]    cutoff_int1 : square integer cutoff distance
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_table_general_pme_linear(ene_info, cutoff2, &
                                        cutoff_int, cutoff_int1, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    real(wp),                intent(in)    :: cutoff2
    integer,                 intent(in)    :: cutoff_int
    integer,                 intent(in)    :: cutoff_int1
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    real(wp)                 :: switchdist, cutoff
    real(wp)                 :: switchdist2
    real(wp)                 :: density
    integer                  :: i, switch_int

    real(wp),        pointer :: table_ene(:), table_grad(:)
    real(wp),        pointer :: table_ecor(:), table_decor(:)


    table_ene   => enefunc%table%table_ene
    table_ecor  => enefunc%table%table_ecor
    table_decor => enefunc%table%table_decor
    table_grad  => enefunc%table%table_grad

    density     = enefunc%table%density

    ! set cutoff and switch distance etc.
    !
    switchdist  = ene_info%switchdist
    cutoff      = ene_info%cutoffdist

    ! set lennard-jones parameters
    !
    switchdist2 = switchdist * switchdist
    switch_int  = int(1.0_wp/switchdist2*cutoff2*density)

    table_ene  (1:3*cutoff_int1) = 0.0_wp
    table_ecor (1:cutoff_int1)   = 0.0_wp
    table_decor(1:cutoff_int1)   = 0.0_wp
    table_grad (1:3*cutoff_int1) = 0.0_wp

    if (ene_info%vdw_force_switch) then
      call table_pme_linear_fswitch(switch_int, cutoff_int, density,      &
                                   cutoff2, switchdist2, el_fact, alpha,     &
                                   alpha2sp, alpha2m,                        &
                                   table_ene, table_grad, table_ecor,        &
                                   table_decor)

    else if (enefunc%forcefield == ForcefieldAMBER) then

      call table_pme_linear_noswitch(switch_int, cutoff_int, density,        &
                                   cutoff2, switchdist2, el_fact, alpha,     &
                                   alpha2sp, alpha2m,                        &
                                   table_ene, table_grad, table_ecor,        &
                                   table_decor)

    else if (enefunc%forcefield == ForcefieldGROAMBER .or.                 &
            enefunc%forcefield == ForcefieldGROMARTINI) then
!     if (abs(switchdist) < EPS) then
      if (ene_info%vdw_shift) then
        call table_pme_linear_groshift(switch_int, cutoff_int, density,      &
                                   cutoff, switchdist, el_fact, alpha,       &
                                   alpha2sp, alpha2m,                        &
                                   table_ene, table_grad, table_ecor,        &
                                   table_decor)

      else if (abs(switchdist-cutoff) < EPS) then
        call table_pme_linear_noswitch(switch_int, cutoff_int, density,      &
                                   cutoff2, switchdist2, el_fact, alpha,     &
                                   alpha2sp, alpha2m,                        &
                                   table_ene, table_grad, table_ecor,        &
                                   table_decor)
      else
        call table_pme_linear_groswitch(switch_int, cutoff_int, density,     &
                                   cutoff, switchdist, el_fact, alpha,       &
                                   alpha2sp, alpha2m,                        &
                                   table_ene, table_grad, table_ecor,        &
                                   table_decor, enefunc%eswitch,             &
                                   enefunc%vswitch)

      endif
    else 
      call table_pme_linear_switch(switch_int, cutoff_int, density,          &
                                   cutoff2, switchdist2, el_fact, alpha,     &
                                   alpha2sp, alpha2m,                        &
                                   table_ene, table_grad, table_ecor,        &
                                   table_decor)

    end if

    return

  end subroutine setup_table_general_pme_linear

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_table_general_cutoff_cubic
  !> @brief        define lookup table for cutoff
  !! @authors      JJ, CK
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[in]    cutoff2     : square cutoff distance
  !! @param[in]    cutoff_int  : integer cutoff distance
  !! @param[in]    cutoff_int1 : square integer cutoff distance
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_table_general_cutoff_cubic(ene_info, cutoff2,  &
                                        cutoff_int, cutoff_int1, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    real(wp),                intent(in)    :: cutoff2
    integer,                 intent(in)    :: cutoff_int
    integer,                 intent(in)    :: cutoff_int1
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    real(wp)                 :: switchdist
    real(wp)                 :: switchdist2, cutoff2_inv, cutoff
    real(wp)                 :: density
    integer                  :: i, switch_int

    real(wp),        pointer :: table_ene(:)
    real(wp),        pointer :: table_grad(:)

    density     = enefunc%table%density
    table_ene  => enefunc%table%table_ene
    table_grad => enefunc%table%table_grad

    cutoff2_inv = 1.0_wp/cutoff2
    el_fact = ELECOEF/enefunc%dielec_const 

    ! set switch distance etc.
    !
    switchdist  = ene_info%switchdist
    cutoff      = ene_info%cutoffdist

    ! set lennard-jones parameters
    !
    switchdist2 = switchdist * switchdist
    switch_int  = int(switchdist2*density)

    table_ene(1:6*cutoff_int1)  = 0.0_wp
    table_grad(1:6*cutoff_int1) = 0.0_wp

    if (ene_info%vdw_force_switch) then

      call table_cutoff_cubic_fswitch(switch_int, cutoff_int, cutoff_int1,  &
                                     density, el_fact,                      &
                                     cutoff2, switchdist2, table_ene,       &
                                     table_grad)


    else if (enefunc%forcefield == ForcefieldAMBER) then

      call table_cutoff_cubic_noswitch(switch_int, cutoff_int, cutoff_int1, &
                                     density, el_fact,                      &
                                     cutoff2, switchdist2, table_ene,       &
                                     table_grad)

    else if (enefunc%forcefield == ForcefieldGROAMBER .or.                &
            enefunc%forcefield == ForcefieldGROMARTINI) then

      if (ene_info%vdw_shift) then

        if (enefunc%forcefield == ForcefieldGROAMBER) then
          call table_cutoff_cubic_groshift(switch_int, cutoff_int,          &
                                         cutoff_int1, density, el_fact,     &
                                         cutoff, switchdist, table_ene,     &
                                         table_grad)
        else
          call table_cutoff_cubic_grodoubleshift(switch_int, cutoff_int,    &
                                         cutoff_int1, density, el_fact,     &
                                         cutoff, switchdist, table_ene,     &
                                         table_grad)
        endif


      else if (abs(switchdist-cutoff) < EPS) then 

        call table_cutoff_cubic_noswitch(switch_int, cutoff_int, cutoff_int1,&
                                     density, el_fact,                       &
                                     cutoff2, switchdist2, table_ene,        &
                                     table_grad)
      else

        call table_cutoff_cubic_groswitch(switch_int, cutoff_int, cutoff_int1,&
                                         density, el_fact,                   &
                                         cutoff, switchdist, table_ene,      &
                                         table_grad)
      endif

    else 

      call table_cutoff_cubic_switch(switch_int, cutoff_int, cutoff_int1,    &
                                     density, el_fact,                       &
                                     cutoff2, switchdist2, table_ene,        &
                                     table_grad)

    end if

    return

  end subroutine setup_table_general_cutoff_cubic

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_table_water_pme_linear
  !> @brief        define lookup table of water (linear interpolation)
  !! @authors      JJ, CK
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[in]    cutoff2     : square cutoff distance
  !! @param[in]    cutoff_int  : integer cutoff distance
  !! @param[in]    cutoff_int1 : square integer cutoff distance
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_table_water_pme_linear(ene_info, cutoff2, &
                                      cutoff_int, cutoff_int1, enefunc)

    ! formal arguments
    type(s_ene_info),         intent(in)    :: ene_info
    real(wp),                 intent(in)    :: cutoff2
    integer,                  intent(in)    :: cutoff_int
    integer,                  intent(in)    :: cutoff_int1
    type(s_enefunc),  target, intent(inout) :: enefunc

    ! local variables
    real(wp)                  :: lj6_OO, lj6_OH, lj6_HH
    real(wp)                  :: lj12_OO, lj12_OH, lj12_HH
    real(wp)                  :: cc_OO, cc_OH, cc_HH, charge_O, charge_H
    real(wp)                  :: switchdist, cutoff
    real(wp)                  :: switchdist2
    real(wp)                  :: density
    integer                   :: i, switch_int, atmcls_O, atmcls_H

    real(wp),         pointer :: table_de_WW(:,:)
    real(wp),         pointer :: table_ene_WW(:,:)


    table_ene_WW => enefunc%table%table_ene_WW
    table_de_WW  => enefunc%table%table_de_WW

    density     = enefunc%table%density

    ! set switch distance etc.
    !
    cutoff      = ene_info%cutoffdist
    switchdist  = ene_info%switchdist

    ! set lennard-jones parameters
    !
    switchdist2 = switchdist * switchdist
    switch_int = int(1.0_wp/switchdist2*cutoff2*density)

    table_ene_WW(1:6*cutoff_int1,1:3) = 0.0_wp
    table_de_WW(1:2*cutoff_int1,1:3)  = 0.0_wp

    ! atom class and charge of water oxygen and hydrogen
    !
    atmcls_O = enefunc%table%atom_cls_no_O
    atmcls_H = enefunc%table%atom_cls_no_H
    charge_O = enefunc%table%charge_O
    charge_H = enefunc%table%charge_H

    lj6_OO  = enefunc%nonb_lj6(atmcls_O,atmcls_O)
    lj6_OH  = enefunc%nonb_lj6(atmcls_O,atmcls_H)
    lj6_HH  = enefunc%nonb_lj6(atmcls_H,atmcls_H)

    lj12_OO = enefunc%nonb_lj12(atmcls_O,atmcls_O)
    lj12_OH = enefunc%nonb_lj12(atmcls_O,atmcls_H)
    lj12_HH = enefunc%nonb_lj12(atmcls_H,atmcls_H)

    cc_OO = charge_O * charge_O * el_fact
    cc_OH = charge_O * charge_H * el_fact
    cc_HH = charge_H * charge_H * el_fact

    if (ene_info%vdw_force_switch) then


      call table_water_pme_linear_fswitch(switch_int, cutoff_int, density,     &
                                     cutoff2, switchdist2, el_fact, alpha,     &
                                     lj6_OO,  lj6_OH,  lj6_HH,                 &
                                     lj12_OO, lj12_OH, lj12_HH,                &
                                     cc_OO, cc_OH, cc_HH,                      &
                                     alpha2sp, alpha2m,                        &
                                     table_ene_WW, table_de_WW)

    else if (enefunc%forcefield == ForcefieldAMBER) then

      call table_water_pme_linear_noswitch(switch_int, cutoff_int, density,    &
                                     cutoff2, switchdist2, el_fact, alpha,     &
                                     lj6_OO,  lj6_OH,  lj6_HH,                 &
                                     lj12_OO, lj12_OH, lj12_HH,                &
                                     cc_OO, cc_OH, cc_HH,                      &
                                     alpha2sp, alpha2m,                        &
                                     table_ene_WW, table_de_WW)

    else if (enefunc%forcefield == ForcefieldGROAMBER .or.                   &
            enefunc%forcefield == ForcefieldGROMARTINI) then
      if (ene_info%vdw_shift) then
        call table_water_pme_linear_groshift(switch_int, cutoff_int, density,  &
                                     cutoff, switchdist, el_fact, alpha,       &
                                     lj6_OO,  lj6_OH,  lj6_HH,                 &
                                     lj12_OO, lj12_OH, lj12_HH,                &
                                     cc_OO, cc_OH, cc_HH,                      &
                                     alpha2sp, alpha2m,                        &
                                     table_ene_WW, table_de_WW)
     else if (abs(switchdist-cutoff) < EPS) then
      call table_water_pme_linear_noswitch(switch_int, cutoff_int, density,    &
                                     cutoff2, switchdist2, el_fact, alpha,     &
                                     lj6_OO,  lj6_OH,  lj6_HH,                 &
                                     lj12_OO, lj12_OH, lj12_HH,                &
                                     cc_OO, cc_OH, cc_HH,                      &
                                     alpha2sp, alpha2m,                        &
                                     table_ene_WW, table_de_WW)
     else
        call table_water_pme_linear_groswitch(switch_int, cutoff_int, density, &
                                     cutoff, switchdist, el_fact, alpha,       &
                                     lj6_OO,  lj6_OH,  lj6_HH,                 &
                                     lj12_OO, lj12_OH, lj12_HH,                &
                                     cc_OO, cc_OH, cc_HH,                      &
                                     alpha2sp, alpha2m,                        &
                                     table_ene_WW, table_de_WW)

     endif

    else

      call table_water_pme_linear_switch(switch_int, cutoff_int, density,      &
                                     cutoff2, switchdist2, el_fact, alpha,     &
                                     lj6_OO,  lj6_OH,  lj6_HH,                 &
                                     lj12_OO, lj12_OH, lj12_HH,                &
                                     cc_OO, cc_OH, cc_HH,                      &
                                     alpha2sp, alpha2m,                        &
                                     table_ene_WW, table_de_WW)

    end if

    return

  end subroutine setup_table_water_pme_linear

end module sp_enefunc_table_mod
