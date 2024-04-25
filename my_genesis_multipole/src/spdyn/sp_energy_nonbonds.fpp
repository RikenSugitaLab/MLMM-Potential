!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_energy_nonbonds_mod
!> @brief   calculate nonbond energy
!! @authors Jaewoon Jung(JJ)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_energy_nonbonds_mod

  use sp_energy_pme_mod
  use sp_energy_table_cubic_mod
  use sp_energy_table_linear_bondcorr_mod
  use sp_energy_table_linear_nowater_mod
  use sp_energy_table_linear_water_mod
  use sp_energy_table_linear_mod
  use sp_pairlist_str_mod
  use sp_boundary_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public  :: compute_energy_nonbond_cutoff
  public  :: compute_energy_nonbond_pme
  public  :: compute_energy_nonbond_pme_short
  public  :: compute_energy_nonbond_pme_long
  !FEP
  public  :: compute_energy_nonbond_cutoff_fep
  public  :: compute_energy_nonbond_pme_fep
  public  :: compute_energy_nonbond_pme_short_fep

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_cutoff
  !> @brief        compute nonbond energy with cutoff 
  !! @authors      JJ
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    pairlist : pair-list information
  !! @param[in]    nonb_ene : flag for calculate nonbonded energy
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : lennard-jones energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_cutoff(domain, enefunc, pairlist, &
                                           nonb_ene, force, virial, eelec, evdw)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    logical,                 intent(in)    :: nonb_ene
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: virial(3,3,nthread)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw(nthread)


    call timer(TimerNonBond, TimerOn)

    if (enefunc%vacuum) then
      ! Without lookup table

      call compute_energy_nonbond14_notable( &
        domain, enefunc, &
        force, virial, eelec, evdw)

      if (nonb_ene) then

        call compute_energy_nonbond_notable( &
          domain, enefunc, pairlist, &
          force, virial, eelec, evdw)

      else

        call compute_force_nonbond_notable( &
          domain, enefunc, pairlist, &
          force, virial)

      end if

    else

      call compute_energy_nonbond14_table( &
                                    domain, enefunc, &
                                    force, virial, eelec, evdw)
   
      ! ==> Type 4 and 10
      if (nonb_ene) then

        call compute_energy_nonbond_table( &
                                    domain, enefunc, pairlist, &
                                    force, virial, eelec, evdw)
      
      else

        call compute_force_nonbond_table( &
                                    domain, enefunc, pairlist, &
                                    force, virial)
      
      end if

    end if

    call timer(TimerNonBond, TimerOff)
 
    return

  end subroutine compute_energy_nonbond_cutoff

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_pme
  !> @brief        Calculate nonbond energy by PME
  !! @authors      JJ, CK
  !! @param[in]    domain       : domain information
  !! @param[in]    enefunc      : potential energy functions information
  !! @param[in]    pairlist     : pair-list information
  !! @param[in]    boundary     : boundary information
  !! @param[in]    npt          : flag for NPT or not
  !! @param[in]    nonb_ene     : flag for calculate nonbonded energy
  !! @param[in]    nonb_limiter : flag for contact checker
  !! @param[inout] coord_pbc    : !TODO
  !! @param[inout] force        : forces of target systems
  !! @param[inout] force_pbc    : !TODO
  !! @param[inout] virial_cell  : !TODO
  !! @param[inout] virial       : virial term of target systems
  !! @param[inout] eelec        : electrostatic energy of target systems
  !! @param[inout] evdw         : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_pme(domain, enefunc, pairlist, boundary, &
                                        npt, nonb_ene, nonb_limiter,         &
                                        coord_pbc,                           &
                                        force_long, force, force_pbc,  &
                                        virial_cell, virial,      &
                                        eelec, evdw)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: nonb_ene
    logical,                 intent(in)    :: nonb_limiter
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(dp),                intent(inout) :: force_long(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3,nthread)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw(nthread)

    ! local
    real(dp) :: ene_virial(1:5)

    call timer(TimerNonBond, TimerOn)
 

    ! Calculate PME real part
    !
    if (real_calc) then

      call timer(TimerPmeReal, TimerOn)

#ifndef USE_GPU
      if (nonb_limiter) then
        call compute_energy_nonbond_table_linear_check( &
                                    domain, enefunc, pairlist, &
                                    coord_pbc, force_pbc, virial_cell, &
                                    eelec, evdw)
      else

        if (nonb_ene) then
         
          call compute_energy_nonbond_table_linear( &
                                      domain, enefunc, pairlist, &
                                      coord_pbc, force_pbc, virial_cell, &
                                      eelec, evdw)
         
        else
         
          call compute_force_nonbond_table_linear( &
                                      domain, enefunc, pairlist, &
                                      coord_pbc, force_pbc, virial_cell)
         
        end if
      end if

#else
      if (nonb_ene) then
        call compute_energy_nonbond_table_linear_gpu( &
                                  domain, enefunc, pairlist, &
                                  npt, coord_pbc, force_pbc, virial, &
                                  eelec, evdw, ene_virial)
      else
        call compute_force_nonbond_table_linear_gpu( &
                                  domain, enefunc, pairlist, npt, &
                                  .false., coord_pbc, force,      &
                                  force_pbc, virial, ene_virial)
      end if
#endif

#ifndef USE_GPU
      call timer(TimerPmeReal, TimerOff)
#endif

    end if


    ! Calculate PME reciprocal part
    !
    if (reciprocal_calc) then
 
      call timer(TimerPmeRecip, TimerOn)
 
      if (npt) call pme_pre(domain, boundary)

      call pme_recip(domain, force_long, virial, eelec)

      call timer(TimerPmeRecip, TimerOff)

   end if

  ! Add self energy
  !
    eelec(1) = eelec(1) + u_self

#ifndef USE_GPU
    call timer(TimerNonBond, TimerOff)
#endif
 
    return

  end subroutine compute_energy_nonbond_pme

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_pme_short
  !> @brief        Calculate nonbond energy (real part) by PME
  !! @authors      JJ
  !! @param[in]    domain      : domain information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[in]    pairlist    : pair-list information
  !! @param[in]    boundary    : boundary information
  !! @param[inout] coord_pbc   : !TODO
  !! @param[inout] force       : forces of target systems
  !! @param[inout] force_pbc   : !TODO
  !! @param[inout] virial_cell : !TODO
  !! @param[inout] virial      : virial term of target systems
  !! @param[inout] eelec       : electrostatic energy of target systems
  !! @param[inout] evdw        : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_pme_short(domain, enefunc, pairlist, &
                                              npt, nonb_ene, coord_pbc, &
                                              force, force_pbc,    &
                                              virial_cell, virial, &
                                              eelec, evdw)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: nonb_ene
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3,nthread)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw(nthread)

    ! local
    real(dp)    :: ene_virial(1:5)

    call timer(TimerNonBond, TimerOn)
 
    ! Calculate PME real part
    !
    if (real_calc) then

      call timer(TimerPmeReal, TimerOn)

#ifndef USE_GPU
      call compute_force_nonbond_table_linear( &
                                  domain, enefunc, pairlist, &
                                  coord_pbc, force_pbc, virial_cell)

      call timer(TimerPmeReal, TimerOff)
#else
      call compute_force_nonbond_table_linear_gpu( &
                                domain, enefunc, pairlist,  &
                                npt, .false., coord_pbc, &
                                force, force_pbc, virial,   &
                                ene_virial)
#endif

    end if

#ifndef USE_GPU
    call timer(TimerNonBond, TimerOff)
#endif
 
    return

  end subroutine compute_energy_nonbond_pme_short

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_pme_long
  !> @brief        Calculate nonbond energy (reciprocal part) by PME
  !! @authors      JJ
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    boundary : boundary information
  !! @param[in]    npt      : flag for NPT or not
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial term of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_pme_long(domain, enefunc, boundary, &
                                             npt, force, virial, eelec)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_boundary),        intent(in)    :: boundary
    logical,                 intent(in)    :: npt
    real(dp),                intent(inout) :: force(:,:,:)
    real(dp),                intent(inout) :: virial(3,3,nthread)
    real(dp),                intent(inout) :: eelec(nthread)


    call timer(TimerNonBond, TimerOn)

    ! Calculate PME reciprocal part
    !
    if (reciprocal_calc) then

      call timer(TimerPmeRecip, TimerOn)

      if (npt) call pme_pre(domain, boundary)

      call pme_recip(domain, force, virial, eelec)

      call timer(TimerPmeRecip, TimerOff)

    end if

    ! Add self energy
    !
    eelec(1) = eelec(1) + u_self

    call timer(TimerNonBond, TimerOff)

    return

  end subroutine compute_energy_nonbond_pme_long

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_cutoff_fep
  !> @brief        compute nonbond energy with cutoff for FEP
  !! @authors      HO
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    pairlist : pair-list information
  !! @param[in]    nonb_ene : flag for calculate nonbonded energy
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : lennard-jones energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_cutoff_fep(domain, enefunc, pairlist, &
      nonb_ene, force, virial, eelec, evdw)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    logical,                 intent(in)    :: nonb_ene
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: virial(3,3,nthread)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw(nthread)

    call timer(TimerNonBond, TimerOn)

    if (enefunc%vacuum) then
      ! Without lookup table

      call compute_energy_nonbond14_notable_fep( &
        domain, enefunc, &
        force, virial, eelec, evdw)

      if (nonb_ene) then

        call compute_energy_nonbond_notable_fep( &
          domain, enefunc, pairlist, &
          force, virial, eelec, evdw)

      else

        call compute_force_nonbond_notable_fep( &
          domain, enefunc, pairlist, &
          force, virial)

      end if

    else

      call compute_energy_nonbond14_table_fep( &
        domain, enefunc, &
        force, virial, eelec, evdw)
   
      ! ==> Type 4 and 10
      if (nonb_ene) then

        call compute_energy_nonbond_table_fep( &
          domain, enefunc, pairlist, &
          force, virial, eelec, evdw)
      
      else

        call compute_force_nonbond_table_fep( &
          domain, enefunc, pairlist, &
          force, virial)
      
      end if

    end if

    call timer(TimerNonBond, TimerOff)
 
    return

  end subroutine compute_energy_nonbond_cutoff_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_pme_fep
  !> @brief        Calculate nonbond energy by PME for FEP calculation
  !! @authors      NK
  !! @param[in]    domain       : domain information
  !! @param[in]    enefunc      : potential energy functions information
  !! @param[in]    pairlist     : pair-list information
  !! @param[in]    boundary     : boundary information
  !! @param[in]    npt          : flag for NPT or not
  !! @param[in]    nonb_ene     : flag for calculate nonbonded energy
  !! @param[in]    nonb_limiter : flag for contact checker
  !! @param[inout] coord_pbc    : !TODO
  !! @param[inout] force        : forces of target systems
  !! @param[inout] force_pbc    : !TODO
  !! @param[inout] virial_cell  : !TODO
  !! @param[inout] virial       : virial term of target systems
  !! @param[inout] eelec        : electrostatic energy of target systems
  !! @param[inout] evdw         : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_pme_fep(domain, enefunc, pairlist, boundary, &
                                        npt, nonb_ene, nonb_limiter,    &
                                        coord_pbc,                           &
                                        force_long, force, force_pbc,  &
                                        virial_cell, virial,      &
                                        eelec, evdw)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    type(s_boundary),        intent(in)    :: boundary
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: nonb_ene
    logical,                 intent(in)    :: nonb_limiter
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(dp),                intent(inout) :: force_long(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3,nthread)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw(nthread)

    ! local
    real(dp) :: ene_virial(1:5)
    ! FEP
    integer  :: i
    logical  :: flag_pme
    integer  :: nstate

    nstate = enefunc%num_fep_neighbor + 1

    call timer(TimerNonBond, TimerOn)


    ! Calculate PME real part
    !
    if (real_calc) then

      call timer(TimerPmeReal, TimerOn)

#ifndef USE_GPU
      if (nonb_limiter) then

        call compute_energy_nonbond_table_linear_check_fep( &
          domain, enefunc, pairlist, &
          coord_pbc, force_pbc, virial_cell, &
          eelec, evdw)

      else

        if (nonb_ene) then

          call compute_energy_nonbond_table_linear_fep( &
            domain, enefunc, pairlist, &
            coord_pbc, force_pbc, virial_cell, &
            eelec, evdw)

        else

          call compute_force_nonbond_table_linear_fep( &
            domain, enefunc, pairlist, &
            coord_pbc, force_pbc, virial_cell)

        end if

      end if

#else
      if (nonb_ene) then
        call compute_energy_nonbond_table_linear_gpu_fep( &
                                  domain, enefunc, pairlist, &
                                  npt, coord_pbc, force_pbc, virial, &
                                  eelec, evdw, ene_virial)
      else
        call compute_force_nonbond_table_linear_gpu_fep( &
                                  domain, enefunc, pairlist, npt, &
                                  .false., coord_pbc, force,      &
                                  force_pbc, virial, ene_virial)
      end if
#endif

#ifndef USE_GPU
      call timer(TimerPmeReal, TimerOff)
#endif

    end if


    ! Calculate PME reciprocal part
    !
    if (reciprocal_calc) then

      call timer(TimerPmeRecip, TimerOn)

      if (npt) call pme_pre_fep(domain, boundary)

      ! PME reciprocal for PRESERVE
      !
      if ((enefunc%lambelA + enefunc%lambelB) /= 1.0) &
        call pme_recip_fep(domain, enefunc, FEP_PRESERVE, force_long, &
                           virial, eelec)

      ! PME reciprocal for APPEAR
      !
      if (enefunc%lambelB /= 0.0) &
        call pme_recip_fep(domain, enefunc, FEP_APPEAR, force_long, virial, &
                         eelec)

      ! PME reciprocal for VANISH
      !
      if (enefunc%lambelA /= 0.0) &
        call pme_recip_fep(domain, enefunc, FEP_VANISH, force_long, virial, &
                         eelec)

      call timer(TimerPmeRecip, TimerOff)

   end if

   ! Add self energy
   !
   eelec(1) = eelec(1) &
     + (1.0_wp - enefunc%lambelB - enefunc%lambelA)*u_self_preserve &
     + enefunc%lambelB*u_self_appear &
     + enefunc%lambelA*u_self_vanish

#ifndef USE_GPU
    call timer(TimerNonBond, TimerOff)
#endif

    return

  end subroutine compute_energy_nonbond_pme_fep

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond_pme_short_fep
  !> @brief        Calculate nonbond energy (real part) by PME for FEP
  !! @authors      HO
  !! @param[in]    domain      : domain information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[in]    pairlist    : pair-list information
  !! @param[in]    boundary    : boundary information
  !! @param[inout] coord_pbc   : !TODO
  !! @param[inout] force       : forces of target systems
  !! @param[inout] force_pbc   : !TODO
  !! @param[inout] virial_cell : !TODO
  !! @param[inout] virial      : virial term of target systems
  !! @param[inout] eelec       : electrostatic energy of target systems
  !! @param[inout] evdw        : van der Waals energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_pme_short_fep(domain, enefunc, pairlist, &
                                              npt, nonb_ene, coord_pbc, &
                                              force, force_pbc,    &
                                              virial_cell, virial, &
                                              eelec, evdw)

    ! formal arguments
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_pairlist),        intent(in)    :: pairlist
    logical,                 intent(in)    :: npt
    logical,                 intent(in)    :: nonb_ene
    real(wp),                intent(inout) :: coord_pbc(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(wp),                intent(inout) :: force_pbc(:,:,:,:)
    real(dp),                intent(inout) :: virial_cell(:,:)
    real(dp),                intent(inout) :: virial(3,3,nthread)
    real(dp),                intent(inout) :: eelec(nthread)
    real(dp),                intent(inout) :: evdw(nthread)

    ! local
    real(dp)    :: ene_virial(1:5)

    call timer(TimerNonBond, TimerOn)
 
    ! Calculate PME real part
    !
    if (real_calc) then

      call timer(TimerPmeReal, TimerOn)

#ifndef USE_GPU
      call compute_force_nonbond_table_linear_fep( &
        domain, enefunc, pairlist, &
        coord_pbc, force_pbc, virial_cell)

      call timer(TimerPmeReal, TimerOff)
#else
      call compute_force_nonbond_table_linear_gpu_fep( &
                                domain, enefunc, pairlist,  &
                                npt, .false., coord_pbc, &
                                force, force_pbc, virial,   &
                                ene_virial)
#endif

    end if

#ifndef USE_GPU
    call timer(TimerNonBond, TimerOff)
#endif
 
    return

  end subroutine compute_energy_nonbond_pme_short_fep

end module sp_energy_nonbonds_mod
