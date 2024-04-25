!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_energy_contacts_mod
!> @brief   calculate native contact energy
!! @authors Jaewoon Jung (JJ)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_energy_go_mod

  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use timers_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public :: compute_energy_contact_126
  public :: compute_energy_noncontact_nobc

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_contact_126
  !> @brief        calculate contact energy
  !! @authors      JJ
  !! @param[in]    domain  : domain information
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[inout] force   : forces of target systems
  !! @param[inout] econtact: contact energy of target systems
  !! @param[inout] enoncontact: non-contact energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_contact_126(domain, enefunc, coord, force, &
                                        econtact, enoncontact)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(in)    :: enefunc
    real(dp),                intent(in)    :: coord(:,:,:)
    real(wp),                intent(inout) :: force(:,:,:,:)
    real(dp),                intent(inout) :: econtact(nthread)
    real(dp),                intent(inout) :: enoncontact(nthread)

    ! local variables
    real(wp)                 :: dij(1:3), rij2, inv_rij2, inv_rij6, inv_rij12
    real(wp)                 :: lj6, lj12, term_lj12, term_lj6
    real(wp)                 :: econt, encont, work(3), coef, cutoff2
    integer                  :: list(2)
    integer                  :: i, j, ix, icel1, icel2, i1, i2
    integer                  :: iatmcls, jatmcls
    integer                  :: id, omp_get_thread_num

    real(wp),        pointer :: contact_lj12(:,:), contact_lj6(:,:)
    real(wp),        pointer :: nonb_lj12(:,:)
    integer,         pointer :: ncontact(:), contactlist(:,:,:)
    integer,         pointer :: atmcls(:,:)
    integer,         pointer :: ncell_local
    integer,         pointer :: id_g2l(:,:)


    ncell_local => domain%num_cell_local
    id_g2l      => domain%id_g2l
    atmcls      => domain%atom_cls_no

    ncontact    => enefunc%num_contact
    contactlist => enefunc%contact_list
    contact_lj12=> enefunc%contact_lj12
    contact_lj6 => enefunc%contact_lj6
    nonb_lj12   => enefunc%nonb_lj12

    cutoff2     = enefunc%cutoffdist * enefunc%cutoffdist

    ! calculate bond energy
    !
    !$omp parallel default(shared) &
    !$omp private(id, i, ix, icel1, i1, icel2, i2, lj6, lj12, dij, rij2, &
    !$omp         inv_rij2, inv_rij6, inv_rij12, term_lj12, term_lj6,    &
    !$omp         econt, encont, iatmcls, jatmcls, coef, work, list)
    !
#ifdef OMP
    id  = omp_get_thread_num()
#else
    id  = 0
#endif

    do i = id+1, ncell_local, nthread

      econt  = 0.0_wp
      encont = 0.0_wp

      do ix = 1, ncontact(i)

        list(1:2) = contactlist(1:2,ix,i)
        icel1 = id_g2l(1,list(1))
        i1    = id_g2l(2,list(1))
        icel2 = id_g2l(1,list(2))
        i2    = id_g2l(2,list(2))
        dij(1:3)  = coord(1:3,i1,icel1) - coord(1:3,i2,icel2)

        ! contact energy
        !
        lj6       = contact_lj6(ix,i)
        lj12      = contact_lj12(ix,i)
        rij2      = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        inv_rij2  = 1.0_wp / rij2
        inv_rij6  = inv_rij2 * inv_rij2 * inv_rij2
        inv_rij12 = inv_rij6 * inv_rij6

        term_lj12 = lj12 * inv_rij12
        term_lj6  = lj6 * inv_rij6
        econt     = econt + term_lj12 - term_lj6
        coef      = 12.0_wp*term_lj12 - 6.0_wp*term_lj6

        ! noncontact energy
        !
        if (rij2 < cutoff2) then
          iatmcls   = atmcls(i1,icel1)
          jatmcls   = atmcls(i2,icel2)
          lj12      = nonb_lj12(iatmcls, jatmcls)
          term_lj12 = lj12 * inv_rij12
          encont    = encont - term_lj12
          coef      = coef - 12.0_wp*term_lj12
        end if

        ! gradient: dE/dX
        !
        coef      = - inv_rij2 * coef
        work(1:3) = coef * dij(1:3)

        ! store force: F=-dE/dX
        !
        force(1:3,i1,icel1,id+1) = force(1:3,i1,icel1,id+1) - work(1:3)
        force(1:3,i2,icel2,id+1) = force(1:3,i2,icel2,id+1) + work(1:3)

      end do

      econtact(id+1) = econtact(id+1) + econt
      enoncontact(id+1) = enoncontact(id+1) + encont

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_contact_126

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_noncontact_nobc
  !> @brief        calculate nonbonded energy
  !! @authors      JJ
  !! @param[in]    domain   : domain information
  !! @param[in]    enefunc  : potential energy functions
  !! @param[in]    pairlist : interaction list in each domain
  !! @param[inout] force    : forces for each cell
  !! @param[inout] enoncontact : non-native contact energy
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_noncontact_nobc(domain, enefunc, pairlist, &
                                            coord, force, enoncontact)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(dp),                 intent(in)    :: coord(:,:,:)
    real(wp),                 intent(inout) :: force(:,:,:,:)
    real(dp),                 intent(inout) :: enoncontact(nthread)

    ! local variables
    real(wp)                  :: dij(1:3), rij2, inv_rij2, inv_rij6, inv_rij12
    real(wp)                  :: lj12, term_lj12, cutoff2
    real(wp)                  :: grad_coef, work(1:3)
    real(wp)                  :: rtmp(1:3)
    real(wp)                  :: encont
    real(wp)                  :: force_local(3)
    integer                   :: i, ix, iy, j, k
    integer                   :: num_nb15
    integer                   :: id, omp_get_thread_num, j_list(MaxAtom)
    integer                   :: ncell, ncell_local
    integer                   :: iatmcls, jatmcls

    real(wp),         pointer :: nonb_lj12(:,:)
    real(wp),         pointer :: cutoff
    real(wp),         pointer :: table_ene(:), table_grad(:)
    integer,          pointer :: cell_pairlist(:,:), atmcls(:,:)
    integer,          pointer :: natom(:), nsolute(:), nwater(:)
    integer,          pointer :: nb15_list(:,:)
    integer,          pointer :: nb15_cell(:)
    integer,          pointer :: num_nb15_nobc(:,:)
    integer,          pointer :: nb15_calc_list(:,:,:,:)


    atmcls          => domain%atom_cls_no
    natom           => domain%num_atom

    cutoff          => enefunc%cutoffdist
    nonb_lj12       => enefunc%nonb_lj12

    num_nb15_nobc   => pairlist%num_nb15_nobc
    nb15_calc_list  => pairlist%nb15_calc_list_nobc

    ncell           =  domain%num_cell_local + domain%num_cell_boundary
    cutoff2         =  cutoff * cutoff

    !$omp parallel default(shared)                             &
    !$omp private(id, i, ix, j, iy, k, rtmp, num_nb15, encont, &
    !$omp         dij, rij2, inv_rij2, inv_rij6, inv_rij12,    &
    !$omp         term_lj12, grad_coef, work, force_local,     &
    !$omp         iatmcls, jatmcls, lj12)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do i = id+1, ncell, nthread

      do ix = 1, natom(i)

        rtmp(1:3) = coord(1:3,ix,i)

        num_nb15 = num_nb15_nobc(ix,i)
        iatmcls  = atmcls(ix,i)
        encont   = 0.0_wp
        force_local(1:3) = 0.0_wp

        do k = 1, num_nb15

          j  = nb15_calc_list(1,k,ix,i)
          iy = nb15_calc_list(2,k,ix,i)

          ! compute distance
          !
          dij(1:3) = rtmp(1:3) - coord(1:3,iy,j)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)

          if (rij2 < cutoff2) then

            jatmcls = atmcls(iy,j)
            lj12 = nonb_lj12(iatmcls, jatmcls)

            inv_rij2  = 1.0_wp / rij2
            inv_rij6  = inv_rij2 * inv_rij2 * inv_rij2
            inv_rij12 = inv_rij6 * inv_rij6
            term_lj12 = lj12 * inv_rij12

            encont = encont + term_lj12

            grad_coef   = - inv_rij2 * 12.0_wp * term_lj12
            work(1:3) = grad_coef * dij(1:3)

            ! store force
            !
            force_local(1:3) = force_local(1:3) - work(1:3)
            force(1:3,iy,j,id+1)  = force(1:3,iy,j,id+1) + work(1:3)

          end if

        end do

        force(1:3,ix,i,id+1) = force(1:3,ix,i,id+1) + force_local(1:3)
        enoncontact(id+1) = enoncontact(id+1) + encont

      end do

    end do

    !$omp end parallel

    return

  end subroutine compute_energy_noncontact_nobc

end module sp_energy_go_mod

