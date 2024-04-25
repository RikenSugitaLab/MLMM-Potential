!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_enefunc_charmm_mod
!> @brief   define potential energy functions
!! @authors Jaewoon Jung (JJ), Takeshi Imai (TI), Chigusa Kobayashi (CK),
!!          Takaharu Mori (TM), Yuji Sugita (YS)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_enefunc_charmm_mod

  use sp_enefunc_localres_mod
  use sp_enefunc_restraints_mod
  use sp_enefunc_fit_mod
  use sp_enefunc_table_mod
  use sp_energy_mod
  use sp_restraints_str_mod
  use sp_constraints_str_mod
  use sp_enefunc_str_mod
  use sp_energy_str_mod
  use sp_domain_str_mod
  use dihedral_libs_mod
  use molecules_str_mod
  use fileio_par_mod
  use fileio_localres_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use fileio_pdb_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: define_enefunc_charmm
  private :: setup_enefunc_bond
  private :: setup_enefunc_bond_constraint
  private :: setup_enefunc_angl
  private :: setup_enefunc_angl_constraint
  private :: setup_enefunc_dihe
  private :: setup_enefunc_impr
  private :: setup_enefunc_cmap
  private :: setup_enefunc_nonb
  public  :: count_nonb_excl
  public  :: count_nonb_excl_rest
  public  :: count_nonb_excl_solute
  public  :: count_nonb_excl_solute_rest
  public  :: count_nonb_excl_constraint
  public  :: count_nonb_excl_constraint_rest
  public  :: pack_array_4to3
  public  :: pack_array_3to2

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_enefunc_charmm
  !> @brief        a driver subroutine for defining potential energy functions
  !! @authors      YS, TI, JJ, CK
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[in]    par         : CHARMM PAR information
  !! @param[in]    localres    : local restraint information
  !! @param[in]    molecule    : molecule information
  !! @param[inout] constraints : constraints information
  !! @param[in]    restraints  : restraints information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_enefunc_charmm(ene_info, par, localres, molecule, &
                                   constraints, restraints, domain, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_par),             intent(in)    :: par
    type(s_localres),        intent(in)    :: localres
    type(s_molecule),        intent(in)    :: molecule
    type(s_constraints),     intent(inout) :: constraints
    type(s_restraints),      intent(in)    :: restraints
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: ncel, ncelb


    ! base
    !
    ncel  = domain%num_cell_local
    ncelb = domain%num_cell_local + domain%num_cell_boundary

    call alloc_enefunc(enefunc, EneFuncBase, ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncBond, ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncAngl, ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncDihe, ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncImpr, ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncBondCell, ncel, ncelb)

    if (.not. constraints%rigid_bond) then

      ! bond
      !
      call setup_enefunc_bond(par, molecule, domain, enefunc)

      ! angle
      !
      call setup_enefunc_angl(par, molecule, domain, enefunc)

    else

      ! bond
      !
      call setup_enefunc_bond_constraint( &
                              par, molecule, domain, constraints, enefunc)

      ! angle
      !
      call setup_enefunc_angl_constraint( &
                              par, molecule, domain, constraints, enefunc)

    end if

    ! dihedral
    !
    call setup_enefunc_dihe(par, molecule, domain, enefunc)

    ! improper
    !
    call setup_enefunc_impr(par, molecule, domain, enefunc)

    ! cmap
    !
    call setup_enefunc_cmap(ene_info, par, molecule, domain, enefunc)

    ! nonbonded
    !
    call setup_enefunc_nonb(par, molecule, constraints, domain, enefunc)

    ! lookup table
    !
    if (ene_info%table) then
      if (.not. enefunc%vacuum) then
        call setup_enefunc_table(ene_info, enefunc)
      end if
    end if

    ! restraint
    !
    call setup_enefunc_restraints(molecule, restraints, domain, enefunc)

    call setup_enefunc_localres(localres, domain, enefunc)

    ! write summary of energy function
    !
    if (main_rank) then
      write(MsgOut,'(A)') &
           'Define_Enefunc_Charmm> Number of Interactions in Each Term'
      write(MsgOut,'(A20,I10,A20,I10)')                  &
           '  bond_ene        = ', enefunc%num_bond_all, &
           '  angle_ene       = ', enefunc%num_angl_all
      write(MsgOut,'(A20,I10,A20,I10)')                  &
           '  torsion_ene     = ', enefunc%num_dihe_all, &
           '  improper_ene    = ', enefunc%num_impr_all
      write(MsgOut,'(A20,I10)')                          &
           '  cmap_ene        = ', enefunc%num_cmap_all
      write(MsgOut,'(A20,I10,A20,I10)')                  &
           '  nb_exclusions   = ', enefunc%num_excl_all, &
           '  nb14_calc       = ', enefunc%num_nb14_all
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine define_enefunc_charmm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_bond
  !> @brief        define BOND term for each cell in potential energy function
  !! @authors      YS, JJ, TM
  !! @param[in]    par      : CHARMM PAR information
  !! @param[in]    molecule : molecule information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_bond(par, molecule, domain, enefunc)

    ! formal arguments
    type(s_par),             intent(in)    :: par
    type(s_molecule),        intent(in)    :: molecule
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variable
    integer                  :: i, j, icel_local
    integer                  :: icel1, icel2
    integer                  :: nbond, nbond_p, i1, i2, found
    character(6)             :: ci1, ci2

    real(wp),        pointer :: force(:,:), dist(:,:)
    integer,         pointer :: bond(:), list(:,:,:)
    integer,         pointer :: ncel, cell_pair(:,:)
    integer,         pointer :: id_g2l(:,:)
    integer                  :: nbond_fep, ii1, ii2

    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l

    bond      => enefunc%num_bond
    list      => enefunc%bond_list
    force     => enefunc%bond_force_const
    dist      => enefunc%bond_dist_min

    nbond     = molecule%num_bonds
    nbond_p   = par%num_bonds

    do i = 1, nbond

      ci1 = molecule%atom_cls_name(molecule%bond_list(1,i))
      ci2 = molecule%atom_cls_name(molecule%bond_list(2,i))
      i1  = molecule%bond_list(1,i)
      i2  = molecule%bond_list(2,i)

      if (domain%fep_use) then
        ! FEP: If the bond are not set to any group of FEP, exclude this bond.
        ii1 = molecule%fepgrp(i1)
        ii2 = molecule%fepgrp(i2)
        if (molecule%fepgrp_bond(ii1,ii2) == 0) cycle
      end if

      icel1 = id_g2l(1,i1)
      icel2 = id_g2l(1,i2)

      ! Check if it is in my domain
      !
      if (icel1 /= 0 .and. icel2 /= 0) then

        icel_local = cell_pair(icel1,icel2)

        if (icel_local > 0 .and. icel_local <= ncel) then

          do j = 1, nbond_p
            if ((ci1 == par%bond_atom_cls(1,j) .and.  &
                 ci2 == par%bond_atom_cls(2,j)) .or.  &
                (ci1 == par%bond_atom_cls(2,j) .and.  &
                 ci2 == par%bond_atom_cls(1,j))) then

              bond (icel_local) = bond(icel_local) + 1
              list (1,bond(icel_local),icel_local) = i1
              list (2,bond(icel_local),icel_local) = i2
              force(bond(icel_local),icel_local)   = par%bond_force_const(j)
              dist (bond(icel_local),icel_local)   = par%bond_dist_min(j)
              exit

            end if
          end do

          if (j == nbond_p + 1) &
            write(MsgOut,*) &
              'Setup_Enefunc_Bond> not found BOND: [', &
              ci1, ']-[', ci2, '] in parameter file. (ERROR)'

        end if

      end if

    end do

    found = 0
    do i = 1, ncel
      found = found + bond(i)
      if (bond(i) > MaxBond) &
        call error_msg('Setup_Enefunc_Bond> Too many bonds.')
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_bond_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_bond_all = found
#endif

    if (domain%fep_use) then
      nbond_fep = 0
      do i = 1, 5
        nbond_fep = nbond_fep + molecule%num_bonds_fep(i)
      end do
      if (enefunc%num_bond_all /= (nbond_fep)) &
        call error_msg('Setup_Enefunc_Bond> Some bond paremeters are missing.')
    else
      if (enefunc%num_bond_all /= nbond) &
        call error_msg('Setup_Enefunc_Bond> Some bond paremeters are missing.')
    end if

    return

  end subroutine setup_enefunc_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_bond_constraint
  !> @brief        define BOND term between heavy atoms
  !! @authors      JJ
  !! @param[in]    par         : CHARMM PAR information
  !! @param[in]    molecule    : molecule information
  !! @param[in]    domain      : domain information
  !! @param[inout] constraints : constraints information
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_bond_constraint(par, molecule, domain, &
                                           constraints, enefunc)

    ! formal arguments
    type(s_par),                 intent(in)    :: par
    type(s_molecule),            intent(in)    :: molecule
    type(s_domain),      target, intent(in)    :: domain
    type(s_constraints), target, intent(inout) :: constraints
    type(s_enefunc),     target, intent(inout) :: enefunc

    ! local variable
    integer                      :: i, j, k, l, m, ih, icel_local, connect
    integer                      :: i1, i2, ih1, ih2, icel1, icel2, icel
    integer                      :: nbond, nbond_p, nbond_a, nbond_c
    integer                      :: wat_bonds, tmp_mole_no, mole_no, wat_found
    character(6)                 :: ci1, ci2
    character(6)                 :: ri1, ri2
    character(4)                 :: atom1, atom2
    logical                      :: mi1, mi2
    logical                      :: cl1, cl2

    real(wp),            pointer :: force(:,:), dist(:,:)
    real(dp),            pointer :: HGr_bond_dist(:,:,:,:)
    integer,             pointer :: bond(:), list(:,:,:), num_water, ncel
    integer,             pointer :: cell_pair(:,:), id_g2l(:,:), id_l2g(:,:)
    integer,             pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)
    integer                      :: nbond_fep, ii1, ii2

    ncel          => domain%num_cell_local
    cell_pair     => domain%cell_pair
    id_g2l        => domain%id_g2l
    id_l2g        => domain%id_l2g

    HGr_local     => constraints%HGr_local
    HGr_bond_list => constraints%HGr_bond_list
    HGr_bond_dist => constraints%HGr_bond_dist

    bond          => enefunc%num_bond
    list          => enefunc%bond_list
    force         => enefunc%bond_force_const
    dist          => enefunc%bond_dist_min
    num_water     => enefunc%table%num_water

    nbond         = molecule%num_bonds
    nbond_p       = par%num_bonds

    connect       =  constraints%connect

    nbond_a       = 0
    nbond_c       = 0

    do i = 1, nbond

      i1  = molecule%bond_list(1,i)
      i2  = molecule%bond_list(2,i)

      if (domain%fep_use) then
        ! FEP: If the bond are not set to any group of FEP, exclude this bond.
        ii1 = molecule%fepgrp(i1)
        ii2 = molecule%fepgrp(i2)
        if (molecule%fepgrp_bond(ii1,ii2) == 0) cycle
      end if

      ci1 = molecule%atom_cls_name(i1)
      ci2 = molecule%atom_cls_name(i2)
      mi1 = molecule%light_atom_mass(i1)
      mi2 = molecule%light_atom_mass(i2)
      cl1 = molecule%light_atom_name(i1)
      cl2 = molecule%light_atom_name(i2)
      if (constraints%hydrogen_type == ConstraintAtomMass) then
        cl1 = mi1
        cl2 = mi2
      else if (constraints%hydrogen_type == ConstraintAtomBoth) then
        cl1 = (cl1 .or. mi1)
        cl2 = (cl2 .or. mi2)
      endif

      if (.not. (cl1 .or.  cl2)) then
        icel1 = id_g2l(1,i1)
        icel2 = id_g2l(1,i2)

        ! Check if it is in my domain
        !
        if (icel1 /= 0 .and. icel2 /= 0) then

          icel_local = cell_pair(icel1,icel2)

          if (icel_local > 0 .and. icel_local <= ncel) then

            do j = 1, nbond_p
              if ((ci1 == par%bond_atom_cls(1, j) .and.  &
                   ci2 == par%bond_atom_cls(2, j)) .or.  &
                  (ci1 == par%bond_atom_cls(2, j) .and.  &
                   ci2 == par%bond_atom_cls(1, j))) then

                nbond_a = nbond_a + 1
                bond (icel_local) = bond(icel_local) + 1
                list (1,bond(icel_local),icel_local) = i1
                list (2,bond(icel_local),icel_local) = i2
                force(bond(icel_local),icel_local) = par%bond_force_const(j)
                dist (bond(icel_local),icel_local) = par%bond_dist_min(j)
                exit

              end if
            end do

            if (j == nbond_p + 1) &
              write(MsgOut,*) &
                'Setup_Enefunc_Bond_Constraint> not found BOND: [', &
                ci1, ']-[', ci2, '] in parameter file. (ERROR)'

          end if

        end if

      else

        icel1 = id_g2l(1,i1)
        icel2 = id_g2l(1,i2)

        if (icel1 /= 0 .and. icel2 /= 0) then

          icel = cell_pair(icel1,icel2)
          if (icel > 0 .and. icel <= ncel) then

            do j = 1, connect

              do k = 1, HGr_local(j,icel)
                ih1 = id_l2g(HGr_bond_list(1,k,j,icel),icel)
                do ih = 1, j
                  ih2 = id_l2g(HGr_bond_list(ih+1,k,j,icel),icel)

                  if (domain%fep_use) then
                    ! FEP: Hydrogen rewiring
                    ! In singleA-dualB bonds, one end has atom index of singleA, 
                    ! and so par%bond_dist_min becomes singleA-dualB value.
                    ! However, in real, the bonds should be considered as
                    ! singleB-dualB. The parameter of singleB-dualB should be
                    ! used. To avoid this problem, for singleA-dualB bond
                    ! including hydrogen, the atom index of singleA is replaced
                    ! with the corresponding atom index of singleB.
                    if ((int(molecule%fepgrp(ih1)) == 1) .and. &
                        (int(molecule%fepgrp(ih2)) == 4)) then
                      do l = 1, molecule%num_atoms_fep(1)
                        if (molecule%id_singleA(l) == ih1) then
                          ih1 = molecule%id_singleB(l)
                          exit
                        end if
                      end do
                    else if ((int(molecule%fepgrp(ih1)) == 4) .and. &
                             (int(molecule%fepgrp(ih2)) == 1)) then
                      do l = 1, molecule%num_atoms_fep(2)
                        if (molecule%id_singleB(l) == ih2) then
                          ih2 = molecule%id_singleA(l)
                          exit
                        end if
                      end do
                    end if
                  end if

                  if (ih1 == i1 .and. ih2 == i2 .or. &
                      ih2 == i1 .and. ih1 == i2) then

                    do m = 1, nbond_p
                      if ((ci1 == par%bond_atom_cls(1, m) .and.  &
                           ci2 == par%bond_atom_cls(2, m)) .or.  &
                          (ci1 == par%bond_atom_cls(2, m) .and.  &
                           ci2 == par%bond_atom_cls(1, m))) then

                        nbond_c = nbond_c + 1
                        HGr_bond_dist(ih+1,k,j,icel) = par%bond_dist_min(m)
                        exit
                      end if
                    end do

                  end if

                end do
              end do

            end do
          end if
        end if

      end if

    end do


    ! for water molecule
    !
    wat_bonds = 0
    wat_found = 0

    if (constraints%fast_water) then

      tmp_mole_no = -1

      do i = 1, nbond

        ri1 = molecule%residue_name(molecule%bond_list(1,i))
        ri2 = molecule%residue_name(molecule%bond_list(2,i))

        if (ri1 == constraints%water_model .and. &
            ri2 == constraints%water_model) then

          wat_bonds=wat_bonds+1
          mole_no = molecule%molecule_no(molecule%bond_list(1,i))

          if (mole_no /= tmp_mole_no) then
            wat_found = wat_found +1
            tmp_mole_no = mole_no
          end if

        end if

      end do

      if (wat_found /= num_water) &
        call error_msg( &
          'Setup_Enefunc_Bond_Constraint> # of water is incorrect')

    end if

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(nbond_a, enefunc%num_bond_all,  1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)

    call mpi_allreduce(nbond_c, constraints%num_bonds, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_bond_all  = nbond_a
    constraints%num_bonds = nbond_c
#endif

    if (domain%fep_use) then
      nbond_fep = 0
      do i = 1, 5
        nbond_fep = nbond_fep + molecule%num_bonds_fep(i)
      end do
      if (enefunc%num_bond_all /= (nbond_fep-constraints%num_bonds-wat_bonds)) then
        call error_msg( &
          'Setup_Enefunc_Bond_Constraint> Some bond paremeters are missing.')
      end if
    else
      if (enefunc%num_bond_all /= (nbond-constraints%num_bonds-wat_bonds)) then
        call error_msg( &
          'Setup_Enefunc_Bond_Constraint> Some bond paremeters are missing.')
      end if
    end if

    return

  end subroutine setup_enefunc_bond_constraint

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_angl
  !> @brief        define ANGLE term for each cell in potential energy function
  !! @authors      YS, JJ, TM
  !! @param[in]    par      : CHARMM PAR information
  !! @param[in]    molecule : molecule information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_angl(par, molecule, domain, enefunc)

    ! formal arguments
    type(s_par),     target, intent(in)    :: par
    type(s_molecule),        intent(in)    :: molecule
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, icel_local
    integer                  :: icel1, icel2
    integer                  :: nangl, nangl_p, found
    integer                  :: list(3)
    character(6)             :: ci1, ci2, ci3

    real(wp),        pointer :: force(:,:), theta(:,:)
    real(wp),        pointer :: ubforce(:,:), ubrmin(:,:)
    integer,         pointer :: angle(:), alist(:,:,:)
    integer,         pointer :: ncel, cell_pair(:,:), id_g2l(:,:)
    integer                  :: nangl_fep, i1, i2, i3

    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l

    angle     => enefunc%num_angle
    alist     => enefunc%angle_list
    force     => enefunc%angle_force_const
    theta     => enefunc%angle_theta_min
    ubforce   => enefunc%urey_force_const
    ubrmin    => enefunc%urey_rmin

    nangl     = molecule%num_angles
    nangl_p   = par%num_angles

    do i = 1, nangl

      ci1 = molecule%atom_cls_name(molecule%angl_list(1,i))
      ci2 = molecule%atom_cls_name(molecule%angl_list(2,i))
      ci3 = molecule%atom_cls_name(molecule%angl_list(3,i))

      list(1) = molecule%angl_list(1,i)
      list(2) = molecule%angl_list(2,i)
      list(3) = molecule%angl_list(3,i)

      icel1 = id_g2l(1,list(1))
      icel2 = id_g2l(1,list(3))

      if (domain%fep_use) then
        ! FEP: If the angle are not set to any group of FEP, exclude this angle.
        i1 = molecule%fepgrp(list(1))
        i2 = molecule%fepgrp(list(2))
        i3 = molecule%fepgrp(list(3))
        if (molecule%fepgrp_angl(i1,i2,i3) == 0) cycle
      end if

      if (icel1 /= 0 .and. icel2 /= 0) then

        icel_local = cell_pair(icel1,icel2)

        if (icel_local >= 1 .and. icel_local <= ncel) then

          do j = 1, nangl_p
            if ((ci1 == par%angl_atom_cls(1,j) .and. &
                 ci2 == par%angl_atom_cls(2,j) .and. &
                 ci3 == par%angl_atom_cls(3,j)) .or. &
                (ci1 == par%angl_atom_cls(3,j) .and. &
                 ci2 == par%angl_atom_cls(2,j) .and. &
                 ci3 == par%angl_atom_cls(1,j))) then

              angle(icel_local) = angle(icel_local) + 1
              alist(1:3,angle(icel_local),icel_local) = list(1:3)

              force(angle(icel_local),icel_local)   = &
                   par%angl_force_const(j)
              theta(angle(icel_local),icel_local)   = &
                   par%angl_theta_min(j)*RAD
              ubforce(angle(icel_local),icel_local) = &
                   par%angl_ub_force_const(j)
              ubrmin(angle(icel_local),icel_local)  = &
                   par%angl_ub_rmin(j)
              exit

            end if
          end do

          if (j == nangl_p + 1) &
            write(MsgOut,*) &
              'Setup_Enefunc_Angl> not found ANGL: [',&
              ci1, ']-[', ci2, ']-[', ci3, '] in parameter file. (ERROR)'

        end if

      end if

    end do

    found = 0
    do i = 1, ncel
      found = found + angle(i)
      if (angle(i) > MaxAngle) &
        call error_msg('Setup_Enefunc_Angl> Too many angles.')
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_angl_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_angl_all = found
#endif

    if (domain%fep_use) then
      nangl_fep = 0
      do i = 1, 5
        nangl_fep = nangl_fep + molecule%num_angles_fep(i)
      end do
      if (enefunc%num_angl_all /= nangl_fep) &
        call error_msg('Setup_Enefunc_Angl> Some angle paremeters are missing.')
    else
      if (enefunc%num_angl_all /= nangl) &
        call error_msg('Setup_Enefunc_Angl> Some angle paremeters are missing.')
    end if

    return

  end subroutine setup_enefunc_angl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_angl_constraint
  !> @brief        define ANGLE term for each cell in potential energy function
  !                with SETTLE constraint
  !! @authors      JJ
  !! @param[in]    par         : CHARMM PAR information
  !! @param[in]    molecule    : molecule information
  !! @param[in]    domain      : domain information
  !! @param[in]    constraints : constraints information
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_angl_constraint(par, molecule, domain, &
                                           constraints, enefunc)

    ! formal arguments
    type(s_par),                 intent(in)    :: par
    type(s_molecule),            intent(in)    :: molecule
    type(s_domain),      target, intent(in)    :: domain
    type(s_constraints), target, intent(in)    :: constraints
    type(s_enefunc),     target, intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, icel_local
    integer                  :: icel1, icel2
    integer                  :: nangl, nangl_p, found, nangl_per_water
    integer                  :: list(3)
    character(6)             :: ci1, ci2, ci3
    character(6)             :: ri1, ri2, ri3

    real(wp),        pointer :: force(:,:), theta(:,:)
    real(wp),        pointer :: ubforce(:,:), ubrmin(:,:)
    integer,         pointer :: angle(:), alist(:,:,:), num_water
    integer,         pointer :: ncel, cell_pair(:,:), id_g2l(:,:)
    integer                  :: nangl_fep, i1, i2, i3

    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l

    angle     => enefunc%num_angle
    alist     => enefunc%angle_list
    force     => enefunc%angle_force_const
    theta     => enefunc%angle_theta_min
    ubforce   => enefunc%urey_force_const
    ubrmin    => enefunc%urey_rmin
    num_water => enefunc%table%num_water

    nangl     = molecule%num_angles
    nangl_p   = par%num_angles

    nangl_per_water = 0

    if (num_water > 0) then

      do i = 1, nangl

        list(1:3) = molecule%angl_list(1:3,i)
        ri1 = molecule%residue_name(list(1))
        ri2 = molecule%residue_name(list(2))
        ri3 = molecule%residue_name(list(3))

        if (ri1 == constraints%water_model .and. &
            ri2 == constraints%water_model .and. &
            ri3 == constraints%water_model) then

          nangl_per_water = nangl_per_water + 1

        end if

      end do

      if (mod(nangl_per_water,num_water) /= 0) then
        write(MsgOut,*) &
             'Setup_Enefunc_Angl_Constraint> invalid ANGL count: ', &
             'number of angle terms in a water molecule is not integer.'
        call error_msg('Setup_Enefunc_Angl_Constraint> abort')
      end if
      nangl_per_water = nangl_per_water / num_water

    end if

    do i = 1, nangl

      ci1 = molecule%atom_cls_name(molecule%angl_list(1,i))
      ci2 = molecule%atom_cls_name(molecule%angl_list(2,i))
      ci3 = molecule%atom_cls_name(molecule%angl_list(3,i))
      ri1 = molecule%residue_name (molecule%angl_list(1,i))
      ri2 = molecule%residue_name (molecule%angl_list(2,i))
      ri3 = molecule%residue_name (molecule%angl_list(3,i))

      if (ri1(1:4) /= constraints%water_model .and. &
          ri2(1:4) /= constraints%water_model .and. &
          ri3(1:4) /= constraints%water_model) then

        list(1) = molecule%angl_list(1,i)
        list(2) = molecule%angl_list(2,i)
        list(3) = molecule%angl_list(3,i)

        icel1 = id_g2l(1,list(1))
        icel2 = id_g2l(1,list(3))

        if (domain%fep_use) then
          ! FEP: If the angle are not set to any group of FEP, exclude this angle.
          i1 = molecule%fepgrp(list(1))
          i2 = molecule%fepgrp(list(2))
          i3 = molecule%fepgrp(list(3))
          if (molecule%fepgrp_angl(i1,i2,i3) == 0) cycle
        end if

        if (icel1 /= 0 .and. icel2 /= 0) then

          icel_local = cell_pair(icel1,icel2)

          if (icel_local >= 1 .and. icel_local <= ncel) then

            do j = 1, nangl_p
              if ((ci1 == par%angl_atom_cls(1,j) .and. &
                   ci2 == par%angl_atom_cls(2,j) .and. &
                   ci3 == par%angl_atom_cls(3,j)) .or. &
                  (ci1 == par%angl_atom_cls(3,j) .and. &
                   ci2 == par%angl_atom_cls(2,j) .and. &
                   ci3 == par%angl_atom_cls(1,j))) then

                angle(icel_local) = angle(icel_local) + 1
                alist(1:3,angle(icel_local),icel_local) = list(1:3)

                force(angle(icel_local),icel_local)   = &
                     par%angl_force_const(j)
                theta(angle(icel_local),icel_local)   = &
                     par%angl_theta_min(j)*RAD
                ubforce(angle(icel_local),icel_local) = &
                     par%angl_ub_force_const(j)
                ubrmin(angle(icel_local),icel_local)  = &
                     par%angl_ub_rmin(j)
                exit

              end if
            end do

            if (j == nangl_p + 1) &
              write(MsgOut,*) &
                'Setup_Enefunc_Angl_Constraint> not found ANGL: [', &
                ci1, ']-[', ci2, ']-[', ci3, '] in parameter file. (ERROR)'

          end if

        end if

      end if

    end do

    found = 0
    do i = 1, ncel
      found = found + angle(i)
      if (angle(i) > MaxAngle) &
        call error_msg('Setup_Enefunc_Angl_Constraint> Too many angles.')
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_angl_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_angl_all = found
#endif

    if (domain%fep_use) then
      nangl_fep = 0
      do i = 1, 5
        nangl_fep = nangl_fep + molecule%num_angles_fep(i)
      end do
      if (enefunc%num_angl_all /= (nangl_fep - nangl_per_water*num_water)) &
        call error_msg('Setup_Enefunc_Angl> Some angle paremeters are missing.')
    else
      if (enefunc%num_angl_all /= (nangl - nangl_per_water*num_water)) &
        call error_msg( &
          'Setup_Enefunc_Angl_Constraint> Some angle paremeters are missing.')
    end if

    return

  end subroutine setup_enefunc_angl_constraint

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_dihe
  !> @brief        define DIHEDRAL term in potential energy function
  !! @authors      YS, JJ, TM
  !! @param[in]    par      : CHARMM PAR information
  !! @param[in]    molecule : molecule information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_dihe(par, molecule, domain, enefunc)

    ! formal arguments
    type(s_par),      target, intent(in)    :: par
    type(s_molecule), target, intent(in)    :: molecule
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(inout) :: enefunc

    ! local variables
    integer                   :: ndihe, ndihe_p
    integer                   :: i, j, icel_local
    integer                   :: icel1, icel2
    integer                   :: found, nw_found
    integer                   :: list(4)
    character(6)              :: ci1, ci2, ci3, ci4

    real(wp),         pointer :: force(:,:), phase(:,:)
    integer,          pointer :: dihedral(:), dlist(:,:,:), period(:,:)
    integer,          pointer :: ncel, cell_pair(:,:), id_g2l(:,:)
    logical,      allocatable :: no_wild(:)
    integer,          pointer :: notation
    integer                   :: ndihe_fep, i1, i2, i3, i4

    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l

    dihedral  => enefunc%num_dihedral
    dlist     => enefunc%dihe_list
    force     => enefunc%dihe_force_const
    period    => enefunc%dihe_periodicity
    phase     => enefunc%dihe_phase
    notation  => enefunc%notation_14types
    notation = 100

    ndihe     = molecule%num_dihedrals
    ndihe_p   = par%num_dihedrals

    ! check usage of wild card
    !
    allocate(no_wild(ndihe_p))

    do i = 1, ndihe_p

      if ((par%dihe_atom_cls(1,i) /= WildCard) .and. &
          (par%dihe_atom_cls(4,i) /= WildCard)) then
        ! A-B-C-D type
        no_wild(i) = .true.
      else
        ! X-B-C-D type
        no_wild(i) = .false.
      end if

    end do

    ! find number of interactions
    !
    do i = 1, ndihe

      ci1 = molecule%atom_cls_name(molecule%dihe_list(1,i))
      ci2 = molecule%atom_cls_name(molecule%dihe_list(2,i))
      ci3 = molecule%atom_cls_name(molecule%dihe_list(3,i))
      ci4 = molecule%atom_cls_name(molecule%dihe_list(4,i))

      list(1) = molecule%dihe_list(1,i)
      list(2) = molecule%dihe_list(2,i)
      list(3) = molecule%dihe_list(3,i)
      list(4) = molecule%dihe_list(4,i)

      icel1 = id_g2l(1,list(1))
      icel2 = id_g2l(1,list(4))
      
      if (domain%fep_use) then
        ! FEP: If the dihedral are not set to any group of FEP, exclude this dihedral.
        i1 = molecule%fepgrp(list(1))
        i2 = molecule%fepgrp(list(2))
        i3 = molecule%fepgrp(list(3))
        i4 = molecule%fepgrp(list(4))
        if (molecule%fepgrp_dihe(i1,i2,i3,i4) == 0) cycle
      end if

      if (icel1 /= 0 .and. icel2 /= 0) then

        icel_local = cell_pair(icel1,icel2)

        if (icel_local >= 1 .and. icel_local <= ncel) then

          nw_found = 0
          do j = 1, ndihe_p
            if (no_wild(j)) then
              if (((ci1 == par%dihe_atom_cls(1,j)) .and. &
                   (ci2 == par%dihe_atom_cls(2,j)) .and. &
                   (ci3 == par%dihe_atom_cls(3,j)) .and. &
                   (ci4 == par%dihe_atom_cls(4,j))) .or. &
                  ((ci1 == par%dihe_atom_cls(4,j)) .and. &
                   (ci2 == par%dihe_atom_cls(3,j)) .and. &
                   (ci3 == par%dihe_atom_cls(2,j)) .and. &
                   (ci4 == par%dihe_atom_cls(1,j)))) then

                nw_found = nw_found + 1
                dihedral(icel_local) = dihedral(icel_local) + 1
                dlist(1:4,dihedral(icel_local),icel_local) = list(1:4)

                force (dihedral(icel_local),icel_local) = &
                     par%dihe_force_const(j)
                period(dihedral(icel_local),icel_local) = &
                     par%dihe_periodicity(j)
                phase (dihedral(icel_local),icel_local) = &
                     par%dihe_phase(j) * RAD

                if (period(dihedral(icel_local),icel_local) >  &
                enefunc%notation_14types) &
                call error_msg('Setup_Enefunc_Dihe> Too many periodicity.')

              end if
            end if
          end do

          if (nw_found == 0) then
            do j = 1, ndihe_p
              if (.not.no_wild(j)) then
                if (((ci2 == par%dihe_atom_cls(2,j)) .and. &
                     (ci3 == par%dihe_atom_cls(3,j))) .or. &
                    ((ci2 == par%dihe_atom_cls(3,j)) .and. &
                     (ci3 == par%dihe_atom_cls(2,j)))) then

                  dihedral(icel_local) = dihedral(icel_local) + 1
                  dlist(1:4,dihedral(icel_local),icel_local) = list(1:4)

                  force (dihedral(icel_local),icel_local) = &
                       par%dihe_force_const(j)
                  period(dihedral(icel_local),icel_local) = &
                       par%dihe_periodicity(j)
                  phase (dihedral(icel_local),icel_local) = &
                       par%dihe_phase(j) * RAD
                  if (period(dihedral(icel_local),icel_local) >  &
                  enefunc%notation_14types) &
                  call error_msg('Setup_Enefunc_Dihe> Too many periodicity.')

                end if
              end if
            end do
          end if

        end if

      end if

    end do

    deallocate(no_wild)

    found = 0
    do i = 1, ncel
      found = found + dihedral(i)
      if (dihedral(i) > MaxDihe) &
        call error_msg('Setup_Enefunc_Dihe> Too many dihedral angles.')
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_dihe_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_dihe_all = found
#endif

    if (domain%fep_use) then
      ndihe_fep = 0
      do i = 1, 5
        ndihe_fep = ndihe_fep + molecule%num_dihedrals_fep(i)
      end do
      if (enefunc%num_dihe_all < ndihe_fep) &
        call error_msg( &
           'Setup_Enefunc_Dihe> Some dihedral paremeters are missing.')
    else
      if (enefunc%num_dihe_all < ndihe) &
        call error_msg( &
           'Setup_Enefunc_Dihe> Some dihedral paremeters are missing.')
    end if

    return

  end subroutine setup_enefunc_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_impr
  !> @brief        define IMPROPER term in potential energy function
  !! @authors      YS,JJ
  !! @param[in]    par      : CHARMM PAR information
  !! @param[in]    molecule : molecule information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_impr(par, molecule, domain, enefunc)

    ! formal variables
    type(s_par),      target, intent(in)    :: par
    type(s_molecule), target, intent(in)    :: molecule
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(inout) :: enefunc

    ! local variables
    integer                   :: nimpr, nimpr_p
    integer                   :: i, j, icel_local
    integer                   :: icel1, icel2
    integer                   :: found
    integer                   :: list(4)
    character(6)              :: ci1, ci2, ci3, ci4

    real(wp),         pointer :: force(:,:), phase(:,:)
    integer,          pointer :: improper(:), ilist(:,:,:)
    integer,          pointer :: ncel, cell_pair(:,:), id_g2l(:,:)
    integer,      allocatable :: wc_type(:)
    logical,      allocatable :: no_wild(:)
    integer                   :: nimpr_fep, i1, i2, i3, i4

    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l

    improper  => enefunc%num_improper
    ilist     => enefunc%impr_list
    force     => enefunc%impr_force_const
    phase     => enefunc%impr_phase

    nimpr     = molecule%num_impropers
    nimpr_p   = par%num_impropers

    ! check usage of wild card
    !
    allocate(wc_type(nimpr_p), no_wild(nimpr))

    do i = 1, nimpr_p

      if ((par%impr_atom_cls(1,i) /= WildCard) .and. &
          (par%impr_atom_cls(2,i) /= WildCard) .and. &
          (par%impr_atom_cls(3,i) /= WildCard) .and. &
          (par%impr_atom_cls(4,i) /= WildCard)) then

        ! A-B-C-D type
        wc_type(i) = 0

      else if ((par%impr_atom_cls(1,i) == WildCard) .and. &
               (par%impr_atom_cls(2,i) /= WildCard) .and. &
               (par%impr_atom_cls(3,i) /= WildCard) .and. &
               (par%impr_atom_cls(4,i) /= WildCard)) then

        ! X-B-C-D type
        wc_type(i) = 1

      else if ((par%impr_atom_cls(1,i) == WildCard) .and. &
               (par%impr_atom_cls(2,i) == WildCard) .and. &
               (par%impr_atom_cls(3,i) /= WildCard) .and. &
               (par%impr_atom_cls(4,i) /= WildCard)) then

        ! X-X-C-D type
        wc_type(i) = 2

      else if ((par%impr_atom_cls(1,i) /= WildCard) .and. &
               (par%impr_atom_cls(2,i) == WildCard) .and. &
               (par%impr_atom_cls(3,i) == WildCard) .and. &
               (par%impr_atom_cls(4,i) /= WildCard)) then

        ! A-X-X-D type
        wc_type(i) = 3

      else
        call error_msg('Setup_Enefunc_Impr> Undefined Wild Card')

      end if

    end do

    ! setup parameters
    !
    do i = 1, nimpr

      no_wild(i) = .false.

      ci1 = molecule%atom_cls_name(molecule%impr_list(1,i))
      ci2 = molecule%atom_cls_name(molecule%impr_list(2,i))
      ci3 = molecule%atom_cls_name(molecule%impr_list(3,i))
      ci4 = molecule%atom_cls_name(molecule%impr_list(4,i))
      list(1:4) = molecule%impr_list(1:4,i)

      icel1 = id_g2l(1,list(1))
      icel2 = id_g2l(1,list(4))

      if (domain%fep_use) then
        ! FEP: If the improper are not set to any group of FEP, exclude this improper.
        i1 = molecule%fepgrp(list(1))
        i2 = molecule%fepgrp(list(2))
        i3 = molecule%fepgrp(list(3))
        i4 = molecule%fepgrp(list(4))
        if (molecule%fepgrp_dihe(i1,i2,i3,i4) == 0) cycle
      end if

      if (icel1 /= 0 .and. icel2 /= 0) then

        icel_local = cell_pair(icel1,icel2)

        if (icel_local >= 1 .and. icel_local <= ncel) then

          ! A-B-C-D type
          !
          do j = 1, nimpr_p
            if (wc_type(j) == 0) then
              if (((ci1 == par%impr_atom_cls(1,j)) .and. &
                   (ci2 == par%impr_atom_cls(2,j)) .and. &
                   (ci3 == par%impr_atom_cls(3,j)) .and. &
                   (ci4 == par%impr_atom_cls(4,j))) .or. &
                  ((ci1 == par%impr_atom_cls(4,j)) .and. &
                   (ci2 == par%impr_atom_cls(3,j)) .and. &
                   (ci3 == par%impr_atom_cls(2,j)) .and. &
                   (ci4 == par%impr_atom_cls(1,j)))) then

                improper(icel_local) = improper(icel_local) + 1
                ilist(1:4,improper(icel_local),icel_local) = list(1:4)

                force(improper(icel_local),icel_local) = &
                     par%impr_force_const(j)
                phase(improper(icel_local),icel_local) = &
                     par%impr_phase(j) * RAD
                no_wild(i) = .true.
                exit

              end if
            end if
          end do

          ! X-B-C-D type
          !
          if (.not.no_wild(i)) then
            do j = 1, nimpr_p
              if (wc_type(j) == 1) then
                if (((ci2 == par%impr_atom_cls(2,j)) .and. &
                     (ci3 == par%impr_atom_cls(3,j)) .and. &
                     (ci4 == par%impr_atom_cls(4,j))) .or. &
                    ((ci2 == par%impr_atom_cls(4,j)) .and. &
                     (ci3 == par%impr_atom_cls(3,j)) .and. &
                     (ci4 == par%impr_atom_cls(2,j)))) then

                  improper(icel_local) = improper(icel_local) + 1
                  ilist(1:4,improper(icel_local),icel_local) = list(1:4)

                  force(improper(icel_local),icel_local) = &
                       par%impr_force_const(j)
                  phase(improper(icel_local),icel_local) = &
                       par%impr_phase(j) * RAD
                  no_wild(i) = .true.
                  exit

                end if
              end if
            end do
          end if

          ! X-X-C-D type
          !
          if (.not.no_wild(i)) then
            do j = 1, nimpr_p
              if (wc_type(j) == 2) then
                if (((ci3 == par%impr_atom_cls(3,j)) .and. &
                     (ci4 == par%impr_atom_cls(4,j))) .or. &
                    ((ci3 == par%impr_atom_cls(4,j)) .and. &
                     (ci4 == par%impr_atom_cls(3,j)))) then

                  improper(icel_local) = improper(icel_local) + 1
                  ilist(1:4,improper(icel_local),icel_local) = list(1:4)

                  force(improper(icel_local),icel_local) = &
                       par%impr_force_const(j)
                  phase(improper(icel_local),icel_local) = &
                       par%impr_phase(j) * RAD
                  no_wild(i) = .true.
                  exit

                end if
              end if
            end do
          end if

          ! A-X-X-D type
          !
          if (.not.no_wild(i)) then
            do j = 1, nimpr_p
              if (wc_type(j) == 3) then
                if (((ci1 == par%impr_atom_cls(1,j)) .and. &
                     (ci4 == par%impr_atom_cls(4,j))) .or. &
                    ((ci1 == par%impr_atom_cls(4,j)) .and. &
                     (ci4 == par%impr_atom_cls(1,j)))) then

                  improper(icel_local) = improper(icel_local) + 1
                  ilist(1:4,improper(icel_local),icel_local) = list(1:4)

                  force(improper(icel_local),icel_local) = &
                       par%impr_force_const(j)
                  phase(improper(icel_local),icel_local) = &
                       par%impr_phase(j) * RAD
                  no_wild(i) = .true.
                  exit

                end if
              end if
            end do
          end if

          if (.not.no_wild(i)) &
            write(MsgOut,*) &
              'Setup_Enefunc_Impr> Unknown IMPR type. [', &
              ci1, ']-[', ci2, ']-[', ci3, ']-[', ci4, '] (ERROR)'

        end if
      end if

    end do

    deallocate(wc_type, no_wild)

    found = 0
    do i = 1, ncel
      found = found + improper(i)
      if (improper(i) > MaxImpr) &
        call error_msg('Setup_Enefunc_Impr> Too many improper dihedral angles')
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_impr_all, 1, mpi_integer, mpi_sum, &
                       mpi_comm_country, ierror)
#else
    enefunc%num_impr_all = found
#endif

    if (domain%fep_use) then
      nimpr_fep = 0
      do i = 1, 5
        nimpr_fep = nimpr_fep + molecule%num_impropers_fep(i)
      end do
      if (enefunc%num_impr_all < nimpr_fep) &
        call error_msg( &
          'Setup_Enefunc_Impr> Some improper paremeters are missing.')
    else
      if (enefunc%num_impr_all < nimpr) &
        call error_msg( &
          'Setup_Enefunc_Impr> Some improper paremeters are missing.')
    end if

    return

  end subroutine setup_enefunc_impr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_cmap
  !> @brief        define cmap term in potential energy function with DD
  !! @authors      TY, TM
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[in]    par         : CHARMM PAR information
  !! @param[in]    molecule    : molecule information
  !! @param[in]    domain      : domain information
  !! @param[inout] enefunc     : energy potential functions informationn
  !!
  !! @note       In str "par", following variables have been defined:
  !!   cmap_atom_cls(8,imap) (char) : 8 atom classes (4 for phi and 4 for psi)
  !!   cmap_resolution(imap) (int ) : = 24 (for all imap) = 360degree/15degree
  !!   cmap_data(i,j,imap)   (real) : Ecmap for grid points. i and j are the
  !!                                  1st (psi?) and the 2nd (phi?) grid IDs.
  !!                                  1<=i,j<=24.
  !!   Where imap is ID of cmap type (1 <= imap <= 6).
  !!
  !! @note       TY determined to use natural spline (periodic = .false.)
  !!             because it is similar to the CHARMM way.
  !!             However, I notice that force will be more accurately continuous
  !!             at phi (or psi) = +-180 degree if periodic spline was used.
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_cmap(ene_info, par, molecule, domain, enefunc)

    ! formal variables
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_par),             intent(in)    :: par
    type(s_molecule),        intent(in)    :: molecule
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k, l, ityp
    integer                  :: ncmap_p, found, ngrid0
    integer                  :: list(2), icel1, icel2, icel_local
    integer                  :: flag_cmap_type, alloc_stat, dealloc_stat
    character(6)             :: ci1, ci2, ci3, ci4, ci5, ci6, ci7, ci8
    logical                  :: periodic

    integer,         pointer :: ncel, cell_pair(:,:), id_g2l(:,:)

    real(wp),    allocatable :: c_ij(:,:,:,:) ! cmap coeffs

    integer                  :: ncmap_fep, idx
    integer                  :: i1, i2, i3, i4, i5, i6, i7, i8

    ncel            => domain%num_cell_local
    cell_pair       => domain%cell_pair
    id_g2l          => domain%id_g2l


    ! If 'periodic' is .true.,
    !   then cubic spline with periodic (in dihedral-angle space) boundary
    !   will be applied.
    ! If 'periodic' is .false.,
    !   then natural cubic spline without periodic boudnary
    !   will be applied to expanded cross-term maps.
    !   This is similar with CHARMM's source code.
    !
    periodic = ene_info%cmap_pspline
    ncmap_p  = par%num_cmaps

    ! begin
    !
    ngrid0 = 0
    do i = 1, ncmap_p
      ngrid0 = max(ngrid0, par%cmap_resolution(i))
    end do

    call alloc_enefunc(enefunc, EneFuncCmap, ncel, ngrid0, ncmap_p)

    alloc_stat = 0
    allocate(c_ij(4,4,ngrid0,ngrid0), stat = alloc_stat)
    if (alloc_stat /= 0) &
      call error_msg_alloc

    do i = 1, ncmap_p
      enefunc%cmap_resolution(i) = par%cmap_resolution(i)
    end do

    ! derive cmap coefficients by bicubic interpolation
    !
    do ityp = 1, ncmap_p

      if (periodic) then
        call derive_cmap_coefficients_p(ityp, par, c_ij)
      else
        call derive_cmap_coefficients_np(ityp, par, c_ij)
      end if

      do l = 1, ngrid0
        do k = 1, ngrid0
          do j = 1, 4
            do i = 1, 4
              enefunc%cmap_coef(i,j,k,l,ityp) = c_ij(i,j,k,l)
            end do
          end do
        end do
      end do
    end do

    enefunc%num_cmap(1:ncel) = 0
    do i = 1, molecule%num_cmaps

      list(1) = molecule%cmap_list(1,i)
      list(2) = molecule%cmap_list(8,i)

      icel1 = id_g2l(1,list(1))
      icel2 = id_g2l(1,list(2))

      if (domain%fep_use) then
        ! FEP: If the cmap are not set to any group of FEP, exclude this cmap.
        i1 = molecule%fepgrp(molecule%cmap_list(1,i))
        i2 = molecule%fepgrp(molecule%cmap_list(2,i))
        i3 = molecule%fepgrp(molecule%cmap_list(3,i))
        i4 = molecule%fepgrp(molecule%cmap_list(4,i))
        i5 = molecule%fepgrp(molecule%cmap_list(5,i))
        i6 = molecule%fepgrp(molecule%cmap_list(6,i))
        i7 = molecule%fepgrp(molecule%cmap_list(7,i))
        i8 = molecule%fepgrp(molecule%cmap_list(8,i))
        idx = i1 + 5*(i2-1 + 5*(i3-1 + 5*(i4-1 + 5*(i5-1 + 5*(i6-1 + 5*(i7-1 + 5*(i8-1)))))))
        if (molecule%fepgrp_cmap(idx) == 0) cycle
      end if

      if (icel1 /= 0 .and. icel2 /= 0) then

        icel_local = cell_pair(icel1,icel2)

        if (icel_local >= 1 .and. icel_local <= ncel) then

          ! ci* will be atom-type strings
          !
          ci1 = molecule%atom_cls_name(molecule%cmap_list(1,i))
          ci2 = molecule%atom_cls_name(molecule%cmap_list(2,i))
          ci3 = molecule%atom_cls_name(molecule%cmap_list(3,i))
          ci4 = molecule%atom_cls_name(molecule%cmap_list(4,i))
          ci5 = molecule%atom_cls_name(molecule%cmap_list(5,i))
          ci6 = molecule%atom_cls_name(molecule%cmap_list(6,i))
          ci7 = molecule%atom_cls_name(molecule%cmap_list(7,i))
          ci8 = molecule%atom_cls_name(molecule%cmap_list(8,i))
          flag_cmap_type = -1

          ! assign cmap type ID to each (psi,phi) pair
          !
          do j = 1, ncmap_p
            if (ci1 == par%cmap_atom_cls(1, j) .and. &
                ci2 == par%cmap_atom_cls(2, j) .and.   &
                ci3 == par%cmap_atom_cls(3, j) .and.   &
                ci4 == par%cmap_atom_cls(4, j) .and.   &
                ci5 == par%cmap_atom_cls(5, j) .and.   &
                ci6 == par%cmap_atom_cls(6, j) .and.   &
                ci7 == par%cmap_atom_cls(7, j) .and.   &
                ci8 == par%cmap_atom_cls(8, j)   ) then

              enefunc%num_cmap(icel_local) = enefunc%num_cmap(icel_local) + 1
              enefunc%cmap_list(1:8,enefunc%num_cmap(icel_local),icel_local) &
                   = molecule%cmap_list(1:8,i)
              enefunc%cmap_type(enefunc%num_cmap(icel_local),icel_local) = j
              flag_cmap_type = j
              exit

            end if
          end do

          ! if not found, print detail about the error.
          !

          if (flag_cmap_type <= 0) then
            write(MsgOut,*) 'Setup_Enefunc_Cmap> not found CMAP: '
            write(MsgOut,*) ' [',ci1,']-[',ci2,']-[',ci3,']-[',ci4,']-'
            write(MsgOut,*) ' [',ci5,']-[',ci6,']-[',ci7,']-[',ci8,'] '
            write(MsgOut,*) '  in parameter file. (ERROR)'
          end if

        end if
      end if

    end do

    deallocate(c_ij, stat=dealloc_stat)
    if (dealloc_stat /= 0) &
      call error_msg_dealloc

    ! write summary
    !
    if (main_rank) then
      if (periodic) then
        write(MsgOut,'(A)') &
    'Setup_Enefunc_Cmap> Periodic-boundary spline is used to derive cmap coefs.'
        write(MsgOut,'(A)') ''
      else
        write(MsgOut,'(A)') &
            'Setup_Enefunc_Cmap> Natural spline is used to derive cmap coefs.'
        write(MsgOut,'(A)') ''
      end if
    end if

    ! stop if parameter is missing
    !
    found = 0
    do i = 1, ncel
      found = found + enefunc%num_cmap(i)
      if (enefunc%num_cmap(i) > MaxCmap) &
        call error_msg('Setup_Enefunc_Cmap> Too many cmaps.')
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_cmap_all, 1, mpi_integer, mpi_sum, &
                       mpi_comm_country, ierror)
#else
    enefunc%num_cmap_all = found
#endif

    if (domain%fep_use) then
      ncmap_fep = 0
      do i = 1, 5
        ncmap_fep = ncmap_fep + molecule%num_cmaps_fep(i)
      end do
      if (enefunc%num_cmap_all /= ncmap_fep) &
        call error_msg('Setup_Enefunc_Cmap> Some cmap parameters are missing.')
    else
      if (enefunc%num_cmap_all /=  molecule%num_cmaps) &
        call error_msg('Setup_Enefunc_Cmap> Some cmap parameters are missing.')
    end if

    return

  end subroutine setup_enefunc_cmap

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_nonb
  !> @brief        define NON-BOND term in potential energy function
  !! @authors      YS, JJ, TI
  !! @param[in]    par         : CHARMM PAR information
  !! @param[in]    molecule    : molecule information
  !! @param[in]    constraints : constraints information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_nonb(par, molecule, constraints, domain, enefunc)

    ! formal arguments
    type(s_par),             intent(in)    :: par
    type(s_molecule),        intent(in)    :: molecule
    type(s_constraints),     intent(in)    :: constraints
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    real(wp)                 :: eps14, rmin14, eps, rmin
    integer                  :: i, j, k, ix, jx, kx
    integer                  :: identical_sign
    integer                  :: nonb_p, nbfx_p, cls_local, ncel
    character(6)             :: ci1, ci2

    integer,  allocatable    :: nonb_atom_cls(:), check_cls(:)
    integer,  allocatable    :: atmcls_map_g2l(:), atmcls_map_l2g(:)
    real(wp), allocatable    :: nb14_lj6(:,:), nb14_lj12(:,:)
    real(wp), allocatable    :: nonb_lj6(:,:), nonb_lj12(:,:)

    enefunc%num_atom_cls = par%num_atom_cls

    ELECOEF          = ELECOEF_CHARMM

    ! set lennard-jones parameters
    !
    nonb_p = enefunc%num_atom_cls

    allocate(nonb_atom_cls(nonb_p),     &
             check_cls(nonb_p),         &
             atmcls_map_g2l(nonb_p),    &
             atmcls_map_l2g(nonb_p),    &
             nb14_lj6 (nonb_p, nonb_p), &
             nb14_lj12(nonb_p, nonb_p), &
             nonb_lj6 (nonb_p, nonb_p), &
             nonb_lj12(nonb_p, nonb_p))

    nonb_atom_cls(1:nonb_p)           = 0
    check_cls(1:nonb_p)               = 0
    nb14_lj6     (1:nonb_p, 1:nonb_p) = 0.0_wp
    nb14_lj12    (1:nonb_p, 1:nonb_p) = 0.0_wp
    nonb_lj6     (1:nonb_p, 1:nonb_p) = 0.0_wp
    nonb_lj12    (1:nonb_p, 1:nonb_p) = 0.0_wp

    do i = 1, nonb_p

      lj_coef(1,i) = par%nonb_eps(i)
      lj_coef(2,i) = par%nonb_rmin(i)

      do j = 1, nonb_p

        ! combination rule
        eps14  = sqrt(par%nonb_eps_14(i) * par%nonb_eps_14(j))
        rmin14 = par%nonb_rmin_14(i) + par%nonb_rmin_14(j)
        eps    = sqrt(par%nonb_eps(i) * par%nonb_eps(j))
        rmin   = par%nonb_rmin(i) + par%nonb_rmin(j)

        ! set parameters
        nb14_lj12(i,j) = eps14 * (rmin14 ** 12)
        nb14_lj6(i,j)  = 2.0_wp * eps14 * (rmin14 ** 6)
        nonb_lj12(i,j) = eps * (rmin ** 12)
        nonb_lj6(i,j)  = 2.0_wp * eps * (rmin ** 6)

      end do
    end do

    ! overwrite lennard-jones parameters by NBFIX parameters
    !
    nbfx_p = par%num_nbfix

    do k = 1, nbfx_p
      ci1 = par%nbfi_atom_cls(1,k)
      ci2 = par%nbfi_atom_cls(2,k)
      do i = 1, nonb_p
        do j = 1, nonb_p
          if ((ci1 == par%nonb_atom_cls(i)  .and. &
               ci2 == par%nonb_atom_cls(j)) .or.  &
              (ci2 == par%nonb_atom_cls(i) .and. &
               ci1 == par%nonb_atom_cls(j))) then

            ! combination rule
            !
            eps14  = abs(par%nbfi_eps_14 (k)) !TODO CHECK
            rmin14 = par%nbfi_rmin_14(k)
            eps    = abs(par%nbfi_eps    (k))
            rmin   = par%nbfi_rmin   (k)

            ! set parameters
            !
            nb14_lj12(i,j) = eps14 * (rmin14 ** 12)
            nb14_lj6 (i,j) = 2.0_wp * eps14 * (rmin14 ** 6)
            nonb_lj12(i,j) = eps * (rmin ** 12)
            nonb_lj6 (i,j) = 2.0_wp * eps * (rmin ** 6)
          end if
        end do
      end do
    end do

    ! check the usage of atom class
    !
    do i = 1, molecule%num_atoms
      k = molecule%atom_cls_no(i)
      if (k < 1) then
        call error_msg( &
        'Setup_Enefunc_Nonb> atom class is not defined: "'&
        //trim(molecule%atom_cls_name(i))//'"')
      endif

      check_cls(k) = 1
    end do

    k = 0
    do i = 1, nonb_p
      if (check_cls(i) == 1) then
        k = k + 1
        atmcls_map_g2l(i) = k
        atmcls_map_l2g(k) = i
      end if
    end do
    cls_local = k
    max_class = cls_local

    call alloc_enefunc(enefunc, EneFuncNbon, cls_local)

    do i = 1, cls_local
      ix = atmcls_map_l2g(i)
      do j = 1, cls_local
        jx = atmcls_map_l2g(j)
        enefunc%nb14_lj12(i,j) = nb14_lj12(ix,jx)
        enefunc%nb14_lj6 (i,j) = nb14_lj6 (ix,jx)
        enefunc%nonb_lj12(i,j) = nonb_lj12(ix,jx)
        enefunc%nonb_lj6 (i,j) = nonb_lj6 (ix,jx)
      end do
    end do

    ! update domain information
    !
    do i = 1, domain%num_cell_local+domain%num_cell_boundary
      do ix = 1, domain%num_atom(i)
        domain%atom_cls_no(ix,i) = atmcls_map_g2l(domain%atom_cls_no(ix,i))
      end do
    end do

    if (enefunc%table%num_water > 0) then
      domain%water%atom_cls_no(1:3) =  &
          atmcls_map_g2l(domain%water%atom_cls_no(1:3))
      enefunc%table%atom_cls_no_O   =  &
         atmcls_map_g2l(enefunc%table%atom_cls_no_O)
      enefunc%table%atom_cls_no_H   =  &
         atmcls_map_g2l(enefunc%table%atom_cls_no_H)

      if (constraints%tip4 .or. enefunc%table%tip4) then
        domain%water%atom_cls_no(4) =  &
           atmcls_map_g2l(domain%water%atom_cls_no(4))
        enefunc%table%atom_cls_no_D =  &
           atmcls_map_g2l(enefunc%table%atom_cls_no_D)
      end if
    end if

    deallocate(nonb_atom_cls,  &
               check_cls,      &
               atmcls_map_g2l, &
               atmcls_map_l2g, &
               nb14_lj6,       &
               nb14_lj12,      &
               nonb_lj6,       &
               nonb_lj12)

    enefunc%num_atom_cls = cls_local

    ! treatment for 1-2, 1-3, 1-4 interactions
    !
    ncel   = domain%num_cell_local

    call alloc_enefunc(enefunc, EneFuncNonb,     ncel, maxcell_near)
    call alloc_enefunc(enefunc, EneFuncNonbList, ncel, maxcell_near)

    if (constraints%rigid_bond) then
      call count_nonb_excl_constraint(.true., .false., constraints, &
                                      domain, enefunc)

    else
      call count_nonb_excl(.true., domain, enefunc)

    end if

    return

  end subroutine setup_enefunc_nonb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    count_nonb_excl
  !> @brief        exclude 1-2, 1-3 interactions
  !! @authors      JJ
  !! @param[in]    first   : flag for first call or not
  !! @param[inout] domain  : structure of domain
  !! @param[inout] enefunc : structure of enefunc
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine count_nonb_excl(first, domain, enefunc)

    ! formal arguments
    logical,                 intent(in)    :: first
    type(s_domain),  target, intent(inout) :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: ncell, i, ii, ix, k, i1, i2
    integer                  :: index(4)
    integer                  :: ij ,j
    integer                  :: icel, icel1, icel2, jcel
    integer                  :: ini_nb14, fin_nb14
    integer                  :: num_excl, num_nb14
    integer                  :: found1, found2
    logical                  :: duplicate

    integer,         pointer :: ncell_local, max_atom, natom(:), id_g2l(:,:)
    integer,         pointer :: nwater(:), water_list(:,:,:)
    integer,         pointer :: cell_pair(:,:)
    integer,         pointer :: num_nonb_excl(:,:), num_nonb_excl1(:,:)
    integer,         pointer :: num_nb14_calc(:,:), num_nb14_calc1(:,:)
    integer,         pointer :: nbond(:), bondlist(:,:,:)
    integer,         pointer :: nangle(:), anglelist(:,:,:)
    integer,         pointer :: ndihedral(:), dihelist(:,:,:), pairlist(:,:)
    integer,         pointer :: nonb_excl_list(:,:,:), nonb_excl_list1(:,:,:)
    integer,         pointer :: nb14_calc_list(:,:,:), nb14_calc_list1(:,:,:)
    integer,         pointer :: sc_calc_list(:,:,:), sc_calc_list1(:,:,:)
    integer,         pointer :: num_excl_total(:), num_nb14_total(:)
    integer,         pointer :: num_excl_total1(:), num_nb14_total1(:)


    ncell_local     => domain%num_cell_local
    max_atom        => domain%max_num_atom
    natom           => domain%num_atom
    nwater          => domain%num_water
    water_list      => domain%water_list
    id_g2l          => domain%id_g2l
    cell_pair       => domain%cell_pairlist1
    pairlist        => domain%cell_pairlist2

    nbond           => enefunc%num_bond
    bondlist        => enefunc%bond_list
    nangle          => enefunc%num_angle
    anglelist       => enefunc%angle_list
    ndihedral       => enefunc%num_dihedral
    dihelist        => enefunc%dihe_list
    nonb_excl_list  => enefunc%nonb_list
    nonb_excl_list1 => enefunc%nonb_list1
    nb14_calc_list  => enefunc%nb14_list
    nb14_calc_list1 => enefunc%nb14_list1
    sc_calc_list    => enefunc%sc_list
    sc_calc_list1   => enefunc%sc_list1
    num_excl_total  => enefunc%num_excl_total
    num_excl_total1 => enefunc%num_excl_total1
    num_nb14_total  => enefunc%num_nb14_total
    num_nb14_total1 => enefunc%num_nb14_total1
    num_nonb_excl   => enefunc%num_nonb_excl
    num_nonb_excl1  => enefunc%num_nonb_excl1
    num_nb14_calc   => enefunc%num_nb14_calc
    num_nb14_calc1  => enefunc%num_nb14_calc1

    ncell = domain%num_cell_local + domain%num_cell_boundary

    max_atom = 0
    do i = 1, ncell
      max_atom = max(max_atom,natom(i))
    end do

    ! initialization
    !
    num_nonb_excl  (1:max_atom,1:maxcell_near)= 0
    num_nonb_excl1 (1:max_atom,1:ncell_local) = 0
    num_nb14_calc  (1:max_atom,1:maxcell_near)= 0
    num_nb14_calc1 (1:max_atom,1:ncell_local) = 0
    num_excl_total (1:maxcell_near)           = 0
    num_excl_total1(1:ncell_local)            = 0
    num_nb14_total (1:maxcell_near)           = 0
    num_nb14_total1(1:ncell_local)            = 0

    ! exclude 1-2 interaction
    !
    if (enefunc%excl_level > 0) then

      do i = 1, ncell_local
        do ix = 1, nbond(i)

          icel1 = id_g2l(1,bondlist(1,ix,i))
          i1    = id_g2l(2,bondlist(1,ix,i))
          icel2 = id_g2l(1,bondlist(2,ix,i))
          i2    = id_g2l(2,bondlist(2,ix,i))

          if (icel1 < icel2) then
            num_excl = num_nonb_excl(i1,pairlist(icel2,icel1)) + 1
            num_nonb_excl(i1,pairlist(icel2,icel1)) = num_excl
            nonb_excl_list(num_excl,i1,pairlist(icel2,icel1)) = i2
            num_excl_total(pairlist(icel2,icel1)) = &
                 num_excl_total(pairlist(icel2,icel1)) + 1

          else if (icel1 > icel2) then
            num_excl = num_nonb_excl(i2,pairlist(icel1,icel2)) + 1
            num_nonb_excl(i2,pairlist(icel1,icel2)) = num_excl
            nonb_excl_list(num_excl,i2,pairlist(icel1,icel2)) = i1
            num_excl_total(pairlist(icel1,icel2)) = &
                 num_excl_total(pairlist(icel1,icel2)) + 1

          else if (i1 < i2) then
            num_excl = num_nonb_excl1(i1,i) + 1
            num_nonb_excl1(i1,i) = num_excl
            nonb_excl_list1(num_excl,i1,i) = i2
            num_excl_total1(i) = num_excl_total1(i) + 1

          else if (i1 > i2) then
            num_excl = num_nonb_excl1(i2,i) + 1
            num_nonb_excl1(i2,i) = num_excl
            nonb_excl_list1(num_excl,i2,i) = i1
            num_excl_total1(i) = num_excl_total1(i) + 1

          end if
        end do
      end do

    end if

    ! exclude 1-3 interaction
    !
    if (enefunc%excl_level > 1) then

      do i = 1, ncell_local
        do ix = 1, nangle(i)

          icel1 = id_g2l(1,anglelist(1,ix,i))
          i1    = id_g2l(2,anglelist(1,ix,i))
          icel2 = id_g2l(1,anglelist(3,ix,i))
          i2    = id_g2l(2,anglelist(3,ix,i))

          if (icel1 < icel2) then

            num_excl = num_nonb_excl(i1,pairlist(icel2,icel1))
            duplicate = .false.

            do k = 1, num_excl
              if (i2 == nonb_excl_list(k,i1,pairlist(icel2,icel1))) &
                duplicate = .true.
            end do

            if (.not. duplicate) then
              num_excl = num_excl + 1
              num_nonb_excl(i1,pairlist(icel2,icel1)) = num_excl
              nonb_excl_list(num_excl,i1,pairlist(icel2,icel1)) = i2
              num_excl_total(pairlist(icel2,icel1)) = &
                   num_excl_total(pairlist(icel2,icel1)) + 1
            end if

          else if (icel1 > icel2) then

            num_excl = num_nonb_excl(i2,pairlist(icel1,icel2))
            duplicate = .false.

            do k = 1, num_excl
              if (i1 == nonb_excl_list(k,i2,pairlist(icel1,icel2))) &
                duplicate = .true.
            end do

            if (.not.duplicate) then
              num_excl = num_excl + 1
              num_nonb_excl(i2,pairlist(icel1,icel2)) = num_excl
              nonb_excl_list(num_excl,i2,pairlist(icel1,icel2)) = i1
              num_excl_total(pairlist(icel1,icel2)) = &
                   num_excl_total(pairlist(icel1,icel2)) + 1
            end if

          else if (i1 < i2) then

            num_excl = num_nonb_excl1(i1,i)
            duplicate = .false.

            do k = 1, num_excl
              if (i2 == nonb_excl_list1(k,i1,i)) &
                duplicate = .true.
            end do

            if (.not.duplicate) then
              num_excl = num_excl + 1
              num_nonb_excl1(i1,i) = num_excl
              nonb_excl_list1(num_excl,i1,i) = i2
              num_excl_total1(i) = num_excl_total1(i) + 1
            end if

          else if (i1 > i2) then

            num_excl = num_nonb_excl1(i2,i)
            duplicate = .false.

            do k = 1, num_excl
              if (i1 == nonb_excl_list1(k,i2,i)) &
                duplicate = .true.
            end do

            if (.not.duplicate) then
              num_excl = num_excl + 1
              num_nonb_excl1(i2,i) = num_excl
              nonb_excl_list1(num_excl,i2,i) = i1
              num_excl_total1(i) = num_excl_total1(i) + 1
            end if

          end if

        end do
      end do

    end if

    ! exclude water
    !
    if (enefunc%excl_level > 1) then

      if (enefunc%table%tip4) then

        do i = 1, ncell_local
          do ix = 1, nwater(i)

            index(1:4) = water_list(1:4,ix,i)

            do i1 = 1, 3
              do i2 = i1+1, 4

                num_excl = num_nonb_excl1(index(i1),i)
                duplicate = .false.

                do k = 1, num_excl
                  if (index(i2) == nonb_excl_list1(k,index(i1),i)) &
                    duplicate = .true.
                end do

                if (.not. duplicate) then
                  num_excl = num_excl + 1
                  num_nonb_excl1(index(i1),i) = num_excl
                  nonb_excl_list1(num_excl,index(i1),i) = index(i2)
                  num_excl_total1(i) = num_excl_total1(i) + 1
                end if

              end do
            end do

          end do
        end do

      end if
    end if

    ! count 1-4 interaction
    !
    if (enefunc%excl_level > 2) then

      do ii = 1, 2

        if (ii == 2) then
          ndihedral => enefunc%num_rb_dihedral
          dihelist  => enefunc%rb_dihe_list
        end if

        do i = 1, ncell_local
          do ix = 1, ndihedral(i)

            icel1 = id_g2l(1,dihelist(1,ix,i))
            i1    = id_g2l(2,dihelist(1,ix,i))
            icel2 = id_g2l(1,dihelist(4,ix,i))
            i2    = id_g2l(2,dihelist(4,ix,i))

            if (icel1 < icel2) then

              num_nb14 = num_nb14_calc(i1,pairlist(icel2,icel1))
              duplicate = .false.

              do k = 1, num_nonb_excl(i1,pairlist(icel2,icel1))
                if (i2 == nonb_excl_list(k,i1,pairlist(icel2,icel1))) &
                  duplicate = .true.
              end do

              do k = 1, num_nb14
                if (i2 == nb14_calc_list(k,i1,pairlist(icel2,icel1))) &
                  duplicate = .true.
              end do

              if (.not.duplicate) then
                num_nb14 = num_nb14 + 1
                num_nb14_calc(i1,pairlist(icel2,icel1)) = num_nb14
                nb14_calc_list(num_nb14,i1,pairlist(icel2,icel1)) = i2
                num_nb14_total(pairlist(icel2,icel1)) = &
                     num_nb14_total(pairlist(icel2,icel1)) + 1
                sc_calc_list(num_nb14,i1,pairlist(icel2,icel1)) = &
                   int(enefunc%dihe_periodicity(ix,i)/enefunc%notation_14types)
              end if

            else if (icel1 > icel2) then

              num_nb14 = num_nb14_calc(i2,pairlist(icel1,icel2))
              duplicate = .false.

              do k = 1, num_nonb_excl(i2,pairlist(icel1,icel2))
                if (i1 == nonb_excl_list(k,i2,pairlist(icel1,icel2))) &
                  duplicate = .true.
              end do

              do k = 1, num_nb14
                if (i1 == nb14_calc_list(k,i2,pairlist(icel1,icel2))) &
                  duplicate = .true.
              end do

              if (.not.duplicate) then
                num_nb14 = num_nb14 + 1
                num_nb14_calc(i2,pairlist(icel1,icel2)) = num_nb14
                nb14_calc_list(num_nb14,i2,pairlist(icel1,icel2)) = i1
                num_nb14_total(pairlist(icel1,icel2)) = &
                     num_nb14_total(pairlist(icel1,icel2)) + 1
                sc_calc_list(num_nb14,i2,pairlist(icel1,icel2)) = &
                  int(enefunc%dihe_periodicity(ix,i)/enefunc%notation_14types)
              end if

            else if (i1 < i2) then

              num_nb14 = num_nb14_calc1(i1,i)
              duplicate = .false.

              do k = 1, num_nb14
                if (i2 == nb14_calc_list1(k,i1,i)) &
                  duplicate = .true.
              end do

              do k = 1, num_nonb_excl1(i1,i)
                if (i2 == nonb_excl_list1(k,i1,i)) &
                  duplicate = .true.
              end do

              if (.not.duplicate) then
                num_nb14 = num_nb14 + 1
                num_nb14_calc1(i1,i) = num_nb14
                nb14_calc_list1(num_nb14,i1,i) = i2
                num_nb14_total1(i) = num_nb14_total1(i) + 1
                sc_calc_list1(num_nb14,i1,i) = &
                  int(enefunc%dihe_periodicity(ix,i)/enefunc%notation_14types)
              end if

            else if (i1 > i2) then

              num_nb14 = num_nb14_calc1(i2,i)
              duplicate = .false.

              do k = 1, num_nb14
                if (i1 == nb14_calc_list1(k,i2,i)) &
                  duplicate = .true.
              end do

              do k = 1, num_nonb_excl1(i2,i)
                if (i1 == nonb_excl_list1(k,i2,i)) &
                  duplicate = .true.
              end do

              if (.not.duplicate) then
                num_nb14 = num_nb14 + 1
                num_nb14_calc1(i2,i) = num_nb14
                nb14_calc_list1(num_nb14,i2,i) = i1
                num_nb14_total1(i) = num_nb14_total1(i) + 1
                sc_calc_list1(num_nb14,i2,i) = &
                   int(enefunc%dihe_periodicity(ix,i)/enefunc%notation_14types)
              end if

            end if

          end do

        end do

      end do

    end if

    ! Check the total number of exclusion list
    !
    if (first) then

      found1 = 0
      found2 = 0

      do icel = 1, maxcell_near
        found1 = found1 + num_excl_total(icel)
        found2 = found2 + num_nb14_total(icel)
      end do

      do icel = 1, ncell_local
        found1 = found1 + num_excl_total1(icel)
        found2 = found2 + num_nb14_total1(icel)
      end do
#ifdef HAVE_MPI_GENESIS
      call mpi_reduce(found1, enefunc%num_excl_all, 1, mpi_integer, mpi_sum, &
                      0, mpi_comm_country, ierror)
      call mpi_reduce(found2, enefunc%num_nb14_all, 1, mpi_integer, mpi_sum, &
                      0, mpi_comm_country, ierror)
#else
    enefunc%num_excl_all = found1
    enefunc%num_nb14_all = found2
#endif
    end if

    ! pack into small dimension
    !
    call pack_array_4to3(natom, domain%cell_pairlist1, num_nonb_excl, &
                     nonb_excl_list, enefunc%nonb_excl_list)

    call pack_array_4to3(natom, domain%cell_pairlist1, num_nb14_calc, &
                     nb14_calc_list, enefunc%nb14_calc_list)

    call pack_array_4to3(natom, domain%cell_pairlist1, num_nb14_calc, &
                     sc_calc_list, enefunc%sc_calc_list)

    call pack_array_3to2(natom, num_nonb_excl1, ncell_local,          &
                     nonb_excl_list1, enefunc%nonb_excl_list1)

    call pack_array_3to2(natom, num_nb14_calc1, ncell_local,          &
                     nb14_calc_list1, enefunc%nb14_calc_list1)

    call pack_array_3to2(natom, num_nb14_calc1, ncell_local,          &
                     sc_calc_list1, enefunc%sc_calc_list1)

   ! scnb/fudge_lj & scee/fudge_qq
   !
   if (enefunc%forcefield == ForcefieldAMBER .or. &
       enefunc%forcefield == ForcefieldGROAMBER .or. &
       enefunc%forcefield == ForcefieldGROMARTINI) then
      do i = 1, ncell_local
        num_nb14 = 0
        do ix = 1, natom(i) - 1
          ini_nb14  = num_nb14 + 1
          fin_nb14  = num_nb14 + enefunc%num_nb14_calc1(ix,i)
          num_nb14  = fin_nb14
          do k = ini_nb14, fin_nb14
            if (enefunc%forcefield == ForcefieldAMBER) then
              enefunc%nb14_lj_scale1(k, i) =  &
                        enefunc%dihe_scnb(enefunc%sc_calc_list1(k,i))
              enefunc%nb14_qq_scale1(k, i) =  &
                        enefunc%dihe_scee(enefunc%sc_calc_list1(k,i))
            else
              enefunc%nb14_lj_scale1(k, i) = enefunc%fudge_lj
              enefunc%nb14_qq_scale1(k, i) = enefunc%fudge_qq
            endif
          end do
        end do
      end do
      do ij = 1, maxcell_near
        i = cell_pair(1,ij)
        j = cell_pair(2,ij)
        num_nb14 = 0
        do ix = 1, natom(i)
          ini_nb14  = num_nb14 + 1
          fin_nb14  = num_nb14 + enefunc%num_nb14_calc(ix,ij)
          num_nb14  = fin_nb14
          do k = ini_nb14, fin_nb14
            if (enefunc%forcefield == ForcefieldAMBER) then
              enefunc%nb14_lj_scale(k, ij) =  &
                   enefunc%dihe_scnb(enefunc%sc_calc_list(k,ij))
              enefunc%nb14_qq_scale(k, ij) =  &
                   enefunc%dihe_scee(enefunc%sc_calc_list(k,ij))
            else
              enefunc%nb14_lj_scale(k, ij) = enefunc%fudge_lj
              enefunc%nb14_qq_scale(k, ij) = enefunc%fudge_qq
            endif
          end do
        end do
      end do
    end if

    return

  end subroutine count_nonb_excl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    count_nonb_excl_rest
  !> @brief        exclude 1-2, 1-3 interactions w/o local restraint
  !! @authors      JJ, CK
  !! @param[in]    first   : flag for first call or not
  !! @param[inout] domain  : structure of domain
  !! @param[inout] enefunc : structure of enefunc
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine count_nonb_excl_rest(first, domain, enefunc)

    ! formal arguments
    logical,                 intent(in)    :: first
    type(s_domain),  target, intent(inout) :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: ncell, i, j, ii, ij, ix, k, i1, i2
    integer                  :: index(4)
    integer                  :: icel, icel1, icel2, jcel
    integer                  :: num_excl, num_nb14
    integer                  :: ini_nb14, fin_nb14
    integer                  :: found1, found2
    integer                  :: fkind
    logical                  :: duplicate

    integer,         pointer :: ncell_local, max_atom, natom(:), id_g2l(:,:)
    integer,         pointer :: nwater(:), water_list(:,:,:)
    integer,         pointer :: cell_pair(:,:)
    integer,         pointer :: num_nonb_excl(:,:), num_nonb_excl1(:,:)
    integer,         pointer :: num_nb14_calc(:,:), num_nb14_calc1(:,:)
    integer,         pointer :: nbond(:), bondlist(:,:,:)
    integer,         pointer :: nangle(:), anglelist(:,:,:)
    integer,         pointer :: ndihedral(:), dihelist(:,:,:), pairlist(:,:)
    integer,         pointer :: nonb_excl_list(:,:,:), nonb_excl_list1(:,:,:)
    integer,         pointer :: nb14_calc_list(:,:,:), nb14_calc_list1(:,:,:)
    integer,         pointer :: num_excl_total(:), num_nb14_total(:)
    integer,         pointer :: num_excl_total1(:), num_nb14_total1(:)


    ncell_local     => domain%num_cell_local
    max_atom        => domain%max_num_atom
    natom           => domain%num_atom
    nwater          => domain%num_water
    water_list      => domain%water_list
    id_g2l          => domain%id_g2l
    cell_pair       => domain%cell_pair
    pairlist        => domain%cell_pairlist2

    nbond           => enefunc%num_bond
    bondlist        => enefunc%bond_list
    nangle          => enefunc%num_angle
    anglelist       => enefunc%angle_list
    ndihedral       => enefunc%num_dihedral
    dihelist        => enefunc%dihe_list
    nonb_excl_list  => enefunc%nonb_list
    nonb_excl_list1 => enefunc%nonb_list1
    nb14_calc_list  => enefunc%nb14_list
    nb14_calc_list1 => enefunc%nb14_list1
    num_excl_total  => enefunc%num_excl_total
    num_excl_total1 => enefunc%num_excl_total1
    num_nb14_total  => enefunc%num_nb14_total
    num_nb14_total1 => enefunc%num_nb14_total1
    num_nonb_excl   => enefunc%num_nonb_excl
    num_nonb_excl1  => enefunc%num_nonb_excl1
    num_nb14_calc   => enefunc%num_nb14_calc
    num_nb14_calc1  => enefunc%num_nb14_calc1

    ncell = domain%num_cell_local + domain%num_cell_boundary

    max_atom = 0
    do i = 1, ncell
      max_atom = max(max_atom,natom(i))
    end do

    ! initialization
    !
    num_nonb_excl  (1:max_atom,1:maxcell_near)= 0
    num_nonb_excl1 (1:max_atom,1:ncell_local) = 0
    num_nb14_calc  (1:max_atom,1:maxcell_near)= 0
    num_nb14_calc1 (1:max_atom,1:ncell_local) = 0
    num_excl_total (1:maxcell_near)           = 0
    num_excl_total1(1:ncell_local)            = 0
    num_nb14_total (1:maxcell_near)           = 0
    num_nb14_total1(1:ncell_local)            = 0

    ! exclude 1-2 interaction
    !
    if (enefunc%excl_level > 0) then

      do i = 1, ncell_local
        do ix = 1, nbond(i)

          icel1 = id_g2l(1,bondlist(1,ix,i))
          i1    = id_g2l(2,bondlist(1,ix,i))
          icel2 = id_g2l(1,bondlist(2,ix,i))
          i2    = id_g2l(2,bondlist(2,ix,i))
          fkind = enefunc%bond_kind(ix,i)

          if (fkind == 0) then

            if (icel1 < icel2) then
              num_excl = num_nonb_excl(i1,pairlist(icel2,icel1)) + 1
              num_nonb_excl(i1,pairlist(icel2,icel1)) = num_excl
              nonb_excl_list(num_excl,i1,pairlist(icel2,icel1)) = i2
              num_excl_total(pairlist(icel2,icel1)) = &
                   num_excl_total(pairlist(icel2,icel1)) + 1

            else if (icel1 > icel2) then
              num_excl = num_nonb_excl(i2,pairlist(icel1,icel2)) + 1
              num_nonb_excl(i2,pairlist(icel1,icel2)) = num_excl
              nonb_excl_list(num_excl,i2,pairlist(icel1,icel2)) = i1
              num_excl_total(pairlist(icel1,icel2)) = &
                   num_excl_total(pairlist(icel1,icel2)) + 1

            else if (i1 < i2) then
              num_excl = num_nonb_excl1(i1,i) + 1
              num_nonb_excl1(i1,i) = num_excl
              nonb_excl_list1(num_excl,i1,i) = i2
              num_excl_total1(i) = num_excl_total1(i) + 1

            else if (i1 > i2) then
              num_excl = num_nonb_excl1(i2,i) + 1
              num_nonb_excl1(i2,i) = num_excl
              nonb_excl_list1(num_excl,i2,i) = i1
              num_excl_total1(i) = num_excl_total1(i) + 1

            end if

          end if
        end do
      end do

    end if

    ! exclude 1-3 interaction
    !
    if (enefunc%excl_level > 1) then

      do i = 1, ncell_local
        do ix = 1, nangle(i)

          icel1 = id_g2l(1,anglelist(1,ix,i))
          i1    = id_g2l(2,anglelist(1,ix,i))
          icel2 = id_g2l(1,anglelist(3,ix,i))
          i2    = id_g2l(2,anglelist(3,ix,i))
          fkind = enefunc%angle_kind(ix,i)

          if (fkind == 0) then

            if (icel1 < icel2) then

              num_excl = num_nonb_excl(i1,pairlist(icel2,icel1))
              duplicate = .false.

              do k = 1, num_excl
                if (i2 == nonb_excl_list(k,i1,pairlist(icel2,icel1))) &
                  duplicate = .true.
              end do

              if (.not. duplicate) then
                num_excl = num_excl + 1
                num_nonb_excl(i1,pairlist(icel2,icel1)) = num_excl
                nonb_excl_list(num_excl,i1,pairlist(icel2,icel1)) = i2
                num_excl_total(pairlist(icel2,icel1)) = &
                     num_excl_total(pairlist(icel2,icel1)) + 1
              end if

            else if (icel1 > icel2) then

              num_excl = num_nonb_excl(i2,pairlist(icel1,icel2))
              duplicate = .false.

              do k = 1, num_excl
                if (i1 == nonb_excl_list(k,i2,pairlist(icel1,icel2))) &
                  duplicate = .true.
              end do

              if (.not.duplicate) then
                num_excl = num_excl + 1
                num_nonb_excl(i2,pairlist(icel1,icel2)) = num_excl
                nonb_excl_list(num_excl,i2,pairlist(icel1,icel2)) = i1
                num_excl_total(pairlist(icel1,icel2)) = &
                     num_excl_total(pairlist(icel1,icel2)) + 1
              end if

            else if (i1 < i2) then

              num_excl = num_nonb_excl1(i1,i)
              duplicate = .false.

              do k = 1, num_excl
                if (i2 == nonb_excl_list1(k,i1,i)) &
                  duplicate = .true.
              end do

              if (.not.duplicate) then
                num_excl = num_excl + 1
                num_nonb_excl1(i1,i) = num_excl
                nonb_excl_list1(num_excl,i1,i) = i2
                num_excl_total1(i) = num_excl_total1(i) + 1
              end if

            else if (i1 > i2) then

              num_excl = num_nonb_excl1(i2,i)
              duplicate = .false.

              do k = 1, num_excl
                if (i1 == nonb_excl_list1(k,i2,i)) &
                  duplicate = .true.
              end do

              if (.not.duplicate) then
                num_excl = num_excl + 1
                num_nonb_excl1(i2,i) = num_excl
                nonb_excl_list1(num_excl,i2,i) = i1
                num_excl_total1(i) = num_excl_total1(i) + 1
              end if

            end if

          end if

        end do
      end do

    end if

    ! exclude water
    !
    if (enefunc%excl_level > 1) then

      if (enefunc%table%tip4) then

        do i = 1, ncell_local
          do ix = 1, nwater(i)

            index(1:4) = water_list(1:4,ix,i)

            do i1 = 1, 3
              do i2 = i1+1, 4

                num_excl = num_nonb_excl1(index(i1),i)
                duplicate = .false.

                do k = 1, num_excl
                  if (index(i2) == nonb_excl_list1(k,index(i1),i)) &
                    duplicate = .true.
                end do

                if (.not. duplicate) then
                  num_excl = num_excl + 1
                  num_nonb_excl1(index(i1),i) = num_excl
                  nonb_excl_list1(num_excl,index(i1),i) = index(i2)
                  num_excl_total1(i) = num_excl_total1(i) + 1
                end if

              end do
            end do

          end do
        end do

      end if
    end if

    ! count 1-4 interaction
    !
    if (enefunc%excl_level > 2) then

      do ii = 1, 2

        if (ii == 2) then
          ndihedral => enefunc%num_rb_dihedral
          dihelist  => enefunc%rb_dihe_list
        end if

        do i = 1, ncell_local
          do ix = 1, ndihedral(i)

            icel1 = id_g2l(1,dihelist(1,ix,i))
            i1    = id_g2l(2,dihelist(1,ix,i))
            icel2 = id_g2l(1,dihelist(4,ix,i))
            i2    = id_g2l(2,dihelist(4,ix,i))
            fkind = enefunc%dihe_kind(ix,i)

            if (fkind == 0) then

              if (icel1 < icel2) then

                num_nb14 = num_nb14_calc(i1,pairlist(icel2,icel1))
                duplicate = .false.

                do k = 1, num_nonb_excl(i1,pairlist(icel2,icel1))
                  if (i2 == nonb_excl_list(k,i1,pairlist(icel2,icel1))) &
                    duplicate = .true.
                end do

                do k = 1, num_nb14
                  if (i2 == nb14_calc_list(k,i1,pairlist(icel2,icel1))) &
                    duplicate = .true.
                end do

                if (.not.duplicate) then
                  num_nb14 = num_nb14 + 1
                  num_nb14_calc(i1,pairlist(icel2,icel1)) = num_nb14
                  nb14_calc_list(num_nb14,i1,pairlist(icel2,icel1)) = i2
                  num_nb14_total(pairlist(icel2,icel1)) = &
                       num_nb14_total(pairlist(icel2,icel1)) + 1
                end if

              else if (icel1 > icel2) then

                num_nb14 = num_nb14_calc(i2,pairlist(icel1,icel2))
                duplicate = .false.

                do k = 1, num_nonb_excl(i2,pairlist(icel1,icel2))
                  if (i1 == nonb_excl_list(k,i2,pairlist(icel1,icel2))) &
                    duplicate = .true.
                end do

                do k = 1, num_nb14
                  if (i1 == nb14_calc_list(k,i2,pairlist(icel1,icel2))) &
                    duplicate = .true.
                end do

                if (.not.duplicate) then
                  num_nb14 = num_nb14 + 1
                  num_nb14_calc(i2,pairlist(icel1,icel2)) = num_nb14
                  nb14_calc_list(num_nb14,i2,pairlist(icel1,icel2)) = i1
                  num_nb14_total(pairlist(icel1,icel2)) = &
                       num_nb14_total(pairlist(icel1,icel2)) + 1
                end if

              else if (i1 < i2) then

                num_nb14 = num_nb14_calc1(i1,i)
                duplicate = .false.

                do k = 1, num_nb14
                  if (i2 == nb14_calc_list1(k,i1,i)) &
                    duplicate = .true.
                end do

                do k = 1, num_nonb_excl1(i1,i)
                  if (i2 == nonb_excl_list1(k,i1,i)) &
                    duplicate = .true.
                end do

                if (.not.duplicate) then
                  num_nb14 = num_nb14 + 1
                  num_nb14_calc1(i1,i) = num_nb14
                  nb14_calc_list1(num_nb14,i1,i) = i2
                  num_nb14_total1(i) = num_nb14_total1(i) + 1
                end if

              else if (i1 > i2) then

                num_nb14 = num_nb14_calc1(i2,i)
                duplicate = .false.

                do k = 1, num_nb14
                  if (i1 == nb14_calc_list1(k,i2,i)) &
                    duplicate = .true.
                end do

                do k = 1, num_nonb_excl1(i2,i)
                  if (i1 == nonb_excl_list1(k,i2,i)) &
                    duplicate = .true.
                end do

                if (.not.duplicate) then
                  num_nb14 = num_nb14 + 1
                  num_nb14_calc1(i2,i) = num_nb14
                  nb14_calc_list1(num_nb14,i2,i) = i1
                  num_nb14_total1(i) = num_nb14_total1(i) + 1
                end if

              end if

            end if

          end do
        end do
      end do

    end if

    ! Check the total number of exclusion list
    !
    if (first) then

      found1 = 0
      found2 = 0

      do icel = 1, maxcell_near
        found1 = found1 + num_excl_total(icel)
        found2 = found2 + num_nb14_total(icel)
      end do

      do icel = 1, ncell_local
        found1 = found1 + num_excl_total1(icel)
        found2 = found2 + num_nb14_total1(icel)
      end do
#ifdef HAVE_MPI_GENESIS
      call mpi_reduce(found1, enefunc%num_excl_all, 1, mpi_integer, mpi_sum, &
                      0, mpi_comm_country, ierror)
      call mpi_reduce(found2, enefunc%num_nb14_all, 1, mpi_integer, mpi_sum, &
                      0, mpi_comm_country, ierror)
#else
    enefunc%num_excl_all = found1
    enefunc%num_nb14_all = found2
#endif
    end if

    ! pack into small dimension
    !
    call pack_array_4to3(natom, domain%cell_pairlist1, num_nonb_excl, &
                     nonb_excl_list, enefunc%nonb_excl_list)

    call pack_array_4to3(natom, domain%cell_pairlist1, num_nb14_calc, &
                     nb14_calc_list, enefunc%nb14_calc_list)

    call pack_array_3to2(natom, num_nonb_excl1, ncell_local,          &
                     nonb_excl_list1, enefunc%nonb_excl_list1)

    call pack_array_3to2(natom, num_nb14_calc1, ncell_local,          &
                     nb14_calc_list1, enefunc%nb14_calc_list1)

   ! scnb/fudge_lj & scee/fudge_qq
   !
   if (enefunc%forcefield == ForcefieldAMBER .or. &
       enefunc%forcefield == ForcefieldGROAMBER .or. &
       enefunc%forcefield == ForcefieldGROMARTINI) then
      do i = 1, ncell_local
        num_nb14 = 0
        do ix = 1, natom(i) - 1
          ini_nb14  = num_nb14 + 1
          fin_nb14  = num_nb14 + enefunc%num_nb14_calc1(ix,i)
          num_nb14  = fin_nb14
          do k = ini_nb14, fin_nb14
            if (enefunc%forcefield == ForcefieldAMBER) then
              enefunc%nb14_lj_scale1(k, i) =  &
                        enefunc%dihe_scnb(enefunc%sc_calc_list1(k,i))
              enefunc%nb14_qq_scale1(k, i) =  &
                        enefunc%dihe_scee(enefunc%sc_calc_list1(k,i))
            else
              enefunc%nb14_lj_scale1(k, i) = enefunc%fudge_lj
              enefunc%nb14_qq_scale1(k, i) = enefunc%fudge_qq
            endif
          end do
        end do
      end do
      do ij = 1, maxcell_near
        i = cell_pair(1,ij)
        j = cell_pair(2,ij)
        num_nb14 = 0
        do ix = 1, natom(i)
          ini_nb14  = num_nb14 + 1
          fin_nb14  = num_nb14 + enefunc%num_nb14_calc(ix,ij)
          num_nb14  = fin_nb14
          do k = ini_nb14, fin_nb14
            if (enefunc%forcefield == ForcefieldAMBER) then
              enefunc%nb14_lj_scale(k, ij) =  &
                   enefunc%dihe_scnb(enefunc%sc_calc_list(k,ij))
              enefunc%nb14_qq_scale(k, ij) =  &
                   enefunc%dihe_scee(enefunc%sc_calc_list(k,ij))
            else
              enefunc%nb14_lj_scale(k, ij) = enefunc%fudge_lj
              enefunc%nb14_qq_scale(k, ij) = enefunc%fudge_qq
            endif
          end do
        end do
      end do
    end if

    return

  end subroutine count_nonb_excl_rest

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    count_nonb_excl_solute
  !> @brief        exclude 1-2, 1-3 interactions
  !! @authors      JJ
  !! @param[in]    first   : flag for first call or not
  !! @param[inout] domain  : structure of domain
  !! @param[inout] enefunc : structure of enefunc
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine count_nonb_excl_solute(first, domain, enefunc)

    ! formal arguments
    logical,                 intent(in)    :: first
    type(s_domain),  target, intent(inout) :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: ncell, ncell_local, i, ii, ix, k, i1, i2
    integer                  :: ij, j
    integer                  :: icel, icel1, icel2, jcel
    integer                  :: ini_nb14, fin_nb14
    integer                  :: num_excl, num_nb14, id, omp_get_thread_num
    integer                  :: found1, found2
    logical                  :: duplicate

    integer,         pointer :: nsolute(:), id_g2l(:,:), natom(:)
    integer,         pointer :: nbond(:), bond_list(:,:,:)
    integer,         pointer :: nangle(:), angl_list(:,:,:)
    integer,         pointer :: ndihedral(:), dihe_list(:,:,:)
    integer,         pointer :: num_nonb_excl(:,:), num_nonb_excl1(:,:)
    integer,         pointer :: num_nb14_calc(:,:), num_nb14_calc1(:,:)
    integer,         pointer :: cell_pair(:,:), cell_pairlist2(:,:)
    integer,         pointer :: nonb_excl_list(:,:,:), nonb_excl_list1(:,:,:)
    integer,         pointer :: nb14_calc_list(:,:,:), nb14_calc_list1(:,:,:)
    integer,         pointer :: sc_calc_list(:,:,:), sc_calc_list1(:,:,:)
    integer,         pointer :: num_excl_total(:), num_nb14_total(:)
    integer,         pointer :: num_excl_total1(:), num_nb14_total1(:)


    nsolute         => domain%num_solute
    natom           => domain%num_atom
    id_g2l          => domain%id_g2l
    cell_pair       => domain%cell_pairlist1
    cell_pairlist2  => domain%cell_pairlist2

    nbond           => enefunc%num_bond
    nangle          => enefunc%num_angle
    ndihedral       => enefunc%num_dihedral
    bond_list       => enefunc%bond_list
    angl_list       => enefunc%angle_list
    dihe_list       => enefunc%dihe_list
    num_nonb_excl   => enefunc%num_nonb_excl
    num_nonb_excl1  => enefunc%num_nonb_excl1
    num_nb14_calc   => enefunc%num_nb14_calc
    num_nb14_calc1  => enefunc%num_nb14_calc1
    nonb_excl_list  => enefunc%nonb_list
    nonb_excl_list1 => enefunc%nonb_list1
    nb14_calc_list  => enefunc%nb14_list
    nb14_calc_list1 => enefunc%nb14_list1
    sc_calc_list    => enefunc%sc_list
    sc_calc_list1   => enefunc%sc_list1
    num_excl_total  => enefunc%num_excl_total
    num_excl_total1 => enefunc%num_excl_total1
    num_nb14_total  => enefunc%num_nb14_total
    num_nb14_total1 => enefunc%num_nb14_total1

    ncell_local = domain%num_cell_local
    ncell       = domain%num_cell_local + domain%num_cell_boundary

    domain%max_num_atom = 0
    do i = 1, ncell
      domain%max_num_atom = max(domain%max_num_atom,domain%num_atom(i))
    end do

    ! initialization
    !
    num_nonb_excl  (1:domain%max_num_atom,1:maxcell_near)= 0
    num_nonb_excl1 (1:domain%max_num_atom,1:ncell_local) = 0
    num_nb14_calc  (1:domain%max_num_atom,1:maxcell_near)= 0
    num_nb14_calc1 (1:domain%max_num_atom,1:ncell_local) = 0
    num_excl_total (1:maxcell_near)                      = 0
    num_excl_total1(1:ncell_local)                       = 0
    num_nb14_total (1:maxcell_near)                      = 0
    num_nb14_total1(1:ncell_local)                       = 0

    ! exclude 1-2 interaction
    !
    !$omp parallel default(shared)                                             &
    !$omp private(id, i, ix, icel1, icel2, icel, i1, i2, num_excl, duplicate,  &
    !$omp          k, num_nb14)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    if (enefunc%excl_level > 0) then

      do i = 1, ncell_local
        do ix = 1, nbond(i)

          icel1 = id_g2l(1,bond_list(1,ix,i))
          icel2 = id_g2l(1,bond_list(2,ix,i))

          if (icel1 /= icel2) then

            icel  = cell_pairlist2(icel1,icel2)
            if (mod(icel,nthread) == id) then

              i1 = id_g2l(2,bond_list(1,ix,i))
              i2 = id_g2l(2,bond_list(2,ix,i))

              if (i1 <= nsolute(icel1) .and. i2 <= nsolute(icel2)) then
                if (icel1 < icel2) then
                  num_excl = num_nonb_excl(i1,icel) + 1
                  num_nonb_excl(i1,icel) = num_excl
                  nonb_excl_list(num_excl,i1,icel) = i2
                  num_excl_total(icel) = num_excl_total(icel) + 1

                else if (icel1 > icel2) then
                  num_excl = num_nonb_excl(i2,icel) + 1
                  num_nonb_excl(i2,icel) = num_excl
                  nonb_excl_list(num_excl,i2,icel) = i1
                  num_excl_total(icel) = num_excl_total(icel) + 1

                end if
              end if

            end if

          else

            if (mod(i,nthread) == id) then

              i1    = id_g2l(2,bond_list(1,ix,i))
              i2    = id_g2l(2,bond_list(2,ix,i))

              if (i1 <= nsolute(icel1) .and. i2 <= nsolute(icel2)) then
                if (i1 < i2) then

                  num_excl = num_nonb_excl1(i1,i) + 1
                  num_nonb_excl1(i1,i) = num_excl
                  nonb_excl_list1(num_excl,i1,i) = i2
                  num_excl_total1(i) = num_excl_total1(i) + 1

                else if (i1 > i2) then

                  num_excl = num_nonb_excl1(i2,i) + 1
                  num_nonb_excl1(i2,i) = num_excl
                  nonb_excl_list1(num_excl,i2,i) = i1
                  num_excl_total1(i) = num_excl_total1(i) + 1

                end if
              end if

            end if

          end if

        end do
      end do

    end if

    ! exclude 1-3 interaction
    !
    if (enefunc%excl_level > 1) then

      do i = 1, ncell_local
        do ix = 1, nangle(i)

          icel1 = id_g2l(1,angl_list(1,ix,i))
          icel2 = id_g2l(1,angl_list(3,ix,i))

          if (icel1 /= icel2) then

            icel  = cell_pairlist2(icel1,icel2)
            if (mod(icel,nthread) == id) then

              i1 = id_g2l(2,angl_list(1,ix,i))
              i2 = id_g2l(2,angl_list(3,ix,i))

              if (i1 <= nsolute(icel1) .and. i2 <= nsolute(icel2)) then
                if (icel1 < icel2) then

                  num_excl = num_nonb_excl(i1,icel)
                  duplicate = .false.

                  do k = 1, num_excl
                    if (i2 == nonb_excl_list(k,i1,icel)) &
                      duplicate = .true.
                  end do

                  if (.not. duplicate) then
                    num_excl = num_excl + 1
                    num_nonb_excl(i1,icel) = num_excl
                    nonb_excl_list(num_excl,i1,icel) = i2
                    num_excl_total(icel) =  num_excl_total(icel) + 1
                  end if

                else if (icel1 > icel2) then

                  num_excl = num_nonb_excl(i2,icel)
                  duplicate = .false.

                  do k = 1, num_excl
                    if (i1 == nonb_excl_list(k,i2,icel)) &
                      duplicate = .true.
                  end do

                  if (.not.duplicate) then
                    num_excl = num_excl + 1
                    num_nonb_excl(i2,icel) = num_excl
                    nonb_excl_list(num_excl,i2,icel) = i1
                    num_excl_total(icel) = num_excl_total(icel) + 1
                  end if

                end if
              end if

            end if

          else

            if (mod(i,nthread) == id) then

              i1    = id_g2l(2,angl_list(1,ix,i))
              i2    = id_g2l(2,angl_list(3,ix,i))

              if (i1 <= nsolute(icel1) .and. i2 <= nsolute(icel2)) then
                if (i1 < i2) then

                  num_excl = num_nonb_excl1(i1,i)
                  duplicate = .false.

                  do k = 1, num_excl
                    if (i2 == nonb_excl_list1(k,i1,i)) &
                      duplicate = .true.
                  end do

                  if (.not.duplicate) then
                    num_excl = num_excl + 1
                    num_nonb_excl1(i1,i) = num_excl
                    nonb_excl_list1(num_excl,i1,i) = i2
                    num_excl_total1(i) = num_excl_total1(i) + 1
                  end if

                else if (i1 > i2) then

                  num_excl = num_nonb_excl1(i2,i)
                  duplicate = .false.

                  do k = 1, num_excl
                    if (i1 == nonb_excl_list1(k,i2,i)) &
                      duplicate = .true.
                  end do

                  if (.not.duplicate) then
                    num_excl = num_excl + 1
                    num_nonb_excl1(i2,i) = num_excl
                    nonb_excl_list1(num_excl,i2,i) = i1
                    num_excl_total1(i) = num_excl_total1(i) + 1
                  end if

                end if
              end if

            end if
          end if

        end do
      end do

    end if

    !$omp end parallel

    ! count 1-4 interaction
    !
    if (enefunc%excl_level > 2) then

      do ii = 1, 2

        if (ii == 2) then
          ndihedral => enefunc%num_rb_dihedral
          dihe_list => enefunc%rb_dihe_list
        end if

        !$omp parallel default(shared)                       &
        !$omp private(id, i, ix, icel1, icel2, icel, i1, i2, &
        !$omp          duplicate, k, num_nb14)
        !

#ifdef OMP
        id = omp_get_thread_num()
#else
        id = 0
#endif

        do i = 1, ncell_local
          do ix = 1, ndihedral(i)

            icel1 = id_g2l(1,dihe_list(1,ix,i))
            icel2 = id_g2l(1,dihe_list(4,ix,i))

            if (icel1 /= icel2) then

              icel  = cell_pairlist2(icel1,icel2)
              if (mod(icel,nthread) == id) then

                i1 = id_g2l(2,dihe_list(1,ix,i))
                i2 = id_g2l(2,dihe_list(4,ix,i))

                if (icel1 < icel2) then

                  num_nb14 = num_nb14_calc(i1,icel)
                  duplicate = .false.

                  do k = 1, num_nonb_excl(i1,icel)
                    if (i2 == nonb_excl_list(k,i1,icel)) &
                      duplicate = .true.
                  end do

                  do k = 1, num_nb14
                    if (i2 == nb14_calc_list(k,i1,icel)) &
                      duplicate = .true.
                  end do

                  if (.not.duplicate) then
                    num_nb14 = num_nb14 + 1
                    num_nb14_calc(i1,icel) = num_nb14
                    nb14_calc_list(num_nb14,i1,icel) = i2
                    num_nb14_total(icel) = num_nb14_total(icel) + 1
                    sc_calc_list(num_nb14,i1,icel) = &
                   int(enefunc%dihe_periodicity(ix,i)/enefunc%notation_14types)
                  end if

                else if (icel1 > icel2) then

                  num_nb14 = num_nb14_calc(i2,icel)
                  duplicate = .false.

                  do k = 1, num_nonb_excl(i2,icel)
                    if (i1 == nonb_excl_list(k,i2,icel)) &
                      duplicate = .true.
                  end do

                  do k = 1, num_nb14
                    if (i1 == nb14_calc_list(k,i2,icel)) &
                      duplicate = .true.
                  end do

                  if (.not.duplicate) then
                    num_nb14 = num_nb14 + 1
                    num_nb14_calc(i2,icel) = num_nb14
                    nb14_calc_list(num_nb14,i2,icel) = i1
                    num_nb14_total(icel) = num_nb14_total(icel) + 1
                    sc_calc_list(num_nb14,i2,icel) = &
                    int(enefunc%dihe_periodicity(ix,i)/enefunc%notation_14types)
                  end if

                end if
              end if

            else

              if (mod(i,nthread) == id) then

                i1    = id_g2l(2,dihe_list(1,ix,i))
                i2    = id_g2l(2,dihe_list(4,ix,i))

                if (i1 < i2) then

                  num_nb14 = num_nb14_calc1(i1,i)
                  duplicate = .false.

                  do k = 1, num_nb14
                    if (i2 == nb14_calc_list1(k,i1,i)) &
                      duplicate = .true.
                  end do

                  do k = 1, num_nonb_excl1(i1,i)
                    if (i2 == nonb_excl_list1(k,i1,i)) &
                      duplicate = .true.
                  end do

                  if (.not.duplicate) then
                    num_nb14 = num_nb14 + 1
                    num_nb14_calc1(i1,i) = num_nb14
                    nb14_calc_list1(num_nb14,i1,i) = i2
                    num_nb14_total1(i) = num_nb14_total1(i) + 1
                    sc_calc_list1(num_nb14,i1,i) = &
                  int(enefunc%dihe_periodicity(ix,i)/enefunc%notation_14types)
                  end if

                else if (i1 > i2) then

                  num_nb14 = num_nb14_calc1(i2,i)
                  duplicate = .false.

                  do k = 1, num_nb14
                    if (i1 == nb14_calc_list1(k,i2,i)) &
                      duplicate = .true.
                  end do

                  do k = 1, num_nonb_excl1(i2,i)
                    if (i1 == nonb_excl_list1(k,i2,i)) &
                      duplicate = .true.
                  end do

                  if (.not.duplicate) then
                    num_nb14 = num_nb14 + 1
                    num_nb14_calc1(i2,i) = num_nb14
                    nb14_calc_list1(num_nb14,i2,i) = i1
                    num_nb14_total1(i) = num_nb14_total1(i) + 1
                    sc_calc_list1(num_nb14,i2,i) = &
                   int(enefunc%dihe_periodicity(ix,i)/enefunc%notation_14types)
                  end if

                end if
              end if

            end if

          end do
        end do

        !$omp end parallel
  
      end do

    end if

    ! Check the total number of exclusion list
    !
    if (first) then

      found1 = 0
      found2 = 0

      do icel = 1, maxcell_near
        found1 = found1 + num_excl_total(icel)
        found2 = found2 + num_nb14_total(icel)
      end do

      do icel = 1, ncell_local
        found1 = found1 + num_excl_total1(icel)
        found2 = found2 + num_nb14_total1(icel)
      end do

#ifdef HAVE_MPI_GENESIS
      call mpi_reduce(found1, enefunc%num_excl_all, 1, mpi_integer, mpi_sum, &
                      0, mpi_comm_country, ierror)
      call mpi_reduce(found2, enefunc%num_nb14_all, 1, mpi_integer, mpi_sum, &
                      0, mpi_comm_country, ierror)
#else
    enefunc%num_excl_all = found1
    enefunc%num_nb14_all = found2
#endif
    end if

    ! pack into small dimension
    !
    call pack_array_4to3(natom, domain%cell_pairlist1,       &
                     num_nonb_excl, nonb_excl_list, enefunc%nonb_excl_list)

    call pack_array_4to3(natom, domain%cell_pairlist1,       &
                     num_nb14_calc, nb14_calc_list, enefunc%nb14_calc_list)

    call pack_array_4to3(natom, domain%cell_pairlist1, num_nb14_calc,&
                     sc_calc_list, enefunc%sc_calc_list)

    call pack_array_3to2(natom, num_nonb_excl1, ncell_local, &
                     nonb_excl_list1, enefunc%nonb_excl_list1)

    call pack_array_3to2(natom, num_nb14_calc1, ncell_local, &
                     nb14_calc_list1, enefunc%nb14_calc_list1)

    call pack_array_3to2(natom, num_nb14_calc1, ncell_local,       &
                     sc_calc_list1, enefunc%sc_calc_list1)

   ! scnb/fudge_lj & scee/fudge_qq
   !
   if (enefunc%forcefield == ForcefieldAMBER .or. &
       enefunc%forcefield == ForcefieldGROAMBER .or. &
       enefunc%forcefield == ForcefieldGROMARTINI) then
      do i = 1, ncell_local
        num_nb14 = 0
        do ix = 1, natom(i) - 1
          ini_nb14  = num_nb14 + 1
          fin_nb14  = num_nb14 + enefunc%num_nb14_calc1(ix,i)
          num_nb14  = fin_nb14
          do k = ini_nb14, fin_nb14
            if (enefunc%forcefield == ForcefieldAMBER) then
              enefunc%nb14_lj_scale1(k, i) =  &
                        enefunc%dihe_scnb(enefunc%sc_calc_list1(k,i))
              enefunc%nb14_qq_scale1(k, i) =  &
                        enefunc%dihe_scee(enefunc%sc_calc_list1(k,i))
            else
              enefunc%nb14_lj_scale1(k, i) = enefunc%fudge_lj
              enefunc%nb14_qq_scale1(k, i) = enefunc%fudge_qq
            endif
          end do
        end do
      end do
      do ij = 1, maxcell_near
        i = cell_pair(1,ij)
        j = cell_pair(2,ij)
        num_nb14 = 0
        do ix = 1, natom(i)
          ini_nb14  = num_nb14 + 1
          fin_nb14  = num_nb14 + enefunc%num_nb14_calc(ix,ij)
          num_nb14  = fin_nb14
          do k = ini_nb14, fin_nb14
            if (enefunc%forcefield == ForcefieldAMBER) then
              enefunc%nb14_lj_scale(k, ij) =  &
                   enefunc%dihe_scnb(enefunc%sc_calc_list(k,ij))
              enefunc%nb14_qq_scale(k, ij) =  &
                   enefunc%dihe_scee(enefunc%sc_calc_list(k,ij))
            else
              enefunc%nb14_lj_scale(k, ij) = enefunc%fudge_lj
              enefunc%nb14_qq_scale(k, ij) = enefunc%fudge_qq
            endif
          end do
        end do
      end do
    end if

    return

  end subroutine count_nonb_excl_solute

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    count_nonb_excl_solute_rest
  !> @brief        exclude 1-2, 1-3 interactions w/o local restraint
  !! @authors      JJ
  !! @param[in]    first   : flag for first call or not
  !! @param[inout] domain  : structure of domain
  !! @param[inout] enefunc : structure of enefunc
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine count_nonb_excl_solute_rest(first, domain, enefunc)

    ! formal arguments
    logical,                 intent(in)    :: first
    type(s_domain),  target, intent(inout) :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: ncell, ncell_local, i, j, ii, ij, ix, k, i1, i2
    integer                  :: icel, icel1, icel2, jcel
    integer                  :: num_excl, num_nb14, id, omp_get_thread_num
    integer                  :: ini_nb14, fin_nb14
    integer                  :: found1, found2
    integer                  :: fkind
    logical                  :: duplicate

    integer,         pointer :: nsolute(:), id_g2l(:,:), natom(:)
    integer,         pointer :: nbond(:), bond_list(:,:,:)
    integer,         pointer :: nangle(:), angl_list(:,:,:)
    integer,         pointer :: ndihedral(:), dihe_list(:,:,:)
    integer,         pointer :: num_nonb_excl(:,:), num_nonb_excl1(:,:)
    integer,         pointer :: num_nb14_calc(:,:), num_nb14_calc1(:,:)
    integer,         pointer :: cell_pair(:,:), cell_pairlist2(:,:)
    integer,         pointer :: nonb_excl_list(:,:,:), nonb_excl_list1(:,:,:)
    integer,         pointer :: nb14_calc_list(:,:,:), nb14_calc_list1(:,:,:)
    integer,         pointer :: num_excl_total(:), num_nb14_total(:)
    integer,         pointer :: num_excl_total1(:), num_nb14_total1(:)


    nsolute         => domain%num_solute
    natom           => domain%num_atom
    id_g2l          => domain%id_g2l
    cell_pair       => domain%cell_pair
    cell_pairlist2  => domain%cell_pairlist2

    nbond           => enefunc%num_bond
    nangle          => enefunc%num_angle
    ndihedral       => enefunc%num_dihedral
    bond_list       => enefunc%bond_list
    angl_list       => enefunc%angle_list
    dihe_list       => enefunc%dihe_list
    num_nonb_excl   => enefunc%num_nonb_excl
    num_nonb_excl1  => enefunc%num_nonb_excl1
    num_nb14_calc   => enefunc%num_nb14_calc
    num_nb14_calc1  => enefunc%num_nb14_calc1
    nonb_excl_list  => enefunc%nonb_list
    nonb_excl_list1 => enefunc%nonb_list1
    nb14_calc_list  => enefunc%nb14_list
    nb14_calc_list1 => enefunc%nb14_list1
    num_excl_total  => enefunc%num_excl_total
    num_excl_total1 => enefunc%num_excl_total1
    num_nb14_total  => enefunc%num_nb14_total
    num_nb14_total1 => enefunc%num_nb14_total1

    ncell_local = domain%num_cell_local
    ncell       = domain%num_cell_local + domain%num_cell_boundary

    domain%max_num_atom = 0
    do i = 1, ncell
      domain%max_num_atom = max(domain%max_num_atom,domain%num_atom(i))
    end do

    ! initialization
    !
    num_nonb_excl  (1:domain%max_num_atom,1:maxcell_near)= 0
    num_nonb_excl1 (1:domain%max_num_atom,1:ncell_local) = 0
    num_nb14_calc  (1:domain%max_num_atom,1:maxcell_near)= 0
    num_nb14_calc1 (1:domain%max_num_atom,1:ncell_local) = 0
    num_excl_total (1:maxcell_near)                      = 0
    num_excl_total1(1:ncell_local)                       = 0
    num_nb14_total (1:maxcell_near)                      = 0
    num_nb14_total1(1:ncell_local)                       = 0

    ! exclude 1-2 interaction
    !
    !$omp parallel default(shared)                                             &
    !$omp private(id, i, ix, icel1, icel2, icel, i1, i2, num_excl, duplicate,  &
    !$omp          k, num_nb14, fkind)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    if (enefunc%excl_level > 0) then

      do i = 1, ncell_local
        do ix = 1, nbond(i)

          icel1 = id_g2l(1,bond_list(1,ix,i))
          icel2 = id_g2l(1,bond_list(2,ix,i))
          fkind = enefunc%bond_kind(ix,i)

          if (fkind == 0) then

            if (icel1 /= icel2) then

              icel  = cell_pairlist2(icel1,icel2)
              if (mod(icel,nthread) == id) then

                i1 = id_g2l(2,bond_list(1,ix,i))
                i2 = id_g2l(2,bond_list(2,ix,i))

                if (i1 <= nsolute(icel1) .and. i2 <= nsolute(icel2)) then
                  if (icel1 < icel2) then
                    num_excl = num_nonb_excl(i1,icel) + 1
                    num_nonb_excl(i1,icel) = num_excl
                    nonb_excl_list(num_excl,i1,icel) = i2
                    num_excl_total(icel) = num_excl_total(icel) + 1

                  else if (icel1 > icel2) then
                    num_excl = num_nonb_excl(i2,icel) + 1
                    num_nonb_excl(i2,icel) = num_excl
                    nonb_excl_list(num_excl,i2,icel) = i1
                    num_excl_total(icel) = num_excl_total(icel) + 1

                  end if
                end if

              end if

            else

              if (mod(i,nthread) == id) then

                i1    = id_g2l(2,bond_list(1,ix,i))
                i2    = id_g2l(2,bond_list(2,ix,i))

                if (i1 <= nsolute(icel1) .and. i2 <= nsolute(icel2)) then
                  if (i1 < i2) then

                    num_excl = num_nonb_excl1(i1,i) + 1
                    num_nonb_excl1(i1,i) = num_excl
                    nonb_excl_list1(num_excl,i1,i) = i2
                    num_excl_total1(i) = num_excl_total1(i) + 1

                  else if (i1 > i2) then

                    num_excl = num_nonb_excl1(i2,i) + 1
                    num_nonb_excl1(i2,i) = num_excl
                    nonb_excl_list1(num_excl,i2,i) = i1
                    num_excl_total1(i) = num_excl_total1(i) + 1

                  end if
                end if

              end if

            end if

          end if

        end do
      end do

    end if

    ! exclude 1-3 interaction
    !
    if (enefunc%excl_level > 1) then

      do i = 1, ncell_local
        do ix = 1, nangle(i)

          icel1 = id_g2l(1,angl_list(1,ix,i))
          icel2 = id_g2l(1,angl_list(3,ix,i))
          fkind = enefunc%angle_kind(ix,i)

          if (fkind == 0) then

            if (icel1 /= icel2) then

              icel  = cell_pairlist2(icel1,icel2)
              if (mod(icel,nthread) == id) then

                i1 = id_g2l(2,angl_list(1,ix,i))
                i2 = id_g2l(2,angl_list(3,ix,i))

                if (i1 <= nsolute(icel1) .and. i2 <= nsolute(icel2)) then
                  if (icel1 < icel2) then

                    num_excl = num_nonb_excl(i1,icel)
                    duplicate = .false.

                    do k = 1, num_excl
                      if (i2 == nonb_excl_list(k,i1,icel)) &
                        duplicate = .true.
                    end do

                    if (.not. duplicate) then
                      num_excl = num_excl + 1
                      num_nonb_excl(i1,icel) = num_excl
                      nonb_excl_list(num_excl,i1,icel) = i2
                      num_excl_total(icel) =  num_excl_total(icel) + 1
                    end if

                  else if (icel1 > icel2) then

                    num_excl = num_nonb_excl(i2,icel)
                    duplicate = .false.

                    do k = 1, num_excl
                      if (i1 == nonb_excl_list(k,i2,icel)) &
                        duplicate = .true.
                    end do

                    if (.not.duplicate) then
                      num_excl = num_excl + 1
                      num_nonb_excl(i2,icel) = num_excl
                      nonb_excl_list(num_excl,i2,icel) = i1
                      num_excl_total(icel) = num_excl_total(icel) + 1
                    end if

                  end if
                end if

              end if

            else

              if (mod(i,nthread) == id) then

                i1    = id_g2l(2,angl_list(1,ix,i))
                i2    = id_g2l(2,angl_list(3,ix,i))

                if (i1 <= nsolute(icel1) .and. i2 <= nsolute(icel2)) then
                  if (i1 < i2) then

                    num_excl = num_nonb_excl1(i1,i)
                    duplicate = .false.

                    do k = 1, num_excl
                      if (i2 == nonb_excl_list1(k,i1,i)) &
                        duplicate = .true.
                    end do

                    if (.not.duplicate) then
                      num_excl = num_excl + 1
                      num_nonb_excl1(i1,i) = num_excl
                      nonb_excl_list1(num_excl,i1,i) = i2
                      num_excl_total1(i) = num_excl_total1(i) + 1
                    end if

                  else if (i1 > i2) then

                    num_excl = num_nonb_excl1(i2,i)
                    duplicate = .false.

                    do k = 1, num_excl
                      if (i1 == nonb_excl_list1(k,i2,i)) &
                        duplicate = .true.
                    end do

                    if (.not.duplicate) then
                      num_excl = num_excl + 1
                      num_nonb_excl1(i2,i) = num_excl
                      nonb_excl_list1(num_excl,i2,i) = i1
                      num_excl_total1(i) = num_excl_total1(i) + 1
                    end if

                  end if
                end if

              end if
            end if

          end if

        end do
      end do

    end if

    !$omp end parallel

    ! count 1-4 interaction
    !
    if (enefunc%excl_level > 2) then

      do ii = 1, 2

        if (ii == 2) then
          ndihedral => enefunc%num_rb_dihedral
          dihe_list => enefunc%rb_dihe_list
        end if

        !$omp parallel default(shared)                       &
        !$omp private(id, i, ix, icel1, icel2, icel, i1, i2, &
        !$omp          duplicate, k, num_nb14, fkind)
        !

#ifdef OMP
        id = omp_get_thread_num()
#else
        id = 0
#endif

        do i = 1, ncell_local
          do ix = 1, ndihedral(i)

            icel1 = id_g2l(1,dihe_list(1,ix,i))
            icel2 = id_g2l(1,dihe_list(4,ix,i))
            fkind = enefunc%dihe_kind(ix,i)

            if (fkind == 0) then
              if (icel1 /= icel2) then

                icel  = cell_pairlist2(icel1,icel2)
                if (mod(icel,nthread) == id) then

                  i1 = id_g2l(2,dihe_list(1,ix,i))
                  i2 = id_g2l(2,dihe_list(4,ix,i))

                  if (icel1 < icel2) then

                    num_nb14 = num_nb14_calc(i1,icel)
                    duplicate = .false.

                    do k = 1, num_nonb_excl(i1,icel)
                      if (i2 == nonb_excl_list(k,i1,icel)) &
                        duplicate = .true.
                    end do

                    do k = 1, num_nb14
                      if (i2 == nb14_calc_list(k,i1,icel)) &
                        duplicate = .true.
                    end do

                    if (.not.duplicate) then
                      num_nb14 = num_nb14 + 1
                      num_nb14_calc(i1,icel) = num_nb14
                      nb14_calc_list(num_nb14,i1,icel) = i2
                      num_nb14_total(icel) = num_nb14_total(icel) + 1
                    end if

                  else if (icel1 > icel2) then

                    num_nb14 = num_nb14_calc(i2,icel)
                    duplicate = .false.

                    do k = 1, num_nonb_excl(i2,icel)
                      if (i1 == nonb_excl_list(k,i2,icel)) &
                        duplicate = .true.
                    end do

                    do k = 1, num_nb14
                      if (i1 == nb14_calc_list(k,i2,icel)) &
                        duplicate = .true.
                    end do

                    if (.not.duplicate) then
                      num_nb14 = num_nb14 + 1
                      num_nb14_calc(i2,icel) = num_nb14
                      nb14_calc_list(num_nb14,i2,icel) = i1
                      num_nb14_total(icel) = num_nb14_total(icel) + 1
                    end if

                  end if
                end if

              else

                if (mod(i,nthread) == id) then

                  i1    = id_g2l(2,dihe_list(1,ix,i))
                  i2    = id_g2l(2,dihe_list(4,ix,i))

                  if (i1 < i2) then

                    num_nb14 = num_nb14_calc1(i1,i)
                    duplicate = .false.

                    do k = 1, num_nb14
                      if (i2 == nb14_calc_list1(k,i1,i)) &
                        duplicate = .true.
                    end do

                    do k = 1, num_nonb_excl1(i1,i)
                      if (i2 == nonb_excl_list1(k,i1,i)) &
                        duplicate = .true.
                    end do

                    if (.not.duplicate) then
                      num_nb14 = num_nb14 + 1
                      num_nb14_calc1(i1,i) = num_nb14
                      nb14_calc_list1(num_nb14,i1,i) = i2
                      num_nb14_total1(i) = num_nb14_total1(i) + 1
                    end if

                  else if (i1 > i2) then

                    num_nb14 = num_nb14_calc1(i2,i)
                    duplicate = .false.

                    do k = 1, num_nb14
                      if (i1 == nb14_calc_list1(k,i2,i)) &
                        duplicate = .true.
                    end do

                    do k = 1, num_nonb_excl1(i2,i)
                      if (i1 == nonb_excl_list1(k,i2,i)) &
                        duplicate = .true.
                    end do

                    if (.not.duplicate) then
                      num_nb14 = num_nb14 + 1
                      num_nb14_calc1(i2,i) = num_nb14
                      nb14_calc_list1(num_nb14,i2,i) = i1
                      num_nb14_total1(i) = num_nb14_total1(i) + 1
                    end if

                  end if
                end if

              end if

            end if

          end do
        end do

        !$omp end parallel

      end do

    end if


    ! Check the total number of exclusion list
    !
    if (first) then

      found1 = 0
      found2 = 0

      do icel = 1, maxcell_near
        found1 = found1 + num_excl_total(icel)
        found2 = found2 + num_nb14_total(icel)
      end do

      do icel = 1, ncell_local
        found1 = found1 + num_excl_total1(icel)
        found2 = found2 + num_nb14_total1(icel)
      end do

#ifdef HAVE_MPI_GENESIS
      call mpi_reduce(found1, enefunc%num_excl_all, 1, mpi_integer, mpi_sum, &
                      0, mpi_comm_country, ierror)
      call mpi_reduce(found2, enefunc%num_nb14_all, 1, mpi_integer, mpi_sum, &
                      0, mpi_comm_country, ierror)
#else
    enefunc%num_excl_all = found1
    enefunc%num_nb14_all = found2
#endif
    end if

    ! pack into small dimension
    !
    call pack_array_4to3(natom, domain%cell_pairlist1,       &
                     num_nonb_excl, nonb_excl_list, enefunc%nonb_excl_list)

    call pack_array_4to3(natom, domain%cell_pairlist1,       &
                     num_nb14_calc, nb14_calc_list, enefunc%nb14_calc_list)

    call pack_array_3to2(natom, num_nonb_excl1, ncell_local, &
                     nonb_excl_list1, enefunc%nonb_excl_list1)

    call pack_array_3to2(natom, num_nb14_calc1, ncell_local, &
                     nb14_calc_list1, enefunc%nb14_calc_list1)

   ! scnb/fudge_lj & scee/fudge_qq
   !
   if (enefunc%forcefield == ForcefieldAMBER .or. &
       enefunc%forcefield == ForcefieldGROAMBER .or. &
       enefunc%forcefield == ForcefieldGROMARTINI) then
      do i = 1, ncell_local
        num_nb14 = 0
        do ix = 1, natom(i) - 1
          ini_nb14  = num_nb14 + 1
          fin_nb14  = num_nb14 + enefunc%num_nb14_calc1(ix,i)
          num_nb14  = fin_nb14
          do k = ini_nb14, fin_nb14
            if (enefunc%forcefield == ForcefieldAMBER) then
              enefunc%nb14_lj_scale1(k, i) =  &
                        enefunc%dihe_scnb(enefunc%sc_calc_list1(k,i))
              enefunc%nb14_qq_scale1(k, i) =  &
                        enefunc%dihe_scee(enefunc%sc_calc_list1(k,i))
            else
              enefunc%nb14_lj_scale1(k, i) = enefunc%fudge_lj
              enefunc%nb14_qq_scale1(k, i) = enefunc%fudge_qq
            endif
          end do
        end do
      end do
      do ij = 1, maxcell_near
        i = cell_pair(1,ij)
        j = cell_pair(2,ij)
        num_nb14 = 0
        do ix = 1, natom(i)
          ini_nb14  = num_nb14 + 1
          fin_nb14  = num_nb14 + enefunc%num_nb14_calc(ix,ij)
          num_nb14  = fin_nb14
          do k = ini_nb14, fin_nb14
            if (enefunc%forcefield == ForcefieldAMBER) then
              enefunc%nb14_lj_scale(k, ij) =  &
                   enefunc%dihe_scnb(enefunc%sc_calc_list(k,ij))
              enefunc%nb14_qq_scale(k, ij) =  &
                   enefunc%dihe_scee(enefunc%sc_calc_list(k,ij))
            else
              enefunc%nb14_lj_scale(k, ij) = enefunc%fudge_lj
              enefunc%nb14_qq_scale(k, ij) = enefunc%fudge_qq
            endif
          end do
        end do
      end do
    end if

    return

  end subroutine count_nonb_excl_solute_rest

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    count_nonb_excl_constraint
  !> @brief        exclude 1-2, 1-3 interactions and constraints
  !! @authors      JJ
  !! @param[in]    first       : flag for first call or not
  !! @param[in]    water_table : flag for water table
  !! @param[in]    constraints : constraints information
  !! @param[inout] domain      : structure of domain
  !! @param[inout] enefunc     : structure of enefunc
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine count_nonb_excl_constraint(first, water_table, constraints, &
                                        domain, enefunc)

    ! formal arguments
    logical,                     intent(in)    :: first
    logical,                     intent(in)    :: water_table
    type(s_constraints), target, intent(in)    :: constraints
    type(s_domain),      target, intent(inout) :: domain
    type(s_enefunc),     target, intent(inout) :: enefunc

    ! local variables
    integer                  :: ncell, ncell_local, i, ii, ix, k, i1, i2, i3
    integer                  :: icel, icel1, icel2, jcel
    integer                  :: ic, j, ih, ij
    integer                  :: ini_nb14, fin_nb14
    integer                  :: num_excl, num_nb14, id, omp_get_thread_num
    integer                  :: index(4)
    integer                  :: found1, found2
    logical                  :: duplicate

    integer,         pointer :: nsolute(:), nwater(:), water_list(:,:,:)
    integer,         pointer :: id_g2l(:,:), natom(:)
    integer,         pointer :: nbond(:), bond_list(:,:,:)
    integer,         pointer :: nangle(:), angl_list(:,:,:)
    integer,         pointer :: ndihedral(:), dihe_list(:,:,:)
    integer,         pointer :: num_nonb_excl(:,:), num_nonb_excl1(:,:)
    integer,         pointer :: num_nb14_calc(:,:), num_nb14_calc1(:,:)
    integer,         pointer :: cell_pair(:,:), cell_pairlist2(:,:)
    integer,         pointer :: nonb_excl_list(:,:,:), nonb_excl_list1(:,:,:)
    integer,         pointer :: nb14_calc_list(:,:,:), nb14_calc_list1(:,:,:)
    integer,         pointer :: sc_calc_list(:,:,:), sc_calc_list1(:,:,:)
    integer,         pointer :: num_excl_total(:), num_nb14_total(:)
    integer,         pointer :: num_excl_total1(:), num_nb14_total1(:)
    integer,         pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)

    ! FEP
    integer                  :: iA, iB, iatomA, iatomB
    integer                  :: fg1, fg2

    HGr_local       => constraints%HGr_local
    HGr_bond_list   => constraints%HGr_bond_list

    nsolute         => domain%num_solute
    nwater          => domain%num_water
    natom           => domain%num_atom
    water_list      => domain%water_list
    id_g2l          => domain%id_g2l
    cell_pair       => domain%cell_pairlist1
    cell_pairlist2  => domain%cell_pairlist2

    nbond           => enefunc%num_bond
    nangle          => enefunc%num_angle
    ndihedral       => enefunc%num_dihedral
    bond_list       => enefunc%bond_list
    angl_list       => enefunc%angle_list
    dihe_list       => enefunc%dihe_list
    num_nonb_excl   => enefunc%num_nonb_excl
    num_nonb_excl1  => enefunc%num_nonb_excl1
    num_nb14_calc   => enefunc%num_nb14_calc
    num_nb14_calc1  => enefunc%num_nb14_calc1
    nonb_excl_list  => enefunc%nonb_list
    nonb_excl_list1 => enefunc%nonb_list1
    nb14_calc_list  => enefunc%nb14_list
    nb14_calc_list1 => enefunc%nb14_list1
    sc_calc_list    => enefunc%sc_list
    sc_calc_list1   => enefunc%sc_list1
    num_excl_total  => enefunc%num_excl_total
    num_excl_total1 => enefunc%num_excl_total1
    num_nb14_total  => enefunc%num_nb14_total
    num_nb14_total1 => enefunc%num_nb14_total1

    ncell_local = domain%num_cell_local
    ncell       = domain%num_cell_local + domain%num_cell_boundary

    domain%max_num_atom = 0
    do i = 1, ncell
      domain%max_num_atom = max(domain%max_num_atom,domain%num_atom(i))
    end do

    ! initialization
    !
    num_nonb_excl  (1:domain%max_num_atom,1:maxcell_near)= 0
    num_nonb_excl1 (1:domain%max_num_atom,1:ncell_local) = 0
    num_nb14_calc  (1:domain%max_num_atom,1:maxcell_near)= 0
    num_nb14_calc1 (1:domain%max_num_atom,1:ncell_local) = 0
    num_excl_total (1:maxcell_near)                      = 0
    num_excl_total1(1:ncell_local)                       = 0
    num_nb14_total (1:maxcell_near)                      = 0
    num_nb14_total1(1:ncell_local)                       = 0

    ! exclude 1-2 interaction
    !
    !$omp parallel default(shared)                                     &
    !$omp private(id, i, ix, icel1, icel2, icel, i1, i2, i3, num_excl, &
    !$omp         duplicate, k, num_nb14, ic, j, ih, index,            &
    !$omp         iA, iB, iatomA, iatomB, fg1, fg2)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    if (enefunc%excl_level > 0) then

      do i = 1, ncell_local
        do ix = 1, nbond(i)

          icel1 = id_g2l(1,bond_list(1,ix,i))
          icel2 = id_g2l(1,bond_list(2,ix,i))

          if (icel1 /= icel2) then

            icel  = cell_pairlist2(icel1,icel2)
            if (mod(icel-1,nthread) == id) then

              i1 = id_g2l(2,bond_list(1,ix,i))
              i2 = id_g2l(2,bond_list(2,ix,i))

              if (i1 <= nsolute(icel1) .and. i2 <= nsolute(icel2)) then
                if (icel1 < icel2) then
                  num_excl = num_nonb_excl(i1,icel) + 1
                  num_nonb_excl(i1,icel) = num_excl
                  nonb_excl_list(num_excl,i1,icel) = i2
                  num_excl_total(icel) = num_excl_total(icel) + 1

                else if (icel1 > icel2) then
                  num_excl = num_nonb_excl(i2,icel) + 1
                  num_nonb_excl(i2,icel) = num_excl
                  nonb_excl_list(num_excl,i2,icel) = i1
                  num_excl_total(icel) = num_excl_total(icel) + 1

                end if
              end if
            end if

          else

            if (mod(i-1,nthread) == id) then

              i1 = id_g2l(2,bond_list(1,ix,i))
              i2 = id_g2l(2,bond_list(2,ix,i))

              if (i1 <= nsolute(icel1) .and. i2 <= nsolute(icel2)) then
                if (i1 < i2) then
                  num_excl = num_nonb_excl1(i1,i) + 1
                  num_nonb_excl1(i1,i) = num_excl
                  nonb_excl_list1(num_excl,i1,i) = i2
                  num_excl_total1(i) = num_excl_total1(i) + 1

                else if (i1 > i2) then
                  num_excl = num_nonb_excl1(i2,i) + 1
                  num_nonb_excl1(i2,i) = num_excl
                  nonb_excl_list1(num_excl,i2,i) = i1
                  num_excl_total1(i) = num_excl_total1(i) + 1

                end if
              end if

            end if

          end if

        end do
      end do

      ! exclude constraint
      !
      do icel = id+1, ncell_local, nthread
        do ic = 1, constraints%connect
          do j = 1, HGr_local(ic,icel)

            i1 = HGr_bond_list(1,j,ic,icel)
            do ih = 1, ic
              i2 = HGr_bond_list(ih+1,j,ic,icel)

              if (domain%fep_use) then
                ! FEP: Hydrogen rewiring
                ! To avoid SHAKE problem in FEP, singleB-dualB bonds is virtually
                ! rewired as singleA-dualB. In this case, singleA-dualB bonds are
                ! excluded, but singleB-dualB bonds are not excluded.
                ! To exclude singleB-dualB bonds, for singleA-dualB bond
                ! including hydrogen, the atom index of singleA is replaced
                ! with the corresponding atom index of singleB.
                fg1 = domain%fepgrp(i1,icel)
                fg2 = domain%fepgrp(i2,icel)
                if ((fg1 == 1) .and. (fg2 == 4)) then
                  do i = 1, domain%num_atom_singleB(icel)
                    iA = domain%id_singleA(i, icel, 3)
                    iB = domain%id_singleB(i, icel, 3)
                    iatomA = domain%id_singleA(iA, icel, 1)
                    iatomB = domain%id_singleB(iB, icel, 1)
                    if (iatomA == i1) then
                      i1 = iatomB
                      exit
                    end if
                  end do
                else if ((fg1 == 4) .and. (fg2 == 1)) then
                  do i = 1, domain%num_atom_singleB(icel)
                    iA = domain%id_singleA(i, icel, 3)
                    iB = domain%id_singleB(i, icel, 3)
                    iatomA = domain%id_singleA(iA, icel, 1)
                    iatomB = domain%id_singleB(iB, icel, 1)
                    if (iatomA == i2) then
                      i2 = iatomB
                      exit
                    end if
                  end do
                end if
              end if

              num_excl = num_nonb_excl1(i1,icel) + 1
              num_nonb_excl1(i1,icel) = num_excl
              nonb_excl_list1(num_excl,i1,icel) = i2
              num_excl_total1(icel) = num_excl_total1(icel) + 1

            end do

          end do
        end do
      end do

    end if

    ! exclude water
    !
    if (enefunc%excl_level > 1) then

      if (constraints%tip4) then

        do icel = id+1, ncell_local, nthread
          do ix = 1, nwater(icel)

            index(1:4) = water_list(1:4,ix,icel)

            do i1 = 1, 3
              do i2 = i1+1, 4

                num_excl = num_nonb_excl1(index(i1),icel)
                duplicate = .false.

                do k = 1, num_excl
                  if (index(i2) == nonb_excl_list1(k,index(i1),icel)) &
                    duplicate = .true.
                end do

                if (.not. duplicate) then
                  num_excl = num_excl + 1
                  num_nonb_excl1(index(i1),icel) = num_excl
                  nonb_excl_list1(num_excl,index(i1),icel) = index(i2)
                  num_excl_total1(icel) = num_excl_total1(icel) + 1
                end if

              end do
            end do

          end do
        end do

      else

        do icel = id+1, ncell_local, nthread
          do ic = 1, nwater(icel)

            index(1:3) = water_list(1:3,ic,icel)

            do i1 = 1, 2
              do i2 = i1+1, 3

                num_excl = num_nonb_excl1(index(i1),icel)
                duplicate = .false.

                do k = 1, num_excl
                  if (index(i2) == nonb_excl_list1(k,index(i1),icel)) &
                    duplicate = .true.
                end do

                if (.not. duplicate) then
                  num_excl = num_excl + 1
                  num_nonb_excl1(index(i1),icel) = num_excl
                  nonb_excl_list1(num_excl,index(i1),icel) = index(i2)
                  num_excl_total1(icel) = num_excl_total1(icel) + 1
                end if

              end do
            end do

          end do
        end do

      end if

      !$mpi barrier

      ! exclude 1-3 interaction
      !
      do i = 1, ncell_local
        do ix = 1, nangle(i)

          icel1 = id_g2l(1,angl_list(1,ix,i))
          icel2 = id_g2l(1,angl_list(3,ix,i))

          if (icel1 /= icel2) then

            icel  = cell_pairlist2(icel1,icel2)
            if (mod(icel-1,nthread) == id) then

              i1 = id_g2l(2,angl_list(1,ix,i))
              i2 = id_g2l(2,angl_list(3,ix,i))

              if (i1 <= nsolute(icel1) .and. i2 <= nsolute(icel2)) then
                if (icel1 < icel2) then

                  num_excl = num_nonb_excl(i1,icel)
                  duplicate = .false.

                  do k = 1, num_excl
                    if (i2 == nonb_excl_list(k,i1,icel)) &
                      duplicate = .true.
                  end do

                  if (.not. duplicate) then
                    num_excl = num_excl + 1
                    num_nonb_excl(i1,icel) = num_excl
                    nonb_excl_list(num_excl,i1,icel) = i2
                    num_excl_total(icel) =  num_excl_total(icel) + 1
                  end if

                else if (icel1 > icel2) then

                  num_excl = num_nonb_excl(i2,icel)
                  duplicate = .false.

                  do k = 1, num_excl
                    if (i1 == nonb_excl_list(k,i2,icel)) &
                      duplicate = .true.
                  end do

                  if (.not.duplicate) then
                    num_excl = num_excl + 1
                    num_nonb_excl(i2,icel) = num_excl
                    nonb_excl_list(num_excl,i2,icel) = i1
                    num_excl_total(icel) = num_excl_total(icel) + 1
                  end if

                end if
              end if
            end if

          else

            if (mod(i-1,nthread) == id) then

              i1 = id_g2l(2,angl_list(1,ix,i))
              i2 = id_g2l(2,angl_list(3,ix,i))

              if (i1 <= nsolute(icel1) .and. i2 <= nsolute(icel2)) then
                if (i1 < i2) then

                  num_excl = num_nonb_excl1(i1,i)
                  duplicate = .false.

                  do k = 1, num_excl
                    if (i2 == nonb_excl_list1(k,i1,i)) &
                      duplicate = .true.
                  end do

                  if (.not.duplicate) then
                    num_excl = num_excl + 1
                    num_nonb_excl1(i1,i) = num_excl
                    nonb_excl_list1(num_excl,i1,i) = i2
                    num_excl_total1(i) = num_excl_total1(i) + 1
                  end if

                else if (i1 > i2) then

                  num_excl = num_nonb_excl1(i2,i)
                  duplicate = .false.

                  do k = 1, num_excl
                    if (i1 == nonb_excl_list1(k,i2,i)) &
                      duplicate = .true.
                  end do

                  if (.not.duplicate) then
                    num_excl = num_excl + 1
                    num_nonb_excl1(i2,i) = num_excl
                    nonb_excl_list1(num_excl,i2,i) = i1
                    num_excl_total1(i) = num_excl_total1(i) + 1
                  end if

                end if
              end if

            end if
          end if

        end do
      end do

    end if

    !$omp end parallel

    ! count 1-4 interaction
    !
    if (enefunc%excl_level > 2) then

      do ii = 1, 2

        if (ii == 2) then
          ndihedral => enefunc%num_rb_dihedral
          dihe_list => enefunc%rb_dihe_list
        end if

        !$omp parallel default(shared)                           &
        !$omp private(id, i, ix, icel1, icel2, icel, i1, i2, i3, &
        !$omp          duplicate, k, num_nb14, ic, j, ih)
        !

#ifdef OMP
        id = omp_get_thread_num()
#else
        id = 0
#endif

        do i = 1, ncell_local
          do ix = 1, ndihedral(i)

            icel1 = id_g2l(1,dihe_list(1,ix,i))
            icel2 = id_g2l(1,dihe_list(4,ix,i))

            if (icel1 /= icel2) then

              icel  = cell_pairlist2(icel1,icel2)
              if (mod(icel-1,nthread) == id) then

                i1 = id_g2l(2,dihe_list(1,ix,i))
                i2 = id_g2l(2,dihe_list(4,ix,i))

                if (icel1 < icel2) then

                  num_nb14 = num_nb14_calc(i1,icel)
                  duplicate = .false.

                  do k = 1, num_nonb_excl(i1,icel)
                    if (i2 == nonb_excl_list(k,i1,icel)) &
                      duplicate = .true.
                  end do

                  do k = 1, num_nb14
                    if (i2 == nb14_calc_list(k,i1,icel)) &
                      duplicate = .true.
                  end do

                  if (.not.duplicate) then
                    num_nb14 = num_nb14 + 1
                    num_nb14_calc(i1,icel) = num_nb14
                    nb14_calc_list(num_nb14,i1,icel) = i2
                    num_nb14_total(icel) = num_nb14_total(icel) + 1
                    sc_calc_list(num_nb14,i1,icel) = &
                   int(enefunc%dihe_periodicity(ix,i)/enefunc%notation_14types)
                  end if

                else if (icel1 > icel2) then

                  num_nb14 = num_nb14_calc(i2,icel)
                  duplicate = .false.

                  do k = 1, num_nonb_excl(i2,icel)
                    if (i1 == nonb_excl_list(k,i2,icel)) &
                      duplicate = .true.
                  end do

                  do k = 1, num_nb14
                    if (i1 == nb14_calc_list(k,i2,icel)) &
                      duplicate = .true.
                  end do

                  if (.not.duplicate) then
                    num_nb14 = num_nb14 + 1
                    num_nb14_calc(i2,icel) = num_nb14
                    nb14_calc_list(num_nb14,i2,icel) = i1
                    num_nb14_total(icel) = num_nb14_total(icel) + 1
                    sc_calc_list(num_nb14,i2,icel) = &
                    int(enefunc%dihe_periodicity(ix,i)/enefunc%notation_14types)
                  end if

                end if
              end if

            else

              if (mod(i-1,nthread) == id) then

                i1 = id_g2l(2,dihe_list(1,ix,i))
                i2 = id_g2l(2,dihe_list(4,ix,i))

                if (i1 < i2) then

                  num_nb14 = num_nb14_calc1(i1,i)
                  duplicate = .false.

                  do k = 1, num_nb14
                    if (i2 == nb14_calc_list1(k,i1,i)) &
                      duplicate = .true.
                  end do

                  do k = 1, num_nonb_excl1(i1,i)
                    if (i2 == nonb_excl_list1(k,i1,i)) &
                      duplicate = .true.
                  end do

                  if (.not.duplicate) then
                    num_nb14 = num_nb14 + 1
                    num_nb14_calc1(i1,i) = num_nb14
                    nb14_calc_list1(num_nb14,i1,i) = i2
                    num_nb14_total1(i) = num_nb14_total1(i) + 1
                    sc_calc_list1(num_nb14,i1,i) = &
                  int(enefunc%dihe_periodicity(ix,i)/enefunc%notation_14types)
                  end if

                else if (i1 > i2) then

                  num_nb14 = num_nb14_calc1(i2,i)
                  duplicate = .false.

                  do k = 1, num_nb14
                    if (i1 == nb14_calc_list1(k,i2,i)) &
                      duplicate = .true.
                  end do

                  do k = 1, num_nonb_excl1(i2,i)
                    if (i1 == nonb_excl_list1(k,i2,i)) &
                      duplicate = .true.
                  end do

                  if (.not.duplicate) then
                    num_nb14 = num_nb14 + 1
                    num_nb14_calc1(i2,i) = num_nb14
                    nb14_calc_list1(num_nb14,i2,i) = i1
                    num_nb14_total1(i) = num_nb14_total1(i) + 1
                    sc_calc_list1(num_nb14,i2,i) = &
                   int(enefunc%dihe_periodicity(ix,i)/enefunc%notation_14types)
                  end if

                end if

              end if

            end if

          end do
        end do

        !$omp end parallel

      end do

    end if

    ! Check the total number of exclusion list
    !
    if (first) then

      found1 = 0
      found2 = 0

      do icel = 1, maxcell_near
        found1 = found1 + num_excl_total(icel)
        found2 = found2 + num_nb14_total(icel)
      end do

      do icel = 1, ncell_local
        found1 = found1 + num_excl_total1(icel)
        found2 = found2 + num_nb14_total1(icel)
      end do

#ifdef HAVE_MPI_GENESIS
      call mpi_reduce(found1, enefunc%num_excl_all, 1, mpi_integer, mpi_sum, &
                      0, mpi_comm_country, ierror)
      call mpi_reduce(found2, enefunc%num_nb14_all, 1, mpi_integer, mpi_sum, &
                      0, mpi_comm_country, ierror)
#else
    enefunc%num_excl_all = found1
    enefunc%num_nb14_all = found2
#endif
    end if

    ! pack into small dimension
    !
    call pack_array_4to3(natom, domain%cell_pairlist1, &
                     num_nonb_excl, nonb_excl_list, enefunc%nonb_excl_list)

    call pack_array_4to3(natom, domain%cell_pairlist1, &
                     num_nb14_calc, nb14_calc_list, enefunc%nb14_calc_list)

    call pack_array_4to3(natom, domain%cell_pairlist1, &
                         num_nb14_calc,&
                     sc_calc_list, enefunc%sc_calc_list)

    call pack_array_3to2(natom, num_nonb_excl1, ncell_local, &
                     nonb_excl_list1, enefunc%nonb_excl_list1)

    call pack_array_3to2(natom, num_nb14_calc1, ncell_local, &
                     nb14_calc_list1, enefunc%nb14_calc_list1)

    call pack_array_3to2(natom, num_nb14_calc1, ncell_local,       &
                     sc_calc_list1, enefunc%sc_calc_list1)

   ! scnb/fudge_lj & scee/fudge_qq
   !
   if (enefunc%forcefield == ForcefieldAMBER .or. &
       enefunc%forcefield == ForcefieldGROAMBER .or. &
       enefunc%forcefield == ForcefieldGROMARTINI) then
      do i = 1, ncell_local
        num_nb14 = 0
        do ix = 1, natom(i) - 1
          ini_nb14  = num_nb14 + 1
          fin_nb14  = num_nb14 + enefunc%num_nb14_calc1(ix,i)
          num_nb14  = fin_nb14
          do k = ini_nb14, fin_nb14
            if (enefunc%forcefield == ForcefieldAMBER) then
              enefunc%nb14_lj_scale1(k, i) =  &
                        enefunc%dihe_scnb(enefunc%sc_calc_list1(k,i))
              enefunc%nb14_qq_scale1(k, i) =  &
                        enefunc%dihe_scee(enefunc%sc_calc_list1(k,i))
            else
              enefunc%nb14_lj_scale1(k, i) = enefunc%fudge_lj
              enefunc%nb14_qq_scale1(k, i) = enefunc%fudge_qq
            endif
          end do
        end do
      end do
      do ij = 1, maxcell_near
        i = cell_pair(1,ij)
        j = cell_pair(2,ij)
        num_nb14 = 0
        do ix = 1, natom(i)
          ini_nb14  = num_nb14 + 1
          fin_nb14  = num_nb14 + enefunc%num_nb14_calc(ix,ij)
          num_nb14  = fin_nb14
          do k = ini_nb14, fin_nb14
            if (enefunc%forcefield == ForcefieldAMBER) then
              enefunc%nb14_lj_scale(k, ij) =  &
                   enefunc%dihe_scnb(enefunc%sc_calc_list(k,ij))
              enefunc%nb14_qq_scale(k, ij) =  &
                   enefunc%dihe_scee(enefunc%sc_calc_list(k,ij))
            else
              enefunc%nb14_lj_scale(k, ij) = enefunc%fudge_lj
              enefunc%nb14_qq_scale(k, ij) = enefunc%fudge_qq
            endif
          end do
        end do
      end do
    end if

    return

  end subroutine count_nonb_excl_constraint

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    count_nonb_excl_constraint_rest
  !> @brief        exclude 1-2, 1-3 interactions and constraints w/o local
  !                restraint
  !! @authors      JJ
  !! @param[in]    first       : flag for first call or not
  !! @param[in]    water_table : flag for water table
  !! @param[in]    constraints : constraints information
  !! @param[inout] domain      : structure of domain
  !! @param[inout] enefunc     : structure of enefunc
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine count_nonb_excl_constraint_rest(first, water_table, constraints, &
                                             domain, enefunc)

    ! formal arguments
    logical,                     intent(in)    :: first
    logical,                     intent(in)    :: water_table
    type(s_constraints), target, intent(in)    :: constraints
    type(s_domain),      target, intent(inout) :: domain
    type(s_enefunc),     target, intent(inout) :: enefunc

    ! local variables
    integer                  :: ncell, ncell_local, i, ii, ij, ix, k, i1, i2, i3
    integer                  :: icel, icel1, icel2, jcel
    integer                  :: ic, j, ih
    integer                  :: ini_nb14, fin_nb14
    integer                  :: num_excl, num_nb14, id, omp_get_thread_num
    integer                  :: index(4)
    integer                  :: found1, found2
    integer                  :: fkind
    logical                  :: duplicate

    integer,         pointer :: nsolute(:), nwater(:), water_list(:,:,:)
    integer,         pointer :: id_g2l(:,:), natom(:)
    integer,         pointer :: nbond(:), bond_list(:,:,:)
    integer,         pointer :: nangle(:), angl_list(:,:,:)
    integer,         pointer :: ndihedral(:), dihe_list(:,:,:)
    integer,         pointer :: num_nonb_excl(:,:), num_nonb_excl1(:,:)
    integer,         pointer :: num_nb14_calc(:,:), num_nb14_calc1(:,:)
    integer,         pointer :: cell_pair(:,:), cell_pairlist2(:,:)
    integer,         pointer :: nonb_excl_list(:,:,:), nonb_excl_list1(:,:,:)
    integer,         pointer :: nb14_calc_list(:,:,:), nb14_calc_list1(:,:,:)
    integer,         pointer :: num_excl_total(:), num_nb14_total(:)
    integer,         pointer :: num_excl_total1(:), num_nb14_total1(:)
    integer,         pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)

    ! FEP
    integer                  :: iA, iB, iatomA, iatomB
    integer                  :: fg1, fg2

    HGr_local       => constraints%HGr_local
    HGr_bond_list   => constraints%HGr_bond_list

    nsolute         => domain%num_solute
    nwater          => domain%num_water
    natom           => domain%num_atom
    water_list      => domain%water_list
    id_g2l          => domain%id_g2l
    cell_pair       => domain%cell_pair
    cell_pairlist2  => domain%cell_pairlist2

    nbond           => enefunc%num_bond
    nangle          => enefunc%num_angle
    ndihedral       => enefunc%num_dihedral
    bond_list       => enefunc%bond_list
    angl_list       => enefunc%angle_list
    dihe_list       => enefunc%dihe_list
    num_nonb_excl   => enefunc%num_nonb_excl
    num_nonb_excl1  => enefunc%num_nonb_excl1
    num_nb14_calc   => enefunc%num_nb14_calc
    num_nb14_calc1  => enefunc%num_nb14_calc1
    nonb_excl_list  => enefunc%nonb_list
    nonb_excl_list1 => enefunc%nonb_list1
    nb14_calc_list  => enefunc%nb14_list
    nb14_calc_list1 => enefunc%nb14_list1
    num_excl_total  => enefunc%num_excl_total
    num_excl_total1 => enefunc%num_excl_total1
    num_nb14_total  => enefunc%num_nb14_total
    num_nb14_total1 => enefunc%num_nb14_total1

    ncell_local = domain%num_cell_local
    ncell       = domain%num_cell_local + domain%num_cell_boundary

    domain%max_num_atom = 0
    do i = 1, ncell
      domain%max_num_atom = max(domain%max_num_atom,domain%num_atom(i))
    end do

    ! initialization
    !
    num_nonb_excl  (1:domain%max_num_atom,1:maxcell_near)= 0
    num_nonb_excl1 (1:domain%max_num_atom,1:ncell_local) = 0
    num_nb14_calc  (1:domain%max_num_atom,1:maxcell_near)= 0
    num_nb14_calc1 (1:domain%max_num_atom,1:ncell_local) = 0
    num_excl_total (1:maxcell_near)                      = 0
    num_excl_total1(1:ncell_local)                       = 0
    num_nb14_total (1:maxcell_near)                      = 0
    num_nb14_total1(1:ncell_local)                       = 0

    ! exclude 1-2 interaction
    !
    !$omp parallel default(shared)                                     &
    !$omp private(id, i, ix, icel1, icel2, icel, i1, i2, i3, num_excl, &
    !$omp         duplicate, k, num_nb14, ic, j, ih, fkind,            &
    !$omp         iA, iB, iatomA, iatomB, fg1, fg2)
    !
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    if (enefunc%excl_level > 0) then

      do i = 1, ncell_local
        do ix = 1, nbond(i)

          icel1 = id_g2l(1,bond_list(1,ix,i))
          icel2 = id_g2l(1,bond_list(2,ix,i))
          fkind = enefunc%bond_kind(ix,i)

          if (fkind == 0) then
            if (icel1 /= icel2) then

              icel  = cell_pairlist2(icel1,icel2)
              if (mod(icel-1,nthread) == id) then

                i1 = id_g2l(2,bond_list(1,ix,i))
                i2 = id_g2l(2,bond_list(2,ix,i))

                if (i1 <= nsolute(icel1) .and. i2 <= nsolute(icel2)) then
                  if (icel1 < icel2) then
                    num_excl = num_nonb_excl(i1,icel) + 1
                    num_nonb_excl(i1,icel) = num_excl
                    nonb_excl_list(num_excl,i1,icel) = i2
                    num_excl_total(icel) = num_excl_total(icel) + 1

                  else if (icel1 > icel2) then
                    num_excl = num_nonb_excl(i2,icel) + 1
                    num_nonb_excl(i2,icel) = num_excl
                    nonb_excl_list(num_excl,i2,icel) = i1
                    num_excl_total(icel) = num_excl_total(icel) + 1

                  end if
                end if
              end if

            else

              if (mod(i-1,nthread) == id) then

                i1 = id_g2l(2,bond_list(1,ix,i))
                i2 = id_g2l(2,bond_list(2,ix,i))

                if (i1 <= nsolute(icel1) .and. i2 <= nsolute(icel2)) then
                  if (i1 < i2) then
                    num_excl = num_nonb_excl1(i1,i) + 1
                    num_nonb_excl1(i1,i) = num_excl
                    nonb_excl_list1(num_excl,i1,i) = i2
                    num_excl_total1(i) = num_excl_total1(i) + 1

                  else if (i1 > i2) then
                    num_excl = num_nonb_excl1(i2,i) + 1
                    num_nonb_excl1(i2,i) = num_excl
                    nonb_excl_list1(num_excl,i2,i) = i1
                    num_excl_total1(i) = num_excl_total1(i) + 1

                  end if
                end if

              end if

            end if

          end if

        end do
      end do

      ! exclude constraint
      !
      do icel = id+1, ncell_local, nthread
        do ic = 1, constraints%connect
          do j = 1, HGr_local(ic,icel)

            i1 = HGr_bond_list(1,j,ic,icel)
            do ih = 1, ic
              i2 = HGr_bond_list(ih+1,j,ic,icel)

              if (domain%fep_use) then
                ! FEP: Hydrogen rewiring
                ! To avoid SHAKE problem in FEP, singleB-dualB bonds is virtually
                ! rewired as singleA-dualB. In this case, singleA-dualB bonds are
                ! excluded, but singleB-dualB bonds are not excluded.
                ! To exclude singleB-dualB bonds, for singleA-dualB bond
                ! including hydrogen, the atom index of singleA is replaced
                ! with the corresponding atom index of singleB.
                fg1 = domain%fepgrp(i1,icel)
                fg2 = domain%fepgrp(i2,icel)
                if ((fg1 == 1) .and. (fg2 == 4)) then
                  do i = 1, domain%num_atom_singleB(icel)
                    iA = domain%id_singleA(i, icel, 3)
                    iB = domain%id_singleB(i, icel, 3)
                    iatomA = domain%id_singleA(iA, icel, 1)
                    iatomB = domain%id_singleB(iB, icel, 1)
                    if (iatomA == i1) then
                      i1 = iatomB
                      exit
                    end if
                  end do
                else if ((fg1 == 4) .and. (fg2 == 1)) then
                  do i = 1, domain%num_atom_singleB(icel)
                    iA = domain%id_singleA(i, icel, 3)
                    iB = domain%id_singleB(i, icel, 3)
                    iatomA = domain%id_singleA(iA, icel, 1)
                    iatomB = domain%id_singleB(iB, icel, 1)
                    if (iatomA == i2) then
                      i2 = iatomB
                      exit
                    end if
                  end do
                end if
              end if

              num_excl = num_nonb_excl1(i1,icel) + 1
              num_nonb_excl1(i1,icel) = num_excl
              nonb_excl_list1(num_excl,i1,icel) = i2
              num_excl_total1(icel) = num_excl_total1(icel) + 1

            end do

          end do
        end do
      end do

    end if

    ! exclude water
    !
    if (enefunc%excl_level > 1) then

      if (constraints%tip4) then

        do icel = id+1, ncell_local, nthread
          do ic = 1, nwater(icel)

            index(1:4) = water_list(1:4,ic,icel)

            do i1 = 1, 3

              num_excl = 0

              do i2 = i1+1, 4
                num_excl = num_excl + 1      
                nonb_excl_list1(num_excl,index(i1),icel) = index(i2)
              end do

              num_nonb_excl1(index(i1),icel) = num_excl

            end do

          end do
        end do

      else 

        do icel = id+1, ncell_local, nthread
          do ic = 1, nwater(icel)

            i1 = water_list(1,ic,icel)
            i2 = water_list(2,ic,icel)
            i3 = water_list(3,ic,icel)

            num_nonb_excl1(i1,icel) = 2
            num_nonb_excl1(i2,icel) = 1
            num_nonb_excl1(i3,icel) = 0
            nonb_excl_list1(1,i1,icel) = i2
            nonb_excl_list1(2,i1,icel) = i3
            nonb_excl_list1(1,i2,icel) = i3
            num_excl_total1(icel) = num_excl_total1(icel) + 3

          end do
        end do

      end if

      ! exclude 1-3 interaction
      !
      do i = 1, ncell_local
        do ix = 1, nangle(i)

          icel1 = id_g2l(1,angl_list(1,ix,i))
          icel2 = id_g2l(1,angl_list(3,ix,i))
          fkind = enefunc%angle_kind(ix,i)

          if (fkind == 0) then

            if (icel1 /= icel2) then

              icel  = cell_pairlist2(icel1,icel2)
              if (mod(icel-1,nthread) == id) then

                i1 = id_g2l(2,angl_list(1,ix,i))
                i2 = id_g2l(2,angl_list(3,ix,i))

                if (i1 <= nsolute(icel1) .and. i2 <= nsolute(icel2)) then
                  if (icel1 < icel2) then

                    num_excl = num_nonb_excl(i1,icel)
                    duplicate = .false.

                    do k = 1, num_excl
                      if (i2 == nonb_excl_list(k,i1,icel)) &
                        duplicate = .true.
                    end do

                    if (.not. duplicate) then
                      num_excl = num_excl + 1
                      num_nonb_excl(i1,icel) = num_excl
                      nonb_excl_list(num_excl,i1,icel) = i2
                      num_excl_total(icel) =  num_excl_total(icel) + 1
                    end if

                  else if (icel1 > icel2) then

                    num_excl = num_nonb_excl(i2,icel)
                    duplicate = .false.

                    do k = 1, num_excl
                      if (i1 == nonb_excl_list(k,i2,icel)) &
                        duplicate = .true.
                    end do

                    if (.not.duplicate) then
                      num_excl = num_excl + 1
                      num_nonb_excl(i2,icel) = num_excl
                      nonb_excl_list(num_excl,i2,icel) = i1
                      num_excl_total(icel) = num_excl_total(icel) + 1
                    end if

                  end if
                end if
              end if

            else

              if (mod(i-1,nthread) == id) then

                i1 = id_g2l(2,angl_list(1,ix,i))
                i2 = id_g2l(2,angl_list(3,ix,i))

                if (i1 <= nsolute(icel1) .and. i2 <= nsolute(icel2)) then
                  if (i1 < i2) then

                    num_excl = num_nonb_excl1(i1,i)
                    duplicate = .false.

                    do k = 1, num_excl
                      if (i2 == nonb_excl_list1(k,i1,i)) &
                        duplicate = .true.
                    end do

                    if (.not.duplicate) then
                      num_excl = num_excl + 1
                      num_nonb_excl1(i1,i) = num_excl
                      nonb_excl_list1(num_excl,i1,i) = i2
                      num_excl_total1(i) = num_excl_total1(i) + 1
                    end if

                  else if (i1 > i2) then

                    num_excl = num_nonb_excl1(i2,i)
                    duplicate = .false.

                    do k = 1, num_excl
                      if (i1 == nonb_excl_list1(k,i2,i)) &
                        duplicate = .true.
                    end do

                    if (.not.duplicate) then
                      num_excl = num_excl + 1
                      num_nonb_excl1(i2,i) = num_excl
                      nonb_excl_list1(num_excl,i2,i) = i1
                      num_excl_total1(i) = num_excl_total1(i) + 1
                    end if

                  end if
                end if

              end if
            end if

          end if

        end do
      end do

    end if

    !$omp end parallel

    ! count 1-4 interaction
    !
    if (enefunc%excl_level > 2) then

      do ii = 1, 2

        if (ii == 2) then
          ndihedral => enefunc%num_rb_dihedral
          dihe_list => enefunc%rb_dihe_list
        end if

        !$omp parallel default(shared)                           &
        !$omp private(id, i, ix, icel1, icel2, icel, i1, i2, i3, &
        !$omp          duplicate, k, num_nb14, ic, j, ih, fkind)
        !

#ifdef OMP
        id = omp_get_thread_num()
#else
        id = 0
#endif

        do i = 1, ncell_local
          do ix = 1, ndihedral(i)

            icel1 = id_g2l(1,dihe_list(1,ix,i))
            icel2 = id_g2l(1,dihe_list(4,ix,i))
            fkind = enefunc%dihe_kind(ix,i)

            if (fkind == 0) then

              if (icel1 /= icel2) then

                icel  = cell_pairlist2(icel1,icel2)
                if (mod(icel-1,nthread) == id) then

                  i1 = id_g2l(2,dihe_list(1,ix,i))
                  i2 = id_g2l(2,dihe_list(4,ix,i))

                  if (icel1 < icel2) then

                    num_nb14 = num_nb14_calc(i1,icel)
                    duplicate = .false.

                    do k = 1, num_nonb_excl(i1,icel)
                      if (i2 == nonb_excl_list(k,i1,icel)) &
                        duplicate = .true.
                    end do

                    do k = 1, num_nb14
                      if (i2 == nb14_calc_list(k,i1,icel)) &
                        duplicate = .true.
                    end do

                    if (.not.duplicate) then
                      num_nb14 = num_nb14 + 1
                      num_nb14_calc(i1,icel) = num_nb14
                      nb14_calc_list(num_nb14,i1,icel) = i2
                      num_nb14_total(icel) = num_nb14_total(icel) + 1
                    end if

                  else if (icel1 > icel2) then

                    num_nb14 = num_nb14_calc(i2,icel)
                    duplicate = .false.

                    do k = 1, num_nonb_excl(i2,icel)
                      if (i1 == nonb_excl_list(k,i2,icel)) &
                        duplicate = .true.
                    end do

                    do k = 1, num_nb14
                      if (i1 == nb14_calc_list(k,i2,icel)) &
                        duplicate = .true.
                    end do

                    if (.not.duplicate) then
                      num_nb14 = num_nb14 + 1
                      num_nb14_calc(i2,icel) = num_nb14
                      nb14_calc_list(num_nb14,i2,icel) = i1
                      num_nb14_total(icel) = num_nb14_total(icel) + 1
                    end if

                  end if
                end if

              else

                if (mod(i-1,nthread) == id) then

                  i1 = id_g2l(2,dihe_list(1,ix,i))
                  i2 = id_g2l(2,dihe_list(4,ix,i))

                  if (i1 < i2) then

                    num_nb14 = num_nb14_calc1(i1,i)
                    duplicate = .false.

                    do k = 1, num_nb14
                      if (i2 == nb14_calc_list1(k,i1,i)) &
                        duplicate = .true.
                    end do

                    do k = 1, num_nonb_excl1(i1,i)
                      if (i2 == nonb_excl_list1(k,i1,i)) &
                        duplicate = .true.
                    end do

                    if (.not.duplicate) then
                      num_nb14 = num_nb14 + 1
                      num_nb14_calc1(i1,i) = num_nb14
                      nb14_calc_list1(num_nb14,i1,i) = i2
                      num_nb14_total1(i) = num_nb14_total1(i) + 1
                    end if

                  else if (i1 > i2) then

                    num_nb14 = num_nb14_calc1(i2,i)
                    duplicate = .false.

                    do k = 1, num_nb14
                      if (i1 == nb14_calc_list1(k,i2,i)) &
                        duplicate = .true.
                    end do

                    do k = 1, num_nonb_excl1(i2,i)
                      if (i1 == nonb_excl_list1(k,i2,i)) &
                        duplicate = .true.
                    end do

                    if (.not.duplicate) then
                      num_nb14 = num_nb14 + 1
                      num_nb14_calc1(i2,i) = num_nb14
                      nb14_calc_list1(num_nb14,i2,i) = i1
                      num_nb14_total1(i) = num_nb14_total1(i) + 1
                    end if

                  end if

                end if

              end if
            end if

          end do
        end do

        !$omp end parallel

      end do

    end if

    ! Check the total number of exclusion list
    !
    if (first) then

      found1 = 0
      found2 = 0

      do icel = 1, maxcell_near
        found1 = found1 + num_excl_total(icel)
        found2 = found2 + num_nb14_total(icel)
      end do

      do icel = 1, ncell_local
        found1 = found1 + num_excl_total1(icel)
        found2 = found2 + num_nb14_total1(icel)
      end do

#ifdef HAVE_MPI_GENESIS
      call mpi_reduce(found1, enefunc%num_excl_all, 1, mpi_integer, mpi_sum, &
                      0, mpi_comm_country, ierror)
      call mpi_reduce(found2, enefunc%num_nb14_all, 1, mpi_integer, mpi_sum, &
                      0, mpi_comm_country, ierror)
#else
    enefunc%num_excl_all = found1
    enefunc%num_nb14_all = found2
#endif
    end if

    ! pack into small dimension
    !
    call pack_array_4to3(natom, domain%cell_pairlist1, &
                     num_nonb_excl, nonb_excl_list, enefunc%nonb_excl_list)

    call pack_array_4to3(natom, domain%cell_pairlist1, &
                     num_nb14_calc, nb14_calc_list, enefunc%nb14_calc_list)

    call pack_array_3to2(natom, num_nonb_excl1, ncell_local, &
                     nonb_excl_list1, enefunc%nonb_excl_list1)

    call pack_array_3to2(natom, num_nb14_calc1, ncell_local, &
                     nb14_calc_list1, enefunc%nb14_calc_list1)

   ! scnb/fudge_lj & scee/fudge_qq
   !
   if (enefunc%forcefield == ForcefieldAMBER .or. &
       enefunc%forcefield == ForcefieldGROAMBER .or. &
       enefunc%forcefield == ForcefieldGROMARTINI) then
      do i = 1, ncell_local
        num_nb14 = 0
        do ix = 1, natom(i) - 1
          ini_nb14  = num_nb14 + 1
          fin_nb14  = num_nb14 + enefunc%num_nb14_calc1(ix,i)
          num_nb14  = fin_nb14
          do k = ini_nb14, fin_nb14
            if (enefunc%forcefield == ForcefieldAMBER) then
              enefunc%nb14_lj_scale1(k, i) =  &
                        enefunc%dihe_scnb(enefunc%sc_calc_list1(k,i))
              enefunc%nb14_qq_scale1(k, i) =  &
                        enefunc%dihe_scee(enefunc%sc_calc_list1(k,i))
            else
              enefunc%nb14_lj_scale1(k, i) = enefunc%fudge_lj
              enefunc%nb14_qq_scale1(k, i) = enefunc%fudge_qq
            endif
          end do
        end do
      end do
      do ij = 1, maxcell_near
        i = cell_pair(1,ij)
        j = cell_pair(2,ij)
        num_nb14 = 0
        do ix = 1, natom(i)
          ini_nb14  = num_nb14 + 1
          fin_nb14  = num_nb14 + enefunc%num_nb14_calc(ix,ij)
          num_nb14  = fin_nb14
          do k = ini_nb14, fin_nb14
            if (enefunc%forcefield == ForcefieldAMBER) then
              enefunc%nb14_lj_scale(k, ij) =  &
                   enefunc%dihe_scnb(enefunc%sc_calc_list(k,ij))
              enefunc%nb14_qq_scale(k, ij) =  &
                   enefunc%dihe_scee(enefunc%sc_calc_list(k,ij))
            else
              enefunc%nb14_lj_scale(k, ij) = enefunc%fudge_lj
              enefunc%nb14_qq_scale(k, ij) = enefunc%fudge_qq
            endif
          end do
        end do
      end do
    end if

    return

  end subroutine count_nonb_excl_constraint_rest

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pack_array_4to3
  !> @brief        pack 4D-array into 3D-array
  !! @authors      YS, JJ
  !! @param[in]    num_ele  :
  !! @param[in]    num_list :
  !! @param[in]    list_1d  :
  !! @param[out]   list_3d  :
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pack_array_4to3(num_ele, cell_list, num_list, list_3d, list_2d)

    ! formal arguments
    integer,                 intent(in)    :: num_ele(:)
    integer,                 intent(in)    :: cell_list(:,:)
    integer,                 intent(in)    :: num_list(:,:)
    integer,                 intent(in)    :: list_3d(:,:,:)
    integer,                 intent(inout) :: list_2d(:,:)

    ! local variables
    integer                  :: i, j, k, icel, jcel
    integer                  :: id, omp_get_thread_num


    !$omp parallel default(shared) private(id, jcel, icel, k, i, j)
    !

#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do jcel = id+1, maxcell_near, nthread
      icel = cell_list(1,jcel)
      k = 0
      do i = 1, num_ele(icel)
        do j = 1, num_list(i,jcel)
          k = k + 1
          list_2d(k,jcel) = list_3d(j,i,jcel)
        end do
      end do
    end do

    !$omp end parallel

    return

  end subroutine pack_array_4to3

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pack_array_3to2
  !> @brief        pack 3D-array into 2D-array
  !! @authors      YS
  !! @param[in]    num_ele  :
  !! @param[in]    num_list :
  !! @param[in]    list_1d  :
  !! @param[out]   list_2d  :
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pack_array_3to2(num_ele, num_list, ncell_local, list_2d, list_1d)

    ! formal arguments
    integer,                 intent(in)    :: num_ele(:)
    integer,                 intent(in)    :: num_list(:,:)
    integer,                 intent(in)    :: ncell_local
    integer,                 intent(in)    :: list_2d(:,:,:)
    integer,                 intent(inout) :: list_1d(:,:)

    ! local variables
    integer                  :: i, j, k, icel
    integer                  :: id, omp_get_thread_num


    !$omp parallel default(shared) private(id, icel, k, i, j)
    !

#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif

    do icel = id+1, ncell_local, nthread
      k = 0
      do i = 1, num_ele(icel)
        do j = 1, num_list(i,icel)
          k = k + 1
          list_1d(k,icel) = list_2d(j,i,icel)
        end do
      end do
    end do

    !$omp end parallel

    return

  end subroutine pack_array_3to2

end module sp_enefunc_charmm_mod
