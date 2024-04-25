!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_enefunc_gromacs_mod
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

module sp_enefunc_gromacs_mod

  use sp_enefunc_charmm_mod
  use sp_enefunc_restraints_mod
  use sp_enefunc_table_mod
  use sp_energy_mod
  use sp_restraints_str_mod
  use sp_constraints_str_mod
  use sp_enefunc_str_mod
  use sp_energy_str_mod
  use sp_domain_str_mod
  use molecules_str_mod
  use fileio_grotop_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: define_enefunc_gromacs
  public  :: count_nonb_excl_go
  private :: setup_enefunc_bond
  private :: setup_enefunc_bond_constraint
  private :: setup_enefunc_angl
  private :: setup_enefunc_angl_constraint
  private :: setup_enefunc_dihe
  private :: setup_enefunc_rb_dihe
  private :: setup_enefunc_impr
  private :: setup_enefunc_nonb
  private :: setup_enefunc_gro_restraints

  ! paramters
  integer, parameter :: WaterIdx(2,3) = reshape((/1,2,1,3,2,3/),shape=(/2,3/))

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_enefunc_gromacs
  !> @brief        a driver subroutine for defining potential energy functions
  !! @authors      NT
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[in]    grotop      : GROMACS TOP information
  !! @param[in]    molecule    : molecule information
  !! @param[in]    constraints : constraints information
  !! @param[in]    restraints  : restraints information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_enefunc_gromacs(ene_info, grotop, molecule, &
                                    constraints, restraints, domain, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info 
    type(s_grotop),          intent(in)    :: grotop
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

    call alloc_enefunc(enefunc, EneFuncBase,     ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncBond,     ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncAngl,     ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncDihe,     ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncRBDihe,   ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncImpr,     ncel, ncel)
    call alloc_enefunc(enefunc, EneFuncBondCell, ncel, ncelb)

    if (.not. constraints%rigid_bond) then

      ! bond
      !
      call setup_enefunc_bond(grotop, molecule, domain, constraints, enefunc)

      ! angle
      !
      call setup_enefunc_angl(grotop, molecule, domain, enefunc)

    else

      ! bond
      !
      call setup_enefunc_bond_constraint(grotop, molecule, domain, constraints,&
                                         enefunc)

      ! angle
      !
      call setup_enefunc_angl_constraint(grotop, molecule, domain, constraints,&
                                         enefunc)

    end if

    ! dihedral
    !
    call setup_enefunc_dihe(grotop, molecule, domain, enefunc)

    ! Ryckaert-Bellemans dihedral
    !
    call setup_enefunc_rb_dihe(grotop, molecule, domain, enefunc)

    ! improper
    !
    call setup_enefunc_impr(grotop, molecule, domain, enefunc)

    ! nonbonded
    !
    call setup_enefunc_nonb(grotop, molecule, constraints, domain, enefunc)

    ! lookup table
    !
    call setup_enefunc_table(ene_info, enefunc)

    ! restraints
    !
    call setup_enefunc_restraints(molecule, restraints, domain, enefunc)

    ! restraints (gromacs)
    !
    call setup_enefunc_gro_restraints(molecule, grotop, domain, enefunc)


    ! write summary of energy function
    !
    if (main_rank) then
      write(MsgOut,'(A)') &
           'Define_Enefunc_Gromacs> Number of Interactions in Each Term'
      write(MsgOut,'(A20,I10,A20,I10)')                         &
           '  bond_ene        = ', enefunc%num_bond_all,        &
           '  angle_ene       = ', enefunc%num_angl_all
      write(MsgOut,'(A20,I10,A20,I10)')                         &
           '  torsion_ene     = ', enefunc%num_dihe_all,        &
           '  rb_torsion_ene  = ', enefunc%num_rb_dihe_all
      write(MsgOut,'(A20,I10)')                                 &
           '  improper_ene    = ', enefunc%num_impr_all
      if (.not. ene_info%table) then
        write(MsgOut,'(A20,I10,A20,I10)')                       &
           '  nb_exclusions   = ', enefunc%num_excl_all,        &
           '  nb14_calc       = ', enefunc%num_nb14_all
      end if
      write(MsgOut,'(A20,I10,A20,I10)')                         &
           ' restraint_groups = ', enefunc%num_restraintgroups, &
           ' restraint_funcs  = ', enefunc%num_restraintfuncs
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine define_enefunc_gromacs

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_bond
  !> @brief        define BOND term for each cell in potential energy function
  !! @authors      NT
  !! @param[in]    grotop   : GROMACS TOP information
  !! @param[in]    molecule : molecule information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_bond(grotop, molecule, domain, constraints, enefunc)

    ! formal arguments
    type(s_grotop),          intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
    type(s_domain),  target, intent(in)    :: domain
    type(s_constraints),     intent(in)    :: constraints
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    real(wp)                 :: water_dist(3)
    integer                  :: i, j, k
    integer                  :: ioffset, nbond_a
    integer                  :: idx1, idx2, icel1, icel2, icel_local
    integer                  :: ii1, ii2

    type(s_grotop_mol), pointer :: gromol
    real(wp),           pointer :: force(:,:), dist(:,:)
    integer,            pointer :: bond(:), list(:,:,:)
    integer,            pointer :: ncel, cell_pair(:,:)
    integer,            pointer :: id_g2l(:,:)


    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l

    bond      => enefunc%num_bond
    list      => enefunc%bond_list
    force     => enefunc%bond_force_const
    dist      => enefunc%bond_dist_min

    ioffset   = 0
    nbond_a   = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count

        if (gromol%settles%func == 0) then

          do k = 1, gromol%num_bonds

            idx1 = gromol%bonds(k)%atom_idx1 + ioffset
            idx2 = gromol%bonds(k)%atom_idx2 + ioffset

            if (domain%fep_use) then
              ! FEP: If the bond are not set to any group of FEP, exclude this bond.
              ii1 = molecule%fepgrp(idx1)
              ii2 = molecule%fepgrp(idx2)
              if (molecule%fepgrp_bond(ii1,ii2) == 0) cycle
            end if

            icel1 = id_g2l(1,idx1)
            icel2 = id_g2l(1,idx2)

            if (icel1 /= 0 .and. icel2 /= 0) then

              icel_local = cell_pair(icel1,icel2)

              if (icel_local > 0 .and. icel_local <= ncel) then

                if (bond(icel_local)+1 > MaxBond) &
                  call error_msg('Setup_Enefunc_Bond> Too many bonds.') 

                nbond_a = nbond_a + 1

                bond(icel_local) = bond(icel_local) + 1
                list (1,bond(icel_local),icel_local) = idx1
                list (2,bond(icel_local),icel_local) = idx2
                force(  bond(icel_local),icel_local) = &
                                 gromol%bonds(k)%kb * 0.01_wp * JOU2CAL * 0.5_wp
                dist (  bond(icel_local),icel_local) = &
                                 gromol%bonds(k)%b0 * 10.0_wp
              end if

            end if

          end do

        else

          water_dist = (/gromol%settles%doh * 10.0_wp, &
                         gromol%settles%doh * 10.0_wp, &
                         gromol%settles%dhh * 10.0_wp/)

          do k = 1, 3

            idx1 = WaterIdx(1,k) + ioffset
            idx2 = WaterIdx(2,k) + ioffset

            icel1 = id_g2l(1,idx1)
            icel2 = id_g2l(1,idx2)

            if (icel1 /= 0 .and. icel2 /= 0) then

              icel_local = cell_pair(icel1,icel2)

              if (icel_local > 0 .and. icel_local <= ncel) then

                if (bond(icel_local)+1 > MaxBond) &
                  call error_msg('Setup_Enefunc_Bond> Too many bonds.') 

                nbond_a = nbond_a + 1

                bond(icel_local) = bond(icel_local) + 1
                list (1,bond(icel_local),icel_local) = idx1
                list (2,bond(icel_local),icel_local) = idx2
                force(  bond(icel_local),icel_local) = 0.0_wp
                dist (  bond(icel_local),icel_local) = water_dist(k)
              end if

            end if

          end do

        end if

        ioffset = ioffset + gromol%num_atoms

      end do
    end do

    if (enefunc%table%tip4) then

      ioffset = 0

      do i = 1, grotop%num_molss
        gromol => grotop%molss(i)%moltype%mol
        do j = 1, grotop%molss(i)%count

          if (gromol%settles%func == 1) then

            idx1 = 2 + ioffset
            idx2 = 4 + ioffset
            
            icel1 = id_g2l(1,idx1)
            icel2 = id_g2l(1,idx2)

            if (icel1 /= 0 .and. icel2 /= 0) then

              icel_local = cell_pair(icel1,icel2)

              if (icel_local > 0 .and. icel_local <= ncel) then

                if (bond(icel_local)+1 > MaxBond) &
                  call error_msg('Setup_Enefunc_Bond> Too many bonds.')

                nbond_a = nbond_a + 1

                bond(icel_local) = bond(icel_local) + 1
                list (1,bond(icel_local),icel_local) = idx1
                list (2,bond(icel_local),icel_local) = idx2
                force(  bond(icel_local),icel_local) = 0.0_wp
                dist (  bond(icel_local),icel_local) = constraints%water_rOD
              end if
            end if

            idx1 = 3 + ioffset
            idx2 = 4 + ioffset

            icel1 = id_g2l(1,idx1)
            icel2 = id_g2l(1,idx2)

            if (icel1 /= 0 .and. icel2 /= 0) then

              icel_local = cell_pair(icel1,icel2)

              if (icel_local > 0 .and. icel_local <= ncel) then

                if (bond(icel_local)+1 > MaxBond) &
                  call error_msg('Setup_Enefunc_Bond> Too many bonds.')

                nbond_a = nbond_a + 1

                bond(icel_local) = bond(icel_local) + 1
                list (1,bond(icel_local),icel_local) = idx1
                list (2,bond(icel_local),icel_local) = idx2
                force(  bond(icel_local),icel_local) = 0.0_wp
                dist (  bond(icel_local),icel_local) = constraints%water_rOD
              end if
            end if

          end if

          ioffset = ioffset + gromol%num_atoms

        end do
      end do

    end if
      
#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(nbond_a, enefunc%num_bond_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_bond_all = nbond_a
#endif

    return
    
  end subroutine setup_enefunc_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_bond_constraint
  !> @brief        define BOND term between heavy atoms
  !! @authors      NT
  !! @param[in]    grotop      : CHARMM grotop information
  !! @param[in]    molecule    : molecule information
  !! @param[in]    domain      : domain information
  !! @param[in]    constraints : constraints information
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_bond_constraint(grotop, molecule, domain, &
                                           constraints, enefunc)

    ! formal arguments
    type(s_grotop),              intent(in)    :: grotop
    type(s_molecule),            intent(in)    :: molecule
    type(s_domain),      target, intent(in)    :: domain
    type(s_constraints), target, intent(inout) :: constraints
    type(s_enefunc),     target, intent(inout) :: enefunc

    ! local variables
    real(wp)                 :: water_dist(3)
    integer                  :: i, j, k, l, m, n
    integer                  :: nwat, ioffset, nbond_a, nbond_c
    integer                  :: icel, connect, ih, ih1, ih2
    integer                  :: i1, i2, idx1, idx2, icel1, icel2, icel_local
    integer                  :: wat_bonds, wat_found, nbond_sys
    integer                  :: ii1, ii2
    character(6)             :: res1, res2
    character(4)             :: atm1, atm2
    logical                  :: cl1, cl2

    type(s_grotop_mol), pointer :: gromol
    real(wp),           pointer :: force(:,:), dist(:,:)
    real(dp),           pointer :: HGr_bond_dist(:,:,:,:)
    integer,            pointer :: bond(:), list(:,:,:)
    integer,            pointer :: ncel, cell_pair(:,:)
    integer,            pointer :: id_g2l(:,:), id_l2g(:,:)
    integer,            pointer :: HGr_local(:,:), HGr_bond_list(:,:,:,:)


    ncel          => domain%num_cell_local
    cell_pair     => domain%cell_pair
    id_g2l        => domain%id_g2l
    id_l2g        => domain%id_l2g

    HGr_local     => constraints%HGr_local
    HGr_bond_list => constraints%HGr_bond_list
    HGr_bond_dist => constraints%HGr_bond_dist
    connect       =  constraints%connect

    bond          => enefunc%num_bond
    list          => enefunc%bond_list
    force         => enefunc%bond_force_const
    dist          => enefunc%bond_dist_min
    nwat          =  enefunc%table%num_water

    ioffset   = 0
    nbond_sys = 0

    nbond_a   = 0
    nbond_c   = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count

        if (gromol%settles%func == 0) then

          do k = 1, gromol%num_bonds

            i1 = gromol%bonds(k)%atom_idx1
            i2 = gromol%bonds(k)%atom_idx2

            if (domain%fep_use) then
              ! FEP: If the bond are not set to any group of FEP,
              ! exclude this bond.
              ii1 = molecule%fepgrp(i1)
              ii2 = molecule%fepgrp(i2)
              if (molecule%fepgrp_bond(ii1,ii2) == 0) cycle
            end if

            atm1 = gromol%atoms(i1)%atom_type
            atm2 = gromol%atoms(i2)%atom_type

            cl1 = (atm1(1:1) /= 'H' .and. atm1(1:1) /= 'h')
            cl2 = (atm2(1:1) /= 'H' .and. atm2(1:1) /= 'h')
            if (constraints%hydrogen_type == ConstraintAtomMass) then
              cl1 = (gromol%atoms(i1)%mass > LIGHT_ATOM_MASS_LIMIT) 
              cl2 = (gromol%atoms(i2)%mass > LIGHT_ATOM_MASS_LIMIT) 
            else if (constraints%hydrogen_type == ConstraintAtomBoth) then
              cl1 = (cl1 .and. &
                  gromol%atoms(i1)%mass > LIGHT_ATOM_MASS_LIMIT) 
              cl2 = (cl2 .and. &
                 gromol%atoms(i2)%mass > LIGHT_ATOM_MASS_LIMIT) 
            endif
            idx1 = i1 + ioffset
            idx2 = i2 + ioffset

            if (cl1 .and. cl2) then

              icel1 = id_g2l(1,idx1)
              icel2 = id_g2l(1,idx2)

              if (icel1 /= 0 .and. icel2 /= 0) then

                icel_local = cell_pair(icel1,icel2)

                if (icel_local > 0 .and. icel_local <= ncel) then

                  if (bond(icel_local)+1 > MaxBond) &
               call error_msg('Setup_Enefunc_Bond_Constraint> Too many bonds.') 

                  nbond_a = nbond_a + 1

                  bond(icel_local) = bond(icel_local) + 1
                  list (1,bond(icel_local),icel_local) = idx1
                  list (2,bond(icel_local),icel_local) = idx2
                  force(  bond(icel_local),icel_local) = &
                                 gromol%bonds(k)%kb * 0.01_wp * JOU2CAL * 0.5_wp
                  dist (  bond(icel_local),icel_local) = &
                                 gromol%bonds(k)%b0 * 10.0_wp
                end if

              end if

            else

              icel1 = id_g2l(1,idx1)
              icel2 = id_g2l(1,idx2)

              if (icel1 /= 0 .and. icel2 /= 0) then

                icel = cell_pair(icel1,icel2)

                if (icel > 0 .and. icel <= ncel) then
                
                  do m = 1, connect

                    do n = 1, HGr_local(m,icel)
                      ih1 = id_l2g(HGr_bond_list(1,n,m,icel),icel)
                      do ih = 1, m
                        ih2 = id_l2g(HGr_bond_list(ih+1,n,m,icel),icel)

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

                        if (ih1 == idx1 .and. ih2 == idx2 .or. &
                          ih2 == idx1 .and. ih1 == idx2) then

                          nbond_c = nbond_c + 1
                          HGr_bond_dist(ih+1,n,m,icel) = &
                                     gromol%bonds(k)%b0 * 10.0_wp
  
                          goto 1

                        end if

                      end do
                    end do

                  end do
                end if
              end if
1           continue

            end if

          end do

          nbond_sys = nbond_sys + gromol%num_bonds

        else

          water_dist = (/gromol%settles%doh * 10.0_wp, &
                         gromol%settles%doh * 10.0_wp, &
                         gromol%settles%dhh * 10.0_wp/)

          do k = 1, 3

            idx1 = WaterIdx(1,k) + ioffset
            idx2 = WaterIdx(2,k) + ioffset

            icel1 = id_g2l(1,idx1)
            icel2 = id_g2l(1,idx2)

            if (icel1 /= 0 .and. icel2 /= 0) then

              icel = cell_pair(icel1,icel2)

              if (icel > 0 .and. icel <= ncel) then

                do m = 1, connect

                  do n = 1, HGr_local(m,icel)
                    ih1 = id_l2g(HGr_bond_list(1,n,m,icel),icel)
                    do ih = 1, m
                      ih2 = id_l2g(HGr_bond_list(ih+1,n,m,icel),icel)

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

                      if (ih1 == idx1 .and. ih2 == idx2 .or. &
                          ih2 == idx1 .and. ih1 == idx2) then

                        nbond_c = nbond_c + 1
                        HGr_bond_dist(ih+1,n,m,icel) = water_dist(k)
  
                        goto 2

                      end if

                    end do
                  end do

                end do
              end if
            end if
2           continue

          end do

          nbond_sys = nbond_sys + 3

        end if

        ioffset   = ioffset   + gromol%num_atoms

      end do
    end do

    ! for water molecule
    !
    wat_bonds = 0
    wat_found = 0

    if (constraints%fast_water) then

      do i = 1, grotop%num_molss
        gromol => grotop%molss(i)%moltype%mol
        do j = 1, grotop%molss(i)%count

          if (gromol%settles%func == 0) then

            do k = 1, gromol%num_bonds

              i1 = gromol%bonds(k)%atom_idx1
              i2 = gromol%bonds(k)%atom_idx2

              res1 = gromol%atoms(i1)%residue_name
              res2 = gromol%atoms(i2)%residue_name

              if (res1 == constraints%water_model .and. &
                  res2 == constraints%water_model) then

                wat_bonds = wat_bonds+1

                if (k == 1) &
                wat_found = wat_found + 1

              end if

            end do

          else

            wat_bonds = wat_bonds + 3
            wat_found = wat_found + 1

          end if

        end do
      end do

      if (wat_found /= nwat) &
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

    if (constraints%fast_water) then

      if (enefunc%num_bond_all /= &
          (nbond_sys-constraints%num_bonds-wat_bonds)) then
        call error_msg( &
          'Setup_Enefunc_Bond_Constraint> Some bond paremeters are missing.')
      end if

    else

      if (enefunc%num_bond_all /= &
          (nbond_sys-constraints%num_bonds-3*nwat))then
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
  !! @authors      NT
  !! @param[in]    grotop   : GROMACS topology informaiton
  !! @param[in]    molecule : molecule information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_angl(grotop, molecule, domain, enefunc)

    ! formal arguments
    type(s_grotop),  target, intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k
    integer                  :: ioffset, nangl_a
    integer                  :: idx1, idx2, idx3, icel1, icel2, icel_local
    integer                  :: i1, i2, i3

    type(s_grotop_mol), pointer :: gromol
    real(wp),           pointer :: force(:,:), theta(:,:)
    integer,            pointer :: angle(:), list(:,:,:)
    integer,            pointer :: ncel, cell_pair(:,:)
    integer,            pointer :: id_g2l(:,:)


    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l

    angle     => enefunc%num_angle
    list      => enefunc%angle_list
    force     => enefunc%angle_force_const
    theta     => enefunc%angle_theta_min

    ioffset   = 0
    nangl_a   = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count

        if (gromol%settles%func == 0) then

          do k = 1, gromol%num_angls

            idx1 = gromol%angls(k)%atom_idx1 + ioffset
            idx2 = gromol%angls(k)%atom_idx2 + ioffset
            idx3 = gromol%angls(k)%atom_idx3 + ioffset

            if (domain%fep_use) then
              ! FEP: If the angle are not set to any group of FEP, exclude this angle.
              i1 = molecule%fepgrp(idx1)
              i2 = molecule%fepgrp(idx2)
              i3 = molecule%fepgrp(idx3)
              if (molecule%fepgrp_angl(i1,i2,i3) == 0) cycle
            end if

            icel1 = id_g2l(1,idx1)
            icel2 = id_g2l(1,idx3)

            if (icel1 /= 0 .and. icel2 /= 0) then

              icel_local = cell_pair(icel1,icel2)

              if (icel_local > 0 .and. icel_local <= ncel) then

                if (angle(icel_local)+1 > MaxAngle) &
                  call error_msg('Setup_Enefunc_Angl> Too many angles.') 

                nangl_a = nangl_a + 1

                angle(icel_local) = angle(icel_local) + 1
                list (1:3,angle(icel_local),icel_local) = (/idx1, idx2, idx3/)
                force(    angle(icel_local),icel_local) = &
                                    gromol%angls(k)%kt * JOU2CAL * 0.5_wp
                theta(    angle(icel_local),icel_local) = &
                                    gromol%angls(k)%theta_0 * RAD
              end if

            end if

          end do

        else

          idx1 = 2 + ioffset
          idx2 = 1 + ioffset
          idx3 = 3 + ioffset

          icel1 = id_g2l(1,idx1)
          icel2 = id_g2l(1,idx3)

          if (icel1 /= 0 .and. icel2 /= 0) then

            icel_local = cell_pair(icel1,icel2)

            if (icel_local > 0 .and. icel_local <= ncel) then

              if (angle(icel_local)+1 > MaxAngle) &
                call error_msg('Setup_Enefunc_Angl> Too many angles.') 

              nangl_a = nangl_a + 1

              angle(icel_local) = angle(icel_local) + 1
              list (1:3,angle(icel_local),icel_local) = (/idx1, idx2, idx3/)
              force(    angle(icel_local),icel_local) = 0.0_wp
              theta(    angle(icel_local),icel_local) = 0.0_wp
            end if

          end if

        end if

        ioffset = ioffset + gromol%num_atoms

      end do
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(nangl_a, enefunc%num_angl_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_angl_all = nangl_a
#endif

    return

  end subroutine setup_enefunc_angl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_angl_constraint
  !> @brief        define ANGLE term for each cell in potential energy function
  !                with SETTLE constraint
  !! @authors      NT
  !! @param[in]    grotop      : GROMACS topology information
  !! @param[in]    domain      : domain information
  !! @param[in]    constraints : constraints information
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_angl_constraint(grotop, molecule, domain, &
                                           constraints, enefunc)

    ! formal arguments
    type(s_grotop),              intent(in)    :: grotop
    type(s_molecule),            intent(in)    :: molecule
    type(s_domain),      target, intent(in)    :: domain
    type(s_constraints), target, intent(in)    :: constraints
    type(s_enefunc),     target, intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k
    integer                  :: nwat, ioffset, icel_local, nangl_a
    integer                  :: i1, i2, i3, idx1, idx2, idx3, icel1, icel2
    integer                  :: ii1, ii2, ii3
    integer                  :: nangl_sys
    character(6)             :: res1, res2, res3

    type(s_grotop_mol), pointer :: gromol
    real(wp),           pointer :: force(:,:), theta(:,:)
    integer,            pointer :: angle(:), list(:,:,:)
    integer,            pointer :: ncel, cell_pair(:,:)
    integer,            pointer :: id_g2l(:,:)


    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l

    angle     => enefunc%num_angle
    list      => enefunc%angle_list
    force     => enefunc%angle_force_const
    theta     => enefunc%angle_theta_min
    nwat      =  enefunc%table%num_water

    ioffset   = 0
    nangl_sys = 0
    nangl_a   = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count

        if (gromol%settles%func == 0) then

          do k = 1, gromol%num_angls

            i1 = gromol%angls(k)%atom_idx1
            i2 = gromol%angls(k)%atom_idx2
            i3 = gromol%angls(k)%atom_idx3

            if (domain%fep_use) then
              ! FEP: If the angle is not set to any group of FEP, exclude this angle.
              ii1 = molecule%fepgrp(i1)
              ii2 = molecule%fepgrp(i2)
              ii3 = molecule%fepgrp(i3)
              if (molecule%fepgrp_angl(ii1,ii2,ii3) == 0) cycle
            end if

            res1 = gromol%atoms(i1)%residue_name
            res2 = gromol%atoms(i2)%residue_name
            res3 = gromol%atoms(i3)%residue_name

            if (res1(1:4) /= constraints%water_model .and. &
                res2(1:4) /= constraints%water_model .and. &
                res3(1:4) /= constraints%water_model) then

              idx1 = i1 + ioffset
              idx2 = i2 + ioffset
              idx3 = i3 + ioffset

              icel1 = id_g2l(1,idx1)
              icel2 = id_g2l(1,idx3)

              if (icel1 /= 0 .and. icel2 /= 0) then

                icel_local = cell_pair(icel1,icel2)

                if (icel_local > 0 .and. icel_local <= ncel) then

                  if (angle(icel_local)+1 > MaxAngle) &
              call error_msg('Setup_Enefunc_Angl_Constraint> Too many angles.')

                  nangl_a = nangl_a + 1

                  angle(icel_local) = angle(icel_local) + 1
                  list (1:3,angle(icel_local),icel_local) = (/idx1,idx2,idx3/)
                  force(    angle(icel_local),icel_local) = &
                                      gromol%angls(k)%kt * JOU2CAL * 0.5_wp
                  theta(    angle(icel_local),icel_local) = &
                                      gromol%angls(k)%theta_0 * RAD
                end if

              end if

            end if

          end do

          nangl_sys = nangl_sys + gromol%num_angls

        else

          nangl_sys = nangl_sys + 1

        end if

        ioffset   = ioffset   + gromol%num_atoms
        
      end do
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(nangl_a, enefunc%num_angl_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_angl_all = nangl_a
#endif

    if (enefunc%num_angl_all /= (nangl_sys - nwat)) &
      call error_msg( &
        'Setup_Enefunc_Angl_Constraint> Some angle paremeters are missing.')


    return

  end subroutine setup_enefunc_angl_constraint

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_dihe
  !> @brief        define DIHEDRAL term in potential energy function
  !! @authors      NT
  !! @param[in]    grotop   : GROMACS TOP information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_dihe(grotop, molecule, domain, enefunc)

    ! formal arguments
    type(s_grotop),  target, intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k
    integer                  :: ioffset, found
    integer                  :: idx1, idx2, idx3, idx4, icel1, icel2, icel_local
    integer                  :: ii1, ii2, ii3, ii4

    type(s_grotop_mol), pointer :: gromol
    real(wp),           pointer :: force(:,:), phase(:,:)
    integer,            pointer :: dihe(:), list(:,:,:), period(:,:)
    integer,            pointer :: ncel, cell_pair(:,:)
    integer,            pointer :: id_g2l(:,:)
    integer,            pointer :: ndihe
    integer,            pointer :: notation


    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l

    dihe      => enefunc%num_dihedral
    list      => enefunc%dihe_list
    force     => enefunc%dihe_force_const
    phase     => enefunc%dihe_phase
    period    => enefunc%dihe_periodicity
    notation  => enefunc%notation_14types
    notation  = 100

    ioffset   = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count

        do k = 1, gromol%num_dihes

          if (gromol%dihes(k)%func /= 1 .and. gromol%dihes(k)%func /= 4 &
              .and. gromol%dihes(k)%func /= 9) &
            cycle

          idx1 = gromol%dihes(k)%atom_idx1 + ioffset
          idx2 = gromol%dihes(k)%atom_idx2 + ioffset
          idx3 = gromol%dihes(k)%atom_idx3 + ioffset
          idx4 = gromol%dihes(k)%atom_idx4 + ioffset

          if (domain%fep_use) then
            ! FEP: If the dihe is not set to any group of FEP, exclude this dihe.
            ii1 = molecule%fepgrp(idx1)
            ii2 = molecule%fepgrp(idx2)
            ii3 = molecule%fepgrp(idx3)
            ii4 = molecule%fepgrp(idx4)
            if (molecule%fepgrp_dihe(ii1,ii2,ii3,ii4) == 0) cycle
          end if

          icel1 = id_g2l(1,idx1)
          icel2 = id_g2l(1,idx4)

          if (icel1 /= 0 .and. icel2 /= 0) then

            icel_local = cell_pair(icel1,icel2)

            if (icel_local > 0 .and. icel_local <= ncel) then

              ndihe => dihe(icel_local)
              ndihe = ndihe + 1

              if (ndihe > MaxDihe) &
                call error_msg('Setup_Enefunc_Dihe> Too many dihedrals.') 

              list (1:4,ndihe,icel_local) = (/idx1, idx2, idx3, idx4/)
              force (   ndihe,icel_local) = gromol%dihes(k)%kp * JOU2CAL
              phase (   ndihe,icel_local) = gromol%dihes(k)%ps * RAD
              period(   ndihe,icel_local) = gromol%dihes(k)%multiplicity
              if (period(ndihe,icel_local) >  enefunc%notation_14types) &
              call error_msg('Setup_Enefunc_Dihe> Too many periodicity.')
            end if

          end if

        end do

        ioffset = ioffset + gromol%num_atoms

      end do
    end do

    found = 0 
    do i = 1, ncel
      found = found + dihe(i)
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_dihe_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_dihe_all = found
#endif

    return

  end subroutine setup_enefunc_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_rb_dihe
  !> @brief        define Ryckaert-Bellemans DIHEDRAL term in potential energy
  !                function
  !! @authors      NT
  !! @param[in]    grotop   : GROMACS TOP information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_rb_dihe(grotop, molecule, domain, enefunc)

    ! formal arguments
    type(s_grotop),   target, intent(in)    :: grotop
    type(s_molecule),         intent(in)    :: molecule
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k
    integer                  :: ioffset, found
    integer                  :: idx1, idx2, idx3, idx4, icel1, icel2, icel_local
    integer                  :: ii1, ii2, ii3, ii4

    type(s_grotop_mol), pointer :: gromol
    real(wp),           pointer :: c(:,:,:)
    integer,            pointer :: dihe(:), list(:,:,:)
    integer,            pointer :: ncel, cell_pair(:,:)
    integer,            pointer :: id_g2l(:,:)
    integer,            pointer :: ndihe


    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l

    dihe      => enefunc%num_rb_dihedral
    list      => enefunc%rb_dihe_list
    c         => enefunc%rb_dihe_c

    ioffset   = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count

        do k = 1, gromol%num_dihes

          if (gromol%dihes(k)%func /= 3) &
            cycle

          idx1 = gromol%dihes(k)%atom_idx1 + ioffset
          idx2 = gromol%dihes(k)%atom_idx2 + ioffset
          idx3 = gromol%dihes(k)%atom_idx3 + ioffset
          idx4 = gromol%dihes(k)%atom_idx4 + ioffset

          if (domain%fep_use) then
            ! FEP: If the dihe is not set to any group of FEP, exclude this dihe.
            ii1 = molecule%fepgrp(idx1)
            ii2 = molecule%fepgrp(idx2)
            ii3 = molecule%fepgrp(idx3)
            ii4 = molecule%fepgrp(idx4)
            if (molecule%fepgrp_dihe(ii1,ii2,ii3,ii4) == 0) cycle
          end if

          icel1 = id_g2l(1,idx1)
          icel2 = id_g2l(1,idx4)

          if (icel1 /= 0 .and. icel2 /= 0) then

            icel_local = cell_pair(icel1,icel2)

            if (icel_local > 0 .and. icel_local <= ncel) then

              ndihe => dihe(icel_local)
              ndihe = ndihe + 1

              if (ndihe > MaxDihe) &
                call error_msg('Setup_Enefunc_RB_Dihe> Too many dihedrals.') 

              list (1:4,ndihe,icel_local) = (/idx1, idx2, idx3, idx4/)
              c    (1:6,ndihe,icel_local) = gromol%dihes(k)%c(1:6) * JOU2CAL
            end if

          end if

        end do

        ioffset = ioffset + gromol%num_atoms

      end do
    end do

    found = 0 
    do i = 1, ncel
      found = found + dihe(i)
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_rb_dihe_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_rb_dihe_all = found
#endif

    return

  end subroutine setup_enefunc_rb_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_impr
  !> @brief        define improper DIHEDRAL term in potential energy function
  !! @authors      NT
  !! @param[in]    grotop   : GROMACS TOP information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_impr(grotop, molecule, domain, enefunc)

    ! formal arguments
    type(s_grotop),  target, intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k
    integer                  :: ioffset, found
    integer                  :: idx1, idx2, idx3, idx4, icel1, icel2, icel_local
    integer                  :: ii1, ii2, ii3, ii4

    type(s_grotop_mol), pointer :: gromol
    real(wp),           pointer :: force(:,:), phase(:,:)
    integer,            pointer :: impr(:), list(:,:,:), period(:,:)
    integer,            pointer :: ncel, cell_pair(:,:)
    integer,            pointer :: id_g2l(:,:)
    integer,            pointer :: nimpr


    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l

    impr      => enefunc%num_improper
    list      => enefunc%impr_list
    force     => enefunc%impr_force_const
    phase     => enefunc%impr_phase
    period    => enefunc%impr_periodicity

    ioffset   = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count

        do k = 1, gromol%num_dihes

          if (gromol%dihes(k)%func /= 2) &
            cycle

          idx1 = gromol%dihes(k)%atom_idx1 + ioffset
          idx2 = gromol%dihes(k)%atom_idx2 + ioffset
          idx3 = gromol%dihes(k)%atom_idx3 + ioffset
          idx4 = gromol%dihes(k)%atom_idx4 + ioffset

          if (domain%fep_use) then
            ! FEP: If the dihe is not set to any group of FEP, exclude this dihe.
            ii1 = molecule%fepgrp(idx1)
            ii2 = molecule%fepgrp(idx2)
            ii3 = molecule%fepgrp(idx3)
            ii4 = molecule%fepgrp(idx4)
            if (molecule%fepgrp_dihe(ii1,ii2,ii3,ii4) == 0) cycle
          end if

          icel1 = id_g2l(1,idx1)
          icel2 = id_g2l(1,idx4)

          if (icel1 /= 0 .and. icel2 /= 0) then

            icel_local = cell_pair(icel1,icel2)

            if (icel_local > 0 .and. icel_local <= ncel) then

              nimpr => impr(icel_local)
              nimpr = nimpr + 1

              if (nimpr > MaxImpr) &
                call error_msg('Setup_Enefunc_Impr> Too many impropers.') 

              list (1:4,nimpr,icel_local) = (/idx1, idx2, idx3, idx4/)
              force (   nimpr,icel_local) = gromol%dihes(k)%kp & 
                                            * JOU2CAL * 0.5_wp
              phase (   nimpr,icel_local) = gromol%dihes(k)%ps * RAD
              period(   nimpr,icel_local) = gromol%dihes(k)%multiplicity
            end if

          end if

        end do

        ioffset = ioffset + gromol%num_atoms

      end do
    end do

    found = 0 
    do i = 1, ncel
      found = found + impr(i)
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_impr_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_impr_all = found
#endif

    return

  end subroutine setup_enefunc_impr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_nonb
  !> @brief        define NON-BOND term in potential energy function
  !! @authors      NT
  !! @param[in]    grotop      : GROMACS information
  !! @param[in]    constraints : constraints information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_nonb(grotop, molecule, constraints, domain, enefunc)

    ! formal arguments
    type(s_grotop),          intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
    type(s_constraints),     intent(in)    :: constraints
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    real(wp)                 :: eps, sig, ei, ej, si, sj
    real(wp)                 :: c6i, c6j, c12i, c12j, c6, c12
    real(wp)                 :: vi, vj, wi, wj, vij, wij
    integer                  :: nnonb, ncel, i, j, k, excl_level
    integer                  :: ix, jx, kx
    integer                  :: cls_local

    integer,    allocatable  :: check_cls(:)
    integer,    allocatable  :: atmcls_map_g2l(:), atmcls_map_l2g(:)
    real(wp),   allocatable  :: nb14_lj6(:,:), nb14_lj12(:,:)
    real(wp),   allocatable  :: nonb_lj6(:,:), nonb_lj12(:,:)

    enefunc%num_atom_cls = grotop%num_atomtypes
    enefunc%fudge_lj     = grotop%defaults%fudge_lj
    enefunc%fudge_qq     = grotop%defaults%fudge_qq

    ncel                 = domain%num_cell_local
    ELECOEF              = ELECOEF_GROMACS

    ! set lennard-jones parameters
    !
    nnonb = enefunc%num_atom_cls

    allocate(check_cls(nnonb),        &
             atmcls_map_g2l(nnonb),   &
             atmcls_map_l2g(nnonb),   &
             nb14_lj6 (nnonb, nnonb), &
             nb14_lj12(nnonb, nnonb), &
             nonb_lj6 (nnonb, nnonb), &
             nonb_lj12(nnonb, nnonb))

    check_cls(1:nnonb)          = 0
    nb14_lj6 (1:nnonb, 1:nnonb) = 0.0_wp
    nb14_lj12(1:nnonb, 1:nnonb) = 0.0_wp
    nonb_lj6 (1:nnonb, 1:nnonb) = 0.0_wp
    nonb_lj12(1:nnonb, 1:nnonb) = 0.0_wp

    do i = 1, nnonb

      lj_coef(1,i) = grotop%atomtypes(i)%v * 10.0_wp
      lj_coef(2,i) = grotop%atomtypes(i)%w * JOU2CAL

      do j = 1, nnonb

        ! combination rule
        !
        if (grotop%num_nbonparms == 0) then

          vi = grotop%atomtypes(i)%v
          vj = grotop%atomtypes(j)%v
          wi = grotop%atomtypes(i)%w
          wj = grotop%atomtypes(j)%w

          if (grotop%defaults%combi_rule == 2) then

            si = vi * 10.0_wp
            sj = vj * 10.0_wp

            ei = wi * JOU2CAL
            ej = wj * JOU2CAL

            sig = (si + sj) * 0.5_wp
            eps = sqrt(ei * ej)

            c6  = 4.0_wp * eps * (sig ** 6)
            c12 = 4.0_wp * eps * (sig ** 12)

          else ! combi_rule == 1 or 3

            c6i  = vi * 1000000.0_wp * JOU2CAL
            c6j  = vj * 1000000.0_wp * JOU2CAL

            c12i = wi * 1000000.0_wp * 1000000.0_wp * JOU2CAL
            c12j = wj * 1000000.0_wp * 1000000.0_wp * JOU2CAL

            c6  = sqrt(c6i  * c6j)
            c12 = sqrt(c12i * c12j)

          end if

        else

          vij = 0.0_wp
          wij = 0.0_wp

          do k = 1, grotop%num_nbonparms
            if (grotop%atomtypes(i)%type_name == &
                  grotop%nbonparms(k)%atom_type1 .and. &
                grotop%atomtypes(j)%type_name == &
                  grotop%nbonparms(k)%atom_type2 .or.  &
                grotop%atomtypes(j)%type_name == &
                  grotop%nbonparms(k)%atom_type1 .and. &
                grotop%atomtypes(i)%type_name == &
                  grotop%nbonparms(k)%atom_type2) then

              vij = grotop%nbonparms(k)%v
              wij = grotop%nbonparms(k)%w

              exit
            end if
          end do

          if (vij == 0.0_wp .or. wij == 0.0_wp) &
            call error_msg('Setup_Enefunc_Nonb> combination is not found.')

          if (grotop%defaults%combi_rule == 2) then

            sig = vij * 10.0_wp
            eps = wij * JOU2CAL

            c6  = 4.0_wp * eps * (sig ** 6)
            c12 = 4.0_wp * eps * (sig ** 12)

          else ! combi_rule = 1 or 3

            c6  = vij * 1000000.0_wp * JOU2CAL
            c12 = wij * 1000000.0_wp * 1000000.0_wp * JOU2CAL

          end if

        end if

        ! set parameters
        !
        nb14_lj12(i,j) = c12
        nb14_lj6 (i,j) = c6

        nonb_lj12(i,j) = c12
        nonb_lj6 (i,j) = c6

      end do
    end do

    ! create native contact list
    if (enefunc%forcefield == ForcefieldAAGO) then
      call setup_enefunc_contact(grotop, domain, enefunc)
    end if

    ! check # of exclusion level
    !
    if (enefunc%forcefield == ForcefieldGROMARTINI) then

      !TODO

      enefunc%excl_level = -1

      do i = 1, grotop%num_molss
        excl_level = grotop%molss(i)%moltype%exclude_nbon
        if (enefunc%excl_level == -1) then
          enefunc%excl_level = excl_level
        else if (enefunc%excl_level /= excl_level) then
          call error_msg( &
               'Setup_Enefunc_Nonb> multiple "exclude_nbon" is not supported.')
        end if

      end do

    end if

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
    do i = 1, nnonb
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
      if (constraints%rigid_bond .or. enefunc%table%water_table) then
        domain%water%atom_cls_no(1:3)  = &
           atmcls_map_g2l(domain%water%atom_cls_no(1:3))
        enefunc%table%atom_cls_no_O    = &
           atmcls_map_g2l(enefunc%table%atom_cls_no_O)
        enefunc%table%atom_cls_no_H    = &
           atmcls_map_g2l(enefunc%table%atom_cls_no_H)
      end if
    end if

    deallocate(check_cls,      &
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

    if (enefunc%forcefield == ForcefieldAAGO) then

      call count_nonb_excl_go(.true., domain, enefunc)

    else 

      if (constraints%rigid_bond) then
        call count_nonb_excl_constraint(.true., .false., constraints, &
                                        domain, enefunc)
      else
        call count_nonb_excl(.true., domain, enefunc)
      end if

    end if

    return

  end subroutine setup_enefunc_nonb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_gro_restraints
  !> @brief        setup restraints from GROTOP information
  !! @authors      NT
  !! @param[in]    molecule : molecule information
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_gro_restraints(molecule, grotop, domain, enefunc)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_grotop),          intent(in)    :: grotop
    type(s_domain),          intent(in)    :: domain
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    type(s_enefunc)          :: ef0
    real(wp)                 :: kx, ky, kz
    integer                  :: nposres, npr_atom, max_pr_atom, natom, ioffset
    integer                  :: i, j, k, n, n2, n3, group0, func0, istart, iend
    integer                  :: ix, iatm, icel

    type(s_grotop_mol), pointer :: gromol


    ! count # of position restraints
    !
    nposres     = 0
    max_pr_atom = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol

      if (gromol%num_posress > 0) &
        nposres = nposres + 1

      npr_atom = grotop%molss(i)%count * gromol%num_posress
      max_pr_atom = max(max_pr_atom,npr_atom)
    end do

    if (nposres == 0) &
      return

    enefunc%restraint = .true.


    ! setup EneFuncRefc
    !
    if (.not. allocated(enefunc%restraint_refcoord)) then
      n = molecule%num_atoms
      call alloc_enefunc(enefunc, EneFuncRefc, n)
      enefunc%restraint_refcoord(1:3,1:n) = molecule%atom_coord(1:3,1:n)
    end if


    ! setup EneFuncRefg
    !
    if (allocated(enefunc%restraint_numatoms)) then

      n  = size(enefunc%restraint_numatoms)
      n2 = size(enefunc%restraint_atomlist(:,1))
      call alloc_enefunc(ef0, EneFuncRefg, n, n2)
      ef0%restraint_numatoms     (1:n) = enefunc%restraint_numatoms     (1:n)
      ef0%restraint_atomlist(1:n2,1:n) = enefunc%restraint_atomlist(1:n2,1:n)
      ef0%restraint_masscoef(1:n2,1:n) = enefunc%restraint_masscoef(1:n2,1:n)

      call alloc_enefunc(enefunc, EneFuncRefg, n+nposres, max(n2,max_pr_atom))
      enefunc%restraint_numatoms     (1:n) = ef0%restraint_numatoms     (1:n)
      enefunc%restraint_atomlist(1:n2,1:n) = ef0%restraint_atomlist(1:n2,1:n)
      enefunc%restraint_masscoef(1:n2,1:n) = ef0%restraint_masscoef(1:n2,1:n)

      call dealloc_enefunc(ef0, EneFuncRefg)

    else

      n = 0
      call alloc_enefunc(enefunc, EneFuncRefg, nposres, max_pr_atom)

    end if


    group0 = n

    natom   = 0
    nposres = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol

      if (gromol%num_posress > 0) &
        nposres = nposres + 1

      npr_atom = 0
      do j = 1, grotop%molss(i)%count
        ioffset = natom
        natom   = natom + gromol%num_atoms

        do k = 1, gromol%num_posress
          npr_atom = npr_atom + 1
          enefunc%restraint_atomlist(npr_atom, group0+nposres) = &
               gromol%posress(k)%atom_idx + ioffset

          if (k > 1 .and. main_rank) then
            if (gromol%posress(k-1)%kx /= gromol%posress(k)%kx .or. &
                gromol%posress(k-1)%ky /= gromol%posress(k)%ky .or. &
                gromol%posress(k-1)%kz /= gromol%posress(k)%kz) then
              write(MsgOut,'(a)') 'Setup_Enefunc_Gro_Restraints> WARNING:'
              write(MsgOut,'(a,a,a)') &
  '   different restraint constant between foreach atoms is not supported. [', &
  trim(grotop%molss(i)%moltype%name), ']'
              write(MsgOut,'(a)') ' '
            end if
          end if

        end do

      end do

      if (npr_atom > 0) &
        enefunc%restraint_numatoms(group0+nposres) = npr_atom

    end do


    ! setup EneFuncReff
    !
    if (allocated(enefunc%restraint_kind)) then

      n  = size(enefunc%restraint_kind)
      n2 = size(enefunc%restraint_grouplist(:,1))
      n3 = max(int(n2/2),1)
      call alloc_enefunc(ef0, EneFuncReff, n, n2)
      ef0%restraint_kind          (1:n) = enefunc%restraint_kind          (1:n)
      ef0%restraint_grouplist(1:n2,1:n) = enefunc%restraint_grouplist(1:n2,1:n)
      ef0%restraint_const    (1:4, 1:n) = enefunc%restraint_const    (1:4, 1:n)
      ef0%restraint_ref      (1:2, 1:n) = enefunc%restraint_ref      (1:2, 1:n)
      ef0%restraint_funcgrp       (1:n) = enefunc%restraint_funcgrp       (1:n)
      ef0%restraint_exponent_func (1:n) = enefunc%restraint_exponent_func (1:n)
      ef0%restraint_exponent_dist (1:n3,1:n) &
                                   = enefunc%restraint_exponent_dist (1:n3,1:n)
      ef0%restraint_weight_dist   (1:n3,1:n) &
                                   = enefunc%restraint_weight_dist   (1:n3,1:n)

      call alloc_enefunc(enefunc, EneFuncReff, n+nposres, max(n2,1))
      enefunc%restraint_kind          (1:n) = ef0%restraint_kind          (1:n)
      enefunc%restraint_grouplist(1:n2,1:n) = ef0%restraint_grouplist(1:n2,1:n)
      enefunc%restraint_const    (1:4, 1:n) = ef0%restraint_const    (1:4, 1:n)
      enefunc%restraint_ref      (1:2, 1:n) = ef0%restraint_ref      (1:2, 1:n)
      enefunc%restraint_funcgrp       (1:n) = ef0%restraint_funcgrp       (1:n)
      enefunc%restraint_exponent_func (1:n) = ef0%restraint_exponent_func (1:n)
      enefunc%restraint_exponent_dist (1:n3,1:n) &
                                       = ef0%restraint_exponent_dist (1:n3,1:n)
      enefunc%restraint_weight_dist   (1:n3,1:n) &
                                       = ef0%restraint_weight_dist   (1:n3,1:n)

      call dealloc_enefunc(ef0, EneFuncReff)

    else

      n = 0
      call alloc_enefunc(enefunc, EneFuncReff, nposres, 1)

    end if


    func0 = n

    nposres = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol

      if (gromol%num_posress == 0) &
        cycle

      nposres = nposres + 1

      kx = gromol%posress(1)%kx * 0.01_wp * JOU2CAL * 0.5_wp
      ky = gromol%posress(1)%ky * 0.01_wp * JOU2CAL * 0.5_wp
      kz = gromol%posress(1)%kz * 0.01_wp * JOU2CAL * 0.5_wp

      enefunc%restraint_kind         (func0+nposres) = RestraintsFuncPOSI
      enefunc%restraint_funcgrp      (func0+nposres) = 1
      enefunc%restraint_grouplist  (1,func0+nposres) = group0 + nposres
      enefunc%restraint_const      (1,func0+nposres) = 1.0_wp
      enefunc%restraint_const      (2,func0+nposres) = kx
      enefunc%restraint_const      (3,func0+nposres) = ky
      enefunc%restraint_const      (4,func0+nposres) = kz
      enefunc%restraint_exponent_func(func0+nposres) = 2

    end do

    enefunc%num_restraintgroups = group0+nposres
    enefunc%num_restraintfuncs  = func0+nposres


    ! setup domain restraints
    !
    call alloc_enefunc(enefunc, EneFuncRest, domain%num_cell_local, 0)

    do i = 1, nposres

      do ix = 1, enefunc%restraint_numatoms( &
                    enefunc%restraint_grouplist(1,func0+i))

        iatm = enefunc%restraint_atomlist( &
                ix, enefunc%restraint_grouplist(1,func0+i))

        icel = domain%id_g2l(1,iatm)

        if (icel > 0 .and. icel <= domain%num_cell_local) then
          enefunc%num_restraint(icel) = enefunc%num_restraint(icel) + 1

          n = enefunc%num_restraint(icel)
          enefunc%restraint_atom (    n,icel) = iatm
          enefunc%restraint_force(1:4,n,icel) = enefunc%restraint_const(1:4,i)
          enefunc%restraint_coord(1:3,n,icel) = &
                                          enefunc%restraint_refcoord(1:3,iatm)
        end if

      end do

    end do


    ! summary of setup enefunc_gro_restraints
    !
    if (main_rank) then

      write(MsgOut,'(A)')&
           'Setup_Enefunc_Gro_Restraints> Setup restraint functions'

      do i = func0+1, enefunc%num_restraintfuncs

        write(MsgOut,'(A,I5,A,I5)') &
          ' func  = ', i, ' kind  = ', enefunc%restraint_kind(i)

        ! summary for positional restraint
        !
        write(MsgOut,'(A,4F8.3)') &
             ' const(total, x, y, z) = ', enefunc%restraint_const(1:4,i) 
        write(MsgOut,'(A,I5)') ' exponend of function = ', &
                                          enefunc%restraint_exponent_func(i)

        write(MsgOut,'(A,I5)') ' # of groups  = ', enefunc%restraint_funcgrp(i)
        write(MsgOut,'(" grouplist: ",$)') 
        do j = 1, enefunc%restraint_funcgrp(i) 
          write(MsgOut,'(i3,$)') enefunc%restraint_grouplist(j,i)
          if (mod(j,20) == 0 .and. j /= enefunc%restraint_funcgrp(i) )  &
            write(MsgOut,'(A)') ''
        end do
        write(MsgOut,'(A)') ''

      end do

      write(MsgOut,*) ''

    end if

    return

  end subroutine setup_enefunc_gro_restraints

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    count_nonb_excl_go
  !> @brief        exclude 1-2, 1-3 interactions with Go potential
  !! @authors      JJ
  !! @param[in]    first   : flag for first call or not
  !! @param[inout] domain  : structure of domain
  !! @param[inout] enefunc : structure of enefunc
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine count_nonb_excl_go(first, domain, enefunc)

    ! formal arguments
    logical,                 intent(in)    :: first
    type(s_domain),  target, intent(inout) :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: ncell, i, ii, ix, k, i1, i2
    integer                  :: ij ,j
    integer                  :: icel, icel1, icel2, jcel
    integer                  :: ini_nb14, fin_nb14
    integer                  :: num_excl, num_nb14
    integer                  :: found1, found2
    integer                  :: list1, list2
    integer                  :: id, omp_get_thread_num
    logical                  :: duplicate

    integer,         pointer :: ncell_local, max_atom, natom(:), id_g2l(:,:)
    integer,         pointer :: cell_pair(:,:)
    integer,         pointer :: num_nonb_excl(:,:), num_nb14_calc(:,:)
    integer,         pointer :: num_nonb_excl1(:,:), num_nb14_calc1(:,:)
    integer,         pointer :: nbond(:), bondlist(:,:,:)
    integer,         pointer :: nangle(:), anglelist(:,:,:)
    integer,         pointer :: ndihedral(:), dihelist(:,:,:), pairlist(:,:)
    integer,         pointer :: nimpr(:), imprlist(:,:,:)
    integer,         pointer :: nonb_excl_list(:,:,:), nb14_calc_list(:,:,:)
    integer,         pointer :: nonb_excl_list1(:,:,:), nb14_calc_list1(:,:,:)
    integer(1),      pointer :: exclusion_mask1(:,:,:), exclusion_mask(:,:,:)


    ncell_local     => domain%num_cell_local
    max_atom        => domain%max_num_atom
    natom           => domain%num_atom
    id_g2l          => domain%id_g2l
    cell_pair       => domain%cell_pairlist1
    pairlist        => domain%cell_pairlist2

    nbond           => enefunc%num_bond
    bondlist        => enefunc%bond_list
    nangle          => enefunc%num_angle
    anglelist       => enefunc%angle_list
    ndihedral       => enefunc%num_dihedral
    dihelist        => enefunc%dihe_list
    nimpr           => enefunc%num_improper
    imprlist        => enefunc%impr_list
    num_nonb_excl   => enefunc%num_nonb_excl
    num_nonb_excl1  => enefunc%num_nonb_excl1
    num_nb14_calc   => enefunc%num_nb14_calc
    num_nb14_calc1  => enefunc%num_nb14_calc1
    nonb_excl_list  => enefunc%nonb_list
    nonb_excl_list1 => enefunc%nonb_list1
    nb14_calc_list  => enefunc%nb14_list
    nb14_calc_list1 => enefunc%nb14_list1
    exclusion_mask  => enefunc%exclusion_mask
    exclusion_mask1 => enefunc%exclusion_mask1

    ncell = domain%num_cell_local + domain%num_cell_boundary

    max_atom = 0
    do i = 1, ncell
      max_atom = max(max_atom,natom(i))
    end do

    ! initialization
    !
    num_nonb_excl1 (1:max_atom,1:ncell_local) = 0
    num_nonb_excl  (1:max_atom,1:maxcell_near) = 0
    num_nb14_calc1 (1:max_atom,1:ncell_local) = 0
    num_nb14_calc  (1:max_atom,1:maxcell_near) = 0

    do i = 1, ncell_local
      k = natom(i)
      exclusion_mask1(1:k,1:k,i) = 0
      do i1 = 1, k
        exclusion_mask1(i1,i1,i) = 1
      end do
    end do
    do ij = 1, maxcell_near
      i = cell_pair(1,ij)
      j = cell_pair(2,ij)
      i1 = max(natom(i),natom(j))
      exclusion_mask(1:i1,1:i1,ij) = 0
    end do

    ! exclude 1-2 interaction
    !
    if (enefunc%excl_level > 0) then

      do i = 1, ncell_local
        do ix = 1, nbond(i)

          list1 = bondlist(1,ix,i)
          list2 = bondlist(2,ix,i)
          icel1 = id_g2l(1,list1)
          icel2 = id_g2l(1,list2)
          i1    = id_g2l(2,list1)
          i2    = id_g2l(2,list2)

          if (icel1 < icel2) then
            icel = pairlist(icel2,icel1)
            num_excl = num_nonb_excl(i1,icel) + 1
            num_nonb_excl(i1,icel) = num_excl
            nonb_excl_list(num_excl,i1,icel) = i2
            exclusion_mask(i2,i1,icel) = 1

          else if (icel1 > icel2) then
            icel = pairlist(icel1,icel2)
            num_excl = num_nonb_excl(i2,icel) + 1
            nonb_excl_list(num_excl,i2,icel) = i1
            exclusion_mask(i1,i2,icel) = 1

          else if (i1 < i2) then
            num_excl = num_nonb_excl1(i1,i) + 1
            num_nonb_excl1(i1,i) = num_excl
            nonb_excl_list1(num_excl,i1,i) = i2
            exclusion_mask1(i2,i1,i) = 1
          
          else if (i1 > i2) then
            num_excl = num_nonb_excl1(i2,i) + 1
            num_nonb_excl1(i2,i) = num_excl
            nonb_excl_list1(num_excl,i2,i) = i1
            exclusion_mask1(i1,i2,i) = 1
          end if

        end do
      end do
    end if

    ! exclude 1-3 interaction
    !
    if (enefunc%excl_level > 1) then

      do i = 1, ncell_local
        do ix = 1, nangle(i)

          list1 = anglelist(1,ix,i)
          list2 = anglelist(3,ix,i)
          icel1 = id_g2l(1,list1)
          icel2 = id_g2l(1,list2)
          i1    = id_g2l(2,list1)
          i2    = id_g2l(2,list2)

          if (icel1 < icel2) then

            icel = pairlist(icel2,icel1)
            num_excl = num_nonb_excl(i1,icel)
            duplicate = .false.
            do k = 1, num_excl
              if (i2 == nonb_excl_list(k,i1,icel)) duplicate = .true.
            end do
            if (.not. duplicate) then
              num_excl = num_excl + 1
              num_nonb_excl(i1,icel) = num_excl
              nonb_excl_list(num_excl,i1,icel) = i2
              exclusion_mask(i2,i1,icel) = 1
            end if

          else if (icel1 > icel2) then

            icel = pairlist(icel1,icel2)
            num_excl = num_nonb_excl(i2,icel)
            duplicate = .false.
            do k = 1, num_excl
              if (i1 == nonb_excl_list(k,i2,icel)) duplicate = .true.
            end do
            if (.not. duplicate) then
              num_excl = num_excl + 1
              num_nonb_excl(i2,icel) = num_excl
              nonb_excl_list(num_excl,i2,icel) = i1
              exclusion_mask(i1,i2,icel) = 1
            end if

          else if (i1 < i2) then
 
            num_excl = num_nonb_excl1(i1,i)
            duplicate = .false.
            do k = 1, num_excl
              if (i2 == nonb_excl_list1(k,i1,i)) duplicate = .true.
            end do
            if (.not. duplicate) then
              num_excl = num_excl + 1
              nonb_excl_list1(num_excl,i1,i) = i2
              exclusion_mask1(i2,i1,i) = 1
            end if

          else if (i1 > i2) then

            num_excl = num_nonb_excl1(i2,i)
            duplicate = .false.
            do k = 1, num_excl
              if (i1 == nonb_excl_list1(k,i2,i)) duplicate = .true.
            end do
            if (.not. duplicate) then
              num_excl = num_excl + 1
              nonb_excl_list1(num_excl,i2,i) = i1
              exclusion_mask1(i1,i2,i) = 1
            end if

          end if
        end do
      end do

    end if

    ! count 1-4 interaction
    !
    if (enefunc%excl_level > 2) then

      do ii = 1, 2

        if (ii == 2) then
          ndihedral => enefunc%num_rb_dihedral
          dihelist => enefunc%rb_dihe_list
        end if

        do i = 1, ncell_local
          do ix = 1, ndihedral(i)

            list1 = dihelist(1,ix,i)
            list2 = dihelist(4,ix,i)
            icel1 = id_g2l(1,list1)
            icel2 = id_g2l(1,list2)
            i1    = id_g2l(2,list1)
            i2    = id_g2l(2,list2)

            if (icel1 < icel2) then

              icel = pairlist(icel2,icel1)
              num_nb14 = num_nb14_calc(i1,icel)
              duplicate = .false.
              do k = 1, num_nonb_excl(i1,icel)
                if (i2 == nonb_excl_list(k,i1,icel)) duplicate = .true.
              end do
              do k = 1, num_nb14
                if (i2 == nb14_calc_list(k,i1,icel)) duplicate = .true.
              end do
              if (.not. duplicate) then
                num_nb14 = num_nb14 + 1
                num_nb14_calc(i1,icel) = num_nb14
                nb14_calc_list(num_nb14,i1,icel) = i2
                exclusion_mask(i2,i1,icel) = 1
              end if

            else if (icel1 > icel2) then

              icel = pairlist(icel1,icel2)
              num_nb14 = num_nb14_calc(i2,icel)
              duplicate = .false.
              do k = 1, num_nonb_excl(i2,icel)
                if (i1 == nonb_excl_list(k,i2,icel)) duplicate = .true.
              end do
              do k = 1, num_nb14
                if (i1 == nb14_calc_list(k,i2,icel)) duplicate = .true.
              end do
              if (.not. duplicate) then
                num_nb14 = num_nb14 + 1
                num_nb14_calc(i2,icel) = num_nb14
                nb14_calc_list(num_nb14,i2,icel) = i1
                exclusion_mask(i1,i2,icel) = 1
              end if

            else if (i1 < i2) then

              num_nb14 = num_nb14_calc1(i1,i)
              duplicate = .false.
              do k = 1, num_nb14
                if (i2 == nb14_calc_list1(k,i1,i)) duplicate = .true.
              end do
              do k = 1, num_nonb_excl1(i1,i)
                if (i2 == nonb_excl_list1(k,i1,i)) duplicate = .true.
              end do
              if (.not. duplicate) then
                num_nb14 = num_nb14 + 1
                num_nb14_calc1(i1,i) = num_nb14
                nb14_calc_list1(num_nb14,i1,i) = i2
                exclusion_mask1(i2,i1,i) = 1
              end if

            else if (i1 > i2) then

              num_nb14 = num_nb14_calc1(i2,i)
              duplicate = .false.
              do k = 1, num_nb14
                if (i1 == nb14_calc_list1(k,i2,i)) duplicate = .true.
              end do
              do k = 1, num_nonb_excl1(i2,i)
                if (i1 == nonb_excl_list1(k,i2,i)) duplicate = .true.
              end do
              if (.not. duplicate) then
                num_nb14 = num_nb14 + 1
                num_nb14_calc1(i2,i) = num_nb14
                nb14_calc_list1(num_nb14,i2,i) = i1
                exclusion_mask1(i1,i2,i) = 1
              end if

            end if 

          end do
        end do

      end do

      do i = 1, ncell_local
        do ix = 1, nimpr(i)

          list1 = imprlist(1,ix,i)
          list2 = imprlist(4,ix,i)
          icel1 = id_g2l(1,list1)
          i1    = id_g2l(2,list1)
          icel2 = id_g2l(1,list2)
          i2    = id_g2l(2,list2)

          if (icel1 < icel2) then

            icel = pairlist(icel2,icel1)
            num_nb14 = num_nb14_calc(i1,icel)
            duplicate = .false.
            do k = 1, num_nonb_excl(i1,icel)
              if (i2 == nonb_excl_list(k,i1,icel)) duplicate = .true.
            end do
            do k = 1, num_nb14
              if (i2 == nb14_calc_list(k,i1,icel)) duplicate = .true.
            end do
            if (.not. duplicate) then
              num_nb14 = num_nb14 + 1
              num_nb14_calc(i1,icel) = num_nb14
              nb14_calc_list(num_nb14,i1,icel) = i2
              exclusion_mask(i2,i1,icel) = 1
            end if

          else if (icel1 > icel2) then

            icel = pairlist(icel1,icel2)
            num_nb14 = num_nb14_calc(i2,icel)
            duplicate = .false.
            do k = 1, num_nonb_excl(i2,icel)
              if (i1 == nonb_excl_list(k,i2,icel)) duplicate = .true.
            end do
            do k = 1, num_nb14
              if (i1 == nb14_calc_list(k,i2,icel)) duplicate = .true.
            end do
            if (.not. duplicate) then
              num_nb14 = num_nb14 + 1
              num_nb14_calc(i2,icel) = num_nb14
              nb14_calc_list(num_nb14,i2,icel) = i1
              exclusion_mask(i1,i2,icel) = 1
            end if

          else if (i1 < i2) then

            num_nb14 = num_nb14_calc1(i1,i)
            duplicate = .false.
            do k = 1, num_nb14
              if (i2 == nb14_calc_list1(k,i1,i)) duplicate = .true.
            end do
            do k = 1, num_nonb_excl1(i1,i)
              if (i2 == nonb_excl_list1(k,i1,i)) duplicate = .true.
            end do
            if (.not. duplicate) then
              num_nb14 = num_nb14 + 1
              num_nb14_calc1(i1,i) = num_nb14
              nb14_calc_list1(num_nb14,i1,i) = i2
              exclusion_mask1(i2,i1,i) = 1
            end if

          else if (i1 > i2) then

            num_nb14 = num_nb14_calc1(i2,i)
            duplicate = .false.
            do k = 1, num_nb14
              if (i1 == nb14_calc_list1(k,i2,i)) duplicate = .true.
            end do
            do k = 1, num_nonb_excl1(i2,i)
              if (i1 == nonb_excl_list1(k,i2,i)) duplicate = .true.
            end do
            if (.not. duplicate) then
              num_nb14 = num_nb14 + 1
              num_nb14_calc1(i2,i) = num_nb14
              nb14_calc_list1(num_nb14,i2,i) = i1
              exclusion_mask1(i1,i2,i) = 1
            end if
 
          end if

        end do
      end do

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

  end subroutine count_nonb_excl_go

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_contact
  !> @brief        define Contact term for each cell
  !! @authors      JJ
  !! @param[in]    grotop   : GROMACS TOP information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_contact(grotop, domain, enefunc)

    ! formal arguments
    type(s_grotop),          intent(in)    :: grotop
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    real(wp)                 :: sig, eps
    integer                  :: step
    integer                  :: i, j, k
    integer                  :: ioffset, ncontact, pcontact
    integer                  :: idx1, idx2, icel1, icel2, icel_local

    type(s_grotop_mol), pointer :: gromol
    real(wp),           pointer :: lj12(:,:), lj6(:,:)
    integer,            pointer :: contact(:), list(:,:,:)
    integer,            pointer :: ncel, cell_pair(:,:)
    integer,            pointer :: id_g2l(:,:)


    ncel      => domain%num_cell_local
    cell_pair => domain%cell_pair
    id_g2l    => domain%id_g2l

    contact   => enefunc%num_contact
    list      => enefunc%contact_list
    lj12      => enefunc%contact_lj12
    lj6       => enefunc%contact_lj6

    do step = 1, 2

      ioffset   = 0
      ncontact  = 0
      contact(1:ncel) = 0

      do i = 1, grotop%num_molss
        gromol => grotop%molss(i)%moltype%mol
        do j = 1, grotop%molss(i)%count

          do k = 1, gromol%num_pairs

            idx1 = gromol%pairs(k)%atom_idx1 + ioffset
            idx2 = gromol%pairs(k)%atom_idx2 + ioffset

            icel1 = id_g2l(1,idx1)
            icel2 = id_g2l(1,idx2)

            if (icel1 /= 0 .and. icel2 /= 0) then

              icel_local = cell_pair(icel1,icel2)

              if (icel_local > 0 .and. icel_local <= ncel) then

                ncontact = ncontact + 1
                pcontact = contact(icel_local)
                pcontact = pcontact + 1
                contact(icel_local) = pcontact

                if (step == 2) then
                  enefunc%contact_list(1,contact(icel_local),icel_local) = idx1
                  enefunc%contact_list(2,contact(icel_local),icel_local) = idx2
                  if (grotop%defaults%combi_rule == 2) then
                    sig = 10.0_wp * gromol%pairs(k)%v
                    eps = JOU2CAL * gromol%pairs(k)%w
                    enefunc%contact_lj12(pcontact,icel_local) = eps * (sig**12)
                    enefunc%contact_lj6 (pcontact,icel_local) = &
                                         2.0_wp * eps * (sig** 6)
                  else
                    enefunc%contact_lj12(pcontact,icel_local) = &
                                         gromol%pairs(k)%w *  1.0E12_wp * JOU2CAL
                    enefunc%contact_lj6 (pcontact,icel_local) = &
                                         gromol%pairs(k)%v *  1.0E6_wp * JOU2CAL
                  end if
                end if

              end if
            end if

          end do

          ioffset = ioffset + gromol%num_atoms

        end do
      end do

      if (step == 1) then

        MaxContact = 0
        do i = 1, ncel
          MaxContact = max(MaxContact,contact(i))
        end do
#ifdef HAVE_MPI_GENESIS
        call mpi_allreduce(mpi_in_place, MaxContact, 1, mpi_integer, &
                           mpi_max, mpi_comm_country, ierror)
#endif
        MaxContact = MaxContact * 3 / 2
        ContactMove  = MaxContact / 2
        call alloc_enefunc(enefunc, EneFuncContact,  ncel, ncel)

      end if

    end do
#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(ncontact, enefunc%num_contact_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_contact_all = ncontact
#endif

    return
  end subroutine setup_enefunc_contact

end module sp_enefunc_gromacs_mod
