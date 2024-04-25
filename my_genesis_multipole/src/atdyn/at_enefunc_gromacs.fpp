!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_enefunc_gromacs_mod
!> @brief   define potential energy functions
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_enefunc_gromacs_mod

  use at_enefunc_restraints_mod
  use at_enefunc_table_mod
  use at_enefunc_pme_mod
  use at_energy_mod
  use at_boundary_str_mod
  use at_restraints_str_mod
  use at_enefunc_str_mod
  use at_energy_str_mod
  use nbond_list_mod
  use math_libs_mod
  use molecules_str_mod
  use fileio_grotop_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! constants
  integer,   parameter :: MAXWRK   = 12
  integer,   parameter :: MAX_EXCL = 16 ! = 4(1-2) + 4x3 (1-3)
  integer,   parameter :: MAX_NB14 = 36 ! = 3x4x3 (1-4) 

  ! subroutines
  public  :: define_enefunc_gromacs
  private :: setup_enefunc_bond
  private :: setup_enefunc_angl
  private :: setup_enefunc_dihe
  private :: setup_enefunc_rb_dihe
  private :: setup_enefunc_impr
  private :: setup_enefunc_nonb
  private :: setup_enefunc_gro_restraints
  private :: count_nonb_excl
  private :: count_nonb_excl_go

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_enefunc_gromacs
  !> @brief        a driver subroutine for defining potential energy functions
  !! @authors      NT
  !! @param[in]    ene_info   : ENERGY section control parameters information
  !! @param[in]    boundary   : boundary condition
  !! @param[in]    grotop     : GROMACS parameter topology information
  !! @param[in]    molecule   : molecule information
  !! @param[in]    restraints : restraints information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_enefunc_gromacs(ene_info, boundary, grotop, &
                                    molecule, restraints, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info 
    type(s_boundary),        intent(in)    :: boundary
    type(s_grotop),          intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
    type(s_restraints),      intent(in)    :: restraints
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: alloc_stat
    integer                  :: i, num_nonb_excl, num_nb14_calc
    integer                  :: nwork


    ! bond
    !
    call setup_enefunc_bond(grotop, enefunc)

    ! angle
    !
    call setup_enefunc_angl(grotop, enefunc)

    ! dihedral
    !
    call setup_enefunc_dihe(grotop, enefunc)

    ! Ryckaert-Bellemans dihedral
    !
    call setup_enefunc_rb_dihe(grotop, enefunc)

    ! improper
    !
    call setup_enefunc_impr(grotop, enefunc)

    ! nonbonded
    !
    call setup_enefunc_nonb(ene_info, grotop, molecule, enefunc)

    ! lookup table
    !
    call setup_enefunc_table(ene_info, molecule, enefunc)

    ! PME
    !
    call define_enefunc_pme(ene_info, boundary, molecule, enefunc)

    ! restraints
    !
    call setup_enefunc_restraints(molecule, restraints, enefunc)

    ! restraints (gromacs)
    !
    call setup_enefunc_gro_restraints(molecule, grotop, enefunc)

    ! allocate working array
    !
    nwork = max(molecule%num_atoms, &
                enefunc%num_bonds,  &
                enefunc%num_angles, &
                enefunc%num_ureys,  &
                enefunc%num_dihedrals, &
                enefunc%num_rb_dihedrals, &
                enefunc%num_impropers, &
                enefunc%num_cmaps*2,   &
                enefunc%num_restraintfuncs, &
                enefunc%num_restraintgroups)

    allocate(enefunc%work(MAXWRK,nwork),stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    ! write summary of energy function
    !
    if (main_rank) then

      if (.not. ene_info%table) then

        num_nonb_excl = 0
        num_nb14_calc = 0

        do i = 1, molecule%num_atoms
          num_nonb_excl = num_nonb_excl + enefunc%num_nonb_excl(i)
          num_nb14_calc = num_nb14_calc + enefunc%num_nb14_calc(i)
        end do

      end if

      write(MsgOut,'(A)') &
           'Define_Enefunc_Gromacs> Number of Interactions in Each Term'
      write(MsgOut,'(A20,I10,A20,I10)')                         &
           '  bond_ene        = ', enefunc%num_bonds,           &
           '  angle_ene       = ', enefunc%num_angles
      write(MsgOut,'(A20,I10,A20,I10)')                         &
           '  torsion_ene     = ', enefunc%num_dihedrals,       &
           '  rb_torsion_ene  = ', enefunc%num_rb_dihedrals
      write(MsgOut,'(A20,I10)')                                 &
           '  improper_ene    = ', enefunc%num_impropers
      if (.not. ene_info%table) then
        write(MsgOut,'(A20,I10,A20,I10)')                       &
           '  nb_exclusions   = ', num_nonb_excl,               &
           '  nb14_calc       = ', num_nb14_calc
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
  !> @brief        define BOND term in potential energy function
  !! @authors      NT
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_bond(grotop, enefunc)

    ! formal arguments
    type(s_grotop),          intent(in)    :: grotop
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k
    integer                  :: natom, nbond, ioffset
    integer                  :: istart, iend

    type(s_grotop_mol), pointer :: gromol


    nbond = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count

        if (gromol%settles%func == 0) then
          nbond = nbond + gromol%num_bonds
        else
          nbond = nbond + 3
        end if

      end do
    end do

    call alloc_enefunc(enefunc, EneFuncBond, nbond)

    natom = 0
    nbond = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        ioffset = natom
        natom   = natom + gromol%num_atoms

        if (gromol%settles%func == 0) then

          do k = 1, gromol%num_bonds
            nbond = nbond + 1
            enefunc%bond_list(1, nbond) = gromol%bonds(k)%atom_idx1 + ioffset
            enefunc%bond_list(2, nbond) = gromol%bonds(k)%atom_idx2 + ioffset
            enefunc%bond_force_const(nbond) = &
                                gromol%bonds(k)%kb * 0.01_wp * JOU2CAL * 0.5_wp
            enefunc%bond_dist_min   (nbond) = &
                                gromol%bonds(k)%b0 * 10.0_wp
          end do

        else

          ! settle O-H (1)
          nbond = nbond + 1
          enefunc%bond_list(1, nbond)     = 1 + ioffset
          enefunc%bond_list(2, nbond)     = 2 + ioffset
          enefunc%bond_force_const(nbond) = 0.0_wp
          enefunc%bond_dist_min(nbond)    = gromol%settles%doh * 10.0_wp

          ! settle O-H (2)
          nbond = nbond + 1
          enefunc%bond_list(1, nbond)     = 1 + ioffset
          enefunc%bond_list(2, nbond)     = 3 + ioffset
          enefunc%bond_force_const(nbond) = 0.0_wp
          enefunc%bond_dist_min(nbond)    = gromol%settles%doh * 10.0_wp

          ! settle H-H
          nbond = nbond + 1
          enefunc%bond_list(1, nbond)     = 2 + ioffset
          enefunc%bond_list(2, nbond)     = 3 + ioffset
          enefunc%bond_force_const(nbond) = 0.0_wp
          enefunc%bond_dist_min(nbond)    = gromol%settles%dhh * 10.0_wp

        end if

      end do
    end do

    enefunc%num_bonds = nbond


    call get_loop_index(enefunc%num_bonds, istart, iend)

    enefunc%istart_bond = istart
    enefunc%iend_bond   = iend

    return

  end subroutine setup_enefunc_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_angl
  !> @brief        define ANGLE term in potential energy function
  !! @authors      NT
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_angl(grotop, enefunc)

    ! formal arguments
    type(s_grotop),          intent(in)    :: grotop
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k
    integer                  :: natom, nangl, ioffset
    integer                  :: istart, iend

    type(s_grotop_mol), pointer :: gromol


    nangl = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count

        if (gromol%settles%func == 0) then
          nangl = nangl + gromol%num_angls
        else
          nangl = nangl + 1
        end if

      end do
    end do

    call alloc_enefunc(enefunc, EneFuncAngl, nangl)

    natom = 0
    nangl = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        ioffset = natom
        natom   = natom + gromol%num_atoms

        if (gromol%settles%func == 0) then

          do k = 1, gromol%num_angls
            nangl = nangl + 1
            enefunc%angl_list(1, nangl) = gromol%angls(k)%atom_idx1 + ioffset
            enefunc%angl_list(2, nangl) = gromol%angls(k)%atom_idx2 + ioffset
            enefunc%angl_list(3, nangl) = gromol%angls(k)%atom_idx3 + ioffset
            enefunc%angl_force_const(nangl) = &
                                          gromol%angls(k)%kt * JOU2CAL * 0.5_wp
            enefunc%angl_theta_min  (nangl) = &
                                          gromol%angls(k)%theta_0 * RAD
          end do

        else

          nangl = nangl + 1
          enefunc%angl_list(1, nangl) = 2 + ioffset
          enefunc%angl_list(2, nangl) = 1 + ioffset
          enefunc%angl_list(3, nangl) = 3 + ioffset
          enefunc%angl_force_const(nangl) = 0.0_wp
          enefunc%angl_theta_min  (nangl) = 0.0_wp

        end if

      end do
    end do

    enefunc%num_angles = nangl


    call get_loop_index(enefunc%num_angles, istart, iend)

    enefunc%istart_angle = istart
    enefunc%iend_angle   = iend

    return

  end subroutine setup_enefunc_angl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_dihe
  !> @brief        define proper DIHEDRAL term in potential energy function
  !! @authors      NT
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_dihe(grotop, enefunc)

    ! formal arguments
    type(s_grotop),          intent(in)    :: grotop
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k
    integer                  :: natom, ndihe, ioffset
    integer                  :: istart, iend

    type(s_grotop_mol), pointer :: gromol


    ndihe = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        do k = 1, gromol%num_dihes
          if (gromol%dihes(k)%func /= 1 .and. gromol%dihes(k)%func /= 4 &
          .and. gromol%dihes(k)%func /= 9 ) &
            cycle
          ndihe = ndihe + 1
        end do
      end do
    end do

    call alloc_enefunc(enefunc, EneFuncDihe, ndihe)
    
    natom = 0
    ndihe = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        ioffset = natom
        natom   = natom + gromol%num_atoms
        do k = 1, gromol%num_dihes
          if (gromol%dihes(k)%func /= 1 .and. gromol%dihes(k)%func /= 4 &
          .and. gromol%dihes(k)%func /= 9 ) &
            cycle
          ndihe = ndihe + 1
          enefunc%dihe_list    (1, ndihe) = gromol%dihes(k)%atom_idx1 + ioffset
          enefunc%dihe_list    (2, ndihe) = gromol%dihes(k)%atom_idx2 + ioffset
          enefunc%dihe_list    (3, ndihe) = gromol%dihes(k)%atom_idx3 + ioffset
          enefunc%dihe_list    (4, ndihe) = gromol%dihes(k)%atom_idx4 + ioffset
          enefunc%dihe_force_const(ndihe) = gromol%dihes(k)%kp * JOU2CAL
          enefunc%dihe_phase      (ndihe) = gromol%dihes(k)%ps * RAD
          enefunc%dihe_periodicity(ndihe) = gromol%dihes(k)%multiplicity
        end do
      end do
    end do

    enefunc%num_dihedrals = ndihe


    call get_loop_index(enefunc%num_dihedrals, istart, iend)

    enefunc%istart_dihedral = istart
    enefunc%iend_dihedral   = iend

    return

  end subroutine setup_enefunc_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_rb_dihe
  !> @brief        define Ryckaert-Bellemans DIHEDRAL term in potential energy
  !!               function
  !! @authors      NT
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_rb_dihe(grotop, enefunc)

    ! formal arguments
    type(s_grotop),          intent(in)    :: grotop
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k
    integer                  :: natom, ndihe, ioffset
    integer                  :: istart, iend

    type(s_grotop_mol), pointer :: gromol


    ndihe = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        do k = 1, gromol%num_dihes
          if (gromol%dihes(k)%func /= 3) &
            cycle
          ndihe = ndihe + 1
        end do
      end do
    end do

    call alloc_enefunc(enefunc, EneFuncRBDihe, ndihe)
    
    natom = 0
    ndihe = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        ioffset = natom
        natom   = natom + gromol%num_atoms
        do k = 1, gromol%num_dihes
          if (gromol%dihes(k)%func /= 3) &
            cycle
          ndihe = ndihe + 1
          enefunc%rb_dihe_list(1, ndihe) = gromol%dihes(k)%atom_idx1 + ioffset
          enefunc%rb_dihe_list(2, ndihe) = gromol%dihes(k)%atom_idx2 + ioffset
          enefunc%rb_dihe_list(3, ndihe) = gromol%dihes(k)%atom_idx3 + ioffset
          enefunc%rb_dihe_list(4, ndihe) = gromol%dihes(k)%atom_idx4 + ioffset
          enefunc%rb_dihe_c (1:6, ndihe) = gromol%dihes(k)%c(1:6) * JOU2CAL
        end do
      end do
    end do

    enefunc%num_rb_dihedrals = ndihe


    call get_loop_index(enefunc%num_rb_dihedrals, istart, iend)

    enefunc%istart_rb_dihed = istart
    enefunc%iend_rb_dihed   = iend

    return

  end subroutine setup_enefunc_rb_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_impr
  !> @brief        define improper DIHEDRAL term in potential energy function
  !! @authors      NT
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_impr(grotop, enefunc)

    ! formal arguments
    type(s_grotop),          intent(in)    :: grotop
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k
    integer                  :: natom, nimpr, ioffset
    integer                  :: istart, iend

    type(s_grotop_mol), pointer :: gromol


    nimpr = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        do k = 1, gromol%num_dihes
          if (gromol%dihes(k)%func /= 2) &
            cycle
          nimpr = nimpr + 1
        end do
      end do
    end do

    call alloc_enefunc(enefunc, EneFuncImpr, nimpr)
    
    natom = 0
    nimpr = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        ioffset = natom
        natom   = natom + gromol%num_atoms
        do k = 1, gromol%num_dihes
          if (gromol%dihes(k)%func /= 2) &
            cycle
          nimpr = nimpr + 1
          enefunc%impr_list    (1, nimpr) = gromol%dihes(k)%atom_idx1 + ioffset
          enefunc%impr_list    (2, nimpr) = gromol%dihes(k)%atom_idx2 + ioffset
          enefunc%impr_list    (3, nimpr) = gromol%dihes(k)%atom_idx3 + ioffset
          enefunc%impr_list    (4, nimpr) = gromol%dihes(k)%atom_idx4 + ioffset
          enefunc%impr_phase      (nimpr) = gromol%dihes(k)%ps * RAD
          enefunc%impr_force_const(nimpr) = gromol%dihes(k)%kp &
                                            * JOU2CAL * 0.5_wp
          enefunc%impr_periodicity(nimpr) = 0
        end do
      end do
    end do

    enefunc%num_impropers = nimpr


    call get_loop_index(enefunc%num_impropers, istart, iend)

    enefunc%istart_improper = istart
    enefunc%iend_improper   = iend

    return

  end subroutine setup_enefunc_impr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_nonb
  !> @brief        define NON-BOND term in potential energy function
  !! @authors      NT
  !! @param[in]    ene_info : ENERGY section control parameters information
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_nonb(ene_info, grotop, molecule, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_grotop),          intent(in)    :: grotop
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    real(wp)                 :: eps, sig, ei, ej, si, sj
    real(wp)                 :: c6i, c6j, c12i, c12j, c6, c12
    real(wp)                 :: vi, vj, wi, wj, vij, wij
    integer                  :: nnonb, i, j, k
    integer                  :: nexcl, natom, npairs, ioffset, excl_level
    integer                  :: istart, iend

    integer,        allocatable :: excls(:,:), nrexc(:)
    type(s_grotop_mol), pointer :: gromol


    enefunc%num_atom_cls = grotop%num_atomtypes
    enefunc%fudge_lj     = grotop%defaults%fudge_lj
    enefunc%fudge_qq     = grotop%defaults%fudge_qq

    ELECOEF = ELECOEF_GROMACS

    ! set lennard-jones parameters
    !
    nnonb = enefunc%num_atom_cls

    call alloc_enefunc(enefunc, EneFuncNbon, nnonb)

    do i = 1, nnonb
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

            c6i  = vi * 1.0E6_wp * JOU2CAL
            c6j  = vj * 1.0E6_wp * JOU2CAL

            c12i = wi * 1.0E12_wp * JOU2CAL
            c12j = wj * 1.0E12_wp * JOU2CAL

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
        enefunc%nb14_lj12(i,j) = c12
        enefunc%nb14_lj10(i,j) = 0.0_wp
        enefunc%nb14_lj6 (i,j) = c6

        enefunc%nonb_lj12(i,j) = c12
        enefunc%nonb_lj10(i,j) = 0.0_wp
        enefunc%nonb_lj6 (i,j) = c6

      end do
    end do

    ! create native contact list
    if (enefunc%forcefield == ForcefieldKBGO .or. &
        enefunc%forcefield == ForcefieldAAGO .or. &
        enefunc%forcefield == ForcefieldCAGO ) then

      npairs = 0
      do i = 1, grotop%num_molss
        gromol => grotop%molss(i)%moltype%mol
        do j = 1, grotop%molss(i)%count
          npairs  = npairs + gromol%num_pairs
        end do
      end do

      enefunc%num_contacts = npairs
      call alloc_enefunc(enefunc, EneFuncCntc, npairs)

      natom  = 0
      npairs = 0
      do i = 1, grotop%num_molss
        gromol => grotop%molss(i)%moltype%mol
        do j = 1, grotop%molss(i)%count
          ioffset = natom
          natom   = natom + gromol%num_atoms
          do k = 1, gromol%num_pairs
            npairs = npairs + 1

            enefunc%contact_list(1,npairs) = gromol%pairs(k)%atom_idx1  &
                                             + ioffset
            enefunc%contact_list(2,npairs) = gromol%pairs(k)%atom_idx2  &
                                             + ioffset

            if (grotop%defaults%combi_rule == 2) then
              sig = 10.0_wp * gromol%pairs(k)%v
              eps  = JOU2CAL * gromol%pairs(k)%w
           
              if (enefunc%forcefield == ForcefieldKBGO) then
                enefunc%contact_lj12(npairs) = 13.0_wp * eps * (sig ** 12)
                enefunc%contact_lj10(npairs) = 18.0_wp * eps * (sig ** 10)
                enefunc%contact_lj6 (npairs) =  4.0_wp * eps * (sig ** 6)
              else if (enefunc%forcefield == ForcefieldCAGO) then
                enefunc%contact_lj12(npairs) = 5.0_wp * eps * (sig ** 12)
                enefunc%contact_lj10(npairs) = 6.0_wp * eps * (sig ** 10)
                enefunc%contact_lj6 (npairs) = 0.0_wp
              else 
                enefunc%contact_lj12(npairs) = eps * (sig ** 12)
                enefunc%contact_lj10(npairs) = 0.0_wp
                enefunc%contact_lj6 (npairs) = 2.0_wp * eps * (sig ** 6)
              endif

            else ! combi_rule == 1 or 3
              if (enefunc%forcefield == ForcefieldCAGO) then
                enefunc%contact_lj12(npairs) = gromol%pairs(k)%w *  &
                                               1.0E12_wp * JOU2CAL
                enefunc%contact_lj10(npairs) = gromol%pairs(k)%v *  &
                                               1.0E10_wp * JOU2CAL
                enefunc%contact_lj6 (npairs) = 0.0_wp
              else
                enefunc%contact_lj12(npairs) = gromol%pairs(k)%w *  &
                                               1.0E12_wp * JOU2CAL
                enefunc%contact_lj10(npairs) = 0.0_wp
                enefunc%contact_lj6 (npairs) = gromol%pairs(k)%v *  &
                                               1.0E6_wp * JOU2CAL
              endif
            endif
          end do
        end do
      end do
      call get_loop_index(enefunc%num_contacts, istart, iend)
      enefunc%istart_contact = istart
      enefunc%iend_contact   = iend
    endif
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

    ! create exclusion list & exclude neighbours Nx
    !
    nexcl = 0
    natom = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        nexcl  = nexcl + gromol%num_excls
        natom  = natom + gromol%num_atoms
      end do
    end do

    allocate(excls(2,nexcl), nrexc(natom))

    nexcl = 0
    natom = 0

    do i = 1, grotop%num_molss
      gromol => grotop%molss(i)%moltype%mol
      do j = 1, grotop%molss(i)%count
        ioffset = natom
        natom   = natom + gromol%num_atoms

        ! exclusion list
        do k = 1, gromol%num_excls
          
          nexcl = nexcl + 1
          excls(1,nexcl) = gromol%excls(k)%atom_idx1 + ioffset
          excls(2,nexcl) = gromol%excls(k)%atom_idx2 + ioffset
        end do

        ! exclude neighbours Nx
        nrexc(1+ioffset:gromol%num_atoms+ioffset) = &
             grotop%molss(i)%moltype%exclude_nbon

      end do
    end do


    ! treatment for 1-2, 1-3, 1-4 interactions
    !
    if (enefunc%forcefield == ForcefieldKBGO .or. &
        enefunc%forcefield == ForcefieldAAGO .or. &
        enefunc%forcefield == ForcefieldCAGO ) then
       call count_nonb_excl_go(molecule, enefunc, excls, nrexc)
    else if (.not.ene_info%table) then
       call count_nonb_excl(molecule, enefunc, excls, nrexc)
    endif


    deallocate(excls)

    return

  end subroutine setup_enefunc_nonb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_nonb
  !> @brief        define NON-BOND term in potential energy function
  !! @authors      NT
  !! @param[in]    molecule : molecule information
  !! @param[in]    grotop   : GROMACS parameter topology information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_gro_restraints(molecule, grotop, enefunc)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_grotop),          intent(in)    :: grotop
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    type(s_enefunc)          :: ef0
    real(wp)                 :: kx, ky, kz
    integer                  :: nposres, npr_atom, max_pr_atom, natom, ioffset
    integer                  :: i, j, k, n, n2, n3, group0, func0, istart, iend

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

    enefunc%restraint_flag = .true.


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
   '   different restraint constant between foreach atoms is not supported. [',&
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


    ! define restraint functions for each processor
    !
    call get_loop_index(enefunc%num_restraintfuncs, istart, iend)

    enefunc%istart_restraint = istart
    enefunc%iend_restraint   = iend


    ! summary of setup enefunc_gro_restraint
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
  !  Subroutine    count_nonb_excl
  !> @brief        exclude 1-2, 1-3 interactions
  !! @authors      NT
  !! @param[in]    molecule : molecule information
  !! @param[out]   enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine count_nonb_excl(molecule, enefunc, excls, nrexc)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    integer,                 intent(in)    :: excls(:,:)
    integer,                 intent(in)    :: nrexc(:)

    ! local variables
    integer                  :: i, k, k1, k2, k3, k4
    integer                  :: num_excl, num_excl_total
    integer                  :: num_nb14, num_nb14_total
    integer                  :: max_nb14_num
    integer                  :: alloc_stat, dealloc_stat
    logical                  :: duplicate

    integer,     allocatable :: nonb_excl_list(:,:)
    integer,     allocatable :: nb14_calc_list(:,:)
    real(wp),    allocatable :: nb14_qq_scale (:,:)
    real(wp),    allocatable :: nb14_lj_scale (:,:)


    ! allocate nonbonded exclusion list and 1-4 interaction list
    !
    allocate(enefunc%num_nonb_excl(molecule%num_atoms), stat=alloc_stat)
    if (alloc_stat /=0) &
      call error_msg_alloc

    allocate(nonb_excl_list(max_excl,molecule%num_atoms),stat=alloc_stat)
    if (alloc_stat /=0) &
      call error_msg_alloc

    allocate(enefunc%num_nb14_calc(molecule%num_atoms), stat=alloc_stat)
    if (alloc_stat /=0) &
      call error_msg_alloc

    allocate(nb14_calc_list(MAX_NB14,molecule%num_atoms), &
             nb14_qq_scale (MAX_NB14,molecule%num_atoms), &
             nb14_lj_scale (MAX_NB14,molecule%num_atoms),stat=alloc_stat)
    if (alloc_stat /=0) &
      call error_msg_alloc

    enefunc%num_nonb_excl(1:molecule%num_atoms) = 0
    enefunc%num_nb14_calc(1:molecule%num_atoms) = 0

    ! exclude 1-2 interaction
    !
    num_excl_total = 0
    do k = 1, molecule%num_bonds
      k1 = molecule%bond_list(1,k)
      k2 = molecule%bond_list(2,k)

      if (nrexc(k1) < 1) then
        cycle
      end if

      if (k1 < k2) then
        num_excl = enefunc%num_nonb_excl(k1) + 1
        enefunc%num_nonb_excl(k1) = num_excl
        nonb_excl_list(num_excl,k1) = k2
        num_excl_total = num_excl_total + 1
      else
        num_excl = enefunc%num_nonb_excl(k2) + 1
        enefunc%num_nonb_excl(k2) = num_excl
        nonb_excl_list(num_excl,k2) = k1
        num_excl_total = num_excl_total + 1
      end if
    end do

    ! exclude 1-3 interaction
    !
    do k = 1, molecule%num_angles

      k1 = molecule%angl_list(1,k)
      k3 = molecule%angl_list(3,k)

      if (nrexc(k1) < 2) then
        cycle
      end if

      if (k1 < k3) then
        num_excl = enefunc%num_nonb_excl(k1)
        duplicate = .false.
        do i = 1, enefunc%num_nonb_excl(k1)
          if (k3 == nonb_excl_list(i,k1)) duplicate = .true.
        end do
        if (.not. duplicate) then
          num_excl = num_excl + 1
          enefunc%num_nonb_excl(k1) = num_excl
          nonb_excl_list(num_excl,k1) = k3
          num_excl_total = num_excl_total + 1
        end if
      else
        num_excl = enefunc%num_nonb_excl(k3)
        duplicate = .false.
        do i = 1, enefunc%num_nonb_excl(k3)
          if (k1 == nonb_excl_list(i,k3)) duplicate = .true.
        end do
        if (.not. duplicate) then
          num_excl = num_excl + 1
          enefunc%num_nonb_excl(k3) = num_excl
          nonb_excl_list(num_excl,k3) = k1
          num_excl_total = num_excl_total + 1
        end if
      end if
    end do

    ! count 1-4 interaction
    !
    num_nb14_total = 0
    do k = 1, molecule%num_dihedrals
      k1 = molecule%dihe_list(1,k)
      k4 = molecule%dihe_list(4,k)

      if (nrexc(k1) < 3) then
        cycle
      end if

      if (k1 < k4) then
        num_nb14 = enefunc%num_nb14_calc(k1)
        duplicate = .false.
        do i = 1, enefunc%num_nonb_excl(k1)
          if (k4 == nonb_excl_list(i,k1)) duplicate = .true.
        end do
        do i = 1, num_nb14
          if (k4 == nb14_calc_list(i,k1)) duplicate = .true.
        end do
        if (.not. duplicate) then
          num_nb14 = num_nb14 + 1
          enefunc%num_nb14_calc(k1) = num_nb14
          nb14_calc_list(num_nb14,k1) = k4
          nb14_qq_scale (num_nb14,k1) = enefunc%fudge_qq
          nb14_lj_scale (num_nb14,k1) = enefunc%fudge_lj
          num_nb14_total = num_nb14_total + 1
        end if
      else
        num_nb14 = enefunc%num_nb14_calc(k4)
        duplicate = .false.
        do i = 1, enefunc%num_nonb_excl(k4)
          if (k1 == nonb_excl_list(i,k4)) duplicate = .true.
        end do
        do i = 1, num_nb14
          if (k1 == nb14_calc_list(i,k4)) duplicate = .true.
        end do
        if (.not. duplicate) then
          num_nb14 = num_nb14 + 1
          enefunc%num_nb14_calc(k4) = num_nb14
          nb14_calc_list(num_nb14,k4) = k1
          nb14_qq_scale (num_nb14,k4) = enefunc%fudge_qq
          nb14_lj_scale (num_nb14,k4) = enefunc%fudge_lj
          num_nb14_total = num_nb14_total + 1
        end if
      end if
      if (k1 < k4) then
        num_nb14 = enefunc%num_nb14_calc(k1)
        duplicate = .false.
        do i = 1, enefunc%num_nonb_excl(k1)
          if (k4 == nonb_excl_list(i,k1)) duplicate = .true.
        end do
        do i = 1, num_nb14
          if (k4 == nb14_calc_list(i,k1)) duplicate = .true.
        end do
        if (.not. duplicate) then
          num_nb14 = num_nb14 + 1
          enefunc%num_nb14_calc(k1) = num_nb14
          nb14_calc_list(num_nb14,k1) = k4
          nb14_qq_scale (num_nb14,k1) = enefunc%fudge_qq
          nb14_lj_scale (num_nb14,k1) = enefunc%fudge_lj
          num_nb14_total = num_nb14_total + 1
        end if
      else
        num_nb14 = enefunc%num_nb14_calc(k4)
        duplicate = .false.
        do i = 1, enefunc%num_nonb_excl(k4)
          if (k1 == nonb_excl_list(i,k4)) duplicate = .true.
        end do
        do i = 1, num_nb14
          if (k1 == nb14_calc_list(i,k4)) duplicate = .true.
        end do
        if (.not. duplicate) then
          num_nb14 = num_nb14 + 1
          enefunc%num_nb14_calc(k4) = num_nb14
          nb14_calc_list(num_nb14,k4) = k1
          nb14_qq_scale (num_nb14,k4) = enefunc%fudge_qq
          nb14_lj_scale (num_nb14,k4) = enefunc%fudge_lj
          num_nb14_total = num_nb14_total + 1
        end if
      end if
    end do

    ! exclude exclusions
    !
    do k = 1, size(excls(1,:))
      k1 = excls(1,k)
      k2 = excls(2,k)
      if (k1 < k2) then
        duplicate = .false.
        do i = 1, enefunc%num_nonb_excl(k1)
          if (k2 == nonb_excl_list(i,k1)) duplicate = .true.
        end do

        if (.not. duplicate) then
          num_excl = enefunc%num_nonb_excl(k1) + 1
          enefunc%num_nonb_excl(k1) = num_excl
          nonb_excl_list(num_excl,k1) = k2
          num_excl_total = num_excl_total + 1
        end if
      else
        duplicate = .false.
        do i = 1, enefunc%num_nonb_excl(k2)
          if (k1 == nonb_excl_list(i,k2)) duplicate = .true.
        end do

        if (.not. duplicate) then
          num_excl = enefunc%num_nonb_excl(k2) + 1
          enefunc%num_nonb_excl(k2) = num_excl
          nonb_excl_list(num_excl,k2) = k1
          num_excl_total = num_excl_total + 1
        end if
      end if
    end do


    ! pack 2D-array into 1D-array
    !
    allocate(enefunc%nonb_excl_list(num_excl_total), stat = alloc_stat)
      if (alloc_stat /=0) call error_msg_alloc

    call pack_array_i(molecule%num_atoms, enefunc%num_nonb_excl,  &
                      nonb_excl_list, enefunc%nonb_excl_list)
      deallocate(nonb_excl_list, stat = dealloc_stat)

    !allocate(enefunc%nb14_calc_list(num_nb14_total), &
    !         enefunc%nb14_qq_scale (num_nb14_total), &
    !         enefunc%nb14_lj_scale (num_nb14_total), stat = alloc_stat)
    max_nb14_num = max(1,maxval(enefunc%num_nb14_calc(1:molecule%num_atoms)))
    allocate(enefunc%nb14_calc_list(max_nb14_num,molecule%num_atoms), &
             enefunc%nb14_qq_scale (max_nb14_num,molecule%num_atoms), &
             enefunc%nb14_lj_scale (max_nb14_num,molecule%num_atoms), &
             stat = alloc_stat)
      if (alloc_stat /=0) call error_msg_alloc

    enefunc%nb14_calc_list(1:max_nb14_num,1:molecule%num_atoms) = 0
    enefunc%nb14_qq_scale (1:max_nb14_num,1:molecule%num_atoms) = 0.0_wp
    enefunc%nb14_lj_scale (1:max_nb14_num,1:molecule%num_atoms) = 0.0_wp
    do i = 1, molecule%num_atoms
      enefunc%nb14_calc_list(1:enefunc%num_nb14_calc(i),i) = &
             nb14_calc_list(1:enefunc%num_nb14_calc(i),i)
      enefunc%nb14_qq_scale (1:enefunc%num_nb14_calc(i),i) = &
             nb14_qq_scale (1:enefunc%num_nb14_calc(i),i)
      enefunc%nb14_lj_scale (1:enefunc%num_nb14_calc(i),i) = &
             nb14_lj_scale (1:enefunc%num_nb14_calc(i),i)
    end do
    !call pack_array_i(molecule%num_atoms, enefunc%num_nb14_calc,  &
    !                  nb14_calc_list, enefunc%nb14_calc_list)

    !call pack_array_r(molecule%num_atoms, enefunc%num_nb14_calc,  &
    !                  nb14_qq_scale, enefunc%nb14_qq_scale)

    !call pack_array_r(molecule%num_atoms, enefunc%num_nb14_calc,  &
    !                  nb14_lj_scale, enefunc%nb14_lj_scale)

    deallocate(nb14_calc_list, &
               nb14_qq_scale,  &
               nb14_lj_scale, stat = dealloc_stat)

    return

  end subroutine count_nonb_excl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    count_nonb_excl_go
  !> @brief        exclude 1-2, 1-3, 1-4 interactions
  !! @authors      NT
  !! @param[in]    molecule : molecule information
  !! @param[out]   enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine count_nonb_excl_go(molecule, enefunc, excls, nrexc)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc
    integer,                 intent(in)    :: excls(:,:)
    integer,                 intent(in)    :: nrexc(:)

    ! local variables
    integer                  :: i, k, k1, k2, k3, k4
    integer                  :: max_contact
    integer                  :: num_excl, num_excl_total
    integer                  :: num_nb14, num_nb14_total
    integer                  :: max_nb14_num
    integer                  :: alloc_stat, dealloc_stat
    logical                  :: duplicate

    integer,     allocatable :: nonb_excl_list(:,:)
    integer,     allocatable :: nb14_calc_list(:,:)
    real(wp),    allocatable :: nb14_qq_scale (:,:)
    real(wp),    allocatable :: nb14_lj_scale (:,:)
    integer,     allocatable :: temporary_list(:,:)


    ! allocate nonbonded exclusion list and 1-4 interaction list
    !
    allocate(enefunc%num_nonb_excl(molecule%num_atoms), stat=alloc_stat)
    if (alloc_stat /=0) &
      call error_msg_alloc

    enefunc%num_nonb_excl(1:molecule%num_atoms) = 0
    do k = 1, size(excls(1,:))
      k1 = excls(1,k)
      k2 = excls(2,k)
      if (k1 < k2) then
        num_excl = enefunc%num_nonb_excl(k1) + 1
        enefunc%num_nonb_excl(k1) = num_excl
      else
        num_excl = enefunc%num_nonb_excl(k2) + 1
        enefunc%num_nonb_excl(k2) = num_excl
      end if
    end do
    max_contact = 0
    do k = 1, molecule%num_atoms
      if (max_contact < enefunc%num_nonb_excl(k)) then
        max_contact = enefunc%num_nonb_excl(k)
      endif
    end do
    enefunc%num_nonb_excl(1:molecule%num_atoms) = 0

    max_contact = max_excl+max_contact
    allocate(nonb_excl_list(max_contact,molecule%num_atoms),stat=alloc_stat)
    if (alloc_stat /=0) &
      call error_msg_alloc

    allocate(enefunc%num_nb14_calc(molecule%num_atoms), stat=alloc_stat)
    if (alloc_stat /=0) &
      call error_msg_alloc

    allocate(nb14_calc_list(MAX_NB14,molecule%num_atoms), &
             nb14_qq_scale (MAX_NB14,molecule%num_atoms), &
             nb14_lj_scale (MAX_NB14,molecule%num_atoms),stat=alloc_stat)
    if (alloc_stat /=0) &
      call error_msg_alloc

    enefunc%num_nb14_calc(1:molecule%num_atoms) = 0

    ! exclude 1-2 interaction
    !
    num_excl_total = 0
    do k = 1, molecule%num_bonds
      k1 = molecule%bond_list(1,k)
      k2 = molecule%bond_list(2,k)

      if (nrexc(k1) < 1) then
        cycle
      end if

      if (k1 < k2) then
        num_excl = enefunc%num_nonb_excl(k1) + 1
        enefunc%num_nonb_excl(k1) = num_excl
        nonb_excl_list(num_excl,k1) = k2
        num_excl_total = num_excl_total + 1
      else
        num_excl = enefunc%num_nonb_excl(k2) + 1
        enefunc%num_nonb_excl(k2) = num_excl
        nonb_excl_list(num_excl,k2) = k1
        num_excl_total = num_excl_total + 1
      end if
    end do

    ! exclude 1-3 interaction
    !

    call create_nbond_list(molecule%bond_list, 3, temporary_list)

    do k = 1, size(temporary_list(1,:))

      k1 = temporary_list(1,k)
      k3 = temporary_list(3,k)

      if (nrexc(k1) < 2) then
        cycle
      end if

      if (k1 < k3) then
        num_excl = enefunc%num_nonb_excl(k1)
        duplicate = .false.
        do i = 1, enefunc%num_nonb_excl(k1)
          if (k3 == nonb_excl_list(i,k1)) duplicate = .true.
        end do
        if (.not. duplicate) then
          num_excl = num_excl + 1
          enefunc%num_nonb_excl(k1) = num_excl
          nonb_excl_list(num_excl,k1) = k3
          num_excl_total = num_excl_total + 1
        end if
      else
        num_excl = enefunc%num_nonb_excl(k3)
        duplicate = .false.
        do i = 1, enefunc%num_nonb_excl(k3)
          if (k1 == nonb_excl_list(i,k3)) duplicate = .true.
        end do
        if (.not. duplicate) then
          num_excl = num_excl + 1
          enefunc%num_nonb_excl(k3) = num_excl
          nonb_excl_list(num_excl,k3) = k1
          num_excl_total = num_excl_total + 1
        end if
      end if
    end do
    deallocate(temporary_list)

    ! count 1-4 interaction
    !
    call create_nbond_list(molecule%bond_list, 4, temporary_list)
    num_nb14_total = 0
    do k = 1, size(temporary_list(1,:))
      k1 = temporary_list(1,k)
      k4 = temporary_list(4,k)

      if (nrexc(k1) < 3) then
        cycle
      end if

      if (k1 < k4) then
        num_nb14 = enefunc%num_nb14_calc(k1)
        duplicate = .false.
        do i = 1, enefunc%num_nonb_excl(k1)
          if (k4 == nonb_excl_list(i,k1)) duplicate = .true.
        end do
        do i = 1, num_nb14
          if (k4 == nb14_calc_list(i,k1)) duplicate = .true.
        end do
        if (.not. duplicate) then
          num_nb14 = num_nb14 + 1
          enefunc%num_nb14_calc(k1) = num_nb14
          nb14_calc_list(num_nb14,k1) = k4
          nb14_qq_scale (num_nb14,k1) = enefunc%fudge_qq
          nb14_lj_scale (num_nb14,k1) = enefunc%fudge_lj
          num_nb14_total = num_nb14_total + 1
        end if
      else
        num_nb14 = enefunc%num_nb14_calc(k4)
        duplicate = .false.
        do i = 1, enefunc%num_nonb_excl(k4)
          if (k1 == nonb_excl_list(i,k4)) duplicate = .true.
        end do
        do i = 1, num_nb14
          if (k1 == nb14_calc_list(i,k4)) duplicate = .true.
        end do
        if (.not. duplicate) then
          num_nb14 = num_nb14 + 1
          enefunc%num_nb14_calc(k4) = num_nb14
          nb14_calc_list(num_nb14,k4) = k1
          nb14_qq_scale (num_nb14,k4) = enefunc%fudge_qq
          nb14_lj_scale (num_nb14,k4) = enefunc%fudge_lj
          num_nb14_total = num_nb14_total + 1
        end if
      end if
      if (k1 < k4) then
        num_nb14 = enefunc%num_nb14_calc(k1)
        duplicate = .false.
        do i = 1, enefunc%num_nonb_excl(k1)
          if (k4 == nonb_excl_list(i,k1)) duplicate = .true.
        end do
        do i = 1, num_nb14
          if (k4 == nb14_calc_list(i,k1)) duplicate = .true.
        end do
        if (.not. duplicate) then
          num_nb14 = num_nb14 + 1
          enefunc%num_nb14_calc(k1) = num_nb14
          nb14_calc_list(num_nb14,k1) = k4
          nb14_qq_scale (num_nb14,k1) = enefunc%fudge_qq
          nb14_lj_scale (num_nb14,k1) = enefunc%fudge_lj
          num_nb14_total = num_nb14_total + 1
        end if
      else
        num_nb14 = enefunc%num_nb14_calc(k4)
        duplicate = .false.
        do i = 1, enefunc%num_nonb_excl(k4)
          if (k1 == nonb_excl_list(i,k4)) duplicate = .true.
        end do
        do i = 1, num_nb14
          if (k1 == nb14_calc_list(i,k4)) duplicate = .true.
        end do
        if (.not. duplicate) then
          num_nb14 = num_nb14 + 1
          enefunc%num_nb14_calc(k4) = num_nb14
          nb14_calc_list(num_nb14,k4) = k1
          nb14_qq_scale (num_nb14,k4) = enefunc%fudge_qq
          nb14_lj_scale (num_nb14,k4) = enefunc%fudge_lj
          num_nb14_total = num_nb14_total + 1
        end if
      end if
    end do
    deallocate(temporary_list)


    ! exclude exclusions
    !
    do k = 1, size(excls(1,:))
      k1 = excls(1,k)
      k2 = excls(2,k)
      if (k1 < k2) then
        duplicate = .false.
        do i = 1, enefunc%num_nonb_excl(k1)
          if (k2 == nonb_excl_list(i,k1)) duplicate = .true.
        end do

        if (.not. duplicate) then
          num_excl = enefunc%num_nonb_excl(k1) + 1
          enefunc%num_nonb_excl(k1) = num_excl
          nonb_excl_list(num_excl,k1) = k2
          num_excl_total = num_excl_total + 1
        end if
      else
        duplicate = .false.
        do i = 1, enefunc%num_nonb_excl(k2)
          if (k1 == nonb_excl_list(i,k2)) duplicate = .true.
        end do

        if (.not. duplicate) then
          num_excl = enefunc%num_nonb_excl(k2) + 1
          enefunc%num_nonb_excl(k2) = num_excl
          nonb_excl_list(num_excl,k2) = k1
          num_excl_total = num_excl_total + 1
        end if
      end if
    end do


    ! pack 2D-array into 1D-array
    !
    allocate(enefunc%nonb_excl_list(num_excl_total), stat = alloc_stat)
      if (alloc_stat /=0) call error_msg_alloc

    call pack_array_i(molecule%num_atoms, enefunc%num_nonb_excl,  &
                      nonb_excl_list, enefunc%nonb_excl_list)
      deallocate(nonb_excl_list, stat = dealloc_stat)

    max_nb14_num = max(1,maxval(enefunc%num_nb14_calc(1:molecule%num_atoms)))
    allocate(enefunc%nb14_calc_list(max_nb14_num,molecule%num_atoms), &
             enefunc%nb14_qq_scale (max_nb14_num,molecule%num_atoms), &
             enefunc%nb14_lj_scale (max_nb14_num,molecule%num_atoms), &
             stat = alloc_stat)
      if (alloc_stat /=0) call error_msg_alloc

    do i = 1, molecule%num_atoms
      enefunc%nb14_calc_list(1:enefunc%num_nb14_calc(i),i) = &
             nb14_calc_list(1:enefunc%num_nb14_calc(i),i)
      enefunc%nb14_qq_scale (1:enefunc%num_nb14_calc(i),i) = &
             nb14_qq_scale (1:enefunc%num_nb14_calc(i),i)
      enefunc%nb14_lj_scale (1:enefunc%num_nb14_calc(i),i) = &
             nb14_lj_scale (1:enefunc%num_nb14_calc(i),i)
    end do


    deallocate(nb14_calc_list, &
               nb14_qq_scale,  &
               nb14_lj_scale, stat = dealloc_stat)

    return

  end subroutine count_nonb_excl_go

end module at_enefunc_gromacs_mod
