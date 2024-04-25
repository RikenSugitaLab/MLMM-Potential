!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_enefunc_mod
!> @brief   define potential energy functions in each domain
!! @authors Jaewoon Jung (JJ), Yuji Sugita (YS), Chigusa Kobayashi (CK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_enefunc_mod

  use sp_enefunc_gromacs_mod
  use sp_enefunc_amber_mod
  use sp_enefunc_charmm_mod
  use sp_enefunc_localres_mod
  use sp_enefunc_table_mod
  use sp_communicate_mod
  use sp_migration_mod
  use sp_energy_mod
  use sp_restraints_str_mod
  use sp_constraints_str_mod
  use sp_enefunc_str_mod
  use sp_energy_str_mod
  use sp_domain_str_mod
  use molecules_str_mod
  use fileio_localres_mod
  use fileio_grotop_mod
  use fileio_prmtop_mod
  use fileio_par_mod
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
  public  :: define_enefunc
  public  :: define_enefunc_pio
  public  :: update_enefunc
  public  :: update_enefunc_contact
  private :: setup_enefunc_bond_pio
  private :: setup_enefunc_bond_constraint_pio
  private :: setup_enefunc_angl_pio
  private :: setup_enefunc_dihe_pio
  private :: setup_enefunc_rb_dihe_pio
  private :: setup_enefunc_impr_pio
  private :: setup_enefunc_cmap_pio
  private :: setup_enefunc_nonb_pio
  private :: setup_enefunc_dispcorr
  private :: check_bonding
  ! FEP
  private :: setup_enefunc_dispcorr_fep

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_enefunc
  !> @brief        a driver subroutine for defining potential energy functions
  !! @authors      YS, JJ, CK
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[in]    par         : CHARMM PAR information
  !! @param[in]    prmtop      : AMBER parameter topology information
  !! @param[in]    grotop      : GROMACS parameter topology information
  !! @param[in]    localres    : local restraint information
  !! @param[in]    molecule    : molecule information
  !! @param[in]    constraints : constraints information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_enefunc(ene_info, par, prmtop, grotop, localres, &
                            molecule, constraints, restraints, domain, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_par),             intent(in)    :: par
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_grotop),          intent(in)    :: grotop
    type(s_localres),        intent(in)    :: localres
    type(s_molecule),        intent(inout) :: molecule
    type(s_constraints),     intent(inout) :: constraints
    type(s_restraints),      intent(in)    :: restraints
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc


    enefunc%forcefield        = ene_info%forcefield
    enefunc%output_style      = ene_info%output_style
    enefunc%table%water_table = ene_info%table .and. &
                                (ene_info%water_model(1:4) /= 'NONE' .and. &
                                (.not. ene_info%nonb_limiter))

    enefunc%switchdist        = ene_info%switchdist
    enefunc%cutoffdist        = ene_info%cutoffdist
    enefunc%pairlistdist      = ene_info%pairlistdist
    enefunc%dielec_const      = ene_info%dielec_const
    enefunc%force_switch      = ene_info%vdw_force_switch
    enefunc%vdw_shift         = ene_info%vdw_shift

    enefunc%pme_use           = ene_info%electrostatic == ElectrostaticPME
    enefunc%pme_alpha         = ene_info%pme_alpha
    enefunc%pme_ngrid_x       = ene_info%pme_ngrid_x
    enefunc%pme_ngrid_y       = ene_info%pme_ngrid_y
    enefunc%pme_ngrid_z       = ene_info%pme_ngrid_z
    enefunc%pme_nspline       = ene_info%pme_nspline
    enefunc%fft_scheme        = ene_info%fft_scheme
    enefunc%pme_max_spacing   = ene_info%pme_max_spacing
    enefunc%dispersion_corr   = ene_info%dispersion_corr
    enefunc%contact_check     = ene_info%contact_check
    enefunc%nonb_limiter      = ene_info%nonb_limiter
    enefunc%minimum_contact   = ene_info%minimum_contact
    enefunc%err_minimum_contact = ene_info%err_minimum_contact

    ! FEP
    if (domain%fep_use) &
      enefunc%fep_topology = molecule%fep_topology

    enefunc%vacuum = ene_info%vacuum

    if (ene_info%structure_check == StructureCheckDomain) then
      enefunc%pairlist_check = .true.
      enefunc%bonding_check  = .true.
    endif


    ! charmm
    !
    if (par%num_bonds > 0) then

      if (ene_info%forcefield /= ForcefieldCHARMM) &
        call error_msg('Define_Enefunc> ERROR: Bad combination between the input files and force field type (see Chapter "INPUT SECTION" in the user manual).')

      call define_enefunc_charmm(ene_info, par, localres, molecule, &
                                 constraints, restraints, domain, enefunc)

    ! amber
    !
    else if (prmtop%num_atoms > 0) then

      if (ene_info%forcefield /= ForcefieldAMBER) &
        call error_msg('Define_Enefunc> ERROR: Bad combination between the input files and force field type (see Chapter "INPUT SECTION" in the user manual).')

      call define_enefunc_amber (ene_info, prmtop, molecule, &
                                 constraints, restraints, domain, enefunc)

    ! gromacs
    !
    else if (grotop%num_atomtypes > 0) then

      if (domain%fep_use) &
        call error_msg('Define_Enefunc> FEP supports only CHARMM and AMBER forcefields.')
!      if (domain%fep_use) then
!        if (ene_info%forcefield == ForcefieldGROMARTINI .or. &
!            ene_info%forcefield == ForcefieldAAGO) &
!          call error_msg('Define_Enefunc> FEP does not support Martini and AAGo.')
!      end if

      if (ene_info%forcefield /= ForcefieldGROAMBER .and.   &
          ene_info%forcefield /= ForcefieldGROMARTINI .and. &
          ene_info%forcefield /= ForcefieldAAGO) &
        call error_msg('Define_Enefunc> ERROR: Bad combination between the input files and force field type (see Chapter "INPUT SECTION" in the user manual).')

      call define_enefunc_gromacs(ene_info, grotop, molecule, &
                                 constraints, restraints, domain, enefunc)

    end if

    ! dispersion correction
    !
    if (domain%fep_use) then
      ! FEP
      call setup_enefunc_dispcorr_fep(ene_info, domain, enefunc)
    else
      call setup_enefunc_dispcorr(ene_info, domain, enefunc)
    end if

    ! bonding_checker
    !
    if (ene_info%structure_check /= StructureCheckNone)  &
      call check_bonding(enefunc, domain)


    return

  end subroutine define_enefunc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_enefunc_pio
  !> @brief        a driver subroutine for defining potential energy functions
  !! @authors      JJ
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[in]    localres    : local restraint information
  !! @param[inout] constraints : constraints information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_enefunc_pio(ene_info, localres, constraints, &
                                domain,enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_localres),        intent(in)    :: localres
    type(s_constraints),     intent(inout) :: constraints
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: ncel, ncelb

    enefunc%forcefield        = ene_info%forcefield
    enefunc%output_style      = ene_info%output_style
    enefunc%table%water_table = ene_info%table .and. &
                                (ene_info%water_model(1:4) /= 'NONE')

    enefunc%switchdist        = ene_info%switchdist
    enefunc%cutoffdist        = ene_info%cutoffdist
    enefunc%pairlistdist      = ene_info%pairlistdist
    enefunc%dielec_const      = ene_info%dielec_const
    enefunc%force_switch      = ene_info%vdw_force_switch
    enefunc%vdw_shift         = ene_info%vdw_shift

    enefunc%pme_use           = ene_info%electrostatic == ElectrostaticPME
    enefunc%pme_alpha         = ene_info%pme_alpha
    enefunc%pme_ngrid_x       = ene_info%pme_ngrid_x
    enefunc%pme_ngrid_y       = ene_info%pme_ngrid_y
    enefunc%pme_ngrid_z       = ene_info%pme_ngrid_z
    enefunc%pme_nspline       = ene_info%pme_nspline
    enefunc%pme_max_spacing   = ene_info%pme_max_spacing
    enefunc%fft_scheme        = ene_info%fft_scheme
    enefunc%dispersion_corr   = ene_info%dispersion_corr
    enefunc%contact_check     = ene_info%contact_check
    enefunc%nonb_limiter      = ene_info%nonb_limiter
    enefunc%minimum_contact   = ene_info%minimum_contact
    enefunc%err_minimum_contact = ene_info%err_minimum_contact

    ! base
    !
    ncel  = domain%num_cell_local
    ncelb = domain%num_cell_local + domain%num_cell_boundary

    call alloc_enefunc  (enefunc, EneFuncBondCell, ncel, ncelb)

    if (.not. constraints%rigid_bond) then

      ! bond
      !
      call setup_enefunc_bond_pio(domain, enefunc)

      ! angle
      !
      call setup_enefunc_angl_pio(domain, enefunc)

    else

      ! bond
      !
      call setup_enefunc_bond_constraint_pio(domain, constraints, enefunc)

      ! angle
      !
      call setup_enefunc_angl_pio(domain, enefunc)

    end if

    ! dihedral
    !
    call setup_enefunc_dihe_pio(domain, enefunc)

    ! Ryckaert-Bellemans dihedral
    !
    call setup_enefunc_rb_dihe_pio(domain, enefunc)

    ! improper
    !
    call setup_enefunc_impr_pio(domain, enefunc)

    ! cmap
    !
    call setup_enefunc_cmap_pio(ene_info, domain, enefunc)

    ! nonbonded
    !
    call setup_enefunc_nonb_pio(ene_info, constraints, domain, enefunc)

    ! dispersion correction
    !
    call setup_enefunc_dispcorr(ene_info, domain, enefunc)

    ! lookup table
    !
    if(ene_info%table) &
    call setup_enefunc_table(ene_info, enefunc)

    ! restraint
    !
    call setup_enefunc_localres(localres, domain, enefunc)

    if (ene_info%structure_check /= StructureCheckNone)  &
      call check_bonding(enefunc, domain)

    ! write summary of energy function
    !
    if (main_rank) then
      write(MsgOut,'(A)') &
           'Define_Enefunc_Pio> Number of Interactions in Each Term'
      write(MsgOut,'(A20,I10,A20,I10)')                  &
           '  bond_ene        = ', enefunc%num_bond_all, &
           '  angle_ene       = ', enefunc%num_angl_all
      write(MsgOut,'(A20,I10,A20,I10)')                  &
           '  torsion_ene     = ', enefunc%num_dihe_all, &
           '  rb_torsion_ene  = ', enefunc%num_rb_dihe_all
      write(MsgOut,'(A20,I10,A20,I10)')                  &
           '  improper_ene    = ', enefunc%num_impr_all, &
           '  cmap_ene        = ', enefunc%num_cmap_all
      write(MsgOut,'(A20,I10,A20,I10)')                  &
           '  nb_exclusions   = ', enefunc%num_excl_all, &
           '  nb14_calc       = ', enefunc%num_nb14_all
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine define_enefunc_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_enefunc
  !> @brief        a driver subroutine for updating potential energy functions
  !! @authors      JJ
  !! @param[in]    table       : flag for table or not
  !! @param[inout] domain      : domain information
  !! @param[inout] comm        : communication information
  !! @param[inout] enefunc     : energy potential functions information
  !! @param[inout] constraints : constraints information [optional]
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_enefunc(table, domain, comm, enefunc, constraints)

    ! formal arguments
    logical,                       intent(in)    :: table
    type(s_domain),                intent(inout) :: domain
    type(s_comm),                  intent(inout) :: comm
    type(s_enefunc),               intent(inout) :: enefunc
    type(s_constraints), optional, intent(inout) :: constraints

    ! local variables
    logical                        :: first


    ! sending the bonding information to other domain
    !

    ! bond
    !
    call update_outgoing_enefunc_bond(domain, enefunc)

    ! angle
    !
    call update_outgoing_enefunc_angl(domain, enefunc)

    ! dihedral
    !
    call update_outgoing_enefunc_dihe(domain, enefunc)

    ! Ryckaert-Bellemans dihedral
    !
    call update_outgoing_enefunc_rb_dihe(domain, enefunc)

    ! improper dihedral
    !
    call update_outgoing_enefunc_impr(domain, enefunc)

    ! cmap
    !
    call update_outgoing_enefunc_cmap(domain, enefunc)

    ! restraint
    !
    call update_outgoing_enefunc_restraint(domain, enefunc)

    ! fitting
    !
    call update_outgoing_enefunc_fitting(domain, enefunc)

    ! communicate neighbour domain
    !
    call communicate_bond(domain, comm, enefunc)


    ! bond
    !
    call update_incoming_enefunc_bond(domain, enefunc)

    ! angle
    !
    call update_incoming_enefunc_angl(domain, enefunc)

    ! dihedral
    !
    call update_incoming_enefunc_dihe(domain, enefunc)

    ! Ryckaert-Bellemans dihedral
    !
    call update_incoming_enefunc_rb_dihe(domain, enefunc)

    ! improper dihedral
    !
    call update_incoming_enefunc_impr(domain, enefunc)

    ! cmap
    !
    call update_incoming_enefunc_cmap(domain, enefunc)

    ! restraint
    !
    call update_incoming_enefunc_restraint(domain, enefunc)

    ! fitting
    !
    call update_incoming_enefunc_fitting(domain, enefunc)

    ! re-count nonbond exclusion list
    !

    first = .false.

    if (present(constraints)) then

      if (enefunc%local_restraint) then
        call count_nonb_excl_constraint_rest(first, table, constraints, &
                                             domain, enefunc)

      else
        call count_nonb_excl_constraint(first, table, constraints, &
                                             domain, enefunc)

      end if

    else

      if (enefunc%local_restraint) then
        if (table) then
          call count_nonb_excl_solute_rest(first, domain, enefunc)

        else
          call count_nonb_excl_rest(first, domain, enefunc)

        end if
      else
        if (table) then
          if (domain%fep_use) then
            call error_msg('Update_Enefunc> table is not available without constraints in FEP.')
          else
            call count_nonb_excl_solute(first, domain, enefunc)
          end if
        else
          call count_nonb_excl(first, domain, enefunc)

        end if
      end if

    end if

    if (enefunc%bonding_check) call check_bonding(enefunc, domain)

    return

  end subroutine update_enefunc

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    update_enefunc_contact
  !> @brief        a driver subroutine for updating potential energy functions
  !! @authors      JJ
  !! @param[in]    table       : flag for table or not
  !! @param[inout] domain      : domain information
  !! @param[inout] comm        : communication information
  !! @param[inout] enefunc     : energy potential functions information
  !! @param[inout] constraints : constraints information [optional]
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine update_enefunc_contact(domain, comm, enefunc)

    ! formal arguments
    type(s_domain),                intent(inout) :: domain
    type(s_comm),                  intent(inout) :: comm
    type(s_enefunc),               intent(inout) :: enefunc

    ! local variables
    logical                        :: first
    integer                        :: ncell

    ! sending the bonding information to other domain
    !

    ! bond
    !
    call update_outgoing_enefunc_bond(domain, enefunc)

    ! contact
    !
    call update_outgoing_enefunc_contact(domain, enefunc)

    ! angle
    !
    call update_outgoing_enefunc_angl(domain, enefunc)

    ! dihedral
    !
    call update_outgoing_enefunc_dihe(domain, enefunc)

    ! Ryckaert-Bellemans dihedral
    !
    call update_outgoing_enefunc_rb_dihe(domain, enefunc)

    ! improper dihedral
    !
    call update_outgoing_enefunc_impr(domain, enefunc)

    ! restraint
    !
    call update_outgoing_enefunc_restraint(domain, enefunc)

    ! communicate neighbour domain
    !
    call communicate_contact(domain, comm, enefunc)

    ! bond
    !
    call update_incoming_enefunc_bond(domain, enefunc)

    ! contact
    !
    call update_incoming_enefunc_contact(domain, enefunc)

    ! angle
    !
    call update_incoming_enefunc_angl(domain, enefunc)

    ! dihedral
    !
    call update_incoming_enefunc_dihe(domain, enefunc)

    ! Ryckaert-Bellemans dihedral
    !
    call update_incoming_enefunc_rb_dihe(domain, enefunc)

    ! improper dihedral
    !
    call update_incoming_enefunc_impr(domain, enefunc)

    ! restraint
    !
    call update_incoming_enefunc_restraint(domain, enefunc)

    ! re-count nonbond exclusion list
    !
    first = .false.
    call count_nonb_excl_go(first, domain, enefunc)

    return

  end subroutine update_enefunc_contact

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_bond_pio
  !> @brief        define BOND term for each cell in potential energy function
  !! @authors      NT
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_bond_pio(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variable
    integer                  :: i, found
    integer,         pointer :: bond(:)


    bond      => enefunc%num_bond

    found = 0
    do i = 1, domain%num_cell_local
      found = found + bond(i)
      if (bond(i) > MaxBond) &
        call error_msg('Setup_Enefunc_Bond_Pio> Too many bonds.')
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_bond_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_bond_all = found
#endif

    return

  end subroutine setup_enefunc_bond_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_bond_constraint_pio
  !> @brief        define BOND term for each cell in potential energy function
  !! @authors      NT
  !! @param[in]    domain      : domain information
  !! @param[in]    constraints : constraints information
  !! @param[inout] enefunc     : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_bond_constraint_pio(domain, constraints, enefunc)

    ! formal arguments
    type(s_domain),      target, intent(in)    :: domain
    type(s_constraints), target, intent(inout) :: constraints
    type(s_enefunc),     target, intent(inout) :: enefunc

    ! local variable
    integer                  :: i, j, k, ih, found, ncel, connect
    integer,         pointer :: bond(:)
    integer,         pointer :: HGr_local(:,:)


    ! count # of bonds
    bond      => enefunc%num_bond

    found = 0
    do i = 1, domain%num_cell_local
      found = found + bond(i)
      if (bond(i) > MaxBond) &
        call error_msg('Setup_Enefunc_Bond_Constraint_Pio> Too many bonds.')
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_bond_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_bond_all = found
#endif


    ! count # of constraints bonds
    HGr_local => constraints%HGr_local

    ncel    = domain%num_cell_local
    connect = constraints%connect

    found = 0

    do i = 1, ncel
      do j = 1, connect
        do k = 1, HGr_local(j,i)
          do ih = 1, j
            found = found + 1
          end do
        end do
      end do
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, constraints%num_bonds, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    constraints%num_bonds = found
#endif


    return

  end subroutine setup_enefunc_bond_constraint_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_angl_pio
  !> @brief        define ANGLE term for each cell in potential energy function
  !! @authors      NT
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_angl_pio(domain, enefunc)

    ! formal arguments
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: i, found
    integer,         pointer :: angle(:)


    angle     => enefunc%num_angle

    found = 0
    do i = 1, domain%num_cell_local
      found = found + angle(i)
      if (angle(i) > MaxAngle) &
        call error_msg('Setup_Enefunc_Angl_Pio> Too many angles.')
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_angl_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_angl_all = found
#endif

    return

  end subroutine setup_enefunc_angl_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_dihe_pio
  !> @brief        define DIHEDRAL term in potential energy function
  !! @authors      NT
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_dihe_pio(domain, enefunc)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(inout) :: enefunc

    ! local variables
    integer                   :: i, found
    integer,          pointer :: dihedral(:)


    dihedral  => enefunc%num_dihedral

    found = 0
    do i = 1, domain%num_cell_local
      found = found + dihedral(i)
      if (dihedral(i) > MaxDihe) &
        call error_msg('Setup_Enefunc_Dihe_Pio> Too many dihedral angles.')
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_dihe_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_dihe_all = found
#endif

    return

  end subroutine setup_enefunc_dihe_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_rb_dihe_pio
  !> @brief        define Ryckaert-Bellemans DIHEDRAL term in potential energy
  !!               function
  !! @authors      NT
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_rb_dihe_pio(domain, enefunc)

    ! formal arguments
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(inout) :: enefunc

    ! local variables
    integer                   :: i, found
    integer,          pointer :: dihedral(:)


    dihedral  => enefunc%num_rb_dihedral

    found = 0
    do i = 1, domain%num_cell_local
      found = found + dihedral(i)
      if (dihedral(i) > MaxDihe) &
        call error_msg('Setup_Enefunc_RB_Dihe_Pio> Too many dihedral angles.')
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_rb_dihe_all, 1, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#else
    enefunc%num_rb_dihe_all = found
#endif

    return

  end subroutine setup_enefunc_rb_dihe_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_impr_pio
  !> @brief        define IMPROPER term in potential energy function
  !! @authors      NT
  !! @param[in]    domain  : domain information
  !! @param[inout] enefunc : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_impr_pio(domain, enefunc)

    ! formal variables
    type(s_domain),   target, intent(in)    :: domain
    type(s_enefunc),  target, intent(inout) :: enefunc

    ! local variables
    integer                   :: i, found
    integer,          pointer :: improper(:)


    improper  => enefunc%num_improper

    found = 0
    do i = 1, domain%num_cell_local
      found = found + improper(i)
      if (improper(i) > MaxImpr) &
        call error_msg( &
            'Setup_Enefunc_Impr_Pio> Too many improper dihedral angles')
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_impr_all, 1, mpi_integer, mpi_sum, &
                       mpi_comm_country, ierror)
#else
    enefunc%num_impr_all = found
#endif

    return

  end subroutine setup_enefunc_impr_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_cmap_pio
  !> @brief        define cmap term in potential energy function with DD
  !! @authors      NT
  !! @param[in]    ene_info : ENERGY section control parameters information
  !! @param[in]    domain   : domain information
  !! @param[inout] enefunc  : energy potential functions informationn
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_cmap_pio(ene_info, domain, enefunc)

    ! formal variables
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_domain),  target, intent(in)    :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: i, found


    found = 0
    do i = 1, domain%num_cell_local
      found = found + enefunc%num_cmap(i)
      if (enefunc%num_cmap(i) > MaxCmap) &
        call error_msg('Setup_Enefunc_Cmap_Pio> Too many cmaps.')
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(found, enefunc%num_cmap_all, 1, mpi_integer, mpi_sum, &
                       mpi_comm_country, ierror)
#else
    enefunc%num_cmap_all = found
#endif

    return

  end subroutine setup_enefunc_cmap_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_nonb_pio
  !> @brief        define NON-BOND term in potential energy function
  !! @authors      NT
  !! @param[in]    ene_info    : ENERGY section control parameters information
  !! @param[in]    molecule    : molecule information
  !! @param[in]    constraints : constraints information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_nonb_pio(ene_info, constraints, domain, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_constraints),     intent(in)    :: constraints
    type(s_domain),          intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: ncel


    ! treatment for 1-2, 1-3, 1-4 interactions
    !
    ncel   = domain%num_cell_local

    call alloc_enefunc(enefunc, EneFuncNonb,     ncel, maxcell)
    call alloc_enefunc(enefunc, EneFuncNonbList, ncel, maxcell)

    if (enefunc%table%water_table) then

      if (constraints%rigid_bond) then
        call count_nonb_excl_constraint(.true., .true., constraints, &
                                        domain, enefunc)

      else
        call count_nonb_excl_solute(.true., domain, enefunc)

      end if

    else

      if (constraints%rigid_bond) then
        call count_nonb_excl_constraint(.true., .false., constraints, &
                                        domain, enefunc)

      else
        call count_nonb_excl(.true., domain, enefunc)

      end if

    end if

    return

  end subroutine setup_enefunc_nonb_pio

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_dispcorr
  !> @brief        define dispersion correction term
  !! @authors      CK
  !! @param[in]    ene_info : ENERGY section control parameters information
  !! @param[inout] domain   : domain information
  !! @param[inout] enefunc  : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  subroutine setup_enefunc_dispcorr(ene_info, domain, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_domain),  target, intent(inout) :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, iatmcls, ntypes
    integer                  :: icel1, icel2, icel3, icel4
    integer                  :: i1, i2, i3, i4
    integer                  :: num_all_atoms, natom2, nexpair
    real(dp)                 :: lj6_tot, lj6_diff, lj6_ex
    real(dp)                 :: factor, rpair
    real(dp)                 :: diff_cs, diff_cs2, diff_cs3, diff_cs4
    real(dp)                 :: cutoff , cutoff2, cutoff3, cutoff4
    real(dp)                 :: cutoff5, cutoff6, cutoff7, cutoff8
    real(dp)                 :: cutoff14
    real(dp)                 :: inv_cutoff3, inv_cutoff6, inv_cutoff12
    real(dp)                 :: switchdist , switchdist2, switchdist3
    real(dp)                 :: switchdist4, switchdist5
    real(dp)                 :: switchdist6, switchdist7, switchdist8
    real(dp)                 :: shift_a, shift_b, shift_c
    real(dp)                 :: vswitch, eswitch, vlong

    integer,         pointer :: id_g2l(:,:)
    integer,         pointer :: bondlist(:,:,:),anglelist(:,:,:)
    integer,         pointer :: dihelist(:,:,:),rb_dihelist(:,:,:)
    integer,         pointer :: atmcls(:,:),imprlist(:,:,:)
    integer,     allocatable :: atype(:)


    if (ene_info%dispersion_corr == Disp_corr_NONE) return

    bondlist    => enefunc%bond_list
    anglelist   => enefunc%angle_list
    dihelist    => enefunc%dihe_list
    rb_dihelist => enefunc%rb_dihe_list
    imprlist    => enefunc%impr_list
    atmcls      => domain%atom_cls_no
    id_g2l      => domain%id_g2l

    ntypes = enefunc%num_atom_cls
    allocate(atype(1:ntypes))

    atype(1:ntypes) = 0
    num_all_atoms   = 0

    do i = 1, domain%num_cell_local
      do j = 1, domain%num_atom(i)
        iatmcls = atmcls(j,i)
        atype(iatmcls) = atype(iatmcls)+1
      end do
      num_all_atoms = num_all_atoms + domain%num_atom(i)
    end do

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, atype, ntypes, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#endif

    lj6_tot = 0.0_dp
    do i = 1, ntypes
      do j = 1, ntypes
        lj6_tot = lj6_tot + enefunc%nonb_lj6(i,j)*atype(i)*atype(j)
      end do
    end do
    deallocate(atype)

    cutoff       = enefunc%cutoffdist
    cutoff2      = cutoff*cutoff
    cutoff3      = cutoff2*cutoff
    inv_cutoff3  = 1.0_dp/cutoff3

    eswitch = 0.0_dp
    vswitch = 0.0_dp
    vlong   = inv_cutoff3/3.0_dp

    if (enefunc%forcefield == ForcefieldAMBER ) then

      factor       = 2.0_dp*PI*lj6_tot
      enefunc%dispersion_energy = -factor*vlong
      enefunc%dispersion_virial = -2.0_dp*factor*vlong

    else if (enefunc%forcefield == ForcefieldGROAMBER .or.  &
             enefunc%forcefield == ForcefieldGROMARTINI) then
      !
      ! remove exclusion
      !
      lj6_ex = 0.0_dp
      nexpair = 0
      do i = 1, domain%num_cell_local
        ! self
        do j = 1, domain%num_atom(i)
          iatmcls = atmcls(j,i)
          lj6_ex  = lj6_ex + enefunc%nonb_lj6(iatmcls,iatmcls)
        end do

        ! bonds
        do j = 1, enefunc%num_bond(i)
          icel1 = id_g2l(1,bondlist(1,j,i))
          i1    = id_g2l(2,bondlist(1,j,i))
          icel2 = id_g2l(1,bondlist(2,j,i))
          i2    = id_g2l(2,bondlist(2,j,i))
          lj6_ex= lj6_ex + enefunc%nb14_lj6(atmcls(i1,icel1),atmcls(i2,icel2))
        end do

        ! angles
        do j = 1, enefunc%num_angle(i)
          icel1 = id_g2l(1,anglelist(1,j,i))
          i1    = id_g2l(2,anglelist(1,j,i))
          icel3 = id_g2l(1,anglelist(3,j,i))
          i3    = id_g2l(2,anglelist(3,j,i))
          lj6_ex= lj6_ex + enefunc%nb14_lj6(atmcls(i1,icel1),atmcls(i3,icel3))
        end do

        ! dihedral
        do j = 1, enefunc%num_dihedral(i)
          icel1 = id_g2l(1,dihelist(1,j,i))
          i1    = id_g2l(2,dihelist(1,j,i))
          icel4 = id_g2l(1,dihelist(4,j,i))
          i4    = id_g2l(2,dihelist(4,j,i))
          lj6_ex= lj6_ex + enefunc%nb14_lj6(atmcls(i1,icel1),atmcls(i4,icel4))
        end do

        ! RB dihedral
        do j = 1, enefunc%num_rb_dihedral(i)
          icel1 = id_g2l(1,rb_dihelist(1,j,i))
          i1    = id_g2l(2,rb_dihelist(1,j,i))
          icel4 = id_g2l(1,rb_dihelist(4,j,i))
          i4    = id_g2l(2,rb_dihelist(4,j,i))
          lj6_ex= lj6_ex + enefunc%nb14_lj6(atmcls(i1,icel1),atmcls(i4,icel4))
        end do

        ! improper
        do j = 1, enefunc%num_improper(i)
          icel1 = id_g2l(1,imprlist(1,j,i))
          i1    = id_g2l(2,imprlist(1,j,i))
          icel4 = id_g2l(1,imprlist(4,j,i))
          i4    = id_g2l(2,imprlist(4,j,i))
          lj6_ex= lj6_ex + enefunc%nb14_lj6(atmcls(i1,icel1),atmcls(i4,icel4))
        end do

        nexpair = nexpair + domain%num_atom(i)        &
                          + enefunc%num_bond(i)        &
                          + enefunc%num_angle(i)       &
                          + enefunc%num_dihedral(i)    &
                          + enefunc%num_rb_dihedral(i) &
                          + enefunc%num_improper(i)
      end do
#ifdef HAVE_MPI_GENESIS
      call mpi_allreduce(mpi_in_place, num_all_atoms, 1, mpi_integer, &
                         mpi_sum, mpi_comm_country, ierror)
      call mpi_allreduce(mpi_in_place, nexpair, 1, mpi_integer, &
                         mpi_sum, mpi_comm_country, ierror)
      call mpi_allreduce(mpi_in_place, lj6_ex, 1, mpi_real8, &
                         mpi_sum, mpi_comm_country, ierror)
#endif
      lj6_diff = (lj6_tot - lj6_ex)

      natom2 = num_all_atoms*num_all_atoms
      rpair  = real(natom2/(natom2-nexpair),dp)
      factor       = 2.0_dp*PI*rpair*lj6_diff

      switchdist   = enefunc%switchdist
      diff_cs      = (cutoff - switchdist)

      if (diff_cs > EPS) then

        if (enefunc%vdw_shift) then
          cutoff4      = cutoff3*cutoff
          cutoff5      = cutoff4*cutoff
          cutoff6      = cutoff5*cutoff
          cutoff7      = cutoff6*cutoff
          cutoff8      = cutoff7*cutoff
          cutoff14     = cutoff7*cutoff7
          inv_cutoff6  = inv_cutoff3*inv_cutoff3
          inv_cutoff12 = inv_cutoff6*inv_cutoff6
  
          diff_cs2     = diff_cs*diff_cs
          diff_cs3     = diff_cs2*diff_cs
          diff_cs4     = diff_cs3*diff_cs
  
          switchdist2  = switchdist*switchdist
          switchdist3  = switchdist2*switchdist
          switchdist4  = switchdist3*switchdist
          switchdist5  = switchdist4*switchdist
          switchdist6  = switchdist5*switchdist
          switchdist7  = switchdist6*switchdist
          switchdist8  = switchdist7*switchdist
  
          ! LJ6
          !
          shift_a = -(10.0_dp*cutoff - 7.0_dp*switchdist)/(cutoff8*diff_cs2)
          shift_b =  ( 9.0_dp*cutoff - 7.0_dp*switchdist)/(cutoff8*diff_cs3)
  
          shift_c = inv_cutoff6 - 2.0_dp * shift_a * diff_cs3  &
                    - 1.5_dp * shift_b * diff_cs4
  
          eswitch = -2.0_dp * shift_a * ((1.0_dp/6.0_dp)*cutoff6                &
                                        -(3.0_dp/5.0_dp)*cutoff5*switchdist     &
                                        +(3.0_dp/4.0_dp)*cutoff4*switchdist2    &
                                        -(1.0_dp/3.0_dp)*cutoff3*switchdist3    &
                                        +(1.0_dp/6.0e1_dp)*switchdist6)         &
                    -1.5_dp * shift_b * ((1.0_dp/7.0_dp)*cutoff7                &
                                        -(2.0_dp/3.0_dp)*cutoff6*switchdist     &
                                        +(6.0_dp/5.0_dp)*cutoff5*switchdist2    &
                                        -                cutoff4*switchdist3    &
                                        +(1.0_dp/3.0_dp)*cutoff3*switchdist4    &
                                        -(1.0_dp/1.05e2_dp)*switchdist7)        &
                    -(1.0_dp/3.0_dp) * shift_c * (cutoff3)
    
          ! LJ12
          !
          shift_a = -(16.0_dp*cutoff - 13.0_dp*switchdist)/(cutoff14*diff_cs2)
          shift_b =  (15.0_dp*cutoff - 13.0_dp*switchdist)/(cutoff14*diff_cs3)
          shift_c = inv_cutoff12 - 2.0_dp * shift_a * diff_cs3  &
                    - 1.5_dp * shift_b * diff_cs4
  
 
          shift_a = -(10.0_dp*cutoff - 7.0_dp*switchdist)/(cutoff8*diff_cs2)
          shift_b =  ( 9.0_dp*cutoff - 7.0_dp*switchdist)/(cutoff8*diff_cs3)
 
          vswitch = shift_a * ( (1.0_dp/6.0_dp)*cutoff6                         &
                               -(2.0_dp/5.0_dp)*cutoff5*switchdist              &
                               +(1.0_dp/4.0_dp)*cutoff4*switchdist2             &
                               -(1.0_dp/6.0e1_dp)*switchdist6)                  &
                   +shift_b * ( (1.0_dp/7.0_dp)*cutoff7                         &
                               -(1.0_dp/2.0_dp)*cutoff6*switchdist              &
                               +(3.0_dp/5.0_dp)*cutoff5*switchdist2             &
                               -(1.0_dp/4.0_dp)*cutoff4*switchdist3             &
                               +(1.0_dp/1.4e2_dp)*switchdist7)
        enefunc%dispersion_energy = factor*(eswitch-vlong)
        enefunc%dispersion_virial = -2.0_dp*factor*(-vswitch+vlong)

        else

          eswitch = enefunc%eswitch
          vswitch = enefunc%vswitch
          enefunc%dispersion_energy = factor*(eswitch-vlong)
          enefunc%dispersion_virial = -factor*(vswitch+vlong)

        end if

      else 

        enefunc%dispersion_energy = factor*(eswitch-vlong)
        enefunc%dispersion_virial = -2.0_dp*factor*(-vswitch+vlong)

      end if

    else
      call error_msg('Setup_Enefunc_DispCorr> This force field is not allowed')
    end if

!   enefunc%dispersion_energy = factor*(eswitch-vlong)
!   enefunc%dispersion_virial = -2.0_wp*factor*(-vswitch+vlong)

  end subroutine setup_enefunc_dispcorr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    check_bonding
  !> @brief        check bonds
  !! @authors      CK
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    domain   : domain information
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_bonding(enefunc, domain)

    ! formal arguments
    type(s_enefunc), target, intent(in)    :: enefunc
    type(s_domain),  target, intent(in)    :: domain

    ! local variables
    real(wp)                 :: d12(1:3), r12, r_dif
    integer                  :: i, j, ix, icel1, icel2, i1, i2
    integer                  :: icel3, i3,icel4, i4
    integer                  :: id, my_id, omp_get_thread_num
    real(wp), parameter      :: maxdistance = 0.5_wp
    real(wp)                 :: maxcell_size 

    real(wp),        pointer :: r0(:,:)
    integer,         pointer :: nbond(:), bondlist(:,:,:)
    integer,         pointer :: nangle(:), anglelist(:,:,:)
    integer,         pointer :: ndihe(:),  dihelist(:,:,:)
    integer,         pointer :: nrbdihe(:),  rb_dihelist(:,:,:)
    integer,         pointer :: nimpr(:),  imprlist(:,:,:)
    integer,         pointer :: ncell_local
    integer,         pointer :: id_g2l(:,:)
    real(dp),        pointer :: coord(:,:,:)

    ncell_local => domain%num_cell_local
    id_g2l      => domain%id_g2l
    coord       => domain%coord

    maxcell_size = max(domain%cell_size(1),  &
                       domain%cell_size(2),  &
                       domain%cell_size(3))

    nbond       => enefunc%num_bond
    bondlist    => enefunc%bond_list
    r0          => enefunc%bond_dist_min

    nangle      => enefunc%num_angle
    anglelist   => enefunc%angle_list

    ndihe       => enefunc%num_dihedral
    dihelist    => enefunc%dihe_list

    nrbdihe     => enefunc%num_rb_dihedral
    rb_dihelist => enefunc%rb_dihe_list

    nimpr       => enefunc%num_improper
    imprlist    => enefunc%impr_list

    !$omp parallel default(shared)                                     &
    !$omp private(id, i, j, ix, icel1, i1, icel2, i2, d12, r12, r_dif, &
    !$omp         my_id, icel3, i3, icel4, i4)
    !
#ifdef OMP
    id  = omp_get_thread_num()
#else
    id  = 0
#endif
    my_id = id

    do i = my_id+1, ncell_local, nthread

      do ix = 1, nbond(i)

        icel1 = id_g2l(1,bondlist(1,ix,i))
        i1    = id_g2l(2,bondlist(1,ix,i))
        icel2 = id_g2l(1,bondlist(2,ix,i))
        i2    = id_g2l(2,bondlist(2,ix,i))

        d12(1:3) = coord(1:3,i1,icel1) - coord(1:3,i2,icel2)
        r12   = sqrt( d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3) )
        r_dif = r12 - r0(ix,i)
        if (r_dif > maxdistance) &
           write(MsgOut,'(A,I10,I10,F10.5)') &
          'WARNING: too long bond:',bondlist(1,ix,i),bondlist(2,ix,i),r12
        if (r_dif < -maxdistance) &
           write(MsgOut,'(A,I10,I10,F10.5)') &
          'WARNING: too short bond:',bondlist(1,ix,i),bondlist(2,ix,i),r12
        if (r12 > maxcell_size) then
           write(MsgOut,'(A,2I10,F10.5)') &
          'Check_bonding> distance is grater than cellsize:', &
           bondlist(1,ix,i),bondlist(2,ix,i),r12
           call error_msg('')
        endif

      end do

      do ix = 1, nangle(i)

        icel1 = id_g2l(1,anglelist(1,ix,i))
        i1    = id_g2l(2,anglelist(1,ix,i))
        icel3 = id_g2l(1,anglelist(3,ix,i))
        i3    = id_g2l(2,anglelist(3,ix,i))

        d12(1:3) = coord(1:3,i1,icel1) - coord(1:3,i3,icel3)
        r12   = sqrt( d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3) )

        if (r12 > maxcell_size) then
           write(MsgOut,'(A,2I10,F10.5)') &
           'Check_bonding> distance in angle is grater than cellsize:', &
           anglelist(1,ix,i),anglelist(3,ix,i),r12
           call error_msg('')
        endif

      end do

      do ix = 1, ndihe(i)

        icel1 = id_g2l(1,dihelist(1,ix,i))
        i1    = id_g2l(2,dihelist(1,ix,i))
        icel4 = id_g2l(1,dihelist(4,ix,i))
        i4    = id_g2l(2,dihelist(4,ix,i))

        d12(1:3) = coord(1:3,i1,icel1) - coord(1:3,i4,icel4)
        r12   = sqrt( d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3) )

        if (r12 > maxcell_size) then
           write(MsgOut,'(A,2I10,F10.5)') &
           'Check_bonding> distance in dihedral is grater than cellsize:', &
           dihelist(1,ix,i),dihelist(4,ix,i),r12
           call error_msg('')
        endif

      end do

      if (nrbdihe(i) > 0) then
        do ix = 1, nrbdihe(i)
       
          icel1 = id_g2l(1,rb_dihelist(1,ix,i))
          i1    = id_g2l(2,rb_dihelist(1,ix,i))
          icel4 = id_g2l(1,rb_dihelist(4,ix,i))
          i4    = id_g2l(2,rb_dihelist(4,ix,i))
       
          d12(1:3) = coord(1:3,i1,icel1) - coord(1:3,i4,icel4)
          r12   = sqrt( d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3) )
          if (r12 > maxcell_size) then
            write(MsgOut,'(A,2I10,F10.5)') &
           'Check_bonding> distance in rb dihedral is grater than cellsize:', &
            rb_dihelist(1,ix,i), rb_dihelist(4,ix,i),r12
            call error_msg('')
          endif
       
        end do
      endif

      do ix = 1, nimpr(i)

        icel1 = id_g2l(1,imprlist(1,ix,i))
        i1    = id_g2l(2,imprlist(1,ix,i))
        icel4 = id_g2l(1,imprlist(4,ix,i))
        i4    = id_g2l(2,imprlist(4,ix,i))

        d12(1:3) = coord(1:3,i1,icel1) - coord(1:3,i4,icel4)
        r12   = sqrt( d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3) )

        if (r12 > maxcell_size) then
          write(MsgOut,'(A,2I10,F10.5)') &
      'Check_bonding> distance in improper dihedral is grater than cellsize:', &
          imprlist(1,ix,i), imprlist(4,ix,i),r12
          call error_msg('')
        endif
      end do

    end do

    !$omp end parallel 

    return

  end subroutine check_bonding

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_dispcorr_fep
  !> @brief        define dispersion correction term for FEP
  !! @authors      HO
  !! @param[in]    ene_info : ENERGY section control parameters information
  !! @param[inout] domain   : domain information
  !! @param[inout] enefunc  : energy potential functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_dispcorr_fep(ene_info, domain, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_domain),  target, intent(inout) :: domain
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, iatmcls, ntypes
    integer                  :: icel1, icel2, icel3, icel4
    integer                  :: i1, i2, i3, i4
!    integer                  :: num_all_atoms, natom2, nexpair
!    real(dp)                 :: lj6_tot, lj6_diff, lj6_ex
!    real(dp)                 :: factor, rpair
    real(dp)                 :: diff_cs, diff_cs2, diff_cs3, diff_cs4
    real(dp)                 :: cutoff , cutoff2, cutoff3, cutoff4
    real(dp)                 :: cutoff5, cutoff6, cutoff7, cutoff8
    real(dp)                 :: cutoff14
    real(dp)                 :: inv_cutoff3, inv_cutoff6, inv_cutoff12
    real(dp)                 :: switchdist , switchdist2, switchdist3
    real(dp)                 :: switchdist4, switchdist5
    real(dp)                 :: switchdist6, switchdist7, switchdist8
    real(dp)                 :: shift_a, shift_b, shift_c
    real(dp)                 :: vswitch, eswitch, vlong

    integer,         pointer :: id_g2l(:,:)
    integer,         pointer :: bondlist(:,:,:),anglelist(:,:,:)
    integer,         pointer :: dihelist(:,:,:),rb_dihelist(:,:,:)
    integer,         pointer :: atmcls(:,:),imprlist(:,:,:)
!    integer,     allocatable :: atype(:)

    ! FEP
    integer                  :: num_all_atoms_preserve, natom_preserve2
    integer                  :: num_all_atoms_vanish, natom_vanish2
    integer                  :: num_all_atoms_appear, natom_appear2
    integer                  :: nexpair_preserve, nexpair_vanish, nexpair_appear
    integer                  :: fg1, fg2, pert_flag
    integer                  :: table_pert(5,5)
    real(dp)                 :: lj6_tot_preserve, lj6_tot_vanish, lj6_tot_appear
    real(dp)                 :: lj6_diff_preserve, lj6_diff_vanish, lj6_diff_appear
    real(dp)                 :: lj6_ex_preserve, lj6_ex_vanish, lj6_ex_appear
    real(dp)                 :: factor_preserve, factor_vanish, factor_appear
    real(dp)                 :: rpair_preserve, rpair_vanish, rpair_appear
    integer,     allocatable :: atype_preserve(:), atype_vanish(:), atype_appear(:)

    if (ene_info%dispersion_corr == Disp_corr_NONE) return

    bondlist    => enefunc%bond_list
    anglelist   => enefunc%angle_list
    dihelist    => enefunc%dihe_list
    rb_dihelist => enefunc%rb_dihe_list
    imprlist    => enefunc%impr_list
    atmcls      => domain%atom_cls_no
    id_g2l      => domain%id_g2l

    ntypes = enefunc%num_atom_cls
!    allocate(atype(1:ntypes))
    allocate(atype_preserve(1:ntypes), atype_vanish(1:ntypes), atype_appear(1:ntypes))

!    atype(1:ntypes) = 0
    atype_preserve(1:ntypes) = 0
    atype_vanish(1:ntypes)   = 0
    atype_appear(1:ntypes)   = 0
!    num_all_atoms            = 0
    num_all_atoms_preserve   = 0
    num_all_atoms_vanish     = 0
    num_all_atoms_appear     = 0

    do i = 1, domain%num_cell_local
      do j = 1, domain%num_atom(i)
        iatmcls = atmcls(j,i)
        if(enefunc%fep_topology == 2) then
          if (domain%fepgrp(j,i) == 5) then
            atype_preserve(iatmcls) = atype_preserve(iatmcls)+1
            num_all_atoms_preserve = num_all_atoms_preserve + 1
          else if (domain%fepgrp(j,i) == 1) then
            atype_preserve(iatmcls) = atype_preserve(iatmcls)+1
            num_all_atoms_preserve = num_all_atoms_preserve + 1
          else if (domain%fepgrp(j,i) == 2) then
            atype_preserve(iatmcls) = atype_preserve(iatmcls)+1
            num_all_atoms_preserve = num_all_atoms_preserve + 1
          else if (domain%fepgrp(j,i) == 3) then
            atype_vanish(iatmcls) = atype_vanish(iatmcls)+1
            num_all_atoms_vanish = num_all_atoms_vanish + 1
          else if (domain%fepgrp(j,i) == 4) then
            atype_appear(iatmcls) = atype_appear(iatmcls)+1
            num_all_atoms_appear = num_all_atoms_appear + 1
          end if
        else
          if (domain%fepgrp(j,i) == 5) then
            atype_preserve(iatmcls) = atype_preserve(iatmcls)+1
            num_all_atoms_preserve = num_all_atoms_preserve + 1
          else if (domain%fepgrp(j,i) == 1) then
            atype_vanish(iatmcls) = atype_vanish(iatmcls)+1
            num_all_atoms_vanish = num_all_atoms_vanish + 1
          else if (domain%fepgrp(j,i) == 2) then
            atype_appear(iatmcls) = atype_appear(iatmcls)+1
            num_all_atoms_appear = num_all_atoms_appear + 1
          else if (domain%fepgrp(j,i) == 3) then
            atype_vanish(iatmcls) = atype_vanish(iatmcls)+1
            num_all_atoms_vanish = num_all_atoms_vanish + 1
          else if (domain%fepgrp(j,i) == 4) then
            atype_appear(iatmcls) = atype_appear(iatmcls)+1
            num_all_atoms_appear = num_all_atoms_appear + 1
          end if
        end if
      end do
!      num_all_atoms = num_all_atoms + domain%num_atom(i)
    end do

#ifdef HAVE_MPI_GENESIS
!    call mpi_allreduce(mpi_in_place, atype, ntypes, mpi_integer, &
!                       mpi_sum, mpi_comm_country, ierror)
    call mpi_allreduce(mpi_in_place, atype_preserve, ntypes, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
    call mpi_allreduce(mpi_in_place, atype_vanish, ntypes, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
    call mpi_allreduce(mpi_in_place, atype_appear, ntypes, mpi_integer, &
                       mpi_sum, mpi_comm_country, ierror)
#endif

    lj6_tot_preserve = 0.0_dp
    lj6_tot_vanish   = 0.0_dp
    lj6_tot_appear   = 0.0_dp
    do i = 1, ntypes
      do j = 1, ntypes
!        lj6_tot = lj6_tot + enefunc%nonb_lj6(i,j)*atype(i)*atype(j)
        lj6_tot_preserve = lj6_tot_preserve + enefunc%nonb_lj6(i,j)*atype_preserve(i)*atype_preserve(j)
        lj6_tot_vanish   = lj6_tot_vanish   + enefunc%nonb_lj6(i,j)*atype_vanish(i)*atype_preserve(j)
        lj6_tot_vanish   = lj6_tot_vanish   + enefunc%nonb_lj6(i,j)*atype_preserve(i)*atype_vanish(j)
        lj6_tot_vanish   = lj6_tot_vanish   + enefunc%nonb_lj6(i,j)*atype_vanish(i)*atype_vanish(j)
        lj6_tot_appear   = lj6_tot_appear   + enefunc%nonb_lj6(i,j)*atype_appear(i)*atype_preserve(j)
        lj6_tot_appear   = lj6_tot_appear   + enefunc%nonb_lj6(i,j)*atype_preserve(i)*atype_appear(j)
        lj6_tot_appear   = lj6_tot_appear   + enefunc%nonb_lj6(i,j)*atype_appear(i)*atype_appear(j)
      end do
    end do

!    deallocate(atype)
    deallocate(atype_preserve, atype_vanish, atype_appear)

    ! FEP: Make a table for perturbation
    do fg1 = 1, 5
      do fg2 = 1, 5
        if (enefunc%fep_topology == 2) then
          if (((fg1==5).and.(fg2==5)) .or. &
            ((fg1==1).and.(fg2==1)) .or. &
            ((fg1==1).and.(fg2==5)) .or. &
            ((fg1==5).and.(fg2==1)) .or. &
            ((fg1==2).and.(fg2==2)) .or. &
            ((fg1==2).and.(fg2==5)) .or. &
            ((fg1==5).and.(fg2==2))) then
            ! preserved-preserved
            table_pert(fg1,fg2) = 5
          else if (((fg1==3).and.(fg2==3)) .or. &
            ((fg1==3).and.(fg2==5)) .or. &
            ((fg1==5).and.(fg2==3)) .or. &
            ((fg1==3).and.(fg2==1)) .or. &
            ((fg1==1).and.(fg2==3))) then
            ! dualA-dualA and dualA-other
            table_pert(fg1,fg2) = 3
          else if (((fg1==4).and.(fg2==4)) .or. &
            ((fg1==4).and.(fg2==5)) .or. &
            ((fg1==5).and.(fg2==4)) .or. &
            ((fg1==4).and.(fg2==2)) .or. &
            ((fg1==2).and.(fg2==4))) then
            ! dualB-dualB and dualB-other
            table_pert(fg1,fg2) = 4
          end if
        else
          if ((fg1==5).and.(fg2==5)) then
            ! preserved-preserved
            table_pert(fg1,fg2) = 5
          else if (((fg1==1).and.(fg2==1)) .or. &
            ((fg1==1).and.(fg2==5)) .or. &
            ((fg1==5).and.(fg2==1))) then
            ! singleA-singleA and singleA-preserved
            table_pert(fg1,fg2) = 1
          else if (((fg1==2).and.(fg2==2)) .or. &
            ((fg1==2).and.(fg2==5)) .or. &
            ((fg1==5).and.(fg2==2))) then
            ! singleB-singleB and singleB-preserved
            table_pert(fg1,fg2) = 2
          else if (((fg1==3).and.(fg2==3)) .or. &
            ((fg1==3).and.(fg2==5)) .or. &
            ((fg1==5).and.(fg2==3)) .or. &
            ((fg1==3).and.(fg2==1)) .or. &
            ((fg1==1).and.(fg2==3))) then
            ! dualA-dualA and dualA-other
            table_pert(fg1,fg2) = 3
          else if (((fg1==4).and.(fg2==4)) .or. &
            ((fg1==4).and.(fg2==5)) .or. &
            ((fg1==5).and.(fg2==4)) .or. &
            ((fg1==4).and.(fg2==2)) .or. &
            ((fg1==2).and.(fg2==4))) then
            ! dualB-dualB and dualB-other
            table_pert(fg1,fg2) = 4
          end if
        end if
      end do
    end do

    cutoff       = enefunc%cutoffdist
    cutoff2      = cutoff*cutoff
    cutoff3      = cutoff2*cutoff
    inv_cutoff3  = 1.0_dp/cutoff3

    eswitch = 0.0_dp
    vswitch = 0.0_dp
    vlong   = inv_cutoff3/3.0_dp

    if (enefunc%forcefield == ForcefieldAMBER ) then

      factor_preserve = 2.0_dp*PI*lj6_tot_preserve
      enefunc%dispersion_energy_preserve = -factor_preserve*vlong
      enefunc%dispersion_virial_preserve = -2.0_dp*factor_preserve*vlong

      factor_vanish = 2.0_dp*PI*lj6_tot_vanish
      enefunc%dispersion_energy_vanish = -factor_vanish*vlong
      enefunc%dispersion_virial_vanish = -2.0_dp*factor_vanish*vlong

      factor_appear = 2.0_dp*PI*lj6_tot_appear
      enefunc%dispersion_energy_appear = -factor_appear*vlong
      enefunc%dispersion_virial_appear = -2.0_dp*factor_appear*vlong

    else if (enefunc%forcefield == ForcefieldGROAMBER) then
      !
      ! remove exclusion
      !
!      lj6_ex = 0.0_dp
      lj6_ex_preserve = 0.0_dp
      lj6_ex_vanish = 0.0_dp
      lj6_ex_appear = 0.0_dp
!      nexpair = 0
      nexpair_preserve = 0
      nexpair_vanish = 0
      nexpair_appear = 0
      do i = 1, domain%num_cell_local
        ! self
        do j = 1, domain%num_atom(i)
          iatmcls = atmcls(j,i)
!          lj6_ex  = lj6_ex + enefunc%nonb_lj6(iatmcls,iatmcls)

          ! FEP
          fg1 = domain%fepgrp(j,i) 
          pert_flag = table_pert(fg1,fg1)
          if (pert_flag == 5) then
            lj6_ex_preserve = lj6_ex_preserve + &
              enefunc%nonb_lj6(iatmcls,iatmcls)
            nexpair_preserve = nexpair_preserve + 1
          else if ((pert_flag == 1) .or. (pert_flag == 3)) then
            lj6_ex_vanish = lj6_ex_vanish + &
              enefunc%nonb_lj6(iatmcls,iatmcls)
            nexpair_vanish = nexpair_vanish + 1
          else if ((pert_flag == 2) .or. (pert_flag == 4)) then
            lj6_ex_appear = lj6_ex_appear + &
              enefunc%nonb_lj6(iatmcls,iatmcls)
            nexpair_appear = nexpair_appear + 1
          end if

        end do

        ! bonds
        do j = 1, enefunc%num_bond(i)
          icel1 = id_g2l(1,bondlist(1,j,i))
          i1    = id_g2l(2,bondlist(1,j,i))
          icel2 = id_g2l(1,bondlist(2,j,i))
          i2    = id_g2l(2,bondlist(2,j,i))
!          lj6_ex= lj6_ex + enefunc%nb14_lj6(atmcls(i1,icel1),atmcls(i2,icel2))

          ! FEP
          fg1 = domain%fepgrp(i1,icel1) 
          fg2 = domain%fepgrp(i2,icel2) 
          pert_flag = table_pert(fg1,fg2)
          if (pert_flag == 5) then
            lj6_ex_preserve = lj6_ex_preserve + &
              enefunc%nb14_lj6(atmcls(i1,icel1),atmcls(i2,icel2))
            nexpair_preserve = nexpair_preserve + 1
          else if ((pert_flag == 1) .or. (pert_flag == 3)) then
            lj6_ex_vanish = lj6_ex_vanish + &
              enefunc%nb14_lj6(atmcls(i1,icel1),atmcls(i2,icel2))
            nexpair_vanish = nexpair_vanish + 1
          else if ((pert_flag == 2) .or. (pert_flag == 4)) then
            lj6_ex_appear = lj6_ex_appear + &
              enefunc%nb14_lj6(atmcls(i1,icel1),atmcls(i2,icel2))
            nexpair_appear = nexpair_appear + 1
          end if

        end do

        ! angles
        do j = 1, enefunc%num_angle(i)
          icel1 = id_g2l(1,anglelist(1,j,i))
          i1    = id_g2l(2,anglelist(1,j,i))
          icel3 = id_g2l(1,anglelist(3,j,i))
          i3    = id_g2l(2,anglelist(3,j,i))
!          lj6_ex= lj6_ex + enefunc%nb14_lj6(atmcls(i1,icel1),atmcls(i3,icel3))

          ! FEP
          fg1 = domain%fepgrp(i1,icel1) 
          fg2 = domain%fepgrp(i3,icel3) 
          pert_flag = table_pert(fg1,fg2)
          if (pert_flag == 5) then
            lj6_ex_preserve = lj6_ex_preserve + &
              enefunc%nb14_lj6(atmcls(i1,icel1),atmcls(i3,icel3))
            nexpair_preserve = nexpair_preserve + 1
          else if ((pert_flag == 1) .or. (pert_flag == 3)) then
            lj6_ex_vanish = lj6_ex_vanish + &
              enefunc%nb14_lj6(atmcls(i1,icel1),atmcls(i3,icel3))
            nexpair_vanish = nexpair_vanish + 1
          else if ((pert_flag == 2) .or. (pert_flag == 4)) then
            lj6_ex_appear = lj6_ex_appear + &
              enefunc%nb14_lj6(atmcls(i1,icel1),atmcls(i3,icel3))
            nexpair_appear = nexpair_appear + 1
          end if

        end do

        ! dihedral
        do j = 1, enefunc%num_dihedral(i)
          icel1 = id_g2l(1,dihelist(1,j,i))
          i1    = id_g2l(2,dihelist(1,j,i))
          icel4 = id_g2l(1,dihelist(4,j,i))
          i4    = id_g2l(2,dihelist(4,j,i))
!          lj6_ex= lj6_ex + enefunc%nb14_lj6(atmcls(i1,icel1),atmcls(i4,icel4))

          ! FEP
          fg1 = domain%fepgrp(i1,icel1) 
          fg2 = domain%fepgrp(i4,icel4) 
          pert_flag = table_pert(fg1,fg2)
          if (pert_flag == 5) then
            lj6_ex_preserve = lj6_ex_preserve + &
              enefunc%nb14_lj6(atmcls(i1,icel1),atmcls(i4,icel4))
            nexpair_preserve = nexpair_preserve + 1
          else if ((pert_flag == 1) .or. (pert_flag == 3)) then
            lj6_ex_vanish = lj6_ex_vanish + &
              enefunc%nb14_lj6(atmcls(i1,icel1),atmcls(i4,icel4))
            nexpair_vanish = nexpair_vanish + 1
          else if ((pert_flag == 2) .or. (pert_flag == 4)) then
            lj6_ex_appear = lj6_ex_appear + &
              enefunc%nb14_lj6(atmcls(i1,icel1),atmcls(i4,icel4))
            nexpair_appear = nexpair_appear + 1
          end if

        end do

        ! RB dihedral
        do j = 1, enefunc%num_rb_dihedral(i)
          icel1 = id_g2l(1,rb_dihelist(1,j,i))
          i1    = id_g2l(2,rb_dihelist(1,j,i))
          icel4 = id_g2l(1,rb_dihelist(4,j,i))
          i4    = id_g2l(2,rb_dihelist(4,j,i))
!          lj6_ex= lj6_ex + enefunc%nb14_lj6(atmcls(i1,icel1),atmcls(i4,icel4))

          ! FEP
          fg1 = domain%fepgrp(i1,icel1) 
          fg2 = domain%fepgrp(i4,icel4) 
          pert_flag = table_pert(fg1,fg2)
          if (pert_flag == 5) then
            lj6_ex_preserve = lj6_ex_preserve + &
              enefunc%nb14_lj6(atmcls(i1,icel1),atmcls(i4,icel4))
            nexpair_preserve = nexpair_preserve + 1
          else if ((pert_flag == 1) .or. (pert_flag == 3)) then
            lj6_ex_vanish = lj6_ex_vanish + &
              enefunc%nb14_lj6(atmcls(i1,icel1),atmcls(i4,icel4))
            nexpair_vanish = nexpair_vanish + 1
          else if ((pert_flag == 2) .or. (pert_flag == 4)) then
            lj6_ex_appear = lj6_ex_appear + &
              enefunc%nb14_lj6(atmcls(i1,icel1),atmcls(i4,icel4))
            nexpair_appear = nexpair_appear + 1
          end if

        end do

        ! improper
        do j = 1, enefunc%num_improper(i)
          icel1 = id_g2l(1,imprlist(1,j,i))
          i1    = id_g2l(2,imprlist(1,j,i))
          icel4 = id_g2l(1,imprlist(4,j,i))
          i4    = id_g2l(2,imprlist(4,j,i))
!          lj6_ex= lj6_ex + enefunc%nb14_lj6(atmcls(i1,icel1),atmcls(i4,icel4))

          ! FEP
          fg1 = domain%fepgrp(i1,icel1) 
          fg2 = domain%fepgrp(i4,icel4) 
          pert_flag = table_pert(fg1,fg2)
          if (pert_flag == 5) then
            lj6_ex_preserve = lj6_ex_preserve + &
              enefunc%nb14_lj6(atmcls(i1,icel1),atmcls(i4,icel4))
            nexpair_preserve = nexpair_preserve + 1
          else if ((pert_flag == 1) .or. (pert_flag == 3)) then
            lj6_ex_vanish = lj6_ex_vanish + &
              enefunc%nb14_lj6(atmcls(i1,icel1),atmcls(i4,icel4))
            nexpair_vanish = nexpair_vanish + 1
          else if ((pert_flag == 2) .or. (pert_flag == 4)) then
            lj6_ex_appear = lj6_ex_appear + &
              enefunc%nb14_lj6(atmcls(i1,icel1),atmcls(i4,icel4))
            nexpair_appear = nexpair_appear + 1
          end if

        end do

!        nexpair = nexpair + domain%num_atom(i)        &
!                          + enefunc%num_bond(i)        &
!                          + enefunc%num_angle(i)       &
!                          + enefunc%num_dihedral(i)    &
!                          + enefunc%num_rb_dihedral(i) &
!                          + enefunc%num_improper(i)
      end do
#ifdef HAVE_MPI_GENESIS
!      call mpi_allreduce(mpi_in_place, num_all_atoms, 1, mpi_integer, &
!                         mpi_sum, mpi_comm_country, ierror)
      call mpi_allreduce(mpi_in_place, num_all_atoms_preserve, 1, mpi_integer, &
                         mpi_sum, mpi_comm_country, ierror)
      call mpi_allreduce(mpi_in_place, num_all_atoms_vanish, 1, mpi_integer, &
                         mpi_sum, mpi_comm_country, ierror)
      call mpi_allreduce(mpi_in_place, num_all_atoms_appear, 1, mpi_integer, &
                         mpi_sum, mpi_comm_country, ierror)

!      call mpi_allreduce(mpi_in_place, nexpair, 1, mpi_integer, &
!                         mpi_sum, mpi_comm_country, ierror)
      call mpi_allreduce(mpi_in_place, nexpair_preserve, 1, mpi_integer, &
                         mpi_sum, mpi_comm_country, ierror)
      call mpi_allreduce(mpi_in_place, nexpair_vanish, 1, mpi_integer, &
                         mpi_sum, mpi_comm_country, ierror)
      call mpi_allreduce(mpi_in_place, nexpair_appear, 1, mpi_integer, &
                         mpi_sum, mpi_comm_country, ierror)

!      call mpi_allreduce(mpi_in_place, lj6_ex, 1, mpi_real8, &
!                         mpi_sum, mpi_comm_country, ierror)
      call mpi_allreduce(mpi_in_place, lj6_ex_preserve, 1, mpi_real8, &
                         mpi_sum, mpi_comm_country, ierror)
      call mpi_allreduce(mpi_in_place, lj6_ex_vanish, 1, mpi_real8, &
                         mpi_sum, mpi_comm_country, ierror)
      call mpi_allreduce(mpi_in_place, lj6_ex_appear, 1, mpi_real8, &
                         mpi_sum, mpi_comm_country, ierror)
#endif
!      lj6_diff = (lj6_tot - lj6_ex)
      lj6_diff_preserve = lj6_tot_preserve - lj6_ex_preserve
      lj6_diff_vanish = lj6_tot_vanish - lj6_ex_vanish
      lj6_diff_appear = lj6_tot_appear - lj6_ex_appear

!      natom2 = num_all_atoms*num_all_atoms
      natom_preserve2 = num_all_atoms_preserve*num_all_atoms_preserve
      natom_vanish2 = num_all_atoms_preserve*num_all_atoms_vanish &
        + num_all_atoms_vanish*num_all_atoms_vanish
      natom_appear2 = num_all_atoms_preserve*num_all_atoms_appear &
        + num_all_atoms_appear*num_all_atoms_appear
!      rpair  = real(natom2/(natom2-nexpair),dp)
      if (natom_preserve2 /= nexpair_preserve) then
        rpair_preserve  &
          = real(natom_preserve2/(natom_preserve2-nexpair_preserve),dp)
      else
        rpair_preserve = 0.0_dp
      end if
      if (natom_vanish2 /= nexpair_vanish) then
        rpair_vanish  &
          = real(natom_vanish2/(natom_vanish2-nexpair_vanish),dp)
      else
        rpair_vanish = 0.0_dp
      end if
      if (natom_appear2 /= nexpair_appear) then
        rpair_appear  &
          = real(natom_appear2/(natom_appear2-nexpair_appear),dp)
      else
        rpair_appear = 0.0_dp
      end if

!      factor       = 2.0_dp*PI*rpair*lj6_diff
      factor_preserve = 2.0_dp*PI*rpair_preserve*lj6_diff_preserve
      factor_vanish = 2.0_dp*PI*rpair_vanish*lj6_diff_vanish
      factor_appear = 2.0_dp*PI*rpair_appear*lj6_diff_appear

      switchdist   = enefunc%switchdist
      diff_cs      = (cutoff - switchdist)

      if (diff_cs > EPS) then

        if (enefunc%vdw_shift) then
          cutoff4      = cutoff3*cutoff
          cutoff5      = cutoff4*cutoff
          cutoff6      = cutoff5*cutoff
          cutoff7      = cutoff6*cutoff
          cutoff8      = cutoff7*cutoff
          cutoff14     = cutoff7*cutoff7
          inv_cutoff6  = inv_cutoff3*inv_cutoff3
          inv_cutoff12 = inv_cutoff6*inv_cutoff6
  
          diff_cs2     = diff_cs*diff_cs
          diff_cs3     = diff_cs2*diff_cs
          diff_cs4     = diff_cs3*diff_cs
  
          switchdist2  = switchdist*switchdist
          switchdist3  = switchdist2*switchdist
          switchdist4  = switchdist3*switchdist
          switchdist5  = switchdist4*switchdist
          switchdist6  = switchdist5*switchdist
          switchdist7  = switchdist6*switchdist
          switchdist8  = switchdist7*switchdist
  
          ! LJ6
          !
          shift_a = -(10.0_dp*cutoff - 7.0_dp*switchdist)/(cutoff8*diff_cs2)
          shift_b =  ( 9.0_dp*cutoff - 7.0_dp*switchdist)/(cutoff8*diff_cs3)
  
          shift_c = inv_cutoff6 - 2.0_dp * shift_a * diff_cs3  &
                    - 1.5_dp * shift_b * diff_cs4
  
          eswitch = -2.0_dp * shift_a * ((1.0_dp/6.0_dp)*cutoff6                &
                                        -(3.0_dp/5.0_dp)*cutoff5*switchdist     &
                                        +(3.0_dp/4.0_dp)*cutoff4*switchdist2    &
                                        -(1.0_dp/3.0_dp)*cutoff3*switchdist3    &
                                        +(1.0_dp/6.0e1_dp)*switchdist6)         &
                    -1.5_dp * shift_b * ((1.0_dp/7.0_dp)*cutoff7                &
                                        -(2.0_dp/3.0_dp)*cutoff6*switchdist     &
                                        +(6.0_dp/5.0_dp)*cutoff5*switchdist2    &
                                        -                cutoff4*switchdist3    &
                                        +(1.0_dp/3.0_dp)*cutoff3*switchdist4    &
                                        -(1.0_dp/1.05e2_dp)*switchdist7)        &
                    -(1.0_dp/3.0_dp) * shift_c * (cutoff3)
    
          ! LJ12
          !
          shift_a = -(16.0_dp*cutoff - 13.0_dp*switchdist)/(cutoff14*diff_cs2)
          shift_b =  (15.0_dp*cutoff - 13.0_dp*switchdist)/(cutoff14*diff_cs3)
          shift_c = inv_cutoff12 - 2.0_dp * shift_a * diff_cs3  &
                    - 1.5_dp * shift_b * diff_cs4
  
 
          shift_a = -(10.0_dp*cutoff - 7.0_dp*switchdist)/(cutoff8*diff_cs2)
          shift_b =  ( 9.0_dp*cutoff - 7.0_dp*switchdist)/(cutoff8*diff_cs3)
 
          vswitch = shift_a * ( (1.0_dp/6.0_dp)*cutoff6                         &
                               -(2.0_dp/5.0_dp)*cutoff5*switchdist              &
                               +(1.0_dp/4.0_dp)*cutoff4*switchdist2             &
                               -(1.0_dp/6.0e1_dp)*switchdist6)                  &
                   +shift_b * ( (1.0_dp/7.0_dp)*cutoff7                         &
                               -(1.0_dp/2.0_dp)*cutoff6*switchdist              &
                               +(3.0_dp/5.0_dp)*cutoff5*switchdist2             &
                               -(1.0_dp/4.0_dp)*cutoff4*switchdist3             &
                               +(1.0_dp/1.4e2_dp)*switchdist7)
!          enefunc%dispersion_energy = factor*(eswitch-vlong)
!          enefunc%dispersion_virial = -2.0_dp*factor*(-vswitch+vlong)
          enefunc%dispersion_energy_preserve = factor_preserve*(eswitch-vlong)
          enefunc%dispersion_virial_preserve &
            = -2.0_dp*factor_preserve*(-vswitch+vlong)
          enefunc%dispersion_energy_vanish = factor_vanish*(eswitch-vlong)
          enefunc%dispersion_virial_vanish &
            = -2.0_dp*factor_vanish*(-vswitch+vlong)
          enefunc%dispersion_energy_appear = factor_appear*(eswitch-vlong)
          enefunc%dispersion_virial_appear &
            = -2.0_dp*factor_appear*(-vswitch+vlong)

        else

          eswitch = enefunc%eswitch
          vswitch = enefunc%vswitch
!          enefunc%dispersion_energy = factor*(eswitch-vlong)
!          enefunc%dispersion_virial = -factor*(vswitch+vlong)
          enefunc%dispersion_energy_preserve = factor_preserve*(eswitch-vlong)
          enefunc%dispersion_virial_preserve = -factor_preserve*(vswitch+vlong)
          enefunc%dispersion_energy_vanish = factor_vanish*(eswitch-vlong)
          enefunc%dispersion_virial_vanish = -factor_vanish*(vswitch+vlong)
          enefunc%dispersion_energy_appear = factor_appear*(eswitch-vlong)
          enefunc%dispersion_virial_appear = -factor_appear*(vswitch+vlong)

        end if

      else 

!        enefunc%dispersion_energy = factor*(eswitch-vlong)
!        enefunc%dispersion_virial = -2.0_dp*factor*(-vswitch+vlong)
        enefunc%dispersion_energy_preserve = factor_preserve*(eswitch-vlong)
        enefunc%dispersion_virial_preserve &
          = -2.0_dp*factor_preserve*(-vswitch+vlong)
        enefunc%dispersion_energy_vanish = factor_vanish*(eswitch-vlong)
        enefunc%dispersion_virial_vanish &
          = -2.0_dp*factor_vanish*(-vswitch+vlong)
        enefunc%dispersion_energy_appear = factor_appear*(eswitch-vlong)
        enefunc%dispersion_virial_appear &
          = -2.0_dp*factor_appear*(-vswitch+vlong)

      end if

    else

      call error_msg('Setup_Enefunc_DispCorr> This force field is not allowed in FEP')

    end if

  end subroutine setup_enefunc_dispcorr_fep

end module sp_enefunc_mod
