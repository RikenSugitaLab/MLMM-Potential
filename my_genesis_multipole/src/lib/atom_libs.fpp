!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   atom_lib_mod
!> @brief   utilities of atom libraries
!! @authors Daisuke Matsuoka (DM), Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module atom_libs_mod

  use messages_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public :: atomic_number
  public :: atomic_number_by_name

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    get_atomic_number
  !> @brief        define atomic number from molecule
  !! @authors      TM
  !! @param[in]    mass      : mass
  !! @param[in]    atom_type : atom class name
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function atomic_number(mass, atom_type)

    ! formal arguments
    real(wp),              intent(in)    :: mass
    character(6),          intent(inout) :: atom_type

    ! local variables
    integer                              :: atomic_number
    integer                              :: z, int_mass

    z = 0
    int_mass = nint(mass)

    if      (int_mass == 1  ) then  ! H
      z = 1
    else if (int_mass == 7  ) then  ! Li
      z = 3
    else if (int_mass == 12 ) then  ! C
      z = 6
    else if (int_mass == 14 ) then  ! N
      z = 7
    else if (int_mass == 16 ) then  ! O
      z = 8
    else if (int_mass == 19 ) then  ! F
      z = 9
    else if (int_mass == 23 ) then  ! Na
      z = 11
    else if (int_mass == 24 ) then  ! Mg
      z = 12
    else if (int_mass == 27 ) then  ! Al
      z = 13
    else if (int_mass == 31 ) then  ! P
      z = 15
    else if (int_mass == 32 ) then  ! S
      z = 16
    else if (int_mass == 35 ) then  ! Cl
      z = 17
    else if (int_mass == 39 ) then  ! K
      z = 19
    else if (int_mass == 40 ) then  ! Ca
      z = 20
    else if (int_mass == 56 ) then  ! Fe
      z = 26
    else if (int_mass == 64 ) then  ! Cu
      z = 29
    else if (int_mass == 65 ) then  ! Zn
      z = 30
    else if (int_mass == 79 ) then  ! Se
      z = 34
    else if (int_mass == 80 ) then  ! Br
      z = 35
    else if (int_mass == 85 ) then  ! Rb
      z = 37
    else if (int_mass == 112) then  ! Cd
      z = 48
    else if (int_mass == 127) then  ! I
      z = 53
    else if (int_mass == 133) then  ! Cs
      z = 55
    else if (int_mass == 137) then  ! Ba
      z = 56
    else
      z = atomic_number_by_name(atom_type)
    end if

    atomic_number = z

    return

  end function atomic_number

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      atomic_number_by_name
  !> @brief        get atomic number of an atom
  !! @authors      DM
  !! @return       atomic number
  !! @param[in]    atom_type   : atom type
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
              
  function atomic_number_by_name(atom_type)

    ! formal variables
    character(6), intent(in) :: atom_type

    ! local variable
    integer                  :: atomic_number_by_name
    integer                  :: z
    character(3)             :: atom_element


    atom_element = adjustl(atom_type(1:3))

    if      (atom_element(1:1) == 'H'  .or. atom_element(1:1) == 'h'  .or. &
             atom_element(1:2) == 'TH') then       ! amber param19_ipq.dat
      z = 1
  
    else if (atom_element(1:2) == 'Na' .or. atom_element(1:3) == 'SOD') then
      z = 11

    else if (atom_element(1:2) == 'MG') then
      z = 12

    else if (atom_element(1:1) == 'K'  .or. atom_element(1:3) == 'POT') then
      z = 19
  
    else if (atom_element(1:2) == 'Cl' .or.   &     ! parm10.dat
             atom_element(1:2) == 'cl' .or.   &     ! gaff.dat
             atom_element(1:3) == 'CLA') then
      z = 17

    else if (atom_element(1:2) == 'Cu'  .or.  &
             atom_element(1:2) == 'CU') then
      z = 29

    else if (atom_element(1:2) == 'C0'   .or. &     ! amber param10.dat
             atom_element(1:3) == 'CAL') then       ! CHARMM
      z = 20

    else if (atom_element(1:2) == 'FE') then
      z = 26

    else if (atom_element(1:2) == 'ZN' .or. atom_element(1:2) == 'Zn') then
      z = 30

    else if (atom_element(1:1) == 'C'  .or. atom_element(1:1) == 'c') then
      z = 6

    else if (atom_element(1:2) == '2C' .or. atom_element(1:2) == '3C' .or. &
             atom_element(1:2) == 'TG' .or. atom_element(1:2) == 'TJ' .or. &
             atom_element(1:2) == 'TP' .or. atom_element(1:2) == 'TM' .or. &
             atom_element(1:2) == 'TA') then    ! amber param19_ipq.dat
      z = 6
  
    else if (atom_element(1:1) == 'N'  .or. atom_element(1:1) == 'n') then
      z = 7
  
    else if (atom_element(1:2) == 'TN') then    ! amber param14ipq
      z = 7

    else if (atom_element(1:1) == 'O' .or. atom_element(1:1) == 'o') then
      z = 8

    else if (atom_element(1:1) == 'F' .or. atom_element(1:1) == 'f') then
      z = 9
  
    else if (atom_element(1:1) == 'S' .or. atom_element(1:1) == 's') then
      z = 16
  
    else if (atom_element(1:1) == 'P' .or. atom_element(1:1) == 'p') then
      z = 15

    ! lone pair of TIP4P/TiP5P
    else if (atom_element(1:2) == 'LP' .or.  &   ! charmm
             atom_element(1:2) == 'EP' .or.  &   ! amber
             atom_element(1:2) == 'IW') then     ! gromos
      z = 0

    else
      call error_msg('Atomic_Number_By_Name> Unknown atom element ['//trim(atom_type)//']')
  
    end if

    atomic_number_by_name = z

    return

  end function atomic_number_by_name

end module atom_libs_mod
