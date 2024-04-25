!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pc_option_str_mod
!> @brief   definition of option
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module pc_option_str_mod

  use select_atoms_str_mod
  use string_mod
  use messages_mod

  implicit none
  private

  ! structures
  type, public :: s_option

    logical             :: check_only
    integer             :: trjout_format
    integer             :: trjout_type
    character(MaxLine)  :: trjout_atom_exp
    type(s_selatoms)    :: trjout_atom
    type(s_selatoms)    :: trjout_atom_trj
    integer             :: pbcc_mode

  end type s_option

  ! subroutines
  public  :: dealloc_option

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_option
  !> @brief        deallocate option
  !! @authors      NT
  !! @param[inout] option : information of option
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine dealloc_option(option)

    ! formal argments
    type(s_option),          intent(inout) :: option


    call dealloc_selatoms(option%trjout_atom)
    call dealloc_selatoms(option%trjout_atom_trj)

    return

  end subroutine dealloc_option

end module pc_option_str_mod
