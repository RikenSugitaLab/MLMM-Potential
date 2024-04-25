!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   rg_option_str_mod
!> @brief   structure of option information
!! @authors Motoshi Kamiya (MK), Takaharu Mori (TM)
! 
!  (c) Copyright 2016 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module rg_option_str_mod

  use select_atoms_str_mod
  use constants_mod
  use string_mod

  implicit none
  private

  ! structures
  type, public :: s_option

    logical            :: check_only
    logical            :: mass_weighted
    character(MaxLine) :: analysis_atom_exp
    type(s_selatoms)   :: analysis_atom

  end type s_option

  public :: dealloc_option

  contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_option
  !> @brief        deallocate option
  !! @authors      MK, TM
  !! @param[inout] option : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_option(option)

    ! formal arguments
    type(s_option),   intent(inout) :: option

    call dealloc_selatoms(option%analysis_atom)

    return

  end subroutine dealloc_option

end module rg_option_str_mod
