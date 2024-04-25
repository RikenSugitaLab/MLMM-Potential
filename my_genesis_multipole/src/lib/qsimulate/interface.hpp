#ifdef QSIMULATE
  interface
    subroutine qsimulate_interface(mpicomm, myrank, input, output, natoms, atoms, coord, charges, &
                                   energy, force, dipole, qmcharges, retry, error) &
                                   bind(C)
      use constants_mod
      use iso_c_binding

      integer :: mpicomm, myrank
      character(kind=c_char) :: input(*)
      character(kind=c_char) :: output(*)
      integer(c_int), value :: natoms
      type(c_ptr), value :: atoms
      type(c_ptr), value :: coord
      type(c_ptr), value :: charges
      real(c_double) :: energy
      type(c_ptr), value :: force
      real(c_double) :: dipole(:)
      type(c_ptr), value :: qmcharges
      logical  :: retry
      logical  :: error
    end subroutine qsimulate_interface
  end interface
#endif
