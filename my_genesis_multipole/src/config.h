/* src/config.h.  Generated from config.h.in by configure.  */
/* src/config.h.in.  Generated from configure.ac by autoheader.  */

/* c compiler version */
#define COMPILE_CC_VER "icc (ICC) 2021.4.0 20210910"

/* c flags */
#define COMPILE_CFLAGS "-O3 -ip -axCORE-AVX2  -qopenmp"

/* defined variables */
#define COMPILE_DEFINED_VARIABLES " -DHAVE_MPI_GENESIS -DOMP -DFFTE -DLAPACK -DDSFMT_MEXP=19937 -DINTEL"

/* fortran flags */
#define COMPILE_FCFLAGS "-xHost -O3 -ip -mkl=parallel  -assume byterecl -qopenmp "

/* fortran compiler version */
#define COMPILE_FC_VER "ifort (IFORT) 2021.4.0 20210910"

/* genesis version */
#define COMPILE_GENESIS_VERSION "1.6.1 [2021-07-27 01:53:56 +0900]"

/* hostname */
#define COMPILE_HOST "kelp"

/* ld flags */
#define COMPILE_LDFLAGS " -assume byterecl -qopenmp  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_lapack95_lp64 "

/* cuda version */
/* #undef COMPILE_NVCC_VER */

/* username */
#define COMPILE_USER "yklei"

/* defined if cuda_gpu is used. */
/* #undef CUDAGPU */

/* defined if Debug is used. */
/* #undef DEBUG */

/* defined always. */
#define DSFMT_MEXP 19937

/* defined if FFTE is used. */
#define FFTE 1

/* Define to 1 if you have the <bagel.h> header file. */
/* #undef HAVE_BAGEL_H */

/* Define to 1 if you have the <inttypes.h> header file. */
/* #undef HAVE_INTTYPES_H */

/* Define to 1 if you have the <memory.h> header file. */
/* #undef HAVE_MEMORY_H */

/* defined if MPI is used. */
#define HAVE_MPI_GENESIS 1

/* defined if MPI is used. */
/* #undef HAVE_MPI_H */

/* Define to 1 if you have the <stdint.h> header file. */
/* #undef HAVE_STDINT_H */

/* Define to 1 if you have the <stdlib.h> header file. */
/* #undef HAVE_STDLIB_H */

/* Define to 1 if you have the <strings.h> header file. */
/* #undef HAVE_STRINGS_H */

/* Define to 1 if you have the <string.h> header file. */
/* #undef HAVE_STRING_H */

/* Define to 1 if you have the <sys/stat.h> header file. */
/* #undef HAVE_SYS_STAT_H */

/* Define to 1 if you have the <sys/types.h> header file. */
/* #undef HAVE_SYS_TYPES_H */

/* Define to 1 if you have the <unistd.h> header file. */
/* #undef HAVE_UNISTD_H */

/* defined if HM_DISK is used. */
/* #undef HM_DISK */

/* defined if Intel compiler is used. */
#define INTEL 1

/* defined if K-computer compiler is used. */
/* #undef KCOMP */

/* defined if LAPACK is used. */
#define LAPACK 1

/* defined if OpenMP is used. */
#define OMP 1

/* Name of package */
#define PACKAGE "genesis"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "genesis@riken.jp"

/* Define to the full name of this package. */
#define PACKAGE_NAME "genesis"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "genesis 1.6.1"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "genesis"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "1.6.1"

/* defined if pgi and cuda are used. */
/* #undef PGICUDA */

/* defined if QSimulate is used. */
/* #undef QSIMULATE */

/* define if platform is RICC. */
/* #undef RICC */

/* Define to 1 if you have the ANSI C header files. */
/* #undef STDC_HEADERS */

/* defined if gpu is used. */
/* #undef USE_GPU */

/* Version number of package */
#define VERSION "1.6.1"

/* defined if _SINGLE is used. */
/* #undef _SINGLE */

/* defined if GCC gfortran compiler is used. */
/* #undef __GFORTRAN__ */

/* defined if pgi compiler is used. */
/* #undef __PGI */
