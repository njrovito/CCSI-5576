
# -----------------------------------------------------
#
#               Makefile for lookup
#
# Run this build process by entering "make" at the
# command prompt.
#
# -----------------------------------------------------
#
# -----------------------------------------------------
#               N O T E S
# -----------------------------------------------------
#
# -----------------------------------------------------
# On local Ubuntu install MPI wrappers, if needed
# -----------------------------------------------------
#
# Be sure to use MPI compiler wrappers.
# Whenever compiling an MPI program, you should use the MPI wrappers:
#
#    C - mpicc
#    C++ - mpiCC, mpicxx, mpic++
#    FORTRAN - mpifort, mpif77, mpif90
#
# These wrappers do all of the dirty work for you of making sure that
# all of the appropriate compiler flags, libraries, include directories,
# library directories, etc. are included when you compile your program.
#
# -----------------------------------------------------
# On Alpine, be sure to...
# -----------------------------------------------------
#
# module load intel
# module load impi
#
# -----------------------------------------------------
# It might be handy to set the compiler based on
# platform
# -----------------------------------------------------

# Default compiler is Intel MPI compiler

COMPILER := mpicc
UBUNTU   := ubuntu

# Change if on my local system.  Change "scott" to your
# user name on your own computer.

ifeq ($(strip $(LOGNAME)),scott)
COMPILER := mpicxx
endif

# -----------------------------------------------------
# Make lookup
# -----------------------------------------------------

lookup: lookup.cpp 
	$(COMPILER) lookup.cpp -o lookup

