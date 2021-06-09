#
# On alcyone.grid.helsinki.fi, first of all:
#
# module load openmpi/1.6-gcc
#
# On louhi and vuori, first:
#
# module swap PrgEnv-pgi PrgEnv-gnu
#
# On sisu:
#
# module swap PrgEnv-cray PrgEnv-gnu
#
# On taito:
#
# module swap intel gcc/4.7.2
# module swap intelmpi mvapich2/1.8.1
# module load openblas/0.2.6
#
# then:
#
# CC = CC
#
# or:
#
# CC = mpiCC
#
# On other machines:
#
CC=g++
#
# or:
#
# CC = mpiCC

MAIN=sun.cc

SOURCE_FILES= $(MAIN) \
	check_unitarity.cc \
	compute_action.cc \
	compute_action_from_first_boundary.cc \
	compute_action_from_second_boundary.cc \
	deallocate_arrays.cc \
	gauge_transform.cc \
	generate_SUN_matrix_near_identity.cc \
	heat_bath.cc \
	heat_bath_for_one_link.cc \
	initialize.cc \
	locked_gauge_transform.cc \
	measurements.cc \
	measure_correlators_from_averaged_timelike_blocks.cc \
	measure_Jarzynski_SF.cc \
	measure_Polyakov_loops.cc \
	measure_Wilson_loops.cc \
	metropolis.cc \
	multiplication_routines_for_SU*.cc \
	norm_su*.cc \
	nrectangularstaple.cc \
	nstaple.cc \
	open_Polyakov_loop.cc \
	overrelaxation.cc \
	overrelaxation_for_one_link.cc \
	plaquette.cc \
	prectangularstaple.cc \
	pstaple.cc \
	read_last_conf.cc \
	reunitarize.cc \
	save_last_conf.cc \
	simulation.cc \
	slab_average.cc \
	smear.cc \
	space_plaq.cc \
	terminate_run.cc \
	thermalize.cc \
	time_plaq.cc \
	time_to_stop.cc \
	update.cc \
	update_first_boundary.cc \
	update_second_boundary.cc

HEADER_FILES=headers.h \
	more_headers.h \
	parameters.h \
	included_routines.h \
	MersenneTwister.h

# FFLAGS = -O3 -Wall -fstrict-aliasing -ffast-math -fexpensive-optimizations -finline-functions -funroll-loops -fbounds-check

FFLAGS = -O3 -Wall -ffast-math -fexpensive-optimizations -finline-functions -funroll-loops -fbounds-check

INCLUDE =



# On the CSC machines: first:
#
# module swap PrgEnv-pgi PrgEnv-gnu
# then compile (on louhi) with:
# LIB =
#
# or (on vuori) with:
# LIB = -lacml


# On brutus:
# GOTOBLAS2=/cluster/apps/gotoblas2/1.06_par/x86_64/gcc_4.1.2
# LIB= -L$(GOTOBLAS2) -lgoto2 -lgfortran

# On marconi: first:
#
# module load autoload
#
# module load gnu/6.1.0
#
# and for parallel code also:
#
# module load openmpi/1.10.3-threadmultiple--gnu--6.1.0
#
# Then:
#
# module load intel/pe-xe-2017--binary
# module load lapack/3.6.1--intel--pe-xe-2017--binary
# module load blas/3.6.0--intel--pe-xe-2017--binary
#
# then compile with:
# LIB = -L$(LAPACK_LIB) -llapack -L$(BLAS_LIB) -lblas -lm

# On other machines:
LIB = -L/usr/lib/ #-llapack -lm

scalar_main = main_scalar.cc
scalar_sim_files = $(scalar_main) $(SOURCE_FILES)

scalar : $(scalar_sim_files) \
	$(HEADER_FILES) \
	Makefile 
	$(CC) $(scalar_main) $(FFLAGS) $(INCLUDE) $(LIB) -o $@



# EXECUTABLE=run.x

# $(EXECUTABLE): $(SOURCE_FILES) \
# 	$(HEADER_FILES) \
# 	Makefile
# 	$(CC) $(MAIN) $(FFLAGS) $(INCLUDE) $(LIB) -o $@

clean:
	rm -f $(EXECUTABLE)
	rm -f scalar

# To run:
#
#  mpirun ${executable}
#
# (assuming that openmpi is installed).
# or:
# aprun ${executable}
#
# To run a parallel mpi job:
# salloc -n 2 --ntasks-per-node=2 --mem-per-cpu=1000 -t 00:03:00 -p test srun ${executable}
