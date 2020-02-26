# List of source code file
obj = m_logger.o m_exception.o mt19937.o c_utils.o utils.o \
	  crop_and_extract.o cgal_delaunay.o cgal_delaunay_wrapper.o \
	  kdtree2.o settings.o likelihood_settings.o \
	  fastMarching_wrapper.o run_time_info.o \
	  likelihood_surf.o likelihood_body.o likelihood.o \
	  read_write.o slice_sampling2.o hmc.o mcmc_loc2.o hmcmc.o

# Likelihood program

# Whether to run in debugging mode
DEBUG =     

# Openmp support
OMP = 1  

# MPI support
MPI = 1     

# netcdf support
NETCDF =

# Set the compiler
COMPILER_TYPE = gnu

# Set the specific compiler to be that defined in COMPILER_TYPE
ifeq ($(COMPILER_TYPE),intel)
F90 = ifort
CC = icc
endif

ifeq ($(COMPILER_TYPE),gnu)
F90 = gfortran
CC = g++
endif

ifdef MPI
F90 = mpif90
FFLAGS += -DMPI
endif

# Archive toll
AR = ar r

# gfortran settings
# =================
ifeq ($(COMPILER_TYPE),gnu)
#FFLAGS += -ffree-line-length-0 -cpp 
FFLAGS += -ffree-line-length-0 -ffixed-line-length-0 -cpp
#FFLAGS += -fopenmp

ifdef DEBUG
FFLAGS += -O0 -g -Wall -Wtabs -Wextra -Wconversion -fimplicit-none -fbacktrace -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan 
F77FLAGS = -O0 -g -ffree-line-length-0 -ffixed-line-length-0 -cpp
CFLAGS += -O0 -g -frounding-math -Wall -Wextra -ansi -pedantic -Wc++0x-compat -Wlong-long -std=c++11 
else
# gnu -Ofast sometimes fails to produce the right code. 
# If speed is needed, try -Ofast. If it fails, use -O3. 
FFLAGS += -O3 
CFLAGS += -O3 -std=c++11 -frounding-math
endif

ifdef OMP
CFLAGS += -fopenmp
FFLAGS += -fopenmp
endif

endif

# ifort settings
# =================
ifeq ($(COMPILER_TYPE),intel)
FFLAGS += -extend-source -fpp

ifdef DEBUG
FFLAGS += -O3 -g -warn all -traceback -check all,noarg_temp_created -debug all -gen-interfaces -warn-interfaces -assume realloc-lhs
#F77FLAGS = -O0 -g 
CFLAGS += -O0 -g -std=c++11 -traceback -debug all
else
FFLAGS += -O3 -g
CFLAGS += -O3 -std=c++11 
endif

ifdef OMP
CFLAGS += -qopenmp
FFLAGS += -qopenmp
endif

endif

# Where bianries are stored
BIN_DIR = ../bin

# Likelihood program
LIKE3D = fmm3d
LIKE2D = fmm2d
SURF = surfmodes

# Where forward modeling code are stored
LIKE2D_DIR = ../fm2d
LIKE3D_DIR = ../fastmarching
SURF_DIR = ../surfmodes

# Add forward modeling library to the program
INCFLAGS += -I$(LIKE2D_DIR) -I$(LIKE3D_DIR) -I$(SURF_DIR)
LDFLAGS += -L$(LIKE2D_DIR) -L$(LIKE3D_DIR) -L$(SURF_DIR)
# Add cgal related library to the program (boost, cgal)
INCFLAGS += -I/home/s1575033/Program/include
#INCFLAGS += -I/exports/igmm/software/pkg/el7/libs/boost/1.59.0/include -I/exports/csce/eddie/geos/groups/eddie_geos_eps_eip/Xin/usr/include
#LDFLAGS += -L/exports/igmm/software/pkg/el7/libs/boost/1.59.0/lib -L/exports/csce/eddie/geos/groups/eddie_geos_eps_eip/Xin/usr/lib64 
LDLIBS += -l$(LIKE2D) -l$(LIKE3D) -l$(SURF) -lboost_system -lCGAL -lgmp -lmpfr -lm -lstdc++
#LDLIBS += -l$(LIKE2D) -l$(LIKE3D) -l$(SURF)  -lgmp -lmpfr -lCGAL -lm -lstdc++

# likelihood library
LIKE2D_LIB = $(patsubst %,lib%.a,$(LIKE2D))
LIKE3D_LIB = $(patsubst %,lib%.a,$(LIKE3D))
SURF_LIB =$(patsubst %,lib%.a,$(SURF))

# NETCDF
ifdef NETCDF
obj += netcdf_read_write.o initialise.o
FFLAGS += -DNETCDF
INCFLAGS += -I/exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/usr/include
LDFLAGS += -L/exports/csce/datasotre/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/usr/lib
#INCFLAGS += -I/exports/csce/eddie/geos/groups/eddie_geos_eps_eip/Xin/usr/include
#LDFLAGS += -L/exports/csce/eddie/geos/groups/eddie_geos_eps_eip/Xin/usr/lib
LDLIBS += -lnetcdff -lnetcdf
else
obj += initialise.o
endif

# Export variables
export DEBUG F90 CC FFLAGS F77FLAGS AR CFLAGS  LIKE2D LIKE2D_LIB LIKE3D LIKE3D_LIB SURF SURF_LIB


#------------------------------------------------------------------------
all :   $(LIKE2D_LIB) $(LIKE3D_LIB) $(SURF_LIB)  MCTomo sample 

MCTomo  :: $(obj) MCTomo.o
	$(F90) $(FFLAGS) $(LDFLAGS) $^  $(LDLIBS) -o $(BIN_DIR)/$@

sample    :: $(obj) sample.o
	$(F90) $(FFLAGS) $(LDFLAGS) $^  $(LDLIBS) -o $(BIN_DIR)/$@

sample_test    :: $(obj) sample_test.o
	$(F90) $(FFLAGS) $(LDFLAGS) $^  $(LDLIBS) -o $(BIN_DIR)/$@

%.o : %.F90
	$(F90) $(INCLIKE) $(INCFLAGS) $(FFLAGS) -c -o $@ $<

%.o : %.f90
	$(F90) $(INCLIKE) $(INCFLAGS) $(FFLAGS) -c -o $@ $<

%.o : %.f
	$(F90) $(F77FLAGS) -c -o $@ $<

%.o : %.cpp
	$(CC) $(CFLAGS) $(INCFLAGS)  -c -o $@ $< 
	
$(LIKE2D_LIB) : 
	$(MAKE) -C $(LIKE2D_DIR) $@

$(LIKE3D_LIB) : 
	$(MAKE) -C $(LIKE3D_DIR) $@

$(SURF_LIB) : 
	$(MAKE) -C $(SURF_DIR) $@

clean: 
	rm -f *.o *.mod
	$(MAKE) -C $(LIKE2D_DIR) clean
	$(MAKE) -C $(LIKE3D_DIR) clean
	$(MAKE) -C $(SURF_DIR) clean

cleanall: clean 
	$(MAKE) -C $(LIKE2D_DIR) cleanall
	$(MAKE) -C $(LIKE3D_DIR) cleanall
	rm -f $(BIN_DIR)/MCTomo
