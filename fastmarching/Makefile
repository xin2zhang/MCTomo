# List of source file
#CC = g++
#FMM3D_LIB = $(patsubst %,lib%.a,fmm3d)
#CFLAGS = -std=c++11 -g -O3
AR = ar r
OBJ = fastMarching.o console/console.o ndgridmap/cell.o fmm/fmdata/fmcell.o

INCFLAGS += -Iconsole
INCFLAGS += -Indgridmap
INCFLAGS += -Ifmm
INCFLAGS += -Ifmm/fmdata
INCFLAGS += -Iio
INCFLAGS += -Igradientdescent

# rule for building foward modeling code
$(LIKE3D_LIB): $(OBJ)
	$(AR) $@ $^

# rule for building object file
%.o: %.cpp
	$(CC) $(CFLAGS) $(INCFLAGS) -c $< -o $@

.PHONY: clean cleanall

clean:
	rm -f console/*.o ndgridmap/*.o fmm/fmdata/*.o *.o *.mod

cleanall: clean
	rm -f *.a
