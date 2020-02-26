# List of source code file
obj = GRT.o Rayleigh.o util.o Love.o C_interval.o C_interval_L.o SearchRayleigh.o SearchLove.o surfdisp96.o surfmodes.o

# rule for building surface modes code
$(SURF_LIB): $(obj)
	$(AR) $@ $^

# rule for building object file
%.o : %.f90
	$(F90) $(FFLAGS) -c -o $@ $<

%.o : %.f
	$(F90) $(F77FLAGS) -c -o $@ $<

.PHONY: clean cleanall

clean: 
	rm -f *.o *.mod

cleanall: clean
	rm -f *.a
