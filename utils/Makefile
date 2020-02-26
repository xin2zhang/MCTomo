# List of source file
OBJ = m_logger.o m_exception.o mt19937.o c_utils.o utils.o kdtree2.o

# rule for building foward modeling code
$(UTILS_LIB): $(OBJ)
	$(AR) $@ $^

# rule for building object file
%.o: %.f90
	$(F90) $(FFLAGS) -c $< -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: clean cleanall

clean:
	rm -f *.o *.mod

cleanall: clean
	rm -f *.a
