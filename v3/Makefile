VPATH = %.inc xml-fortran
VPATH = %.mod xml-fortran
F90 = gfortran
F90FLAGS = -cpp -fbacktrace -g

xml_fortran_objects = xml-fortran/iso_varying_string_1.3-F.o xml-fortran/read_xml_primitives.o xml-fortran/write_xml_primitives.o \
  xml-fortran/xmlparse.o

inc = read_from_buffer.inc read_xml_array.inc read_xml_scalar.inc read_xml_word.inc

objects = constants.o dict_header.o list_header.o geometry_header.o tally_header.o \
					mesh_header.o global.o tally_initialize.o error.o string.o xmlparse.o \
	  read_xml_primitives.o write_xml_primitives.o xml_data_geometry_t.o xml_data_tallies_t.o

main: $(xml_fortran_objects) $(objects) main.f90 
	$(F90) $(F90FLAGS) -I xml-fortran -J xml-fortran $^ -o $@

xml_data_tallies_t.o: tallies_t.f90 xml-fortran/xmlreader
	$(F90) $(F90FLAGS) -I xml-fortran -c $< -o $@

tallies_t.f90: xml-fortran/xmlreader tallies_t.xml
	$< tallies_t


main: $(xml_fortran_objects) $(objects) main.f90 
	$(F90) $(F90FLAGS) -I xml-fortran -J xml-fortran $^ -o $@

xml_data_geometry_t.o: geometry_t.f90 xml-fortran/xmlreader
	$(F90) $(F90FLAGS) -I xml-fortran -c $< -o $@

geometry_t.f90: xml-fortran/xmlreader geometry_t.xml
	$< geometry_t

xml-fortran/xmlreader $(xml_fortran_objects): xml-fortran/Makefile
	cd xml-fortran; make F90=$(F90) F90FLAGS="$(F90FLAGS)"


%.o : %.f90
	$(F90) -c $(F90FLAGS) $< -o $@

clean:
	cd xml-fortran; make clean
	rm *.o *.mod main *_t.f90

FORCE:
