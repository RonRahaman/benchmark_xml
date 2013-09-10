VPATH = %.inc xml-fortran
VPATH = %.mod xml-fortran
F90 = gfortran
F90FLAGS = -cpp -fbacktrace -g

xml_fortran_objects = xml-fortran/iso_varying_string_1.3-F.o xml-fortran/read_xml_primitives.o xml-fortran/write_xml_primitives.o \
  xml-fortran/xmlparse.o

inc = read_from_buffer.inc read_xml_array.inc read_xml_scalar.inc read_xml_word.inc

objects = constants.o dict_header.o list_header.o geometry_header.o global.o error.o \
	  string.o xmlparse.o \
	  read_xml_primitives.o write_xml_primitives.o xml_data_geometry_t.o

src = constants.f90 dict_header.f90 list_header.f90 geometry_header.f90 global.f90 error.f90 \
	  string.f90 xmlparse.f90 \
	  read_xml_primitives.f90 write_xml_primitives.f90 geometry_t_new.f90 \
	  geometry_t_old.f90 xml_interface.f90


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
	rm *.o *.mod main benchmark_xml.txt geometry.xml

FORCE:
