VPATH = %.inc 
VPATH = %.mod xml/lib
F90 = gfortran
FC = gfortran
FFLAGS = -cpp -fbacktrace
F90FLAGS = -cpp -fbacktrace

# For Fox DOM
xml_lib = -Lxml/lib -lxml

inc = read_from_buffer.inc read_xml_array.inc read_xml_scalar.inc read_xml_word.inc

objects = constants.o dict_header.o list_header.o geometry_header.o global.o error.o \
	  string.o xml_interface.o

%.o : %.f90
	$(FC) -c $(FFLAGS) $(xml_lib) -Ixml/include $< -o $@

main: main.f90 xml $(objects)
	$(FC) $(FFLAGS) $(objects) $(xml_lib) -Ixml/include $< -o $@

xml: FORCE
	cd xml; make F90=$(F90) F90FLAGS="$(F90FLAGS)"

clean:
	cd xml; make clean
	rm *.o *.mod main benchmark_xml.txt geometry.xml

FORCE:
