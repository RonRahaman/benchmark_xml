all: fox v0 v1 v3 v4

fox:
	cd fox; make

v0:
	cd v0; make

v1: 
	cd v1; make

v3:
	cd v3; make

v4: 
	cd v4; make

clean:
	cd fox; make clean
	cd v0; make clean
	cd v1; make clean
	cd v3; make clean
	cd v4; make clean
