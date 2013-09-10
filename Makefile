all: v0_v1_fox v3 v4

v0_v1_fox:
	cd v0_v1_fox; make

v3:
	cd v3; make

v4: 
	cd v4; make

clean:
	cd v0_v1_fox; make clean
	cd v3; make clean
	cd v4; make clean
