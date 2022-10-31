#!/bin/sh

DEBUG=0

if [ $DEBUG = 0 ] 
then

	gfortran -c "parametri.f90" "list_mod.f90" "physics.f90"
	gfortran -c "initialization.f90" "integratore.f90" "measures.f90" "Statistica.f90" 
	gfortran -c "Dinamica Molecolare.f90"
	gfortran "parametri.o" "list_mod.o" "physics.o" "initialization.o" "Statistica.o" "integratore.o" "measures.o" "Dinamica Molecolare.o" -o "Dinamica Molecolare"
	time ./"Dinamica Molecolare"

else

	gfortran -c "parametri.f90" -ggdb -O0
	gfortran -c "list_mod.f90" -ggdb -O0
	gfortran -c "physics.f90" -ggdb -O0
	
	gfortran -c "initialization.f90" -ggdb -O0
	gfortran -c "Statistica.f90" -ggdb -O0
	gfortran -c "integratore.f90" -ggdb -O0
	gfortran -c "measures.f90" -ggdb -O0
	
	gfortran -c "Dinamica Molecolare.f90" -ggdb -O0
	
	gfortran "parametri.o" "list_mod.o" "physics.o" "initialization.o" "Statistica.o"  "integratore.o" "measures.o" "Dinamica Molecolare.o" -o "Dinamica Molecolare"
	nemiver ./"Dinamica Molecolare"
	
fi
