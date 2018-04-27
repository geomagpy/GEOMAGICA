#!/bin/bash
DBFILE="dB_2017-09-07T23:25:00.txt"
MODEL1D=39
MODEL2D="aerosq1000"
P=300

# Compile the FORTRAN code (-std=legacy may not always be needed):
gfortran -std=legacy -g thinsheet/nonuni_short_2d.f thinsheet/anomal.f thinsheet/seidel_new.f90 thinsheet/newmodel.f -o thinsheet.exe

# Run FORTRAN code:
./thinsheet.exe $DBFILE $P $MODEL1D ts_$MODEL2D.txt

# Compute GIC with output E-field file:
#python GIC_Model_Horton.py -d $DATESTR -p $P -l $MODEL1D -s $MODEL2D
