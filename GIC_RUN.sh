#!/bin/bash
DATESTR="2017-09-07T23:25:00"
MODEL1D=39
MODEL2D="aerosq1000"
P=300

# Compile the FORTRAN code (-std=legacy may not always be needed):
gfortran -std=legacy -g thinsheet/nonuni_short_2d.f thinsheet/anomal.f thinsheet/seidel_new.f90 thinsheet/newmodel.f -o thinsheet.exe

# Run FORTRAN code:
DBFILE="dB_${DATESTR}.txt"
./thinsheet.exe $DBFILE $P $MODEL1D ts_$MODEL2D.txt

# Compute GIC with output E-field file:
EFILE="Efiles/E_${MODEL1D}_${DATESTR}.txt"
python GIC_Model_Horton.py -e $EFILE
