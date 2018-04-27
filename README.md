# Geomagnetically induced currents in Austria (GEOMAGICA)

Project GEOMAGICA has been an FFG-funded project conducted in Austria with contributions from international partners. The aim was to model the nature of geomagnetically induced currents (GICs) in Austria and provide predictions or warnings of areas in the power grid particularly at risk. Details on the Austrian power grid cannot be provided, so the modelling code has instead been applied here to an example grid taken from the literature. Geoelectric fields induced by geomagnetic variations are modelled using a thinsheet approximation that includes information on the lateral surface and layered subsurface resistivity. From a given geoelectric field, the GICs in the network are computed using the Lehtinen-Pirola method.

## Basic usage

Requires some standard scientific Python packages and a FORTRAN compiler.

To run the most basic example using the Horton (2012) example grid and a 1 V/km geoelectric field in the N and E components, then print the results:

```
python GIC_Model_Horton.py
```

Instructions on how to run the scripts for more complex models are included in GUN_RUN.sh. The only non-intuitive part is the compilation of the thinsheet FORTRAN code, which may need an option for legacy code:

```
gfortran -std=legacy -g thinsheet/nonuni_short_2d.f thinsheet/anomal.f thinsheet/seidel_new.f90 thinsheet/newmodel.f -o thinsheet.exe
```

The thinsheet code should be run with input files specified:

```
./thinsheet.exe <DBDT-FILE> <PERIOD> <1DMODEL-CODE> <2DMODEL-FILE>
OR
./thinsheet.exe dB_2017-09-07T23:25:00.txt 300 39 ts_aerosq1000.txt
```

## Modelling GIC - GIC_Model_Horton.py

This script uses the original method devised by Lehtinen and Pirjola (see 1985 paper) for computing the amount of GIC flowing between a grounded conductive network (such as a power grid) and the Earth. 

In this example, the example network provided in the Horton et al. (2012) paper is used for testing purposes.

### Input

* **Efiles/E_\<1DMODEL\>_\<DATETIME\>.txt** - Coming soon.

## Thinsheet code - thinsheet.exe

This program takes a file describing the geomagnetic variations in a region, a 1D subsurface resistivity layer model and a 2D surface thinsheet model to compute the induced geoelectric field in the region. It is based off of the formulae derived in the Vasseur and Weidelt (1977) paper and this particular code was developed by P. Weidelt and later adapted by British Geological Survey researchers for modelling regions in Britain. This has been adapted to work in Austria.

### Code

* **nonuni_short_2d.f** is where the main computation happens and contains most of the algorithm. Geomagnetic variation data is read directly from the folder dBfiles/dB_\<DATETIME\>.txt, and output geoelectric field is written to folder Efiles/E_\<1DMODEL\>_\<DATETIME\>.txt.
* **anomal.f** reads in the thinsheet surface conductivity model stored under condmodels/ts_\<2DMODEL\>.txt and returns it to the main code.
* **newmodel.f** reads in the layered subsurface resistivity model stored under condmodels/layers_\<1DMODEL\>.txt and returns it to the main code.
* **seidel_new.f90** determines the solution of the integral equations (see Vasseur and Weidelt 1977 paper) using the Gauss-Seidel iteration method. (Old version **seidel.f** is also included.)

### Setting parameters

The parameters required for determining how thinsheet.exe executes are listed in the file *thinsheet_Parameters.txt*. Below is a description of each parameter:

```
 10000.0                # DX, dimension of square cells in the thinsheet in m
 42                     # NX, number of cells in the x (North-South) direction
 65                     # NY, number of cells in the y (East-West) direction
 5.37                   # X_INC, cell spacing in x direction in arcmin
 8.08                   # Y_INC, cell spacing in y direction in arcmin
 0.0                    # ZA, depth of thinsheet, 0 for sheet at surface
 0.0                    # ZD, depth of data level, 0 for surface
 49.81                  # NBOUND, northern boundary of data in lat
 46.05                  # SBOUND, southern boundary of data in lat
 17.85                  # EBOUND, eastern boundary of data in lon
 9.1                    # WBOUND, western boundary of data in lon
 6                      # NITER, iteration count, usually converges after 5 in this case
 10                     # NL, number of layers in resistivity layers
 20.0                   # TAUN, conductance of thinsheet outside the anomalous domain
```

Note that the last two parameters are read in the *newmodel.f* subroutine. The others are read in the main code *nonuni_short_2d.f*.

### Input files

* **dBfiles/dB_\<DATETIME\>.txt** - should contain *NX x NY* lines for the geomagnetic variations (in the northward and eastward directions) over each cell. These variations should be in nT per period, with the period *P* defined when calling *thinsheet.exe*. In the example file *dB_2017-09-07T23:25:00.txt*, variations from September 7 2017 were used as these were larger than usual. The original data was filtered down to 300 s to calculate the dB/dt value over the time window. The result is 24.5 nT/300s in X and 31.9 nT/300s in Y, which is assumed to be constant across the whole region. The grid used is *42 x 65*, roughly 10 x 10 km, resulting in 2730 cells.
* **condmodels/ts_\<2DMODEL\>.txt** - this contains the "thinsheet" in conductance S. The example provided, *ts_aerosq1000.txt*, is the inhomogeneous model based off aerogeomagnetic measurements conducted in Austria of the surface conductivity projected onto a hydrogeological map of Austria (see Schattauer et al. 2017). The same shape constraints apply as with the dB/dt values, and the grid of *NX x NY* cells must be provided in the same format.
* **condmodels/layers_\<1DMODEL\>.txt** - contains 1D information on the resistivity of the Earth in layers. The values for Austria in *layers_austria.txt* were taken from the European Rho Model EURHOM (see Adam et al. 2012). For additional tests, a model of Quebec is also included in *layers_quebec.txt*. Here, the layers are defined as *D* and *RHO*, so the first column (after the layer numbers) is thickness in m followed by resistivity in Ohm. EURHOM model 39 (third set of columns), for example, has a top layer thickness of 55 km with a resistivity of 1000 Ohm. Below that is a layer with thickness 45 km and a resistivity of 300 Ohm, and so on...

### Output files

* **Efiles/E_\<1DMODEL\>_\<DATETIME\>.txt** - the output is written in the same grid format as the input dB/dt and thinsheet files, assuming that the cell numbers and cell spacings were defined correctly in the thinsheet parameters.

## Resources

1. Adam, A., E. Pracser, and V. Wesztergom (2012), Estimation of the electric resistivity distribution (EURHOM) in the European lithosphere in the frame of the EURISGIC WP2 project, Acta Geodaetica et Geophysica Hungarica, 47 (4), 377–387, doi:10.1556/AGeod.47.2012.4.1.

2. Horton, R., D. Boteler, T. J. Overbye, R. Pirjola, and R. C. Dugan (2012), A test case for the calculation of geomagnetically induced currents, IEEE Transactions on Power Delivery, 27 (4), 2368–2373.

3. Lehtinen, M., and R. Pirjola (1985), Currents produced in earthed conductor networks by geomagnetically-induced electric fields, Annales Geophysicae, 3, 479–484.

4. Schattauer, I., A. Römer, R. Bailey, R. Leonhardt, K. Motschka, R. Supper, and A. Schiller (2017), Lateral conductivity variations within Austria and its surroundings by extrapolating airborne electromagnetic data, 2nd European Airborne Electromagnetics Conference 2017, Held at Near Surface Geoscience Conference and Exhibition 2017, pp. 136–140.

5. Vasseur, G., and P. Weidelt (1977), Bimodal electromagnetic induction in non-uniform thin sheets with an application to the northern Pyrenean induction anomaly, Geophysical Journal International, 51 (3), 669–690.

## Authors

* **P. Weidelt** - developer of the thinsheet approximation algorithm
* **C. Beggan, K. Turnbull, A. Mckay** (BGS) - adapted and provided the original code
* **R. Bailey** - adapted said code to region of Austria

<!---## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details--->



