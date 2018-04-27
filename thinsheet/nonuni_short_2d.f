      PROGRAM nonuni_short_2d

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
C                      'S H O R T'                                  C
C     PROGRAM FOR MODELLING LATERALLY NON-UNIFORM THIN SHEETS       C
C     AUTHOR: P. WEIDELT, INSTITUT FUER GEOPHYSIK UND METEOROLOGIE, C
C     ------ MENDELSSOHNSTRASSE 3, D-38106 BRAUNSCHWEIG, FRG        C
C                                                                   C
C     PARTICULAR FEATURES OF THE PROGRAM:                           C
C     ===================================                           C
C     1) INTEGRAL EQUATION TECHNIQUE (cf. VASSEUR & WEIDELT, GEO-   C
C        PHYS. J. R. ASTR. SOC, v.51, p.669-690, 1977)              C
C     2) THE THIN SHEET MAY BE POSITIONED AT ANY DEPTH ZA IN A LAY- C
C        ERED HALFSPACE.                                            C
C     3) THE DATA CAN BE COMPUTED AT ANY DATA LEVEL ZD IN THE HALF- C
C        SPACE.                                                     C
C     4) INDUCTION ARROWS CAN BE CALCULATED FOR EACH POINT IN THE   C
C        ANOMALOUS DOMAIN.                                          C
C     5) ALL QUANTITIES ARE GIVEN IN SI-UNITS (e.g. B IN T, LENGTH  C
C        IN m, FREQUENCY IN Hz, CONDUCTANCE IN S, etc.              C
C     6) ON ALL DISPLAYS THE X-DIRECTION SHOWS 'UPWARDS', THE Y-    C
C        AXIS POINTS TO THE RIGHT AND Z SHOWS INTO THE PAPER  (OR   C
C        SCREEN).                                                   C
C     7) 2D ARRAYS IN THE (X,Y)-PLANE ARE GENERALLY STORED AS 1D    C
C        ARRAYS, WHERE THE NUMBERING PROCEEDS first ALONG X-DIREC-  C
C        TION.                                -----                 C
C                                                                   C
C     INPUT PARAMETERS:                                             C
C     =================                                             C
C     NORMAL CONDUCTIVITY STRUCTURE (without THIN SHEET)            C
C     --------------------------------------------------            C
C     NL:      NUMBER OF LAYERS (NL=1 FOR UNIFORM HALF-SPACE)       C
C          ****RESTRICTION: NL.LE.10                                C
C     RHO:     1D ARRAY OF DIMENSION NL FOR LAYER RESISTIVITIES     C 
C     D:       1D ARRAY OF DIMENSION NL FOR LAYER THICKNESSES,      C
C              D(NL) ARBITRARY                                      C
C                                                                   C
C     THIN SHEET                                                    C
C     ----------                                                    C
C     ZA:      DEPTH OF THIN SHEET (.GE.0, ARBITRARY ELSE, cf. RHONEW
C     TAUN:    CONDUCTANCE OF THIN SHEET OUTSIDE THE RECTANGULAR    C
C              ANOMALOUS DOMAIN                                     C
C     DX:      DIMENSION OF SQUARE CELLS IN THE ANOMALOUS DOMAIN    C
C     NX:      NUMBER OF CELLS IN X-DIRECTION                       C
C     NY:      NUMBER OF CELLS IN Y-DIRECTION                       C
C     TAU:     1D ARRAY OF DIMENSION NX*NY CONTAINING THE CONDUC-   
C              TANCES IN THE ANOMALOUS DOMAIN, cf. SUBROUTINE ANO-  C
C              MAL.                                                 C
C                                                                   C
C     DEFINITION OF DATA DOMAIN                                     C
C     -------------------------                                     C
C     THE FIELD IS CALCULATED IN THE CENTER OF SQUARE CELLS OF DI-  C
C     MENSION DX (cf. ABOVE) ON A RECTANGULAR GRID. THIS GRID CAN   C
C     BE SHIFTED WITH RESPECT TO THE ANOMALOUS DOMAIN BY INTEGER    C
C     GRID WIDTHS IN X- AND Y-DIRECTION.                            C
C     ZD:      DEPTH OF DATA LEVEL (SURFACE: ZD=0)                  C
C     NX2:     NUMBER OF SQUARE CELLS IN X-DIRECTION                C
C     NY2:     NUMBER OF SQUARE CELLS IN Y-DIRECTION                C
C     NX0:     X-ORIGIN OF THE DATA DOMAIN IN UNITS OF DX           C
C     NY0:     Y-ORIGIN OF THE DATA DOMAIN IN UNITS OF DX           C
C ****IF ZD COINCIDES WITH ZA OR A LAYER BOUNDARY, THE FIELD VALUES C
C     DISPLAYED REFER TO ZD-0.                                      C
C                                                                   C
C     FREQUENCY                                                     C
C     ---------                                                     C
C     F:       FREQUENCY                                            C
C                                                                   C
C     SOURCE PARAMETERS                                             C
C     -----------------                                             C
C     B0:      NORMAL QUASI-UNIFORM MAGNETIC FIELD (T)              C
C     THETA:   AZIMUTH OF B0 (deg), COUNTED POSITIVE FROM X TO Y.   C
C                                                                   C
C     ITERATION PARAMETERS                                          C
C     --------------------                                          C
C     NITER:   MAXIMUM NUMBER OF ITERATIONS FOR GAUSS-SEIDEL        C
C     ITERD:   FOR IWS=1 OR IWK=1 INTERMEDIATE RESULTS ARE PRINTED  C
C              AFTER ITERD ITERATIONS.                              C
C     SOR:     SUCCESSIVE OVERRELAXATION PARAMETER (GAUSS-SEIDEL)   C
C     EPS:     ITERATION STOPS, IF RELATIVE DEVIATION BETWEEN SUC-  C
C              CESSIVE ITERATIONS IS LESS THAN THIS THRESHOLD.      C
C                                                                   C
C     PRINTING CONTROL                                              C
C     ----------------                                              C
C     NH:      NUMBER OF COMPLEX FIELD VALUES IN HORIZONTAL LINE    C
C              (NH.LE.3 FOR SCREEN DISPLAY, NH.LE.5 FOR ORDINARY    C
C              PRINTER DISPLAY)                                     C
C     IWG1:    =1 FOR DISPLAY OF INTEGRAL EQUATION KENELS (GREEN1)  C
C     IWG2:    =1 FOR DISPLAY OF DATA KERNELS (GREEN2)              C
C     IWN1:    =1 FOR DISPLAY OF NORMAL E-FIELD AT ZA (NORMAL1)     C
C     IWN2:    =1 FOR DISPLAY OF NORMAL FIELD AT ZD (NORMAL2)       C
C     IWS:     =1 FOR DISPLAY OF RESULTS FROM GAUSS-SEIDEL ITERATIONC
C     IWD:     =1 FOR DISPLAY OF RESULTS FROM FRECHET DERIVATIVES   C
C     IWF:     =1 FOR DISPLAY OF MAGNETIC FIELD COMPONENTS AT ZD    C
C     IWZ:     =1 FOR DISPLAY OF INDUCTION ARROWS                   C
C     IWQ:     =1 FOR DISPLAY OF SOURCE INFORMATION                 C
C                                                                   C
C     MINIMUM DIMENSIONS REQUIRED:                                  C
C     ============================                                  C
C     DEFINITIONS:                                                  C
C     N1=MAX(NX,NY), N11=N1*N1                                      C
C     N2=MAX(NX2+NX0,NY2+NY0,NX-NX0,NY-NY0), N22=N2*N2              C
C     NXY=NX*NY, NXY2=NX2*NY2                                       C
C                                                                   C
C     COMPLEX GXX(N11),GXY(N11),ENX(NXY),ENY(NXY),ESX(NXY),ESY(NXY) C
C     COMPLEX BXX(N22),BXY(N22),BXZ(N22)                            C
C     COMPLEX BX(NXY2),BY(NXY2),BZ(NXY2)                            C
C     COMPLEX ZZ(6,NXY2),ES(4,NXY)                                  C
C     DIMENSION TAU(NXY),TAUA(NXY)                                  C
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
!       PARAMETER (NXD=NX,NYD=NY,N1D=NYD,N11D=N1D*N1D)
!       PARAMETER (NX2D=200,NY2D=300,N2D=400,N22D=N2D*N2D)
!       PARAMETER (NXYD=NXD*NYD,NXY2D=NX2D*NY2D)
      INTEGER NXD, NYD, N1D, N11D, NXYD, NXY2D
      PARAMETER (NX2D=200,NY2D=300,N2D=400,N22D=N2D*N2D)
!       COMPLEX GXX(N11D),GXY(N11D)
!       COMPLEX ENX(NXYD),ENY(NXYD),ESX(NXYD),ESY(NXYD)
!       COMPLEX BXX(N22D),BXY(N22D),BXZ(N22D)
!       COMPLEX BX(NXY2D),BY(NXY2D),BZ(NXY2D)
!       COMPLEX ZZ(6,NXY2D),ES(4,NXYD),
!      2     EXSUM(43460,2),EYSUM(43460,2)
      COMPLEX, ALLOCATABLE :: GXX(:),GXY(:)
      COMPLEX, ALLOCATABLE :: ENX(:),ENY(:),ESX(:),ESY(:)
      COMPLEX, ALLOCATABLE :: BXX(:),BXY(:),BXZ(:)
      COMPLEX, ALLOCATABLE :: BX(:),BY(:),BZ(:)
      COMPLEX, ALLOCATABLE :: ZZ(:,:),ES(:,:)
      COMPLEX EXSUM(43460,2),EYSUM(43460,2)

!      complex trial3(NXD,NYD,NXD,NYD)

      REAL, ALLOCATABLE :: dBx(:)
      REAL west,east,south,north,inc,RLAT,RLON,
     1     BGRADX(43460),BGRADY(43460),
     1     EXSUMR,EXSUMI,EYSUMR,EYSUMI,XIN,YIN

!       DIMENSION TAU(NXYD),TAUA(NXYD)
      REAL, ALLOCATABLE :: TAU(:),TAUA(:)
      COMMON /RESIST/NL,RHO(20),D(20),NL1,RHO1(22),D1(22),
     1     NLA,NLD,TAUN
      COMMON /PARAM/F,U,DX,NITER,ITERD,EPS,SOR,NH,ZA,ZD
      COMMON /HANKEL/NC,NC0,HR(2,100)
      COMMON /NUMBER/NX,NY,NX2,NY2,NX0,NY0,N1,N2
      COMMON /SOURCE/B0,THETA,BGRAD(43460)

      DATA IWG1,IWG2,IWN1,IWN2,IWS,IWD,IWF,IWZ,IWQ/0,0,0,0,1,0,1,0,1/

      CHARACTER*2 modstr
      CHARACTER*6 filex
      CHARACTER*12 cdffile
      CHARACTER*14 efilex,efiley
      CHARACTER*32 filein, Fstring, cmodelfile
      CHARACTER*19 datestring
      CHARACTER*46 Efieldout
      CHARACTER*32 arg
c
c      OPEN (6,FILE='SHORT.RES')
c      open(100,FILE='efield.dat')
     
!     READ IN PARAMETERS (RLB added 2016):
!     ------------------------------------
      open(33,FILE='thinsheet_Parameters.txt',STATUS='OLD')
      read(33,*) DX
      read(33,*) NX
      read(33,*) NY
      read(33,*) x_inc
      read(33,*) y_inc
      read(33,*) ZA
      read(33,*) ZD
      read(33,*) north
      read(33,*) south
      read(33,*) east
      read(33,*) west
      read(33,*) NITER
      
      write(*,*) NITER
      
      NXD=NX
      NYD=NY
      N1D=NYD
      N11D=N1D*N1D
      NXYD=NXD*NYD
      NXY2D=NX2D*NY2D
      
      allocate( GXX(N11D) )
      allocate( GXY(N11D) )
      allocate( ENX(NXYD) )
      allocate( ENY(NXYD) )
      allocate( ESX(NXYD) )
      allocate( ESY(NXYD) )
      allocate( BXX(N22D) )
      allocate( BXY(N22D) )
      allocate( BXZ(N22D) )
      allocate( BX(NXY2D) )
      allocate( BY(NXY2D) )
      allocate( BZ(NXY2D) )
      allocate( ZZ(6,NXY2D) )
      allocate( ES(4,NXYD) )     
      allocate( dBx(NXYD) )
      allocate( TAU(NXYD) )
      allocate( TAUA(NXYD) )

C     DEFINITION OF ANOMALOUS DOMAIN AND DATA DOMAIN:
C     -----------------------------------------------
      
      NX0=-2    ! for positioning field arrows we need a slight offset. Not for E-field
      NY0=-2

      NX2=NX+2*ABS(NX0)
      NY2=NY+2*ABS(NY0)

!       ZA=0.       ! thin sheet depth is 0.0
!       ZD=0.

      N1=MAX0(NX,NY)
      N2=MAX0(NX2+NX0,NY2+NY0,NX-NX0,NY-NY0)
      NXY=NX*NY
      NXY2=NX2*NY2


c     DEFINITIONS FOR GRID FILES:
c     ---------------------------

      x_inc = (1./60.)*x_inc      ! spacing of GMT grid in arcmin
      y_inc = (1./60.)*y_inc      ! spacing of GMT grid in arcmin
! 

c **TEST DATA**
C
C Define THETA (degrees east from true north), and driving frequency F (Hz).
C These numbers need to be read in together with the BGRAD field variations in nT.
C
c      THETA= 0.0 ! Specify 0.0=Northward perturbation; 90.0=eastward perturbation, now set within program
c      F=1./120.0  ! 2 minutes period assumed for external magnetic field
c      B0=1.0E-09  ! Value for scaling B-field perturbation, here converting nT to T
C
c     Define test functions of a north->south and west->east magnetic field gradient 
C
c      do i =1,NX  ! south->north gradient between 0 and 1 nT
c         dBx(i) = 1. * i/NX
c      end do
c      open(100,file='testB2.txt') ! save values to new file secs.txt, for checking      
c      do K=1,NY       ! EAST-WEST &
c          DO I=1,NX   ! NORTH-SOUTH gradient field
c            IK=(K-1)*NX+I
c            BGRADX(IK) = 0.+100.0*dBx(I)                                  ! Assume *1000nT field variation
c     *                   1.0*(1.0+dBx(I))                          ! Bx gradient in S->N   direction
c     1                  0.5*(dBx(I)+(1.0-float(K-1)/float(NY-1)))  ! Bx gradient in SE->NW direction
c            BGRADY(IK) = 0.+100.0*dBx(I)                                  ! Assume *1000nT field variation
c     *                         (dBx(I))                           ! By gradient in S->N   direction
c     1                    0.5*(dBx(I) + (float(K-1)/float(NY-1))) ! By gradient in SW->NE direction
c            print *, BGRADX(IK), BGRADY(IK)
c            write(100,'(4E16.7)')
c     1      (south+(I-1)*x_inc),(west+(K-1)*y_inc),
c     2       BGRADX(IK),BGRADY(IK)
c          END DO
c      END DO
c      close(100)
c
c **END TEST DATA**


C Here is a FORTRAN read statement, for reading from an external array "secs.txt", 
C at 5 arc-min spacing in both X (north) and Y (east), in the geographic array bounded by 
C the area "south/north/west/east" defined above, i.e. using same SW origin of coordinates.
C
C There needs to be 2 loops within this programme. 
C
C First with the grid of just the components in the North direction (THETA=0), to get
C an Ex, Ey pair. Then a second pass with the grid of B components just in the East direction (THETA=90), 
C to get a second pair of Ex, Ey. 
C
C Then the program sums these 2 pairs of Ex, Ey to get the full E-field 
C vector appropriate for your initial grid in X and Y directions of field perturbations BGRAD,i.e.
C Ex(i,j)=Ex(i,j,theta=0)+Ex(i,j,theta=90), and same for Ey. NB, Ex, Ey are complex fields.
C
C Here the BGRAD (X,Y) field is read in, assuming it is in nT.
C Therefore to get correct units, scale by a fixed B0=1.0E-09.
C
c      write(*,*) ' Period (seconds)?'
c      read(*,*) F
c      write(*,*) ' Input filename (xxxxxxxx.out only)?'
c      read(*,'(a)') filein

c     PROVIDE INPUT FILES WHEN CALLING (RLB added 2016):
c     --------------------------------------------------

      call get_command_argument(1, filein)
      call get_command_argument(2, Fstring)
      read(Fstring,*)P
      call get_command_argument(3, modstr)
      read(modstr,*)layermodel
      call get_command_argument(4, cmodelfile)
      datestring = filein(4:24)
      write(*,*) datestring

      open(100,file='dBfiles/'//filein,status='old')  ! i.e. file already exists
C
      F=1.0/P   ! Frequency 1/P
C      F=1.0/360.0 ! Or whatever you want to set it to. THETA is now implicit and B0 is fixed=1.E-09
      B0=1.0E-09  ! Value for scaling B-field perturbation, here converting nT to T

      do K=1,NY         ! loop over EAST-WEST
          DO I=1,NX,1   ! loop over NORTH-SOUTH first
            READ(100,'(4E16.7)') RLAT,RLON,XIN,YIN              ! BGRAD is assumed to be read in as nT
c            write(*,*) K, I, RLAT,RLON,XIN,YIN 
C            READ(100,*) RLAT,RLON,BGRADX(IK),BGRADY(IK) ! BGRAD is assumed to be read in as nT
            IF(I.GE.1) THEN
             IK=(K-1)*NX+I
             BGRADX(IK)=XIN
             BGRADY(IK)=YIN
             if(((i.eq.NX).and.(k.eq.NY)).or.
     1          ((i.eq.1).and.(k.eq.1))) 
     1       write(*,*) I,K,IK,B0,F,RLAT,RLON,BGRADX(IK),BGRADY(IK)
c             write(*,*) I,K,IK,B0,F,RLAT,RLON,BGRADX(IK),BGRADY(IK)
c             print *, BGRADX(IK), BGRADY(IK)
             BGRADX(IK)=B0*BGRADX(IK)                   ! BGRADX is now in T.
             BGRADY(IK)=B0*BGRADY(IK)                   ! BGRADY is now in T.
            END IF
          END DO
      END DO
      close(100)
C
C     TEST OF DIMENSIONS:
C     -------------------
      IERROR=0
      IF (NXY.LE.NXYD) GOTO 10
      WRITE (*,*) "IERROR #1", NXYD,NXY
      IERROR=IERROR+1
   10 IF (NXY2.LE.NXY2D) GOTO 20
      WRITE (*,*) "IERROR #2", NXY2D,NXY2
      IERROR=IERROR+1
   20 IF (N1.LE.N1D) GOTO 30
      WRITE (*,*) "IERROR #3", N1D,N1
      IERROR=IERROR+1
   30 IF (N2.LE.N2D) GOTO 40
      WRITE (*,*) "IERROR #4", N2D,N2
      IERROR=IERROR+1
   40 IF (IERROR.GT.0) STOP


c     AUSTRIAN MODELS:
c     ----------------
c      call NEWMODEL(3)   !  "newmodel.f": EURHOM model for Vienna basin (easternmost Austria) region
c      call NEWMODEL(16)   !  "newmodel.f": EURHOM model for sub-alpine (northernmost Austria) region
c      call NEWMODEL(39)   !  "newmodel.f": EURHOM model for sub-alpine (central/Vienna) region
c      call NEWMODEL(55)   !  "newmodel.f": EURHOM model for alpine region
      call NEWMODEL(layermodel)

      
C     ITERATION PARAMETERS:
C     ---------------------
!       NITER=3
      ITERD=50
      SOR=1.0
      EPS=1.E-2


C     PRINTING PARAMETER:
C     -------------------
      NH=3


C     PREPARATORY SOURCE-FIELD INDEPENDENT CALCULATIONS:
C     --------------------------------------------------
      CALL RHONEW

C ****RHONEW HAS TO BE CALLED AFTER EACH CHANGE OF ZA AND/OR ZD
c      CALL PRINTRHO
      
      CALL ANOMAL(TAUN,TAU,TAUA,cmodelfile) ! this is where 2D conductivity model is specified - see "anomal.f"

C     CREATE netCDF file
      
c      CALL netcdf_create(cdffile)
c      CALL write_params(cdffile)

      CALL FILTER
      CALL PRINTN
      CALL GREEN1(GXX,GXY,IWG1)   
      CALL GREEN2(BXX,BXY,BXZ,IWG2)


c     print greens functions as a function of distance from source cell
c      open (20,FILE='/scratch/allan/thin/modelruns/gyx.dat')
c      DO i =1,NY
c         WRITE(20,*) i,REAL(GXY(i)),REAL(GXY(NX+i)),REAL(GXY((2*NX)+i))
c      END DO
c      close(20)   
c         PAUSE


C     MAGNETIC FIELD AND INDUCTION ARROWS (TIPPER):
C     --------------------------------------------
C      B0=1. ! test value

      DO 100 IPOL=1,2  ! IPOL=1=>Theta=0; IPOL=2=>THETA=90

      IF(IPOL.EQ.1) THEN
       THETA= 0.0
       do K=1,NY       ! EAST-WEST 
          DO I=1,NX    ! NORTH-SOUTH
            IK=(K-1)*NX+I
            BGRAD(IK)=BGRADX(IK)
          end do
       end do
      ELSE
       THETA=90.0
       do K=1,NY       ! EAST-WEST 
          DO I=1,NX    ! NORTH-SOUTH
            IK=(K-1)*NX+I
            BGRAD(IK)=BGRADY(IK)
          end do
       end do
      END IF

      IF (IWQ.EQ.1) CALL PRINTQ
      
      if (IPOL.eq.1) write(*,'(/,a)') ' Bx-> Ex and Ey'
      if (IPOL.eq.2) write(*,'(/,a,/)') ' By-> Ex and Ey'

      write(*,*)' MAIN PROG: into NORMAL1 '
      CALL NORMAL1(ENX,ENY,IWN1) ! This is used for E field

      write(*,*)' MAIN PROG: into SEIDEL with 1/F, NMO = ',1./f,NMO
      CALL SEIDEL_new(TAUA,GXX,GXY,ENX,ENY,ESX,ESY,IWS,NXD,NYD)
c      CALL SEIDEL(TAUA,GXX,GXY,ENX,ENY,ESX,ESY,IWS,NXD,NYD)

c      CALL cgsolve(TAUA,GXX,GXY,ENX,ENY,ESX,ESY,IWS,NXD,NYD)

cc      CALL write_params(ESX,ESY,IPOL)

      write(*,*)' MAIN PROG: into NORMAL2'
      CALL NORMAL2(BX,BY,BZ,IWN2) ! This is also used for E field
      
      write(*,*)' MAIN PROG: into BE'
c      CALL BE(TAUA,BXX,BXY,BXZ,BX,BY,BZ,ESX,ESY,IWF,IPOL,cdffile) 
      CALL BE(TAUA,BXX,BXY,BXZ,BX,BY,BZ,ESX,ESY,IWF,IPOL) 
c
c     added IPOL in call for writing netcdf files

C     STORE MAGNETIC FIELD ZZ:
      DO 80 IK=1,NXY2
      DO 80 IC=1,3
      IDX=3*(IPOL-1)+IC
      IF (IC.EQ.1) ZZ(IDX,IK)=BX(IK)
      IF (IC.EQ.2) ZZ(IDX,IK)=BY(IK)
   80 IF (IC.EQ.3) ZZ(IDX,IK)=BZ(IK)

C     STORE SHEET ELECTRIC FIELD IN ARRAY ES:
      DO 90 IK=1,NXY
      DO 90 IC=1,2
      IDX=2*(IPOL-1)+IC
      IF (IC.EQ.1) ES(IDX,IK)=ESX(IK)
   90 IF (IC.EQ.2) ES(IDX,IK)=ESY(IK)

C
C Save the electric field vectors.

      do K=1,NY       ! EAST-WEST 
          DO I=1,NX   ! NORTH-SOUTH
            IK=(K-1)*NX+I
            EXSUM(IK,IPOL)=ESX(IK)
            EYSUM(IK,IPOL)=ESY(IK)
          end do
      end do
      
  100 CONTINUE


C
C *NEW* write statement to O/P E-field as simple X,Y arrays. NB: Units output are V/km.
C Output real and imag parts of E-field in X and Y direction, as well as the magnitude, as check
C
c      Efieldout='E_XXXXXXXxxxx_XXX_7.txt'
c      write(Efieldout(3:13),'(a11)') filein(1:11)
c      write(Efieldout(15:17),'(I3.3)') NINT(1.0/F)
      Efieldout='Efiles/E_'//modstr//'_'//datestring//'.txt'
      write(*,*) ' Write combined Ex, Ey fields to ',Efieldout
c      open(150,file='boutput/'//Efieldout)
      open(150,file=Efieldout)
      do K=1,NY       ! EAST-WEST 
          DO I=1,NX   ! NORTH-SOUTH
            IK=(K-1)*NX+I
            EXSUMR=REAL(EXSUM(IK,1))+REAL(EXSUM(IK,2))
            EXSUMI=AIMAG(EXSUM(IK,1))+AIMAG(EXSUM(IK,2))
            EYSUMR=REAL(EYSUM(IK,1))+REAL(EYSUM(IK,2))
            EYSUMI=AIMAG(EYSUM(IK,1))+AIMAG(EYSUM(IK,2))   
            WRITE(150,'(8F11.4)') 
     1       (south+(I-1)*x_inc),(west+(K-1)*y_inc), ! lat, lon position
     2        1000.0*EXSUMR,1000.0*EXSUMI,           ! real, imag parts of Ex in v/km
     3        1000.0*SQRT(EXSUMR**2 + EXSUMI**2),    ! scalar magnitude of Ex in v/km
     2        1000.0*EYSUMR,1000.0*EYSUMI,           ! real, imag parts of Ey in v/km
     3        1000.0*SQRT(EYSUMR**2 + EYSUMI**2)     ! scalar magnitude of Ey in v/km
          end do
      end do
      close(150)
      
c      

      close(6)
      close(100)
c
c     The end of Freak
c      END DO
c     THe end of NMOD - models runs
c      END DO

c     close netcdf file
c      call netcdf_error(nf90_close(ncFileID))
      STOP
 1000 FORMAT ('*** CHANGE PARAMETER NXYD FROM ',I3,' TO AT LEAST',I3)
 1010 FORMAT ('*** CHANGE PARAMETER NXY2D FROM ',I3,' TO AT LEAST',I3)
 1020 FORMAT ('*** CHANGE PARAMETER N1D FROM ',I3,' TO AT LEAST',I3)
 1030 FORMAT ('*** CHANGE PARAMETER N2D FROM ',I3,' TO AT LEAST',I3)

 1040 write(*,*) 'End.'

      END
      

C     ***************************************************************

      SUBROUTINE RHONEW
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
C     THIS SUBROUTINE DETERMINES A NEW NORMAL RESISTIVITY MODEL,    C
C     IN WHICH THE POSITION OF THE THIN SHEET ZA AND THE DATA LEVEL C
C     ZD COINCIDE WITH (ARTIFICIAL) LAYER BOUNDARIES. IF ZA OR ZD   C
C     DO NOT AGREE WITH AN ORIGINAL BOUNDARY, THE NUMBER OF LAYERS  C
C     INCREASES CORRESPONDINGLY. THE NEW MODEL IS GIVEN BY NL1, D1, C
C     AND RHO1. THE THIN SHEET LIES ON THE TOP OF LAYER NLA AND THE C
C     DATA LEVEL LIES ON THE TOP OF LAYER NLD. - BECAUSE  OF FLOAT- C
C     ING POINT  IMPRECISION  ZA  AND ZD  ARE ASSUMED TO COINCIDE   C
C     WITH A LAYER  BOUNDARY IF THEY  DEVIATE  BY A FRACTION  LESS  C
C     THAN EPS FROM THE RESPECTIVE DEPTH.                           C
C ****RHONEW HAS TO BE CALLED AFTER EACH CHANGE OF ZA and/or ZD**** C
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      COMMON /RESIST/NL,RHO(20),D(20),NL1,RHO1(22),D1(22),NLA,NLD,TAUN
      COMMON /PARAM/F,U,DX,NITER,ITERD,EPS,SOR,NH,ZA,ZD
      DIMENSION Z(2),N(2)
      ZA=AMAX1(ZA,0.)
      ZD=AMAX1(ZD,0.)
      Z(1)=ZA
      Z(2)=ZD
      D(NL)=777.777

      DO 10 I=1,NL
      RHO1(I)=RHO(I)
   10 D1(I)=D(I)
      NL1=NL
C
      DO 70 IZ=1,2
      H=0.
      NST=NL1
      DO 60 I=1,NST
      IF (   (Z(IZ)-H).LT.-EPS*H) GOTO 20
      IF (ABS(Z(IZ)-H).LE.+EPS*H) GOTO 40
      IF (   (Z(IZ)-H).GT.+EPS*H) GOTO 50
   20 N(IZ)=I
      IF (IZ.EQ.2.AND.N(1).GE.N(2)) N(1)=N(1)+1
      NL1=NST+1
      DO 30 II=NST,I,-1
      RHO1(II+1)=RHO1(II)
   30 D1(II+1)=D1(II)
      RHO1(I)=RHO1(I-1)
      D1(I)=H-Z(IZ)
      D1(I-1)=D1(I-1)-D1(I)
      GOTO 70
   40 N(IZ)=I
      NL1=NST
      GOTO 70
   50 IF (I.LT.NST) GOTO 60
      N(IZ)=NST+1
      NL1=NST+1
      RHO1(NL1)=RHO1(NST)
      D1(NL1)=D1(NST)
      D1(NST)=Z(IZ)-H
      GOTO 70
   60 H=H+D1(I)
   70 CONTINUE
      NLA=N(1)
      NLD=N(2)
      RETURN
      END
C
C     ***************************************************************
C
      SUBROUTINE GREEN1(GXX,GXY,IWG1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
C     THIS SUBROUTINE CALCULATES THE GREEN'S KERNELS REQUIRED IN    C
C     THE SOLUTION OF THE THIN SHEET INTEGRAL EQUATION. BECAUSE OF  C
C     THE SQUARE CELLS IN THE ANOMALOUS DOMAIN, COMPUTATION CAN BE  C
C     RESTRICTED TO GXX AND GXY. THEREFORE IT IS ONLY NECESSARY TO  C
C     PLACE AN ELECTRIC DIPOLE IN X-DIRECTION IN THE CENTRE OF THE  C
C     LOWER LEFT CORNER OF A SQUARED ANOMALOUS DOMAIN OF SIZE N1*N1,C
C     WHERE N1=MAX(NX,NY), AND TO CALCULATE THE HORIZONTAL ELEC-    C
C     TRIC  FIELD  IN  THE  CENTERS  OF  ALL  OTHER  CELLS.  THE    C
C     OTHER KERNEL ELEMENTS ARE DERIVED FROM THESE NUMBERS BY SYM-  C
C     METRY CONSIDERATIONS (cf. SUBROUTINE SEIDEL). THE DISCRETIZA- C
C     TION OF THE INTEGRAL EQUATON REQUIRES THE INTEGRATIION OVER   C
C     THE SOURCE CELL. FOR POINTS OF OBSERVATION OUTSIDE THE SINGU- C
C     LAR CELL THIS IS PERFORMED BY AN OPTIMIZED 9-POINT INTEGRA-   C
C     TION RULE (cf. ABRAMOWITZ & STEGUN, p.892/893). THE CONTRI-   C
C     BUTION OF THE SINGULAR CELL IS OBTAINED BY REPLACING THE      C
C     SQUARE CELL BY A CIRCLE OF EQUAL AREA. - THE WAVENUMBER IN-   C
C     TEGRATION IS PERFORMED BY A FAST HANKEL TRANSFORM, THE INTE-  C
C     GRATION KERNEL IS PROVIDED BY SUBROUTINE KERNLG1.             C
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      COMPLEX GXX(N1*N1),GXY(N1*N1)
      COMPLEX FT,FPD,AUX,ER,EP,C,CC,CD,CF
      COMMON /RESIST/NL,RHO(20),D(20),NL1,RHO1(22),D1(22),NLA,NLD,TAUN
      COMMON /NUMBER/NX,NY,NX2,NY2,NX0,NY0,N1,N2
      COMMON /PARAM/F,U,DX,NITER,ITERD,EPS,SOR,NH,ZA,ZD
      COMMON /HANKEL/NC,NC0,HR(2,100)
      COMMON /INTPOL/RR(28),CC(3,28,10),CD(10),CF(28,10)
      REAL X0(9),Y0(9),W0(9),XS(9),YS(9),WS(9)
      DATA PI,Q/3.141592654,1.258925412/
      DATA X0/-1.,0.,1.,-1.,0.,1.,-1.,0.,1./
      DATA Y0/3*-1,3*0.,3*1./
      DATA W0/25.,40.,25.,40.,64.,40.,25.,40.,25./
C
C     WEIGHTS FOR INTEGRATION OVER SOURCE CELL:
C     -----------------------------------------
      C1=0.5*DX*SQRT(3./5.)
      C2=DX*DX/324.
      DO 10 L=1,9
      XS(L)=X0(L)*C1
      YS(L)=Y0(L)*C1
   10 WS(L)=W0(L)*C2
C
      N11=N1*N1
      C=(0.,8.E-7)*PI*PI*F
      RMIN=DX/SQRT(PI)
      RMAX=N1*SQRT(2.)*DX
      NR=2.5+10.*ALOG10(RMAX/RMIN)
      NR=MAX0(NR,7)
      IF (NR.GT.21) WRITE (6,110) NR
      NR=MIN0(NR,21)
      NCNR=NC+NR-1
C
C     HANKEL TRANSFORM OF AUXILIARY FUNCTIONS:
C     ----------------------------------------
      DO 20 I=1,NR
      RR(I)=RMIN*Q**(I-1)
      DO 20 M=1,4
   20 CF(I,M)=0.
      DO 30 NN=1,NCNR
      N=-NC+NC0+NN
      U=Q**(1-N)/RMIN
      RD=EXPC(-0.005*U*RMIN)
C     REDUCTION FACTOR 'RD' FORCES CONVERGENCE
      CALL KERNLG1(FT,FPD)
      FT=FT*RD
      FPD=FPD*RD
      CD(1)=U*FT
      CD(2)=U*FPD
      CD(3)=FT-FPD
      CD(4)=FT+FPD
      I1=MAX0(1,NN-NC+1)
      I2=MIN0(NR,NN)
      DO 30 I=I1,I2
      IN=NC-NN+I
      CF(I,1)=CF(I,1)+CD(1)*HR(1,IN)
      CF(I,2)=CF(I,2)+CD(2)*HR(1,IN)
      CF(I,3)=CF(I,3)+CD(3)*HR(2,IN)
   30 CF(I,4)=CF(I,4)+CD(4)*HR(2,IN)
      DO 40 I=1,NR
      AUX=1./(2.*PI*RR(I))
      DO 35 M=1,4
   35 CF(I,M)=CF(I,M)*AUX
      CF(I,3)=CF(I,3)/RR(I)
   40 CF(I,4)=CF(I,4)/RR(I)
C
C     INTERPOLATION COEFFICIENTS FROM SPLINE1:
C     ----------------------------------------
      DO 50 M=1,3
   50 CALL SPLINE1(NR,CF(1,M),CC(1,1,M))
      RL1=ALOG10(RR(1))
      RL2=ALOG10(RR(NR))
C
      DO 60 IK=1,N11
      GXX(IK)=0.
   60 GXY(IK)=0.
C
C     CONTRIBUTION FROM SINGULAR CELL:
C     --------------------------------
      GXX(1)=DX*DX*CF(1,4)
C
C     INTERPOLATION AND INTEGRATION OVER SOURCE CELL:
C     -----------------------------------------------
      DO 80 L=1,9
      DO 80 I=1,N1
      X=(I-1)*DX-XS(L)
      DO 80 K=1,N1
      IK=(K-1)*N1+I
      IF (IK.EQ.1) GOTO 80
      Y=(K-1)*DX-YS(L)
      R=SQRT(X*X+Y*Y)
      CP=X/R
      SP=Y/R
      RL=ALOG10(R)
      DO 70 M=1,3
   70 CALL SPLINE2(NR,RL,RL1,RL2,CC(1,1,M),CD(M))
      ER=(CD(3)+CD(2))*CP
      EP=(CD(3)-CD(1))*SP
      GXX(IK)=GXX(IK)+WS(L)*(ER*CP-EP*SP)
      GXY(IK)=GXY(IK)+WS(L)*(EP*CP+ER*SP)
   80 CONTINUE
C
      DO 90 I=1,N1
      K=(I-1)*N1+1
      GXY(I)=0.
   90 GXY(K)=0.
C
C     OPTIONAL PRINT:
C     ---------------
      IF (IWG1.NE.1) RETURN
      WRITE (6,100)
      CALL PRINTF(GXX,N1,N1,1)
      CALL PRINTF(GXY,N1,N1,2)
      RETURN
  100 FORMAT (/'  GREEN"S KERNELS FOR INTEGRAL EQUATION'/2X,37(1H=))
  110 FORMAT (/'  ***REPLACE DIMENSION 21 IN /INTPOL/ BY',I3,'***'/)
      END
C
C     ***************************************************************
C
      SUBROUTINE KERNLG1(FT,FPD)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
C     THIS SUBROUTINE COMPUTES THE TOROIDAL AND POLOIDAL TRANSFER   C
C     FUNCTIONS REQUIRED FOR THE CALCULATION OF THE TANGENTIAL      C
C     ELECTRIC FIELD AT THE LEVEL NLA CAUSED BY A HORIZONTAL ELEC-  C
C     TRIC DIPOLE AT THE SAME LEVEL. THIS SUBROUTINE IS REQUIRED    C
C     IN GREEN1.                                                    C
C     FT=FT(ZA), FPD=FP'(ZD)                                        C
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      COMMON /RESIST/NL,RHO(20),D(20),NL1,RHO1(22),D1(22),NLA,NLD,TAUN
      COMMON /PARAM/F,U,DX,NITER,ITERD,EPS,SOR,NH,ZA,ZD
      COMPLEX C,CTH,FT,FPD,CEXPC,CSQRT
      COMPLEX A(12),BTD(12),BTU(12),BPD(12),BPU(12)
      DIMENSION B(12)
      DATA PI/3.141592654/
C
      C=(0.,8.E-7)*PI*PI*F
      U2=U*U
C
C     UPWARD TRANSFER FUNCTIONS:
C     --------------------------
      BTU(1)=U
      BPU(1)=U
      IF (NLA.EQ.1) GOTO 20
      DO 10 I=1,NLA-1
      A(I)=CSQRT(U2+C/RHO1(I))
      CTH=CEXPC(-2.*A(I)*D1(I))
      CTH=(1.-CTH)/(1.+CTH)
      B(I)=0.
      IF (I.GT.1) B(I)=RHO1(I)/RHO1(I-1)
      BTU(I+1)=A(I)*(BTU(I)+A(I)*CTH)/(A(I)+BTU(I)*CTH)
   10 BPU(I+1)=A(I)*(BPU(I)+A(I)*B(I)*CTH)/(A(I)*B(I)+BPU(I)*CTH)
C
C     DOWNWARD TRANSFER FUNCTIONS:
C     ----------------------------
   20 A(NL1)=CSQRT(U2+C/RHO1(NL1))
      BTD(NL1)=A(NL1)
      BPD(NL1)=A(NL1)
      IF (NLA.EQ.NL1) GOTO 40
      DO 30 I=NL1-1,NLA,-1
      A(I)=CSQRT(U2+C/RHO1(I))
      CTH=CEXPC(-2.*A(I)*D1(I))
      CTH=(1.-CTH)/(1.+CTH)
      B(I)=RHO1(I)/RHO1(I+1)
      BTD(I)=A(I)*(BTD(I+1)+A(I)*CTH)/(A(I)+BTD(I+1)*CTH)
   30 BPD(I)=A(I)*(BPD(I+1)+A(I)*B(I)*CTH)/(A(I)*B(I)+BPD(I+1)*CTH)
C
   40 FT=1./(BTD(NLA)+BTU(NLA)+C*TAUN)
      FPD=C*(TAUN+1./RHO1(NLA)/BPD(NLA))
      IF (NLA.EQ.1) GOTO 50
      FPD=FPD+C/RHO1(NLA-1)/BPU(NLA)
   50 FPD=1./FPD
      RETURN
      END
C
C     ***************************************************************
C
      SUBROUTINE GREEN2(BXX,BXY,BXZ,IWG2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
C     THIS SUBROUTINE CALCULATES THE GREEN'S KERNELS FOR THE MAPP-  C
C     ING OF THE E-FIELD IN THE ANOMALOUS DOMAIN ONTO THE ELECTRO-  C
C     MAGNETIC FIELD AT ZD. BECAUSE OF IDENTICAL CELL SIZE AND      C
C     TRANSLATION OF THE GRID BY INTEGER NUMBER OF CELLS (NX0,NY0), C
C     THE COMPUTATION REDUCES TO PLACE AN ELECTRIC DIPOLE IN X-DI-  C
C     RECTION IN THE CENTRE OF THE LOWER LEFT CORNER OF A SQUARED   C
C     DOMAIN OF SIZE N2*N2, WHERE                                   C
C                 N2=MAX(NX2+NX0,NY2+NY0,NX-NX0,NY-NY0)             C
C     AND TO DETERMINE THE MAGNETIC AND ELECTRIC FIELD IN ALL CELLS C
C     OF THE SQUARED DOMAIN, ASSUMING THE DIPOLE AT THE DEPTH Z=ZA  C
C     AND THE POINTS OF OBSERVATION AT ZD.  ALL REQUIRED KERNEL     C
C     ELEMENTS ARE DERIVED FROM THESE NUMBERS BY SYMMETRY CONSIDE-  C
C     RATIONS (cf. SUBROUTINE BE). - THE DISCRETIZATION OF THE IN-  C
C     TEGRAL  EQUATION  REQUIRES  THE INTEGRATIION  OVER THE SOURCE C
C     CELL. FOR POINTS OF OBSERVATION OUTSIDE THE SINGULAR CELL     C
C     CELL  THIS IS  PERFORMED  BY  AN OPTIMIZED 9-POINT INTEGRA-   C
C     TION RULE (cf. ABRAMOWITZ & STEGUN, p.892/893). THE CONTRI-   C
C     BUTION OF THE SINGULAR CELL IS OBTAINED BY REPLACING THE      C
C     SQUARE CELL BY A CIRCLE OF EQUAL AREA. - THE WAVENUMBER IN-   C
C     TEGRATION IS PERFORMED BY A FAST HANKEL TRANSFORM, THE INTE-  C
C     GRATION KERNEL IS PROVIDED BY SUBROUTINE KERNLG2 .            C
C     FOR IWG2=1 THE RESULTS ARE DISPLAYED.                         C
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      COMPLEX BXX(N2*N2),BXY(N2*N2),BXZ(N2*N2)
      COMPLEX FT,FP,FTD,FPD,FEZ,C,CC,CD,CF,AUX
      COMPLEX BR,BP,BZ
      COMMON /RESIST/NL,RHO(20),D(20),NL1,RHO1(22),D1(22),NLA,NLD,TAUN
      COMMON /PARAM/F,U,DX,NITER,ITERD,EPS,SOR,NH,ZA,ZD
      COMMON /HANKEL/NC,NC0,HR(2,100)
      COMMON /INTPOL/RR(28),CC(3,28,10),CD(10),CF(28,10)
      COMMON /NUMBER/NX,NY,NX2,NY2,NX0,NY0,N1,N2
      REAL X0(9),Y0(9),W0(9),XS(9),YS(9),WS(9)
      DATA PI,Q/3.141592654,1.258925412/
      DATA X0/-1.,0.,1.,-1.,0.,1.,-1.,0.,1./
      DATA Y0/3*-1,3*0.,3*1./
      DATA W0/25.,40.,25.,40.,64.,40.,25.,40.,25./
C
C     WEIGHTS FOR INTEGRATION OVER SOURCE CELL:
C     -----------------------------------------
      C1=0.5*DX*SQRT(3./5.)
      C2=DX*DX/324.
      DO 10 L=1,9
      XS(L)=X0(L)*C1
      YS(L)=Y0(L)*C1
   10 WS(L)=W0(L)*C2
C
      N22=N2*N2
      C=(0.,8.E-7)*PI*PI*F
C     RMIN=RADIUS OF CIRCLE WITH SAME AREA AS SHEET ELEMENT:
      RMIN=DX/SQRT(PI)
      RMAX=N2*SQRT(2.)*DX
      NR=2.5+10.*ALOG10(RMAX/RMIN)
      NR=MAX0(NR,7)
      IF (NR.GT.21) WRITE (6,1010) NR
      NR=MIN0(NR,21)
      NCNR=NC+NR-1
C
C     HANKEL TRANSFORM OF AUXILIARY FUNCTIONS:
C     ----------------------------------------
      DO 20 I=1,NR
      RR(I)=RMIN*Q**(I-1)
      DO 20 K=1,5
   20 CF(I,K)=0.
C
      DO 30 NN=1,NCNR
      N=-NC+NC0+NN
      U=Q**(1-N)/RMIN
      RD=EXPC(-0.005*U*RMIN)
C     REDUCTION FACTOR 'RD' FORCES CONVERGENCE
      CALL KERNLG2(FT,FP,FTD,FPD)
      FEZ=FPD/U
      IF (NLD.GT.1) FEZ=FP*RHO1(NLD-1)
      FT=FT*RD
      FP=FP*RD
      FTD=FTD*RD
      FPD=FPD*RD
      FEZ=FEZ*RD
      CD(1)=U*FTD
      CD(2)=U*C*FP
      CD(3)=C*FP-FTD
      CD(4)=U*U*FT
      CD(5)=C*FP+FTD
      I1=MAX0(1,NN-NC+1)
      I2=MIN0(NR,NN)
      DO 30 I=I1,I2
      IN=NC-NN+I
      DO 25 M=1,2
   25 CF(I,M)=CF(I,M)+CD(M)*HR(1,IN)
      DO 30 M=3,5
   30 CF(I,M)=CF(I,M)+CD(M)*HR(2,IN)
      DO 40 I=1,NR
      AUX=1./(2.*PI*RR(I))
      DO 35 M=1,5
   35 CF(I,M)=CF(I,M)*AUX
      CF(I,3)=CF(I,3)/RR(I)
   40 CF(I,5)=CF(I,5)/RR(I)
C
C     INTERPOLATION COEFFICIENTS FROM SPLINE1:
C     ---------------------------------------
      DO 50 M=1,4
   50 CALL SPLINE1(NR,CF(1,M),CC(1,1,M))
      RL1=ALOG10(RR(1))
      RL2=ALOG10(RR(NR))
C
      DO 60 IK=1,N22
      BXX(IK)=0.
      BXY(IK)=0.
   60 BXZ(IK)=0.
C
C     CONTRIBUTION FROM SINGULAR CELL:
C     --------------------------------
      BXY(1)=DX*DX*CF(1,5)
C
C     INTERPOLATION AND INTEGRATION OVER SOURCE CELL:
C     -----------------------------------------------
      DO 80 L=1,9
      DO 80 I=1,N2
      X=(I-1)*DX-XS(L)
      DO 80 K=1,N2
      IK=(K-1)*N2+I
      IF (IK.EQ.1) GOTO 80
      Y=(K-1)*DX-YS(L)
      R=SQRT(X*X+Y*Y)
      CP=X/R
      SP=Y/R
      RL=ALOG10(R)
      DO 70 M=1,4
   70 CALL SPLINE2(NR,RL,RL1,RL2,CC(1,1,M),CD(M))
      BR=(+CD(3)+CD(1))*SP
      BP=(-CD(3)+CD(2))*CP
      BZ=CD(4)*SP
      BXX(IK)=BXX(IK)+WS(L)*(BR*CP-BP*SP)
      BXY(IK)=BXY(IK)+WS(L)*(BP*CP+BR*SP)
      BXZ(IK)=BXZ(IK)+WS(L)*BZ
   80 CONTINUE
      AUX=-1./((0.,2.)*PI*F)
      DO 90 IK=1,N22
      BXX(IK)=BXX(IK)*AUX
      BXY(IK)=BXY(IK)*AUX
   90 BXZ(IK)=BXZ(IK)*AUX
C
      DO 100 IK=1,N2
      BXX(IK)=0.
  100 BXZ(IK)=0.
      DO 110 K=1,N2
      IK=(K-1)*N2+1
  110 BXX(IK)=0.
C
C     OPTIONAL PRINT:
C     ---------------
      IF (IWG2.NE.1) RETURN
      WRITE (6,1000)
      CALL PRINTF(BXX,N2,N2,7)
      CALL PRINTF(BXY,N2,N2,8)
      CALL PRINTF(BXZ,N2,N2,9)
      RETURN
 1000 FORMAT (/'  GREEN"S KERNELS FOR EM FIELD AT ZD'/2X,34(1H=))
 1010 FORMAT (/'  ***REPLACE DIMENSION 21 IN /INTPOL/ BY',I3,'***'/)
      END
C
C     ***************************************************************
C
      SUBROUTINE KERNLG2(FT,FP,FTD,FPD)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
C     THIS SUBROUTINE COMPUTES THE TOROIDAL AND POLOIDAL FIELDS     C
C     CAUSED BY A HORIZONTAL ELECTRIC DIPOLE IN A THIN SHEET        C
C     PROPAGATING FROM THE LEVEL NLA OF THE THIN SHEET TO THE LEVEL C
C     NLD OF THE POINT OF OBSERVATION. THIS SUBROUTINE IS REQUIRED  C
C     IN GREEN2.                                                    C
C     FT=FT(ZD), FP=SIGMA(ZD)*FP(ZD), FTD=FT'(ZD-0), FPD=FP'(ZD)    C
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      COMMON /RESIST/NL,RHO(20),D(20),NL1,RHO1(22),D1(22),NLA,NLD,TAUN
      COMMON /PARAM/F,U,DX,NITER,ITERD,EPS,SOR,NH,ZA,ZD
      COMPLEX C,CTH,FT,FP,FTD,FPD,E,CEXPC,CSQRT
      COMPLEX A(12),BTD(12),BTU(12),BPD(12),BPU(12)
      DIMENSION B(12)
      DATA PI/3.141592654/
C
      C=(0.,8.E-7)*PI*PI*F
      U2=U*U
C
C     UPWARD TRANSFER FUNCTIONS:
C     --------------------------
      BTU(1)=U
      BPU(1)=U
      IF (NLA.EQ.1) GOTO 20
      DO 10 I=1,NLA-1
      A(I)=CSQRT(U2+C/RHO1(I))
      CTH=CEXPC(-2.*A(I)*D1(I))
      CTH=(1.-CTH)/(1.+CTH)
      B(I)=0.
      IF (I.GT.1) B(I)=RHO1(I)/RHO1(I-1)
      BTU(I+1)=A(I)*(BTU(I)+A(I)*CTH)/(A(I)+BTU(I)*CTH)
   10 BPU(I+1)=A(I)*(BPU(I)+A(I)*B(I)*CTH)/(A(I)*B(I)+BPU(I)*CTH)
C
C     DOWNWARD TRANSFER FUNCTIONS:
C     ----------------------------
   20 A(NL1)=CSQRT(U2+C/RHO1(NL1))
      BTD(NL1)=A(NL1)
      BPD(NL1)=A(NL1)
      IF (NLA.EQ.NL1) GOTO 40
      DO 30 I=NL1-1,NLA,-1
      A(I)=CSQRT(U2+C/RHO1(I))
      CTH=CEXPC(-2.*A(I)*D1(I))
      CTH=(1.-CTH)/(1.+CTH)
      B(I)=RHO1(I)/RHO1(I+1)
      BTD(I)=A(I)*(BTD(I+1)+A(I)*CTH)/(A(I)+BTD(I+1)*CTH)
   30 BPD(I)=A(I)*(BPD(I+1)+A(I)*B(I)*CTH)/(A(I)*B(I)+BPD(I+1)*CTH)
C
   40 FT=1./(BTD(NLA)+BTU(NLA)+C*TAUN)
      FPD=C*(TAUN+1./RHO1(NLA)/BPD(NLA))
      IF (NLA.EQ.1) GOTO 50
      FPD=FPD+C/RHO1(NLA-1)/BPU(NLA)
   50 FPD=1./FPD
      E=0.
      IF (NLA-NLD) 60,80,80
C
C     DATA LEVEL BELOW NON-UNIFORM SHEET:
C     -----------------------------------
   60 DO 70 I=NLA,NLD-1
      FT=FT*(A(I)+BTD(I))/(A(I)+BTD(I+1))
      FPD=FPD*(1.+A(I)/BPD(I))/(1.+A(I)*B(I)/BPD(I+1))
   70 E=E+A(I)*D1(I)
      E=CEXPC(-E)
      FT=FT*E
      FPD=FPD*E
      FTD=-FT*BTD(NLD)
      FP=-FPD/BPD(NLD)/RHO1(NLD)
      RETURN
C
C     DATA LEVEL ABOVE NON-UNIFORM SHEET:
C     -----------------------------------
   80 FP=0.
      IF (NLA.GT.1) FP=FPD/BPU(NLA)/RHO1(NLA-1)
      FTD=FT*BTU(NLA)
      IF (NLA.EQ.NLD) RETURN
      DO 90 I=NLA-1,NLD,-1
      FT=FT*(A(I)+BTU(I+1))/(A(I)+BTU(I))
      FPD=FPD*(1.+A(I)/BPU(I+1))/(1.+A(I)*B(I)/BPU(I))
   90 E=E+A(I)*D1(I)
      E=CEXPC(-E)
      FT=FT*E
      FPD=FPD*E
      FTD=FT*BTU(NLD)
      FP=0.
      IF (NLD.GT.1) FP=FPD/BPU(NLD)/RHO1(NLD-1)
      RETURN
      END
C
C     ***************************************************************
C
      SUBROUTINE NORMAL1(ENX,ENY,IWN1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
C     THIS SUBROUTINE CALCULATES THE TANGENTIAL ELECTRIC FIELD IN   C
C     THE ANOMALOUS DOMAIN OF THE THIN SHEET FOR THE NORMAL CONDUC- C
C     TIVITY STRUCTURE. THESE COMPONENTS ARE OUTPUT IN THE ONE-     C
C     DIMENSIONAL ARRAYS ENX AND ENY. THEY ARE PRINTED FOR IWN1=1.  C
C     THE INTEGRATION KERNEL IS PROVIDED BY SUBROUTINE KERNLN1.     C
C     THREE SOURCE FIELDS CAN BE USED:                              C
C                                                                   C
C     ITYPE=1: UNIFORM MAGNETIC SURFACE FIELD B0 (T) WITH  AZIMUTH  C
C              THETA (deg), COUNTED FROM X IN DIRECTION TO Y        C
C     ITYPE=2: VERTICAL MAGNETIC DIPOLE OF MOMENT TM (A*m**2) AT    C
C              (XD,YD), WHERE THE LOWER LEFT CORNER OF THE ANOMA-   C
C              LOUS DOMAIN IS THE ORIGIN OF COORDINATES.            C
C     ITYPE=3: HORIZONTAL ELECTRIC DIPOLE OF CURRENT MOMENT TD (A*m)C
C              AT (XD,YD) WITH THE SAME ORIGIN AS ITYPE=2. THE DI-  C
C              POLE HAS THE AZIMUTH THETA (deg) (cf. ITYPE=1).      C
C                                                                   C
C ****IF A DIPOLE SOURCE LIES over THE ANOMALOUS DOMAIN, IT SHOULD  C
C     BE POSITIONED OVER THE centre OF A CELL. COMPUTATIONALLY IT   C
C     WILL BE 'SMEARED OUT' OVER THAT CELL.                         C
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      COMPLEX C,ENX(NX*NY),ENY(NX*NY),FT
      complex STFTC,CTFTZ
      
      COMMON /RESIST/NL,RHO(20),D(20),NL1,RHO1(22),D1(22),NLA,NLD,TAUN
      COMMON /PARAM/F,U,DX,NITER,ITERD,EPS,SOR,NH,ZA,ZD
      COMMON /SOURCE/B0,THETA,BGRAD(43460)
      COMMON /NUMBER/NX,NY,NX2,NY2,NX0,NY0,N1,N2
      DATA PI/3.141592654/
C
C     UNIFORM SURFACE MAGNETIC FIELD B0 WITH AZIMUTH THETA
C     ----------------------------------------------------
C     THETA IS COUNTED POSITIVE FROM NORTH (X) TO EAST (Y):
      U=0.
C
      CT=COS(THETA*PI/180.)
      ST=SIN(THETA*PI/180.)
      C=(0.,8.E-7)*PI*PI*F

      CALL KERNLN1(FT)

      STFTC = ST*FT*(0.,2.)*PI*F
      CTFTZ = CT*FT*(0.,2.)*PI*F

c     change 2.12.03 
c      ENX(1)=+ST*FT*B0*(0.,2.)*PI*F
c      ENY(1)=-CT*FT*B0*(0.,2.)*PI*F

      NXY=NX*NY
c     change 2.12.03 
c      DO 10 IK=1,NXY  ! used for a uniform field everywhere
c      ENX(IK)=ENX(1)
c   10 ENY(IK)=ENY(1)
      do IK=1,NXY
         ENX(IK) =  BGRAD(IK)*STFTC ! for a non-uniform magnetic field: BGRAD, in T
         ENY(IK) = -BGRAD(IK)*CTFTZ
      end do

      write(*,*) 
     1 ' MAIN PROG/SUBR NORMAL: Normal field Ex,Ey',ENX(1),ENY(1)


C
C     OPTIONAL PRINT:
C     ---------------
      IF (IWN1.NE.1) RETURN
      WRITE (6,100)
      CALL PRINTF(ENX,NX,NY,3)
      CALL PRINTF(ENY,NX,NY,4)
      RETURN
  100 FORMAT (/'  NORMAL TANGENTIAL E-FIELD IN THE ANOMALOUS DOMAIN'/
     *2X,49(1H=)/)
      END
C
C     ***************************************************************
C
      SUBROUTINE KERNLN1(FT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
C     THIS SUBROUTINE CALCULATES KERNELS IN THE WAVENUMBER DOMAIN   C
C     FOR THE POLOIDAL AND TOROIDAL CURRENT MODE, WHICH RELATE      C
C     FIELDS AT THE SURFACE Z=0 TO FIELDS AT THE DEPTH OF THE ANO-  C
C     MALOUS SHEET AT Z = ZA (LEVEL NLA). THEY ARE USED BOTH FOR    C
C     THE CALCULATION OF THE NORMAL FIELD AT Z = ZA (cf. NORMAL1)   C
C     AND FOR THE GREEN'S KERNELS MAPPING THE E-FIELD AT Z=ZA TO    C
C     THE SURFACE (cf. GREEN2). THE KERNELS CALCULATED ARE:         C
C                                                                   C
C     FT=FT(ZA)/FT(0)/(U+BT), BT=-FT'(0)/FT(0)                      C
C                                                                   C
C     FT(Z) IS A DOWNWARD DIFFUSING FIELD.                          C
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
c     i.e. find the surface response for a given conductivity model.

      COMPLEX A(22),BT(22),C,FT,CTH,E,CSQRT,CEXPC
      COMMON /RESIST/NL,RHO(20),D(20),NL1,RHO1(22),D1(22),NLA,NLD,TAUN
      COMMON /PARAM/F,U,DX,NITER,ITERD,EPS,SOR,NH,ZA,ZD
      DATA PI/3.141592654/
C
C     DOWNWARD TRANSFER FUNCTIONS:
C     ----------------------------
c     this the coding of the recurrence relation on p75 of Weaver. 
c     The complex hyperbolic functions in the recurrence relation
c     replaced by complex exponentials. 

      U2=U*U
      C=(0.,8.E-7)*PI*PI*F
      A(NL1)=CSQRT(U2+C/RHO1(NL1))
      BT(NL1)=A(NL1)
      IF (NL1.NE.NLA) GOTO 10
      BT(NL1)=BT(NL1)+C*TAUN
   10 IF (NL1.EQ.1) GOTO 30
C
      DO 20 I=NL1-1,1,-1
      A(I)=CSQRT(U2+C/RHO1(I))
      CTH=CEXPC(-2.*A(I)*D1(I))
      CTH=(1.-CTH)/(1.+CTH)
      BT(I)=A(I)*(BT(I+1)+A(I)*CTH)/(A(I)+BT(I+1)*CTH)
      IF (I.NE.NLA) GOTO 20
      BT(I)=BT(I)+C*TAUN
   20 CONTINUE
C
C     DOWNWARD CONTINUATION FROM THE SURFACE TO THE LEVEL ZA:
C     -------------------------------------------------------
   30 FT=1./(U+BT(1))
      IF (NLA.EQ.1) RETURN
      E=0.
      DO 40 I=1,NLA-1
      E=E+A(I)*D1(I)
   40 FT=FT*(A(I)+BT(I))/(A(I)+BT(I+1))
      E=CEXPC(-E)
      FT=FT*E
      RETURN
      END
C
C     ***************************************************************
C
      SUBROUTINE NORMAL2(BX,BY,BZ,IWN2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
C     THIS SUBROUTINE CALCULATES THE SIX COMPONENTS OF THE NORMAL   C
C     ELECTROMAGNETIC FIELD  AT THE DATA LEVEL ZD FOR THE THREE     C
C     SOURCE FIELDS DESCRIBED IN THE CAPTION TO NORMAL1. FOR IWN2=1 C
C     THE FIELD IS DISPLAYED. - THE INTEGRATION KERNEL IS PROVIDED  C
C     BY SUBROUTINE KERNLN2.                                        C
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      COMPLEX C,C0,FT,FTD
      complex CTFTD,STFTD

      COMPLEX BX(NX2*NY2),BY(NX2*NY2),BZ(NX2*NY2)
      COMMON /RESIST/NL,RHO(20),D(20),NL1,RHO1(22),D1(22),NLA,NLD,TAUN
      COMMON /PARAM/F,U,DX,NITER,ITERD,EPS,SOR,NH,ZA,ZD
      COMMON /SOURCE/B0,THETA,BGRAD(43460)
      COMMON /NUMBER/NX,NY,NX2,NY2,NX0,NY0,N1,N2
      DATA PI/3.141592654/
C
C     UNIFORM SURFACE MAGNETIC FIELD B0 WITH AZIMUTH THETA
C     ----------------------------------------------------
C     THETA IS COUNTED POSITIVE FROM NORTH (X) TO EAST (Y):
      U=0.
C
      CT=COS(THETA*PI/180.)
      ST=SIN(THETA*PI/180.)
      C=(0.,8.E-7)*PI*PI*F
C
      NXY2=NX2*NY2
      DO 10 IK=1,NXY2
      BX(IK)=0.
      BY(IK)=0.
   10 BZ(IK)=0.
C
      CALL KERNLN2(FT,FTD,C0)
c     added 2.12.03
      CTFTD = CT*FTD
      STFTD = ST*FTD
c      change 2.12.03
c      BX(1)=-B0*CT*FTD
c      BY(1)=-B0*ST*FTD
c      DO 20 IK=1,NXY2
c      BX(IK)=BX(1) ! used for a uniform field
c   20 BY(IK)=BY(1)

c     added 2.12.03
      do IK=1,NXY2
         BX(IK) = -CTFTD*BGRAD(IK) ! for a non-uniform magnetic field: BGRAD, in T
         BY(IK) = -STFTD*BGRAD(IK)
      end do
      

C
C     OPTIONAL PRINT:
C     ---------------
      IF (IWN2.NE.1) RETURN
      WRITE (6,1000) ZD
      CALL PRINTF(BX,NX2,NY2,13)
      CALL PRINTF(BY,NX2,NY2,14)
      CALL PRINTF(BZ,NX2,NY2,15)
      RETURN
 1000 FORMAT (/'  NORMAL EM FIELD IN THE DATA DOMAIN (ZD=',F8.1,' m)'
     */2X,50(1H=)/)
      END
C
C     ***************************************************************
C
      SUBROUTINE KERNLN2(FT,FTD,C0)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
C     THIS SUBROUTINE CALCULATES KERNELS IN THE WAVENUMBER DOMAIN   C
C     FOR THE POLOIDAL AND TOROIDAL CURRENT MODE, WHICH RELATE      C
C     FIELDS AT THE SURFACE Z=0 TO FIELDS AT THE DATA LEVEL ZD      C
C     (LEVEL NLD). THEY ARE USED FOR THE CALCULATION OF THE NORMAL  C
C     FIELD AT THE DATA LEVEL ZD, cf. NORMAL2.                      C
C     THE KERNELS CALCULATED ARE:                                   C
C                                                                   C
C     FTD=FT'(ZD)/FT(0)/(U+BT)                                      C
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
c      COMPLEX A(12),BT(12),C,C0,FT,FTD,CTH,E,CSQRT,CEXPC
      COMPLEX A(22),BT(22),C,C0,FT,FTD,CTH,E,CSQRT,CEXPC
      COMMON /RESIST/NL,RHO(20),D(20),NL1,RHO1(22),D1(22),NLA,NLD,TAUN
      COMMON /PARAM/F,U,DX,NITER,ITERD,EPS,SOR,NH,ZA,ZD
      DATA PI/3.141592654/
C
C     DOWNWARD TRANSFER FUNCTIONS:
C     ----------------------------
      U2=U*U
      C=(0.,8.E-7)*PI*PI*F
      A(NL1)=CSQRT(U2+C/RHO1(NL1))
      BT(NL1)=A(NL1)
      IF (NL1.NE.NLA) GOTO 10
      BT(NL1)=BT(NL1)+C*TAUN
   10 IF (NL1.EQ.1) GOTO 30
C
c     from lowest layer to the top..
      DO 20 I=NL1-1,1,-1
      A(I)=CSQRT(U2+C/RHO1(I))
      CTH=CEXPC(-2.*A(I)*D1(I))
      CTH=(1.-CTH)/(1.+CTH)
      BT(I)=A(I)*(BT(I+1)+A(I)*CTH)/(A(I)+BT(I+1)*CTH)
      IF (I.NE.NLA) GOTO 20
      BT(I)=BT(I)+C*TAUN
   20 CONTINUE
C
C     DOWNWARD CONTINUATION FROM THE SURFACE TO THE LEVEL ZD:
C     -------------------------------------------------------
c     for a thin sheet at the surface NLD=1 and this part is skipped
   30 FT=1./(U+BT(1))
      IF (NLD.EQ.1) GOTO 60
      E=0.
      DO 50 I=1,NLD-1
      IF (I.NE.NLA) GOTO 40
      BT(I)=BT(I)-C*TAUN
   40 FT=FT*(A(I)+BT(I))/(A(I)+BT(I+1))
   50 E=E+A(I)*D1(I)
      E=CEXPC(-E)
      FT=FT*E
      write(*,*) 'FT----', FT
      write(*,*) 'BT----', BT
      write(*,*) 'NLD----', NLD
   60 FTD=-FT*BT(NLD)
      C0=1./BT(1)
      RETURN
      END
C
C     ***************************************************************
C

C
      SUBROUTINE BE(TAUA,BXX,BXY,BXZ,BX,BY,BZ,ESX,ESY,IWF,IPOL)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
C     THIS SUBROUTINE CALCULATES THE TOTAL MAGNETIC SURFACE         C
C     FIELD FROM THE ELECTRIC FIELD IN THE ANOMALOUS DOMAIN. IN     C
C     ADDITION TO THE E-FIELD FROM SEIDEL IT REQUIRES THE KERNELS   C
C     GREEN2 AND THE NORMAL SURFACE FIELD FROM NORMAL2.             C
C     FOR IWF=1 THE FIELD COMPONENTS ARE DISPLAYED.                 C
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      COMPLEX BXX(N2*N2),BXY(N2*N2),BXZ(N2*N2)
      COMPLEX BX(NX2*NY2),BY(NX2*NY2),BZ(NX2*NY2)
      COMPLEX ESX(1),ESY(1)
      COMPLEX C,CTAU
      COMPLEX FXX,FXY,FXZ,FYX,FYY,FYZ
      REAL TAUA(1)
      COMMON /NUMBER/NX,NY,NX2,NY2,NX0,NY0,N1,N2
      COMMON /PARAM/F,U,DX,NITER,ITERD,EPS,SOR,NH,ZA,ZD
      DATA PI/3.141592654/
c      character*12 file
C
      C=(0.,8.E-7)*PI*PI*F
C
      DO 20 I=1,NX
      DO 20 K=1,NY
      IK=(K-1)*NX+I
      IF (ABS(TAUA(IK)).LT.1.E-5) GOTO 20
      CTAU=C*TAUA(IK)
      DO 10 L=1,NX2
      DO 10 M=1,NY2
      LM=(M-1)*NX2+L
      LI=L-I+NX0
      MK=M-K+NY0
      IDX=IABS(MK)*N2+IABS(LI)+1
      IDY=IABS(LI)*N2+IABS(MK)+1
      SX=ISIGN(1,LI)
      SY=ISIGN(1,MK)
      SXY=SX*SY
      FXX=BXX(IDX)*SXY
      FXY=BXY(IDX)
      FXZ=BXZ(IDX)*SY
      FYX=-BXY(IDY)
      FYY=-BXX(IDY)*SXY
      FYZ=-BXZ(IDY)*SX
      BX(LM)=BX(LM)-CTAU*(ESX(IK)*FXX+ESY(IK)*FYX)
      BY(LM)=BY(LM)-CTAU*(ESX(IK)*FXY+ESY(IK)*FYY)
      BZ(LM)=BZ(LM)-CTAU*(ESX(IK)*FXZ+ESY(IK)*FYZ)
   10 CONTINUE
   20 CONTINUE
C
C     OPTIONAL PRINT:
C     ---------------
      IF (IWF.NE.1) RETURN
      WRITE (6,100) ZD
c      CALL PRINTF(BX,NX2,NY2,13)
c      CALL PRINTF(BY,NX2,NY2,14)
c      CALL PRINTF(BZ,NX2,NY2,15)

c      call write_magfield(file,BX,BY,BZ,IPOL)
      

      RETURN
  100 FORMAT (/'  TOTAL MAGNETIC FIELD IN THE DATA DOMAIN (ZD=',F8.1,
     *' m)'/2X,55(1H=)/)
      END
C
C     ***************************************************************
C
      SUBROUTINE TIPPER(ZZ,IWZ)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
C     THIS SUBROUTINE CALCULATES THE INDUCTION ARROWS (TIPPER)      C
C     FROM THE DATA FOR TWO LINEARLY INDEPENDENT POLARISATIONS OF A C
C     UNIFORM SOURCE FIELD. FOR A POINT IK IN THE MEASURING AREA    C
C     THE ELEMENTS OF THE OUTPUT ARRAY ZZ HAVE THE FOLLOWING MEAN-  C
C     ING MEANING                                                   C
C          ZZ(1,IK):  A               ZZ(2,IK):  B                  C
C     WITH IK=(K-1)*NX2+I, I=1,..., NX2, K=1,..., NY2 AND           C
C                        BZ=A*BX+B*BY                               C
C                                                                   C
C     FOR EACH SQUARE IN THE DATA DOMAIN THE FOLLOWING QUANTITIES   C
C     ARE COMPUTED AND PRINTED (IWZ=1):                             C
C      A : REAL AND IMAGINARY X-COMPONENT OF INDUCTION ARROW        C
C      B : REAL AND IMAGINARY Y-COMPONENT OF INDUCTION ARROW        C
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      COMPLEX ZZ(6,1),DET
      COMPLEX BX1,BX2,BY1,BY2,BZ1,BZ2
      COMMON /PARAM/F,U,DX,NITER,ITERD,EPS,SOR,NH,ZA,ZD
      COMMON /NUMBER/NX,NY,NX2,NY2,NX0,NY0,N1,N2
      COMMON /RESIST/NL,RHO(20),D(20),NL1,RHO1(22),D1(22),NLA,NLD,TAUN
      REAL NAME(2)
      DATA NAME/3H A ,3H B /
C
      NXY2=NX2*NY2
C
      DO 10 IK=1,NXY2
      BX1=ZZ(1,IK)
      BY1=ZZ(2,IK)
      BZ1=ZZ(3,IK)
      BX2=ZZ(4,IK)
      BY2=ZZ(5,IK)
      BZ2=ZZ(6,IK)
C
      DET=BX1*BY2-BX2*BY1
      ZZ(1,IK)=(BZ1*BY2-BZ2*BY1)/DET
   10 ZZ(2,IK)=(BZ2*BX1-BZ1*BX2)/DET

      !CALL write_tipper(ZZ)
      IF (IWZ.NE.1) RETURN
C
      WRITE (6,140)
      NP=1+(NY2-1)/NH
      DO 80 L=1,NP
      K1=1+(L-1)*NH
      K2=MIN0(NY2,L*NH)
      WRITE (6,100) (K, K=K1,K2)
      DO 80 I=NX2,1,-1
      IK1=(K1-1)*NX2+I
      IK2=(K2-1)*NX2+I
      DO 70 M=1,2
      IF (M.EQ.1) WRITE (6,110) I,NAME(M),(ZZ(M,IK),IK=IK1,IK2,NX2)
   70 IF (M.NE.1) WRITE (6,120)   NAME(M),(ZZ(M,IK),IK=IK1,IK2,NX2)
   80 WRITE (6,130)
      RETURN
  100 FORMAT(I20,4I22)
  110 FORMAT (1X,I3,1X,A3,5(1PE12.3,E10.3))
  120 FORMAT (4X,A3,5(1PE12.3,E10.3))
  130 FORMAT ()
  140 FORMAT (/'  TIPPER'/2X,6(1H=)//,
     #'  FOR EACH SQUARE IN THE DATA DOMAIN THE FOLLOWING QUANTI- '/
     #'  TIES ARE PRINTED:'/
     #'    A : REAL AND IMAGINARY X-COMPONENT OF INDUCTION ARROW'/
     #'    B : REAL AND IMAGINARY Y-COMPONENT OF INDUCTION ARROW'/)
      END
C
C     ***************************************************************
C
      SUBROUTINE PRINTRHO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
C     THIS SUBROUTINE PRINTS THE NORMAL CONDUCTIVITY DISTRIBUTION   C
C     (INCLUDING THE THIN SHEET OF CONDUCTANCE TAUN AT Z=ZA).       C
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      USE netcdf
c      USE netcdf_varid


      COMMON /RESIST/NL,RHO(20),D(20),NL1,RHO1(22),D1(22),NLA,NLD,TAUN
      COMMON /PARAM/F,U,DX,NITER,ITERD,EPS,SOR,NH,ZA,ZD
C
      WRITE (6,100) ZD

C     Netcdf call to write Data Level

c      call netcdf_error(nf90_put_var(ncfileID, sheetdepthID, ZA))

      A1=0.
      IF (NLD.EQ.1) WRITE (6,140)
      IF (NLA.EQ.1) WRITE (6,130) TAUN
      IF (NL1.EQ.1) GOTO 20
      DO 10 IL=1,NL1-1
      A2=A1+D1(IL)
      WRITE (6,110) A1,A2,RHO1(IL),D1(IL)
      IF (NLD.EQ.IL+1) WRITE (6,140)
      IF (NLA.EQ.IL+1) WRITE (6,130) TAUN
   10 A1=A2
   20 WRITE (6,120) A1,RHO1(NL1)
      RETURN
  100 FORMAT (/'  NORMAL CONDUCTIVITY STRUCTURE WITH DATA LEVEL AT',
     *' ZD =',F7.0,' m'/2X,62(1H=)/)
  110 FORMAT ('  FROM  ',F7.0,' m TO ',F7.0,' m: RHO=',F7.1,
     *' Ohm*m (THICKNESS',F7.0, ' m)')
  120 FORMAT ('  BELOW ',F7.0,' m:              RHO=',F7.1,' Ohm*m')
  130 FORMAT ('  THIN SHEET OF CONDUCTANCE ',F7.1,' S')
  140 FORMAT ('  DATA LEVEL')
      END
C
C     ***************************************************************
C
      SUBROUTINE PRINTF(A,NX,NY,L)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
C     THIS SUBROUTINE PRINTS A TWO-DIMENSIONAL COMPLEX ARRAY A(I,K) C
C     WITH I=1,..,NX (INCREASING FROM BOTTOM TO TOP) AND K=1,..,NY  C
C     (INCREASING FROM LEFT TO RIGHT). THE ELEMENTS ARE STORED IN   C
C     THE ONE-DIMENSIONAL ARRAY A(IK), IK=(K-1)*NX+I.- NAME(L) IS A C
C     THREE CHARACTER WORD PRINTED AS HEADING. - FOR SCREEN DISPLAY C
C     NH=3 COMPLEX ELEMENTS ARE PRINTED IN HORIZONTAL DIRECTON. FOR C
C     PRINTER OUTPUT NH=5 IS APPROPRIATE.                           C
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      COMPLEX A(NX*NY)
      COMMON /PARAM/F,U,DX,NITER,ITERD,EPS,SOR,NH,ZA,ZD
      REAL NAME(22)
      DATA NAME/3HGXX,3HGXY,3HENX,3HENY,3HESX,3HESY,
     *          3HBXX,3HBXY,3HBXZ,3HEXX,3HEXY,3HEXZ,
     *          3H BX,3H BY,3H BZ,3H EX,3H EY,3H EZ,
     *          3H FX,3H FY,3H DA,3H DB/
C
      WRITE (6,100) NAME(L)
c      WRITE (*,*) NAME(L)
      IF ((L.EQ.5).OR.(L.EQ.6)) THEN
         goto 30
      end if
         NP=1+(NY-1)/NH
         DO 20 J=1,NP
            K1=1+(J-1)*NH
            K2=MIN0(NY,J*NH)
            IF (J.EQ.1) WRITE (6,110) (K,K=K1,K2)
            IF (J.NE.1) WRITE (6,120) (K,K=K1,K2)
            DO 10 I=NX,1,-1
               IK1=(K1-1)*NX+I
               IK2=(K2-1)*NX+I
 10            WRITE (6,130) I,(A(IK),IK=IK1,IK2,NX)


 20            WRITE (6,140)


c     write electric field to a seperate file
 30     write(*,*) ' MAIN PROG/ SUB PRINTF: ASCII efield write' 
c      IF ((L.EQ.5).OR.(L.EQ.6)) THEN

c         write(100,100) NAME(L)
c         DO j=1,NY
c            DO i=NX,1,-1
c               IJ=(j-1)*NX+i
c               write(100,150) i,j, A(IJ)
c            END DO
c         END DO
c      END IF


      RETURN
  100 FORMAT (2X,A3)
  110 FORMAT (2X,3(1H-),I12,4I22)
  120 FORMAT (I17,4I22)
  130 FORMAT (2X,I2,5(1PE12.3,E10.3))
  140 FORMAT ()
 150  FORMAT (2I5,1X,5(1PE10.3,1X,E10.3))
      END
C
C     ***************************************************************
C
      SUBROUTINE PRINTH(TAU)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
C     THIS SUBROUTINE PRINTS MAJOR PARTS OF THE HEADING, INCLUDING  C
C     THE CONDUCTANCE DISTRIBUTION IN THE ANOMALOUS DOMAIN AND THE  C
C     DEFINITION OF THE DATA DOMAIN.                                C
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DIMENSION TAU(NX*NY),IDX(100)
      COMMON /PARAM/F,U,DX,NITER,ITERD,EPS,SOR,NH,ZA,ZD
      COMMON /NUMBER/NX,NY,NX2,NY2,NX0,NY0,N1,N2
      NH2=2*NH
      WRITE (6,1060)
      WRITE (6,1000) ZA,DX
      WRITE (6,1010) NX,NY
      WRITE (6,1020)
      NP=1+(NY-1)/NH2
      DO 20 L=1,NP
      K1=1+(L-1)*NH2
      K2=MIN0(NY,L*NH2)
      WRITE (6,1030) (K,K=K1,K2)
      DO 10 I=NX,1,-1
      IK1=(K1-1)*NX+I
      IK2=(K2-1)*NX+I
   10 WRITE (6,1040) I,(TAU(IK),IK=IK1,IK2,NX)
   20 WRITE (6,1050)
C
      WRITE (6,1070)
      WRITE (6,1080) NX2,NY2,NX0,NY0
C
      IMAX=MAX0(NX,NX2,NX-NX0,NX2+NX0)
      KMAX=MAX0(NY,NY2,NY-NY0,NY2+NY0)
      IF (KMAX.LE.100) GOTO 30
      WRITE (6,1110) KMAX
      RETURN
   30 I1=MAX0(1,1-NX0)
      I2=I1+NX-1
      I3=MAX0(1,1+NX0)
      I4=I3+NX2-1
      K1=MAX0(1,1-NY0)
      K2=K1+NY-1
      K3=MAX0(1,1+NY0)
      K4=K3+NY2-1
      NH2=10*NH
      NP=1+(KMAX-1)/NH2
C
      DO 90 I=IMAX,1,-1
      DO 40 K=1,KMAX
   40 IDX(K)=0
      IF (I.LT.I1.OR.I.GT.I2) GOTO 60
      DO 50 K=K1,K2
   50 IDX(K)=1
   60 IF (I.LT.I3.OR.I.GT.I4) GOTO 80
      DO 70 K=K3,K4
   70 IDX(K)=3-IDX(K)
C
   80 DO 90 L=1,NP
      IK1=(L-1)*NH2+1
      IK2=MIN0(KMAX,L*NH2)
   90 WRITE (6,1090) (IDX(IK),IK=IK1,IK2)
      WRITE (6,1100)
      RETURN
 1000 FORMAT('  THE ANOMALOUS THIN SHEET LIES AT THE DEPTH ZA=',F9.1,
     *' m'/'  IT CONSISTS OF SQUARE CELLS OF WIDTH DX=',F9.1,' m')
 1010 FORMAT ('  NUMBER OF CELLS IN X-DIRECTION: NX=',I2/
     *        '  NUMBER OF CELLS IN Y-DIRECTION: NY=',I2/)
 1020 FORMAT ('  CONDUCTANCE TAU (S) IN THE ANOMALOUS DOMAIN:'/)
 1030 FORMAT (3X,10I9)
 1040 FORMAT (2X,I2,10F9.1)
 1050 FORMAT ()
 1060 FORMAT (/'  ANOMALOUS CONDUCTIVITY STRUCTURE'/2X,32(1H=)/)
 1070 FORMAT (/'  DEFINITION OF THE DATA DOMAIN'
     */2X,29(1H=)//
     *'  THE FIELD AT DATA LEVEL ZD IS DEFINED IN THE CENTER OF'/
     *'  SQUARE CELLS OF THE SAME SIZE AS IN THE ANOMALOUS DOMAIN.'/
     /'  THE CELLS FORM A RECTANGLE PARALLEL TO THE ANOMALY.')
 1080 FORMAT ('  NX2=',I3,' CELLS LIE IN X-DIRECTION AND',
     *        ' NY2=',I3,' CELLS IN Y-DIRECTION.'/
     *'  THE RECTANGLE IS SHIFTED BY NX0=',I3,' CELLS IN X-DIRECTION'
     */'  AND NY0=',I3,' CELLS IN Y-DIRECTION. SCHEMATICALLY:'/)
 1090 FORMAT (12X,50I2)
 1100 FORMAT (
     */'  1: ONLY ANOMALY         2:  ANOMALY AND DATA DOMAIN'
     */'  3: ONLY DATA DOMAIN     0:  NEITHER NOR'/)
 1110 FORMAT (/' ***DIMENSION IDX(',I3,') REQUIRED***'/)
      END
C
C     ***************************************************************
C
      SUBROUTINE PRINTN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
C     SUBROUTINE COMPUTES AND PRINTS NORMAL PLANE WAVE RESPONSE     C
C     AT THE DATA LEVEL ZD AND AT THE SURFACE                       C
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      COMPLEX CC,C0,FT,FTD,CLOG
      COMMON /PARAM/F,U,DX,NITER,ITERD,EPS,SOR,NH,ZA,ZD
      DATA PI/3.141592654/
      WRITE (6,100)
      WRITE (6,110) F
      U=0.
      CALL KERNLN2(FT,FTD,C0)
C
C     PLANE WAVE RESPONSE AT DATA LEVEL ZD:
C     -------------------------------------
      WRITE (6,120) ZD
      CC=-FT/FTD
      RA=8.E-7*PI*PI*F*CABS(CC)**2
      PH=90.+AIMAG(CLOG(CC))*180./PI
      WRITE (6,130) CC
      WRITE (6,140) RA
      WRITE (6,150) PH
      IF (ZD.EQ.0.) RETURN
C
C     PLANE WAVE RESPONSE AT THE SURFACE:
C     -----------------------------------
      WRITE (6,160)
      RA=8.E-7*PI*PI*F*CABS(C0)**2
      PH=90.+AIMAG(CLOG(C0))*180./PI
      WRITE (6,130) C0
      WRITE (6,140) RA
      WRITE (6,150) PH
      RETURN
  100 FORMAT ('  FREQUENCY AND NORMAL PLANE WAVE RESPONSE'/2X,40(1H=))
  110 FORMAT ('  FREQUENCY F= ',1PE11.4,' Hz')
  120 FORMAT (/'  NORMAL PLANE WAVE RESPONSE AT THE DATA LEVEL ZD=',
     *F8.1,' m:')
  130 FORMAT (4X,'  MODIFIED PENETRATION DEPTH C = (',1P2E10.3,') m')
  140 FORMAT (4X,'  APPARENT RESISTIVITY RA = ',1PE10.3,' Ohm*m')
  150 FORMAT (4X,'  PHASE PH = ', F7.1,' deg')
  160 FORMAT (/'  NORMAL PLANE WAVE RESPONSE AT THE SURFACE:')
      END
C
C     ***************************************************************
C
      SUBROUTINE PRINTQ
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
C     SUBROUTINE PRINTS PARAMETERS OF THE SOURCE FIELD              C
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      COMMON /SOURCE/B0,THETA,BGRAD(43460)
      WRITE (6,100)
      WRITE (6,110) B0,THETA/10
      RETURN
  100 FORMAT (/'  SOURCE FIELD'/2X,12(1H=))
  110 FORMAT ('  QUASI-UNIFORM MAGNETIC SURFACE FIELD WITH B0 =',
     *1PE8.1,' TESLA'/'  ALONG AZIMUTH THETA = ',F4.0,' deg')
      END
C
C     ***************************************************************
C
      SUBROUTINE SPLINE1(N,Y,C)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
C     CUBIC SPLINE INTERPOLATION FOR N EQUIDISTANT ORDINATES Y(1),  C
C     ...,Y(N) WITH BOUNDARY CONDITIONS Y"(1)=Y"(N)=0 AT THE END-   C
C     POINTS. THE OUTPUT ARRAY C(3,N) CONTAINS THE INTERPOLATION    C
C     COEFFICIENTS TO BE USED IN THE SUBROUTINE SPLINE2.            C
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      COMPLEX Y(N),C(3,N),P
      N1=N-1
      DO 10 I=2,N1
   10 C(1,I)=Y(I+1)-2.*Y(I)+Y(I-1)
C
      C(2,1)=0.
      C(3,1)=0.
      DO 20 I=2,N1
      P=4.+C(2,I-1)
      C(2,I)=-1./P
   20 C(3,I)=(C(1,I)-C(3,I-1))/P
C
      C(1,N)=0.
      DO 30 II=2,N1
      I=N+1-II
   30 C(1,I)=C(2,I)*C(1,I+1)+C(3,I)
      C(1,1)=0.
C
      DO 40 I=1,N1
      C(2,I)=Y(I+1)-Y(I)-C(1,I+1)+C(1,I)
   40 C(3,I)=Y(I)-C(1,I)
      C(3,N)=Y(N)
      RETURN
      END
C
C     ***************************************************************
C
      SUBROUTINE SPLINE2(N,XINT,X1,X2,C,YINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
C     THIS SUBROUTINE CALCULATES FOR THE ABSCISSA XINT THE FUNC-    C
C     TION YINT FROM THE COEFFICIENTS C(3,N) DETERMINED IN SPLINE1. C
C     X1 AND X2 ARE THE ABSCISSAE OF Y(1) AND Y(N), RESPECTIVELY.   C
C     OUTSIDE THIS RANGE YINT IS LINEARLY EXTRAPOLATED.             C
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      COMPLEX C(3,N),YINT
      H=(X2-X1)/FLOAT(N-1)
      IF (XINT.LT.X1) GOTO 10
      IF (XINT.GE.X2) GOTO 20
C
      U=(XINT-X1)/H
      I=1+INT(U)
      P=U-I+1
      Q=1.-P
      YINT=C(1,I)*Q**3+C(1,I+1)*P**3+C(2,I)*P+C(3,I)
      RETURN
C
   10 P=(XINT-X1)/H
      YINT=C(2,1)*P+C(3,1)
      RETURN
C
   20 P=(XINT-X2)/H
      YINT=C(2,N-1)*P+C(3,N)
      RETURN
      END
C
C     ***************************************************************
C
      COMPLEX FUNCTION CEXPC(Z)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
C     THE FUNCTIONS CEXPC AND EXPC CORRESPOND TO THE ORDINARY EX-   C
C     PONENTIAL FUNCTIONS CEXP AND EXP, BUT AVOID UNDERFLOW AND     C
C     OVERFLOW BY LIMITING THE (REAL PART OF THE) ARGUMENT WITHIN   C
C     -XL < X < XL, WHERE XL IS PRESCRIBED.                         C
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      COMPLEX Z,Z1,CEXP
      DATA XL/50./
      Z1=Z
      X1=REAL(Z1)
      IF (ABS(X1).GT.XL) Z1=SIGN(XL,X1)
      CEXPC=CEXP(Z1)
      RETURN
      END
      REAL FUNCTION EXPC(X)
      DATA XL/50./
      X1=X
      IF (ABS(X1).GT.XL) X1=SIGN(XL,X1)
      EXPC=EXP(X1)
      RETURN
      END
C
C     ****************************************************************
C
      SUBROUTINE FILTER
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
C     FILTER COEFFICIENTS HR(K,I) FOR FAST HANKEL TRANSFORM         C
C     K=1: J0, K=2: J1                                              C
C     NC=NUMBER OF COEFFICIENTS                                     C
C     COEFFICIENT HR(K,NC0) REFERS TO ZERO ARGUMENT                 C
C     10 DATA POINTS PER DECADE                                     C
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DIMENSION HJ01(48),HJ02(52),HJ11(48),HJ12(52)
      COMMON /HANKEL/NC,NC0,HR(2,100)
      DATA HJ01/
     * 2.89878288E-07, 3.64935144E-07, 4.59426126E-07, 5.78383226E-07,
     * 7.28141338E-07, 9.16675639E-07, 1.15402625E-06, 1.45283298E-06,
     * 1.82900834E-06, 2.30258511E-06, 2.89878286E-06, 3.64935148E-06,
     * 4.59426119E-06, 5.78383236E-06, 7.28141322E-06, 9.16675664E-06,
     * 1.15402621E-05, 1.45283305E-05, 1.82900824E-05, 2.30258527E-05,
     * 2.89878259E-05, 3.64935186E-05, 4.59426051E-05, 5.78383329E-05,
     * 7.28141144E-05, 9.16675882E-05, 1.15402573E-04, 1.45283354E-04,
     * 1.82900694E-04, 2.30258630E-04, 2.89877891E-04, 3.64935362E-04,
     * 4.59424960E-04, 5.78383437E-04, 7.28137738E-04, 9.16674828E-04,
     * 1.15401453E-03, 1.45282561E-03, 1.82896826E-03, 2.30254535E-03,
     * 2.89863979E-03, 3.64916703E-03, 4.59373308E-03, 5.78303238E-03,
     * 7.27941497E-03, 9.16340705E-03, 1.15325691E-02, 1.45145832E-02/
      DATA HJ02/
     * 1.82601199E-02, 2.29701042E-02, 2.88702619E-02, 3.62691810E-02,
     * 4.54794031E-02, 5.69408192E-02, 7.09873072E-02, 8.80995426E-02,
     * 1.08223889E-01, 1.31250483E-01, 1.55055715E-01, 1.76371506E-01,
     * 1.85627738E-01, 1.69778044E-01, 1.03405245E-01,-3.02583233E-02,
     *-2.27574393E-01,-3.62173217E-01,-2.05500446E-01, 3.37394873E-01,
     * 3.17689897E-01,-5.13762160E-01, 3.09130264E-01,-1.26757592E-01,
     * 4.61967890E-02,-1.80968674E-02, 8.35426050E-03,-4.47368304E-03,
     * 2.61974783E-03,-1.60171357E-03, 9.97717882E-04,-6.26275815E-04,
     * 3.94338818E-04,-2.48606354E-04, 1.56808604E-04,-9.89266288E-05,
     * 6.24152398E-05,-3.93805393E-05, 2.48472358E-05,-1.56774945E-05,
     * 9.89181741E-06,-6.24131160E-06, 3.93800058E-06,-2.48471018E-06,
     * 1.56774609E-06,-9.89180896E-07, 6.24130948E-07,-3.93800005E-07,
     * 2.48471005E-07,-1.56774605E-07, 9.89180888E-08,-6.24130946E-08/
      DATA HJ11/
     * 1.84909557E-13, 2.85321327E-13, 4.64471808E-13, 7.16694771E-13,
     * 1.16670043E-12, 1.80025587E-12, 2.93061898E-12, 4.52203829E-12,
     * 7.36138206E-12, 1.13588466E-11, 1.84909557E-11, 2.85321327E-11,
     * 4.64471808E-11, 7.16694771E-11, 1.16670043E-10, 1.80025587E-10,
     * 2.93061898E-10, 4.52203829E-10, 7.36138206E-10, 1.13588466E-09,
     * 1.84909557E-09, 2.85321326E-09, 4.64471806E-09, 7.16694765E-09,
     * 1.16670042E-08, 1.80025583E-08, 2.93061889E-08, 4.52203807E-08,
     * 7.36138149E-08, 1.13588452E-07, 1.84909521E-07, 2.85321237E-07,
     * 4.64471580E-07, 7.16694198E-07, 1.16669899E-06, 1.80025226E-06,
     * 2.93060990E-06, 4.52201549E-06, 7.36132477E-06, 1.13587027E-05,
     * 1.84905942E-05, 2.85312247E-05, 4.64449000E-05, 7.16637480E-05,
     * 1.16655653E-04, 1.79989440E-04, 2.92971106E-04, 4.51975783E-04/
      DATA HJ12/
     * 7.35565435E-04, 1.13444615E-03, 1.84548306E-03, 2.84414257E-03,
     * 4.62194743E-03, 7.10980590E-03, 1.15236911E-02, 1.76434485E-02,
     * 2.84076233E-02, 4.29770596E-02, 6.80332569E-02, 9.97845929E-02,
     * 1.51070544E-01, 2.03540581E-01, 2.71235377E-01, 2.76073871E-01,
     * 2.16691977E-01,-7.83723737E-02,-3.40675627E-01,-3.60693673E-01,
     * 5.13024526E-01,-5.94724729E-02,-1.95117123E-01, 1.99235600E-01,
     *-1.38521553E-01, 8.79320859E-02,-5.50697146E-02, 3.45637848E-02,
     *-2.17527180E-02, 1.37100291E-02,-8.64656417E-03, 5.45462758E-03,
     *-3.44138864E-03, 2.17130686E-03,-1.36998628E-03, 8.64398952E-04,
     *-5.45397874E-04, 3.44122545E-04,-2.17126585E-04, 1.36997597E-04,
     *-8.64396364E-05, 5.45397224E-05,-3.44122382E-05, 2.17126544E-05,
     *-1.36997587E-05, 8.64396338E-06,-5.45397218E-06, 3.44122380E-06,
     *-2.17126543E-06, 1.36997587E-06,-8.64396337E-07, 5.45397218E-07/
C
      NC=100
      NC0=60
      DO 10 I=1,48
      HR(1,I)=HJ01(I)
   10 HR(2,I)=HJ11(I)
      DO 20 I=1,52
      HR(1,I+48)=HJ02(I)
   20 HR(2,I+48)=HJ12(I)
      RETURN
      END
C
C     ***************************************************************
C
      SUBROUTINE FRECH(IX0,IY0,TAUA,GXX,GXY,BXX,BXY,BXZ,ES,ZZ,IWD)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
C     THIS SUBROUTINE DETERMINES THE MAGNETIC FIELD FRECHET DERIVA- C
C     TIVES BY SOLVING FOR EACH COMPONENT AN INTEGRAL EQUATION SIMI-C
C     LAR TO THE ELECTRIC FIELD INTEGRAL EQUATION IN SUBROUTINE SEI-C
C     DEL. FOR THE POINT OF OBSERVATION (IX0,IY0), 1.LE.IX0.LE.NX2, C
C     1.LE.IY0.LE.NY2, THE DERIVATIVES ARE OBTAINED FOR EACH CELL OFC
C     THE ANOMALOUS DOMAIN AND ARE STORED IN THE ARRAY DB(IK,IC),   C
C     IK=1,..., NX*NY, IC=1,...,3 (IC=1: BX, IC=2: BY, IC=3: BZ).   C
C     FOR IWD=1 INTERMEDIATE ITERATIONS OF THE ITERATIVE INTEGRAL   C
C     EQUATION SOLUTION AND THE FINAL SOLUTION ARE DISPLAYED (cf.   C
C     SUBROUTINE SEIDEL).                                           C
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      PARAMETER (NXYD=200)
      REAL TAUA(1)
      COMPLEX GXX(1),GXY(1),BXX(1),BXY(1),BXZ(1),ES(4,1),ZZ(6,1)
      COMPLEX FXX,FXY,FYX,FYY,FXZ,FYZ
      COMPLEX FX(NXYD),FY(NXYD),FNX(NXYD),FNY(NXYD),DAB(6,NXYD)
      COMPLEX AXY,CX,CY,SX,SY,C,CTAU
      COMPLEX BX1,BY1,BZ1,BX2,BY2,BZ2,DETX,DETY,DETZ
      COMPLEX DBX1,DBY1,DBZ1,DBX2,DBY2,DBZ2,DA,DB,A,B,R1,R2
      COMMON /PARAM/F,U,DX,NITER,ITERD,EPS,SOR,NH,ZA,ZD
      COMMON /NUMBER/NX,NY,NX2,NY2,NX0,NY0,N1,N2
      DATA PI/3.141592654/
C
      NXY=NX*NY
      C=(0.,8.E-7)*PI*PI*F
C
C     TEST OF DIMENSIONS:
      IF (NXY.LE.NXYD) GOTO 5
      WRITE (6,1040) NXYD,NXY
      STOP
    5 CONTINUE

C
      DO 90 IC=1,3
      DO 10 I=1,NX
      DO 10 K=1,NY
      IK=(K-1)*NX+I
      CTAU=C*TAUA(IK)
      L=IX0
      M=IY0
      LM=(M-1)*NX2+L
      LI=L-I+NX0
      MK=M-K+NY0
      IDX=IABS(MK)*N2+IABS(LI)+1
      IDY=IABS(LI)*N2+IABS(MK)+1
      SX=ISIGN(1,LI)
      SY=ISIGN(1,MK)
      SXY=SX*SY
      FXX=BXX(IDX)*SXY
      FXY=BXY(IDX)
      FXZ=BXZ(IDX)*SY
      FYX=-BXY(IDY)
      FYY=-BXX(IDY)*SXY
      FYZ=-BXZ(IDY)*SX
      IF (IC.EQ.1) FNX(IK)=FXX
      IF (IC.EQ.1) FNY(IK)=FYX
      IF (IC.EQ.2) FNX(IK)=FXY
      IF (IC.EQ.2) FNY(IK)=FYY
      IF (IC.EQ.3) FNX(IK)=FXZ
      IF (IC.EQ.3) FNY(IK)=FYZ
      FX(IK)=FNX(IK)
   10 FY(IK)=FNY(IK)
C
C     START OF ITERATIONS:
C     -------------------
      DO 40 ITER=1,NITER
      S1=0.
      S2=0.
      DO 30 I=1,NX
      DO 30 K=1,NY
      IK=(K-1)*NX+I
      CTAU=C*TAUA(IK)
      AXY=1.+CTAU*GXX(1)
      CX=FX(IK)-FNX(IK)
      CY=FY(IK)-FNY(IK)
C
      DO 20 L=1,NX
      DO 20 M=1,NY
      LM=(M-1)*NX+L
      IF (ABS(TAUA(LM)).LT.1.E-5) GOTO 20
      IL=IABS(I-L)
      KM=IABS(K-M)
      IDX=KM*N1+IL+1
      IDY=IL*N1+KM+1
      FXX=GXX(IDX)
      FXY=GXY(IDX)*ISIGN(1,(I-L)*(K-M))
      FYX=FXY
      FYY=GXX(IDY)
      CTAU=C*TAUA(LM)
      CX=CX+CTAU*(FX(LM)*FXX+FY(LM)*FYX)
      CY=CY+CTAU*(FX(LM)*FXY+FY(LM)*FYY)
   20 CONTINUE
      SX=-SOR
      SY=-SOR*CY/AXY
      FX(IK)=FX(IK)+SX
      FY(IK)=FY(IK)+SY
      S1=S1+CABS(SX)**2+CABS(SY)**2
   30 S2=S2+CABS(FX(IK))**2+CABS(FY(IK))**2
      EPS1=SQRT(S1/S2)
      IF (EPS1.LT.EPS) GOTO 50
C
C     PRINTOUT AFTER ITERD ITERATIONS:
C     -------------------------------
      IF (MOD(ITER,ITERD).NE.0.OR.IWD.EQ.0) GOTO 40
      WRITE (6,1010) ITER,EPS1
      CALL PRINTF(FX,NX,NY,19)
      CALL PRINTF(FY,NX,NY,20)
   40 CONTINUE
C     END OF ITERATIONS
C
C     FAILURE EXIT: NO CONVERGENCE WITHIN NITER ITERATIONS
C     ----------------------------------------------------
      WRITE (6,1020) EPS,NITER,EPS1
      GOTO 60
C
C     REGULAR EXIT: CONVERGENCE AFTER NOT MORE THAN NITER ITERATIONS
C     --------------------------------------------------------------
   50 CONTINUE
C     WRITE (6,1030) ITER,EPS
C
C     OPTIONAL PRINT:
C     ---------------
   60 IF (IWD.NE.1) GOTO 70
      CALL PRINTF(FX,NX,NY,19)
      CALL PRINTF(FY,NX,NY,20)
C
C     CALCULATE FRECHET DERIVATIVES OF ALL THREE MAGNETIC FIELD COMPONENTS
C     --------------------------------------------------------------------
C     FOR BOTH POLARIZATIONS AT POSITION (IX0,IY0):
C     ---------------------------------------------
   70 DO 80 IK=1,NXY
      DO 80 IPOL=1,2
      I2=2*(IPOL-1)
      I3=3*(IPOL-1)+IC
   80 DAB(I3,IK)=-C*(ES(I2+1,IK)*FX(IK)+ES(I2+2,IK)*FY(IK))
   90 CONTINUE
C
      IK0=(IY0-1)*NX2+IX0
      BX1=ZZ(1,IK0)
      BY1=ZZ(2,IK0)
      BZ1=ZZ(3,IK0)
      BX2=ZZ(4,IK0)
      BY2=ZZ(5,IK0)
      BZ2=ZZ(6,IK0)
      DETX=BY1*BZ2-BY2*BZ1
      DETY=BZ1*BX2-BZ2*BX1
      DETZ=BX1*BY2-BX2*BY1
      A=-DETX/DETZ
      B=-DETY/DETZ
      DO 100 IK=1,NXY
      DBX1=DAB(1,IK)
      DBY1=DAB(2,IK)
      DBZ1=DAB(3,IK)
      DBX2=DAB(4,IK)
      DBY2=DAB(5,IK)
      DBZ2=DAB(6,IK)
      R1=DBZ1-A*DBX1-B*DBY1
      R2=DBZ2-A*DBX2-B*DBY2
      DA=(R1*BY2-R2*BY1)/DETZ
      DB=(R2*BX1-R1*BX2)/DETZ
      DAB(1,IK)=DA
  100 DAB(2,IK)=DB
C
      WRITE (6,1000) IX0,IY0
      DO 120 IC=1,2
      DO 110 IK=1,NXY
  110 FNX(IK)=DAB(IC,IK)
  120 CALL PRINTF(FNX,NX,NY,20+IC)
      RETURN
C
 1000 FORMAT(/'  SENSITIVITY OF THE TIPPER COMPONENTS A AND B',
     #/'  AT THE POINT IX0=',I2,', IY0=',I2,' OF THE DATA DOMAIN'
     #/'  TO CONDUCTANCE CHANGES AT EACH OF THE NX*NY POINTS'
     #/'  IN THE ANOMALOUS DOMAIN. GIVEN IS DA:=dA[IX0,IY0]/dS[i,k] '
     #/'  AND DB:=dB[IX0,IY0]/dS[I,K] IN UNITS OF 1/S (= Ohm).'/)

 1010 FORMAT (/'  ITERATION  ',I3/
     *'  RELATIVE MEAN DEVIATION OF SUCCESSIVE ITERATIONS:',1PE10.3/)
 1020 FORMAT (//'  ***THRESHOLD EPS =',1PE10.3,' NOT REACHED AFTER '
     *,I3,' ITERATIONS'/'  ACTUAL VALUE: ',1PE10.3/)
 1030 FORMAT(/'  CONVERGENCE AFTER ',I3,' ITERATIONS (EPS=',
     *1PE10.3,')'/)
 1040 FORMAT (/'  INCREASE PARAMETER NXYD IN SUBROUTINE FRECH FROM ',
     *I3,' TO AT LEAST ',I4,'***STOP')
      END
C
C     ********************************************************



