SUBROUTINE seidel_new(TAUA,GXX,GXY,ENX,ENY,ESX,ESY,IWS,NXD,NYD)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C                                                                   C
!C     THIS SUBROUTINE DETERMINES THE SOLUTION OF THE INTEGRAL EQUA- C
!C     TION BY GAUSS SEIDEL ITERATION.                               C
!C                                                                   C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C

      REAL TAUA(NX*NY)
      COMPLEX GXX(N1*N1),GXY(N1*N1),ENX(NX*NY),ENY(NX*NY)
      COMPLEX ESX(NXD*NYD),ESY(NXD*NYD)
      COMPLEX FXX,FXY,FYX,FYY
      COMPLEX AXY,CX,CY,SX,SY,C,CTAU
      REAL,SAVE,DIMENSION(43460) :: n_zerotaua
      REAL,SAVE,DIMENSION(43460) :: x_cord
      REAL,SAVE,DIMENSION(43460) :: y_cord
      REAL :: rmax

      COMMON /PARAM/F,U,DX,NITER,ITERD,EPS,SOR,NH,ZA,ZD
      COMMON /NUMBER/NX,NY,NX2,NY2,NX0,NY0,N1,N2
      
      DATA PI/3.141592654/
      DATA rmax/20/
!
      IF (IWS.EQ.1) WRITE (6,1000)
!     iwmu_0
      C=(0.,8.E-7)*PI*PI*F
      NXY=NX*NY
      
!     set electric field to normal field-the first 'guess'
      DO 10 IK=1,NXY
      ESX(IK)=ENX(IK)
 10   ESY(IK) = ENY(IK)

!     Alternatively Open old electrci field file 
!     and read in the soluttion
!     read in a previously solution to speed up the 
!      process
!      write(*,*) 'Reading in electric field'
!      open(50,FILE='/scratch/allan/thin/modelruns/esx_0.dat')  
!      open(60,FILE='/scratch/allan/thin/modelruns/esy_0.dat')  
!      do k=1,NY
!          DO i=1,NX
!            IK=(K-1)*NX+I
!                READ(50,*) ESX(IK)
!                READ(60,*) ESY(IK)
!          END DO
!       END DO
!       close(50)
!       close(60)
!       write(*,*) 'Finished reading in electric field'
      
      j=1
      DO i=1,NX
         DO k=1,NY
            ik = (k-1)*NX+i
            IF (ABS(TAUA(ik)).GT.1.E-5) THEN
               n_zerotaua(j)=TAUA(ik)
               x_cord(j)=i
               y_cord(j)=k
               j=j+1
            END IF
         END DO
      END DO
      WRITE(*,*) 'Size of Conductance array', j
      jmax=j-1


!     START OF ITERATIONS:
      DO 40 ITER=1,NITER
         write(*,*) ITER,NITER
      S1=0.
      S2=0.

!     for each grid point (Do 30)- sum contrimbution from all other nodes (DO 20)

      DO 30 I=1,NX
      DO 30 K=1,NY
         
         IK=(K-1)*NX+I
         CTAU=C*TAUA(IK)
         AXY=1.+CTAU*GXX(1)

         !     Want to solve for the Anomalous Field =total efield - normal field
         CX=ESX(IK)-ENX(IK)
         CY=ESY(IK)-ENY(IK)

         !     Integrations to compute distorted E-field at the K,I grid point
         !call mxprod(.TRUE.,rmax,f,TAUA,GXX,GXY,ESX,ESY,NXD,NYD,NX,NY,N1,cx,cy)

         ! just include cells which have non-zero conductance
         ! and where rmax < |r-r_o] the distance between source and observer
         DO j=1,jmax
            l=x_cord(j)
            m=y_cord(j)
            LM=(m-1)*NX+l
            ! Distance between 'source point' and 'grid point' in x and y
            IL=IABS(I-L)
            KM=IABS(K-M)
            r_mag = SQRT(REAL((IL)**2+(KM)**2))
            IF (r_mag.LT.rmax) THEN 
               IDX=KM*N1+IL+1
               IDY=IL*N1+KM+1
               ! the appropriate greens vector
               FXX=GXX(IDX)
               FXY=GXY(IDX)*ISIGN(1,(I-L)*(K-M))
               FYX=FXY
               FYY=GXX(IDY)
               CTAU=C*TAUA(LM)
               CX=CX+CTAU*(ESX(LM)*FXX+ESY(LM)*FYX)
               CY=CY+CTAU*(ESX(LM)*FXY+ESY(LM)*FYY)
            ELSE
               CYCLE
            END IF
         END DO

         SX=-SOR*CX/AXY
         SY=-SOR*CY/AXY
         ESX(IK)=ESX(IK)+SX
         ESY(IK)=ESY(IK)+SY
         S1=S1+CABS(SX)**2+CABS(SY)**2
30       S2=S2+CABS(ESX(IK))**2+CABS(ESY(IK))**2
         EPS1=SQRT(S1/S2)
         IF (EPS1.LT.EPS) GOTO 50
!
!     PRINTOUT AFTER ITERD ITERATIONS:
      IF (MOD(ITER,ITERD).NE.0.OR.IWS.EQ.0) GOTO 40
      WRITE (6,1010) ITER,EPS1
      CALL PRINTF(ESX,NX,NY,5)
      CALL PRINTF(ESY,NX,NY,6)
   40 CONTINUE
!     END OF ITERATIONS
!
!     FAILURE EXIT: NO CONVERGENCE WITHIN NITER ITERATIONS
!     ----------------------------------------------------
      WRITE (6,1020) EPS,NITER,EPS1
      GOTO 60
!
!     REGULAR EXIT: CONVERGENCE AFTER NOT MORE THAN NITER ITERATIONS
!     --------------------------------------------------------------
   50 WRITE (6,1030) ITER,EPS
!
!     OPTIONAL PRINT:
!     ---------------
   60 IF (IWS.NE.1) RETURN
      CALL PRINTF(ESX,NX,NY,5)
      CALL PRINTF(ESY,NX,NY,6)
      RETURN
 1000 FORMAT(/'  ELECTRIC FIELD BY GAUSS-SEIDEL ITERATION'/2X,40(1H=))
 1010 FORMAT (/'  ITERATION  ',I3/'  RELATIVE MEAN DEVIATION OF SUCCESSIVE ITERATIONS:',1PE10.3/)
 1020 FORMAT (//'  ***THRESHOLD EPS =',1PE10.3,' NOT REACHED AFTER ',I3,' ITERATIONS'/'  ACTUAL VALUE: ',1PE10.3/)
 1030 FORMAT(/'  SEIDEL: CONVERGENCE AFTER ',I3,' ITERATIONS (EPS=',1PE10.3,')'/)

END SUBROUTINE seidel_new
!C
!C     ***************************************************************
