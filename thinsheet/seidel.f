      SUBROUTINE SEIDEL(TAUA,GXX,GXY,ENX,ENY,ESX,ESY,IWS,NXD,NYD)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
C     THIS SUBROUTINE DETERMINES THE SOLUTION OF THE INTEGRAL EQUA- C
C     TION BY GAUSS-SEIDEL ITERATION. THE CONVERGENCE MIGHT BE AC-  C
C     CELERATED BY SUCCESSIVE OVERRELAXATION, USING A PRESCRIBED    C
C     OVERRELAXATION PARAMETER SOR (1<SOR<1.5, POSSIBLY ALSO COM-   C
C     PLEX).A MAXIMUM OF NITER ITERATIONS IS USED. THE INTERMEDIATE C
C     RESULT IS PRINTED AFTER ITERD ITERATIONS. THE ITERATION STOPS C
C     IF THE RELATIVE DEVIATION OF TWO SUCCESSIVE ITERATIONS IS BE- C
C     LOW A GIVEN THRESHOLD EPS (TYPICALLY 1.E-5). THE ELEMENTS  OF C
C     THE SYSTEM MATRIX ARE CONSTRUCTED AT EACH ITERATION FROM THE  C
C     ELEMENTS OF GXX AND GXY, DETERMINED IN GREEN1.                C
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      REAL TAUA(NX*NY)
      COMPLEX GXX(N1*N1),GXY(N1*N1),ENX(NX*NY),ENY(NX*NY)
      COMPLEX ESX(NXD*NYD),ESY(NXD*NYD)
      COMPLEX FXX,FXY,FYX,FYY
      COMPLEX AXY,CX,CY,SX,SY,C,CTAU
      COMMON /PARAM/F,U,DX,NITER,ITERD,EPS,SOR,NH,ZA,ZD
      COMMON /NUMBER/NX,NY,NX2,NY2,NX0,NY0,N1,N2
      
      DATA PI/3.141592654/
C
      IF (IWS.EQ.1) WRITE (6,1000)
c     iwmu_0
      C=(0.,8.E-7)*PI*PI*F
      NXY=NX*NY
      

c     set electric field to normal field-the first 'guess'
      DO 10 IK=1,NXY
      ESX(IK)=ENX(IK)
 10   ESY(IK) = ENY(IK)

c     Alternatively Open old electrci field file 
c     and read in the soluttion
c     read in a previously solution to speed up the 
c      process
c      write(*,*) 'Reading in electric field'
c      open(50,FILE='/scratch/allan/thin/modelruns/esx_0.dat')  
c      open(60,FILE='/scratch/allan/thin/modelruns/esy_0.dat')  
c      do k=1,NY
c          DO i=1,NX
c            IK=(K-1)*NX+I
c                READ(50,*) ESX(IK)
c                READ(60,*) ESY(IK)
c          END DO
c       END DO
c       close(50)
c       close(60)
c       write(*,*) 'Finished reading in electric field'

C     START OF ITERATIONS:
      DO 40 ITER=1,NITER
         write(*,*) 'Iteration, max interations, eps1:',ITER,NITER,EPS1
      S1=0.
      S2=0.

c     for each grid point (Do 30)- sum contrimbution from all other nodes (DO 20)

      DO 30 I=1,NX
      DO 30 K=1,NY
      IK=(K-1)*NX+I
      CTAU=C*TAUA(IK)
      AXY=1.+CTAU*GXX(1)

c     Want to solve for the Anomalous Field =total efield - normal field
      CX=ESX(IK)-ENX(IK)
      CY=ESY(IK)-ENY(IK)

c     Integrations to compute distorted E-field at the K,I grid point
      DO 20 L=1,NX
      DO 20 M=1,NY
      LM=(M-1)*NX+L
      IF (ABS(TAUA(LM)).LT.1.E-1) GOTO 20
c     Distance between 'source point' and 'grid point' in x and y
      IL=IABS(I-L)
      KM=IABS(K-M)
      IDX=KM*N1+IL+1
      IDY=IL*N1+KM+1
c     the appropriate greens vector-surely could do this once for the whole grid and store
      FXX=GXX(IDX)
      FXY=GXY(IDX)*ISIGN(1,(I-L)*(K-M))
      FYX=FXY
      FYY=GXX(IDY)
      CTAU=C*TAUA(LM)
      CX=CX+CTAU*(ESX(LM)*FXX+ESY(LM)*FYX)
      CY=CY+CTAU*(ESX(LM)*FXY+ESY(LM)*FYY)
   20 CONTINUE

      SX=-SOR*CX/AXY
      SY=-SOR*CY/AXY
      ESX(IK)=ESX(IK)+SX
      ESY(IK)=ESY(IK)+SY
      S1=S1+CABS(SX)**2+CABS(SY)**2
   30 S2=S2+CABS(ESX(IK))**2+CABS(ESY(IK))**2
      EPS1=SQRT(S1/S2)
      IF (EPS1.LT.EPS) GOTO 50
C
C     PRINTOUT AFTER ITERD ITERATIONS:
      IF (MOD(ITER,ITERD).NE.0.OR.IWS.EQ.0) GOTO 40
      WRITE (6,1010) ITER,EPS1
      CALL PRINTF(ESX,NX,NY,5)
      CALL PRINTF(ESY,NX,NY,6)
   40 CONTINUE
c     END OF ITERATIONS
C
C     FAILURE EXIT: NO CONVERGENCE WITHIN NITER ITERATIONS
C     ----------------------------------------------------
      WRITE (6,1020) EPS,NITER,EPS1
      GOTO 60
C
C     REGULAR EXIT: CONVERGENCE AFTER NOT MORE THAN NITER ITERATIONS
C     --------------------------------------------------------------
   50 WRITE (6,1030) ITER,EPS
C
C     OPTIONAL PRINT:
C     ---------------
   60 IF (IWS.NE.1) RETURN
      CALL PRINTF(ESX,NX,NY,5)
      CALL PRINTF(ESY,NX,NY,6)
      RETURN
 1000 FORMAT(/'  ELECTRIC FIELD BY GAUSS-SEIDEL ITERATION'/
     *2X,40(1H=))
 1010 FORMAT (/'  ITERATION  ',I3/
     *'  RELATIVE MEAN DEVIATION OF SUCCESSIVE ITERATIONS:',1PE10.3/)
 1020 FORMAT (//'  ***THRESHOLD EPS =',1PE10.3,' NOT REACHED AFTER '
     *,I3,' ITERATIONS'/'  ACTUAL VALUE: ',1PE10.3/)
 1030 FORMAT(/'  SEIDEL: CONVERGENCE AFTER ',I3,' ITERATIONS (EPS=',
     *1PE10.3,')'/)
      END
C
C     ***************************************************************
