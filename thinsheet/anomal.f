      SUBROUTINE ANOMAL(TAUN,TAU,TAUA,CMODELFILE)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
C     THIS SUBROUTINE DEFINES THE CONDUCTANCE TAU (S) IN THE ANO-   C
C     MALOUS DOMAIN. IT IS STORED AS ONE-DIMENSIONAL ARRAY TAU(IK), C
C     IK=(K-1)*NX+I, I=1,..,NX, K=1,..,NY.  ALSO DETERMINED IS THE  C
C     ANOMALOUS CONDUCTANCE TAUA(IK)=TAU(IK)-TAUN.                  C
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DIMENSION TAU(NX*NY),TAUA(NX*NY)
      CHARACTER*32 CMODELFILE
      COMMON /NUMBER/NX,NY,NX2,NY2,NX0,NY0,N1,N2
C
C     DEFINITION OF CONDUCTANCE IN THE ANOMALOUS DOMAIN
C     -------------------------------------------------
c      TAU1=100.
c
      NXY=NX*NY
c
c      DO 10 IK=1,NXY
c   10 TAU(IK)=TAUN
c      TAU(6)=TAU1
c      TAU(7)=TAU1
c      TAU(13)=TAU1
c      TAU(14)=TAU1
c      DO 20 I=1,4
c      DO 20 K=2,3
c      IK=(K-1)*NX+I
c   20 TAU(IK)=TAU1
c      TAU(22)=TAU1
C


C  See /users/storm/work/thinsheet/taumod for alternative models
C  The Allan McKay model 3: PhD thesis Figures 6.9-6.11 for info
C
       OPEN(50,FILE='condmodels/'//CMODELFILE) 

c      DO i = 1,NXY
c         READ(50,*) dummy, dummy, TAU(i)
c      END DO
c      CLOSE(50)

c     see notes 2/08/01
c     increasing the size of the anomalous domain from that contained in the model file
c     NXMODEL=121; NYNODEL=169 in model file (NX geograph-north)
c     Increase size of model 10 Grid points in each direction: NX=131, NY=179

       NXMODEL=NX+NX*0.5
       NYMODEL=NY+NX*0.5
       DEX=(NX-NXMODEL)/2
       DEY=(NY-NYMODEL)/2
       NXT=NX-DEX
       NYT=NY-DEY
       TAUN=0.1

       do k=1,NY
          DO i=1,NX
            IK=(K-1)*NX+I
            IF((k.GT.DEY).AND.(k.LE.NYT).AND.
     1        (i.GT.DEX).AND.(i.LE.NXT)) THEN
                READ(50,*) dummy, dummy, TAU(IK)
c                write(*,*) "1", k,i,ik,TAU(IK)
            ELSE
c                write(*,*) "2", k,i,ik,TAU(IK)
                TAU(IK)=TAUN ! =600, see edge of the grid that is read in
            END IF
          END DO
       END DO

       CLOSE(50)
c       do k=1,NY
c          DO i=1,NX
c            IK=(K-1)*NX+I
c            write(60,*) k,i,TAU(IK)
c         END Do
c      END DO           

c       CLOSE(60)

      write(*,*) ' PROG ANOMAL: NX*NY = ',NXY 
c--------------------------------
Cc this is the small square model
C       DO i=1,NXY
C       TAU(i)=600
C       END DO
C
Cc y-component
C           do k=25,35
Cc x-component
C             do i=30,40 
C
C              IK=(K-1)*NX+I
C              TAU(IK)=50
C           end do
C       end do
c--------------------------------

       write(*,*) 
     1 ' PROG ANOMAL: Normal conduct = ',TAUN

C     ANOMALOUS CONDUCTANCE:
C     ----------------------
      DO 30 IK=1,NXY
   30 TAUA(IK)=TAU(IK)-TAUN
      RETURN
      END
