
      subroutine NEWMODEL(NMOD)

c     open input file containing resist model
c     read it in
      
      COMMON /RESIST/NL,RHO(20),D(20),NL1,RHO1(22),D1(22),NLA,NLD,
     1     TAUN

      read(33,*) NL
      read(33,*) TAUN

      IF (NMOD.NE.99) THEN
         OPEN(2,FILE='condmodels/layers_austria.txt',STATUS='OLD')
      ELSE IF (NMOD.EQ.99) THEN
         OPEN(2,FILE='condmodels/layers_quebec.txt',STATUS='OLD')
      END IF

c     --------------------------------
c     Adapted script for Austria, RLB 2015-10-14 (source Adam et al. (2012))
      IF (NMOD.EQ.3) THEN
      write(*,*) 'EURHOM conductivity model #3 (Vienna basin)...'
      DO i=1,4
         READ(2,*) rdummy, D(i), RHO(i)
      END DO
           
      DO i=5,NL-1
         READ(2,*) rdummy, D(i), RHO(i)
      END DO

      READ(2,*) rdummy, rdummy, RHO(NL)
      END IF
c     -----------

      IF (NMOD.EQ.16) THEN
      write(*,*) 'EURHOM conductivity model #16 (sub-alps)...'
      DO i=1,4
         READ(2,*) rdummy, rdummy, rdummy, D(i), RHO(i)
      END DO
           
      DO i=5,NL-1
         READ(2,*) rdummy, D(i), RHO(i)
      END DO

      READ(2,*) rdummy, rdummy, RHO(NL)
      END IF
c     -----------

      IF (NMOD.EQ.39) THEN
      write(*,*) 'EURHOM conductivity model #39 (plains)...'
      DO i=1,4
         READ(2,*) rdummy, rdummy, rdummy, rdummy, rdummy, D(i), RHO(i)
      END DO
           
      DO i=5,NL-1
         READ(2,*) rdummy, D(i), RHO(i)
      END DO

      READ(2,*) rdummy, rdummy, RHO(NL)
      END IF
c     -----------

      IF (NMOD.EQ.55) THEN
      write(*,*) 'EURHOM conductivity model #55 (alps)...'
      DO i=1,4
         READ(2,*) rdummy, rdummy, rdummy, rdummy, rdummy, rdummy, 
     1     rdummy, D(i), RHO(i)
      END DO
           
      DO i=5,NL-1
         READ(2,*) rdummy, D(i), RHO(i)
      END DO

      READ(2,*) rdummy, rdummy, RHO(NL)
      END IF
c     -----------

      IF (NMOD.EQ.99) THEN
      write(*,*) 'NERC Quebec conductivity model (#99)...'
      DO i=1,4
         READ(2,*) rdummy, D(i), RHO(i)
      END DO
           
      DO i=5,NL-1
         READ(2,*) rdummy, D(i), RHO(i)
      END DO

      READ(2,*) rdummy, rdummy, RHO(NL)
      END IF
c     --------------------------------


      RETURN
      END subroutine NEWMODEL
