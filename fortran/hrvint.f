c
C     George C. Wallick 1970, https://doi.org/10.1145/355598.362781
C     ALGORITHM 400 COLLECTED ALGORITHMS FROM ACM.
C     ALGORITHM APPEARED IN COMM. ACM, VOL. 13, NO. 10,
C     P. 622.
C     Original FORTRAN IV modified to make it compatible
C     with gnu/gfortran compilator
C     
C     cnaw@as.arizona.edu
C     2020-09-04
C
      DOUBLE PRECISION FUNCTION HRVINT(F, A, B, MAX, ACC, FAC, MFIN)
C HAVIE INTEGRATION WITH AN EXPANDED RUTISHAUSER-
C TYPE SUMMATION PROCEDURE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION T(17), U(17), TPREV(17), UPREV(17)
      EXTERNAL F
C TEST FOR MAX GREATER THAN 16
      MUX=MAX
      if(max.gt.16) MUX=16
C INITIALIZATION
 10   ENPT=0.5D0*(F(A)+F(B))
      SUMT=0.0D0
      MFIN=1
      N=1
      H=B-A
      SH=H
C BEGIN REPETITIVE LOOP FROM ORDER 1 TO ORDER MAX
 15   T(1)=H*(ENPT+SUMT)
      SUM=0.D0
      NN=N+N
      EN=NN
      EM=SH/EN
C BEGIN RUTISHAUSER EVALUATION OF RECTANGULAR SUMS
C INITIALIZATION
      IF(NN .le. 16) then
         NZ = NN
      else
         NZ=16
      end if
C
      IF(NN.le.256 ) then
         NA = NN
      else
         NA=256
      end if
C
      IF(NN .le. 4096 ) then
         NB=NN
      else
         NB=4096
      endif
C     DEVELOPMENT OF RECTANGULAR SUMS
      do kc = 1, nn, 4096
         sumb = 0.0d0
         kk = kc + nb -1
         do kb = kc, kk, 256
            suma = 0.0d0
            kkk = kb + na -1
            do ka = kb, kkk, 16
               sumz = 0.0d0
               kfr = ka + nz -1
               do kz = ka, kfr, 2
                  zkz = kz
                  sumz = sumz + f(a+zkz*em)
               end do
               suma = sumz + suma
            end do
            sumb = suma + sumb
         end do
         sum = sumb + sum
      end do
C
C     original code
C
C         IF(NN-16) 20, 20, 25
C   20    NZ=NN
C         GO TO 30
C   25    NZ=16
C         IF(NN-256) 30, 30, 35
C   30    NA=NN
C         GO TO 40
C   35    NA=256
C         IF(NN-4096) 40, 40, 45
C   40    NB=NN
C         GO TO 50
C   45    NB=4096
c 50   DO 70 KC=1,NN,4096
c         SUMB=0.0D0
c         KK=KC+NB-1
c         DO 65 KB=KC,KK,256
c            SUMA=0.0D0
c            KKK=KB+NA-1
c            DO 60 KA=KB,KKK,16
c               SUMZ=0.0D0
c               KFR=KA+NZ-1
c               DO 55 KZ=KA,KFR,2
c                  ZKZ=KZ
c 55               SUMZ=SUMZ+F(A+ZKZ*EM)
c 60               SUMA=SUMZ+SUMA
c 65               SUMB=SUMA+SUMB
c 70               SUM=SUMB+SUM
C END OF RUTISHAUSER PROCEDURE
      U(1)=H*SUM
      K=1
C Original code:
C
C BEGIN EXTRAPOLATION LOOP
C 75   FAC=DABS(T(K)-U(K))
CC      IF(T(K)) 80, 85, 80
CC TEST FOR RELATIVE ACCURACY
C 80   CONTINUE
C      IF(FAC-DABS(ACC*T(K))) 90, 90, 100
CC TEST FOR ABSOLUTE ACCURACY WHEN T(K)=0
C 85   continue
C      IF(FAC-DABS(ACC)) 95, 95, 100
C 90   continue
C      FAC=FAC/DABS(T(K))
CC INTEGRAL EVALUATION BEFORE EXIT
C 95   continue
C      HRVINT=0.5D0*(T(K)+U(K))
C      RETURN
C
 75   FAC = DABS(T(K)-U(K))
c      print *, 'k, t(k), u(k)',k, t(k), u(k)
      if(t(k).eq.0.0d0 .AND.fac-dabs(acc).LE.0.D0) then
         hrvint=0.5d0*(t(k)+u(k))
         return
      end if
C
      IF(t(k).ne.0.0d0 .AND.FAC-DABS(ACC*T(K)).LE. 0.0D0) THEN
         FAC = FAC/DABS(T(K))
         hrvint=0.5d0*(t(k)+u(k))
         return
      end if
 100  continue
c      IF(K-MFIN) 105, 115, 115
      IF(K-MFIN .ge.0 ) go to 115
 105  AK=K+K
      D=2.D0**AK
      DMA=D-1.D0
C     BEGIN EXTRAPOLATION
      T(K+1)=(D*T(K)-TPREV(K))/DMA
      TPREV(K)=T(K)
      U(K+1)=(D*U(K)-UPREV(K))/DMA
      UPREV(K)=U(K)
C     END EXTRAPOLATION
      K=K+1
c      IF(K-MUX) 75, 110, 110
      IF(K-MUX .lt.0 ) go to 75
C     END EXTRAPOLATION LOOP
 110  FAC=ABS(T(K)-U(K))
c      IF(T(K)) 90, 95, 90
      if(t(k).eq.0.0d0) then
         hrvint=0.5d0*(t(k)+u(k))
         RETURN
      else
         FAC = FAC/DABS(T(K))
         hrvint=0.5d0*(t(k)+u(k))
         return
      end if
C     ORDER IS INCREASED BY ONE
 115  H=0.5D0*H
      SUMT=SUMT+SUM
      TPREV(K)=T(K)
      UPREV(K)=U(K)
      MFIN=MFIN+1
      N=NN
      GO TO 15
C     RETURN FOR NEXT ORDER EXTRAPOLATION
      END
