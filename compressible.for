      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

      INCLUDE 'ABA_PARAM.INC'

      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)

      DOUBLE PRECISION B(6), DET
      PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, FOUR=4.D0,
     1 SIX=6.D0)

      mu = PROPS(1)
      lam = PROPS(2)

      DET=DFGRD1(1, 1)*DFGRD1(2, 2)*DFGRD1(3, 3)
     1 -DFGRD1(1, 2)*DFGRD1(2, 1)*DFGRD1(3, 3)
      IF(NSHR.EQ.3) THEN
      DET=DET+DFGRD1(1, 2)*DFGRD1(2, 3)*DFGRD1(3, 1)
     1 +DFGRD1(1, 3)*DFGRD1(3, 2)*DFGRD1(2, 1)
     2 -DFGRD1(1, 3)*DFGRD1(3,1)*DFGRD1(2, 2)
     3 -DFGRD1(2, 3)*DFGRD1(3, 2)*DFGRD1(1, 1)
      END IF

      B(1)=DFGRD1(1, 1)**2+DFGRD1(1, 2)**2+DFGRD1(1, 3)**2
      B(2)=DFGRD1(2, 1)**2+DFGRD1(2, 2)**2+DFGRD1(2, 3)**2
      B(3)=DFGRD1(3, 3)**2+DFGRD1(3, 1)**2+DFGRD1(3, 2)**2
      B(4)=DFGRD1(1, 1)*DFGRD1(2, 1)+DFGRD1(1, 2)*DFGRD1(2, 2)
     1 +DFGRD1(1, 3)*DFGRD1(2, 3)
      IF(NSHR.EQ.3) THEN
          B(5)=DFGRD1(1, 1)*DFGRD1(3, 1)+DFGRD1(1, 2)*DFGRD1(3, 2)
     1 +DFGRD1(1, 3)*DFGRD1(3, 3)
          B(6)=DFGRD1(2, 1)*DFGRD1(3, 1)+DFGRD1(2, 2)*DFGRD1(3, 2)
     1 +DFGRD1(2, 3)*DFGRD1(3, 3)
      END IF
C
C CALCULATE THE STRESS
C
      tp1 = (lam*LOG(det)-mu)/det
      tp2 = mu/det
      DO i = 1,NDI
          STRESS(i) = tp1 + tp2*B(i)
      END DO
      DO i = NDI+1, NDI+NSHR
          STRESS(i) = tp2*B(i)
      END DO

      tp3 = lam/det
      DDSDDE(1, 1) = tp3 + tp2*two*B(1)
      DDSDDE(2, 2) = tp3 + tp2*two*B(2)
      DDSDDE(3, 3) = tp3 + tp2*two*B(3)
      DDSDDE(1, 2) = tp3
      DDSDDE(1, 3) = tp3
      DDSDDE(2, 3) = tp3
      DDSDDE(1, 4) = tp2*B(4)
      DDSDDE(2, 4) = tp2*B(4)
      DDSDDE(3, 4) = 0.D0
      DDSDDE(4, 4) = tp2/two*(B(1)+B(2))
      IF(NSHR.EQ.3) THEN
          DDSDDE(1, 5) = tp2*B(5)
          DDSDDE(2, 5) = 0.D0
          DDSDDE(3, 5) = tp2*B(5)
          DDSDDE(1, 6) = 0.D0
          DDSDDE(2, 6) = tp2*B(6)
          DDSDDE(3, 6) = tp2*B(6)
          DDSDDE(5, 5) = tp2/two*(B(1)+B(3))
          DDSDDE(6, 6) = tp2/two*(B(2)+B(3))
          DDSDDE(4, 5) = tp2/two*B(6)
          DDSDDE(4, 6) = tp2/two*B(5)
          DDSDDE(5, 6) = tp2/two*B(4)
      END IF
      DO i =1, NTENS
          DO j = 1, i-1
              DDSDDE(i, j) = DDSDDE(j, i)
          END DO
      END DO

      STATEV(1) = DTIME

      RETURN
      END

      SUBROUTINE SDVINI(STATEV,COORDS,NSTATV,NCRDS,NOEL,NPT,
     1 LAYER,KSPT)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION STATEV(NSTATV),COORDS(NCRDS)

      STATEV(1) = 1.0D0 ! theta

      RETURN
      END
