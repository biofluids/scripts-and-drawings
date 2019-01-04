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

      parameter(ZERO=0.D0, ONE=1.D0, TWO=2.D0, TOL=1.D-4)
      double precision det, mu, lam, tp1, tp2, tp3, tp4,
     1 theta, theta_pre
      integer i, j, cnt
      logical flag
      double precision Fe(3,3), Be(6), Ce(6)
      double precision Je, trce, trme, dr1, dr2, kg, local_res, local_K
C
C Note theta is the growth factor at the current time step (increment)
C Don't want to update it unless this is the last iteration of the current time step (increment)
C According to abaqus, statev is updated only at the end of a time step
C
      theta_pre = statev(1)
      theta = statev(1)
      cnt = 0
      flag = .false.
      mu = props(1)
      lam = props(2)
      mecr = props(3) ! growth criterion
      tmax = props(4) ! theta_max
      gmax = props(5) ! positive growth exponent
      amax = props(6) ! adaption speed (alpha = 1/tau) for positive growth
C      tmin = props(7) ! theta_min
C      gmin = props(8) ! negative growth exponent
C      amin = props(9) ! adaption speed for negative growth
      local_res = 1.D0

      call get_elastic(DFGRD1,theta,mu,lam,Fe,Be,Ce,Je,trce,trme)
C Local Newton iteration to calculate theta, currently only positive growth is allowed
      do while (kstep .eq. 3 .and. local_res > tol
     1  .and. dabs(tmax - theta) > tol .and. trme - mecr > tol)
          kg = amax*((tmax - theta_pre)/(tmax - 1.D0))**gmax
          local_res = theta - theta_pre - kg*(trme - mecr)*dtime
          dr1 = -(9.D0*lam + 2.D0*mu*trce)/theta ! partial derivative of trme w.r.t. theta
          dr2 = -gmax*kg/(tmax - theta) ! partial derivative of kg w.r.t. theta
          local_K = 1.D0 - (kg*dr1 + (trme - mecr)*dr2)*dtime
          theta = theta - local_res/local_K
          flag = .true.
          call get_elastic(DFGRD1,theta,mu,lam,Fe,Be,Ce,Je,trce,trme)
          cnt = cnt + 1
          if (theta - tmax > tol) then
              write(*,*) 'Maximum value of theta exceeded!'
          else if (cnt .eq. 500) then
              write(*,*) 'Maximum number of theta iteration exceeded!'
          end if
      end do

C
C Total Jacobian
C
      det = Je*theta**3
C
C CALCULATE THE STRESS
C
      tp1 = (lam*LOG(Je)-mu)/det
      tp2 = mu/det
      DO i = 1,NDI
          STRESS(i) = tp1 + tp2*Be(i)
      END DO
      DO i = NDI+1, NDI+NSHR
          STRESS(i) = tp2*Be(i)
      END DO

C Calculate the tangent with kinematic terms required by Abaqus
C
      tp3 = lam/det
      DDSDDE(1, 1) = tp3 + tp2*two*Be(1)
      DDSDDE(2, 2) = tp3 + tp2*two*Be(2)
      DDSDDE(3, 3) = tp3 + tp2*two*Be(3)
      DDSDDE(1, 2) = tp3
      DDSDDE(1, 3) = tp3
      DDSDDE(2, 3) = tp3
      DDSDDE(1, 4) = tp2*Be(4)
      DDSDDE(2, 4) = tp2*Be(4)
      DDSDDE(3, 4) = 0.D0
      DDSDDE(4, 4) = tp2/two*(Be(1)+Be(2))
      IF(NSHR.EQ.3) THEN
          DDSDDE(1, 5) = tp2*Be(5)
          DDSDDE(2, 5) = 0.D0
          DDSDDE(3, 5) = tp2*Be(5)
          DDSDDE(1, 6) = 0.D0
          DDSDDE(2, 6) = tp2*Be(6)
          DDSDDE(3, 6) = tp2*Be(6)
          DDSDDE(5, 5) = tp2/two*(Be(1)+Be(3))
          DDSDDE(6, 6) = tp2/two*(Be(2)+Be(3))
          DDSDDE(4, 5) = tp2/two*Be(6)
          DDSDDE(4, 6) = tp2/two*Be(5)
          DDSDDE(5, 6) = tp2/two*Be(4)
      END IF
C
C Add growth term to tangent
C
      if (flag) then
          tp4 = -kg*dtime/(det*local_K*theta)
          DDSDDE(1, 1) = DDSDDE(1, 1) + tp4*(3.D0*lam+2.D0*mu*Be(1))**2
          DDSDDE(2, 2) = DDSDDE(2, 2) + tp4*(3.D0*lam+2.D0*mu*Be(2))**2
          DDSDDE(3, 3) = DDSDDE(3, 3) + tp4*(3.D0*lam+2.D0*mu*Be(3))**2
          DDSDDE(1, 2) = DDSDDE(1, 2)
     1 + tp4*(3.D0*lam+2.D0*mu*Be(1))*(3.D0*lam+2.D0*mu*Be(2))
          DDSDDE(1, 3) = DDSDDE(1, 3)
     1 + tp4*(3.D0*lam+2.D0*mu*Be(1))*(3.D0*lam+2.D0*mu*Be(3))
          DDSDDE(2, 3) = DDSDDE(2, 3)
     1 + tp4*(3.D0*lam+2.D0*mu*Be(2))*(3.D0*lam+2.D0*mu*Be(3))
          DDSDDE(1, 4) = DDSDDE(1, 4)
     1 + tp4*(3.D0*lam+2.D0*mu*Be(1))*(2.D0*mu*Be(4))
          DDSDDE(2, 4) = DDSDDE(2, 4)
     1 + tp4*(3.D0*lam+2.D0*mu*Be(2))*(2.D0*mu*Be(4))
          DDSDDE(3, 4) = DDSDDE(3, 4)
     1 + tp4*(3.D0*lam+2.D0*mu*Be(3))*(2.D0*mu*Be(4))
          DDSDDE(4, 4) = DDSDDE(4, 4) + tp4*(2.D0*mu*Be(4))**2
          IF(NSHR.EQ.3) THEN
              DDSDDE(1, 5) = DDSDDE(1, 5)
     1 + tp4*(3.D0*lam+2.D0*mu*Be(1))*(2.D0*mu*Be(5))
              DDSDDE(2, 5) = DDSDDE(2, 5)
     1 + tp4*(3.D0*lam+2.D0*mu*Be(2))*(2.D0*mu*Be(5))
              DDSDDE(3, 5) = DDSDDE(3, 5)
     1 + tp4*(3.D0*lam+2.D0*mu*Be(3))*(2.D0*mu*Be(5))
              DDSDDE(1, 6) = DDSDDE(1, 6)
     1 + tp4*(3.D0*lam+2.D0*mu*Be(1))*(2.D0*mu*Be(6))
              DDSDDE(2, 6) = DDSDDE(2, 6)
     1 + tp4*(3.D0*lam+2.D0*mu*Be(2))*(2.D0*mu*Be(6))
              DDSDDE(3, 6) = DDSDDE(3, 6)
     1 + tp4*(3.D0*lam+2.D0*mu*Be(3))*(2.D0*mu*Be(6))
              DDSDDE(5, 5) = DDSDDE(5, 5) + tp4*(2.D0*mu*Be(5))**2
              DDSDDE(6, 6) = DDSDDE(6, 6) + tp4*(2.D0*mu*Be(6))**2
              DDSDDE(4, 5) = DDSDDE(4, 5) 
     1 + tp4*(2.D0*mu*Be(4))*(2.D0*mu*Be(5))
              DDSDDE(4, 6) = DDSDDE(4, 6)
     1 + tp4*(2.D0*mu*Be(4))*(2.D0*mu*Be(6))
              DDSDDE(5, 6) = DDSDDE(5, 6)
     1 + tp4*(2.D0*mu*Be(5))*(2.D0*mu*Be(6))
          END IF
      end if
      
      DO i =1, NTENS
          DO j = 1, i-1
              DDSDDE(i, j) = DDSDDE(j, i)
          END DO
      END DO
      
      STATEV(1) = theta
      STATEV(2) = trme
      RETURN
      END

      subroutine get_elastic(F,theta,mu,lam,Fe,Be,Ce,Je,trce,trme)
      double precision F(3,3), Fe(3,3), Be(6), Ce(6)
      double precision theta, mu, lam, trce, trme, Je

      Fe = F/theta

      Je = Fe(1, 1)*Fe(2, 2)*Fe(3, 3)
     1 -Fe(1, 2)*Fe(2, 1)*Fe(3, 3)
      Je = Je + Fe(1, 2)*Fe(2, 3)*Fe(3, 1)
     1 +Fe(1, 3)*Fe(3, 2)*Fe(2, 1)
     2 -Fe(1, 3)*Fe(3,1)*Fe(2, 2)
     3 -Fe(2, 3)*Fe(3, 2)*Fe(1, 1)

      Be(1)=Fe(1, 1)**2+Fe(1, 2)**2+Fe(1, 3)**2
      Be(2)=Fe(2, 1)**2+Fe(2, 2)**2+Fe(2, 3)**2
      Be(3)=Fe(3, 3)**2+Fe(3, 1)**2+Fe(3, 2)**2
      Be(4)=Fe(1, 1)*Fe(2, 1)+Fe(1, 2)*Fe(2, 2)
     1 +Fe(1, 3)*Fe(2, 3)
      Be(5)=Fe(1, 1)*Fe(3, 1)+Fe(1, 2)*Fe(3, 2)
     1 +Fe(1, 3)*Fe(3, 3)
      Be(6)=Fe(2, 1)*Fe(3, 1)+Fe(2, 2)*Fe(3, 2)
     1 +Fe(2, 3)*Fe(3, 3)

      Ce(1)=Fe(1, 1)**2+Fe(2, 1)**2+Fe(3, 1)**2
      Ce(2)=Fe(1, 2)**2+Fe(2, 2)**2+Fe(3, 2)**2
      Ce(3)=Fe(1, 3)**2+Fe(2, 3)**2+Fe(3, 3)**2
      Ce(4)=Fe(1, 1)*Fe(1, 2)+Fe(2, 1)*Fe(2, 2)
     1 +Fe(3, 1)*Fe(3, 2)
      Ce(5)=Fe(1, 1)*Fe(1, 3)+Fe(2, 1)*Fe(2, 3)
     1 +Fe(3, 1)*Fe(3, 3)
      Ce(6)=Fe(1, 2)*Fe(1, 3)+Fe(2, 2)*Fe(2, 3)
     1 +Fe(3, 2)*Fe(3, 3)

      trce = Ce(1) + Ce(2) + Ce(3)
      trme = (lam*LOG(Je)-mu)*3.D0 + mu*trce
      return
      end

      SUBROUTINE SDVINI(STATEV,COORDS,NSTATV,NCRDS,NOEL,NPT,
     1 LAYER,KSPT)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION STATEV(NSTATV),COORDS(NCRDS)


      STATEV(1) = 1.0D0 ! theta
      STATEV(2) = 0.0D0 ! trace of Me

      RETURN
      END
      
      SUBROUTINE  DISP(U,KSTEP,KINC,TIME,NODE,NOEL,JDOF,COORDS)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION U(3),TIME(2),COORDS(3)
C
      parameter(pi = 3.1415926d0,phi=0.52359878d0,rin0=1.7467112d-3,
     1 rout0=2.6782905d-3, tol=1.d-6, lam=1.13d0, dr = 0.0003)
      double precision r0, theta0, rout, rin, r, theta
      
C phi = 45/30 degrees (0.7853982d0/0.52359878d0),rin0=1.5/1.81659021246/1.72336879 mm,rout0=2.3/2.78543832577/2.64249882 mm
C phi = 30 degree (0.52359878d0) lambda = 1.2 rout0 = 2.76, rin0 = 1.8
C Assuming dr/r0 = const, dtheta/theta0 = [-1,0], dz/z0 = lam-1
C (pi-phi)/pi*(rout0**2-rin0**2) = lam*(rout**2-rin**2)
C (rout-rout0)/rout0 = (rin - rin0)/rin0
C And don't forget to change Line 44 and other hard-coded parameters
      r0 = sqrt(coords(1)**2 + coords(2)**2)
      theta0 = asin(coords(2)/r0)
      if (coords(1) < 0.D0) then
          theta0 = pi - theta0
      end if
      rout = rout0*sqrt((pi-phi)/pi/lam)
      rin = rin0*sqrt((pi-phi)/pi/lam)
      r = r0*sqrt((pi-phi)/pi/lam) ! new r
      rout = 2.3d-3
      rin = 1.5d-3
      r = r0*rout/rout0
      theta = theta0 - (theta0-pi)/(phi-pi)*phi ! new theta
      if (kstep .eq. 1) then
          if (JDOF .eq. 1) then
              U(1) = r*cos(theta) - coords(1)
          else if (JDOF .eq. 2) then
              U(1) = r*sin(theta) - coords(2)
          else if (JDOF .eq. 3) then
              U(1) = 0.13*coords(3)
          end if
      else if (kstep .eq. 2) then
          if (JDOF .eq. 1) then
              U(1) =dr*cos(theta0)+coords(1)
     1 - cos((5*theta0+pi)/6)*rin0
C     - cos(0.75*theta0+pi/4)*rin0
          else if (JDOF .eq. 2) then
              U(1) =dr*sin(theta0)+coords(2)
     1 - sin((5*theta0+pi)/6)*rin0
C     - sin(0.75*theta0+pi/4)*rin0
          end if
      end if
      RETURN
      END
