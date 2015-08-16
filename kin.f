*...last modified on 25.01.2009, ETA is redefined.
*-------------------------- KIN.F ----------------------------------*
       SUBROUTINE PINIT
       IMPLICIT REAL *8(A-H,M,O-Z)
       COMMON /PARTCL/P(5,30),M(2,30)
       DO L=1,30
       DO K=1,5
       P(K,L)=0.d0
       ENDDO   
       M(1,L)=0.d0
       M(2,L)=0.d0
       ENDDO   
       RETURN
       END
*-------------------------------------------*
       SUBROUTINE P_L(I,PL)
       IMPLICIT REAL *8(A-H,M,O-Z)
       COMMON /PARTCL/P(5,30),M(2,30)

       PL=P(3,I)
       RETURN
       END
*-------------------------------------------*
       SUBROUTINE E_L(I,EL)
       IMPLICIT REAL *8(A-H,M,O-Z)
       COMMON /PARTCL/P(5,30),M(2,30)
       EL=DSQRT(M(2,I)+P(3,I)**2)
       RETURN
       END
*-------------------------------------------*
       SUBROUTINE DOT_4(I,J,DOT4)
       IMPLICIT REAL *8(A-H,M,O-Z)
       COMMON /PARTCL/P(5,30),M(2,30)
       CALL DOT_3(I,J,D3)
       DOT4=P(4,I)*P(4,J)-D3
       RETURN
       END
*-------------------------------------------*
       SUBROUTINE SPEED(I,SI,IFL)
       IMPLICIT REAL *8(A-H,M,O-Z)
       COMMON /PARTCL/P(5,30),M(2,30)
       IF(P(4,I).EQ.0.d0)THEN
       SI=0.d0
       IFL=1
       ELSE
       SI=P(5,I)/P(4,I)  
       IFL=0
       ENDIF
       RETURN
       END
*-------------------------------------------*
       SUBROUTINE ANGLE(I,J,A,IFL)
       IMPLICIT REAL *8(A-H,M,O-Z)
       COMMON /PARTCL/P(5,30),M(2,30)
       IF(P(5,I)*P(5,J).EQ.0.d0)THEN
       A=0.d0
       IFL=1
       RETURN
       ENDIF
       CALL DOT_3(I,J,D3)
       COSA=D3/P(5,I)/P(5,J)  
       IF(DABS(COSA).GT. 1.d0)THEN
       A=0.d0
       IFL=1
       ELSE
       A=DACOS(COSA)
       IFL=0
       ENDIF
       RETURN
       END
*-------------------------------------------*
       SUBROUTINE M_INV(I,J,MINV,IFL)
       IMPLICIT REAL *8(A-H,M,O-Z)
       COMMON /PARTCL/P(5,30),M(2,30)
       PX=P(1,I)+P(1,J)
       PY=P(2,I)+P(2,J)
       PZ=P(3,I)+P(3,J)
       PE=P(4,I)+P(4,J)
       M2=PE*PE-PX*PX-PY*PY-PZ*PZ
       IF (M2.LT.0.d0)THEN
       MINV=0.0d0
       IFL=1
       ELSE
       MINV=DSQRT(M2)
       IFL=0
       ENDIF
       RETURN
       END
*-------------------------------------------*
       SUBROUTINE M_INV4(I,J,I1,J1,MINV4,IFL)
       IMPLICIT REAL *8(A-H,M,O-Z)
       COMMON /PARTCL/P(5,30),M(2,30)
       PX=P(1,I)+P(1,J)+P(1,I1)+P(1,J1)
       PY=P(2,I)+P(2,J)+P(2,I1)+P(2,J1)
       PZ=P(3,I)+P(3,J)+P(3,I1)+P(3,J1)
       PE=P(4,I)+P(4,J)+P(4,I1)+P(4,J1)
       M2=PE*PE-PX*PX-PY*PY-PZ*PZ
       IF (M2.LT.0.d0)THEN
       MINV4=0.0d0
       IFL=1
       ELSE
       MINV4=DSQRT(M2)
       IFL=0
       ENDIF
       RETURN
       END
*-------------------------------------------*
       SUBROUTINE M_T(I,J,MT,IFL)
       IMPLICIT REAL *8(A-H,M,O-Z)
       COMMON /PARTCL/P(5,30),M(2,30)
       CALL P_T(I,PTI)
       CALL P_T(J,PTJ)
       CALL E_T(I,ETI)
       CALL E_T(J,ETJ)
       CALL DOT_2(I,J,D2)
       M2=(ETI+ETJ)**2-PTI**2-PTJ**2-2.*D2
       IF (M2.LT.0.d0)THEN
       MT=0.0d0
       IFL=1
       ELSE
       MT=DSQRT(M2)
       IFL=0
       ENDIF
       RETURN
       END
*-------------------------------------------*
       SUBROUTINE RAPDT (I,YI,IFL)
       IMPLICIT REAL *8(A-H,M,O-Z)
       COMMON /PARTCL/P(5,30),M(2,30)
       ARG=(P(4,I)+P(3,I))/(P(4,I)-P(3,I))
       IF(ARG.LT.0.d0)THEN
       YI=1.d6
       IFL=1
       ELSE
*      YI=DABS(DLOG(ARG)) /2.d0
       YI=DLOG(ARG) /2.d0
       IFL=0
       ENDIF
       RETURN
       END
*-------------------------------------------*
       SUBROUTINE DELT_R(I,J,DELTR,IFL)
       IMPLICIT REAL *8(A-H,M,O-Z)
       COMMON /PARTCL/P(5,30),M(2,30)
       PI = DACOS(-1.0d0)
       CALL AZIMUTH(I,PHI) 
       CALL AZIMUTH(J,PHJ) 
       DPHI=(PHI-PHJ)
       IF (DPHI .GT. PI)DPHI = 2.d0*PI - DPHI
       IF (DPHI .LT. -PI)DPHI = 2.d0*PI + DPHI
       DPHI2 = DPHI**2
                           IFLI = 0
                           IFLJ = 0
       CALL PSRAPDT(I,ETAI,IFLI)
       CALL PSRAPDT(J,ETAJ,IFLJ)
       IF (IFLI.EQ.1.OR.IFLJ.EQ.1) THEN
       DELTR = 1.D6
       IFL = 1
       RETURN
       ELSE 
       DETA2=(ETAI-ETAJ)**2
       DELTR=DSQRT(DPHI2+DETA2)
       IFL = 0
       ENDIF
       RETURN
       END
*-------------------------------------------*
       SUBROUTINE JETMERGE(I,J)
       IMPLICIT REAL *8(A-H,M,O-Z)
       COMMON /PARTCL/P(5,30),M(2,30)
       IF(P(5,I).EQ.0.d0.OR.P(5,J).EQ.0.d0)RETURN
       DO K=1,4
       P(K,I)=P(K,I)+P(K,J)
       P(K,J)=0.0d0                   
       ENDDO
       CALL P_MAG(I,P(5,I))
       P(5,J)=0.0d0
       RETURN
       END
*-------------------------------------------*
       SUBROUTINE BOOST(N,I,V,IFL)
       IMPLICIT REAL *8(A-H,M,O-Z)
       COMMON /PARTCL/P(5,30),M(2,30)
       IF(IFL.NE.0) RETURN
       IF(V.GE.1.d0) THEN
       IFL=1
       RETURN
       ENDIF
       GAM=1.d0/DSQRT(1.d0-V**2)
       PP1=P(N,I)
       PP2=P(4,I)
       P(N,I)=GAM*(PP1+V*PP2)
       P(4,I)=GAM*(PP2+V*PP1)
       CALL P_MAG(I,P(5,I))  
       RETURN
       END
*-------------------------------------------*
       SUBROUTINE ROTATE(N,I,ANGL)
       IMPLICIT REAL *8(A-H,M,O-Z)
       COMMON /PARTCL/P(5,30),M(2,30)
       J=N+1
       K=N+2
       IF(J.GT.3)J=J-3
       IF(K.GT.3)K=K-3
       CTH=DCOS(ANGL)         
       STH=DSIN(ANGL)         
       PP1=P(J,I)
       PP2=P(K,I)
       P(J,I)=CTH*PP1+STH*PP2 
       P(K,I)=-STH*PP1+CTH*PP2 
       RETURN
       END
*-------------------------------------------*
       SUBROUTINE P_MAG(I,PMAG)
       IMPLICIT REAL*8(A-H,M,O-Z)
       COMMON /PARTCL/P(5,30),M(2,30)
       PMAG=DSQRT(P(1,I)**2+P(2,I)**2+P(3,I)**2)
       RETURN
       END
*--------------------------------------------*
       SUBROUTINE P_E(I,PE)
       IMPLICIT REAL*8(A-H,M,O-Z)
       COMMON /PARTCL/P(5,30),M(2,30)
       CALL P_MAG(I,PI) 
       PE=DSQRT(M(2,I)+PI**2)
       RETURN
       END
*-------------------------------------------*
       SUBROUTINE P_T(I,PT)
       IMPLICIT REAL*8(A-H,M,O-Z)
       COMMON /PARTCL/P(5,30),M(2,30)
       PT=DSQRT(P(1,I)**2+P(2,I)**2)
       RETURN
       END
*-------------------------------------------*
       SUBROUTINE E_T(I,ET)
       IMPLICIT REAL*8(A-H,M,O-Z)
       COMMON /PARTCL/P(5,30),M(2,30)
       CALL P_T(I,PTI) 
       ET=DSQRT(M(2,I)+PTI**2)
       RETURN
       END
*-------------------------------------------*
       SUBROUTINE DOT_3(I,J,DOT3)
       IMPLICIT REAL*8(A-H,M,O-Z)
       COMMON /PARTCL/P(5,30),M(2,30)
       DOT3=P(1,I)*P(1,J)+P(2,I)*P(2,J)+P(3,I)*P(3,J)
       RETURN
       END
*-------------------------------------------*
       SUBROUTINE DOT_2(I,J,DOT2)
       IMPLICIT REAL*8(A-H,M,O-Z)
       COMMON /PARTCL/P(5,30),M(2,30)
       DOT2=P(1,I)*P(1,J)+P(2,I)*P(2,J)
       RETURN
       END
*-------------------------------------------*
       SUBROUTINE COLAT(I,THETA,IFL)
       IMPLICIT REAL*8(A-H,M,O-Z)
       COMMON /PARTCL/P(5,30),M(2,30)
	PI = DACOS(-1.D0)
       IF(P(5,I).EQ.0.d0) THEN
       THETA=0.d0
       IFL=1
       ELSE
       IF(P(3,I).EQ.0.d0)THEN
       THETA=DASIN(1.d0)
       IFL=0
       ELSE
       CSTH=P(3,I)/P(5,I)
*--------Theta in Rad.:
       THETA= DACOS(CSTH)
       ENDIF
       ENDIF
       RETURN
       END
*-------------------------------------------*
       SUBROUTINE AZIMUTH(I,PHI)
       IMPLICIT REAL*8(A-H,M,O-Z)
       COMMON /PARTCL/P(5,30),M(2,30)
       IF(P(2,I).EQ.0.d0) THEN
       PHI=0.d0
       ELSE
       IF(P(1,I).EQ.0.d0)THEN
       PHI=DASIN(1.d0)
       ELSE
       TANPHI=P(2,I)/P(1,I)
       PHI=DATAN (TANPHI)
       IF (P(1,I) .LT. 0.d0) THEN
                           PHI = PHI + DACOS(-1.d0)
       ELSE
           IF (P(2,I) .LT. 0.d0) PHI = PHI + 2.d0 * DACOS(-1.d0)
       ENDIF
       ENDIF
       ENDIF
       RETURN
       END
*-------------------------------------------*
       SUBROUTINE PSRAPDT(I,ETA,IFL)
       IMPLICIT REAL*8(A-H,M,O-Z)
       COMMON /PARTCL/P(5,30),M(2,30)
c
       PTI = dsqrt(p(1,i)**2 + p(2,i)**2)
       arg = (p(4,i) + p(3,i))/pti
       if (arg .lt. 1D-13 .or. arg .gt. 1D13)then
         eta = 100D0
         ifl = 1
        else
       eta = Dlog(arg)
        ifl = 0
       endif
c      eta = dsqrt(p(1,i)**2 + p(2,i)**2 + p(3,i)**2 + m(2,i))
c      eta = (eta + p(3,i))/(eta - p(3,i))
c      if (eta .lt. 1D-13) then
c      eta = 100D0
c      ifl = 1
c      else
c      eta = 0.5D0 * DLOG(eta)
c      endif
       RETURN
       END
*-------------------------------------------*
       FUNCTION KALLEN(X,Y,Z)
       IMPLICIT REAL*8(A-H,M,O-Z)
       KALLEN=(X+Y-Z)**2-4.d0*X*Y
       RETURN
       END
*-------------------------------------------*
       FUNCTION  Y_JADE(I,J,S,IFL)
       IMPLICIT REAL*8(A-H,M,O-Z)
       COMMON /PARTCL/P(5,30),M(2,30)
       CALL ANGLE(I,J,A,IFL)
       IF (IFL.EQ.1) THEN
       Y_JADE=0.d0
       RETURN
       ENDIF
       Y_JADE=2.d0*P(4,I)*P(4,J)*(1.d0-COS(A))/ S
       RETURN
       END
*-------------------------------------------*
       FUNCTION  Y_DURM(I,J,S,IFL)
       IMPLICIT REAL*8(A-H,M,O-Z)
       COMMON /PARTCL/P(5,30),M(2,30)
       CALL ANGLE(I,J,A,IFL)
       IF (IFL.EQ.1) THEN
       Y_DURM=0.d0
       RETURN
       ENDIF
       PMAX = P(4,I)
       IF(P(4,J).GT.P(4,I))PMAX=P(4,J)
       Y_DURM=2.d0*PMAX**2 *(1.d0-COS(A))/S
       RETURN
       END
*--------------------------------------------*
       SUBROUTINE EVENT_LOSS
       IMPLICIT REAL*8(A-H,M,O-Z)
       COMMON /SHOTS/ NPT, ITN
       COMMON /CUTS/ NCUT(30), NERR(30)
       NERR_TOT = 0
       DO I = 1,30
       NERR_TOT = NERR_TOT + NERR(I)
       ENDDO
       IF(NERR_TOT .NE. 0) THEN
       PCTOT = FLOAT(NERR_TOT)/NPT*100
       WRITE(4,*)
     > 'There were freak configurations in ',PCTOT,
     > ' % of the generated events:'
       DO I = 1,30
       IF (NERR(I) .NE. 0) THEN
       PCERR = FLOAT(NERR(I))/NPT*100
       WRITE(4,*)'Error:', I,' affected',PCERR,' %'
       ENDIF
       ENDDO
       ELSE
       WRITE(4,*) 'No freak configurations were detected'
       ENDIF
       NCUT_TOT = 0
       DO I = 1,30
       NCUT_TOT = NCUT_TOT + NCUT(I)
       ENDDO
       IF(NCUT_TOT .NE. 0) THEN
       PCTOT = FLOAT(NCUT_TOT)/NPT*100
       WRITE(4,*) 
     > 'Kinematic cuts removed     ',PCTOT,
     > ' % of the generated events:'
       DO I = 1,30
       IF (NCUT(I) .NE. 0) THEN
       PCCUT = FLOAT(NCUT(I))/NPT*100
       WRITE(4,*)'Cut:  ', I,' removed',PCCUT,' %'
       ENDIF
       ENDDO
       ELSE
       WRITE(4,*) 'No events were lost due to kinematic cuts'
       ENDIF
       WRITE(4,200)
 200   FORMAT(  1X, 75(1H=))
       RETURN
       END
*-------------------------------------------------------------------*
       subroutine r2(i,j,deltar,ifl)
       IMPLICIT REAL*8(A-H,M,O-Z)
       COMMON /PARTCL/P(5,30),M(2,30)
c----calculate the jets separation between p(i) and p(j)
               ifl = 0
      ei=dsqrt(p(1,i)**2+p(2,i)**2+p(3,i)**2)
      ej=dsqrt(p(1,j)**2+p(2,j)**2+p(3,j)**2)

      r1= (ei+p(3,i))*(ej-p(3,j))/
     .     ((ej+p(3,j))*(ei-p(3,i)))

        if ( r1 .lt. 0.d0) then
                           deltar = 1.d06
                           ifl = 1
                           return
        endif
      dely=0.5d0*dlog(r1)
      rr= (p(1,i)*p(1,j)+p(2,i)*p(2,j))
     .     /dsqrt((p(1,i)**2+p(2,i)**2)*(p(1,j)**2+p(2,j)**2))
      if (rr .gt. +0.9999999D0) rr=+1D0
      if (rr .lt. -0.9999999D0) rr=-1D0
      delphi=dacos(rr)

       deltar =dsqrt(dely**2+delphi**2)
      
      return
      end
***********************************************************************
           subroutine sphercity(ipar,npar,sT)
c
c          Transverse one
c
           IMPLICIT REAL*8(A-H,M,O-Z)
           real*8 lambda_1, lambda_2
           COMMON /PARTCL/P(5,30),M(2,30)
           dimension ipar(10)
            AA = 0.d0
            BB = 0.d0
            CC = 0.d0
*...Construction of Sphericity tensor S:
              do ii = 1,npar
               AA =  AA + p(1,ipar(ii))* p(1,ipar(ii))
               BB =  BB + p(1,ipar(ii))* p(2,ipar(ii))
               CC =  CC + p(2,ipar(ii))* p(2,ipar(ii))
              enddo
*...Eigenvalues :
            DISC  = dsqrt( (AA+CC)**2 - 4.d0 * (AA*CC - BB*BB) )
            lambda_1 = 0.5d0*((AA+CC) + disc )
            lambda_2 = 0.5d0*((AA+CC) - disc )

*...Transverse Sphericity:

           sT = 2.d0 * lambda_2 /( lambda_1 + lambda_2)
 
           return
           end
********************************************************************
         subroutine sort(n,x)
          IMPLICIT REAL*8(A-H,M,O-Z)
           dimension x(n)

              do ii =1,n-1
                 jj = n - ii
                   do 1110 ll  =1,jj
                    if(x(ll).le.x(ll+1))goto 1110 
                    temp = x(ll)
                    x(ll) = x(ll+1)
                    x(ll+1) = temp
1110                 continue
            enddo                
                    
            return  
            end
*********************************************************************
