**************************************************************************
      SUBROUTINE DK4(IY,I1,I2,I3,I4,WPS)                                    
      IMPLICIT REAL*8(A-H,M,O-Z)
      PARAMETER (PI = 3.141592653d0)                                            
      DIMENSION YY(8),IFLAG(7)
      COMMON/PARTCL/P(5,30),M(2,30)                                          
      COMMON/Dk4CMN/YY,IFLAG
C **********************************************************************
C                    DOUBLE PRECISION VERSION                            
C **********************************************************************
C                                                     
C     FOUR BODY PHASE SPACE DECAY :  IY --> I1 I2 I3 I4 
C     WPS: RETURNED PHASE SPACE WEIGHT
C     Multiply the phase space by the normalization factor: 
C                   1/(2*PI)**8
C     YY array of 8 random numbers supplied in main routine 
C     INDEX: FIRST OF FIVE CONSECUTIVE INDICES FOR HTUPLE                     
C     /PARTCL/:                                                               
C     P(I,*),I=1,3    THREE-MOMENTUM COMPONENTS                               
C     P(4,*)          ENERGY                                                  
C     P(5,*)          MOMENTUM MAGNITUDE                                      
C     M(1,*),M(2,*)   MASS MASS**2
C     ARRAYS M and P(IY,*) HAVE TO BE INITIALISED IN THE CALLING 
C     ROUTINE                                            
C   notation same as Barger and Phillips ................
C ******************************************************************
                                     
      WPS=0.d0                                                              

      DO I = 1,7
      IFLAG(I)=0
      ENDDO

      MY=M(1,IY)                                                              
      MYSQ=M(2,IY)                                                            
      M1SQ=M(2,I1)                                                            
      M2SQ=M(2,I2)                                                            
      M3SQ=M(2,I3)
      M4SQ= M(2,I4)                                                            
      MX1MIN=(M(1,I1)+M(1,I2))**2 
      MX2MIN=(M(1,I3)+M(1,I4))**2                                              
      MX1MAX=(MY-DSQRT(MX2MIN))*(MY-DSQRT(MX2MIN))                                               
      IF(MX1MAX.LT.MX1MIN) THEN
                           IFLAG(1)=1   
                           RETURN 
      ELSE
                        MX1SQ = MX1MIN + (MX1MAX-MX1MIN)*YY(1)
      ENDIF
      WPS = (MX1MAX - MX1MIN) 
      MX1 = DSQRT(MX1SQ)      
      MX2MAX = (MY - MX1)**2
      MX2SQ = MX2MIN + (MX2MAX-MX2MIN)*YY(2)
      MX2 = DSQRT(MX2SQ)
      WPS = WPS*(MX2MAX-MX2MIN)
      XLA1 = ((MYSQ-MX1SQ-MX2SQ)**2-4.d0*MX1SQ*MX2SQ) 
      IF(XLA1.LT.0.d0) THEN       
                        IFLAG(2)=1
                        RETURN                                     
      else                        
                       XLA1=DSQRT(XLA1)
      ENDIF
      WPS = WPS*XLA1/8.0d0/MYSQ
      XLA2=(MX1SQ-M1SQ-M2SQ)**2-4.d0*M1SQ*M2SQ
      IF(XLA2.LT.0.d0) THEN
                       IFLAG(3)=1 
                       RETURN     
      ELSE                                                              
                      XLA2=DSQRT(XLA2)     
      ENDIF                                                    
      WPS=WPS*XLA2/8.0d0/MX1SQ
      XLA3 = (MX2SQ-M3SQ-M4SQ)**2-4.*M3SQ*M4SQ
      IF(XLA3.LT.0.0d0) THEN
                         IFLAG(4) = 1
                         RETURN
      ELSE 
                         XLA3 = DSQRT(XLA3)
      ENDIF
      WPS = WPS*XLA3/8.0d0/MX2SQ
            
C                                                                              
C                                                                              
C     CALCULATE EX1,EX2 ETC IN THE REST FRAME OF Y.
C    THEN CALCULATE THE DECAY OF X1--> 1+2, X2 --> 3+4
C     USE EX1,....TO CALUCLATE    
C     ROTATION AND BOOST TO GET THE MMTA.IN THE REST FRAME OF Y.          
C                                                   
      PX1CM=XLA1/(2.*MY)                                        
      EX1CM=(MYSQ+MX1SQ-MX2SQ)*0.5d0/MY                                       
      EX2CM=(MYSQ+MX2SQ-MX1SQ)*0.5d0/MY
      P1CM = XLA2/(2.d0*MX1)
      P3CM = XLA3/(2.d0*MX2)
      P(1,I1)=0.d0                                                               
      P(2,I1)=0.d0                                                               
      P(3,I1)=P1CM
      P(4,I1)=(MX1SQ+M1SQ-M2SQ)/(2.d0*MX1)
      P(1,I2)=0.d0 
      P(2,I2)=0.d0 
      P(3,I2)=-P1CM
      P(4,I2)=(MX1SQ+M2SQ-M1SQ)/(2.d0*MX1)
      P(1,I3)=0.0d0
      P(2,I3)=0.d0
      P(3,I3)=P3CM 
      P(4,I3)=(MX2SQ+M3SQ-M4SQ)/(2.d0*MX2)
      P(1,I4)=0.0d0
      P(2,I4)=0.d0
      P(3,I4)=-P3CM 
      P(4,I4)=(MX2SQ+M4SQ-M3SQ)/(2.d0*MX2)  
c     write(*,*)'I1',(p(i,i1),i=1,4)

C     *************************************************
C     HERE WE PICK UP POINTS IN THETA AND PHI FOR DECAYS OF
C     X1.
C    ***************************************************
                                     
      CV1 = -1.d0+ 2.0d0*YY(3)  
      SV1=DSQRT(DABS(1.0d0-CV1*CV1))  
      WPS = WPS*2.0d0                                               
      CALL ROT12(I2,1,3,SV1,CV1)                                               
      CALL ROT12(I1,1,3,SV1,CV1)                                               
      PH1 = 2.0d0*PI*YY(4)
      WPS = WPS*2.0d0*PI                                                       
      SPH1 =SIN(PH1)                                                          
      CPH1=COS(PH1)                                                           
      CALL ROT12(I2,2,1,SPH1,CPH1)                                             
      CALL ROT12(I1,2,1,SPH1,CPH1)               

C************************************************************************
C          NOW HERE DO THE SAME FOR THE DECAY X2 -- 3 +4
C***********************************************************************
              
      CV3 = -1.d0+ 2.0d0*YY(5)  
      SV3=DSQRT(ABS(1.0d0-CV3*CV3))  
      WPS = WPS*2.0d0                                               
      CALL ROT12(I3,1,3,SV3,CV3)                                               
      CALL ROT12(I4,1,3,SV3,CV3)                                               
      PH2 = 2.0d0*PI*YY(6)                                                         
      SPH2 =SIN(PH2)                                                          
      CPH2=COS(PH2)                                                           
      CALL ROT12(I3,2,1,SPH2,CPH2)                                             
      CALL ROT12(I4,2,1,SPH2,CPH2)               
      WPS = WPS*2.0d0*PI                                                     

C*********************************************************************
C    CALCULATE NOW THE BOOST FACTORS TO GET THE BOOSTED FOUR MOMENTA
C    FIRST FOR X1 DECAY PRODUCTS AND THEN FOR X2 DECAY PRODUCTS.
C
C********************************************************************
                                                              
      GA1=EX1CM/MX1                                                            
      ETA1=PX1CM/MX1   
      CALL BOOSTZ(I1,GA1,ETA1)   
      CALL BOOSTZ(I2,GA1,ETA1) 

C**********************************************
C     FOR X2 HERE
C   *******************************************
      GA2=EX2CM/MX2                                                            
      ETA2=-PX1CM/MX2                                                          
      CALL BOOSTZ(I3,GA2,ETA2)   
      CALL BOOSTZ(I4,GA2,ETA2)    
                                              
C                                                                              
C     COMPLETE DECAY CONFIGURATION IN Y REST FRAME                             
C                                   
C                                                                              
                                          
      CVY = -1.0d0 + 2.0d0*YY(7)                                                   
      SVY=DSQRT(DABS(1.0d0-CVY*CVY))
      CALL ROT12(I1,1,3,SVY,CVY)
      CALL ROT12(I2,1,3,SVY,CVY)
      CALL ROT12(I3,1,3,SVY,CVY)
      CALL ROT12(I4,1,3,SVY,CVY)
      PHY = 2.0d0*PI*YY(8)
      WPS = WPS*4.0d0*PI
      SPHY=SIN(PHY)
      CPHY=COS(PHY)
      CALL ROT12(I1,2,1,SPHY,CPHY)
      CALL ROT12(I2,2,1,SPHY,CPHY)
      CALL ROT12(I3,2,1,SPHY,CPHY)
      CALL ROT12(I4,2,1,SPHY,CPHY)
C                                               
C     BOOST/ROTATE TO LAB FRAME                                                
C                                                                              
      ETAY=P(5,IY)/MY                                                          
      IF(ETAY.LE.1.D-4)GO TO 200                                               
      GAMMAY=P(4,IY)/MY                                                        
      CALL BOOSTZ(I1,GAMMAY,ETAY)                                              
      CALL BOOSTZ(I2,GAMMAY,ETAY)                                              
      CALL BOOSTZ(I3,GAMMAY,ETAY)  
      CALL BOOSTZ(I4,GAMMAY,ETAY)                                            
      CVLAB=P(3,IY)/P(5,IY)                                                   
      SVLAB=1.0d0-CVLAB*CVLAB
      IF(SVLAB.LT.-.001D0) THEN                               
                             WPS=0.0d0                      
                             IFLAG(5)=1  
                             RETURN 
      ELSE
                        SVLAB=DSQRT(ABS(SVLAB))                           
      ENDIF                                                              
      CALL ROT12(I1,1,3,SVLAB,CVLAB)  
      CALL ROT12(I2,1,3,SVLAB,CVLAB)                                          
      CALL ROT12(I3,1,3,SVLAB,CVLAB)
      CALL ROT12(I4,1,3,SVLAB,CVLAB)               
C*******************************************************
C       CALCULATE PHI FOR Y
C*******************************************************                      
      PTY=DSQRT(P(1,IY)**2+P(2,IY)**2)                                         
      IF(PTY.GE.1.D-4) THEN                                                
                        CPHLB=P(1,IY)/PTY                             
                        SPHLB=P(2,IY)/PTY            
                        CALL ROT12(I1,2,1,SPHLB,CPHLB)        
                        CALL ROT12(I2,2,1,SPHLB,CPHLB)        
                        CALL ROT12(I3,2,1,SPHLB,CPHLB)
                        CALL ROT12(I4,2,1,SPHLB,CPHLB)
      ENDIF     
200   CALL PSET(I1)                                                            
      CALL PSET(I2)                                                            
      CALL PSET(I3)        
      CALL PSET(I4)
c     pe = p(4,i1)+p(4,i2)
c     px = p(1,i1)+p(1,i2)
c     py = p(2,i1)+p(2,i2)
c     pz = p(3,i1)+p(3,i2)
c     pbb2 = pe**2 - (px+py+pz)**2
c     if (pbb2 .lt. 0.d0)then 
c                        WPS=0.0d0                      
c                        IFLAG(6)=1  
c                        RETURN 
c          else
c        pbb = dsqrt(pbb2)
c     endif 
c     mh = 120.d0
c     mhcut = 40.d0
c     if (dabs (pbb-mH) .gt. mHcut)then
c                        WPS=0.0d0                      
c                        IFLAG(7)=1  
c                        RETURN 
c     endif                             
             

c     write(*,*)'III1',(p(i,i1),i=1,4)
c     write(*,*)'I2',(p(i,i2),i=1,4)
c     write(*,*)'I3',(p(i,i3),i=1,4)
c     write(*,*)'I4',(p(i,i4),i=1,4)
c     write(*,*)'======================'
      RETURN                                                                   
      END                                                                      
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
