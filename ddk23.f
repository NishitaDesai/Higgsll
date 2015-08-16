      SUBROUTINE DK3(IY,I1,I2,I3,WPS)                                         
      IMPLICIT DOUBLE PRECISION(A-H,M,O-Z)                                               
      PARAMETER (PI = 3.1415926536d0)                                             
      COMMON /PARTCL/P(5,30),M(2,30)                                           
      COMMON /DK3CMN/S1,S2,S3,YY(5),IFLAG                                      
C                                                                               
C     THREE BODY PHASE SPACE DECAY :  IY --> I1 I2 I3                           
C     WPS: RETURNED PHASE SPACE WEIGHT                                          
C     INDEX: FIRST OF FIVE CONSECUTIVE INDICES FOR HTUPLE                       
C     /PARTCL/:                                                                 
C     P(I,*),I=1,3    THREE-MOMENTUM COMPONENTS                                 
C     P(4,*)          ENERGY                                                    
C     P(5,*)          MOMENTUM MAGNITUDE                                        
C     M(1,*),M(2,*)   MASS MASS**2                                              
C                                                                               
      WPS=0.0D0                                                                
      IFLAG=0                                                                  
      MY=M(1,IY)                                                               
      MYSQ=M(2,IY)                                                             
      M1SQ=M(2,I1)                                                             
      M2SQ=M(2,I2)                                                             
      M3SQ=M(2,I3)                                                             
      MXMIN=(M(1,I2)+M(1,I3))**2                                               
      MXMAX=(MY-M(1,I1))**2     
      IF(MXMAX.GT.MXMIN)GO TO 100                                              
      IFLAG=1                                                                  
      RETURN                                                                   
100   MXSQ = MXMIN + (MXMAX-MXMIN)*YY(1)                                       
      MX = DSQRT(MXSQ)                                                         
      XLA=((MYSQ-MXSQ-M1SQ)**2-4.d0*M1SQ*MXSQ)                                   
      IF(XLA.GT.0.0D0)GO TO 110                                                
      IFLAG=2                                                                  
      RETURN                                                                   
110   XLA=DSQRT(XLA)                                                           
      XLB=(MXSQ-M2SQ-M3SQ)**2-4.d0*M2SQ*M3SQ                                     
      IF(XLB.GT.0.D0)GO TO 120                                                   
      IFLAG=3                                                                  
      RETURN                                                                   
120   XLB=DSQRT(XLB)                                                           
      WPS=(MXMAX-MXMIN)*(XLB/MXSQ)*(XLA/MYSQ)                                  

C                                                                              
C     REST FRAME DECAY OF X  :  X-->2+3                                        
C     ROTATION AND BOOST                                                       
C                                                                              
      PX=XLB/(2.d0*MX)                                              
      P(1,I2)=0.d0                                                              
      P(2,I2)=0.d0                                                             
      P(3,I2)=PX                                                               
      P(4,I2)=(MXSQ+M2SQ-M3SQ)/(2.d0*MX)                                        
      P(1,I3)=0.d0                                                       
      P(2,I3)=0.d0                                                             
      P(3,I3)=-PX                                                              
      P(4,I3)=(MXSQ+M3SQ-M2SQ)/(2.d0*MX) 
      CV = -1.d0+ 2.d0*YY(2)                                                  
      SV=DSQRT(DABS(1.d0-CV*CV))                                                 
      CALL ROT12(I2,1,3,SV,CV)                                                 
      CALL ROT12(I3,1,3,SV,CV)                                                 
      PH = 2.d0*PI*YY(3)                                                       
      SPH =DSIN(PH)                                                            
      CPH=DCOS(PH)                                                             
      CALL ROT12(I2,2,1,SPH,CPH)                                               
      CALL ROT12(I3,2,1,SPH,CPH)                                               
      EXCM=(MYSQ+MXSQ-M1SQ)/(2.d0*MY)                                          
      PXCM=XLA/(2.d0*MY)                                                       
      GA=EXCM/MX                                                               
      ETA=PXCM/MX                                                              
      CALL BOOSTZ(I2,GA,ETA)                                                   
      CALL BOOSTZ(I3,GA,ETA)                                                   
C                                                                              
C     COMPLETE DECAY CONFIGURATION IN Y REST FRAME                             
C     COMPUTE/STORE INVARIANTS                                                 
C                                                                              
      P(1,I1)=0.d0                                                              
      P(2,I1)=0.d0                                                              
      P(3,I1)=-PXCM                                                           
      P(4,I1)=(MYSQ+M1SQ-MXSQ)/(2.*MY)                                         
      CV = -1.d0 + 2.d0*YY(4)                                                   
      SV=DSQRT(DABS(1.d0-CV*CV))                                                
      CALL ROT12(I1,1,3,SV,CV)                                                
      CALL ROT12(I2,1,3,SV,CV)                                                
      CALL ROT12(I3,1,3,SV,CV)                                                
      PH = 2.d0*PI*YY(5)                                                       
      SPH=DSIN(PH)                                                            
      CPH=DCOS(PH)                                                            
      CALL ROT12(I1,2,1,SPH,CPH)                                              
      CALL ROT12(I2,2,1,SPH,CPH)                                              
      CALL ROT12(I3,2,1,SPH,CPH)                                              
      S1=MYSQ+M1SQ-2.d0*MY*P(4,I1)                                              
      S2=MYSQ+M2SQ-2.d0*MY*P(4,I2)                                              
      S3=MYSQ+M3SQ-2.d0*MY*P(4,I3)                                              
C                                                                            
C     BOOST/ROTATE TO LAB FRAME                                                
C                                                                              
      ETAY=P(5,IY)/MY                                                          
      IF(ETAY.LE.1.D-4)GO TO 200                                              
      GAMMAY=P(4,IY)/MY                                                       
      CALL BOOSTZ(I1,GAMMAY,ETAY)                                             
      CALL BOOSTZ(I2,GAMMAY,ETAY)                                              
      CALL BOOSTZ(I3,GAMMAY,ETAY)                                              
      CV=P(3,IY)/P(5,IY)                                                       
      SV=1.d0-CV*CV                                                            
      IF(SV.GT.-.001D0)GO TO 130                                               
      WPS=0.D0                                                                 
      IFLAG=4                                                                  
      RETURN                                                                   
130   SV=DSQRT(DABS(SV))                                                       
      CALL ROT12(I1,1,3,SV,CV)                                                 
      CALL ROT12(I2,1,3,SV,CV)                                                 
      CALL ROT12(I3,1,3,SV,CV)                                                
      PTY=DSQRT(P(1,IY)**2+P(2,IY)**2)                                         
      IF(PTY.LE.1.D-4)GO TO 200                                                
      CPH=P(1,IY)/PTY                                                          
      SPH=P(2,IY)/PTY                                                          
      CALL ROT12(I1,2,1,SPH,CPH)                                               
      CALL ROT12(I2,2,1,SPH,CPH)                                               
      CALL ROT12(I3,2,1,SPH,CPH)                                               
200   CALL PSET(I1)                                                            
      CALL PSET(I2)                                                            
      CALL PSET(I3)                                                            
      RETURN                                                                   
      END                                                                      
                                                                               
C-------------------------------------------------------------------------

      SUBROUTINE DK2(IY,I1,I2,WPS)
      IMPLICIT REAL*8(A-H,M,O-Z)
      COMMON /PARTCL/P(5,30),M(2,30)
      COMMON/DK2CMN/YY(2),IFLAG
C
C     TWO BODY PHASE SPACE DECAY :  IY --> I1 I2
C     WPS: RETURNED PHASE SPACE WEIGHT
C     /PARTCL/:
C     P(I,*),I=1,3    THREE-MOMENTUM COMPONENTS
C     P(4,*)          ENERGY
C     P(5,*)          MOMENTUM MAGNITUDE
C     M(1,*),M(2,*)   MASS MASS**2
C
       PI = 3.1415926536D0
       WPS = 0.0D0
       IFLAG = 0
       MY = M(1,IY)
       MYSQ = M(2,IY)
       M1SQ = M(2,I1)
       M2SQ = M(2,I2)
       XLA = ((MYSQ-M2SQ-M1SQ) * (MYSQ-M2SQ-M1SQ)-4. * M1SQ * M2SQ)
C       WRITE(*,*) 'IY,I1,I2',IY,I1,I2
       IF(XLA.GT.0)GO TO 110
       IFLAG = 2
       RETURN
110    XLA = DSQRT(DABS(XLA))
       WPS = 1.d0
       EXCM = (MYSQ + M2SQ-M1SQ)/(2.d0 * MY)
       PXCM = XLA/(2.d0 * MY)
C     COMPLETE DECAY CONFIGURATION IN Y REST FRAME
C     COMPUTE/STORE INVARIANTS
C
       P(1,I2) = 0.d0
       P(2,I2) = 0.d0
       P(3,I2) = PXCM
       P(4,I2) = (MYSQ + M2SQ-M1SQ)/(2.d0 * MY)
       P(1,I1) = 0.d0
       P(2,I1) = 0.d0
       P(3,I1) = -PXCM
       P(4,I1) = (MYSQ + M1SQ-M2SQ)/(2.d0 * MY)
       CV = -1.d0 + 2.d0 * YY(1)
       SV = DSQRT(DABS(1.d0-CV * CV))
       CALL ROT12(I1,1,3,SV,CV)
       CALL ROT12(I2,1,3,SV,CV)
       PH = 2.d0 * PI * YY(2)
       SPH = DSIN(PH)
       CPH = DCOS(PH)
       CALL ROT12(I1,2,1,SPH,CPH)
       CALL ROT12(I2,2,1,SPH,CPH)
       S1 = MYSQ + M1SQ-2.d0 * MY * P(4,I1)
       S2 = MYSQ + M2SQ-2.d0 * MY * P(4,I2)
C
C     BOOST/ROTATE TO LAB FRAME
C
       ETAY = P(5,IY)/MY
       IF(ETAY.LE.1.D-4)GO TO 200
       GAMMAY = P(4,IY)/MY
       CALL BOOSTZ(I1,GAMMAY,ETAY)
       CALL BOOSTZ(I2,GAMMAY,ETAY)

       CV = P(3,IY)/P(5,IY)
       SV = 1.d0-CV * CV
       IF(SV.GT.-.001d0)GO TO 130
       WPS = 0.E0
       IFLAG = 4
       RETURN
130    SV = DSQRT(DABS(SV))
       CALL ROT12(I1,1,3,SV,CV)
       CALL ROT12(I2,1,3,SV,CV)
       PTY = DSQRT(DABS(P(1,IY) * P(1,IY) + P(2,IY) * P(2,IY)))
       IF(PTY.LE.1.D-4)GO TO 200
       CPH = P(1,IY)/PTY
       SPH = P(2,IY)/PTY
       CALL ROT12(I1,2,1,SPH,CPH)
       CALL ROT12(I2,2,1,SPH,CPH)
c       WRITE(*,*)'Y', p(1,iy),p(2,iy),p(3,iy),p(4,iy)
c       WRITE(*,*)'1', p(1,i1),p(2,i1),p(3,i1),p(4,i1)
c       WRITE(*,*)'2', p(1,i2),p(2,i2),p(3,i2),p(4,i2)
c       WRITE(*,*)'=========================='
200    CALL PSET(I1)
       CALL PSET(I2)

       RETURN
       END

      SUBROUTINE BOOSTZ(I,GAMMA,ETA)                                           
      IMPLICIT REAL*8(A-H,M,O-Z)                                               
      COMMON /PARTCL/ P(5,30),M(2,30)                                          
C                                                                              
C     LORENTZ BOOST ALONG Z AXIS                                               
C                                                                              
      TMP=GAMMA*P(3,I)+ETA*P(4,I)                                              
      P(4,I)=GAMMA*P(4,I)+ETA*P(3,I)                                           
      P(3,I)=TMP                                                               
      P(5,I)=DSQRT(P(1,I)**2+P(2,I)**2+P(3,I)**2)                              
      RETURN                                                                   
      END                                                                      
C                                                                              
C     LORENTZ BOOST ALONG X AXIS                                               
C                                                                              
      SUBROUTINE BOOSTXX(I,GAMMA,ETA)                                           
      IMPLICIT REAL*8(A-H,M,O-Z)                                               
      COMMON/PARTCL/ P(5,30),M(2,30)                                           
      TMP=GAMMA*P(1,I)+ETA*P(4,I)                                              
      P(4,I)=GAMMA*P(4,I)+ETA*P(1,I)                                           
      P(1,I)=TMP                                                               
      P(5,I)=DSQRT(P(1,I)*P(1,I)+P(2,I)*P(2,I)+P(3,I)*P(3,I))                  
      RETURN                                                                   
      END                                                                      
C                                                                              
C     LORENTZ BOOST ALONG Y AXIS                                               
C                                                                              
      SUBROUTINE BOOSTY(I,GAMMA,ETA)                                           
      IMPLICIT REAL*8(A-H,M,O-Z)                                               
      COMMON/PARTCL/P(5,30),M(2,30)                                            
      TMP=GAMMA*P(2,I)+ETA*P(4,I)                                              
      P(4,I)=GAMMA*P(4,I)+ETA*P(2,I)                                           
      P(2,I)=TMP                                                               
      P(5,I)=DSQRT(P(1,I)*P(1,I)+P(2,I)*P(2,I)+P(3,I)*P(3,I))                  
      RETURN                                                                   
      END                                                                      
C                                                                              
C     THREE SPACE ROTATION                                                     
C                                                                              
      SUBROUTINE ROT12(I,K1,K2,S,C)                                            
      IMPLICIT REAL*8(A-H,M,O-Z)                                               
      COMMON/PARTCL/P(5,30),M(2,30)                                            
      TMP=C*P(K1,I)+S*P(K2,I)                                                  
      P(K2,I)=C*P(K2,I)-S*P(K1,I)                                              
      P(K1,I)=TMP                                                              
      RETURN                                                                   
      END                                                                      
C                                                                              
C     STORE MAGNITUDE OF MOMENTUM VECTOR                                       
C                                                                              
      SUBROUTINE PSET(I)                                                       
      IMPLICIT REAL*8(A-H,M,O-Z)                                               
      COMMON/PARTCL/P(5,30),M(2,30)                                            
      TMP=P(1,I)**2+P(2,I)**2+P(3,I)**2                                        
      P(5,I)=DSQRT(TMP)                                                        
      RETURN                                                                   
      END                                                                      
