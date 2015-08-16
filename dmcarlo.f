*****************************************************************************
*************************|-------------------|*******************************
*************************| Program: DMONTES.F|*******************************
*************************|-------------------|*******************************
*****************************************************************************
***********************|--------------------------|**************************
***********************|  MONTE CARLO INTEGRATOR  |**************************
***********************|--------------------------|**************************
*****************************************************************************
**|    Last modified on October 15, 2004 at Wuerzburg                    |***
**|----------------------------------------------------------------------|***
**|  Now in its current form it gives normalized distribution :          |***
**|     1 (d\sigma/dX)                                                   |***
**|     ---                                                              |***
**|    \sigma                                                            |***
**|----------------------------------------------------------------------|***
**|                         How to call MONTES:                          |***
**| CALL MONTES (CROSS,NDIM,NG,IPNT,ITN,XS,E,.FALSE.,0,.FALSE.,.FALSE.)  |***
**|----------------------------------------------------------------------|***
*****************************************************************************
***|---------------------------------------------------------------------|***
***|            Dummy argument list and explanations.                    |***
***|            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                     |***
***| FUNC     real function returning the value of the integrand         |***
***|          multiplied by the hypervolume of integration. input is an  |***
***|          array of dimension NDIM  containing numbers in the range 0 |***
***|          to 1.  The actual argument values are computed in FUNC.    |***
***|          The vector X may be reset. function must not change the    |***
***|          values in the array. The function must be declared         |***
***|          EXTERNAL in the calling routine.                           |***
***|---------------------------------------------------------------------|***
***| NDIM     dimension of the integral (NDIM .LE. MAXDIM)               |***
***| NG       number of grid points (NG .LE. MAXGRD)                     |***
***| MANY     total number of function evaluations                       |***
***| ITN      number of monte-carlo iterations. Suggested value is 4.    |***
***| TINT     contains the value of the integral on output               |***
***| TERR     contains an estimate of the error on output                |***
***| AVRAGE   logical variable specifying whether average values of each |***
***|          integration variable are to be computed or not. Enables the|***
***|          computation when the value is set to .TRUE.                |***
***|..NOTE... when set to .TRUE., "FUNC" has to reset the  input         |***
***|          variables X(J) to the actual value used in the             |***
***|          function evaluation.                                       |***
***|          plots are made only when AVRAGE is set to .TRUE.           |***
***| INTER    specifies whether DINIT uses interactive                   |***
***|          initialisation or not (=1 means interactive).              |***
***|          if not, then unit 8 is used as input to DINIT.             |***
***| DUMP     logical variable specifying whether importance             |***
***|          info is to be dumped (.true. forces output on              |***
***|          unit 9)                                                    |***
***| GET      logical variable specifying whether importance             |***
***|          info is to be read (.true. forces read from                |***
***|          unit 10)                                                   |***
***|---------------------------------------------------------------------|***
***|                 Unit numbers for files                              |***
***|                 ~~~~~~~~~~~~~~~~~~~~~~                              |***
***| INPUT  (8)  for DINIT, only when AVRAGE=.TRUE. and INTER=0          |***
***| INPUT  (10) importance information read only when GET=.TRUE.        |***
***| OUTPUT (9)  all output from this subroutine                         |***
***| OUTPUT (11) importance information dumped only when DUMP=.TRUE.     |***
***| OUTPUT (13) differential information dumped only when AVRAGE=.TRUE. |***
***|---------------------------------------------------------------------|***
***|                 Internal variables and usage                        |***
***|                 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~                        |***
***| MAXDIM =16   is the maximum dimension that can be handled           |***
***| MAXGRD =64   is the maximum number of grid points allowed per axis  |***
***|              must be a power of 2 because of use in routine points  |***
***| KNT          is the number of function evaluations per iteration    |***
***| PINT         is the value of the integral in the current iteration  |***
***| PERR         is the error estimate in the current iteration         |***
***| CHI2         is a measure of the mutual agreement of the estimates  |***
***| X(1:NDIM)    is the sample point in a unit hypercube, used as input |***
***|              to the function returning the integrand.               |***
***| WGT          is the previous importance information                 |***
***| WGRD(*,*,1)  is the histogram of importance information on each axis|***
***| WGRD(*,*,2)  is the grid position                                   |***
***| AVRG(J)      finally contains the average value of the j-th variable|***
***| AERR(J)      finally contains an estimate of the error in AVRG(J)   |***
***| CHI(J)       finally contains the standard deviation in successive  |***
***|              estimates of AVRG(J)                                   |***
***| AVST(J)      are the successive estimates of AVRG(J)                |***
***| AVER(J)      are the estimated errors in AVST(J)                    |***
***|---------------------------------------------------------------------|***
*****************************************************************************
      SUBROUTINE MONTES(FUNC,NDIM,NG,MANY,ITN,TINT,TERR,       
     >                               AVRAGE,INTER,DUMP,GET,CHI2)
*---------------------------------------------------------------------------*
      implicit double precision (a-h,o-z) 
C      implicit real*8(a-h,o-z) 
      PARAMETER (MAXDIM=16, MAXGRD=64,MDIM=16,MBIN=2000,MCOR=5)       
      PARAMETER (EPSLON=1.D-30)              
      DIMENSION X(MAXDIM),DUM(MAXGRD)        
      DIMENSION AVRG(MAXDIM), AERR(MAXDIM), CHI(MAXDIM),      
     >    AVST(MAXDIM), AVER(MAXDIM)         
      COMMON /TUPL/ YTPL(MAXDIM),IPTR(MAXDIM)
      COMMON /WORK/ WGRD(MAXGRD,MAXDIM,2)    
      COMMON /AVTL/ AVRG, AERR            
      COMMON /BINS/ LBIN(MDIM), LCOR(2,MCOR), XL(MDIM), XU(MDIM), 
     >              DX(MDIM), Y(MBIN,MDIM), COR(MBIN,MBIN,MCOR), 
     >              LOFL(MDIM), LUFL(MDIM), NPLOT, NCOR, NHIT, KALL 
      COMMON /NAME/ ID(MDIM)  
      COMMON /STYL/ LOGY(MDIM)
      common/dist_index/idmtt
      COMMON/EVENT/IT,NITER

      LOGICAL AVRAGE, DUMP, GET    
      CHARACTER*10 ID 
      DOUBLE PRECISION YTPL
      EXTERNAL FUNC 
*--------------------------------------------------------------------------*
* Check whether input variables are sensible              
*--------------------------------------------------------------------------*
      IF((NDIM.LE.0).OR.(NDIM.GT.MAXDIM)) STOP '# DIM OUT OF RANGE'
      IF((NG.LE.0).OR.(NG.GT.MAXGRD)) STOP '# GRID PTS OUT OF RANGE'
      IF(ITN .LE.0) STOP '# ITERATIONS .LE. ZERO'             
      IF(MANY.LE.ITN) STOP 'LESS THAN 1 FN EVALUATION PER ITERATION' 
*     WRITE (6,*) NDIM, ITN, MANY, NG        
*-------------------------------------------------------------------*
* Initiation of grids  
*-------------------------------------------------------------------*
C -- ND++
      NITER=ITN
C -- ND--      
      KNT=(MANY-1)/ITN
      KNTOT=0         
      XNG=NG          
      VOL=1.d0/MANY    
      TINT=0.d0         
      TERR=0.d0         
      CHI2=0.d0         
      DO 10 I=1,MAXDIM
         AVRG(I)=0.d0   
         AERR(I)=0.d0   
         CHI(I)=0.d0    
 10      YTPL(I)=0.D0 
      DO 11 I=1,NG    
      DO 11 J=1,NDIM  
         WGRD(I,J,1) =0    
 11      WGRD(I,J,2)=I/XNG 
      IF (GET) THEN
         READ (10,*) NTG, NTDIM 
         READ (10,*) ((WGRD(I,J,2),I=1,NTG),J=1,NTDIM)
c        write(*,*) ntg,ntdim
c        write(*,*) ((wgrd(i,j,2),i=1,ntg),j=1,ntdim)
      END IF  
      NPLOT=0
      KALL=0
      IF(AVRAGE) CALL DINIT(INTER,NPLOT)
          WRITE(9,3000) ITN,KNT,MANY,NG          
*-------------------------------------------------------------------*
* Compute the integral of the function   
*-------------------------------------------------------------------*
      DO 30 IT=1,ITN  
         PRINT*,'STARTING ITERATION',IT
         PINT=0.d0      
         PERR=0.d0      
         DO 12 J=1,NPLOT   
            AVST(J)=0.d0
 12         AVER(J)=0.d0
         DO 31 KN=1,KNT    
*-------------------------------------------------------------------*
* Locate sample point in the hypercube              
*-------------------------------------------------------------------*
            WGT=VOL
            KNTOT=KNTOT+1  
            CALL POINTS(NDIM,KNTOT,XNG,WGT,X)
*-------------------------------------------------------------------*
* Get new weight by function evaluation             
*-------------------------------------------------------------------*
            WGT=WGT*FUNC(X,NDIM)
*-------------------------------------------------------------------*
* Dump previous importance information in histogram 
*-------------------------------------------------------------------*
            PINT=PINT+WGT  
            PERR=PERR+WGT*WGT 
            DO 33 J=1,NDIM 
 33            WGRD(IPTR(J),J,1)=WGRD(IPTR(J),J,1)+DABS(WGT)   
            IF(.NOT.AVRAGE) GO TO 31         
            CALL DSTUFF(X,WGT)
            DO 34 J=1,NPLOT
               AVST(J)=AVST(J)+WGT*X(J)      
 34            AVER(J)=AVER(J)+WGT*WGT*X(J)*X(J)              
 31         CONTINUE  
*-------------------------------------------------------------------*
* Update the values of PINT, PERR, TINT, TERR and CHI2 
*-------------------------------------------------------------------*
         PINT=PINT/VOL/KNT 
c	    write(*,*)'pint,wgt=',pint,wgt
         PERR=(PERR/VOL/VOL-KNT*PINT*PINT)/(KNT-1)/KNT        
         IF (PERR.EQ.0.) THEN 
            TINT=PINT 
           WRITE(9,*) 'UNIFORM DISTRIBUTION IN INTEGRAND'    
            RETURN    
         END IF       
         TINT=TINT+PINT*(PINT*PINT/PERR)     
         TERR=TERR+(PINT*PINT/PERR)          
         CHI2=CHI2+PINT*PINT*(PINT*PINT/PERR)
c       write(*,*) it, pint
         WRITE(9,3001) IT,PINT,DSQRT(PERR), (ITN*AVST(J)/PINT,     
     >    ITN*DSQRT(AVER(J)+PERR*(AVST(J)/PINT)**2)/PINT,J=1,NPLOT)
         IF(.NOT.AVRAGE) GO TO 42            
         DO 36 J=1,NPLOT   
            AVST(J)=AVST(J)/VOL/KNT          
            AVER(J)=(AVER(J)/VOL/VOL-KNT*AVST(J)*AVST(J))/(KNT-1)/KNT           
            AVRG(J)=AVRG(J)+AVST(J)*(AVST(J)*AVST(J)/AVER(J)) 
            AERR(J)=AERR(J)+(AVST(J)*AVST(J)/AVER(J))         
 36         CHI(J)=CHI(J)+AVST(J)*AVST(J)*(AVST(J)*AVST(J)/AVER(J))             
 42      CONTINUE     
*--------------------------------------------------------------------------*
* Refine grids  --- not done on last iteration         
*--------------------------------------------------------------------------*
         IF(IT.EQ.ITN) GO TO 30              
         DO 29 J=1,NDIM    
*--------------------------------------------------------------------------*
* Update importance information in WGRD(*,*,1)      
*--------------------------------------------------------------------------*
* Weights in 3 adjacent bins are added and the new weights normalized to 1.*
* these weights are transformed by a function                              *
*                 F(W) = ( (W-1)/ALOG(W) )**(3/2)                          * 
* and the average of the transformed weights is found. The grid is reset   *
* so that each bin in the new grid is expected to contain a weight equal   *
* to this average.                                                         *
*--------------------------------------------------------------------------*
            RC=0.d0     
            XN=0.d0     
            DR=0.d0     
            K=0       
            DUM(1)=0.5d0*(WGRD(1,J,1)+WGRD(2,J,1))
            DUM(NG)=0.5d0*(WGRD(NG-1,J,1)+WGRD(NG,J,1))          
            X(J)=DUM(1)+DUM(NG)              
            DO 21 I=2,NG-1 
               DUM(I)=(WGRD(I-1,J,1)+WGRD(I,J,1)+WGRD(I+1,J,1))/3.d0 
 21            X(J)=X(J)+DUM(I)              
            DO 22 I=1,NG   
               XO=EPSLON+DUM(I)/X(J)         
               WGRD(I,J,1)=((XO-1)/DLOG(XO))**1.5d0
 22            RC=RC+WGRD(I,J,1)             
            RC=RC/XNG 
*--------------------------------------------------------------------------*
            DO 26 I=1,NG-1 
 24            IF(DR.GT.RC) GO TO 25         
               K=K+1    
               DR=DR+WGRD(K,J,1)          
               XO=XN    
               XN=WGRD(K,J,2)             
               GO TO 24 
 25            DR=DR-RC    
 26            DUM(I)=XN-(XN-XO)*DR/WGRD(K,J,1)
            XO=0.     
            DO 27 I=1,NG   
 27            WGRD(I,J,1)=0.d0 
            WGRD(NG,J,2)=1.d0
            DO 28 I=1,NG-1 
 28            WGRD(I,J,2)=DUM(I)            
 29         CONTINUE  
 30      CONTINUE     
*--------------------------------------------------------------------------*
* Compute the value of the integral      
*--------------------------------------------------------------------------*
      TINT=TINT/TERR  
c     cs = tint
      CHI2=(CHI2/TINT/TINT-TERR)*(ITN-1)/(ITN-0.999d0)**2       
      TERR=TINT/DSQRT(TERR) 
       WRITE(9,3002) TINT,TERR,CHI2           
cc       WRITE(*,*) 'Chi Squared =', CHI2           
      IF (DUMP) THEN
         WRITE (11,*) NG, NDIM
         WRITE (11,*) ((WGRD(I,J,2),I=1,NG),J=1,NDIM)         
      END IF          
      IF(.NOT.AVRAGE) RETURN  
      DO 45 J=1,NPLOT 
         AVRG(J)=AVRG(J)/AERR(J)             
         CHI(J)=(CHI(J)/AVRG(J)/AVRG(J)-AERR(J))*(ITN-1)
     &                    /(ITN-.999d0)**2
         AERR(J)=1/DSQRT(AERR(J))             
         AVRG(J)=AVRG(J)/TINT 
         AERR(J)=DABS(AVRG(J))*DSQRT(AERR(J)*AERR(J)+(TERR/TINT)**2)    
        WRITE(9,3003) J,AVRG(J),AERR(J),CHI(J)
45      CONTINUE
        CALL DHIST(MANY)
      RETURN          
*--------------------------------------------------------------------------*
* Format statements    
*-------------------------------------------------------------------*
 3000 FORMAT('S'/'1 MONTE CARLO: ',I2,' CYCLES,',I5,          
     > ' POINTS/CYCLE,',I6,' TOTAL POINTS,',I3,'GRIDS.'/,     
     > '0',8X,'INTEG',5X,'ERROR')            
 3001 FORMAT(2X,I2,2X,2G10.3/(26X,5X,2G12.4))
 3002 FORMAT('0 AVERAGES:      DIM',4X,'MEAN',9X,'ERR',8X,'CHI SQR'
     > /,16X,'(INT)',1X,G12.4,2X,G12.4,2X,G10.4)              
 3003 FORMAT(16X,I3,3X,G12.4,2X,G12.4,2X,G10.4)
*--------------------------------------------------------------------------*
      END             
*==========================================================================*
      SUBROUTINE POINTS(ND,N,XU,WGT,X)       
      implicit real*8(a-h,o-z)
      PARAMETER (MAXDIM=16, MAXGRD=64, NPOWER=6,MAXDI1=MAXDIM+1) 
      INTEGER PRIME(1:MAXDIM), ADD(1:MAXDI1), DELTA(1:381)  
      DOUBLE PRECISION Y   
      REAL*8 X(1:MAXDIM)
      COMMON /TUPL/ Y(1:MAXDIM), IPTR(1:MAXDIM)
      common/dist_index/idmtt
      COMMON /WORK/ WGRD(1:MAXGRD,1:MAXDIM,1:2)
      DATA PRIME/2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53/ 
      DATA ADD/1,3,6,11,18,29,42,59,78,101,130,161,198,239,282,329,382/ 
      DATA (DELTA(I),I=1,366)/-1,1,-2,1,1,-2,3,-2,3,-2,-3,4,-2,4,-5,
     >4,-2,-4,5,3,-6,8,-7,3,-5,8,-2,-3,-7,6,4,-8,6,-4,8,-11,8,-4,6,-8,4,
     >-9,8,5,-10,8,-6,11,-15,9,-3,7,-10,8,-10,13,-9,3,        
     >-13,9,5,-11,14,-11,5,-10,14,-8,5,-8,14,-10,-6,14,-6,-5,8,-12,11,6,
     >-13,16,-13,6,-11,20,-13,6,-10,13,-17,13,-4,11,-15,10,-13,16,-11,4,
     >-21,15,-8,17,-13,9,-18,25,-18,9,-14,18,-9,13,-21,11,-6,13,-22,18, 
     >9,-22,8,3,-14,22,-13,-4,13,-13,15,8,-18,22,-18,9,-16,27,-17,8,-13,
     >18,-14,6,-14,27,-16,8,-21,20,-13,18,-16,6,12,-24,15,-13,18,-11,   
     >-21,18,10,-22,17,-12,23,-31,22,-11,17,-23,12,16,-35,15,11,-17,12, 
     > -9,19,-28,25,-12,-10,28,-16,-17,24,-14,18,-21,15,-9,18,-28,+16,  
     >-19,20,11,-24,19,-14,26,-35,20,11,-20,3,12,-24,35,-30,14,-23,34,  
     >-17,10,-19,24,-18,6,-17,33,-24,17,-22,31,-17,-20,25,-11,16,-21,14,
     >-19,30,-17,-16,21,11,-25,31,-25,12,-22,32,-18,11,-18,31,-36,18,7,-
     >15,22,-36,18,14,-22,15,16,-34,10,11,-25,35,-25,8,12,-28,18,-12,28,
     > -38,29,-11,7,-18,27,-20, -23,24,-12,27,-33,27,-13,24,-41, 26,-13,
     >20,-26,32,-20,-14,23,-5,-12,32,-45,34,-17,10,-23,35,-21,18,-26,   
     >14,18,-39,26,-15,19,-25,36,-24,-19,30,-15,24,-28,14,-20,31,-15,   
     >-25,26,14,-31,24,-17,33,-45,32,-15,24,-33,17,-23,45,-28,15,-24,29,
     >-42,29,-11,28,-37,24,-10,18,-39,24,25,-37,3,21,-32,39,-15,-20,24/ 
      DATA (DELTA(I),I=367,381)/-15,28,-46,26,13,-33,14,28,-37,19,-15,
     >27,-39,32,-12/
      I=ND-1          
      NPTR=ADD(ND+1)  
 10   IP=1            
         M=N          
         NPTR=NPTR-PRIME(I+1) 
 20      CONTINUE     
            NDIV=M/PRIME(I+1) 
            NREM=M-NDIV*PRIME(I+1)           
            M=NDIV    
            IP=IP*PRIME(I+1)  
            Y(I+1)=Y(I+1)+DBLE(DELTA(NREM+NPTR))/IP           
            IF(NREM.EQ.0) GO TO 20           
         XINT=XU*Y(I+1)    
         IM=XINT      
         IPTR(I+1)=IM+1    
         XO=WGRD(IPTR(I+1),I+1,2)            
         IF(IM.EQ.0) GO TO 30 
         XO=XO-WGRD(IM,I+1,2) 
 30      CONTINUE     
         X(I+1)=WGRD(IPTR(I+1),I+1,2)-XO*(IPTR(I+1)-XINT)     
         WGT=WGT*XO*XU
         I=I-1        
         IF(I.GE.0) GO TO 10  
      RETURN          
      END             
*****************************************************************************
****|--------------------------------------------------------------------|***
****| Set of routines for plotting distributions generated by the program|***
****| MONTE. histograms are generated from the points & weights generated|***
****|--------------------------------------------------------------------|***
****| The various sections are-                                          |***
****| DINIT     initializes all variables for the histogram.             |***
****| DSTUFF    stuffs the bins in accordance with the data from monte.  |***
****| DHIST     prints the histogram of distribution.                    |***
****|--------------------------------------------------------------------|***
****| Program variables-                                                 |***
****| BIN(MDIM)    number of bins (must be .LE. MBIN)                    |***
****| COR(MBIN,MBIN,MCOR)                                                |***
****|              holds the correlation information.                    |***
****| DX(MDIM)     width of each bin in the histogram.                   |***
****| HIST(MBIN)   used for printing the histogram.                      |***
****| ID(MDIM)     contains a character variable for histogram           |***
****|              identification                                        |***
****| INTER        =1  if initialisation is to be done interactively     |***
****| LCOR(2,MCOR) holds pointers to the variables to be correlated.     |***
****|              LCOR(1,*) and LCOR(2,*) hold the pointers to be       |***
****|              correlated for the *-th plot.                         |***
****| NHIT         the number of times a non-zero weight is obtained     |***
****| OFL(MDIM)    number of times a hit is scored beyond the upper cut  |***
****| UFL(MDIM)    number of times a hit is scored below the lower cut X |***
****| XL(MDIM)     lower limit for the histogram axis.                   |***
****| XU(MDIM)     upper limit for the histogram axis.                   |***
****| Y(MBIN,*)    contains the weight in each bin.                      |***
****|--------------------------------------------------------------------|***
****| Fixed parameters-                                                  |***
****| MBIN=50      maximum number of bins in the plot.                   |***
****| MCOR =5      maximum number of correlation plots allowed.          |***
****| MCONT =8     number of distinct contours plotted.                  |***
****| MDIM =16     maximum of plots that can be generated; equals the    |***
****|              maximum dimensionality that  MONTE  can handle.       |***
****| MLEV =30     the number of printed lines between top and bottom    |***
****|              of the histogram.                                     |***
****|--------------------------------------------------------------------|***
*****************************************************************************
      SUBROUTINE DINIT(INTER,NDIM)          
      implicit real*8 (a-h,o-z)
      PARAMETER (MBIN=2000, MCOR=5, MDIM=16)   
      CHARACTER*10 ID 
      INTEGER BIN, LCOR, OFL, UFL            
      COMMON /BINS/ BIN(MDIM), LCOR(2,MCOR), XL(MDIM), XU(MDIM), 
     >              DX(MDIM), Y(MBIN,MDIM), COR(MBIN,MBIN,MCOR), 
     >              OFL(MDIM), UFL(MDIM), NNDIM, NCOR, NHIT, KALL 
      COMMON /NAME/ ID(MDIM)  
      COMMON /STYL/ LOGY(MDIM)
      common/dist_index/idmtt
*---------------------------------------------------------------------------*
* Initializations for histogram plots    
*---------------------------------------------------------------------------*
      KALL=KALL+1     
      NHIT=0          
*---------------------------------------------------------------------------*
* Initialisations on first call          
*---------------------------------------------------------------------------*
      IF (KALL.EQ.1) THEN  
         IF (INTER.EQ.1) THEN 
*---------------------------------------------------------------------------*
* Histograms
*---------------------------------------------------------------------------*
            WRITE(*,*) 'NUMBER OF PLOTS (.LE.',MDIM,')'       
            READ(*,*) NDIM 
            IF (NDIM.GT.MDIM) STOP 'NPLOT>MDIM'
            DO 100 J=1,NDIM
               WRITE(*,*) 'AXIS ',J,' TO BE LABELLED'         
               READ(*,'(A10)') ID(J)         
               WRITE(*,*) 'LO CUT, UP CUT, NO OF BINS, LOGY'  
               READ(*,*) XL(J), XU(J), BIN(J), LOGY(J)        
  100          BIN(J)=MIN(BIN(J),MBIN)       
*---------------------------------------------------------------------------*
* Correlation plots 
*---------------------------------------------------------------------------*
            WRITE(*,*) 'NUMBER OF CORRELATIONS (.LE.',MCOR,')'
            READ(*,*) NCOR 
            IF (NCOR.GT.MCOR) STOP 'NCOR>MCOR' 
            IF (NCOR.GT.0) THEN              
            WRITE(*,*) 'PAIRS FOR CORRELATION? (ENTER NUMBERS)'              
            WRITE (*,'(5X,I2,5X,A10)') (J,ID(J),J=1,NDIM)  
            READ (*,*) (LCOR(1,J),LCOR(2,J),J=1,NCOR)      
            END IF    
         ELSE         
*---------------------------------------------------------------------------*
* Histograms
*---------------------------------------------------------------------------*
                 READ(8,*) NDIM 
            IF (NDIM.GT.MDIM) STOP 'NPLOT>MDIM'
            DO 110 J=1,NDIM
               READ(8,'(A10)') ID(J)         
               READ(8,*) XL(J), XU(J), BIN(J), LOGY(J)        
  110          BIN(J)=MIN(BIN(J),MBIN)       
*---------------------------------------------------------------------------*
* Correlation plots 
*---------------------------------------------------------------------------*
            READ(8,*) NCOR 
            IF (NCOR.GT.MCOR) STOP 'NCOR>MCOR' 
            IF (NCOR.GT.0) READ (8,*) (LCOR(1,J),LCOR(2,J),J=1,NCOR)
         END IF       
      END IF          
*---------------------------------------------------------------------------*
* Initialise accumulators to zero        
*---------------------------------------------------------------------------*
      DO 120 JDIM=1,NDIM   
         UFL(JDIM)=0  
         OFL(JDIM)=0  
         DX(JDIM)=(XU(JDIM)-XL(JDIM))/BIN(JDIM)
         DO 120 K=1,BIN(JDIM) 
 120        Y(K,JDIM)=0.d0   
      DO 121 JCOR=1,MCOR   
      DO 121 K=1,MBIN 
      DO 121 J=1,MBIN 
 121     COR(J,K,JCOR)=0.d0  
      RETURN          
      END             
*==========================================================================*
      SUBROUTINE DSTUFF(X,WEIGHT)            
      implicit real*8 (a-h,o-z)
      PARAMETER (MBIN=2000, MCOR=5, MDIM=16)   
      DIMENSION X(MDIM)    
      INTEGER BIN, LCOR, OFL, UFL, K(MDIM)   
      common/dist_index/idmtt
      COMMON /BINS/ BIN(MDIM), LCOR(2,MCOR), XL(MDIM), XU(MDIM), 
     >              DX(MDIM), Y(MBIN,MDIM), COR(MBIN,MBIN,MCOR), 
     >              OFL(MDIM), UFL(MDIM), NDIM, NCOR, NHIT, KALL 
*---------------------------------------------------------------------------*
* Stuffing action only if weight is non-zero.             
*---------------------------------------------------------------------------*
      IF(WEIGHT.LE.0.d0) RETURN 
      NHIT=NHIT+1     
*---------------------------------------------------------------------------*
* Histograms      
*---------------------------------------------------------------------------*
      DO 200 JDIM=1,NDIM   
         IF(X(JDIM).LT.XL(JDIM)) THEN        
            UFL(JDIM)=UFL(JDIM)+1            
            K(JDIM)=0 
         ELSE IF(X(JDIM).GT.XU(JDIM)) THEN   
            OFL(JDIM)=OFL(JDIM)+1            
            K(JDIM)=0 
         ELSE         
            K(JDIM)=1+(X(JDIM)-XL(JDIM))/DX(JDIM)             
            Y(K(JDIM),JDIM)=Y(K(JDIM),JDIM)+WEIGHT            
         END IF       
 200     CONTINUE     
*---------------------------------------------------------------------------*
* Correlations    
*---------------------------------------------------------------------------*
      DO 210 JCOR=1,NCOR   
         I=LCOR(1,JCOR)    
         J=LCOR(2,JCOR)    
         IF ((K(I).GT.0).AND.(K(J).GT.0))    
     >      COR(K(I),K(J),JCOR)=COR(K(I),K(J),JCOR)+WEIGHT    
 210     CONTINUE     
      RETURN          
      END             
*==========================================================================*
      SUBROUTINE DHIST(MANY)  
      implicit real*8 (a-h,o-z)
      PARAMETER (MBIN=2000, MCONT=8, MCOR=5, MDIM=16, MLEV=10)  
      PARAMETER (NFILL=MCONT+1)              
      CHARACTER*1 HIST(MBIN), FILL(NFILL)    
      CHARACTER*10 ID 
      INTEGER BIN, LCOR, OFL, UFL, BINL(MDIM), BINU(MDIM)     
      LOGICAL MINIM, DRAW  
      COMMON /BINS/ BIN(MDIM), LCOR(2,MCOR), XL(MDIM), XU(MDIM), 
     >              DX(MDIM), Y(MBIN,MDIM), COR(MBIN,MBIN,MCOR), 
     >              OFL(MDIM), UFL(MDIM), NDIM, NCOR, NHIT, KALL 
      COMMON /NAME/ ID(MDIM)  
      common/dist_index/idmtt
      COMMON /STYL/ LOGY(MDIM)
      DATA FILL/'.','1','2','3','4','5','6','7','8'/          
*---------------------------------------------------------------------------*
* Check for empty distribution           
*---------------------------------------------------------------------------*
*     WRITE(13,1300) KALL, MANY              
      IF(NHIT.EQ.0) THEN   
          WRITE(9,9000) ID(NDIM)              
         RETURN       
      END IF          
*---------------------------------------------------------------------------*
* Construct the histogram 
*---------------------------------------------------------------------------*
      DO 310 JDIM=1,NDIM   
*---------------------------------------------------------------------------*
* Figure out the maximum in a bin and the limits on variables            
*---------------------------------------------------------------------------*
         DELTA=0.d0     
         MINIM=.FALSE.
         DO 311 K=1,BIN(JDIM) 
            IF (.NOT.MINIM.AND.(Y(K,JDIM).GT.0.d0)) THEN        
            MINIM=.TRUE.
            BINL(JDIM)=K
            END IF    
            IF (Y(K,JDIM).GT.0) BINU(JDIM)=K 
 311        DELTA=MAX(DELTA,Y(K,JDIM))       
            WRITE(9,9001) ID(JDIM), BIN(JDIM), XL(JDIM), XU(JDIM),  
     >      XL(JDIM)+(BINL(JDIM)-1)*DX(JDIM), XL(JDIM)+DX(JDIM)* 
     >      (BINU(JDIM)+1), MANY, NHIT, UFL(JDIM), OFL(JDIM), DELTA             
         IF (LOGY(JDIM).EQ.0) THEN           
            DELTA=DELTA/MLEV  
         ELSE         
            ZERO=LOG10(DELTA)-LOGY(JDIM)     
            DELTA=FLOAT(LOGY(JDIM))/MLEV     
         END IF       
*---------------------------------------------------------------------------*
* Write out the distribution in the histogram          
*---------------------------------------------------------------------------*
*.....Distribution in : d\sigma form (un-normalized):
*---------------------------------------------------------------------------*
         WRITE(34,1301) ID(JDIM), XL(JDIM)+DX(JDIM)*(BINL(JDIM)-0.5d0), 
     >      XL(JDIM)+DX(JDIM)*(BINU(JDIM)-0.5d0), DX(JDIM)      
         WRITE(34,1302) (XL(JDIM) + DX(JDIM)*(K - 0.5d0),Y(K,JDIM),
     >                        K=BINL(JDIM),BINU(JDIM))   
*---------------------------------------------------------------------------*
*.....Un-normalized Distribution in : (d\sigma/dx) form:
*---------------------------------------------------------------------------*
c      if (idmtt .eq. 1) then
         WRITE(45,1301) ID(JDIM), XL(JDIM)+DX(JDIM)*(BINL(JDIM)-0.5d0), 
     >      XL(JDIM)+DX(JDIM)*(BINU(JDIM)-0.5d0), DX(JDIM)      
         WRITE(45,1302) (XL(JDIM) + DX(JDIM)*(K - 0.5d0),
     >                  Y(K,JDIM)/DX(JDIM),
     >                        K=BINL(JDIM),BINU(JDIM))   
c       elseif (idmtt .eq. 2)then
c        WRITE(46,1301) ID(JDIM), XL(JDIM)+DX(JDIM)*(BINL(JDIM)-0.5d0), 
c    >      XL(JDIM)+DX(JDIM)*(BINU(JDIM)-0.5d0), DX(JDIM)      
c        WRITE(46,1302) (XL(JDIM) + DX(JDIM)*(K - 0.5d0),
c    >                  Y(K,JDIM)/DX(JDIM),
c    >                        K=BINL(JDIM),BINU(JDIM))   
c       elseif (idmtt .eq. 3)then
c        WRITE(47,1301) ID(JDIM), XL(JDIM)+DX(JDIM)*(BINL(JDIM)-0.5d0), 
c    >      XL(JDIM)+DX(JDIM)*(BINU(JDIM)-0.5d0), DX(JDIM)      
c        WRITE(47,1302) (XL(JDIM) + DX(JDIM)*(K - 0.5d0),
c    >                  Y(K,JDIM)/DX(JDIM),
c    >                        K=BINL(JDIM),BINU(JDIM))   
c       endif
*---------------------------------------------------------------------------*
* Print the histogram  : un-normalized: d\sigma 
*---------------------------------------------------------------------------*
           DO 312 J=1,MBIN   
312         HIST(J)=' '    
         HIST(BIN(JDIM)+1)='.'
         WRITE(9,9002) ('.',K=1,BIN(JDIM)+1) 
         DO 313 LEVEL=MLEV,1,-1              
            VTEST=(LEVEL+.5)*DELTA           
           DO 314 K=1,BIN(JDIM)             
             IF (LOGY(JDIM).EQ.0) THEN     
                DRAW=Y(K,JDIM).GT.VTEST    
               ELSE   
                IF(Y(K,JDIM).GT.0) THEN    
                  DRAW=(DLOG10(Y(K,JDIM))-ZERO).GT.VTEST    
                 ELSE
                   DRAW=.FALSE.            
                 END IF   
                 END IF 
                 IF(DRAW) THEN  
                    HIST(K)='#' 
                 ELSE   
                    HIST(K)=' ' 
                 END IF 
 314           CONTINUE    
            WRITE(9,9002) (HIST(K),K=1,BIN(JDIM)+1) 
            WRITE(9,9002) (HIST(K),K=1,MBIN) 
 313        CONTINUE  
             WRITE(9,9002) ('.',K=1,BIN(JDIM)+1) 
 310     CONTINUE     
*---------------------------------------------------------------------------*
* Construct the correlation plot         
*---------------------------------------------------------------------------*
      DO 320 JCOR=1,NCOR   
         IX=LCOR(1,JCOR)   
         IY=LCOR(2,JCOR)   
*---------------------------------------------------------------------------*
* The maximum weight
*---------------------------------------------------------------------------*
         DELTA=0.d0     
         DO 321 K=1,BIN(IY)
         DO 321 J=1,BIN(IX)
  321       DELTA=MAX(DELTA,COR(J,K,JCOR))   
          WRITE(9,9005) ID(IX), ID(IY), DELTA, MCONT           
         DELTA=DELTA/MCONT 
*---------------------------------------------------------------------------*
* Print the contour map of correlations 
*---------------------------------------------------------------------------*
         DO 322 J=1,MBIN   
 322        HIST(J)=' '    
         HIST(BIN(IX)+1)='.'  
          WRITE(9,9002)('.',J=1,BIN(IX)+1)    
         DO 323 K=BIN(IY),1,-1
            DO 324 J=1,BIN(IX)
                IF (COR(J,K,JCOR).GT.0.) THEN 
                   IH=NINT(COR(J,K,JCOR)/DELTA)+1              
                   HIST(J)=FILL(IH)           
                ELSE   
                   HIST(J)=' ' 
                END IF 
 324           CONTINUE    
323        WRITE(9,9002) (HIST(J),J=1,MBIN) 
320        WRITE(9,9002)('.',J=1,BIN(IX)+1)    
C323        CONTINUE
C320        CONTINUE
*---------------------------------------------------------------------------*
* Print out the double differential cross sections        
*---------------------------------------------------------------------------*
      DO 340 JCOR=1,NCOR   
         IX=LCOR(1,JCOR)   
         IY=LCOR(2,JCOR)   
         WRITE(13,1305) ID(IX), XL(IX)+DX(IX)*(BINL(IX)-0.5d0), 
     >      XL(IX)+DX(IX)*(BINU(IX)-0.5d0), DX(IX), ID(IY)      
         DO 340 K=BINL(IY),BINU(IY)          
            WRITE(13,1306) ID(IY),XL(IY)+DX(IY)*(K-0.5d0)       
340         WRITE(13,1307)(COR(J,K,JCOR),J=BINL(IX),BINU(IX))
*---------------------------------------------------------------------------*
      RETURN          
*---------------------------------------------------------------------------*
 9000 FORMAT(//1X,19(1H*),'DISTRIBUTION ',A10,'IS EMPTY'/)    
 9001 FORMAT(5X,'DIFFERENTIAL CROSS SECTION IN ',A/           
     >       11X,I2,' BINS BETWEEN ',E12.4,' AND ',E12.4/11X, 
     >       ' CROSS SECTION NON-ZERO BETWEEN ',E12.4,' AND ',E12.4/ 
     >       11X,I6,' POINTS, ',I6,' HITS ',I6,' UNDERFLOWS ', I6,   
     >       ' OVERFLOWS'/11X,'MAXIMUM OF DISTRIBUTION=',E12.4)  
 9002 FORMAT(10X,'.',60A1) 
 9005 FORMAT(5X,'CORRELATIONS BETWEEN ',A,' (X-AXIS) AND ',A, 
     > ' (Y-AXIS)'/T11,'ALL INFORMATION SAME AS IN HISTOGRAMS'/  
     > T11,'MAXIMUM IN BIN=',E12.4,' CONTOURS IN STEPS OF ',I2)  
 1300 FORMAT('#',5x, I2,'-TH CALL TO MONTE. NUMBER OF POINTS=',I6) 
 1301 FORMAT(/'#%HISTOGRAM'/'#',A10,T21,E12.4,',',E12.4,',',E12.4)  
 1302 FORMAT(2E14.6)  
 1305 FORMAT('#%CORRELATION'/A10,T21,E12.4,',',E12.4,',',E12.4,','/ 
     >  T21,A,' AS SHOWN') 
 1306 FORMAT('#',A,'=',E12.4)
 1307 FORMAT(5E12.4)
*---------------------------------------------------------------------------*
      END             
*****************************************************************************
