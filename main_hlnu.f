************************************************************************
********************|-------------------------------|*******************
********************| Program: main_higgs_4f.f      |*******************
********************|-------------------------------|*******************
************************************************************************
******************|----------------------------------|******************
******************| Last modified: 13.04.2010        |******************
******************|----------------------------------|******************
******************|     IACS                         |******************
******************|----------------------------------|******************
************************************************************************
************|--------------------------------------------------|********
************| Parton-level Monte Carlo event generator for a   |********
************| p-pbar collider. This main program reads in the  |********
************| parameters, cuts, etc. and then calls the actual |********
************| integration package. It also flashes various     |********
************| information on the screen if required to do so.  |********
************|--------------------------------------------------|********
************| This is a main program which calls the function  |********
************| ppglgl.f where the integrand is computed.        |********
************|--------------------------------------------------|********
************| Process studied:                                 |********
************|                      _    _                      |********
************| H  --> W* W* --->  f f' f f'                     |********
************|  in SM.                                          |********
************|--------------------------------------------------|********
************|  both t ---> b W                                 |********
************|  both W ---> q qb'                               |********
************| ==> 4b +  4j       + mpT                         |********
************************************************************************
*--   declarations--------------------------------------------------------*
      IMPLICIT REAL*8 (a-h,M,o-z)

      CHARACTER*20:: PARM(20)
      DOUBLE PRECISION::VALUE(20),ALPHASPDF
*......................................................................*
      COMMON /PARTCL/P(5,30),M(2,30)
      common/distrib/var(30)
      COMMON/PDFINFO/ALPHAS,MH,nneg,nev 
*......................................................................*
      EXTERNAL PPHLNU
*......................................................................*
      include 'const.h'
*--   I/O devices---------------------------------------------------------*
      OPEN (UNIT = 1, STATUS = 'OLD', FILE = 'events.in')
      OPEN (UNIT = 8, STATUS = 'OLD', FILE = 'hist.in')
*......................................................................*
      OPEN (UNIT = 2, STATUS ='UNKNOWN', FILE = 'out' )
      OPEN (UNIT = 45, STATUS ='UNKNOWN', FILE = 'dist.out' )
*--   reading in data-----------------------------------------------------*
      READ (1,*) NPT, ITN, IHIST
      READ (1,*) MH,MW
*--   initialising masses, momenta, etc.----------------------------------*
      CALL PINIT

C      q(1)q'(2) -> W*(3) -> l(4) nu (5) h(6)
      
      M(1,6)=MH 

      DO I = 1,30
         M(2,I) = M(1,I)**2
      ENDDO
*----------------Initialise LHAPDF -------------- 

c$$$      PARM(1)='DEFAULT'
c$$$      VALUE(1)=10042

      PARM(1)='NPTYPE'
      VALUE(1)=1
      PARM(2)='NGROUP'
      VALUE(2)=4
      PARM(3)='NSET'
      VALUE(3)=46

      CALL PDFSET(PARM,VALUE)

      PRINT*,'FINISHED CALLING PDF'
      nneg=0
      nev=0

*------------------! Actual Calculation Done Here !-----------------*
      NDIM = 9
      IF(IHIST .EQ. 0 )THEN
         CALL MONTES
     >        (PPHLNU,ndim,40,NPT,ITN,CS,ECS,
     >        .FALSE.,0,.FALSE.,.FALSE.,CHI2)
      ENDIF
      IF(IHIST .NE. 0 )THEN
         CALL MONTES
     >        (PPHLNU,ndim,40,NPT,ITN,CS,ECS,
     >        .TRUE.,0,.FALSE.,.FALSE.,CHI2)
      ENDIF

      CALL PDFSTA

*-------------------------------------------------------------------*
      write(*,*)cs,ecs
      write(2,*)cs,ecs
      print*,'Percent of negative points',nneg*100.0/nev
*--   program ends--------------------------------------------------------*
      STOP
      END
************************************************************************
