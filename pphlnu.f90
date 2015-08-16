! File: pphlnu.f
! Author: Nishita Desai
! Description: Full program for production of higgsvis W-associated
! production and decay into l nu l nu 
! Revision: $Revision: 1.12 $
! Date: $Date: 2011-04-05 06:27:59 $

! Documentation:
! Asymmetry definitions
! 1. cos(angle) LR
! 2. cos(angle) C
! 3. dPhi LR
! 4. dPhi C
! 5.  d[abs(eta)] LR
! 6. d[abs(eta)] C

program pphlnu

  implicit none

  double precision::p(5,30),m(2,30),phi1,phi2,sphi
  common /partcl/p,m

  double precision::mw,mz,gamw,swsq,ee,mh,s,gamma
  integer::nerr
  logical::dohist
  common/input/mw,mz,gamw,swsq,ee,mh,s,gamma,nerr,dohist

  character*20:: parm(20)
  double precision::value(20)
  double precision::xmin,xmax,q2min,q2max
  COMMON/W50513/XMIN,XMAX,Q2MIN,Q2MAX 

  integer,parameter::nhist=20
  double precision::hist(nhist*100),hmin(nhist),hmax(nhist),totwgt(nhist),wgtmax
  common/histo/hist,hmin,hmax,totwgt,wgtmax
  integer::iof(nhist),iuf(nhist),ndata(nhist)
  common/stat/iof,iuf,ndata

  logical::iSaveWgt,doLHEF,iProd
  common/LHEF/iSaveWgt,doLHEF,iProd

  double precision::wcut(20),wnocut
  common/cutflow/wcut,wnocut

  double precision::hist2(100,100), totwgt2
  common/hist2d/hist2, totwgt2
  integer::i1,i2

  double precision::asym(10),asymnocut(10)
  integer::igen,ievent
  common/asymm/asym,asymnocut,igen
  common/iev/ievent

  double precision::b_in,c_in
  common/input/b_in,c_in

  integer::nev,iev,niter,ndim,nevunw
  double precision::crossx,errorx,mesq,chi2,gev2fb
  integer::i, ip

  external mesq

  ! --- Read Input ---
!!$  read(*,*)b_in,c_in
!!$  print*,b_in," ",c_in

  ! --- Flags ---
  dohist = .true.   ! Histograms
  doLHEF = .false.  ! LHEF files
  iProd  = .false.  ! Production only

  ! -------- input params ----------

  mz = 91.187
  ee = 0.31223 !0.31344 PDG
  swsq = 0.481d0**2.0d0

  mw = mz*dsqrt(1-swsq)
  gamw = 2.089 ! From PDG
  mh = 150.0d0
  gamma = 1.73e-2  ! mh = 150
  !gamma =  4.87e-3 ! mh = 130

  s = 14000**2
  m(:,:)=0.0d0
  wgtmax = 0.0d0
  wcut(:) = 0
  wnocut = 0
  asym(:)= 0
  asymnocut(:)= 0
  igen=0
  ievent=0

  nev = 1000000
  niter = 10
  nevunw = 1000000000

  if(iProd) then
     ndim = 7  ! Production only
     dohist = .false.
     doLHEF = .false. 
  else
     ndim = 15 ! Full run 
  end if

  ! -------------Initialise PDF-------------- 

  ip = 4
  call SetCtq6 (ip) 

  if(dohist) then
  ! --- Book histograms ---
     call inithist(1,dble(0.),dble(200.))   ! pT of 1st Lepton
     call inithist(2,dble(0.),dble(200.))   ! pT of 2nd Lepton
     call inithist(3,dble(0.),dble(200.))   ! pT(1)-pT(2)
     call inithist(4,dble(0.),dble(600.))   ! pT(1)+pT(2)
     call inithist(5,dble(0.),dble(1.0))    ! (pT(1)-pT(2))/(pT(1)+pT(2))
     call inithist(6,dble(0.),dble(600.))   ! Inv. mass of ssd
     call inithist(7,dble(0.),dble(600.0))  ! MET (SSD)
     call inithist(8,dble(0.),dble(1000.0)) ! sumpt+MET (SSD)
     call inithist(9,dble(-1.0),dble(1.0))  ! cos(angle between SSD)
     call inithist(10,dble(-3.15),dble(3.15))! delta-phi
     call inithist(11,dble(-5.0),dble(5.0)) ! Eta of 1st Lepton
     call inithist(12,dble(-5.0),dble(5.0)) ! Eta of 2nd Lepton
     call inithist(13,dble(-5.0),dble(5.0)) ! Delta-abs(Eta )
  end if

  ! -- Integrate --
  p(:,:)=0.0d0
  m(1,5) = mh
  m(2,5) = mh**2

  ! --- First iteration ---
  iSaveWgt = .true.

  call montes(mesq,ndim,50,nev,niter,crossx,errorx, &
       & .false.,0,.false.,.false.,chi2)
  
  ! -- Conversions --
  gev2fb = 3.89d11 !in fb
  
  crossx = crossx * gev2fb 
  errorx = errorx * gev2fb 

  ! --- Normalisations of asymmetries ---
  asymnocut(:)=asymnocut(:)/wnocut
  asym(:)=asym(:)/wcut(3)

  ! -- Print on terminal ---
  print*,'Higgs mass:',mh ,'crossx:',crossx,'+/-',errorx
  if (nerr>0) print*,'Percent of error:',nerr*100.0/nev
  print*,'Maximum weight',wgtmax
  if(.not. iProd) then
     write(*,'("SSD before:",f10.4)')56.0*crossx
     write(*,'("SSD after cuts and multiplicity:",f10.4)')wcut(3)*56.0*crossx/wnocut
     do i=1,10
        write(*,'("Asym(",i3,")=",f10.5)')i,asym(i)
     end do
  end if

  ! -- Write to file ---
  write(998,*) 'Higgs mass:',mh ,'crossx:',crossx,'+/-',errorx
  if (nerr>0) write(998,*) 'Percent of errors:',nerr*100.0/nev
  write(998,*) 'Maximum weight',wgtmax
  if(.not. iProd) then
     write(998,*)'Cut flow table'
     write(998,*)'Total weight:', wnocut
     write(998,*)'No of errors:', nerr
     do i=1,10
        write(998,'("CUT:",i5,f10.5)')i,wcut(i)/wnocut
     end do
     write(998,'("SSD before:",f10.4)')56.0*crossx
     write(998,'("SSD after cuts and multiplicity:",f10.4)')wcut(3)*56.0*crossx/wnocut
     do i=1,10
        write(998,'("Asym(",i3,")=",f10.5,"   Asym_nocut(",i3,")=",f10.5)')i,asym(i),i,asymnocut(i)
     end do
  end if

  if(dohist) then
     do i=1,20
        call histnorm(i,crossx)
        call printhist(1000+i,i)
        call histnorm(i,-1.0d0)
        call printhist(2000+i,i)
     enddo
  end if

  if(doLHEF) then
     dohist = .false.

     ! --- Write LHEF Header ---
     write(997,*)'<LesHouchesEvents version="1.0">'
     write(997,*)'<!-- Generated for pp to hlnu '
     write(997,*)' Author: Nishita Desai -->'
     write(997,*)'<init>'
     write(997,*)'2212   2212   7.0E+3  7.0E+3  0  0  10042 10042  3  1'
     write(997,*)crossx, errorx, wgtmax, "1001"
     write(997,*)'</init>'

     ! --- Second iteration ---
     iSaveWgt = .false.
     
     nev = nevunw
     niter = 1
     do i=1,10
        call montes(mesq,ndim,50,nev,niter,crossx,errorx, &
             & .false.,0,.false.,.false.,chi2)
     end do
     write(997,*)'</LesHouchesEvents>'

  end if

end program pphlnu

!!$-------------------------------------------------------------------

double precision function mesq(rf,ndim)

  implicit none
  integer,intent(in)::ndim
  double precision,intent(in)::rf(16)

  double precision::p(5,30),m(2,30)
  common /partcl/p,m

  double precision::s1,s2,s3,zz
  integer::iflag
  common /dk3cmn/s1,s2,s3,zz(5),iflag                                      

  double precision::yy(2),wgtdk21,wgtdk22
  integer::iflag2
  COMMON/DK2CMN/YY,IFLAG2

  double precision::xx(8)
  integer::iflag4(7)
  common/dk4cmn/xx,iflag4

  double precision::mw,mz,gamw,swsq,ee,mh,s,gamma
  integer::nerr
  logical::dohist
  common/input/mw,mz,gamw,swsq,ee,mh,s,gamma,nerr,dohist

  integer::itn,niter
  common/event/itn,niter

  double precision::xmin,xmax,q2min,q2max
  COMMON/W50513/XMIN,XMAX,Q2MIN,Q2MAX 

  logical::iSaveWgt,doLHEF,iProd
  common/LHEF/iSaveWgt,doLHEF,iProd

  integer,parameter::nhist=20
  double precision::hist(nhist*100),hmin(nhist),hmax(nhist),totwgt(nhist),wgtmax
  common/histo/hist,hmin,hmax,totwgt,wgtmax
  integer::iof(nhist),iuf(nhist),ndata(nhist)
  common/stat/iof,iuf,ndata

  double precision::wcut(20),wnocut
  common/cutflow/wcut,wnocut

  double precision::hist2(100,100), totwgt2
  common/hist2d/hist2,totwgt2
  integer::i1,i2

  double precision::asym(10),asymnocut(10)
  integer::igen,ievent
  common/asymm/asym,asymnocut,igen
  common/iev/ievent


  integer::i,j
  double precision:: x1,x2,scl,e,mesq1,mesq2,x1min,x2min
  double precision:: u1,d1,ub1,db1,sb1,ch1,bt1,tp1,gl1
  double precision:: u2,d2,ub2,db2,sb2,ch2,bt2,tp2,gl2
  double precision:: p1_p3,p2_p4,p1_p4,p2_p3,p1_p2,p3_p4, denom, eps, ep
  double precision:: p1(5),p2(5),p3(5),p4(5),p6(5),p5(5),a,b,c,p3sq,p4sq,p5sq
  double precision:: p8_p9,p8_p10,p8_p11,p9_p10,p9_p11,p10_p11,p8(5),p9(5),p10(5),p11(5)
  double precision:: pi,wps3,g,flux,norm,x,const,pdf1,pdf2
  double precision:: mdecy,wps4,decaydenom,dknorm,dkconst
  double precision:: pt(2),eta(2),phi(2),dphi,dabseta,ptasym
  double precision:: pt1,pt2,phi1,phi2,eta1,eta2,p3_p8,angle
  double precision:: rj1,rj2,phij1,phij2,etaj1,etaj2,rj1min,rj2min
  double precision:: psum(5),sl1l2,pnorm(2),sigma, tmp
  double precision:: metssd,ran,term2,term1,wgt,dphi1,dphi2
  ! External functions
  double precision:: sphi, Ctq6Pdf

  logical:: dodecay,allowssd,swapUD,isplus

  igen = igen+1

  dodecay = .false.
  swapUD = .false.
  isplus = .true.

  if(.not. iProd)  then
     dodecay = .true.
     allowssd = .true.
  end if
  ! -- begin: Initialisation --
  pi = atan(1.0)*4.0d0
  mesq1 = 0.0d0
  mesq2 = 0.0d0
  mesq  = 0.0d0

  b = 0.20d0/(mw**2)
  c = 0.0d0/(mw**2)
  e = dsqrt(s)/2.0d0
  g = ee/dsqrt(swsq)

  ! -- begin: PDFs --
  x1min = mh**2/s
  x1 = x1min+(1.0-x1min) * rf(1)
  x2min = mh**2/(s*x1)
  x2 = x2min + (1.0-x2min)* rf(2)

  scl = sqrt(s*x1*x2)

  p(1,1)= 0.0d0; p(2,1)= 0.0d0; p(3,1)= e*x1; p(4,1)= e*x1; p(5,1)= p(4,1)
  p(1,2)= 0.0d0; p(2,2)= 0.0d0; p(3,2)=-e*x2; p(4,2)= e*x2; p(5,2)= p(4,2)
  p(:,6)= p(:,1)+p(:,2)
  p(5,6)= abs(p(3,6))
  call dot_4(6,6, m(2,6))
  m(1,6) = dsqrt(m(2,6))
  call dot_4(7,7, m(2,7))
  m(1,7) = dsqrt(m(2,7))


  ! -- begin: 3-body decay of W* --

  do i=1,5 ! uses upto rf(7)
     zz(i)=rf(i+2) ! 
  enddo

  ! W*->h(W->l nu)
  call dk3(6,5,3,4,wps3)

  p1(:)=p(:,1);p2(:)=p(:,2);p3(:)=p(:,3);p4(:)=p(:,4)
  p5(:)=p(:,5);p6(:)=p(:,6)

  p(:,7) = p(:,3)+p(:,4)
  p(5,7) = sqrt(p(1,7)**2+p(2,7)**2+p(3,7)**2)

  if(iflag/=0) then
     mesq=0.0d0
     print*,'Point discarded',igen,ievent
     print*, wps3, (rf(i),i=1,15)
     stop
     return
  end if

  call dot_4(1,3,p1_p3)
  call dot_4(2,4,p2_p4)
  call dot_4(1,4,p1_p4)
  call dot_4(2,3,p2_p3)
  call dot_4(1,2,p1_p2)
  call dot_4(3,4,p3_p4)

  ! Let the u come from either proton with 50 percent probability
  ran = rand(0)
  if(ran < 0.5) then
     swapUD = .true.
     tmp = x1
     x1 = x2
     x2 = tmp

     ! Exchange 1<->2
     call dot_4(1,3,p2_p3)
     call dot_4(2,4,p1_p4)
     call dot_4(1,4,p2_p4)
     call dot_4(2,3,p1_p3)
  end if

  u1  = Ctq6Pdf ( 1, x1, scl)
  ub1 = Ctq6Pdf (-1, x1, scl)
  db2 = Ctq6Pdf (-2, x2, scl)
  d2  = Ctq6Pdf ( 2, x2, scl)


  ! For native CTEQ6L1; no adding of sea required.
  pdf1 = u1 * db2 ! W+
  pdf2 = ub1 * d2 ! W-

  denom = ((2.0*p1_p2-mw**2)**2 + (mw*gamw)**2) &
       &* ((2.0*p3_p4-mw**2)**2 + (mw*gamw)**2)
  const = g**6 * mw**2 /(4.0*3.0*16.0)  ! 1/4=> spin avg, 3/9=>color avg
  flux = 0.50d0/m(2,6)                  ! 1/(2 s^hat)
  norm = 1.0d0/(128.0d0 * pi**3)        ! = 1/(2pi)^6 * (pi^3)/2
  const = const*(1.0-x2min)*(1.0-x1min) ! Jacobian factor

  ! SM only
!!$   mesq = 64.0d0 * p1_p4 * p2_p3 * pdf1
!!$   mesq = const * mesq * wps3 * flux * norm / denom 
!!$   return

  eps = ep(p1,p2,p3,p4)

  ! --- Matrix element ---
  ! For W+ 
  mesq1 = + 128*p1_p2*p1_p3*p2_p4*p3_p4*c**2 &
       &  + 128*p1_p2*p1_p3*p2_p4*p3_p4*b**2 &
       &  + 64*p1_p2*p1_p4*p3_p4*c &
       &  + 64*p1_p2*p1_p4**2*p3_p4*c**2 &
       &  + 64*p1_p2*p1_p4**2*p3_p4*b**2 &
       &  - 64*p1_p2*p2_p3*p3_p4*c &
       &  + 64*p1_p2*p2_p3**2*p3_p4*c**2 &
       &  + 64*p1_p2*p2_p3**2*p3_p4*b**2 &
       &  - 64*p1_p2**2*p3_p4**2*c**2 &
       &  - 64*p1_p2**2*p3_p4**2*b**2 &
       &  + 128*p1_p3*p1_p4*p2_p3*p2_p4*c**2 &
       &  + 128*p1_p3*p1_p4*p2_p3*p2_p4*b**2 &
       &  - 64*p1_p3*p1_p4*p2_p4*c &
       &  + 64*p1_p3*p2_p3*p2_p4*c &
       &  - 64*p1_p3**2*p2_p4**2*c**2

  mesq1 = mesq1 - 64*p1_p3**2*p2_p4**2*b**2 &
       &  + 64*p1_p4*p2_p3 &
       &  - 64*p1_p4*p2_p3**2*c &
       &  + 64*p1_p4**2*p2_p3*c &
       &  - 64*p1_p4**2*p2_p3**2*c**2 &
       &  - 64*p1_p4**2*p2_p3**2*b**2 &
       &  - 64*eps*p1_p4*b &
       &  - 64*eps*p2_p3*b

  ! For W-

  mesq2 = + 128*p1_p2*p1_p3*p2_p4*p3_p4*c**2 &
       &  + 128*p1_p2*p1_p3*p2_p4*p3_p4*b**2 &
       &  - 64*p1_p2*p1_p4*p3_p4*c &
       &  + 64*p1_p2*p1_p4**2*p3_p4*c**2 &
       &  + 64*p1_p2*p1_p4**2*p3_p4*b**2 &
       &  + 64*p1_p2*p2_p3*p3_p4*c &
       &  + 64*p1_p2*p2_p3**2*p3_p4*c**2 &
       &  + 64*p1_p2*p2_p3**2*p3_p4*b**2 &
       &  - 64*p1_p2**2*p3_p4**2*c**2 &
       &  - 64*p1_p2**2*p3_p4**2*b**2 &
       &  + 128*p1_p3*p1_p4*p2_p3*p2_p4*c**2 &
       &  + 128*p1_p3*p1_p4*p2_p3*p2_p4*b**2 &
       &  + 64*p1_p3*p1_p4*p2_p4*c &
       &  - 64*p1_p3*p2_p3*p2_p4*c &
       &  - 64*p1_p3**2*p2_p4**2*c**2

  mesq2 = mesq2 - 64*p1_p3**2*p2_p4**2*b**2 &
       &  + 64*p1_p4*p2_p3 &
       &  + 64*p1_p4*p2_p3**2*c &
       &  - 64*p1_p4**2*p2_p3*c &
       &  - 64*p1_p4**2*p2_p3**2*c**2 &
       &  - 64*p1_p4**2*p2_p3**2*b**2 &
       &  + 64*eps*p1_p4*b &
       &  + 64*eps*p2_p3*b


  if(isplus) then
     mesq = mesq1 * pdf1 ! For W+; LHEF
  else
     mesq = mesq2 * pdf2 ! For W-; LHEF
  end if

  mesq = const * mesq * wps3 * flux * norm / denom

  !  -- begin: decay ---

  ! For decay check
  ! p(1,5)=0.0d0;p(2,5)=0.0d0;p(3,5)=0.0d0;p(4,5)=m(1,5);p(5,5)=0.0d0;

  if(dodecay) then
     do i=1,8 !uses upto rf(15)
        xx(i)=rf(i+7) !Change to i for decay only
     enddo

     call dk4(5,8,9,10,11,wps4)

     call dot_4(8,9,p8_p9)
     call dot_4(8,10,p8_p10)
     call dot_4(8,11,p8_p11)
     call dot_4(9,10,p9_p10)
     call dot_4(9,11,p9_p11)
     call dot_4(10,11,p10_p11)

     p8(:)=p(:,8);p9(:)=p(:,9);p10(:)=p(:,10);p11(:)=p(:,11)

     mdecy = + 128*p8_p9*p8_p10*p9_p11*p10_p11*c**2&
          &  + 128*p8_p9*p8_p10*p9_p11*p10_p11*b**2&
          &  - 64*p8_p9*p8_p11*p10_p11*c&
          &  + 64*p8_p9*p8_p11**2*p10_p11*c**2&
          &  + 64*p8_p9*p8_p11**2*p10_p11*b**2&
          &  + 64*p8_p9*p9_p10*p10_p11*c&
          &  + 64*p8_p9*p9_p10**2*p10_p11*c**2&
          &  + 64*p8_p9*p9_p10**2*p10_p11*b**2&
          &  - 64*p8_p9**2*p10_p11**2*c**2&
          &  - 64*p8_p9**2*p10_p11**2*b**2&
          &  + 128*p8_p10*p8_p11*p9_p10*p9_p11*c**2&
          &  + 128*p8_p10*p8_p11*p9_p10*p9_p11*b**2&
          &  + 64*p8_p10*p8_p11*p9_p11*c&
          &  - 64*p8_p10*p9_p10*p9_p11*c&
          &  - 64*p8_p10**2*p9_p11**2*c**2

     mdecy = mdecy - 64*p8_p10**2*p9_p11**2*b**2&
          &  + 64*p8_p11*p9_p10&
          &  + 64*p8_p11*p9_p10**2*c&
          &  - 64*p8_p11**2*p9_p10*c&
          &  - 64*p8_p11**2*p9_p10**2*c**2&
          &  - 64*p8_p11**2*p9_p10**2*b**2&
          &  + 64*ep(p8,p9,p10,p11)*p8_p11*b&
          &  + 64*ep(p8,p9,p10,p11)*p9_p10*b


     decaydenom = ((2.0*p8_p9-mw**2)**2 + (mw*gamw)**2) &
          &     * ((2.0*p10_p11-mw**2)**2 + (mw*gamw)**2)

     dknorm = 1.0d0/((2.0d0*pi)**8) /(m(1,5)*2.0d0)
     dkconst = g**6 * mw**2 / 16.0d0 
     mdecy = dkconst * mdecy * wps4 * dknorm /decaydenom

     ! Fold in the decay part
     mesq = mesq * mdecy / gamma

     ! Decay only
     ! mesq = mdecy

  end if ! decay

  if(.not. mesq > 0.0d0) then
     if(pdf1>0 .and. pdf2>0) then
        if (mesq1<0) print*,'Mesq1<0'
        if (mesq2<0) print*,'Mesq2<0'
        if (mdecy<0) print*,'Mdecy<0'
     end if
     nerr=nerr+1
     mesq=0.0d0
     return
  endif

  if(dodecay) then

     wnocut = wnocut + mesq

     !Smearing for detector effects
     sigma  = 0.05 * dsqrt(p3(4)) + 0.02 * p3(4) 
     call gaussno(pnorm,2)
     p3(1) = p3(1) + pnorm(1) * sigma
     p3(2) = p3(2) + pnorm(2) * sigma

     sigma  = 0.05 * dsqrt(p8(4)) + 0.02 * p8(4) 
     call gaussno(pnorm,2)
     p8(1) = p8(1) + pnorm(1) * sigma
     p8(2) = p8(2) + pnorm(2) * sigma

     ! -- pT(1,2,3)
     pt1 = dsqrt(p3(1)**2 + p3(2)**2)   !l+
     pt2 = dsqrt(p8(1)**2 + p8(2)**2)   !l+

     eta1 = 0.5d0 * log((p(4,3)+p(3,3))  / (p(4,3)-p(3,3)))
     eta2 = 0.5d0 * log((p(4,8)+p(3,8))  / (p(4,8)-p(3,8)))

     ! pT-ordering
     if(pt1>pt2) then
        pt(1) = pt1
        pt(2) = pt2
        eta(1) = eta1
        eta(2) = eta2
     else
        pt(1) = pt2
        pt(2) = pt1
        eta(1) = eta2
        eta(2) = eta1
     end if

     phi1 = sphi(p3(1),p3(2))
     phi2 = sphi(p8(1),p8(2))

     ! --- Order by picking the one closest to the jets as the one
     ! --- from higgs decay.

     rj1min = 10; rj2min=10
     do i=10,11
        etaj1 = 0.5d0 * log((p(4,i)+p(3,i))  / (p(4,i)-p(3,i)))
        phij1 = sphi(p(1,i),p(2,i))
        etaj2 = 0.5d0 * log((p(4,i)+p(3,i))  / (p(4,i)-p(3,i)))
        phij2 = sphi(p(1,i),p(2,i))
        dphi1 = phij1-phi1
        dphi2 = phij2-phi2
        if(abs(dphi1)>pi) dphi1 = 2*pi-abs(dphi1)
        if(abs(dphi2)>pi) dphi2 = 2*pi-abs(dphi2)
        rj1 = dsqrt(dphi1**2+(etaj1-eta1)**2)
        rj2 = dsqrt(dphi2**2+(etaj1-eta2)**2)
        if (rj1<rj1min) rj1min = rj1
        if (rj2<rj2min) rj2min = rj2
     end do

     if(rj1min>rj2min) then
        phi(1) = phi1
        phi(2) = phi2
     else
        phi(1) = phi2
        phi(2) = phi1
     end if

     ! === Variables ===
     ! -- 3: Delta(pT) --
     pTasym = (pt(1)-pt(2))/(pt(1)+pt(2)) 

     ! -- 6: Minv --
     psum(:)=p3(:)+p8(:)
     sl1l2 = psum(4)**2 - psum(1)**2 - psum(2)**2 - psum(3)**2
     if(sl1l2>0) then
        sl1l2 = dsqrt(sl1l2)
     else
        sl1l2=0.d0
     end if

     ! -- 7: MET --
     metssd = dsqrt(p4(1)**2 + p4(2)**2) + dsqrt(p9(1)**2 + p9(2)**2)

     ! -- 9:cos(theta_SSD)
     call dot_3(3,8,p3_p8)
     angle = p3_p8/(p(5,3)*p(5,8))

     ! -- 10: Delta-phi
     dphi = phi(1)-phi(2)
     if(dphi > pi) dphi = dphi - 2.0 * pi
     if(dphi <-pi) dphi = 2.0 * pi + dphi

     ! -- 13:Delta eta
     dabseta = abs(eta(1))-abs(eta(2))

     ! === Plots before cuts ====

     if(dohist) then
        ! ---- 1,2: pT of leptons ---

        call histfill(1,pt(1),mesq)
        call histfill(2,pt(2),mesq)

        ! ---- 7:MET SSD -----
        call histfill(7,metssd,mesq)

        ! ---- 11, 12: eta of leptons ---
        call histfill(11,eta1,mesq)
        call histfill(12,eta2,mesq)

     endif

     ! === Asymmetries ===

     ! Angle: L-R
     if(angle<0) asymnocut(1)=asymnocut(1)-mesq
     if(angle>0) asymnocut(1)=asymnocut(1)+mesq

     ! Angle: C-NC
     if(abs(angle)<0.5) then
        asymnocut(2)=asymnocut(2)+mesq
     else
        asymnocut(2)=asymnocut(2)-mesq
     end if

     ! Delta-Phi: L-R
     if(dphi<0) asymnocut(3)=asymnocut(3)-mesq
     if(dphi>0) asymnocut(3)=asymnocut(3)+mesq

     ! Delta-Phi: C-NC
     if(abs(dphi)<pi/2) then
        asymnocut(4)=asymnocut(4)+mesq
     else
        asymnocut(4)=asymnocut(4)-mesq
     end if

     ! 
     if(dabseta<0) asymnocut(5)=asymnocut(5)-mesq
     if(dabseta>0) asymnocut(5)=asymnocut(5)+mesq

     if(abs(dabseta)<1.25) then
        asymnocut(6)=asymnocut(6)+mesq
     else
        asymnocut(6)=asymnocut(6)-mesq
     end if

     ! === Cuts ===
     allowssd=.false.
     if (abs(eta1)<2.5 .and. abs(eta2)<2.5) then
        wcut(1)=wcut(1)+mesq
        if(pt(1)>40.0 .and. pt(2)>30.0) then
           wcut(2)=wcut(2)+mesq
           if(metssd>30.0) then
              wcut(3)=wcut(3)+mesq
              allowssd=.true.
           end if
        end if
     end if

     if(allowssd) then
        if(dohist) then
           ! === Fill histograms ===
           call histfill(3,(pt(1)-pt(2)),mesq)  ! Delta pT
           call histfill(4,(pt(1)+pt(2)),mesq)  ! Sum pT
           call histfill(5,pTasym,mesq)         ! pT-asymmetry
           call histfill(6,sl1l2,mesq)          ! Minv
           call histfill(8,(pt1+pt2+metssd),mesq) ! Sum pT + MET
           call histfill(9,angle,mesq)          ! cos(theta)
           call histfill(10,dphi,mesq)          ! Delta-phi
           call histfill(13,dabseta,mesq)       ! Delta-|eta|
        end if

        ! -- Asymmetries
        ! Angle: L-R
        if(angle<0) asym(1)=asym(1)-mesq
        if(angle>0) asym(1)=asym(1)+mesq

        ! Angle: C
        if(abs(angle)<0.5) then
           asym(2)=asym(2)+mesq
        else
           asym(2)=asym(2)-mesq
        end if


        ! Delta-Phi: L-R
        if(dphi<0) asym(3)=asym(3)-mesq
        if(dphi>0) asym(3)=asym(3)+mesq

        ! Delta-Phi: C-NC
        if(abs(dphi)<pi/2) then
           asym(4)=asym(4)+mesq
        else
           asym(4)=asym(4)-mesq
        end if

        ! Delta-|eta|: L-R
        if(dabseta<0) asym(5)=asym(5)-mesq
        if(dabseta>0) asym(5)=asym(5)+mesq

        ! Delta-|eta|: C
        if(abs(dabseta)<1.25) then
           asym(6)=asym(6)+mesq
        else
           asym(6)=asym(6)-mesq
        end if
     end if
  endif

  ! === LHEF file data ===
  ! Update maximum weight
  if(iSaveWgt) then
     if(mesq > wgtmax) wgtmax = mesq
     return
  end if

  ran = dble(rand(0))
  if(doLHEF .and. mesq/wgtmax > ran .and. ran>1.0d-10) then
     print*,"Adding event", mesq/wgtmax, ran
     write(997,*)'<event>'
     write(997,*)"9   1001   1.00E+0", scl, " 0.7816531E-02   1.219584e-01"
     if(isplus) then
        if(swapUD) then
           write(997,*)" -1   -1   0   0     0  101",(p(i,1),i=1,4),m(1,1), " 0.0  9"
           write(997,*)"  2   -1   0   0   101    0",(p(i,2),i=1,4),m(1,2), " 0.0  9"
        else
           write(997,*)"  2   -1   0   0   101    0",(p(i,1),i=1,4),m(1,1), " 0.0  9"
           write(997,*)" -1   -1   0   0     0  101",(p(i,2),i=1,4),m(1,2), " 0.0  9"
        end if
        write(997,*)"-13    1   1   2     0    0",(p(i,3),i=1,4),m(1,3), " 0.0  9"
        write(997,*)" 14    1   1   2     0    0",(p(i,4),i=1,4),m(1,4), " 0.0  9"
        write(997,*)" 25    2   1   2     0    0",(p(i,5),i=1,4),m(1,5), " 0.0  9"
        write(997,*)"-13    1   5   5     0    0",(p(i,8),i=1,4),m(1,8), " 0.0  9"
        write(997,*)" 14    1   5   5     0    0",(p(i,9),i=1,4),m(1,9), " 0.0  9"
        write(997,*)"  1    1   5   5   102    0",(p(i,10),i=1,4),m(1,10), " 0.0  9"
        write(997,*)" -2    1   5   5     0  102",(p(i,11),i=1,4),m(1,11), " 0.0  9"
     else
        if(swapUD) then
           write(997,*)"  1   -1   0   0   101    0",(p(i,1),i=1,4),m(1,1), " 0.0  9"
           write(997,*)" -2   -1   0   0     0  101",(p(i,2),i=1,4),m(1,2), " 0.0  9"
        else
           write(997,*)" -2   -1   0   0     0  101",(p(i,1),i=1,4),m(1,1), " 0.0  9"
           write(997,*)"  1   -1   0   0   101    0",(p(i,2),i=1,4),m(1,2), " 0.0  9"
        end if
        write(997,*)" 13    1   1   2     0    0",(p(i,3),i=1,4),m(1,3), " 0.0  9"
        write(997,*)"-14    1   1   2     0    0",(p(i,4),i=1,4),m(1,4), " 0.0  9"
        write(997,*)" 25    2   1   2     0    0",(p(i,5),i=1,4),m(1,5), " 0.0  9"
        write(997,*)" 13    1   5   5     0    0",(p(i,8),i=1,4),m(1,8), " 0.0  9"
        write(997,*)"-14    1   5   5     0    0",(p(i,9),i=1,4),m(1,9), " 0.0  9"
        write(997,*)" -1    1   5   5     0  102",(p(i,10),i=1,4),m(1,10), " 0.0  9"
        write(997,*)"  2    1   5   5   102    0",(p(i,11),i=1,4),m(1,11), " 0.0  9"
     endif
     write(997,*)'</event>'
  end if

end function mesq

!!$===================================================================

function ep(p1,p2,p3,p4)
  implicit none
  DOUBLE PRECISION::P1(4),P2(4),p3(4),p4(4),exp, ep

  exp = p1(1)*(p2(2)*(p3(3)*p4(4)-p3(4)*p4(3))  &
       &     - p2(3)*(p3(2)*p4(4)-p3(4)*p4(2)) &
       &     + p2(4)*(p3(2)*p4(3)-p3(3)*p4(2)))

  exp = exp - p1(2)*(p2(1)*(p3(3)*p4(4)-p3(4)*p4(3))  &
       &           - p2(3)*(p3(1)*p4(4)-p3(4)*p4(1)) &
       &           + p2(4)*(p3(1)*p4(3)-p3(3)*p4(1)))

  exp = exp + p1(3)*(p2(1)*(p3(2)*p4(4)-p3(4)*p4(2))  &
       &           - p2(2)*(p3(1)*p4(4)-p3(4)*p4(1)) &
       &           + p2(4)*(p3(1)*p4(2)-p3(2)*p4(1)))

  exp = exp - p1(4)*(p2(1)*(p3(2)*p4(3)-p3(3)*p4(2))  &
       &           - p2(2)*(p3(1)*p4(3)-p3(3)*p4(1)) &
       &           + p2(3)*(p3(1)*p4(2)-p3(2)*p4(1)))

  ep = exp

  return
end function ep

!!$===================================================================

function sphi(x,y)
  implicit none
  double precision::x,y,sphi,phi,pi
  pi = 4.0*atan(1.0)
  phi=atan(abs(y)/abs(x))
  if(x<0 .and. y<0) phi=phi+pi
  if(x<0 .and. y>0) phi=pi-phi
  if(x>0 .and. y<0) phi=2.0d0*pi-phi
  sphi=phi
  return
end function sphi

!!$===================================================================

subroutine inithist(ihist,hminin,hmaxin)

  implicit none
  double precision::hminin,hmaxin
  integer::ihist
  integer,parameter::nhist=20
  double precision::hist(nhist*100),hmin(nhist),hmax(nhist),totwgt(nhist)
  common/histo/hist,hmin,hmax,totwgt
  integer::iof(nhist),iuf(nhist),ndata(nhist)
  common/stat/iof,iuf,ndata

  hmin(ihist) = hminin
  hmax(ihist) = hmaxin
  totwgt(ihist)=0.0d0

  iof(ihist)=0
  iuf(ihist)=0
  ndata(ihist)=0
  hist(:)=0

  WRITE(999,*)'Initialising histogram',ihist,'with ',hmin(ihist),hmax(ihist)

end subroutine inithist

!!$-------------------------------------------------------------------

subroutine histfill(ihist,data,wgt)
  implicit none
  double precision::data,wgt
  integer::ihist
  integer,parameter::nhist=20
  double precision::hist(nhist*100),hmin(nhist),hmax(nhist),totwgt(nhist)
  common/histo/hist,hmin,hmax,totwgt
  integer::iof(nhist),iuf(nhist),ndata(nhist)
  common/stat/iof,iuf,ndata

  double precision::delta
  integer:: num

  delta = (hmax(ihist)-hmin(ihist))/100.0 
  ndata(ihist)=ndata(ihist)+1
  totwgt(ihist)=totwgt(ihist)+wgt

  if (data < hmax(ihist) .and. data >= hmin(ihist)) then
     num = (ihist-1)*100 + 1 + floor((data-hmin(ihist))/delta)
     hist(num) = hist(num) + wgt
  else if(data > hmax(ihist)) then
     iof(ihist)=iof(ihist)+1
     hist(ihist*100)=hist(ihist*100)+wgt
  else 
     iuf(ihist)=iuf(ihist)+1
     hist((ihist-1)*100+1)=hist((ihist-1)*100+1)+wgt
  end if

end subroutine histfill

!!$-------------------------------------------------------------------

subroutine printhist(unit,ihist)

  implicit none
  integer::ihist,i,unit
  integer,parameter::nhist=20
  double precision::hist(nhist*100),hmin(nhist),hmax(nhist),totwgt(nhist)
  common/histo/hist,hmin,hmax,totwgt
  integer::iof(nhist),iuf(nhist),ndata(nhist)
  common/stat/iof,iuf,ndata
  double precision::delta

  delta = (hmax(ihist)-hmin(ihist))/100.0 

  do i=1,100
     write(unit,*)hmin(ihist)+(i-1)*delta,hist(i+(ihist-1)*100)
  end do

  WRITE(999,*)'Histogram no:',ihist,hmin(ihist),hmax(ihist)
  WRITE(999,*)'Total data:',ndata(ihist)
  WRITE(999,*)'Overflow:',iof(ihist)
  WRITE(999,*)'Underflow:',iuf(ihist)

end subroutine printhist

!!$-------------------------------------------------------------------

subroutine histnorm(ihist,norm)

  !norm<0 => normalise to unity

  implicit none
  integer::ihist,i,ndata2
  double precision::norm,twgt
  integer,parameter::nhist=20
  double precision::hist(nhist*100),hmin(nhist),hmax(nhist),totwgt(nhist)
  common/histo/hist,hmin,hmax,totwgt
  integer::iof(nhist),iuf(nhist),ndata(nhist)
  common/stat/iof,iuf,ndata

  ndata2=ndata(ihist)-iof(ihist)-iuf(ihist)
  if(ndata2>0 .and. totwgt(ihist)>0) then
     if(norm>0.0) then
        do i=(ihist-1)*100+1,ihist*100
           hist(i)=hist(i)*norm/totwgt(ihist)
        enddo
     else
        twgt=0
        do i=(ihist-1)*100+1,ihist*100
           twgt = twgt + hist(i)
        enddo
        do i=(ihist-1)*100+1,ihist*100
           hist(i)=hist(i)/twgt
        enddo
     endif
  end if



end subroutine histnorm

SUBROUTINE GAUSSNO(X,NUM)

  implicit none
  DOUBLE PRECISION, INTENT(INOUT)::X(NUM)
  INTEGER :: NUM,i
  DOUBLE PRECISION::RAN,RAN1,RAN2

  DO I=1,NUM,2
     CHECK: DO
        RAN1 = 2*RAND(0)-1
        RAN2 = 2*RAND(0)-1
        RAN = RAN1**2+RAN2**2
        IF(RAN<1.0 .AND. RAN>0) EXIT CHECK
     END DO CHECK
     X(I) = RAN1*SQRT(-2D0*LOG(MAX(1D-10,RAN))/RAN)
     IF(NUM/=1) X(I+1) = RAN2*SQRT(-2D0*LOG(MAX(1D-10,RAN))/RAN)
  END DO

END SUBROUTINE GAUSSNO
