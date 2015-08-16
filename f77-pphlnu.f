      program pphlnu

      implicit none

      integer nev,iev,niter,ndim
      double precision mh,mw,swsq,mz,s,k,pdgcross,pi,alpha,e
      double precision crossx,errorx,mesq,chi2,gev2pb
      integer i,j

      double precision p(5,30),m(2,30)
      common /partcl/p,m

      common/input/mw,mz,swsq,mh,s

      character*20 parm(20)
      double precision value(20)

      external mesq

! -------- input params ----------
      nev = 1000000
      niter=20
      iev=0
      mw=80.0
      mz=91.0
      swsq=0.23
      pi=datan(1.0d0)*4.0d0
      s=1.4d4**2
      mh=150.0d0

      do i = 1,30
         do j=1,5
            p(j,i)=0.0d0
            m(j,i)=0.0d0
         enddo
      enddo

      m(1,6) = mh
      m(2,6) = mh**2

      ndim=9

! -------------Initialise CERNLIB/LHAPDF -------------- 

      parm(1)='NPTYPE'
      value(1)=1
      parm(2)='NGROUP'
      value(2)=4
      parm(3)='NSET'
      value(3)=46

      call pdfset(parm,value)

! ------------- Integrate ------------------------
      call montes(mesq,ndim,40,nev,niter,crossx,errorx,
     & .false.,0,.false.,.false.,chi2)

      gev2pb = 3.89d8

      crossx = crossx * gev2pb 
      errorx = errorx * gev2pb 

! ----------- Results and statistics -------------
      call pdfsta
      
      print*,'---------------------------------'
      print*,'The cross section is',crossx,'+/-', errorx
      print*,'higgs mass',mh,'sqrt(s)',sqrt(s)

      end program pphlnu

      function mesq(rf,ndim)

      implicit none
      integer ndim
      double precision rf(ndim)
      double precision mesq

      double precision p(5,30),m(2,30)
      common /partcl/p,m

      double precision s1,s2,s3,zz
      integer iflag
      common /dk3cmn/s1,s2,s3,zz(5),iflag                                      

      double precision mw,mz,swsq,mh,s
      common/input/mw,mz,swsq,mh,s

      integer i
      double precision pi, wps3, g, flux, norm, x, cosqw, const
      double precision p1_p3,p2_p4,p1_p4,p2_p3,p1_p2,p3_p4, denom, gamw
      double precision  x1,x2,scl,e,pdf1,pdf2
      double precision  u1,d1,ub1,db1,sb1,ch1,bt1,tp1,gl1
      double precision  u2,d2,ub2,db2,sb2,ch2,bt2,tp2,gl2

      pi=atan(1.0)*4.0d0

      e=dsqrt(s)/2.0d0
      x1=rf(1)
      x2=rf(2)

!     ---- check that on-shell higgs is possible, otherwise dk3 gives error    
      If(x1*x2*S.lt.mh**2) then
         mesq=0.0
         return
      endif

      p(1,1)=0.0d0
      p(2,1)=0.0d0
      p(3,1)=e*x1
      p(4,1)=e*x1
      p(5,1)=p(4,1)

      p(1,2)=0.0d0
      p(2,2)=0.0d0
      p(3,2)=-e*x2
      p(4,2)=e*x2
      p(5,2)=p(4,2)
      
      do i=1,5
         p(i,3)=p(i,1)+p(i,2)
      enddo

      m(2,3)=s*x1*x2
      m(1,3)=dsqrt(m(2,3))
      scl = m(2,3)

      call structm(x1,scl,u1,d1,ub1,db1,sb1,ch1,bt1,tp1,gl1)
      call structm(x2,scl,u2,d2,ub2,db2,sb2,ch2,bt2,tp2,gl2)

      ! u(ubar) <=> x1
      pdf1 = (u1*db2+ch1*sb2+ub1*d2+ch1*sb2) /(x1*x2) 
      ! u(ubar) <=> x2
      pdf2 = (u2*db1+ch2*sb1+ub2*d1+ch2*sb1) /(x1*x2) 

      do i=1,5
         zz(i)=rf(i)
      enddo

      call dk3(3,4,5,6,wps3)

      if(iflag/=0) then
         mesq=0.0d0
         print*,'Point discarded'
         return
      end if

      call dot_4(1,4,p1_p3)
      call dot_4(2,5,p2_p4)
      call dot_4(1,5,p1_p4)
      call dot_4(2,4,p2_p3)
      call dot_4(1,2,p1_p2)
      call dot_4(4,5,p3_p4)


      mesq = 64.0d0 * (p1_p4*p2_p3*pdf1 + p1_p3*p2_p4*pdf2)

      g = 0.65d0

      gamw = 2.141d0*0.107d0    ! From PDG for leptonic decay only

      const = g**6 * mw**2 / 16.0d0

      denom = ((2.0*p1_p2-mw**2)**2 + (mw*gamw)**2) 
     & * ((2.0*p3_p4-mw**2)**2 + (mw*gamw)**2)

      flux = 0.50d0/m(2,3)      ! 1/(2 s^hat)

!     print*,mesq*wps3,denom,pdf

      norm = 1.0d0/(128.0d0 * pi**3) ! = 1/(2pi)^6 * pi^3/2
      mesq = const * mesq * wps3 * flux * norm / denom

      return

      end function mesq

  
