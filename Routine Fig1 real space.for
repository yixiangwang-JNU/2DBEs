*************************************************
*     Calculate the energy levels of s states
*************************************************
      module value
      implicit double precision (a-h,o-z)      
      real(8), parameter :: Pi=3.14159 
      real(8), parameter :: xme=0.29               ! electron mass, in unit of m0
      real(8), parameter :: xmh=0.64               ! hole mass  
      real(8), parameter :: xmr=0.2                ! reduced mass
      real(8), parameter :: r0=5                   ! screening length, in unit of nm
      real(8), parameter :: epsilonV=3.97          ! relative dielectric constant  
      real(8), parameter :: enV=2261.69/r0         ! energy scale of V   
      real(8), parameter :: deltar=0.05            ! lattice spacing
cc
      integer, parameter :: Ncut=1200              ! cutoff 
      integer, parameter :: Ndim=2*Ncut            ! dimension
      integer, parameter :: LDA=Ndim
      integer, parameter :: LDB=Ndim
      integer, parameter :: LDVL=Ndim
      integer, parameter :: LDVR=Ndim
      integer, parameter :: LWORK=8*Ndim   
      character(1), parameter :: JOBVL='N'
      character(1), parameter :: JOBVR='N'
cc
      real(8) Amat(Ndim,Ndim),Bmat(Ndim,Ndim),
     &        VL(LDVL,Ndim),VR(LDVR,Ndim),
     &        ALPHAR(Ndim),ALPHAI(Ndim),BETA(Ndim),WORK(LWORK),
     &        EVAL(Ndim)
      common l
cc      
      end module value
cc      
**********************
*     main program 
**********************
      program main
      use value
      implicit double precision (a-h,o-z)  
      character*100 filename
cc       
      l=0 
      write(filename,'(a3,I2,a4)') 'l=',l,'.dat'
      open(10,file=filename)  
cc 
      do 100 xB=-65,65.01,0.5                       ! magnetic field         
      write(*,*) xB  
      Amat=0
      do 201 i=1,Ncut-1                             ! Schrodinger equation
         r1=(i-1)*deltar
         r2=r1+deltar
         r=(r1+r2)/2
         Ci=-(0.5/r-1./deltar)*76.272/xmr           ! Kinetic energy
         Di=-(0.5/r+1./deltar)*76.272/xmr 
         Amat(i,Ncut+i)=Ci
         Amat(i,Ncut+i+1)=Di       
cc 
         F1i=38.136*l*l/(xmr*r*r) 
         F2i=2.198D-5*xB*xB*r*r/xmr 
         F3i=l*(1./xme-1./xmh)*0.05791*xB          
         call Keldysh(r,VKi)                        ! Keldsh potential
         Fi=F1i+F2i+F3i+VKi
         Amat(i,i)=Fi
         Amat(i,i+1)=Fi
201   continue
      Amat(Ncut,Ncut)=1                             ! boundary condition1  
      Amat(Ncut+1,Ncut+1)=1                         ! boundary condition2  
cc
      do i=1,Ncut-1                                 ! another equations   
         Amat(Ncut+i+1,i)=2
         Amat(Ncut+i+1,i+1)=-2
cc
         Amat(Ncut+i+1,Ncut+i)=deltar
         Amat(Ncut+i+1,Ncut+i+1)=deltar
      end do
cc
      Bmat=0
      do i=1,Ncut-1
         Bmat(i,i)=1
         Bmat(i,i+1)=1 
      end do
cc 
      call dggev(JOBVL,JOBVR,Ndim,Amat,LDA,Bmat,LDB,ALPHAR,ALPHAI,BETA,
     &           VL,LDVL,VR,LDVR,WORK,LWORK,INFO) 
cc
cc    sort the (Ncut-1) eigenvalues         
      do i=1,Ncut-1 
         EVAL(i)=ALPHAR(i)/BETA(i)    
      end do
cc      
      do i=1,Ncut-2
      do j=i+1,Ncut-1
         if(EVAL(i)>EVAL(j)) then
            temp=EVAL(i)
            EVAL(i)=EVAL(j)
            EVAL(j)=temp       
         end if 
      end do
      end do
      write(10,'(100f16.8)') xB,(EVAL(i),i=1,4)
cc      
100   continue
cc 
      end program main     
cc   
***************************************
*     subroutine: Keldysh potential
***************************************  
      subroutine Keldysh(r,VK)
      use value
      implicit double precision (a-h,o-z) 
cc
      call Stvh0(epsilonV*r/r0,sh0)
      call Jy01a(epsilonV*r/r0,by0)  
      VK=-enV*(sh0-by0)    
cc    
      end subroutine Keldysh
cc          
************************************
*     subroutine: Stuve function
************************************
      subroutine Stvh0(x,sh0)
!     Input: x, the argument 
!     Output: sh0, the value of H0(x) 
!
      implicit double precision (a-h,o-z)
      pi=3.141592653589793D+00
      s=1.0D+00
      r=1.0D+00
cc
      if(x<=20.0D+00) then
      a0=2.0D+00*x/pi
      do k=1,60
      r=-r*x/(2.0D+00*k+1.0D+00)*x/(2.0D+00*k+1.0D+00)
      s=s+r
      if(abs(r)<abs(s)*1.0D-12) then
        exit
      end if
      end do
cc
      sh0=a0*s
cc
      else 
cc      
      if(x<50.0D+00) then
      km=int(0.5D+00*(x+1.0D+00))
      else
      km=25
      end if
cc
      do k=1,km
      r=-r*((2.0D+00*k-1.0D+00)/x)**2
      s=s+r
      if(abs(r)<abs(s)*1.0D-12) then
        exit
      end if
      end do
cc
      t=4.0D+00/x
      t2=t*t
cc
      p0=((((- 0.37043D-05*t2   
     &    +0.173565D-04)*t2 
     &    -0.487613D-04 )*t2 
     &    +0.17343D-03 )*t2 
     &    -0.1753062D-02)*t2 
     &    +0.3989422793D+00
cc
      q0=t*(((((0.32312D-05*t2      
     &    -0.142078D-04)*t2
     &    +0.342468D-04)*t2 
     &    -0.869791D-04)*t2 
     &    +0.4564324D-03)*t2 
     &    -0.0124669441D+00)
cc
      ta0=x-0.25D+00*pi
      by0=2.0D+00/sqrt(x)*(p0*sin(ta0)+q0*cos(ta0))
      sh0=2.0D+00/(pi*x)*s+by0
cc 
      end if
cc
      return
      end subroutine Stvh0      
cc
********************************************************
*     subroutine: Bessel function of the second kind
********************************************************
      subroutine Jy01a(x,by0)
!     Input: x, the argument
!     Output: BY0, the value of Y0(x)
!
      implicit double precision (a-h,o-z)
      real(8) :: a(12)=(/ 
     & -0.7031250000000000D-01, 0.1121520996093750D+00, 
     & -0.5725014209747314D+00, 0.6074042001273483D+01, 
     & -0.1100171402692467D+03, 0.3038090510922384D+04, 
     & -0.1188384262567832D+06, 0.6252951493434797D+07, 
     & -0.4259392165047669D+09, 0.3646840080706556D+11, 
     & -0.3833534661393944D+13, 0.4854014686852901D+15 /)
      real(8) :: b(12)=(/ 
     & 0.7324218750000000D-01, -0.2271080017089844D+00, 
     & 0.1727727502584457D+01, -0.2438052969955606D+02, 
     & 0.5513358961220206D+03, -0.1825775547429318D+05, 
     & 0.8328593040162893D+06, -0.5006958953198893D+08, 
     & 0.3836255180230433D+10, -0.3649010818849833D+12, 
     & 0.4218971570284096D+14, -0.5827244631566907D+16 /)
cc
      pi=3.141592653589793D+00
      rp2=0.63661977236758D+00
      x2=x*x
      if (x==0.0D+00) then
      by0=-1.0D+300
      return
      end if
cc
      if(x<=12.0D+00) then
cc
      bj0=1.0D+00
      r=1.0D+00
      do k=1,30
      r=-0.25D+00*r*x2/(k*k)
      bj0=bj0+r
      if(abs(r)<abs(bj0)*1.0D-15) then
        exit
      end if
      end do
cc
      ec=log(x/2.0D+00)+0.5772156649015329D+00
      cs0=0.0D+00
      w0=0.0D+00
      r0=1.0D+00
      do k=1,30
      w0=w0+1.0D+00/k
      r0=-0.25D+00*r0/(k*k)*x2
      r=r0*w0
      cs0=cs0+r
      if(abs(r)<abs(cs0)*1.0D-15) then
        exit
      end if
      end do
cc 
      by0=rp2*(ec*bj0-cs0)
      cs1=1.0D+00
      w1=0.0D+00
      r1=1.0D+00
      do k=1,30
      w1=w1+1.0D+00/k
      r1=-0.25D+00*r1/(k*(k+1))*x2
      r=r1*(2.0D+00*w1+1.0D+00/(k+1.0D+00))
      cs1=cs1+r
      if(abs(r)<abs(cs1)*1.0D-15) then
        exit
      end if
      end do
cc
      else
cc
      if(x<35.0D+00) then
      k0=12
      else if(x<50.0D+00) then
      k0=10
      else
      k0=8
      end if
cc
      t1=x-0.25D+00*pi
      p0=1.0D+00
      q0=-0.125D+00/x
      do k=1,k0
      p0=p0+a(k)*x**(-2*k)
      q0=q0+b(k)*x**(-2*k-1)
      end do
      cu=sqrt(rp2/x)
      by0=cu*(p0*sin(t1)+q0*cos(t1))
cc
      end if
cc 
      return
      end subroutine Jy01a
