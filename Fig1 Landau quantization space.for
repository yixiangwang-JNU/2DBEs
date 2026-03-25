*****************************************************************
*     Calculate the exciton energy levels within the LL basis
*****************************************************************
      module value
      implicit double precision (a-h,o-z)      
      real(8), parameter :: Pi=3.14159 
      real(8), parameter :: xme=0.29             ! electron mass, in unit of m0
      real(8), parameter :: xmh=0.64             ! hole mass  
      real(8), parameter :: xmr=0.2              ! reduced mass
      real(8), parameter :: r0=5                 ! screening length, in unit of nm
      real(8), parameter :: epsilonV=3.97        ! relative dielectric constant  
      real(8), parameter :: Alower=0             ! lower limit of the integral       
cc
      integer, parameter :: l=0                  ! angular quantum number
      integer, parameter :: Ncut=300             ! cutoff 
      integer, parameter :: Ndim=Ncut            ! dimension
      integer, parameter :: LDA=Ndim
      integer, parameter :: LWORK=3*Ndim   
      character(1), parameter :: JOBZ='N'
      character(1), parameter :: UPLO='U'  
cc
      real(8) H(Ndim,Ndim),EVAL(Ndim),WORK(LWORK)   
      common xB,xlB,nmin,nmax,nabs    
cc
      end module value
cc      
**********************
*     main program
**********************
      program main
      use value
      use QDAG_INT  
      include 'link_fnl_shared.h'       
      implicit double precision (a-h,o-z)
      external V1
      character(100) filename
      write(filename,'(a3,I2,a4)') 'l=',l,'.dat'
      open(10,file='10.dat')  
      open(20,file=filename)
cc
      do 100 xB=1,65,1       ! magnetic field
      write(*,*) xB
      xlB=25.6/sqrt(xB)      ! magnetic length 
      omegae=0.11578*xB/xme
      omegah=0.11578*xB/xmh
cc
      H=0
      do n=1,Ndim            ! n is the hole index that begins from 0  
         nmin=n-1            ! note the difference between the LL index and matrix index  
         nmax=n-1
         nabs=0
cc
         Ten=omegae*(n-1+l+0.5)+omegah*(n-1+0.5) 
cc       Bupper=nmax*0.6+20                      ! set the upper limit 
         Bupper=nmax*0.2+24       
         call QDAG(V1,Alower,Bupper,Res) 
         Ven=-Res*1439.83                        ! potential energy  
cc 
         H(n,n)=Ten+Ven  
cc       write(10,'(I4,100f16.8)') n,Ten,Ven 
      end do
cc
      do n1=1,Ndim-1
      do n2=n1+1,Ndim
         nmin=n1-1                           ! note the difference between the LL index and matrix index 
         nmax=n2-1
         nabs=n2-n1     
cc
cc       Bupper=nmax*0.6+20                  ! set the upper limit 
         Bupper=nmax*0.2+24         
         call QDAG(V1,Alower,Bupper,Res) 
         Ven=-Res*1439.83                    ! potential energy  
cc
         H(n1,n2)=Ven
         H(n2,n1)=Ven
cc       write(10,'(2I4,f16.8)') n1,n2,Ven        
      end do
      end do
cc       
      call DSYEV(JOBZ,UPLO,Ndim,H,LDA,EVAL,WORK,LWORK,INFO) 
      write(20,'(100f16.8)') xB,(EVAL(i),i=1,5) 
cc      
100   continue      
cc 
      end program main 
cc       
********************************************
*     integrand function of the integral
********************************************
      function V1(x)
      use value
      implicit double precision (a-h,o-z)
cc
      x1=x*x/2
      call TLag_nm(nmin+l,nabs,x1,Res1) 
      call TLag_nm(nmin,nabs,x1,Res2) 
      V1=Res1*Res2/(xlB*epsilonV+x*r0) 
cc       
      end function V1
cc
***************************
*     calculate TLnm(x)
*************************** 
      subroutine TLag_nm(n,m,x,Res)
      implicit double precision (a-h,o-z)
cc   
      Resa=exp(-x/2.)                       ! n=0, TLag_0m
      do i=1,m
         Resa=Resa*sqrt(x)/sqrt(real(i))
      end do   
cc 
      Resb=(1+m-x)*exp(-x/2.)               ! n=1, TLag_1m 
      do i=1,m
         Resb=Resb*sqrt(x)/sqrt(real(i+1))
      end do        
cc
      if(n==0) Res=Resa 
      if(n==1) Res=Resb 
cc      
      if(n>1) then
      do i=2,n 
         Res=Resb*(2*i+m-x-1)/sqrt(real(i*(i+m)))
     &       -Resa*sqrt(real((i-1)*(i+m-1)))/sqrt(real(i*(i+m)))   
         Resa=Resb
         Resb=Res 
      end do   
      end if 
cc
      end subroutine TLag_nm
cc