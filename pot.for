      program pot
      
      implicit real*8(a-h,o-z)
      dimension vn(1000), vc(1000)
      common/ppot/ v0,a,r0,rc0,nz,na
      
      open(unit = 1, file = 'pot_huy.out')
      open(unit = 2, file = 'pot.out')
      open( unit = 3, file = 'cpot.out')
      print*, 'Which type of potential do you choose?'
      print*, '[1=squared; 2=harmonic; 3=Gaussian; 4=Woods-Xason]'
      read*, kpot
      print*, 'The largest distance of rmax ='
      read*, rmax
      print*, 'The step of points h ='
      read*, h
      print*, 'The number of proton Nz='
      read*, nz
      print*, 'The number of mass A ='
      read*, na
      print*, 'The nuclear radius of RC0 = '
      read*, rc0
!Number of points
      n = int(rmax/h)              
      
      write(2,*) 'NUCLEAR CENTRAL POTENTIALS'
      write(2,*) ' '
!Squared potential
      if (kpot.eq.1) then
       print*, 'The squared potential V(x) = -V0 with 0 <= x <= a'
       print*, 'The depth of potential V0 ='
       read*, v0
       print*, 'The width of potential a ='
       read*, a
       print*, 'r (fm)               V(MeV)'
       do i = 1, n
        r = i*h
        v = vsquare(r)
        vn(i) = v
        print*, r, v
        write(1,*) r, v
       end do
       call simpson(vsquare,0.d0,rmax,tp,n)
       write(1,*) 'Potential volume =', tp
      end if

!Parabolic potential
      if (kpot.eq.2) then
       print*, 'The parabolic potential V(x) = -V0 + ax^2 < 0'
       print*, 'The depth of potential V0 ='
       read*, v0
       print*, 'a ='
       read*, a
       print*, 'r (fm)               V(MeV)'
       do i = 1, n
         r = i*h
         v = vpara(r)
         vn(i) = v
       print*, r, v
       write(1,*) r, v
       end do
       call simpson(vpara,0.d0,rmax,tp,n)
       write(1,*) 'Potential volume =', tp
      end if
      
!Gaussian potential
      if(kpot.eq.3) then
       print*, 'The Gaussian potential V(x) = -V0*exp(-ax^2) <0'
       print*, 'The depth of potential V0='
       read*, v0
       print*, 'a='
       read*, a
       print*, 'r (fm)               V(MeV)'
       do i = 1, n
        r = i*h
        v = vgauss(r)
        vn(i) = v
        print*, r, v
        write(1,*) r, v
       end do
       call simpson(vgauss,0.d0,rmax,tp,n)
       write(1,*) 'Potential volume =', tp
      end if
      
!Woods-saxon potential
      if(kpot.eq.4) then
       print*, 'The Woods-saxon potential V(x) = -V0/(1+exp((x-x0)/a)) '
       print*, 'The depth of potential V0='
       read*, v0
       print*, 'a='
       read*, a
       print*, 'x0='
       read*, r0
       print*, 'r (fm)               V(MeV)'
       do i = 1, n
        r = i*h
        v = vws(r)
        vn(i) = v
        print*, r, v
        write(1,*) r, v
       end do
       call simpson(vws,0.d0,rmax,tp,n)
       write(1,*) 'Potential volume =', tp
      end if
      
      write(2,'(10E14.7)')(vn(i),i=1,n)
      
!Compute coulomb

      print*, ' r(fm)       V(MeV)'
       do i = 1, n
        r = i*h
        v1 = vcou(r)
        vc(i) = v1
        print*, r, vcou(r)
       write(3,*) r,v1
      end do
      call simpson(vcou,0.d0,rmax,tp,n)
      write(3,*) 'Potential volume =' , tp
      
!Write Coulomb
      write(2,*) 'COULOMB POTENTIALS'
      write(2,*) ' '
      write(2,'(10E14.7)')(vc(i),i=1,n)

      end program

!Squared potential
      function vsquare(x)
      implicit real*8(a-h,o-z)
      common/ppot/ v0,a,r0,rc0,nz,na
      if (x.ge.0.and.x.le.a.or.x.eq.0.or.x.eq.a) then
       vsquare = -v0
      else 
       vsquare = 0d0
      end if
      end function 

!Parabolic potential
      function vpara(x)
      implicit real*8(a-h,o-z)
      common/ppot/ v0,a,r0,rc0,nz,na
      xmax = sqrt(v0/a)
      if (x.ge.0.and.x.le.xmax.or.x.eq.0.or.x.eq.xmax) then
       vpara = -v0 + a*x**2
      else 
       vpara = 0d0
      end if
      end function 

!Gaussian potential  
      function vgauss(x)
      implicit real*8 (a-h,o-z)
      common/ppot/ v0,a,r0,rc0,nz,na
      vgauss = -v0*exp(-a*x**2)
      end function
      
!Woods-saxon potential
      function vws(x)
      implicit real*8 (a-h,o-z)
      common/ppot/ v0,a,r0,rc0,nz,na
      vws = -v0/((1+exp((x-r0)/a)))
      end function
      
!Coulomb potential    
      function  vcou(r)
      implicit real*8 (a-h,o-z)
      common/ppot/ v0,a,r0,rc0,nz,na
      rc = rc0*dfloat(na)**(1.d0/3.d0)
      if (r.ge.rc)  then
       vcou = (dfloat(nz)*1.44)/r
      else
       vcou = ((dfloat(nz)*1.44)/(2*rc))*(3-r**2/rc**2)
      end if
      end function
      
      
       Subroutine simpson(f,a,b,tp,n)
!==========================================================
! Integration of f(x) on [a,b]
! Method: Simpson rule for n intervals  
! written by: Alex Godunov (October 2009)
!----------------------------------------------------------
! IN:
! f   - Function to integrate (supplied by a user)
! a	  - Lower limit of integration
! b	  - Upper limit of integration
! n   - number of intervals
! OUT:
! integral - Result of integration
!==========================================================
      implicit none
      double precision f, a, b, tp,s
      double precision h, x
      integer n, i

! if n is odd we add +1 to make it even
      if((n/2)*2.ne.n) n=n+1

! loop over n (number of intervals)
       s = 0.0
       h = (b-a)/dfloat(n)
       do i=2, n-2, 2
        x   = a+dfloat(i)*h
        s = s + 2.0*f(x) + 4.0*f(x+h)
       end do
      tp = (s + f(a) + f(b) + 4.0*f(a+h))*h/3.0
      return
      end subroutine simpson
      
