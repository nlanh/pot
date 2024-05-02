      program pot
      
      implicit real*8(a-h,o-z)
      dimension vn(1000)
      
      open(unit = 1, file = 'potin.dat')
      open(unit = 2, file = 'potout.dat')      
      read(1,*)h,e0,rmax
      n = int((rmax-e0)/h+1e-8)              

       do i = 1, n
        write(1,*) r, vn(i)
       end do
       write(2,'(10E14.7)')(vn(i),i=1,n)

      end program
