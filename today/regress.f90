program regress; implicit none

integer, parameter :: n = 1500

real ::  x(n), y(n), a, b, sx=0.0,sy=0.0,sxy=0.0,sxx=0.0, sxsy=0.0, sxxsy=0.0, sxsxy=0.0

integer i,j, status

! read data from stdin, at most n points
do i = 1,n
   read (*,*,iostat=status) x(i), y(i)
   if (status < 0) exit
end do


do i=1,n
   sx = sx+x(i)
   sy = sy+y(i)
   sxx = sxx+x(i)*x(i)
   sxy = sxy+x(i)*y(i)
end do

do i=1,n
   sxsy = sxsy + x(i)*sy
   sxxsy = sxxsy + x(i)*x(i)*sy
   sxsxy = sxsxy + x(i)*sxy
end do


a=(n*sxy - sxsy)/(n*sxx-sx*sx)
b=-1.0*(sxsxy - sxxsy)/(n*sxx-sx*sx)

write(*,*) a,b

end program regress
