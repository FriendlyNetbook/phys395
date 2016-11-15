! Finds the roots of x^3-x+1/4 to machine precision
program q1
  integer, parameter :: nMax = 100
  real x
  integer i
  write(*,*) "Roots of x^3 - x + 1/4"
  write(*,*) newton(-1.0,nMax)
  write(*,*) newton(0.0,nMax)
  write(*,*) newton(1.0,nMax)

contains
  function newton(a,maxN)
    real a, x, dx, newton
    integer i, maxN
    x = a
    do i = 1,maxN
       dx=-1*f(x)/df(x) !Iterate x
       x = x + dx
       if (dx == 0.0) then !If x has stopped changing, dx must be below machine epsilon 
          exit
       end if
    end do
    newton = x
  end function newton
  
  elemental function f(x) !function
    real f, x; intent(in) x
    f = x**3-x+0.25
  end function f
  
  elemental function df(x) !function derivative
    real df, x; intent(in) x
    df = 3.0*(x**2) - 1.0
  end function df
end program q1
