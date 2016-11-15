program bisect; implicit none
  real :: a0=-5.0, b0=15.0
  write(*,*) bisection(a0,b0)
  
contains
  function bisection(a, b) ! a < b 
    real a, b, fa, fb, fc, c, bisection!; intent(in) a, b
    real, parameter :: eps = 1.0e-22
    fa=f(a); fb = f(b)
    fc = f(c)
    do while (b-a > eps)
       if (fa*fb > 0) call abort
       c=(a+b)/2.0
       fc=f(c)
       if (fc == 0.0) exit
       if (SIGN(1.0,fa)*SIGN(1.0,fc) < 0.0) then
          b=c; fb=fc
       end if
       if (SIGN(1.0,fc)*SIGN(1.0,fb) < 0.0) then
          a=c; fa=fc
       end if
    end do
    bisection = c
  
  end function bisection
  
  elemental function f(x)
    real f, x; intent(in) x
    f = exp(x)*x - 5.0
  end function f
end program bisect
