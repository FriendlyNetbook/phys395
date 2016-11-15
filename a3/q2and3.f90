! Uses the Golden-Sections Method to find local minima of (x^2-1)^2+x
program q2and3
  real, parameter :: mrange = -10, prange = 10 ! Range to search for minima
  real, parameter :: dx = 0.00001 !h to use for derivative approximation for minimum-checking
  integer, parameter :: n=50 ! Number of blocks to check for minima
  real x, max, min, fx, fmax, fmin
  integer i
  write(*,*) "Local Minima of (x^2-1)^2+x"
  write(*,*) "                   x                      f(x)"
  do i=0,n-1
     min = mrange + i*(prange-mrange)/n
     max = mrange +(i+1.0)*(prange-mrange)/n
     fmax = f(max) !storing to save function evaluations as these are used in goldenSection() as fa, fb
     fmin = f(min)
     if ((f(min+dx) .LE. fmin) .and. (f(max-dx) .LE. fmax)) then !Checks that (min,max) brackets a minima to find
        x = goldenSection(min, max, fmax, fmin)
        write (*,*) x, f(x)
     end if
  end do
contains
  function goldenSection(a0, b0, fmax, fmin)
    real a0, b0, goldenSection, fmax, fmin
    real a, b, c,d, fa, fb, fc,fd
    real, parameter :: eps = 1.0e-16, phi = (1+sqrt(5.0))/2.0
    a=a0; b=b0; c=(b+a*phi)/(1+phi); d= (b+c*phi)/(1+phi)
    fa=fmax;fb=fmin;fc=f(c);fd=f(d)
    do while (abs(a-b) > eps)
       if (fc < fd) then !if minimum bracketed by (a,d)
          b=d; fb=fd
          d=c; fd=fc
          c=(b+a*phi)/(1+phi); fc=f(c) 
       else if (fc > fd) then !if minimum bracketed by (c,b)
          a=c; fa=fc
          c=d; fc=fd
          d=(b+c*phi)/(1+phi); fd=f(d)
       else
          exit !if changes in f are no longer distinguishable
       end if
    end do
   goldenSection = (c+d)/2
  end function goldenSection
  
  elemental function f(x)
    real f, x; intent(in) x
    f = x*(x*(x**2-2)+1)+1
  end function f
  
end program q2and3
