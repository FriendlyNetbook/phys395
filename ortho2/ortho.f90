program ortho; implicit none
  integer, parameter :: n = 150, k=2
  real, parameter :: pi= 3.1415925435Q0
  real x(n), y(n)
  integer i

  do i = 1,n
     x(i) = (2*i-1.0-n)/(n-1.0)
     y(i) = f(x(i))

     write (*,*) x(i), y(i), sin(2*pi*k*x(i)), ChebyshevTk(k, x(i)), LegendreP(k,x(i))
  end do
  
contains

  pure function f(x)
    real f,x; intent(in) x
    f = 1.0/(1.0+10.0*x*x)
  end function f

  pure function ChebyshevTk(k, x)
    real ChebyshevTk, x; integer k; intent(in) k, x
    ChebyshevTk = cos(k*acos(x))
  end function ChebyshevTk

  pure function LegendreP(k, x)
    real LegendreP, x; integer k; intent(in) k, x
    real P(0:k); integer n
    
    P(0) = 1.0; if (k==0) then; LegendreP = P(0); return; end if
    P(1) = x; if (k==1) then; LegendreP = P(1); return; end if
    
    do n = 1,k-1
       P(n+1) = ((2*n+1)*x*P(n) - n*P(n-1)/(n+1))
    end do

    LegendreP = P(k)
  end function LegendreP
  
end program ortho
