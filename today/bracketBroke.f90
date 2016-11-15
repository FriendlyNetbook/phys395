program bracketing; implicit none

real x

x = bracket(1.0, 2.0)

!write (*,*) x

contains

function bracket(a0, c0)
  real a0, c0, bracket
  real a, b, c, fa, fb, fc, ma, mc
  real, parameter :: eps = 1.0e-16
   
  a = a0; fa = f(a)
  b = (a0 + c0)/2.0; fb = f(b) 
  c = c0; fc = f(c)
  
  !if ((fb-fa)*(fc-fb) > 0) call abort
  !if (fb >= fa .or. fb >= fc) call abort
  
  do while (abs(c-a) > eps)
     write(*,*) a,b,c, 'f', fa, fb,fc
     ! test if ab or bc bracket the minimum
     ma=f(((a+b)/2))
     mc=f(((c+b)/2))
     if (ma < fb .AND. ma < fa) then
        c = b; b=(a+b)/2; fc=fb; fb=ma
     end if
     
     if (mc < fc .AND. mc < fb) then
        a = b; b=(b+c)/2; fa=fb; fb=mc
     end if
     
  end do
        
  bracket = b
end function bracket

elemental function f(x)
  real f, x; intent(in) x
  
  f = (x*x - 1.0)**2 + x
end function f

end program bracketing
    
