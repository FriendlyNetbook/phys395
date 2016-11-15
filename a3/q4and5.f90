! Computes the minima of (x^2+y^2-1)^2+mu x (mu = 0.01)
! using the undamped Levenberg-Marquardt method
program q4and5; implicit none
  real, parameter :: mu = 0.01
  real, parameter :: lambda = 1.0
  real, parameter :: eps = 1.0e-12
  real fx, x(2), grad(2), dx(2)
  integer i
  
  x = [0.0,3.0]
  write(*,*) 0, x 
  do i = 1,100 !max 100 iterations
     grad = df(x)
     dx = lambda * matmul(inverse(ddf(x)), grad)
     x=x-dx
     write (*,*) i, x
     if (maxval(abs(dx)) .LE. eps) then
        exit
     end if
  end do
  open(unit = 1, file="results.txt")
  write(1,*) "Question 5"
  write(1,*) "Using undamped Levenberg-Marquardt,"
  write(1,*) "Number of Iterations for accuracy of 1.0e-12: ", i
  close(1)
  
  
contains
  function f(x)
    real x(2), f
    associate ( x => x(1), y => x(2) )
      f = x*(x*((x**2)+y*(y*2.0)-2.0)+mu)+(y**2)*(y**2-2.0)+1.0
    end associate
  end function f
  
  function df(x)
    real x(2), df(2)
    associate ( x => x(1), y => x(2) )
      df(1) = 4.0*(x*x + y*y - 1.0) * x + mu
      df(2) = 4.0*(x*x + y*y - 1.0) * y
    end associate
  end function df
  
  function ddf(x)
    real x(2), ddf(2,2)
    associate ( x => x(1), y => x(2) )
      ddf(1,1) = 12.0*x*x + 4.0*y*y - 4.0
      ddf(2,1) = 8.0*x*y
      ddf(1,2) = 8.0*x*y
      ddf(2,2) = 12.0*y*y + 4.0*x*x - 4.0
    end associate
  end function ddf
  
  function inverse(M)
    real M(2,2), inverse(2,2), det
    det = M(1,1)*M(2,2) - M(1,2)*M(2,1)
    inverse(1,1) =  M(2,2)/det
    inverse(1,2) = -M(2,1)/det
    inverse(2,1) = -M(1,2)/det
    inverse(2,2) =  M(1,1)/det
  end function inverse
end program q4and5
