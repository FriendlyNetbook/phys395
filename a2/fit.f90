program expansion; implicit none
integer, parameter :: m = 3001, bases = 4 !m is the resolution of the fit evaluation, bases is the number of coefficients to evaluate - correspods to "n"+1 from the assignment
real, parameter :: pi = 3.141592653589793238462643383279502884197169399375Q0
integer :: allocations = 1, status, n
integer i
real, allocatable :: x(:), y(:), xT(:), yT(:), B(:,:)

call readInAndAllocate !Read in the file "Assignment #2.dat", and allocate B

do i=1,n !Make basis array for 'bases' basis functions evaluated at n points
   B(i,:)=basis(bases,x(i))
end do

!solve the matrix Bx=y, storing x in y(1:bases)
call msolve(n, B, y, 0.0) !Do SVD to find fit
call dump(m) !return the fit results

contains
  subroutine readInAndAllocate !Read in "Assignment #2.dat" in 3000-line chunks into x, y
    allocate(x(3000),y(3000))
    i=0
    do
       i=i+1
       read (*,*,iostat=status) x(i), y(i) !read a line
       if (status < 0) exit !check that line is not EOF, otherwise end loop
       if(mod(i,3000) == 0) then !if that was the end of a chunk
          allocations = allocations+1 
          allocate(xT (i), yT (i)) !copy of the contents of x,y
          xT = x
          yT = y
          deallocate(x,y)
          allocate(x (allocations*3000), y (allocations*3000)) !extend x,y another 3000
          x(:)=xT !copy back contents of x,y
          y(:)=yT
          deallocate(xT, yT)
       end if
    end do
    n=i-1 !set n based on # of iterations
    allocate(yT (n)) !copy y for calulation of Chi^2 later
    yT=y
    allocate(B(n,bases))
  end subroutine readInAndAllocate
  
  ! function to be approximated
  elemental function f(x)
    real f, x; intent(in) x
    f = 1.0/(1.0 + 10.0*x*x)
  end function f
  
  ! evaluate basis functions b_i(x) at point x
  pure function basis(t, x)
    integer k, t; real x, basis(t); intent(in) t, x
    forall (k=0:t-1) basis(k+1) = cos(k*2*pi*x)
  end function basis
  
  ! uniformly spaced grid of n points covering interval [min,max]
  pure function uniform(n, min, max)
    integer i, n; real uniform(n), min, max; intent(in) n, min, max
    forall (i=1:n) uniform(i) = min + i*abs(min-max)/(n*1.0)
  end function uniform

  !Minimize |Mx-d| using LAPACK SVD
  subroutine msolve(n, M, d, rcond)
    integer n; real M(n,bases), d(bases), S(n), W(6*n), rcond
    integer rank, status
    status = 0
    select case (kind(M))
    case(4); call sgelss(n, bases, 1, M, n, d, n, S, rcond, rank, W, 6*n, status)
    case(8); call dgelss(n, bases, 1, M, n, d, n, S, rcond, rank, W, 6*n, status)
    case default; call abort
    end select
    ! bail at first sign of trouble
    if (status /= 0) call abort
  end subroutine msolve
  
  subroutine dump(m)
    integer i, m; real xEval(m), chi2
    xEval = uniform(m, minval(x), maxval(x)) !Make xEval on the same range as input data x
    open(unit = 1, file="FitEvaluation")
    do i = 1,m !Write the fit evaluation to a file for plotting
       write(1,*) xEval(i), sum(y(1:bases)*basis(bases,xEval(i)))
    end do
    close(1)
    do i = 1,n !Calculate Chi^2
       chi2 = chi2 + (yT(i)-sum(y(1:bases)*basis(bases,x(i))))**2
    end do
    chi2 = chi2 /(n*0.5)
    write(*,*) "Coefficients: ", y(1:bases)
    write(*,*) "Chi^2 = ", chi2 
  end subroutine dump
end program
