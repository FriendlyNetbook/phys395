!Part 4: Calculate Error of Chebyshev Interpolant and Interpolant Derivative
program expansion; implicit none

integer, parameter :: n = 100, m = 3001
real, parameter :: pi = 3.141592653589793238462643383279502884197169399375Q0

real x(n)	! function is evaluated at these (grid) points
real v(n)	! values of the function f(x) at these points
real B(n,n)	! B(i,j) contains value of j-th basis function at point i

integer i

write(*,*) 'N =', n
write(*,*) 'Uniform Grid'
write(*,*) '------------'
! initialize grid and evaluate function
x = uniform(n); v = f(x)

! evaluate the basis on the grid points
do i = 1,n; B(i,:) = basis(n, x(i)); end do

! find coefficients of polynomial expansion
call lsolve(n, B, v)

! dump polynomial expansion at a finer grid
call dump(m)

contains

! function to be approximated
elemental function f(x)
	real f, x; intent(in) x
	
	f = 1.0/(1.0 + 10.0*x*x)
end function f

elemental function fP(x)
  real fP, x; intent(in) x
  fP = -20.0*x/((10*x*x+1)*(10*x*x+1))
end function fP

! evaluate basis Derivative functions b_i(x) at point x
pure function basisPrime(n, x)
	integer k, n; real x, basisPrime(n); intent(in) n, x
	forall (k=1:n) basisPrime(k) = (k-1.0)*sin((k-1.0)*acos(x))/sqrt(1.0-x*x)
end function

! evaluate basis functions b_i(x) at point x
pure function basis(n, x)
	integer k, n; real x, basis(n); intent(in) n, x
	
	! Chebyshev polynomials
	forall (k=1:n) basis(k) = cos((k-1)*acos(x))
end function

! uniformly spaced grid of n points covering interval [-1,1]
pure function uniform(n)
	integer i, n; real uniform(n); intent(in) n
	
	! uniform grid
	forall (i=1:n) uniform(i) = (2*i-n-1.0)/(n-1.0)
end function

! Gauss-Lobatto grid of n points covering interval [-1,1]
pure function lobatto(n)
	integer i, n; real lobatto(n); intent(in) n
	
	! Gauss-Lobatto grid for Chebyshev polynomials
	forall (i=1:n) lobatto(i) = cos(pi*(n-i)/(n-1.0))
end function

subroutine lsolve(n, A, b)
  integer n
  real A(n,n), b(n),temp, multiplier, temporary(n)
  integer i, j, temp2
  do i=1, n !index row of matrix
     temp = A(i,i) !find the ith pivot
     temp2 = i
     do j=i+1, n 
        if (A(i,j) .ge. temp) then
           temp2 = j
        end if
     end do
     temporary = A(temp2,:) !copy the pivot row
     A(temp2,:) = A(i,:) !swap it with the ith row
     A(i,:) = temporary
     temp = b(temp2)
     b(temp2)=b(i)
     b(i)=temp
     do j=i+1, n !index column of matrix
        multiplier = A(j,i)/A(i,i) !compute multiplier based on pivot element
        A(j,:) = A(j,:)-multiplier*A(i,:) !subtract
        A(j,i) = 0 !prevent rounding error here
        b(j) = b(j)-multiplier*b(i) !
     end do
  end do
  b(n) = b(n)/A(n,n)
  do i=1,n-1 ! index row where i is the n-i'th row, i.e. solving the n-i'th eqn
     b(n-i) = b(n-i)
     do j=0,i-1 !index column where j is the n-j'th col, i.e. eliminating the n-j'th coeff.
        b(n-i)= b(n-i) - A(n-i,n-j)*b(n-j)
     end do
     b(n-i) = b(n-i) / A(n-i,n-i)
  end do
  
end subroutine lsolve

! evaluate polynomial expansion
subroutine dump(m)
  integer i, m, tempIndP, tempInd; real x(m), vPmax, vMax, temp, v2(m), vP(m)
  vPmax=0.0
  vMax =0.0
  x = uniform(m)
  v2 = f(x)
  vP = fP(x)
  tempInd = 1
  tempIndP = 1
  do i = 2,m-1 !Find the max |interpolant error|, |interpolant derivative error|
     temp=abs(sum(v*basis(n,x(i))) - v2(i))
     if (temp .GE. vMax) then
        vMax = temp
        tempInd =i
     end if
     temp=abs(sum(v*basisPrime(n,x(i))) - vP(i))
     if (temp .GE. vPmax) then
        vPmax = temp
        tempIndP = i
     end if
  end do
  write(*,*) "max error in P(x)", vMax, "occurs at", x(tempInd)
  write(*,*) "max error in P'(x)", vPmax, "occurrs at", x(tempIndP)
end subroutine

end program
