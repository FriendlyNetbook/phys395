! gaussj.f90  -  starting point for Gauss-Jordan elimination
! compile with: gfortran -O3 -fdefault-real-8 gaussj.f90

program gaussj; implicit none
integer, parameter :: n = 2
real A(n,n), B(n)

A(1,1) = 1.0
A(1,2) = 7.0
A(2,1) = 5.0
A(2,2) = 1.0

B(1) = 1.0
B(2) = 3.0

call lsolve(n, A, B)
write (*,*) B

contains
! solve A.x = B using Gauss-Jordan elimination
! A gets destroyed, answer is returned in B
  subroutine lsolve(n, A, b)
    integer n
    real A(n,n), b(n), temp, multiplier, temporary(n)
    integer i, j, temp2
    do i=1, n !index row of matrix
       temp = A(i,i) !find the ith pivot
       temp2 = i !store the pivot element number
       do j=i+1, n 
          if (A(i,j) .ge. temp) then
             temp2 = j !iterate through elements below a pivot to find the max
          end if
       end do
       temporary = A(temp2,:) !copy the pivot row
       A(temp2,:) = A(i,:) !swap it with the ith row
       A(i,:) = temporary 
       temp = b(temp2)!do the same operation on the b vector
       b(temp2)=b(i)
       b(i)=temp
       do j=i+1, n !index column of matrix
          multiplier = A(j,i)/A(i,i) !compute multiplier based on pivot element
          A(j,:) = A(j,:)-multiplier*A(i,:) !subtract
          A(j,i) = 0 !set elements below pivot to zero
          b(j) = b(j)-multiplier*b(i) !subtract b as well
       end do
    end do
    b(n) = b(n)/A(n,n)
    do i=1,n-1 ! index row where i is the n-i'th row, i.e. solving the n-i'th eqn
       b(n-i) = b(n-i)
       do j=0,i-1 !index column where j is the n-j'th col, i.e. eliminating the n-j'th coeff.
          b(n-i)= b(n-i) - A(n-i,n-j)*b(n-j) !back-substitute
       end do
       b(n-i) = b(n-i) / A(n-i,n-i)
    end do
  end subroutine lsolve
  
end program 
