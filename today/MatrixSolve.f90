program MatrixSolve; implicit none
  integer, parameter :: n =2 
  real :: A(n,n) = reshape([1, 5, 7, 1],[n,n]) !Matrix
  real :: b(n) = [1,3], x(n)
  real :: pivot(n,n) = 0!Keep Track of pivoting
  real :: temporary(n), multiplier
  integer :: i = 0, j = 0, temp2 = 0, t
  real :: temp = 0.0

  do i=1, n !Make Pivot-tracking matrix initially I
     pivot(i,i) = 1
  end do
  do i=1, n !index row of matrix
     temp = A(i,i) !find the ith pivot
     temp2 = i
     do j=i+1, n 
        if (A(i,j) > temp) then
           temp2 = j
        end if
     end do
     temporary = A(temp2,:) !copy the pivot row
     A(temp2,:) = A(i,:) !swap it with the ith row
     A(i,:) = temporary
     temporary = pivot(temp2,:) ! do the same for the pivot matrix
     pivot(temp2,:) = pivot(i,:)
     pivot(i,:) = temporary
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
  !x(n) = b(n)/A(n,n)
  do i=1,n-1 ! index row where i is the n-i'th row, i.e. solving the n-i'th eqn
     x(n-i) = b(n-i)
     do j=0,i-1 !index column where j is the n-j'th col, i.e. eliminating the n-j'th coeff.
        x(n-i)= x(n-i) - A(n-i,n-j)*x(n-j)
     end do
     x(n-i) = x(n-i) / A(n-i,n-i)
  end do
  write(*,*) x
end program MatrixSolve
