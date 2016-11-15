! solve A.x = B using Gauss-Jordan elimination
! A gets destroyed, answer is returned in B

subroutine gaussj(n, A, B)
  integer n,i,j,k,p(n,n),maxElement; real A(n,n), B(n), M(n), maxElementValue, tempRow, temp
  p = 0 !q is the permutation matrix
  do i=1,n
     p(i,i)=1
  end do
   
  do i=1,n
     do j=1,n
        if (A(i,j) < maxElementValue)
        maxElementValue = A(i,j)
        maxElement = j
        end if
     end do
     tempRow = A(:,i) !pivot A
     A(:,i) = A(:,maxElement)
     A(:,maxElement) = tempRow
     tempRow = p(:,i) !pivot p
     p(:,i) = p(:,maxElement)
     p(:,maxElement) = tempRow
     temp = B(i) !pivot B
     B(i) = B(maxElement)
     B(maxElement) = temp
     A(i,1:j-1)=0
     do j=i+1,n
        j=-1*B(i,j)/A(i,i)
        do k=j,n
           A(i,k)=A(i,k)+A(i,j)*m
           B(k)=B(k)+B(j)*m
        end do
  end do

  do i=nm1
     do j=n,i+1
        C(i)=C(i)+A(i,j)*C(j)
     end do
     print (*,*) C(i) !Need to multiply c by p-transpose to reorder
  end do
  
end subroutine
