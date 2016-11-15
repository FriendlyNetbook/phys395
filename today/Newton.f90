program Newton; implicit none
  M=3  ! M is the size of the wings
contains
  
  function MPtMin(x,M)
    integer, parameter :: N = 90; intent(in) N
    integer, parameter :: M; intent(in) M
    real x(2*M+1)

    do i = M+1, N-M
       x = lobatto(2*M+1)
       

    end do
    
  end function 3PtMin

  elemental function f(x)
    real f, x; intent(in) x
    f= x*x-4
  end function f
  
    ! Gauss-Lobatto grid of n points covering interval [-1,1]
pure function lobatto(n)
        integer i, n; real lobatto(n); intent(in) n
        ! Gauss-Lobatto grid for Chebyshev polynomials
        forall (i=1:n) lobatto(i) = cos(pi*(n-i)/(n-1.0))
      end function lobatto

      !evaluate basis functions b_i(x) at point x
pure function basis(n, x)
        integer k, n; real x, basis(n); intent(in) n, x

        ! Chebyshev polynomials
        forall (k=1:n) basis(k) = cos((k-1)*acos(x))
end function


  ! solve A.x = B using LAPACK canned routine (A gets destroyed, the answer is returned in B)
  subroutine lsolve(n, A, B)
    integer n; real A(n,n), B(n)
    
    integer status, pivot(n)
    
    ! find static solution by direct inversion
    status = 0
    select case (kind(A))
    case(4); call sgesv(n, 1, A, n, pivot, B, n, status)
    case(8); call dgesv(n, 1, A, n, pivot, B, n, status)
    case default; call abort
    end select
 
    ! bail at first sign of trouble
    if (status /= 0) call abort
end subroutine

end program Newton
