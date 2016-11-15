program spectral; implicit none

  real, parameter :: pi = 3.141592653589793238462643383279502884197169399375Q0
  integer, parameter :: nn = 150

  real x(nn), phi(nn), H(nn,nn), lam

  integer i

  ! initialize grid and Laplacian operator matric
  call initg(x)
  call initH()

  ! initialize initial solution guess
  call evalb(3,x,phi)

  ! try to relax using Newton's iteration
  do i = 1,50
     phi = itPhi(phi,lambda(phi,H))
  end do

contains

  function lambda(psi, H)
    real psi(nn), H(nn,nn), BOM(nn)
    real a,lambda

    a=0
    do i=1,nn
       a=a+psi(i)**2
    end do

    lambda=0
    BOM = matmul(H,psi)
    do i=1,nn
       lambda=lambda+psi(i)*BOM(i)
    end do

  end function lambda
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! derivative of the potential
  elemental function DV(phi)
    real DV, phi; intent(in) phi

    DV = (phi**2 - 1.0) * phi
  end function DV

  ! second derivative of the potential
  elemental function DDV(phi)
    real DDV, phi; intent(in) phi

    DDV = 3.0*phi**2 - 1.0
  end function DDV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! solve linear Laplace problem [L - m_eff^2] phi = RHS
  function itPhi(psi, lamb)
    real psi(nn)
    real itPhi(nn)
    real A(nn,nn)
    real lamb
    real pivot(nn)
    integer status

    A = H
    do i = 1,nn
       A(i,i) = A(i,i) - lamb
    end do

    ! find static solution by direct inversion
    status = 0
    select case (kind(A))
    case(4); call sgesv(nn, 1, A, nn, pivot, psi, nn, status)
    case(8); call dgesv(nn, 1, A, nn, pivot, psi, nn, status)
    case default; call abort
    end select

    ! bail at first sign of trouble
    if (status /= 0) call abort

    ! return solution
    itPhi = psi
  end function itPhi


  ! initialize collocation grid
  subroutine initg(x)
    integer i; real x(nn)

    forall(i=1:nn) x(i) = atanh(cos(pi*(nn-i+0.5)/nn))
  end subroutine initg


  ! evaluate basis function and its derivatives on a grid
  subroutine evalb(n, x, Tn, Tnxx)
    integer n; real, dimension(nn) :: x, theta, thetaO,Tn, Tnxx; optional Tnxx
    if (mod(n,2) ==0) then
       theta = n*acos(tanh(x)); Tn = cos(theta)-1.0
    else
       thetaO = acos(tanh(x))
       theta = n*acos(tanh(x)); Tn = cos(theta)-cos(thetaO)
    end if
    if (present(Tnxx)) then
       Tnxx = -n*(sinh(x)*sin(theta) + n*cos(theta))/cosh(x)**2
    end if
  end subroutine evalb


  ! ! evaluate basis function and its derivatives on a grid
  ! subroutine evalb(n, x, Tn, Tnxx)
  !   integer n; real, dimension(nn) :: x, theta, Tn, Tnxx; optional Tnxx

  !   theta = n*acos(tanh(x)); Tn = cos(theta)
  !   if (present(Tnxx)) Tnxx = -n*(sinh(x)*sin(theta) + n*cos(theta))/cosh(x)**2
  ! end subroutine evalb

  ! initialize Laplacian operator matrix
  subroutine initH()
    integer i, pivot(nn), status; real A(nn,nn), B(nn,nn)

    ! evaluate basis and differential operator values on a grid
    do i = 1,nn
       call evalb(i-1, x, A(i,:), B(i,:))
    end do

    ! find linear differentiation matrix
    status = 0
    select case (kind(A))
    case(4); call sgesv(nn, nn, A, nn, pivot, B, nn, status)
    case(8); call dgesv(nn, nn, A, nn, pivot, B, nn, status)
    case default; call abort
    end select

    ! bail at first sign of trouble
    if (status /= 0) call abort

    ! to evaluate Laplacian of function f, simply do matmul(L,f)
    H = transpose(B)
  end subroutine initH


end program spectral
