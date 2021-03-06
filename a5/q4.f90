program spectral; implicit none
  real, parameter :: pi = 3.141592653589793238462643383279502884197169399375Q0
  integer, parameter :: nn = 200
  real x(nn), psi(nn), h, L(nn,nn), Hamiltonian(nn,nn), identity(nn,nn), potential(nn), lambda, matinv(nn,nn)
  integer i

  call initgrid(x, potential)
  call initl()
  call inith()

  identity = 0
  do i=1,nn
     identity(i,i) = 1
  end do

  ! initialize initial solution guess
  call evalb(2,x,psi)

  open(unit=1, file='RConv.dat')
  
  do i = 1,20
     lambda = dot_product(psi,(matmul(Hamiltonian, psi)))/dot_product(psi,psi)
     write(1,*) i, log10(abs(lambda-0.5))
     matinv = Hamiltonian - lambda*identity
     psi = lsolve(matinv, psi)
     psi = psi/sqrt(sum(psi**2)*h)
  end do

contains
  ! ititialize grid
  subroutine initgrid(x, potential)
    integer i; real x(nn), potential(nn)
    ! Initialize grid
    forall(i=1:nn) x(i) = atanh(cos(pi*(nn-i+0.5)/nn))
    ! initialize potential
    forall (i=1:nn) potential(i) = V(x(i))
    h=(x(nn)-x(1))/nn
  end subroutine initgrid


  ! initialize Laplacian operator matrix
  subroutine initl()
    integer i, pivot(nn), status; real A(nn,nn), B(nn,nn)

    ! evaluate basis and differential operator values on a grid
    do i = 1,nn
       call evalb(i+1, x, A(i,:), B(i,:))
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
    L = transpose(B)
  end subroutine initl

  subroutine inith()
    Hamiltonian = -0.5*L
    do i=1,nn
       Hamiltonian(i,i) = Hamiltonian(i,i) + potential(i)
       !write(*,*) L(i,:)
       !write(*,*) Hamiltonian(i,:)
    end do
  end subroutine inith

  ! Define potential
  pure function V(x)
    real V, x
    intent(in) x
    V = 0.5*x**2        
  end function V

  ! evaluate basis functions and derivatives
  subroutine evalb(n, x, Tn, Tnxx)
    integer k, n 
    real theta(nn), x(nn), Tn(nn), Tnxx(nn); optional Tnxx

    theta = n*acos(tanh(x))

    ! Modified Chebyshev polynomials
    Tn = cos(theta)
    if (n>1) then
       if (mod(n,2) == 0) then
          Tn = Tn - 1
       else if (mod(n,2) == 1) then
          Tn = Tn - x
       end if
    end if
    ! Second derivative of Chebyshev polynomials
    if (present(Tnxx)) Tnxx = -n*(sinh(x)*sin(theta) + n*cos(theta))/cosh(x)**2
  end subroutine evalb


  ! solve linear problem Ax = B
  function lsolve(A, B)
    real lsolve(nn), psi(nn), rhs(nn), A(nn,nn), B(nn)
    integer i, pivot(nn), status

    ! find static solution by direct inversion
    status = 0
    select case (kind(A))
    case(4); call sgesv(nn, 1, A, nn, pivot, B, nn, status)
    case(8); call dgesv(nn, 1, A, nn, pivot, B, nn, status)
    case default; call abort
    end select

    ! bail at first sign of trouble
    if (status /= 0) call abort

    ! return solution
    lsolve = B
  end function lsolve

end program spectral
