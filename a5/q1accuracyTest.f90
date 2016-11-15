program spectral; implicit none
  real, parameter :: pi = 3.141592653589793238462643383279502884197169399375Q0, hbar = 1, m=1
  integer, parameter :: nn = 200
  real, parameter :: h = 0.001
  real x(nn), psi(nn), L(nn,nn), Hamiltonian(nn,nn), identity(nn,nn), potential(nn), lambda, matinv(nn,nn)
  integer i
  real :: E=1.5

  call initgrid(x, potential)
  call initl()
  call inith()

  identity = 0
  do i=1,nn
     identity(i,i) = 1
  end do

  call gl8(psi,h,E)
  write(*,*) abs(matmul(Hamiltonian,psi)-E*psi)
  
  

contains

  ! Calculates the derivative of psi
  subroutine evalf(psi, psiP, E)
    real :: psi(2), psiP(2),  E
    psiP(1) = psi(2)
    psiP(2) = psi(1)*(V(x) - E)*2*m/(hbar**2)
  end subroutine evalf

  ! 8th order implicit Gauss-Legendre integrator
  subroutine gl8(psi,h,E)
    integer, parameter :: s = 4, n = 2
    real psi(n), g(n,s), h, E; integer i, k

    ! Butcher tableau for 8th order Gauss-Legendre method
    real, parameter :: a(s,s) = reshape((/ &
         0.869637112843634643432659873054998518Q-1, -0.266041800849987933133851304769531093Q-1, &
         0.126274626894047245150568805746180936Q-1, -0.355514968579568315691098184956958860Q-2, &
         0.188118117499868071650685545087171160Q0,   0.163036288715636535656734012694500148Q0,  &
         -0.278804286024708952241511064189974107Q-1,  0.673550059453815551539866908570375889Q-2, &
         0.167191921974188773171133305525295945Q0,   0.353953006033743966537619131807997707Q0,  &
         0.163036288715636535656734012694500148Q0,  -0.141906949311411429641535704761714564Q-1, &
         0.177482572254522611843442956460569292Q0,   0.313445114741868346798411144814382203Q0,  &
         0.352676757516271864626853155865953406Q0,   0.869637112843634643432659873054998518Q-1 /), (/s,s/))
    real, parameter ::   b(s) = (/ &
         0.173927422568726928686531974610999704Q0,   0.326072577431273071313468025389000296Q0,  &
         0.326072577431273071313468025389000296Q0,   0.173927422568726928686531974610999704Q0  /)

    ! iterate trial steps
    g = 0.0
    do k = 1,16
       g = matmul(g,a)
       do i = 1,s
          call evalf(psi + g(:,i)*h, g(:,i), E)
       end do
    end do

    ! update the solution
    psi = psi + matmul(g,b)*h
  end subroutine gl8

  ! ititialize grid
  subroutine initgrid(x, potential)
    integer i; real x(nn), potential(nn)
    ! Initialize grid
    forall(i=1:nn) x(i) = atanh(cos(pi*(nn-i+0.5)/nn))
    ! initialize potential
    forall (i=1:nn) potential(i) = V(x(i))
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
