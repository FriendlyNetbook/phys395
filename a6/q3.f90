program q3; implicit none
  real, parameter :: pi = 3.141592653589793238462643383279502884197169399375Q0
  real, parameter :: dt = 0.01
  integer, parameter :: nn = 250, tmax = 40
  integer, parameter :: timesteps = tmax/dt

  integer i, j
  real dx, x(nn), L(nn,nn), V(nn), II(nn,nn), H(nn,nn)
  complex psi(nn), psixx(nn), Qinv(nn,nn), S(nn,nn)

  integer ipiv, info
  complex work(nn,nn)

  ! initialize grid, V, L, I, avg space step
  call init(x, V, II); call initl(); dx = (x(nn)-x(1))/nn

  ! Calculate H
  H = -0.5*L
  do i=1,nn
     H(i,i) = H(i,i) + V(i)
  end do

  ! Calculate Q
  Qinv = dcmplx(II, 0.5*H*dt)
  call zgetrf( nn, nn, Qinv, nn, ipiv, info)
  call zgetri( nn, Qinv, nn, ipiv, work, nn, info)
  S = matmul(Qinv,dcmplx(II, -0.5*dt*H))

  ! Initial Conditions
  psi = exp(-0.5*(x-1)**2)
  psi = psi/sqrt(dot_product(psi,psi)*dx)

  open(unit = 1, file = 'q3data.dat')

  ! Time-Evolve
  do i=0,timesteps
     if (mod(i,15) == 0) then
        write(1,*)
        write(1,*)
        do j=1,nn
           write(1,*) x(j), realpart(psi(j)), imagpart(psi(j)), realpart(psi(j)*conjg(psi(j)))
        end do
     end if
     psi = matmul(S, psi)
  end do

contains

  ! Init grid, V, identity
  subroutine init(x, V, II)
    real x(nn), V(nn), II(nn,nn)
    integer k

    II = 0.0
    forall(k=1:nn) x(k) = atanh(cos(pi*(nn-k+0.5)/nn)) ! Gridpoints
    forall(k=1:nn) V(k) = 0.5*x(k)**2 ! Harmonic Potential
    forall(k=1:nn) II(k,k) = 1.0 !I
  end subroutine init

  ! initialize Laplacian operator matrix
  subroutine initl()
    integer i, pivot(nn), status; real A(nn,nn), B(nn,nn)

    do i = 1,nn
       call evalb(i-1, x, A(i,:), B(i,:))
    end do
    status = 0
    select case (kind(A))
    case(4); call sgesv(nn, nn, A, nn, pivot, B, nn, status)
    case(8); call dgesv(nn, nn, A, nn, pivot, B, nn, status)
    case default; call abort
    end select
    if (status /= 0) call abort
    L = transpose(B)
  end subroutine initl
  
  ! Evaluate Basis Functions
  subroutine evalb(n, x, psi, psixx)
    integer n
    real theta(nn), x(nn),  psi(nn), psixx(nn); optional psixx

    theta = n*acos(tanh(x)); psi = cos(theta)
    if (n>1) then
       if (mod(n,2)==0) then; psi = psi - 1
       else if (mod(n,2)==1) then; psi = psi - x
       end if
    end if
    if (present(psixx)) psixx = -n*(sinh(x)*sin(theta) + n*cos(theta))/cosh(x)**2
  end subroutine evalb

end program q3
