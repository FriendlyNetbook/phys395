program putz
  integer, parameter :: nn =150
  real, parameter :: pi = 3.141592653521
  real x(nn), phi(nn)
  integer i
  call initg(x)
  call evalb(1,x,phi,phi)
  do i=1,nn
     write(*,*) i, phi(i)
  end do
  
contains
  ! Initialize grid
  subroutine initg(x)
    integer i
    real x(nn)
    do i=1,nn
       x(i) = atanh(cos(pi*(nn-i+0.5)/nn))
    end do
  end subroutine initg
  

  ! evaluate basis function & derivatives
  subroutine evalb(n,x,Tn,Tnxx)
    integer n
    real,dimension(n) :: x, Tn,theta, Tnxx; optional Tnxx
    theta=n*acos(tanh(x))
    Tn=cos(theta)
    if(present(Tnxx)) then
       Tnxx = -n*(sinh(x)*sin(theta) + n*cos(theta))/cosh(x)**2
    end if
    
  end subroutine evalb
  
end program putz
