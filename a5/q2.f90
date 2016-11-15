program q2
  real, parameter :: hBar=1, w = 1, m=1
  real, parameter :: h = 1E-1
  real, parameter :: init(2) =[1.0,0.0]
  real, parameter :: bullseye = 0
  real, parameter :: xmax = 3
  integer i

  do i=1,10
     write(*,*) bisect(i*1.0,i+1.0,1.0)
  end do

contains
  function bisect(a0, b0, odd)
    real a0, b0, bisect
    real a, b, c, fa, fb, fc, odd
    real, parameter :: eps = 1.0e-16

    a = a0; fa = boundary(a,odd)-bullseye
    b = b0; fb = boundary(b,odd)-bullseye

    if (fa*fb > 0) call abort

    do while (abs(b-a) > eps*10)
       c = (a+b)/2.0; fc = boundary(c,odd)-bullseye; if (fc == 0.0) exit
       if (fa*fc < 0.0) then
          b = c; fb = fc;
       else if (fc*fb < 0.0) then
          a = c; fa = fc
       else
          exit
       end if
    end do
    bisect = c
  end function bisect
  
  function boundary(E, odd)
    real boundary
    real psi(2), psiMinus(2)
    real x, odd
    intent(in) E
    x=0
    psi = init
    psiMinus = -1*init
    do while (x .LE. xmax)
       call gl8(psi,h,x, E)
       call gl8(psiMinus,h,-1*x, E)
       x=x+h
    end do
    boundary = psi(1) + odd*psiMinus(1)
  end function boundary
  
  function V(u)
    real :: u, V
    V = 0.5*m*(u*w)**2
  end function V

  subroutine evalf(psi, psiP,x, E)
    real :: psi(2), psiP(2), x, E
    psiP(1) = psi(2)
    psiP(2) = psi(1)*(V(x) - E)*2*m/(hbar**2)
  end subroutine evalf

  ! 8th order implicit Gauss-Legendre integrator
  subroutine gl8(psi,h, x, E)
    integer, parameter :: s = 4, n = 2
    real psi(n), g(n,s), h,x ,E; integer i, k

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
          call evalf(psi + g(:,i)*h, g(:,i), x, E)
       end do
    end do

    ! update the solution
    psi = psi + matmul(g,b)*h
  end subroutine gl8
end program q2
