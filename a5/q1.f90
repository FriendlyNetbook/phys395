program q1
  real, parameter :: hBar=1, w = 1, m=1
  real, parameter :: h = 1E-4


  real :: init(2) =[0.0,1.0]
  call integrate(0.0)
  call integrate(0.2)
  call integrate(0.4)
  call integrate(0.6)
  call integrate(0.8)
  call integrate(1.0)
  call integrate(1.2)
  call integrate(1.4)
  call integrate(1.6)
  call integrate(1.8)
contains
  subroutine integrate(E)
    real psi(2)
    real minusPsi(2)
    real x
    real E
    character (len=20) :: filename
    intent(in) E
    write(filename, fmt = '(a, F3.1, a)') "Eis",E, "Minus.dat"
    open(unit = 1, file = filename)
    filename =''
    write(filename, fmt = '(a, F3.1, a)') "Eis",E, "Plus.dat"
    open(unit = 2, file = filename)
    x=0
    psi = init
    minusPsi= -1*init
    do while  ( abs(psi(1)) < 4.0 .and. ( abs(psi(1)) > 0.1 .or. abs(psi(2)) > 0.1 ) )       
       call gl8(psi,h,x,E) 
       call gl8(minusPsi,h,-1*x,E)
       write(1,*) x, (psi-minusPsi)/2
       write(2,*) x, (psi+minusPsi)/2
       x=x+h
    end do
    close(1)
    close(2)
  end subroutine integrate
  
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
  subroutine gl8(psi,h, x,E)
    integer, parameter :: s = 4, n = 2
    real psi(n), g(n,s), h,x, E; integer i, k

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
end program q1
