program q1
  real, parameter :: hBar=1, w = 1, m=1
  real, parameter :: h = 1E-4

  real :: init(2) =[0.0,1.0]
  integer j

  do j=0,9
     call integrate(j*0.2, init)
  end do

  init = [1.0, 0.0]
  do j=0,9
     call integrate(j*0.2,init)
  end do
  
contains
  subroutine integrate(E, init)
    real psi(2) ! Wavefunction
    real minusPsi(2) ! Wavefunction given by psi(-x)
    real x, init(2)
    real E
    character (len=40) :: filename, filename2
    write(filename, '(a, F3.1, a, F2.0,F2.0,a)') "Eis",E,"init",init(1),init(2), "Odd.dat"
    open(unit = 1, file = filename)

    write(filename2, fmt = '(a, F3.1, a, F2.0, F2.0, a)') "Eis",E,"init",init(1),init(2), "Even.dat"
    open(unit = 2, file = filename2)

    x=0
    psi = init
    minusPsi= init
    do while  ( abs(psi(1)) < 4.0 .and. ( abs(psi(1)) > 0.1 .or. abs(psi(2)) > 0.1 ) )
       call gl8(psi,h,x,E)
       call gl8(minusPsi,-1*h,-1*x,E)

       write(2,*) x, (psi+minusPsi)/2 ! write the even
       write(2,*) -1*x, (psi+minusPsi)/2 ! write the even
       write(1,*) x, (psi-minusPsi)/2 ! and odd components
       write(1,*) -1*x, -1*(psi-minusPsi)/2 ! and odd components
       x=x+h
    end do

    close(1)
    close(2)
  end subroutine integrate
  
  function V(u)
    real :: u, V
    V = 0.5*m*(u*w)**2
  end function V

  ! Calculates the derivative of psi
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
