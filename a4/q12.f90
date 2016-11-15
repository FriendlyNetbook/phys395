program pend
  !  Solving the equations for momentum conjugates p1, p2
  real(8), parameter :: dt =  0.0005, tMax = 100*sqrt(1/9.8)

  real(8) ::  y(4), yp(4), m=1, l=1, g=9.8,  En, Energy
  real(8) ::  t = 0
  integer :: energyOverflows = 0
  y=[3.14156/3.0, 0.0, -3.14156/3.0 ,0.0]
  yp = [3.0, 0.0, 0.0, 0.0]
  call evalf(y, yp)
  En = Energy

  open(unit = 1, file = "energyError")	
  open(unit = 2, file = "traceData")
  open(unit = 3, file = 'animateData')

  do while (t .LE. tMax)
     write(1,*) abs(Energy/En -1)
     write(3,*) 0, 0
     write(3,*) l*sin(y(1)), -l*cos(y(1))
     write(3,*) l*(sin(y(1))+sin(y(3))), -l*(cos(y(1))+cos(y(3)))
     write(3,*) " "
     write(3,*) " "
     write(2,*)  l*(sin(y(1))+sin(y(3))), -l*(cos(y(1))+cos(y(3)))
     if (abs(Energy/En -1) .GE. 1e-12) then
        write(*,*) "Energy Fluctutation of", abs(Energy/En -1)
        energyOverflows = 1
     end if

     call gl8(y, dt)
     t=t+dt
  end do
  if (energyOverflows == 0) then
     write(*,*) "The energy fluctuations stayed below 1E-12 with a timestep of ", dt
  end if

contains 
  subroutine evalf(y, dydx)
    real(8) y(4), dydx(4), kinetic, potential; intent(in) y
    dydx(1) = 6.0/(l**2) * (2.0*y(2) - 3.0*y(4)*cos(y(1) - y(3)))/(16.0 - 9.0*cos(y(1)-y(3))*cos(y(1) - y(3)))
    dydx(3) = 6.0/(l**2) * (8.0*y(4) - 3.0*y(2)*cos(y(1) - y(3)))/(16.0 - 9.0*cos(y(1)-y(3))*cos(y(1) - y(3)))
    dydx(2) = -0.5 * l**2 * (dydx(1) * dydx(3) * sin(y(1) - y(3)) + 3.0*(g/l)*sin(y(1)))
    dydx(4) = -0.5 * l**2 * (-dydx(1) * dydx(3) * sin(y(1) - y(3)) + (g/l)*sin(y(3)))

    kinetic = (1/6.0)*m*l*l*(dydx(3)**2 + 4*dydx(1)*dydx(1) + 3*dydx(1)*dydx(3)*cos(y(1)-y(3)))
    potential = (1/2.0)*m*g*l*(3*cos(y(1)) + cos(y(3)))
    Energy = kinetic - potential


  end subroutine evalf

  ! 8th order implicit Gauss-Legendre integrator
  subroutine gl8(y, dt)
    integer, parameter :: s = 4, n = 4
    real(8) y(n), temp(n), g(n,s), dt; integer i, k

    ! Butcher tableau for 8th order Gauss-Legendre method
    real(8), parameter :: a(s,s) = reshape((/ &
         0.869637112843634643432659873054998518Q-1, -0.266041800849987933133851304769531093Q-1, &
         0.126274626894047245150568805746180936Q-1, -0.355514968579568315691098184956958860Q-2, &
         0.188118117499868071650685545087171160Q0,   0.163036288715636535656734012694500148Q0,  &
         -0.278804286024708952241511064189974107Q-1,  0.673550059453815551539866908570375889Q-2, &
         0.167191921974188773171133305525295945Q0,   0.353953006033743966537619131807997707Q0,  &
         0.163036288715636535656734012694500148Q0,  -0.141906949311411429641535704761714564Q-1, &
         0.177482572254522611843442956460569292Q0,   0.313445114741868346798411144814382203Q0,  &
         0.352676757516271864626853155865953406Q0,   0.869637112843634643432659873054998518Q-1 /), (/s,s/))
    real(8), parameter ::   b(s) = (/ &
         0.173927422568726928686531974610999704Q0,   0.326072577431273071313468025389000296Q0,  &
         0.326072577431273071313468025389000296Q0,   0.173927422568726928686531974610999704Q0  /)

    ! iterate trial steps
    g = 0.0
    do k = 1,16
       g = matmul(g,a)
       do i = 1,s
          call evalf(y + g(:,i)*dt , g(:,i))
       end do
    end do

    ! update the solution
    y = y+dt* matmul(g,b)
 
  end subroutine gl8

end program
