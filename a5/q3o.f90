program q3; implicit none

  real, parameter :: dx = 1e-4           
  integer, parameter :: n = 10
  integer i, flag, flagup
  real y(2), points(n,2), start(2), delta
  delta = 0.5
  i=0
  flag = 0
  start = 0

  do while (i < n)
     flag = 0
     flagup = 0
     points(i,:)= bisect(start(1), start(1)+delta,1, flag)
     if(flag == 1) then
        start(1) = start(1)+delta
     else
        start(1) = points(i,1)+0.00001
        i=i+1
        write(*,*) start(1)
     end if

     points(i,:)= bisect(start(2), start(2)+delta,-1,flagup)
     if (flagup == 1) then
        flag = 1
        start(2) = start(2)+delta
     else
        start(2) = points(i,1)+0.00001
        i=i+1
        write(*,*) start(2)
     end if
     write(*,*) i, start, delta
  end do

  write(*,*) points(:,1)
 
contains
  subroutine normalizedDist(energy, norm, init)
    real energy, norm, x, init(2), y(2)
    character (len = 40) :: filename
    x = 0.0

    write(filename, fmt = '(a, F4.1, a)') "Norm_Anharmonic_E",energy,".dat"
    open(unit = 1, file = filename)

    y = init     
    write(1,*) x, y(1)**2/(2*norm)
    do while ( abs(y(1)) < 4.0 .and. ( abs(y(1)) > 0.1 .or. abs(y(2)) > 0.1 ) )       
       call gl8(x, y, dx, energy)
       x = x + dx
       write(1,*) x, y(1)**2/(2*norm)
    end do
    close(1)
  end subroutine normalizedDist


  function calc(energy,even)
    real x, y(2)
    real energy, calc(2), norm
    integer even

    x = 0.0
    if(even == 1) then
       y = (/ 1.0, 0.0 /) ! Even Initial Condition
    else
       y = (/ 0.0, 1.0 /) ! Odd initial condition
    end if

    norm = y(1)**2

    do while ( abs(y(1)) < 4.0 .and. ( abs(y(1)) > 0.1 .or. abs(y(2)) > 0.1 ) )       
       call gl8(x, y, dx, energy)
       x = x + dx
       norm = norm + dx*(y(1)**2) ! Continuing to integrate the normalization factor
    end do
    calc = [y(1),norm]      ! Tells us whether we blow up to positive or negative infinity
  end function calc

  function bisect(Emin, Emax, even, flag)
    real Emin, Emax, bisect(2)
    real a, b, c, fa(2), fb(2), fc(2)
    real, parameter :: eps = 1.0e-12
    integer even, flag
    flag = 0

    a = Emin; fa = calc(a,even)
    b = Emax; fb = calc(b,even)

    if (fb(1)*fa(1) > 0) then
       flag = 1
       b=a
    end if

    do while (abs(b-a) > eps)
       c = (a+b)/2.0; fc = calc(c,even)
       if (abs(fc(1)) < 0.001) exit
       if (fa(1)*fc(1) < 0.0) then
          b = c; fb = fc
       end if
       if (fc(1)*fb(1) < 0.0) then
          a = c; fa = fc
       end if
    end do

    bisect = [c, fc(2)]
  end function bisect

  pure function V(x)
    real V, x
    intent(in) x
    V = 0.5*x**4        
  end function V

  subroutine evalf(x, y, dydx, E)
    real x, y(2), dydx(2), E

    dydx(1) = y(2)
    dydx(2) = 2*(V(x) - E)*y(1)
  end subroutine evalf

  ! 8th order implicit Gauss-Legendre integrator
  subroutine gl8(x, y, dx, E)
    integer, parameter :: s = 4, n = 2
    real x, y(n), g(n,s), E, dx; integer i, k

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
    do k = 1,12
       g = matmul(g,a)
       do i = 1,s
          call evalf(x, y + g(:,i)*dx, g(:,i), E)
       end do
    end do
    
    ! update the solution
    y = y + matmul(g,b)*dx
  end subroutine gl8
  
end program q3

