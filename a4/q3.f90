program pend
  !  Solving the equations for momentum conjugates p1, p2
  real, parameter :: dt = 0.05
  integer, parameter :: N=8
  real ::  y(4), m=1, l=1, g=9.8, xx, yy
  real ::  t = 0
  integer p,q 
  integer, parameter :: nx = 2*N, ny = 2*N
  real(8) :: xrange(2) = [-3.14156,3.14156], yrange(2) = [-3.14156, 3.14156]
  real(4), allocatable :: image(:,:,:)
  allocate(image(1,nx,ny))

  do p=0,2*N-1
     do q=0,2*N-1
        t=0
        y = [3.14156*(p-N)/(N-1), 0.0, 3.14156*(q-N)/(N-1), 0.0]
        do while ((abs(y(3)) .LE. 3.14156) .AND. (t .LE. 1000/sqrt(g)))
           call gl8(y, dt)
           t=t+dt
        end do
        image(1,p+1,q+1) = t
     end do
  end do
  call write2fits('fractal1.fit', image, xrange, yrange, ['N'])
  write(*,*)
  do p=0,2*N-1
     do q=0,2*N-1
        t=0
        y = [-1.5 + 3.14156*(p-N)/(3*N-1), 0.0, -3.14/4 + 3.14156*(q-N)/(3*N-1), 0.0]
        do while ((abs(y(3)) .LE. 3.14156) .AND. (t .LE. 1000/sqrt(g)))
           call gl8(y, dt)
           t=t+dt
        end do
        image(1,p+1,q+1) = t
     end do
  end do
  call write2fits('fractal2.fit', image, xrange, yrange, ['N'])

  do p=0,2*N-1
     do q=0,2*N-1
        t=0
        y = [-1.5 + 3.14156*(p-N)/(9*N-1), 0.0, -3.14/4 + 3.14156*(q-N)/(9*N-1), 0.0]
        do while ((abs(y(3)) .LE. 3.14156) .AND. (t .LE. 1000/sqrt(g)))
           call gl8(y, dt)
           t=t+dt
        end do
        image(1,p+1,q+1) = t
     end do
  end do
  call write2fits('fractal3.fit', image, xrange, yrange, ['N'])

  do p=0,2*N-1
     do q=0,2*N-1
        t=0
        y = [-1.5 + 3.14156*(p-N)/(27*N-1), 0.0, -3.14/4 + 3.14156*(q-N)/(27*N-1), 0.0]
        do while ((abs(y(3)) .LE. 3.14156) .AND. (t .LE. 1000/sqrt(g)))
           call gl8(y, dt)
           t=t+dt
        end do
        image(1,p+1,q+1) = t
     end do
  end do
  call write2fits('fractal4.fit', image, xrange, yrange, ['N'])



contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! write array data into FITS file as sequence of image extensions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine write2fits(file, array, xx, yy, vars, coords)
    character(len=*) file, vars(:), coords
    real(4) array(:,:,:); real(8) xx(2), yy(2)
    optional xx, yy, vars, coords

    integer i, j, status, unit
    integer :: hdus, naxis = 2, n(2), npix
    integer :: bitpix = -32, group = 1, blocksize = -1

    ! data dimansions
    hdus = size(array,1)
    n(1) = size(array,2)
    n(2) = size(array,3)
    npix = n(1)*n(2)

    ! delete file if it already exists
    open(unit=1234, iostat=status, file=file, status='old')
    if (status == 0) close(1234, status='delete'); status = 0

    ! initialize FITS file
    call ftgiou(unit, status)
    call ftinit(unit, file, blocksize, status)

    ! write image extensions
    do i = 1,hdus
       call ftiimg(unit, bitpix, naxis, n, status)
       call ftppre(unit, group, 1, npix, array(i,:,:), status)

       if (present(vars)) then
          if (present(coords)) then
             call ftpkys(unit, 'EXTNAME', vars(i)//coords, 'variable stored in extension', status)
          else
             call ftpkys(unit, 'EXTNAME', vars(i), 'variable stored in extension', status)
          end if
       end if
       if (present(xx)) then
          call ftpkyj(unit, 'CRPIX1', 1, 'x-axis origin pixel', status)
          call ftpkyd(unit, 'CRVAL1', xx(1), 14, 'x-axis origin coordinate', status)
          call ftpkyd(unit, 'CDELT1', (xx(2)-xx(1))/n(1), 14, 'x-axis increment', status)
       end if
       if (present(yy)) then
          call ftpkyj(unit, 'CRPIX2', 1, 'y-axis origin pixel', status)
          call ftpkyd(unit, 'CRVAL2', yy(1), 14, 'y-axis origin coordinate', status)
          call ftpkyd(unit, 'CDELT2', (yy(2)-yy(1))/n(2), 14, 'y-axis increment', status)
       end if
    end do

    ! clean up
    call ftclos(unit, status)
    call ftfiou(unit, status)
  end subroutine write2fits

  subroutine evalf(y, dydx)
    real y(4), dydx(4), kinetic, potential; intent(in) y
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
    real y(n), temp(n), g(n,s), dt; integer i, k

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
          call evalf(y + g(:,i)*dt , g(:,i))
       end do
    end do

    ! update the solution
    y = y+dt* matmul(g,b)
 
  end subroutine gl8

end program
