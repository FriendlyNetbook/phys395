program a0; implicit none

  integer, parameter :: n = 50 !Number of Fibonnacci Numbers to Print
  real :: x = 1
  real :: lp0 = 1, lp1 = x, temp !Starting 
  integer i

  write (*,*) fn0
  write (*,*) fn1

  do i=3, n
    temp = fn0
    fn0 = fn1
    fn1 = fn1 + temp
    write (*,*) fn1
  end do

end program
