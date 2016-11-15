program a0; implicit none

  integer, parameter :: n = 50 !Number of Fibonnacci Numbers to Print
  integer*8 :: fn0 = 0, fn1 = 1, temp !Starting F1 and F2, alternate convention is fn0=1, fn1=1
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
