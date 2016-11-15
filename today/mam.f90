program mam; implicit none
  real, parameter :: dt = 0.1, tDomain(2) = [0,10], xDomain(2) = [-10,10]
  integer, parameter :: N = 2**10, Nt = (tDomain(2)-tDomain(1))/dt
  integer, parameter :: fields = 2
  integer l
  !real :: t(Nt) = (/( (tDomain(1)+i*dt),i=1,Nt)/)
  real samples(fields,N,3)

  call init(samples(:,:,1), samples(:,:,2))

  ! Evaluate loop reusing slices in a ring buffer
  do l=1,Nt
     if (mod(l-1,8)==0 ) then
        call dump(samples(:,:,2))
     end if
     call step(samples(:,:,1), samples(:,:,2), samples(:,:,3))
     call step(samples(:,:,2), samples(:,:,3), samples(:,:,1))
     call step(samples(:,:,3), samples(:,:,1), samples(:,:,2))
  end do

contains
  subroutine step()
    
  end subroutine step
  subroutine init(dn, hr)
    real, dimension(fields,n) :: dn, hr
    intent(in) dn; intent(out) hr

    dn=0.0; up = 0.0
  end subroutine init
  subroutine dump(hr)
    real,dimension(fieldsNt) :: hr
    integer i
    intent(in) hr
    do i=1,N
       write(*,*) t,x,hr(:,i*Nt/N)
    end do
    
  end subroutine dump
end program mam
