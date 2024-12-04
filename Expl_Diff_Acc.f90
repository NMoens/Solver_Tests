program Explicit_Diffusion_Serial
use openacc

!> Fortran shenanigans
implicit none

!> Declare variables
integer, parameter :: N = 16384
double precision :: x(N), y(N)
double precision :: rho(N,N), rho_old(N,N), Diff

double precision, parameter :: tmax = 1.d0
double precision :: time

double precision :: t1, t2, t3

integer :: ii, jj, it
integer, parameter :: it_max = 1024 !131072

double precision :: dx, dt

!> Check/set number of threads/gpu device stuff
print*, 'Nr of devices ', acc_get_num_devices(acc_get_device_type()), 'device type ', acc_get_device_type()

!> Start timer
call cpu_time(t1)

!> Setup grid
dx = 2.d0/(N-1)

x(1) = -1.d0
y(1) = -1.d0
do ii = 2,N
  x(ii) = x(ii-1) + dx
  y(ii) = y(ii-1) + dx
end do

!> Initial conditions
Diff = 1.d0
do ii = 1, N
  do jj = 1, N
    rho(ii,jj) = dexp(-(x(ii)**2 + y(ii)**2)/4.d0)
  end do
end do

!> Time integration
time = 0.d0
it = 0
dt = 0.8 * dx**2/(4*Diff)

call cpu_time(t2)
!$acc data copy(rho, time) copyin(dt, dx, diff) create(rho_old)
do while (time < tmax .and. it < it_max)

  !> Save old state
  !$acc parallel loop collapse(2)
  do ii = 1,N
    do jj = 1,N
      rho_old(ii,jj) = rho(ii,jj)
    end do
  end do
  !$acc end parallel loop

  !> Loop over grid
  !$acc parallel loop collapse(2)
  do ii = 2, N-1
    do jj = 2, N-1
      rho(ii,jj) = rho_old(ii,jj) &
                 + dt*diff/dx**2 &
                 *(rho_old(ii-1,jj) + rho_old(ii+1,jj) + rho_old(ii,jj-1) + rho_old(ii,jj+1) - 4*rho_old(ii,jj))
    end do
  end do
  !$acc end parallel loop
  time = time + dt
  it = it + 1

end do
!$acc end data
call cpu_time(t3)

print*, t1,t2,t3
print*, 'Setup took ', t2 - t1, ' s'
print*, 'Integration took :', t3 - t2, ' s'

end program Explicit_Diffusion_Serial
                                            
