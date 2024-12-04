program Explicit_Diffusion_Serial
use mod_output

!> Fortran shenanigans
implicit none

!> Declare variables
integer, parameter :: N = 100
double precision :: x(N), y(N)
double precision :: rho(N,N), rho_old(N,N), Diff

double precision, parameter :: tmax = 1.d-1
double precision :: time
integer :: it
double precision :: t1, t2, t3

integer :: ii, jj
double precision :: dx, dt

!> Start timer
call cpu_time(t1)

!> Setup grid
dx = 10.d0/(N-1)

x(1) = -5.d0
y(1) = -5.d0
do ii = 2,N
  x(ii) = x(ii-1) + dx
  y(ii) = y(ii-1) + dx
end do

!> Initial conditions
Diff = 1.d0
do ii = 1, N
  do jj = 1, N
    rho(ii,jj) = dexp(-(x(ii)**2 + y(jj)**2))
  end do
end do

!> Time integration
time = 0.d0
it = 0
dt = 0.8 * dx**2/(4*Diff)

call cpu_time(t2)
do while (time < tmax)

  !> Save old state
  rho_old = rho

  !> Loop over grid
  do ii = 2, N-1
    do jj = 2, N-1
      rho(ii,jj) = rho_old(ii,jj) &
                 + dt*diff/dx**2 &
                 *(rho_old(ii-1,jj) + rho_old(ii+1,jj) + rho_old(ii,jj-1) + rho_old(ii,jj+1) - 4*rho_old(ii,jj))   
    end do
  end do

  !> Update time
  time = time + dt
  it = it + 1

end do
call cpu_time(t3)

print*, t1,t2,t3
print*, 'Setup took ', t2 - t1, ' s'
print*, 'Integration took :', t3 - t2, ' s' 

!> Write output
call write_2D_slice(rho,N,'Serial_out')

end program Explicit_Diffusion_Serial
