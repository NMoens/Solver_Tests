program Explicit_Diffusion_Serial

use OMP_lib

!> Fortran shenanigans
implicit none

!> Declare variables
integer, parameter :: N = 100
double precision :: x(N), y(N)
double precision :: rho(N,N), rho_old(N,N), Diff

double precision, parameter :: tmax = 1.d0
double precision :: time

double precision :: t1, t2, t3

integer :: ii, jj
double precision :: dx, dt

!> Check/set number of available threads:
call omp_set_num_threads( 16 )
print*, 'Running on: ', omp_get_max_threads(), ' threads'

!> Start timer
t1 = omp_get_wtime()

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
dt = 0.8 * dx**2/(4*Diff)

t2 = omp_get_wtime()
do while (time < tmax)

  !> Save old state
  rho_old = rho

  !> Loop over grid
  !$omp parallel do
  do ii = 2, N-1
    do jj = 2, N-1
      rho(ii,jj) = rho_old(ii,jj) &
                 + dt*diff/dx**2 &
                 *(rho_old(ii-1,jj) + rho_old(ii+1,jj) + rho_old(ii,jj-1) + rho_old(ii,jj+1) - 4*rho_old(ii,jj))   
    end do
  end do
  !$omp end parallel do

  time = time + dt

end do
t3 = omp_get_wtime()

print*, t1,t2,t3
print*, 'Setup took ', t2 - t1, ' s'
print*, 'Integration took :', t3 - t2, ' s' 

end program Explicit_Diffusion_Serial
