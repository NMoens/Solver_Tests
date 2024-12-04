program Explicit_Diffusion_Serial
use DiffSolver


!> Fortran shenanigans
implicit none

!> Declare variables
integer, parameter :: N = 256
double precision :: x(N), y(N)
double precision :: rho(N,N), rho_old(N,N), Diff

double precision, parameter :: tmax = 1.d0
double precision :: time

double precision :: t1, t2, t3

integer :: ii, jj
double precision :: dx, dt

!> Device declaration (use _d naming convention)
double precision, device :: rho_d(N,N), rho_old_d(N,N)

!> grid and tBlock variables
type(dim3) :: grid, tBlock

tBlock = dim3(N,N,1)
grid = dim3(ceiling(real(N)/tBlock%x),ceiling(real(N)/tBlock%y),1)


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
dt = 0.8 * dx**2/(4*Diff)

!> transfer host data to device
rho_d = rho
rho_old_d = rho

call cpu_time(t2)

do while (time < tmax)

  call FTCS<<<grid,tBlock>>>(rho_d,rho_old_d,dt*diff/dx**2 )
  time = time + dt

end do
call cpu_time(t3)

!> transfer device data to host
rho = rho_d


print*, t1,t2,t3
print*, 'Setup took ', t2 - t1, ' s'
print*, 'Integration took :', t3 - t2, ' s' 

end program Explicit_Diffusion_Serial


module DiffSolver
  implicit none
contains
  
  attributes(global) subroutine FTCS(rho, alpha)

  double precision :: rho(:,:), rho_old(:,:)
  double precision, value :: alpha

  integer :: ii, jj

  !> Get cell index from thread id
  ii = blockDim%x * (blockIdx%x - 1) + threadIdx%x
  jj = blockDim%y * (blockIdx%y - 1) + threadIdx%y

  !> Save old state
  rho_old(ii,jj) = rho(ii,jj)

  !> Cell update  
  rho(ii,jj) = rho_old(ii,jj) &
                 + alpha &
                 *(rho_old(ii-1,jj) + rho_old(ii+1,jj) + rho_old(ii,jj-1) + rho_old(ii,jj+1) - 4*rho_old(ii,jj))
  
  end subroutine FTCS

end module DiffSolver
