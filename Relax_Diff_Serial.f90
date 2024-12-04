module mod_relaxation

implicit none

contains

subroutine LinOp(Lu, u, N, dt, dx, Diff)

  integer, intent(in) :: N
  double precision, intent(out) :: Lu(N,N)
  double precision, intent(in) :: u(N,N)
  double precision, intent(in) :: dt, dx, Diff

  integer :: ii,jj

  !> Loop over grid
  Lu = 0.d0
  do ii = 2, N-1
    do jj = 2, N-1
      Lu(ii,jj) = diff/dx**2 *(u(ii-1,jj) + u(ii+1,jj) + u(ii,jj-1) + u(ii,jj+1) - 4*u(ii,jj)) - u(ii,jj)/dt
    end do
  end do

  !> BC?
  do jj = 2,N-1
    Lu(1,jj) = diff/dx**2 *(u(1,jj) + u(1+1,jj) + u(1,jj-1) + u(1,jj+1) - 4*u(1,jj)) - u(1,jj)/dt
    Lu(N,jj) = diff/dx**2 *(u(N-1,jj) + u(N,jj) + u(N,jj-1) + u(N,jj+1) - 4*u(N,jj)) - u(N,jj)/dt
    ii = jj
    Lu(ii,1) = diff/dx**2 *(u(ii-1,1) + u(ii+1,1) + u(ii,1) + u(ii,1+1) - 4*u(ii,1)) - u(ii,1)/dt
    Lu(ii,N) = diff/dx**2 *(u(ii-1,N) + u(ii+1,N) + u(ii,N-1) + u(ii,N) - 4*u(ii,N)) - u(ii,N)/dt
  end do

  !> Corners
  Lu(1,1) = diff/dx**2 *(u(1,jj) + u(1+1,1) + u(1,1) + u(1,1+1) - 4*u(1,1)) - u(1,1)/dt
  Lu(N,1) = diff/dx**2 *(u(N-1,1) + u(N,1) + u(N,1) + u(N,1+1) - 4*u(N,1)) - u(N,1)/dt
  Lu(1,N) = diff/dx**2 *(u(1,N) + u(1+1,N) + u(1,N-1) + u(1,N) - 4*u(1,N)) - u(1,N)/dt
  Lu(N,N) = diff/dx**2 *(u(N-1,N) + u(N,jj) + u(N,N-1) + u(N,N) - 4*u(N,N)) - u(N,N)/dt

end subroutine LinOp

subroutine jacobi(rho, N, dt, dx, Diff)
use mod_output
  integer, intent(in) :: N
  double precision, intent(inout) :: rho(N,N)
  double precision, intent(in) :: dt, dx, Diff

  double precision :: Lu(N,N), f(N,N)
  double precision :: u(N,N), u_old(N,N)
  double precision :: dtau, err

  integer :: ii,jj,it

  print*, 'timestep', dt

  !> Init guess for u
  u = 0.d0
  f = -rho/dt

  !> Max stable ftcs timestep
  dtau = 0.8d0 * dx**2/(4*Diff)
  err = 1.d99
  it  = 1

  !> iteration pseudotime
  do while (err > 1.d-8)
    !> Update u
    u_old = u
    call LinOp(Lu, u, N, dt, dx, Diff)
    u(2:N-1,2:N-1) = u(2:N-1,2:N-1) + dtau*(Lu(2:N-1,2:N-1) - f(2:N-1,2:N-1))

    !> calc error
    err = maxval(abs((u-u_old)/u_old))
    print*, it, err
    it = it+1 
  
  end do

  rho = u

end subroutine jacobi

end module mod_Relaxation


program Explicit_Diffusion_Serial
use mod_output
use mod_relaxation

!> Fortran shenanigans
implicit none

!> Declare variables
integer, parameter :: N = 16
double precision :: x(N), y(N)
double precision :: rho(N,N), rho_old(N,N), Diff

!> Analytical solution
double precision :: sol(N,N)

!> New jacobi variables
double precision :: u(N,N), u_old(N,N), dtau, err

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
    rho(ii,jj) = 1.d0/(4*3.14*Diff)*dexp(-(x(ii)**2 + y(jj)**2)/4*Diff)
  end do
end do

!> Time integration
time = 0.d0
it = 0
dt = 10.d0 * dx**2/(4*Diff)

call cpu_time(t2)
do while (time < tmax)

  !> Execute jacobi solver
  call jacobi(rho, N, dt, dx, Diff)

  !> Update time
  time = time + dt
  it = it + 1

  !> Analytical solution
  do ii = 1, N
    do jj = 1, N
      sol(ii,jj) = 1.d0/(4*3.14*Diff*(1+time)) * dexp(-(x(ii)**2 + y(jj)**2)/(4*Diff*(1+time)))
    end do
  end do
 
  print*, time

  !> Write output
  call write_2D_slice(rho,N,'test')
  call write_2D_slice(sol,N,'test2')
  stop

end do
call cpu_time(t3)

!> Write output
call write_2D_slice(rho-sol,N,'test')

print*, t1,t2,t3
print*, 'Setup took ', t2 - t1, ' s'
print*, 'Integration took :', t3 - t2, ' s' 

end program Explicit_Diffusion_Serial

