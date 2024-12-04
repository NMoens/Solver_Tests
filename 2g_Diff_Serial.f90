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

end subroutine LinOp

!> Solves Lu = f for u
subroutine jacobi(u, f, N, dt, dx, Diff, j_it)
use mod_output
  integer, intent(in) :: N
  double precision, intent(in) :: f(N,N)
  double precision, intent(in) :: dt, dx, Diff
  integer, intent(in) :: j_it

  double precision :: Lu(N,N)
  double precision :: u(N,N), u_old(N,N)
  double precision :: dtau, err

  integer :: ii,jj,it

  !> Max stable ftcs timestep
  dtau = 0.8 * dx**2/(4*Diff)
  err = 1.d99
  it  = 1

  !> iteration pseudotime
  do while (err > 1.d-3 .and. it < j_it)
    !> Update u
    u_old = u
    call LinOp(Lu, u, N, dt, dx, Diff)
    u(2:N-1,2:N-1) = u(2:N-1,2:N-1) + dtau*(Lu(2:N-1,2:N-1) - f(2:N-1,2:N-1))

    !> calc error
    err = maxval(abs((u-u_old)/u_old))
    print*, it, err
    it = it+1 
  
  end do


end subroutine jacobi

subroutine two_cycle(rho, N, dt, dx, Diff)
use mod_output

  integer, intent(in) :: N
  double precision, intent(inout) :: rho(N,N)
  double precision, intent(in) :: dt, dx, Diff

  !> Fine grid
  double precision :: u_old(N,N), u(N,N), d(N,N), f(N,N), Lu(N,N)

  !> Coarse grid
  double precision :: u2(N/2,N/2), d2(N/2,N/2)

  integer :: ii,jj,it

  !> init guess for u
  u = rho

  !> Right hand side
  f = -rho/dt

  !> pre-smoothing 
  !> The main goal is to smoothen the solution a little, to prevent big errors on the coarsening step
  u_old = u
  call jacobi(u, f, N, dt, dx, Diff, 4)

  !> Calculate defect
  call LinOp(Lu, u, N, dt, dx, Diff)
  d(2:N-1,2:N-1) = Lu(2:N-1,2:N-1) - f(2:N-1,2:N-1)

  !> Coarsen, restrict
  do ii = 1,N/2
  do jj = 1,N/2
  d2(ii,jj) = (d(2*ii-1,2*jj-1) + d(2*ii,2*jj-1) + d(2*ii-1,2*jj) + d(2*ii,2*jj))/4
  end do
  end do

  !> Main relaxation
  !> Solve Lu = -d
  u2 = d2 !> Save input args
  call jacobi(u2, -d2, N/2, dt, 2*dx, Diff, 5000)

  !> Prolong
  do ii = 1,N/2
  do jj = 1,N/2
  u(1+(ii-1)*2,1+(jj-1)*2) = u2(ii,jj)
  u(1+(ii-1)*2+1,1+(jj-1)*2) = u2(ii,jj)
  u(1+(ii-1)*2,1+(jj-1)*2+1) = u2(ii,jj)
  u(1+(ii-1)*2+1,1+(jj-1)*2+1) = u2(ii,jj)
  !> rho has not changed
  end do
  end do

  !> Correct
  u(2:N-1,2:N-1) = u_old(2:N-1,2:N-1) + u(2:N-1,2:N-1)

  !> post-smoothing
  call jacobi(u, f, N, dt, dx, Diff, 4)

  !> Update rho
  rho(2:N-1,2:N-1) = u(2:N-1,2:N-1)

end subroutine two_cycle


end module mod_Relaxation


program Explicit_Diffusion_Serial
use mod_output
use mod_relaxation

!> Fortran shenanigans
implicit none

!> Declare variables
integer, parameter :: N = 100
double precision :: x(N), y(N)
double precision :: rho(N,N), rho_old(N,N), Diff

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
    rho(ii,jj) = dexp(-(x(ii)**2 + y(jj)**2))
  end do
end do

!> Time integration
time = 0.d0
it = 0
dt = 10 * dx**2/(4*Diff)

call cpu_time(t2)
do while (time < tmax)

  !> Execute jacobi solver
  call two_cycle(rho, N, dt, dx, Diff)

  !> Update time
  time = time + dt
  it = it + 1

end do
call cpu_time(t3)

!> Write output
call write_2D_slice(rho,N,'test')

print*, t1,t2,t3
print*, 'Setup took ', t2 - t1, ' s'
print*, 'Integration took :', t3 - t2, ' s' 

end program Explicit_Diffusion_Serial

