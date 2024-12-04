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

!> Solves Lu = f for u
subroutine jacobi(u, f, N, dt, dx, Diff, j_it)
use mod_output
  integer, intent(in) :: N
  double precision, intent(inout) :: u(N,N)
  double precision, intent(in) :: f(N,N)
  double precision, intent(in) :: dt, dx, Diff
  integer, intent(in) :: j_it

  double precision :: Lu(N,N)
  double precision :: u_old(N,N)
  double precision :: dtau, err

  double precision :: dtau_fac

  integer :: ii,jj,it

  print*, 'timestep', dt

  !> Max stable ftcs timestep
  dtau_fac = 0.5d0 !dt/(dx**2/(4*Diff)) !> coarser grids have trouble bc timestep is too small, adapt for this
  dtau = dtau_fac  * dx**2/(4*Diff)
  err = 1.d99
  it  = 1

  !> iteration pseudotime
  do while (err > 1.d-8 .and. it <= j_it)
    !> Update u
    u_old = u
    call LinOp(Lu, u, N, dt, dx, Diff)
    u = u + dtau*(Lu - f)

    !> calc error
    err = maxval(abs((u-u_old)/u))!/dtau_fac
    print*, it, err, N
    it = it+1 
 
    err = 1.d0

  end do

end subroutine jacobi

subroutine coarsen(u, u2, N)
  integer, intent(in) :: N
  double precision, intent(in) :: u(N,N)
  double precision, intent(out) :: u2(N/2,N/2)

  integer :: ii, jj

  !> Coarsen, restrict
  do ii = 1,N/2
  do jj = 1,N/2
    u2(ii,jj) = (u(2*ii-1,2*jj-1) + u(2*ii,2*jj-1) + u(2*ii-1,2*jj) + u(2*ii,2*jj))/4
  end do
  end do

end subroutine coarsen

subroutine prolong(u,u2,N)
  integer, intent(in) :: N
  double precision, intent(in) :: u2(N/2,N/2)
  double precision, intent(out) :: u(N,N)

  integer :: ii, jj

  do ii = 1,N/2
  do jj = 1,N/2
  u(1+(ii-1)*2,1+(jj-1)*2) = u2(ii,jj)
  u(1+(ii-1)*2+1,1+(jj-1)*2) = u2(ii,jj)
  u(1+(ii-1)*2,1+(jj-1)*2+1) = u2(ii,jj)
  u(1+(ii-1)*2+1,1+(jj-1)*2+1) = u2(ii,jj)
  end do
  end do

end subroutine prolong


recursive subroutine V_cycle(u, f, N, dt, dx, Diff)
use mod_output

  integer, intent(in) :: N
  double precision, intent(inout) :: u(N,N)
  double precision, intent(in) :: f(N,N)
  double precision, intent(in) :: dt, dx, Diff

  !> Fine grid
  double precision :: u_old(N,N), d(N,N), Lu(N,N)
  double precision :: u_bc(N,N), corr(N,N)

  !> Coarse grid
  double precision :: u2(N/2,N/2), d2(N/2,N/2), u2_old(N/2,N/2), f2(N/2,N/2), Lu2(N/2,N/2), corr2(N/2,N/2)

  integer :: ii,jj,it

  print*, 'start Vcycle', N

  !> Save old state
  u_old = u

  !> pre-smoothing 
  !> The main goal is to smoothen the solution a little, to prevent big errors on the coarsening step
  call jacobi(u, f, N, dt, dx, Diff, 2)

  !> Calculate defect
  call LinOp(Lu, u, N, dt, dx, Diff)
  d = f - Lu

  !> Coarsen, restrict
  call coarsen(u,u2,N)
  call coarsen(d,d2,N)
  u2_old = u2

  !> Update f on coarse grid
  call LinOp(Lu2, u2, N/2, dt, 2*dx, Diff)
  f2 = d2 + Lu2

  !> Go deeper or back
  if (N/2 > 16) then
    call V_cycle(u2, f2, N/2, dt, 2*dx, Diff)
  else
    !> Main relaxation
    print*, 'Exact solution', N/2

    !> Solve Lu = f
    call jacobi(u2, f2, N/2, dt, 2*dx, Diff, 15)

    print*, 'Exact solution done'
  end if

  !> Correction
  corr2 = u2 - u2_old

  !> Prolong
  call prolong(corr, corr2, N)

  !> Correct.
  u = u + corr

  !> post-smoothing
  call jacobi(u, f, N, dt, dx, Diff, 4)

  print*, 'end Vcycle', N

end subroutine V_cycle

!> Solves Lrho = rho
subroutine Multi_Grid(rho, N, dt, dx, Diff)
use mod_output

  integer, intent(in) :: N
  double precision, intent(inout) :: rho(N,N)
  double precision, intent(in) :: dt, dx, Diff

  !> Fine grid
  double precision :: u(N,N), d(N,N), f(N,N), Lu(N,N)

  double precision :: maxd

  integer :: ii,jj,it

  !> init guess for u
  u = rho

  !> Right hand side
  f = -rho/dt

  !> Repeat V-cycles until converged
  maxd = 1.d0
  do while (maxd > 1.d-3)
    !> Start V cycle
    call V_cycle(u, f, N, dt, dx, Diff)

    !> Calculate defect
    call LinOp(Lu, u, N, dt, dx, Diff)
    d = Lu - f

    print*, 'max defect', maxval(abs(d))
    maxd = maxval(abs(d))
  end do

  !> Update rho
  rho(2:N-1,2:N-1) = u(2:N-1,2:N-1)

end subroutine Multi_Grid


end module mod_Relaxation


program Explicit_Diffusion_Serial
use mod_output
use mod_relaxation

!> Fortran shenanigans
implicit none

!> Declare variables
integer, parameter :: N = 4096
double precision :: x(N), y(N)
double precision :: rho(N,N), rho_old(N,N), Diff

!> New jacobi variables
double precision :: u(N,N), u_old(N,N), dtau, err

!> Analytic solution
double precision :: sol(N,N)

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
    rho(ii,jj) = 1.d0/(4*3.14*Diff) * dexp(-(x(ii)**2 + y(jj)**2)/(4*Diff))
  end do
end do

!> Time integration
time = 0.d0
it = 0
dt = 1.d0 !1.d3 * dx**2/(4*Diff)

call cpu_time(t2)
do while (time < tmax)

  !> Execute jacobi solver
  call Multi_Grid(rho, N, dt, dx, Diff)

  !> Update time
  time = time + dt
  it = it + 1

  !> Analytical solution
  do ii = 1, N
    do jj = 1, N
      sol(ii,jj) = 1.d0/(4*3.14*Diff*(1+time)) * dexp(-(x(ii)**2 + y(jj)**2)/(4*Diff*(1+time)))
    end do
  end do

  call write_2D_slice(rho,N,'test')
  call write_2D_slice(sol,N,'test2')
  stop

end do
call cpu_time(t3)

!> Write output
call write_2D_slice(rho,N,'test')

print*, t1,t2,t3
print*, 'Setup took ', t2 - t1, ' s'
print*, 'Integration took :', t3 - t2, ' s' 

end program Explicit_Diffusion_Serial

