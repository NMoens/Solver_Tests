module mod_output

implicit none

contains

subroutine write_2D_slice(rho,N,fname)
integer, intent(in) :: N
double precision, intent(in) :: rho(N,N)
character(len=*), intent(in) :: fname
integer :: ii

open(unit=1,file=fname)

do ii = 1,N
  write(1,*) rho(ii,:)
end do

close(1)

        
end subroutine write_2D_slice


end module mod_output
