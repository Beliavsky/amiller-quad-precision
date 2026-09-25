program test_transcendentals
   use quadruple_precision
   implicit none
   type(quad) :: x, y
   real(qp), parameter :: tol = 2.0e-30_qp

   x = 2.0_dp
   y = exp(log(x)) - x
   call require(abs(qval(y)) < tol, 'exp(log(x))')

   y = sin(piby6) - 0.5_dp
   call require(abs(qval(y)) < tol, 'sin(pi/6)')

   y = cos(piby3) - 0.5_dp
   call require(abs(qval(y)) < tol, 'cos(pi/3)')

   y = atan2(quad(1.0_dp, 0.0_dp), quad(1.0_dp, 0.0_dp)) - piby4
   call require(abs(qval(y)) < tol, 'atan2')

   print '(a)', 'test_transcendentals: PASS'
contains
   subroutine require(ok, label)
      logical, intent(in) :: ok
      character(len=*), intent(in) :: label
      if (.not. ok) then
         print '(a,1x,a)', 'FAIL:', label
         error stop 1
      end if
   end subroutine require
end program test_transcendentals
