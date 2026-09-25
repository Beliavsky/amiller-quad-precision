program test_arithmetic
   use quadruple_precision
   implicit none
   type(quad) :: a, b, c
   real(qp), parameter :: tol = 1.0e-30_qp

   a = 1.0_dp
   a = a / 3
   b = a * 3
   c = b - 1
   call require(abs(qval(c)) < tol, 'division/multiplication')

   a = sqrt(quad(2.0_dp, 0.0_dp))
   c = a*a - 2
   call require(abs(qval(c)) < tol, 'sqrt')

   c = (quad(1.0_dp, 0.0_dp) + epsilon(a)) - quad(1.0_dp, 0.0_dp)
   call require(qval(c) > 0.0_qp, 'quad epsilon')

   print '(a)', 'test_arithmetic: PASS'
contains
   subroutine require(ok, label)
      logical, intent(in) :: ok
      character(len=*), intent(in) :: label
      if (.not. ok) then
         print '(a,1x,a)', 'FAIL:', label
         error stop 1
      end if
   end subroutine require
end program test_arithmetic
