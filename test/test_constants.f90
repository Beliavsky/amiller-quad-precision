program test_constants
   use quadruple_precision
   implicit none
   real(qp), parameter :: tol = 2.0e-30_qp

   call require(abs(qval(pi) - acos(-1.0_qp)) < tol, 'pi')
   call require(abs(qval(ln2) - log(2.0_qp)) < tol, 'ln2')
   call require(abs(qval(sqrt2) - sqrt(2.0_qp)) < tol, 'sqrt2')

   print '(a)', 'test_constants: PASS'
contains
   subroutine require(ok, label)
      logical, intent(in) :: ok
      character(len=*), intent(in) :: label
      if (.not. ok) then
         print '(a,1x,a)', 'FAIL:', label
         error stop 1
      end if
   end subroutine require
end program test_constants
