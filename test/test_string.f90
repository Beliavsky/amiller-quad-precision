program test_string
   use quadruple_precision
   implicit none
   type(quad) :: x, y
   character(len=48) :: text
   real(qp), parameter :: tol = 1.0e-30_qp

   x = str_quad('3.141592653589793238462643383279')
   text = quad_str(x)
   y = str_quad(text)
   call require(abs(qval(x) - qval(y)) < tol, 'string round trip')
   call require(abs(qval(x) - acos(-1.0_qp)) < 1.0e-29_qp, 'string parse')

   print '(a)', 'test_string: PASS'
contains
   subroutine require(ok, label)
      logical, intent(in) :: ok
      character(len=*), intent(in) :: label
      if (.not. ok) then
         print '(a,1x,a)', 'FAIL:', label
         error stop 1
      end if
   end subroutine require
end program test_string
