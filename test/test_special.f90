program test_special
   use quadruple_precision
   use quad_special_functions
   implicit none
   type(quad) :: x, e1, e2, lg
   real(qp), parameter :: tol = 2.0e-30_qp

   x = 1.0_dp
   call q_erf_eval(x, e1, e2)
   call require(abs(qval(e1) - erf(1.0_qp)) < tol, 'erf')
   call require(abs(qval(e2) - erfc(1.0_qp)) < tol, 'erfc')

   x = 5.0_dp
   lg = q_lngm_eval(x)
   call require(abs(qval(lg) - log(24.0_qp)) < tol, 'log gamma')

   print '(a)', 'test_special: PASS'
contains
   subroutine require(ok, label)
      logical, intent(in) :: ok
      character(len=*), intent(in) :: label
      if (.not. ok) then
         print '(a,1x,a)', 'FAIL:', label
         error stop 1
      end if
   end subroutine require
end program test_special
