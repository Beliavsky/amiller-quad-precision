program test_legacy_entry_points
   use quadruple_precision
   implicit none
   type(quad) :: x, ev, ec, lg
   real(qp), parameter :: tol = 2.0e-30_qp

   interface
      subroutine q_erf(x, erf, erfc)
         use quadruple_precision, only : quad
         type(quad), intent(in) :: x
         type(quad), intent(out), optional :: erf
         type(quad), intent(out), optional :: erfc
      end subroutine q_erf
      function q_lngm(x) result(b)
         use quadruple_precision, only : quad
         type(quad), intent(in) :: x
         type(quad) :: b
      end function q_lngm
   end interface

   x = 1.0_dp
   call q_erf(x, ev, ec)
   call require(abs(qval(ev) - erf(1.0_qp)) < tol, 'external q_erf')
   call require(abs(qval(ec) - erfc(1.0_qp)) < tol, 'external q_erf erfc')

   x = 5.0_dp
   lg = q_lngm(x)
   call require(abs(qval(lg) - log(24.0_qp)) < tol, 'external q_lngm')

   print '(a)', 'test_legacy_entry_points: PASS'
contains
   subroutine require(ok, label)
      logical, intent(in) :: ok
      character(len=*), intent(in) :: label
      if (.not. ok) then
         print '(a,1x,a)', 'FAIL:', label
         error stop 1
      end if
   end subroutine require
end program test_legacy_entry_points
