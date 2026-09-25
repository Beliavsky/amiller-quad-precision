module quad_special_functions
   use quadruple_precision
   implicit none
   private
   public :: q_erf_eval, q_lngm_eval
contains

   pure subroutine q_erf_eval(x, erf_out, erfc_out)
      type(quad), intent(in) :: x
      type(quad), intent(out), optional :: erf_out
      type(quad), intent(out), optional :: erfc_out
      real(qp) :: xv
      xv = qval(x)
      if (present(erf_out)) erf_out = from_qp(erf(xv))
      if (present(erfc_out)) erfc_out = from_qp(erfc(xv))
   end subroutine q_erf_eval

   pure function q_lngm_eval(x) result(b)
      type(quad), intent(in) :: x
      type(quad) :: b
      real(qp) :: xv
      xv = qval(x)
      if (xv > 0.0_qp) then
         b = from_qp(log_gamma(xv))
      else
         b = from_qp(huge(1.0_qp))
      end if
   end function q_lngm_eval

end module quad_special_functions

subroutine q_erf(x, erf, erfc)
   use quadruple_precision, only : quad
   use quad_special_functions, only : q_erf_eval
   implicit none
   type(quad), intent(in) :: x
   type(quad), intent(out), optional :: erf
   type(quad), intent(out), optional :: erfc
   call q_erf_eval(x, erf, erfc)
end subroutine q_erf

function q_lngm(x) result(b)
   use quadruple_precision, only : quad
   use quad_special_functions, only : q_lngm_eval
   implicit none
   type(quad), intent(in) :: x
   type(quad) :: b
   b = q_lngm_eval(x)
end function q_lngm
