program test_complex
   use quadruple_precision
   use quad_prec_complex
   implicit none
   type(qc) :: z, w
   type(quad) :: zero, one
   real(qp), parameter :: tol = 2.0e-30_qp

   zero = 0.0_dp
   one = 1.0_dp
   z = qc(-one, zero)
   w = sqrt(z)
   call require(abs(qval(w%qr)) < tol, 'sqrt(-1) real')
   call require(abs(abs(qval(w%qi)) - 1.0_qp) < tol, 'sqrt(-1) imag')

   z = qc(one, one)
   w = exp(log(z)) - z
   call require(abs(qval(w%qr)) < 1.0e-29_qp, 'complex exp/log real')
   call require(abs(qval(w%qi)) < 1.0e-29_qp, 'complex exp/log imag')

   print '(a)', 'test_complex: PASS'
contains
   subroutine require(ok, label)
      logical, intent(in) :: ok
      character(len=*), intent(in) :: label
      if (.not. ok) then
         print '(a,1x,a)', 'FAIL:', label
         error stop 1
      end if
   end subroutine require
end program test_complex
