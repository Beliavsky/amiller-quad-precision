program test_arrays
   use quadruple_precision
   implicit none
   type(quad) :: a(2), b(2), m(2,2), v(2), s
   real(qp), parameter :: tol = 1.0e-30_qp

   a = [quad(1.0_dp, 0.0_dp), quad(2.0_dp, 0.0_dp)]
   b = [quad(3.0_dp, 0.0_dp), quad(4.0_dp, 0.0_dp)]
   s = dot_product(a, b)
   call require(abs(qval(s) - 11.0_qp) < tol, 'dot product')

   m(1,1) = 1.0_dp
   m(1,2) = 2.0_dp
   m(2,1) = 3.0_dp
   m(2,2) = 4.0_dp
   v = matmul(m, a)
   call require(abs(qval(v(1)) - 5.0_qp) < tol, 'matmul 1')
   call require(abs(qval(v(2)) - 11.0_qp) < tol, 'matmul 2')

   print '(a)', 'test_arrays: PASS'
contains
   subroutine require(ok, label)
      logical, intent(in) :: ok
      character(len=*), intent(in) :: label
      if (.not. ok) then
         print '(a,1x,a)', 'FAIL:', label
         error stop 1
      end if
   end subroutine require
end program test_arrays
