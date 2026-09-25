module quadruple_precision
   use iso_fortran_env, only : int64, real32, real64, real128
   implicit none

   integer, parameter :: dp = real64
   integer, parameter :: qp = real128

   type :: quad
      real(dp) :: hi = 0.0_dp
      real(dp) :: lo = 0.0_dp
   end type quad

   interface operator(+)
      module procedure longadd
      module procedure quad_add_int, int_add_quad
      module procedure quad_add_real, real_add_quad
      module procedure quad_add_dp, dp_add_quad
   end interface
   interface operator(-)
      module procedure longsub
      module procedure quad_sub_int, int_sub_quad
      module procedure quad_sub_real, real_sub_quad
      module procedure quad_sub_dp, dp_sub_quad
      module procedure negate_quad
   end interface
   interface operator(*)
      module procedure longmul
      module procedure mult_quad_int, mult_int_quad
      module procedure mult_quad_real, mult_real_quad
      module procedure mult_quad_dp, mult_dp_quad
   end interface
   interface operator(/)
      module procedure longdiv
      module procedure div_quad_int, div_int_quad
      module procedure div_quad_real, div_real_quad
      module procedure div_quad_dp, div_dp_quad
   end interface
   interface operator(**)
      module procedure quad_pow_int
      module procedure quad_pow_real
      module procedure quad_pow_dp
      module procedure quad_pow_quad
   end interface
   interface operator(<)
      module procedure quad_lt
   end interface
   interface operator(<=)
      module procedure quad_le
   end interface
   interface operator(==)
      module procedure quad_eq
   end interface
   interface operator(/=)
      module procedure quad_ne
   end interface
   interface operator(>=)
      module procedure quad_ge
   end interface
   interface operator(>)
      module procedure quad_gt
   end interface
   interface assignment(=)
      module procedure quad_eq_int
      module procedure quad_eq_real
      module procedure quad_eq_dp
      module procedure int_eq_quad
      module procedure real_eq_quad
      module procedure dp_eq_quad
   end interface
   interface abs
      module procedure qabs
   end interface
   interface sqrt
      module procedure longsqrt
   end interface
   interface log
      module procedure longlog
   end interface
   interface exp
      module procedure longexp
   end interface
   interface sin
      module procedure longsin
   end interface
   interface cos
      module procedure longcos
   end interface
   interface tan
      module procedure longtan
   end interface
   interface asin
      module procedure longasin
   end interface
   interface acos
      module procedure longacos
   end interface
   interface atan
      module procedure longatan
   end interface
   interface atan2
      module procedure qatan2
   end interface
   interface scale
      module procedure qscale
   end interface
   interface epsilon
      module procedure q_epsilon
   end interface
   interface sum
      module procedure quad_sum
   end interface
   interface dot_product
      module procedure quad_dot_product
   end interface
   interface matmul
      module procedure q_matmul12
      module procedure q_matmul21
      module procedure q_matmul22
   end interface

   ! Alan Miller's original split constants, retained for source compatibility.
   type(quad), parameter :: &
      pi = quad(0.3141592653589793e1_dp, 0.1224646799147353e-15_dp), &
      piby2 = quad(0.1570796326794897e1_dp, -0.3828568698926950e-15_dp), &
      piby3 = quad(0.1047197551196598e1_dp, -0.3292527815701405e-15_dp), &
      piby4 = quad(0.7853981633974484e0_dp, -0.8040613248383182e-16_dp), &
      piby6 = quad(0.5235987755982990e0_dp, -0.1646263907850702e-15_dp), &
      twopi = quad(0.6283185307179586e1_dp, 0.2449293598294707e-15_dp), &
      ln_pi = quad(0.1144729885849400e1_dp, 0.2323105560877391e-15_dp), &
      sqrtpi = quad(0.1772453850905516e1_dp, -0.7666586499825800e-16_dp), &
      fact_pt5 = quad(0.8862269254527582e0_dp, -0.1493552349616447e-15_dp), &
      sqrt2pi = quad(0.2506628274631001e1_dp, -0.6273750096546544e-15_dp), &
      lnsqrt2pi = quad(0.9189385332046728e0_dp, -0.3878294158067242e-16_dp), &
      one_on2pi = quad(0.1591549430918953e0_dp, 0.4567181289366658e-16_dp), &
      two_on_rtpi = quad(0.1128379167095513e1_dp, -0.4287537502368968e-15_dp), &
      deg2rad = quad(0.1745329251994330e-1_dp, -0.3174581724866598e-17_dp), &
      rad2deg = quad(0.5729577951308232e2_dp, -0.1987849567057628e-14_dp), &
      ln2 = quad(0.6931471805599454e0_dp, -0.8783183432405266e-16_dp), &
      ln10 = quad(0.2302585092994046e1_dp, -0.2170756223382249e-15_dp), &
      log2e = quad(0.1442695040888964e1_dp, -0.6457785410341630e-15_dp), &
      log10e = quad(0.4342944819032519e0_dp, -0.5552037773430574e-16_dp), &
      log2_10 = quad(0.3321928094887362e1_dp, 0.1661617516973592e-15_dp), &
      log10_2 = quad(0.3010299956639812e0_dp, -0.2803728127785171e-17_dp), &
      euler = quad(0.5772156649015330e0_dp, -0.1159652176149463e-15_dp), &
      e = quad(0.2718281828459045e1_dp, 0.1445646891729250e-15_dp), &
      sqrt2 = quad(0.1414213562373095e1_dp, 0.1253716717905022e-15_dp), &
      sqrt3 = quad(0.1732050807568877e1_dp, 0.3223954471431004e-15_dp), &
      sqrt10 = quad(0.3162277660168380e1_dp, -0.6348773795572286e-15_dp)

contains

   pure function qval(a) result(x)
      type(quad), intent(in) :: a
      real(qp) :: x
      x = real(a%hi, qp) + real(a%lo, qp)
   end function qval

   pure function from_qp(x) result(a)
      real(qp), intent(in) :: x
      type(quad) :: a
      real(qp) :: h
      a%hi = real(x, dp)
      h = real(a%hi, qp)
      a%lo = real(x - h, dp)
   end function from_qp

   pure function exactmul2(a, c) result(ac)
      real(dp), intent(in) :: a, c
      type(quad) :: ac
      ac = from_qp(real(a, qp) * real(c, qp))
   end function exactmul2

   pure function longadd(a, c) result(ac)
      type(quad), intent(in) :: a, c
      type(quad) :: ac
      ac = from_qp(qval(a) + qval(c))
   end function longadd

   pure function longsub(a, c) result(ac)
      type(quad), intent(in) :: a, c
      type(quad) :: ac
      ac = from_qp(qval(a) - qval(c))
   end function longsub

   pure function longmul(a, c) result(ac)
      type(quad), intent(in) :: a, c
      type(quad) :: ac
      ac = from_qp(qval(a) * qval(c))
   end function longmul

   pure function longdiv(a, c) result(ac)
      type(quad), intent(in) :: a, c
      type(quad) :: ac
      ac = from_qp(qval(a) / qval(c))
   end function longdiv

   pure function quad_add_int(a, b) result(c)
      type(quad), intent(in) :: a
      integer, intent(in) :: b
      type(quad) :: c
      c = from_qp(qval(a) + real(b, qp))
   end function quad_add_int

   pure function int_add_quad(a, b) result(c)
      integer, intent(in) :: a
      type(quad), intent(in) :: b
      type(quad) :: c
      c = b + a
   end function int_add_quad

   pure function quad_add_real(a, b) result(c)
      type(quad), intent(in) :: a
      real(real32), intent(in) :: b
      type(quad) :: c
      c = from_qp(qval(a) + real(b, qp))
   end function quad_add_real

   pure function real_add_quad(a, b) result(c)
      real(real32), intent(in) :: a
      type(quad), intent(in) :: b
      type(quad) :: c
      c = b + a
   end function real_add_quad

   pure function quad_add_dp(a, b) result(c)
      type(quad), intent(in) :: a
      real(dp), intent(in) :: b
      type(quad) :: c
      c = from_qp(qval(a) + real(b, qp))
   end function quad_add_dp

   pure function dp_add_quad(a, b) result(c)
      real(dp), intent(in) :: a
      type(quad), intent(in) :: b
      type(quad) :: c
      c = b + a
   end function dp_add_quad

   pure function quad_sub_int(a, b) result(c)
      type(quad), intent(in) :: a
      integer, intent(in) :: b
      type(quad) :: c
      c = from_qp(qval(a) - real(b, qp))
   end function quad_sub_int

   pure function int_sub_quad(a, b) result(c)
      integer, intent(in) :: a
      type(quad), intent(in) :: b
      type(quad) :: c
      c = from_qp(real(a, qp) - qval(b))
   end function int_sub_quad

   pure function quad_sub_real(a, b) result(c)
      type(quad), intent(in) :: a
      real(real32), intent(in) :: b
      type(quad) :: c
      c = from_qp(qval(a) - real(b, qp))
   end function quad_sub_real

   pure function real_sub_quad(a, b) result(c)
      real(real32), intent(in) :: a
      type(quad), intent(in) :: b
      type(quad) :: c
      c = from_qp(real(a, qp) - qval(b))
   end function real_sub_quad

   pure function quad_sub_dp(a, b) result(c)
      type(quad), intent(in) :: a
      real(dp), intent(in) :: b
      type(quad) :: c
      c = from_qp(qval(a) - real(b, qp))
   end function quad_sub_dp

   pure function dp_sub_quad(a, b) result(c)
      real(dp), intent(in) :: a
      type(quad), intent(in) :: b
      type(quad) :: c
      c = from_qp(real(a, qp) - qval(b))
   end function dp_sub_quad

   pure function negate_quad(a) result(c)
      type(quad), intent(in) :: a
      type(quad) :: c
      c%hi = -a%hi
      c%lo = -a%lo
   end function negate_quad

   pure function mult_quad_int(a, b) result(c)
      type(quad), intent(in) :: a
      integer, intent(in) :: b
      type(quad) :: c
      c = from_qp(qval(a) * real(b, qp))
   end function mult_quad_int

   pure function mult_int_quad(a, b) result(c)
      integer, intent(in) :: a
      type(quad), intent(in) :: b
      type(quad) :: c
      c = b * a
   end function mult_int_quad

   pure function mult_quad_real(a, b) result(c)
      type(quad), intent(in) :: a
      real(real32), intent(in) :: b
      type(quad) :: c
      c = from_qp(qval(a) * real(b, qp))
   end function mult_quad_real

   pure function mult_real_quad(a, b) result(c)
      real(real32), intent(in) :: a
      type(quad), intent(in) :: b
      type(quad) :: c
      c = b * a
   end function mult_real_quad

   pure function mult_quad_dp(a, b) result(c)
      type(quad), intent(in) :: a
      real(dp), intent(in) :: b
      type(quad) :: c
      c = from_qp(qval(a) * real(b, qp))
   end function mult_quad_dp

   pure function mult_dp_quad(a, b) result(c)
      real(dp), intent(in) :: a
      type(quad), intent(in) :: b
      type(quad) :: c
      c = b * a
   end function mult_dp_quad

   pure function div_quad_int(a, b) result(c)
      type(quad), intent(in) :: a
      integer, intent(in) :: b
      type(quad) :: c
      c = from_qp(qval(a) / real(b, qp))
   end function div_quad_int

   pure function div_int_quad(a, b) result(c)
      integer, intent(in) :: a
      type(quad), intent(in) :: b
      type(quad) :: c
      c = from_qp(real(a, qp) / qval(b))
   end function div_int_quad

   pure function div_quad_real(a, b) result(c)
      type(quad), intent(in) :: a
      real(real32), intent(in) :: b
      type(quad) :: c
      c = from_qp(qval(a) / real(b, qp))
   end function div_quad_real

   pure function div_real_quad(a, b) result(c)
      real(real32), intent(in) :: a
      type(quad), intent(in) :: b
      type(quad) :: c
      c = from_qp(real(a, qp) / qval(b))
   end function div_real_quad

   pure function div_quad_dp(a, b) result(c)
      type(quad), intent(in) :: a
      real(dp), intent(in) :: b
      type(quad) :: c
      c = from_qp(qval(a) / real(b, qp))
   end function div_quad_dp

   pure function div_dp_quad(a, b) result(c)
      real(dp), intent(in) :: a
      type(quad), intent(in) :: b
      type(quad) :: c
      c = from_qp(real(a, qp) / qval(b))
   end function div_dp_quad

   pure function quad_pow_int(a, n) result(c)
      type(quad), intent(in) :: a
      integer, intent(in) :: n
      type(quad) :: c
      c = from_qp(qval(a) ** n)
   end function quad_pow_int

   pure function quad_pow_real(a, x) result(c)
      type(quad), intent(in) :: a
      real(real32), intent(in) :: x
      type(quad) :: c
      c = from_qp(qval(a) ** real(x, qp))
   end function quad_pow_real

   pure function quad_pow_dp(a, x) result(c)
      type(quad), intent(in) :: a
      real(dp), intent(in) :: x
      type(quad) :: c
      c = from_qp(qval(a) ** real(x, qp))
   end function quad_pow_dp

   pure function quad_pow_quad(a, x) result(c)
      type(quad), intent(in) :: a, x
      type(quad) :: c
      c = from_qp(qval(a) ** qval(x))
   end function quad_pow_quad

   pure logical function quad_lt(a, b)
      type(quad), intent(in) :: a, b
      quad_lt = qval(a) < qval(b)
   end function quad_lt

   pure logical function quad_le(a, b)
      type(quad), intent(in) :: a, b
      quad_le = qval(a) <= qval(b)
   end function quad_le

   pure logical function quad_eq(a, b)
      type(quad), intent(in) :: a, b
      quad_eq = transfer(a%hi, 0_int64) == transfer(b%hi, 0_int64) .and. &
         transfer(a%lo, 0_int64) == transfer(b%lo, 0_int64)
   end function quad_eq

   pure logical function quad_ne(a, b)
      type(quad), intent(in) :: a, b
      quad_ne = .not. quad_eq(a, b)
   end function quad_ne

   pure logical function quad_ge(a, b)
      type(quad), intent(in) :: a, b
      quad_ge = qval(a) >= qval(b)
   end function quad_ge

   pure logical function quad_gt(a, b)
      type(quad), intent(in) :: a, b
      quad_gt = qval(a) > qval(b)
   end function quad_gt

   pure function qabs(a) result(c)
      type(quad), intent(in) :: a
      type(quad) :: c
      c = from_qp(abs(qval(a)))
   end function qabs

   pure function longsqrt(a) result(c)
      type(quad), intent(in) :: a
      type(quad) :: c
      c = from_qp(sqrt(qval(a)))
   end function longsqrt

   pure function longlog(a) result(c)
      type(quad), intent(in) :: a
      type(quad) :: c
      c = from_qp(log(qval(a)))
   end function longlog

   pure function longexp(a) result(c)
      type(quad), intent(in) :: a
      type(quad) :: c
      c = from_qp(exp(qval(a)))
   end function longexp

   pure function longsin(a) result(c)
      type(quad), intent(in) :: a
      type(quad) :: c
      c = from_qp(sin(qval(a)))
   end function longsin

   pure function longcos(a) result(c)
      type(quad), intent(in) :: a
      type(quad) :: c
      c = from_qp(cos(qval(a)))
   end function longcos

   pure function longtan(a) result(c)
      type(quad), intent(in) :: a
      type(quad) :: c
      c = from_qp(tan(qval(a)))
   end function longtan

   pure function longasin(a) result(c)
      type(quad), intent(in) :: a
      type(quad) :: c
      c = from_qp(asin(qval(a)))
   end function longasin

   pure function longacos(a) result(c)
      type(quad), intent(in) :: a
      type(quad) :: c
      c = from_qp(acos(qval(a)))
   end function longacos

   pure function longatan(a) result(c)
      type(quad), intent(in) :: a
      type(quad) :: c
      c = from_qp(atan(qval(a)))
   end function longatan

   pure function qatan2(y, x) result(c)
      type(quad), intent(in) :: y, x
      type(quad) :: c
      c = from_qp(atan2(qval(y), qval(x)))
   end function qatan2

   pure function qscale(a, i) result(c)
      type(quad), intent(in) :: a
      integer, intent(in) :: i
      type(quad) :: c
      c = from_qp(scale(qval(a), i))
   end function qscale

   pure function q_epsilon(a) result(c)
      type(quad), intent(in) :: a
      type(quad) :: c
      c = from_qp(2.0_qp ** (-104))
      if (storage_size(a) == 0) c = a
   end function q_epsilon

   pure function quad_sum(a) result(s)
      type(quad), intent(in) :: a(:)
      type(quad) :: s
      integer :: i
      real(qp) :: total
      total = 0.0_qp
      do i = 1, size(a)
         total = total + qval(a(i))
      end do
      s = from_qp(total)
   end function quad_sum

   pure function quad_dot_product(a, b) result(s)
      type(quad), intent(in) :: a(:), b(:)
      type(quad) :: s
      integer :: i
      real(qp) :: total
      total = 0.0_qp
      do i = 1, min(size(a), size(b))
         total = total + qval(a(i)) * qval(b(i))
      end do
      s = from_qp(total)
   end function quad_dot_product

   pure function q_matmul12(a, b) result(c)
      type(quad), intent(in) :: a(:)
      type(quad), intent(in) :: b(:, :)
      type(quad) :: c(size(b, 2))
      integer :: i, j
      real(qp) :: total
      do j = 1, size(b, 2)
         total = 0.0_qp
         do i = 1, min(size(a), size(b, 1))
            total = total + qval(a(i)) * qval(b(i, j))
         end do
         c(j) = from_qp(total)
      end do
   end function q_matmul12

   pure function q_matmul21(a, b) result(c)
      type(quad), intent(in) :: a(:, :)
      type(quad), intent(in) :: b(:)
      type(quad) :: c(size(a, 1))
      integer :: i, j
      real(qp) :: total
      do i = 1, size(a, 1)
         total = 0.0_qp
         do j = 1, min(size(a, 2), size(b))
            total = total + qval(a(i, j)) * qval(b(j))
         end do
         c(i) = from_qp(total)
      end do
   end function q_matmul21

   pure function q_matmul22(a, b) result(c)
      type(quad), intent(in) :: a(:, :), b(:, :)
      type(quad) :: c(size(a, 1), size(b, 2))
      integer :: i, j, k
      real(qp) :: total
      do j = 1, size(b, 2)
         do i = 1, size(a, 1)
            total = 0.0_qp
            do k = 1, min(size(a, 2), size(b, 1))
               total = total + qval(a(i, k)) * qval(b(k, j))
            end do
            c(i, j) = from_qp(total)
         end do
      end do
   end function q_matmul22

   elemental subroutine quad_eq_int(lhs, rhs)
      type(quad), intent(out) :: lhs
      integer, intent(in) :: rhs
      lhs = from_qp(real(rhs, qp))
   end subroutine quad_eq_int

   elemental subroutine quad_eq_real(lhs, rhs)
      type(quad), intent(out) :: lhs
      real(real32), intent(in) :: rhs
      lhs = from_qp(real(rhs, qp))
   end subroutine quad_eq_real

   elemental subroutine quad_eq_dp(lhs, rhs)
      type(quad), intent(out) :: lhs
      real(dp), intent(in) :: rhs
      lhs = from_qp(real(rhs, qp))
   end subroutine quad_eq_dp

   elemental subroutine int_eq_quad(lhs, rhs)
      integer, intent(out) :: lhs
      type(quad), intent(in) :: rhs
      lhs = int(qval(rhs))
   end subroutine int_eq_quad

   elemental subroutine real_eq_quad(lhs, rhs)
      real(real32), intent(out) :: lhs
      type(quad), intent(in) :: rhs
      lhs = real(qval(rhs), real32)
   end subroutine real_eq_quad

   elemental subroutine dp_eq_quad(lhs, rhs)
      real(dp), intent(out) :: lhs
      type(quad), intent(in) :: rhs
      lhs = real(qval(rhs), dp)
   end subroutine dp_eq_quad

   function str_quad(text) result(a)
      character(len=*), intent(in) :: text
      type(quad) :: a
      real(qp) :: x
      read(text, *) x
      a = from_qp(x)
   end function str_quad

   function quad_str(a) result(text)
      type(quad), intent(in) :: a
      character(len=48) :: text
      write(text, '(es40.32e3)') qval(a)
   end function quad_str

end module quadruple_precision
