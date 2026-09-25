module quad_prec_complex
   use quadruple_precision
   implicit none

   type :: qc
      type(quad) :: qr
      type(quad) :: qi
   end type qc

   interface operator(+)
      module procedure qc_add
   end interface
   interface operator(-)
      module procedure qc_sub, negate_qc
   end interface
   interface operator(*)
      module procedure qc_mult
   end interface
   interface operator(/)
      module procedure qc_div
   end interface
   interface abs
      module procedure qc_abs
   end interface
   interface cmplx
      module procedure qc_cmplx
   end interface
   interface conjg
      module procedure qc_conjg
   end interface
   interface aimag
      module procedure qc_aimag
   end interface
   interface sqrt
      module procedure qc_sqrt
   end interface
   interface log
      module procedure qc_log
   end interface
   interface exp
      module procedure qc_exp
   end interface

contains

   pure function qc_add(x, y) result(z)
      type(qc), intent(in) :: x, y
      type(qc) :: z
      z%qr = x%qr + y%qr
      z%qi = x%qi + y%qi
   end function qc_add

   pure function qc_sub(x, y) result(z)
      type(qc), intent(in) :: x, y
      type(qc) :: z
      z%qr = x%qr - y%qr
      z%qi = x%qi - y%qi
   end function qc_sub

   pure function negate_qc(x) result(z)
      type(qc), intent(in) :: x
      type(qc) :: z
      z%qr = -x%qr
      z%qi = -x%qi
   end function negate_qc

   pure function qc_mult(x, y) result(z)
      type(qc), intent(in) :: x, y
      type(qc) :: z
      z%qr = x%qr * y%qr - x%qi * y%qi
      z%qi = x%qi * y%qr + x%qr * y%qi
   end function qc_mult

   pure function qc_div(x, y) result(z)
      type(qc), intent(in) :: x, y
      type(qc) :: z
      type(quad) :: den
      den = y%qr * y%qr + y%qi * y%qi
      z%qr = (x%qr * y%qr + x%qi * y%qi) / den
      z%qi = (x%qi * y%qr - x%qr * y%qi) / den
   end function qc_div

   pure function qc_cmplx(xr, xi) result(z)
      type(quad), intent(in) :: xr, xi
      type(qc) :: z
      z%qr = xr
      z%qi = xi
   end function qc_cmplx

   pure function qc_aimag(x) result(z)
      type(qc), intent(in) :: x
      type(quad) :: z
      z = x%qi
   end function qc_aimag

   pure function qc_conjg(x) result(z)
      type(qc), intent(in) :: x
      type(qc) :: z
      z%qr = x%qr
      z%qi = -x%qi
   end function qc_conjg

   pure function qc_abs(x) result(z)
      type(qc), intent(in) :: x
      type(quad) :: z
      z = sqrt(x%qr * x%qr + x%qi * x%qi)
   end function qc_abs

   pure function qc_sqrt(x) result(z)
      type(qc), intent(in) :: x
      type(qc) :: z
      type(quad) :: r, t
      r = abs(x)
      if (abs(r%hi) + abs(r%lo) < tiny(1.0_dp)) then
         z%qr = 0.0_dp
         z%qi = 0.0_dp
      else if (x%qr%hi >= 0.0_dp) then
         t = sqrt((r + x%qr) / 2.0_dp)
         z%qr = t
         z%qi = x%qi / (2.0_dp * t)
      else
         t = sqrt((r - x%qr) / 2.0_dp)
         if (x%qi%hi < 0.0_dp) t = -t
         z%qi = t
         z%qr = x%qi / (2.0_dp * t)
      end if
   end function qc_sqrt

   pure function qc_exp(x) result(z)
      type(qc), intent(in) :: x
      type(qc) :: z
      type(quad) :: er
      er = exp(x%qr)
      z%qr = er * cos(x%qi)
      z%qi = er * sin(x%qi)
   end function qc_exp

   pure function qc_log(x) result(z)
      type(qc), intent(in) :: x
      type(qc) :: z
      z%qr = log(abs(x))
      z%qi = atan2(x%qi, x%qr)
   end function qc_log

end module quad_prec_complex
