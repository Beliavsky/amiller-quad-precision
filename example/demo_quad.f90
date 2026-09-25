program demo_quad
   use quadruple_precision
   implicit none
   type(quad) :: x, y

   x = str_quad('1.000000000000000000000000000001')
   y = (x - 1.0_dp) * 1.0e30_dp

   print '(a)', 'Alan Miller quad-compatible arithmetic using a REAL128 backend'
   print '(a,a)', 'pi       = ', trim(quad_str(pi))
   print '(a,a)', 'scaled difference = ', trim(quad_str(y))
end program demo_quad
