from math import exp

print "BC_x_w(1) = 0.0"
print "BC_x_e(1) = 200.0"
print "BC_y_s(1) = 0.0"
print "BC_y_n(1) = 0.0"
print "BC_TYPE(1) = 'MASS_INFLOW'"
print "BC_P_g(1) = 101325.0"
print "\n"

print "BC_x_w(2) = 0.0"
print "BC_x_e(2) = 7000.0"
print "BC_y_s(2) = 7000.0"
print "BC_y_n(2) = 7000.0"
print "BC_TYPE(2) = 'P_OUTFLOW'"
print "BC_P_g(2) = 44906.587482293"
print "\n"

n = 3
for bc in range(0, 7000,100):
   print "BC_x_e(%d) = %d" % (n, 7000)
   print "BC_x_w(%d) = %d" % (n, 7000)
   print "BC_y_n(%d) = %d" % (n, bc+100)
   print "BC_y_s(%d) = %d" % (n, bc)
   print "BC_TYPE(%d) = 'P_OUTFLOW'" % (n)
   print "BC_P_g(%d) = %f" % (n, 101325.0*exp(-(9.8*28.42*(bc+100))/(8314.56*288.15)))
   print "\n"
   n = n+1
