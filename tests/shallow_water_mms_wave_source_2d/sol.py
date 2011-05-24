# This file was *autogenerated* from the file sol.sage.
from sage.all_cmdline import *   # import sage library
_sage_const_2 = Integer(2); _sage_const_1 = Integer(1); _sage_const_0 = Integer(0)
t = var('t')
x = var('x')
y = var('y')
k = var('k')
g = var('g') # gravity
H = var('H') # depth
a = _sage_const_0 
b = _sage_const_1 

eta = -sin(_sage_const_2 *pi*x+t)
u_x = sin(_sage_const_2 *pi*x+t)
u_y = cos(_sage_const_2 *pi*y+t)

print "u_x: ", u_x
print "u_y: ", u_y
print "eta: ", eta
print "u_x source: ", g*diff(eta, x) + diff(u_x, t)
print "u_y source: ", g*diff(eta, y) + diff(u_y, t)
print "cont source: ", diff(eta, t) + diff(H*u_x, x) + diff(H*u_y, y)
