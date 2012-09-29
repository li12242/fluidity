# This file was *autogenerated* from the file mphase_mms_rans_p2p1_keps.sage.
from sage.all_cmdline import *   # import sage library
_sage_const_2 = Integer(2); _sage_const_1 = Integer(1); _sage_const_10p0 = RealNumber('10.0'); _sage_const_3p = RealNumber('3.'); _sage_const_8p2 = RealNumber('8.2'); _sage_const_1p8 = RealNumber('1.8'); _sage_const_1p0 = RealNumber('1.0'); _sage_const_1p3 = RealNumber('1.3'); _sage_const_1p7 = RealNumber('1.7'); _sage_const_0p0 = RealNumber('0.0'); _sage_const_3p8 = RealNumber('3.8'); _sage_const_0p4 = RealNumber('0.4'); _sage_const_0p7 = RealNumber('0.7'); _sage_const_0p6 = RealNumber('0.6'); _sage_const_0p9 = RealNumber('0.9'); _sage_const_0p8 = RealNumber('0.8'); _sage_const_3p7 = RealNumber('3.7'); _sage_const_2p1 = RealNumber('2.1'); _sage_const_2p5 = RealNumber('2.5'); _sage_const_5p2 = RealNumber('5.2'); _sage_const_2p = RealNumber('2.'); _sage_const_4p3 = RealNumber('4.3')
y = var('y')

def function(phi_0, phi_x, phi_y, phi_xy, 
             f_sin_x, f_cos_x, f_sin_y, f_cos_y, f_sin_xy, f_cos_xy, 
             alpha_x, alpha_y, alpha_xy):
    
    f_0 = phi_0 
    f_x = phi_x*(f_sin_x*sin(alpha_x*x) + f_cos_x*cos(alpha_x*x)) 
    f_y = phi_y*(f_sin_y*sin(alpha_y*y) + f_cos_y*cos(alpha_y*y)) 
    f_xy = phi_xy*(f_sin_xy*sin(alpha_xy*x*y/pi) + f_cos_xy*cos(alpha_xy*x*y/pi)) 
    f = f_0 + f_x + f_y + f_xy
    return f

p = function(-_sage_const_1p0 , _sage_const_1p0 , _sage_const_1p0 , _sage_const_1p0 ,
             _sage_const_1p0 , _sage_const_0p0 , _sage_const_0p0 , _sage_const_1p0 , _sage_const_1p0 , _sage_const_0p0 ,
             _sage_const_1p0 , _sage_const_1p0 , _sage_const_1p0 )
       
#### Phase 1
u1 = function(_sage_const_2p5 , _sage_const_1p0 , _sage_const_0p6 , _sage_const_0p0 , 
             _sage_const_0p0 , _sage_const_1p0 , _sage_const_1p0 , _sage_const_0p0 , _sage_const_1p0 , _sage_const_0p0 ,
             _sage_const_1p0 , _sage_const_1p0 , _sage_const_1p0 )
v1 = integral(-diff(u1,x),y)  # divergence free
temperature1 = function(_sage_const_5p2 , -_sage_const_1p8 , -_sage_const_1p3 , _sage_const_3p7 , 
               _sage_const_1p0 , _sage_const_0p0 , _sage_const_0p0 , _sage_const_1p0 , _sage_const_1p0 , _sage_const_0p0 ,
               _sage_const_1p7 , _sage_const_2p1 , _sage_const_1p3 )
ke1 = function(_sage_const_0p9 , _sage_const_0p9 , _sage_const_0p6 , _sage_const_0p4 , 
             _sage_const_0p0 , _sage_const_1p0 , _sage_const_1p0 , _sage_const_0p0 , _sage_const_0p0 , _sage_const_1p0 ,
             _sage_const_0p6 , _sage_const_0p7 , _sage_const_0p8 )
eps1 = function(_sage_const_8p2 , -_sage_const_3p8 , _sage_const_4p3 , _sage_const_1p7 , 
             _sage_const_1p0 , _sage_const_0p0 , _sage_const_0p0 , _sage_const_1p0 , _sage_const_1p0 , _sage_const_0p0 ,
             _sage_const_0p7 , _sage_const_0p8 , _sage_const_0p6 )
rho1 = _sage_const_10p0 *temperature1 + _sage_const_10p0 
nu_T1 = rho1*(ke1**_sage_const_2 )/eps1
nu1 = _sage_const_0p7 
vfrac1 = _sage_const_0p7 

tau_xx1 = _sage_const_2 *nu1*diff(u1,x)            
tau_xy1 = nu1*(diff(u1,y) + diff(v1,x))
tau_yy1 = _sage_const_2 *nu1*diff(v1,y)            
tau_yx1 = nu1*(diff(u1,y) + diff(v1,x))  
tau_xx_R1 = _sage_const_2 *nu_T1*diff(u1,x)
tau_xy_R1 = nu_T1*(diff(u1,y) + diff(v1,x))
tau_yy_R1 = _sage_const_2 *nu_T1*diff(v1,y)
tau_yx_R1 = nu_T1*(diff(u1,y) + diff(v1,x))

Su1 = vfrac1*rho1*u1*diff(u1,x) + vfrac1*rho1*v1*diff(u1,y) - diff(vfrac1*tau_xx1, x) - diff(vfrac1*tau_xy1, y) - diff(vfrac1*tau_xx_R1, x) - diff(vfrac1*tau_xy_R1, y) - vfrac1*rho1 + vfrac1*diff(p,x)  
Sv1 = vfrac1*rho1*u1*diff(v1,x) + vfrac1*rho1*v1*diff(v1,y) - diff(vfrac1*tau_yx1, x) - diff(vfrac1*tau_yy1, y) - diff(vfrac1*tau_yx_R1, x) - diff(vfrac1*tau_yy_R1, y) - vfrac1*rho1 + vfrac1*diff(p,y)  

Stemperature1 = u1*diff(temperature1,x) + v1*diff(temperature1,y) - (_sage_const_1p0  + nu_T1)*(diff(temperature1, x, x) + diff(temperature1, y, y)) - diff(nu_T1, x)*diff(temperature1, x) -  diff(nu_T1, y)*diff(temperature1, y)

P1 = nu_T1*(_sage_const_2 *(diff(u1,x)**_sage_const_2  + diff(v1,y)**_sage_const_2  + diff(u1,y)*diff(v1,x)) + diff(u1,y)**_sage_const_2  + diff(v1,x)**_sage_const_2 ) - (_sage_const_2p /_sage_const_3p )*ke1*rho1*(diff(u1,x) + diff(v1,y))

C3 = _sage_const_1p0 
g_x = _sage_const_1p0 
g_y = _sage_const_1p0 
B1 = -nu_T1*(g_x*diff(rho1,x) + g_y*diff(rho1,y))/rho1

pr = _sage_const_1 
ab = _sage_const_1 
bo = _sage_const_1 

Ske1 = vfrac1*rho1*u1*diff(ke1,x) + vfrac1*rho1*v1*diff(ke1,y) - vfrac1*(nu1 + nu_T1)*(diff(ke1, x, x) + diff(ke1, y, y)) - diff(vfrac1*nu_T1, x)*diff(ke1, x) -  diff(vfrac1*nu_T1, y)*diff(ke1, y) - pr*vfrac1*P1 + ab*vfrac1*rho1*eps1 - bo*vfrac1*B1
Seps1 = vfrac1*rho1*u1*diff(eps1,x) + vfrac1*rho1*v1*diff(eps1,y) - vfrac1*(nu1 + nu_T1)*(diff(eps1, x, x) + diff(eps1, y, y)) - diff(vfrac1*nu_T1, x)*diff(eps1, x) -  diff(vfrac1*nu_T1, y)*diff(eps1, y) - pr*vfrac1*(eps1/ke1)*P1 + ab*vfrac1*rho1*(eps1**_sage_const_2 /ke1) - bo*vfrac1*C3*(eps1/ke1)*B1
  
  
##### Phase 2 #####
u2 = function(_sage_const_2p5 , _sage_const_1p0 , _sage_const_0p6 , _sage_const_0p0 , 
             _sage_const_0p0 , _sage_const_1p0 , _sage_const_1p0 , _sage_const_0p0 , _sage_const_1p0 , _sage_const_0p0 ,
             _sage_const_1p0 , _sage_const_1p0 , _sage_const_1p0 )
v2 = integral(-diff(u2,x),y)  # divergence free

temperature2 = function(_sage_const_5p2 , -_sage_const_1p8 , -_sage_const_1p3 , _sage_const_3p7 , 
               _sage_const_1p0 , _sage_const_0p0 , _sage_const_0p0 , _sage_const_1p0 , _sage_const_1p0 , _sage_const_0p0 ,
               _sage_const_1p7 , _sage_const_2p1 , _sage_const_1p3 )
ke2 = function(_sage_const_0p9 , _sage_const_0p9 , _sage_const_0p6 , _sage_const_0p4 , 
             _sage_const_0p0 , _sage_const_1p0 , _sage_const_1p0 , _sage_const_0p0 , _sage_const_0p0 , _sage_const_1p0 ,
             _sage_const_0p6 , _sage_const_0p7 , _sage_const_0p8 )
eps2 = function(_sage_const_8p2 , -_sage_const_3p8 , _sage_const_4p3 , _sage_const_1p7 , 
             _sage_const_1p0 , _sage_const_0p0 , _sage_const_0p0 , _sage_const_1p0 , _sage_const_1p0 , _sage_const_0p0 ,
             _sage_const_0p7 , _sage_const_0p8 , _sage_const_0p6 )
rho2 = _sage_const_10p0 *temperature2 + _sage_const_10p0 
nu_T2 = rho2*(ke2**_sage_const_2 )/eps2
nu2 = _sage_const_0p7 
vfrac2 = _sage_const_1p0  - vfrac1

tau_xx2 = _sage_const_2 *nu2*diff(u2,x)            
tau_xy2 = nu2*(diff(u2,y) + diff(v2,x))
tau_yy2 = _sage_const_2 *nu2*diff(v2,y)            
tau_yx2 = nu2*(diff(u2,y) + diff(v2,x))  
tau_xx_R2 = _sage_const_2 *nu_T2*diff(u2,x)
tau_xy_R2 = nu_T2*(diff(u2,y) + diff(v2,x))
tau_yy_R2 = _sage_const_2 *nu_T2*diff(v2,y)
tau_yx_R2 = nu_T2*(diff(u2,y) + diff(v2,x))

Su2 = vfrac2*rho2*u2*diff(u2,x) + vfrac2*rho2*v2*diff(u2,y) - diff(vfrac2*tau_xx2, x) - diff(vfrac2*tau_xy2, y) - diff(vfrac2*tau_xx_R2, x) - diff(vfrac2*tau_xy_R2, y) - vfrac2*rho2 + vfrac2*diff(p,x)  
Sv2 = vfrac2*rho2*u2*diff(v2,x) + vfrac2*rho2*v2*diff(v2,y) - diff(vfrac2*tau_yx2, x) - diff(vfrac2*tau_yy2, y) - diff(vfrac2*tau_yx_R2, x) - diff(vfrac2*tau_yy_R2, y) - vfrac2*rho2 + vfrac2*diff(p,y)  

Stemperature2 = u2*diff(temperature2,x) + v2*diff(temperature2,y) - (_sage_const_1p0  + nu_T2)*(diff(temperature2, x, x) + diff(temperature2, y, y)) - diff(nu_T2, x)*diff(temperature2, x) -  diff(nu_T2, y)*diff(temperature2, y)

P2 = nu_T2*(_sage_const_2 *(diff(u2,x)**_sage_const_2  + diff(v2,y)**_sage_const_2  + diff(u2,y)*diff(v2,x)) + diff(u2,y)**_sage_const_2  + diff(v2,x)**_sage_const_2 ) - (_sage_const_2p /_sage_const_3p )*ke2*rho2*(diff(u2,x) + diff(v2,y))

C3 = _sage_const_1p0 
g_x = _sage_const_1p0 
g_y = _sage_const_1p0 
B2 = -nu_T2*(g_x*diff(rho2,x) + g_y*diff(rho2,y))/rho2

pr = _sage_const_1 
ab = _sage_const_1 
bo = _sage_const_1 

Ske2 = vfrac2*rho2*u2*diff(ke2,x) + vfrac2*rho2*v2*diff(ke2,y) - vfrac2*(nu2 + nu_T2)*(diff(ke2, x, x) + diff(ke2, y, y)) - diff(vfrac2*nu_T2, x)*diff(ke2, x) -  diff(vfrac2*nu_T2, y)*diff(ke2, y) - pr*vfrac2*P2 + ab*vfrac2*rho2*eps2 - bo*vfrac2*B2
Seps2 = vfrac2*rho2*u2*diff(eps2,x) + vfrac2*rho2*v2*diff(eps2,y) - vfrac2*(nu2 + nu_T2)*(diff(eps2, x, x) + diff(eps2, y, y)) - diff(vfrac2*nu_T2, x)*diff(eps2, x) -  diff(vfrac2*nu_T2, y)*diff(eps2, y) - pr*vfrac2*(eps2/ke2)*P2 + ab*vfrac2*rho2*(eps2**_sage_const_2 /ke2) - bo*vfrac2*C3*(eps2/ke2)*B2

print 'from math import sin, cos, tanh, pi'
print ''
print 'def u1(X):'
print '    return', str(u1).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def v1(X):'
print '    return', str(v1).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print '' 
print 'def u2(X):'
print '    return', str(u2).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def v2(X):'
print '    return', str(v2).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print '' 
print 'def p(X):'
print '    return', str(p).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''  
print 'def temperature1(X):'
print '    return', str(temperature1).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''  
print 'def temperature2(X):'
print '    return', str(temperature2).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''  
print 'def rho1(X):'
print '    return', str(rho1).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''  
print 'def rho2(X):'
print '    return', str(rho2).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''  
print 'def ke1(X):'
print '    return', str(ke1).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def eps1(X):'
print '    return', str(eps1).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def ke2(X):'
print '    return', str(ke2).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def eps2(X):'
print '    return', str(eps2).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_u1(X):'
print '    return', str(Su1).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_v1(X):'
print '    return', str(Sv1).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_u2(X):'
print '    return', str(Su2).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_v2(X):'
print '    return', str(Sv2).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_temperature1(X):'
print '    return', str(Stemperature1).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_temperature2(X):'
print '    return', str(Stemperature2).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_ke1(X):'
print '    return', str(Ske1).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_eps1(X):'
print '    return', str(Seps1).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_ke2(X):'
print '    return', str(Ske2).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def forcing_eps2(X):'
print '    return', str(Seps2).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def P1(X):'
print '    return', str(P1).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def P2(X):'
print '    return', str(P2).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def P_eps1(X):'
print '    return', str((eps1/ke1)*P1).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def P_eps2(X):'
print '    return', str((eps2/ke2)*P2).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def A1(X):'
print '    return', str(rho1*eps1/ke1).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def A2(X):'
print '    return', str(rho2*eps2/ke2).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def B1(X):'
print '    return', str(B1).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def B2(X):'
print '    return', str(B2).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def B_eps1(X):'
print '    return', str(C3*(eps1/ke1)*B1).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def B_eps2(X):'
print '    return', str(C3*(eps2/ke2)*B2).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def EV1(X):'
print '    return', str(nu_T1).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def EV2(X):'
print '    return', str(nu_T2).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'X[0]').replace('y', 'X[1]')
print ''
print 'def velocity1(X):'
print '   return [u1(X), v1(X)]'
print ''
print 'def forcing_velocity1(X):'
print '   return [forcing_u1(X), forcing_v1(X)]'
print ''
print 'def velocity2(X):'
print '   return [u2(X), v2(X)]'
print ''
print 'def forcing_velocity2(X):'
print '   return [forcing_u2(X), forcing_v2(X)]'
print ''
print 'def A_ke1(X):'
print '   return A1([X[0],X[1]])*ke1([X[0],X[1]])'
print ''
print 'def A_eps1(X):'
print '   return A1([X[0],X[1]])*eps1([X[0],X[1]])'
print ''
print 'def A_ke2(X):'
print '   return A2([X[0],X[1]])*ke2([X[0],X[1]])'
print ''
print 'def A_eps2(X):'
print '   return A2([X[0],X[1]])*eps2([X[0],X[1]])'
