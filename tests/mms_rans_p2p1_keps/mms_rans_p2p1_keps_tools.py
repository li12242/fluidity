from math import sin, cos, tanh, pi

def u(X):
    return 0.600*sin(X[1]) + cos(X[0]) + 2.50

def v(X):
    return X[1]*sin(X[0])

def p(X):
    return sin(X[0]*X[1]/pi) + sin(X[0]) + cos(X[1]) - 1.00

def rho(X):
    return 3.70*sin(1.30*X[0]*X[1]/pi) - 1.80*sin(1.70*X[0]) - 1.30*cos(2.10*X[1]) + 5.20

def ke(X):
    return 0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900

def eps(X):
    return 1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20

def forcing_u(X):
    return 0.600*X[1]*sin(X[0])*cos(X[1]) - (-0.600*sin(X[1]) + cos(X[0]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) - (0.600*sin(X[1]) + cos(X[0]) + 2.50)*sin(X[0]) - 2*(1.02*X[1]*cos(0.600*X[0]*X[1]/pi)/pi - 2.66*cos(0.700*X[0]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*sin(X[0])/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)**2 + 4*(-0.320*X[1]*sin(0.800*X[0]*X[1]/pi)/pi - 0.540*sin(0.600*X[0]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)*sin(X[0])/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) + (1.02*X[0]*cos(0.600*X[0]*X[1]/pi)/pi - 3.44*sin(0.800*X[1]))*(X[1]*cos(X[0]) + 0.600*cos(X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)**2 - 2*(-0.320*X[0]*sin(0.800*X[0]*X[1]/pi)/pi + 0.420*cos(0.700*X[1]))*(X[1]*cos(X[0]) + 0.600*cos(X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) + 2*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*cos(X[0])/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) - 0.213333333333333*X[1]*sin(0.800*X[0]*X[1]/pi)/pi + X[1]*cos(X[0]*X[1]/pi)/pi - 3.70*sin(1.30*X[0]*X[1]/pi) - 0.360*sin(0.600*X[0]) + 1.80*sin(1.70*X[0]) + 0.600*sin(X[1]) + 2.00*cos(X[0]) + 1.30*cos(2.10*X[1]) - 5.20

def forcing_v(X):
    return (0.600*sin(X[1]) + cos(X[0]) + 2.50)*X[1]*cos(X[0]) + (0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*X[1]*sin(X[0])/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) + X[1]*sin(X[0])**2 + (1.02*X[1]*cos(0.600*X[0]*X[1]/pi)/pi - 2.66*cos(0.700*X[0]))*(X[1]*cos(X[0]) + 0.600*cos(X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)**2 - 2*(-0.320*X[1]*sin(0.800*X[0]*X[1]/pi)/pi - 0.540*sin(0.600*X[0]))*(X[1]*cos(X[0]) + 0.600*cos(X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) + 2*(1.02*X[0]*cos(0.600*X[0]*X[1]/pi)/pi - 3.44*sin(0.800*X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*sin(X[0])/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)**2 - 4*(-0.320*X[0]*sin(0.800*X[0]*X[1]/pi)/pi + 0.420*cos(0.700*X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)*sin(X[0])/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) + X[1]*sin(X[0]) - 0.213333333333333*X[0]*sin(0.800*X[0]*X[1]/pi)/pi + X[0]*cos(X[0]*X[1]/pi)/pi - 3.70*sin(1.30*X[0]*X[1]/pi) + 1.80*sin(1.70*X[0]) - sin(X[1]) + 0.280*cos(0.700*X[1]) + 1.30*cos(2.10*X[1]) - 5.20

def forcing_rho(X):
    return (4.81*X[0]*cos(1.30*X[0]*X[1]/pi)/pi + 2.73*sin(2.10*X[1]))*X[1]*sin(X[0]) - ((0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) + 1.00)*(-6.25300000000000*X[0]**2*sin(1.30*X[0]*X[1]/pi)/pi**2 - 6.25300000000000*X[1]**2*sin(1.30*X[0]*X[1]/pi)/pi**2 + 5.20200000000000*sin(1.70*X[0]) + 5.73300000000000*cos(2.10*X[1])) + (0.600*sin(X[1]) + cos(X[0]) + 2.50)*(4.81*X[1]*cos(1.30*X[0]*X[1]/pi)/pi - 3.06*cos(1.70*X[0])) + (4.81*X[1]*cos(1.30*X[0]*X[1]/pi)/pi - 3.06*cos(1.70*X[0]))*((1.02*X[1]*cos(0.600*X[0]*X[1]/pi)/pi - 2.66*cos(0.700*X[0]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)**2 - 2*(-0.320*X[1]*sin(0.800*X[0]*X[1]/pi)/pi - 0.540*sin(0.600*X[0]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)) + (4.81*X[0]*cos(1.30*X[0]*X[1]/pi)/pi + 2.73*sin(2.10*X[1]))*((1.02*X[0]*cos(0.600*X[0]*X[1]/pi)/pi - 3.44*sin(0.800*X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)**2 - 2*(-0.320*X[0]*sin(0.800*X[0]*X[1]/pi)/pi + 0.420*cos(0.700*X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20))

def forcing_ke(X):
    return (-0.320*X[0]*sin(0.800*X[0]*X[1]/pi)/pi + 0.420*cos(0.700*X[1]))*X[1]*sin(X[0]) + (0.600*sin(X[1]) + cos(X[0]) + 2.50)*(-0.320*X[1]*sin(0.800*X[0]*X[1]/pi)/pi - 0.540*sin(0.600*X[0])) + (-0.320*X[1]*sin(0.800*X[0]*X[1]/pi)/pi - 0.540*sin(0.600*X[0]))*((1.02*X[1]*cos(0.600*X[0]*X[1]/pi)/pi - 2.66*cos(0.700*X[0]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)**2 - 2*(-0.320*X[1]*sin(0.800*X[0]*X[1]/pi)/pi - 0.540*sin(0.600*X[0]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)) + (-0.320*X[0]*sin(0.800*X[0]*X[1]/pi)/pi + 0.420*cos(0.700*X[1]))*((1.02*X[0]*cos(0.600*X[0]*X[1]/pi)/pi - 3.44*sin(0.800*X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)**2 - 2*(-0.320*X[0]*sin(0.800*X[0]*X[1]/pi)/pi + 0.420*cos(0.700*X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)) - (0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(-0.256*X[0]**2*cos(0.800*X[0]*X[1]/pi)/pi**2 - 0.256*X[1]**2*cos(0.800*X[0]*X[1]/pi)/pi**2 - 0.294*sin(0.700*X[1]) - 0.324*cos(0.600*X[0]))/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) + (0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(4.81*X[0]*cos(1.30*X[0]*X[1]/pi)/pi + 4.81*X[1]*cos(1.30*X[0]*X[1]/pi)/pi + 2.73*sin(2.10*X[1]) - 3.06*cos(1.70*X[0]))/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) - (0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(X[1]**2*cos(X[0])**2 + 1.20*X[1]*cos(X[0])*cos(X[1]) + 4*sin(X[0])**2 + 0.360*cos(X[1])**2)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) + 1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20

def forcing_eps(X):
    return (1.02*X[0]*cos(0.600*X[0]*X[1]/pi)/pi - 3.44*sin(0.800*X[1]))*X[1]*sin(X[0]) + (0.600*sin(X[1]) + cos(X[0]) + 2.50)*(1.02*X[1]*cos(0.600*X[0]*X[1]/pi)/pi - 2.66*cos(0.700*X[0])) + (1.02*X[1]*cos(0.600*X[0]*X[1]/pi)/pi - 2.66*cos(0.700*X[0]))*((1.02*X[1]*cos(0.600*X[0]*X[1]/pi)/pi - 2.66*cos(0.700*X[0]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)**2 - 2*(-0.320*X[1]*sin(0.800*X[0]*X[1]/pi)/pi - 0.540*sin(0.600*X[0]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)) + (1.02*X[0]*cos(0.600*X[0]*X[1]/pi)/pi - 3.44*sin(0.800*X[1]))*((1.02*X[0]*cos(0.600*X[0]*X[1]/pi)/pi - 3.44*sin(0.800*X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)**2 - 2*(-0.320*X[0]*sin(0.800*X[0]*X[1]/pi)/pi + 0.420*cos(0.700*X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)) - (0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(-0.612*X[0]**2*sin(0.600*X[0]*X[1]/pi)/pi**2 - 0.612*X[1]**2*sin(0.600*X[0]*X[1]/pi)/pi**2 + 1.86200000000000*sin(0.700*X[0]) - 2.75200000000000*cos(0.800*X[1]))/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) + (0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)*(4.81*X[0]*cos(1.30*X[0]*X[1]/pi)/pi + 4.81*X[1]*cos(1.30*X[0]*X[1]/pi)/pi + 2.73*sin(2.10*X[1]) - 3.06*cos(1.70*X[0])) - (0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)*(X[1]**2*cos(X[0])**2 + 1.20*X[1]*cos(X[0])*cos(X[1]) + 4*sin(X[0])**2 + 0.360*cos(X[1])**2) + (1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)**2/(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)

def P(X):
    return (0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(X[1]**2*cos(X[0])**2 + 1.20*X[1]*cos(X[0])*cos(X[1]) + 4*sin(X[0])**2 + 0.360*cos(X[1])**2)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)

def P_eps(X):
    return (0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)*(X[1]**2*cos(X[0])**2 + 1.20*X[1]*cos(X[0])*cos(X[1]) + 4*sin(X[0])**2 + 0.360*cos(X[1])**2)

def A(X):
    return (1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)/(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)

def B(X):
    return -(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(4.81*X[0]*cos(1.30*X[0]*X[1]/pi)/pi + 4.81*X[1]*cos(1.30*X[0]*X[1]/pi)/pi + 2.73*sin(2.10*X[1]) - 3.06*cos(1.70*X[0]))/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)

def B_eps(X):
    return -(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)*(4.81*X[0]*cos(1.30*X[0]*X[1]/pi)/pi + 4.81*X[1]*cos(1.30*X[0]*X[1]/pi)/pi + 2.73*sin(2.10*X[1]) - 3.06*cos(1.70*X[0]))

def EV(X):
    return (0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)

def velocity(X):
   return [u(X), v(X)]

def forcing_velocity(X):
   return [forcing_u(X), forcing_v(X)]

def A_ke(X):
   return A([X[0],X[1]])*ke([X[0],X[1]])

def A_eps(X):
   return A([X[0],X[1]])*eps([X[0],X[1]])
