from math import sin, cos, tanh, pi

def u1(X):
    return 0.600*sin(X[1]) + cos(X[0]) + 2.50

def v1(X):
    return X[1]*sin(X[0])

def u2(X):
    return 0.600*sin(X[1]) + cos(X[0]) + 2.50

def v2(X):
    return X[1]*sin(X[0])

def p(X):
    return sin(X[0]*X[1]/pi) + sin(X[0]) + cos(X[1]) - 1.00

def temperature1(X):
    return 3.70*sin(1.30*X[0]*X[1]/pi) - 1.80*sin(1.70*X[0]) - 1.30*cos(2.10*X[1]) + 5.20

def temperature2(X):
    return 3.70*sin(1.30*X[0]*X[1]/pi) - 1.80*sin(1.70*X[0]) - 1.30*cos(2.10*X[1]) + 5.20

def rho1(X):
    return 37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0

def rho2(X):
    return 37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0

def ke1(X):
    return 0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900

def eps1(X):
    return 1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20

def ke2(X):
    return 0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900

def eps2(X):
    return 1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20

def forcing_u1(X):
    return 0.600*(25.9*sin(1.30*X[0]*X[1]/pi) - 12.6*sin(1.70*X[0]) - 9.10*cos(2.10*X[1]) + 43.4)*X[1]*sin(X[0])*cos(X[1]) - 0.700*(-0.600*sin(X[1]) + cos(X[0]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) - (0.600*sin(X[1]) + cos(X[0]) + 2.50)*(25.9*sin(1.30*X[0]*X[1]/pi) - 12.6*sin(1.70*X[0]) - 9.10*cos(2.10*X[1]) + 43.4)*sin(X[0]) + 1.40*(48.1*X[1]*cos(1.30*X[0]*X[1]/pi)/pi - 30.6*cos(1.70*X[0]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*sin(X[0])/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) - 1.40*(1.02*X[1]*cos(0.600*X[0]*X[1]/pi)/pi - 2.66*cos(0.700*X[0]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)*sin(X[0])/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)**2 + 2.80*(-0.320*X[1]*sin(0.800*X[0]*X[1]/pi)/pi - 0.540*sin(0.600*X[0]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)*sin(X[0])/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) - 0.700*(48.1*X[0]*cos(1.30*X[0]*X[1]/pi)/pi + 27.3*sin(2.10*X[1]))*(X[1]*cos(X[0]) + 0.600*cos(X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) + 0.700*(1.02*X[0]*cos(0.600*X[0]*X[1]/pi)/pi - 3.44*sin(0.800*X[1]))*(X[1]*cos(X[0]) + 0.600*cos(X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)**2 - 1.40*(-0.320*X[0]*sin(0.800*X[0]*X[1]/pi)/pi + 0.420*cos(0.700*X[1]))*(X[1]*cos(X[0]) + 0.600*cos(X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) + 1.40*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)*cos(X[0])/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) + 0.666666666666667*(48.1*X[1]*cos(1.30*X[0]*X[1]/pi)/pi - 30.6*cos(1.70*X[0]))*(0.420*sin(0.700*X[1]) + 0.280*cos(0.800*X[0]*X[1]/pi) + 0.630*cos(0.600*X[0]) + 0.630) + 0.666666666666667*(-0.224*X[1]*sin(0.800*X[0]*X[1]/pi)/pi - 0.378*sin(0.600*X[0]))*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0) + 0.700*X[1]*cos(X[0]*X[1]/pi)/pi - 25.9*sin(1.30*X[0]*X[1]/pi) + 12.6*sin(1.70*X[0]) + 0.294*sin(X[1]) + 1.19*cos(X[0]) + 9.10*cos(2.10*X[1]) - 43.4

def forcing_v1(X):
    return (0.600*sin(X[1]) + cos(X[0]) + 2.50)*(25.9*sin(1.30*X[0]*X[1]/pi) - 12.6*sin(1.70*X[0]) - 9.10*cos(2.10*X[1]) + 43.4)*X[1]*cos(X[0]) + 0.700*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)*X[1]*sin(X[0])/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) + (25.9*sin(1.30*X[0]*X[1]/pi) - 12.6*sin(1.70*X[0]) - 9.10*cos(2.10*X[1]) + 43.4)*X[1]*sin(X[0])**2 - 0.700*(48.1*X[1]*cos(1.30*X[0]*X[1]/pi)/pi - 30.6*cos(1.70*X[0]))*(X[1]*cos(X[0]) + 0.600*cos(X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) + 0.700*(1.02*X[1]*cos(0.600*X[0]*X[1]/pi)/pi - 2.66*cos(0.700*X[0]))*(X[1]*cos(X[0]) + 0.600*cos(X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)**2 - 1.40*(-0.320*X[1]*sin(0.800*X[0]*X[1]/pi)/pi - 0.540*sin(0.600*X[0]))*(X[1]*cos(X[0]) + 0.600*cos(X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) - 1.40*(48.1*X[0]*cos(1.30*X[0]*X[1]/pi)/pi + 27.3*sin(2.10*X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*sin(X[0])/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) + 1.40*(1.02*X[0]*cos(0.600*X[0]*X[1]/pi)/pi - 3.44*sin(0.800*X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)*sin(X[0])/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)**2 - 2.80*(-0.320*X[0]*sin(0.800*X[0]*X[1]/pi)/pi + 0.420*cos(0.700*X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)*sin(X[0])/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) + 0.666666666666667*(48.1*X[0]*cos(1.30*X[0]*X[1]/pi)/pi + 27.3*sin(2.10*X[1]))*(0.420*sin(0.700*X[1]) + 0.280*cos(0.800*X[0]*X[1]/pi) + 0.630*cos(0.600*X[0]) + 0.630) + 0.666666666666667*(-0.224*X[0]*sin(0.800*X[0]*X[1]/pi)/pi + 0.294*cos(0.700*X[1]))*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0) + 0.490*X[1]*sin(X[0]) + 0.700*X[0]*cos(X[0]*X[1]/pi)/pi - 25.9*sin(1.30*X[0]*X[1]/pi) + 12.6*sin(1.70*X[0]) - 0.700*sin(X[1]) + 9.10*cos(2.10*X[1]) - 43.4

def forcing_u2(X):
    return 0.600*(11.1*sin(1.30*X[0]*X[1]/pi) - 5.40*sin(1.70*X[0]) - 3.90*cos(2.10*X[1]) + 18.6)*X[1]*sin(X[0])*cos(X[1]) - 0.300*(-0.600*sin(X[1]) + cos(X[0]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) - (0.600*sin(X[1]) + cos(X[0]) + 2.50)*(11.1*sin(1.30*X[0]*X[1]/pi) - 5.40*sin(1.70*X[0]) - 3.90*cos(2.10*X[1]) + 18.6)*sin(X[0]) + 0.600*(48.1*X[1]*cos(1.30*X[0]*X[1]/pi)/pi - 30.6*cos(1.70*X[0]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*sin(X[0])/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) - 0.600*(1.02*X[1]*cos(0.600*X[0]*X[1]/pi)/pi - 2.66*cos(0.700*X[0]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)*sin(X[0])/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)**2 + 1.20*(-0.320*X[1]*sin(0.800*X[0]*X[1]/pi)/pi - 0.540*sin(0.600*X[0]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)*sin(X[0])/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) - 0.300*(48.1*X[0]*cos(1.30*X[0]*X[1]/pi)/pi + 27.3*sin(2.10*X[1]))*(X[1]*cos(X[0]) + 0.600*cos(X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) + 0.300*(1.02*X[0]*cos(0.600*X[0]*X[1]/pi)/pi - 3.44*sin(0.800*X[1]))*(X[1]*cos(X[0]) + 0.600*cos(X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)**2 - 0.600*(-0.320*X[0]*sin(0.800*X[0]*X[1]/pi)/pi + 0.420*cos(0.700*X[1]))*(X[1]*cos(X[0]) + 0.600*cos(X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) + 0.600*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)*cos(X[0])/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) + 0.666666666666667*(48.1*X[1]*cos(1.30*X[0]*X[1]/pi)/pi - 30.6*cos(1.70*X[0]))*(0.180*sin(0.700*X[1]) + 0.120*cos(0.800*X[0]*X[1]/pi) + 0.270*cos(0.600*X[0]) + 0.270) + 0.666666666666667*(-0.0960*X[1]*sin(0.800*X[0]*X[1]/pi)/pi - 0.162*sin(0.600*X[0]))*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0) + 0.300*X[1]*cos(X[0]*X[1]/pi)/pi - 11.1*sin(1.30*X[0]*X[1]/pi) + 5.40*sin(1.70*X[0]) + 0.126*sin(X[1]) + 0.510*cos(X[0]) + 3.90*cos(2.10*X[1]) - 18.6

def forcing_v2(X):
    return (0.600*sin(X[1]) + cos(X[0]) + 2.50)*(11.1*sin(1.30*X[0]*X[1]/pi) - 5.40*sin(1.70*X[0]) - 3.90*cos(2.10*X[1]) + 18.6)*X[1]*cos(X[0]) + 0.300*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)*X[1]*sin(X[0])/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) + (11.1*sin(1.30*X[0]*X[1]/pi) - 5.40*sin(1.70*X[0]) - 3.90*cos(2.10*X[1]) + 18.6)*X[1]*sin(X[0])**2 - 0.300*(48.1*X[1]*cos(1.30*X[0]*X[1]/pi)/pi - 30.6*cos(1.70*X[0]))*(X[1]*cos(X[0]) + 0.600*cos(X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) + 0.300*(1.02*X[1]*cos(0.600*X[0]*X[1]/pi)/pi - 2.66*cos(0.700*X[0]))*(X[1]*cos(X[0]) + 0.600*cos(X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)**2 - 0.600*(-0.320*X[1]*sin(0.800*X[0]*X[1]/pi)/pi - 0.540*sin(0.600*X[0]))*(X[1]*cos(X[0]) + 0.600*cos(X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) - 0.600*(48.1*X[0]*cos(1.30*X[0]*X[1]/pi)/pi + 27.3*sin(2.10*X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*sin(X[0])/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) + 0.600*(1.02*X[0]*cos(0.600*X[0]*X[1]/pi)/pi - 3.44*sin(0.800*X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)*sin(X[0])/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)**2 - 1.20*(-0.320*X[0]*sin(0.800*X[0]*X[1]/pi)/pi + 0.420*cos(0.700*X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)*sin(X[0])/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) + 0.666666666666667*(48.1*X[0]*cos(1.30*X[0]*X[1]/pi)/pi + 27.3*sin(2.10*X[1]))*(0.180*sin(0.700*X[1]) + 0.120*cos(0.800*X[0]*X[1]/pi) + 0.270*cos(0.600*X[0]) + 0.270) + 0.666666666666667*(-0.0960*X[0]*sin(0.800*X[0]*X[1]/pi)/pi + 0.126*cos(0.700*X[1]))*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0) + 0.210*X[1]*sin(X[0]) + 0.300*X[0]*cos(X[0]*X[1]/pi)/pi - 11.1*sin(1.30*X[0]*X[1]/pi) + 5.40*sin(1.70*X[0]) - 0.300*sin(X[1]) + 3.90*cos(2.10*X[1]) - 18.6

def forcing_temperature1(X):
    return (4.81*X[0]*cos(1.30*X[0]*X[1]/pi)/pi + 2.73*sin(2.10*X[1]))*X[1]*sin(X[0]) - ((0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) + 1.00)*(-6.25300000000000*X[0]**2*sin(1.30*X[0]*X[1]/pi)/pi**2 - 6.25300000000000*X[1]**2*sin(1.30*X[0]*X[1]/pi)/pi**2 + 5.20200000000000*sin(1.70*X[0]) + 5.73300000000000*cos(2.10*X[1])) + (0.600*sin(X[1]) + cos(X[0]) + 2.50)*(4.81*X[1]*cos(1.30*X[0]*X[1]/pi)/pi - 3.06*cos(1.70*X[0])) - (4.81*X[1]*cos(1.30*X[0]*X[1]/pi)/pi - 3.06*cos(1.70*X[0]))*((48.1*X[1]*cos(1.30*X[0]*X[1]/pi)/pi - 30.6*cos(1.70*X[0]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) - (1.02*X[1]*cos(0.600*X[0]*X[1]/pi)/pi - 2.66*cos(0.700*X[0]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)**2 + 2*(-0.320*X[1]*sin(0.800*X[0]*X[1]/pi)/pi - 0.540*sin(0.600*X[0]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)) - (4.81*X[0]*cos(1.30*X[0]*X[1]/pi)/pi + 2.73*sin(2.10*X[1]))*((48.1*X[0]*cos(1.30*X[0]*X[1]/pi)/pi + 27.3*sin(2.10*X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) - (1.02*X[0]*cos(0.600*X[0]*X[1]/pi)/pi - 3.44*sin(0.800*X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)**2 + 2*(-0.320*X[0]*sin(0.800*X[0]*X[1]/pi)/pi + 0.420*cos(0.700*X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20))

def forcing_temperature2(X):
    return (4.81*X[0]*cos(1.30*X[0]*X[1]/pi)/pi + 2.73*sin(2.10*X[1]))*X[1]*sin(X[0]) - ((0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) + 1.00)*(-6.25300000000000*X[0]**2*sin(1.30*X[0]*X[1]/pi)/pi**2 - 6.25300000000000*X[1]**2*sin(1.30*X[0]*X[1]/pi)/pi**2 + 5.20200000000000*sin(1.70*X[0]) + 5.73300000000000*cos(2.10*X[1])) + (0.600*sin(X[1]) + cos(X[0]) + 2.50)*(4.81*X[1]*cos(1.30*X[0]*X[1]/pi)/pi - 3.06*cos(1.70*X[0])) - (4.81*X[1]*cos(1.30*X[0]*X[1]/pi)/pi - 3.06*cos(1.70*X[0]))*((48.1*X[1]*cos(1.30*X[0]*X[1]/pi)/pi - 30.6*cos(1.70*X[0]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) - (1.02*X[1]*cos(0.600*X[0]*X[1]/pi)/pi - 2.66*cos(0.700*X[0]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)**2 + 2*(-0.320*X[1]*sin(0.800*X[0]*X[1]/pi)/pi - 0.540*sin(0.600*X[0]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)) - (4.81*X[0]*cos(1.30*X[0]*X[1]/pi)/pi + 2.73*sin(2.10*X[1]))*((48.1*X[0]*cos(1.30*X[0]*X[1]/pi)/pi + 27.3*sin(2.10*X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) - (1.02*X[0]*cos(0.600*X[0]*X[1]/pi)/pi - 3.44*sin(0.800*X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)**2 + 2*(-0.320*X[0]*sin(0.800*X[0]*X[1]/pi)/pi + 0.420*cos(0.700*X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20))

def forcing_ke1(X):
    return (-0.320*X[0]*sin(0.800*X[0]*X[1]/pi)/pi + 0.420*cos(0.700*X[1]))*(25.9*sin(1.30*X[0]*X[1]/pi) - 12.6*sin(1.70*X[0]) - 9.10*cos(2.10*X[1]) + 43.4)*X[1]*sin(X[0]) + (0.600*sin(X[1]) + cos(X[0]) + 2.50)*(-0.320*X[1]*sin(0.800*X[0]*X[1]/pi)/pi - 0.540*sin(0.600*X[0]))*(25.9*sin(1.30*X[0]*X[1]/pi) - 12.6*sin(1.70*X[0]) - 9.10*cos(2.10*X[1]) + 43.4) - 0.700*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)*(X[1]**2*cos(X[0])**2 + 1.20*X[1]*cos(X[0])*cos(X[1]) + 4*sin(X[0])**2 + 0.360*cos(X[1])**2)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) - (0.700*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) + 0.490)*(-0.256*X[0]**2*cos(0.800*X[0]*X[1]/pi)/pi**2 - 0.256*X[1]**2*cos(0.800*X[0]*X[1]/pi)/pi**2 - 0.294*sin(0.700*X[1]) - 0.324*cos(0.600*X[0])) - (-0.320*X[1]*sin(0.800*X[0]*X[1]/pi)/pi - 0.540*sin(0.600*X[0]))*(0.700*(48.1*X[1]*cos(1.30*X[0]*X[1]/pi)/pi - 30.6*cos(1.70*X[0]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) - 0.700*(1.02*X[1]*cos(0.600*X[0]*X[1]/pi)/pi - 2.66*cos(0.700*X[0]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)**2 + 1.40*(-0.320*X[1]*sin(0.800*X[0]*X[1]/pi)/pi - 0.540*sin(0.600*X[0]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)) - (-0.320*X[0]*sin(0.800*X[0]*X[1]/pi)/pi + 0.420*cos(0.700*X[1]))*(0.700*(48.1*X[0]*cos(1.30*X[0]*X[1]/pi)/pi + 27.3*sin(2.10*X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) - 0.700*(1.02*X[0]*cos(0.600*X[0]*X[1]/pi)/pi - 3.44*sin(0.800*X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)**2 + 1.40*(-0.320*X[0]*sin(0.800*X[0]*X[1]/pi)/pi + 0.420*cos(0.700*X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)) + 0.700*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(48.1*X[0]*cos(1.30*X[0]*X[1]/pi)/pi + 48.1*X[1]*cos(1.30*X[0]*X[1]/pi)/pi + 27.3*sin(2.10*X[1]) - 30.6*cos(1.70*X[0]))/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) + (25.9*sin(1.30*X[0]*X[1]/pi) - 12.6*sin(1.70*X[0]) - 9.10*cos(2.10*X[1]) + 43.4)*(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)

def forcing_eps1(X):
    return (1.02*X[0]*cos(0.600*X[0]*X[1]/pi)/pi - 3.44*sin(0.800*X[1]))*(25.9*sin(1.30*X[0]*X[1]/pi) - 12.6*sin(1.70*X[0]) - 9.10*cos(2.10*X[1]) + 43.4)*X[1]*sin(X[0]) + (0.600*sin(X[1]) + cos(X[0]) + 2.50)*(1.02*X[1]*cos(0.600*X[0]*X[1]/pi)/pi - 2.66*cos(0.700*X[0]))*(25.9*sin(1.30*X[0]*X[1]/pi) - 12.6*sin(1.70*X[0]) - 9.10*cos(2.10*X[1]) + 43.4) - 0.700*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)*(X[1]**2*cos(X[0])**2 + 1.20*X[1]*cos(X[0])*cos(X[1]) + 4*sin(X[0])**2 + 0.360*cos(X[1])**2) - (0.700*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) + 0.490)*(-0.612*X[0]**2*sin(0.600*X[0]*X[1]/pi)/pi**2 - 0.612*X[1]**2*sin(0.600*X[0]*X[1]/pi)/pi**2 + 1.86200000000000*sin(0.700*X[0]) - 2.75200000000000*cos(0.800*X[1])) - (1.02*X[1]*cos(0.600*X[0]*X[1]/pi)/pi - 2.66*cos(0.700*X[0]))*(0.700*(48.1*X[1]*cos(1.30*X[0]*X[1]/pi)/pi - 30.6*cos(1.70*X[0]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) - 0.700*(1.02*X[1]*cos(0.600*X[0]*X[1]/pi)/pi - 2.66*cos(0.700*X[0]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)**2 + 1.40*(-0.320*X[1]*sin(0.800*X[0]*X[1]/pi)/pi - 0.540*sin(0.600*X[0]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)) - (1.02*X[0]*cos(0.600*X[0]*X[1]/pi)/pi - 3.44*sin(0.800*X[1]))*(0.700*(48.1*X[0]*cos(1.30*X[0]*X[1]/pi)/pi + 27.3*sin(2.10*X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) - 0.700*(1.02*X[0]*cos(0.600*X[0]*X[1]/pi)/pi - 3.44*sin(0.800*X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)**2 + 1.40*(-0.320*X[0]*sin(0.800*X[0]*X[1]/pi)/pi + 0.420*cos(0.700*X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)) + 0.700*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)*(48.1*X[0]*cos(1.30*X[0]*X[1]/pi)/pi + 48.1*X[1]*cos(1.30*X[0]*X[1]/pi)/pi + 27.3*sin(2.10*X[1]) - 30.6*cos(1.70*X[0])) + (25.9*sin(1.30*X[0]*X[1]/pi) - 12.6*sin(1.70*X[0]) - 9.10*cos(2.10*X[1]) + 43.4)*(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)**2/(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)

def forcing_ke2(X):
    return (-0.320*X[0]*sin(0.800*X[0]*X[1]/pi)/pi + 0.420*cos(0.700*X[1]))*(11.1*sin(1.30*X[0]*X[1]/pi) - 5.40*sin(1.70*X[0]) - 3.90*cos(2.10*X[1]) + 18.6)*X[1]*sin(X[0]) + (0.600*sin(X[1]) + cos(X[0]) + 2.50)*(-0.320*X[1]*sin(0.800*X[0]*X[1]/pi)/pi - 0.540*sin(0.600*X[0]))*(11.1*sin(1.30*X[0]*X[1]/pi) - 5.40*sin(1.70*X[0]) - 3.90*cos(2.10*X[1]) + 18.6) - 0.300*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)*(X[1]**2*cos(X[0])**2 + 1.20*X[1]*cos(X[0])*cos(X[1]) + 4*sin(X[0])**2 + 0.360*cos(X[1])**2)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) - (0.300*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) + 0.210)*(-0.256*X[0]**2*cos(0.800*X[0]*X[1]/pi)/pi**2 - 0.256*X[1]**2*cos(0.800*X[0]*X[1]/pi)/pi**2 - 0.294*sin(0.700*X[1]) - 0.324*cos(0.600*X[0])) - (-0.320*X[1]*sin(0.800*X[0]*X[1]/pi)/pi - 0.540*sin(0.600*X[0]))*(0.300*(48.1*X[1]*cos(1.30*X[0]*X[1]/pi)/pi - 30.6*cos(1.70*X[0]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) - 0.300*(1.02*X[1]*cos(0.600*X[0]*X[1]/pi)/pi - 2.66*cos(0.700*X[0]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)**2 + 0.600*(-0.320*X[1]*sin(0.800*X[0]*X[1]/pi)/pi - 0.540*sin(0.600*X[0]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)) - (-0.320*X[0]*sin(0.800*X[0]*X[1]/pi)/pi + 0.420*cos(0.700*X[1]))*(0.300*(48.1*X[0]*cos(1.30*X[0]*X[1]/pi)/pi + 27.3*sin(2.10*X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) - 0.300*(1.02*X[0]*cos(0.600*X[0]*X[1]/pi)/pi - 3.44*sin(0.800*X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)**2 + 0.600*(-0.320*X[0]*sin(0.800*X[0]*X[1]/pi)/pi + 0.420*cos(0.700*X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)) + 0.300*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(48.1*X[0]*cos(1.30*X[0]*X[1]/pi)/pi + 48.1*X[1]*cos(1.30*X[0]*X[1]/pi)/pi + 27.3*sin(2.10*X[1]) - 30.6*cos(1.70*X[0]))/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) + (11.1*sin(1.30*X[0]*X[1]/pi) - 5.40*sin(1.70*X[0]) - 3.90*cos(2.10*X[1]) + 18.6)*(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)

def forcing_eps2(X):
    return (1.02*X[0]*cos(0.600*X[0]*X[1]/pi)/pi - 3.44*sin(0.800*X[1]))*(11.1*sin(1.30*X[0]*X[1]/pi) - 5.40*sin(1.70*X[0]) - 3.90*cos(2.10*X[1]) + 18.6)*X[1]*sin(X[0]) + (0.600*sin(X[1]) + cos(X[0]) + 2.50)*(1.02*X[1]*cos(0.600*X[0]*X[1]/pi)/pi - 2.66*cos(0.700*X[0]))*(11.1*sin(1.30*X[0]*X[1]/pi) - 5.40*sin(1.70*X[0]) - 3.90*cos(2.10*X[1]) + 18.6) - 0.300*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)*(X[1]**2*cos(X[0])**2 + 1.20*X[1]*cos(X[0])*cos(X[1]) + 4*sin(X[0])**2 + 0.360*cos(X[1])**2) - (0.300*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) + 0.210)*(-0.612*X[0]**2*sin(0.600*X[0]*X[1]/pi)/pi**2 - 0.612*X[1]**2*sin(0.600*X[0]*X[1]/pi)/pi**2 + 1.86200000000000*sin(0.700*X[0]) - 2.75200000000000*cos(0.800*X[1])) - (1.02*X[1]*cos(0.600*X[0]*X[1]/pi)/pi - 2.66*cos(0.700*X[0]))*(0.300*(48.1*X[1]*cos(1.30*X[0]*X[1]/pi)/pi - 30.6*cos(1.70*X[0]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) - 0.300*(1.02*X[1]*cos(0.600*X[0]*X[1]/pi)/pi - 2.66*cos(0.700*X[0]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)**2 + 0.600*(-0.320*X[1]*sin(0.800*X[0]*X[1]/pi)/pi - 0.540*sin(0.600*X[0]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)) - (1.02*X[0]*cos(0.600*X[0]*X[1]/pi)/pi - 3.44*sin(0.800*X[1]))*(0.300*(48.1*X[0]*cos(1.30*X[0]*X[1]/pi)/pi + 27.3*sin(2.10*X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20) - 0.300*(1.02*X[0]*cos(0.600*X[0]*X[1]/pi)/pi - 3.44*sin(0.800*X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)**2 + 0.600*(-0.320*X[0]*sin(0.800*X[0]*X[1]/pi)/pi + 0.420*cos(0.700*X[1]))*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)) + 0.300*(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)*(48.1*X[0]*cos(1.30*X[0]*X[1]/pi)/pi + 48.1*X[1]*cos(1.30*X[0]*X[1]/pi)/pi + 27.3*sin(2.10*X[1]) - 30.6*cos(1.70*X[0])) + (11.1*sin(1.30*X[0]*X[1]/pi) - 5.40*sin(1.70*X[0]) - 3.90*cos(2.10*X[1]) + 18.6)*(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)**2/(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)

def P1(X):
    return (0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)*(X[1]**2*cos(X[0])**2 + 1.20*X[1]*cos(X[0])*cos(X[1]) + 4*sin(X[0])**2 + 0.360*cos(X[1])**2)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)

def P2(X):
    return (0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)*(X[1]**2*cos(X[0])**2 + 1.20*X[1]*cos(X[0])*cos(X[1]) + 4*sin(X[0])**2 + 0.360*cos(X[1])**2)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)

def P_eps1(X):
    return (0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)*(X[1]**2*cos(X[0])**2 + 1.20*X[1]*cos(X[0])*cos(X[1]) + 4*sin(X[0])**2 + 0.360*cos(X[1])**2)

def P_eps2(X):
    return (0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)*(X[1]**2*cos(X[0])**2 + 1.20*X[1]*cos(X[0])*cos(X[1]) + 4*sin(X[0])**2 + 0.360*cos(X[1])**2)

def A1(X):
    return (37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)*(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)/(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)

def A2(X):
    return (37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)*(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)/(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)

def B1(X):
    return -(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(48.1*X[0]*cos(1.30*X[0]*X[1]/pi)/pi + 48.1*X[1]*cos(1.30*X[0]*X[1]/pi)/pi + 27.3*sin(2.10*X[1]) - 30.6*cos(1.70*X[0]))/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)

def B2(X):
    return -(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(48.1*X[0]*cos(1.30*X[0]*X[1]/pi)/pi + 48.1*X[1]*cos(1.30*X[0]*X[1]/pi)/pi + 27.3*sin(2.10*X[1]) - 30.6*cos(1.70*X[0]))/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)

def B_eps1(X):
    return -(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)*(48.1*X[0]*cos(1.30*X[0]*X[1]/pi)/pi + 48.1*X[1]*cos(1.30*X[0]*X[1]/pi)/pi + 27.3*sin(2.10*X[1]) - 30.6*cos(1.70*X[0]))

def B_eps2(X):
    return -(0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)*(48.1*X[0]*cos(1.30*X[0]*X[1]/pi)/pi + 48.1*X[1]*cos(1.30*X[0]*X[1]/pi)/pi + 27.3*sin(2.10*X[1]) - 30.6*cos(1.70*X[0]))

def EV1(X):
    return (0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)

def EV2(X):
    return (0.600*sin(0.700*X[1]) + 0.400*cos(0.800*X[0]*X[1]/pi) + 0.900*cos(0.600*X[0]) + 0.900)**2*(37.0*sin(1.30*X[0]*X[1]/pi) - 18.0*sin(1.70*X[0]) - 13.0*cos(2.10*X[1]) + 62.0)/(1.70*sin(0.600*X[0]*X[1]/pi) - 3.80*sin(0.700*X[0]) + 4.30*cos(0.800*X[1]) + 8.20)

def velocity1(X):
   return [u1(X), v1(X)]

def forcing_velocity1(X):
   return [forcing_u1(X), forcing_v1(X)]

def velocity2(X):
   return [u2(X), v2(X)]

def forcing_velocity2(X):
   return [forcing_u2(X), forcing_v2(X)]

def A_ke1(X):
   return A1([X[0],X[1]])*ke1([X[0],X[1]])

def A_eps1(X):
   return A1([X[0],X[1]])*eps1([X[0],X[1]])

def A_ke2(X):
   return A2([X[0],X[1]])*ke2([X[0],X[1]])

def A_eps2(X):
   return A2([X[0],X[1]])*eps2([X[0],X[1]])
