from math import sin, cos, tanh, pi, sqrt

def u(X):
    return 0.600*sin(X[1]) + cos(X[0]) + 3.00

def v(X):
    return X[1]*sin(X[0])

def p(X):
    return sin(0.31830988618379069*X[0]*X[1]) + sin(X[0]) + cos(X[1]) - 1.00

def rho(X):
    return 2.50

def forcing_u(X):
    return 1.50*X[1]*sin(X[0])*cos(X[1]) - 0.0154212568767200*(-0.300*sin(X[1]) + 0.500*cos(X[0]))*(0.500*X[1]*cos(X[0]) + 0.300*cos(X[1]))*(X[1]*cos(X[0]) + 0.600*cos(X[1]))/sqrt(4*(0.500*X[1]*cos(X[0]) + 0.300*cos(X[1]))**2 + 4*sin(X[0])**2) - (0.00385531421918000*sqrt(4*(0.500*X[1]*cos(X[0]) + 0.300*cos(X[1]))**2 + 4*sin(X[0])**2) + 0.700)*(-0.600*sin(X[1]) + cos(X[0])) + (0.00771062843836000*sqrt(4*(0.500*X[1]*cos(X[0]) + 0.300*cos(X[1]))**2 + 4*sin(X[0])**2) + 1.40)*cos(X[0]) - (1.50*sin(X[1]) + 2.50*cos(X[0]) + 7.50)*sin(X[0]) + 0.31830988618379069*X[1]*cos(0.31830988618379069*X[0]*X[1]) - 0.0154212568767200*((0.500*X[1]*cos(X[0]) + 0.300*cos(X[1]))*X[1]*sin(X[0]) - 2*sin(X[0])*cos(X[0]))*sin(X[0])/sqrt(4*(0.500*X[1]*cos(X[0]) + 0.300*cos(X[1]))**2 + 4*sin(X[0])**2) + cos(X[0])

def forcing_v(X):
    return (0.00385531421918000*sqrt(4*(0.500*X[1]*cos(X[0]) + 0.300*cos(X[1]))**2 + 4*sin(X[0])**2) + 0.700)*X[1]*sin(X[0]) + (1.50*sin(X[1]) + 2.50*cos(X[0]) + 7.50)*X[1]*cos(X[0]) + 2.50*X[1]*sin(X[0])**2 - 0.0308425137534400*(-0.300*sin(X[1]) + 0.500*cos(X[0]))*(0.500*X[1]*cos(X[0]) + 0.300*cos(X[1]))*sin(X[0])/sqrt(4*(0.500*X[1]*cos(X[0]) + 0.300*cos(X[1]))**2 + 4*sin(X[0])**2) + 0.31830988618379069*X[0]*cos(0.31830988618379069*X[0]*X[1]) + 0.00771062843836000*(X[1]*cos(X[0]) + 0.600*cos(X[1]))*((0.500*X[1]*cos(X[0]) + 0.300*cos(X[1]))*X[1]*sin(X[0]) - 2*sin(X[0])*cos(X[0]))/sqrt(4*(0.500*X[1]*cos(X[0]) + 0.300*cos(X[1]))**2 + 4*sin(X[0])**2) - sin(X[1])

def velocity(X):
   return [u(X), v(X)]

def forcing_velocity(X):
   return [forcing_u(X), forcing_v(X)]

