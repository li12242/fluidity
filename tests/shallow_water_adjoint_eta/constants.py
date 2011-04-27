from math import sin, cos, pi
import numpy

d0 = 1.0
theta = 0.5
h = 1.0
functional = 0.5
dfunctional = 1

def eta_src(x, t):
  return 2*pi*h*cos(2*pi*x)

def u_src(x, t):
  return numpy.array([-2*pi*h*sin(2*pi*x),0,0])

def u_exact(x, t):
  return numpy.array([h*sin(2*pi*x), 0.0, 0.0])

def eta_exact(x, t):
  return h*cos(2*pi*x)
