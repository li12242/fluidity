t = var('t')

u = 1/sqrt(2*pi) * exp(-x**2 / 2)
v = 1.0

lmbda = 0.5 * sin(3*pi/20 * (x+10))

j = diff(lmbda, t) - v * diff(lmbda, [x, x])
s = str(j).replace('e^', 'exp').replace('^', '**')
print s

endpoints = (-10.0, 10.0)
bcs = [lmbda(x=y).n() for y in endpoints]
print bcs

print "j * u: ", j * u
true_J = integrate(j * u, x, -10.0, 10.0)
print true_J
# Alas, a bug in sage stops it from being able to
# print out float(true_J). However,
# a numerical value can be extracted
# with wolfram alpha.
