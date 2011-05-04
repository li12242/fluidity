t = var('t')
h = var('h')

# --------------------------------------------------

u = (t+1)*h*sin(2*pi*x)
eta = (t+1)*h*cos(2*pi*x)

u_src = diff(u, t) + diff(eta, x)
eta_src = diff(eta, t) + diff(u, x)

print "u: ", u
print "eta: ", eta

print "u_src: ", u_src
print "eta_src: ", eta_src

J = integrate((eta.subs(t=1))**2, x, 0, 1)
print "J(t=1): ", J.subs(h=1)
print "diff(J(t=1), h).subs(h=1): ", diff(J, h).subs(h=1)
J = integrate((eta.subs(t=0))**2, x, 0, 1)
print "J(t=0): ", J.subs(h=1)
print "diff(J(t=0), h).subs(h=1): ", diff(J, h).subs(h=1)
