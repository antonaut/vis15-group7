import operator
import matplotlib.pyplot as plt
import math

# Helpers


def norm(tup):
    length = math.sqrt(float(tup[0])*float(tup[0]) + float(tup[1])*float(tup[1]))
    return (float(tup[0])/length, float(tup[1])/length)

def smul(scalar, tup):
    return tuple([x * scalar for x in tup])

def add_tup(a, b):
    return tuple(map(operator.add, a, b))

def plot_lists(tuplist):
    xx = []
    yy = []
    for p in tuplist:
        xx.append(p[0])
        yy.append(p[1])
    return (xx, yy)

# Fun stuff

def integrator(field, start, s, nsteps, method):
    x_i = norm(start)
    stream_line = [x_i]
    for i in range(nsteps):
        xp = method(x_i, field, s)
        stream_line.append(xp)
        x_i = xp
    return stream_line

def euler(x_i, field, s):
    return add_tup(x_i, smul(s, field(x_i)))

def rk4(x_i, field, s):
        v1 = smul(1/6.0, field(x_i))
        v2 = smul(1/3.0, field(add_tup(x_i, smul(s/2.0, v1))))
        v3 = smul(1/3.0, field(add_tup(x_i, smul(s/2.0, v2))))
        v4 = smul(1/6.0, field(add_tup(x_i, smul(s, v3))))
        svec = add_tup(add_tup(v1, v2), add_tup(v3, v4))
        xp = add_tup(x_i, smul(s, svec))
        return xp

# Testing

def V(xvec):
    x = xvec[0]
    y = xvec[1]
    return (-y, x/2.0)
    #return ((x/2.0)*(x/2.0), (-y)*(-y))

S = (1, 0)


def pe():
    xx, yy = plot_lists(integrator(V, S, 0.1, 100, euler))
    plt.plot(xx, yy)

def prk():
    xx, yy = plot_lists(integrator(V, S, 0.3, 30, rk4))
    plt.plot(xx, yy)

def c():
    plt.close()
