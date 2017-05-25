# -*- coding: utf-8 -*-
"""
Created on Wed May 17 17:51:26 2017
@author: ANSH
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation


#N discretizes x

N = 2**11


time = 0.0

#set up x as a numpy array
dx = 0.1
x = dx*(np.arange(N)-0.5*N)

dt = 0.1


#analytical solution for wavefunction

def wavefunction(x, time, width, position, mass, hbar, momentum):
    c = 1.0/width
    b = position
    h = hbar
    m = mass
    t = time
    k = momentum/hbar
    
    coeff1 = (2*c/np.pi)**(1/4)
    coeff2 = np.exp(1.0j*k*b-(k**2)/(4*c))
    
    gamma = np.sqrt(1.0 + (2.0j*h*c*t)/m)
    
    coeff = coeff1*coeff2/gamma
    
    exponum = c*(1.0j*(b-x) - k/(2*c))**2
    
    expo = exponum * gamma**(-2)
    
    psi = coeff * np.exp(expo)
    
    
    return psi




#Setting up plot

fig = plt.figure()

xlim = (np.amin(x), np.amax(x))
ymin = -1
ymax = 1
ax1 = fig.add_subplot(111, xlim=xlim,
                      ylim=(ymin - 0.2 * (ymax - ymin),
                            ymax + 0.2 * (ymax - ymin)))
psi_x_line, = ax1.plot([], [], c='r', label=r'$|\psi(x)|$')
#psi_y_line, = ax1.plot([], [], c='b', label='$|\psi(x)|$')
#psi_abs_line, = ax1.plot([], [], c='g', label='$|\psi(x)|$')

title = ax1.set_title("")
ax1.legend(prop=dict(size=12))
ax1.set_xlabel('$x$')
ax1.set_ylabel(r'$|\psi(x)|$')
#ax1.set_ylabel('$|\psi(x)|$')
#ax1.set_ylabel('$|\psi(x)|$')


#initialising lines
def init():
    psi_x_line.set_data([], [])
#    psi_y_line.set_data([],[])
#    psi_abs_line.set_data([],[])

    return psi_x_line, #psi_y_line,



#Setting up animation. showing real and imaginary parts is optional, feel free to remove it

def animate(i):
    
    global time
    global dt
    
    w = 100
    pos = -100
    energy = 5
    mass = 1
    hbar = 1
    forwards = True
    
    p = (2*mass*energy)**(1/2)
    
    if not forwards:
        p = -1*p
    
    psi = wavefunction(x, time, w, pos, mass, hbar, p)
    
    
    
    psi_x_line.set_data(x, psi.real)
    
#    psi_y_line.set_data(x, psi.imag)
#    
#    psi_abs_line.set_data(x, abs(psi))
    
    time += dt
    
    return psi_x_line, #psi_y_line, psi_abs_line,


#using animation function
anim = animation.FuncAnimation(fig, animate, init_func=init,
                              frames=12000, interval=10, blit=True)


plt.show()

