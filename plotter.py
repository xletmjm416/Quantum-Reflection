# -*- coding: utf-8 -*-
"""
Created on Mon May 29 21:58:43 2017

@author: mjm416
"""
import numpy as np
from scipy import stats
from matplotlib import pyplot as plt
def str2complex(str):
    import re
    numeric_const_pattern = r"""
    [-+]? # optional sign
    (?:
    (?: \d* \. \d+ ) # .1 .12 .123 etc 9.1 etc 98.1 etc
    |
    (?: \d+ \.? ) # 1. 12. 123. etc 1 12 123 etc
    )
    # followed by optional exponent part if desired
    (?: [Ee] [+-]? \d+ ) ?"""
    rx = re.compile(numeric_const_pattern, re.VERBOSE)
    return rx.findall(str)

def list2complex(ls):
    return complex(float(ls[0]),float(ls[1]))

def parse_complex(data):
    data = data.split(';')
    retrieved = [list2complex(str2complex(s)) for s in data]
    return retrieved
    
def extract_from_file(filename):
    with open(filename, 'r+') as f:
        params = f.readline()
        data = []
        for line in f:
            data.append(np.array(parse_complex(line)))
    return params.split(','), np.array(data)

amplitude_arr = []
position_arr = []
spread_arr = []
filename = raw_input("Name of the file to analyse: ")
params, data = extract_from_file(filename)

x_step = float(params[0])
t_step = float(params[1])
x_size = float(params[2])
t_size = float(params[3])
init_pos = float(params[4]) #x_0
momentum = float(params[5]) #k
init_spread = float(params[6])   #sigma
N_space = x_size/x_step
N_time = t_size/t_step

x = np.linspace(0,x_size,N_space)
t = np.linspace(0,t_size,N_time)
i=0
for step in data:
    
    amplitude = np.abs(step)**2         #probability amplitude
    norm = np.sum(amplitude)
    position = np.average(x, weights=amplitude)
    spread = (np.average((x-position)**2, weights=amplitude))**0.5
    
    amplitude_arr.append(amplitude)
    position_arr.append(position)
    spread_arr.append(spread)
    
    if (i%10 == 0):
        pass
        plt.plot(x, amplitude, 'r-')
        plt.show()
    
    i+=1

analytic_speed = init_pos + momentum*t
plt.plot(t, position_arr)
plt.grid()
velocity, intercept = np.polyfit(t[0:5], position_arr[0:5], 1)
print "initial velocity", velocity
plt.xlabel('time')
plt.ylabel('position')
plt.show()
# from http://people.physics.tamu.edu/valery/Gaussian%20wave%20packet.pdf
analytic_spread = np.sqrt(init_spread**2+(1.0/(4.0*(init_spread**2))*t**2))
estimate = 1/(2*init_spread)*t
meas = plt.plot(t, spread_arr, 'r-', label='measured')
anal = plt.plot(t, analytic_spread, 'g-', label='analytic')
est = plt.plot(t, estimate, 'b--', label='HUP estimate')
plt.xlabel('time')
plt.ylabel('spread')
ax = plt.gca()
legend = ax.legend(loc='upper left', shadow=True)
plt.grid()
plt.show()