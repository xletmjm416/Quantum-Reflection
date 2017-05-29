# -*- coding: utf-8 -*-
"""
Created on Mon May 29 21:58:43 2017

@author: mjm416
"""
import numpy as np
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
    
filename = raw_input("Name of the file to analyse: ")
f = open(filename, 'r+')
data = f.readline()
data = np.array(parse_complex(data))
x = np.linspace(0,10,128)
plt.plot(x,np.real(data), 'r-')

for i in range(50):
    data = f.readline()
    
data = np.array(parse_complex(data))
x = np.linspace(0,10,128)
plt.plot(x,np.real(data), 'b-')


f.close()