'''
Ok so this program will be deleted and incorporated into "stfmr.py" later once I build that
This is solely to test the vector magnet and control it using callibration files
'''

'''
First function should be to read the callibration files and give us some method to interpolate values for an arbitrary field value
we will use pandas to read the file
'''

import pandas as pd
import numpy as np
from numpy.polynomial import Polynomial
import matplotlib.pyplot as plt

def read_callibration_file(filepath):
    """
    Reads the given filepath as a callibration file where the left column should be in volts (V)
    and the right column in magnetic field (Oe)
       
    args:
        filepath: The path to the callibration file (str)

    returns:
        voltage_input: numpy array of voltage values (ndarray)
        field_output: numpy array of field values (ndarray)
    """    
    df = pd.read_csv(filepath, sep="\t", names=['volts', 'field'])
    voltage_input = df['volts'].values
    field_output = df['field'].values
    return voltage_input, field_output

df_x = pd.read_csv("utils\\1to3_callibration.txt", sep="\t") 
df_y = pd.read_csv("utils\\2to4_callibration.txt", sep="\t")
# display DataFrame 

'''
Now I want a function that I can input whatever field strength I want and it will calculate what voltage I need
make sure to throw erros if out of range
'''

def field_to_volts(field, callibration_file):
    voltage_in, field_out = read_callibration_file(callibration_file)
    #get best fit line
    a, b = Polynomial.fit(field_out, voltage_in, 1)
    volts = a*field + b
    return volts
'''
meow = field_to_volts(500, "utils\\1to3_callibration.txt")
print('herro', meow)
'''
voltage_in, field_out = read_callibration_file("utils\\1to3_callibration.txt")
print(voltage_in)
print(field_out)
voltage_range = np.linspace(0,2, 500)
c,b,a = Polynomial.fit(voltage_in, field_out, 2)
print(b,a)
plt.scatter(voltage_in, field_out)
plt.plot(voltage_in, a*voltage_in**2 + b*voltage_in + c)
plt.show()

#i should just use scipy curve fit and do log fit