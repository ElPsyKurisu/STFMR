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
from ekpy.control.instruments import digilent
import time

def read_callibration_file(filepath):
    """
    Reads the given filepath as a callibration file where the left column should be in volts (V)
    and the right column in magnetic field (Oe). Note that the callibration file should bound all
    desired field values as this program will not extrapolate field values.
       
    args:
        filepath: The path to the callibration file (str)

    returns:
        voltage_input: numpy array of voltage values (ndarray)
        field_output: numpy array of field values (ndarray)
    """    
    df = pd.read_csv(filepath, delimiter=r"\s+", names=['volts', 'field'])
    voltage_input = df['volts'].values
    field_output = df['field'].values
    return voltage_input, field_output

'''
Now I want a function that I can input whatever field strength I want and it will calculate what voltage I need
make sure to throw erros if out of range
'''


def field_to_volts(field, callibration_file, poll_high=True) -> float:
    '''
    Reads the given filepath of the callibration file using read_callibration_file and performs
    a linear interpolation on the data by finding the two nearest points in the callibration file
    and obtaining a line from them to interpolate the desired field. This is being edited to not use 
    numpy search sorted and do it manually to allow for polling
       
    args:
        callibration_file: The path to the callibration file (str)
        field: The desired field (float)

    returns:
        volts: calculated input voltage to reach desired field in volts (float)
        IndexError: returns none if given field is out of range (None)
    '''
    try:
        voltage_in, field_out = read_callibration_file(callibration_file)
        if not poll_high:
            #reverses lists
            voltage_in = voltage_in[::-1]
            field_out = field_out[::-1]
        higher_index = 0
        if field >= 0:
            for i in range(len(field_out)):
                current = field_out[i]
                next = field_out[i+1]
                if next < field and current > field:
                    higher_index = i
                    break
        else:
            for i in range(len(field_out)):
                current = field_out[i]
                next = field_out[i+1]
                if next > field and current < field:
                    higher_index = i
                    break
        lower_index = higher_index - 1
        #this is needed to ensure that we dont take the value of arr[-1] which returns the last element
        if lower_index < 0:
            raise IndexError
        slope = (voltage_in[higher_index]-voltage_in[lower_index])/(field_out[higher_index]-field_out[lower_index])
        intercept = voltage_in[higher_index]-slope*field_out[higher_index]
        volts = slope*field + intercept
        return volts
    except IndexError:
        #change "raise ValueError" to print to remove the error and simply print out the message
        raise IndexError('The given field is out of range of the callibration file')

'''
Now the next function will allow you to set an angle u want with a field strength that you want
by using the two coils
'''

def vectorized_magnetic_field(x_callibration_file, y_callibration_file, angle, magnitude, poll_high=True) -> float:
    '''
    Calculates the parameterized form of the given angle and magnitude and calls the function
    field_to_volts to convert to the required ouput voltages on the x_coil and y_coil.
       
    args:
        x_callibration_file: The path to the callibration file for coil x (str)
        y_callibration_file: The path to the callibration file for coil y (str)
        angle: The desired angle (float)
        magnitude: The desired field (float)

    returns:
        volts: calculated input voltage to reach desired field in volts (float)
        IndexError: returns none if given field is out of range (None)
    '''
    x_magnitude = magnitude*np.cos(np.deg2rad(angle))
    y_magnitude = magnitude*np.sin(np.deg2rad(angle))
    x_volts = field_to_volts(x_magnitude, x_callibration_file, poll_high)
    y_volts = field_to_volts(y_magnitude, y_callibration_file, poll_high)
    return x_volts, y_volts

'''
Our next step is to call the functions that actually output stuff directly here using the digilient instrument library
in the ekpy package and ensures you cycle the field by first polling the field and then bringing it along a certain curve
could add an optional arg that is poll high or poll low, aka uses diff graphs
'''
def output_vector_field(x_callibration_file, y_callibration_file, angle, magnitude, poll_high=True):

    x_volts, y_volts = vectorized_magnetic_field(x_callibration_file, y_callibration_file, angle, magnitude, poll_high)
    if poll_high:
        digilent.v_out(0,0,10)
        digilent.v_out(0,1,10)
        time.sleep(1)
    else:
        digilent.v_out(0,0,-10)
        digilent.v_out(0,1,-10)
        time.sleep(1)
    digilent.v_out(0,0,x_volts)
    digilent.v_out(0,1,y_volts)


#test = vectorized_magnetic_field("utils\\1to3_callibration.txt", "utils\\2to4_callibration.txt", 0, 0, False)
#print(test)
