#script to generate input dimension values [r_b, r_t, h] for 'data' file from theta_max in degrees

#!/usr/bin/env python3

#angle = float(input('Cone Angle theta_max in degrees: '))

#import ~/usr/include/numpy as np

from math import pi
from math import tan
from math import acos
from math import sqrt

def datafileinputvals(thetamax, rt, area):
    #area = 292.1681168 #pi*rb*h-(pi*rt*(hc-h))
    #area = 1800


    #convert thetamax to radians
    thetamax = thetamax*(pi/180)
    
    tantheta = tan(acos(thetamax / (2*pi)))
    rb = sqrt(rt**2 + area / (pi * tantheta))
    h = (rb - rt) * tantheta
    
    vals = [rb, rt, h]
    return vals

#print(datafileinputvals(angle))
