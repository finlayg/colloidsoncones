#Calculates and prints the truncation radius required to fit a given number of equilibrium distances around the circle at the truncation

import numpy as np
import math

circumference = float(input("Circumference of polygon at peak in equilibrium distances: "))

sides = int(math.floor(circumference))
print("Sides: "+str(sides))
sidelength = circumference/sides
print("Side length: "+str(sidelength))

truncrad = sidelength/(2*np.sin((180/sides)*(np.pi/180)))

print("Truncation Radius = "+str(truncrad))
