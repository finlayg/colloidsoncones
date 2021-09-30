#script to generate input 'data' file for cone for initial coords generation or for runs - make sure to check file before use

import numpy as np

outputfile = ''

#INPUTS
coordgenerator = input("Is this for an initial coord generation? Y/N: ")
if coordgenerator is 'N':
	GTKC = input("GTHOMKEEPCOORDS? Y/N: ")
	range = input("Short (rho = 18.0) or long (rho = 6.0) ranged interactions? S/L: ")
elif coordgenerator is 'Y':
	GTKC = 'N'
else:
	print("Error: not recognised")

checkd = input("Checking analytical vs numerical gradients? Y/N: ")
surfno = input("Number of surface? (8 for cylindrical, 9 for cyclical coords): ")

whichangle = input("Input angle: theta_max or theta_p? p/m: ")

if whichangle is 'p':
	anglep = float(input("Cone angle theta_p in degrees: "))
	angle = (2*np.pi*np.sin((anglep*(np.pi/180))/2))/(np.pi/180)
elif whichangle is 'm':
	angle = float(input("Cone angle theta_max in degrees: "))
else:
	print('Error: not recognised')

areaset = input("Use default area of 292.1681168? Y/N: ")

if areaset is 'N':
	area = float(input("Available area to use? "))
elif areaset is 'Y':
	area = float(292.1681168)
else:
	print('Error: Choose Y or N')

rt = float(input("Radius of truncation: "))
if coordgenerator is 'N':
	steps = input("Number of steps: ")
elif coordgenerator is 'Y':
	steps = '10'

#GENERATING FILE LINES
if checkd is 'Y':
	outputfile = outputfile + 'CHECKD \n'
if GTKC is 'Y':
	outputfile = outputfile + 'GTHOMKEEPCOORDS \n'

from dimensiongenerator import datafileinputvals
vals = datafileinputvals(angle, rt, area)
outputfile = outputfile + 'GTHOMSON ' + str(surfno) + ' ' + str(vals[0]) + ' '  + str(vals[1]) + ' ' + str(vals[2]) + '\n' 

#left generic GTHOMSON LINE: outputfile = outputfile + 'GTHOMSON n rb rt h \n'

if coordgenerator is 'Y':
	outputfile = outputfile + 'GTHOMSONPOT 6 1.0 4.0 \n'
elif coordgenerator is 'N':
	if range is 'L':
		outputfile = outputfile + 'GTHOMSONPOT 6 1.0 6.0 \n'
	elif range is 'S':
		outputfile = outputfile + 'GTHOMSONPOT 6 1.0 18.0 \n'

outputfile = outputfile + 'SAVE 100 \n'

outputfile = outputfile + 'SORT \n'

outputfile = outputfile + 'SLOPPYCONV 0.00001 \n'

outputfile = outputfile + 'TIGHTCONV  0.0000001 \n'


outputfile = outputfile + 'STEPS ' + steps + ' 1.0 \n'

outputfile = outputfile + 'MAXIT 10000 10000 \n'

outputfile = outputfile + 'STEP 0.3 \n'

outputfile = outputfile + 'TEMPERATURE 1.0 \n'

outputfile = outputfile + 'MAXBFGS 0.01 \n'

outputfile = outputfile + 'MAXERISE 1.0D-8'

converted_file = open("data", "w+")
converted_file.write(outputfile)

print('Check file before use.')
