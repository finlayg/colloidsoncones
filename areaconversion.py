#Script to convert coordinates of crystal on cone of given angle, rt, area to new cone with same angle and different rt, area
#Takes lowest output and turns to coords file for input

import numpy as np

inputfile = input('lowest Input Filename: ')
inputminima = float(input('Minima coordinates to convert: '))
numparts = int(input('Number of particles: '))

whichangle = input('Input angle of cone as theta_p or theta_max? p/m ')
if whichangle is 'p':
	anglep = float(input('Angle theta_p of cone? '))*(np.pi/180) #angle of cone converted to radians
	angle = (2*np.pi*np.sin(anglep/2))#/(np.pi/180)
elif whichangle is 'm':
	angle = float(input('Angle theta_max of cone? '))*(np.pi/180) #angle of net, converted to radians
	anglep = 2*np.arcsin(angle/(2*np.pi)) 
else:
	print('Error: input not recognised')

inputrt = float(input('r_t of input cone: '))

inputareadef = input('Input area of 292.1681168? Y/N ')
if inputareadef is 'Y':
	inputarea = 292.1681168
elif inputareadef is 'N':
	inputarea = float(input('Area of output cone: '))
else:
	print('Error: input not recognised')

outputrt = float(input('r_t of output cone: '))

outputareadef = input('Output area of 292.1681168? Y/N ')
if outputareadef is 'Y':
	outputarea = 292.1681168
elif outputareadef is 'N':
	outputarea = float(input('Area of output cone: '))
else:
	print('Error: input not recognised')

print('----------')
print('Angle of cone net in degrees: ', angle/(np.pi/180))
print('Input dimensions: r_t: ',inputrt,', area: ',inputarea)
print('Output dimensions: r_t: ',outputrt,', area: ',outputarea)
print('----------')

#zb is vertical dist from top of cone to bottom of available area
#zt is vertical dist from top of cone to top of available area
#input is first area, output is second area
#zinput is absolute location of z origin on area input
#zoutput is absolute location of z origin on area output

ztinput = inputrt / (np.tan(anglep/2))
zbinput = np.sqrt(((inputarea/np.pi) + (inputrt**2)/(np.cos(anglep/2)) ) * (np.cos(anglep/2)) ) / (np.tan(anglep/2))

zinput = (zbinput - ztinput)/2

ztoutput = outputrt / (np.tan(anglep/2))
zboutput = np.sqrt(((outputarea/np.pi) + (outputrt**2)/(np.cos(anglep/2)) ) * (np.cos(anglep/2)) ) / (np.tan(anglep/2))

zoutput = (zboutput - ztoutput)/2

tranval = zoutput - zinput

coordslist = []
n = 0

with open(inputfile, 'r') as stream1:
	for line in stream1:
		if n < (float(inputminima)*2*numparts) and 'C     ' in line:
			splitline = line.rstrip().split(' ')
			new = [item for item in splitline if item!='']
			coordslist.append([new[1], new[2], new[3]])
			n += 1
#second set of points for each minima in lowest is the net of cone - don't want this
wantedcoordslist = coordslist[(-2*numparts):(-1*numparts)]

for i in range(0, numparts):
	coordslist[i][0] = wantedcoordslist[i][0]
	coordslist[i][1] = wantedcoordslist[i][1]
	coordslist[i][2] = float(wantedcoordslist[i][2]) + float(tranval)

#Constructing output coords file
outputfile = ''

for i in range(0, numparts):
        output = '   ' + str(coordslist[i][0]) + '   ' + str(coordslist[i][1]) + '   ' + str(coordslist[i][2]) + '\n'
        outputfile = outputfile + output

print("coordsareaconv file produced")

converted_file = open('coordsareaconv', "w+")
converted_file.write(outputfile)
