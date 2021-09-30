#Script to convert a given minima from the current angle to a slightly higher or lower angle
#Input a lowest file, outputs a coordsstretch file
#Redundant as just used coords from lower angle with new data file instead- this manually changes the coordinates whereas GMIN will naturally project them onto the new surface.

import numpy as np

inputfile = input('Lowest Input Filename: ')
datafile = input('Data filename for dimensions: ')
inputminima = float(input('Minima coordinates to convert: '))
numparts = int(input('Number of particles: '))

whichangle = input('Input angles theta_p or theta_max? p/m ')
if whichangle is 'p':
	anglep = float(input('Angle theta_p of cone? '))*(np.pi/180) #angle of cone converted to radians
	angle = (2*np.pi*np.sin(anglep/2))#/(np.pi/180)
	angleoutp = float(input('Angle theta_p of cone to convert to? '))*(np.pi/180) #output angle converted to radians
	angleout = (2*np.pi*np.sin(angleoutp/2))#/(np.pi/180)
elif whichangle is 'm':
	angle = float(input('Angle theta_max of cone? '))*(np.pi/180) #angle of net, converted to radians
	angleout = float(input('Angle theta_max of cone to convert to? '))*(np.pi/180)
else:
	print('Error: input not recognised')

print('----------')
with open(datafile, 'r') as stream2:
	for line in stream2:
		if 'GTHOMSON ' in line:
			splitline = line.rstrip().split(' ')
			new = [item for item in splitline if item!=' ']
			rb = float(new[2])
			rt = float(new[3])
			h = float(new[4])
			print('Dimensions: r_b = '+str(rb)+', rt = '+str(rt)+', h = '+str(h))

L = np.sqrt(h**2 + (rb-rt)**2)
l = (L*rt)/(rb-rt)
halfconeheight = h/2 

print('Warning: Using '+inputfile+', coneheight = '+str(halfconeheight*2)+', theta_p = '+str(anglep))
print('Output cone: coneheight = '+str(halfconeheight*2)+', theta_p = '+str(angleoutp))
print('----------')


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


#convert 3D Cartesian coords to 3D polars
polarcoordslist = []
#convert from Cartesians in 3D to polars in 3D - [x, y, z] to [phi, z, 0]
for particle in range(0, len(wantedcoordslist)):
	#phi = np.arctan(float(wantedcoordslist[particle][1])/float(wantedcoordslist[particle][0]))
	z = wantedcoordslist[particle][2]
	if ((float(wantedcoordslist[particle][0]) >= 0) and (float(wantedcoordslist[particle][1]) >= 0)): #top right quadrant
		phi = np.arctan(float(wantedcoordslist[particle][1])/float(wantedcoordslist[particle][0]))
	elif ((float(wantedcoordslist[particle][0]) < 0) and (float(wantedcoordslist[particle][1]) >= 0)): #top left quadrant
		phi = (np.arctan(float(wantedcoordslist[particle][1])/float(wantedcoordslist[particle][0])) + (np.pi))
	elif ((float(wantedcoordslist[particle][0]) < 0) and (float(wantedcoordslist[particle][1]) < 0)): #bottom left quadrant
		phi = np.arctan(float(wantedcoordslist[particle][1])/float(wantedcoordslist[particle][0])) + (np.pi)
	elif ((float(wantedcoordslist[particle][0]) >= 0) and (float(wantedcoordslist[particle][1]) < 0)): #bottom right quadrant
		phi = np.arctan(float(wantedcoordslist[particle][1])/float(wantedcoordslist[particle][0])) 
	polarcoordslist.append([float(phi), float(z), 0])
	

#convert 3D polars back to Cartesians but with new cone angle
outputcoordslist = []
#convert from polars in 3D to cartesians in 3D - [phi, z, 0] to [x, y, z]
for particle in range(0, len(polarcoordslist)):
	R = rb + (rt-rb)*(polarcoordslist[particle][1]/h + 0.5)
	Rnew = (R * np.tan(angleoutp/2))/(np.tan(anglep/2))
	#print(R, Rnew)
	x = Rnew * np.cos(polarcoordslist[particle][0])
	y = Rnew * np.sin(polarcoordslist[particle][0])
	z = polarcoordslist[particle][1]
	outputcoordslist.append([x, y, z])

#Constructing output file
outputfile = ''

for i in range(0, len(outputcoordslist)):
        output = '   ' + str(round(outputcoordslist[i][0], 10)) + '   ' + str(round(outputcoordslist[i][1], 10)) + '   ' + str(round(outputcoordslist[i][2], 10)) + '\n'
        outputfile = outputfile + output

print("coordsstretch file produced")

converted_file = open('coordsstretch', "w+")
converted_file.write(outputfile)
