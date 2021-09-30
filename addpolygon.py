#takes lowest file and prints coordinates of vertices to use in polygon potential

import numpy as np

inputfile = input('Lowest input filename: ')
datafile = input('Data filename for dimensions: ')
#inputminima = float(input('Minima coordinates to convert: '))
#numparts = int(input('Number of particles: '))
polygonvert = int(input('Number of vertices on polygon potential: '))


coordslist = []
n = 0

#pulling dimensions from data file
with open(datafile, 'r') as stream2:
	for line in stream2:
		if 'GTHOMSON ' in line:
			splitline = line.rstrip().split(' ')
			new = [item for item in splitline if item!=' ']
			rb = float(new[2])
			rt = float(new[3])
			h = float(new[4])
			
			L = np.sqrt(h**2 + (rb-rt)**2)
			l = (L*rt)/(rb-rt)
		
			print('Dimensions: r_b = '+str(rb)+', r_t = '+str(rt)+', h = '+str(h))

#Setting coordinates of polygpot
polygcart = []
polygnet = []

alpha = np.arccos(1)
for vert in range(0, polygonvert):
	polygcart.append([0, 0, 0])
	polygcart[vert][0] = rt * np.cos(((2*np.pi) / polygonvert) * vert)
	polygcart[vert][1] = rt * np.sin(((2*np.pi) / polygonvert) * vert)
	polygcart[vert][2] = h/2
	polygnet.append([0, 0, 0])

#	r=l+L/2(1-cosalpha), theta=phi (rb-rt)/L
#	x=rcostheta y=rsintheta
#	so 
#	x=(l+L/2(1-cosalpha))*cos(phi(rb-rt/L))
#	y=(l+L/2(1-cosalpha))*sin(phi(rb-rt/L))
#	where
#	phi=((2*np.pi)/polygonvert)*vert
#	alpha = arccos(2*z/H)

	polygnet[vert][0] = (l + (L/2)*(1-np.cos(alpha))) * np.cos((((2*np.pi)/polygonvert)*vert)*((rb-rt)/L))
	polygnet[vert][1] = (l + (L/2)*(1-np.cos(alpha))) * np.sin((((2*np.pi)/polygonvert)*vert)*((rb-rt)/L))
	polygnet[vert][2] = 0.0

linelist = []
with open(inputfile, 'r') as stream1:
	for line in stream1:
		linelist.append(str(line))

numparts = int(linelist[0])
newnumparts = numparts+polygonvert

for i in range(0, len(linelist)):
	if ('         '+str(numparts)+'\n') in linelist[i]:
		linelist[i] = '         '+str(newnumparts)+'\n'
	elif 'Energy ' in linelist[i]:
		for j in range(0, polygonvert):
			linelist.insert(i+1, 'C         '+str(polygcart[j][0])+'   '+str(polygcart[j][1])+'   '+str(polygcart[j][2])+'\n')
	elif 'previous minimum transformed' in linelist[i]:
		for j in range(0, polygonvert):
			linelist.insert(i+1, 'C          '+str(polygnet[j][0])+'   '+str(polygnet[j][1])+'   '+str(polygnet[j][2])+'\n')

#second set of points for each minima in lowest is the net of cone - don't want this
wantedcoordslist = coordslist[(-2*numparts):(-1*numparts)]

#Constructing output file
outputfile = ''

for i in range(0, len(linelist)):
        output = linelist[i]# + '\n'
        outputfile = outputfile + output

print("lowestwithpolygon file produced")

converted_file = open('lowestwithpolygon', "w+")
converted_file.write(outputfile)
