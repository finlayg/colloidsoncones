#script to identify type of defect from a given lowest input file

from scipy.spatial import Voronoi, voronoi_plot_2d
from matplotlib import pyplot as plt
import numpy as np
from collections import Counter
from defectidentifierfunction import defectcategorisation
from defectidentification import defectidentification

standardnames = input('Files are called lowest and data? Y/N ')
if standardnames is 'N':
	inputfile = input('Lowest Filename: ')
	datafile = input('Data filename for dimensions: ')
elif standardnames is 'Y':
	inputfile = 'lowest'
	datafile = 'data'

filetree = input('Location of lowest and data: ')

datafile = filetree+'/'+datafile
inputfile = filetree+'/'+inputfile

categoriseall = input('Categorise multiple minima in lowest and produce a frequency vs energy histogram? Y/N ')
if categoriseall is 'N': 
	inputminimum = int(input('Single minimum to categorise: '))
	inputminimarange = inputminimum
elif categoriseall is 'Y':
	inputminimum = 1
	inputminimarange = int(input('Categorise and plot energy vs frequency histograms using minima up to number: '))
else:
	print('Error: input not recognised')


numparts = int(input('Number of particles (remember to add polygon if relevant): '))

whichangle = input('Input angle theta_p or theta_max? p/m ')
if whichangle is 'p':
	anglep = float(input('Angle theta_p of cone? '))#*(np.pi/180) 	#angle p of cone
	anglepr = anglep*(np.pi/180)					#converted to radians
	angle = (2*np.pi*np.sin(anglepr/2))#/(np.pi/180)
elif whichangle is 'm':
	angle = float(input('Angle theta_max of cone? '))*(np.pi/180) #angle of net, converted to radians
else:
	print('Error: input not recognised')

dimensionlist = []
with open(datafile, 'r') as stream2:
	for line in stream2:
		if 'GTHOMSON ' in line:
			splitline = line.rstrip().split(' ')
			new = [item for item in splitline if item!=' ']
			rb = float(new[2])
			rt = float(new[3])
			h = float(new[4])
			print('Dimensions: rb = '+str(rb)+', rt = '+str(rt)+' , h = '+str(h)+' , Angle of net in degrees = '+str(angle/(np.pi/180)))

#first finding L and l to define sides of net 
L = np.sqrt(h**2 + (rb-rt)**2)
l = (L*rt)/(rb-rt)

rotate = input('Rotate net by a given angle? Y/N ')
if rotate is 'Y':
	angofrot = float(input('Angle of rotation in degrees: '))*(np.pi/180) #angle of rotation in +ve angle direction converted to radians

extend = 'Y' #input('Extend net for visualisation? Y/N ')
colour = 'Y' #input('Add colour? Y/N ')
lengthtolerance = float(0.3) #float(input('Tolerance on length for colour: '))

coordslist = []
n = 0
lowestfile = open(inputfile, 'r')
for line in lowestfile:
	if n < (float(inputminimarange)*2*numparts) and 'C     ' in line:
		splitline = line.rstrip().split(' ')
		new = [item for item in splitline if item!='']
		coordslist.append([new[1], new[2], new[3]])
		n += 1
lowestfile.close()
#print(len(coordslist))

minimalist = []

#loop to generate list of minima coords to plot voronoi diagrams of
for minimum in range(inputminimum, inputminimarange+1):
	#second set of points for each minima in lowest is the net of cone

	wantedcoordslist = []

	for i in range(((numparts*minimum*2) - (numparts)), (numparts*minimum*2)):
		wantedcoordslist.append([coordslist[i][0], coordslist[i][1]])
		#wantedcoordslist = coordslist[(-numparts):]
	#print((wantedcoordslist))

	#Rotating net by a given input angle if desired
	if rotate is 'Y':
		#convert coords to polars
		polarcoordslist = []
		for particle in wantedcoordslist:
			r = np.sqrt(float(particle[0])**2+float(particle[1])**2)
			if ((float(particle[0]) >= 0) and (float(particle[1]) >= 0)): #top right quadrant
				theta = np.arctan(float(particle[1])/float(particle[0]))
			elif ((float(particle[0]) < 0) and (float(particle[1]) > 0)): #top left quadrant
				theta = (np.arctan(float(particle[1])/float(particle[0])) + (np.pi))
			elif ((float(particle[0]) < 0) and (float(particle[1]) < 0)): #bottom left quadrant
				theta = np.arctan(float(particle[1])/float(particle[0])) + (np.pi)
			elif ((float(particle[0]) < 0) and (float(particle[1]) < 0)): #bottom right quadrant
				theta = np.arctan(float(particle[1])/float(particle[0])) 
			polarcoordslist.append([r, theta])

		#rotate all particles by a given angle
		rotatedcoordslist = []
		for particle in polarcoordslist:
			if (particle[1]+angofrot > angle):
				rotatedcoordslist.append([particle[0], (particle[1]+angofrot-angle)])
			elif (particle[1]+angofrot < angle):
				rotatedcoordslist.append([particle[0], (particle[1]+angofrot)])
		#convert back to cartesians
		wantedcoordslist = []
		for particle in rotatedcoordslist:
			x = particle[0]*np.cos(particle[1])
			y = particle[0]*np.sin(particle[1])
			wantedcoordslist.append([x, y])
	
	#Extending net periodically
	if extend is 'Y':
		#convert coords to polars
		polarcoordslist = []
		for particle in wantedcoordslist:
			r = np.sqrt(float(particle[0])**2+float(particle[1])**2)
			if ((float(particle[0]) >= 0) and (float(particle[1]) >= 0)): #top right quadrant
				theta = np.arctan(float(particle[1])/float(particle[0]))
			elif ((float(particle[0]) < 0) and (float(particle[1]) >= 0)): #top left quadrant
				theta = (np.arctan(float(particle[1])/float(particle[0])) + (np.pi))
			elif ((float(particle[0]) >= 0) and (float(particle[1]) < 0)): #bottom left quadrant
				theta = np.arctan(float(particle[1])/float(particle[0])) + (np.pi)
			elif ((float(particle[0]) < 0) and (float(particle[1]) < 0)): #bottom right quadrant
				theta = np.arctan(float(particle[1])/float(particle[0])) 
			polarcoordslist.append([r, theta])
	
		#copy portion of net in both directions
		if angle < np.pi:
			netportion = 1
		if angle >= np.pi:
			netportion = (((2*np.pi) - angle)/2)/(2*np.pi)
			print('Angle > 180, fix this bit of the code')
		extpolarcoordslist = polarcoordslist.copy()
		for particle in polarcoordslist:						#adding top copy of net
			if (particle[1] < (angle*netportion)):
				extpolarcoordslist.append([particle[0], (particle[1]+angle)])
		for particle in polarcoordslist:						#adding bottom copy of net
			if (particle[1] > (angle*(1-netportion))):
				extpolarcoordslist.append([particle[0], (particle[1]-angle)])
	
		#convert back to cartesians
		extendedcoordslist = []
		for particle in extpolarcoordslist:
			x = particle[0]*np.cos(particle[1])
			y = particle[0]*np.sin(particle[1])
			extendedcoordslist.append([x, y])
	else:
		extendedcoordslist = wantedcoordslist.copy()

	minimalist.append(extendedcoordslist)
	#print('extcoordslist length:'+str(len(extendedcoordslist)))

defectidvariablelist = defectidentification(minimalist, l, L, angle, colour, lengthtolerance, anglep, categoriseall, inputminimum)
