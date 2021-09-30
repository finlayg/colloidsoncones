#Produces voronoi tesselation diagram for given minimum coordinates (takes net form)

from scipy.spatial import Voronoi, voronoi_plot_2d
from matplotlib import pyplot as plt
import numpy as np

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

printmin = input('Print just one minimum? Y/N ')
if printmin is 'Y':
	inputminimum = int(input('Minimum: '))
	inputminimarange = inputminimum
elif printmin is 'N':
	inputminimum = 1
	inputminimarange = int(input('Print minima up to number: '))

numparts = int(input('Number of particles (remember to add polygon if relevant): '))

whichangle = input('Input angle theta_p or theta_max? p/m ')
if whichangle is 'p':
	anglep = float(input('Angle theta_p of cone? '))*(np.pi/180) #angle of cone converted to radians
	angle = (2*np.pi*np.sin(anglep/2))#/(np.pi/180)
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

filename = input('Save as: ')

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
#	print(minimum, angle)
#	coordslist = []
#	n = 0
#	lowestfile = open(inputfile, 'r')
#		#numparts = int(stream1.readlines()[0])
#		#print(numparts, type(numparts))
#	for line in lowestfile:
#		if n < (float(minimum)*2*numparts) and 'C     ' in line:
#			splitline = line.rstrip().split(' ')
#			new = [item for item in splitline if item!='']
#			coordslist.append([new[1], new[2], new[3]])
#			n += 1
#	lowestfile.close()
#	print(len(coordslist))

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
			elif ((float(particle[0]) < 0) and (float(particle[1]) >= 0)): #top left quadrant
				theta = (np.arctan(float(particle[1])/float(particle[0])) + (np.pi))
			elif ((float(particle[0]) < 0) and (float(particle[1]) < 0)): #bottom left quadrant
				theta = np.arctan(float(particle[1])/float(particle[0])) + (np.pi)
			elif ((float(particle[0]) >= 0) and (float(particle[1]) < 0)): #bottom right quadrant
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
			elif ((float(particle[0]) < 0) and (float(particle[1]) < 0)): #bottom left quadrant
				theta = np.arctan(float(particle[1])/float(particle[0])) + (np.pi)
			elif ((float(particle[0]) >= 0) and (float(particle[1]) < 0)): #bottom right quadrant
				theta = np.arctan(float(particle[1])/float(particle[0])) 
			polarcoordslist.append([r, theta])
	
		#copy portion of net in both directions
		if angle < np.pi:
			netportion = 0.5
		if angle >= np.pi:
			netportion = (((2*np.pi) - angle)/2)/(2*np.pi)
			print('Angle > 180, fix this bit of the code')
		extpolarcoordslist = polarcoordslist.copy()
		for particle in polarcoordslist:
			if (particle[1] < (angle*netportion)):
				extpolarcoordslist.append([particle[0], (particle[1]+angle)])
			elif (particle[1] > (angle*(1-netportion))):
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
	#print(extendedcoordslist)

#loop to plot voronoi diagram of each minima in list of minima generated above
for minimum in range(0, len(minimalist)):
	vordiag = None
	vordiag = Voronoi(minimalist[minimum])
	
	#print(angle, l, L)

	fig, ax = plt.subplots()
	
	#plt.gca().clf()
	voronoi_plot_2d(vordiag, ax=ax, show_vertices=False, point_size=5) #ax=ax,  show_vertices=False, point_size = 5)
	ax.axis('off')
	
	if angle < (np.pi/12):
		#ax.set_xlim([(l-(L*(2/3))),(l+L+1)])
		#ax.set_xlim([(l-(L*(1/5))),(l+L-(L/4))])
		#ax.set_xlim([(l-(L*(1/3))),(l+L+1)])
		ax.set_xlim([(l-1),(l+L+1)])
		#ax.set_xlim([(l-(L*(2/3)))+15,(l+L+1)-20])
		#ax.set_xlim([0, 3*(L+l)])
		ax.set_ylim([(-1*(l+L+1)*np.sin(0.5*angle)),((l+L+1)*np.sin(1.5*angle))])
	elif angle < (np.pi/2):
		ax.set_xlim([(l*np.cos(angle)*(2/3))-1,(l+L+1)])
		ax.set_ylim([(-1*(l+L+1)*np.sin(0.5*angle)),((l+L+1)*np.sin(1.5*angle))])
	elif angle >= (np.pi/2):
		ax.set_xlim([(-1*(l+L+1)*np.sin(0.5*angle)),(l+L+1)])
		ax.set_ylim([(-1*(l+L+1)*np.sin(0.5*angle)),((l+L+1))])#*np.sin(angle))])
	
	if (h > 40) and (angle > (np.pi/12)):
		#ax.set_xlim([0,((l+L+1)*(2/3))])
		ax.set_xlim([0,(l+L+1)])
		ax.set_ylim([(-1*(l+L+1)*np.sin(0.5*angle)),(((l+L+1)*np.sin(1.5*angle))*(2/3))])
	
	#ax.set_xlim(0, 100)
	#ax.set_ylim(-50, 100)
	
	#-----------------------------------------------------------
	
	lw = '1' #linewidth
	lc = 'k' #colour of net
	
	#Adding solid boundary of cone net at theta = 0
	rrange = np.arange(l, (l+L), 0.01)
	thetarange = rrange*0
	xrange = []
	yrange = []
	for i in range(0, len(rrange)):
		xrange.append(rrange[i]*np.cos(thetarange[i]))
		yrange.append(rrange[i]*np.sin(thetarange[i]))
	ax.plot(xrange, yrange, color=lc, linewidth = lw, label="Outline of net of cone")
	
	#Adding solid boundary of cone net at theta = coneangle
	rrange2 = np.arange(l, (l+L), 0.01)
	thetarange2 = (rrange2*0)+angle
	xrange2 = []
	yrange2 = []
	for j in range(0, len(rrange2)):
		xrange2.append(rrange2[j]*np.cos(thetarange2[j]))
		yrange2.append(rrange2[j]*np.sin(thetarange2[j]))
	ax.plot(xrange2, yrange2, linewidth = lw, color=lc)
	
	#Adding dotted boundary of cone net at theta = 0
	rrange3 = np.arange(0, l, 0.01)
	thetarange3 = rrange3*0
	xrange3 = []
	yrange3 = []
	for i in range(0, len(rrange3)):
		xrange3.append(rrange3[i]*np.cos(thetarange3[i]))
		yrange3.append(rrange3[i]*np.sin(thetarange3[i]))
	ax.plot(xrange3, yrange3, ":", color=lc)
	
	#Adding dotted boundary of cone net at theta = coneangle 
	rrange4 = np.arange(0, l, 0.01)
	thetarange4 = (rrange4*0)+angle
	xrange4 = []
	yrange4 = []
	for i in range(0, len(rrange4)):
		xrange4.append(rrange4[i]*np.cos(thetarange4[i]))
		yrange4.append(rrange4[i]*np.sin(thetarange4[i]))
	ax.plot(xrange4, yrange4, ":", color=lc)
	
	#Adding net arc at bottom of cone
	thetarange5 = np.arange(0, angle+0.001, 0.001)
	rrange5 = (thetarange5*0)+L+l
	xrange5 = []
	yrange5 = []
	for i in range(0, len(rrange5)):
		xrange5.append(rrange5[i]*np.cos(thetarange5[i]))
		yrange5.append(rrange5[i]*np.sin(thetarange5[i]))
	ax.plot(xrange5, yrange5, linewidth = lw, color=lc)
	
	#Adding net arc at top of cone
	thetarange6 = np.arange(0, angle+0.001, 0.001)
	rrange6 = (thetarange6*0)+l
	xrange6 = []
	yrange6 = []
	for i in range(0, len(rrange6)):
		xrange6.append(rrange6[i]*np.cos(thetarange6[i]))
		yrange6.append(rrange6[i]*np.sin(thetarange6[i]))
	ax.plot(xrange6, yrange6, linewidth = lw, color=lc)
	
	ax.set_aspect('equal')
	#set_aspect('equal')
	#-----------------------------------------------------------
	
	#Adding colour
	
	#irreglist = []
	#reglist = []
	
	if colour is 'Y': #if colour is to be added
		n=0
		m=0
		q=0
		k=0
		g=0
		for region in range(0, len(vordiag.regions)): #for each polygon/region on the diagram
			regionlengths = []
			distortedpolygon = False
			notlargepolygon = True
			if not -1 in vordiag.regions[region]: #if the region is bounded
	
	#generating length list containing lengths of each side in given polygon
				for polygonindex in range(0, len(vordiag.regions[region])-1):                                           #for each point in the region (not including the last one)
					for polygonindex2 in range(polygonindex+1, len(vordiag.regions[region])):                       #and for every other point after it in the region...
						for lineindex in range(0, len(vordiag.ridge_vertices)):                #for each line in the voronoi diagram...
							if (vordiag.regions[region][polygonindex] in vordiag.ridge_vertices[lineindex]) and (vordiag.regions[region][polygonindex2] in vordiag.ridge_vertices[lineindex]):                                                                                                                      #...find if that line contains these two vertices...
									                                               #...and if so then find length of the line and store it
								sidelength = np.sqrt((vordiag.vertices[vordiag.ridge_vertices[lineindex][1]][0] - \
	                                                                              vordiag.vertices[vordiag.ridge_vertices[lineindex][0]][0])**2 + \
	                                                                             (vordiag.vertices[vordiag.ridge_vertices[lineindex][1]][1] - \
	                                                                              vordiag.vertices[vordiag.ridge_vertices[lineindex][0]][1])**2)
								regionlengths.append(sidelength)
	
	
				if (len(regionlengths) > 2): 
					thresholdlength = (min(regionlengths))*lengthtolerance
					for length in range(1, len(regionlengths)):                                    #for each side in the polygon
						if (abs(regionlengths[length] - regionlengths[0]) > 1.6):
							notlargepolygon = False
	
	#colouring diagram based on lengths in length list
	#colouring 6 sided polygons, regular and irregular differently
				if (regionlengths != []) and (len(vordiag.regions[region]) == 6):
					thresholdlength = (min(regionlengths))*lengthtolerance
					for length in range(1, len(regionlengths)):                                    #for each side in the polygon
						if (abs(regionlengths[length] - regionlengths[0]) > thresholdlength):  #if this side is not same length as the first within the given tolerance
							distortedpolygon = True
		
					if notlargepolygon:
						if distortedpolygon:
							polygon = [vordiag.vertices[i] for i in vordiag.regions[region]]
							if m<1:
								ax.fill(*zip(*polygon), "gold", label="Irregular hexagonal environment")
								m=1
							else:
								ax.fill(*zip(*polygon), "gold")
						else:	
							polygon = [vordiag.vertices[i] for i in vordiag.regions[region]]
							#reglist.append(np.where(vordiag.point_region == region)[0][0])
							if n<1:
								ax.fill(*zip(*polygon), "greenyellow", label="Regular hexagonal environment")
								n=1
							else:
								ax.fill(*zip(*polygon), "greenyellow")
	
	#colouring polygons of 5 and 7 and other length sides
				elif (regionlengths != []) and (len(vordiag.regions[region]) == 5) and (notlargepolygon):
					polygon = [vordiag.vertices[i] for i in vordiag.regions[region]]
					if k<1:
						ax.fill(*zip(*polygon), "r", label="Pentagonal environment")
						k = 1
					else:
						ax.fill(*zip(*polygon), "r")
				elif (regionlengths != []) and (len(vordiag.regions[region]) == 7) and (notlargepolygon):
					polygon = [vordiag.vertices[i] for i in vordiag.regions[region]]
					if q<1:
						ax.fill(*zip(*polygon), "b", label="Heptagonal environment")
						q=1
					else:
						ax.fill(*zip(*polygon), "b")
				elif (regionlengths != []) and (len(vordiag.regions[region]) in [3, 4, 8, 9]) and (notlargepolygon):
					polygon = [vordiag.vertices[i] for i in vordiag.regions[region]]
					if g<1:
						ax.fill(*zip(*polygon), "orange", label="3/4/5/9 sided environment")
						g=1
					else:
						ax.fill(*zip(*polygon), "orange")
	
	#if angle < (np.pi/12):
	#	plt.legend(fontsize='x-small', framealpha = 1, loc = "upper left")	
	#else:
	#	plt.legend(fontsize='x-small', framealpha = 1)#loc="lower left")
	
	#print(len(reglist), reglist)
	
	fig.savefig(filename+'_min'+str(minimum+1)+'.png', dpi = 300)
	print('Minimum '+str(minimum+1)+' saved')
	fig.clf()
	plt.close()
