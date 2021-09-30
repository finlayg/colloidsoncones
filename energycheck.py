#script to check energy of system (using morse and boundary potentials) of output configs in lowest

import numpy as np
from matplotlib import pyplot as plt
from scipy.spatial import Voronoi

#inputfile = input('Lowest filename: ') 
#datafile = input('Data filename for dimensions: ')
#inputminima = int(input('Which minima config to check? '))
#numparts = int(input('Number of particles: '))


#for debugging speed
inputfile = 'lowest'
datafile = 'data'
numparts = 200
inputminima = 1
bptest = 'N'
particlecalc = 'N'
plothist = 'Y'
anglep = 5.715*(np.pi/180)
angle = (2*np.pi*np.sin(anglep/2))
findlineandboundary = 'Y' 
lengthtolerance = float(0.3)

coordslist = []
n = 0
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

#defining potential parameters
sigma = 1.0
rho = 18 
zcutoff = 1
halfconeheight = h/2 
reppower = -12
print('Warning: Using '+inputfile+', rho = '+str(rho)+', sigma = '+str(sigma)+', coneheight = '+str(halfconeheight*2)+', theta_p = '+str(anglep))
print('----------')

with open(inputfile, 'r') as stream1:
	for line in stream1:
                if n < (inputminima*2*numparts) and 'C     ' in line:
                        splitline = line.rstrip().split(' ')
                        new = [item for item in splitline if item!='']
                        coordslist.append([float(new[1]), float(new[2]), float(new[3])])
                        n += 1
#second set of points for each minima in lowest is the net of cone - don't want this

netcoordslist = []

for i in range(-numparts, 0):
	netcoordslist.append([coordslist[i][0], coordslist[i][1]])

wantedcoordslist = coordslist[(-2*numparts):(-1*numparts)]

#print((netcoordslist))

print('Using newbp')

#-------------------------------
#FOR BPTEST USING Z VALUES TO FIND RANGE OF ZS FROM TOP TO BOTTOM OF CONE AND WHERE BP ACTS
#-------------------------------
#bptest = input('Find range of z that boundary potential acts over? Y/N ')
if bptest is 'Y':
	numvals = 1000
	zlist = []
	intervalno = 0
	rangeinterval = (halfconeheight*2)/numvals
	for i in range(0, numvals):
		zlist.append((-halfconeheight + (intervalno*rangeinterval)+0.01))
		#print(intervalno)
		intervalno = intervalno + 1
	#print(zlist)

	#adding boundary E if required
	for j in range(0, len(zlist)):
		if (zlist[j] > (halfconeheight - zcutoff)): #top of cone 
			boundaryE = (0.001*(-1*zlist[j] + halfconeheight)**(reppower) - 0.001*(-1*zcutoff + halfconeheight)**(reppower) - \
			    	    (-1*zlist[j] + halfconeheight + zcutoff)*(-1*reppower*0.001*(-1*zcutoff + halfconeheight)**(reppower-1)))
		elif (zlist[j] < (zcutoff - halfconeheight)): #bottom of cone
			boundaryE = (0.001*(   zlist[j] + halfconeheight)**(reppower) - 0.001*(   zcutoff + halfconeheight)**(reppower) - \
				    (   zlist[j] + halfconeheight - zcutoff)*(   reppower*0.001*(   zcutoff + halfconeheight)**(reppower-1)))
		else:
			boundaryE = 0
	
		if boundaryE != 0:
			print('Adding bp for z value '+str(j)+' at '+str(zlist[j]))
			print(str(boundaryE))
	print('zcutoff is above '+str(halfconeheight-zcutoff)+' and below '+str(zcutoff-halfconeheight))
	exit()

#--------------------------------
#CALCULATING OVERALL ENERGY OF SYSTEM
#--------------------------------
totalE = 0
for  p1index in range(0, len(wantedcoordslist)):
	#for all particles adding boundary energy and (for up to last) Morse energy with all subsequent particles
#	print('Particle '+str(p1index))	
	p1coords = wantedcoordslist[p1index]	
	
		
	#adding boundary E if required
	if (p1coords[2] > (halfconeheight - zcutoff)): #top of cone 
		boundaryE = (0.001*(-1*p1coords[2] + halfconeheight)**(reppower) - 0.001*(-1*zcutoff + halfconeheight)**(reppower) - \
			    (-1*p1coords[2] + halfconeheight + zcutoff)*(-1*reppower*0.001*(-1*zcutoff + halfconeheight)**(reppower-1)))
	elif (p1coords[2] < (zcutoff - halfconeheight)): #bottom of cone
		boundaryE = (0.001*(   p1coords[2] + halfconeheight)**(reppower) - 0.001*(   zcutoff + halfconeheight)**(reppower) - \
			    (   p1coords[2] + halfconeheight - zcutoff)*(   reppower*0.001*(   zcutoff + halfconeheight)**(reppower-1)))
	else:
		boundaryE = 0
	
	#if boundaryE != 0:
		#print('Adding bp for particle '+str(p1index))
		#print(str(boundaryE))

	
	#adding Morse energy for all particles except the last (all contributions already calculated)
	if p1index != (len(wantedcoordslist)-1):
		totalmorseE = 0
#		print('pair '+str(p1index))	
		for p2index in range (p1index+1, len(wantedcoordslist)):
			p2coords = wantedcoordslist[p2index]
#			print(p2index)

			#1 find distance between p1 and p2
			dist = np.sqrt((p2coords[2]-p1coords[2])**2 + (p2coords[1]-p1coords[1])**2 + (p2coords[0]-p1coords[0])**2)
		
			#2 calculate morse E given by presence of p2
			morseE = ((1 - np.exp(rho*(sigma-dist)))**2 - 1)
			totalmorseE = totalmorseE + morseE
		
	totalE = totalE + totalmorseE + boundaryE

print('Total energy of this configuration: ' + str(totalE))

#-------------------------------------
#FINDING ENERGIES OF INDIVIDUAL PARTICLES AND PLOTTING HISTOGRAM
#-------------------------------------
#particlecalc = input('Print file of individual particle energies? Y/N ')
#plothist = input('Plot histogram of particle energies? Y/N ')

outputfile = ''
particleenergies = []

for p1index in range(0, len(wantedcoordslist)):
	#print(p1index)
	p1coords = wantedcoordslist[p1index]	
	
	#adding boundary E if required
	if (p1coords[2] > (halfconeheight - zcutoff)): #top of cone 
		boundaryE = (0.001*(-1*p1coords[2] + halfconeheight)**(reppower) - 0.001*(-1*zcutoff + halfconeheight)**(reppower) - \
			    (-1*p1coords[2] + halfconeheight + zcutoff)*(-1*reppower*0.001*(-1*zcutoff + halfconeheight)**(reppower-1)))
	elif (p1coords[2] < (zcutoff - halfconeheight)): #bottom of cone
		boundaryE = (0.001*(   p1coords[2] + halfconeheight)**(reppower) - 0.001*(   zcutoff + halfconeheight)**(reppower) - \
			    (   p1coords[2] + halfconeheight - zcutoff)*(   reppower*0.001*(   zcutoff + halfconeheight)**(reppower-1)))
	else:
		boundaryE = 0
	
	#adding Morse energy for all particles except the last (all contributions already calculated)
	totalmorseE = 0
	for p2index in range(0, len(wantedcoordslist)):
		if p2index != p1index:
			#print(p2index)
			p2coords = wantedcoordslist[p2index]

			#1 find distance between p1 and p2
			dist = np.sqrt((p2coords[2]-p1coords[2])**2 + (p2coords[1]-p1coords[1])**2 + (p2coords[0]-p1coords[0])**2)
	
			#2 calculate morse E given by presence of p2
			morseE = ((1 - np.exp(rho*(sigma-dist)))**2 - 1)
			totalmorseE = totalmorseE + morseE
		
	particleenergy = totalmorseE + boundaryE
	particleenergies.append(particleenergy)
#	output = str(p1index)+' , '+str(particleenergy)+'\n'
#	outputfile = outputfile + output


if particlecalc is 'Y': 
#	print("particleenergies.txt file produced")
	converted_file = open('particleenergies.txt', "w+")
	converted_file.write(outputfile)

if plothist is 'Y':
	newparticleenergies = []
	for i in particleenergies:
		if i < 0:
			newparticleenergies.append(i)
	fig, ax = plt.subplots()
	numbins = 24 
	ax.hist(newparticleenergies, range = [-6, 0], bins=numbins)
	plt.xlabel('Particle Energy / $\epsilon$')
	plt.ylabel('Number of Particles')
	plt.ylim(0, 160)

	filename = 'partenergieshistogram'
	fig.savefig(filename+'.png', dpi = 300)

#----------------------------------------
#TAKING ANGLE OF CONE INPUT AND EXTENDING NET OF CONE UP AND DOWN
#----------------------------------------

#whichangle = input('Input angle theta_p or theta_max? p/m ')
#if whichangle is 'p':
#	anglep = float(input('Angle theta_p of cone? '))*(np.pi/180) #angle of cone converted to radians
#	angle = (2*np.pi*np.sin(anglep/2))#/(np.pi/180)
#elif whichangle is 'm':
#	angle = float(input('Angle theta_max of cone? '))*(np.pi/180) #angle of net, converted to radians
#else:
#	print('Error: input not recognised')

#Extending net periodically
#convert coords to polars
polarcoordslist = []
for particle in netcoordslist:
	r = np.sqrt(float(particle[0])**2+float(particle[1])**2)
	if ((float(particle[0]) > 0) and (float(particle[1]) > 0)): #top right quadrant
		theta = np.arctan(float(particle[1])/float(particle[0]))
	elif ((float(particle[0]) < 0) and (float(particle[1]) > 0)): #top left quadrant
		theta = (np.arctan(float(particle[1])/float(particle[0])) + (np.pi))
#		print(particle[0], particle[1], 'adding 90 - >', r, theta)
	elif ((float(particle[0]) < 0) and (float(particle[1]) < 0)): #bottom left quadrant
		theta = np.arctan(float(particle[1])/float(particle[0])) + (np.pi)
	elif ((float(particle[0]) < 0) and (float(particle[1]) < 0)): #bottom right quadrant
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

extendedcoordsxlist = []
extendedcoordsylist = []
for i in range(0, len(extendedcoordslist)):
	extendedcoordsxlist.append(extendedcoordslist[i][0])
	extendedcoordsylist.append(extendedcoordslist[i][1])

#---------------------------------------------------------------
#IDENTIFYING DEFECT AND BOUNDARY PARTICLES USING VORONOI REGIONS: Irregular regions are defect and boundary particles
#---------------------------------------------------------------

vordiag = Voronoi(extendedcoordslist)

#List which will contain the indices of all particles with irregular polygons.
irregpolygpartlist = []
regpolygpartlist = []

#findlineandboundary = input('Find line and boundary particles? Y/N ')

if findlineandboundary is 'Y': #if colour is to be added
#	lengthtolerance = float(input('Tolerance on length for irregular hexagons: '))
	for region in range(0, len(vordiag.regions)): #for each polygon/region on the diagram
		regionlengths = []
		distortedpolygon = False
		if -1 in vordiag.regions[region]: #if region is unbounded
			irregpolygpartlist.append(np.where(vordiag.point_region == region)[0][0])
		elif not -1 in vordiag.regions[region]: #if the region is bounded

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

#Assigning different types of polygon

#6 sided polygons, regular and irregular
			if (regionlengths != []) and (len(vordiag.regions[region]) == 6):
				thresholdlength = (min(regionlengths))*lengthtolerance
				for length in range(1, len(regionlengths)):                                    #for each side in the polygon
					if (abs(regionlengths[length] - regionlengths[0]) > thresholdlength):  #if this side is not same length as the first within the given tolerance
						distortedpolygon = True
				if distortedpolygon:
					polygon = [vordiag.vertices[i] for i in vordiag.regions[region]]
					#plt.fill(*zip(*polygon), "gold", label="Irregular hexagonal environment")
					irregpolygpartlist.append(np.where(vordiag.point_region == region)[0][0])
				else:	
					polygon = [vordiag.vertices[i] for i in vordiag.regions[region]]
					#plt.fill(*zip(*polygon), "greenyellow", label="Regular hexagonal environment")
					regpolygpartlist.append(np.where(vordiag.point_region == region)[0][0])

#polygons of 5 and 7 and other length sides
			elif (regionlengths != []) and (len(vordiag.regions[region]) == 5):
				#polygon = [vordiag.vertices[i] for i in vordiag.regions[region]]
				#plt.fill(*zip(*polygon), "r", label="Pentagonal environment")
				irregpolygpartlist.append(np.where(vordiag.point_region == region)[0][0])
			
			elif (regionlengths != []) and (len(vordiag.regions[region]) == 7):
				#polygon = [vordiag.vertices[i] for i in vordiag.regions[region]]
				#plt.fill(*zip(*polygon), "b", label="Heptagonal environment")
				irregpolygpartlist.append(np.where(vordiag.point_region == region)[0][0])
			
			elif (regionlengths != []) and (len(vordiag.regions[region]) in [3, 4, 8, 9]):
				#polygon = [vordiag.vertices[i] for i in vordiag.regions[region]]
				#plt.fill(*zip(*polygon), "orange", label="3/4/5/9 sided environment")
				irregpolygpartlist.append(np.where(vordiag.point_region == region)[0][0])

#print(len(irregpolygpartlist), irregpolygpartlist)
#print(len(regpolygpartlist), regpolygpartlist)
fulllist = irregpolygpartlist + regpolygpartlist

#print(vordiag.points[irregpolygpartlist[0]])

xlist = []
ylist = []
edgexlist = []
edgeylist = []
defectxlist = []
defectylist = []
allpartsxlist = []
allpartsylist = []
boundedgelistx = []
boundedgelisty = []
edgeindlist = []
defindlist = []

for i in range(0, numparts):
	allpartsxlist.append(vordiag.points[i][0])
	allpartsylist.append(vordiag.points[i][1])

for i in range(0, len(irregpolygpartlist)):
	boundedgelistx.append(vordiag.points[irregpolygpartlist[i]][0])
	boundedgelisty.append(vordiag.points[irregpolygpartlist[i]][1])

allpartslist = []
for i in range(0, len(allpartsxlist)):
	allpartslist.append([allpartsxlist[i], allpartsylist[i]])

#---------------------------------------------
#DIFFERENTIATING BETWEEN EDGE AND DEFECT PARTICLES
#---------------------------------------------

#1) choose a particle in edge or boundary list

for partindind1 in range(0, len(irregpolygpartlist)):
	x1 = vordiag.points[irregpolygpartlist[partindind1]][0]
	y1 = vordiag.points[irregpolygpartlist[partindind1]][1]
	count = 0

#2) cycle through all particles on edge or boundary, find those within certain distance (that border it)
	
	for partindind2 in range(0, len(irregpolygpartlist)):
		x2 = vordiag.points[irregpolygpartlist[partindind2]][0]
		y2 = vordiag.points[irregpolygpartlist[partindind2]][1]
		dist12 = np.sqrt((x1-x2)**2 + (y1-y2)**2)
		if ((dist12 < 1.1) and (dist12 > 0.001)): #then we have a pair
#			print(x1, y1, x2, y2)
			#find two positions at either side of the pair, p3 and p4

			pairm = (y2-y1)/(x2-x1)
			pairc = (y1 - (((y2-y1)/(x2-x1))*x1))
			midpointx = x1+0.5*(x2-x1)
			midpointy = y1+0.5*(y2-y1)
			perpm = -1*(1/pairm)
			perpc = midpointy - (perpm*midpointx)
			vecx = midpointx
			vecy = midpointy - perpc
			veclength = np.sqrt((midpointy - perpc)**2+(midpointx**2))

			x3 = midpointx + 0.85*(vecx/veclength)
			y3 = midpointy + 0.85*(vecy/veclength)
			x4 = midpointx - 0.85*(vecx/veclength)
			y4 = midpointy - 0.85*(vecy/veclength)
			
			xlist.append(x3)
			xlist.append(x4)
			ylist.append(y3)
			ylist.append(y4)
			
			for bulkpind in range(0, len(extendedcoordslist)): #for all particles in extended bulk
				x5 = extendedcoordslist[bulkpind][0]
				y5 = extendedcoordslist[bulkpind][1]
				if ((x1 != x5) and (x2 != x5) and (y1 != y5) and (y2 != y5)):
					dist35 = np.sqrt((x3-x5)**2 + (y3-y5)**2)
					dist45 = np.sqrt((x4-x5)**2 + (y4-y5)**2)
					if ((dist35 < 0.9) or (dist45 < 0.9)):
						count += 1
#						print('Adding to count')
		if count is 1: #then edge particle
			edgexlist.append(x1)#vordiag.points[irregpolygpartlist[partindind1]][0])
			edgeylist.append(y1)#vordiag.points[irregpolygpartlist[partindind1]][1])
			if (irregpolygpartlist[partindind1]) not in edgeindlist:
				edgeindlist.append(irregpolygpartlist[partindind1])


#3) choose each pair of particles and cycle through all particles on cone; if one is within certain distance of both then boundary, two then defect


#Next : find indices of edge and boundary parts in two lists. Use these to match up with individual particle energies
# then average over each (weighted for amount of each per length) to get energy per length

for i in irregpolygpartlist:
	if i not in edgeindlist:
		defindlist.append(i)

edgepartenergies = []
defpartenergies = []
edgelistx = []
edgelisty = []
deflistx = []
deflisty = []

for i in range(0, len(edgeindlist)):
	if edgeindlist[i] < numparts:
		edgepartenergies.append(particleenergies[edgeindlist[i]])
		edgelistx.append(vordiag.points[edgeindlist[i]][0])
		edgelisty.append(vordiag.points[edgeindlist[i]][1])

for i in range(0, len(defindlist)):
	if defindlist[i] < numparts:
		defpartenergies.append(particleenergies[defindlist[i]])
		deflistx.append(vordiag.points[defindlist[i]][0])
		deflisty.append(vordiag.points[defindlist[i]][1])

fig1, ax1 = plt.subplots()
ax1.plot(allpartsxlist, allpartsylist, marker = 'o', linestyle = 'none')
#ax1.plot(extendedcoordsxlist, extendedcoordsylist, marker = '^', linestyle = 'none')
ax1.plot(deflistx, deflisty, marker = '*', linestyle = 'none')
ax1.plot(edgelistx, edgelisty, marker = '.', linestyle = 'none')
ax1.plot(xlist, ylist, marker = ',', linestyle = 'none')
plt.xlabel('x')
plt.ylabel('y')
ax1.set_aspect('equal')

filename = 'testvis'
fig1.savefig(filename+'.png')


#print(len(allpartsxlist), len(edgelistx), len(deflistx), (edgepartenergies), (defpartenergies))
#print(allpartsxlist, allpartsylist)

#---------------------------------
#FINDING LENGTH OF EDGE AND DEFECT TO ALLOW ENERGY PER LENGTH CALCULATIONS
#---------------------------------

#find halfway point in x on net
midx = (min(allpartsxlist)+max(allpartsxlist))/2
#print(midx)

#DEFECT FIRST
#find coords of max x and min x particles in defect
defectL = deflistx.index(min(deflistx))
defectR = deflistx.index(max(deflistx))

linedeflength = np.sqrt((deflistx[defectL] - deflistx[defectR])**2 + (deflisty[defectL] - deflisty[defectR])**2)
#print(linedeflength)

defeneperpart = sum(defpartenergies) / len(defpartenergies)
defectenergy = sum(defpartenergies) / (linedeflength)
print('-----DEFECT-----')
print('Defect Energy per Length: '+str(defectenergy)+', per Particle: '+str(defeneperpart),', gain per Particle: ',6+defeneperpart)	
print('Defect Particles: ',len(defpartenergies),', Calculated line defect length: ',linedeflength,', Particles per unit length of defect: ', len(defpartenergies)/linedeflength)

#EDGES NEXT
#find coords of particles at top and bottom of edges to find distance between them (approx of length)
#first separating edge particles into near edge and far edge (won't work for high angles, when in/outside not separable in x)
nearedgeind = []
faredgeind = []
nearedgey = []
nearedgex = []
faredgey = []
faredgex = []
for x in range(0, len(edgelistx)):
	if edgelistx[x] < midx:
		nearedgeind.append(x)
		nearedgey.append(edgelisty[x])
		nearedgex.append(edgelistx[x])
	else:
		faredgeind.append(x)
		faredgey.append(edgelisty[x])
		faredgex.append(edgelistx[x])

#NEAR EDGE LENGTH
#topnearedgeind = nearedgeind[nearedgey.index(max(nearedgey))]
#bottomnearedgeind = nearedgeind[nearedgey.index(min(nearedgey))]

nearedgez = []
for i in range(0, len(nearedgey)):
	nearedger = np.sqrt(nearedgey[i]**2 + nearedgex[i]**2)
	#nearedgetheta = np.tan(nearedgey[i]/nearedgex[i])
	nearedgez.append((l*np.cos(anglep/2)) + halfconeheight - ((h/L)*(l+(L/2)-nearedger)))
avgnearedgez = sum(nearedgez)/len(nearedgez)
nearedgelength = 2*np.pi*avgnearedgez*np.tan(anglep/2) #if input is thetam then need to work thetap out properly

#nearedgelength = np.sqrt((edgelistx[bottomnearedgeind] - edgelistx[topnearedgeind])**2 + (edgelisty[bottomnearedgeind] - edgelisty[topnearedgeind])**2)


#FAR EDGE LENGTH
#topfaredgeind = faredgeind[faredgey.index(max(faredgey))]
#bottomfaredgeind = faredgeind[faredgey.index(min(faredgey))]

faredgez = []
for i in range(0, len(faredgey)):
	faredger = np.sqrt(faredgey[i]**2 + faredgex[i]**2)
	#faredgetheta = np.tan(faredgey[i]/faredgex[i])
	faredgez.append((l*np.cos(anglep/2)) + halfconeheight - ((h/L)*(l+(L/2)-faredger)))
avgfaredgez = sum(faredgez)/len(faredgez)
faredgelength = 2*np.pi*avgfaredgez*np.tan(anglep/2) #if input is thetam then need to work thetap out properly

#faredgelength = np.sqrt((edgelistx[bottomfaredgeind] - edgelistx[topfaredgeind])**2 + (edgelisty[bottomfaredgeind] - edgelisty[topfaredgeind])**2)

modellinedeflength = (avgfaredgez/(np.cos(anglep/2))) - (avgnearedgez/(np.cos(anglep/2)))
print('Model line defect length: ', modellinedeflength, 'Parts per unit length in defect using model: ', len(defpartenergies)/modellinedeflength)

edgeeneperpart = sum(edgepartenergies)/len(edgepartenergies)
edgeenergy = sum(edgepartenergies) / (nearedgelength+faredgelength)
print('----EDGES-----')
print('Edge Energy per Length: '+str(edgeenergy)+', per Particle: '+str(edgeeneperpart),' gain per Particle: ',6+edgeeneperpart)
print('Edge Particles: ',len(edgepartenergies),', Near edge length: ',nearedgelength,', Far edge length: ',faredgelength,', Total edge length: ',nearedgelength+faredgelength)
print('Particles per unit Length of Edges: ',len(edgepartenergies)/(nearedgelength+faredgelength))

print('Top of band is ', avgnearedgez, ' from the top of the cone.')

#adjedgeparts = int(input('Adjusted number of edge particles: '))
#print('Particles per unit Length of Edges: ',adjedgeparts/(nearedgelength+faredgelength))

