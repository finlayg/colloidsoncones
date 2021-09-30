#script to identify type of defect from given lowest input files and plot phase diagram of defect type on angle vs range plot
#Requires file tree to be: rho > angle > lowest and data

from scipy.spatial import Voronoi, voronoi_plot_2d
from matplotlib import pyplot as plt
import numpy as np
from collections import Counter
from pathlib import Path
from defectidentifierfunction import defectcategorisation

inputfilename = 'lowest'
datafilename = 'data'

rhostouse = input('What values of rho to use? (Should be same as folder names, input in ascending order separated by commas and with no spaces) ')
splitline = rhostouse.rstrip().split(',')
newrhostouse = [float(item) for item in splitline if item!=',']

inputminimarange = 1
inputminimum = 1

whichangle = 'p'
#whichangle = input('Input angle theta_p or theta_max? p/m ')

loangle = float(input('Lower limit of angle range (inclusive): '))
hiangle = float(input('Upper limit of angle range (inclusive): '))
interval = 0.1 
anglelistz = np.arange(loangle, hiangle+interval, interval)
print(anglelistz)

categoriseall = 'Y'

numparts = int(200)
#numparts = int(input('Number of particles (remember to add polygon if relevant): '))

phasediaglist = []

#First loop over values of rho
for rhoval in range(0, len(newrhostouse)):
	#print(rhoval)
	rholist = []

	#Second loop over values of angle
	for angleval in range(0, len(anglelistz)):
		#print(angleval)
		inputfile = str(int(newrhostouse[rhoval]))+'/'+str(round(anglelistz[angleval], 1))+'/'+inputfilename
		datafile = str(int(newrhostouse[rhoval]))+'/'+str(round(anglelistz[angleval], 1))+'/'+datafilename
		print('---------------------')
		
		if not Path(inputfile).is_file():
			print(inputfile+' does not exist')
			print('---------------------')
			continue 

		print('Categorising: '+inputfile)
	
		if whichangle is 'p':
			anglep = angleval #float(input('Angle theta_p of cone? '))#*(np.pi/180) 	#angle p of cone
			anglepr = anglep*(np.pi/180)					#converted to radians
			angle = (2*np.pi*np.sin(anglepr/2))#/(np.pi/180)
		elif whichangle is 'm':
			angle = angleval #float(input('Angle theta_max of cone? '))*(np.pi/180) #angle of net, converted to radians
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

		rotate = 'N'
		#rotate = input('Rotate net by a given angle? Y/N ')
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


		#loop to plot voronoi diagram of each minima in list of minima generated above
		for minimum in range(0, len(minimalist)):
			vordiag = None
			vordiag = Voronoi(minimalist[minimum])

			fig, ax = plt.subplots()
	
			#plt.gca().clf()
			voronoi_plot_2d(vordiag, ax=ax, show_vertices=False, point_size=5) #ax=ax,  show_vertices=False, point_size = 5)
			ax.axis('off')
	
			#-----------------------------------------------------------
		
			#Adding net of cone to plot

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
			
			#-----------------------------------------------------------
			
			#Adding colour to Voronoi cells based on shape
	
			nongreenindexlist = []
			bluelist = []
			redlist = []	#list of indices of particles in red and blue polygons for use in analysis
			
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
										nongreenindexlist.append(region)
									else:
										ax.fill(*zip(*polygon), "gold")
										nongreenindexlist.append(region)
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
								nongreenindexlist.append(region)
								redlist.append(region)
							else:
								ax.fill(*zip(*polygon), "r")
								nongreenindexlist.append(region)
								redlist.append(region)
						elif (regionlengths != []) and (len(vordiag.regions[region]) == 7) and (notlargepolygon):
							polygon = [vordiag.vertices[i] for i in vordiag.regions[region]]
							if q<1:
								ax.fill(*zip(*polygon), "b", label="Heptagonal environment")
								q=1
								nongreenindexlist.append(region)
								bluelist.append(region)
							else:
								ax.fill(*zip(*polygon), "b")
								nongreenindexlist.append(region)
								bluelist.append(region)
						elif (regionlengths != []) and (len(vordiag.regions[region]) in [3, 4, 8, 9]) and (notlargepolygon):
							polygon = [vordiag.vertices[i] for i in vordiag.regions[region]]
							if g<1:
								ax.fill(*zip(*polygon), "orange", label="3/4/5/9 sided environment")
								g=1
								nongreenindexlist.append(region)
							else:
								ax.fill(*zip(*polygon), "orange")
								nongreenindexlist.append(region)
		
			#-----------------------------------------------------------

			#Identifying largest defect on net and creating list of indices of its polygons
		
			#generating auxlist
			auxlist = []
			for i in range(0, len(nongreenindexlist)):		#making list of labels from 0 to list length - will be changed to defect labels
				auxlist.append(i)


			#print('Auxlistlength: '+str(len(auxlist)))

			for part1 in range(0, len(nongreenindexlist)):
				for part2 in range(0, len(nongreenindexlist)):					#for every subsequent non green particle (not yet looked at)
					adj = 0										#start adjacency count
					if (part2 is not part1):	#make sure not comparing particle with itself and it's on original net
						for i in vordiag.regions[nongreenindexlist[part1]]:			#if the two non green particles share a vertex
							if i in vordiag.regions[nongreenindexlist[part2]]:		#then
								adj = adj+1						#add one to adjacency count
							if adj is 2:							#if share two vertices
								if np.where(vordiag.point_region == nongreenindexlist[part2])[0][0] > 399:	#if neighbour is on lower copy
									
									nbourregindex = (np.where(vordiag.point_region == nongreenindexlist[part2])[0][0])			#index of part in neighbouring region
									ognetnbourregindex = (np.where(vordiag.point_region == (nongreenindexlist[part2]))[0][0]) - 400 	#same but on og net
									if vordiag.point_region[ognetnbourregindex] in nongreenindexlist:
										ognetnbournongreenindex = np.where(nongreenindexlist == vordiag.point_region[ognetnbourregindex])[0][0] #its index in nongreenindexlist
									
										if ognetnbournongreenindex > part1:
											for part in range(0, len(auxlist)):
												if auxlist[part] is auxlist[ognetnbournongreenindex]:
													auxlist[part] = auxlist[part1]
										if part1 > ognetnbournongreenindex:
											for part in range(0, len(auxlist)):
												if auxlist[part] is auxlist[part1]:
													auxlist[part] = auxlist[ognetnbournongreenindex]
			
								elif (np.where(vordiag.point_region == nongreenindexlist[part2])[0][0] > 199) and (np.where(vordiag.point_region == nongreenindexlist[part2])[0][0] < 400):													#if neighbour is on upper copy
									nbourregindex = (np.where(vordiag.point_region == nongreenindexlist[part2])[0][0]) 			#index of particle in neighbouring region
									ognetnbourregindex = (np.where(vordiag.point_region == (nongreenindexlist[part2]))[0][0]) - 200		#same but on og net
									if vordiag.point_region[ognetnbourregindex] in nongreenindexlist:		
										ognetnbournongreenindex = np.where(nongreenindexlist == vordiag.point_region[ognetnbourregindex])[0][0]	#its index in nongreenindexlist
			
										if ognetnbournongreenindex > part1:
											for part in range(0, len(auxlist)):
												if auxlist[part] is auxlist[ognetnbournongreenindex]:
													auxlist[part] = auxlist[part1]
										if part1 > ognetnbournongreenindex:
											for part in range(0, len(auxlist)):
												if auxlist[part] is auxlist[part1]:
													auxlist[part] = auxlist[ognetnbournongreenindex]
			
								elif np.where(vordiag.point_region == nongreenindexlist[part2])[0][0] < 200:	#if neighbour is on original
			
									if part2 > part1:
										for part in range(0, len(auxlist)):			#find each particle in the same defect as defect 2
											if auxlist[part] is auxlist[part2]:		
												auxlist[part] = auxlist[part1]		#and make them all part of defect 1
									if part1 > part2:
										for part in range(0, len(auxlist)):			#find each particle in the same defect as defect 2
											if auxlist[part] is auxlist[part1]:	
												auxlist[part] = auxlist[part2]		#and make them all part of defect 1
			
			#finding modal label in auxlist
			data = Counter(auxlist)
			biggestdefect = data.most_common(1)
			
			biggestdefectlist = []
			for pg in range(0, len(auxlist)):
				if (auxlist[pg] is biggestdefect[0][0]) and (np.where(vordiag.point_region == nongreenindexlist[pg])[0][0] < 200):
					#add index of this particle to a list
					biggestdefectlist.append(nongreenindexlist[pg])
			
			#converting red and blue polys in original net to list of associated particles
			reddeflist = []
			bluedeflist = []
			for pg in range(0, len(redlist)):
				if (redlist[pg] in biggestdefectlist) and (np.where(vordiag.point_region == redlist[pg])[0][0] < 200):
					#add index of this particle to a list
					reddeflist.append(redlist[pg])
			for pg in range(0, len(bluelist)):
				if (bluelist[pg] in biggestdefectlist) and (np.where(vordiag.point_region == bluelist[pg])[0][0] < 200):
					#add index of this particle to a list
					bluedeflist.append(bluelist[pg])
			
			
			#print(biggestdefect[0][0], biggestdefectlist)
				
			for polygonss in biggestdefectlist:
				polygona = [vordiag.vertices[i] for i in vordiag.regions[polygonss]]
				ax.fill(*zip(*polygona), "purple", label="3/4/5/9 sided environment")
			
			anglename = ''
			for char in str(anglep):
				if char is not '.':
					anglename = anglename+char
			
			#----------------------------------

			#Now Analysing Identified Defect - First Finding Usable Coords of Biggest Defect
	
			#biggestdefectlist contains list of polygon indices. Want list of particle indices
			biggestdefectpartlist = []
			xlist = []
			ylist = []
			biggestdefectpartxylist = []
			
			for polyg in range(0, len(biggestdefectlist)):
				part = np.where(vordiag.point_region == biggestdefectlist[polyg])[0][0]
				biggestdefectpartlist.append(part)
				xlist.append(minimalist[minimum][part][0])
				ylist.append(minimalist[minimum][part][1])
				biggestdefectpartxylist.append([minimalist[minimum][part][0],minimalist[minimum][part][1]])
			
			redpartxylist = []
			bluepartxylist = []
			redxlist = []
			redylist = []
			bluexlist = []
			blueylist = []
			redpartlist = []
			bluepartlist = [] 
			
			for polyg in range(0, len(reddeflist)):
				part = np.where(vordiag.point_region == reddeflist[polyg])[0][0]
				redpartlist.append(part)
				redxlist.append(minimalist[minimum][part][0])
				redylist.append(minimalist[minimum][part][1])
				redpartxylist.append([minimalist[minimum][part][0],minimalist[minimum][part][1]])
			for polyg in range(0, len(bluedeflist)):
				part = np.where(vordiag.point_region == bluedeflist[polyg])[0][0]
				bluepartlist.append(part)
				bluexlist.append(minimalist[minimum][part][0])
				blueylist.append(minimalist[minimum][part][1])
				bluepartxylist.append([minimalist[minimum][part][0],minimalist[minimum][part][1]])
			
			#plt.scatter(redxlist, redylist, marker ='*', zorder = 100)
			#plt.scatter(bluexlist, blueylist, marker ='p', zorder = 100)
			
			#now generating full defect from original net by repeating it by the number of distinct, non-connected defect sections there are in the individual net
			#First find number of distinct defect sections on original net
	
			#generating auxlist
			auxlist = []
			for i in range(0, len(biggestdefectlist)):		#making list of labels from 0 to list length - will be changed to defect labels
				auxlist.append(i)
			
			#finding separated defect sections and assigning each polygon to one of them 
			for part1 in range(0, len(biggestdefectlist)):
				for part2 in range(0, len(biggestdefectlist)):							#for every subsequent non green particle (not yet looked at)
					adj = 0											#start adjacency count
					if (part2 is not part1):								#make sure not comparing particle with itself and it's on original net
						for i in vordiag.regions[biggestdefectlist[part1]]:				#if the two non green particles share a vertex
							if i in vordiag.regions[biggestdefectlist[part2]]:			#then
								adj = adj+1							#add one to adjacency count
							if adj is 2:								#if share two vertices
								if part2 > part1:
									for part in range(0, len(auxlist)):			#find each particle in the same defect as defect 2
										if auxlist[part] is auxlist[part2]:		
											auxlist[part] = auxlist[part1]		#and make them all part of defect 1
								if part1 > part2:
									for part in range(0, len(auxlist)):			#find each particle in the same defect as defect 2
										if auxlist[part] is auxlist[part1]:	
											auxlist[part] = auxlist[part2]		#and make them all part of defect 1
			newlist = []
			for i in auxlist:
				if i not in newlist:
					newlist.append(i)
			#print(len(newlist), newlist)
			
			numberdefectsections = len(newlist)
			#print(numberdefectsections)
			
			#now repeat defect that many times and find the biggest defect
			if numberdefectsections > 1:
					
				#convert coords to polars - defect first then red and blue particles
				polarbiggestdefectpartlist = []
				for particle in biggestdefectpartxylist:
					r = np.sqrt(float(particle[0])**2+float(particle[1])**2)
					if ((float(particle[0]) >= 0) and (float(particle[1]) >= 0)): #top right quadrant
						theta = np.arctan(float(particle[1])/float(particle[0]))
					elif ((float(particle[0]) < 0) and (float(particle[1]) >= 0)): #top left quadrant
						theta = (np.arctan(float(particle[1])/float(particle[0])) + (np.pi))
					elif ((float(particle[0]) >= 0) and (float(particle[1]) < 0)): #bottom left quadrant
						theta = np.arctan(float(particle[1])/float(particle[0])) + (np.pi)
					elif ((float(particle[0]) < 0) and (float(particle[1]) < 0)): #bottom right quadrant
						theta = np.arctan(float(particle[1])/float(particle[0])) 
					polarbiggestdefectpartlist.append([r, theta])
			
				polarredpartlist = []
				for particle in redpartxylist:
					r = np.sqrt(float(particle[0])**2+float(particle[1])**2)
					if ((float(particle[0]) >= 0) and (float(particle[1]) >= 0)): #top right quadrant
						theta = np.arctan(float(particle[1])/float(particle[0]))
					elif ((float(particle[0]) < 0) and (float(particle[1]) >= 0)): #top left quadrant
						theta = (np.arctan(float(particle[1])/float(particle[0])) + (np.pi))
					elif ((float(particle[0]) >= 0) and (float(particle[1]) < 0)): #bottom left quadrant
						theta = np.arctan(float(particle[1])/float(particle[0])) + (np.pi)
					elif ((float(particle[0]) < 0) and (float(particle[1]) < 0)): #bottom right quadrant
						theta = np.arctan(float(particle[1])/float(particle[0])) 
					polarredpartlist.append([r, theta])
				polarbluepartlist = []
				for particle in bluepartxylist:
					r = np.sqrt(float(particle[0])**2+float(particle[1])**2)
					if ((float(particle[0]) >= 0) and (float(particle[1]) >= 0)): #top right quadrant
						theta = np.arctan(float(particle[1])/float(particle[0]))
					elif ((float(particle[0]) < 0) and (float(particle[1]) >= 0)): #top left quadrant
						theta = (np.arctan(float(particle[1])/float(particle[0])) + (np.pi))
					elif ((float(particle[0]) >= 0) and (float(particle[1]) < 0)): #bottom left quadrant
						theta = np.arctan(float(particle[1])/float(particle[0])) + (np.pi)
					elif ((float(particle[0]) < 0) and (float(particle[1]) < 0)): #bottom right quadrant
						theta = np.arctan(float(particle[1])/float(particle[0])) 
					polarbluepartlist.append([r, theta])
			
				#copy portion of net in both directions - not robust if total angle exceeds 360
			
				extpolarbiggestdefectpartlist = polarbiggestdefectpartlist.copy()
				extpolarredpartlist = polarredpartlist.copy()
				extpolarbluepartlist = polarbluepartlist.copy()
			
				for number in range(2, numberdefectsections+1):
					for particle in polarbiggestdefectpartlist:						#adding top copy of net
						extpolarbiggestdefectpartlist.append([particle[0], (particle[1]+((number-1)*angle))])
					for particle in polarredpartlist:						#adding top copy of net
						extpolarredpartlist.append([particle[0], (particle[1]+((number-1)*angle))])
					for particle in polarbluepartlist:						#adding top copy of net
						extpolarbluepartlist.append([particle[0], (particle[1]+((number-1)*angle))])
			
				#convert back to cartesians
				extbiggestdefectpartlist = []
				extbiggestdefectpartxlist = []
				extbiggestdefectpartylist = []
			
				extredpartlist = []
				extredpartxlist = []
				extredpartylist = []
				extbluepartlist = []
				extbluepartxlist = []
				extbluepartylist = []
			
				for particle in extpolarbiggestdefectpartlist:
					x = particle[0]*np.cos(particle[1])
					y = particle[0]*np.sin(particle[1])
					extbiggestdefectpartlist.append([x, y])
					extbiggestdefectpartxlist.append(x)
					extbiggestdefectpartylist.append(y)
		
				for particle in extpolarredpartlist:
					x = particle[0]*np.cos(particle[1])
					y = particle[0]*np.sin(particle[1])
					extredpartlist.append([x, y])
					extredpartxlist.append(x)
					extredpartylist.append(y)
				for particle in extpolarbluepartlist:
					x = particle[0]*np.cos(particle[1])
					y = particle[0]*np.sin(particle[1])
					extbluepartlist.append([x, y])
					extbluepartxlist.append(x)
					extbluepartylist.append(y)
		
			elif numberdefectsections is 1:
				extbiggestdefectpartlist = biggestdefectpartxylist.copy()
				extredpartlist = redpartxylist.copy()
				extbluepartlist = bluepartxylist.copy()
		
			#plt.scatter(extbiggestdefectpartxlist, extbiggestdefectpartylist, marker ='*', zorder = 100)
		
			#Now find biggest group of neighbours by using proximity - if closer than 1.5 then neighbours
		
			#generating auxlist
			auxlist = []
			for i in range(0, len(extbiggestdefectpartlist)):		#making list of labels from 0 to list length - will be changed to defect labels
				auxlist.append(i)
		
			#finding separated defect sections and assigning each polygon to one of them 
			for part1 in range(0, len(extbiggestdefectpartlist)):
				for part2 in range(0, len(extbiggestdefectpartlist)):							#for every subsequent non green particle (not yet looked at)
					adj = 0											#start adjacency count
					if (part2 is not part1):								#make sure not comparing particle with itself and it's on original net
						dist = np.sqrt((extbiggestdefectpartlist[part1][0]-extbiggestdefectpartlist[part2][0])**2+(extbiggestdefectpartlist[part1][1]-extbiggestdefectpartlist[part2][1])**2)#distance between particles 1 and 2
						if dist < 1.5:
							if part2 > part1:
								for part in range(0, len(auxlist)):			#find each particle in the same defect as defect 2
									if auxlist[part] is auxlist[part2]:		
										auxlist[part] = auxlist[part1]		#and make them all part of defect 1
							if part1 > part2:
								for part in range(0, len(auxlist)):			#find each particle in the same defect as defect 2
									if auxlist[part] is auxlist[part1]:	
										auxlist[part] = auxlist[part2]		#and make them all part of defect 1
			#finding modal label in auxlist (corresponds to biggest defect)
			if len(auxlist) > 0:
				data2 = Counter(auxlist)
				biggestdefectlabel = data2.most_common(1)[0][0]
				
				defectcoordsforanalysis = []
				defectxcoordsforanalysis = []
				defectycoordsforanalysis = []
				for i in range(0, len(auxlist)):
					if auxlist[i] is biggestdefectlabel:
						defectcoordsforanalysis.append(extbiggestdefectpartlist[i])
						defectxcoordsforanalysis.append(extbiggestdefectpartlist[i][0])
						defectycoordsforanalysis.append(extbiggestdefectpartlist[i][1])
				
				redcoordsforanalysis = []
				redxcoordsforanalysis = []
				redycoordsforanalysis = []
				bluecoordsforanalysis = []
				bluexcoordsforanalysis = []
				blueycoordsforanalysis = []
				
				for i in range(0, len(extredpartlist)):
					for j in range(0, len(defectcoordsforanalysis)):
						if ((extredpartlist[i][0] - defectcoordsforanalysis[j][0])**2 < 0.01) and ((extredpartlist[i][1] - defectcoordsforanalysis[j][1])**2 < 0.01): 
							redcoordsforanalysis.append(extredpartlist[i])
							redxcoordsforanalysis.append(extredpartlist[i][0])
							redycoordsforanalysis.append(extredpartlist[i][1])
				for i in range(0, len(extbluepartlist)):
					for j in range(0, len(defectcoordsforanalysis)):
						if ((extbluepartlist[i][0] - defectcoordsforanalysis[j][0])**2 < 0.01) and ((extbluepartlist[i][1] - defectcoordsforanalysis[j][1])**2 < 0.01): 
							bluecoordsforanalysis.append(extbluepartlist[i])
							bluexcoordsforanalysis.append(extbluepartlist[i][0])
							blueycoordsforanalysis.append(extbluepartlist[i][1])
				
				#plt.scatter(defectxcoordsforanalysis, defectycoordsforanalysis, marker ='*', zorder = 100)
				
				#generating and saving plot
				if categoriseall is 'N':
					fig.savefig('defectid'+str(anglename)+'_min'+str(inputminimum)+'.png')
					print('Defect ID - Minimum '+str(inputminimum)+' saved: '+'defectid'+str(anglename)+'_min'+str(inputminimum)+'.png')
#				elif categoriseall is 'Y':
#					fig.savefig('defectid'+str(anglename)+'_min'+str(minimum+1)+'.png')
#					print('Defect ID - Minimum '+str(minimum+1)+' saved: '+'defectid'+str(anglename)+'_min'+str(minimum+1)+'.png')
				fig.clf()
				plt.close()
	
				#---------------------------------
	
				#Plot Biggest Defect
	
				fig, ax = plt.subplots()
				plt.scatter(defectxcoordsforanalysis, defectycoordsforanalysis, marker ='o', zorder = 50, color = 'gold')
				plt.ylim(min(defectycoordsforanalysis)-1, max(defectycoordsforanalysis)+1)
				if ((max(defectxcoordsforanalysis)+1) - (min(defectxcoordsforanalysis)-1)) > 9:
					plt.xlim(min(defectxcoordsforanalysis)-1, max(defectxcoordsforanalysis)+1)
				else:
					centre = (max(defectxcoordsforanalysis) - min(defectxcoordsforanalysis))/2 + min(defectxcoordsforanalysis)
					plt.xlim(centre-5, centre+5)
			
				#ANALYSIS ON COORDINATES OF BIGGEST DEFECT
	
				#1 - Find CoM and work out if within defect or not, also calculate proportional inertia - m*r**2 scaled to number of parts in defect
				#Finding CoM
				CoMx = np.mean(defectxcoordsforanalysis)
				CoMy = np.mean(defectycoordsforanalysis)
				plt.scatter(CoMx, CoMy, marker = '*', s = 200, zorder = 200, color = 'k')
			
				inertialist = []
				for part in range(0, len(defectcoordsforanalysis)):
					r = np.sqrt((defectcoordsforanalysis[part][0]-CoMx)**2+(defectcoordsforanalysis[part][1]-CoMy)**2)
					inertialist.append(r**2)
				norminert = np.sum(inertialist)/len(defectcoordsforanalysis)
				if categoriseall is 'N':
					print('Total inertia: '+str(np.sum(inertialist)))
					print('Total normalised inertia: '+str(norminert))

				#2 - Vector Histogram? Or Use whether vectors are consistent to define line slip - will be constant + also uniform rotation
				anglelist = []
				for part1 in range(0, len(defectcoordsforanalysis)):
					for part2 in range(0, len(defectcoordsforanalysis)):
						if part1 is not part2 :
							dist = np.sqrt((defectcoordsforanalysis[part1][0]-defectcoordsforanalysis[part2][0])**2+(defectcoordsforanalysis[part1][1]-defectcoordsforanalysis[part2][1])**2)
							if dist < 1.5:
								vecx = defectcoordsforanalysis[part2][0]-defectcoordsforanalysis[part1][0]
								vecy = defectcoordsforanalysis[part2][1]-defectcoordsforanalysis[part1][1]
								vector = np.array([vecx, vecy])
								origin = np.array([defectcoordsforanalysis[part1][0],defectcoordsforanalysis[part1][1]])
								plt.quiver(*origin, vector[0], vector[1], color = 'orange', zorder = 5, scale = 18)
				
								#First converting vector to correct angle from 0-360
								if (vector[0] >= 0) and (vector[1] >= 0): #top right quadrant
									value = (np.arctan(vector[1]/vector[0])*(180/np.pi))
								elif (vector[0] < 0) and (vector[1] >= 0): #top left quadrant
									value = (np.arctan(vector[1]/vector[0]) + (np.pi)) * (180/np.pi)
								elif (vector[0] < 0) and (vector[1] < 0): #bottom left quadrant
									value = (np.arctan(vector[1]/vector[0]) + (np.pi)) * (180/np.pi)
								elif (vector[0] >= 0) and (vector[1] < 0): #bottom right quadrant
									value = np.arctan(vector[1]/vector[0])*(180/np.pi)
								#rotating angles so all from 0-360							
								if value < 0:
									value = value + 360
								if (value <= 360) and (value > 180):
									value = value - 180
							
								anglelist.append(value)
			
				#3 - Find longest chain of particles in defect
			
				#pick one particle (furthest from com maybe?)
				#find vectors to all of its neighbours 
				#continue the vectors to the next neighbour (within certain angle tolerance)
				#if there's another neighbour then make the vector the vector from the first part to this particle instead
				#repeat until no more neighbours lie on same line as vector chain
				#store overall vector, list of parts in chain, vectors between each pair of particles
				#particle1: origin of chain
				#particle2: penultimate particle in current chain
				#particle3: current end of chain
				#particle4: new particle to be added to chain - neighbour of particle2
		
		
				longestvectorlist = []
				vectorsinchainlist = []
				vectorchainlist = []
				vecchainanglelist = []
				vectorchainxlist = []
				vectorchainylist = []
			
				finalchainlist = []
				anglethreshold = 4 	#upper threshold for difference squared between two angles for angles to be similar enough
				
				for part1 in range(0, len(defectcoordsforanalysis)):
					particle1index = part1
					particle2index = 0
					particle3index = 0
					particle4index = 0
					for part2 in range(0, len(defectcoordsforanalysis)):
						if part2 is not particle1index :
							dist = np.sqrt((defectcoordsforanalysis[particle1index][0]-defectcoordsforanalysis[part2][0])**2+(defectcoordsforanalysis[particle1index][1]-defectcoordsforanalysis[part2][1])**2)
							if dist < 1.5: #defining vectors + angles for all pairs of confirmed neighbours of particle 1
								currentchainanglelist = []
								currentchainpartlist = []
								particle2index = part2
								vecx1 = defectcoordsforanalysis[part2][0]-defectcoordsforanalysis[particle1index][0]
								vecy1 = defectcoordsforanalysis[part2][1]-defectcoordsforanalysis[particle1index][1]
								vector1 = np.array([vecx1, vecy1])
								origin = np.array([defectcoordsforanalysis[particle1index][0],defectcoordsforanalysis[particle1index][1]])
								if (vector1[0] >= 0) and (vector1[1] >= 0): #top right quadrant
									angle1 = (np.arctan(vector1[1]/vector1[0])*(180/np.pi))
								elif (vector1[0] < 0) and (vector1[1] >= 0): #top left quadrant
									angle1 = (np.arctan(vector1[1]/vector1[0]) + (np.pi)) * (180/np.pi)
								elif (vector1[0] < 0) and (vector1[1] < 0): #bottom left quadrant
									angle1 = (np.arctan(vector1[1]/vector1[0]) + (np.pi)) * (180/np.pi)
								elif (vector1[0] >= 0) and (vector1[1] < 0): #bottom right quadrant
									angle1 = np.arctan(vector1[1]/vector1[0])*(180/np.pi)
				
								for part3 in range(0, len(defectcoordsforanalysis)): #looking for next particle and see if add to vector chain
									if part3 is not particle2index and part3 is not particle1index:	#particle3 is not 1 or 2
										dist2 = np.sqrt((defectcoordsforanalysis[part3][0]-defectcoordsforanalysis[particle2index][0])**2+(defectcoordsforanalysis[part3][1]-defectcoordsforanalysis[particle2index][1])**2)
										if dist2 < 1.5: #for all confirmed neighbours of particle 2
											vecx2 = defectcoordsforanalysis[part3][0]-defectcoordsforanalysis[particle2index][0]
											vecy2 = defectcoordsforanalysis[part3][1]-defectcoordsforanalysis[particle2index][1]
											vector2 = np.array([vecx2, vecy2])
											if (vector2[0] >= 0) and (vector2[1] >= 0): #top right quadrant
												angle2 = (np.arctan(vector2[1]/vector2[0])*(180/np.pi))
											elif (vector2[0] < 0) and (vector2[1] >= 0): #top left quadrant
												angle2 = (np.arctan(vector2[1]/vector2[0]) + (np.pi)) * (180/np.pi)
											elif (vector2[0] < 0) and (vector2[1] < 0): #bottom left quadrant
												angle2 = (np.arctan(vector2[1]/vector2[0]) + (np.pi)) * (180/np.pi)
											elif (vector2[0] >= 0) and (vector2[1] < 0): #bottom right quadrant
												angle2 = np.arctan(vector2[1]/vector2[0])*(180/np.pi)
				
											if ((angle1-angle2)**2 < anglethreshold): 	#confirming p3-p2 angle is similar enough to p2-p1
												particle3index = part3
												currentchainanglelist.append(angle1)	#adding p2-p1 angle (between first two parts in chain)
												currentchainanglelist.append(angle2)	#adding p3-p2 angle (between second and third particle in chain)
												currentchainpartlist.append(particle1index)	#adding p1 p2 and p3 to list
												currentchainpartlist.append(particle2index)
												currentchainpartlist.append(particle3index)
												#print('P1: ', particle1index, defectcoordsforanalysis[particle1index][0], defectcoordsforanalysis[particle1index][1])
												#print('P2: ', particle2index, defectcoordsforanalysis[particle2index][0], defectcoordsforanalysis[particle2index][1])
												#print('V1: ', vector1, 'A1: ', angle1, 'V2: ', vector2, 'A2: ', angle2, currentchainpartlist, currentchainanglelist)
				
												particle4index = None
												startover = True
												while startover:
													startover = False
													for part4 in range(0, len(defectcoordsforanalysis)):
														particle3index = currentchainpartlist[-1]
														particle2index = currentchainpartlist[-2]
														if part4 not in currentchainpartlist:
															dist3 = np.sqrt((defectcoordsforanalysis[part4][0]-defectcoordsforanalysis[particle3index][0])**2+(defectcoordsforanalysis[part4][1]-defectcoordsforanalysis[particle3index][1])**2)
															if dist3 < 1.5:
																vecx3 = defectcoordsforanalysis[part4][0]-defectcoordsforanalysis[particle3index][0]
																vecy3 = defectcoordsforanalysis[part4][1]-defectcoordsforanalysis[particle3index][1]
																vector3 = np.array([vecx3, vecy3])
																if (vector3[0] >= 0) and (vector3[1] >= 0): #top right quadrant
																	angle3 = (np.arctan(vector3[1]/vector3[0])*(180/np.pi))
																elif (vector3[0] < 0) and (vector3[1] >= 0): #top left quadrant
																	angle3 = (np.arctan(vector3[1]/vector3[0]) + (np.pi)) * (180/np.pi)
																elif (vector3[0] < 0) and (vector3[1] < 0): #bottom left quadrant
																	angle3 = (np.arctan(vector3[1]/vector3[0]) + (np.pi)) * (180/np.pi)
																elif (vector3[0] >= 0) and (vector3[1] < 0): #bottom right quadrant
																	angle3 = np.arctan(vector3[1]/vector3[0])*(180/np.pi)
					
																if ((angle2-angle3)**2 < anglethreshold):#confirming p4-p3 angle is similar enough to p3-p2
																	currentchainanglelist.append(angle3)
																	currentchainpartlist.append(part4)
																	angle2 = angle3
																	particle4index = False			#have added a particle4
																	startover = True
																	break
				
												else:
													if particle4index is None: 	#chain is three particles long
														overallvecx = defectcoordsforanalysis[particle3index][0]-defectcoordsforanalysis[particle1index][0]
														overallvecy = defectcoordsforanalysis[particle3index][1]-defectcoordsforanalysis[particle1index][1]
														overallvector = np.array([overallvecx, overallvecy])
														if (overallvector[0] >= 0) and (overallvector[1] >= 0): #top right quadrant
															overallangle = (np.arctan(overallvector[1]/overallvector[0])*(180/np.pi))
														elif (overallvector[0] < 0) and (overallvector[1] >= 0): #top left quadrant
															overallangle = (np.arctan(overallvector[1]/overallvector[0]) + (np.pi)) * (180/np.pi)
														elif (overallvector[0] < 0) and (overallvector[1] < 0): #bottom left quadrant
															overallangle = (np.arctan(overallvector[1]/overallvector[0]) + (np.pi)) * (180/np.pi)
														elif (overallvector[0] >= 0) and (overallvector[1] < 0): #bottom right quadrant
															overallangle = np.arctan(overallvector[1]/overallvector[0])*(180/np.pi)
														finalchainlist.append([currentchainpartlist, currentchainanglelist, overallangle, overallvector])
														currentchainpartlist = []
					
													elif particle4index is False:	#chain more than three long but no more particles to add
														particle4index = currentchainpartlist[-1]
														overallvecx = defectcoordsforanalysis[particle4index][0]-defectcoordsforanalysis[particle1index][0]
														overallvecy = defectcoordsforanalysis[particle4index][1]-defectcoordsforanalysis[particle1index][1]
														overallvector = np.array([overallvecx, overallvecy])
														if (overallvector[0] >= 0) and (overallvector[1] >= 0): #top right quadrant
															overallangle = (np.arctan(overallvector[1]/overallvector[0])*(180/np.pi))
														elif (overallvector[0] < 0) and (overallvector[1] >= 0): #top left quadrant
															overallangle = (np.arctan(overallvector[1]/overallvector[0]) + (np.pi)) * (180/np.pi)
														elif (overallvector[0] < 0) and (overallvector[1] < 0): #bottom left quadrant
															overallangle = (np.arctan(overallvector[1]/overallvector[0]) + (np.pi)) * (180/np.pi)
														elif (overallvector[0] >= 0) and (overallvector[1] < 0): #bottom right quadrant
															overallangle = np.arctan(overallvector[1]/overallvector[0])*(180/np.pi)
				
														if overallangle < 0:
															overallangle = overallangle + 180
														elif (overallangle > 180) and (overallangle < 360):
															overallangle = overallangle - 180
				
														finalchainlist.append([currentchainpartlist, currentchainanglelist, overallangle, overallvector])
														currentchainpartlist = []
					
				
				#Removing duplicate chains sharing two particles, keeping longer one
				startagain = True
				while startagain:
					startagain = False
					#print('started while loop')
					for chain in range(0, len(finalchainlist)):
						for chain2 in range(chain, len(finalchainlist)):
							count = 0
							if chain is not chain2:
								#print(finalchainlist[chain][0])
								#print(finalchainlist[chain2][0])
								for chainpart in range(0, len(finalchainlist[chain][0])):
									#print(chainpart)
									if ((finalchainlist[chain][0][chainpart]) in (finalchainlist[chain2][0])):
										#print('adding to count')
										count += 1
								if count > 1:
									#print(len(finalchainlist))
									if len(finalchainlist[chain][0]) >= len(finalchainlist[chain2][0]):
										del finalchainlist[chain2]
									elif len(finalchainlist[chain][0]) < len(finalchainlist[chain2][0]):
										del finalchainlist[chain]
									#print('DEL, ', len(finalchainlist))
									startagain = True
									break
							if startagain:
								break			
						if startagain:
							break	
				
				if categoriseall is 'N':
					print('Final Chain List contains: ',len(finalchainlist),' chains')#, finalchainlist)
				
				#Plotting points and vectors in 4 longest Chains
				if len(finalchainlist) > 0:
					whichchain = 0
					chainxplot = []
					chainyplot = []
					for i in range(0, len(finalchainlist[whichchain][0])):
						chainxplot.append(defectxcoordsforanalysis[finalchainlist[whichchain][0][i]])
						chainyplot.append(defectycoordsforanalysis[finalchainlist[whichchain][0][i]])
					plt.scatter(chainxplot, chainyplot, marker = '1', s = 150, zorder = 300, color = 'm', label='Chain 0')
					origin = np.array([defectxcoordsforanalysis[finalchainlist[whichchain][0][0]],defectycoordsforanalysis[finalchainlist[whichchain][0][0]]])
					plt.quiver(*origin, 2*finalchainlist[whichchain][3][0], 2*finalchainlist[whichchain][3][1], color = 'g', zorder = 5, scale = 18, headaxislength=0,headwidth=1,headlength=0)
				if len(finalchainlist) > 1:
					whichchain = 1
					chainxplot = []
					chainyplot = []
					for i in range(0, len(finalchainlist[whichchain][0])):
						chainxplot.append(defectxcoordsforanalysis[finalchainlist[whichchain][0][i]])
						chainyplot.append(defectycoordsforanalysis[finalchainlist[whichchain][0][i]])
					plt.scatter(chainxplot, chainyplot, marker = '2', s = 150, zorder = 300, color = 'c', label='Chain 1')
					origin = np.array([defectxcoordsforanalysis[finalchainlist[whichchain][0][0]],defectycoordsforanalysis[finalchainlist[whichchain][0][0]]])
					plt.quiver(*origin, 2*finalchainlist[whichchain][3][0], 2*finalchainlist[whichchain][3][1], color = 'g', zorder = 5, scale = 18, headaxislength=0,headwidth=1,headlength=0)
				if len(finalchainlist) > 2:
					whichchain = 2
					chainxplot = []
					chainyplot = []
					for i in range(0, len(finalchainlist[whichchain][0])):
						chainxplot.append(defectxcoordsforanalysis[finalchainlist[whichchain][0][i]])
						chainyplot.append(defectycoordsforanalysis[finalchainlist[whichchain][0][i]])
					plt.scatter(chainxplot, chainyplot, marker = '3', s = 150, zorder = 300, color = 'tan', label='Chain 2')
					origin = np.array([defectxcoordsforanalysis[finalchainlist[whichchain][0][0]],defectycoordsforanalysis[finalchainlist[whichchain][0][0]]])
					plt.quiver(*origin, 2*finalchainlist[whichchain][3][0], 2*finalchainlist[whichchain][3][1], color = 'g', zorder = 5, scale = 18, headaxislength=0,headwidth=1,headlength=0)
				if len(finalchainlist) > 3:
					whichchain = 3
					chainxplot = []
					chainyplot = []
					for i in range(0, len(finalchainlist[whichchain][0])):
						chainxplot.append(defectxcoordsforanalysis[finalchainlist[whichchain][0][i]])
						chainyplot.append(defectycoordsforanalysis[finalchainlist[whichchain][0][i]])
					plt.scatter(chainxplot, chainyplot, marker = '4', s = 150, zorder = 300, color = 'lightcoral', label='Chain 3')
					origin = np.array([defectxcoordsforanalysis[finalchainlist[whichchain][0][0]],defectycoordsforanalysis[finalchainlist[whichchain][0][0]]])
					plt.quiver(*origin, 2*finalchainlist[whichchain][3][0], 2*finalchainlist[whichchain][3][1], color = 'g', zorder = 5, scale = 18, headaxislength=0,headwidth=1,headlength=0)
				if len(finalchainlist) > 4:
					whichchain = 4
					chainxplot = []
					chainyplot = []
					for i in range(0, len(finalchainlist[whichchain][0])):
						chainxplot.append(defectxcoordsforanalysis[finalchainlist[whichchain][0][i]])
						chainyplot.append(defectycoordsforanalysis[finalchainlist[whichchain][0][i]])
					plt.scatter(chainxplot, chainyplot, marker = 'x', s = 150, zorder = 300, color = 'lightgreen', label='Chain 4')
					origin = np.array([defectxcoordsforanalysis[finalchainlist[whichchain][0][0]],defectycoordsforanalysis[finalchainlist[whichchain][0][0]]])
					plt.quiver(*origin, 2*finalchainlist[whichchain][3][0], 2*finalchainlist[whichchain][3][1], color = 'g', zorder = 5, scale = 18, headaxislength=0,headwidth=1,headlength=0)
				if len(finalchainlist) > 5:
					whichchain = 5
					chainxplot = []
					chainyplot = []
					for i in range(0, len(finalchainlist[whichchain][0])):
						chainxplot.append(defectxcoordsforanalysis[finalchainlist[whichchain][0][i]])
						chainyplot.append(defectycoordsforanalysis[finalchainlist[whichchain][0][i]])
					plt.scatter(chainxplot, chainyplot, marker = '+', s = 150, zorder = 300, color = 'saddlebrown', label='Chain 5')
					origin = np.array([defectxcoordsforanalysis[finalchainlist[whichchain][0][0]],defectycoordsforanalysis[finalchainlist[whichchain][0][0]]])
					plt.quiver(*origin, 2*finalchainlist[whichchain][3][0], 2*finalchainlist[whichchain][3][1], color = 'g', zorder = 5, scale = 18, headaxislength=0,headwidth=1,headlength=0)
				if len(finalchainlist) > 6:
					whichchain = 6
					chainxplot = []
					chainyplot = []
					for i in range(0, len(finalchainlist[whichchain][0])):
						chainxplot.append(defectxcoordsforanalysis[finalchainlist[whichchain][0][i]])
						chainyplot.append(defectycoordsforanalysis[finalchainlist[whichchain][0][i]])
					plt.scatter(chainxplot, chainyplot, marker = '1', s = 150, zorder = 300, color = 'gray', label='Chain 6')
					origin = np.array([defectxcoordsforanalysis[finalchainlist[whichchain][0][0]],defectycoordsforanalysis[finalchainlist[whichchain][0][0]]])
					plt.quiver(*origin, 2*finalchainlist[whichchain][3][0], 2*finalchainlist[whichchain][3][1], color = 'g', zorder = 5, scale = 18, headaxislength=0,headwidth=1,headlength=0)
				
				
				#4 - Red/Blue pair or not - difficult as voronoi dependent
		
				plt.scatter(redxcoordsforanalysis, redycoordsforanalysis, marker = 'h', s = 150, zorder = 200, color = 'r')
				plt.scatter(bluexcoordsforanalysis, blueycoordsforanalysis, marker = 'h', s = 150, zorder = 200, color = 'b')
				
				#generating and saving defect visualisation plot
				#plt.legend()
				ax.set_aspect('equal')
				if categoriseall is 'N':
					fig.savefig('defectidjustdef'+str(anglename)+'_min'+str(inputminimum)+'.png', dpi = 300)
					print('Just Defect - Minimum '+str(inputminimum)+' saved: '+'defectidjustdef'+str(anglename)+'_min'+str(inputminimum)+'.png')
#				elif categoriseall is 'Y':
#					fig.savefig('defectidjustdef'+str(anglename)+'_min'+str(minimum+1)+'.png', dpi = 300)
#					print('Just Defect - Minimum '+str(minimum+1)+' saved: '+'defectidjustdef'+str(anglename)+'_min'+str(minimum+1)+'.png')
				fig.clf()
				plt.close()
	
	
				#5 - Directional vectors - find an opening wedge or not - angles between longest chains?
			
				#Calculating difference between every pair of chains		
				chainangles = []
				for chain in range(0, len(finalchainlist)):
					for chain2 in range(chain+1, len(finalchainlist)):
						#condition checking if chains neighbour each other
						neighbours = 0
						for particle in (finalchainlist[chain][0]):
							for particle2 in (finalchainlist[chain2][0]):
								if particle is not particle2:
									distances = np.sqrt((defectcoordsforanalysis[particle2][0] - defectcoordsforanalysis[particle][0])**2 + (defectcoordsforanalysis[particle2][1] - defectcoordsforanalysis[particle][1])**2)
									if distances < 1.70:
										neighbours += 1
						if neighbours > 3:
							if finalchainlist[chain][2] > finalchainlist[chain2][2]:
								anglebetweenchains = finalchainlist[chain][2]-finalchainlist[chain2][2]
								if anglebetweenchains < 45:
									chainangles.append([anglebetweenchains, chain, chain2])
								elif anglebetweenchains >= 135:
									chainangles.append([180-anglebetweenchains, chain, chain2])
							elif finalchainlist[chain][2] < finalchainlist[chain2][2]:
								anglebetweenchains = finalchainlist[chain2][2] - finalchainlist[chain][2]
								if anglebetweenchains < 45:
									chainangles.append([anglebetweenchains, chain2, chain])
								elif anglebetweenchains >= 135:
									chainangles.append([180-anglebetweenchains, chain2, chain])
				if categoriseall is 'N':
					print('Angles between neighbouring chains (and constituent chains): ', chainangles)
				
				#plotting defect lattice vector histogram
				fig, ax = plt.subplots()
				n, bins, patches = plt.hist(anglelist, bins=90, range = [0, 180])
				ax.set_xlabel('Angle in degrees')
				ax.set_ylabel('Frequency')
				plt.xticks(rotation = 45) #rotates xtick labels by 45 degrees
				plt.tight_layout() #stops x label being cut off
				ax.xaxis.set_ticks(np.arange(0, 180, 10))
				if categoriseall is 'N':
					fig.savefig('defectidvecthist'+str(anglename)+'_min'+str(inputminimum)+'.png')
					print('Histogram of defect lattice vectors saved: '+'defectidvecthist'+str(anglename)+'_min'+str(inputminimum)+'.png')
				fig.clf()
				plt.close()


			print(redxcoordsforanalysis, redycoordsforanalysis, bluexcoordsforanalysis, blueycoordsforanalysis)
			#-----------------------------------------------------------------------
			#print('---------------------')
			if categoriseall is 'N':
				print('DEFECT CATEGORISATION - minimum '+str(inputminimum))
			elif categoriseall is 'Y':
				print('DEFECT CATEGORISATION - minimum '+str(minimum+1))
		
			#types: helical (complete and starting at edge), normal wedge (across crystal or up/down crystal), crack through whole crystal (across and up/down separately), wedge and line joined together (L shape), lightning shape/bent wedge, no defect
		
			rholist.append([newrhostouse[rhoval], round(angleval/10, 1), defectcategorisation(auxlist, chainangles, yrange5, yrange6, redxcoordsforanalysis, bluexcoordsforanalysis, redycoordsforanalysis, blueycoordsforanalysis, defectxcoordsforanalysis, defectycoordsforanalysis, norminert, finalchainlist)])
		
			print('---------------------')

	phasediaglist.append(rholist)

#print(phasediaglist)

fig, ax = plt.subplots()

labelledlist = []
for i in range(0, len(newrhostouse)):
	anglelist = []
	rangelist = []
	defcatlist = []
	markertype = []
	colourtype = []
	labeltype = []
	for points in range(0, len(phasediaglist[i])):
		#print(phasediaglist[i][points])
		anglelist.append(phasediaglist[i][points][1]+0.1)
		rangelist.append(phasediaglist[i][points][0])
		defcatlist.append(phasediaglist[i][points][2])
		if phasediaglist[i][points][2] is 0: 	#unknown defect
			markertype.append('$?$')
			colourtype.append('red')
			labeltype.append('Unknown defect')
		elif phasediaglist[i][points][2] is 1:	#defectfree
			markertype.append('o')
			colourtype.append('g')
			labeltype.append('Defect free crystal')
		elif phasediaglist[i][points][2] is 6:	#complete helical defect
			markertype.append('X')
			colourtype.append('y')
			labeltype.append('Complete helical defect')
		elif phasediaglist[i][points][2] is 3:	#complete crack up and down cone
			markertype.append('P')
			colourtype.append('r')
			labeltype.append('Complete crack up and down cone')
		else:					#all other types
			markertype.append('*')
			colourtype.append('k')
			labeltype.append('Other Types')
	for values in range(0, len(anglelist)):
		if labeltype[values] not in labelledlist:
			plt.scatter(rangelist[values], anglelist[values], marker = markertype[values], s = 100, zorder = 300, color = colourtype[values], label = labeltype[values])
			labelledlist.append(labeltype[values])
		else:
			plt.scatter(rangelist[values], anglelist[values], marker = markertype[values], s = 100, zorder = 300, color = colourtype[values])

ax.legend(bbox_to_anchor=(1.05, 1))
plt.tight_layout()
fig.savefig('defectcatphasediag.png')
print('Angle vs. Range Defect Categorisation Phase Diagram saved: '+'defectcatphasediag.png')
