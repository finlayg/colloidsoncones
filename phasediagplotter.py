#script to identify type of defect from given lowest input files and plot phase diagram of defect type on angle vs range plot
#Requires file tree to be: rho > angle > lowest and data

from scipy.spatial import Voronoi, voronoi_plot_2d
from matplotlib import pyplot as plt
from matplotlib import markers
import numpy as np
from collections import Counter
from pathlib import Path
from defectidentifierfunction import defectcategorisation
from defectidentification import defectidentification

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
#print(anglelistz)

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

			defectidentifiervariablelist = defectidentification(minimalist, l, L, angle, colour, lengthtolerance, anglep, categoriseall, inputminimum)
			auxlist = defectidentifiervariablelist[0]
			chainangles = defectidentifiervariablelist[1]
			yrange5 = defectidentifiervariablelist[2]
			yrange6 = defectidentifiervariablelist[3]
			redxcoordsforanalysis = defectidentifiervariablelist[4]
			bluexcoordsforanalysis = defectidentifiervariablelist[5]
			redycoordsforanalysis = defectidentifiervariablelist[6]
			blueycoordsforanalysis = defectidentifiervariablelist[7]
			defectxcoordsforanalysis = defectidentifiervariablelist[8]
			defectycoordsforanalysis = defectidentifiervariablelist[9]
			norminert = defectidentifiervariablelist[10]
			finalchainlist = defectidentifiervariablelist[11]

			#print(redxcoordsforanalysis, redycoordsforanalysis, bluexcoordsforanalysis, blueycoordsforanalysis)
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
	rholist = []	

#manualphasediagram
phasediaglist = [ [ [6 ,0.1,1 ],[6 ,0.2,1 ],[6 ,0.3,1 ],[6 ,0.4,1 ],[6 ,0.5,6 ],[6 ,0.6,6 ],[6 ,0.7,6 ],[6 ,0.8,8 ],[6 ,0.9,6 ],[6 ,1.0,6 ],    
                    [6 ,1.1,6 ],[6 ,1.2,25],[6 ,1.3,25],[6 ,1.4,25],[6 ,1.5,17],[6 ,1.6,16],[6 ,1.7,16],[6 ,1.8,16],[6 ,1.9,16],[6 ,2.0,15], 
                    [6 ,2.1,15],[6 ,2.2,15],[6 ,2.3,15],[6 ,2.4,15],[6 ,2.5,15],[6 ,2.6,15],[6 ,2.7,22],[6 ,2.8,22],[6 ,2.9,22],[6 ,3.0,22] ], 
                  [ [9 ,0.1,1 ],[9 ,0.2,1 ],[9 ,0.3,1 ],[9 ,0.4,5 ],[9 ,0.5,6 ],[9 ,0.6,8 ],[9 ,0.7,6 ],[9 ,0.8,6 ],[9 ,0.9,6 ],[9 ,1.0,6 ],    
                    [9 ,1.1,6 ],[9 ,1.2,6 ],[9 ,1.3,6 ],[9 ,1.4,17],[9 ,1.5,17],[9 ,1.6,23],[9 ,1.7,23],[9 ,1.8,19],[9 ,1.9,17],[9 ,2.0,23], 
                    [9 ,2.1,19],[9 ,2.2,18],[9 ,2.3,19],[9 ,2.4,19],[9 ,2.5,24],[9 ,2.6,24],[9 ,2.7,24],[9 ,2.8,24],[9 ,2.9,24],[9 ,3.0,24] ], #note: wedges much shorter
                  [ [12,0.1,1 ],[12,0.2,1 ],[12,0.3,1 ],[12,0.4,5 ],[12,0.5,5 ],[12,0.6,13],[12,0.7,6 ],[12,0.8,6 ],[12,0.9,6 ],[12,1.0,13],
                    [12,1.1,6 ],[12,1.2,6 ],[12,1.3,10],[12,1.4,17],[12,1.5,18],[12,1.6,19],[12,1.7,19],[12,1.8,19],[12,1.9,19],[12,2.0,19],
                    [12,2.1,19],[12,2.2,19],[12,2.3,24],[12,2.4,19],[12,2.5,19],[12,2.6,19],[12,2.7,19],[12,2.8,21],[12,2.9,21],[12,3.0,21] ],
                  [ [15,0.1,1 ],[15,0.2,1 ],[15,0.3,7 ],[15,0.4,5 ],[15,0.5,5 ],[15,0.6,8 ],[15,0.7,6 ],[15,0.8,6 ],[15,0.9,6 ],[15,1.0,17],
                    [15,1.1,6 ],[15,1.2,9 ],[15,1.3,10],[15,1.4,18],[15,1.5,23],[15,1.6,19],[15,1.7,19],[15,1.8,19],[15,1.9,19],[15,2.0,19],
                    [15,2.1,10],[15,2.2,19],[15,2.3,19],[15,2.4,19],[15,2.5,19],[15,2.6,21],[15,2.7,21],[15,2.8,21],[15,2.9,21],[15,3.0,21] ],
                  [ [18,0.1,1 ],[18,0.2,1 ],[18,0.3,7 ],[18,0.4,5 ],[18,0.5,8 ],[18,0.6,8 ],[18,0.7,6 ],[18,0.8,6 ],[18,0.9,17],[18,1.0,17],
                    [18,1.1,17],[18,1.2,6 ],[18,1.3,10],[18,1.4,23],[18,1.5,12],[18,1.6,19],[18,1.7,19],[18,1.8,19],[18,1.9,21],[18,2.0,19],
                    [18,2.1,19],[18,2.2,21],[18,2.3,21],[18,2.4,19],[18,2.5,21],[18,2.6,21],[18,2.7,21],[18,2.8,3 ],[18,2.9,3 ],[18,3.0,3 ]  ] ]
#                  [ [18,0.1,1 ],[18,0.2,7 ],[18,0.3,7 ],[18,0.4,0 ],[18,0.5,0 ],[18,0.6,6 ],[18,0.7,6 ],[18,0.8,6 ],[18,0.9,0 ]            , #these are repeat2, above is compiled from sdi
#                    [18,1.1,0 ],[18,1.2,0 ],[18,1.3,9 ],[18,1.4,0 ],[18,1.5,0 ],[18,1.6,0 ],[18,1.7,19],[18,1.8,0 ],[18,1.9,0 ],[18,2.0,19],
#                    [18,2.1,21],[18,2.2,21],[18,2.3,21],[18,2.4,19],[18,2.5,21],[18,2.6,21]                                                 ] ]
newrhostouse = [6,9,12,15,18]

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
		anglelist.append(phasediaglist[i][points][1])#+0.1)
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
		elif (phasediaglist[i][points][2] is 4) or (phasediaglist[i][points][2] is 5):	#simple wedge
			markertype.append('$<$')
			colourtype.append('y')
			labeltype.append('Simple wedge defect')
		elif phasediaglist[i][points][2] is 10:	#compound double wedge
			markertype.append('$w$')
			colourtype.append('y')
			labeltype.append('Compound double wedge defect')
		elif phasediaglist[i][points][2] is 6:	#complete helical defect
			markertype.append('$c$')
			colourtype.append('y')
			labeltype.append('Complete helical defect')
		elif phasediaglist[i][points][2] is 7:	#edge-terminating helical defect
			markertype.append('$u$')
			colourtype.append('y')
			labeltype.append('Edge-terminating helical defect')
		elif phasediaglist[i][points][2] is 2:	#ccomplete crack across cone
			markertype.append('|')
			colourtype.append('y')
			labeltype.append('Complete crack across cone')
		elif phasediaglist[i][points][2] is 3:	#complete crack up and down cone
			markertype.append('_')
			colourtype.append('r')
			labeltype.append('Complete crack up and down cone')
		elif phasediaglist[i][points][2] is 21:	#verging on complete crack up and down cone
			markertype.append('_')
			colourtype.append('orange')
			labeltype.append('Almost complete crack up and down cone')
		elif phasediaglist[i][points][2] is 15:	#2 pentagonal defects
			markertype.append('$dd$')
			colourtype.append('y')
			labeltype.append(0)
#			labeltype.append('Dislocation defects')
		elif phasediaglist[i][points][2] is 22: #3 pentagonal defects
			markertype.append('$ddd$')
			colourtype.append('y')
			labeltype.append(0)
#			labeltype.append('Three dislocation defects')
		elif phasediaglist[i][points][2] is 17:	#multiple helical defects
			markertype.append('$cc$')
			colourtype.append('y')
			labeltype.append(0)
#			labeltype.append('Multiple complete helical defects')
		elif phasediaglist[i][points][2] is 16:	#complete helical and pentagonal defects coexisting
			markertype.append('$cd$')
			colourtype.append('y')
			labeltype.append(0)
#			labeltype.append('Complete helical/dislocation defects coexisting')
		elif phasediaglist[i][points][2] is 18:	#complete helical and double wedge defects coexisting
			markertype.append('$cw$')
			colourtype.append('y')
			labeltype.append(0)
#			labeltype.append('Complete helical/double wedge defects coexisting')
		elif phasediaglist[i][points][2] is 19:	#2 - multiple wedge defects 
			markertype.append('$<<$')
			colourtype.append('y')
			labeltype.append(0)
#			labeltype.append('Two wedge defects')
		elif phasediaglist[i][points][2] is 20:	#Wedge and pentagons coexisting
			markertype.append('$<d$')
			colourtype.append('y')
			labeltype.append(0)
#			labeltype.append('Wedge/dislocation defects coexisting')
		elif phasediaglist[i][points][2] is 23:	#complete helical and simple wedge defects coexisting
			markertype.append('$c<$')
			colourtype.append('y')
			labeltype.append(0)
#			labeltype.append('Complete helical/simple wedge defects coexisting')
		elif phasediaglist[i][points][2] is 24:	#3 wedge defects coexisting
			markertype.append('$<<<$')
			colourtype.append('y')
			labeltype.append(0)
#			labeltype.append('Three wedge defects')
		elif phasediaglist[i][points][2] is 25:	#1 pentagonal dislocation defect
			markertype.append('$d$')
			colourtype.append('y')
#			labeltype.append(0)
			labeltype.append('Dislocation defect')
		elif (phasediaglist[i][points][2] is 8 or 9 or 11 or 12):	#hybrid L defect
			markertype.append('$L$')
			colourtype.append('y')
			labeltype.append('Hybrid L-shaped defect')
		else:					#all other types
			markertype.append('*')
			colourtype.append('k')
			labeltype.append('Other Types')
	for values in range(0, len(anglelist)):
		if (labeltype[values] not in labelledlist) and (labeltype[values] is not 0):
			plt.scatter(rangelist[values], anglelist[values], marker = markertype[values], s = 100, zorder = 300, color = colourtype[values], label = labeltype[values])
			labelledlist.append(labeltype[values])
		else:
			plt.scatter(rangelist[values], anglelist[values], marker = markertype[values], s = 100, zorder = 300, color = colourtype[values])


plt.xlabel(r'Potential range, $\rho$')
plt.ylabel(r'Angle $\theta_{p}$ in degrees')
plt.xlim([min(newrhostouse)-2, max(newrhostouse)+2])
plt.xticks(newrhostouse)
ax.legend(bbox_to_anchor=(1.05, 1))
plt.tight_layout()
fig.savefig('defectcatphasediag.png', dpi = 300)
print('Angle vs. Range Defect Categorisation Phase Diagram saved: '+'defectcatphasediag.png')
