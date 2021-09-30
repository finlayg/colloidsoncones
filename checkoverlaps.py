#checks if any particles in coords file are overlapping

inputfile = input('Coords Input Filename: ')

from math import sqrt

#generating particle coordinate list
particlelist = []
with open(inputfile, 'r') as stream1:
        for line in stream1:
                values = line.split(' ')
                splitline = line.rstrip().split(' ')
                new = [item for item in splitline if item!='']
                #print(new)
                particlelist.append([float(new[0]),float(new[1]),float(new[2])])

#testing distances between particles in list
for p1index in range(0, len(particlelist)):
	for p2index in range((p1index+1), len(particlelist)):
		distance = sqrt((particlelist[p1index][0]-particlelist[p2index][0])**2 + (particlelist[p1index][1]-particlelist[p2index][1])**2 + (particlelist[p1index][2]-particlelist[p2index][2])**2)
		#print(distance)
		if distance < 2:
			print("Warning: Overlap!")
