#takes lowest file and outputs specified minimum as coords-format file

inputfile = input('Lowest Input Filename: ')
inputminima = float(input('Minima coordinates to convert: '))
numparts = int(input('Number of particles: '))

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

#Constructing output file
outputfile = ''

for i in range(0, len(wantedcoordslist)):
        output = '   ' + str(wantedcoordslist[i][0]) + '   ' + str(wantedcoordslist[i][1]) + '   ' + str(wantedcoordslist[i][2]) + '\n'
        outputfile = outputfile + output

print("coordsconv file produced")

converted_file = open('coordsconv', "w+")
converted_file.write(outputfile)
