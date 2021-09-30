#Script to take coords input file and convert to .xyz file for visualisation of initial coordinates

inputfile = input('Coords Input Filename: ')

coordslist = []

with open(inputfile, 'r') as stream1:
	for line in stream1:
		values = line.split(' ')
		splitline = line.rstrip().split(' ')
		new = [item for item in splitline if item!='']
		#print(new)
		coordslist.append([new[0],new[1],new[2]])

lines = len(open(inputfile).readlines(  ))

outputfile = str(lines)+' \n'+'text \n'

#Constructing output file

for i in range(0, len(coordslist)):
	output = 'C ' + str(coordslist[i][0]) + ' ' + str(coordslist[i][1]) + ' ' + str(coordslist[i][2]) + '\n' 
	outputfile = outputfile + output

print(inputfile+"_CONV.xyz file produced")

converted_file = open(inputfile+"_CONV.xyz", "w+")
converted_file.write(outputfile)
