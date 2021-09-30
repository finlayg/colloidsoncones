#script to take csv output file from stress.out and convert to xyz with rgb file for reading by xmakemol using given colour scheme

inputfile = input('.csv Input Filename: ')

intensitylist = []
coordslist = []

with open(inputfile, 'r') as stream1:
	for line in stream1:
		values = line.split(',')
		if len(values) is 3:
			is2D = True
		elif len(values) is 4:
			is2D = False
		else :
			is2D = None

with open(inputfile, 'r') as stream1:
	for line in stream1:
		values = line.split(',')
		intensity = values[-1]
		intensitylist.append(intensity)
		if is2D :
			coordslist.append([float(values[0]), values[1], 0.00])
		elif is2D is False :
                	coordslist.append([values[0],values[1],values[2]])

lines = len(open(inputfile).readlines(  ))

outputfile = str(lines)+' \n'+'text \n'

#converting intensity list to rgb

rgblist = []
highint = float(max(intensitylist))
lowint = float(min(intensitylist))

for j in range(0, len(intensitylist)):
	frac = (float(intensitylist[j]) - lowint)/(highint - lowint)
#	if frac < 0.0:
#		frac = 0.0
#	if frac > 1.0:
#		frac = 1.0
	if frac < 0.25:
		r = 0.0
		g = frac
		b = 1.0
	elif frac < 0.5:
		r = 0.0
		g = 1.0
		b = 1.0 - (frac - 0.25) * 4.0
	elif frac < 0.75:
		r = (frac - 0.5) * 4.0
		g = 1.0
		b = 0.0
	else :
		r = 1.0
		g = 1.0  - (frac - 0.75) * 4.0
		b = 0.0
	rgbvals = [r, g, b]
	rgblist.append(rgbvals)

#Constructing output file depending on whether coords file is net or cone

if is2D is False :
	for i in range(0, len(intensitylist)):
		output = 'C ' + str(coordslist[i][0]) + ' ' + str(coordslist[i][1]) + ' ' + str(coordslist[i][2]) + ' atom_rgb ' + str(rgblist[i][0]) + ' ' + str(rgblist[i][1]) + ' ' + str(rgblist[i][2]) + ' \n'
		outputfile = outputfile + output

elif is2D :
	for i in range(0, len(intensitylist)):
		output = 'C ' + str(coordslist[i][0]) + ' ' + str(coordslist[i][1]) + ' ' + str(coordslist[i][2]) + ' atom_rgb ' + str(rgblist[i][0]) + ' ' + str(rgblist[i][1]) + ' ' + str(rgblist[i][2]) + ' \n'
		outputfile = outputfile + output

converted_file = open(inputfile[:-4]+"_CONV.xyz", "w+")
converted_file.write(outputfile)
