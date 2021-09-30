#takes output file and calculates statistics on b-h steps

import numpy as np

inputfile = 'output' #input('Lowest Input Filename: ')

nostepslist = []

with open(inputfile, 'r') as stream1:
	for line in stream1:
		if 'Qu     ' in line:
			splitline = line.rstrip().split(' ')
			new = [item for item in splitline if item!='']
			if len(new) is 13:
				nostepslist.append(int(new[5]))
			elif len(new) is 12:
				#print(new)
				splitnew = new[4].rstrip().split('=')
				#print(int(splitnew[1]))
				nostepslist.append(int(splitnew[1]))
print(np.mean(nostepslist))


