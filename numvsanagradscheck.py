#takes output file where numerical and analytical gradients are calculated by checkd and 
#outputs some important information

import statistics as stat

filename = input('Output Filename: ')
baderrorthreshold = input('Bad % error threshold: ')
printgrads = input('Print bad error gradients and % errors? Y/N ')

totalgrads = 0
warningcount = 0
baderrorcount = 0
graderrorlist = []
gradpercerrorlist = []
baderrorlist = []

with open(filename, 'r') as stream1:
	for line in stream1:
		if 'Gradient numerical' in line:
			totalgrads += 1
		elif 'WARNING' in line:
			warningcount += 1
			splitline = line.rstrip().split('  ')
			new = [item for item in splitline if item!='']
			abserror = float((new[-1]))
			if (abs(float(new[2]))) == 0:
				print('WARNING: Index ', new[1], ' has zero numerical gradient.')
#			elif (abs(float(new[2]))) < 1E-01:

			else:
				percerror = (abserror/abs(float(new[2])))*100
				if percerror > float(baderrorthreshold):
					baderrorcount += 1
					baderrorlist.append(new[1])
					if printgrads is 'Y':
						print('Grad ', new[1], ' : ', float(new[2]))
						print(percerror, '%')
				graderrorlist.append(float(splitline[-1]))
				gradpercerrorlist.append(percerror)
	
graderroravg = stat.mean(graderrorlist)
graderrorsd = stat.pstdev(graderrorlist)
graderrormin = min(graderrorlist)
graderrormax = max(graderrorlist)


gradpercerroravg = stat.mean(gradpercerrorlist)
gradpercerrorsd = stat.pstdev(gradpercerrorlist)
gradpercerrormin = min(gradpercerrorlist)
gradpercerrormax = max(gradpercerrorlist)


#print('Average error = ', graderroravg)
#print('S.D. error = ', graderrorsd)
#print('Min error = ', graderrormin)
#print('Max error = ', graderrormax)
print('Average % error = ', gradpercerroravg)
print('S.D. % error = ', gradpercerrorsd)
print('Min % error = ', gradpercerrormin)
print('Max % error = ', gradpercerrormax)
print('Warning count = ', warningcount, ' out of ', totalgrads)
print('Errors over threshold count = ', baderrorcount, ' out of ', totalgrads)
if printgrads is 'Y':
	print('Indices of bad errors: ', baderrorlist)
