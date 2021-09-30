#To plot energy vs angle from file containing first two lines of rigid from every angle constructed using 'head -2 */rigid.xyz > energyvsangle.dat'

from matplotlib import pyplot as plt
import numpy as np


multiplelines = input('Plot multiple ranges? Y/N ')
filenameslist = []

if multiplelines is 'N':
	inputfile = input('Energies input Filename: ')
	filenameslist.append(inputfile)

elif multiplelines is 'Y':
	listofranges = input('List the ranges to plot: ')
	splitline = listofranges.rstrip().split(',')
	new1 = [item for item in splitline if item!='']
	print(new1)
	inputfile = input('Energies input Filename root: ')
	for i in range(0, len(new1)):
		filenameslist.append(inputfile+new1[i])

	

energieslistlist = []

for j in range(0, len(filenameslist)):
	energieslist = []
	with open(filenameslist[j], 'r') as stream1:
		for line in stream1:
			if 'Energy of minimum' in line:
				splitline = line.rstrip().split(' ')
				new = [item for item in splitline if item!='']
				energieslist.append(float(new[4]))
	energieslistlist.append(energieslist)

#print(energieslist)

print('Number of ranges = ', len(energieslistlist))
print('Number of energies = ', len(energieslistlist[0]))

rangelower = float(input('Lower limit of angle range: '))
rangeupper = round(rangelower+(int(len(energieslist)-1)*0.1), 2)

print('Range of angles: ', rangelower, ' to ', rangeupper)

fig, ax = plt.subplots()

for k in range(0, len(energieslistlist)):
	linelabel = r'$\rho = $'+new1[k]
	if k is 0:
		linelstyle = 'solid'
	elif k is 1:
		linelstyle = 'dotted'
	elif k is 2:
		linelstyle = 'dashed'
	elif k is 3:
		linelstyle = 'dashdot'
	elif k is 4:
		linelstyle = (0, (3, 1, 1, 1, 1, 1))
	
	
	ax.plot(np.linspace(rangelower, rangeupper, num = len(energieslistlist[k]), endpoint = True), energieslistlist[k], marker = ',', linestyle = linelstyle, label = linelabel)

ax.legend()
ax.set_xlabel(r'Angle $\theta_p$ in degrees')
ax.set_ylabel(r'Energy / $\epsilon$')
ax.xaxis.set_ticks(np.arange(0, round(rangeupper, 1)+0.5, 0.5))
fig.savefig('energyvsangleplot.png', dpi=300)
fig.clf()
plt.close()

print('Saved energyvsangleplot.png')
