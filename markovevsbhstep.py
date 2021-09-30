#takes output file and plots Markov E vs b-h step


from matplotlib import pyplot as plt
import numpy as np

noruns = input('How many runs? 1/2/3 ')

if noruns is '1':
	inputfile  = 'output'

if noruns is '2':
	inputfile  = 'output'
	inputfile2 = 'firstrun/output'


if noruns is '3':
	inputfile  = 'secondrun/firstrun/output'
	inputfile2 = 'secondrun/output' 
	inputfile3 = 'output'


stepslist = []
markovElist = []
if noruns is '1':
	with open(inputfile, 'r') as stream1:
		for line in stream1:
			if 'Qu     ' in line:
				splitline = line.rstrip().split(' ')
				new = [item for item in splitline if item!='']
				if len(new) is 13:
					#print(new)
					stepslist.append(int(new[1]))
					markovElist.append(float(new[10]))
				elif len(new) is 12:
					#print(new)
					stepslist.append(int(new[1]))
					markovElist.append(float(new[9]))

if noruns is '2':
	with open(inputfile, 'r') as stream1:
		for line in stream1:
			if 'Qu     ' in line:
				splitline = line.rstrip().split(' ')
				new = [item for item in splitline if item!='']
				if len(new) is 13:
					#print(new)
					stepslist.append(int(new[1]))
					markovElist.append(float(new[10]))
				elif len(new) is 12:
					#print(new)
					stepslist.append(int(new[1]))
					markovElist.append(float(new[9]))

	add = stepslist[-1]	
	with open(inputfile1, 'r') as stream2:
		for line in stream2:
			if 'Qu     ' in line:
				splitline = line.rstrip().split(' ')
				new = [item for item in splitline if item!='']
				if len(new) is 13:
					#print(new)
					stepslist.append(add + int(new[1]))
					markovElist.append(float(new[10]))
				elif len(new) is 12:
					#print(new)
					stepslist.append(add + int(new[1]))
					markovElist.append(float(new[9]))

if noruns is '3':
	with open(inputfile, 'r') as stream1:
		for line in stream1:
			if 'Qu     ' in line:
				splitline = line.rstrip().split(' ')
				new = [item for item in splitline if item!='']
				if len(new) is 13:
					#print(new)
					stepslist.append(int(new[1]))
					markovElist.append(float(new[10]))
				elif len(new) is 12:
					#print(new)
					stepslist.append(int(new[1]))
					markovElist.append(float(new[9]))
	add = stepslist[-1]	
	with open(inputfile2, 'r') as stream2:
		for line in stream2:
			if 'Qu     ' in line:
				splitline = line.rstrip().split(' ')
				new = [item for item in splitline if item!='']
				if len(new) is 13:
					#print(new)
					stepslist.append(add + int(new[1]))
					markovElist.append(float(new[10]))
				elif len(new) is 12:
					#print(new)
					stepslist.append(add + int(new[1]))
					markovElist.append(float(new[9]))
	add = stepslist[-1]
	with open(inputfile3, 'r') as stream3:
		for line in stream3:
			if 'Qu     ' in line:
				splitline = line.rstrip().split(' ')
				new = [item for item in splitline if item!='']
				if len(new) is 13:
					#print(new)
					stepslist.append(add + int(new[1]))
					markovElist.append(float(new[10]))
				elif len(new) is 12:
					#print(new)
					stepslist.append(add + int(new[1]))
					markovElist.append(float(new[9]))

#print(stepslist[:100])
#print(markovElist[:100])

fig, ax = plt.subplots()

plt.ylim(-600, 0)

ax.plot(stepslist, markovElist)
fig.savefig('markovEplot.png', dpi = 300)


