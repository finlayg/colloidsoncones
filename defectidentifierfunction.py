#defecttypes:
#0: unknown defect
#1: defect free
#2: complete crack separating distinct crystals
#3: complete crack up and down cone
#4: simple wedge up and down cone
#5: simple wedge across cone
#6: complete helical defect
#7: edge-terminating helical defect
#8: simple L-shaped hybrid defect - wedge across, helical up/down
#9: simple L-shaped hybrid defect - helical across, wedge up/down
#10: compound double wedge defect
#11: complex L-shaped hybrid defect - double wedge across, helical up/down
#12: complex L-shaped hybrid defect - helical across, double wedge up/down
#13: indeterminate L-shaped hybrid defect - helical up/down, unknown across
#14: amorphous localised defect

def defectcategorisation(auxlist, chainangles, yrange5, yrange6, redxcoordsforanalysis, bluexcoordsforanalysis, redycoordsforanalysis, blueycoordsforanalysis, defectxcoordsforanalysis, defectycoordsforanalysis, norminert, finalchainlist): #using function so can exit from if statements
	typeofdef = 0
	if len(auxlist) is 0:
		print('->No defect bigger than one particle.')
		print('*DEFECT TYPE: DEFECT-FREE CRYSTAL')
		typeofdef = 1
		return typeofdef
		
	else:
		if len(chainangles) is 0:
			print('->No neighbouring chains.')
			print('*DEFECT TYPE: DEFECT-FREE CRYSTAL')
			typeofdef = 1
			return typeofdef

		elif len(chainangles) is 1:
			print('->One pair of neighbouring chains.')
		
			if (len(redxcoordsforanalysis) > 2) or (len(bluexcoordsforanalysis) > 2):
				print('-->Defect contains more than two 5/7 side polygon pairs.')
				updownrangetop = max(yrange5) - min(yrange5)		
				updownrangebottom = max(yrange6) - min(yrange6)		
				yrangemid = (updownrangebottom + updownrangetop) / 2
				defectrange = max(defectycoordsforanalysis) - min(defectycoordsforanalysis)
				if (yrangemid - defectrange)**2 < 0.5:
					print('--->Particles in defect span vast majority of space across cone.')
					print('*DEFECT TYPE: COMPLETE CRACK ACROSS CONE SEPARATING DISTINCT CRYSTALS')
					typeofdef = 2
					return typeofdef
				else:
					print('--->Particles in defect do not span vast majority of space across cone.')
					if norminert > 10:
						print('---->Defect has normalised inertia of greater than 10')
						print('*DEFECT TYPE: COMPLETE CRACK UP AND DOWN CONE')
						typeofdef = 3
						return typeofdef
					else:
						print('---->Defect has normalised inertia of less than 10')
						print('*DEFECT TYPE: AMORPHOUS LOCALISED DEFECT')
						typeofdef = 14
						return typeofdef
	
			if chainangles[0][0] > 1:
				print('-->Neighbouring chains have opening angle of more than 1 degree.')
				chainwedge1 = chainangles[0][1]
				chainwedge2 = chainangles[0][2]
				if (finalchainlist[chainwedge1][2] > 135) or (finalchainlist[chainwedge1][2] < 45):
					print('--->Chains point up and down cone.')
					print('*DEFECT TYPE: SIMPLE WEDGE - UP AND DOWN CONE')
					typeofdef = 4
					return typeofdef
				else:
					print('--->Chains point across cone.')
					print('*DEFECT TYPE: SIMPLE WEDGE - ACROSS CONE')
					typeofdef = 5
					return typeofdef
			elif chainangles[0][0] < 1:
				print('-->Neighbouring chains have opening angle of less than 1 degree.')
				if (len(redxcoordsforanalysis) > 0) and (len(bluexcoordsforanalysis) > 0):
					print('--->Defect contains a 5/7 side polygon pair.')
					if len(finalchainlist) is 2:
						print('---->Defect contains the 2 chains that make the helix only.')
						print('*DEFECT TYPE: COMPLETE HELICAL DEFECT')
						typeofdef = 6
						return typeofdef
					elif len(finalchainlist) > 2:
						print('---->Defect contains more chains than the 2 that make the helix.')
						print('*DEFECT TYPE: INDETERMINATE L-SHAPED HYBRID DEFECT - HELIX UP AND DOWN CONE AND UNKNOWN ACROSS')
						typeofdef = 13
						return typeofdef
				else:
					print('--->Defect does not contain a 5/7 side polygon pair.')
					print('*DEFECT TYPE: EDGE-TERMINATING HELICAL DEFECT')
					typeofdef = 7
					return typeofdef

		elif len(chainangles) > 1:
			print('->Two or more pairs of neighbouring chains.')	
			if (len(redxcoordsforanalysis) > 2) or (len(bluexcoordsforanalysis) > 2):
				print('-->Defect contains more than two 5/7 side polygon pairs.')
				updownrangetop = max(yrange5) - min(yrange5)		
				updownrangebottom = max(yrange6) - min(yrange6)		
				yrangemid = (updownrangebottom + updownrangetop) / 2
				defectrange = max(defectycoordsforanalysis) - min(defectycoordsforanalysis)
				if (yrangemid - defectrange)**2 < 0.35:
					print('--->Particles in defect span majority of space across cone.')
					print('*DEFECT TYPE: COMPLETE CRACK ACROSS CONE SEPARATING DISTINCT CRYSTALS')
					typeofdef = 2
					return typeofdef
			if (len(chainangles) is 2):
				print('-->Two pairs of neighbouring chains.')
				chainangleslist = [chainangles[0][0], chainangles[1][0]]
				if (max(chainangleslist) < 10) and (max(chainangleslist) > 1) and (min(chainangleslist) < 1):
					print('--->Largest angle between neighbouring chains is between 1 and 10 degrees and smallest is less than 1 degree.')
					for chainangle in range(0, len(chainangles)):
						if chainangles[chainangle][0] < 1:
							chainwedge1 = chainangles[chainangle][1] #taking one of the chains in helical part of defect and seeing if point up or across cone
							if (finalchainlist[chainwedge1][2] > 135) or (finalchainlist[chainwedge1][2] < 45):
								print('---->Helical end of defect points up and down cone.')
								print('*DEFECT TYPE: SIMPLE L-SHAPED HYBRID DEFECT - WEDGE ACROSS CONE AND HELICAL UP AND DOWN CONE')
								typeofdef = 8
								return typeofdef
							else: 
								print('---->Helical end of defect points across cone.')
								print('*DEFECT TYPE: SIMPLE L-SHAPED HYBRID DEFECT - WEDGE UP AND DOWN CONE AND HELICAL ACROSS CONE')
								typeofdef = 9
								return typeofdef
				elif (max(chainangleslist) < 10) and (max(chainangleslist) > 1) and (min(chainangleslist) < 10) and (min(chainangleslist) > 1):
					print('--->Both angles between neighbouring chains are between 1 and 10 degrees.')
					print('*DEFECT TYPE: COMPOUND DOUBLE WEDGE DEFECT')
					typeofdef = 10
					return typeofdef

			elif (len(chainangles) is 3):
				print('-->Three pairs of neighbouring chains.')
				chainangleslist = [chainangles[0][0], chainangles[1][0], chainangles[2][0]]
				chainangleslist.sort()
				if (max(chainangleslist) < 10) and (max(chainangleslist) > 1) and (min(chainangleslist) < 1) and ((chainangleslist[1] < 10) and (chainangleslist[1] > 1)):
					print('--->Largest and second largest angles between neighbouring chains is between 1 and 10 and smallest is less than 1.')
					for chainangle in range(0, len(chainangles)):
						if chainangles[chainangle][0] < 1:
							chainwedge1 = chainangles[chainangle][1] #taking one of the chains in helical part of defect and seeing if point up or across cone
							if (finalchainlist[chainwedge1][2] > 135) or (finalchainlist[chainwedge1][2] < 45):
								print('---->Helical end of defect points up and down cone.')
								print('*DEFECT TYPE: COMPLEX L-SHAPED HYBRID DEFECT - DOUBLE WEDGE ACROSS CONE AND HELICAL UP AND DOWN CONE')
								typeofdef = 11
								return typeofdef
							else: 
								print('---->Helical end of defect points across cone.')
								print('*DEFECT TYPE: COMPLEX L-SHAPED HYBRID DEFECT - DOUBLE WEDGE UP AND DOWN CONE AND HELICAL ACROSS CONE')
								typeofdef = 12
								return typeofdef
			if (len(redxcoordsforanalysis) > 2) or (len(bluexcoordsforanalysis) > 2):
				print('--->Particles in defect do not span majority of space across cone.')
				if norminert > 10:
					print('---->Defect has normalised inertia of greater than 10')
					print('*DEFECT TYPE: COMPLETE CRACK UP AND DOWN CONE')
					typeofdef = 3
					return typeofdef
				else:
					print('---->Defect has normalised inertia of less than 10')
					print('*DEFECT TYPE: AMORPHOUS LOCALISED DEFECT')
					typeofdef = 14
					return typeofdef
		print('*DEFECT TYPE UNKNOWN')
		return typeofdef
