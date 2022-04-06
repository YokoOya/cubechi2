import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from math import *


saveDir = 'chi2plot_varRin-limI60-120-ev10/'

beamsize = 0.832 * 0.488 / 0.01**2 # FWHM_major (as) * FWHM_minor (as) / pixsize^2
velwidth = 1.0 / 0.2 # FWHM (km/s) / pixsize
num_indepPix = beamsize / 4. * velwidth / 2.

f_prodIndepPix = False	#True
f_colorLog = True


def colorArray_log(temp_max, temp_min, chi2pDOFmax, chi2pDOFmin):
	base = 10
	while 1:
		digit_min = floor(log(temp_min) / log(base))
		digit_max = ceil(log(temp_max) / log(base))
		if digit_max - digit_min < 20:
			base = sqrt(base)
		elif digit_max - digit_min > 100:
			base = base**2
		else:
			break
	
	return np.logspace(digit_min, digit_max, digit_max - digit_min + 1, base = base)
	

def colorArray_lin(temp_max, temp_min, chi2pDOFmax, chi2pDOFmin):
	ndigit = floor(log(temp_max - temp_min) / log(10))
	dcolor = 10**(ndigit - 1)

	numContour = (temp_max - temp_min) / dcolor
	if numContour > 100:	#50:
		dcolor = dcolor * 10
	elif numContour > 40:	#20:
		dcolor = dcolor * 5
	elif numContour > 20:	#10:
		dcolor = dcolor * 2
	
	return np.arange(floor(chi2pDOFmin / dcolor) * dcolor, chi2pDOFmax, dcolor)


def calcColorbarArray(temp_max, temp_min, chi2pDOFmax, chi2pDOFmin):
	if f_colorLog:
		return colorArray_log(temp_max, temp_min, chi2pDOFmax, chi2pDOFmin)
	return colorArray_lin(temp_max, temp_min, chi2pDOFmax, chi2pDOFmin)


try:
	os.makedirs(saveDir)
except:
	pass


for i in range(4):
	case = i
	
	if case == 0:
		infilename = 'Elias29_SO_onlyIRE_varRin-limI60-120-ev10.out'
		MList = ['0.2', '0.4', '0.6', '0.8', '1.0', '1.2', '1.4', '1.6', '1.8', '2.0', '3.0', '4.0', '5.0', '6.0', '8.0', '10.0']
		CBList = ['1.', '5.', '10.', '20.', '30.', '40.']	#['1.', '10.', '20.', '30.', '40.', '50.', '60.', '70.', '80.', '90.']
		IList = ['60', '70', '80', '90', '100', '110', '120']
		RinList = ['CB']
	elif case == 1:
		continue
		infilename = 'Elias29_SO_onlyKep-Rot1_varRin-limI60-120-ev10.out'
		MList = ['0.2', '0.4', '0.6', '0.8', '1.0', '1.2', '1.4', '1.6', '1.8', '2.0', '3.0', '4.0', '5.0', '6.0', '8.0', '10.0']
		CBList = ['100.']
		IList = ['60', '70', '80', '90', '100', '110', '120']
		RinList = ['1.', '5.', '10.', '20.', '30.', '40.']	#['0', '10.', '20.', '30.', '40.', '50.', '60.', '70.', '80.', '90.', '100.'] #['CB'] #
	elif case == 2:
		infilename = 'Elias29_SO_onlyKep-Rot-1_varRin-limI60-120-ev10.out'
		MList = ['0.2', '0.4', '0.6', '0.8', '1.0', '1.2', '1.4', '1.6', '1.8', '2.0', '3.0', '4.0', '5.0', '6.0', '8.0', '10.0']
		CBList = ['100.']
		IList = ['60', '70', '80', '90', '100', '110', '120']
		RinList = ['1.', '5.', '10.', '20.', '30.', '40.']	#['0', '10.', '20.', '30.', '40.', '50.', '60.', '70.', '80.', '90.', '100.'] #['CB'] #
	elif case == 3:
		infilename = 'Elias29_SO_IRE-Kep_varRin-limI60-120-ev10.out'
		MList = ['0.2', '0.4', '0.6', '0.8', '1.0', '1.2', '1.4', '1.6', '1.8', '2.0', '3.0', '4.0', '5.0', '6.0', '8.0', '10.0']
		CBList = ['1.', '5.', '10.', '20.', '30.', '40.']	#['1.', '10.', '20.', '30.', '40.', '50.', '60.', '70.', '80.', '90.']
		IList = ['60', '70', '80', '90', '100', '110', '120']
		RinList = ['1.']
	else :
		print "Choose '.out' file to be input."
		sys.exit(1)
	
	print "\nProcessing: %s\n" % infilename
	saveepsname = saveDir + infilename + '.'


	#############################


	nM = len(MList)
	nCB = len(CBList)
	nI = len(IList)
	nRin = len(RinList)

	fin = open(infilename, 'r')

	chi2regionfltTable = []
	chi2pDOFregionfltTable = []
	chi2thresholdfltTable = []
	chi2pDOFthresholdfltTable = []


	fileList = []


	searchStr_region = "In the specified region      #Cell = "
	searchStr_threshold = "Larger than the threshold:   #Cell = "
	searchStr1 = "  Sum of (Cube1 - Cube2)^2 / rms^2 =      "
	searchStr2 = "                              /DOF =      "
	searchStr_filename = "Saved: "

	fregion = 1
	fthreshold = 2



	flag_regionthreshold = 0

	for line in fin:
		if line[:len(searchStr_region)] == searchStr_region:
			flag_regionthreshold = fregion
		if line[:len(searchStr_threshold)] == searchStr_threshold:
			flag_regionthreshold = fthreshold
		
		if flag_regionthreshold == fregion:
			if line[:len(searchStr1)] == searchStr1:
				val = (float)(line[len(searchStr1):])
				chi2regionfltTable.append(val)
			if line[:len(searchStr2)] == searchStr2:
				val = (float)(line[len(searchStr2):])
				chi2pDOFregionfltTable.append(val)
		
		if flag_regionthreshold == fthreshold:
			if line[:len(searchStr1)] == searchStr1:
				val = (float)(line[len(searchStr1):])
				chi2thresholdfltTable.append(val)
			if line[:len(searchStr2)] == searchStr2:
				val = (float)(line[len(searchStr2):])
				chi2pDOFthresholdfltTable.append(val)
		
		if flag_regionthreshold == fregion and line[:len(searchStr_filename)] == searchStr_filename:
			fileList.append(line[len(searchStr_filename):])


	chi2regionfltTable = np.array(chi2regionfltTable)
	chi2regionfltTable = chi2regionfltTable.reshape(nM, nCB, nI, nRin)
	chi2pDOFregionfltTable = np.array(chi2pDOFregionfltTable)
	chi2pDOFregionfltTable = chi2pDOFregionfltTable.reshape(nM, nCB, nI, nRin)

	chi2thresholdfltTable = np.array(chi2thresholdfltTable)
	chi2thresholdfltTable = chi2thresholdfltTable.reshape(nM, nCB, nI, nRin)
	chi2pDOFthresholdfltTable = np.array(chi2pDOFthresholdfltTable)
	chi2pDOFthresholdfltTable = chi2pDOFthresholdfltTable.reshape(nM, nCB, nI, nRin)




	chi2pDOFfltTable = chi2pDOFregionfltTable			#plot the chi-squared value in the specified region
	#chi2pDOFfltTable = chi2pDOFthresholdfltTable		#plot the chi-squared value for the data points with a intensity larger than the threshold
	
	if f_prodIndepPix:
		chi2pDOFfltTable = chi2pDOFfltTable * num_indepPix
		saveepsname = saveepsname + 'indepPixconsidered.'
	if f_colorLog:
		saveepsname = saveepsname + 'logColor.'
	
	
	### plot ###
	"""
	chi2pDOFmin = 1e8
	chi2pDOFmax = -1e8
	for iRin in range(nRin):
		for iCB in range(nCB):
			for iI in range(nI):
				for iM in range(nM):
					val = chi2pDOFfltTable[iM][iCB][iI][iRin]
					if chi2pDOFmin > val:
						chi2pDOFmin = val
					if chi2pDOFmax < val:
						chi2pDOFmax = val
	"""
	
	chi2pDOFmin = np.nanmin(chi2pDOFfltTable)
	chi2pDOFmax = np.nanmax(chi2pDOFfltTable)


	MList_flt = []
	CBList_flt = []
	IList_flt = []
	RinList_flt = []

	for M in MList:
		MList_flt.append(float(M))
	for CB in CBList:
		CBList_flt.append(float(CB))
	for I in IList:
		IList_flt.append(float(I))
	
	if RinList[0] == 'CB':
		for CB in CBList:
			RinList_flt.append(float(CB))
	else:
		for Rin in RinList:
			RinList_flt.append(float(Rin))

	Marray = np.array(MList_flt)
	CBarray = np.array(CBList_flt)
	Iarray = np.array(IList_flt)
	Rinarray = np.array(RinList_flt)



	### color bar ###
	"""
	ndigit = floor(log(chi2pDOFmax - chi2pDOFmin) / log(10))
	dcolor = 10**(ndigit - 1)

	numContour = (chi2pDOFmax - chi2pDOFmin) / dcolor
	if numContour > 50:
		dcolor = dcolor * 5
	elif numContour > 20:
		dcolor = dcolor * 2
	"""
	
	### for M ###
	
	if nCB > 1 and nI > 1:
		X, Y = np.meshgrid(CBarray, Iarray)

		Zmin = np.full((nI, nCB), np.nan)
		
		for iM in range(nM):
			ZList = []
			for iI in range(nI):
				ZList.append([])
				for iCB in range(nCB):
					val = np.nanmin(chi2pDOFfltTable[iM][iCB][iI])
					ZList[iI].append(val)
					if np.isnan(Zmin[iI][iCB]) or Zmin[iI][iCB] > val:
						Zmin[iI][iCB] = val

			Z = np.array(ZList)
			
			temp_max = np.nanmax(Z)
			temp_min = np.nanmin(Z)
			
			"""
			dcolor = calcDcolor(temp_max, temp_min)
			if np.isnan(dcolor):
				continue
			
			colorbarArray = np.arange(floor(chi2pDOFmin / dcolor) * dcolor, chi2pDOFmax, dcolor)
			"""
			colorbarArray = calcColorbarArray(temp_max, temp_min, chi2pDOFmax, chi2pDOFmin)
			
			fig = plt.figure()
			ax = fig.add_subplot(1, 1, 1)
			cont = ax.contour(X, Y, Z, colorbarArray)
			
			cont.clabel(fmt = '%1.2f', fontsize = 12)
			ax.set_title('Chi-Squared Value (Central Mass ' + MList[iM] + ' Msun; Rin from ' + RinList[0] + ' to ' + RinList[-1] + ' au)')
			plt.xlabel('Radius of the Centrifugal Barrier (au)', fontsize = 14)
			plt.ylabel('Inclination Angle (degree)', fontsize = 14)
			
			plt.colorbar(cont)
			
			plt.savefig(saveepsname + 'chi2_CBvsI_allRin_M' + MList[iM] + '.eps')
			plt.close(fig)
		
		temp_max = np.nanmax(Zmin)
		temp_min = np.nanmin(Zmin)
		"""
		dcolor = calcDcolor(temp_max, temp_min)
		colorbarArray = np.arange(floor(chi2pDOFmin / dcolor) * dcolor, chi2pDOFmax, dcolor)
		"""
		colorbarArray = calcColorbarArray(temp_max, temp_min, chi2pDOFmax, chi2pDOFmin)
		
		fig = plt.figure()
		ax = fig.add_subplot(1, 1, 1)
		cont = ax.contour(X, Y, Zmin, colorbarArray)
		cont.clabel(fmt = '%1.2f', fontsize = 12)
		ax.set_title('Chi-Squared Value (Central Mass from ' + MList[0] + ' to ' + MList[-1] + ' Msun; Rin from ' + RinList[0] + ' to ' + RinList[-1] + ' au)')
		plt.xlabel('Radius of the Centrifugal Barrier (au)', fontsize = 14)
		plt.ylabel('Inclination Angle (degree)', fontsize = 14)
		plt.colorbar(cont)

		plt.savefig(saveepsname + 'chi2_CBvsI_allRin_allM.eps')
		plt.close(fig)
	

	### for CB ###

	if nM > 1 and nI > 1:
		X, Y = np.meshgrid(Marray, Iarray)

		Zmin = np.full((nI, nM), np.nan)

		for iCB in range(nCB):
			ZList = []
			for iI in range(nI):
				ZList.append([])
				for iM in range(nM):
					val = np.nanmin(chi2pDOFfltTable[iM][iCB][iI])
					ZList[iI].append(val)
					if np.isnan(Zmin[iI][iM]) or Zmin[iI][iM] > val:
						Zmin[iI][iM] = val
			Z = np.array(ZList)
			
			temp_max = np.nanmax(Z)
			temp_min = np.nanmin(Z)
			"""
			dcolor = calcDcolor(temp_max, temp_min)
			if isnan(dcolor):
				continue
			
			colorbarArray = np.arange(floor(chi2pDOFmin / dcolor) * dcolor, chi2pDOFmax, dcolor)
			"""
			colorbarArray = calcColorbarArray(temp_max, temp_min, chi2pDOFmax, chi2pDOFmin)
			
			fig = plt.figure()
			ax = fig.add_subplot(1, 1, 1)
			cont = ax.contour(X, Y, Z, colorbarArray)
			
			cont.clabel(fmt = '%1.2f', fontsize = 12)
			ax.set_title('Chi-Squared Value (Centrifugal Barrier ' + CBList[iCB] + ' au; Rin from ' + RinList[0] + ' to ' + RinList[-1] + ' au)')
			plt.xlabel('Central Mass (Msun)', fontsize = 14)
			plt.ylabel('Inclination Angle (degree)', fontsize = 14)
			plt.colorbar(cont)
			
			plt.savefig(saveepsname + 'chi2_MvsI_allRin_CB' + CBList[iCB] + '.eps')
			plt.close(fig)

		temp_max = np.nanmax(Zmin)
		temp_min = np.nanmin(Zmin)
		"""
		dcolor = calcDcolor(temp_max, temp_min)
		colorbarArray = np.arange(floor(chi2pDOFmin / dcolor) * dcolor, chi2pDOFmax, dcolor)
		"""
		colorbarArray = calcColorbarArray(temp_max, temp_min, chi2pDOFmax, chi2pDOFmin)

		fig = plt.figure()
		ax = fig.add_subplot(1, 1, 1)
		cont = ax.contour(X, Y, Zmin, colorbarArray)
		cont.clabel(fmt = '%1.2f', fontsize = 12)
		ax.set_title('Chi-Squared Value (Centrifugal Barrier from ' + CBList[0] + ' to ' + CBList[-1] + ' au; Rin from ' + RinList[0] + ' to ' + RinList[-1] + ' au)')
		plt.xlabel('Central Mass (Msun)', fontsize = 14)
		plt.ylabel('Inclination Angle (degree)', fontsize = 14)
		plt.colorbar(cont)

		plt.savefig(saveepsname + 'chi2_MvsI_allRin_allCB.eps')
		plt.close(fig)

	### for I ###

	if nCB > 1 and nM > 1:
		X, Y = np.meshgrid(CBarray, Marray)

		Zmin = np.full((nM, nCB), np.nan)

		for iI in range(nI):
			ZList = []
			for iM in range(nM):
				ZList.append([])
				for iCB in range(nCB):
					val = np.nanmin(chi2pDOFfltTable[iM][iCB][iI])
					ZList[iM].append(val)
					if np.isnan(Zmin[iM][iCB]) or Zmin[iM][iCB] > val:
						Zmin[iM][iCB] = val
			Z = np.array(ZList)
			
			temp_max = np.nanmax(Z)
			temp_min = np.nanmin(Z)
			"""
			dcolor = calcDcolor(temp_max, temp_min)
			if isnan(dcolor):
				continue
				
			colorbarArray = np.arange(floor(chi2pDOFmin / dcolor) * dcolor, chi2pDOFmax, dcolor)
			"""
			colorbarArray = calcColorbarArray(temp_max, temp_min, chi2pDOFmax, chi2pDOFmin)
			
			fig = plt.figure()
			ax = fig.add_subplot(1, 1, 1)
			cont = ax.contour(X, Y, Z, colorbarArray)
			
			cont.clabel(fmt = '%1.2f', fontsize = 12)
			ax.set_title('Chi-Squared Value (Inclination Angle ' + IList[iI] + ' degree; Rin from ' + RinList[0] + ' to ' + RinList[-1] + ' au)')
			plt.xlabel('Radius of the Centrifugal Barrier (au)', fontsize = 14)
			plt.ylabel('Central Mass (Msun)', fontsize = 14)
			plt.colorbar(cont)
			
			plt.savefig(saveepsname + 'chi2_CBvsM_allRin_I' + IList[iI] + '.eps')
			plt.close(fig)

		temp_max = np.nanmax(Zmin)
		temp_min = np.nanmin(Zmin)
		"""
		dcolor = calcDcolor(temp_max, temp_min)
		colorbarArray = np.arange(floor(chi2pDOFmin / dcolor) * dcolor, chi2pDOFmax, dcolor)
		"""
		colorbarArray = calcColorbarArray(temp_max, temp_min, chi2pDOFmax, chi2pDOFmin)
		
		fig = plt.figure()
		ax = fig.add_subplot(1, 1, 1)
		cont = ax.contour(X, Y, Zmin, colorbarArray)
		cont.clabel(fmt = '%1.2f', fontsize = 12)
		ax.set_title('Chi-Squared Value (Inclination Angle from ' + IList[0] + ' to ' + IList[-1] + ' degree; Rin from ' + RinList[0] + ' to ' + RinList[-1] + ' au)')
		plt.xlabel('Radius of the Centrifugal Barrier (au)', fontsize = 14)
		plt.ylabel('Central Mass (Msun)', fontsize = 14)
		plt.colorbar(cont)

		plt.savefig(saveepsname + 'chi2_CBvsM_allRin_allI.eps')
		plt.close(fig)
	
	
	### for Rin ###
	
	if nRin > 1:
		if nI > 1:
			X, Y = np.meshgrid(Rinarray, Iarray)

			Zmin = np.full((nI, nRin), np.nan)
			
			for iM in range(nM):
				ZList = []
				for iI in range(nI):
					ZList.append([])
					for iRin in range(nRin):
						val = chi2pDOFfltTable[iM][0][iI][iRin]
						for iCB in range(nCB):
							temp = chi2pDOFfltTable[iM][iCB][iI][iRin]
							if val > temp:
								val = temp
						ZList[iI].append(val)
						if np.isnan(Zmin[iI][iRin]) or Zmin[iI][iRin] > val:
							Zmin[iI][iRin] = val
				
				Z = np.array(ZList)
				
				temp_max = np.nanmax(Z)
				temp_min = np.nanmin(Z)
				"""
				dcolor = calcDcolor(temp_max, temp_min)
				if isnan(dcolor):
					continue

				colorbarArray = np.arange(floor(chi2pDOFmin / dcolor) * dcolor, chi2pDOFmax, dcolor)
				"""
				colorbarArray = calcColorbarArray(temp_max, temp_min, chi2pDOFmax, chi2pDOFmin)
				
				fig = plt.figure()
				ax = fig.add_subplot(1, 1, 1)
				cont = ax.contour(X, Y, Z, colorbarArray)
				
				cont.clabel(fmt = '%1.2f', fontsize = 12)
				ax.set_title('Chi-Squared Value (Central Mass ' + MList[iM] + ' Msun; CB from ' + CBList[0] + ' to ' + CBList[-1] + ' au)')
				plt.xlabel('Inner Radius (au)', fontsize = 14)
				plt.ylabel('Inclination Angle (degree)', fontsize = 14)
				
				plt.colorbar(cont)
				
				plt.savefig(saveepsname + 'chi2_RinvsI_allCB_M' + MList[iM] + '.eps')
				plt.close(fig)
			
			temp_max = np.nanmax(Zmin)
			temp_min = np.nanmin(Zmin)
			"""
			dcolor = calcDcolor(temp_max, temp_min)
			colorbarArray = np.arange(floor(chi2pDOFmin / dcolor) * dcolor, chi2pDOFmax, dcolor)
			"""
			colorbarArray = calcColorbarArray(temp_max, temp_min, chi2pDOFmax, chi2pDOFmin)

			fig = plt.figure()
			ax = fig.add_subplot(1, 1, 1)
			cont = ax.contour(X, Y, Zmin, colorbarArray)
			cont.clabel(fmt = '%1.2f', fontsize = 12)
			ax.set_title('Chi-Squared Value (Central Mass from ' + MList[0] + ' to ' + MList[-1] + ' Msun; CB from ' + CBList[0] + ' to ' + CBList[-1] + ' au)')
			plt.xlabel('Inner Radius (au)', fontsize = 14)
			plt.ylabel('Inclination Angle (degree)', fontsize = 14)
			plt.colorbar(cont)

			plt.savefig(saveepsname + 'chi2_RinvsI_allCB_allM.eps')
			plt.close(fig)
		
		if nM > 1:
			X, Y = np.meshgrid(Rinarray, Marray)

			Zmin = np.full((nM, nRin), np.nan)

			for iI in range(nI):
				ZList = []
				for iM in range(nM):
					ZList.append([])
					for iRin in range(nRin):
						val = chi2pDOFfltTable[iM][0][iI][iRin]
						for iCB in range(nCB):
							temp = chi2pDOFfltTable[iM][iCB][iI][iRin]
							if val > temp:
								val = temp
						ZList[iM].append(val)
						if np.isnan(Zmin[iM][iRin]) or Zmin[iM][iRin] > val:
							Zmin[iM][iRin] = val
				Z = np.array(ZList)
				
				temp_max = np.nanmax(Z)
				temp_min = np.nanmin(Z)
				"""
				dcolor = calcDcolor(temp_max, temp_min)
				if isnan(dcolor):
					continue
					
				colorbarArray = np.arange(floor(chi2pDOFmin / dcolor) * dcolor, chi2pDOFmax, dcolor)
				"""
				colorbarArray = calcColorbarArray(temp_max, temp_min, chi2pDOFmax, chi2pDOFmin)
				
				fig = plt.figure()
				ax = fig.add_subplot(1, 1, 1)
				cont = ax.contour(X, Y, Z, colorbarArray)
				
				cont.clabel(fmt = '%1.2f', fontsize = 12)
				ax.set_title('Chi-Squared Value (Inclination Angle ' + IList[iI] + ' degree; CB from ' + CBList[0] + ' to ' + CBList[-1] + ' au)')
				plt.xlabel('Inner Radius (au)', fontsize = 14)
				plt.ylabel('Central Mass (Msun)', fontsize = 14)
				plt.colorbar(cont)
				
				plt.savefig(saveepsname + 'chi2_RinvsM_allCB_I' + IList[iI] + '.eps')
				plt.close(fig)

			temp_max = np.nanmax(Zmin)
			temp_min = np.nanmin(Zmin)
			"""
			dcolor = calcDcolor(temp_max, temp_min)
			colorbarArray = np.arange(floor(chi2pDOFmin / dcolor) * dcolor, chi2pDOFmax, dcolor)
			"""
			colorbarArray = calcColorbarArray(temp_max, temp_min, chi2pDOFmax, chi2pDOFmin)
			
			fig = plt.figure()
			ax = fig.add_subplot(1, 1, 1)
			cont = ax.contour(X, Y, Zmin, colorbarArray)
			cont.clabel(fmt = '%1.2f', fontsize = 12)
			ax.set_title('Chi-Squared Value (Inclination Angle from ' + IList[0] + ' to ' + IList[-1] + ' degree; CB from ' + CBList[0] + ' to ' + CBList[-1] + ' au)')
			plt.xlabel('Inner Radius (au)', fontsize = 14)
			plt.ylabel('Central Mass (Msun)', fontsize = 14)
			plt.colorbar(cont)

			plt.savefig(saveepsname + 'chi2_RinvsM_allCB_allI.eps')
			plt.close(fig)


