import os
import sys
import numpy as np

saveDir = 'chi2plot_varRin-limI60-120-ev10/'

beamsize = 0.832 * 0.488 / 0.01**2 # FWHM_major (as) * FWHM_minor (as) / pixsize^2
velwidth = 1.0 / 0.2 # FWHM (km/s) / pixsize
num_indepPix = beamsize / 4. * velwidth / 2.

f_prodIndepPix = False	#True

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
	saveTableFilename = saveDir + infilename + '.table.tex'
	
	
	#############################
	
	axisID_M = 0
	axisID_CB = 1
	axisID_I = 2
	axisID_Rin = 3


	nM = len(MList)
	nCB = len(CBList)
	nI = len(IList)
	nRin = len(RinList)

	fin = open(infilename, 'r')
	#fout = open(saveTableFilename, 'w')

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
		saveTableFilename = saveDir + infilename + '.indepPixconsidered.table.tex'
	
	fout = open(saveTableFilename, 'w')
	
	###### Table ######
		
	fout.write("\\documentclass[manuscript]{aastex}\n\n")
	fout.write("\\newcommand{\Msun}{$M_\\odot$}\n\n")
	
	fout.write("\\begin{document}\n\n")
	
	
	### min value ###
	min = 1e100
	minID = [-1, -1, -1, -1]
	for iM in range(nM):
		for iCB in range(nCB):
			for iI in range(nI):
				for iRin in range(nRin):
					val = chi2pDOFfltTable[iM][iCB][iI][iRin]
					if min > val:
						min = val
						minID[0] = iM
						minID[1] = iCB
						minID[2] = iI
						minID[3] = iRin
	
	fout.write("min value = {} at {}\n\n".format(min, minID))
	fout.write("\\quad M = {}\n\n".format(MList[minID[0]]))
	fout.write("\\quad CB = {}\n\n".format(CBList[minID[1]]))
	fout.write("\\quad I = {}\n\n".format(IList[minID[2]]))
	fout.write("\\quad Rin = {}\n\n\n\n".format(RinList[minID[3]]))

	fout.write("Parameters:\n\n")
	fout.write("\\quad M from {} to {} Msun\n\n".format(MList[0], MList[-1]))
	fout.write("\\quad CB from {} to {} au\n\n".format(CBList[0], CBList[-1]))
	fout.write("\\quad I from {} to {} degree\n\n".format(IList[0], IList[-1]))
	fout.write("\\quad Rin from {} to {} au\n\n".format(RinList[0], RinList[-1]))
	fout.write("\clearpage\n\n")
	
	### for M ###
	
	temp_chi2pDOFfltTable = np.nanmin(chi2pDOFfltTable, axis = axisID_Rin)
	for iM in range(nM):
		fout.write("\\begin{table}\n\\begin{center}\n")
		
		strCaption = "Reduced $\\chi^2$ Test on the *** Model with the Central Mass of %s \\Msun\\ with the inner radius from %s to %s au" % (MList[iM], RinList[0], RinList[-1])
		fout.write("\\caption{%s\n\t\t\label{tb:%s_iM%d}}\n" % (strCaption, infilename, iM))
		
		strTemp = "c"
		for iCB in range(nCB):
			strTemp = strTemp + "c"
		fout.write("\\begin{tabular}{%s}\n" % strTemp)
		
		fout.write("\\hline \n& \\multicolumn{%d}{c}{Radius of the Centrifugal Barrier (au)} \\\\ \n" % nCB)
		strTemp = "Inclination Angle (\\degr)"
		for iCB in range(nCB):
			strTemp = strTemp + (" & %s" % CBList[iCB])
		fout.write(strTemp + " \\\\ \\hline \\hline \n")
		
		for iiI in range(nI):
			iI = nI - 1 - iiI
			chi2pDOFstr = str(IList[iI])
			for iCB in range(nCB):
				chi2pDOFstr = chi2pDOFstr + (" & %.3lf" % temp_chi2pDOFfltTable[iM][iCB][iI])
			fout.write(chi2pDOFstr + " \\\\ \n")
		
		fout.write("\\hline\n\\end{tabular}\n\\end{center}\n\\end{table}\n\n\\clearpage\n\n")
			

	### for CB ###
	
	for iCB in range(nCB):
		fout.write("\\begin{table}\n\\begin{center}\n")
		
		strCaption = "Reduced $\\chi^2$ Test on the *** Model with the radius of the centrifugal barrier of %s au with the inner radius from %s to %s au" % (CBList[iCB], RinList[0], RinList[-1])
		fout.write("\\caption{%s\n\t\t\label{tb:%s_iCB%d}}\n" % (strCaption, infilename, iCB))
		
		strTemp = "c"
		for iM in range(nM):
			strTemp = strTemp + "c"
		fout.write("\\begin{tabular}{%s}\n" % strTemp)
		
		fout.write("\\hline \n& \\multicolumn{%d}{c}{Central Mass (\\Msun)} \\\\ \n" % nM)
		strTemp = "Inclination Angle (\\degr)"
		for iM in range(nM):
			strTemp = strTemp + (" & %s" % MList[iM])
		fout.write(strTemp + " \\\\ \\hline \\hline \n")
		
		for iiI in range(nI):
			iI = nI - 1 - iiI
			chi2pDOFstr = str(IList[iI])
			for iM in range(nM):
				chi2pDOFstr = chi2pDOFstr + (" & %.3lf" % temp_chi2pDOFfltTable[iM][iCB][iI])
			fout.write(chi2pDOFstr + " \\\\ \n")
		
		fout.write("\\hline\n\\end{tabular}\n\\end{center}\n\\end{table}\n\n\\clearpage\n\n")


	### for I ###

	for iI in range(nI):
		fout.write("\\begin{table}\n\\begin{center}\n")
		
		strCaption = "Reduced $\\chi^2$ Test on the *** Model with the Inclination Angle of %s\\degr\\ with the inner radius from %s to %s au" % (IList[iI], RinList[0], RinList[-1])
		fout.write("\\caption{%s\n\t\t\label{tb:%s_iI%d}}\n" % (strCaption, infilename, iI))
		
		strTemp = "c"
		for iCB in range(nCB):
			strTemp = strTemp + "c"
		fout.write("\\begin{tabular}{%s}\n" % strTemp)
		
		fout.write("\\hline \n& \\multicolumn{%d}{c}{Radius of the Centrifugal Barrier (au)} \\\\ \n" % nCB)
		strTemp = "Central Mass (\\Msun)"
		for iCB in range(nCB):
			strTemp = strTemp + (" & %s" % CBList[iCB])
		fout.write(strTemp + " \\\\ \\hline \\hline \n")
		
		for iiM in range(nM):
			iM = nM - 1 - iiM
			chi2pDOFstr = str(MList[iM])
			for iCB in range(nCB):
				chi2pDOFstr = chi2pDOFstr + (" & %.3lf" % temp_chi2pDOFfltTable[iM][iCB][iI])
			fout.write(chi2pDOFstr + " \\\\ \n")
		
		fout.write("\\hline\n\\end{tabular}\n\\end{center}\n\\end{table}\n\n\\clearpage\n\n")
	
	
	### for Rin ###
	
	temp_chi2pDOFfltTable = np.amin(chi2pDOFfltTable, axis = axisID_CB)
	
	for iRin in range(nRin):
		fout.write("\\begin{table}\n\\begin{center}\n")
		
		strCaption = "Reduced $\\chi^2$ Test on the *** Model with the inner radius of %s au with the radius of the centrifugal barrier from %s to %s au" % (RinList[iRin], CBList[0], CBList[-1])
		fout.write("\\caption{%s\n\t\t\label{tb:%s_iRin%d}}\n" % (strCaption, infilename, iRin))
		
		strTemp = "c"
		for iM in range(nM):
			strTemp = strTemp + "c"
		fout.write("\\begin{tabular}{%s}\n" % strTemp)
		
		fout.write("\\hline \n& \\multicolumn{%d}{c}{Central Mass (\\Msun)} \\\\ \n" % nM)
		strTemp = "Inclination Angle (\\degr)"
		for iM in range(nM):
			strTemp = strTemp + (" & %s" % MList[iM])
		fout.write(strTemp + " \\\\ \\hline \\hline \n")
		
		for iiI in range(nI):
			iI = nI - 1 - iiI
			chi2pDOFstr = str(IList[iI])
			for iM in range(nM):
				chi2pDOFstr = chi2pDOFstr + (" & %.3lf" % temp_chi2pDOFfltTable[iM][iI][iRin])
			fout.write(chi2pDOFstr + " \\\\ \n")
		
		fout.write("\\hline\n\\end{tabular}\n\\end{center}\n\\end{table}\n\n\\clearpage\n\n")
	
	
	fout.write("\n\n\\end{document}\n\n")

