import os
import sys
import glob
import subprocess


obsFitsname = "SO.67-56.autocleaned.robust0.5.niter5000.threshold18mJy.notaper.image.pbcor.regrid.fits"


for i in range(4):
	case = i
	
	if case == 0:
		outfilename = 'Elias29_SO_onlyIRE_varRin.in'
		MList = ['0.2', '0.4', '0.6', '0.8', '1.0', '1.2', '1.4', '1.6', '1.8', '2.0', '3.0', '4.0', '5.0', '6.0', '8.0', '10.0']
		CBList = ['1.', '5.', '10.', '20.', '30.', '40.']
		IList = ['0', '30', '60', '90', '120', '150', '180', '210', '240', '270', '300', '330']
		RotList = ['-1']
		RinList = ['CB']
	elif case == 1:
		outfilename = 'Elias29_SO_onlyKep-Rot1_varRin.in'
		MList = ['0.2', '0.4', '0.6', '0.8', '1.0', '1.2', '1.4', '1.6', '1.8', '2.0', '3.0', '4.0', '5.0', '6.0', '8.0', '10.0']
		CBList = ['100.']
		IList = ['0', '30', '60', '90', '120', '150', '180', '210', '240', '270', '300', '330']
		RotList = ['1']
		RinList = ['1.', '5.', '10.', '20.', '30.', '40.']
	elif case == 2:
		outfilename = 'Elias29_SO_onlyKep-Rot-1_varRin.in'
		MList = ['0.2', '0.4', '0.6', '0.8', '1.0', '1.2', '1.4', '1.6', '1.8', '2.0', '3.0', '4.0', '5.0', '6.0', '8.0', '10.0']
		CBList = ['100.']
		IList = ['0', '30', '60', '90', '120', '150', '180', '210', '240', '270', '300', '330']
		RotList = ['-1']
		RinList = ['1.', '5.', '10.', '20.', '30.', '40.']
	elif case == 3:
		outfilename = 'Elias29_SO_IRE-Kep_varRin.in'
		MList = ['0.2', '0.4', '0.6', '0.8', '1.0', '1.2', '1.4', '1.6', '1.8', '2.0', '3.0', '4.0', '5.0', '6.0', '8.0', '10.0']
		CBList = ['5.', '10.', '20.', '30.', '40.']
		IList = ['0', '30', '60', '90', '120', '150', '180', '210', '240', '270', '300', '330']	#
		RotList = ['-1']
		RinList = ['0']
	else:
		print("Choose '.in' file to be created.")
		sys.exit(1)




	VsysOffset_mps = 0
	rms_jypb = 0.007
	threshold_jypb = 0.021
	Ramin_deg = 0 #274.3
	Ramax_deg = 360 #274.4
	Decmin_deg = -180 #-4.7
	Decmax_deg = 180 #-4.6
	Vmin_mps = -12000 #Vsys+-16km/s
	Vmax_mps = 20000 #Vsys+-16km/s
	weight = 0	#0 to calculate the weight. specified value for a fixed weight.

	Obj = 'Elias29'
	Ra = '16h27m09.4358s'
	Dec = '-24d37m19.286s'
	Vsys = '4.0'
	
	line = 'SO67-56'
	restfreq = '261.8437210'
	
	pixSize = '0.01'
	velRes = '0.2'
	
	
	DList = ['137']
	PAList = ['30']
	RoutList = ['50']
	
	ireHeightList = ['0']
	ireFlareList = ['30']
	ireNProfList = ['-1.5']
	ireTProfList = ['0.0']
	
	kepHeightList = ['0']
	kepFlareList = ['30']
	kepNProfList = ['-1.5']
	kepTProfList = ['0.0']
	
	
	CBNList = ['1e8']
	fracList = ['1e-10']
	CBTList = ['10']
	
	LWList = ['1.0']
	BmajList = ['0.832']
	BminList = ['0.488']
	BpaList = ['-85.833']
	
	
	
	### DO NOT EDIT BELOW ###


	fpout = open(outfilename, 'w')


	def writeParams(D, M, CB, I, PA, Rot, Rout, Rin, ireHeight, ireFlare, ireNProf, ireTProf, kepHeight, kepFlare, kepNProf, kepTProf, CBN, frac, CBT, LW, Bmaj, Bmin, Bpa):
		if Rin == 'CB':
			Rin = CB
		modelFitsname = Obj + "-Vsys" + Vsys + "_Line" + line + "_Pix" + pixSize + "as" + velRes + "kmps_D" + D + "M" + M + "CB" + CB + "I" + I + "PA" + PA + "Rot" + Rot + "Rout" + Rout + "Rin" + Rin + "_IRE-T" + ireHeight + "Flare" + ireFlare + "Nprof" + ireNProf + "Tprof" + ireTProf + "_Kep-T" + kepHeight + "Flare" + kepFlare + "Nprof" + kepNProf + "Tprof" + kepTProf + "_LW" + LW + "_Beam" + Bmaj + "x" + Bmin + "PA" + Bpa
		fpout.write("cubechi2." + modelFitsname + ".out\n")
		fpout.write(obsFitsname + "\n")
		fpout.write(modelFitsname + ".fits\n")
		fpout.write(str(VsysOffset_mps) + "\t# Offset for Systemic Velocity (m/s)\n")
		fpout.write(str(rms_jypb) + " " + str(threshold_jypb) + "\t# rms (Jy/beam), threshold (Jy/beam)\n")
		fpout.write(str(Ramin_deg) + " " + str(Ramax_deg) + " " + str(Decmin_deg) + " " + str(Decmax_deg) + " " + str(Vmin_mps) + " " + str(Vmax_mps) + "\t# Ramin, Ramax, Decmin, Decmax, Vmin, Vmax\n")
		fpout.write(str(weight) + "\t# Weight // If this value is <=0, weight is calculated.\n")
		fpout.write("\n")





	nCube = 1
	nPV = 1
	iCube = 0
	iPV = 0
	nCubeidList = [len(DList), len(MList), len(CBList), len(IList), len(PAList), len(RotList), len(RoutList), len(RinList), len(ireHeightList), len(ireFlareList), len(ireNProfList), len(ireTProfList), len(kepHeightList), len(kepFlareList), len(kepNProfList), len(kepTProfList), len(CBNList), len(fracList), len(CBTList), len(LWList), len(BmajList), len(BminList), len(BpaList)]

	for item in nCubeidList:
		nCube = nCube * item

	Cubeid = [0 for i in range(len(nCubeidList))]

	for iCube in range(nCube):
		writeParams(DList[Cubeid[0]], MList[Cubeid[1]], CBList[Cubeid[2]], IList[Cubeid[3]], PAList[Cubeid[4]], RotList[Cubeid[5]], RoutList[Cubeid[6]], RinList[Cubeid[7]], ireHeightList[Cubeid[8]], ireFlareList[Cubeid[9]], ireNProfList[Cubeid[10]], ireTProfList[Cubeid[11]], kepHeightList[Cubeid[12]], kepFlareList[Cubeid[13]], kepNProfList[Cubeid[14]], kepTProfList[Cubeid[15]], CBNList[Cubeid[16]], fracList[Cubeid[17]], CBTList[Cubeid[18]], LWList[Cubeid[19]], BmajList[Cubeid[20]], BminList[Cubeid[21]], BpaList[Cubeid[22]])
		
		Cubeid[-1] = Cubeid[-1] + 1
		for i in range(1, len(nCubeidList), 1):
			if Cubeid[-i] == nCubeidList[-i]:
				Cubeid[-i] = 0
				Cubeid[-(i + 1)] = Cubeid[-(i + 1)] + 1



	fpout.close()

