#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "fitsio.h"

using namespace std;

#define EPS 1e-308
#define INF 1e+308

#define naxis 3
#define nCubes 2
#define delim " "
#define noval 0.0
#define calcRatioLim 100
#define lenfilename 512
#define lenline 512
#define DataDir "fits/"
#define OutDir "out/"

const bool f_skipNegative = true;
const bool logger = true;
const bool saveRegrid = false;
const double db_erratio = 5e-15;


FILE *fp_out;

void printerror(int status) {
	if (status) {
		fits_report_error(stderr, status);
		exit(status);
	}
	return;
}

//////////
class Cubedata {
public:
	long npix[naxis], npixels;
	char filename[lenfilename];
	fitsfile *fptr;
	double *valData;
	double valMax, valMin;
	double crpix[naxis], crval[naxis], cdelt[naxis];
	
	
	Cubedata();
	//Cubedata(const Cubedata &);
	~Cubedata();
	
	void dataInit();
	void pixInit(long *npix_in, double *crpix_in, double *crval_in, double *cdelt_in);
	void readFits(char* filename_in);
	
	int IP2Pix(double *ip);
	void Pix2IP(int ipix, double *ip);
	void setMaxMin();
	void getMaxMin(double *max, double *min);
	
	bool skipCell(double *ip, bool f_threshold, double threshold4skip);
	void setVal(double *ip, double val);
	void getVal(double *ip, double *val);
	
	void getP(double *ip, double *p);
	void getIP(double *ip, double *p);
};


Cubedata::Cubedata() {
}

Cubedata::~Cubedata() {
	free(valData);
	valData = NULL;
}

void Cubedata::dataInit() {
	
	valData = (double*) malloc(sizeof(double) * npixels);
	
	for (int i = 0; i < npixels; ++i) valData[i] = 0.;
	
	return ;
}

void Cubedata::pixInit(long *npix_in, double *crpix_in, double *crval_in, double *cdelt_in) {
	npixels = 1;
	for (int i = 0; i < naxis; ++i) {
		npix[i] = npix_in[i];
		crpix[i] = crpix_in[i];
		crval[i] = crval_in[i];
		cdelt[i] = cdelt_in[i];
		npixels *= npix[i];
	}
	
	dataInit();
	return ;
}

void Cubedata::readFits(char* filename_in) {
	sscanf(filename_in, "%s", filename);
	
	int status = 0, nfound, dumint, firstelem = 1;
	double nullval = 0;
	
	
		
	if (fits_open_file(&fptr, filename, READONLY, &status)) printerror(status);
	
	if (fits_read_keys_lng(fptr, "NAXIS", 1, naxis, npix, &nfound, &status)) printerror(status);
	
	if (fits_read_keys_dbl(fptr, "CRPIX", 1, naxis, crpix, &dumint, &status)) printerror(status);
	if (fits_read_keys_dbl(fptr, "CRVAL", 1, naxis, crval, &dumint, &status)) printerror(status);
	if (fits_read_keys_dbl(fptr, "CDELT", 1, naxis, cdelt, &dumint, &status)) printerror(status);
	
	npixels = 1;
	for (int i = 0; i < naxis; ++i) {
		npixels *= npix[i];
		crpix[i] -= 1; // 1-origin -> 0-origin
		
		/*
		if (i == naxis - 1) {
			crval[i] *= 1.0e-3;		// m/s -> km/s
			cdelt[i] *= 1.0e-3;
		} else {
			crval[i] *= 60.*60.;	//radian -> arcsec
			cdelt[i] *= 60.*60.;
		}
		//*/
	}
	
	

	dataInit();
	
	if (fits_read_img(fptr, TDOUBLE, firstelem, npixels, &nullval, valData, &dumint, &status)) printerror(status);
	
	if (fits_close_file(fptr, &status)) printerror(status);
	return ;
}

int Cubedata::IP2Pix(double *ip) {
	int ans;
	int ip_dum[naxis];
	
	for (int iaxis = 0; iaxis < naxis; ++iaxis) ip_dum[iaxis] = floor(ip[iaxis] + 0.5);
	ans = ip_dum[naxis - 1];
	
	for (int iaxis = naxis - 2; iaxis >= 0; --iaxis) ans = ans * npix[iaxis] + ip_dum[iaxis];
	
	return ans;
}

void Cubedata::Pix2IP(int ipix, double *ip) {
	for (int iaxis = 0; iaxis < naxis; ++iaxis) {
		ip[iaxis] = double(ipix % npix[iaxis]);
		ipix = ipix / npix[iaxis];
	}
	return ;
}

void Cubedata::setMaxMin() {
	valMax = -INF;
	valMin = INF;
	double val;
	
	for (int i = 0; i < npixels; ++i) {
		val = valData[i];
		if (valMax < val) valMax = val;
		if (valMin > val) valMin = val;
	}
	
	return ;
}

void Cubedata::getMaxMin(double *max, double *min) {
	setMaxMin();
	if (max) *max = valMax;
	if (min) *min = valMin;
	return ;
}


bool Cubedata::skipCell(double *ip, bool f_threshold = false, double threshold4skip = 0.) {
	double val = noval;
	for (int i = 0; i < naxis; ++i) {
		if (ip[i] < 0 || ip[i] >= npix[i]) return true;
	}
	
	if (f_threshold) {
		getVal(ip, &val);
		if (abs(val) < threshold4skip) return true;
		if (f_skipNegative && val < threshold4skip) return true;
	}
	return false;
}


void Cubedata::setVal(double *ip, double val) {
	if (skipCell(ip)) {
		valData[IP2Pix(ip)] = noval;
	} else {
		valData[IP2Pix(ip)] = val;
	}
	return ;
}

void Cubedata::getVal(double *ip, double *val) {
	if (skipCell(ip)) {
		if (val) *val = noval;
	} else {
		if (val) *val = valData[IP2Pix(ip)];
	}
	return ;
}

void Cubedata::getP(double *ip, double *p) {
	for (int i = 0; i < naxis; ++i) p[i] = crval[i] + (ip[i] - crpix[i]) * cdelt[i];
	return ;
}

void Cubedata::getIP(double *ip, double *p) {
	for (int i = 0; i < naxis; ++i) ip[i] = crpix[i] + (p[i] - crval[i]) / cdelt[i];
	return ;
}



//////////
class Cubepair {
public:
	Cubedata *Cube1, *Cube2, *Cube2regrid;
	char filename1[lenfilename], filename2[lenfilename];

	double VsysOffset, weight4Cube2regrid, rms_Cube1, threshold_Cube1;
	
	double limP[2][naxis], rP[naxis][2];
	int rIP[naxis][2], nregion;
	
	
	Cubepair();
	~Cubepair();
	
	void init(Cubedata *Cube1_in, Cubedata *Cube2_in, Cubedata *Cube2_regrid, char* filename1_in, char* filename2_in, double weight_in, double rms_Cube1_in, double threshold_Cube1_in, double VsysOffset_in);
	bool skipCell(double *ip);
	bool skipCell_andThreshold(double *ip, double threshold4skip_Cube1_in, double threshold4skip_Cube2regrid_in);
	bool skipCell_orThreshold(double *ip, double threshold4skip_Cube1_in, double threshold4skip_Cube2regrid_in);
	bool skipCell_region4comparing(double *ip, bool f_threshold, double threshold4skip_Cube1_in, double threshold4skip_Cube2regrid_in);
	void regrid();
	void setRegion(double *limP_in);
	
	void calcRatio();
	void setRatio();
	void setWeight();
	void compare();
};

Cubepair::Cubepair() {
}

Cubepair::~Cubepair() {
}

 
void Cubepair::init(Cubedata *Cube1_in, Cubedata *Cube2_in, Cubedata *Cube2regrid_in, char* filename1_in, char* filename2_in, double weight_in, double rms_Cube1_in, double threshold_Cube1_in, double VsysOffset_in = 0.) {
	Cube1 = Cube1_in;
	Cube2 = Cube2_in;
	Cube2regrid = Cube2regrid_in;
	
	sscanf(filename1_in, "%s", filename1);
	sscanf(filename2_in, "%s", filename2);
	
	weight4Cube2regrid = weight_in;
	rms_Cube1 = rms_Cube1_in;
	threshold_Cube1 = threshold_Cube1_in;
	VsysOffset = VsysOffset_in;
	
	
	Cube1->readFits(filename1);
	Cube2->readFits(filename2);
	regrid();
}

bool Cubepair::skipCell(double *ip) {
	return Cube1->skipCell(ip) || Cube2regrid->skipCell(ip);
}

bool Cubepair::skipCell_andThreshold(double *ip, double threshold4skip_Cube1_in, double threshold4skip_Cube2regrid_in) {
	return skipCell(ip) || (Cube1->skipCell(ip, true, threshold4skip_Cube1_in) && Cube2regrid->skipCell(ip, true, threshold4skip_Cube2regrid_in));
}

bool Cubepair::skipCell_orThreshold(double *ip, double threshold4skip_Cube1_in, double threshold4skip_Cube2regrid_in) {
	return skipCell(ip) || (Cube1->skipCell(ip, true, threshold4skip_Cube1_in) || Cube2regrid->skipCell(ip, true, threshold4skip_Cube2regrid_in));
}

bool Cubepair::skipCell_region4comparing(double *ip, bool f_threshold = false, double threshold4skip_Cube1_in = 0., double threshold4skip_Cube2regrid_in = 0.) {
	for (int i = 0; i < naxis; ++i) {
		if (ip[i] < rIP[i][0] || ip[i] >= rIP[i][1]) return true;
	}
	
	if (f_threshold) {
		return skipCell_orThreshold(ip, threshold4skip_Cube1_in, threshold4skip_Cube2regrid_in);
	} else {
		return skipCell(ip);
	}
}


void Cubepair::regrid() {
	Cube2regrid->pixInit(Cube1->npix, Cube1->crpix, Cube1->crval, Cube1->cdelt);
	
	double p_regrid[naxis], ip[naxis], ip_regrid[naxis], ipLU[naxis][2], dip[naxis], ip_dum[naxis];
	double val, valLU[2 << naxis];
	bool f_continue = false;
	
	for (int ipix = 0; ipix < Cube2regrid->npixels; ++ipix) {
		f_continue = false;
		Cube2regrid->Pix2IP(ipix, ip_regrid);
		Cube2regrid->getP(ip_regrid, p_regrid);
		p_regrid[naxis - 1] -= VsysOffset;
		
		Cube2->getIP(ip, p_regrid);
		for (int iaxis = 0; iaxis < naxis; ++ iaxis) {
			ipLU[iaxis][0] = floor(ip[iaxis] * (1. + db_erratio));
			ipLU[iaxis][1] = ceil(ip[iaxis] * (1. - db_erratio));
			dip[iaxis] = ip[iaxis] - ipLU[iaxis][0];
			
			if (ipLU[iaxis][0] < 0 || ipLU[iaxis][1] >= Cube2->npix[iaxis]) f_continue = true;
		}
		if (f_continue) continue;
		
		
		val = 0.;
		for (int iLU = 0; iLU < (2 << naxis); ++iLU) {
			for (int iaxis = 0; iaxis < naxis; ++iaxis) {
				ip_dum[iaxis] = ipLU[iaxis][(iLU / (2 << (naxis - 1 - iaxis))) % 2];
			}
			Cube2->getVal(ip_dum, &valLU[iLU]);
			for (int iaxis = 0; iaxis < naxis; ++iaxis) {
				if ((iLU / (2 << (naxis - 1 - iaxis))) % 2) {
					valLU[iLU] *= dip[iaxis];
				} else {
					valLU[iLU] *= (1. - dip[iaxis]);
				}
			}
			val += valLU[iLU];
		}
		
		Cube2regrid->setVal(ip_regrid, val);
	}
	
	Cube2regrid->setMaxMin();
	return ;
}


void Cubepair::setRegion(double *limP_in) {
	double ip_dum[naxis], p_dum[naxis], p_dum2[naxis], dum;
	char str_region[lenline], str_dum[lenline];
	

	for (int iaxis = 0; iaxis < naxis; ++iaxis) for (int i = 0; i < 2; ++i) limP[i][iaxis] = limP_in[iaxis * 2 + i];
	
	for (int i = 0; i < 2; ++i) {
		for (int iaxis = 0; iaxis < naxis; ++iaxis) p_dum[iaxis] = limP[i][iaxis];
		Cube1->getIP(ip_dum, p_dum);

		for (int iaxis = 0; iaxis < naxis; ++iaxis) {
			rIP[iaxis][i] = ip_dum[iaxis];
			if (rIP[iaxis][i] < 0) rIP[iaxis][i] = 0;
			if (rIP[iaxis][i] >= Cube1->npix[iaxis]) rIP[iaxis][i] = Cube1->npix[iaxis] - 1;
		}
		
		for (int iaxis = 0; iaxis < naxis; ++iaxis) ip_dum[iaxis] = rIP[iaxis][i];
		Cube1->getP(ip_dum, p_dum);
		for (int iaxis = 0; iaxis < naxis; ++iaxis) rP[iaxis][i] = p_dum[iaxis];
	}
	for (int iaxis = 0; iaxis < naxis; ++iaxis) {
		if (rIP[iaxis][0] > rIP[iaxis][1]) {
			dum = rIP[iaxis][0];
			rIP[iaxis][0] = rIP[iaxis][1];
			rIP[iaxis][1] = dum;
			dum = rP[iaxis][0];
			rP[iaxis][0] = rP[iaxis][1];
			rP[iaxis][1] = dum;
		}
	}
	

	nregion = 1;
	for (int iaxis = 0; iaxis < naxis; ++iaxis) nregion *= (rIP[iaxis][1] - rIP[iaxis][0] + 1);
	
	
	for (int iaxis = 0; iaxis < naxis; ++iaxis) ip_dum[iaxis] = 0;
	Cube1->getP(ip_dum, p_dum);
	for (int iaxis = 0; iaxis < naxis; ++iaxis) ip_dum[iaxis] = Cube1->npix[iaxis] - 1;
	Cube1->getP(ip_dum, p_dum2);
	
	sprintf(str_region, "");
	for (int iaxis = 0; iaxis < naxis; ++iaxis) {
		sprintf(str_dum, "\taxis %d [0, %ld] (pix) <-> [%.4f, %.4f]", iaxis, Cube1->npix[iaxis], p_dum[iaxis], p_dum2[iaxis]);
		if (iaxis == naxis - 1) {
			strcat(str_dum, " (m/s), ");
		} else {
			strcat(str_dum, " (deg), ");
		}
		strcat(str_region, str_dum);
	}
	if (saveRegrid) fprintf(fp_out, "!\tWhole region in Cube 1         : %s\n", str_region);
	
	sprintf(str_region, "");
	for (int iaxis = 0; iaxis < naxis; ++iaxis) {
		sprintf(str_dum, "\taxis %d [%d, %d] (pix) <-> [%.4f, %.4f]", iaxis, rIP[iaxis][0], rIP[iaxis][1], rP[iaxis][0], rP[iaxis][1]);
		if (iaxis == naxis - 1) {
			strcat(str_dum, " (m/s), ");
		} else {
			strcat(str_dum, " (deg), ");
		}
		strcat(str_region, str_dum);
	}
	if (saveRegrid) fprintf(fp_out, "!\tRegion for Comparison in Cube 1: %s\n", str_region);
	if (saveRegrid) fprintf(fp_out, "!\t#Cell = %d\n", nregion);
	
	return ;
}


//////////

void Cubepair::calcRatio() {
	double val_Cube1, val_Cube2regrid;
	double ip[naxis];
	int ndatapoints4weight = 0;
	bool f_continue = false;
	
	double sxx = 0., sxy = 0.;
	
	
	for (int ipix = 0; ipix < Cube1->npixels; ++ipix) {
		f_continue = false;
		Cube1->Pix2IP(ipix, ip);
		
		for (int iaxis = 0; iaxis < naxis; ++iaxis) {
			if (ip[iaxis] < rIP[iaxis][0] || ip[iaxis] >= rIP[iaxis][1]) f_continue = true;
		}
		if (f_continue) continue;
		if (skipCell_orThreshold(ip, threshold_Cube1, threshold_Cube1 / weight4Cube2regrid)) continue;
		
		++ndatapoints4weight;
		
		Cube1->getVal(ip, &val_Cube1);
		Cube2regrid->getVal(ip, &val_Cube2regrid);
		sxx += val_Cube1 * val_Cube1;
		sxy += val_Cube1 * val_Cube2regrid;
	}
	
	
	if (ndatapoints4weight == 0) {
		printf("\n\tERROR :: There is no data point for calculating the normalization constant.  Check the threshold level. \n\n");
		//exit(1);
	}
	
	weight4Cube2regrid = sxx / sxy;
	
	if (abs(weight4Cube2regrid) < EPS) {
		printf("\n\tWARN :: The normalization constant is 0. \n");
		printf("\t\t\tweight = %le\n\n", weight4Cube2regrid);
		//exit(1);
	}
	
	if (weight4Cube2regrid < 0.) {
		printf("\n\tWARN :: The normalization constant is negative. \n");
		printf("\t\t\tweight = %le\n", weight4Cube2regrid);
	}
	
	return ;
}



void Cubepair::setRatio() {
	double maxval_Cube1, maxval_Cube2regrid, dum_dbl;

	if (abs(weight4Cube2regrid) > 0.) return;
	
	
	
	Cube1->getMaxMin(&maxval_Cube1, &dum_dbl);
	Cube2->getMaxMin(&maxval_Cube2regrid, &dum_dbl);
	
	weight4Cube2regrid = maxval_Cube1 / maxval_Cube2regrid;

	dum_dbl = weight4Cube2regrid;
	for (int i = 0; i < calcRatioLim; ++i) {
		calcRatio();
		if (abs(weight4Cube2regrid - dum_dbl) < EPS) break;
		dum_dbl = weight4Cube2regrid;
	}
}



void Cubepair::setWeight() {
	double ip[naxis];
	double val;
	
	for (int ipix = 0; ipix < Cube2regrid->npixels; ++ipix) {
		Cube2regrid->Pix2IP(ipix, ip);
		Cube2regrid->getVal(ip, &val);
		Cube2regrid->setVal(ip, val * weight4Cube2regrid);
	}
	
	Cube2regrid->setMaxMin();
	return ;
}




void Cubepair::compare() {
	int ncount_region = 0, ncount_threshold = 0, f_used = 0;
	double ip[naxis], p[naxis];
	double val1, val2regrid, sumres2_region = 0., sumres2_threshold = 0.;
	
	printf("\nNormalization constant for regridded Cube data = %.4e\n", weight4Cube2regrid);
	if (saveRegrid) fprintf(fp_out, "!Normalization constant for regridded Cube data = %.4e\n", weight4Cube2regrid);
	if (saveRegrid) fprintf(fp_out, "!\n");
	
    if (saveRegrid) fprintf(fp_out, "!  \t  ");
	if (saveRegrid) fprintf(fp_out, "\n!");
	for (int iaxis = 0; iaxis < naxis; ++iaxis) if (saveRegrid) fprintf(fp_out, "|axis[%d] pix\t\tval\t", iaxis);
	
	if (saveRegrid) fprintf(fp_out, "|Cube1\t\tCube2regrid\t\tused? (0 for unused, 1 for 'in Region' but 'less than threshold', 2 for used)\n");
	
	for (int ipix = 0; ipix < Cube1->npixels; ++ipix) {
		Cube1->Pix2IP(ipix, ip);
		Cube1->getP(ip, p);
		for (int iaxis = 0; iaxis < naxis; ++iaxis) if (saveRegrid) fprintf(fp_out, "\t%lf\t%.8lf", ip[iaxis], p[iaxis]);
		
		if (skipCell(ip)) {
			val1 = noval;
			val2regrid = noval;
		} else {
			Cube1->getVal(ip, &val1);
			Cube2regrid->getVal(ip, &val2regrid);
		}
		
		if (saveRegrid) fprintf(fp_out, "\t%.4e\t%.4e", val1, val2regrid);
		
		if (skipCell(ip) || skipCell_region4comparing(ip)) {
			f_used = 0;
		} else if (skipCell_region4comparing(ip, true, threshold_Cube1, threshold_Cube1)) {
			f_used = 1;
			++ncount_region;
			sumres2_region += (val1 - val2regrid) * (val1 - val2regrid);
		} else {
			f_used = 2;
			++ncount_threshold;
			sumres2_threshold += (val1 - val2regrid) * (val1 - val2regrid);
		}
		
		if (saveRegrid) fprintf(fp_out, "\t%d\n", f_used);
	}
	
	ncount_region += ncount_threshold;
	sumres2_region += sumres2_threshold;
	
	sumres2_region /= rms_Cube1 * rms_Cube1;
	sumres2_threshold /= rms_Cube1 * rms_Cube1;
	
	
	if (saveRegrid) fprintf(fp_out, "!In the specified region      #Cell = %d\n", ncount_region);
	if (saveRegrid) fprintf(fp_out, "!  Sum of (Cube1 - Cube2)^2 / rms^2 = %15.4e\n", sumres2_region);
	if (saveRegrid) fprintf(fp_out, "!                              /DOF = %15.4e\n", sumres2_region / ncount_region);
	if (saveRegrid) fprintf(fp_out, "!\n");
	if (saveRegrid) fprintf(fp_out, "!Larger than the threshold:   #Cell = %d\n", ncount_threshold);
	if (saveRegrid) fprintf(fp_out, "!  Sum of (Cube1 - Cube2)^2 / rms^2 = %15.4e\n", sumres2_threshold);
	if (saveRegrid) fprintf(fp_out, "!                              /DOF = %15.4e\n", sumres2_threshold / ncount_threshold);
	if (saveRegrid) fprintf(fp_out, "\n\n");
	
	if (logger) {
		puts("");
		printf("In the specified region      #Cell = %d\n", ncount_region);
		printf("  Sum of (Cube1 - Cube2)^2 / rms^2 = %15.4e\n", sumres2_region);
		printf("                              /DOF = %15.4e\n", sumres2_region / ncount_region);
		printf("\n");
		printf("Larger than the threshold:   #Cell = %d\n", ncount_threshold);
		printf("  Sum of (Cube1 - Cube2)^2 / rms^2 = %15.4e\n", sumres2_threshold);
		printf("                              /DOF = %15.4e\n", sumres2_threshold / ncount_threshold);
		puts("");
	}
	
	return ;
	
}





//////////


int main() {
	
	char tempfilename[lenfilename];
	
	while (~scanf("%s ", tempfilename)) {
		puts("---------------------------------------------\n");
		char paramline[lenline], *paramline_delim;
		char outfilename[lenfilename], Cubefilename[nCubes][lenfilename];
		double Vsys, rms_Cube1, threshold, weight;
		double limP[naxis * 2];
		
        sscanf(OutDir, "%s", outfilename);
        strcat(outfilename, tempfilename);
		
		for (int i = 0; i < nCubes; ++i) {
			fgets(tempfilename, lenline, stdin);
			sscanf(tempfilename, "%s", tempfilename);
			sscanf(DataDir, "%s", Cubefilename[i]);
			strcat(Cubefilename[i], tempfilename);
		}
		
        
		fgets(paramline, lenline, stdin);
		sscanf(paramline, "%lf ", &Vsys);
		
        fgets(paramline, lenline, stdin);
        sscanf(paramline, "%lf %lf ", &rms_Cube1, &threshold);
		
        fgets(paramline, lenline, stdin);
		
		
		paramline_delim = strtok(paramline, delim);
		for (int i = 0; i < naxis * 2; ++i) {
			sscanf(paramline_delim, "%lf", &limP[i]);
			paramline_delim = strtok(NULL, delim);
		}
		
		
        fgets(paramline, lenline, stdin);
		sscanf(paramline, "%lf ", &weight);
		
		fp_out = fopen(outfilename, "w");
		
        if (saveRegrid) fprintf(fp_out, "!Systemic Velocity: %.4f km/s\n", Vsys);
        if (saveRegrid) fprintf(fp_out, "!rms in Cube1: %.4f Jy/beam\n", rms_Cube1);
		if (saveRegrid) fprintf(fp_out, "!Region for comparison: threshold = %lf Jy/beam\n", threshold);
		for (int iaxis = 0; iaxis < naxis; ++iaxis) if (saveRegrid) fprintf(fp_out, "!                   range for axis %d = [%lf, %lf]\n", iaxis, limP[iaxis * 2], limP[iaxis * 2 + 1]);
		
		
		Cubedata CubedataList[nCubes + 1];
		Cubepair Cubepair_in;
		
		
		Cubepair_in.init(&CubedataList[0], &CubedataList[1], &CubedataList[2], Cubefilename[0], Cubefilename[1], weight, rms_Cube1, threshold, Vsys);

		Cubepair_in.setRegion(limP);
		Cubepair_in.setRatio();
		Cubepair_in.setWeight();
		
		if (saveRegrid) fprintf(fp_out, "!Weight for Cube2: %.4e\n", Cubepair_in.weight4Cube2regrid);
		
		if (saveRegrid) fprintf(fp_out, "!----------\n");
		for (int i = 0; i < nCubes; ++i) {
			if (saveRegrid) fprintf(fp_out, "!Cube %d: %s\n", i, Cubefilename[i]);
			if (logger) printf("Cube %d: %s\n", i, Cubefilename[i]);
			
			
			for (int iaxis = 0; iaxis < naxis; ++iaxis) {
				if (saveRegrid) fprintf(fp_out, "!\t\t\taxis[%d]: crpix = %.4e, crval = %.4e, cdelt = %.4e\t", iaxis, Cubepair_in.Cube1->crpix[iaxis], Cubepair_in.Cube1->crval[iaxis], Cubepair_in.Cube1->cdelt[iaxis]);
				if (iaxis == naxis - 1) {
					if (saveRegrid) fprintf(fp_out, "(unit: m/s)\n");
				} else {
					if (saveRegrid) fprintf(fp_out, "(unit: deg)\n");
				}
			}
		}
		if (saveRegrid) fprintf(fp_out, "!----------\n");
		
		
		Cubepair_in.compare();
        
        fclose(fp_out);
        if (saveRegrid) printf("\nSaved: %s\n\n", outfilename);
		
    }
	
    
    return 0;
}


