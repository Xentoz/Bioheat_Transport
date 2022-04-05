//---------------------------------------------------------------------------------------------
// Program Name: 2D2M_CN_v2.0.2.cpp
// 
// Purpose:      Uses the crank-nicollson method to solve the heat equation in two dimensions dispersing
//				 through multiple mediums and with a heat source. Creates a .dat file which can be 
//				 inputted into gnuplot for a graphical representation. 
// 
// Input data:   initial_conditions_cn.txt
// Output data:	 2D2M_CN_v2.0.2.dat
//---------------------------------------------------------------------------------------------
// Original Author:           Callum Mann
// Date of initial Creation:  13Mar2022
// Previous Version:          2.0.1
// Current Version:           2.0.2
//---------------------------------------------------------------------------------------------
// Updates from previous version:
// Date        Initials   Reason/s
// DDMMMYYYY   XXX        N)___
// 24Mar2022   CJM        1) Test by having a single nanoparticle of size 20nm in a 200nm x 200nm
//							 square with isolating BC's
// 29Mar2022   CJM        2) Add density to heat source equation
//---------------------------------------------------------------------------------------------


////////////////////
// START PROGRAM  //
////////////////////

//Include mathmatical constants
#define _USE_MATH_DEFINES
//For matmatical constants
#include <cmath>
//For outputting the data to a file
#include <iostream>
//For using in-built mathematics from the standard library
#include <cmath>
//This is the standard library
#include <cstdlib>
//For writing files
#include <fstream>
//To operate on strings
#include <sstream>
//For arrays
#include <vector>
//Needed for std::setw and std::setfill
#include <iomanip>
//Needed to calculate runtime
#include <time.h>


//This is for printing with a standard width
#define FIXOUT(a,b) a.width(75);a << std::left << b;


//Initalise Constants
//---------------------------------------------------------------------------------------------
//Intial space
//Units - m
const double Ix = 0;
//---------------------------------------------------------------------------------------------
//Inital time
//Units - s
const double It = 0;
//---------------------------------------------------------------------------------------------
//Boundry temprature
//Units - Degrees Celcius
const double boundt = 37;
//---------------------------------------------------------------------------------------------

//Initailise functions
//---------------------------------------------------------------------------------------------
//Alpha constant
//Combines all constant values into one to simplify calculations
void ALPHA(std::vector< std::vector<double> > K, std::vector < std::vector<double> > P, std::vector < std::vector<double> > C, std::vector < std::vector<double> >& A, double NX) {

	//X Position Loop
	for (double x = 1; x < NX; x++) {

		//Y Position Loop
		for (double y = 1; y < NX; y++) {

			A[x][y] = K[x][y] / (P[x][y] * C[x][y]);

		}
	}
}
//---------------------------------------------------------------------------------------------
//Stability
//Calcualtes the R constant in the scheme (equivalent to stabilty in explicit scheme) 
void STABLE(std::vector< std::vector<double> > K, std::vector < std::vector<double> > P, std::vector < std::vector<double> > C, std::vector < std::vector<double> >& A, double NX, double DT, double DX, std::vector < std::vector<double> >& R) {

	ALPHA(K, P, C, A, NX);

	//X Position Loop
	for (double x = 1; x < NX; x++) {

		//Y Position Loop
		for (double y = 1; y < NX; y++) {

			R[x][y] = (A[x][y] * DT) / pow(DX, 2);

		}
	}
}
//---------------------------------------------------------------------------------------------
//Heat Source
//Determines amount and location of energy thats inputted into system from particle qi is for time t, qi2 is for time t + deltat
void SOURCE(double NANOT, std::vector< std::vector<double> > Q, std::vector< std::vector<double> >& QI, double NX, double TIME, std::vector< std::vector<double> > P, std::vector< std::vector<double> > C) {

	//X Position Loop
	for (double x = 1; x < NX; x++) {

		//Y Position Loop
		for (double y = 1; y < NX; y++) {

			double a1 = NANOT; //q at time t
			double a2 = NANOT; //q at time t + deltat

			QI[x][y] = ((a1 + a2) / (2.0 * P[x][y] * C[x][y])) * Q[x][y] * P[x][y];

		}
	}

}
//---------------------------------------------------------------------------------------------
//Bij
//The beta value in the Crank - Nicholson scheme
void BEE(std::vector< std::vector<double> > K, std::vector < std::vector<double> > P, std::vector < std::vector<double> > C, std::vector < std::vector<double> >& A, double NX, double DT, double DX, std::vector < std::vector<double> >& R, std::vector< std::vector<double> > NTEMP, std::vector< std::vector<double> >& BIJ, double NANOT, std::vector < std::vector<double> >& QI, std::vector < std::vector<double> >& Q, double TIME) {
	ALPHA(K, P, C, A, NX);
	STABLE(K, P, C, A, NX, DT, DX, R);
	SOURCE(NANOT, Q, QI, NX, TIME, P, C);

	//X Position Loop
	for (double x = 1; x < NX; x++) {

		//Y Position Loop
		for (double y = 1; y < NX; y++) {

			BIJ[x][y] = ((R[x][y] / 2) * (NTEMP[x + 1][y] + NTEMP[x - 1][y] + NTEMP[x][y - 1] + NTEMP[x][y + 1])) + ((1 - (2 * R[x][y])) * NTEMP[x][y]) + (QI[x][y] * DT);
		}
	}
}
//---------------------------------------------------------------------------------------------
//OMEGA
//The overcorrection value in the SOR scheme. Used to reduce itterations.
double OMEGA(double NX) {
	double omega = 4 / (2 + sqrt(4 - (2 * cos(M_PI / NX))));
	return omega;
}
//---------------------------------------------------------------------------------------------
//SOR
double SOR(double RC, double DX, double DT, double TU, double TD, double TL, double TR, double TC, double NX, double bee) {
	double omega = OMEGA(NX);
	double sor = ((1 - omega) * TC) + ((omega / (1 + (2 * RC))) * (((RC / 2) * (TU + TD + TL + TR)) + bee));
	return sor;
}
//---------------------------------------------------------------------------------------------
//Boundries
//Sets boundry conditions for situation
void TSIN(double DX, double DT, double FX, double NX, std::vector< std::vector<double> >& NTEMP) {

	for (int i = 0; i < NX; i++)
	{
		NTEMP[0][i] = NTEMP[1][i];
		NTEMP[NX - 1][i] = NTEMP[NX - 2][i];
		NTEMP[i][0] = NTEMP[i][1];
		NTEMP[i][NX] = NTEMP[i][NX - 1];
	}
}
//---------------------------------------------------------------------------------------------
//Itteration
void ITR(std::vector< std::vector<double> >& R, double DX, double DT, double NX, std::vector< std::vector<double> >& NTEMP, double NI, std::vector< std::vector<double> > BIJ, std::vector< std::vector<double> >& OTEMP, std::vector< std::vector<double> >& ETEMP, double FX, double NANOT, std::vector< std::vector<double> > Q, std::vector< std::vector<double> > QI) {

	//Itteration loop
	for (int ix = 0; ix <= NI; ix++) {

		//X Position Loop
		for (double x = 1; x < NX; x++) {

			//Y Position Loop
			for (double y = 1; y < NX; y++) {

				//Calculate new temprature
				NTEMP[x][y] = SOR(R[x][y], DX, DT, NTEMP[x + 1][y], NTEMP[x - 1][y], NTEMP[x][y - 1], NTEMP[x][y + 1], NTEMP[x][y], NX, BIJ[x][y]);

				//Calculate Error due to itterations
				if (ix == NI - 1) {
					OTEMP[x][y] = NTEMP[x][y];
				}

				if (ix == NI) {
					ETEMP[x][y] = OTEMP[x][y] - NTEMP[x][y];
				}

				//Reapply boundry conditions to temp arrays
				TSIN(DX, DT, FX, NX, NTEMP);

			}
		}
	}
}
//---------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------
//Start of main code
//---------------------------------------------------------------------------------------------

int main() {

	//Arrays for 1D Variables
	std::vector<double> time, spacex, spacey, initial;

	//Arrays for 2D variables
	std::vector< std::vector<double> >  newtemp, bij, oldtemp, error, kay, roh, cee, stable, alpha, q, Qi;

	//Open inital conditions file and check to see that the file was opened correctly:
	std::ifstream ifile("initial_conditions_CN.txt", std::ios::in);
	if (!ifile.is_open()) {
		std::cerr << "There was a problem opening the input file!\n";
		exit(1);//exit or do additional error checking
	}

	double num = 0.0;
	//keep storing values from the text file so long as data exists:
	while (ifile >> num) {
		initial.push_back(num);
	}

	//Set initial conditions as constants
	//-----------------------------------
	//Time step
	const double deltat = initial[0];
	//-----------------------------------
	//Final time
	const double Ft = initial[1];
	//-----------------------------------
	//Space step
	const double deltax = initial[2];
	//-----------------------------------
	//Final space
	const double Fx = initial[3];
	//-----------------------------------
	//Nanoparticle SAR
	const double nanot = initial[4];
	//-----------------------------------
	//Number of itterations
	const double Ni = initial[5];
	//-----------------------------------

	// Calculate number of steps from IC's
	//-----------------------------------
	//Number of time steps
	const int Nt = Ft / deltat;
	//-----------------------------------
	//Number of space steps
	const int Nx = (Fx / deltax) + 1;
	//-----------------------------------

	//Resize and intialise the array for time
	time.resize(Nt);
	std::fill(time.begin(), time.end(), 0.0);

	//Resize and intialise the array for space in x and y directions
	spacex.resize(Nx);
	std::fill(spacex.begin(), spacex.end(), 0.0);

	spacey.resize(Nx);
	std::fill(spacey.begin(), spacey.end(), 0.0);

	//Resize and intialise the array for new/old temp, Error, and bij
	newtemp.resize(Nx + 1, std::vector<double>(Nx + 1, 37));
	bij.resize(Nx + 1, std::vector<double>(Nx + 1, 0));
	oldtemp.resize(Nx + 1, std::vector<double>(Nx + 1, 37));
	error.resize(Nx + 1, std::vector<double>(Nx + 1, 0));
	kay.resize(Nx + 1, std::vector<double>(Nx + 1, 0.4692));
	roh.resize(Nx + 1, std::vector<double>(Nx + 1, 1160));
	cee.resize(Nx + 1, std::vector<double>(Nx + 1, 3600));
	alpha.resize(Nx + 1, std::vector<double>(Nx + 1, 0));
	stable.resize(Nx + 1, std::vector<double>(Nx + 1, 0));
	q.resize(Nx + 1, std::vector<double>(Nx + 1, 0));
	Qi.resize(Nx + 1, std::vector<double>(Nx + 1, 0));

	// Set location and parameters of nanoparticle/s

	for (double x = 6; x < 16; x++) {

		for (double y = 6; y < 16; y++) {

			q[x][y] = 1;
			kay[x][y] = 0.4692;
			roh[x][y] = 1160;
			cee[x][y] = 3600;

		}
	}


	// Apply boundry conditions to temp arrays
	TSIN(deltax, deltat, Fx, Nx, newtemp);

	//for outputting the SAR value as diameter varies
	std::ofstream outpst("2D2M_CN_v2.0.2_dx" + std::to_string(deltax) + "_fx" + std::to_string(Fx) + "_dt" + std::to_string(deltat) + "_ft" + std::to_string(Ft) + "_ni" + std::to_string(Ni) + "_sar" + std::to_string(nanot) + ".dat");
	//check the file is open for writing to
	if (!outpst.is_open())
	{
		std::cerr << "Could not open file, exiting..." << std::endl;
		exit(0);
	}
	else
	{
		outpst << "#time (s) - x position(m) - y position - temprature(C) - Error(C)" << std::endl;
	}

	//two lines to make indices in gnuplot
	outpst << std::endl;
	outpst << std::endl;

	//print the system information at the top of the output file
	FIXOUT(std::cout, "# Space step                 : " << deltax << " [...]" << std::endl);
	FIXOUT(std::cout, "# Final space                : " << Fx << " [...]" << std::endl);
	FIXOUT(std::cout, "# Time step                  : " << deltat << " [...]" << std::endl);
	FIXOUT(std::cout, "# Final time                 : " << Ft << " [...]" << std::endl);
	FIXOUT(std::cout, "# Number of itterations      : " << Ni << " [...]" << std::endl);
	FIXOUT(std::cout, "# SAR Input                  : " << nanot << " [...]" << std::endl);
	FIXOUT(std::cout, "# Outputting data to         : " << "2D2M_CN_v2.0.2.dat" << std::endl);


	//-----------------------------------------------------------------------------------------
	// Start Calculations
	//-----------------------------------------------------------------------------------------

	//Start timer
	clock_t tStart = clock();

	//Time Loop
	for (double t = 0; t < Nt; t++) {

		//Add time to array
		time[t] = deltat * t;

		//Populate space arrays
		for (double x = 0; x < Nx; x++) {

			spacex[x] = x * deltax;

			for (double y = 0; y < Nx; y++) {

				spacey[y] = y * deltax;

			}
		}

		//Calculate new temprature at time
		BEE(kay, roh, cee, alpha, Nx, deltat, deltax, stable, newtemp, bij, nanot, Qi, q, Nt);
		ITR(stable, deltax, deltat, Nx, newtemp, Ni, bij, oldtemp, error, Fx, nanot, q, Qi);

		//Add to output file
		for (double x = 0; x < Nx; x++) {

			for (double y = 0; y < Nx; y++) {

				outpst << time[t] << "\t" << spacex[x] << "\t" << spacey[y] << "\t" << newtemp[x][y] << "\t" << error[x][y] << "\t" << std::endl;

			}

			outpst << std::endl;

		}

		//two lines to make indices in gnuplot
		//outpst << std::endl;

	}

	//Output runtime
	printf("# Time taken                 : %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

	//close the output files
	outpst.close();
	//check the output file closed successfully
	if (outpst.is_open())
	{
		std::cerr << "Warning: could not close output file, continuing..." << std::endl;
	}
	outpst.close();

	time.clear();
	spacex.clear();
	spacey.clear();
	newtemp.clear();
	initial.clear();
	oldtemp.clear();
	kay.clear();
	roh.clear();
	cee.clear();
	alpha.clear();
	stable.clear();
	q.clear();
	Qi.clear();

	//end of the program will return safely if the program is executed successfully
	return(0);
}
