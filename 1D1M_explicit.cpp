//---------------------------------------------------------------------------------------------
// Program Name: 1D1M_explicit.cpp
// 
// Purpose:      Uses the explict method to solve the heat equation in one dimension dispersing
//				 through one medium. Creates a .dat file which can be inputted into gnuplot for 
//				 a graphical representation.
// 
// Input data:   initial_conditions.txt
// Output data:	 1D1M_explicit.dat
//---------------------------------------------------------------------------------------------
// Author:       Callum Mann
// Date created: 01Dec2021
// Version:      1
//---------------------------------------------------------------------------------------------
// Updates:
// Date        Initials   Reason/s
// DDMMMYYYY   XXX        N)___
// 
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


//This is for printing with a standard width
#define FIXOUT(a,b) a.width(75);a << std::left << b;

//Initalise Constants
//---------------------------------------------------------------------------------------------
//Thermal Conductivity
//Units - 
const double kay = 0.5;
//---------------------------------------------------------------------------------------------
//Density
//Units - kg/m^3
const double roh = 1000;
//---------------------------------------------------------------------------------------------
//Specific heat
//Units - 
const double cee = 4000;
//---------------------------------------------------------------------------------------------
//Space step
//Units - m
const double deltax = 1e-9;
//---------------------------------------------------------------------------------------------
//Time step
//Units - s
const double deltat = 0.1;
//---------------------------------------------------------------------------------------------
//Number of space steps
//Units - m
const int Nx = 10000;
//---------------------------------------------------------------------------------------------
//Number of time steps
//Units - s
const double Nt = 10000;
//---------------------------------------------------------------------------------------------
//Intial space
//Units - m
const double Ix = 0;
//---------------------------------------------------------------------------------------------
//Inital time
//Units - s
const double It = 0;
//---------------------------------------------------------------------------------------------
//Final space
//Units - m
const double Fx = 0.01;
//---------------------------------------------------------------------------------------------
//Final time
//Units - s
const double Ft = 1000;
//---------------------------------------------------------------------------------------------
//Temprature of nanoparticle
//Units - Degrees Celcius
const double nanot = 45;
//---------------------------------------------------------------------------------------------
//Boundry temprature
//Units - Degrees Celcius
const double boundt = 0;
//---------------------------------------------------------------------------------------------
 

//Initailise functions
//---------------------------------------------------------------------------------------------
//Alpha constant
//Units - 
double ALPHA(double K, double P, double C) {
    double alpha = K / (P * C);
    return alpha;
}
//---------------------------------------------------------------------------------------------
//Stability
double STABLE(double K, double P, double C, double DX, double DT) {
    double alpha = ALPHA(K, P, C);
    double stable = (alpha * DT) / pow(DX,2);
    return stable;
}
//---------------------------------------------------------------------------------------------
//Explicit Formula
double EXPLICIT(double K, double P, double C, double DX, double DT, double BT, double PT, double AT) {
	double alpha = ALPHA(K, P, C);
	double R = (alpha * DT) / pow(DX, 2);
	double ntemp = (R*BT) + ((1-(2*R))*PT) + (R*AT);
	return ntemp;
}
//---------------------------------------------------------------------------------------------

int main()
{	
	//Arrays for Variables
	std::vector<double> time, space, oldtemp, newtemp, initial, error;

	//Open inital conditions file and check to see that the file was opened correctly:
	std::ifstream ifile("initial_conditions.txt", std::ios::in);
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
		
	//Time step
	const double deltat2 = initial[0];
	//Final time
	const double Ft2 = initial[1];
	//Space step
	const double deltax2 = initial[2];
	//Final space
	const double Fx2 = initial[3];
	//Nanoparticle temp
	const double nanot2 = initial[4];

	//Number of time steps
	const int Nt2= Ft2 / deltat2;
	//Number of space steps
	const int Nx2 = Fx2 / deltax2;

	//Resize and intialise the array for time and space
	time.resize(Nt2*Nx2);
	std::fill(time.begin(), time.end(), 0.0);

	space.resize(Nt2 * Nx2);
	std::fill(space.begin(), space.end(), 0.0);

	//Resize and intialise the array for temp
	oldtemp.resize(Nt2 * Nx2);
	newtemp.resize(Nt2 * Nx2);
	std::fill(oldtemp.begin(), oldtemp.end(), 0.0);
	std::fill(newtemp.begin(), newtemp.end(), 0.0);

	//Resize and intialise the array for error on temp
	error.resize(Nt2 * Nx2);
	std::fill(time.begin(), time.end(), 0.0);

	//Boundry Conditions
	newtemp[Nx2 / 2] = nanot2;
	newtemp[0] = 37;
	newtemp[Nx2+1] = 37;

	//for outputting the SAR value as diameter varies
	std::ofstream outpst("1D1M_explicit.dat");
	//check the file is open for writing to
	if (!outpst.is_open())
	{
		std::cerr << "Could not open file, exiting..." << std::endl;
		exit(0);
	}
	else
	{
		outpst << "#time (s) - postion(m) - temprature(C) - Error(C)" << std::endl;
	}

	//two lines to make indices in gnuplot
	outpst << std::endl;
	outpst << std::endl;

	//print the system information at the top of the output file
	FIXOUT(std::cout, "# Thermal Conductivity       : " << kay << " [...]" << std::endl);
	FIXOUT(std::cout, "# Density                    : " << roh << " [...]" << std::endl);
	FIXOUT(std::cout, "# Specific heat              : " << cee << " [...]" << std::endl);
	FIXOUT(std::cout, "# Space step                 : " << deltax2 << " [...]" << std::endl);
	FIXOUT(std::cout, "# Final space                : " << Fx2 << " [...]" << std::endl);
	FIXOUT(std::cout, "# Time step                  : " << deltat2 << " [...]" << std::endl);
	FIXOUT(std::cout, "# Final time                 : " << Ft2 << " [...]" << std::endl);
	FIXOUT(std::cout, "# Temprature of nanoparticle : " << nanot2 << " [...]" << std::endl);
	FIXOUT(std::cout, "# Outputting data to         : " << "1D1M_explicit.dat" << std::endl);

	//Stability of Model Calculation
	const double stable = STABLE(kay, roh, cee, deltax2, deltat2);

	if (stable > 0.5)
	{
		std::cerr << "----------------------------------------------------------------------------------------------" << std::endl;
		std::cerr << " Warning: Model unstable, amend Space step and Time step in initail_conditions.txt and re-run"  << std::endl;
		std::cerr << "          Current Stability: "<< stable                               << std::endl;
		std::cerr << "          This value needs to be less than 0.5"                       << std::endl;
		std::cerr << "----------------------------------------------------------------------------------------------" << std::endl;
		return(0);
	}
	else
	{
		std::cerr << "----------------------------------------------------------------------------------------------" << std::endl;
		std::cerr << " Model is stable" << std::endl;
		std::cerr << " Current Stability: " << stable << std::endl;
		std::cerr << "----------------------------------------------------------------------------------------------" << std::endl;
	}


	//counter for outputting every N time-steps
	int it_count = 0;

	//Time loop
	for (double t = It; t < Ft2; t += deltat2)
	{
		time[it_count] = t;

		//counter for postion
		int it_counp = 0;

		//set new temps
		oldtemp = newtemp;

		//Position Loop
		for (double x = Ix; x <= Fx2; x += deltax2)
		{
			//Move to new position
			it_counp = it_counp + 1;

			//Calculate new temprature
			newtemp[it_counp] = EXPLICIT(kay, roh, cee, deltax2, deltat2, oldtemp[it_counp - 1], oldtemp[it_counp], oldtemp[it_counp + 1]);

			//Error on newtemp
			error[it_counp] = newtemp[it_counp] * (deltat2 / 2);

			//Reapply Boundry conditions
			newtemp[Nx2 / 2] = nanot2;
			newtemp[0] = 37;
			newtemp[Nx2 + 1] = 37;

			error[Nx2 / 2] = 0;
			error[0] = 0;
			error[Nx2 + 1] = 0;

			//Add postion to array
			space[it_counp] = x;

			//Add result to output file
			outpst << time[it_count] << "\t" << space[it_counp] << "\t" << newtemp[it_counp] << "\t" << error[it_counp] << std::endl;
		}

		//two lines to make indices in gnuplot
		outpst << std::endl;
		outpst << std::endl;

		//Increase time iteration counter
		it_count = it_count + 1;
	}

	//close the output files
	outpst.close();
	//check the output file closed successfully
	if (outpst.is_open())
	{
		std::cerr << "Warning: could not close output file, continuing..." << std::endl;
	}
	outpst.close();

	time.clear();
	space.clear();
	oldtemp.clear();
	newtemp.clear();
	initial.clear();

	//end of the program will return safely if the program is executed successfully
	return(0);
}

///////////////////
//  END PROGRAM  //
///////////////////