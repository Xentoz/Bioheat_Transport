//---------------------------------------------------------------------------------------------
// Program Name: magnetic_susceptibility.cpp
// 
// Purpose:      Output .dat file of the susceptibility components as the frequency of an   
//               applied magnetic field varies. The output will be in a form usable by gnuplot
//				 so it can be plotted as a scatter graph.
// 
// Input data:   _null_
// Output data:	 magnetic_susceptability.dat
//---------------------------------------------------------------------------------------------
// Author:       Callum Mann
// Date created: 08Nov2021
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
//Thickness of the surfactant layer
//Units - m
const double delta = 1e-12;
//---------------------------------------------------------------------------------------------
//Fluid viscosity
//Units - 
const double eta = 0.1;
//---------------------------------------------------------------------------------------------
//Boltzmann constant
//Units - J K^1
const double k_b = 1.38e-23;
//---------------------------------------------------------------------------------------------
//Thermodynamic temperature
//Units - K
const double temp = 310.15;
//---------------------------------------------------------------------------------------------
//Radius
//Units - m
const double rad = 4.45e-9;
//---------------------------------------------------------------------------------------------
//Frequency step
//Units - Hz
const double dfreq = 1e2;
//---------------------------------------------------------------------------------------------
//Initial Frequency
//Units - Hz
const double freqinit = 0;
//---------------------------------------------------------------------------------------------
//Overal change in frequency
//Units - Hz
const double fdfreq = 3e6;
//---------------------------------------------------------------------------------------------
//Number of points
//Units - No units
const double Nz = 300000;
//---------------------------------------------------------------------------------------------
//Saturation magnetisation
//Units - A/m
const double m_s = 400e3;
//---------------------------------------------------------------------------------------------
//Magnetic anisotropy
//Units - J/m^3
const double K = 1e4;
//---------------------------------------------------------------------------------------------
//Magnetic Permeability
//Units - H/m
const double mu_0 = (4 * M_PI) * pow(10,-7);
//---------------------------------------------------------------------------------------------
//Attempt time
//Units - s
const double tau_0 = 1e-9;
//---------------------------------------------------------------------------------------------
//Volume
//Units - m^3
const double V = (4/3) * M_PI * pow(rad,3);
//---------------------------------------------------------------------------------------------

int main()
{
	//Arrays for Variables dependant on frequency
	std::vector<double> freq, chiprime, chidubprime, chiprime_chi_0, chidubprime_chi_0, omega;

	//Resize the array for the Frequency
	freq.resize(Nz);
	//Resize the arrays for Real and Imaginary parts of the susceptability, respectively
	chiprime.resize(Nz);
	chidubprime.resize(Nz);
	//Resize the arrays for Real and Imaginary parts of the susceptability over equilibrium, respectively
	chiprime_chi_0.resize(Nz);
	chidubprime_chi_0.resize(Nz);
	//Resize the array for omega
	omega.resize(Nz);


	//Initialise the arrays
	std::fill(freq.begin(), freq.end(), freqinit);
	std::fill(chiprime.begin(), chiprime.end(), 0.0);
	std::fill(chidubprime.begin(), chidubprime.end(), 0.0);
	std::fill(omega.begin(), omega.end(), 0.0);

	//for outputting the Susceptability components as the field varies
	std::ofstream outpst("magnetic_susceptability.dat");
	//check the file is open for writing to
	if (!outpst.is_open())
	{
		std::cerr << "Could not open file, exiting..." << std::endl;
		exit(0);
	}
	else
	{
		outpst << "#Frequency (Hz) - chiprime / chi_0 - chidubprime / chi_0" << std::endl;
	}

	//two lines to make indices in gnuplot
	outpst << std::endl;
	outpst << std::endl;

	//print the system information at the top of the output file
	FIXOUT(std::cout, "# Fluid viscosity                          : " << eta << " [units]" << std::endl);
	FIXOUT(std::cout, "# Boltzmann constant                       : " << k_b << " [J K^-1]" << std::endl);
	FIXOUT(std::cout, "# Thermodynamic temperature                : " << temp << " [K]" << std::endl);
	FIXOUT(std::cout, "# Thickness of the surfactant layer        : " << delta << " [m]" << std::endl);
	FIXOUT(std::cout, "# Saturation magnetisation                 : " << m_s << " [A/m]" << std::endl);
	FIXOUT(std::cout, "# Magnetic anisotropy                      : " << K << " [J/m^3]" << std::endl);
	FIXOUT(std::cout, "# Magnetic Permeability                    : " << mu_0 << " [H/m]" << std::endl);
	FIXOUT(std::cout, "# Radius                                   : " << rad << " [m]" << std::endl);
	FIXOUT(std::cout, "# Initial frequency                        : " << freqinit << " [Hz]" << std::endl);
	FIXOUT(std::cout, "# Final frequency                          : " << fdfreq << " [Hz]" << std::endl);
	FIXOUT(std::cout, "# Frequency step                           : " << dfreq << " [Hz]" << std::endl);
	FIXOUT(std::cout, "# Outputting relaxation time by radius to  : " << "magnetic_susceptability.dat" << std::endl);


	//----------------------------------------------------------------------------------------
	// Caluculate relaxation time
	//----------------------------------------------------------------------------------------
	
	//Magnetic Volume
	const double v_m = (4 / 3) * M_PI * pow(rad, 3);

	//Hydrodynamic Volume
	const double v_h = v_m * pow(1 + (delta / rad), 3);

	//Brownian Relaxation Time
	const double tau_b = (3 * eta * v_h) / (k_b * temp);

	//Capgamma
	const double capgamma = (K * v_m) / (k_b * temp);

	//Neel Relaxation Time 
	const double tau_n = (tau_0 / 2) * exp(capgamma);

	//Relaxation Time
	const double tau = 1 / ((1 / tau_b) + (1 / tau_n));

	//----------------------------------------------------------------------------------------
	// Caluculate susceptability
	//----------------------------------------------------------------------------------------

	//Equilibrium Susceptibility
	const double chi_0 = (mu_0 * pow(m_s, 2) * V) / (3 * k_b * temp);

	//counter for outputting every N frequency-steps
	int it_count = 0;

	//Loop over frequency
	for (double F = freqinit; F < fdfreq; F += dfreq)
	{
		std::stringstream filesstr;

		//Add frequency to array
		freq[it_count] = F;

		//Omega
		omega[it_count] = 2 * M_PI * freq[it_count];

		//Susceptability (Real)
		chiprime[it_count] = chi_0 / (1 + pow((omega[it_count] * tau),2));

		//Susceptability (Imaginary)
		chidubprime[it_count] = (chi_0 * omega[it_count] * tau) / (1 + pow((omega[it_count] * tau), 2));

		//Susceptability / equilibrium (Real)
		chiprime_chi_0[it_count] = chiprime[it_count] / chi_0;

		//Susceptability / equilibrium (Imaginary)
		chidubprime_chi_0[it_count] = chidubprime[it_count] / chi_0;

		//Add 1 to count
		it_count++;
	}

	//output full profile
	for (unsigned int i = 0; i < it_count; i++)
	{
		outpst << freq[i] << "\t" << chiprime_chi_0[i] << "\t" << chidubprime_chi_0[i] << std::endl;
	}

	//close the output files
	outpst.close();
	//check the output file closed successfully
	if (outpst.is_open())
	{
		std::cerr << "Warning: could not close output file, continuing..." << std::endl;
	}
	outpst.close();

	freq.clear();
	chiprime.clear();
	chidubprime.clear();
	chiprime_chi_0.clear();
	chidubprime_chi_0.clear();
	omega.clear();

	//end of the program will return safely if the program is executed successfully
	return(0);
}

///////////////////
//  END PROGRAM  //
///////////////////