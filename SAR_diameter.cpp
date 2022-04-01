//---------------------------------------------------------------------------------------------
// Program Name: SAR_diameter.cpp
// 
// Purpose:      Output .dat file of the SAR value generated by varing sized nonopartilces.
//			     The output will be in a form usable by gnuplot so it can be plotted as a scatter
//				 graph.
// 
// Input data:   _null_
// Output data:	 SAR_diameter.dat
//---------------------------------------------------------------------------------------------
// Author:       Callum Mann
// Date created: 15Nov2021
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
const double eta = 0.00235;
//---------------------------------------------------------------------------------------------
//Boltzmann constant
//Units - J K^1
const double k_b = 1.38e-23;
//---------------------------------------------------------------------------------------------
//Thermodynamic temperature
//Units - K
const double temp = 300;
//---------------------------------------------------------------------------------------------
//Inital Radius
//Units - m
const double radinit = 5e-9;
//---------------------------------------------------------------------------------------------
//Radius step
//Units - Hz
const double drad = 0.5e-12;
//---------------------------------------------------------------------------------------------
//Final radius
//Units - Hz
const double fdrad = 12.5e-9;
//---------------------------------------------------------------------------------------------
//Frequency
//Units - Hz
const double freq = 100e3;
//---------------------------------------------------------------------------------------------
//Saturation magnetisation
//Units - A/m
const double m_s = 400e3;
//---------------------------------------------------------------------------------------------
//Magnetic anisotropy
//Units - J/m^3
const double kay = 1e4;
//---------------------------------------------------------------------------------------------
//Magnetic Permeability
//Units - H/m
const double mu_0 = (4 * M_PI) * pow(10, -7);
//---------------------------------------------------------------------------------------------
//Attempt time
//Units - s
const double tau_0 = 1e-9;
//---------------------------------------------------------------------------------------------
//Applied Magnetic Feild Strength
//Units - A/m
const double H_app = 100 * (pow(10,3) / (4 * M_PI));
//---------------------------------------------------------------------------------------------
//Number of points
//Units - No units
double Nz = 30000;
//---------------------------------------------------------------------------------------------
//Density of Magnetic Material
//Units - Kg/m^3
double roh = 5200;
//---------------------------------------------------------------------------------------------


//Initailise functions
//---------------------------------------------------------------------------------------------
//Volume
//Units - m^3
double VOLUME(double R) {
    double V = (4 / 3) * M_PI * pow(R, 3);
    return V;
}
//---------------------------------------------------------------------------------------------
//Neels Relaxation
double TAU_NEEL(double R, double K, double T) {
    double Vol = VOLUME(R);
    double capgamma = (kay * Vol) / (K * T);
    double tau_n = (tau_0 / 2) * exp(capgamma);
    return tau_n;
}
//---------------------------------------------------------------------------------------------
//Brownian Relaxation
double TAU_BROWN(double R, double E, double T, double K) {
    double Vol = VOLUME(R);
    double v_h = Vol * pow(1 + (delta / R), 3);
    double tau_b = (3 * E * v_h) / (K * T);
    return tau_b;
}
//---------------------------------------------------------------------------------------------
//Time of Relaxation
double TAU(double R, double K, double T, double E) {
    double taub = TAU_BROWN(R,E,T,K);
    double taun = TAU_NEEL(R,K,T);
    double tau = 1 / ((1 / taub) + (1 / taun));
    return tau;
}
//---------------------------------------------------------------------------------------------
//Omega
double OMEGA(double F) {
    double omega = 2 * M_PI * F;
    return omega;
}
//---------------------------------------------------------------------------------------------
//Equilibrium Susceptibility
double CHI_0(double R, double K, double T) {
    double Vol = VOLUME(R);
    double chi_0 = (mu_0 * pow(m_s, 2) * Vol) / (3 * K * T);
    return chi_0;
}
//---------------------------------------------------------------------------------------------
//Susceptability (Real)
double CHIPRIME(double R, double K, double T,double F, double E) {
    double tau = TAU(R, K, T, E);
    double omega = OMEGA(F);
    double chi_0 = CHI_0(R, K, T);
    double chiprime = chi_0 / (1 + pow((omega * tau),2));
    return chiprime;
}
//---------------------------------------------------------------------------------------------
//Susceptability (Imaginary)
double CHIDUBPRIME(double R, double K, double T, double F, double E) {
    double tau = TAU(R, K, T, E);
    double omega = OMEGA(F);
    double chi_0 = CHI_0(R, K, T);
    double chidubprime = (chi_0 * omega * tau) / (1 + pow((omega * tau),2));
    return chidubprime;
}
//---------------------------------------------------------------------------------------------
//Power
double POWER(double R, double K, double T, double F, double E) {
    double chidubprime = CHIDUBPRIME(R, K, T, F, E);
    double power = M_PI * mu_0 * chidubprime * pow(H_app,2) * F;
    return power;
}
//---------------------------------------------------------------------------------------------
//Diameter
double DIAMETER(double R) {
    double diameter = 2 * R;
    return diameter;
}
//---------------------------------------------------------------------------------------------


int main()
{
	//Arrays for Variables dependant on frequency
	std::vector<double> dia, sar, b_tau, n_tau, r_tau;

	//Resize the array for the Diameter
	dia.resize(Nz);
	//Resize the array for the SAR value
	sar.resize(Nz);
	//Resize the array for tau values
	b_tau.resize(Nz);
	n_tau.resize(Nz);
	r_tau.resize(Nz);

	//Initialise the arrays
	std::fill(dia.begin(), dia.end(), (radinit*2));
	std::fill(sar.begin(), sar.end(), 0.0);

	//for outputting the SAR value as diameter varies
	std::ofstream outpst("SAR_diameter.dat");
	//check the file is open for writing to
	if (!outpst.is_open())
	{
		std::cerr << "Could not open file, exiting..." << std::endl;
		exit(0);
	}
	else
	{
		outpst << "#Diameter (m) - SAR (W/g) - Brownian Relaxation (s) - Neel Relaxation (s) - Relaxation constant (s)" << std::endl;
	}

	//two lines to make indices in gnuplot
	outpst << std::endl;
	outpst << std::endl;

	//print the system information at the top of the output file
	FIXOUT(std::cout, "# Applied Magnetic Feild Strength          : " << H_app << " [A/m]" << std::endl);
	FIXOUT(std::cout, "# Density of Magnetic Material             : " << roh << " [A/m]" << std::endl);
	FIXOUT(std::cout, "# Fluid viscosity                          : " << eta << " [units]" << std::endl);
	FIXOUT(std::cout, "# Boltzmann constant                       : " << k_b << " [J K^-1]" << std::endl);
	FIXOUT(std::cout, "# Thermodynamic temperature                : " << temp << " [K]" << std::endl);
	FIXOUT(std::cout, "# Thickness of the surfactant layer        : " << delta << " [m]" << std::endl);
	FIXOUT(std::cout, "# Saturation magnetisation                 : " << m_s << " [A/m]" << std::endl);
	FIXOUT(std::cout, "# Magnetic anisotropy                      : " << kay << " [J/m^3]" << std::endl);
	FIXOUT(std::cout, "# Magnetic Permeability                    : " << mu_0 << " [H/m]" << std::endl);
	FIXOUT(std::cout, "# Frequency                                : " << freq << " [Hz]" << std::endl);
	FIXOUT(std::cout, "# Initial Radius                           : " << radinit << " [m]" << std::endl);
	FIXOUT(std::cout, "# Final radius                             : " << fdrad << " [m]" << std::endl);
	FIXOUT(std::cout, "# Radius step                              : " << drad << " [m]" << std::endl);
	FIXOUT(std::cout, "# Outputting relaxation time by radius to  : " << "SAR_diameter.dat" << std::endl);

	//counter for outputting every N frequency-steps
	int it_count = 0;

	//Loop over radius
	for (double rad = radinit; rad < fdrad; rad += drad)
	{
		std::stringstream filesstr;

		//Add Diameter to array
		dia[it_count] = DIAMETER(rad);

		//Power
		const double power = POWER(rad, k_b, temp, freq, eta);

		//SAR
		sar[it_count] = power / roh;

		//Tau's
		n_tau[it_count] = TAU_NEEL(rad, k_b, temp);

		b_tau[it_count] = TAU_BROWN(rad, eta, temp, k_b);

		r_tau[it_count] = TAU(rad, k_b, temp, eta);


		//Add 1 to count
		it_count++;
	}

	//output full profile
	for (unsigned int i = 0; i < it_count; i++)
	{
		outpst << dia[i] << "\t" << sar[i] << "\t" << b_tau[i] << "\t" << n_tau[i] << "\t" << r_tau[i] << std::endl;
	}

	//close the output files
	outpst.close();
	//check the output file closed successfully
	if (outpst.is_open())
	{
		std::cerr << "Warning: could not close output file, continuing..." << std::endl;
	}
	outpst.close();

	dia.clear();
	sar.clear();
	n_tau.clear();
	b_tau.clear();
	r_tau.clear();

	//end of the program will return safely if the program is executed successfully
	return(0);
}

///////////////////
//  END PROGRAM  //
///////////////////