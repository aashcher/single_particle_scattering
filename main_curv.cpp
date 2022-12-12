/*

Copyright © 2022 Alexey A. Shcherbakov. All rights reserved.

This file is part of the SingleParticleScattering project.

This is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

This sortware is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

*/

#include "./matrix.h"
#include "./particle.h"
#include "./gsm_sphere.h"

#include <math.h>
#include <cmath>
#include <complex>
#include <fstream>

#define USE_MATH_DEFINES

using namespace std;

int main(int argc, char* argv[]) {
	cout.precision(16);

	int N = 15; // maximum degree = N-1
	int NN = N*N;
	int NN2 = 2*NN;
	wl = 1.; wv = 2.*M_PI/wl; // wavelength and wavevector
	R = 0.4; kR = wv*R; // radius
	
	ec = Complex(1.5,0.); ec *= ec; // particle refractive index
	es = 1.; es *= es; // surrounding medium refractive index
	eb = 1.; // basis permittivity for the method
	
	ellipsoid_ab PE(1.*kR,1.2*kR,ec);
	Method3D_GSMSC MGC(N);

	double th0 = 0., ph0 = 0.; // incidence angles
	Vector Vinc(NN2);
	Vinc = MGC.calc_pw(1.,0.,th0,ph0); // incidence amplitude vector

	int ns = 50; // number of spherical shells
	double ca1 = 0.95, ca2 = 1.05; // coefficients before radii of inscribed and superscribed spheres
	
	Vector Vsca(NN2);
	Vsca = MGC.calc_SVt(PE,VI1,PE.R1*ca1,PE.R2*ca2,ec,eb,es,ns); // scattered amplitude vector
	
	Vsca -= Vinc; // substract incident fiels
	cout << MGC.balance1(Vinc, Vsca) << endl; // power conservation accuracy

		// calculate scattering diagaram
	ofstream out_sd("./sd.dat"); // output file
	for (double th=0.1; th<180.5; th+=1.) {
		out_sd << th;
		VF = MGC.calc_far(Vsca, th*M_PI/180., 0.0*M_PI);
		out_sd << " " << abs(VF(1)*VF(1));
		VF = MGC.calc_far(Vsca, th*M_PI/180., 0.5*M_PI);
		out_sd << " " << abs(VF(0)*VF(0));
		out_sd << endl;
	}

	return 0;
}
