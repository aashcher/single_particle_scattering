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

#ifndef _SPFUNC_H
#define _SPFUNC_H

#define _USE_MATH_DEFINES
#include <math.h>

#include "./matrix.h"

void xyz2rtp(double x, double y, double z, double &r, double &th, double &ph);
Vector cart2sph(const Vector&, double, double);
Vector sph2cart(const Vector&, double, double);

namespace sph_const {
	const double M_SQRT3_2 = sqrt(3.)*M_SQRT1_2;
	const double M_1_2_SQRT3 = 0.5*sqrt(3.);
	const double M_SQRT1_2PI = 1./sqrt(2.*M_PI);
	const double M_SQRT2PI = sqrt(2.*M_PI);
}

using namespace std;
using namespace sph_const;
typedef complex<double> Complex;

	//////////////////////////////

inline double flog(int n) {double tv = 0.; for (int i=2; i<n+1; ++i) tv += log(double(i)); return tv;}
inline double flog2(int n) {double tv = 0.; for (int i=n; i>1; i-=2) tv += log(double(i)); return tv;}
inline double flog(int n1, int n2) {double tv = 0.; for (int i=n1; i<n2+1; ++i) tv += log(double(i)); return tv;}

	// spherical Bessel functions //

inline Complex besj0(Complex z) {if (abs(z) < 1.e-7) return 1.; else return sin(z)/z;};
inline Complex besj1(Complex z) {if (abs(z) < 1.e-7) return z/3.; else return (sin(z)-z*cos(z))/z/z;};
inline Complex besj0d(Complex z) {return -besj1(z);}
inline Complex besj1d(Complex z) {if (abs(z) < 1.e-7) return 1./3.; else return besj0(z)-2.*besj1(z)/z;}
inline Complex besy0(Complex z) {return -cos(z)/z;}
inline Complex besy1(Complex z) {return (-cos(z) - z*sin(z))/z/z;}
inline Complex besy0d(Complex z) {return -besy1(z);}
inline Complex besy1d(Complex z) {return besy0(z) - 2.*besy1(z)/z;}
inline Complex besh10(Complex z) {return -j_*exp(j_*z)/z;}
inline Complex besh11(Complex z) {return (-z-j_)*exp(j_*z)/z/z;}
inline Complex besh10d(Complex z) {return -besh11(z);}
inline Complex besh11d(Complex z) {return besh10(z)-2.*besh11(z)/z;}
inline Complex besh20(Complex z) {return j_*exp(-j_*z)/z;}
inline Complex besh21(Complex z) {return (-z+j_)*exp(-j_*z)/z/z;}
inline Complex besh20d(Complex z) {return -besh21(z);}
inline Complex besh21d(Complex z) {return besh20(z)-2.*besh21(z)/z;}

Complex besj(Complex, int);
Complex besjd(Complex, int);
Complex besy(Complex, int);
Complex besyd(Complex, int);
Complex besh1(Complex, int);
Complex besh2(Complex, int);
Complex besh1d(Complex, int);
Complex besh2d(Complex, int);

inline Complex bes_dzj(Complex z, int n) {return besj(z,n)+z*besjd(z,n);};
inline Complex bes_dzy(Complex z, int n) {return besy(z,n)+z*besyd(z,n);};
inline Complex bes_dzh1(Complex z, int n) {return besh1(z,n)+z*besh1d(z,n);};
inline Complex bes_dzh2(Complex z, int n) {return besh2(z,n)+z*besh2d(z,n);};

	// Legendre and associated Legendre polynomials //

inline double pLeg0(double t) {return 1.;}
inline double pLeg1(double t) {return cos(t);}
inline double paLeg11(double t) {return sin((t));}
inline double paLeg1m1(double t) {return -0.5*sin((t));}

double pLeg(double t, int n);
double paLeg(double t, int n, int m);
double paLegd(double t, int n, int m);

inline double pLegn0(double t) {return M_SQRT1_2;} // 1./sqrt(2.)
inline double pLegnd0(double t) {return 0.;}
inline double pLegn1(double t) {return M_SQRT3_2*cos(t);} // sqrt(1.5)
inline double pLegnd1(double t) {return M_SQRT3_2;} // sqrt(1.5)
inline double paLegn11(double t) {return M_1_2_SQRT3*sin(t);} // 0.5*sqrt(3.)
inline double paLegn1m1(double t) {return -M_1_2_SQRT3*sin(t);} // 0.5*sqrt(3.)

double pLegn(double t, int n);
double pLegnd(double t, int n);
double paLegn(double t, int n, int m);
double paLegnd(double t, int n, int m);
void pLegn(double *PL, double t, int n);
void paLegn(double *PL, double t, int n, int m);

double LPi(double t, int n, int m);
double LTau(double t, int n, int m);
double LPin(double t, int n, int m);
double LTaun(double t, int n, int m);

	// spherical functions //

inline Complex sY(double th, double ph, int n, int m) {return paLeg(th,n,m)*exp(j_*double(m)*ph);}
inline Complex sYn(double th, double ph, int n, int m) {return M_SQRT1_2PI*paLegn(th,n,m)*exp(j_*double(m)*ph);}
inline Complex sfRg(Complex z, double th, double ph, int n, int m)
	{return ((n<0)||((n>=0)&&(abs(m)>n))) ?  0. : M_SQRT1_2PI*besj(z,n)*paLegn(th,n,m)*exp(j_*double(m)*ph);};
inline Complex sf1(Complex z, double th, double ph, int n, int m)
	{return ((n<0)||((n>=0)&&(abs(m)>n))) ?  0. : M_SQRT1_2PI*besh1(z,n)*paLegn(th,n,m)*exp(j_*double(m)*ph);};
inline Complex sf2(Complex z, double th, double ph, int n, int m)
	{return ((n<0)||((n>=0)&&(abs(m)>n))) ?  0. : M_SQRT1_2PI*besh2(z,n)*paLegn(th,n,m)*exp(j_*double(m)*ph);};

	// Wigner d-functions

double wig(double, int, int, int);
double wig1(double, int, int, int);
inline Complex Wig(double a, double b, double c, int n, int m1, int m2) {return wig(b,n,m1,m2)*exp(j_*(m1*a + m2*b));}

	// Clebsch-Gordan coefficients

double ClG(int,int,int,int,int);
double ClG1(int,int,int,int,int);
void ClG(double*,int,int,int,int);

int get_nClG(int, int, int, int);
void ClGv(double*, int, int, int, int);

	// Gaunt coefficients
double gaunt_axial(int, int, int, int);

	// spherical vector functions

Vector svfRgM(Complex, double, double, int, int);
Vector svfM1(Complex, double, double, int, int);
Vector svfM2(Complex, double, double, int, int);
Vector svfRgMx(Complex, double, double, int, int);
Vector svfM1x(Complex, double, double, int, int);
Vector svfM2x(Complex, double, double, int, int);
Vector svfRgN(Complex, double, double, int, int);
Vector svfN1(Complex, double, double, int, int);
Vector svfN2(Complex, double, double, int, int);
Vector svfRgNx(Complex, double, double, int, int);
Vector svfN1x(Complex, double, double, int, int);
Vector svfN2x(Complex, double, double, int, int);


#endif
