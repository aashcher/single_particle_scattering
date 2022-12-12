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

#include "./particle.h"
#include "../sfun/spfunc.h"
#include "../sfun/fun.h"
#include "../sfun/fastgl.h"

#include <math.h>
#define USE_MATH_DEFINES

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Matrix sphere::calc_de_r(double kr, Complex eb, Complex es) const {
	Complex t1, t2, t3, t4;
	t1 = (ep-eb)/ep, t2 = (ep-eb)/eb; //M_SQRT2*
	t3 = (es-eb)/es, t4 = (es-eb)/eb; //M_SQRT2*
	Matrix DE(2,1);
	if (kr < kR) {DE.Data[0] = t1; DE.Data[1] = t2;}
	else {DE.Data[0] = t3; DE.Data[1] = t4;}
	return DE;
}

Matrix sphere::calc_de_r(int NL, double *kr, Complex eb, Complex es) const {
	Complex t1, t2, t3, t4;
	t1 = (ep-eb)/ep, t2 = (ep-eb)/eb; // M_SQRT2*
	t3 = (es-eb)/es, t4 = (es-eb)/eb; // M_SQRT2*
	Matrix DE(2,NL);
	for (int k=0; k<NL; ++k) {
		if (kr[k] < kR) {DE.Data[2*k] = t1; DE.Data[2*k+1] = t2;}
		else {DE.Data[k] = t3; DE.Data[NL+k] = t4;}
	}
	return DE;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// multilayer sphere
Matrix sphere_ml::calc_de_r(int NL, double *kr, Complex eb, Complex es) const {
	int k, ii = 0; double rr;
	Complex t1, t2, ee;
	while ((ii < ns)&&(kr[0] > kR[ii])) ++ii;
	if (ii == ns) {ee = es; rr = 2.*kr[NL-1];} else {ee = ep[ii]; rr = kR[ii];}
	t1 = (ee-eb)/ee, t2 = (ee-eb)/eb;//M_SQRT2*
	Matrix DE(2,NL);
	for (int k=0; k<NL; ++k) {
		if (kr[k] > rr) {
			if (++ii >= ns) {ee = es; rr = 2.*kr[NL-1];} else {ee = ep[ii]; rr = kR[ii];}
			t1 = (ee-eb)/ee, t2 = (ee-eb)/eb; // M_SQRT2*
		}
		DE.Data[k] = t1; DE.Data[NL+k] = t2;
	}
	return DE;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// spheroids
Matrix ellipsoid_ab::calc_de_t(int N, double kr, Complex eb, Complex es) const {
	double th, tv, tvv;
	Complex tc1, tc2, tc3, tc4;
	Matrix DE(2,2*N);
	memset(DE.Data,0,4*N*sizeof(Complex));
	if (kz > kx) {
		tc1 = tc2 = M_SQRT2*(ep - eb); tc1 /= ep; tc2 /= eb;
		tc3 = tc4 = M_SQRT2*(es - eb); tc3 /= es; tc4 /= eb;
	}
	else {
		tc3 = tc4 = M_SQRT2*(ep - eb); tc3 /= ep; tc4 /= eb;
		tc1 = tc2 = M_SQRT2*(es - eb); tc1 /= es; tc2 /= eb;
	}
	if (kr <= R1) {DE.Data[0] = tc1; DE.Data[2*N] = tc2;}
	else if (kr >= R2) {DE.Data[0] = tc3; DE.Data[2*N] = tc4;}
	else {
		//th = acos(kx/kr*sqrt((kr*kr - kz*kz)/(kx*kx - kz*kz))); tv = cos(th);
		th = acos(kz/kr*sqrt((kr*kr - kx*kx)/(kz*kz - kx*kx))); tv = cos(th);
		DE.Data[0] = (1.-tv)*tc1 + tv*tc3; DE.Data[2*N] = (1.-tv)*tc2 + tv*tc4;
		for (int n=1; n<N; ++n) {
			tv = M_SQRT2*( pLegn(th,2*n+1)/sqrt(4*n+3.) - pLegn(th,2*n-1)/sqrt(4*n-1.) )/sqrt(4*n+1.);
			DE.Data[2*n] = (tc3 - tc1)*tv; DE.Data[2*N+2*n] = (tc4 - tc2)*tv;
		}
	}
	return DE;
}

Matrix ellipsoid_ab::calc_de_p(int N, double kr, Complex eb, Complex es) const {
	int NN = N*N; double th, tv, tvv;
	Complex tc1, tc2, tc3, tc4;
	Matrix DE(2,4*NN);
	memset(DE.Data,0,8*NN*sizeof(Complex));
	if (kz > kx) {
		tc1 = tc2 = M_SQRT2*(ep - eb); tc1 /= ep; tc2 /= eb;
		tc3 = tc4 = M_SQRT2*(es - eb); tc3 /= es; tc4 /= eb;
	}
	else {
		tc3 = tc4 = M_SQRT2*(ep - eb); tc3 /= ep; tc4 /= eb;
		tc1 = tc2 = M_SQRT2*(es - eb); tc1 /= es; tc2 /= eb;
	}
	if (kr <= R1) {DE.Data[0] = tc1; DE.Data[4*NN] = tc2;}
	else if (kr >= R2) {DE.Data[0] = tc3; DE.Data[4*NN] = tc4;}
	else {
		//th = acos(kx/kr*sqrt((kr*kr - kz*kz)/(kx*kx - kz*kz))); tv = cos(th);
		th = acos(kz/kr*sqrt((kr*kr - kx*kx)/(kz*kz - kx*kx))); tv = cos(th);
		DE.Data[0] = (1.-tv)*tc1 + tv*tc3; DE.Data[4*NN] = (1.-tv)*tc2 + tv*tc4;
		for (int n=1; n<N; ++n) {
			tv = M_SQRT2*( pLegn(th,2*n+1)/sqrt(4*n+3.) - pLegn(th,2*n-1)/sqrt(4*n-1.) )/sqrt(4*n+1.);
			DE.Data[2*n*(2*n+1)] = (tc3 - tc1)*tv; DE.Data[4*NN+2*n*(2*n+1)] = (tc4 - tc2)*tv;
		}
	}
	return DE;
}

Matrix ellipsoid_ab::calc_de_t(int N, int NL, double *kr, Complex eb, Complex es) const {
	double th, tv, tvv;
	//cout<<kx<<" "<<kz<<endl; cin.get();
	Complex tc1, tc2, tc3, tc4;
	Matrix DE(2*NL,2*N);
	memset(DE.Data,0,4*NL*N*sizeof(Complex));
	if (kz > kx) {
		tc1 = tc2 = M_SQRT2*(ep - eb); tc1 /= ep; tc2 /= eb;
		tc3 = tc4 = M_SQRT2*(es - eb); tc3 /= es; tc4 /= eb;
	}
	else {
		tc3 = tc4 = M_SQRT2*(ep - eb); tc3 /= ep; tc4 /= eb;
		tc1 = tc2 = M_SQRT2*(es - eb); tc1 /= es; tc2 /= eb;
	}
	for (int k=0; k<NL; ++k) {
		if (kr[k] <= R1) {DE.Data[k*2*N] = tc1; DE.Data[(NL+k)*2*N] = tc2;}
		else if (kr[k] >= R2) {DE.Data[k*2*N] = tc3; DE.Data[(NL+k)*2*N] = tc4;}
		else {
			th = acos(kz/kr[k]*sqrt((kr[k]*kr[k] - kx*kx)/(kz*kz - kx*kx))); tv = cos(th);
			DE.Data[k*2*N] = (1.-tv)*tc1 + tv*tc3; DE.Data[(NL+k)*2*N] = (1.-tv)*tc2 + tv*tc4;
			for (int n=1; n<N; ++n) {
				tv = M_SQRT2*( pLegn(th,2*n+1)/sqrt(4*n+3.) - pLegn(th,2*n-1)/sqrt(4*n-1.) )/sqrt(4*n+1.);
				DE.Data[k*2*N+2*n] = (tc3 - tc1)*tv; DE.Data[(NL+k)*2*N+2*n] = (tc4 - tc2)*tv;
			}
		}
	}
	return DE;
}

Matrix ellipsoid_ab::calc_de_p(int N, int NL, double *kr, Complex eb, Complex es) const {
	int NN = N*N; double th, tv, tvv;
	Complex tc1, tc2, tc3, tc4;
	Matrix DE(2*NL,4*NN);
	memset(DE.Data,0,8*NL*NN*sizeof(Complex));
	if (kz > kx) {
		tc1 = tc2 = M_SQRT2*(ep - eb); tc1 /= ep; tc2 /= eb;
		tc3 = tc4 = M_SQRT2*(es - eb); tc3 /= es; tc4 /= eb;
	}
	else {
		tc3 = tc4 = M_SQRT2*(ep - eb); tc3 /= ep; tc4 /= eb;
		tc1 = tc2 = M_SQRT2*(es - eb); tc1 /= es; tc2 /= eb;
	}
	for (int k=0; k<NL; ++k) {
		if (kr[k] <= R1) {DE.Data[k*4*NN] = tc1; DE.Data[(NL+k)*4*NN] = tc2;}
		else if (kr[k] >= R2) {DE.Data[k*4*NN] = tc3; DE.Data[(NL+k)*4*NN] = tc4;}
		else {
			th = acos(kz/kr[k]*sqrt((kr[k]*kr[k] - kx*kx)/(kz*kz - kx*kx))); tv = cos(th);
			DE.Data[k*4*NN] = (1.-tv)*tc1 + tv*tc3; DE.Data[(NL+k)*4*NN] = (1.-tv)*tc2 + tv*tc4;
			for (int n=1; n<N; ++n) {
				tv = M_SQRT2*( pLegn(th,2*n+1)/sqrt(4*n+3.) - pLegn(th,2*n-1)/sqrt(4*n-1.) )/sqrt(4*n+1.);
				DE.Data[k*4*NN+2*n*(2*n+1)] = (tc3 - tc1)*tv; DE.Data[(NL+k)*4*NN+2*n*(2*n+1)] = (tc4 - tc2)*tv;
			}
		}
	}
	return DE;
}

Matrix ellipsoid_ab::calc_de3(int nx, int ny, int nz, double kdx, double kdy, double kdz, Complex eb) const {
	int ix, iy, iz, ii, nn = nx*ny*nz;
	double dx, dy, dz, tv1, tv2, r1, r2, a = -1./5., tt = kz/kx;
	Complex de = ep - eb, dde = eb - eb;
	Matrix M(3,nn);

	dx = 1./double(nx); dy = 1./double(ny); dz = 1./double(nz); tt *= tt;

	tt = kz/kx;
	for (ix=0; ix<nx; ix++) for (iy=0; iy<ny; iy++) for (iz=0; iz<nz; iz++) {
		r1 = r2 = 0;
		tv1 = (ix-nx*0.5)*(ix-nx*0.5)*tt*tt; tv2 = (ix+1-nx*0.5)*(ix+1-nx*0.5)*tt*tt;
		if (tv1 > tv2) {r1 += tv2*dx*dx; r2 += tv1*dx*dx;} else {r1 += tv1*dx*dx; r2 += tv2*dx*dx;}
		tv1 = (iy-ny*0.5)*(iy-ny*0.5)*tt*tt; tv2 = (iy+1-ny*0.5)*(iy+1-ny*0.5)*tt*tt;
		if (tv1 > tv2) {r1 += tv2*dy*dy; r2 += tv1*dy*dy;} else {r1 += tv1*dy*dy; r2 += tv2*dy*dy;}
		tv1 = (iz-nz*0.5)*(iz-nz*0.5); tv2 = (iz+1-nz*0.5)*(iz+1-nz*0.5);
		if (tv1 > tv2) {r1 += tv2*dz*dz; r2 += tv1*dz*dz;} else {r1 += tv1*dz*dz; r2 += tv2*dz*dz;}
		r1 = sqrt(r1); r2 = sqrt(r2); ii = ix*ny*nz+iy*nz+iz;
		if (r2 <= 0.5) M.Data[ii] = M.Data[nn+ii] = M.Data[2*nn+ii] = de;
		else if (r1 >= 0.5) M.Data[ii] = M.Data[nn+ii] = M.Data[2*nn+ii] = dde;
		else {
			tv1 = (0.5-r1)/(r2-r1);
			if (tv1 < 0.5) tv2 = (tv1 < 1./3.) ? (6*tv1*tv1+4*tv1+1)*3*tv1/13 : tv1*(24*tv1+1)/13;
			else tv2 = (tv1 > 2./3.) ? 1-(6*(1-tv1)*(1-tv1)+4*(1-tv1)+1)*3*(1-tv1)/13 : 1-(1-tv1)*(24*(1-tv1)+1)/13;
			M.Data[ii] = M.Data[nn+ii] = M.Data[2*nn+ii] = aver(ep,eb,tv2,a)-eb;
		}
	}

	return M;
/**/
}

//////////////////////

void ellipsoid_ab::get_fgi(int n, double *p, double *ff, int Ni) const {
		// p[0]=b^2; p[1]=a^2-b^2; p[2]=b^4; p[3]=a^4-b^4;
		// p[4]=a'b^3; p[5]=b'a^3-a'b^3; p[6] = (b^2-a^2)^2
	double th, tw, tL0, tL1, tL2;
	double x, x2, xs, xs2;
	double fg[7], ts, tc, ts2, tc2, tvd, tv2, tv32;
	fastgl::QuadPair qp;

	memset(ff,0,8*sizeof(double));

	for (int i=0; i<Ni; ++i) {
		qp = fastgl::GLPair(Ni,i+1);
		x = qp.x(); //th = M_PI*(qp.x()+1.)*0.5;
		tw = qp.weight; //tw = 0.5*M_PI*qp.weight*sin(th);
		th = acos(x);

		//cout<<i<<"\t"<<x<<"   "<<180.*th/M_PI<<endl;

		x2 = x*x; xs2 = x2*(1. - x2); xs = x*sqrt(1. - x2);
		tvd = (p[4] + p[5]*x2) / (p[2] + p[3]*x2);
		tv2 = sqrt(p[0] + p[1]*x2); tv32 = tv2*tv2*tv2;

		fg[0] = tv32*tvd;
		fg[1] = xs*tvd;
		fg[2] = xs*tv2*tvd;
		fg[3] = (2.*tv2 + p[6]*xs2/tv32)*tvd;
		fg[4] = xs2*tvd/tv32;
		fg[5] = xs2*tvd/tv2;
		fg[6] = tv2*tv2*tvd;

		tL0 = pLegn(th,n)*tw;
		tL1 = paLegn(th,n,1)*tw;
		tL2 = paLegn(th,n,2)*tw;

		ff[0] += fg[0]*tL0;
		ff[1] += fg[1]*tL1;
		ff[2] += fg[2]*tL1;
		ff[3] += fg[3]*tL0;
		ff[4] += fg[4]*tL2;
		ff[5] += fg[5]*tL0;
		ff[6] += fg[6]*tL0;
		ff[7] += fg[5]*tL2;
	}
	//cin.get();
}

Matrix ellipsoid_ab::calc_ge_t(int N, int NL, double *kr, double kr1, double kr2, double kr0, Complex rib, Complex es) const {
		// kr - in basis medium; kr0,1,2 - in vacuum; kx, kz - in vacuum
	int i, k, n, nn, ni, nii;
	double ka, kb, th, ga, gb, dga, dgb, p[8], fg[7], tm[8], tw, tL0, tL1, tL2, tn;//kr0, 
	Complex tc1, tc2, tc3, tc4;
	Matrix M(8*NL,2*N);

	memset(M.Data,0,16*NL*N*sizeof(Complex));

	ni = 10*N*N;

	ka = kx; kb = kz;
	//ka = kz; kb = kx;
	for (k=0; k<NL; ++k) {
		if (kr[k] < rib.real()*kr0) {
			ga = (kr[k]*(ka-kr1) - rib.real()*kr1*(ka-kr0))/(kr0-kr1); dga = (ka-kr1)/(kr0-kr1); // in basis medium
			gb = (kr[k]*(kb-kr1) - rib.real()*kr1*(kb-kr0))/(kr0-kr1); dgb = (kb-kr1)/(kr0-kr1); // in basis medium
		}
		else {
			ga = (kr[k]*(kr2-ka) + rib.real()*kr2*(ka-kr0))/(kr2-kr0); dga = (kr2-ka)/(kr2-kr0);
			gb = (kr[k]*(kr2-kb) + rib.real()*kr2*(kb-kr0))/(kr2-kr0); dgb = (kr2-kb)/(kr2-kr0);
		}

		p[1] = ga*ga - (p[0] = gb*gb);
		p[3] = pow(ga,4) - (p[2] = pow(gb,4));
		p[5] = dgb*pow(ga,3) - (p[4] = dga*pow(gb,3));
		p[6] = p[1]*p[1];

		for (n=0; n<N; ++n) {
			nn = 2*n;
			get_fgi(nn,p,tm,ni);
			M.Data[(k+0*NL)*2*N+nn] = M_SQRT1_2*kr[k]*kr[k]/(ga*ga*gb*gb) * tm[0];
			M.Data[(k+1*NL)*2*N+nn] = M_SQRT1_2*0.5*kr[k]*p[1]/(ga*gb) * tm[1];
			M.Data[(k+2*NL)*2*N+nn] = M_SQRT1_2*0.5*kr[k]*p[1]/(ga*ga*gb*gb) * tm[2];
			M.Data[(k+3*NL)*2*N+nn] = 0.25*M_SQRT1_2 * tm[3];
			M.Data[(k+4*NL)*2*N+nn] = 0.25*M_SQRT1_2*p[6] * tm[4];
			M.Data[(k+5*NL)*2*N+nn] = 0.25*M_SQRT1_2*p[6]/(ga*ga*gb*gb) * tm[5];
			M.Data[(k+6*NL)*2*N+nn] = 0.5*M_SQRT1_2/(ga*gb) * tm[6];
			M.Data[(k+7*NL)*2*N+nn] = 0.25*M_SQRT1_2*p[6]/(ga*ga*gb*gb) * tm[7];
		}
	}

	return M;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// simply connected cassini oval
Matrix cassini::calc_de_t(int N, double kr, Complex eb, Complex es) const {
	int n;
	double r2, r4, th, tv, tvv, c4, a2, a4;
	Complex tc1, tc2, tc3, tc4;
	Matrix M(2,2*N);
	memset(M.Data,0,4*N*sizeof(Complex));

	tc1 = tc2 = M_SQRT2*(ep - eb); tc1 /= ep; tc2 /= eb;
	tc3 = tc4 = M_SQRT2*(es - eb); tc3 /= es; tc4 /= eb;
	c4 = pow(kc,4); a2 = ka*ka; a4 = a2*a2;

	r2 = kr*kr; r4 = r2*r2;
	th = acos(tv = 0.5*(1.0 + 0.5*(c4-a4-r4)/(r2*a2)));
	M.Data[0] = (1.-tv)*tc2 + tv*tc4; M.Data[2*N] = (1.-tv)*tc1 + tv*tc3;
	for (n=1; n<N; ++n) {
		tv = M_SQRT2*( pLegn(th,2*n+1)/sqrt(double(4*n+3))
			- pLegn(th,2*n-1)/sqrt(double(4*n-1)) )/sqrt(double(4*n+1));
		M.Data[2*n] = (tc4 - tc2)*tv; M.Data[2*N+2*n] = (tc3 - tc1)*tv;
	}

	return M;
}

Matrix cassini::calc_de_p(int N, double kr, Complex eb, Complex es) const {
	int n, NN = N*N;
	double r2, r4, th, tv, tvv, c4, a2, a4;
	Complex tc1, tc2, tc3, tc4;
	Matrix M(2,4*NN);
	memset(M.Data,0,8*NN*sizeof(Complex));

	tc1 = tc2 = M_SQRT2*(ep - eb); tc1 /= ep; tc2 /= eb;
	tc3 = tc4 = M_SQRT2*(es - eb); tc3 /= es; tc4 /= eb;
	c4 = pow(kc,4); a2 = ka*ka; a4 = a2*a2;

	r2 = kr*kr; r4 = r2*r2;
	th = acos(tv = 0.5*(1.0 + 0.5*(c4-a4-r4)/(r2*a2)));
	M.Data[0] = (1.-tv)*tc2 + tv*tc4; M.Data[4*NN] = (1.-tv)*tc1 + tv*tc3;
	for (n=1; n<N; ++n) {
		tv = M_SQRT2*( pLegn(th,2*n+1)/sqrt(double(4*n+3))
			- pLegn(th,2*n-1)/sqrt(double(4*n-1)) )/sqrt(double(4*n+1));
		M.Data[2*n*(2*n+1)] = (tc4 - tc2)*tv; M.Data[4*NN+2*n*(2*n+1)] = (tc3 - tc1)*tv;
	}

	return M;
}

Matrix cassini::calc_de_t(int N, int NL, double *kr, Complex eb, Complex es) const {
	int n;
	double r2, r4, th, tv, tvv, c4, a2, a4;
	Complex tc1, tc2, tc3, tc4;
	Matrix M(2*NL,2*N);
	memset(M.Data,0,4*NL*N*sizeof(Complex));

	tc1 = tc2 = M_SQRT2*(ep - eb); tc1 /= ep; tc2 /= eb;
	tc3 = tc4 = M_SQRT2*(es - eb); tc3 /= es; tc4 /= eb;
	c4 = pow(kc,4); a2 = ka*ka; a4 = a2*a2;

	for (int k=0; k<NL; ++k) {
		r2 = kr[k]*kr[k]; r4 = r2*r2;
		//th = 0.5*acos(0.5*(cR4-aR4-rr4)/rr2/aR2);
		th = acos(tv = 0.5*(1.0 + 0.5*(c4-a4-r4)/(r2*a2)));
		M.Data[k*2*N] = (1.-tv)*tc2 + tv*tc4; M.Data[(k+NL)*2*N] = (1.-tv)*tc1 + tv*tc3;
		for (n=1; n<N; ++n) {
			tv = M_SQRT2*( pLegn(th,2*n+1)/sqrt(double(4*n+3))
				- pLegn(th,2*n-1)/sqrt(double(4*n-1)) )/sqrt(double(4*n+1));
			M.Data[k*2*N+2*n] = (tc4 - tc2)*tv;
			M.Data[(k+NL)*2*N+2*n] = (tc3 - tc1)*tv;
		}
	}
	return M;
}

Matrix cassini::calc_de_p(int N, int NL, double *kr, Complex eb, Complex es) const {
	int n, NN = N*N;
	double r2, r4, th, tv, tvv, c4, a2, a4;
	Complex tc1, tc2, tc3, tc4;
	Matrix M(2*NL,4*NN);
	memset(M.Data,0,8*NL*NN*sizeof(Complex));

	tc1 = tc2 = M_SQRT2*(ep - eb); tc1 /= ep; tc2 /= eb;
	tc3 = tc4 = M_SQRT2*(es - eb); tc3 /= es; tc4 /= eb;
	c4 = pow(kc,4); a2 = ka*ka; a4 = a2*a2;

	for (int k=0; k<NL; ++k) {
		r2 = kr[k]*kr[k]; r4 = r2*r2;
		th = acos(tv = 0.5*(1.0 + 0.5*(c4-a4-r4)/(r2*a2)));
		M.Data[k*4*NN] = (1.-tv)*tc2 + tv*tc4; M.Data[(k+NL)*4*NN] = (1.-tv)*tc1 + tv*tc3;
		for (n=1; n<N; ++n) {
			tv = M_SQRT2*( pLegn(th,2*n+1)/sqrt(double(4*n+3))
				- pLegn(th,2*n-1)/sqrt(double(4*n-1)) )/sqrt(double(4*n+1));
			M.Data[k*4*NN+2*n*(2*n+1)] = (tc4 - tc2)*tv;
			M.Data[(NL+k)*4*NN+2*n*(2*n+1)] = (tc3 - tc1)*tv;
		}
	}
	return M;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// cylinder
Matrix cylinder::calc_de_t(int N, double kr, Complex eb, Complex es) const {
	double th1, th2, tv, tvv;
	Complex tc1, tc2, tc3, tc4;
	Matrix DE(2,2*N);
	memset(DE.Data,0,4*N*sizeof(Complex));

	tc1 = tc2 = M_SQRT2*(es - eb); tc1 /= es; tc2 /= eb;
	tc3 = tc4 = M_SQRT2*(ep - eb); tc3 /= ep; tc4 /= eb;

	if (kr <= R1) {DE.Data[0] = tc1; DE.Data[2*N] = tc2;}
	else if (kr >= R2) {DE.Data[0] = tc3; DE.Data[2*N] = tc4;}
	else {
		th1 = (kr > 0.5*kH) ? acos(0.5*kH/kr) : 0.;
		th2 = (kr > kR) ? asin(kR/kr) : 0.5*M_PI;
		tv = cos(th1) - cos(th2);
		DE.Data[0] = (1.-tv)*tc1 + tv*tc3; DE.Data[2*N] = (1.-tv)*tc2 + tv*tc4;
		for (int n=1; n<N; ++n) {
			tv = M_SQRT2*( (pLegn(th1,2*n+1) - pLegn(th2,2*n+1))/sqrt(4*n+3.)
				- (pLegn(th1,2*n-1) - pLegn(th2,2*n-1))/sqrt(4*n-1.) )/sqrt(4*n+1.);
			DE.Data[2*n] = (tc3 - tc1)*tv; DE.Data[2*N+2*n] = (tc4 - tc2)*tv;
		}
	}

	return DE;
}

Matrix cylinder::calc_de_p(int N, double kr, Complex eb, Complex es) const {
	int NN = N*N; double th1, th2, tv, tvv;
	Complex tc1, tc2, tc3, tc4;
	Matrix DE(2,4*NN);
	memset(DE.Data,0,8*NN*sizeof(Complex));

	tc1 = tc2 = M_SQRT2*(es - eb); tc1 /= es; tc2 /= eb;
	tc3 = tc4 = M_SQRT2*(ep - eb); tc3 /= ep; tc4 /= eb;

	if (kr <= R1) {DE.Data[0] = tc1; DE.Data[4*NN] = tc2;}
	else if (kr >= R2) {DE.Data[0] = tc3; DE.Data[4*NN] = tc4;}
	else {
		th1 = (kr > 0.5*kH) ? acos(0.5*kH/kr) : 0.;
		th2 = (kr > kR) ? asin(kR/kr) : 0.5*M_PI;
		tv = cos(th1) - cos(th2);
		DE.Data[0] = (1.-tv)*tc1 + tv*tc3; DE.Data[4*NN] = (1.-tv)*tc2 + tv*tc4;
		for (int n=1; n<N; ++n) {
			tv = M_SQRT2*( (pLegn(th1,2*n+1) - pLegn(th2,2*n+1))/sqrt(4*n+3.)
				- (pLegn(th1,2*n-1) - pLegn(th2,2*n-1))/sqrt(4*n-1.) )/sqrt(4*n+1.);
			DE.Data[2*n*(2*n+1)] = (tc3 - tc1)*tv; DE.Data[4*NN+2*n*(2*n+1)] = (tc4 - tc2)*tv;
		}
	}

	return DE;
}

Matrix cylinder::calc_de_t(int N, int NL, double *kr, Complex eb, Complex es) const {
	double th1, th2, tv, tvv;
	Complex tc1, tc2, tc3, tc4;
	Matrix DE(2*NL,2*N);
	memset(DE.Data,0,4*NL*N*sizeof(Complex));

	tc1 = tc2 = M_SQRT2*(es - eb); tc1 /= es; tc2 /= eb;
	tc3 = tc4 = M_SQRT2*(ep - eb); tc3 /= ep; tc4 /= eb;

	for (int k=0; k<NL; ++k) {
		if (kr[k] <= R1) {DE.Data[k*2*N] = tc1; DE.Data[(NL+k)*2*N] = tc2;}
		else if (kr[k] >= R2) {DE.Data[k*2*N] = tc3; DE.Data[(NL+k)*2*N] = tc4;}
		else {
			th1 = (kr[k] > 0.5*kH) ? acos(0.5*kH/kr[k]) : 0.;
			th2 = (kr[k] > kR) ? asin(kR/kr[k]) : 0.5*M_PI;
			tv = cos(th1) - cos(th2);
			DE.Data[k*2*N] = (1.-tv)*tc1 + tv*tc3; DE.Data[(NL+k)*2*N] = (1.-tv)*tc2 + tv*tc4;
			for (int n=1; n<N; ++n) {
				tv = M_SQRT2*( (pLegn(th1,2*n+1) - pLegn(th2,2*n+1))/sqrt(4*n+3.)
					- (pLegn(th1,2*n-1) - pLegn(th2,2*n-1))/sqrt(4*n-1.) )/sqrt(4*n+1.);
				DE.Data[k*2*N+2*n] = (tc3 - tc1)*tv; DE.Data[(NL+k)*2*N+2*n] = (tc4 - tc2)*tv;
			}
		}
	}

	return DE;
}

Matrix cylinder::calc_de_p(int N, int NL, double *kr, Complex eb, Complex es) const {
	int NN = N*N; double th1, th2, tv, tvv;
	Complex tc1, tc2, tc3, tc4;
	Matrix DE(2*NL,4*NN);
	memset(DE.Data,0,8*NL*NN*sizeof(Complex));

	tc1 = tc2 = M_SQRT2*(es - eb); tc1 /= es; tc2 /= eb;
	tc3 = tc4 = M_SQRT2*(ep - eb); tc3 /= ep; tc4 /= eb;

	for (int k=0; k<NL; ++k) {
		if (kr[k] <= R1) {DE.Data[k*4*NN] = tc1; DE.Data[(NL+k)*4*NN] = tc2;}
		else if (kr[k] >= R2) {DE.Data[k*4*NN] = tc3; DE.Data[(NL+k)*4*NN] = tc4;}
		else {
			th1 = (kr[k] > 0.5*kH) ? acos(0.5*kH/kr[k]) : 0.;
			th2 = (kr[k] > kR) ? asin(kR/kr[k]) : 0.5*M_PI;
			tv = cos(th1) - cos(th2);
			DE.Data[k*4*NN] = (1.-tv)*tc1 + tv*tc3; DE.Data[(NL+k)*4*NN] = (1.-tv)*tc2 + tv*tc4;
			for (int n=1; n<N; ++n) {
				tv = M_SQRT2*( (pLegn(th1,2*n+1) - pLegn(th2,2*n+1))/sqrt(4*n+3.)
					- (pLegn(th1,2*n-1) - pLegn(th2,2*n-1))/sqrt(4*n-1.) )/sqrt(4*n+1.);
				DE.Data[k*4*NN+2*n*(2*n+1)] = (tc3 - tc1)*tv; DE.Data[(NL+k)*4*NN+2*n*(2*n+1)] = (tc4 - tc2)*tv;
			}
		}
	}
	return DE;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// parallelepiped
Matrix box::calc_de3(int nx, int ny, int nz, double kdx, double kdy, double kdz, Complex eb) const {
	int ii, nn = nx*ny*nz; Matrix M(3,nn);
	for (int ii=0; ii<nn; ++ii) {M.Data[ii] = epx - eb; M.Data[nn+ii] = epy - eb; M.Data[2*nn+ii] = epz - eb;}
	return M;
}
/**/
Matrix box::calc_de_p(int N, double kr, Complex eb, Complex es) const {
	int k, i, np = 100*N, m, n, NN = N*N;
	double xx, ww, tt, wws, *pl, df;
	Complex t1, t2, tc, ep = epx;
	fastgl::QuadPair qp;
	Matrix DE(2,4*NN);

	memset(DE.Data,0,8*NN*sizeof(Complex));
	pl = new double [2*N];
	t1 = 4./M_PI*(eb/es - eb/ep); t2 = 4./M_PI*(ep - es)/eb;

	if (kr < (0.5*kDx+1.e-12)) {DE.Data[0] = M_SQRT2*(1. - eb/ep); DE.Data[4*NN] = M_SQRT2*(ep/eb - 1.);}
	else {
		DE.Data[0] = M_SQRT2*(1. - eb/es); DE.Data[4*NN] = M_SQRT2*(es/eb - 1.);
		if (kr < 0.5*sqrt(3.)*kDx) {
			for (k=0; k<np/2; ++k) { // symmetry relative to XY plane
				qp = fastgl::GLPair(np,k+1);
				tt = 0.5*M_PI*(qp.x()+1.); ww = 0.5*M_PI*qp.weight;
				if (fabs(kr*cos(tt)) < 0.5*kDx) { // plane intersects the box
						// find delta phi range
					if (kr*sin(tt) > 0.5*kDx) { // circle-box intersecion
						if (kr*sin(tt) < M_SQRT2*0.5*kDx) df = 0.5*M_PI - 2.*acos(0.5*kDx/kr/sin(tt));
						else df = 0.;
					}
					else df = 0.5*M_PI;
					ww *= sin(tt)*df;
					pLegn(pl,tt,2*N-1);
					for (n=0; n<2*N; n+=2) {
						DE.Data[n*(n+1)] += ww*t1*pl[n]; DE.Data[4*NN+n*(n+1)] += ww*t2*pl[n];
					}
					for (m=4; m<2*N; m+=4) {
						wws = ((m/4)%2) ? -ww*sinc(0.5*m*df) : ww*sinc(0.5*m*df);
						paLegn(pl,tt,2*N-1,m);
						for (n=m; n<2*N; n+=2) {
							DE.Data[n*(n+1)+m] += wws*t1*pl[n]; DE.Data[4*NN+n*(n+1)+m] += wws*t2*pl[n];
							DE.Data[n*(n+1)-m] += wws*t1*pl[n]; DE.Data[4*NN+n*(n+1)-m] += wws*t2*pl[n];
						}
					}
				}
			}
		}
	}

	delete [] pl;
	return DE;
}

Matrix box::calc_de_p(int N, int NL, double *kr, Complex eb, Complex es) const {
	int k, i, nm, j, np = 100*N, m, n, NN = N*N;//, np = 100*N
	double xx, ww, tt, wws, *pl, df, tv;
	Complex t1, t2, tc, ep = epx;
	fastgl::QuadPair qp;
	Matrix DE(2*NL,4*NN), DDE(2*NL,4*NN);

	memset(DE.Data,0,8*NL*NN*sizeof(Complex));
	memset(DDE.Data,0,8*NL*NN*sizeof(Complex));
	pl = new double [2*N];
	t1 = 4./M_PI*(eb/es - eb/ep); t2 = 4./M_PI*(ep - es)/eb;

	for (j=0; j<NL; ++j) {
		if (kr[j] < (0.5*kDx+1.e-12)) { // sphere inside the cube
			DE.Data[j*4*NN] = M_SQRT2*(1. - eb/ep); DE.Data[(j+NL)*4*NN] = M_SQRT2*(ep/eb - 1.);
		}
		else {
			DE.Data[j*4*NN] = M_SQRT2*(1. - eb/es); DE.Data[(j+NL)*4*NN] = M_SQRT2*(es/eb - 1.);
			if (kr[j] < 0.5*sqrt(3.)*kDx) { // sphere intersects the cube
				for (k=0; k<np/2; ++k) { // loop over circles || symmetry relative to XY plane
					qp = fastgl::GLPair(np,k+1);
					tt = 0.5*M_PI*(qp.x()+1.); ww = qp.weight*0.5*M_PI;
					if ((fabs(kr[j]*cos(tt)) < 0.5*kDz) && (kr[j]*sin(tt) < M_SQRT2*0.5*kDx)) { // plane intersects the cube //
							// find delta phi range
						if (kr[j]*sin(tt) > 0.5*kDx) { // circle intersects the cube
							//if (kr[j]*sin(tt) < M_SQRT2*0.5*kDx)
								df = 0.5*M_PI - 2.*acos(0.5*kDx/kr[j]/sin(tt));
							//else df = 0.;
						}
						else df = 0.5*M_PI; // circle inside the cube
						ww *= sin(tt)*df;
						pLegn(pl,tt,2*N-1); // m==0
						for (n=0; n<2*N; n+=2) {tv = ww*pl[n]; DE.Data[j*4*NN+n*(n+1)] += t1*tv; DE.Data[(j+NL)*4*NN+n*(n+1)] += t2*tv;}
						for (m=4; m<2*N; m+=4) {
							wws = ((m/4)%2) ? -ww*sinc(0.5*m*df) : ww*sinc(0.5*m*df);
							paLegn(pl,tt,2*N-1,m);
							for (n=m; n<2*N; n+=2) {
								tv = wws*pl[n]; tc = t1*tv; DE.Data[j*4*NN+n*(n+1)+m] += tc; DE.Data[j*4*NN+n*(n+1)-m] += tc;
								tc = t2*tv; DE.Data[(j+NL)*4*NN+n*(n+1)+m] += tc; DE.Data[(j+NL)*4*NN+n*(n+1)-m] += tc;
							}
						}
					}
				}
			}
		}
	}
	delete [] pl;
	return DE;
}
