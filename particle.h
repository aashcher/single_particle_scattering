
#ifndef PARTICLE_H
#define PARTICLE_H

#include <algorithm>
#include "./matrix.h"

struct particle3D {
public:
	double R1, R2;
};

struct sphere : public particle3D {
public:
	double kR;
	Complex ep;

	sphere() {kR = R1 = R2 = 0.; ep = 1.;}
	sphere(double kR_, Complex ep_) {kR = R2 = kR_; R1 = 0.; ep = ep_;}
	sphere& operator = (const sphere &S) {kR = R2 = S.kR; ep = S.ep;}

	Matrix calc_de_r(double kr, Complex eb, Complex es) const;
	Matrix calc_de_r(int NL, double *kr, Complex eb, Complex es) const;
	Matrix calc_de3(int nx, int ny, int nz, double kdx, double kdy, double kdz, Complex eb) const;
};

struct sphere_ml : public particle3D {
public:
	int ns;
	double *kR;
	Complex *ep;

	sphere_ml() {ns = 0;}
	sphere_ml(int ns_, double *kR_, Complex *ep_) {
		if (ns = ns_) {
			kR = new double [ns]; ep = new Complex [ns+1];
			for (int k=0; k<ns; ++k) {kR[k] = kR_[k]; ep[k] = ep_[k];} ep[ns] = ep_[ns];
			R1 = kR[0]; R2 = kR[ns-1];
		}
	}
	~sphere_ml() {if (ns) {delete [] kR; delete [] ep;}}
	sphere_ml& operator = (const sphere_ml &S) {
		if ((ns != S.ns)&&(ns)) {delete [] kR; delete [] ep;}
		if (ns = S.ns) {
			kR = new double [ns]; ep = new Complex [ns];
			for (int k=0; k<ns; ++k) {kR[k] = S.kR[k]; ep[k] = S.ep[k];}
			R1 = kR[0]; R2 = kR[ns-1];
		}
	}

	Matrix calc_de_r(int NL, double *kr, Complex eb, Complex es) const;
	Matrix calc_de3(Complex eb) const;
};

struct sphere_rad : public particle3D {
public:
	int nl;
	double *kr;
	Complex ec;
};

struct ellipsoid_ab : public particle3D {
public:
	double kx, kz;
	Complex ep;

	ellipsoid_ab() {kx = kz = R1 = R2 = 0.; ep = 1.;}
	ellipsoid_ab(double kx_, double kz_, Complex ep_) {
		kx = kx_; kz = kz_; ep = ep_;
		R1 = min(kx,kz); R2 = max(kx,kz);
	}
	ellipsoid_ab& operator = (const ellipsoid_ab &P) {kx = P.kx; kz = P.kz; ep = P.ep; R1 = P.R1; R2 = P.R2;}

	Matrix calc_de_t(int N, double kr, Complex eb, Complex es) const;
	Matrix calc_de_p(int N, double kr, Complex eb, Complex es) const;
	Matrix calc_de_t(int N, int NL, double *kr, Complex eb, Complex es) const;
	Matrix calc_de_p(int N, int NL, double *kr, Complex eb, Complex es) const;

	//Matrix calc_ge_t(int N, double kr, double kr1, double kr2, Complex eb, Complex es) const;
	Matrix calc_ge_t(int N, int NL, double *kr, double kr1, double kr2, double kr0, Complex rib, Complex es) const;

	Matrix calc_de3(int nx, int ny, int nz, double kdx, double kdy, double kdz, Complex eb) const;

private:
	void get_fg(double t, double *p, double *fg) const;
	void get_fgi(int n, double *p, double *ff, int Ni) const;
};

struct cassini : public particle3D {
public:
	double ka, kc;
	Complex ep;

	cassini() {ka = kc = R1 = R2 = 0.; ep = 1.;}
	cassini(double ka_, double kc_, Complex ep_) {
		ka = ka_; kc = kc_; ep = ep_;
		if (ka > (kc-1.e-5)) {cout<<"incorrect cassini oval parameters\n"; return;}
		R1 = sqrt(kc*kc - ka*ka); R2 = sqrt(kc*kc + ka*ka);
	}
	cassini& operator = (const cassini &P) {ka = P.ka; kc = P.kc; ep = P.ep; R1 = P.R1; R2 = P.R2;}

	Matrix calc_de_t(int N, double kr, Complex eb, Complex es) const;
	Matrix calc_de_p(int N, double kr, Complex eb, Complex es) const;
	Matrix calc_de_t(int N, int NL, double *kr, Complex eb, Complex es) const;
	Matrix calc_de_p(int N, int NL, double *kr, Complex eb, Complex es) const;

	Matrix calc_de3(int nx, int ny, int nz, double kdx, double kdy, double kdz, Complex eb) const;
};

struct cylinder : public particle3D {
public:
	double kR, kH;
	Complex ep;

	cylinder() {kR = kH = R1 = R2 = 0.; ep = 1.;}
	cylinder(double kR_, double kH_, Complex ep_) {
		kR = kR_; kH = kH_; ep = ep_;
		R1 = min(kR,0.5*kH); R2 = sqrt(kR*kR + 0.25*kH*kH);
	}
	cylinder& operator = (const cylinder &P) {kR = P.kR; kH = P.kH; ep = P.ep; R1 = P.R1; R2 = P.R2;}

	Matrix calc_de_t(int N, double kr, Complex eb, Complex es) const;
	Matrix calc_de_p(int N, double kr, Complex eb, Complex es) const;
	Matrix calc_de_t(int N, int NL, double *kr, Complex eb, Complex es) const;
	Matrix calc_de_p(int N, int NL, double *kr, Complex eb, Complex es) const;
	Matrix calc_de3(int nx, int ny, int nz, double kdx, double kdy, double kdz, Complex eb) const;
};

struct box : public particle3D {
public:
	double kDx, kDy, kDz;
	Complex epx, epy, epz;

	box() {kDx = kDy = kDz = R1 = R2 = 0.; epx = epy = epz = 1.;}
	box(double kDx_, double kDy_, double kDz_, Complex epx_, Complex epy_, Complex epz_) {
		kDx = kDx_; kDy = kDy_; kDz = kDz_; epx = epx_; epy = epy_; epz = epz_;
		R1 = 0.5*(min(kDx,min(kDy,kDz))); R2 = 0.5*sqrt(kDx*kDx + kDy*kDy + kDz*kDz);
	}
	box& operator = (const box &B) {
		kDx = B.kDx; kDy = B.kDy; kDz = B.kDz; epx = B.epx; epy = B.epy; epz = B.epz;
		R1 = 0.5*(min(kDx, min(kDy, kDz))); R2 = 0.5*sqrt(kDx*kDx + kDy*kDy + kDz*kDz);
	}

	Matrix calc_de3(int nx, int ny, int nz, double kdx, double kdy, double kdz, Complex eb) const;
	Matrix calc_de_p(int N, double kr, Complex eb, Complex es) const;
	Matrix calc_de_p(int N, int NL, double *kr, Complex eb, Complex es) const;
};

#endif