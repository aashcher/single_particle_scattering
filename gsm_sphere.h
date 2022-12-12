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

#ifndef GSM_SPHERE_H
#define GSM_SPHERE_H

#include "./matrix.h"
#include "./spfunc.h"

struct mscat_r {
};

class Method3D {
public:
	int N;

	Method3D(int N_) {N = N_;}

	Vector calc_pw(double as, double ap, double th, double ph); // plane wave decomposition
	Vector calc_pw(double as, double ap, double th, double ph, double kx0, double ky0, double kz0); // plane wave decomposition with additional phase
	Vector calc_edz(double px, double py, double pz, Complex krz, int in); // electric dipole field on axis Z
	Vector calc_edz_1(double px, double py, double pz, Complex krz, int in); // electric dipole field on axis Z, |m|<=1

	Vector calc_far(const Vector &V, double th, double ph);
	Vector calc_far_1(const Vector &V, double th, double ph);// |m|<=1

	Vector mul_Vr(const Vector &VI, const Vector &R);

	double calc_Wsca(const Vector &VS, double tC);
	double calc_Wsca(const Vector &VS, double Z, double k);
	double calc_Csca(const Vector &VS, double Z, double k, double E0);
	Vector calc_Csca_part(const Vector &VS, double k, double E0);
	double calc_Wsca_1(const Vector &VS, double tC);// |m|<=1
	double calc_Wext(const Vector &VI, const Vector &VS, double tC);
	double balance1(const Vector&, const Vector&); // Winc-Wsca
	double balance2(const Vector&, const Vector&); // Wsca-Wext

	double directivity(const Vector &VS, double th, double ph, double tC);
	double directivity_1(const Vector &VS, double th, double ph, double tC);// |m|<=1

		// multiple scattering
};

class Method3D_Mie : public Method3D {
public:
	Method3D_Mie(int N_) : Method3D(N_) {}

	Matrix calc_RTi(double kr, Complex e1, Complex e2, Complex m1, Complex m2); // j <-> h1h2
	Matrix calc_RTii(double kr, Complex e1, Complex e2, Complex m1, Complex m2); // j <-> jh1
	Matrix calc_RTo(double kr, Complex e1, Complex e2, Complex m1, Complex m2); // h1h2 <-> h1h2
	Matrix calc_RTio(double kr, Complex e1, Complex e2, Complex m1, Complex m2); // jh1 <-> h1h2

	Matrix calc_RT(double kr, Complex e1, Complex e2, Complex m1, Complex m2); // jh1 <-> jh1, S-matrix
	Matrix calc_TT(double kr, Complex e1, Complex e2, Complex m1, Complex m2); // jh1 <-> jh1, T-matrix, increment

	Vector calc_R(double kr, Complex e1, Complex e2, Complex m1, Complex m2);
	Vector calc_RR(double kr, Complex e1, Complex e2, Complex m1, Complex m2);
	Vector calc_RL(double *kr, Complex *ep, Complex *mp, int ns); // layered sphere

	Matrix calc_SML(double *kr, Complex *ep, Complex *mp, int ns); // layered sphere, S-matrix algorithm
	Matrix calc_TML(double *kr, Complex *ep, Complex *mp, int ns); // layered sphere, T-matrix algorithm
	Matrix calc_SML(Matrix **ML, int ns); // layered sphere, S-matrix algorithm
	Matrix calc_TML(Matrix **ML, int ns); // layered sphere, T-matrix algorithm
};

class Method3D_SM : public Method3D {
public:
	int NL;
	Complex eb;
	Matrix **mc;//, *mc1, *mc0;

	Method3D_SM(int N_);// : Method3D(N_) {mc0 = new Matrix(2*N,N*N);}
	~Method3D_SM() {for (int i=0; i<N*N; ++i) delete mc[i]; delete [] mc;}

	void calc_mc(void);

	Matrix calc_Qt(const Matrix &DE);
	Matrix calc_Qp(const Matrix &DE);
	Matrix calc_Vr(double kR);
	Matrix calc_Vt(double kR);

	Vector calc_Rin_r(double kr, Complex ec);
	Matrix calc_Rin_t(double kr, Complex ec);
	void calc_RTout_r(double kr, Complex es, Vector &R00, Vector &R11, Vector &T01, Vector &T10);
	void calc_RTout_t(double kr, Complex es, Matrix &R00, Matrix &R11, Matrix &T01, Matrix &T10);
	void calc_RT_r(const Matrix &V, const Matrix &Q, const Complex wt, Vector &R00, Vector &R11, Vector &T01, Vector &T10);
	void calc_RT_t(const Matrix &V, const Matrix &Q, const Complex wt, Matrix &R00, Matrix &R11, Matrix &T01, Matrix &T10);
	void calc_RT_p(const Matrix &V, const Matrix &Q, const Complex wt, Matrix &R00, Matrix &R11, Matrix &T01, Matrix &T10);

	Vector mul_Sr(Vector &R, Vector &R00, Vector &R11, Vector &T01, Vector &T10);
	Matrix mul_St(Matrix &R, Matrix &R00, Matrix &R11, Matrix &T01, Matrix &T10);

	template <class SType> Vector calc_Sr(SType &ST, double kr1, double kr2, Complex ec, Complex eb_, Complex es, int NL_);
	template <class SType> Matrix calc_St(SType &ST, double kr1, double kr2, Complex ec, Complex eb_, Complex es, int NL_);
	template <class SType> Matrix calc_Sp(SType &ST, double kr1, double kr2, Complex ec, Complex eb_, Complex es, int NL_);
};

Vector mul_ASr(const Vector &V, void *P);
Vector mul_ASt(const Vector &V, void *P);
Vector mul_ASp(const Vector &V, void *P);

class Method3D_GSMS : public Method3D {
public:
	int NL;
	double eb, *kR;
	Complex *cce, *cch, *wR;
	Matrix *v,*qm;

	Method3D_GSMS(int N_) : Method3D(N_) {cce = new Complex [8*N]; cch = new Complex [8*N];}
	~Method3D_GSMS() {delete [] cce; delete [] cch;}

	Matrix calc_Qt(const Matrix DE);
	Matrix calc_Qp(const Matrix DE);
	Matrix calc_Vr(void);
	Matrix calc_Vt(void);

	void calc_RT(double kr1, double kr2, Complex ec, Complex es);
	void calc_Rin(double kr, Complex ec, Vector &R11);
	void calc_RTout(double kr, Complex es, Vector &R00, Vector &R11, Vector &T01, Vector &T10);

	Vector calc_inc(const Vector &VA);
	Vector calc_rad_r(const Vector &VA);
	Vector calc_rad_t(const Vector &VA);
	Vector calc_rad_p(const Vector &VA);
	Vector calc_dif(const Vector &VA);
	Vector calc_out(const Vector &VA);
	void add_inc(const Vector &VA, Vector &VB);

	template <class SType> Vector calc_SVr(const SType &ST, const Vector &VA, double kr1, double kr2, Complex ec, double eb_, Complex es, int NL_);
	template <class SType> Vector calc_SVt(const SType &ST, const Vector &VA, double kr1, double kr2, Complex ec, double eb_, Complex es, int NL_);
	template <class SType> Vector calc_SVt_in(const SType &ST, const Vector &VA, double kr1, double kr2, Complex ec, double eb_, Complex es, int NL_);
	template <class SType> Vector calc_SVt(const SType &ST, const Vector &VA, const Vector &W, double kr1, double kr2, Complex ec, double eb_, Complex es, int NL_);
	template <class SType> Vector calc_SVp(const SType &ST, const Vector &VA, double kr1, double kr2, Complex ec, double eb_, Complex es, int NL_);
};

Vector mul_ASct(const Vector &V, void *P);

class Method3D_GSMSC : public Method3D {
public:
	int NL;
	double eb, mb, *kR;
	Complex *cce, *cch, *wR, *cfe, *cfm;
	Matrix **v_reg, **v_out, **R;

	Method3D_GSMSC(int N_) : Method3D(N_) {cce = new Complex [8*N]; cch = new Complex [8*N];}
	~Method3D_GSMSC() {delete [] cce; delete [] cch;}

	//Matrix calc_Qt(const Matrix DE);
	void calc_Qt(const Matrix DE);
	void calc_Vt(void);

	void calc_RT(double kr1, double kr2, Complex ec, Complex es);
	void calc_Rin(double kr, Complex ec, Vector &R11);
	void calc_RTout(double kr, Complex es, Vector &R00, Vector &R11, Vector &T01, Vector &T10);

	Vector calc_inc(const Vector &VA);
	Vector calc_rad_t(const Vector &VA);
	Vector calc_dif(const Vector &VA);
	Vector calc_out(const Vector &VA);
	void add_inc(const Vector &VA, Vector &VB);

	template <class SType> Vector calc_SVt(const SType &ST, const Vector &VA,
		double kr1, double kr2, Complex ec, double eb_, Complex es, int NL_);
};

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

template <class SType> Vector Method3D_SM::calc_Sr(SType &ST, double kr1, double kr2, Complex ec, Complex eb_, Complex es, int NL_) {
	NL = NL_; eb = eb_;
	double kR, kdR, nb = sqrt(eb.real());
	Vector T(2*N), R00(2*N), R11(2*N), T01(2*N), T10(2*N);
	Matrix DE(2,1);

	kdR = nb*(kr2-kr1)/double(NL); kR = nb*kr1 + 0.5*kdR;
	T = calc_Rin_r(kr1,ec);
	for (int k=0; k<NL; ++k,kR+=kdR) {
		DE = ST.calc_de_r(kR/nb,eb,es); //if (k==NL/2) cout<<"de: "<<DE(0,0)<<" "<<DE(1,0)<<endl;
		calc_RT_r(calc_Vr(kR),DE,j_*kR*kR*kdR,R00,R11,T01,T10);
		//if (k==NL/2) for (int n=0; n<N; ++n) cout<<n<<"\t"<<R00(n+N)<<"   "<<R11(n+N)<<endl;
		T = mul_Sr(T,R00,R11,T01,T10);
	}
	calc_RTout_r(kr2,es,R00,R11,T01,T10);
	T = mul_Sr(T,R00,R11,T01,T10);

	return T;
}

template <class SType> Matrix Method3D_SM::calc_St(SType &ST, double kr1, double kr2, Complex ec, Complex eb_, Complex es, int NL_) {
	int NN = N*N, NN2 = 2*NN;
	double kR, kdR, nb;
	NL = NL_; eb = eb_; nb = sqrt(eb.real());
	Matrix T(NN2), R00(NN2), R11(NN2), T01(NN2), T10(NN2), DE(2,2*N), Q(3,NN*N);//, V(6,N);

	kdR = nb*(kr2-kr1)/double(NL); kR = nb*kr1 + 0.5*kdR;
	T = calc_Rin_t(kr1,ec);
	for (int k=0; k<NL; ++k,kR+=kdR) {
		cout<<"layer "<<k<<endl;
		DE = ST.calc_de_t(N,kR/nb,eb,es);
		//cout<<"  D: "<<DE.norm()<<endl;
		Q = calc_Qt(DE);
		//cout<<"  Q: "<<Q.norm()<<endl;
		//V = calc_Vt(kR);
		//cout<<"  V: "<<V.norm()<<endl;
		calc_RT_t(calc_Vt(kR),Q,j_*kR*kR*kdR,R00,R11,T01,T10);
		T = mul_St(T,R00,R11,T01,T10);
	}
	calc_RTout_t(kr2,es,R00,R11,T01,T10);
	T = mul_St(T,R00,R11,T01,T10);

	return T;
}

template <class SType> Matrix Method3D_SM::calc_Sp(SType &ST, double kr1, double kr2, Complex ec, Complex eb_, Complex es, int NL_) {
	int NN = N*N, NN2 = 2*NN;
	double kR, kdR, nb, tm;
	NL = NL_; eb = eb_; nb = sqrt(eb.real());
	Matrix T(NN2), R00(NN2), R11(NN2), T01(NN2), T10(NN2), DE(2,4*NN), Q(3,NN*NN);

	calc_mc();
	kdR = nb*(kr2-kr1)/double(NL); kR = nb*kr1 + 0.5*kdR;
	T = calc_Rin_t(kr1,ec);
	for (int k=0; k<NL; ++k,kR+=kdR) {
		cout<<"layer "<<k<<endl;
		DE = ST.calc_de_p(N,kR/nb,eb,es);
		//cout<<"  de ok\n";
		Q = calc_Qp(DE);
		//cout<<"  q ok\n";
		calc_RT_p(calc_Vt(kR),Q,j_*kR*kR*kdR,R00,R11,T01,T10);
		//cout<<"  rt ok\n";
		T = mul_St(T,R00,R11,T01,T10);
		//cout<<"  mul ok\n";
	}
	calc_RTout_t(kr2,es,R00,R11,T01,T10);
	T = mul_St(T,R00,R11,T01,T10);

	return T;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

template <class SType> Vector Method3D_GSMS::calc_SVr(const SType &ST, const Vector &VA, double kr1, double kr2,
																											Complex ec, double eb_, Complex es, int NL_) {
	int k, NN = N*N; double kdr, nb;
	NL = NL_; nb = sqrt(eb = eb_);
	Vector V(4*NL*NN), VB(2*NN);

	kR = new double [NL]; wR = new Complex [NL];
	kdr = (kr2-kr1)/double(NL); kR[0] = kr1+0.5*kdr;
	for (k=1; k<NL; ++k) kR[k] = kR[k-1] + kdr;
	qm = new Matrix (2,NL); *qm = ST.calc_de_r(NL,kR,eb,es);
	//cout<<qm->norm(); cin.get();
	
	kdr *= nb; for (k=0; k<NL; ++k) {kR[k] *= nb; wR[k] = j_*kR[k]*kR[k]*kdr;}
	v = new Matrix (9*NL,N); *v = calc_Vr();

	calc_RT(kr1,kr2,ec,es);
	V = calc_inc(VA);
	V.GMResAG(reinterpret_cast<void*>(this),&mul_ASr,1000,1.e-5);
	VB = calc_out(calc_rad_r(V));
//	cout<<VB.norm(); cin.get();
	add_inc(VA,VB);

	delete v; delete qm; delete [] kR; delete [] wR;
	return VB;
}

template <class SType> Vector Method3D_GSMS::calc_SVt(const SType &ST, const Vector &VA, double kr1, double kr2,
																											Complex ec, double eb_, Complex es, int NL_) {
	int k, NN = N*N; double kdr, nb;
	NL = NL_; nb = sqrt(eb = eb_);
	Vector V(4*NL*NN), V0(4*NL*NN), VB(2*NN);

	memset(V0.Data,0,4*NL*NN*sizeof(Complex));

	kR = new double [NL]; wR = new Complex [NL];
	kdr = (kr2-kr1)/double(NL); kR[0] = kr1+0.5*kdr;
	for (k=1; k<NL; ++k) kR[k] = kR[k-1] + kdr;
	qm = new Matrix(2*NL,NN*N);
	*qm = calc_Qt(ST.calc_de_t(N,NL,kR,eb,es));

	kdr *= nb; for (k=0; k<NL; ++k) {kR[k] *= nb; wR[k] = j_*kR[k]*kR[k]*kdr;}
	v = new Matrix(6*NL,N); *v = calc_Vt();

	calc_RT(kr1,kr2,ec,es);
	V = calc_inc(VA);
	V.GMResAG(reinterpret_cast<void*>(this),&mul_ASt,2000,1.e-5);
	//V.GMResM(reinterpret_cast<void*>(this),&mul_ASt,V0,20,1.e-8);
	//V.BiCGs(reinterpret_cast<void*>(this),&mul_ASt,1000,1.e-10);
	VB = calc_out(calc_rad_t(V)); add_inc(VA,VB);

	delete v; delete qm; delete [] kR; delete [] wR;
	return VB;
}

template <class SType> Vector Method3D_GSMS::calc_SVt_in(const SType &ST, const Vector &VA, double kr1, double kr2, Complex ec, double eb_, Complex es, int NL_) {
	int k, NN = N*N; double kdr, nb;
	NL = NL_; nb = sqrt(eb = eb_);
	Vector V(4*NL*NN), V0(4*NL*NN);

	memset(V0.Data,0,4*NL*NN*sizeof(Complex));

	kR = new double [NL]; wR = new Complex [NL];
	kdr = (kr2-kr1)/double(NL); kR[0] = kr1+0.5*kdr;
	for (k=1; k<NL; ++k) kR[k] = kR[k-1] + kdr;
	qm = new Matrix(2*NL,NN*N);
	*qm = calc_Qt(ST.calc_de_t(N,NL,kR,eb,es));

	kdr *= nb; for (k=0; k<NL; ++k) {kR[k] *= nb; wR[k] = j_*kR[k]*kR[k]*kdr;}
	v = new Matrix(6*NL,N); *v = calc_Vt();

	calc_RT(kr1,kr2,ec,es);
	V = calc_inc(VA);
	//V.GMResAG(reinterpret_cast<void*>(this),&mul_ASt,1000,1.e-8);
	//V.GMResAG(reinterpret_cast<void*>(this),&mul_ASt,V0,1000,1.e-8);
	V.BiCGs(reinterpret_cast<void*>(this),&mul_ASt,1000,1.e-10);

	delete v; delete qm; delete [] kR; delete [] wR;
	return V;
}

template <class SType> Vector Method3D_GSMS::calc_SVt(const SType &ST, const Vector &VA, const Vector &W, double kr1, double kr2, Complex ec, double eb_, Complex es, int NL_) {
	int k, NN = N*N; double kdr, nb;
	NL = NL_; nb = sqrt(eb = eb_);
	Vector V(4*NL*NN), VB(2*NN);

	kR = new double [NL]; wR = new Complex [NL];
	kdr = (kr2-kr1)/double(NL); kR[0] = kr1+0.5*kdr;
	for (k=1; k<NL; ++k) kR[k] = kR[k-1] + kdr;
	qm = new Matrix(2*NL,NN*N);
	*qm = calc_Qt(ST.calc_de_t(N,NL,kR,eb,es));

	kdr *= nb; for (k=0; k<NL; ++k) {kR[k] *= nb; wR[k] = j_*kR[k]*kR[k]*kdr;}
	v = new Matrix(6*NL,N); *v = calc_Vt();

	calc_RT(kr1,kr2,ec,es);
	V = calc_inc(VA);
	//V.GMResAG(reinterpret_cast<void*>(this),&mul_ASt,W,1000,1.e-8);
	V.BiCGs(reinterpret_cast<void*>(this),&mul_ASt,W,1000,1.e-10);
	VB = calc_out(calc_rad_t(V)); add_inc(VA,VB);

	delete v; delete qm; delete [] kR; delete [] wR;
	return VB;
}


template <class SType> Vector Method3D_GSMS::calc_SVp(const SType &ST, const Vector &VA, double kr1, double kr2,
																											Complex ec, double eb_, Complex es, int NL_) {
	int k, NN = N*N; double kdr, nb;
	NL = NL_; nb = sqrt(eb = eb_);
	Vector V(4*NL*NN), VB(2*NN);

	kR = new double [NL]; wR = new Complex [NL];
	kdr = (kr2-kr1)/double(NL); kR[0] = kr1+0.5*kdr;
	for (k=1; k<NL; ++k) kR[k] = kR[k-1] + kdr;
	qm = new Matrix(2*NL,NN*NN);
	*qm = calc_Qp(ST.calc_de_p(N,NL,kR,eb,es));

	kdr *= nb; for (k=0; k<NL; ++k) {kR[k] *= nb; wR[k] = j_*kR[k]*kR[k]*kdr;}
	v = new Matrix(6*NL,N); *v = calc_Vt();

	calc_RT(kr1,kr2,ec,es);
	V = calc_inc(VA);
	V.GMResAG(reinterpret_cast<void*>(this),&mul_ASp,1000,1.e-5);
	VB = calc_out(calc_rad_p(V)); add_inc(VA,VB);

	delete v; delete qm; delete [] kR; delete [] wR;
	return VB;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

template <class SType> Vector Method3D_GSMSC::calc_SVt(const SType &ST, const Vector &VA, double kr1_, double kr2_,
				Complex ec, double eb_, Complex es, int NL_) {
	int k, NN = N*N;
	double kdr, rib, kr0, kr1, kr2;
	NL = NL_; rib = sqrt(eb = eb_);
	Vector V(4*NL*NN), V0(4*NL*NN), VB(2*NN);

	memset(V0.Data,0,4*NL*NN*sizeof(Complex));
	kR = new double [NL]; wR = new Complex [NL];
	cfe = new Complex [NL]; cfm = new Complex [NL];

	kr1 = kr1_; kr2 = kr2_; // in vacuum

	kr0 = 0.5*(ST.R1 + ST.R2); // in vacuum
	kdr = (kr2-kr1)/double(NL); // in vacuum
	kR[0] = kr1 + 0.5*kdr; cfe[0] = ST.ep/eb; cfm[0] = 1.;
	k = 1;
	while (((kR[k-1]+kdr) < kr0) && (k < NL)) {
		kR[k] = kR[k-1] + kdr; cfe[k] = cfe[k-1];
		cfm[k++] = 1.;
	}
	while (k < NL) {kR[k] = kR[k-1] + kdr; cfe[k] = es/eb; cfm[k++] = 1.;}

	kdr *= rib; for (k=0; k<NL; ++k) {kR[k] *= rib; wR[k] = 0.5*kR[k]*kR[k]*kdr;} // in basis medium

	R = new Matrix*[8];
	for (int i=0; i<8; ++i) R[i] = new Matrix(NL,NN*N);
	calc_Qt(ST.calc_ge_t(N,NL,kR,kr1,kr2,kr0,rib,es));

	//v = new Matrix(6*NL,N);
	//*v = calc_Vt(); // input kR - in basis medium
	v_reg = new Matrix*[3]; v_out = new Matrix*[3];
	for (int i=0; i<3; ++i) {v_reg[i] = new Matrix(NL,N); v_out[i] = new Matrix(NL,N);}
	calc_Vt();

	calc_RT(kr1,kr2,ec,es); // input kr1,2 - in vacuum
	V = calc_inc(VA);
	V.GMResAG(reinterpret_cast<void*>(this),&mul_ASct,2000,1.e-10);
	VB = calc_out(calc_rad_t(V)); add_inc(VA,VB);

	for (int i=0; i<8; ++i) delete R[i]; delete []R;
	for (int i=0; i<3; ++i) {delete v_reg[i]; delete v_out[i];} delete []v_reg; delete []v_out;
	delete [] kR; delete [] wR;
	delete [] cfe; delete [] cfm;
	return VB;
}

#endif
