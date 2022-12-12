
#include "./gsm_sphere.h"
#include <algorithm>

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		// Method3D //

Vector Method3D::calc_pw(double as, double ap, double th, double ph) {
	int n, m; double tv, tp, tt; Complex tc = j_, tcc, te; Vector VA(2*N*N);
	memset(VA.Data,0,2*N*N*sizeof(Complex));
	for (n=1; n<N; ++n) {
		tcc = tc*4.*sqrt(M_PI)/sqrt(2.*n*(n+1.));
		for (m=-n; m<n+1; ++m) {
			te = exp(-j_*double(m)*ph);
			tp = LPin(th,n,m); tt = LTaun(th,n,m);
			VA.Data[n*(n+1)+m] = -tcc*(j_*tp*ap + as*tt)*te;
			VA.Data[N*N+n*(n+1)+m] = -tcc*(j_*tt*ap + as*tp)*te;
		}
		tc *= j_;
	}
	return VA;
}

Vector Method3D::calc_pw(double as, double ap, double th, double ph, double kx0, double ky0, double kz0) {
	int n, m; double tv, tp, tt;
	Complex tc = j_, tcc, te; Vector VA(2*N*N);
	memset(VA.Data,0,2*N*N*sizeof(Complex));
	for (n=1; n<N; ++n) {
		tcc = tc*4.*sqrt(M_PI)/sqrt(2.*n*(n+1.));
		for (m=-n; m<n+1; ++m) {
			te = exp(-j_*double(m)*ph);
			tp = LPin(th,n,m); tt = LTaun(th,n,m);
			VA.Data[n*(n+1)+m] = -tcc*(j_*tp*ap + as*tt)*te;
			VA.Data[N*N+n*(n+1)+m] = -tcc*(j_*tt*ap + as*tp)*te;
		}
		tc *= j_;
	}
	te = exp(j_*(kx0*sin(th)*cos(ph) + ky0*sin(th)*sin(ph) + kz0*cos(th)));
	for (n=0; n<2*N*N; ++n) VA.Data[n] *= te;
	return VA;
}

Vector Method3D::calc_edz(double px, double py, double pz, Complex krz, int in) {
	int n, m;
	double tv = -0.25/sqrt(M_PI), tvn, tv1, tv2;
	Complex pp = Complex(px,py), pm = conj(pp), zf, zfd;
	Vector VA(2*N*N);
	memset(VA.Data,0,2*N*N*sizeof(Complex));

	if (in == 1) { // field inside dipole radius
		for (n=1; n<N; ++n) {
			tvn = tv*sqrt(2*n+1.);
			zf = besh1(krz,n); zfd = besh1d(krz,n);
			VA.Data[n*(n+1)-1] = pm*( VA.Data[n*(n+1)+1] = tvn*zf );
			VA.Data[n*(n+1)+1] *= pp;
			VA.Data[N*N+n*(n+1)+0] = -2.*j_*tvn*pz*sqrt(n*(n+1.))*zf/krz;
			VA.Data[N*N+n*(n+1)+1] = -pp*( VA.Data[N*N+n*(n+1)-1] = j_*tvn*(zfd + zf/krz) );
			VA.Data[N*N+n*(n+1)-1] *= pm;
			tv = -tv;
		}
	}
	else { // field outside dipole radius
		if (abs(krz) > 1.e-10) for (n=1; n<N; ++n) {
			tvn = tv*sqrt(2*n+1.);
			zf = besj(krz,n); zfd = besjd(krz,n);
			VA.Data[n*(n+1)-1] = pm*( VA.Data[n*(n+1)+1] = tvn*zf );
			VA.Data[n*(n+1)+1] *= pp;
			VA.Data[N*N+n*(n+1)+0] = -2.*j_*tvn*pz*sqrt(n*(n+1.))*zf/krz;
			VA.Data[N*N+n*(n+1)+1] = -pp*( VA.Data[N*N+n*(n+1)-1] = j_*tvn*(zfd + zf/krz) );
			VA.Data[N*N+n*(n+1)-1] *= pm;
			tv = -tv;
		}
		else {
			zf = 0.1*j_/sqrt(3.*M_PI);
			VA.Data[N*N+2] = zf*M_SQRT2*pz;
			VA.Data[N*N+3] = zf*pp;
			VA.Data[N*N+1] = -zf*pm;
		}
	}

	return VA;
}

Vector Method3D::calc_edz_1(double px, double py, double pz, Complex krz, int in) {
	int n, m;
	double tv = -0.25/sqrt(M_PI), tvn, tv1, tv2;
	Complex pp = Complex(px,py), pm = conj(pp), zf, zfd;
	Vector VA(5*N);
	memset(VA.Data,0,5*N*sizeof(Complex));

	if (in == 1) { // field inside dipole radius
		for (n=1; n<N; ++n) {
			tvn = tv*sqrt(2*n+1.);
			zf = besh1(krz,n); zfd = besh1d(krz,n);
			VA.Data[2*n] = pm*( VA.Data[2*n+1] = tvn*zf ); // m = -1,1
			VA.Data[2*n+1] *= pp; // m = 1
			VA.Data[2*N+3*n+1] = -2.*j_*tvn*pz*sqrt(n*(n+1.))*zf/krz; // m = 0
			VA.Data[2*N+3*n+2] = -pp*( VA.Data[2*N+3*n] = j_*tvn*(zfd + zf/krz) ); // m = 1,-1
			VA.Data[2*N+2*N+3*n] *= pm; // m = -1
			tv = -tv;
		}
	}
	else { // field outside dipole radius
		for (n=1; n<N; ++n) {
			tvn = tv*sqrt(2*n+1.);
			zf = besj(krz,n); zfd = besjd(krz,n);
			VA.Data[2*n] = pm*( VA.Data[2*n+1] = tvn*zf ); // m = -1,1
			VA.Data[2*n+1] *= pp;
			VA.Data[2*N+3*n+1] = -2.*j_*tvn*pz*sqrt(n*(n+1.))*zf/krz;
			VA.Data[2*N+3*n+2] = -pp*( VA.Data[2*N+3*n] = j_*tvn*(zfd + zf/krz) );
			VA.Data[2*N+3*n] *= pm;
			tv = -tv;
		}
	}

	return VA;
}

Vector Method3D::calc_far(const Vector &V, double th, double ph) {
	int m, n, NN = N*N; double tv;
	Complex tc, tc1, tc2, tc3, tc4, *te;
	Vector VE(2);
	te = new Complex [N]; for (m=0; m<N; ++m) te[m] = exp(j_*double(m)*ph);
	tc = -j_; VE.Data[0] = VE.Data[1] = 0.;
	for (n=1; n<N; ++n) {
		tc1 = V(n*(n+1))*LPin(th,n,0); tc2 = V(NN+n*(n+1))*LTaun(th,n,0);
		tc3 = V(n*(n+1))*LTaun(th,n,0); tc4 = V(NN+n*(n+1))*LPin(th,n,0);
		for (m=1; m<n+1; ++m) {
			tc1 += V(n*(n+1)-m)*LPin(th,n,-m)*conj(te[m]) + V(n*(n+1)+m)*LPin(th,n,m)*te[m]; // ae*pi
			tc2 += V(NN+n*(n+1)-m)*LTaun(th,n,-m)*conj(te[m]) + V(NN+n*(n+1)+m)*LTaun(th,n,m)*te[m]; // ah*tau
			tc3 += V(n*(n+1)-m)*LTaun(th,n,-m)*conj(te[m]) + V(n*(n+1)+m)*LTaun(th,n,m)*te[m]; // ae*tau
			tc4 += V(NN+n*(n+1)-m)*LPin(th,n,-m)*conj(te[m]) + V(NN+n*(n+1)+m)*LPin(th,n,m)*te[m]; // ah*pi
		}
		tv = 1./sqrt(n*(n+1.));
		VE.Data[0] -= tc*tv*(tc1 + tc2);
		tc *= -j_; VE.Data[1] += tc*tv*(tc3 + tc4);
	}
	VE.Data[0] *= M_SQRT1_2PI; VE.Data[1] *= M_SQRT1_2PI;
	delete [] te;
	return VE;
}

Vector Method3D::calc_far_1(const Vector &V, double th, double ph) {
	int n; double tv;
	Complex tc, tc1, tc2, tc3, tc4, te;
	Vector VE(2);
	te = exp(j_*ph);
	tc = -j_; VE.Data[0] = VE.Data[1] = 0.;
	for (n=1; n<N; ++n) {
		tc1 = V(2*n)*LPin(th,n,-1)*conj(te) + V(2*n+1)*LPin(th,n,1)*te; // ae*pi
		tc3 = V(2*n)*LTaun(th,n,-1)*conj(te) + V(2*n+1)*LTaun(th,n,1)*te; // ae*tau
		tc2 = V(2*N+3*n+1)*LTaun(th,n,0) + V(2*N+3*n)*LTaun(th,n,-1)*conj(te) + V(2*N+3*n+2)*LTaun(th,n,1)*te; // ah*tau
		tc4 = V(2*N+3*n+1)*LPin(th,n,0) + V(2*N+3*n)*LPin(th,n,-1)*conj(te) + V(2*N+3*n+2)*LPin(th,n,1)*te; // ah*pi
		
		tv = 1./sqrt(n*(n+1.));
		VE.Data[0] -= tc*tv*(tc1 + tc2);
		tc *= -j_; VE.Data[1] += tc*tv*(tc3 + tc4);
	}
	VE.Data[0] *= M_SQRT1_2PI; VE.Data[1] *= M_SQRT1_2PI;
	return VE;
}

double Method3D::calc_Wsca(const Vector &VS, double tC) {
	int n, NN = N*N; double tv = 0.;
	for (n=1; n<2*NN; ++n) tv += abs(VS(n)*VS(n));// + abs(VS(NN+n)*VS(NN+n));
	return 0.5*tv/tC;
}

double Method3D::calc_Wsca(const Vector &VS, double Z, double k) {
	int n, NN = N*N; double tv = 0.;
	for (n=1; n<2*NN; ++n) tv += abs(VS(n)*VS(n));
	return 0.5*tv/Z/k/k;
}

double Method3D::calc_Csca(const Vector &VS, double Z, double k, double E0) {
	return 2*Z/E0/E0*calc_Wsca(VS,Z,k);
}

Vector Method3D::calc_Csca_part(const Vector &VS, double k, double E0) {
	double tv = 1./k/k/E0/E0;
	Vector VC(VS);
	for (int n=1; n<2*N*N; ++n) VC.Data[n] = tv*abs(VS(n)*VS(n));
	return VC;
}

double Method3D::calc_Wsca_1(const Vector &VS, double tC) {
	double tv = 0.;
	for (int n=0; n<5*N; ++n) tv += abs(VS(n)*VS(n));
	return 0.5*tv/tC;
}

double Method3D::calc_Wext(const Vector &VI, const Vector &VS, double tC) {
	int n, NN = N*N; double tv = 0.;
	for (n=1; n<NN; ++n) tv += (VS(n)*conj(VI(n))).real() + (VS(NN+n)*conj(VI(NN+n))).real();
	return -tv/tC;
}

double Method3D::directivity(const Vector &VS, double th, double ph, double tC) {
	int n, m, nm, NN = N*N;
	double tp, tt;
	Complex tc1, tc2, tc3, tc4, tc = -j_, tcc, tce;
	tc1 = tc2 = tc3 = tc4 = 0.;
	for (n=nm=1; n<N; ++n) {
		tcc = tc/sqrt(double(n*(n+1)));
		for (m=-n; m<n+1; m++,nm++) {
			tp = LPin(th,n,m); tt = LTaun(th,n,m);
			tce = exp(j_*double(m)*ph);
			tc1 += tcc*tce*(VS(nm)*tp + VS(nm+NN)*tt);
			tc2 += tcc*tce*(VS(nm)*tt + VS(nm+NN)*tp);
		}
		tc *= -j_;
	}
	return (tc1*conj(tc1) + tc2*conj(tc2)).real()/calc_Wsca(VS,tC)/tC;
}

double Method3D::directivity_1(const Vector &VS, double th, double ph, double tC) {
	int n, m, nm, NN = N*N;
	double tp1, tt1, tp0, tt0;
	Complex tc1, tc2, tc3, tc4, tc = -j_, tcc, tce;
	tc1 = tc2 = tc3 = tc4 = 0.;
	for (n=1; n<N; ++n) {
		tcc = tc/sqrt(double(n*(n+1)));
		tp0 = LPin(th,n,0); tt0 = LTaun(th,n,0);
		tp1 = LPin(th,n,1); tt1 = LTaun(th,n,1);
		tce = exp(j_*ph);
		tc1 += tcc*((VS(2*n+1)*tce + VS(2*n)*conj(tce))*tp1
			+ (VS(2*N+3*n+2)*tce - VS(2*N+3*n)*conj(tce))*tt1 + tt0*VS(2*N+3*n+1));
		tc2 += tcc*((VS(2*n+1)*tce - VS(2*n)*conj(tce))*tt1
			+ (VS(2*N+3*n+2)*tce + VS(2*N+3*n)*conj(tce))*tp1 + tp0*VS(2*N+3*n+1));
		tc *= -j_;
	}
	return (tc1*conj(tc1) + tc2*conj(tc2)).real()/calc_Wsca_1(VS,tC)/tC;
}

double Method3D::balance1(const Vector &VI, const Vector &VS) {
	int n, NN = N*N; double tv1, tv2;
	tv1 = tv2 = 0.;
	for (n=1; n<NN; ++n) {
		tv1 += abs(VI(n)*VI(n)) + abs(VI(NN+n)*VI(NN+n));
		tv2 += abs(VS(n)*VS(n)) + abs(VS(NN+n)*VS(NN+n));
	}
	return fabs(tv2/tv1-1.);
}

double Method3D::balance2(const Vector &VI, const Vector &VS) {
	int n, NN = N*N; double tv1, tv2;
	tv1 = tv2 = 0.;
	for (n=1; n<NN; ++n) {
		tv1 -= (VS(n)*conj(VI(n))).real() + (VS(NN+n)*conj(VI(NN+n))).real(); // ext cs
		tv2 += abs(VS(n)*VS(n)) + abs(VS(NN+n)*VS(NN+n)); // scat cs
	}
	return fabs(0.5*tv2/tv1-1.);
}

Vector Method3D::mul_Vr(const Vector &VI, const Vector &R) {
	int n, m, nm, NN = N*N;
	Vector VS(VI);
	for (n=nm=0; n<N; ++n) for (m=-n; m<n+1; ++m,++nm) {
		VS.Data[nm] = VI(nm)*R(n); VS.Data[NN+nm] = VI(NN+nm)*R(N+n);
	}
	return VS;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			// Method3D_Mie //
 
Matrix Method3D_Mie::calc_RTi(double kr, Complex e1, Complex e2, Complex m1, Complex m2) {
	Complex kR1, kR2, te, tm, tc;
	kR1 = kr*sqrt(e1*m1); if (arg(kR1) < -1.e-8) kR1 = -kR1;
	kR2 = kr*sqrt(e2*m2); if (arg(kR2) < -1.e-8) kR2 = -kR2;
	te = e2/e1; tm = m2/m1;
	Matrix M(4,2*N);// --, --, 10, 11
	memset(M.Data,0,4*N*sizeof(Complex));
	for (int n=0; n<N; ++n) {
		tc = 1./(besj(kR1,n)*bes_dzh1(kR2,n) - tm*besh1(kR2,n)*bes_dzj(kR1,n));
		M.Data[6*N+n] = tc*(tm*besh2(kR2,n)*bes_dzj(kR1,n) - besj(kR1,n)*bes_dzh2(kR2,n)); // 11e
		M.Data[4*N+n] = tc*2.*j_/kR2; // 10e
		tc = 1./(besj(kR1,n)*bes_dzh1(kR2,n) - te*besh1(kR2,n)*bes_dzj(kR1,n));
		M.Data[7*N+n] = tc*(te*besh2(kR2,n)*bes_dzj(kR1,n) - besj(kR1,n)*bes_dzh2(kR2,n)); // 11h
		M.Data[5*N+n] = tc*2.*j_/kR1/tm; // 10h
	}
	return M;
}

Matrix Method3D_Mie::calc_RTii(double kr, Complex e1, Complex e2, Complex m1, Complex m2) {
	Complex kR1, kR2, te, tm, tc;
	kR1 = kr*sqrt(e1*m1); if (arg(kR1) < -1.e-8) kR1 = -kR1;
	kR2 = kr*sqrt(e2*m2); if (arg(kR2) < -1.e-8) kR2 = -kR2;
	te = e2/e1; tm = m2/m1;
	Matrix M(4,2*N);// --, --, 10, 11
	memset(M.Data,0,4*N*sizeof(Complex));
	for (int n=0; n<N; ++n) {
		tc = 1./(besj(kR1,n)*bes_dzh1(kR2,n) - tm*besh1(kR2,n)*bes_dzj(kR1,n));
		M.Data[6*N+n] = tc*(tm*besj(kR2,n)*bes_dzj(kR1,n) - besj(kR1,n)*bes_dzj(kR2,n)); // 11e
		M.Data[4*N+n] = tc*j_/kR2; // 10e
		tc = 1./(besj(kR1,n)*bes_dzh1(kR2,n) - te*besh1(kR2,n)*bes_dzj(kR1,n));
		M.Data[7*N+n] = tc*(te*besj(kR2,n)*bes_dzj(kR1,n) - besj(kR1,n)*bes_dzj(kR2,n)); // 11h
		M.Data[5*N+n] = tc*j_/kR1/tm; // 10h
	}
	return M;
}

Matrix Method3D_Mie::calc_RTo(double kr, Complex e1, Complex e2, Complex m1, Complex m2) {
	Complex kR1, kR2, te, tm, tc;
	kR1 = kr*sqrt(e1*m1); if (arg(kR1) < -1.e-8) kR1 = -kR1;
	kR2 = kr*sqrt(e2*m2); if (arg(kR2) < -1.e-8) kR2 = -kR2;
	te = e1/e2; tm = m1/m2;
	Matrix M(4,2*N);// 00 01 10 11
	memset(M.Data,0,8*N*sizeof(Complex));
	for (int n=0; n<N; ++n) {
		tc = 1./(tm*besh2(kR1,n)*bes_dzh1(kR2,n) - besh1(kR2,n)*bes_dzh2(kR1,n));
		M.Data[0*N+n] = tc*(besh1(kR2,n)*bes_dzh1(kR1,n) - tm*besh1(kR1,n)*bes_dzh1(kR2,n)); // 00e
		M.Data[2*N+n] = tc*2.*j_/kR1; // 01e
		tc = 1./(te*besh2(kR1,n)*bes_dzh1(kR2,n) - besh1(kR2,n)*bes_dzh2(kR1,n));
		M.Data[1*N+n] = tc*(besh1(kR2,n)*bes_dzh1(kR1,n) - te*besh1(kR1,n)*bes_dzh1(kR2,n)); // 00h
		M.Data[3*N+n] = tc*2.*j_/kR2/tm; // 01h
		tc = 1./(tm*besh2(kR1,n)*bes_dzh1(kR2,n) - besh1(kR2,n)*bes_dzh2(kR1,n));
		M.Data[6*N+n] = tc*(besh2(kR2,n)*bes_dzh2(kR1,n) - tm*besh2(kR1,n)*bes_dzh2(kR2,n)); // 11e
		M.Data[4*N+n] = tc*2.*j_*tm/kR2; // 10e
		tc = 1./(te*besh2(kR1,n)*bes_dzh1(kR2,n) - besh1(kR2,n)*bes_dzh2(kR1,n));
		M.Data[7*N+n] = tc*(besh2(kR2,n)*bes_dzh2(kR1,n) - te*besh2(kR1,n)*bes_dzh2(kR2,n)); // 11h
		M.Data[5*N+n] = tc*2.*j_*tm*te/kR1; // 10h
	}
	return M;
}

Matrix Method3D_Mie::calc_RTio(double kr, Complex e1, Complex e2, Complex m1, Complex m2) {
	Complex kR1, kR2, te, tm, tc;
	kR1 = kr*sqrt(e1*m1); if (arg(kR1) < -1.e-8) kR1 = -kR1;
	kR2 = kr*sqrt(e2*m2); if (arg(kR2) < -1.e-8) kR2 = -kR2;
	te = e1/e2; tm = m1/m2;
	Matrix M(4,2*N);// 00 01 10 11
	memset(M.Data,0,8*N*sizeof(Complex));
	for (int n=0; n<N; ++n) {
		tc = 1./(tm*besj(kR1,n)*bes_dzh1(kR2,n) - besh1(kR2,n)*bes_dzj(kR1,n));
		M.Data[0*N+n] = tc*(besh1(kR2,n)*bes_dzh1(kR1,n) - tm*besh1(kR1,n)*bes_dzh1(kR2,n)); // 00e
		M.Data[2*N+n] = tc*j_/kR1; // 01e
		M.Data[6*N+n] = tc*(besh2(kR2,n)*bes_dzj(kR1,n) - tm*besj(kR1,n)*bes_dzh2(kR2,n)); // 11e
		M.Data[4*N+n] = tc*2.*j_*tm/kR2; // 10e
		tc = 1./(te*besj(kR1,n)*bes_dzh1(kR2,n) - besh1(kR2,n)*bes_dzj(kR1,n));
		M.Data[1*N+n] = tc*(besh1(kR2,n)*bes_dzh1(kR1,n) - te*besh1(kR1,n)*bes_dzh1(kR2,n)); // 00h
		M.Data[3*N+n] = tc*j_/kR2; // 01h
		M.Data[7*N+n] = tc*(besh2(kR2,n)*bes_dzj(kR1,n) - te*besj(kR1,n)*bes_dzh2(kR2,n)); // 11h
		M.Data[5*N+n] = tc*2.*j_*tm*te/kR1; // 10h
	}
	return M;
}

Matrix Method3D_Mie::calc_RT(double kr, Complex e1, Complex e2, Complex m1, Complex m2) {
	Complex kR1, kR2, te, tm, tc, tcc;
	kR1 = kr*sqrt(e1*m1); if (arg(kR1) < -1.e-8) kR1 = -kR1;
	kR2 = kr*sqrt(e2*m2); if (arg(kR2) < -1.e-8) kR2 = -kR2;
	te = e1/e2; tm = 1.;//m1/m2;
	Matrix M(4,2*N);//0 1 2 3 -> 00 01 10 11
	memset(M.Data,0,8*N*sizeof(Complex));
	//cout<<kR1<<" "<<kR2<<endl;
	for (int n=0; n<N; ++n) {
		//cout<<n<<endl;
		//tcc = 0.01+j_*0.3;
		//cout<<"  j: "<<besj(tcc,n)<<"   "<<besjd(tcc,n)<<endl;
		//cout<<"  y: "<<besy(tcc,n)<<"   "<<besyd(tcc,n)<<endl;
		//cout<<"  h: "<<besh1(tcc,n)<<"   "<<besh1d(tcc,n)<<endl;
		//cout<<besj(kR1,n)*bes_dzh1(kR2,n) - besh1(kR2,n)*bes_dzj(kR1,n)<<endl;
		tc = 1./(besj(kR1,n)*bes_dzh1(kR2,n) - besh1(kR2,n)*bes_dzj(kR1,n));
		M.Data[0*N+n] = tc*(besh1(kR2,n)*bes_dzh1(kR1,n) - besh1(kR1,n)*bes_dzh1(kR2,n)); // 00e
		M.Data[2*N+n] = tc*j_/kR1; // 01e
		M.Data[6*N+n] = tc*(besj(kR2,n)*bes_dzj(kR1,n) - besj(kR1,n)*bes_dzj(kR2,n)); // 11e
		M.Data[4*N+n] = tc*j_/kR2; // 10e
		tc = 1./(te*besj(kR1,n)*bes_dzh1(kR2,n) - besh1(kR2,n)*bes_dzj(kR1,n));
		M.Data[1*N+n] = tc*(besh1(kR2,n)*bes_dzh1(kR1,n) - te*besh1(kR1,n)*bes_dzh1(kR2,n)); // 00h
		M.Data[3*N+n] = tc*j_/kR2; // 01h
		M.Data[7*N+n] = tc*(besj(kR2,n)*bes_dzj(kR1,n) - te*besj(kR1,n)*bes_dzj(kR2,n)); // 11h
		M.Data[5*N+n] = tc*j_*te/kR1; // 10h

		//cout<<n<<" "<<M.Data[6*N+n]<<" "<<M.Data[7*N+n]<<endl;
		//cin.get();
	}
	//cin.get();
	return M;
}

Matrix Method3D_Mie::calc_TT(double kr, Complex e1, Complex e2, Complex m1, Complex m2) {
	Complex kR1, kR2, te, tm, tc;
	kR1 = kr*sqrt(e1*m1); if (arg(kR1) < -1.e-8) kR1 = -kR1;
	kR2 = kr*sqrt(e2*m2); if (arg(kR2) < -1.e-8) kR2 = -kR2;
	te = e2/e1; tm = 1.;//m1/m2;
	Matrix M(4,2*N);//0 1 2 3 -> a->a b->a a->b b->b
	memset(M.Data,0,8*N*sizeof(Complex));
	for (int n=0; n<N; ++n) {
		tc = -j_*kR2;
		M.Data[0*N+n] = tc*(besj(kR2,n)*bes_dzh1(kR1,n) - besh1(kR1,n)*bes_dzj(kR2,n)); // a->a,e
		M.Data[2*N+n] = tc*(besj(kR2,n)*bes_dzj(kR1,n) - besj(kR1,n)*bes_dzj(kR2,n)); // b->a, e
		M.Data[4*N+n] = tc*(besh1(kR1,n)*bes_dzh1(kR2,n) - besh1(kR2,n)*bes_dzh1(kR1,n)); // a->b,e
		M.Data[6*N+n] = tc*(besj(kR1,n)*bes_dzh1(kR2,n) - besh1(kR2,n)*bes_dzj(kR1,n)); // b->b,e
		tc = -j_*kR1;
		M.Data[1*N+n] = tc*(te*besj(kR2,n)*bes_dzh1(kR1,n) - besh1(kR1,n)*bes_dzj(kR2,n)); // a->a,h
		M.Data[3*N+n] = tc*(te*besj(kR2,n)*bes_dzj(kR1,n) - besj(kR1,n)*bes_dzj(kR2,n)); // b->a, h
		M.Data[5*N+n] = tc*(besh1(kR1,n)*bes_dzh1(kR2,n) - te*besh1(kR2,n)*bes_dzh1(kR1,n)); // a->b,h
		M.Data[7*N+n] = tc*(besj(kR1,n)*bes_dzh1(kR2,n) - te*besh1(kR2,n)*bes_dzj(kR1,n)); // b->b,h
	}
	return M;
}

Vector Method3D_Mie::calc_R(double kr, Complex e1, Complex e2, Complex m1, Complex m2) {
	int n, k;
	Matrix RT(4,2*N); Vector R(2*N);
	RT = calc_RTi(kr,e1,e2,m1,m2);
	memcpy(R.Data,RT.Data+6*N,2*N*sizeof(Complex));
	return R;
}

Vector Method3D_Mie::calc_RR(double kr, Complex e1, Complex e2, Complex m1, Complex m2) {
	int n, k;
	Matrix RT(4,2*N); Vector R(2*N);
	RT = calc_RTii(kr,e1,e2,m1,m2);
	memcpy(R.Data,RT.Data+6*N,2*N*sizeof(Complex));
	return R;
}

Vector Method3D_Mie::calc_RL(double *kr, Complex *ep, Complex *mp, int ns) {
	int n, k;
	Matrix RT(4,2*N); Vector R(2*N);
	RT = calc_RTi(kr[0],ep[0],ep[1],mp[0],mp[1]);
	memcpy(R.Data,RT.Data+6*N,2*N*sizeof(Complex));
	for (k=1; k<ns; ++k) {
		RT = calc_RTo(kr[k],ep[k],ep[k+1],mp[k],mp[k+1]);
		for (n=1; n<N; ++n) {
			R.Data[n] = RT.Data[6*N+n] + R.Data[n]*RT.Data[2*N+n]*RT.Data[4*N+n]/(1. - R.Data[n]*RT.Data[n]);
			R.Data[N+n] = RT.Data[7*N+n] + R.Data[N+n]*RT.Data[3*N+n]*RT.Data[5*N+n]/(1. - R.Data[N+n]*RT.Data[N+n]);
		}
	}
	return R;
}

Matrix Method3D_Mie::calc_SML(double *kr, Complex *ep, Complex *mp, int ns) { // S-matrix algorithm
	int n, k; Complex tc;
	Matrix SML1(4,2*N), SML2(4,2*N), SML(4,2*N);
	SML2 = calc_RT(kr[0],ep[0],ep[1],mp[0],mp[1]);
	for (k=1; k<ns; ++k) {
		SML1 = SML2;
		SML = calc_RT(kr[k],ep[k],ep[k+1],mp[k],mp[k+1]);
		for (n=0; n<N; ++n) {
			tc = 1./(1. - SML1(3,n)*SML(0,n));
			SML2.Data[n+0*N] = SML1(0,n) + tc*SML1(1,n)*SML1(2,n)*SML(0,n); // 00e
			SML2.Data[n+2*N] = tc*SML1(1,n)*SML(1,n); // 01e
			SML2.Data[n+4*N] = tc*SML1(2,n)*SML(2,n); // 10e
			SML2.Data[n+6*N] = SML(3,n) + tc*SML(1,n)*SML(2,n)*SML1(3,n); // 11e
			tc = 1./(1. - SML1(3,n+N)*SML(0,n+N));
			SML2.Data[n+1*N] = SML1(0,n+N) + tc*SML1(1,n+N)*SML1(2,n+N)*SML(0,n+N); // 00h
			SML2.Data[n+3*N] = tc*SML1(1,n+N)*SML(1,n+N); // 01h
			SML2.Data[n+5*N] = tc*SML1(2,n+N)*SML(2,n+N); // 10h
			SML2.Data[n+7*N] = SML(3,n+N) + tc*SML(1,n+N)*SML(2,n+N)*SML1(3,n+N); // 11h
		}
	}
	return SML2;
}

Matrix Method3D_Mie::calc_TML(double *kr, Complex *ep, Complex *mp, int ns) { // T-matrix algorithm
	int n, k; Complex tc;
	Matrix TML1(4,2*N), TML2(4,2*N), TML(4,2*N);
	TML2 = calc_TT(kr[0],ep[0],ep[1],mp[0],mp[1]);
	for (k=1; k<ns; ++k) {
		TML1 = TML2;
		TML = calc_TT(kr[k],ep[k],ep[k+1],mp[k],mp[k+1]);
		for (n=0; n<N; ++n) {
			TML2.Data[n+0*N] = TML(0,n)*TML1(0,n) + TML(1,n)*TML1(2,n); // a->a,e
			TML2.Data[n+2*N] = TML(0,n)*TML1(1,n) + TML(1,n)*TML1(3,n); // b->a,e
			TML2.Data[n+4*N] = TML(2,n)*TML1(0,n) + TML(3,n)*TML1(2,n); // a->b,e
			TML2.Data[n+6*N] = TML(2,n)*TML1(1,n) + TML(3,n)*TML1(3,n); // b->b,e
			TML2.Data[n+1*N] = TML(0,n+N)*TML1(0,n+N) + TML(1,n+N)*TML1(2,n+N); // a->a,h
			TML2.Data[n+3*N] = TML(0,n+N)*TML1(1,n+N) + TML(1,n+N)*TML1(3,n+N); // b->a,h
			TML2.Data[n+5*N] = TML(2,n+N)*TML1(0,n+N) + TML(3,n+N)*TML1(2,n+N); // a->b,h
			TML2.Data[n+7*N] = TML(2,n+N)*TML1(1,n+N) + TML(3,n+N)*TML1(3,n+N); // b->b,h
		}
	}
	return TML2;
}

Matrix Method3D_Mie::calc_SML(Matrix **SM, int ns) { // S-matrix algorithm
	int n, k; Complex tc;
	Matrix SML1(4,2*N), SML2(4,2*N);
	memcpy(SML2.Data,SM[0]->Data,8*N*sizeof(Complex));
	for (k=1; k<ns; ++k) {
		SML1 = SML2;
		for (n=0; n<N; ++n) {
			tc = 1./(1. - SML1(3,n)*(*SM[k])(0,n));
			SML2.Data[n+0*N] = SML1(0,n) + tc*SML1(1,n)*SML1(2,n)*(*SM[k])(0,n); // 00e
			SML2.Data[n+2*N] = tc*SML1(1,n)*(*SM[k])(1,n); // 01e
			SML2.Data[n+4*N] = tc*SML1(2,n)*(*SM[k])(2,n); // 10e
			SML2.Data[n+6*N] = (*SM[k])(3,n) + tc*(*SM[k])(1,n)*(*SM[k])(2,n)*SML1(3,n); // 11e
			tc = 1./(1. - SML1(3,n+N)*(*SM[k])(0,n+N));
			SML2.Data[n+1*N] = SML1(0,n+N) + tc*SML1(1,n+N)*SML1(2,n+N)*(*SM[k])(0,n+N); // 00h
			SML2.Data[n+3*N] = tc*SML1(1,n+N)*(*SM[k])(1,n+N); // 01h
			SML2.Data[n+5*N] = tc*SML1(2,n+N)*(*SM[k])(2,n+N); // 10h
			SML2.Data[n+7*N] = (*SM[k])(3,n+N) + tc*(*SM[k])(1,n+N)*(*SM[k])(2,n+N)*SML1(3,n+N); // 11h
		}
	}
	return SML2;
}

Matrix Method3D_Mie::calc_TML(Matrix **TM, int ns) { // T-matrix algorithm
	int n, k; Complex tc;
	Matrix TML1(4,2*N), TML2(4,2*N);
	memcpy(TML2.Data,TM[0]->Data,8*N*sizeof(Complex));
	for (k=1; k<ns; ++k) {
		TML1 = TML2;
		for (n=0; n<N; ++n) {
			TML2.Data[n+0*N] = (*TM[k])(0,n)*TML1(0,n) + (*TM[k])(1,n)*TML1(2,n); // a->a,e
			TML2.Data[n+2*N] = (*TM[k])(0,n)*TML1(1,n) + (*TM[k])(1,n)*TML1(3,n); // b->a,e
			TML2.Data[n+4*N] = (*TM[k])(2,n)*TML1(0,n) + (*TM[k])(3,n)*TML1(2,n); // a->b,e
			TML2.Data[n+6*N] = (*TM[k])(2,n)*TML1(1,n) + (*TM[k])(3,n)*TML1(3,n); // b->b,e
			TML2.Data[n+1*N] = (*TM[k])(0,n+N)*TML1(0,n+N) + (*TM[k])(1,n+N)*TML1(2,n+N); // a->a,h
			TML2.Data[n+3*N] = (*TM[k])(0,n+N)*TML1(1,n+N) + (*TM[k])(1,n+N)*TML1(3,n+N); // b->a,h
			TML2.Data[n+5*N] = (*TM[k])(2,n+N)*TML1(0,n+N) + (*TM[k])(3,n+N)*TML1(2,n+N); // a->b,h
			TML2.Data[n+7*N] = (*TM[k])(2,n+N)*TML1(1,n+N) + (*TM[k])(3,n+N)*TML1(3,n+N); // b->b,h
		}
	}
	return TML2;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Method3D_GSMSC //

void Method3D_GSMSC::calc_Qt(const Matrix M) {
	int n, m, nm, p, u, N0, N1, N2, NN = N*N, k;
	double *c, *CGm, *CG00, *CG10, *CG11, *CG1m1, cn, cnp;
	
	for (int i=0; i<8; ++i) memset(R[i]->Data,0,NL*N*NN*sizeof(Complex));

	c = new double [2*N-1];
	CGm = new double [2*N-1]; CG00 = new double [2*N-1]; CG10 = new double [2*N-1];
	CG1m1 = new double [2*N-1]; CG11 = new double [2*N-1];
	for (n=0; n<2*N-1; ++n) c[n] = sqrt(double(2*n+1));
	
	for (n=nm=0; n<N; ++n) {
		cn = (n%2) ? -c[n] : c[n]; // (-1)^m * sqrt(2n+1)
		for (m=-n; m<n+1; ++m,++nm) {
			for (p=abs(m); p<N; ++p) {
				cnp = cn*c[p]; // * sqrt(2p+1)
				N0 = abs(n-p);
				N1 = max(N0,1);
				N2 = max(N0,2);

				memset(CGm,0,(2*N-1)*sizeof(double));
				memset(CG00,0,(2*N-1)*sizeof(double));
				memset(CG10,0,(2*N-1)*sizeof(double));
				memset(CG11,0,(2*N-1)*sizeof(double));
				memset(CG1m1,0,(2*N-1)*sizeof(double));

				ClG(CGm,p,n,m,-m);
				ClG(CG00,p,n,0,0);
				ClG(CG1m1,p,n,1,-1);
				ClG(CG10,p,n,1,0);
				ClG(CG11,p,n,1,1);

				for (u=N0; u<n+p+1; ++u) {
					CGm[u-N0] *= cnp/c[u]; CG00[u-N0] *= CGm[u-N0]; CG1m1[u-N0] *= CGm[u-N0];
				}
				for (u=N1; u<n+p+1; ++u) CG10[u-N1] *= CGm[u-N0];
				for (u=N2; u<n+p+1; ++u) CG11[u-N2] *= CGm[u-N0];

				for (k=0; k<NL; ++k) {
					for (u=N0; u<n+p+1; ++u) {
						R[0]->Data[(k*NN+nm)*N+p] += CG00[u-N0]*M(k+0*NL,u);
						R[3]->Data[(k*NN+nm)*N+p] += CG1m1[u-N0]*M(k+3*NL,u);
						R[5]->Data[(k*NN+nm)*N+p] += CG1m1[u-N0]*M(k+5*NL,u);
						R[6]->Data[(k*NN+nm)*N+p] += CG1m1[u-N0]*M(k+6*NL,u);
					}
					for (u=N1; u<n+p+1; ++u) {
						R[1]->Data[(k*NN+nm)*N+p] += CG10[u-N1]*M(k+1*NL,u);
						R[2]->Data[(k*NN+nm)*N+p] += CG10[u-N1]*M(k+2*NL,u);
					}
					for (u=N2; u<n+p+1; ++u) {
						R[4]->Data[(k*NN+nm)*N+p] += CG11[u-N2]*M(k+4*NL,u);
						R[7]->Data[(k*NN+nm)*N+p] += CG11[u-N2]*M(k+7*NL,u);
					}
				}
			}
			cn = -cn;
		}
	}

	delete [] CGm; delete [] CG00; delete [] CG10; delete [] CG1m1; delete [] CG11;
	delete [] c;
}
/**/
/**
Matrix Method3D_GSMSC::calc_Qt(const Matrix M) {
	int n, m, nm, p, u, N0, N1, N2, NN = N*N, k;
	double *c, *cc, *ccc, *CGm, *CG00, *CG10, *CG11, *CG1m1, cn, cnp;
	
	Matrix Q(8*NL,NN*N);
	memset(Q.Data,0,8*NL*N*NN*sizeof(Complex));

	c = new double [2*N-1]; cc = new double [2*N-1];
	CGm = new double [2*N-1]; CG00 = new double [2*N-1]; CG10 = new double [2*N-1];
	CG1m1 = new double [2*N-1]; CG11 = new double [2*N-1];
	for (n=0; n<2*N-1; ++n) cc[n] = 1./(c[n] = sqrt(double(2*n+1)));
	
	for (n=nm=0; n<N; ++n) {
		cn = (n%2) ? -c[n] : c[n]; // (-1)^m * sqrt(2n+1)
		for (m=-n; m<n+1; ++m,++nm) {
			for (p=abs(m); p<N; ++p) {
				cnp = cn*c[p]; // * sqrt(2p+1)
				N0 = abs(n-p);
				N1 = max(N0,1);
				N2 = max(N0,2);

				memset(CGm,0,(2*N-1)*sizeof(double));
				memset(CG00,0,(2*N-1)*sizeof(double));
				memset(CG10,0,(2*N-1)*sizeof(double));
				memset(CG11,0,(2*N-1)*sizeof(double));
				memset(CG1m1,0,(2*N-1)*sizeof(double));

				ClG(CGm,p,n,m,-m);
				ClG(CG00,p,n,0,0);
				ClG(CG10,p,n,1,0);
				ClG(CG11,p,n,1,1);
				ClG(CG1m1,p,n,1,-1);

				for (u=N0; u<n+p+1; ++u) {
					CGm[u-N0] *= cnp*cc[u]; CG00[u-N0] *= CGm[u-N0]; CG1m1[u-N0] *= CGm[u-N0];
				}
				for (u=N1; u<n+p+1; ++u) CG10[u-N1] *= CGm[u-N0]*cc[u];
				for (u=N2; u<n+p+1; ++u) CG11[u-N2] *= CGm[u-N0]*cc[u];

				for (k=0; k<NL; ++k) {
					for (u=N0; u<n+p+1; ++u) {
						Q.Data[((k+0*NL)*NN+nm)*N+p] += CG00[u-N0]*M(k+0*NL,u);
						Q.Data[((k+3*NL)*NN+nm)*N+p] += CG1m1[u-N0]*M(k+3*NL,u);
						Q.Data[((k+5*NL)*NN+nm)*N+p] += CG1m1[u-N0]*M(k+5*NL,u);
						Q.Data[((k+6*NL)*NN+nm)*N+p] += CG1m1[u-N0]*M(k+6*NL,u);
					}
					for (u=N1; u<n+p+1; ++u) {
						Q.Data[((k+1*NL)*NN+nm)*N+p] += CG10[u-N1]*M(k+1*NL,u);
						Q.Data[((k+2*NL)*NN+nm)*N+p] += CG10[u-N1]*M(k+2*NL,u);
					}
					for (u=N2; u<n+p+1; ++u) {
						Q.Data[((k+4*NL)*NN+nm)*N+p] += CG11[u-N2]*M(k+4*NL,u);
						Q.Data[((k+7*NL)*NN+nm)*N+p] += CG11[u-N2]*M(k+7*NL,u);
					}
				}
			}
			cn = -cn;
		}
	}

	delete [] CGm; delete [] CG00; delete [] CG10; delete [] CG1m1; delete [] CG11;
	delete [] c; delete [] cc;

	return Q;
}
/**/
void Method3D_GSMSC::calc_Vt(void) {
		// kR - in the basis medium
	double tn, kR1; // v0 = z; v1 = (z + r*z')/r; v2 = sqrt(n(n+1))*z/r
	for (int k=0; k<NL; ++k) {
		kR1 = 1./kR[k];
		for (int n=0; n<N; ++n) {
			tn = kR1*sqrt(double(n*(n+1)));
			v_reg[2]->Data[k*N+n] = ( v_reg[0]->Data[k*N+n] = v_reg[1]->Data[k*N+n] = besj(kR[k],n) ) * tn;
			v_reg[1]->Data[k*N+n] = v_reg[1]->Data[k*N+n] * kR1 + besjd(kR[k],n);
			v_out[2]->Data[k*N+n] = ( v_out[0]->Data[k*N+n] = v_out[1]->Data[k*N+n] = besh1(kR[k],n) ) * tn;
			v_out[1]->Data[k*N+n] = v_out[1]->Data[k*N+n] * kR1 + besh1d(kR[k],n);
		}
	}
}

void Method3D_GSMSC::calc_Rin(double kr, Complex ec, Vector &R11) {
	Complex kR1, kR2, te = eb/ec;
	kR1 = kr*sqrt(ec); if (arg(kR1) < -1.e-8) kR1 = -kR1;
	kR2 = kr*sqrt(eb); if (arg(kR2) < -1.e-8) kR2 = -kR2;
	for (int n=0; n<N; ++n) {
		R11.Data[n] = (besj(kR1,n)*bes_dzj(kR2,n)-bes_dzj(kR1,n)*besj(kR2,n))/
			(besh1(kR2,n)*bes_dzj(kR1,n)-besj(kR1,n)*bes_dzh1(kR2,n));
		R11.Data[n+N] = (besj(kR1,n)*bes_dzj(kR2,n)-te*bes_dzj(kR1,n)*besj(kR2,n))/
			(te*besh1(kR2,n)*bes_dzj(kR1,n)-besj(kR1,n)*bes_dzh1(kR2,n));
	}
}

void Method3D_GSMSC::calc_RTout(double kr, Complex es, Vector &R00, Vector &R11, Vector &T01, Vector &T10) {
	Complex kR1, kR2, te = es/eb, tc;
	kR1 = kr*sqrt(eb); if (arg(kR1) < -1.e-8) kR1 = -kR1;
	kR2 = kr*sqrt(es); if (arg(kR2) < -1.e-8) kR2 = -kR2;
	for (int n=0; n<N; ++n) {
		tc = 1./(besj(kR1,n)*bes_dzh1(kR2,n) - bes_dzj(kR1,n)*besh1(kR2,n));
		R00.Data[n] = tc*(besh1(kR2,n)*bes_dzh1(kR1,n) - bes_dzh1(kR2,n)*besh1(kR1,n));
		R11.Data[n] = tc*(bes_dzj(kR1,n)*besh2(kR2,n) - besj(kR1,n)*bes_dzh2(kR2,n));
		T01.Data[n] = tc*j_/kR1; T10.Data[n] = tc*2.*j_/kR2;
		tc = 1./(besj(kR1,n)*bes_dzh1(kR2,n) - te*bes_dzj(kR1,n)*besh1(kR2,n));
		R00.Data[n+N] = tc*(te*besh1(kR2,n)*bes_dzh1(kR1,n) - bes_dzh1(kR2,n)*besh1(kR1,n));
		R11.Data[n+N] = tc*(te*bes_dzj(kR1,n)*besh2(kR2,n) - besj(kR1,n)*bes_dzh2(kR2,n));
		T01.Data[n+N] = tc*te*j_/kR2; T10.Data[n+N] = tc*2.*j_/kR1;
	}
}

void Method3D_GSMSC::calc_RT(double kr1, double kr2, Complex ec, Complex es) {
	Complex te, th;
	Vector R00(2*N), R11(2*N), T01(2*N), T10(2*N), RR(2*N);
	calc_Rin(kr1,ec,RR); calc_RTout(kr2,es,R00,R11,T01,T10);
	for (int n=0; n<N; ++n) {
		te = 1./(1.-RR(n)*R00(n)); th = 1./(1.-RR(n+N)*R00(n+N));
			// inc coeff
		cce[n] = RR(n)*( cce[n+N] = T10(n)*te );
		cch[n] = RR(n+N)*( cch[n+N] = T10(n+N)*th );
			// coeff
		cce[n+2*N] = RR(n)*te; cce[n+3*N] = R00(n)*te; cce[n+4*N] = RR(n)*R00(n)*te;
		cch[n+2*N] = RR(n+N)*th; cch[n+3*N] = R00(n+N)*th; cch[n+4*N] = RR(n+N)*R00(n+N)*th;
			// out coeff
		cce[n+6*N] = RR(n)*( cce[n+5*N] = T01(n)*te );
		cce[n+7*N] = R11(n) + RR(n)*T01(n)*T10(n)*te;
		cch[n+6*N] = RR(n+N)*( cch[n+5*N] = T01(n+N)*th );
		cch[n+7*N] = R11(n+N) + RR(n+N)*T01(n+N)*T10(n+N)*th;
	}
}

Vector Method3D_GSMSC::calc_inc(const Vector &VA) {
	int n, m, mn, NN = N*N; Complex t1, t2;
	Vector V(4*NL*NN);
	for (n=mn=0; n<N; ++n) for (m=-n; m<n+1; ++m,++mn) {
		V.Data[0*NL*NN+mn] = VA(mn)*cce[n]; V.Data[2*NL*NN+mn] = VA(mn)*cce[n+N];
		V.Data[1*NL*NN+mn] = VA(NN+mn)*cch[n]; V.Data[3*NL*NN+mn] = VA(NN+mn)*cch[n+N];
	}
	for (int k=1; k<NL; ++k) {
		memcpy(V.Data+(0*NL+k)*NN,V.Data+0*NL*NN,NN*sizeof(Complex));
		memcpy(V.Data+(1*NL+k)*NN,V.Data+1*NL*NN,NN*sizeof(Complex));
		memcpy(V.Data+(2*NL+k)*NN,V.Data+2*NL*NN,NN*sizeof(Complex));
		memcpy(V.Data+(3*NL+k)*NN,V.Data+3*NL*NN,NN*sizeof(Complex));
	}
	return V;
}

Vector Method3D_GSMSC::calc_rad_t(const Vector &VA) {
	int k, n, m, nm, p, pm, nmm, pmm, NN = N*N, nmp, pmn, nmmp, pmmn;
	Complex ce[3], ch[3], JE[3], JM[3], ae, ah, be, bh;
	Vector VB(VA);
	for (k=0; k<NL; ++k) for (n=0; n<N; ++n) for (m=-n; m<n+1; ++m) {
		nm = n*(n+1)+m;
		nmm = n*(n+1)-m;
		JE[0] = JE[1] = JE[2] = JM[0] = JM[1] = JM[2] = 0.;
		for (p=abs(m); p<N; ++p) {
			pm = p*(p+1)+m;
			pmm = p*(p+1)-m;
			
			nmp = nm*N+p;
			nmmp = nmm*N+p;
			pmn = pm*N+n;
			pmmn = pmm*N+n;

			ae = VA((0*NL+k)*NN+pm);
			ah = VA((1*NL+k)*NN+pm);
			be = VA((2*NL+k)*NN+pm);
			bh = VA((3*NL+k)*NN+pm);
				// modified amplitudes
			ce[0] = -j_*( (*v_out[2])(k,p) * ae + (*v_reg[2])(k,p) * be );
			ch[0] = (*v_out[2])(k,p) * ah + (*v_reg[2])(k,p) * bh;
			ce[1] = j_*( (*v_out[0])(k,p) * ae + (*v_reg[0])(k,p) * be )
								 + (*v_out[1])(k,p) * ah + (*v_reg[1])(k,p) * bh;
			ce[2] = j_*( (*v_out[0])(k,p) * ae + (*v_reg[0])(k,p) * be )
								 - (*v_out[1])(k,p) * ah - (*v_reg[1])(k,p) * bh;
			ch[1] =			 (*v_out[0])(k,p) * ah + (*v_reg[0])(k,p) * bh
						- j_*( (*v_out[1])(k,p) * ae + (*v_reg[1])(k,p) * be );
			ch[1] =			 (*v_out[0])(k,p) * ah + (*v_reg[0])(k,p) * bh
						+ j_*( (*v_out[1])(k,p) * ae + (*v_reg[1])(k,p) * be );
				// amplitudes to sources:
			JE[0] += (1./cfe[k]) * ( (*R[0])(k,nmp) * ch[0] + (*R[2])(k,nmp) * ch[1] + (*R[2])(k,nmmp) * ch[2] )
							- (*R[1])(k,nmp) * ce[1] + (*R[1])(k,nmmp) * ce[2];
			JE[1] += cfe[k] * ( (*R[3])(k,nmp) * ce[1] - (*R[4])(k,nmmp) * ce[2] )
							+ (1./cfm[k]) * ( (*R[2])(k,pmn) * ce[0] + (*R[5])(k,nmp) * ce[1] - (*R[7])(k,nmmp) * ce[2] )
							- (*R[1])(k,pmn) * ch[0] - (*R[6])(k,nmp) * ch[1];
			JE[2] += cfe[k] * ( (*R[3])(k,nmmp) * ce[2] - (*R[4])(k,nmp) * ce[1] )
							+ (1./cfm[k]) * ( (*R[2])(k,pmmn) * ce[0] + (*R[5])(k,nmmp) * ce[2] - (*R[7])(k,nmp) * ce[1] )
							+ (*R[1])(k,pmmn) * ch[0] + (*R[6])(k,nmmp) * ch[2];
			JM[0] += (1./cfm[k]) * ( (*R[0])(k,nmp) * ce[0] - (*R[2])(k,nmp) * ce[1] + (*R[2])(k,nmmp) * ce[2] )
							- (*R[1])(k,nmp) * ch[1] + (*R[1])(k,nmmp) * ch[2];
			JM[1] += cfm[k] * ( (*R[3])(k,nmp) * ch[1] - (*R[4])(k,nmmp) * ch[2] )
							- (1./cfe[k]) * ( (*R[2])(k,pmn) * ch[0] - (*R[5])(k,nmp) * ch[1] + (*R[7])(k,nmmp) * ch[2] )
							- (*R[1])(k,pmn) * ce[0] + (*R[6])(k,nmp) * ce[1];
			JM[2] += cfm[k] * ( (*R[3])(k,nmmp) * ch[2] - (*R[4])(k,nmp) * ch[1] )
							- (1./cfe[k]) * ( (*R[2])(k,pmmn) * ch[0] - (*R[5])(k,nmmp) * ch[2] + (*R[7])(k,nmp) * ch[1] )
							+ (*R[1])(k,pmmn) * ce[0] - (*R[6])(k,nmmp) * ce[2];
				// diagonal terms:
			if (p == n) {
				JE[0] -= ch[0]; JE[1] += 0.5*(ce[1] - ch[1]/kR[k]); JE[2] += 0.5*(ce[2] + ch[2]/kR[k]);
				JM[0] -= ce[0]; JM[1] += 0.5*(ch[1] + ce[1]/kR[k]); JM[2] += 0.5*(ch[2] - ce[2]/kR[k]);
			}
		}
			// sources to amplitudes
		VB.Data[(0*NL+k)*NN+nm] = wR[k]*( -(*v_reg[0])(k,n) * (JE[1] + JE[2]) 
																			+(*v_reg[1])(k,n) * (JM[1] - JM[2]) + (*v_reg[2])(k,n) * JM[0] ); // ae
		VB.Data[(2*NL+k)*NN+nm] = wR[k]*( -(*v_out[0])(k,n) * (JE[1] + JE[2]) 
																			+(*v_out[1])(k,n) * (JM[1] - JM[2]) + (*v_out[2])(k,n) * JM[0] ); // be
		VB.Data[(1*NL+k)*NN+nm] = -j_*wR[k]*( (*v_reg[0])(k,n) * (JM[1] + JM[2])
																				+ (*v_reg[1])(k,n) * (JE[1] - JE[2]) + (*v_reg[2])(k,n) * JE[0] ); // ah
		VB.Data[(3*NL+k)*NN+nm] = -j_*wR[k]*( (*v_out[1])(k,n) * (JM[1] + JM[2])
																				+ (*v_out[2])(k,n) * (JE[1] - JE[2]) + (*v_out[2])(k,n) * JE[0] ); // bh
	}

	return VB;
}

Vector Method3D_GSMSC::calc_dif(const Vector &VA) {
	int n, m, nm, k, NN = N*N;
	Complex ae, ah, be, bh, ta1, ta2, tb1, tb2;
	Vector VB(4*NL*NN);
	memset(VB.Data,0,4*NL*NN*sizeof(Complex));
	for (n=nm=0; n<N; ++n) for (m=-n; m<n+1; ++m,++nm) {
		ae = be = ah = bh = 0.;
		for (k=0; k<NL; ++k) {
			ta1 = VA(k*NN+nm); ta2 = VA((2*NL+NL-1-k)*NN+nm);
			VB.Data[k*NN+nm] += 0.5*ta1 + ae; ae += ta1;
			VB.Data[(2*NL+NL-1-k)*NN+nm] += 0.5*ta2 + be; be += ta2;
			ta1 = VA((NL+k)*NN+nm); ta2 = VA((3*NL+NL-1-k)*NN+nm);
			VB.Data[(NL+k)*NN+nm] = 0.5*ta1 + ah; ah += ta1;
			VB.Data[(3*NL+NL-1-k)*NN+nm] = 0.5*ta2 + bh; bh += ta2;
		}
		ta1 = ae*cce[4*N+n] + be*cce[2*N+n];
		tb1 = ae*cce[3*N+n] + be*cce[4*N+n];
		ta2 = ah*cch[4*N+n] + bh*cch[2*N+n];
		tb2 = ah*cch[3*N+n] + bh*cch[4*N+n];
		for (k=0; k<NL; ++k) {
			VB.Data[(0*NL+k)*NN+nm] += ta1; VB.Data[(2*NL+k)*NN+nm] += tb1;
			VB.Data[(1*NL+k)*NN+nm] += ta2; VB.Data[(3*NL+k)*NN+nm] += tb2;
		}
	}
	return VB;
}

Vector Method3D_GSMSC::calc_out(const Vector &VA) {
	int n, m, nm, k, NN = N*N;
	Complex ae, ah, be, bh, tc1, tc2;
	Vector V(2*NN);
	for (n=nm=0; n<N; ++n) for (m=-n; m<n+1; ++m,++nm) {
		ae = ah = be = bh = 0.;
		for (k=0; k<NL; ++k) {
			ae += VA((0*NL+k)*NN+nm); be += VA((2*NL+k)*NN+nm);
			ah += VA((1*NL+k)*NN+nm); bh += VA((3*NL+k)*NN+nm);
		}
		V.Data[nm] = ae*cce[n+5*N] + be*cce[n+6*N];
		V.Data[NN+nm] = ah*cch[n+5*N] + bh*cch[n+6*N];
	}
	return V;
}

void Method3D_GSMSC::add_inc(const Vector &VA, Vector &VB) {
	int n, m, mn, NN = N*N;
	Complex  te, th;
	for (n=mn=0; n<N; ++n) for (m=-n; m<n+1; ++m,++mn) {
		VB.Data[mn] += VA(mn)*cce[n+7*N];
		VB.Data[NN+mn] += VA(NN+mn)*cch[n+7*N];
	}
}

Vector mul_ASct(const Vector &V, void *P) {
	Vector VV(V);
	Method3D_GSMSC *cm = reinterpret_cast<Method3D_GSMSC*>(P);
	VV = VV - cm->calc_dif(cm->calc_rad_t(V));
	return VV;
}
