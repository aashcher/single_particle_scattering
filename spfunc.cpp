
#include "./spfunc.h"

#include <complex>
#include <algorithm>
#include <math.h>
#define _USE_MATH_DEFINES

#define DG 12

void xyz2rtp(double x, double y, double z, double &r, double &th, double &ph) {
	double tv = x*x + y*y;
	r = sqrt(tv+z*z); tv = sqrt(tv);
	if (r < 1.e-14) {th = ph = 0.; return;}
	th = acos(z/r);
	if (tv < 1.e-14) {ph = 0.; return;}
	ph = (y > 0.) ? acos(x/tv) : (M_PI + acos(x/tv));
}

Vector cart2sph(const Vector &V, double th, double ph) {
	Vector W(V);
	W.Data[0] = W.Data[1] = V(0)*cos(ph) + V(1)*sin(ph);
	W.Data[0] = W(0)*sin(th) + V(2)*cos(th);
	W.Data[1] = W(1)*cos(th) - V(2)*sin(th);
	W.Data[2] = V(1)*cos(ph) - V(0)*sin(ph);
	return W;
}

Vector sph2cart(const Vector &V, double th, double ph) {
	Vector W(V);
	W.Data[0] = W.Data[1] = V(0)*sin(th) + V(1)*cos(th);
	W.Data[0] = W(0)*cos(ph) - V(2)*sin(ph);
	W.Data[1] = W(1)*sin(ph) + V(2)*cos(ph);
	W.Data[2] = V(0)*cos(th) - V(1)*sin(th);
	return W;
}

	// spherical Bessel functions

Complex besj(Complex z, int n) {
	if (n == 0) return besj0(z);
	else if (n == 1) return besj1(z);

	if (abs(z) < 0.6) {
		int i, ip;
		double tvv;
		Complex tc, tcc;

		tc = 1.; for (i=1; i<n+1; i++) tc *= z/double(2*i+1); i = 1; tcc = tc;
		ip = int(0.4343*log(abs(tc))) - DG; tvv = (ip < -100) ? 1.e-100 : exp(2.3*ip);
		do {tcc += ( tc *= -0.5*z*z/double(i*(2*(n+i)+1)) ); i++;} while (abs(tc) > tvv);

		return tcc;
	}
	else if (abs(z)/double(n) < 1.) {
		int nn, i, pw;
		double tv1, tv2, tv3, tx, tx1, tx2, ty1, ty2;
		Complex tc, tc1, tc2, tcc;

		if (n > int(abs(z))) pw = abs(int( (n*log(0.5*M_E*abs(z)/n) - 0.5*log(n*abs(z)) - M_LN2)/M_LN10 )) + DG;
		else pw = DG;

		tv1 = abs(z); tv2 = 2./M_E/tv1; tv3 = pw*M_LN10 - M_LN2;
		tx2 = (tx1 = n) + 2;
		ty1 = tx1*log(tv2*tx1) + 0.5*log(tv1*tx1) - tv3;
		ty2 = tx2*log(tv2*tx2) + 0.5*log(tv1*tx2) - tv3;
		do {
			tx = tx2 - ty2*(tx2-tx1)/(ty2-ty1); tx1 = tx2; tx2 = tx;
			ty1 = ty2; ty2 = tx2*log(tv2*tx2) + 0.5*log(tv1*tx2) - tv3;
		} while (fabs(tx2-tx1) > 0.5);
		nn = int(tx2);

		tc1 = 0.; tc2 = exp(-pw*M_LN10);
		for (i=nn; i>0; i--) {
			tc = (2.*i+1.)/z*tc2 - tc1; tc1 = tc2; tc2 = tc;
			if (i == n+1) tcc = tc2;
		}
		return (abs(besj1(z)) > abs(besj0(z))) ? tcc*besj1(z)/tc1 : tcc*besj0(z)/tc2;
	}
	else {
		int i, tm;
		double tvv;
		Complex tc, tp, tq, pp, qq;
		i = -int(0.4343*log(abs(z)/n)) - DG; tvv = (i < -100) ? 1.e-100 : exp(2.3*i);
		tp = pp = 1.; tm = (2*n+1)*(2*n+1); tq = qq = double(tm-1)*(tc = 0.125/z); tc *= tc; i = 1;
		do {
			pp += ( tp *= -tc*double((tm - (4*i-1)*(4*i-1))*(tm - (4*i-3)*(4*i-3)))/double(2*i*(2*i-1)) );
			qq += ( tq *= -tc*double((tm - (4*i+1)*(4*i+1))*(tm - (4*i-1)*(4*i-1)))/double(2*i*(2*i+1)) );
			i++;
		} while ((abs(tp) > tvv) && (abs(tq) > tvv));
		return (pp*cos(z-M_PI_2*(n+1)) - qq*sin(z-M_PI_2*(n+1)))/z;
	}
}

Complex besjd(Complex z, int n) {
	if (n == 0) return besj0d(z);
	else if (n == 1) return besj1d(z);

	if (abs(z) < 0.6) {
		int i, ip;
		double tv, tvv;
		Complex tc, tcc;

		tc = 1/3.; for (i=2; i<n+1; i++) tc *= z/double(2*i+1); i = 1; tcc = double(n)*tc;
		ip = int(0.4343*log(abs(tc))) - DG; tvv = (ip < -100) ? 1.e-100 : exp(2.3*ip);
		do {tv = abs( tc *= -0.5*z*z/double(i*(2*(n+i)+1)) ); tcc += double(n+2*i)*tc; i++;} while (tv > tvv);

		return tcc;
	}
	else if (abs(z)/n < 1.) {
		int nn, i, pw;
		double tv1, tv2, tv3, tx, tx1, tx2, ty1, ty2;
		Complex tc, tc1, tc2, tcc;

		if (n > int(abs(z))) pw = abs(int( (n*log(0.5*M_E*abs(z)/n) - 0.5*log(n*abs(z)) - M_LN2)/M_LN10 )) + DG;
		else pw = DG;

		tv1 = abs(z); tv2 = 2./M_E/tv1; tv3 = pw*M_LN10 - M_LN2;
		tx2 = (tx1 = n) + 2;
		ty1 = tx1*log(tv2*tx1) + 0.5*log(tv1*tx1) - tv3;
		ty2 = tx2*log(tv2*tx2) + 0.5*log(tv1*tx2) - tv3;
		do {
			tx = tx2 - ty2*(tx2-tx1)/(ty2-ty1); tx1 = tx2; tx2 = tx;
			ty1 = ty2; ty2 = tx2*log(tv2*tx2) + 0.5*log(tv1*tx2) - tv3;
		} while (fabs(tx2-tx1) > 0.5);
		nn = int(tx2);

		tc1 = 0.; tc2 = exp(-pw*M_LN10);
		for (i=nn; i>0; i--) {
			tc = (2.*i+1.)/z*tc2 - tc1; tc1 = tc2; tc2 = tc;
			if (i == n+1) tcc = double(n)/z*tc2 - tc1;
		}
		return (abs(besj1(z)) > abs(besj0(z))) ? tcc*besj1(z)/tc1 : tcc*besj0(z)/tc2;
	}
	else {
		return (double(n)*besj(z,n-1) - double(n+1)*besj(z,n+1))/double(2*n+1);
	}
}

Complex besy(Complex z, int n) {
	if (n == 0) return besy0(z);
	else if (n == 1) return besy1(z);

	if (abs(z) < 0.6) {
		int i, ip;
		double tv, tvv;
		Complex tc, tcc;

		tc = -1./z; for (i=1; i<n+1; i++) tc *= double(2*i-1)/z; i = 1; tcc = tc;
		ip = int(0.4343*log(abs(tc))) - DG; tvv = (ip < -100) ? 1.e-100 : exp(2.3*ip);
		do {tv = abs( tc *= 0.5*z*z/double(2*(n-i+1)-1)/double(i) ); tcc += tc; i++;} while (tv > tvv);

		return tcc;
	}
	else if (abs(z)/n < 1.) {
		Complex tc, tc1 = besy0(z), tc2 = besy1(z);
		for (int i=2; i<n+1; i++) {
			tc = double(2*i-1)*tc2/z - tc1; tc1 = tc2; tc2 = tc;
		}
		return tc2;
	}
	else {
		int i, tm;
		double tvv;
		Complex tc, tp, tq, pp, qq;
		i = -int(0.4343*log(abs(z))) - DG; tvv = (i < -100) ? 1.e-100 : exp(2.3*i);
		tp = pp = 1.; tm = (2*n+1)*(2*n+1); tq = qq = double(tm-1)*(tc = 0.125/z); tc *= tc; i = 1;
		do {
			pp += ( tp *= -tc*double((tm - (4*i-1)*(4*i-1))*(tm - (4*i-3)*(4*i-3)))/double(2*i*(2*i-1)) );
			qq += ( tq *= -tc*double((tm - (4*i+1)*(4*i+1))*(tm - (4*i-1)*(4*i-1)))/double(2*i*(2*i+1)) );
			i++;
		} while ((abs(tp) > tvv) && (abs(tq) > tvv));
		return (pp*sin(z-M_PI_2*(n+1)) + qq*cos(z-M_PI_2*(n+1)))/z;
	}
}

Complex besyd(Complex z, int n) {
	if (n == 0) return besy0d(z);
	else if (n == 1) return besy1d(z);

	if (abs(z) < 0.6) {
		int i, ip;
		double tvv;
		Complex tc, tcc;
		tc = 1./z/z; for (i=1; i<n+1; i++) tc *= double(2*i-1)/z; i = 1; tcc = double(n+1)*tc;
		ip = int(0.4343*log(abs(tc))) - DG; tvv = (ip < -100) ? 1.e-100 : exp(2.3*ip);
		do {tcc -= double(2*i-n-1)*( tc *= 0.5*z*z/double(2*(n-i+1)-1)/double(i) ); i++;} while (abs(tc) > tvv);
		return tcc;
	}
	else if (abs(z)/n < 1.) {
		Complex tc, tc1 = besy0(z), tc2 = besy1(z);
		for (int i=2; i<n+1; i++) {
			tc = double(2*i-1)*tc2/z - tc1; tc1 = tc2; tc2 = tc;
		}
		return (tc1 - double(n+1)/z*tc2);
	}
	else return (double(n)*besy(z,n-1) - double(n+1)*besy(z,n+1))/double(2*n+1);
}

Complex besh1(Complex z, int n) {
	if (n == 0) return besh10(z);
	else if (n == 1) return besh11(z);
	return (besj(z,n) + j_*besy(z,n));

	if (z.imag() > -1.e-16*abs(z)) {
		Complex tc, tc1 = besh10(z), tc2 = besh11(z); int i = 1;
		do {tc = tc2*double(2*i+1)/z-tc1; tc1 = tc2; tc2 = tc;} while (++i < n);
		return tc2;
	}
	else return (2.*besj(z,n) - besh2(z,n));
}

Complex besh2(Complex z, int n) {
	if (n == 0) return besh20(z);
	else if (n == 1) return besh21(z);

	if (z.imag() < 1.e-16*abs(z)) {
		Complex tc, tc1 = besh20(z), tc2 = besh21(z); int i = 1;
		do {tc = tc2*double(2*i+1)/z-tc1; tc1 = tc2; tc2 = tc;} while (++i < n);
		return tc2;
	}
	else return (2.*besj(z,n) + besh1(z,n));
}

Complex besh1d(Complex z, int n) {
	if (n == 0) return besh10d(z);
	else if (n == 1) return besh11d(z);

	if (z.imag() > -1.e-16*abs(z)) {
		Complex tc, tc1 = besh10(z), tc2 = besh11(z); int i = 1;
		do {tc = tc2*double(2*i+1)/z-tc1; tc1 = tc2; tc2 = tc;} while (++i < n);
		return (tc1 - double(n+1)*tc2/z);
	}
	else return (2.*besjd(z,n) - besh2d(z,n));
}

Complex besh2d(Complex z, int n) {
	if (n == 0) return besh20d(z);
	else if (n == 1) return besh21d(z);

	if (z.imag() < 1.e-16*abs(z)) {
		Complex tc, tc1 = besh20(z), tc2 = besh21(z); int i = 1;
		do {tc = tc2*double(2*i+1)/z-tc1; tc1 = tc2; tc2 = tc;} while (++i < n);
		return (tc1 - double(n+1)*tc2/z);
	}
	else return (2.*besjd(z,n) + besh1d(z,n));
}

	// Legendre polynomials

double pLeg(double t, int n) {
	switch (n) {
		case 0: return pLeg0(t);
		case 1: return pLeg1(t);
		default:
			int i = 1;
			double tv1 = pLeg0(t), tv2 = pLeg1(t), tv;
			do {tv = ((i+i+1)*tv2*cos(t)-i*tv1)/(i+1); tv1 = tv2; tv2 = tv;} while (++i < n);
			return tv;
	}
}

double pLegn(double t, int nn) {
	if (nn < 0) return 0;
	else if (fabs(t) < 1.e-14) return sqrt(0.5*(2*nn+1));
	else if (fabs(t-M_PI) < 1.e-14) return (nn%2) ? -sqrt(0.5*(2*nn+1)) : sqrt(0.5*(2*nn+1));
	else switch (nn) {
		case 0: return pLegn0(t);
		case 1: return pLegn1(t);
		default: {
			int n = 1; double tv1 = pLegn0(t), tv2 = pLegn1(t), tv, t2 = sqrt(3.), t1 = 1.;
			do {tv = (t2*tv2*cos(t) - n/t1*tv1)/(n+1.); t1 = t2; t2 = sqrt(2*n+3.); tv1 = tv2; tv2 = tv*t2;} while (++n < nn);
			return tv2;
		}
	}
}

void pLegn(double *PL, double t, int n) { // PL -- array of size n+1
	if (n < 0) return;
	else {
		PL[0] = pLegn0(t);
		if (n > 0) {
			PL[1] = pLegn1(t);
			if (n > 1) {
				int i = 1;
				double tv1 = PL[0], tv2 = PL[1], tv, t2 = sqrt(3.), t1 = 1., tc = cos(t);
				do {
					tv = (t2*tv2*tc - double(i)/t1*tv1)/(i+1.); t1 = t2; t2 = sqrt(2*i+3.);
					tv1 = tv2; PL[i+1] = tv2 = tv*t2;
				} while (++i < n);
			}
		}
	}
}

double pLegnd(double t, int nn) {
	switch (nn) {
		case 0: return pLegnd0(t);
		case 1: return pLegnd1(t);
		default: {
			int n = 1; double tv1 = pLegn0(t), tv2 = pLegn1(t), tv, t2 = sqrt(3.), t1 = 1., t0, td0, td1 = 0., td2 = sqrt(1.5);
			do {
				t0 = t1; t1 = t2; t2 = sqrt(2*n+3.);
				tv = t2*(t1*tv2*cos(t) - n/t0*tv1)/(n+1.); tv1 = tv2; tv2 = tv*t2;
				td0 = td1; td1 = td2; td2 = t2*(t1*tv1 + td0/t0);
			} while (++n < nn);
			return td2;
		}
	}
}

double paLeg(double t, int n, int m) {
	if (abs(m) > n) return 0.;
	if (m == 0) return pLeg(t,n);
	else {
		if (n == 1) {if (m == 1) return sin(t); else return -0.5*sin(t);}
		else if ((fabs(t) < 1.e-14) || (fabs(t-M_PI) < 1.e-14)) return 0.;
		else
			if (m > 0) return ((n+m-1)*paLeg(t,n-1,m-1) - (n-m+1)*cos(t)*paLeg(t,n,m-1))/sin(t);
			else return (cos(t)*paLeg(t,n,m+1) - paLeg(t,n-1,m+1))/sin(t)/double(n-m);
	}
}

double paLegd(double t, int n, int m) {
	if (fabs(t) < 1.e-14) {
		switch (m) {
			case 0: return 0.5*n*(n+1);
			case 1: return std::numeric_limits<double>::infinity();
			case 2: return -0.25*double((n-1)*n*(n+1)*(n+2));
			default: return 0.;
		}
	}
	else if (fabs(t-M_PI) < 1.e-14) {
		switch (m) {
			case 0: return ((n+1)%2) ? -0.5*n*(n+1) : 0.5*n*(n+1);
			case 1: return (n%2) ? -std::numeric_limits<double>::infinity() : std::numeric_limits<double>::infinity();
			case 2: return ((n+1)%2) ? 0.25*double((n-1)*n*(n+1)*(n+2)) : -0.25*double((n-1)*n*(n+1)*(n+2));
			default: return 0.;
		}
	}
	return ((n+1)*cos(t)*paLeg(t,n,m)-(n-m+1)*paLeg(t,n+1,m))/sin(t)/sin(t);
}

double paLegn(double t, int n, int m) {
	if (abs(m) > n) return 0.;
	else if (m == 0) return pLegn(t,n);
	else {
		if (n == 1) return (m > 0) ? paLegn11(t) : paLegn1m1(t);
		else if ((fabs(t) < 1.e-14) || (fabs(t-M_PI) < 1.e-14)) return 0.;
		else if (m < 0) return (abs(m)%2) ? -paLegn(t,n,abs(m)) : paLegn(t,n,abs(m));
		else {
			int i; double tv1, tv2, tv, t1, t2, tc = cos(t);
			tv = log(m+0.5); for (i=2; i<2*m+1; ++i) tv -= log(double(i)); tv2 = exp(0.5*tv);
			tv = 0.; for (i=2; 2*i-1<2*m; ++i) tv += log(double(2*i-1));
			if ((tv1 = sin(t)) > 0.) tv2 *= exp(tv)*exp(m*log(tv1));
			else tv2 = (m%2) ? -tv2*exp(tv)*exp(m*log(fabs(tv1))) : tv2*exp(tv)*exp(m*log(fabs(tv1)));
			if (m == n) return tv2;
			tv1 = t2 = 0.; i = m;
			do {
				t1 = t2; t2 = sqrt((i+m+1.)*(i-m+1.)/(2*i+1.)/(2*i+3));
				tv = (tc*tv2 - t1*tv1)/t2; tv1 = tv2; tv2 = tv;
			} while (++i < n);
			return tv2;
		}
	}
}

double paLegnd(double t, int n, int m) {
	if (m < 0) return (abs(m)%2) ? -paLegnd(t,n,abs(m)) : paLegnd(t,n,abs(m)); 
	if (fabs(t) < 1.e-14) {
		switch (m) {
			case 0: return 0.5*n*(n+1)*sqrt(n+0.5);
			case 1: return -std::numeric_limits<double>::infinity();
			case 2: return -0.25*sqrt((n+0.5)*(n-1)*n*(n+1)*(n+2));
			default: return 0.;
		}
	}
	else if (fabs(t-M_PI) < 1.e-14) {
		switch (m) {
			case 0: return (n%2) ? 0.5*n*(n+1)*sqrt(n+0.5) : -0.5*n*(n+1)*sqrt(0.5*(2*n+1.));
			case 1: return (n%2) ? std::numeric_limits<double>::infinity() : -std::numeric_limits<double>::infinity();
			case 2: return (n%2) ? -0.25*sqrt((n+0.5)*(n-1)*n*(n+1)*(n+2)) : 0.25*sqrt((n+0.5)*(n-1)*n*(n+1)*(n+2));
			default: return 0.;
		}
	}
	return ((n+1)*cos(t)*paLegn(t,n,m) - sqrt((n-m+1.)*(n+m+1.)*(2*n+1.)/(2*n+3.))*paLegn(t,n+1,m))/sin(t)/sin(t);
}

void paLegn(double *PA, double t, int n, int m) { // PA -- array of size n+1
	memset(PA,0,(n+1)*sizeof(double));
	if ((abs(m) > n) || (fabs(t) < 1.e-15) || (fabs(t-M_PI) < 1.e-15)) return;
	else if (m == 0) {pLegn(PA,t,n); return;}
	else {
		if (n == 1) {PA[0] = pLegn0(t); PA[1] = (m > 0) ? paLegn11(t) : paLegn1m1(t); return;}
		else {
			int i, mm = abs(m); double tv1, tv2, tv, t1, t2, tc = cos(t);
			tv = log(mm+0.5); for (i=2; i<2*mm+1; ++i) tv -= log(double(i)); tv2 = exp(0.5*tv);
			tv = 0.; for (i=2; 2*i-1<2*mm; ++i) tv += log(double(2*i-1));
			tv2 *= exp(tv)*exp(mm*log(sin(t)));
			tv1 = t2 = 0.; i = mm; PA[i] = tv2;
			if (i < n) do {
				t1 = t2; t2 = sqrt((i+mm+1.)*(i-mm+1.)/(2*i+1.)/(2*i+3));
				tv = (tc*tv2 - t1*tv1)/t2; tv1 = tv2; tv2 = tv;
				PA[i+1] = tv2;
			} while (++i < n);
			if ((m < 0) && (m%2)) for (i=0; i<n+1; ++i) PA[i] = -PA[i];
		}
	}
}

double LPi(double t, int n, int m) {
	if (m < 0) return (abs(m)%2) ? LPi(t,n,abs(m)) : -LPi(t,n,abs(m)); 
	if (fabs(t) < 1.e-14) {
		switch (m) {
			case 1: return 0.5*n*(n+1);
			case -1: return 0.5;
			default: return 0.;
		}
	}
	else if (fabs(t-M_PI) < 1.e-14) {
		switch (m) {
			case 1: return (n%2) ? 0.5*n*(n+1) : -0.5*n*(n+1);
			case -1: return (n%2) ? 0.5 : -0.5;
			default: return 0.;
		}
	}
	else {return m*paLeg(t,n,m)/sin(t);}
}

double LTau(double t, int n, int m) {
	if (m < 0) return (abs(m)%2) ? -LPi(t,n,abs(m)) : LPi(t,n,abs(m)); 
	if (fabs(t) < 1.e-14) {
		switch (m) {
			case 1: return -0.5*n*(n+1);
			case -1: return 0.5;
			default: return 0.;
		}
	}
	else if (fabs(t-M_PI) < 1.e-14) {
		switch (m) {
			case 1: return (n%2) ? 0.5*n*(n+1) : -0.5*n*(n+1);
			case -1: return (n%2) ? -0.5 : 0.5;
			default: return 0.;
		}
	}
	else {return -sin(t)*paLegd(t,n,m);}
}

double LPin(double t, int n, int m) {
	if (m < 0) return (abs(m)%2) ? LPin(t,n,abs(m)) : -LPin(t,n,abs(m)); 
	if (fabs(t) < 1.e-14) {
		switch (m) {
			case 1: return 0.5*sqrt(n*(n+1)*(n+0.5));
			case -1: return 0.5*sqrt(n*(n+1)*(n+0.5));
			default: return 0.;
		}
	}
	else if (fabs(t-M_PI) < 1.e-14) {
		switch (m) {
			case 1: return (n%2) ? 0.5*sqrt(n*(n+1)*(n+0.5)) : -0.5*sqrt(n*(n+1)*(n+0.5));
			case -1: return (n%2) ? 0.5*sqrt(n*(n+1)*(n+0.5)) : -0.5*sqrt(n*(n+1)*(n+0.5));
			default: return 0.;
		}
	}
	else {return m*paLegn(t,n,m)/sin(t);}
}

double LTaun(double t, int n, int m) {
	if (m < 0) return (abs(m)%2) ? -LTaun(t,n,abs(m)) : LTaun(t,n,abs(m)); 
	if (fabs(t) < 1.e-14) {
		switch (m) {
			case 1: return 0.5*sqrt(n*(n+1)*(n+0.5));
			case -1: return -0.5*sqrt(n*(n+1)*(n+0.5));
			default: return 0.;
		}
	}
	else if (fabs(t-M_PI) < 1.e-14) {
		switch (m) {
			case 1: return (n%2) ? -0.5*sqrt(n*(n+1)*(n+0.5)) : 0.5*sqrt(n*(n+1)*(n+0.5));
			case -1: return (n%2) ? 0.5*sqrt(n*(n+1)*(n+0.5)) : -0.5*sqrt(n*(n+1)*(n+0.5));
			default: return 0.;
		}
	}
	else {return -sin(t)*paLegnd(t,n,m);}
}

//===================================================================================================================

	// Wigner D functions

double wig(double ab, int n, int m1, int m2) {
	if ((abs(m1) > n) || (abs(m2) > n)) return 0.;
	if (m1 < -m2) return (abs(m1-m2)%2) ? -wig(ab,n,-m1,-m2) : wig(ab,n,-m1,-m2);
	if (m1 < m2) return (abs(m1-m2)%2) ? -wig(ab,n,m2,m1) : wig(ab,n,m2,m1);

	if (fabs(ab) < 1.e-8) return (m1 == m2) ? 1. : 0.;
	else if (fabs(ab-M_PI) < 1.e-8) return (m1 == (-m2)) ? (((n-m1)%2) ? -1. : 1.) : 0.;
	else if (ab < 0.) return (abs(m1-m2)%2) ? -wig(fabs(ab),n,m1,m2) : wig(fabs(ab),n,m1,m2);

	if (m1 == 0) return pLeg(ab,n);
	else if (m2 == 0) return sqrt(2./(2*n+1.))*paLegn(ab,n,m1);
	else {
		int ii; double tv, tvv, tv1, tv2, t1, t2;
		tv = flog(2*m1) - 2.*flog(m1-abs(m2)) - flog(m1-abs(m2)+1,m1+abs(m2));
		tvv = cos(ab);
		tv1 = exp(0.5*tv - m1*log(2.) + 0.5*((m1-m2)*log(1.-tvv) + (m1+m2)*log(1.+tvv))); if (n == m1) return tv1;
		tv2 = tv1*((m1+1.)*tvv - m2)*sqrt((2*m1+1.)/(m1+m2+1.)/(m1-m2+1.)); if (n == (m1+1)) return tv2;
		ii = m1+1; t2 = sqrt(double((ii-m1)*(ii+m1)*(ii-m2)*(ii+m2)));
		do {
			ii++; t1 = t2; t2 = sqrt(double((ii-m1)*(ii+m1)*(ii-m2)*(ii+m2)));
			tv = (tv2*(2*ii-1.)*(ii*(ii-1)*tvv - m1*m2) - tv1*t1*ii)/t2/(ii-1.); tv1 = tv2; tv2 = tv;
		} while(ii < n);
		return tv2;
	}
}

double wig1(double ab, int n, int m1, int m2) {
	if ((abs(m1) > n) || (abs(m2) > n)) return 0.;
	if (m1 < -m2) return (abs(m1-m2)%2) ? -wig(ab,n,-m1,-m2) : wig(ab,n,-m1,-m2);
	if (m1 < m2) return (abs(m1-m2)%2) ? -wig(ab,n,m2,m1) : wig(ab,n,m2,m1);

	if (fabs(ab) < 1.e-8) return (m1 == m2) ? 1. : 0.;
	else if (fabs(ab-M_PI) < 1.e-8) return (m1 == (-m2)) ? (((n-m1)%2) ? -1. : 1.) : 0.;
	else if (ab < 0.) return (abs(m1-m2)%2) ? -wig(fabs(ab),n,m1,m2) : wig(fabs(ab),n,m1,m2);

	if (m1 == 0) return pLeg(ab,n);
	else if (m2 == 0) {
		return sqrt(2./(2*n+1.))*paLegn(ab,n,m1);
	}
	else {
		int ii, k, i, m, q; double tv1, tv2, *dd, tcp, tcm, ts;
		dd = new double [2*abs(m2)+1];
		tcp = 0.5*(1.+cos(ab)); tcm = 0.5*(1.-cos(ab)); ts = sin(ab);
		for (i=0,ii=m1-abs(m2); ii<m1+abs(m2)+1; ++ii,++i)
			if (ii <= n) {
				dd[i] = sqrt(2./(2*n+1.))*paLegn(ab,n,ii);
			}
			else dd[i] = 0.;
		if (m2 > 0) for (q=1; q<m2+1; ++q) {
			tv2 = dd[q-1]; k = min(2*(m2-q)+1,n+1);
			for (ii=0; ii<k; ++ii) {
				m = m1 - m2 + q + ii; tv1 = tv2; tv2 = dd[q+ii];
				if (m < n) dd[q+ii] = ( sqrt((n+m)*(n-m+1.))*tcp*tv1 - sqrt((n-m)*(n+m+1.))*tcm*dd[q+ii+1] + ts*m*tv2 )/sqrt((n+q)*(n-q+1.));
				else dd[q+ii] = ( sqrt((n+m)*(n-m+1.))*tcp*tv1 + ts*m*tv2 )/sqrt((n+q)*(n-q+1.));
			}
		}
		else for (q=-1; q>m2-1; --q) {
			tv2 = dd[abs(q)-1]; k = min(2*abs(m2-q)+1,n+1);
			for (ii=0; ii<k; ++ii) {
				m = m1 - abs(m2) + abs(q) + ii; tv1 = tv2; tv2 = dd[abs(q)+ii];
				if (m < n) dd[abs(q)+ii] = ( sqrt((n-m)*(n+m+1.))*tcp*dd[abs(q)+ii+1] - sqrt((n+m)*(n-m+1.))*tcm*tv1 + ts*m*tv2 )/sqrt((n-q)*(n+q+1.));
				else dd[abs(q)+ii] = ( -sqrt((n+m)*(n-m+1.))*tcm*tv1 + ts*m*tv2 )/sqrt((n-q)*(n+q+1.));
			}
		}
		tv1 = dd[abs(m2)]; delete [] dd;
		return tv1;
	}
}

	// Clebsch-Gordan coefficients

double ClG(int n, int n1, int n2, int m1, int m2) {
	if ((n1<0)||(n2<0)||(n<0)) return 0;
	if ((n<abs(n1-n2))||(n>(n1+n2))||(abs(m1)>n1)||(abs(m2)>n2)) return 0.;
	int N1 = max(abs(n1-n2),abs(m1+m2)), N2 = n1+n2; if (n < N1) return 0.;

		int ii = N2-1, m = m1 + m2;
		double tv1 = 0., tv2, tv, E, F, D;
		tv2 = exp(0.5*(flog(2*n1)+flog(2*n2)+flog(n1+n2+m)+flog(n1+n2-m)-flog(2*(n1+n2))-flog(n1+m1)-flog(n1-m1)-flog(n2+m2)-flog(n2-m2)));
		if (n == N2) return tv2;
		do {
			D = 2*(ii+1)*sqrt(double((2*ii+1)*(2*ii+3))/double((ii+m+1)*(ii-m+1)*(ii+n2-n1+1)*(ii+n1-n2+1)*(n1+n2-ii)*(n1+n2+ii+2)));
			E = 0.5*(m1-m2 + m*(n2*(n2+1) - n1*(n1+1))/double((ii+1)*(ii+2)));
			F = 0.5/double(ii+2)*sqrt(double((ii+m+2)*(ii-m+2)*(ii+n2-n1+2.)*(ii+n1-n2+2.)*(n1+n2-ii-1.)*(n1+n2+ii+3.))/double((2*ii+3)*(2*ii+5)));
			tv = D*(E*tv2 - F*tv1); tv1 = tv2; tv2 = tv;
		} while (ii-- > n);
		return tv2;
}

double ClG1(int n, int n1, int n2, int m1, int m2) {
	if ((n1<0)||(n2<0)||(n<0)) return 0;
	if ((n<abs(n1-n2))||(n>(n1+n2))||(abs(m1)>n1)||(abs(m2)>n2)||(n<abs(m1+m2))) return 0.;
	int N1 = max(abs(n1-n2),abs(m1+m2)), N2 = n1+n2; if (n < N1) return 0.;

	int ii = N2, m = m1 + m2;
	double tv1 = 0., tv2 = 1., tv;
	while (ii-- > n) {
		tv = ( tv2*double((2*ii+3)*((m1-m2)*(ii+1)*(ii+2)-(m1+m2)*(n1*(n1+1)-n2*(n2+1))))
			- tv1*double((ii+1)*(ii+n1-n2+2)*(ii-n1+n2+2)*(ii+m1+m2+2)*(ii-m1-m2+2)) )/double((ii+2)*(n1+n2-ii)*(n1+n2+ii+2));
		tv1 = tv2; tv2 = tv;
	}
	tv = exp(0.5*(flog(n+m1+m2)+flog(n-m1-m2)+flog(n+n1-n2)+flog(n+n2-n1)-flog(n1+n2-n+1,n1+n2+n+1)-flog(n1+m1)-flog(n1-m1)-flog(n2+m2)-flog(n2-m2)));
	return sqrt(2*n+1.)*tv*tv2;
}

void ClG(double *C, int n1, int n2, int m1, int m2) {
	int ii, n, N1 = max(abs(n1-n2),abs(m1+m2)), N2 = n1+n2, N = N2-N1+1;
	if (N > 0) {
		if ((n1==0)&&(n2==0)) {if ((m1+m2) != 0) C[0] = 0.; else C[0] = (abs(m1)%2) ? -1. : 1.; return;}
		double tv1 = 0., tv2 = 1., tv;
		n = N2;
		tv = exp(0.5*(flog(n+m1+m2)+flog(n-m1-m2)+flog(n+n1-n2)+flog(n+n2-n1)-flog(n1+n2-n+1,n1+n2+n+1)-flog(n1+m1)-flog(n1-m1)-flog(n2+m2)-flog(n2-m2)));
		C[N-1] = sqrt(2*n+1.)*tv*tv2;
		for (ii=N-2; ii>=0; --ii) {
			n = ii + N1;
			tv = ( tv2*double((2*n+3)*((m1-m2)*(n+1)*(n+2)-(m1+m2)*(n1*(n1+1)-n2*(n2+1))))
				- tv1*double((n+1)*(n+n1-n2+2)*(n-n1+n2+2)*(n+m1+m2+2)*(n-m1-m2+2)) )/double((n+2)*(n1+n2-n)*(n1+n2+n+2));
			tv1 = tv2; tv2 = tv;
			tv = exp(0.5*(flog(n+m1+m2)+flog(n-m1-m2)+flog(n+n1-n2)+flog(n+n2-n1)-flog(n1+n2-n+1,n1+n2+n+1)-flog(n1+m1)-flog(n1-m1)-flog(n2+m2)-flog(n2-m2)));
			C[ii] = sqrt(2*n+1.)*tv*tv2;
		}
	}
}

int get_nClG(int n, int m, int n1, int n2) {
	if ((abs(m) > n) || (n < abs(n1-n2)) || (n > (n1+n2))) return 0;
	return ((m+n1+n2-abs(n1-n2-m))/2 - (m-n1-n2+abs(n1-n2+m))/2 + 1);
}

void ClGv(double *c, int n, int m, int n1, int n2) {
	int N = get_nClG(n,m,n1,n2);
	if (N > 0) {
		if (N == 1) {c[0] = 1.; return;}
		int i, mm1, mm2, i1, i2, ii;
		mm2 = m - (mm1 = (m-n1-n2+abs(n1-n2+m))/2);
		i1 = n1*(n1+1); i2 = n2*(n2+1); ii = i1+i2-n*(n+1);
		c[N-1] = 1.;
		c[N-2] = -(ii+2*(mm1+N-1)*(mm2-N+1))/sqrt(double((i1-(mm1+N-1)*(mm1+N-2))*(i2-(mm2-N+1)*(mm2-N+2))));
		for (i=N-3	; i>=0; --i) c[i] = -(c[i+2]*sqrt(double((i1-(mm1+i+1)*(mm1+i+2))*(i2-(mm2-i-1)*(mm2-i-2)))) + c[i+1]*(ii+2*(mm1+i+1)*(mm2-i-1)))
			/sqrt(double((i1-(mm1+i)*(mm1+i+1))*(i2-(mm2-i)*(mm2-i-1))));
		for (i=1; i<N; ++i) c[N-1] += c[N-1-i]*c[N-1-i]; c[N-1] = 1./sqrt(c[N-1]);
		for (i=1; i<N; ++i) c[N-1-i] *= c[N-1];
	}
}

//================================================================

double gaunt_axial(int n1, int n2, int m, int p) {
	double cg = 0.;
	if ((n1 < 0) || (n2 < 0) || (p < 0) || (abs(m) > n1) || (abs(m) > n2) ||
		(p > (n1+n2)) || ((n1+n2-p)%2)) return 0.;
	else if (p == (n1+n2)) {
		return exp( flog(2*n1) + flog(2*n2) - flog(n1) - flog(n2) + 2*flog(n1+n2) - flog(2*n1+2*n2)
			- 0.5*( flog(n1-m) + flog(n2-m) + flog(n1+m) + flog(n2+m) ) );
	}
	else if (p == (n1+n2-2)) {
		return (2*n1+2*n2-3)*(n1*n2-m*m*(2*n1+2*n2-1))/double(2*n1-1)/double(2*n2-1)/double(n1+n2)*
			exp( flog(2*n1) + flog(2*n2) - flog(n1) - flog(n2) + 2*flog(n1+n2) - flog(2*n1+2*n2)
				- 0.5*( flog(n1-m) + flog(n2-m) + flog(n1+m) + flog(n2+m) ) );
	}
	else {
		int ind = 1, q = 0;
		double cg1, cg2, cc[4];
		cg1 = exp( flog(2*n1) + flog(2*n2) - flog(n1) - flog(n2) + 2*flog(n1+n2) - flog(2*n1+2*n2)
			- 0.5*( flog(n1-m) + flog(n2-m) + flog(n1+m) + flog(n2+m) ) );
		cg = (2*n1+2*n2-3)*(n1*n2-m*m*(2*n1+2*n2-1))/double(2*n1-1)/double(2*n2-1)/double(n1+n2)*cg1;
		ind = 1;
		do {
			cg2 = cg1; cg1 = cg;
			for (int i=0; i<4; ++i) {
				q = n1+n2-2*ind-1+i;
				cc[i] = ((n1+n2+1)*(n1+n2+1) - q*q) * (q*q - (n1-n2)*(n1-n2)) / double(4*q*q-1);
			}
			cg = ((cc[1] + cc[2] - 4.*m*m)*cg1 - cc[3]*cg2) / cc[0];
		} while ((n1+n2-2*(++ind)) != p);
		return cg;
	}
}

//===================================================================================================================
	// spherical vector functions

Vector svfRgM(Complex z, double th, double ph, int n, int m) {
	Vector M(3);
	if (n == 0) {M.Data[0] = M.Data[1] = M.Data[2] = 0.;}
	else {
		M.Data[0] = 0.;
		M.Data[1] = j_*LPin(th,n,m)*( M.Data[2] = besj(z,n)*exp(j_*double(m)*ph)*M_SQRT1_2PI/sqrt(n*(n+1.)) );//double(m)*
		M.Data[2] *= -LTaun(th,n,m);
	}
	return M;
}

Vector svfM1(Complex z, double th, double ph, int n, int m) {
	Vector M(3);
	M.Data[0] = 0.;
	M.Data[1] = j_*LPin(th,n,m)*( M.Data[2] = besh1(z,n)*exp(j_*double(m)*ph)*M_SQRT1_2PI/sqrt(n*(n+1.)) );//double(m)*
	M.Data[2] *= -LTaun(th,n,m);
	return M;
}

Vector svfM2(Complex z, double th, double ph, int n, int m) {
	Vector M(3);
	M.Data[0] = 0.;
	M.Data[1] = j_*LPin(th,n,m)*( M.Data[2] = besh2(z,n)*exp(j_*double(m)*ph)*M_SQRT1_2PI/sqrt(n*(n+1.)) );//double(m)*
	M.Data[2] *= -LTaun(th,n,m);
	return M;
}

Vector svfRgMx(Complex z, double th, double ph, int n, int m) {
	double tv = 0.5/sqrt(n*(n+1.)); Complex tc; Vector M(3);
	memset(M.Data,0,3*sizeof(Complex));
	if ((fabs(m) > n) || (n==0)) return M;
	if (n >= abs(m+1)) {
		M.Data[0] = j_*( M.Data[1] = tv*sqrt((n-m)*(n+m+1))*sfRg(z,th,ph,n,m+1) );
	}
	if (n >= abs(m-1)) {
		tc = tv*sqrt((n+m)*(n-m+1))*sfRg(z,th,ph,n,m-1);
		M.Data[0] += j_*tc; M.Data[1] -= tc;
	}
	M.Data[2] = -2.*j_*tv*double(m)*sfRg(z,th,ph,n,m);
	return M;
}

Vector svfM1x(Complex z, double th, double ph, int n, int m) {
	double tv = 0.5/sqrt(n*(n+1.)); Complex tc; Vector M(3);
	memset(M.Data,0,3*sizeof(Complex));
	if ((fabs(m) > n) || (n==0)) return M;
	if (n >= abs(m+1)) {
		M.Data[0] = j_*( M.Data[1] = tv*sqrt((n-m)*(n+m+1))*sf1(z,th,ph,n,m+1) );
	}
	if (n >= abs(m-1)) {
		tc = tv*sqrt((n+m)*(n-m+1))*sf1(z,th,ph,n,m-1);
		M.Data[0] += j_*tc; M.Data[1] -= tc;
	}
	M.Data[2] = -2.*j_*tv*double(m)*sf1(z,th,ph,n,m);
	return M;
}

Vector svfM2x(Complex z, double th, double ph, int n, int m) {
	double tv = 0.5/sqrt(n*(n+1.)); Complex tc; Vector M(3);
	memset(M.Data,0,3*sizeof(Complex));
	if ((fabs(m) > n) || (n==0)) return M;
	if (n >= abs(m+1)) {
		M.Data[0] = j_*( M.Data[1] = tv*sqrt((n-m)*(n+m+1))*sf2(z,th,ph,n,m+1) );
	}
	if (n >= abs(m-1)) {
		tc = tv*sqrt((n+m)*(n-m+1))*sf2(z,th,ph,n,m-1);
		M.Data[0] += j_*tc; M.Data[1] -= tc;
	}
	M.Data[2] = -2.*j_*tv*double(m)*sf2(z,th,ph,n,m);
	return M;
}

Vector svfRgN(Complex z, double th, double ph, int n, int m) {
	double tv = sqrt(n*(n+1.)); Complex tc, tz; Vector N(3);
	tc = exp(j_*double(m)*ph); tz = besj(z,n)/z;
	N.Data[0] = N.Data[2] = N.Data[1] = 0.;
	if (n > 0) {
		N.Data[0] = tv*M_SQRT1_2PI*tz*paLegn(th,n,m)*tc;
		N.Data[2] = j_*LPin(th,n,m)*( N.Data[1] = M_SQRT1_2PI/tv*tc*(tz + besjd(z,n)) );//double(m)*
		N.Data[1] *= LTaun(th,n,m);
	}
	return N;
}

Vector svfN1(Complex z, double th, double ph, int n, int m) {
	double tv = sqrt(n*(n+1.)); Complex tc, tz; Vector N(3);
	tc = exp(j_*double(m)*ph); tz = besh1(z,n)/z;
	N.Data[0] = N.Data[2] = N.Data[1] = 0.;
	if (n > 0) {
		N.Data[0] = tv*M_SQRT1_2PI*tz*paLegn(th,n,m)*tc;
		N.Data[2] = j_*LPin(th,n,m)*( N.Data[1] = M_SQRT1_2PI/tv*tc*(tz + besh1d(z,n)) );//double(m)*
		N.Data[1] *= LTaun(th,n,m);
	}
	return N;
}

Vector svfN2(Complex z, double th, double ph, int n, int m) {
	double tv = sqrt(n*(n+1.)); Complex tc, tz; Vector N(3);
	tc = exp(j_*double(m)*ph); tz = besh2(z,n)/z;
	N.Data[0] = N.Data[2] = N.Data[1] = 0.;
	if (n > 0) {
		N.Data[0] = tv*M_SQRT1_2PI*tz*paLegn(th,n,m)*tc;
		N.Data[2] = j_*LPin(th,n,m)*( N.Data[1] = M_SQRT1_2PI/tv*tc*(tz + besh2d(z,n)) );//double(m)*
		N.Data[1] *= LTaun(th,n,m);
	}
	return N;
}

Vector svfRgNx(Complex z, double th, double ph, int n, int m) {
	double tv; Complex tc; Vector N(3);
	memset(N.Data,0,3*sizeof(Complex));
	if ((fabs(m) > n) || (n==0)) return N;
	tv = 1./sqrt(n*(n+1.)*(2*n+1.));
	N.Data[0] = N.Data[1] = n*sqrt((n+m+1.)*(n+m+2.)/(2*n+3.))*sfRg(z,th,ph,n+1,m+1);
	if (n > abs(m+1)) {
		tc = -(n+1)*sqrt((n-m)*(n-m-1.)/(2*n-1.))*sfRg(z,th,ph,n-1,m+1);
		N.Data[0] += tc; N.Data[1] += tc;
	}
	tc = -n*sqrt((n-m+1.)*(n-m+2.)/(2*n+3.))*sfRg(z,th,ph,n+1,m-1);
	N.Data[0] += tc; N.Data[1] -= tc;
	if (n > abs(m-1)) {
		tc = (n+1)*sqrt((n+m)*(n+m-1.)/(2*n-1.))*sfRg(z,th,ph,n-1,m-1);
		N.Data[0] += tc; N.Data[1] -= tc;
	}
	N.Data[2] = n*sqrt((n-m+1.)*(n+m+1.)/(2*n+3.))*sfRg(z,th,ph,n+1,m);
	if (n > abs(m)) {
		N.Data[2] += (n+1)*sqrt((n-m)*(n+m)/(2*n-1.))*sfRg(z,th,ph,n-1,m);
	}
	N.Data[0] *= 0.5*tv; N.Data[1] *= -0.5*j_*tv; N.Data[2] *= tv;
	return N;
}

Vector svfN1x(Complex z, double th, double ph, int n, int m) {
	double tv; Complex tc; Vector N(3);
	memset(N.Data,0,3*sizeof(Complex));
	if ((fabs(m) > n) || (n==0)) return N;
	tv = 1./sqrt(n*(n+1.)*(2*n+1.));
	N.Data[0] = N.Data[1] = n*sqrt((n+m+1.)*(n+m+2.)/(2*n+3.))*sf1(z,th,ph,n+1,m+1);
	if (n > abs(m+1)) {
		tc = -(n+1)*sqrt((n-m)*(n-m-1.)/(2*n-1.))*sf1(z,th,ph,n-1,m+1);
		N.Data[0] += tc; N.Data[1] += tc;
	}
	tc = -n*sqrt((n-m+1.)*(n-m+2.)/(2*n+3.))*sf1(z,th,ph,n+1,m-1);
	N.Data[0] += tc; N.Data[1] -= tc;
	if (n > abs(m-1)) {
		tc = (n+1)*sqrt((n+m)*(n+m-1.)/(2*n-1.))*sf1(z,th,ph,n-1,m-1);
		N.Data[0] += tc; N.Data[1] -= tc;
	}
	N.Data[2] = n*sqrt((n-m+1.)*(n+m+1.)/(2*n+3.))*sf1(z,th,ph,n+1,m);
	if (n > abs(m)) {
		N.Data[2] += (n+1)*sqrt((n-m)*(n+m)/(2*n-1.))*sf1(z,th,ph,n-1,m);
	}
	N.Data[0] *= 0.5*tv; N.Data[1] *= -0.5*j_*tv; N.Data[2] *= tv;
	return N;
}

Vector svfN2x(Complex z, double th, double ph, int n, int m) {
	double tv; Complex tc; Vector N(3);
	memset(N.Data,0,3*sizeof(Complex));
	if ((fabs(m) > n) || (n==0)) return N;
	tv = 1./sqrt(n*(n+1.)*(2*n+1.));
	N.Data[0] = N.Data[1] = n*sqrt((n+m+1.)*(n+m+2.)/(2*n+3.))*sf2(z,th,ph,n+1,m+1);
	if (n > abs(m+1)) {
		tc = -(n+1)*sqrt((n-m)*(n-m-1.)/(2*n-1.))*sf2(z,th,ph,n-1,m+1);
		N.Data[0] += tc; N.Data[1] += tc;
	}
	tc = -n*sqrt((n-m+1.)*(n-m+2.)/(2*n+3.))*sf2(z,th,ph,n+1,m-1);
	N.Data[0] += tc; N.Data[1] -= tc;
	if (n > abs(m-1)) {
		tc = (n+1)*sqrt((n+m)*(n+m-1.)/(2*n-1.))*sf2(z,th,ph,n-1,m-1);
		N.Data[0] += tc; N.Data[1] -= tc;
	}
	N.Data[2] = n*sqrt((n-m+1.)*(n+m+1.)/(2*n+3.))*sf2(z,th,ph,n+1,m);
	if (n > abs(m)) {
		N.Data[2] += (n+1)*sqrt((n-m)*(n+m)/(2*n-1.))*sf2(z,th,ph,n-1,m);
	}
	N.Data[0] *= 0.5*tv; N.Data[1] *= -0.5*j_*tv; N.Data[2] *= tv;
	return N;
}

///////////////////
