
#include <iostream>
#include <memory.h>

#include "./matrix.h"

#ifdef _OPENMP
	#include <omp.h>
#endif

const double nz_ = 1.e-15;

////////////////////////////////////////////////////////////////
// implementation of the class Vector //////////////////////////
////////////////////////////////////////////////////////////////

Vector::Vector(unsigned Mown) {
	Nrow = Mown;
	Data = new Complex[Nrow];
}

Vector::Vector(const Vector& V) {
	Nrow = V.Nrow; 
	Data = new Complex[Nrow];
	memcpy(Data,V.Data,Nrow*sizeof(Complex));
}

Vector::~Vector() {delete[] Data;}

Vector& Vector::operator = (const Vector& M) {
	if (Nrow != M.Nrow) {if (Nrow) delete[] Data; Nrow = M.Nrow; Data = new Complex[Nrow];}
	memcpy(Data, M.Data, Nrow*sizeof(Complex));
	return *this;
}

Vector& Vector::operator += (const Vector& M) {
	Complex *D1 = Data;
	Complex *D2 = M.Data;
	Complex *D1end = Data + Nrow;
	if ((*this).Nrow != M.Nrow) {cout<<"Vector operator += : wrong vector length\n"; return *this;}
	while (D1 < D1end) *D1++ += *D2++;
	return *this;
}

Vector& Vector::operator -= (const Vector& M) {
	Complex *D1 = Data;
	Complex *D2 = M.Data;
	Complex *D1end = Data + Nrow;
	if ((*this).Nrow != M.Nrow) {cout<<"Vector operator -= : wrong vector length\n"; return *this;}
	while (D1 < D1end) *D1++ -= *D2++;
	return *this;
}

Vector Vector::operator * (const Complex tc) {
	Vector V(*this);
	for (int i=0; i<Nrow; i++) V.Data[i] *= tc;
	return V;
}

Complex Vector::operator * (const Vector& B) const {
	if (Nrow != B.Nrow) {cout<<"Vector operator * Vector : wrong vector length\n"; return 0.;}
	Complex tc = 0.;
#ifdef _OPENMP
	{
		if (Nrow < 1.e6) {
			double ta, taa, tb, tbb, *va, *vb, *vc, *vd;
			va = reinterpret_cast <double*>(Data); vb = reinterpret_cast <double*>(B.Data);
			vc = reinterpret_cast <double*>(&tc); vd=va+2*Nrow; 
			while (va<vd) {ta = *va++; taa = *va++; tb = *vb++; tbb = *vb++; *vc += ta*tb+taa*tbb; *(vc+1) += ta*tbb-taa*tb;} 
		}
		else {
			double tvr = 0., tvi = 0.;
#pragma omp parallel default(shared) reduction(+:tvr,tvi)
			{
				int nt = omp_get_thread_num(), np = omp_get_num_threads(), nn = Nrow/np;
				double ta, taa, tb, tbb, *va, *vb, *vd;
				va = reinterpret_cast<double*>(Data) + 2*nt*nn; vb = reinterpret_cast<double*>(B.Data) + 2*nt*nn;
				vd = (nt == np-1) ? (reinterpret_cast<double*>(Data) + 2*Nrow) : va + 2*nn;
				while (va<vd) {
					ta = *va++; taa = *va++;
					tb = *vb++; tbb = *vb++;
					tvr += ta*tb + taa*tbb;
					tvi += ta*tbb - taa*tb;
				}
			}
			tc = Complex(tvr,tvi);
		}
	}
#else
	{
		double ta, taa, tb, tbb, *va, *vb, *vc, *vd;
		va = reinterpret_cast <double*>(Data); vb = reinterpret_cast <double*>(B.Data);
		vc = reinterpret_cast <double*>(&tc); vd=va+2*Nrow; 
		while (va<vd) {ta = *va++; taa = *va++; tb = *vb++; tbb = *vb++; *vc += ta*tb+taa*tbb; *(vc+1) += ta*tbb-taa*tb;} 
	}
#endif
	return tc;
}

Vector Vector::operator * (const Matrix& B) {
	Vector C(B.Ncol);
	if (Nrow != B.Nrow) {cout<<"Vector operator * Matrix : wrong vector length\n"; return C;}
	double *va, *vb, *vc, *vaa, *vbb, *vcc, ta, taa, tb, tbb;
	va = reinterpret_cast <double*>(Data); vb = reinterpret_cast <double*>(B.Data); 
	vc = reinterpret_cast <double*>(C.Data); vbb = vb + 2*B.Nrow*B.Ncol;
	for (vcc=vc+2*B.Ncol; vc<vcc; vc+=2) {
		*vc = *(vc+1) = 0.; 
		for (vaa=vb; vaa<vbb; vaa += 2*B.Ncol) 
		{ta = *va++; taa = *va++; tb = *vaa; tbb = *(vaa+1); *vc += ta*tb-taa*tbb; *(vc+1) += ta*tbb+taa*tb;}
		va -= 2*Nrow; vb += 2;
	}
	return C;
}

void Vector::GMResG(const Matrix &M, int kn, int knn, double mmax, bool wr) {
	cout.precision(16);
	double *va, *vb, tvr, tvi, tv1r, tv1i, tv2r, tv2i;
	Vector B(Nrow),  BB(Nrow), R(*this);
	Vector **V, **W; V = new Vector* [knn+2]; W = new Vector* [knn+2];
	int k, kk = 0, i, j;
	double tv, *D, max, lmax; D = new double [knn+2];
	Complex tc;
	for (i=0; i<Nrow; i++) BB.Data[i] = 0.;
lbl2:
	k = 0; V[0] = new Vector(Nrow); W[0] = new Vector(Nrow); *V[0] = R; goto lbl3;
	do {
		tv = 0.; for (i=0; i<Nrow; i++) {tc = W[k]->Data[i]; tv += abs(tc*tc);} D[k] = 1./tv;
		tc = 0.; for (i=0; i<Nrow; i++) tc += R.Data[i]*conj(W[k]->Data[i]); tc *= D[k];
		for (i=0; i<Nrow; i++) BB.Data[i] += V[k]->Data[i]*tc;
		for (i=0; i<Nrow; i++) B.Data[i] = W[k]->Data[i]*tc; R -= B;
		max = 0; for (i=0; i<Nrow; i++) {tv = abs(B.Data[i]); if (tv > max) max = tv;}
		if (k) {
			if (max <= mmax) goto lbl1; 
			if (((max > lmax)&&(k >= kn))||(k > knn)) {for (i=0; i<=k; i++) {delete W[i]; delete V[i];} goto lbl2;} 
		}
		lmax = max; kk++; k++; 
		V[k] = new Vector(Nrow); W[k] = new Vector(Nrow); 
		*V[k] = *W[k-1]; 
		if (wr) cout<<"\tGMResG "<<k<<"\t "<<max<<endl;
	lbl3:
		*W[k] = (*V[k])*M;
		B = *W[k]; 
		for (j=0; j<k; j++) {
			va = reinterpret_cast <double*>(W[j]->Data); vb = reinterpret_cast <double*>(B.Data); tvr = tvi = 0.;
			for (i=0; i<Nrow; i++) 
			{tv1r = *va++; tv2r = *vb++; tv1i = *va++; tv2i = *vb++; tvr += tv1r*tv2r + tv1i*tv2i; tvi += tv2i*tv1r - tv2r*tv1i;}
			tc = Complex(tvr,tvi); tc *= D[j]; tvr = tc.real(); tvi = tc.imag();
			va = reinterpret_cast <double*>(W[j]->Data); vb = reinterpret_cast <double*>(W[k]->Data); 
			for (i=0; i<Nrow; i++) {tv1r = *va++; tv1i = *va++; *vb++ -= tvr*tv1r - tvi*tv1i; *vb++ -= tvr*tv1i + tvi*tv1r;}
			va = reinterpret_cast <double*>(V[j]->Data); vb = reinterpret_cast <double*>(V[k]->Data); 
			for (i=0; i<Nrow; i++) {tv1r = *va++; tv1i = *va++; *vb++ -= tvr*tv1r - tvi*tv1i; *vb++ -= tvr*tv1i + tvi*tv1r;} 
		}
	} while (1);
lbl1:
	for (i=0; i<k; i++) {delete W[i]; delete V[i];}
	delete [] V; delete [] W; delete [] D;
	cout<<kk<<" GMResG error: "<<max<<endl;
	*this = BB;
}

void Vector::GMResG(void *PP, Vector (*mul)(const Vector&, void*), int kn, int knn, double mmax, bool wr) {
	cout.precision(16);
	double *va, *vb, tvr, tvi, tv1r, tv1i, tv2r, tv2i;
	Vector B(Nrow), BB(Nrow), R(*this);
	Vector **V, **W; V = new Vector* [knn+2]; W = new Vector* [knn+2];
	int k, kk = 0, i, j;
	double tv, *D, max, lmax; D = new double [knn+2];
	Complex tc;
	for (i=0; i<Nrow; i++) BB.Data[i] = 0.; 
lbl2:
	k = 0; V[0] = new Vector(Nrow); W[0] = new Vector(Nrow); *V[0] = R;
	goto lbl3;
	do {
		tv = 0.; for (i=0; i<Nrow; i++) {tc = W[k]->Data[i]; tv += abs(tc*tc);} D[k] = 1./tv;
		tc = 0.; for (i=0; i<Nrow; i++) tc += R.Data[i]*conj(W[k]->Data[i]); tc *= D[k];   	   
		for (i=0; i<Nrow; i++) BB.Data[i] += V[k]->Data[i]*tc;
		for (i=0; i<Nrow; i++) B.Data[i] = W[k]->Data[i]*tc; R -= B; 
		max = 0; for (i=0; i<Nrow; i++) {tv = abs(B.Data[i]); if (tv > max) max = tv;}
		if (k) {
			if (max <= mmax) goto lbl1; 
			if (((max > lmax)&&(k >= kn))||(k > knn)) {for (i=0; i<=k; i++) {delete W[i]; delete V[i];} goto lbl2;} 
		}
		lmax = max; kk++; k++; 
		V[k] = new Vector(Nrow); W[k] = new Vector(Nrow); 
		*V[k] = *W[k-1]; 
		if (wr) cout<<"\tGMResG "<<k<<"\t "<<max<<endl;
	lbl3:
		*W[k] = mul(*V[k],PP);
		B = *W[k];
		for (j=0; j<k; j++) {
			va = reinterpret_cast <double*>(W[j]->Data); vb = reinterpret_cast <double*>(B.Data); tvr = tvi = 0.;
			for (i=0; i<Nrow; i++) 
			{tv1r = *va++; tv2r = *vb++; tv1i = *va++; tv2i = *vb++; tvr += tv1r*tv2r + tv1i*tv2i; tvi += tv2i*tv1r - tv2r*tv1i;}
			tc = Complex(tvr,tvi); tc *= D[j]; tvr = tc.real(); tvi = tc.imag();
			va = reinterpret_cast <double*>(W[j]->Data); vb = reinterpret_cast <double*>(W[k]->Data); 
			for (i=0; i<Nrow; i++) {tv1r = *va++; tv1i = *va++; *vb++ -= tvr*tv1r - tvi*tv1i; *vb++ -= tvr*tv1i + tvi*tv1r;}
			va = reinterpret_cast <double*>(V[j]->Data); vb = reinterpret_cast <double*>(V[k]->Data); 
			for (i=0; i<Nrow; i++) {tv1r = *va++; tv1i = *va++; *vb++ -= tvr*tv1r - tvi*tv1i; *vb++ -= tvr*tv1i + tvi*tv1r;} 
		}
	} while (1);
lbl1:
	for (i=0; i<k+1; i++) {delete W[i]; delete V[i];}
	delete [] V; delete [] W; delete [] D;
	if (wr) cout<<" GMResG error: "<<max<<endl;
	*this = BB;
}

void Vector::GMResAG(const Matrix& M, int km, double mmax, bool wr) {
	int i, k, kk;
	double tv, max, lmax;
	double *pv1, *pv2, tvr, tvi, tv1r, tv1i, tv2r, tv2i;
	Complex tc, *ss, *cc, *gg, *y, *h1, *h2;
	Vector Z(*this), Z0(*this), **V, **H;

lbl2:
	V = new Vector* [km+1]; H = new Vector* [km+1];
	ss = new Complex [km]; cc = new Complex [km]; gg = new Complex [km+1]; y = new Complex [km];

	for (i=0; i<km+1; i++) gg[i] = 0.;
	kk = 0; V[0] = new Vector (Nrow);
	gg[0] = tv = Z.normF();
	if (tv>1.e-10) tv = 1./tv;
	for (i=0; i<Nrow; i++) V[0]->Data[i] = tv*Z(i);
	do {
		V[kk+1] = new Vector(Nrow); H[kk] = new Vector(kk+2);

		*V[kk+1] = (*V[kk])*M;
		for (k=0; k<kk+1; k++) {
			H[kk]->Data[k] = tc = (*V[k])*(*V[kk+1]); tvr = tc.real(); tvi = tc.imag();
			pv1 = reinterpret_cast<double*>(V[kk+1]->Data); pv2 = reinterpret_cast<double*>(V[k]->Data);
			for (i=0; i<Nrow; i++) {tv1r = *pv2++; tv1i = *pv2++; *pv1++ -= (tvr*tv1r - tvi*tv1i); *pv1++ -= (tvr*tv1i + tvi*tv1r);}
			//for (i=0; i<Nrow; i++) V[kk+1]->Data[i] -= tc*V[k]->Data[i];
		}

		tv = (*V[kk+1]).normF();
		H[kk]->Data[kk+1] = tv;
		if ( fabs(tv)>1.e-10 ) {tv = 1./tv; for (i=0; i<Nrow; i++) V[kk+1]->Data[i] *= tv;}

		h1 = h2 = H[kk]->Data; h2++;
		for (k=0; k<kk; k++, h1++, h2++) {
			*h1 = cc[k]*(tc = *h1) + ss[k]* *h2;
			*h2 = -conj(ss[k])*tc + cc[k]* *h2;
		}

		h1 = h2 = H[kk]->Data + kk; h2++;
		if (abs(*h2) < 1.e-10) {ss[kk] = 0.; cc[kk] = 1.;}
		else if (abs(*h1) < 1.e-10) {ss[kk] = conj(*h2)/abs(*h2); cc[kk] = 0.;}
		else {
			tv = 1./sqrt( std::norm(*h1) + std::norm(*h2) );
			cc[kk] = tv*abs(*h1); ss[kk] = tv*conj(*h2)*(*h1)/abs(*h1);
		}
		*h1 = cc[kk]*(tc = *h1) + ss[kk]* *h2;
		*h2 = -conj(ss[kk])*tc + cc[kk]* *h2;
		gg[kk+1] = -gg[kk]*conj(ss[kk]); gg[kk] *= cc[kk];

		max = abs(gg[kk+1]);
		if (wr) cout<<"GMResAG "<<kk<<"\t "<<tv<<endl;
		//if ( (kk) && ((kk==knn-1) || ((max > lmax) && (kk > kn))) ) {kk++; goto lbl1;}
		//else {lmax = max; kk++;}
		++kk;
	} while ((max > mmax) && (kk < km));
lbl1:
	for (k=kk-1; k>=0; k--) {y[k] = gg[k]; for (i=kk-1; i>k; i--) y[k] -= y[i]*H[i]->Data[k]; y[k] /= H[k]->Data[k];}
	for (i=0; i<Nrow; i++) Z.Data[i] = 0.;
	for (k=0; k<kk; k++) {
		tvr = y[k].real(); tvi = y[k].imag();
		pv1 = reinterpret_cast<double*>(Z.Data); pv2 = reinterpret_cast<double*>(V[k]->Data);
		for (i=0; i<Nrow; i++) {
			tv2r = *pv2++; tv2i = *pv2++;
			*pv1++ += tvr*tv2r - tvi*tv2i; *pv1++ += tvr*tv2i + tvi*tv2r;
		}
	}
	for (k=0; k<kk; k++) {delete V[k]; delete H[k];} delete V[kk];
	delete [] V; delete [] H;
	delete [] ss; delete [] cc; delete [] gg; delete [] y;
	if (max > mmax) {Z = (Z0 = Z0 - Z*M); goto lbl2;}
	cout<<" GMResAG resudial: "<<max<<endl;
	*this = Z;
}

void Vector::GMResAG(void *PP, Vector (*mul)(const Vector&, void*), int km, double mmax, bool wr) {
	int i, k, kk; long sd = 1;
	double tv, max, lmax;
	double *pv1, *pv2, tvr, tvi, tv1r, tv1i, tv2r, tv2i, res1, res2;
	Complex tc, *ss, *cc, *gg, *y, *h1, *h2;
	Vector Z(*this), Z0(*this), **V, **H, X1(*this), X2(*this);
	cout.precision(10);

lbl2:
	V = new Vector* [km+1]; H = new Vector* [km];
	ss = new Complex [km]; cc = new Complex [km]; gg = new Complex [km+1]; y = new Complex [km];

	for (i=0; i<km+1; i++) gg[i] = 0.;

	kk = 0; V[0] = new Vector (Nrow);
	gg[0] = res2 = tv = Z.normF();
	if (tv>1.e-10) tv = 1./tv;
	for (i=0; i<Nrow; i++) V[0]->Data[i] = tv*Z(i);

	do {
		V[kk+1] = new Vector(Nrow); H[kk] = new Vector(kk+2);
		*V[kk+1] = mul(*V[kk],PP);
		//cout<<V[kk+1]->normF(); cin.get();
		for (k=0; k<kk+1; k++) {
			H[kk]->Data[k] = tc = (*V[k])*(*V[kk+1]); tvr = tc.real(); tvi = tc.imag();
			pv1 = reinterpret_cast<double*>(V[kk+1]->Data); pv2 = reinterpret_cast<double*>(V[k]->Data);
			for (i=0; i<Nrow; i++) {
				tv1r = *pv2++; tv1i = *pv2++;
				*pv1++ -= (tvr*tv1r - tvi*tv1i);
				*pv1++ -= (tvr*tv1i + tvi*tv1r);
			} //for (i=0; i<Nrow; i++) V[kk+1]->Data[i] -= tc*V[k]->Data[i];
		}
		tv = V[kk+1]->normF();
		H[kk]->Data[kk+1] = tv;
		if (fabs(tv) > 1.e-14) {tv = 1./tv; for (i=0; i<Nrow; i++) V[kk+1]->Data[i] *= tv;}
		//else {cout<<"gmres error\n"; cin.get(); return;}
		h1 = h2 = H[kk]->Data; h2++;
		for (k=0; k<kk; k++, h1++, h2++) {
			*h1 = cc[k]*(tc = *h1) + ss[k]* *h2;
			*h2 = -conj(ss[k])*tc + cc[k]* *h2;
		}
		h1 = h2 = H[kk]->Data + kk; h2++;
		if (abs(*h2) < 1.e-10) {ss[kk] = 0.; cc[kk] = 1.;}
		else if (abs(*h1) < 1.e-10) {ss[kk] = conj(*h2)/abs(*h2); cc[kk] = 0.;}
		else {
			tv = 1./sqrt( std::norm(*h1) + std::norm(*h2) );
			cc[kk] = tv*abs(*h1); ss[kk] = tv*conj(*h2)*(*h1)/abs(*h1);
		}
		*h1 = cc[kk]*(tc = *h1) + ss[kk]* *h2;
		*h2 = -conj(ss[kk])*tc + cc[kk]* *h2;
		gg[kk+1] = -gg[kk]*conj(ss[kk]); gg[kk] *= cc[kk];
/**/
			// calculate intermediate solution
		if (1) {
			X1 = X2;
			for (k=kk; k>=0; k--) {y[k] = gg[k]; for (i=kk; i>k; i--) y[k] -= y[i]*H[i]->Data[k]; y[k] /= H[k]->Data[k];}
			for (i=0; i<Nrow; i++) X2.Data[i] = 0.;
			if (kk > 0) for (k=0; k<kk+1; k++) {for (i=0; i<Nrow; i++) X2.Data[i] += y[k]*V[k]->Data[i];}		
		}
/**/
		max = abs(gg[kk+1]); res1 = res2; res2 = max;
		//if (wr) cout<<"\tgmres iteration "<<kk<<" residual: "<<max<<endl; //cin.get();
		if (wr) cout<<"\tgmres iteration "<<kk<<" residual: "<<max<<"   solution: "<<X2.cmp(X1)/X2.normF()<<endl; //cin.get();
	} while ((max > mmax) && (fabs(res2/res1-1.) > 1.e-10) && (++kk < km));
lbl1:
	if ((max > mmax) && (kk == km)) kk--;
	for (k=kk; k>=0; k--) {y[k] = gg[k]; for (i=kk; i>k; i--) y[k] -= y[i]*H[i]->Data[k]; y[k] /= H[k]->Data[k];}
	for (i=0; i<Nrow; i++) Z.Data[i] = 0.;
	if (kk > 0) for (k=0; k<kk+1; k++) {
		//if (wr) {cout<<"  gmres out: "<<k<<" "<<y[k]<<" "<<V[k]->Data[8]<<" "<<V[k]->Data[24]<<endl;}
		tvr = y[k].real(); tvi = y[k].imag();
		pv1 = reinterpret_cast<double*>(Z.Data); pv2 = reinterpret_cast<double*>(V[k]->Data);
		for (i=0; i<Nrow; i++) {
			tv2r = *pv2++; tv2i = *pv2++;
			*pv1++ += tvr*tv2r - tvi*tv2i; *pv1++ += tvr*tv2i + tvi*tv2r;
		}//		for (i=0; i<Nrow; i++) Z.Data[i] += y[k]*V[k]->Data[i];
	}
	else Z = *V[0];
	for (k=0; k<kk+1; k++) {delete V[k]; delete H[k];} delete V[kk+1]; delete [] V; delete [] H;
	delete [] ss; delete [] cc; delete [] gg; delete [] y;
	//if (max > mmax) {Z = mul(Z,PP); Z = (Z0 = Z0 - Z); goto lbl2;}
	cout<<" gmres resudial: "<<max<<"   "<<kk<<endl;
	*this = Z;
}

void Vector::GMResM(void *PP, Vector (*mul)(const Vector&, void*), const Vector &X0, int km, double mmax, bool wr) {
	int i, k, kk;
	double tv, max, lmax;
	double *pv1, *pv2, tvr, tvi, tv1r, tv1i, tv2r, tv2i;
	Complex tc, *ss, *cc, *gg, *y, *h1, *h2;
	Vector Z(Nrow), X(X0), B(*this), **V, **H;
	cout.precision(10);

	do {
		Z = B - mul(X,PP);

		V = new Vector* [km+1]; H = new Vector* [km];
		ss = new Complex [km]; cc = new Complex [km]; gg = new Complex [km+1]; y = new Complex [km];

		for (i=0; i<km+1; i++) gg[i] = 0.;
		kk = 0; V[0] = new Vector (Nrow);
		gg[0] = tv = Z.normF();
		if (tv > 1.e-14) tv = 1./tv;
		for (i=0; i<Nrow; i++) V[0]->Data[i] = tv*Z(i);

		do {
			V[kk+1] = new Vector(Nrow); H[kk] = new Vector(kk+2);
			*V[kk+1] = mul(*V[kk],PP);
			for (k=0; k<kk+1; k++) {
				H[kk]->Data[k] = tc = (*V[k])*(*V[kk+1]);
				for (i=0; i<Nrow; i++) V[kk+1]->Data[i] -= tc*V[k]->Data[i];
			}
			H[kk]->Data[kk+1] = tv = (*V[kk+1]).normF();
			if (fabs(tv) > 1.e-14) {tv = 1./tv; for (i=0; i<Nrow; i++) V[kk+1]->Data[i] *= tv;}

			h1 = h2 = H[kk]->Data; h2++;
			for (k=0; k<kk; k++, h1++, h2++) {
				*h1 = cc[k]*(tc = *h1) + ss[k]* *h2;
				*h2 = -conj(ss[k])*tc + cc[k]* *h2;
			}
			//h1 = h2 = H[kk]->Data + kk; h2++;
			if (abs(*h2) < 1.e-14) {ss[kk] = 0.; cc[kk] = 1.;}
			else if (abs(*h1) < 1.e-14) {ss[kk] = conj(*h2)/abs(*h2); cc[kk] = 0.;}
			else {
				tv = 1./sqrt( std::norm(*h1) + std::norm(*h2) );
				cc[kk] = tv*abs(*h1); ss[kk] = tv*conj(*h2)*(*h1)/abs(*h1);
			}
			*h1 = cc[kk]*(tc = *h1) + ss[kk]* *h2;
			*h2 = -conj(ss[kk])*tc + cc[kk]* *h2;
			gg[kk+1] = -gg[kk]*conj(ss[kk]); gg[kk] *= cc[kk];

			max = abs(gg[kk+1]);
			if (wr) cout<<"\tGMResAG "<<kk<<"\t "<<max<<endl;
		} while ((max > mmax) && (++kk < km));

		for (k=kk-1; k>=0; k--) {y[k] = gg[k]; for (i=kk-1; i>k; i--) y[k] -= y[i]*H[i]->Data[k]; y[k] /= H[k]->Data[k];}
			// retrieve solution
		memset(X.Data,0,Nrow*sizeof(Complex));
		for (k=0; k<kk; k++) X += y[k]* *V[k];
		for (k=0; k<kk; k++) {delete V[k]; delete H[k];} delete V[kk]; delete [] V; delete [] H;
		delete [] ss; delete [] cc; delete [] gg; delete [] y;

		cin.get();
	} while (max > mmax);

	cout<<" GMResAG resudial: "<<max<<"   "<<kk<<endl;
	*this = X;
}

void Vector::GMResAH(const Matrix& M, int kn, int knn, double mmax, bool wr) {
	int i, k, kk;
	double tv, *nv;
	Complex tc, *ss, *cc, *gg, *y, *h1, *h2;
	Vector Z(*this), **V, **H;

	nv = new double [knn+1];
	ss = new Complex [knn]; cc = new Complex [knn];
	gg = new Complex [knn+1]; y = new Complex [knn];
	V = new Vector * [knn+1]; H = new Vector * [knn+1];

	for (i=0; i<knn+1; i++) gg[i] = 0.;

	V[0] = new Vector(Nrow);
	*V[0] = Z; if (abs(Z(0))>1.e-10) V[0]->Data[0] += Z.normF()*Z(0)/abs(Z(0));
	tc = (*V[0])*(*V[0]); nv[0] = 1./tc.real(); tc = (*V[0])*Z; tc *= 2.*nv[0];
	for (i=0; i<Nrow; i++) Z.Data[i] -= tc*V[0]->Data[i];
	gg[0] = Z(0);
	tc = -2.*conj(V[0]->Data[0])*nv[0];
	for (i=0; i<Nrow; i++) Z.Data[i] = tc*V[0]->Data[i]; Z.Data[0] += 1.;
	Z = Z*M;
	tc = (*V[0])*Z; tc *= 2.*nv[0];
	for (i=0; i<Nrow; i++) Z.Data[i] -= tc*V[0]->Data[i];

	kk = 0;
	do {
		V[kk+1] = new Vector(Nrow);
		for (i=0; i<kk+1; i++) V[kk+1]->Data[i] = 0.;
		for (i=kk+1; i<Nrow; i++) V[kk+1]->Data[i] = Z(i);
		if ( abs(Z(kk+1))>1.e-10 ) V[kk+1]->Data[kk+1] += Z.normF(kk+1)*Z(kk+1)/abs(Z(kk+1));
		tc = (*V[kk+1])*(*V[kk+1]); nv[kk+1] = 1./tc.real(); tc = (*V[kk+1])*Z; tc *= 2.*nv[kk+1];
		for (i=0; i<Nrow; i++) Z.Data[i] -= tc*V[kk+1]->Data[i];

		H[kk] = new Vector(kk+2);
		for (i=0; i<kk+2; i++) {H[kk]->Data[i] = Z(i);}// cout<<" H "<<i<<" "<<Z(i)<<endl;}

		tc = -2.*conj(V[kk+1]->Data[kk+1])*nv[kk+1];
		for (i=0; i<Nrow; i++) Z.Data[i] = tc*V[kk+1]->Data[i]; Z.Data[kk+1] += 1.;
		for (k=kk; k>=0; k--) {tc = (*V[k])*Z; tc *= 2.*nv[k]; for (i=0; i<Nrow; i++) Z.Data[i] -= tc*V[k]->Data[i];}
		Z = Z*M;
		for (k=0; k<kk+2; k++) {tc = (*V[k])*Z; tc *= 2.*nv[k]; for (i=0; i<Nrow; i++) Z.Data[i] -= tc*V[k]->Data[i];}

		h1 = h2 = H[kk]->Data; h2++;
		for (k=0; k<kk; k++, h1++, h2++) {
			*h1 = cc[k]*(tc = *h1) + ss[k]* *h2;
			*h2 = -conj(ss[k])*tc + cc[k]* *h2;
		}

		h1 = h2 = H[kk]->Data + kk; h2++;
		if ( abs(*h2)<1.e-10 ) {ss[kk] = 0.; cc[kk] = 1.;}
		else if ( abs(*h1)<1.e-10 ) {ss[kk] = conj(*h2)/abs(*h2); cc[kk] = 0.;}
		else {
			tv = 1./sqrt( std::norm(*h1) + std::norm(*h2) );
			cc[kk] = tv*abs(*h1); ss[kk] = tv*conj(*h2)*(*h1)/abs(*h1);
		}
		*h1 = cc[kk]*(tc = *h1) + ss[kk]* *h2;
		*h2 = -conj(ss[kk])*tc + cc[kk]* *h2;
		gg[kk+1] = -gg[kk]*conj(ss[kk]); gg[kk] *= cc[kk];

		tv = abs(gg[kk+1]);
		if (wr) cout<<"\tGMResAH "<<kk<<"\t "<<tv<<endl;
		kk++;
	} while ( (kk < knn) && (tv > mmax) );
lbl1:
	cout<<" GMResAH resudial: "<<abs(gg[kk+1])<<endl;
	for (k=kk-1; k>=0; k--) {
		y[k] = gg[k];
		for (i=kk-1; i>k; i--) y[k] -= y[i]*H[i]->Data[k];
		y[k] /= H[k]->Data[k];
	}
	for (i=0; i<Nrow; i++) Z.Data[i] = 0.;
	for (k=kk-1; k>=0; k--) {
		Z.Data[k] += y[k];
		tc = (*V[k])*Z; tc *= 2.*nv[k]; for (i=0; i<Nrow; i++) Z.Data[i] -= tc*V[k]->Data[i];
	}

	delete [] nv;
	delete [] ss; delete [] cc;
	delete [] gg; delete [] y;
	for (k=0; k<kk; k++) {delete V[k]; delete H[k];} delete V[kk]; delete [] V; delete [] H;

	*this = Z;
}

void Vector::GMResAH(void *PP, Vector (*mul)(const Vector&, void*), int kn, int knn, double mmax, bool wr) {
	int i, k, kk;
	double tv, *nv;
	double *pv1, *pv2, tvr, tvi, tv1r, tv1i;
	Complex tc, *tcl, *ss, *cc, *gg, *y, *h1, *h2;
	Vector Z(*this), **V, **H;

	nv = new double [knn+1];
	ss = new Complex [knn]; cc = new Complex [knn];
	gg = new Complex [knn+1]; y = new Complex [knn];
	V = new Vector * [knn+1]; H = new Vector * [knn+1];

	for (i=0; i<knn+1; i++) gg[i] = 0.;

	V[0] = new Vector(Z);
	if (abs(Z(0))>1.e-10) V[0]->Data[0] += Z.normF()*Z(0)/abs(Z(0));
	else V[0]->Data[0] += Z.normF();
	tc = (*V[0])*(*V[0]); nv[0] = 1./tc.real(); tc = (*V[0])*Z; tc *= 2.*nv[0];
	for (i=0; i<Nrow; i++) Z.Data[i] -= tc*V[0]->Data[i];
	gg[0] = Z(0);
	tc = -2.*conj(V[0]->Data[0])*nv[0];
	for (i=0; i<Nrow; i++) Z.Data[i] = tc*V[0]->Data[i]; Z.Data[0] += 1.;
	Z = mul(Z,PP);
	tc = (*V[0])*Z; tc *= 2.*nv[0];
	for (i=0; i<Nrow; i++) Z.Data[i] -= tc*V[0]->Data[i];

	kk = 0;
	do {
		V[kk+1] = new Vector(Nrow);
		for (i=0; i<kk+1; i++) V[kk+1]->Data[i] = 0.;
		for (i=kk+1; i<Nrow; i++) V[kk+1]->Data[i] = Z(i);
		if ( abs(Z(kk+1))>1.e-10 ) V[kk+1]->Data[kk+1] += Z.normF(kk+1)*Z(kk+1)/abs(Z(kk+1));
		else V[kk+1]->Data[kk+1] += Z.normF(kk+1);
		tc = (*V[kk+1])*(*V[kk+1]); nv[kk+1] = 1./tc.real(); tc = (*V[kk+1])*Z; tc *= 2.*nv[kk+1];
		pv1 = reinterpret_cast<double*>(Z.Data); pv2 = reinterpret_cast<double*>(V[kk+1]->Data); tvr = tc.real(); tvi = tc.imag();
		for (i=0; i<Nrow; i++) {tv1r = *pv2++; tv1i = *pv2++; *pv1++ -= tvr*tv1r - tvi*tv1i; *pv1++ -= tvr*tv1i + tvi*tv1r;}
		//for (i=0; i<Nrow; i++) Z.Data[i] -= tc*V[kk+1]->Data[i];

		H[kk] = new Vector(kk+2);
		for (i=0; i<kk+2; i++) {H[kk]->Data[i] = Z(i);}

		tc = -2.*conj(V[kk+1]->Data[kk+1])*nv[kk+1];
		pv1 = reinterpret_cast<double*>(Z.Data); pv2 = reinterpret_cast<double*>(V[kk+1]->Data); tvr = tc.real(); tvi = tc.imag();
		for (i=0; i<Nrow; i++) {tv1r = *pv2++; tv1i = *pv2++; *pv1++ = tvr*tv1r - tvi*tv1i; *pv1++ = tvr*tv1i + tvi*tv1r;}
		//for (i=0; i<Nrow; i++) Z.Data[i] = tc*V[kk+1]->Data[i];
		Z.Data[kk+1] += 1.;
		for (k=kk; k>=0; k--) {
			tc = (*V[k])*Z; tc *= 2.*nv[k];
			pv1 = reinterpret_cast<double*>(Z.Data); pv2 = reinterpret_cast<double*>(V[k]->Data); tvr = tc.real(); tvi = tc.imag();
			for (i=0; i<Nrow; i++) {tv1r = *pv2++; tv1i = *pv2++; *pv1++ -= tvr*tv1r - tvi*tv1i; *pv1++ -= tvr*tv1i + tvi*tv1r;}
			//for (i=0; i<Nrow; i++) Z.Data[i] -= tc*V[k]->Data[i];
		}
		Z = mul(Z,PP);
		for (k=0; k<kk+2; k++) {
			tc = (*V[k])*Z; tc *= 2.*nv[k];
			pv1 = reinterpret_cast<double*>(Z.Data); pv2 = reinterpret_cast<double*>(V[k]->Data); tvr = tc.real(); tvi = tc.imag();
			for (i=0; i<Nrow; i++) {tv1r = *pv2++; tv1i = *pv2++; *pv1++ -= tvr*tv1r - tvi*tv1i; *pv1++ -= tvr*tv1i + tvi*tv1r;}
			//for (i=0; i<Nrow; i++) Z.Data[i] -= tc*V[k]->Data[i];
		}

		h1 = h2 = H[kk]->Data; h2++;
		for (k=0; k<kk; k++, h1++, h2++) {
			*h1 = cc[k]*(tc = *h1) + ss[k]* *h2;
			*h2 = -conj(ss[k])*tc + cc[k]* *h2;
		}

		h1 = h2 = H[kk]->Data + kk; h2++;
		if ( abs(*h2)<1.e-10 ) {ss[kk] = 0.; cc[kk] = 1.;}
		else if ( abs(*h1)<1.e-10 ) {ss[kk] = conj(*h2)/abs(*h2); cc[kk] = 0.;}
		else {
			tv = 1./sqrt( std::norm(*h1) + std::norm(*h2) );
			cc[kk] = tv*abs(*h1); ss[kk] = tv*conj(*h2)*(*h1)/abs(*h1);
		}
		*h1 = cc[kk]*(tc = *h1) + ss[kk]* *h2;
		*h2 = -conj(ss[kk])*tc + cc[kk]* *h2;
		gg[kk+1] = -gg[kk]*conj(ss[kk]); gg[kk] *= cc[kk];

		tv = abs(gg[kk+1]);
		if (wr) cout<<"\tGMResAH "<<kk<<"\t "<<tv<<endl;
		kk++;
	} while ( (kk < knn) && (tv > mmax) );
lbl1:
	cout<<" GMResAH resudial: "<<abs(gg[kk+1])<<endl;
	for (k=kk-1; k>=0; k--) {
		y[k] = gg[k];
		for (i=kk-1; i>k; i--) y[k] -= y[i]*H[i]->Data[k];
		y[k] /= H[k]->Data[k];
	}
	for (i=0; i<Nrow; i++) Z.Data[i] = 0.;
	for (k=kk-1; k>=0; k--) {
		Z.Data[k] += y[k]; tc = (*V[k])*Z; tc *= 2.*nv[k];
		pv1 = reinterpret_cast<double*>(Z.Data); pv2 = reinterpret_cast<double*>(V[k]->Data); tvr = tc.real(); tvi = tc.imag();
		for (i=0; i<Nrow; i++) {tv1r = *pv2++; tv1i = *pv2++; *pv1++ -= tvr*tv1r - tvi*tv1i; *pv1++ -= tvr*tv1i + tvi*tv1r;}
		//		for (i=0; i<Nrow; i++) Z.Data[i] -= tc*V[k]->Data[i];
	}

	delete [] nv;
	delete [] ss; delete [] cc;
	delete [] gg; delete [] y;
	for (k=0; k<kk; k++) {delete V[k]; delete H[k];} delete V[kk]; delete [] V; delete [] H;

	*this = Z;
}

void Vector::CG(const Matrix &M, int kn, double mmax, bool wr) {
	int k = 0; double tv; Complex tc, tcc;
	Vector R(*this), P(*this), PP(*this);
	memset(Data,0,Nrow*sizeof(Complex));
	do {
		PP = P*M; tc = (tcc = (R*R))/(PP*P);
		*this += tc*P; R -= tc*PP; tv = abs(tc)*P.normF();
		P = R + ((R*R)/tcc)*P;
		if (wr) cout << "CG " << k << "  " << tv << endl;
	} while ((tv > mmax)&&(++k < kn));
	cout << " CG " << k << "  " << tv << endl;
}

void Vector::BiCG(const Matrix &M, int kn, double mmax, bool wr) {
	int k = 0; double tv; Complex tc, tcc; Matrix MT(M); MT = M.mconj();
	Vector R1(*this), R2(*this), P1(*this), P2(*this), PP(*this);
	memset(Data,0,Nrow*sizeof(Complex));
	do {
		PP = P1*M; tc = (tcc = (R1*R2))/(PP*P2); R1 -= tc*PP;
		*this += tc*P1; tv = abs(tc)*P1.normF();
		PP = P2*MT; R2 -= tc*PP;
		tc = (R1*R2)/tcc; P1 = R1 + tc*P1; P2 = R2 + tc*P2;
		if (wr) cout << "BiCG " << k << "  " << tv << endl;
	} while ((tv > mmax)&&(++k < kn));
	cout << " BiCG " << k << "  " << tv << endl;
}

void Vector::BiCGs(const Matrix &M, int kn, double mmax, bool wr) {
	int k = 0; double tv; Complex ta, tw, tc; Matrix MT(M); MT = M.mconj();
	Vector R(*this), RR(*this), P(*this), SS(*this), PP(*this);
	memset(Data,0,Nrow*sizeof(Complex));
	tc = R*RR;
	do {
		PP = P*M; ta = tc/(PP*RR); *this += ta*P; tv = abs(ta)*P.normF();
		R -= ta*PP; SS = R*M; tw = (SS*R)/(SS*SS); *this += tw*R; tv += abs(tw)*R.normF();
		R -= tw*SS; ta /= tw*tc; ta *= (tc = (R*RR));
		P = R + ta*(P - tw*PP);
		if (wr) cout << "BiCGs " << k << "  " << tv << endl;
	} while ((tv > mmax)&&(++k < kn));
	cout << " BiCGs " << k << "  " << tv << endl;
}

void Vector::BiCGs(void *par, Vector (*mul)(const Vector&, void*), int kn, double mmax, bool wr) {
/**
	int k = 0; double tv; Complex ta, tw, tc;
	Vector R(*this), RR(*this), P(*this), SS(*this), PP(*this);
	memset(Data,0,Nrow*sizeof(Complex));
	tc = R*RR;
	do {
		PP = mul(P,par); //PP = P*M;
		ta = tc/(PP*RR); *this += ta*P; tv = abs(ta)*P.normF();
		R -= ta*PP;
		SS = mul(R,par); //SS = R*M;
		tw = (SS*R)/(SS*SS); *this += tw*R; tv += abs(tw)*R.normF();
		R -= tw*SS; ta /= tw*tc; ta *= (tc = (R*RR));
		P = R + ta*(P - tw*PP);
		if (wr) cout << "BiCGs " << k << "  " << tv << endl;
	} while ((tv > mmax)&&(++k < kn));
	cout << " BiCGs " << k << "  " << tv << endl;
/**/
	int k = 0; double tv; Complex ta, tw, tt, trh1, trh2, tb;
	Vector VP(*this), VR(*this), VRR(*this), VS(*this), VT(*this), VV(*this);
	memset(Data,0,Nrow*sizeof(Complex));
	memset(VV.Data,0,Nrow*sizeof(Complex));
	memset(VP.Data,0,Nrow*sizeof(Complex));
	trh1 = ta = tw = 1.;
	do {
		trh2 = VRR*VR;
		tb = trh2/trh1*ta/tw;
		VP = VR + tb*(VP - tw*VV);
		VV = mul(VP,par);
		ta = trh2/(VRR*VV);
		VS = VR - ta*VV;
		VT = mul(VS,par);
		tw = tt = 0.; for (int i=0; i<Nrow; ++i) {tw += VT(i)*VS(i); tt += VT(i)*VT(i);} tw /= tt;
		*this = *this + tw*VS + ta*VP;
		VR = VS - tw*VT;
		trh1 = trh2; tv = VR.normF();
		if (wr) cout << "BiCGs " << k << "  " << tv << endl;
	} while ((tv > mmax)&&(++k < kn));
}

void Vector::BiCGs(void *par, Vector (*mul)(const Vector&, void*), const Vector &V0, int kn, double mmax, bool wr) {
	int k = 0; double tv; Complex ta, tw, tt, trh1, trh2, tb;
	Vector VP(*this), VR(*this), VRR(*this), VS(*this), VT(*this), VV(*this);
	memset(Data,0,Nrow*sizeof(Complex));
	memset(VV.Data,0,Nrow*sizeof(Complex));
	memset(VP.Data,0,Nrow*sizeof(Complex));
	//cout<<Nrow<<" "<<V0.Nrow<<endl; cin.get();
	VRR = VR = VR - mul(V0,par);
	trh1 = ta = tw = 1.;
	do {
		trh2 = VRR*VR;
		tb = trh2/trh1*ta/tw;
		VP = VR + tb*(VP - tw*VV);
		VV = mul(VP,par);
		ta = trh2/(VRR*VV);
		VS = VR - ta*VV;
		VT = mul(VS,par);
		tw = tt = 0.; for (int i=0; i<Nrow; ++i) {tw += VT(i)*VS(i); tt += VT(i)*VT(i);} tw /= tt;
		*this = *this + tw*VS + ta*VP;
		VR = VS - tw*VT;
		trh1 = trh2; tv = VR.normF();
		if (wr) cout << "BiCGs " << k << "  " << tv << endl;
	} while ((tv > mmax)&&(++k < kn));
}

Vector Vector::calc_LSa(const Matrix &M) {
	if ((Nrow != M.Ncol) || (M.Ncol < M.Nrow)) return *this; // exception
	if (M.Ncol == M.Nrow) return *this/M;

	int i, j, N = M.Nrow;
	Matrix D(N);
	Vector V(N);

	for (i=0; i<N; i++) for (j=0; j<N; j++) {}

	return V;
}

Vector& Vector::operator /= (const Matrix& A) {
	if (A.Ncol != A.Nrow) {cout<<"Vector operator /= : matrix is not square\n"; return *this;}
	if (A.Nrow != Nrow) {cout<<"Vector operator /= : wrong matrix dimension\n"; return *this;}
	int i, ii;
	Matrix AA(Nrow); 
	for (i=0; i<Nrow; i++) for (ii=0; ii<Nrow; ii++) AA.Data[i*Nrow+ii] = A.Data[ii*Nrow+i]; 
	*this %= AA; return *this;
}

Vector& Vector::operator %= (const Matrix& A) {
	if (A.Ncol != A.Nrow) {cout<<"Vector operator %= : matrix is not square\n"; return *this;}
	if (A.Nrow != Nrow) {cout<<"Vector operator %= : wrong matrix dimension\n"; return *this;}

	double MaxMod, NewMod;
	const double NearZero = 1.e-15;
	Matrix R(A);
	int i, ii, j, jj, t = 2*Nrow;
	double *vva, *vvb, *vvc, *vv, va, vb, vc, vd, tv;
	Complex *a = R.Data, *b = Data, *aa, *bb, *aaa = R.Data + Nrow*Nrow - 1, ra; 

	for(i=0; i<Nrow-1; i++,a+=Nrow+1,b++) {
		MaxMod = abs(*a); j = i; vva = reinterpret_cast <double*>(a+Nrow);
		for (ii=i+1; ii<Nrow; ii++,vva+=2*Nrow) {
			va = *vva; vb = *(vva+1); NewMod = sqrt(va*va+vb*vb);
			if (MaxMod < NewMod) {MaxMod = NewMod; j = ii;}
		}
		if (MaxMod <= NearZero) {cout<<"Vector operator %= : matrix is zero\n"; return *this;}
		jj = 2*(Nrow-i);
		if (j > i) {
			vva = reinterpret_cast <double*>(a); vvc = reinterpret_cast <double*>(a+(j-i)*Nrow); vv = vvc+jj;
			while(vvc < vv) {tv = *vva; *vva++ = *vvc; *vvc++ = tv;}
			vva = reinterpret_cast <double*>(b); vvc = reinterpret_cast <double*>(b+(j-i)); 
			tv = *vva; *vva++ = *vvc; *vvc++ = tv; tv = *vva; *vva++ = *vvc; *vvc++ = tv;
		}
		jj -= 2; 
		//ra = -1./(*a);
		va = a->real(); vb = a->imag(); tv = va*va+vb*vb; ra = Complex(-va/tv,vb/tv);
		for (aa=a+Nrow,bb=b+1; aa<aaa; aa+=Nrow,bb++) {
			va = ra.real(); vb = ra.imag(); 
			tv = aa->real(); vd = aa->imag(); vc = va*tv-vb*vd; vd = va*vd+vb*tv;
			vva = reinterpret_cast <double*>(a+1); vvc = reinterpret_cast <double*>(aa+1); vv = vvc+jj;   	
			while(vvc < vv) {va = *vva++; vb = *vva++; *vvc++ += va*vc-vb*vd; *vvc++ += va*vd+vb*vc;}
			vva = reinterpret_cast <double*>(b); vvc = reinterpret_cast <double*>(bb);    	
			va = *vva++; vb = *vva; *vvc++ += va*vc-vb*vd; *vvc += va*vd+vb*vc;
		}
	}
	for (i=Nrow-1; i>=0; i--,a-=Nrow+1,b--) {
		//ra = 1./(*a); vc = ra.real(); vd = ra.imag(); 
		va = a->real(); vb = a->imag(); tv = va*va+vb*vb; vc = va/tv; vd = -vb/tv;
		vvc = reinterpret_cast <double*>(b); 	
		va = *vvc; vb = *(vvc+1); *vvc++ = va*vc-vb*vd; *vvc = va*vd+vb*vc;
		vvb = reinterpret_cast <double*>(R.Data+i);    
		for (bb=Data; bb<b; vvb+=2*Nrow,bb++) {
			vc = *vvb; vd = *(vvb+1); 
			vva = reinterpret_cast <double*>(b); vvc = reinterpret_cast <double*>(bb);   	
			va = *vva++; vb = *vva; *vvc++ -= va*vc-vb*vd; *vvc -= va*vd+vb*vc;
		}
	}
	return *this;
}

double Vector::normF(int nn) const {
	if ( nn >= Nrow ) return 0.;
	double tv = 0., tvr, tvi, *pv;
	pv = reinterpret_cast<double*>(this->Data) + 2*nn;
	for (int i=nn; i<Nrow; i++) {tvr = *pv++; tvi = *pv++; tv += tvr*tvr + tvi*tvi;}
	return sqrt(tv);
}

double Vector::cmp(const Vector &V) {
	if (Nrow != V.Nrow) return -1.;
	double tv = 0., tvv;
	for (int i=0; i<Nrow; i++) if (tv < (tvv = abs(V.Data[i] - Data[i]))) tv = tvv;
	return tv;
}

double Vector::cmpa(const Vector &V) {
	if (Nrow != V.Nrow) return -1.;
	double tv = 0., tvv;
	for (int i=0; i<Nrow; i++) if (tv < (tvv = fabs(abs(V.Data[i]) - abs(Data[i])))) tv = tvv;
	return tv;
}

void qsort(Complex *A, int iLo, int iHi) {
	int Lo, Hi; double Mid; Complex tc;
	Lo = iLo; Hi = iHi; Mid = A[(Lo+Hi)/2].real();
	do {
		while (A[Lo].real() < Mid) Lo++; while (A[Hi].real() > Mid) Hi--;
		if (Lo <= Hi) {tc = A[Lo]; A[Lo++] = A[Hi]; A[Hi--] = tc;}
	} while (Lo <= Hi);
	if (Hi > iLo) qsort(A,iLo,Hi); if (Lo < iHi) qsort(A,Lo,iHi);
}

void qsort(Complex *A, Matrix &ME, int iLo, int iHi) {
	int i, Lo, Hi; double Mid; Complex tc;
	Lo = iLo; Hi = iHi; Mid = A[(Lo+Hi)/2].real();
	do {
		while (A[Lo].real() < Mid) Lo++; while (A[Hi].real() > Mid) Hi--;
		if (Lo <= Hi) {
			tc = A[Lo]; A[Lo] = A[Hi]; A[Hi] = tc;
			for (i=0; i<ME.Nrow; ++i) {tc = ME(i,Lo); ME.Data[i*ME.Nrow+Lo] = ME(i,Hi); ME.Data[i*ME.Nrow+Hi] = tc;}
			Lo++; Hi--;
		}
	} while (Lo <= Hi);
	if (Hi > iLo) qsort(A,ME,iLo,Hi); if (Lo < iHi) qsort(A,ME,Lo,iHi);
}

void Vector::qsort() {
	int Lo, Hi; double Mid; Complex tc;
	Lo = 0; Hi = Nrow-1; Mid = Data[(Lo+Hi)/2].real();
	do {
		while (Data[Lo].real() < Mid) Lo++; while (Data[Hi].real() > Mid) Hi--;
		if (Lo <= Hi) {tc = Data[Lo]; Data[Lo++] = Data[Hi]; Data[Hi--] = tc;}
	} while (Lo <= Hi);
	if (Hi > 0) ::qsort(Data,0,Hi); if (Lo < Nrow-1) ::qsort(Data,Lo,Nrow-1);
}

void Vector::qsorta() {
	int Lo, Hi; double Mid; Complex tc;
	Lo = 0; Hi = Nrow-1; Mid = abs(Data[(Lo+Hi)/2]);
	do {
		while (abs(Data[Lo]) < Mid) Lo++; while (abs(Data[Hi]) > Mid) Hi--;
		if (Lo <= Hi) {tc = Data[Lo]; Data[Lo++] = Data[Hi]; Data[Hi--] = tc;}
	} while (Lo <= Hi);
	if (Hi > 0) ::qsort(Data,0,Hi); if (Lo < Nrow-1) ::qsort(Data,Lo,Nrow-1);
}

void Vector::qsort(Matrix &ME) {
	int i, Lo, Hi; double Mid; Complex tc;
	Lo = 0; Hi = Nrow-1; Mid = Data[(Lo+Hi)/2].real();
	do {
		while (Data[Lo].real() < Mid) Lo++; while (Data[Hi].real() > Mid) Hi--;
		if (Lo <= Hi) {
			tc = Data[Lo]; Data[Lo] = Data[Hi]; Data[Hi] = tc;
			for (i=0; i<Nrow; ++i) {tc = ME(i,Lo); ME.Data[i*Nrow+Lo] = ME(i,Hi); ME.Data[i*Nrow+Hi] = tc;}
			Lo++; Hi--;
		}
	} while (Lo <= Hi);
	if (Hi > 0) ::qsort(Data,ME,0,Hi); if (Lo < Nrow-1) ::qsort(Data,ME,Lo,Nrow-1);
}

Complex Vector::poly(Complex x) {
	Complex s;
	int i;
	if (abs(x) < 1) {s = x; for (i=Nrow-1; i>0; i--) {s += Data[i]; s *= x;} s += Data[0];}
	else {s = Data[0]/x; for (i=1; i<Nrow; i++) {s += Data[i]; s /= x;} s += 1;}
	return s;
}

Vector Vector::pcoeff() {
	int i, j;
	Vector A(Nrow), B(Nrow);
	for (j=1; j<Nrow; j++) A.Data[j] = 0.; A.Data[0] = 1.;
	for (i=0; i<Nrow; i++) {
		B.Data[0] = 0; for (j=1; j<Nrow; j++) B.Data[j] = A.Data[j-1];
		for (j=0; j<=i; j++) B.Data[j] -= A.Data[j]*Data[i];
		for (j=0; j<Nrow; j++) A.Data[j] = B.Data[j];
	}
	return A;
}

Vector Vector::roots() {
	Vector D(Nrow);
	if (Nrow == 1) {D.Data[0] = -Data[0]; return D;}
	int i, j;
	Matrix M(Nrow);
	for (i=0; i<Nrow; i++) for (j=0; j<Nrow; j++) M.Data[i*Nrow+j] = 0.;
	for (i=0; i<Nrow; i++) M.Data[i] = -Data[Nrow-1-i]; for (i=1; i<Nrow; i++) M.Data[i*Nrow+i-1] = 1.;
	D = M.diag(0); D.qsort(); 
	return D;
}

///////////////////////////////////////////////////////////////////////
// implementations of the class Matrix ////////////////////////////////
///////////////////////////////////////////////////////////////////////

Matrix::Matrix(unsigned Mown, unsigned Mext) {
	Nrow = Mown;
	Ncol = (Mext) ? Mext : Mown;
	Data = new Complex[Nrow*Ncol];
}

Matrix::Matrix(const Matrix& M) {
	Nrow = M.Nrow; Ncol = M.Ncol;
	Data = new Complex[Nrow*Ncol];
	memcpy(Data, M.Data, Nrow*Ncol*sizeof(Complex));
}

Matrix::~Matrix() {delete [] Data;}

double Matrix::normF(int n1) const {
	double tv = 0.;
	//for (int i=n1; i<Nrow*Ncol; i++) tv += abs(Data[i]*Data[i]);
	if (n1 > -1) for (int i=0; i<n1; i++) tv += abs(Data[i]*Data[i]);
	else for (int i=0; i<Nrow*Ncol; i++) tv += abs(Data[i]*Data[i]);
	return sqrt(tv);//
}

double Matrix::cmp(const Matrix &M) {
	if ((Nrow != M.Nrow) || (Ncol != M.Ncol)) return -1.;
	double tv = 0., tvv;
	for (int i=0; i<Nrow*Ncol; i++) if (tv < (tvv = abs(M.Data[i] - Data[i]))) tv = tvv;
	return tv;
}

double Matrix::cmpa(const Matrix &M) {
	if ((Nrow != M.Nrow) || (Ncol != M.Ncol)) return -1.;
	double tv = 0., tvv;
	for (int i=0; i<Nrow*Ncol; i++) if (tv < (tvv = fabs(abs(M.Data[i]) - abs(Data[i])))) tv = tvv;
	return tv;
}

void Matrix::print(const int pres) const {
	cout.precision(pres);
	for (int i=0; i<Nrow; ++i) {
		for (int ii=0; ii<Ncol; ++ii) cout << Data[i*Ncol+ii] << " ";
		cout << endl;
	}
}

void Matrix::eye(void) {
	int N = (Nrow < Ncol) ? Nrow : Ncol;
	memset(Data,0,Nrow*Ncol*sizeof(Complex));
	for (int i=0; i<N; ++i) Data[i*Ncol+i] = 1.;
}

Matrix& Matrix::operator = (const Matrix& M) {
	if ((Nrow != M.Nrow)||(Ncol != M.Ncol))
		{if (Nrow*Ncol) delete[] Data; Nrow = M.Nrow; Ncol = M.Ncol; Data = new Complex[Nrow*Ncol];}
	memcpy(Data, M.Data, Nrow*Ncol*sizeof(Complex));
	return *this;
}

Matrix Matrix::transp() const {
	int i, j; Matrix M(Ncol,Nrow);
	for (i=0; i<Nrow; i++) for (j=0; j<Ncol; j++) M.Data[j*Nrow+i] = Data[i*Ncol+j];
	return M;
}

Matrix Matrix::mconj() const {
	int i, j; Matrix M(Ncol,Nrow);
	for (i=0; i<Nrow; i++) for (j=0; j<Ncol; j++) M.Data[j*Nrow+i] = conj(Data[i*Ncol+j]);
	return M;
}

Complex Matrix::trace() const {
	if (Nrow != Ncol) return 0.;
	Complex tc = 0.;
	for (int i=0; i<Nrow; i++) tc += Data[i*Nrow+i];
	return tc;
}

Matrix& Matrix::operator += (const Matrix& M) {
  Complex *D1 = Data;
  Complex *D2 = M.Data;
  Complex *D1end = Data + Nrow*Ncol;
	if ((*this).Nrow != M.Nrow) {cout<<"Matrix operator += : wrong dimension\n"; return *this;}
  if ((*this).Ncol != M.Ncol) {cout<<"Matrix operator += : wrong dimension\n"; return *this;}
  while (D1 < D1end) *D1++ += *D2++;
  return *this;
}

Matrix& Matrix::operator -= (const Matrix& M) {
  Complex *D1 = Data;
  Complex *D2 = M.Data;
  Complex *D1end = Data + Nrow*Ncol;
  if ((*this).Nrow != M.Nrow) {cout<<"Matrix operator -= : wrong dimension\n"; return *this;}
  if ((*this).Ncol != M.Ncol) {cout<<"Matrix operator -= : wrong dimension\n"; return *this;}
  while (D1 < D1end) *D1++ -= *D2++;
  return *this;
}

Vector Matrix::operator * (const Vector& B) const {
	Vector C(Nrow);
	//if (Nrow != Ncol) {cout<<"Matrix operator * Vector : wrong dimension\n"; return C;}
	if (Ncol != B.Nrow) {cout<<"Matrix operator * Vector : wrong dimension\n"; return C;}

	for (int i=0; i<Nrow; ++i) {C.Data[i] = 0.; for (int j=0; j<Ncol; ++j) C.Data[i] += Data[i*Ncol+j]*B(j);}
	return C;

	unsigned tt = 2*Nrow;
	double *va, *vb, *vc, *vbb, *vcc, ta, taa, tb, tbb;
	va = reinterpret_cast <double*>(Data); vb = reinterpret_cast <double*>(B.Data); 
	vc = reinterpret_cast <double*>(C.Data); 
	for (vcc = vc+2*Nrow; vc<vcc; vc+=2) {
		*vc = *(vc+1) = 0.; 
		for (vbb=vb+2*Ncol; vb<vbb;) 
			{ta = *va++; taa = *va++; tb = *vb++; tbb = *vb++; *vc += ta*tb-taa*tbb; *(vc+1) += ta*tbb+taa*tb;}
		vb -= tt;
	}
	return C;
}

Matrix Matrix::operator * (const Matrix& B) const {
	Matrix C(Nrow, B.Ncol);
	if (Ncol != B.Nrow) {cout<<"Matrix operator * Matrix : wrong dimension\n"; return C;}; // exception ???
//if ((Nrow < 0x200)&&(Nrow%2 == 0)&&(Ncol%2 == 0)&&(B.Ncol%2 == 0))
	{
		unsigned tt = 2*Nrow*C.Ncol, tc = 2*C.Ncol, t = 2*Ncol*C.Ncol;
		double *va, *vb, *vc, *vaa, *vbb, *vcc, tva, tvaa, tvb, tvbb;
		va = reinterpret_cast <double*>(Data); vb = reinterpret_cast <double*>(B.Data); vc = reinterpret_cast <double*>(C.Data); 
		for (vcc = vc+tt; vc<vcc;) *vc++ = *vc++ = 0.; vcc = vc-tt+tc;
		for (vbb=va+2*Nrow*Ncol; va<vbb; vcc+=tc, vb-=t) for (vaa=va+2*Ncol; va<vaa;) {
			tva = *va++; tvaa = *va++; vc = vcc-tc; 
			while (vc < vcc) {
				tvb = *vb++; tvbb = *vb++; *vc++ += tva*tvb-tvaa*tvbb; *vc++ += tva*tvbb+tvaa*tvb;
				//*vc++ += tva*(tvb = *vb++)-tvaa*(tvbb = *(++vb)); *vc++ += tva*tvbb+tvaa*tvb;
				//tvb = tva*(tv = *vb++); tvbb = tvaa*(tvv = *vb++); *vc++ += tvb-tvbb; *vc++ += (tv+tva)*(tvv+tvaa)-tvb-tvbb;
			}
		}
	}
	return C;
	//else
	{  // Strassen
		int i, ii, n1 = Nrow/2, n2 = Ncol/2, n3 = B.Ncol/2; 
		Matrix S(n1,n2), SS(n2,n3), M1(n1,n3), M2(n1,n3), M3(n1,n3), M4(n1,n3), M5(n1,n3);
		for (ii=0; ii<n2; ii++) {
			for (i=0; i<n1; i++) S.Data[i*n2+ii] = Data[i*Ncol+ii]-Data[(i+n1)*Ncol+ii];
			for (i=0; i<n3; i++) SS.Data[ii*n3+i] = B.Data[(ii+n2)*B.Ncol+n3+i]-B.Data[ii*B.Ncol+n3+i];
		}     
		M4 = S*SS;
		for (ii=0; ii<n2; ii++) {
			for (i=0; i<n1; i++) S.Data[i*n2+ii] = Data[i*Ncol+n2+ii];
			for (i=0; i<n3; i++) SS.Data[ii*n3+i] = B.Data[(ii+n2)*B.Ncol+i];
		}     
		M3 = S*SS;
		for (ii=0; ii<n2; ii++) {
			for (i=0; i<n1; i++) S.Data[i*n2+ii] = Data[i*Ncol+ii];
			for (i=0; i<n3; i++) SS.Data[ii*n3+i] = B.Data[ii*B.Ncol+i];
		}     
		M2 = S*SS;
		for (ii=0; ii<n2; ii++) {
			for (i=0; i<n1; i++) S.Data[i*n2+ii] = Data[(i+n1)*Ncol+ii]+Data[(i+n1)*Ncol+n2+ii];
			for (i=0; i<n3; i++) SS.Data[ii*n3+i] = B.Data[ii*B.Ncol+n3+i]-B.Data[ii*B.Ncol+i];
		}     
		M5 = S*SS;
		for (ii=0; ii<n2; ii++) {
			for (i=0; i<n1; i++) S.Data[i*n2+ii] -= Data[i*Ncol+ii];
			for (i=0; i<n3; i++) SS.Data[ii*n3+i] = B.Data[(ii+n2)*B.Ncol+n3+i]-SS.Data[ii*n3+i];
		}     
		M1 = S*SS+M2; M4 += M1; M1 += M5;
		for (i=0; i<n1; i++) for (ii=0; ii<n3; ii++) {
			C.Data[i*C.Ncol+ii] = M2.Data[i*n3+ii]+M3.Data[i*n3+ii];
			C.Data[(i+n1)*C.Ncol+n3+ii] = M4.Data[i*n3+ii]+M5.Data[i*n3+ii];
		}
		M5 = SS;
		for (ii=0; ii<n2; ii++) {
			for (i=0; i<n1; i++) S.Data[i*n2+ii] = Data[i*Ncol+n2+ii]-S.Data[i*n2+ii];
			for (i=0; i<n3; i++) SS.Data[ii*n3+i] = B.Data[(ii+n2)*B.Ncol+n3+i];
		}     
		M2 = S*SS;
		for (ii=0; ii<n2; ii++) {
			for (i=0; i<n1; i++) S.Data[i*n2+ii] = Data[(i+n1)*Ncol+n2+ii];
			for (i=0; i<n3; i++) SS.Data[ii*n3+i] = M5.Data[ii*n3+i]-B.Data[(ii+n2)*B.Ncol+i];
		}     
		M3 = S*SS;
		for (i=0; i<n1; i++) for (ii=0; ii<n3; ii++) {
			C.Data[i*C.Ncol+n3+ii] = M1.Data[i*n3+ii]+M2.Data[i*n3+ii];
			C.Data[(i+n1)*C.Ncol+ii] = M4.Data[i*n3+ii]-M3.Data[i*n3+ii];
		}
	}
	return C;
}

Matrix Matrix::operator * (const Complex d) {
	Matrix C(*this);
	for (int i=0; i<Nrow*Ncol; i++) C.Data[i] *= d;
	return C;
}

Matrix& Matrix::operator *= (const Complex d) {
	for (int i=0; i<Nrow*Ncol; i++) this->Data[i] *= d;
	return *this;
}

Matrix& Matrix::operator /= (const Matrix& A) {
	if (this == &A) {cout<<"Matrix operator /= : wrong divider\n"; return *this;}
	if (A.Ncol != A.Nrow) {cout<<"Matrix operator /= : divider is not square\n"; return *this;}
	if (A.Nrow != Ncol) {cout<<"Matrix operator /= : wrong divider dimension\n"; return *this;}
	int i, ii;
	Matrix R(Ncol,Nrow), AA(Ncol);
	for (i=0; i<Ncol; i++) for (ii=0; ii<Nrow; ii++) R.Data[i*Nrow+ii] = Data[ii*Ncol+i];
	for (i=0; i<Ncol; i++) for (ii=0; ii<Ncol; ii++) AA.Data[i*Ncol+ii] = A.Data[ii*Ncol+i]; 
	R %= AA;
	for (i=0; i<Nrow; i++) for (ii=0; ii<Ncol; ii++) Data[i*Ncol+ii] = R.Data[ii*Nrow+i];
	return *this;
}

Matrix& Matrix::operator %= (const Matrix& A) {
  if (this == &A) {cout<<"Matrix operator /= : wrong argument\n"; return *this;}
  if (A.Ncol != A.Nrow) {cout<<"Matrix operator /= : argument is not square\n"; return *this;}
  if (A.Nrow != Nrow) {cout<<"Matrix operator /= : wrong argument dimension\n"; return *this;}
	//if (Nrow < 0x200)
	{
		double MaxMod, NewMod;
		const double NearZero = 1.e-15;
		Matrix R(A);
		int i, ii, j, jj, t = 2*Ncol;
		double *vva, *vvb, *vvc, *vv, va, vb, vc, vd, tv;
		Complex *a = R.Data, *b = Data, *aa, *bb, *aaa = R.Data+Nrow*Nrow-1, ra; 

		for(i=0; i<Nrow-1; i++,a+=Nrow+1,b+=Ncol) {
			MaxMod = abs(*a); j = i; vva = reinterpret_cast <double*>(a+Nrow);
			for (ii=i+1; ii<Nrow; ii++,vva+=2*Nrow) {
				va = *vva; vb = *(vva+1); NewMod = sqrt(va*va+vb*vb);
				if (MaxMod < NewMod) {MaxMod = NewMod; j = ii;}
			}
			if (MaxMod <= NearZero) return *this; // exception ???
			jj = 2*(Nrow-i);
			if (j > i) {
				vva = reinterpret_cast <double*>(a); vvc = reinterpret_cast <double*>(a+(j-i)*Nrow); vv = vvc+jj;
				while(vvc < vv) {tv = *vva; *vva++ = *vvc; *vvc++ = tv;}
				vva = reinterpret_cast <double*>(b); vvc = reinterpret_cast <double*>(b+(j-i)*Ncol); vv = vvc+t;
				while(vvc < vv) {tv = *vva; *vva++ = *vvc; *vvc++ = tv;}
			}
			jj -= 2; 
			//ra = -1./(*a);
			va = a->real(); vb = a->imag(); tv = va*va+vb*vb; ra = Complex(-va/tv,vb/tv);
			for (aa=a+Nrow,bb=b+Ncol; aa<aaa; aa+=Nrow,bb+=Ncol) {
				va = ra.real(); vb = ra.imag(); 
				tv = aa->real(); vd = aa->imag(); vc = va*tv-vb*vd; vd = va*vd+vb*tv;
				vva = reinterpret_cast <double*>(a+1); vvc = reinterpret_cast <double*>(aa+1); vv = vvc+jj;   	
				while(vvc < vv) {va = *vva++; vb = *vva++; *vvc++ += va*vc-vb*vd; *vvc++ += va*vd+vb*vc;}
				vva = reinterpret_cast <double*>(b); vvc = reinterpret_cast <double*>(bb); vv = vvc+t;   	
				while(vvc < vv) {va = *vva++; vb = *vva++; *vvc++ += va*vc-vb*vd; *vvc++ += va*vd+vb*vc;}
			}
		}
		for (i=Nrow-1; i>=0; i--,a-=Nrow+1,b-=Ncol) {
			//ra = 1./(*a); vc = ra.real(); vd = ra.imag(); 
			va = a->real(); vb = a->imag(); tv = va*va+vb*vb; vc = va/tv; vd = -vb/tv;
			vvc = reinterpret_cast <double*>(b); vv = vvc+t;   	
			while(vvc < vv) {va = *vvc; vb = *(vvc+1); *vvc++ = va*vc-vb*vd; *vvc++ = va*vd+vb*vc;}
			vvb = reinterpret_cast <double*>(R.Data+i);    
			for (bb=Data; bb<b; vvb+=2*Nrow,bb+=Ncol) {
				vc = *vvb; vd = *(vvb+1); 
				vva = reinterpret_cast <double*>(b); vvc = reinterpret_cast <double*>(bb); vv = vvc+t;   	
				while(vvc < vv) {va = *vva++; vb = *vva++; *vvc++ -= va*vc-vb*vd; *vvc++ -= va*vd+vb*vc;}
			}
		}
		return *this;
	}
}

Matrix Matrix::calc_QRHess(Vector &VC) { // QR decomposition of Hessenberg matrix
	int i, ii; double tt;
	Complex tc, ts, tv, *h1, *h2, *hh1, *hh2;
	Matrix R(*this);
	if (Nrow != Ncol) return R;
	h2 = (h1 = R.Data) + Ncol;
	for (i=0; i<Nrow-1; ++i,h1+=(Ncol+1),h2+=(Ncol+1)) {
		if (abs(*h2) < 1.e-15) {tc = 1.; ts = 0.;}
		else if (abs(*h1) < 1.e-15) {tc = 0.; ts = conj(*h2)/abs(*h2);}
		else {
			tt = 1./sqrt(abs(*h1 * *h1) + abs(*h2 * *h2));
			tc = tt*abs(*h1); ts = tt*conj(*h2)*(*h1)/abs(*h1);
		}
		VC.Data[2*i] = tc; VC.Data[2*i+1] = ts;
		for (ii=i,hh1=h1,hh2=h2; ii<Ncol; ++ii,hh1++,hh2++) {
			*hh1 = tc*(tv = *hh1) + ts* *hh2;
			*hh2 = -conj(ts)*tv + tc* *hh2;
		}
	}
	return R;
}

Matrix Matrix::calc_SchurHess(Matrix &S) { // Schur desomposition of upper Hessenberg matrix
	if (Nrow != Ncol) return *this;
	int i, k; double tn; Complex tc;
	Vector VG(2*(Nrow-1));
	Matrix Q(*this), R(*this), M(*this), G(*this);
	S.eye();
	do {
		cout<<"QR iteration:\n";
		R = M.calc_QRHess(VG);
/**
		for (i=0; i<Nrow-1; ++i) {
			G.eye();
			G.Data[i*Nrow+i] = G.Data[(i+1)*Nrow+i+1] = VG(2*i);
			G.Data[(i+1)*Nrow+i] = -conj(G.Data[i*Nrow+i+1] = -VG(2*i+1));
			R *= G;
		}
/**/
//		for (i=0; i<Nrow; ++i) {for (int j=0; j<Nrow; ++j) printf("%2.5f ",R(i,j).real()); cout<<endl;} cout<<endl;
/**/
		for (i=0; i<Nrow-1; ++i) {
			for (k=0; k<i+1; ++k) {
				tc = R(k,i);
				R.Data[k*Ncol+i] = tc*VG(2*i) - R(k,i+1)*VG(2*i+1);
				R.Data[k*Ncol+i+1] = tc*conj(VG(2*i+1)) + R(k,i+1)*VG(2*i);
			}
		}
/**/
		M = R;
		for (i=0; i<Nrow; ++i) {for (int j=0; j<Nrow; ++j) printf("%2.5f ",M(i,j).real()); cout<<endl;} cout<<endl;
		tn = 0.; for (i=1; i<Nrow; ++i) for (k=0; k<i; ++k) tn += abs(M(i,k)*M(i,k)); tn = sqrt(tn);
		cout<<tn<<endl; cin.get();
	} while (tn > 1.e-15);
	return M;
}

Vector Matrix::PVector(Complex la) {
	int i;
	Matrix A(*this); 
	Vector V(Nrow);
	double tv;
	for (i=0; i<Nrow; i++) {A.Data[i*(Nrow+1)] -= la; V.Data[i] = 1.;}
	V %= A; tv = 0.;
	for (i=0; i<Nrow; i++) tv += abs(V.Data[i]*V.Data[i]); tv = sqrt(tv);
	for (i=0; i<Nrow; i++) V.Data[i] /= tv;
	return V;
}

Vector Matrix::LVector(Complex la) {
	int i;
	Matrix A(*this); 
	Vector V(Nrow);
	double tv;
	for (i=0; i<Nrow; i++) {A.Data[i*(Nrow+1)] -= la; V.Data[i] = 1.;}
	V /= A; tv = 0.;
	for (i=0; i<Nrow; i++) tv += abs(V.Data[i]*V.Data[i]); tv = sqrt(tv);
	for (i=0; i<Nrow; i++) V.Data[i] /= tv;
	return V;
}

Vector Matrix::diag(Matrix& V, int erm) const {
	Vector D(Nrow);
	Complex *M, *A, *DD;
	M = new Complex[Nrow*Nrow]; A = new Complex[Nrow*Nrow]; 
	int i, j;
	bool rk;
	for (i=0; i<Nrow; i++) for (j=0; j<Nrow; j++) A[i*Nrow+j] = Data[j*Nrow+i]; 
//	rk = 0; j = Nrow; zeigValVecA(D.Data,M,A,&j,&rk);
	for (i=0; i<Nrow; i++) for (j=0; j<Nrow; j++) V.Data[i*Nrow+j] = M[j*Nrow+i];
	delete [] A; delete [] M; 
	return D;
}

Vector Matrix::diag(int erm) const {
	Matrix M(Nrow); Vector D(Nrow);
	int i, j;
	bool rk;
	for (i=0; i<Nrow; i++) for (j=0; j<Nrow; j++) M.Data[i*Nrow+j] = Data[j*Nrow+i]; 
//	rk = 0; j = Nrow; zeigValA(D.Data,M.Data,&j,&rk);
	return D;
}

Matrix Matrix::inv2(int it) {
	int i, ii, n = Nrow/2;
	Matrix A(n), B(n), C(n), D(n), M(*this), N(n);
	
	memset(N.Data,0,n*n*sizeof(Complex));
	for (i=0; i<n; ++i) {
		memcpy(A.Data+n*i,Data+Nrow*i,n*sizeof(Complex));
		memcpy(B.Data+n*i,Data+Nrow*i+n,n*sizeof(Complex));
		memcpy(C.Data+n*i,Data+Nrow*(i+n),n*sizeof(Complex));
		memcpy(D.Data+n*i,Data+Nrow*(i+n)+n,n*sizeof(Complex));
		N.Data[i*n+i] = 1.;
	}
	A = N/A; D -= C*A*B; D = N/D; C = N = D*C*A; N = A*B*N; B = A*B*D; A += N;
	for (i=0; i<n; ++i) {
		for (ii=0; ii<n; ++ii) {B.Data[i*n+ii] = -B.Data[i*n+ii]; C.Data[i*n+ii] = -C.Data[i*n+ii];}
		memcpy(M.Data+Nrow*i,A.Data+n*i,n*sizeof(Complex));
		memcpy(M.Data+Nrow*i+n,B.Data+n*i,n*sizeof(Complex));
		memcpy(M.Data+Nrow*(i+n),C.Data+n*i,n*sizeof(Complex));
		memcpy(M.Data+Nrow*(i+n)+n,D.Data+n*i,n*sizeof(Complex));
	}
	return M;
}

// SMatrix implementations //////////////////////////////////////////////////////////////////////////////////////

SMatrix::SMatrix(unsigned n) {
	Mown = n;
	Data[0][0] = new Matrix(Mown); Data[0][1] = new Matrix(Mown);
	Data[1][0] = new Matrix(Mown); Data[1][1] = new Matrix(Mown);
}

SMatrix::SMatrix(const SMatrix& SM) {
	Mown = SM.Mown;
	Data[0][0] = new Matrix(Mown); *(Data[0][0]) = SM(0,0); Data[0][1] = new Matrix(Mown); *(Data[0][1]) = SM(0,1);
	Data[1][0] = new Matrix(Mown); *(Data[1][0]) = SM(1,0); Data[1][1] = new Matrix(Mown); *(Data[1][1]) = SM(1,1);
}

SMatrix::SMatrix(const SMatrixD& SD) {
	int i, j;
	Mown = SD.Mown;
	Data[0][0] = new Matrix(Mown); Data[0][1] = new Matrix(Mown);
	Data[1][0] = new Matrix(Mown); Data[1][1] = new Matrix(Mown);
	for (i=0; i<Mown; i++) {
		for (j=0; j<Mown; j++)
			Data[0][0]->Data[i*Mown+j] = Data[1][0]->Data[i*Mown+j] =
				Data[0][1]->Data[i*Mown+j] = Data[1][1]->Data[i*Mown+j] = 0;
		Data[0][0]->Data[i*Mown+i] = SD(0,i,0); Data[0][1]->Data[i*Mown+i] = SD(0,i,1);
		Data[1][0]->Data[i*Mown+i] = SD(1,i,0); Data[1][1]->Data[i*Mown+i] = SD(1,i,1);
	}
}

SMatrix::~SMatrix() {
  delete Data[0][0]; delete Data[1][0]; delete Data[0][1]; delete Data[1][1];
}

SMatrix& SMatrix::operator = (const SMatrix& SM) {
	Mown = SM.Mown;
	*(Data[0][0]) = SM(0,0); *(Data[0][1]) = SM(0,1);
	*(Data[1][0]) = SM(1,0); *(Data[1][1]) = SM(1,1);
	return *this;
}

SMatrix& SMatrix::operator = (const SMatrixD& SD) {
	int i, j;
	if (Mown != SD.Mown) {
		Mown = SD.Mown;
		delete Data[0][0]; delete Data[0][1]; delete Data[1][0]; delete Data[1][1];
		Data[0][0] = new Matrix(Mown); Data[0][1] = new Matrix(Mown);
		Data[1][0] = new Matrix(Mown); Data[1][1] = new Matrix(Mown);
	}
	for (i=0; i<Mown; i++) {
		for (j=0; j<Mown; j++)
			Data[0][0]->Data[i*Mown+j] = Data[1][0]->Data[i*Mown+j] =
				Data[0][1]->Data[i*Mown+j] = Data[1][1]->Data[i*Mown+j] = 0;
		Data[0][0]->Data[i*Mown+i] = SD(0,i,0); Data[0][1]->Data[i*Mown+i] = SD(0,i,1);
		Data[1][0]->Data[i*Mown+i] = SD(1,i,0); Data[1][1]->Data[i*Mown+i] = SD(1,i,1);
	}
	return *this;
}

void SMatrix::cut(const SMatrix& B) {
	int i, ii, j, jj, s;
	if (Mown == B.Mown) {*this = B; return;}
	if (Mown > B.Mown) for (j=0; j<2; j++) for (jj=0; jj<2; jj++) {
		s = Mown/2-B.Mown/2;
		for (i=0; i<Mown; i++) for (ii=0; ii<Mown; ii++) Data[j][jj]->Data[i*Mown+ii] = 0;
		for (i=0; i<B.Mown; i++) for (ii=0; ii<B.Mown; ii++) Data[j][jj]->Data[(i+s)*Mown+ii+s] = B.Data[j][jj]->Data[i*B.Mown+ii];
	}
	else for (j=0; j<2; j++) for (jj=0; jj<2; jj++)  {
		s = B.Mown/2-Mown/2;
		for (i=0; i<Mown; i++) for (ii=0; ii<Mown; ii++) Data[j][jj]->Data[i*Mown+ii] = B.Data[j][jj]->Data[(i+s)*B.Mown+ii+s];
	}
}

SMatrix SMatrix::operator * (const SMatrix& B) {
	SMatrix C(Mown);
	if (Mown != B.Mown) return C;
	int i, ii;
	Matrix BB(Mown);
/**
		// matrix by vector multiplication
	BB = B(0,0)*(*this)(1,1);
	for (i=0; i<Mown; ++i) {for (ii=0; ii<Mown; ++ii) BB.Data[i*Mown+ii] = -BB(i,ii); BB.Data[i*Mown+i] += 1.;}
	BB = (*this)(0,1)/BB; *(C.Data[0][1]) = BB*B(0,1); *(C.Data[0][0]) = BB*B(0,0)*(*this)(1,0) + (*this)(0,0);

	BB = (*this)(1,1)*B(0,0);
	for (i=0; i<Mown; ++i) {for (ii=0; ii<Mown; ++ii) BB.Data[i*Mown+ii] = -BB(i,ii); BB.Data[i*Mown+i] += 1.;}
	BB = B(1,0)/BB; *(C.Data[1][0]) = BB*(*this)(1,0); *(C.Data[1][1]) = B(1,1) + BB*(*this)(1,1)*B(0,1);
/**/
		// vector by matrix multiplication
	BB = (*this)(1,1)*B(0,0);
	for (i=0; i<Mown; ++i) {for (ii=0; ii<Mown; ++ii) BB.Data[i*Mown+ii] = -BB(i,ii); BB.Data[i*Mown+i] += 1.;}
	BB = (*this)(1,0)%BB; *(C.Data[1][0]) = B(1,0)*BB; *(C.Data[0][0]) = (*this)(0,0) + (*this)(0,1)*B(0,0)*BB;

	BB = B(0,0)*(*this)(1,1);
	for (i=0; i<Mown; ++i) {for (ii=0; ii<Mown; ++ii) BB.Data[i*Mown+ii] = -BB(i,ii); BB.Data[i*Mown+i] += 1.;}
	BB = B(0,1)%BB; *(C.Data[0][1]) = (*this)(0,1)*BB; *(C.Data[1][1]) = B(1,1) + B(1,0)*(*this)(1,1)*BB;

	return C;
}
	// incidence from outside to amplitudes between layers
SMatrix SMatrix::operator & (const SMatrix& B) {
	if (Mown != B.Mown) return *this;
	int i, j;
	Matrix AA(Mown), E(Mown);
	SMatrix C(Mown);
		// matrix by vector multiplication; (c0+,c1-) -> (+,-)
	AA = B(0,0)*(*this)(1,1);
	for (i=0; i<Mown; i++) {
		for (j=0; j<Mown; j++) {AA.Data[i*Mown+j] = -AA(i,j); E.Data[i*Mown+j] = 0.;}
		AA.Data[i*Mown+i] += 1.; E.Data[i*Mown+i] = 1.;
	}
	*(C.Data[0][0]) = *(C.Data[1][0]) = E/AA; *(C.Data[1][0]) = B(1,0)*(*this)(1,1)*C(1,0); *(C.Data[0][0]) = (*this)(0,1)*C(0,0);

	AA = (*this)(1,1)*B(0,0);
	for (i=0; i<Mown; i++) {
		for (j=0; j<Mown; j++) {AA.Data[i*Mown+j] = -AA(i,j); E.Data[i*Mown+j] = 0.;}
		AA.Data[i*Mown+i] += 1.; E.Data[i*Mown+i] = 1.;
	}
	*(C.Data[0][1]) = *(C.Data[1][1]) = E/AA; *(C.Data[0][1]) = (*this)(0,1)*B(0,0)*C(0,1); *(C.Data[1][1]) = B(1,0)*C(1,1);

	return C;
}
	// effective source amplitudes between layers
SMatrix SMatrix::operator % (const SMatrix& B) {
	if (Mown != B.Mown) return *this;
	int i, j;
	Matrix AA(Mown), E(Mown);
	SMatrix C(Mown);
		// vector by matrix multiplication, 0=="+", 1=="-"
	AA = (*this)(1,1)*B(0,0);
	for (i=0; i<Mown; i++) {
		for (j=0; j<Mown; j++) {AA.Data[i*Mown+j] = -AA(i,j); E.Data[i*Mown+j] = 0.;}
		AA.Data[i*Mown+i] += 1.; E.Data[i*Mown+i] = 1.;
	}
	*(C.Data[0][1]) = B(0,0)*(*(C.Data[1][1]) = E/AA); // C(1,1):r- -> a-; C(0,1): r+ -> a-
	AA = B(0,0)*(*this)(1,1);
	for (i=0; i<Mown; i++) {
		for (j=0; j<Mown; j++) {AA.Data[i*Mown+j] = -AA(i,j); E.Data[i*Mown+j] = 0.;}
		AA.Data[i*Mown+i] += 1.; E.Data[i*Mown+i] = 1.;
	}
	*(C.Data[1][0]) = (*this)(1,1)*(*(C.Data[0][0]) = E/AA); // C(0,0):r+ -> a+; C(0,1): r- -> a+

	return C;
}
	// source between layers emitting to outside
SMatrix SMatrix::operator ^ (const SMatrix& B) {
	if (Mown != B.Mown) return *this;
	int i, j;
	Matrix M(Mown);
	SMatrix S(Mown);
		// vector by matrix multiplication, (r+,r-) -> (a0-,a1+)
	M = (*this)(1,1)*B(0,0);
	for (i=0; i<Mown; i++) {for (j=0; j<Mown; j++) M.Data[i*Mown+j] = -M(i,j); M.Data[i*Mown+i] += 1.;}
	*(S.Data[0][0]) = B(0,0)*(*(S.Data[1][0]) = (*this)(1,0)%M);

	M = B(0,0)*(*this)(1,1);
	for (i=0; i<Mown; i++) {for (j=0; j<Mown; j++) M.Data[i*Mown+j] = -M(i,j); M.Data[i*Mown+i] += 1.;}
	*(S.Data[1][1]) = (*this)(1,1)*(*(S.Data[0][1]) = B(0,1)%M);

	return S;
}

SMatrix SMatrix::operator * (const SMatrixD& B) {
	if (Mown != B.Mown) return *this;
	int i, j;
	SMatrix C(Mown);
	Matrix AA(Mown), BB(Mown);
		// vector by matrix multiplication
	for (i=0; i<Mown; i++) {for (j=0; j<Mown; j++) AA.Data[i*Mown+j] = -(*this)(1,i,1,j)*B(0,j,0); AA.Data[i*Mown+i] += 1.;}
	BB = (*this)(1,0)%AA;
	for (i=0; i<Mown; i++) for (j=0; j<Mown; j++) {
		C.Data[0][0]->Data[i*Mown+j] = B(0,i,0)*BB(i,j); C.Data[1][0]->Data[i*Mown+j] = B(1,i,0)*BB(i,j);
	}
	*(C.Data[0][0]) = (*this)(0,0) + (*this)(0,1)*C(0,0);

	for (i=0; i<Mown; i++) {
		for (j=0; j<Mown; j++) {AA.Data[i*Mown+j] = -B(0,i,0)*(*this)(1,i,1,j); BB.Data[i*Mown+j] = 0.;}
		AA.Data[i*Mown+i] += 1.; BB.Data[i*Mown+i] = B(0,i,1);
	}
	*C.Data[1][1] = (*this)(1,1)*( *C.Data[0][1] = BB%AA ); *C.Data[0][1] = (*this)(0,1)*C(0,1);
	for (i=0; i<Mown; i++) {
		for (j=0; j<Mown; j++) C.Data[1][1]->Data[i*Mown+j] *= B(1,i,0);
		C.Data[1][1]->Data[i*Mown+i] += B(1,i,1);
	}

	return C;
}

SMatrix SMatrix::operator & (const SMatrixD& B) {
	if (Mown != B.Mown) return *this;
	int i, j;
	Matrix AA(Mown), E(Mown); SMatrix C(Mown);

	for (i=0; i<Mown; i++) {
		for (j=0; j<Mown; j++) {AA.Data[i*Mown+j] = -B(0,i,0)*(*this)(1,i,1,j); E.Data[i*Mown+j] = 0.;}
		AA.Data[i*Mown+i] += 1.; E.Data[i*Mown+i] = 1.;
	}
	*C.Data[1][0] = (*this)(1,1)*( *C.Data[0][0] = E/AA ); *C.Data[0][0] = (*this)(0,1)*C(0,0);
	for (i=0; i<Mown; i++) for (j=0; j<Mown; j++) C.Data[1][0]->Data[i*Mown+j] *= B(1,i,0);

	for (i=0; i<Mown; i++) {
		for (j=0; j<Mown; j++) {
			AA.Data[i*Mown+j] = -(*this)(1,i,1,j)*B(0,j,0);
			C.Data[1][1]->Data[i*Mown+j] = C.Data[0][1]->Data[i*Mown+j] = E.Data[i*Mown+j] = 0.;
		}
		AA.Data[i*Mown+i] += 1.; E.Data[i*Mown+i] = 1.;
		C.Data[0][1]->Data[i*Mown+i] = B(0,i,0); C.Data[1][1]->Data[i*Mown+i] = B(1,i,0);
	}
	E /= AA; *C.Data[1][1] *= E; *C.Data[0][1] = (*this)(0,1)*C(0,1)*E;

	return C;
}

SMatrix SMatrix::operator % (const SMatrixD& B) {
	if (Mown != B.Mown) return *this;
	int i, j; SMatrix C(Mown); Matrix AA(Mown), E(Mown);

	for (i=0; i<Mown; i++) {
		for (j=0; j<Mown; j++) {AA.Data[i*Mown+j] = -B(0,i,0)*(*this)(1,i,1,j); E.Data[i*Mown+j] = 0.;}
		AA.Data[i*Mown+i] += 1.; E.Data[i*Mown+i] = 1.;
	}
	*(C.Data[1][0]) = (*this)(1,1)*( *C.Data[0][0] = E/AA );

	for (i=0; i<Mown; i++) {
		for (j=0; j<Mown; j++) {AA.Data[i*Mown+j] = -(*this)(1,i,1,j)*B(0,j,0); E.Data[i*Mown+j] = 0.;}
		AA.Data[i*Mown+i] += 1.; E.Data[i*Mown+i] = 1.;
	}
	*C.Data[0][1] = *C.Data[1][1] = E/AA;
	for (i=0; i<Mown; i++) for (j=0; j<Mown; j++) C.Data[0][1]->Data[i*Mown+j] *= B(0,i,0);

	return C;
}

SMatrix SMatrix::operator ^ (const SMatrixD& B) {
	if (Mown != B.Mown) return *this;
	int i, j;
	SMatrix C(Mown);
	Matrix AA(Mown);

	for (i=0; i<Mown; i++) {for (j=0; j<Mown; j++) AA.Data[i*Mown+j] = -(*this)(1,i,1,j)*B(0,j,0); AA.Data[i*Mown+i] += 1;}
	*(C.Data[0][0]) = *(C.Data[1][0]) = (*this)(1,0)%AA;
	for (i=0; i<Mown; i++) for (j=0; j<Mown; j++) C.Data[0][0]->Data[i*Mown+j] *= B(0,i,0);

	for (i=0; i<Mown; i++) {
		for (j=0; j<Mown; j++) {AA.Data[i*Mown+j] = -B(0,i,0)*(*this)(1,i,1,j); C.Data[0][1]->Data[i*Mown+j] = 0.;}
		AA.Data[i*Mown+i] += 1.; C.Data[0][1]->Data[i*Mown+i] = B(0,i,1);
	}
	*C.Data[1][1] = (*this)(1,1)*( *C.Data[0][1] %= AA );

	return C;
}

SMatrix SMatrix::operator * (const double d) {
  SMatrix C(Mown);
  int i, ii;
  for (i=0; i<Mown; i++)
    for (ii=0; ii<Mown; ii++) {
      C.Data[1][1]->Data[i*Mown+ii] = Data[1][1]->Data[i*Mown+ii]*exp(j_*double(ii-i)*d);
      C.Data[1][0]->Data[i*Mown+ii] = Data[1][0]->Data[i*Mown+ii]*exp(j_*double(ii-i)*d);
      C.Data[0][1]->Data[i*Mown+ii] = Data[0][1]->Data[i*Mown+ii]*exp(j_*double(ii-i)*d);
      C.Data[0][0]->Data[i*Mown+ii] = Data[0][0]->Data[i*Mown+ii]*exp(j_*double(ii-i)*d);
		}
  return C;
}

SMatrix& SMatrix::operator *= (const double d) {
  int i, ii;
  for (i=0; i<Mown; i++)
    for (ii=0; ii<Mown; ii++) {
      Data[1][1]->Data[i*Mown+ii] *= exp(j_*double(ii-i)*d);
      Data[1][0]->Data[i*Mown+ii] *= exp(j_*double(ii-i)*d);
      Data[0][1]->Data[i*Mown+ii] *= exp(j_*double(ii-i)*d);
      Data[0][0]->Data[i*Mown+ii] *= exp(j_*double(ii-i)*d);
		}
  return *this;
}

SMatrix SMatrix::operator * (const SDMatrix& B) {
	if (Mown != 2*B.Mown) return *this;
	int i, j, Mo = B.Mown;
	Matrix M1(2*Mo), M2(2*Mo);
	SMatrix S(Mown);

	for (j=0; j<Mo; j++) {
		for (i=0; i<Mo; i++) {
			M1.Data[i*2*Mo+j] = - B(0,i,0)*(*this)(1,i,1,j) - B(0,i+Mo,0)*(*this)(1,i+Mo,1,j);
			M1.Data[i*2*Mo+j+Mo] = - B(0,i,0)*(*this)(1,i,1,j+Mo) - B(0,i+Mo,0)*(*this)(1,i+Mo,1,j+Mo);
			M1.Data[(i+Mo)*2*Mo+j] = - B(0,i+2*Mo,0)*(*this)(1,i,1,j) - B(0,i+3*Mo,0)*(*this)(1,i+Mo,1,j);
			M1.Data[(i+Mo)*2*Mo+j+Mo] = - B(0,i+2*Mo,0)*(*this)(1,i,1,j+Mo) - B(0,i+3*Mo,0)*(*this)(1,i+Mo,1,j+Mo);
		}
		M1.Data[j*2*Mo+j] += 1.; M1.Data[(j+Mo)*2*Mo+j+Mo] += 1.;
	}
	M2 = (*this)(0,1)/M1;
	for (j=0; j<Mo; j++) for (i=0; i<Mo; i++) {
		S.Data[0][1]->Data[i*2*Mo+j] = M2(i,j)*B(0,j,1) + M2(i,j+Mo)*B(0,j+2*Mo,1);
		S.Data[0][1]->Data[i*2*Mo+j+Mo] = M2(i,j)*B(0,j+Mo,1) + M2(i,j+Mo)*B(0,j+3*Mo,1);
		S.Data[0][1]->Data[(i+Mo)*2*Mo+j] = M2(i+Mo,j)*B(0,j,1) + M2(i+Mo,j+Mo)*B(0,j+2*Mo,1);
		S.Data[0][1]->Data[(i+Mo)*2*Mo+j+Mo] = M2(i+Mo,j)*B(0,j+Mo,1) + M2(i+Mo,j+Mo)*B(0,j+3*Mo,1);
		M1.Data[i*2*Mo+j] = M2(i,j)*B(0,j,0) + M2(i,j+Mo)*B(0,j+2*Mo,0);
		M1.Data[i*2*Mo+j+Mo] = M2(i,j)*B(0,j+Mo,0) + M2(i,j+Mo)*B(0,j+3*Mo,0);
		M1.Data[(i+Mo)*2*Mo+j] = M2(i+Mo,j)*B(0,j,0) + M2(i+Mo,j+Mo)*B(0,j+2*Mo,0);
		M1.Data[(i+Mo)*2*Mo+j+Mo] = M2(i+Mo,j)*B(0,j+Mo,0) + M2(i+Mo,j+Mo)*B(0,j+3*Mo,0);
	}
	*S.Data[0][0] = (*this)(0,0) + (M1*(*this)(1,0));

	memset(M2.Data,0,4*Mo*Mo*sizeof(Complex));
	for (j=0; j<Mo; ++j) {
		for (i=0; i<Mo; ++i) {
			M1.Data[i*2*Mo+j] = - (*this)(1,i,1,j)*B(0,j,0) - (*this)(1,i,1,j+Mo)*B(0,j+2*Mo,0);
			M1.Data[i*2*Mo+j+Mo] = - (*this)(1,i,1,j)*B(0,j+Mo,0) - (*this)(1,i,1,j+Mo)*B(0,j+3*Mo,0);
			M1.Data[(i+Mo)*2*Mo+j] = - (*this)(1,i+Mo,1,j)*B(0,j,0) - (*this)(1,i+Mo,1,j+Mo)*B(0,j+2*Mo,0);
			M1.Data[(i+Mo)*2*Mo+j+Mo] = - (*this)(1,i+Mo,1,j)*B(0,j+Mo,0) - (*this)(1,i+Mo,1,j+Mo)*B(0,j+3*Mo,0);
		}
		M1.Data[j*2*Mo+j] += 1.; M1.Data[(j+Mo)*2*Mo+j+Mo] += 1.;
		M2.Data[j*2*Mo+j] = B(1,j,0); M2.Data[j*2*Mo+j+Mo] = B(1,j+Mo,0);
		M2.Data[(j+Mo)*2*Mo+j] = B(1,j+2*Mo,0); M2.Data[(j+Mo)*2*Mo+j+Mo] = B(1,j+3*Mo,0);
	}
	M2 /= M1; M1 = M2*(*this)(1,1);
	*S.Data[1][0] = M2*(*this)(1,0);
	for (j=0; j<Mo; j++) {
		for (i=0; i<Mo; i++) {
			S.Data[1][1]->Data[i*2*Mo+j] = M1(i,j)*B(0,j,1) + M1(i,j+Mo)*B(0,j+2*Mo,1);
			S.Data[1][1]->Data[i*2*Mo+j+Mo] = M1(i,j)*B(0,j+Mo,1) + M1(i,j+Mo)*B(0,j+3*Mo,1);
			S.Data[1][1]->Data[(i+Mo)*2*Mo+j] = M1(i+Mo,j)*B(0,j,1) + M1(i+Mo,j+Mo)*B(0,j+2*Mo,1);
			S.Data[1][1]->Data[(i+Mo)*2*Mo+j+Mo] = M1(i+Mo,j)*B(0,j+Mo,1) + M1(i+Mo,j+Mo)*B(0,j+3*Mo,1);
		}
		S.Data[1][1]->Data[j*2*Mo+j] += B(1,j,1); S.Data[1][1]->Data[j*2*Mo+j+Mo] += B(1,j+Mo,1);
		S.Data[1][1]->Data[(j+Mo)*2*Mo+j] += B(1,j+2*Mo,1); S.Data[1][1]->Data[(j+Mo)*2*Mo+j+Mo] += B(1,j+3*Mo,1);
	}

	return S;
}

SMatrix SMatrix::ST() {
  SMatrix T(Mown);
  int i, ii;
  for (i=0; i<Mown; i++) {
    for (ii=0; ii<Mown; ii++)
      {T.Data[1][0]->Data[i*Mown+ii] = 0.; T.Data[0][0]->Data[i*Mown+ii] = -Data[0][1]->Data[i*Mown+ii];}
    T.Data[1][0]->Data[i*Mown+i] = 1.;
  }
  *T.Data[1][0] /= (*Data[1][1]);
  *T.Data[1][1] = *T.Data[1][0]*(*Data[1][0]);
  *T.Data[0][1] = *Data[0][0] + *T.Data[0][0]*(*T.Data[1][1]);
  *T.Data[0][0] *= *T.Data[1][0];
  return T;
}

SMatrix SMatrix::STT() {
  SMatrix T(Mown);
  int i, ii;
  for (i=0; i<Mown; i++) {
    for (ii=0; ii<Mown; ii++) T.Data[0][1]->Data[i*Mown+ii] = 0;
    T.Data[0][1]->Data[i*Mown+i] = 1;
  }
  *T.Data[0][1] /= (*Data[0][0]);
  *T.Data[0][0] = *T.Data[0][1]*(*Data[0][0]);
  *T.Data[1][1] = *Data[1][1]*(*T.Data[0][1]);
  for (i=0; i<Mown; i++) for (ii=0; ii<Mown; ii++) T.Data[1][0]->Data[i*Mown+ii] = -T(1,i,1,ii);
  *T.Data[1][0] = *Data[1][0] + *T.Data[1][0]*(*Data[0][0]);
  return T;
}

// implementations of the class SMatrixD

SMatrixD::SMatrixD(unsigned mown) {
	Mown = mown;
	Data[0][0] = new Vector(mown);
	Data[0][1] = new Vector(mown);
	Data[1][0] = new Vector(mown);
	Data[1][1] = new Vector(mown);
}

SMatrixD::SMatrixD(const SMatrixD& SD) {
	Mown = SD.Mown;
	Data[0][0] = new Vector(SD(0,0)); Data[0][1] = new Vector(SD(0,1));
	Data[1][0] = new Vector(SD(1,0)); Data[1][1] = new Vector(SD(1,1));
}

SMatrixD::~SMatrixD() {if (Mown) {delete Data[0][0]; delete Data[0][1]; delete Data[1][0]; delete Data[1][1];}}

SMatrixD& SMatrixD::operator = (const SMatrixD& SD) {
	if (Mown != SD.Mown) {
		Mown = SD.Mown;
		delete Data[0][0]; delete Data[0][1]; delete Data[1][0]; delete Data[1][1];
		Data[0][0] = new Vector(SD(0,0)); Data[0][1] = new Vector(SD(0,1));
		Data[1][0] = new Vector(SD(1,0)); Data[1][1] = new Vector(SD(1,1));
	} 
	else {*(Data[0][0]) = SD(0,0); *(Data[0][1]) = SD(0,1); *(Data[1][0]) = SD(1,0); *(Data[1][1]) = SD(1,1);}
	return *this;
}

SMatrixD SMatrixD::operator * (const SMatrixD& B) {
	if (Mown != B.Mown) return *this;
	SMatrixD C(Mown);
	Complex tc1, tc2;
	for (int i=0; i<Mown; i++) {
		tc1 = 1./(1. - B(0,i,0)*(*this)(1,i,1));
		tc2 = (*this)(0,i,1)*tc1;
		C.Data[0][0]->Data[i] = (*this)(0,i,0) + tc2*B(0,i,0)*(*this)(1,i,0);
		C.Data[0][1]->Data[i] = tc2*B(0,i,1);
		tc2 = B(1,i,0)*tc1;
		C.Data[1][1]->Data[i] = B(1,i,1) + tc2*(*this)(1,i,1)*B(0,i,1);
		C.Data[1][0]->Data[i] = tc2*(*this)(1,i,0);
	}
	return C;
}
	// incidence from outside to amplitudes between layers; (c0+,c1-) ->(b+,b-)
SMatrixD SMatrixD::operator & (const SMatrixD& B) {
	if (Mown != B.Mown) return *this;
	SMatrixD C(Mown);
	Complex tc;
	for (int i=0; i<Mown; i++) {
		tc = 1./(1. - B(0,i,0)*(*this)(1,i,1));
		C.Data[0][1]->Data[i] = ( C.Data[0][0]->Data[i] = (*this)(0,i,1)*tc )*B(0,i,0);
		C.Data[1][0]->Data[i] = ( C.Data[1][1]->Data[i] = B(1,i,0)*tc )*(*this)(1,i,1);
	}
	return C;
}
	// input: source amplitudes between layers (r+,r-); output: effective amplitudes between layers (a+,a-)
SMatrixD SMatrixD::operator % (const SMatrixD& B) {
	if (Mown != B.Mown) return *this;
	SMatrixD C(Mown);
	for (int i=0; i<Mown; i++) {
		C.Data[0][0]->Data[i] = C.Data[0][1]->Data[i]
			= C.Data[1][0]->Data[i] = C.Data[1][1]->Data[i] = 1./(1. - B(0,i,0)*(*this)(1,i,1));
		C.Data[0][1]->Data[i] *= B(0,i,0); C.Data[1][0]->Data[i] *= (*this)(1,i,1);
	}
	return C;
}
	// input: source amplitudes between layers (r+,r-); output: outgoing waves (b-,b+)
SMatrixD SMatrixD::operator ^ (const SMatrixD& B) {
	if (Mown != B.Mown) return *this;
	SMatrixD C(Mown);
	Complex tc;
	for (int i=0; i<Mown; i++) {
		tc = 1./(1. - B(0,i,0)*(*this)(1,i,1));
		C.Data[0][0]->Data[i] = ( C.Data[1][0]->Data[i] = tc*(*this)(1,i,0) )*B(0,i,0);
		C.Data[1][1]->Data[i] = ( C.Data[0][1]->Data[i] = tc*B(0,i,1) )*(*this)(1,i,1);
	}
	return C;
}

SMatrix SMatrixD::operator * (const SMatrix& B) {
	if (Mown != B.Mown) return *this;
	int i, j; SMatrix C(Mown); Matrix AA(Mown), BB(Mown);

	for (i=0; i<Mown; i++) {
		for (j=0; j<Mown; j++) {AA.Data[i*Mown+j] = -(*this)(1,i,1)*B(0,i,0,j); BB.Data[i*Mown+j] = 0.;}
		AA.Data[i*Mown+i] += 1.; BB.Data[i*Mown+i] = (*this)(1,i,0);
	}
	BB %= AA; *C.Data[0][0] = B(0,0)*BB; *C.Data[1][0] = B(1,0)*BB;
	for (i=0; i<Mown; i++) {
		for (j=0; j<Mown; j++) C.Data[0][0]->Data[i*Mown+j] *= (*this)(0,i,1);
		C.Data[0][0]->Data[i*Mown+i] += (*this)(0,i,0);
	}

	for (i=0; i<Mown; i++) {
		for (j=0; j<Mown; j++) {
			AA.Data[i*Mown+j] = -B(0,i,0,j)*(*this)(1,j,1);
			C.Data[0][1]->Data[i*Mown+j] = C.Data[1][1]->Data[i*Mown+j] = 0.;
		}
		AA.Data[i*Mown+i] += 1.;
		C.Data[0][1]->Data[i*Mown+i] = (*this)(0,i,1);
		C.Data[1][1]->Data[i*Mown+i] = (*this)(1,i,1);
	}
	BB = B(0,1)%AA; *C.Data[0][1] *= BB; *C.Data[1][1] = B(1,1) + B(1,0)*C(1,1)*BB;

	return C;
}

SMatrix SMatrixD::operator & (const SMatrix& B) {
	if (Mown != B.Mown) return *this;
	int i, j; Matrix AA(Mown), E(Mown); SMatrix C(Mown);

	for (i=0; i<Mown; i++) {
		for (j=0; j<Mown; j++) {
			AA.Data[i*Mown+j] = -B(0,i,0,j)*(*this)(1,j,1); E.Data[i*Mown+j] = 0.;
			C.Data[0][0]->Data[i*Mown+j] = C.Data[1][0]->Data[i*Mown+j] = 0.;
		}
		AA.Data[i*Mown+i] += 1.; E.Data[i*Mown+i] = 1.;
		C.Data[0][0]->Data[i*Mown+i] = (*this)(0,i,1); C.Data[1][0]->Data[i*Mown+i] = (*this)(1,i,1);
	}
	E /= AA; *C.Data[0][0] *= E; *C.Data[1][0] = B(1,0)*C(1,0)*E;

	for (i=0; i<Mown; i++) {
		for (j=0; j<Mown; j++) {AA.Data[i*Mown+j] = -(*this)(1,i,1)*B(0,i,0,j); E.Data[i*Mown+j] = 0.;}
		AA.Data[i*Mown+i] += 1.; E.Data[i*Mown+i] = 1.;
	}
	E /= AA; *C.Data[0][1] = B(0,0)*E; *C.Data[1][1] = B(1,0)*E;
	for (i=0; i<Mown; i++) for (j=0; j<Mown; j++) C.Data[0][1]->Data[i*Mown+j] *= (*this)(0,i,1);
	
	return C;
}

SMatrix SMatrixD::operator % (const SMatrix& B) {
	if (Mown != B.Mown) return *this;
	int i, j; SMatrix C(Mown); Matrix AA(Mown), E(Mown);

	for (i=0; i<Mown; i++) {
		for (j=0; j<Mown; j++) {
			AA.Data[i*Mown+j] = -B(0,i,0,j)*(*this)(1,j,1); E.Data[i*Mown+j] = 0.;
			C.Data[1][0]->Data[i*Mown+j] = 0.;
		}
		AA.Data[i*Mown+i] += 1.; E.Data[i*Mown+i] = 1.; C.Data[1][0]->Data[i*Mown+i] = (*this)(1,i,1);
	}
	*C.Data[1][0] *= (*C.Data[0][0] = E/AA);

	for (i=0; i<Mown; i++) {
		for (j=0; j<Mown; j++) {AA.Data[i*Mown+j] = -(*this)(1,i,1)*B(0,i,0,j); E.Data[i*Mown+j] = 0.;}
		AA.Data[i*Mown+i] += 1.; E.Data[i*Mown+i] = 1.;
	}
	*C.Data[0][1] = B(0,0)*(*C.Data[1][1] = E/AA);

	return C;
}

SMatrix SMatrixD::operator ^ (const SMatrix& B) {
	if (Mown != B.Mown) return *this;
	int i, j; SMatrix C(Mown); Matrix AA(Mown);

	for (i=0; i<Mown; i++) {
		for (j=0; j<Mown; j++) {AA.Data[i*Mown+j] = -(*this)(1,i,1)*B(0,i,0,j); C.Data[1][0]->Data[i*Mown+j] = 0.;}
		AA.Data[i*Mown+i] += 1.; C.Data[1][0]->Data[i*Mown+i] = (*this)(1,i,0);
	}
	*C.Data[0][0] = B(0,0)*( *C.Data[1][0] %= AA );

	for (i=0; i<Mown; i++) {for (j=0; j<Mown; j++) AA.Data[i*Mown+j] = -B(0,i,0,j)*(*this)(1,j,1); AA.Data[i*Mown+i] += 1.;}
	*C.Data[0][1] = B(0,1)%AA;
	for (i=0; i<Mown; i++) for (j=0; j<Mown; j++) C.Data[1][1]->Data[i*Mown+j] = (*this)(1,i,1)*C(0,i,1,j);

	return C;
}

SDMatrix SMatrixD::operator * (const SDMatrix& B) {
	if (Mown != 2*B.Mown) return *this;
	int i, k, kk, j, Mo = B.Mown; Complex tc;
	Matrix M1(2), M2(2); SDMatrix D(Mo);
	for (i=0; i<Mo; i++) {
		D.Data[0][0]->Data[i] = (*this)(0,i,0);
		D.Data[0][0]->Data[i+3*Mown] = (*this)(0,i+Mo,0);
		for (k=0; k<4; ++k) D.Data[1][1]->Data[i+k*Mo] = B(1,i+k*Mo,1);

		M1.Data[3] = 1. - B(0,i,0)*(*this)(1,i,1); M1.Data[0] = 1. - B(0,i+3*Mo,0)*(*this)(1,i+Mo,1);
		M1.Data[1] = B(0,i+Mo,0)*(*this)(1,i+Mo,1); M1.Data[2] = B(0,i+2*Mo,0)*(*this)(1,i,1);
		tc = 1./(M1(0,0)*M1(1,1) - M1(0,1)*M1(1,0)); for (k=0; k<4; ++k) M1.Data[k] *= tc;

		M2.Data[0] = (*this)(0,i,1)*M1(0,0); M2.Data[1] = (*this)(0,i,1)*M1(0,1);
		M2.Data[2] = (*this)(0,i+Mo,1)*M1(1,0); M2.Data[3] = (*this)(0,i+Mo,1)*M1(1,1);
		D.Data[0][1]->Data[i] = M2(0,0)*B(0,i,1) + M2(0,1)*B(0,i+2*Mo,1);
		D.Data[0][1]->Data[i+Mo] = M2(0,0)*B(0,i+Mo,1) + M2(0,1)*B(0,i+3*Mo,1);
		D.Data[0][1]->Data[i+2*Mo] = M2(1,0)*B(0,i,1) + M2(1,1)*B(0,i+2*Mo,1);
		D.Data[0][1]->Data[i+3*Mo] = M2(1,0)*B(0,i+Mo,1) + M2(1,1)*B(0,i+3*Mo,1);

		D.Data[0][0]->Data[i] += ( M2(0,0)*B(0,i,0) + M2(0,1)*B(0,i+2*Mo,0) )*(*this)(1,i,0);
		D.Data[0][0]->Data[i+Mo] += ( M2(0,0)*B(0,i+Mo,0) + M2(0,1)*B(0,i+3*Mo,0) )*(*this)(1,i+Mo,0);
		D.Data[0][0]->Data[i+2*Mo] += ( M2(1,0)*B(0,i,0) + M2(1,1)*B(0,i+2*Mo,0) )*(*this)(1,i,0);
		D.Data[0][0]->Data[i+3*Mo] += ( M2(1,0)*B(0,i+Mo,0) + M2(1,1)*B(0,i+3*Mo,0) )*(*this)(1,i+Mo,0);

		M1.Data[3] = 1. - (*this)(1,i,1)*B(0,i,0); M1.Data[0] = 1. - (*this)(1,i+Mo,1)*B(0,i+3*Mo,0);
		M1.Data[1] = (*this)(1,i,1)*B(0,i+Mo,0); M1.Data[2] = (*this)(1,i+Mo,1)*B(0,i+2*Mo,0);
		tc = 1./(M1(0,0)*M1(1,1) - M1(0,1)*M1(1,0)); for (k=0; k<4; ++k) M1.Data[k] *= tc;

		M2.Data[0] = B(1,i,0)*M1(0,0) + B(1,i+Mo,0)*M1(1,0);
		M2.Data[1] = B(1,i,0)*M1(0,1) + B(1,i+Mo,0)*M1(1,1);
		M2.Data[2] = B(1,i+2*Mo,0)*M1(0,0) + B(1,i+3*Mo,0)*M1(1,0);
		M2.Data[3] = B(1,i+2*Mo,0)*M1(1,0) + B(1,i+3*Mo,0)*M1(1,1);
		D.Data[1][0]->Data[i] = M2(0,0)*(*this)(1,i,0); D.Data[1][0]->Data[i+Mo] = M2(0,1)*(*this)(1,i+Mo,0);
		D.Data[1][0]->Data[i+2*Mo] = M2(1,0)*(*this)(1,i,0); D.Data[1][0]->Data[i+3*Mo] = M2(1,1)*(*this)(1,i+Mo,0);

		M2.Data[0] *= (*this)(1,i,1); M2.Data[1] *= (*this)(1,i+Mo,1);
		M2.Data[2] *= (*this)(1,i,1); M2.Data[3] *= (*this)(1,i+Mo,1);
		D.Data[1][1]->Data[i] += M2(0,0)*B(0,i,1) + M2(0,1)*B(0,i+2*Mo,1);
		D.Data[1][1]->Data[i+Mo] += M2(0,0)*B(0,i+Mo,1) + M2(0,1)*B(0,i+3*Mo,1);
		D.Data[1][1]->Data[i+2*Mo] += M2(1,0)*B(0,i,1) + M2(1,1)*B(0,i+2*Mo,1);
		D.Data[1][1]->Data[i+3*Mo] += M2(1,0)*B(0,i+Mo,1) + M2(1,1)*B(0,i+3*Mo,1);
	}
	return D;
}

SMatrixD SMatrixD::ST() {
  SMatrixD T(Mown);
  int i;
  Complex tc;
  for (i=0; i<Mown; i++) {
    tc = T.Data[1][1]->Data[i] = 1./(*this)(1,i,1);
    T.Data[1][0]->Data[i] = tc*(*this)(1,i,0);
    T.Data[0][1]->Data[i] = (*this)(0,i,1)*tc;
    T.Data[0][0]->Data[i] = (*this)(0,i,0)-T(0,i,1)*(*this)(1,i,0);
  }
  return T;
}

SMatrixD SMatrixD::STT() {
  SMatrixD T(Mown);
  int i;
  Complex tc;
  for (i=0; i<Mown; i++) {
    tc = T.Data[0][1]->Data[i] = 1./(*this)(0,i,1);
    T.Data[0][0]->Data[i] = tc*(*this)(0,i,0);
    T.Data[1][1]->Data[i] = (*this)(1,i,1)*tc;
    T.Data[1][0]->Data[i] = (*this)(1,i,0)-T(1,i,1)*(*this)(0,i,0);
  }
  return T;
}

// implementations of the class SDMatrix

SDMatrix::SDMatrix(unsigned mown) {
	Mown = mown;
	Data[0][0] = new Vector(4*Mown); Data[0][1] = new Vector(4*Mown);
	Data[1][0] = new Vector(4*Mown); Data[1][1] = new Vector(4*Mown);
}

SDMatrix::SDMatrix(SDMatrix& SD) {
	Mown = SD.Mown;
	Data[0][0] = new Vector(SD(0,0)); Data[0][1] = new Vector(SD(0,1));
	Data[1][0] = new Vector(SD(1,0)); Data[1][1] = new Vector(SD(1,1));
}

SDMatrix::SDMatrix(SMatrixD& SD) {
	Mown = SD.Mown/2;
	Data[0][0] = new Vector(4*Mown); Data[0][1] = new Vector(4*Mown);
	Data[1][0] = new Vector(4*Mown); Data[1][1] = new Vector(4*Mown);
	memset(Data[0][0],0,4*Mown*sizeof(Complex)); memset(Data[0][1],0,4*Mown*sizeof(Complex));
	memset(Data[1][0],0,4*Mown*sizeof(Complex)); memset(Data[1][1],0,4*Mown*sizeof(Complex));
	memcpy(Data[0][0],SD(0,0).Data,Mown*sizeof(Complex));
	memcpy(Data[0][0]+2*Mown,SD(0,0).Data+Mown,Mown*sizeof(Complex));
	memcpy(Data[0][1],SD(0,1).Data,Mown*sizeof(Complex));
	memcpy(Data[0][1]+2*Mown,SD(0,1).Data+Mown,Mown*sizeof(Complex));
	memcpy(Data[1][0],SD(1,0).Data,Mown*sizeof(Complex));
	memcpy(Data[1][0]+2*Mown,SD(1,0).Data+Mown,Mown*sizeof(Complex));
	memcpy(Data[1][1],SD(1,1).Data,Mown*sizeof(Complex));
	memcpy(Data[1][1]+2*Mown,SD(1,1).Data+Mown,Mown*sizeof(Complex));
}

SDMatrix& SDMatrix::operator = (const SDMatrix& SD) {
	if (Mown != SD.Mown) {
		Mown = SD.Mown;
		delete Data[0][0]; delete Data[0][1]; delete Data[1][0]; delete Data[1][1];
		Data[0][0] = new Vector(SD(0,0)); Data[0][1] = new Vector(SD(0,1));
		Data[1][0] = new Vector(SD(1,0)); Data[1][1] = new Vector(SD(1,1));
	} 
	else {*(Data[0][0]) = SD(0,0); *(Data[0][1]) = SD(0,1); *(Data[1][0]) = SD(1,0); *(Data[1][1]) = SD(1,1);}
	return *this;
}

SDMatrix& SDMatrix::operator = (const SMatrixD& SD) {
	if (2*Mown != SD.Mown) {
		Mown = SD.Mown/2;
		delete Data[0][0]; delete Data[0][1]; delete Data[1][0]; delete Data[1][1];
		Data[0][0] = new Vector(4*Mown); Data[0][1] = new Vector(4*Mown);
		Data[1][0] = new Vector(4*Mown); Data[1][1] = new Vector(4*Mown);
	} 
	memset(Data[0][0],0,4*Mown*sizeof(Complex)); memset(Data[0][1],0,4*Mown*sizeof(Complex));
	memset(Data[1][0],0,4*Mown*sizeof(Complex)); memset(Data[1][1],0,4*Mown*sizeof(Complex));
	memcpy(Data[0][0],SD(0,0).Data,Mown*sizeof(Complex)); memcpy(Data[0][0]+2*Mown,SD(0,0).Data+Mown,Mown*sizeof(Complex));
	memcpy(Data[0][1],SD(0,1).Data,Mown*sizeof(Complex)); memcpy(Data[0][1]+2*Mown,SD(0,1).Data+Mown,Mown*sizeof(Complex));
	memcpy(Data[1][0],SD(1,0).Data,Mown*sizeof(Complex)); memcpy(Data[1][0]+2*Mown,SD(1,0).Data+Mown,Mown*sizeof(Complex));
	memcpy(Data[1][1],SD(1,1).Data,Mown*sizeof(Complex)); memcpy(Data[1][1]+2*Mown,SD(1,1).Data+Mown,Mown*sizeof(Complex));
	return *this;
}

SDMatrix SDMatrix::operator * (const SDMatrix& B) {
	if (Mown != B.Mown) return *this;
	int i, k, kk, j; Complex tc;
	Matrix M1(2), M2(2); SDMatrix D(Mown);
	for (i=0; i<Mown; i++) {
		for (k=0; k<4; ++k) {
			D.Data[0][0]->Data[i+k*Mown] = (*this)(0,i+k*Mown,0);
			D.Data[1][1]->Data[i+k*Mown] = B(1,i+k*Mown,1);
		}

		M1.Data[3] = 1. - B(0,i,0)*(*this)(1,i,1) - B(0,i+Mown,0)*(*this)(1,i+2*Mown,1);
		M1.Data[1] = B(0,i,0)*(*this)(1,i+Mown,1) + B(0,i+Mown,0)*(*this)(1,i+3*Mown,1);
		M1.Data[2] = B(0,i+2*Mown,0)*(*this)(1,i,1) + B(0,i+3*Mown,0)*(*this)(1,i+2*Mown,1);
		M1.Data[0] = 1. - B(0,i+2*Mown,0)*(*this)(1,i+Mown,1) - B(0,i+3*Mown,0)*(*this)(1,i+3*Mown,1);
		tc = 1./(M1(0,0)*M1(1,1) - M1(0,1)*M1(1,0)); for (k=0; k<4; ++k) M1.Data[k] *= tc;
		M2.Data[0] = (*this)(0,i,1)*M1(0,0) + (*this)(0,i+Mown,1)*M1(1,0);
		M2.Data[1] = (*this)(0,i,1)*M1(0,1) + (*this)(0,i+Mown,1)*M1(1,1);
		M2.Data[2] = (*this)(0,i+2*Mown,1)*M1(0,0) + (*this)(0,i+3*Mown,1)*M1(1,0);
		M2.Data[3] = (*this)(0,i+2*Mown,1)*M1(0,1) + (*this)(0,i+3*Mown,1)*M1(1,1);
		D.Data[0][1]->Data[i] = M2(0,0)*B(0,i,1) + M2(0,1)*B(0,i+2*Mown,1);
		D.Data[0][1]->Data[i+Mown] = M2(0,0)*B(0,i+Mown,1) + M2(0,1)*B(0,i+3*Mown,1);
		D.Data[0][1]->Data[i+2*Mown] = M2(1,0)*B(0,i,1) + M2(1,1)*B(0,i+2*Mown,1);
		D.Data[0][1]->Data[i+3*Mown] = M2(1,0)*B(0,i+Mown,1) + M2(1,1)*B(0,i+3*Mown,1);

		tc = M2(0,0); M2.Data[0] = tc*B(0,i,0) + M2(0,1)*B(0,i+2*Mown,0);
		M2.Data[1] = tc*B(0,i+Mown,0) + M2(0,1)*B(0,i+3*Mown,0);
		tc = M2(1,0); M2.Data[2] = tc*B(0,i,0) + M2(1,1)*B(0,i+2*Mown,0);
		M2.Data[3] = tc*B(0,i+Mown,0) + M2(1,1)*B(0,i+3*Mown,0);
		D.Data[0][0]->Data[i] += M2(0,0)*(*this)(1,i,0) + M2(0,1)*(*this)(1,i+2*Mown,0);
		D.Data[0][0]->Data[i+Mown] += M2(0,0)*(*this)(1,i+Mown,0) + M2(0,1)*(*this)(1,i+3*Mown,0);
		D.Data[0][0]->Data[i+2*Mown] += M2(1,0)*(*this)(1,i,0) + M2(1,1)*(*this)(1,i+2*Mown,0);
		D.Data[0][0]->Data[i+3*Mown] += M2(1,0)*(*this)(1,i+Mown,0) + M2(1,1)*(*this)(1,i+3*Mown,0);

		M1.Data[3] = 1. - (*this)(1,i,1)*B(0,i,0) - (*this)(1,i+Mown,1)*B(0,i+2*Mown,0);
		M1.Data[0] = 1. - (*this)(1,i+2*Mown,1)*B(0,i+Mown,0) - (*this)(1,i+3*Mown,1)*B(0,i+3*Mown,0);
		M1.Data[1] = (*this)(1,i,1)*B(0,i+Mown,0) + (*this)(1,i+Mown,1)*B(0,i+3*Mown,0);
		M1.Data[2] = (*this)(1,i+2*Mown,1)*B(0,i,0) + (*this)(1,i+3*Mown,1)*B(0,i+2*Mown,0);
		tc = 1./(M1(0,0)*M1(1,1) - M1(0,1)*M1(1,0)); for (k=0; k<4; ++k) M1.Data[k] *= tc;

		M2.Data[0] = B(1,i,0)*M1(0,0) + B(1,i+Mown,0)*M1(1,0);
		M2.Data[1] = B(1,i,0)*M1(0,1) + B(1,i+Mown,0)*M1(1,1);
		M2.Data[2] = B(1,i+2*Mown,0)*M1(0,0) + B(1,i+3*Mown,0)*M1(1,0);
		M2.Data[3] = B(1,i+2*Mown,0)*M1(0,1) + B(1,i+3*Mown,0)*M1(1,1);
		D.Data[1][0]->Data[i] = M2(0,0)*(*this)(1,i,0) + M2(0,1)*(*this)(1,i+2*Mown,0);
		D.Data[1][0]->Data[i+Mown] = M2(0,0)*(*this)(1,i+Mown,0) + M2(0,1)*(*this)(1,i+3*Mown,0);
		D.Data[1][0]->Data[i+2*Mown] = M2(1,0)*(*this)(1,i,0) + M2(1,1)*(*this)(1,i+2*Mown,0);
		D.Data[1][0]->Data[i+3*Mown] = M2(1,0)*(*this)(1,i+Mown,0) + M2(1,1)*(*this)(1,i+3*Mown,0);

		tc = M2(0,0); M2.Data[0] = tc*(*this)(1,i,1) + M2(0,1)*(*this)(1,i+2*Mown,1);
		M2.Data[1] = tc*(*this)(1,i+Mown,1) + M2(0,1)*(*this)(1,i+3*Mown,1);
		tc = M2(1,0); M2.Data[2] = tc*(*this)(1,i,1) + M2(1,1)*(*this)(1,i+2*Mown,1);
		M2.Data[3] = tc*(*this)(1,i+Mown,1) + M2(1,1)*(*this)(1,i+3*Mown,1);
		D.Data[1][1]->Data[i] += M2(0,0)*B(0,i,1) + M2(0,1)*B(0,i+2*Mown,1);
		D.Data[1][1]->Data[i+Mown] += M2(0,0)*B(0,i+Mown,1) + M2(0,1)*B(0,i+3*Mown,1);
		D.Data[1][1]->Data[i+2*Mown] += M2(1,0)*B(0,i,1) + M2(1,1)*B(0,i+2*Mown,1);
		D.Data[1][1]->Data[i+3*Mown] += M2(1,0)*B(0,i+Mown,1) + M2(1,1)*B(0,i+3*Mown,1);
	}
	return D;
}

SDMatrix SDMatrix::operator * (const SMatrixD& B) {
	if ((2*Mown) != B.Mown) return *this;
	int i, k, kk, j, N = B.Mown; Complex tc;
	Matrix M1(2), M2(2); SDMatrix D(Mown);
	for (i=0; i<Mown; i++) {
		for (k=0; k<4; ++k) D.Data[0][0]->Data[i+k*Mown] = (*this)(0,i+k*Mown,0);
		D.Data[1][1]->Data[i] = B(1,i,1); D.Data[1][1]->Data[i+3*Mown] = B(1,i+Mown,1);
		D.Data[1][1]->Data[i+Mown] = D.Data[1][1]->Data[i+2*Mown] = 0.;

		M1.Data[3] = 1. - B(0,i,0)*(*this)(1,i,1); M1.Data[0] = 1. - B(0,i+Mown,0)*(*this)(1,i+3*Mown,1);
		M1.Data[1] = B(0,i,0)*(*this)(1,i+Mown,1); M1.Data[2] = B(0,i+Mown,0)*(*this)(1,i+2*Mown,1);
		tc = 1./(M1(0,0)*M1(1,1) - M1(0,1)*M1(1,0)); for (k=0; k<4; ++k) M1.Data[k] *= tc;

		M2.Data[0] = (*this)(0,i,1)*M1(0,0) + (*this)(0,i+Mown,1)*M1(1,0);
		M2.Data[1] = (*this)(0,i,1)*M1(0,1) + (*this)(0,i+Mown,1)*M1(1,1);
		M2.Data[2] = (*this)(0,i+2*Mown,1)*M1(0,0) + (*this)(0,i+3*Mown,1)*M1(1,0);
		M2.Data[3] = (*this)(0,i+2*Mown,1)*M1(1,0) + (*this)(0,i+3*Mown,1)*M1(1,1);
		D.Data[0][1]->Data[i] = M2(0,0)*B(0,i,1); D.Data[0][1]->Data[i+Mown] = M2(0,1)*B(0,i+Mown,1);
		D.Data[0][1]->Data[i+2*Mown] = M2(1,0)*B(0,i,1); D.Data[0][1]->Data[i+3*Mown] = M2(1,1)*B(0,i+Mown,1);

		M2.Data[0] = M2(0,0)*B(0,i,0); M2.Data[1] = M2(0,1)*B(0,i+Mown,0);
		M2.Data[2] = M2(1,0)*B(0,i,0); M2.Data[3] = M2(1,1)*B(0,i+Mown,0);
		D.Data[0][0]->Data[i] += M2(0,0)*(*this)(1,i,0) + M2(0,1)*(*this)(1,i+2*Mown,0);
		D.Data[0][0]->Data[i+Mown] += M2(0,0)*(*this)(1,i+Mown,0) + M2(0,1)*(*this)(1,i+3*Mown,0);
		D.Data[0][0]->Data[i+2*Mown] += M2(1,0)*(*this)(1,i,0) + M2(1,1)*(*this)(1,i+2*Mown,0);
		D.Data[0][0]->Data[i+3*Mown] += M2(1,0)*(*this)(1,i+Mown,0) + M2(1,1)*(*this)(1,i+3*Mown,0);

		M1.Data[3] = 1. - (*this)(1,i,1)*B(0,i,0); M1.Data[0] = 1. - (*this)(1,i+3*Mown,1)*B(0,i+Mown,0);
		M1.Data[1] = (*this)(1,i+Mown,1)*B(0,i+Mown,0); M1.Data[2] = (*this)(1,i+2*Mown,1)*B(0,i,0);
		tc = 1./(M1(0,0)*M1(1,1) - M1(0,1)*M1(1,0)); for (k=0; k<4; ++k) M1.Data[k] *= tc;

		M2.Data[0] = B(1,i,0)*M1(0,0); M2.Data[1] = B(1,i,0)*M1(0,1);
		M2.Data[2] = B(1,i+Mown,0)*M1(1,0); M2.Data[3] = B(1,i+Mown,0)*M1(1,1);
		D.Data[1][0]->Data[i] = M2(0,0)*(*this)(1,i,0) + M2(0,1)*(*this)(1,i+2*Mown,0);
		D.Data[1][0]->Data[i+Mown] = M2(0,0)*(*this)(1,i+Mown,0) + M2(0,1)*(*this)(1,i+3*Mown,0);
		D.Data[1][0]->Data[i+2*Mown] = M2(1,0)*(*this)(1,i,0) + M2(1,1)*(*this)(1,i+2*Mown,0);
		D.Data[1][0]->Data[i+3*Mown] = M2(1,0)*(*this)(1,i+Mown,0) + M2(1,1)*(*this)(1,i+3*Mown,0);

		tc = M2(0,0); M2.Data[0] = tc*(*this)(1,i,1) + M2(0,1)*(*this)(1,i+2*Mown,1);
		M2.Data[1] = tc*(*this)(1,i+Mown,1) + M2(0,1)*(*this)(1,i+3*Mown,1);
		tc = M2(1,0); M2.Data[2] = tc*(*this)(1,i,1) + M2(1,1)*(*this)(1,i+2*Mown,1);
		M2.Data[3] = tc*(*this)(1,i+Mown,1) + M2(1,1)*(*this)(1,i+3*Mown,1);
		D.Data[1][1]->Data[i] += M2(0,0)*B(0,i,1); D.Data[1][1]->Data[i+Mown] += M2(0,1)*B(0,i+Mown,1);
		D.Data[1][1]->Data[i+2*Mown] += M2(1,0)*B(0,i,1); D.Data[1][1]->Data[i+3*Mown] += M2(1,1)*B(0,i+Mown,1);
	}
	return D;
}

SMatrix SDMatrix::operator * (const SMatrix& B) {
	if (2*Mown != B.Mown) {SMatrix S(B); return S;}
	int i, j; Matrix M1(2*Mown), M2(2*Mown); SMatrix S(2*Mown);

	memset(M2.Data,0,4*Mown*Mown*sizeof(Complex));
	for (j=0; j<Mown; ++j) {
		for (i=0; i<Mown; ++i) {
			M1.Data[i*2*Mown+j] = - B(0,i,0,j)*(*this)(1,j,1) - B(0,i,0,j+Mown)*(*this)(1,j+2*Mown,1);
			M1.Data[i*2*Mown+j+Mown] = - B(0,i,0,j)*(*this)(1,j+Mown,1) - B(0,i,0,j+Mown)*(*this)(1,j+3*Mown,1);
			M1.Data[(i+Mown)*2*Mown+j] = - B(0,i+Mown,0,j)*(*this)(1,j,1) - B(0,i+Mown,0,j+Mown)*(*this)(1,j+2*Mown,1);
			M1.Data[(i+Mown)*2*Mown+j+Mown] = - B(0,i+Mown,0,j)*(*this)(1,j+Mown,1) - B(0,i+Mown,0,j+Mown)*(*this)(1,j+3*Mown,1);
		}
		M1.Data[j*2*Mown+j] += 1.; M1.Data[(j+Mown)*2*Mown+j+Mown] += 1.;
		M2.Data[j*2*Mown+j] = (*this)(0,j,1); M2.Data[j*2*Mown+j+Mown] = (*this)(0,j+Mown,1);
		M2.Data[(j+Mown)*2*Mown+j] = (*this)(0,j+2*Mown,1); M2.Data[(j+Mown)*2*Mown+j+Mown] = (*this)(0,j+3*Mown,1);
	}
	M2 /= M1; M1 = M2*B(0,0);
	*S.Data[0][1] = M2*B(0,1);
	for (j=0; j<Mown; j++) {
		for (i=0; i<Mown; i++) {
			S.Data[0][0]->Data[i*2*Mown+j] = M1(i,j)*(*this)(1,j,0) + M1(i,j+Mown)*(*this)(1,j+2*Mown,0);
			S.Data[0][0]->Data[i*2*Mown+j+Mown] = M1(i,j)*(*this)(1,j+Mown,0) + M1(i,j+Mown)*(*this)(1,j+3*Mown,0);
			S.Data[0][0]->Data[(i+Mown)*2*Mown+j] = M1(i+Mown,j)*(*this)(1,j,0) + M1(i+Mown,j+Mown)*(*this)(1,j+2*Mown,0);
			S.Data[0][0]->Data[(i+Mown)*2*Mown+j+Mown] = M1(i+Mown,j)*(*this)(1,j+Mown,0) + M1(i+Mown,j+Mown)*(*this)(1,j+3*Mown,0);
		}
		S.Data[0][0]->Data[j*2*Mown+j] += (*this)(0,j,0); S.Data[0][0]->Data[j*2*Mown+j+Mown] += (*this)(0,j+Mown,0);
		S.Data[0][0]->Data[(j+Mown)*2*Mown+j] += (*this)(0,j+2*Mown,0); S.Data[0][0]->Data[(j+Mown)*2*Mown+j+Mown] += (*this)(0,j+3*Mown,0);
	}

	for (j=0; j<Mown; j++) {
		for (i=0; i<Mown; i++) {
			M1.Data[i*2*Mown+j] = - (*this)(1,i,1)*B(0,i,0,j) - (*this)(1,i+Mown,1)*B(0,i+Mown,0,j);
			M1.Data[i*2*Mown+j+Mown] = - (*this)(1,i,1)*B(0,i,0,j+Mown) - (*this)(1,i+Mown,1)*B(0,i+Mown,0,j+Mown);
			M1.Data[(i+Mown)*2*Mown+j] = - (*this)(1,i+2*Mown,1)*B(0,i,0,j) - (*this)(1,i+3*Mown,1)*B(0,i+Mown,0,j);
			M1.Data[(i+Mown)*2*Mown+j+Mown] = - (*this)(1,i+2*Mown,1)*B(0,i,0,j+Mown) - (*this)(1,i+3*Mown,1)*B(0,i+Mown,0,j+Mown);
		}
		M1.Data[j*2*Mown+j] += 1.; M1.Data[(j+Mown)*2*Mown+j+Mown] += 1.;
	}
	M2 = B(1,0)/M1;
	for (j=0; j<Mown; j++) for (i=0; i<Mown; i++) {
		S.Data[1][0]->Data[i*2*Mown+j] = M2(i,j)*(*this)(1,j,0) + M2(i,j+Mown)*(*this)(1,j+2*Mown,0);
		S.Data[1][0]->Data[i*2*Mown+j+Mown] = M2(i,j)*(*this)(1,j+Mown,0) + M2(i,j+Mown)*(*this)(1,j+3*Mown,0);
		S.Data[1][0]->Data[(i+Mown)*2*Mown+j] = M2(i+Mown,j)*(*this)(1,j,0) + M2(i+Mown,j+Mown)*(*this)(1,j+2*Mown,0);
		S.Data[1][0]->Data[(i+Mown)*2*Mown+j+Mown] = M2(i+Mown,j)*(*this)(1,j+Mown,0) + M2(i+Mown,j+Mown)*(*this)(1,j+3*Mown,0);
		M1.Data[i*2*Mown+j] = M2(i,j)*(*this)(1,j,1) + M2(i,j+Mown)*(*this)(1,j+2*Mown,1);
		M1.Data[i*2*Mown+j+Mown] = M2(i,j)*(*this)(1,j+Mown,1) + M2(i,j+Mown)*(*this)(1,j+3*Mown,1);
		M1.Data[(i+Mown)*2*Mown+j] = M2(i+Mown,j)*(*this)(1,j,1) + M2(i+Mown,j+Mown)*(*this)(1,j+2*Mown,1);
		M1.Data[(i+Mown)*2*Mown+j+Mown] = M2(i+Mown,j)*(*this)(1,j+Mown,1) + M2(i+Mown,j+Mown)*(*this)(1,j+3*Mown,1);
	}
	*S.Data[1][1] = B(1,1) + (M1*B(0,1));

	return S;
}

// SVector implementations //////////////////////////////////////////////////////////////////////////////

SVector::SVector(unsigned mown) {Mown = mown; Data[0] = new Vector(mown); Data[1] = new Vector(mown);}

SVector::SVector(const SVector& SM) {
	Mown = SM.Mown; Data[0] = new Vector(Mown); *(Data[0]) = SM(0); Data[1] = new Vector(Mown); *(Data[1]) = SM(1);
}

SVector::~SVector() {delete Data[0]; delete Data[1];}

SVector& SVector::operator = (const SVector& SD) {
	if (Mown != SD.Mown) {
		Mown = SD.Mown; delete Data[0]; delete Data[1]; 
		Data[0] = new Vector(SD(0)); Data[1] = new Vector(SD(1));
	}
	else {*(Data[0]) = SD(0); *(Data[1]) = SD(1);}
	return *this;		
}

SVector SVector::operator * (const SMatrix& B) {
	SVector C(Mown);
	if (Mown != B.Mown) return C; // exception ???
	*C.Data[0] = (*Data[0])*B(0,0) + (*Data[1])*B(1,0); *C.Data[1] = (*Data[0])*B(0,1) + (*Data[1])*B(1,1); 
	return C;
}

SVector SVector::operator * (const SMatrixD& B) {
	int i, j; 
	Complex tc;
	SVector C(Mown);
	if (Mown != B.Mown) return C; // exception ???
	for (j=0; j<2; j++) for (i=0; i<Mown; i++)
		C.Data[j]->Data[i] = Data[0]->Data[i]*B(0,i,1-j)+Data[1]->Data[i]*B(1,i,1-j);
	return C;
}

// QMatrix implementations //////////////////////////////////////////////////////////////////////////////

QMatrix& QMatrix::operator = (const QMatrix& Q) {
	Data[0][0] = Q.Data[0][0]; Data[0][1] = Q.Data[0][1];
	Data[1][0] = Q.Data[1][0]; Data[1][1] = Q.Data[1][1];
	return *this;
}

QMatrix QMatrix::operator * (const QMatrix& B) {
	QMatrix M;
	Complex tv;
	tv = 1./(1.-Data[1][1]*B.Data[0][0]);
	M.Data[0][0] = Data[0][0] + Data[0][1]*B.Data[0][0]*Data[1][0]*tv;
	M.Data[0][1] = Data[0][1]*B.Data[0][1]*tv;
	M.Data[1][0] = Data[1][0]*B.Data[1][0]*tv;
	M.Data[1][1] = B.Data[1][1] + Data[1][1]*B.Data[0][1]*B.Data[1][0]*tv;
	return M;
}
	// input: incident amplitudes at structure boundaries (a+,a-); output: amplitudes between structure layers (b+,b-)
QMatrix QMatrix::operator & (const QMatrix& B) {
	QMatrix M;
	Complex tv;
	tv = 1./(1.-Data[1][1]*B.Data[0][0]);
	M.Data[0][1] = B.Data[0][0]*(M.Data[0][0] = Data[0][1]*tv);
	M.Data[1][0] = Data[1][1]*(M.Data[1][1] = B.Data[1][0]*tv);
	return M;
}
	// input: source amplitudes between lyaers (r+,r-); output: effective amplitudes between layers (a+,a-)
QMatrix QMatrix::operator % (const QMatrix& B) {
	QMatrix M;
	M.Data[0][0] = M.Data[1][1] = M.Data[1][0] = M.Data[0][1] = 1./(1. - Data[1][1]*B.Data[0][0]);;
	M.Data[0][1] *= B.Data[0][0]; M.Data[1][0] *= Data[1][1];
  return M;
}
	// input: source amplitudes between lyaers (r+,r-); output: outgoing waves (b-,b+)
QMatrix QMatrix::operator ^ (const QMatrix& B) {
	QMatrix M;
	Complex tv;
	tv = 1./(1.-Data[1][1]*B.Data[0][0]);
	M.Data[0][0] = B.Data[0][0]*( M.Data[1][0] = Data[1][0]*tv );
	M.Data[1][1] = Data[1][1]*( M.Data[0][1] = B.Data[0][1]*tv );
	return M;
}
