
#ifndef _MATRIX_H
#define _MATRIX_H

#include <iostream>
#include <complex>

#define _USE_MATH_DEFINES

using namespace std;

typedef complex<double> Complex;

const Complex j_(0.,1.);
/**
#ifdef _M_X64
	#pragma comment(lib,"LNM_Lapack.lib")
	#define DllExport extern "C" __declspec(dllexport) bool __cdecl  //
	DllExport zXeqAmB(Complex *X, Complex *A, Complex *B, int *ni, int *nj, int *nk, bool *brk);
	DllExport zA_1inplace(Complex *a, int *n, bool *brk);
	DllExport zBeqA_1mB(Complex *A, int *n, Complex *B, int *nx, bool *brk);
	DllExport zeigValVecA(Complex*, Complex*, Complex*, int*, bool*);
	DllExport zeigValA(Complex*, Complex*, int*, bool*);
#else
	#pragma comment(lib,"LAPACK.lib")
	#define DllExport __declspec(dllexport) bool __cdecl // extern "C"
	DllExport zXeqAmB(Complex *X, Complex *A, Complex *B, int *ni, int *nj, int *nk, bool *brk);
	DllExport zA_1inplace(Complex *a, int *n, bool *brk);
	DllExport zBeqA_1mB(Complex *A, int *n, Complex *B, int *nx, bool *brk);
	DllExport zeigValVecA(Complex*, Complex*, Complex*, int*, bool*);
	DllExport zeigValA(Complex*, Complex*, int*, bool*);
#endif
/**/
void qsort(Complex*, int, int);

class Vector;
class Matrix;

class Vector {
public:
	unsigned Nrow;
	Complex* Data;

	Vector(unsigned Mown = 1);
	Vector(const Vector&);
	~Vector();

	Vector& operator = (const Vector&);
	Complex& operator() (unsigned i) const {
		if (i<Nrow) return Data[i];
		else {cout<<"vector: out of boundaries "<<i<<" "<<Nrow<<endl; return Data[0];}
	}

	Vector& operator += (const Vector&);
	Vector& operator -= (const Vector&);

	Vector operator * (const Complex);
	Vector& operator *= (const Complex tc) {(*this) = (*this)*tc; return *this;}

	Complex operator * (const Vector&) const;
	Vector operator * (const Matrix& B);
	Vector& operator *= (const Matrix& B) {*this = (*this)*B; return *this;}
	Vector& operator /= (const Matrix&);
	Vector& operator %= (const Matrix&);
	Vector operator + (const Vector& B) {Vector C(*this); C += B; return C;}
	Vector operator - (const Vector& B) {Vector C(*this); C -= B; return C;}
	Vector operator / (const Matrix& B) {Vector C(*this); C /= B; return C;}
	Vector operator % (const Matrix& B) {Vector C(*this); C %= B; return C;}

	double normF(int nn = 0) const;

	double cmp(const Vector&);
	double cmpa(const Vector&);

	Vector calc_LSa(const Matrix&);

	void GMResG(const Matrix&, int, int, double mmax = 1.e-8, bool wr = true);
	void GMResG(void*, Vector (*mul)(const Vector&, void*), int, int, double mmax = 1.e-8, bool wr = true);
	void GMResAG(const Matrix&, int, double mmax = 1.e-8, bool wr = true);
	void GMResAG(void*, Vector (*mul)(const Vector&, void*), int, double mmax = 1.e-8, bool wr = true);
	void GMResM(void*, Vector (*mul)(const Vector&, void*), const Vector&, int, double mmax = 1.e-8, bool wr = true);
	void GMResDEFL(void*, Vector (*mul)(const Vector&, void*), const Vector&, int, double mmax = 1.e-8, bool wr = true);
	void GMResAH(const Matrix&, int, int, double mmax = 1.e-8, bool wr = true);
	void GMResAH(void*, Vector (*mul)(const Vector&, void*), int, int, double mmax = 1.e-8, bool wr = true);
	void CG(const Matrix&, int, double mmax = 1.e-8, bool wr = true);
	void CG(void*, Vector (*mul)(const Vector&, void*), int, double mmax = 1.e-8, bool wr = true);
	void BiCG(const Matrix&, int, double mmax = 1.e-8, bool wr = true);
	void BiCG(void*, Vector (*mul)(const Vector&, void*), int, double mmax = 1.e-8, bool wr = true);
	void BiCGs(const Matrix&, int, double mmax = 1.e-8, bool wr = true);
	void BiCGs(void*, Vector (*mul)(const Vector&, void*), int, double mmax = 1.e-8, bool wr = true);
	void BiCGs(void*, Vector (*mul)(const Vector&, void*), const Vector&, int, double mmax = 1.e-8, bool wr = true);

	Vector roots();
	Complex poly(Complex);
	Vector pcoeff();
	void qsort();
	void qsorta();
	void qsort(Matrix&);
	Vector polynom(const Vector&);

	friend Vector operator * (Complex tc, const Vector &V) {Vector VV(V); VV *= tc; return VV;};

	Matrix Toeplitz(void);
	void inv_Toeplitz(Vector&, Vector&);
	Vector mul_inv_Toeplitz(Vector&, Vector&);

		// adaptive precision operators
	Complex mul_K1(const Vector&, int K) const;
	void mul_KK(const Vector&, Complex *res, int K) const;
	Vector mulR_K1(const Matrix&, int K) const;
	Vector mulL_K1(const Matrix&, int K) const;
	void mulR_KK(const Matrix&, Vector**, int K) const;
	void mulL_KK(const Matrix&, Vector**, int K) const;
};

class Matrix {
public:
	unsigned Nrow, Ncol;
	Complex* Data;

	Matrix(unsigned Mown = 1, unsigned Mext = 0);
	Matrix(const Matrix&);
	~Matrix();

	void print(const int pres) const;

	Matrix& operator = (const Matrix&);
	Complex& operator () (unsigned i, unsigned j) const {
		if (i<Nrow && j<Ncol) return Data[i*Ncol + j];
		else {cout<<"matrix: out of boundaries "<<i<<" "<<j<<" "<<Nrow<<" "<<Ncol<<endl; return Data[0];}
	}

	double cmp(const Matrix&);
	double cmpa(const Matrix&);
	double normF(int n1=-1) const;

	void eye();

	Matrix& operator += (const Matrix&);
	Matrix& operator -= (const Matrix&);
	Vector operator * (const Vector&) const;
	Matrix operator * (const Matrix&) const;
	Matrix& operator *= (const Matrix& B) {*this = (*this)*B; return *this;}
	Matrix& operator /= (const Matrix&);
	Matrix& operator %= (const Matrix&);
	Matrix operator + (const Matrix& B) {Matrix C(*this); C += B; return C;}
	Matrix operator - (const Matrix& B) {Matrix C(*this); C -= B; return C;}
	Matrix operator / (const Matrix& B) {Matrix C(*this); C /= B; return C;}
	Matrix operator % (const Matrix& B)	{Matrix C(*this); C %= B; return C;}

	Matrix operator * (const Complex);
	Matrix& operator *= (const Complex);
	friend Matrix operator * (const Complex c, const Matrix &B) {Matrix A(B); A *= c; return A;};

	Matrix transp() const;
	Matrix mconj() const;

	Complex trace() const;

	Matrix calc_QRHess(Vector&);
	Matrix calc_SchurHess(Matrix&);

	Complex det();
	Vector PVector(Complex);
	Vector LVector(Complex);

	Vector diag(int erm = 0) const;
	Vector diag(Matrix&, int erm = 0) const;

	Matrix inv2(int);
};

class SMatrixD;
class SDMatrix;

class SMatrix {
public:
	unsigned Mown;
	Matrix *Data[2][2];

	SMatrix(unsigned n = 1);
	SMatrix(const SMatrix&);
	SMatrix(const SMatrixD&);
	~SMatrix();

	SMatrix& operator = (const SMatrix&);
	SMatrix& operator = (const SMatrixD&);
	void cut(const SMatrix&);
	Matrix& operator() (unsigned i, unsigned j) const {return *Data[i][j];}
	Complex operator() (unsigned i, unsigned ii, unsigned j, unsigned jj) const {return (*Data[i][j])(ii,jj);}

	SMatrix operator * (const SMatrix&);
	SMatrix operator & (const SMatrix&);
	SMatrix operator % (const SMatrix&);
	SMatrix operator ^ (const SMatrix&);
	SMatrix& operator *= (const SMatrix& B) {return (*this) = (*this)*B;}
	SMatrix& operator &= (const SMatrix& B) {return (*this) = (*this)&B;}
	SMatrix& operator %= (const SMatrix& B) {return (*this) = (*this)%B;}
	SMatrix& operator ^= (const SMatrix& B) {return (*this) = (*this)^B;}
	SMatrix operator * (const SMatrixD&);
	SMatrix operator & (const SMatrixD&);
	SMatrix operator % (const SMatrixD&);
	SMatrix operator ^ (const SMatrixD&);
	SMatrix& operator *= (const SMatrixD& B) {return (*this) = (*this)*B;}
	SMatrix& operator &= (const SMatrixD& B) {return (*this) = (*this)&B;}
	SMatrix& operator %= (const SMatrixD& B) {return (*this) = (*this)%B;}
	SMatrix& operator ^= (const SMatrixD& B) {return (*this) = (*this)^B;}
	SMatrix operator * (const double);
	SMatrix& operator *= (const double);
	SMatrix operator * (const SDMatrix&);
	SMatrix& operator *= (const SDMatrix& B) {return (*this) = (*this)*B;}

	SMatrix ST(void);
	SMatrix STT(void);
};

class SMatrixD {
public:
	unsigned Mown;
	Vector *Data[2][2];

	SMatrixD(unsigned no = 1);
	SMatrixD(const SMatrixD&);
	~SMatrixD();

	Complex operator() (unsigned i, unsigned ii, unsigned j) const
		{return (*Data[i][j])(ii);}
	Vector operator() (unsigned i, unsigned j) const {return *Data[i][j];}

	SMatrixD& operator = (const SMatrixD&);

	SMatrixD operator * (const SMatrixD&);
	SMatrixD operator & (const SMatrixD&);
	SMatrixD operator % (const SMatrixD&);
	SMatrixD operator ^ (const SMatrixD&);
	SMatrix operator * (const SMatrix&);
	SMatrix operator & (const SMatrix&);
	SMatrix operator % (const SMatrix&);
	SMatrix operator ^ (const SMatrix&);
	SMatrixD& operator *= (const SMatrixD& B) {return (*this) = (*this)*B;}
	SMatrixD& operator &= (const SMatrixD& B) {return (*this) = (*this)&B;}
	SMatrixD& operator ^= (const SMatrixD& B) {return (*this) = (*this)^B;}
	SMatrixD& operator %= (const SMatrixD& B) {return (*this) = (*this)%B;}
	SDMatrix operator * (const SDMatrix&);

	SMatrixD ST(void);
	SMatrixD STT(void);
};

class SDMatrix {
public:
	unsigned Mown;
	Vector *Data[2][2];

	SDMatrix(unsigned no = 1);
	SDMatrix(SDMatrix&);
	SDMatrix(SMatrixD&);

	Complex operator() (unsigned i, unsigned ii, unsigned j) const
		{return (*Data[i][j])(ii);}
	Vector operator() (unsigned i, unsigned j) const {return *Data[i][j];}

	SDMatrix& operator = (const SDMatrix&);
	SDMatrix& operator = (const SMatrixD&);

	SDMatrix operator * (const SDMatrix&);
	SDMatrix operator * (const SMatrixD&);
	SDMatrix& operator *= (const SDMatrix& B) {return (*this) = (*this)*B;}
	SDMatrix& operator *= (const SMatrixD& B) {return (*this) = (*this)*B;}
	SMatrix operator * (const SMatrix&);
};

class SVector {
public:
	unsigned Mown;
	Vector *Data[2];

	SVector(unsigned n = 1);
	SVector(const SVector&);
	~SVector();
    
	Complex operator() (unsigned i, unsigned ii) const {
		return Data[i]->Data[ii];
	}
	Vector operator() (unsigned i) const {return *Data[i];}
	SVector& operator = (const SVector&);
	SVector operator * (const SMatrix&);
	SVector operator * (const SMatrixD&);
	SVector& operator *= (const SMatrix& B) {return (*this) = (*this)*B;}
	SVector& operator *= (const SMatrixD& B) {return (*this) = (*this)*B;}

	SVector& operator += (const SVector &V) {*(this->Data[0]) += *V.Data[0]; *(this->Data[1]) += *V.Data[1]; return *this;}
	SVector operator + (const SVector &V) {SVector W(*this); W += V; return W;}
	SVector& operator -= (const SVector &V) {*(this->Data[0]) -= *V.Data[0]; *(this->Data[1]) -= *V.Data[1]; return *this;}
	SVector operator - (const SVector &V) {SVector W(*this); W -= V; return W;}
};

class QMatrix {
public:
	Complex Data[2][2];

	QMatrix() {}
	QMatrix(QMatrix& Q) {Data[0][0] = Q.Data[0][0]; Data[0][1] = Q.Data[0][1]; Data[1][0] = Q.Data[1][0]; Data[1][1] = Q.Data[1][1]; };    

	Complex operator() (unsigned i, unsigned j) const {return Data[i][j];}

	QMatrix& operator = (const QMatrix&);
	QMatrix operator * (const QMatrix&);
	QMatrix operator & (const QMatrix&);
	QMatrix operator % (const QMatrix&);
	QMatrix operator ^ (const QMatrix&);
	QMatrix& operator *= (const QMatrix& B) {return (*this) = (*this)*B;}
	QMatrix& operator &= (const QMatrix& B) {return (*this) = (*this)&B;}
	QMatrix& operator %= (const QMatrix& B) {return (*this) = (*this)%B;}
	QMatrix& operator ^= (const QMatrix& B) {return (*this) = (*this)^B;}
};

#endif