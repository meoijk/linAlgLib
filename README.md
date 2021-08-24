# linAlgLib

*** Valgrind was used to detect memory management and threading bugs ***

*** ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 0 from 0) ***

*** Linear Algera Library for fun ***
 
*** Linear algebra playground ***


template <class T>
class Matrix

Matrix();

Matrix(int,int);

~Matrix();

T **matrix = NULL;

int rows, cols;

void RandomMatrix();

void NormalizeColumns();

void Ones();

void AllocateMatrix(int, int);

void PrintMatrixrix(void);

void FreeMemory();

void SetElement(int, int, T);

void SetRow(int, T *);

void SetColumn(int, T *);

T GetElement(int, int);

T GetMaxVal();

T Norm();

T GetMinVal();

T *GetMinMax();

T *GetColumn(int);

T *GetRow(int);


// Numerical methods can be implemented in this class
template <class T>

class AlgLib : Matrix<T>

AlgLib();

~AlgLib();

Matrix<T> MultiplyMatrices(Matrix<T>, Matrix<T>);

Matrix<T> SumMatrices(Matrix<T>, Matrix<T>);

Matrix<T> SubtractMatrices(Matrix<T>, Matrix<T>);

Matrix<T> GetDiagonalMatrix(Matrix<T>);

Matrix<T> Transpose(Matrix<T>);

Matrix<T> ScalarMultiplication(T, Matrix<T>);    

Matrix<T> GetLowerMatrix(Matrix<T>);

Matrix<T> GetUpperMatrix(Matrix<T>);

Matrix<T> JacobiMethod(Matrix<T>, Matrix<T>, T);

Matrix<T> SteepestDescent(Matrix<T>, Matrix<T>, T, long int);

void PassingValues(Matrix<T>,Matrix<T>&);


*** How to use it: ***

AlgLib<float> algLib;        

float _A[2][2] = {{2, 1.9},{1.9, 4}};

float _b[2] = {0.2, -4.2};

Matrix<float> A = PassingValues<float,2,2>(_A);    

Matrix<float> b = PassingValues<float,2>(_b);

Matrix<float> x;
    
// D = A * C
Matrix<float> C(A.rows,A.cols);

Matrix<float> D;

// Random mat

C.RandomMatrix(); 

float* minMax = C.GetMinMax();

cout << "C[min,max]: " << "[" << minMax[0] << "," << minMax[1] << "]" << endl;

D = algLib.MultiplyMatrices(A,C);

D.PrintMatrixrix();

// Solving Ax=b

//x = algLib.JacobiMethod(A,b,1e-6);

x = algLib.SteepestDescent(A,b,1e-6,1e4);

x.FreeMemory();

C.FreeMemory();

D.FreeMemory();

b.FreeMemory();

A.FreeMemory();

b.FreeMemory();    

delete[] minMax;
