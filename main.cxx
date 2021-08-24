// Author: Marcelo Elias de Oliveira
// AlgLib playground for fun

#include <iostream>
#include <assert.h>
#include <iomanip>
#include <limits>
#include <math.h>

using namespace std;

template <class T>
class Matrix
{
public:
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

private:
};

// Numerical methods can be implemented in this class
template <class T>
class AlgLib : Matrix<T>
{
public:
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
    

private:
};

template <class T>
Matrix<T>::Matrix() {}

template <class T>
Matrix<T>::Matrix(int r, int c) { AllocateMatrix(r,c); }

template <class T>
Matrix<T>::~Matrix() { }

template <class T>
T *Matrix<T>::GetMinMax()
{
    T *minMax = new T[2];

    minMax[0] = numeric_limits<T>::max();
    minMax[1] = numeric_limits<T>::min();

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            minMax[0] = (matrix[i][j] < minMax[0]) ? matrix[i][j] : minMax[0];
            minMax[1] = (matrix[i][j] > minMax[1]) ? matrix[i][j] : minMax[1];
        }
    }

    return minMax;
}

template <class T>
T Matrix<T>::GetElement(int i, int j)
{
    return matrix[i][j];
}

template <class T>
T Matrix<T>::Norm()
{
    T norm = 0.0;

    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            norm += GetElement(i,j) * GetElement(i,j);

    norm = static_cast<T>( sqrt(norm) );
    return norm;
}

template <class T>
T *Matrix<T>::GetColumn(int id)
{    
    assert(id < cols);

    T *elements = new T[rows];

    for (int i = 0; i < rows; ++i)
        elements[i] = matrix[i][id];    
    return elements;
}

template <class T>
void Matrix<T>::SetElement(int i, int j, T val)
{
    matrix[i][j] = val;
}

template <class T>
T Matrix<T>::GetMaxVal()
{
    T max = numeric_limits<T>::min();
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            max = (matrix[i][j] > max) ? matrix[i][j] : max;
    return max;
}

template <class T>
T Matrix<T>::GetMinVal()
{
    T min = numeric_limits<T>::max();
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            min = (matrix[i][j] < min) ? matrix[i][j] : min;
    return min;
}

template <class T>
void Matrix<T>::AllocateMatrix(int r, int c)
{
    rows = r;
    cols = c;

    matrix = (T **)calloc(rows, sizeof(T *));

    for (int i = 0; i < rows; ++i)
        matrix[i] = (T *)calloc(cols, sizeof(T));

    if (matrix == NULL)
    {
        cerr << "Memory could not be allocated" << endl;
        exit(EXIT_FAILURE);
    }
}

template <class T>
void Matrix<T>::FreeMemory()
{
    if (matrix != NULL)
    {
        for (int i = 0; i < rows; ++i)
            free(matrix[i]);

        free(matrix);
        matrix = NULL;
    }
}

template <class T>
void Matrix<T>::PrintMatrixrix()
{

    if (matrix == NULL)
        return;

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
            cout << std::setprecision(3) << matrix[i][j] << "\t";
        cout << endl;
    }
    cout << endl;
}



template <class T>
void Matrix<T>::RandomMatrix()
{
    if (matrix == NULL)
    {
        cerr << "Memory should be allocated first" << endl;
        exit(EXIT_FAILURE);
    }
    srand(time(NULL));
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            matrix[i][j] = static_cast<T>(rand()) / static_cast<T>(RAND_MAX);
}



template <class T>
void Matrix<T>::NormalizeColumns()
{
    if (matrix == NULL)
    {
        cerr << "Memory should be allocated first" << endl;
        exit(EXIT_FAILURE);
    }

    cout << __func__ << endl;
    for (int j = 0; j < cols; ++j)
    {
        T norm = 0;
        for (int i = 0; i < rows; ++i)        
            norm += matrix[i][j] * matrix[i][j];
        
        norm = static_cast<T>( sqrt(norm)  );
        for (int i = 0; i < rows; ++i) matrix[i][j] /= norm;
    }            
}


template <class T>
void Matrix<T>::Ones()
{
    if (matrix == NULL)
    {
        cerr << "Memory should be allocated first" << endl;
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            matrix[i][j] = 1.0;
}


template <class T>
AlgLib<T>::AlgLib() {}

template <class T>
AlgLib<T>::~AlgLib() {}

template <class T>
Matrix<T> AlgLib<T>::MultiplyMatrices(Matrix<T> A, Matrix<T> B)
{    
    //cout << __func__ << endl;
    assert(A.cols == B.rows);
    Matrix<T> C;    
    C.AllocateMatrix(A.rows, B.cols);

    // C(i,j) = A(i,k) * B(k,j)
    for (int i = 0; i < A.rows; ++i)
    {
        for (int j = 0; j < B.cols; ++j)
        {
            T sum = 0.0;
            for (int k = 0; k < A.cols; ++k)
                sum += A.GetElement(i, k) * B.GetElement(k, j);
            C.SetElement(i, j, sum);
        }
    }
    return C;
}

template <class T>
void AlgLib<T>::PassingValues(Matrix<T> from, Matrix<T>&to)
{
    assert(from.matrix != NULL);

    if (to.matrix == NULL) to.AllocateMatrix(from.rows,from.cols);

    for (int i = 0; i < from.rows; ++i)
        for (int j = 0; j < from.cols; ++j)
            to.SetElement(i,j,from.GetElement(i,j));
}

template <class T>
Matrix<T> AlgLib<T>::SumMatrices(Matrix<T> A, Matrix<T> B)
{
    assert(A.cols == B.cols && A.rows == B.rows);

    Matrix<T> C;    
    C.AllocateMatrix(A.rows, A.cols);

    // C(i,j) = A(i,j) + B(i,j)
    for (int i = 0; i < A.rows; ++i)
        for (int j = 0; j < A.cols; ++j)
            C.SetElement(i, j, A.GetElement(i, j) + B.GetElement(i, j));

    return C;
}

template <class T>
Matrix<T> AlgLib<T>::SubtractMatrices(Matrix<T> A, Matrix<T> B)
{
    Matrix<T> C;
    assert(A.cols == B.cols && A.rows == B.rows);
    C.AllocateMatrix(A.rows, A.cols);

    // C(i,j) = A(i,j) + B(i,j)
    for (int i = 0; i < A.rows; ++i)
        for (int j = 0; j < A.cols; ++j)
            C.SetElement(i, j, A.GetElement(i, j) - B.GetElement(i, j));

    return C;
}


template <class T>
Matrix<T> AlgLib<T>::Transpose(Matrix<T> A)
{
    Matrix<T> Tr;    
    Tr.AllocateMatrix(A.cols, A.rows);

    
    for (int i = 0; i < A.rows; ++i)
        for (int j = 0; j < A.cols; ++j)
            Tr.SetElement(j,i, A.GetElement(i,j));

    return Tr;
}

template <class T>
Matrix<T> AlgLib<T>::ScalarMultiplication(T s, Matrix<T> A)
{       
    assert(A.matrix != NULL);
    Matrix<T> sA;    
    sA.AllocateMatrix(A.rows, A.cols);
    
    for (int i = 0; i < A.rows; ++i)
        for (int j = 0; j < A.cols; ++j)
            sA.SetElement(i, j, s * A.GetElement(i, j));

    return sA;
}



template <class T>
Matrix<T> AlgLib<T>::GetDiagonalMatrix(Matrix<T> A)
{
    Matrix<T> D(A.rows,A.cols);
    
    for (int j = 0; j < A.cols; ++j) 
        D.SetElement(j,j,A.GetElement(j,j));

    return D;
}

template <class T>
Matrix<T> AlgLib<T>::GetLowerMatrix(Matrix<T> A)
{
    Matrix<T> L(A.rows,A.cols);

     for (int k = 0; k < A.cols; ++k)    
        for (int i = 0; i < A.rows; ++i)        
            for (int j = 0; j < A.cols; ++j)            
                if (i > k && j <= k) L.SetElement(i,j,A.GetElement(i,j));                

    return L;
}

template <class T>
Matrix<T> AlgLib<T>::GetUpperMatrix(Matrix<T> A)
{
    Matrix<T> U(A.rows,A.cols);

     for (int k = 0; k < A.cols; ++k)    
        for (int i = 0; i < A.rows; ++i)        
            for (int j = 0; j < A.cols; ++j)            
                if (i < k && j >= k) U.SetElement(i,j,A.GetElement(i,j));                

    return U;
}


template <class T>
Matrix<T> AlgLib<T>::JacobiMethod(Matrix<T> A, Matrix<T> b, T epsilon)
{    
    cout << __func__ << endl;
    // Solving: A.x = b
    assert(A.cols == b.rows && A.rows == A.cols && b.cols == 1);
    Matrix<T> x(A.rows, 1);         
    Matrix<T> y(A.rows, 1);           
    
    T res = numeric_limits<T>::max();    
    while (res > epsilon)            
    {   
        for (int i = 0; i < x.rows; ++i)
        {            
            T sum = 0.0;    
            for (int j = 0; j < A.cols; ++j)
            {
                if (i == j) continue;               
                sum += A.GetElement(i,j) * x.GetElement(j,0);
            }                   
            y.SetElement(i, 0, ( b.GetElement(i,0) - sum ) / A.GetElement(i,i) );                                                            
        }

        Matrix<T> Ay = MultiplyMatrices(A,y);
        Matrix<T> Ay_b = SubtractMatrices( Ay, b );                    
        T _res = Ay_b.Norm();

        res = (_res < res) ? _res : res;            
        //cout << "r:  " << setprecision(6) << _res << endl;        
        PassingValues(y,x);        
        
        Ay.FreeMemory();
        Ay_b.FreeMemory();
    }
    cout << "[A]:" << endl;
    A.PrintMatrixrix();
    cout << "[b]:" << endl;
    b.PrintMatrixrix();
    cout << "[x]:" << endl;
    x.PrintMatrixrix();   

    y.FreeMemory();

    return x;
}

template <class T>
Matrix<T> AlgLib<T>::SteepestDescent(Matrix<T> A, Matrix<T> b, T epsilon, long int maxIter)
{    
    cout << __func__ << endl;

    // Solving: A.x = b
    assert(A.cols == b.rows && b.cols == 1);

    Matrix<T> x(A.rows,1);         
    //x.RandomMatrix();     

    Matrix<T> bt = Transpose(b);        

    long int iter = 0;  
    T supremum = numeric_limits<T>::max();          
    while (supremum > epsilon && iter < maxIter)                
    {   
        Matrix<T> A_x = MultiplyMatrices(A,x);
        
        Matrix<T> g = SumMatrices(A_x,b); 
        Matrix<T> gt =  Transpose(g);        
        Matrix<T> gt_g =  MultiplyMatrices(gt,g);        
        Matrix<T> xt = Transpose(x);
        Matrix<T> A_g = MultiplyMatrices(A,g);
        Matrix<T> gt_A_g = MultiplyMatrices(gt,A_g);

        // alpha: gt.g / gt.A.gn     
        T alpha = -gt_g.GetElement(0,0) / gt_A_g.GetElement(0,0);
        
        // x_{n+1} = x_{n} - alpha . g_{n}
        Matrix<T> alpha_x = ScalarMultiplication(alpha,g);        
        Matrix<T> _x = SumMatrices(x, alpha_x);        
        //gt.PrintMatrixrix();
        PassingValues(_x,x);

        // Stop criterium:  supremum of the gradient vector
        T* minMax = gt.GetMinMax();        
        //gt.PrintMatrixrix();
        supremum = (minMax[1] < supremum) ? minMax[1] : supremum;

        //cout << iter << "\t" << supremum << endl;
        A_x.FreeMemory();
        A_g.FreeMemory();
        gt_A_g.FreeMemory();
        g.FreeMemory();
        gt.FreeMemory();
        gt_g.FreeMemory();
        xt.FreeMemory();
        alpha_x.FreeMemory();
        _x.FreeMemory();
        delete[] minMax;
        iter++;        
    }
    
    cout << "[A]:" << endl;
    A.PrintMatrixrix();
    cout << "[b]:" << endl;
    b.PrintMatrixrix();
    cout << "[x]:" << endl;
    x.PrintMatrixrix();   
    
    bt.FreeMemory();       

    return x;
}

template<typename T, int rows, int cols>
Matrix<T> PassingValues(T array[rows][cols])
{
    Matrix<T> A(rows,cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            A.SetElement(i,j,array[i][j]);
    return A;    
}

template<typename T, int rows>
Matrix<T> PassingValues(T array[rows])
{
    Matrix<T> A(rows,1);
    for (int i = 0; i < rows; ++i)
        A.SetElement(i,0,array[i]);
    return A;    
}

int main(int argc, char *argv[])
{   
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
    
    return EXIT_SUCCESS;     
}

