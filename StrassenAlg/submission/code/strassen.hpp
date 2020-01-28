#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <ctime>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <assert.h>
#include <random>

using namespace std;
/**
 * Checks if the vector/matrix A is n x m.
 * 
 * @param n number of rows
 * @param m number of columns
 * @return 0 if check successful
 *         1 if check unsuccessful
 */
int checkDimensions(vector<double>& A, uint n, uint m);

/**
 * Print matrix M to console
 */
int printMat(vector<double> M, uint n, uint m);

/**
 * Returns element at position (i,j) of matrix.
 * Returns -10000 if (i,j) out of scope
 */
double getElement(vector<double>& A, uint i, uint j);

/**
 * Places the value alpha (default: 1) into the k-th diagonal.
 * 
 * k = 0 is the main diagonal.
 * k > 0 are the subdiagonals below the main diagonal.
 * k < 0 are the subdiagonals above the main diagonal.
 * 
 * Returns 0 if executed successfully.
 */
int MatEye(vector<double>& A, uint rowSize, double alpha=1, uint k=0);

/**
 * Places uniformly distributed random (double) values from 
 * the intervall [min,max) into every element of the vector/matrix A.
 * Seed has to be provided or otherwise stays default (100).
 * 
 * RNG...MT19937_64
 * 
 * Returns 0 if executed successfully.
 */
int MatRand(vector<double>& A, double min=0.0, double max=10.0, uint seed=100);

/**
 * C = A * B
 * 
 * Standard naive matrix-matrix multiplication.
 * C is (n x m)-matrix (vector<double> of size n*m).
 */
int StdMatMult(uint n, uint m, vector<double>& A, vector<double>& B, vector<double>& C);

/**
 * C = A * B
 * 
 * Standard naive matrix-matrix multiplication for submatrices.
 * Submatrices of main matrices (A,B,C) are specified by the top-left
 * most element's index. Size of the submatrices is given by n.
 * C is (n x n)-(sub)matrix (vector<double> of size n*m).
 */
int StdSubMatMult(uint rowSize, uint n, vector<double>& A, uint iaMin, uint jaMin, vector<double>& B, uint ibMin, uint jbMin, vector<double>& C, uint icMin, uint jcMin);

/**
 * Submatrix addition.
 * 
 * C[icMin:+n][jcMin:+m] = A[iaMin:+n][jaMin:+m] + B[ibMin:+n][jbMin:+m]
 */
int AddMat(uint rowSize, vector<double>& A, uint iaMin, uint jaMin, vector<double>& B, uint ibMin, uint jbMin, vector<double>& C, uint icMin, uint jcMin, uint n, uint m);

/**
 * Submatrix subtraction.
 * 
 * C[icMin:+n][jcMin:+m] = A[iaMin:+n][jaMin:+m] - B[ibMin:+n][jbMin:+m]
 */
int SubMat(uint rowSize, vector<double>& A, uint iaMin, uint jaMin, vector<double>& B, uint ibMin, uint jbMin, vector<double>& C, uint icMin, uint jcMin, uint n, uint m);

/**
 * Set submatrix C to content of submatrix A.
 * 
 * C[icMin:+n][jcMin:+m] = A[iaMin:+n][jaMin:+m]
 */
int SetMat(uint rowSize, vector<double>& A, uint iaMin, uint jaMin, vector<double>& C, uint icMin, uint jcMin, uint n, uint m);

/**
 * Strassen algorithm for matrix-matrix multiplication.
 * !!! Change to signature for public use !!!
 */
int StrassenR(uint rowSize, uint n, vector<double>& A, uint iaMin, uint jaMin,
vector<double>& B, uint ibMin, uint jbMin, vector<double>& C, uint icMin, uint jcMin, vector<double>& W, uint iWMin, uint jWMin, uint min_size);

/**
 * Strassen algorithm for matrix-matrix multiplication.
 * !!! Change to signature for recursion use !!!
 */
int Strassen(uint n, uint m, vector<double>& A, vector<double>& B, vector<double>& C, vector<double>& W, uint min_size);