#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <ctime>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <assert.h>

using namespace std;

// Header: Function declarations
int Strassen(uint n, uint m, vector<double>& A, vector<double>& B, vector<double>& C, vector<double>& W);
int StdMatMult(uint n, uint m, vector<double>& A, vector<double>& B, vector<double>& C);
int checkDimensions(vector<double>& A, uint n, uint m);
// Functions
int checkDimensions(vector<double>& A, uint n, uint m){
    if (A.size() != n*m) return 1;
    return 0;
}

int StdMatMult(uint n, uint m, vector<double>& A, vector<double>& B, vector<double>& C){
    if (checkDimensions(A,n,m) + checkDimensions(B,n,m) + checkDimensions(C,n,m) != 0) return 1;
    for (uint i = 0; i < n; i++)
        {
            for (uint k = 0; k < m; k++)
            {
                C[i*n + k] = 0;
                for (uint j = 0; j < n; j++)
                {
                    C[i*n + k] += A[i*n + j] * B[j*n + k];
                }
            }
        }
	return 0;
}

int AddMat(uint n, uint m, vector<double>& A, vector<double>& B, vector<double>& C){
    if (checkDimensions(A,n,m) + checkDimensions(B,n,m) + checkDimensions(C,n,m) != 0) return 1;
    for (uint i = 0; i < n; i++)
    {
        for (uint j = 0; j < m; j++)
        {
            C[i*n + j] = A[i*n + j] * B[i*n + j];
        }
    }
    return 0;
}

int Strassen(uint n, uint m, vector<double>& A, vector<double>& B, vector<double>& C, vector<double>& W){
    if (checkDimensions(A,n,m) + checkDimensions(B,n,m) + checkDimensions(C,n,m) != 0) return 1;
    if (n ==1) StdMatMult(n,n,A,B,C);

    int h = n/2;
    // Load A11+A22 into W11
    for (uint i = 0; i < h; i++)
    {   //W11 = A11 + A22
        for(uint j = 0; j < h; j++)
        {
           // W[i*h+j] = A[i] + A[(n-h)*n + (n-h+i)]; // Tthink about it again...might have mixed up i,j
        }
        
        //W[i][j] = A[0][i] + A[n-h][n-h+i];
        W[i+2] = A[n+i] + A[(n-h+1)*n + (n-h+i)];
        //W[0][i] = A[1][i] + A[n-h+1][n-h+i];
    }
}

int main(void)
{
    /* Testing the Strassen algorithm vs standard matrix multiplication*/
    cout << 'Starting test of Strassen algorithm' << endl;
    uint n = 3;
    vector<double> A(n*n, 1);
    vector<double> B(n*n, 1);
    vector<double> C(n*n, 0);
    vector<double> CRef(n*n, 0);
    vector<double> W(n*n, 0);

    //StdError = StdMatMult(n,n,A,B,C);
    //StrassenError = Strassen(n,n,A,B,C,W);
    //cout << 'Strassen algorithm terminated due to error: ' << StrassenError << endl;
}