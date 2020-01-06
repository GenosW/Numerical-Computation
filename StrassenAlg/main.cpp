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

int Strassen(uint n, uint m, vector<double>& A, vector<double>& B, vector<double>& C, vector<double>& W){
    if (checkDimensions(A,n,n) + checkDimensions(B,n,n) + checkDimensions(C,n,n) != 0) return 1;
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