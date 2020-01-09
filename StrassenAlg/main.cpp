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
int checkDimensions(vector<double>& A, uint n, uint m);
int printMat(vector<double> M, uint n, uint m);
int StdMatMult(uint n, uint m, vector<double>& A, vector<double>& B, vector<double>& C);
int AddMat(uint rowSize, vector<double>& A, uint iaMin, uint jaMin, vector<double>& B, uint ibMin, uint jbMin, vector<double>& C, uint icMin, uint jcMin, uint n, uint m);
int SubMat(uint rowSize, vector<double>& A, uint iaMin, uint jaMin, vector<double>& B, uint ibMin, uint jbMin, vector<double>& C, uint icMin, uint jcMin, uint n, uint m);
int SetMat(uint rowSize, vector<double>& A, uint iaMin, uint jaMin, vector<double>& C, uint icMin, uint jcMin, uint n, uint m);
int Strassen(uint rowSize, uint n, vector<double>& A, uint iaMin, uint jaMin, vector<double>& B, uint ibMin, uint jbMin, vector<double>& C, uint icMin, uint jcMin, vector<double>& W, uint iWMin, uint jWMin);
int StrassenR(uint n, uint m, vector<double>& A, vector<double>& B, vector<double>& C, vector<double>& W);
// Functions
int checkDimensions(vector<double>& A, uint n, uint m){
    if (A.size() != n*m) return 1;
    return 0;
}

int printMat(vector<double> M, uint n, uint m){
    if (checkDimensions(M,n,m)) return 1;
    uint maxIndex = 10;
    cout << "Matrix:";
    for (uint i = 0; i < min(n,maxIndex); i++){
        cout << endl;
        for (uint j = 0; j < min(n,maxIndex); j++)
            cout << M[i*n + j] << '\t';
    }
    cout << endl;
    return 0;
}

int AddMat(uint rowSize, vector<double>& A, uint iaMin, uint jaMin, vector<double>& B, uint ibMin, uint jbMin, vector<double>& C, uint icMin, uint jcMin, uint n, uint m){  
    for (uint i = 0; i < n; i++)
    {
        for (uint j = 0; j < m; j++)
        {
            C[(i+icMin)*rowSize + (j+jcMin)] = A[(i+iaMin)*rowSize + (j+jaMin)] + B[(i+ibMin)*rowSize + (j+jbMin)];
        }
    }
    return 0;
}

int SubMat(uint rowSize, vector<double>& A, uint iaMin, uint jaMin, vector<double>& B, uint ibMin, uint jbMin, vector<double>& C, uint icMin, uint jcMin, uint n, uint m){  
    for (uint i = 0; i < n; i++)
    {
        for (uint j = 0; j < m; j++)
        {
            C[(i+icMin)*rowSize + (j+jcMin)] = A[(i+iaMin)*rowSize + (j+jaMin)] - B[(i+ibMin)*rowSize + (j+jbMin)];
        }
    }
    return 0;
}

int SetMat(uint rowSize, vector<double>& A, uint iaMin, uint jaMin, vector<double>& C, uint icMin, uint jcMin, uint n, uint m){  
    for (uint i = 0; i < n; i++)
    {
        for (uint j = 0; j < m; j++)
        {
            C[(i+icMin)*rowSize + (j+jcMin)] = A[(i+iaMin)*rowSize + (j+jaMin)];
        }
    }
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

int StdSubMatMult(uint rowSize, uint n, vector<double>& A, uint iaMin, uint jaMin, vector<double>& B, uint ibMin, uint jbMin, vector<double>& C, uint icMin, uint jcMin){
    //if (checkDimensions(A,n,m) + checkDimensions(B,n,m) + checkDimensions(C,n,m) != 0) return 1;
    for (uint i = 0; i < n; i++)
        {
            for (uint k = 0; k < n; k++)
            {
                double sum = 0;
                for (uint j = 0; j < n; j++)
                {
                    sum += A.at((i+iaMin)*rowSize + (j+jaMin)) * B.at((j+ibMin)*rowSize + (k+jbMin));
                    // cout << "k= " << j << "-> A[k]= " << A.at((i+iaMin)*rowSize + (j+jaMin)) << endl;
                    // cout << "k= " << j << "-> B[k]= " << B.at((j+ibMin)*rowSize + (k+jbMin)) << endl;
                }
                // cout << "(i,j) = (" << i << "," << k << ")" << endl;
                // cout << "sum = " << sum << endl;
                C[(i+icMin)*rowSize + (k+jcMin)] = sum;
            }
        }
	return 0;
}

int Strassen(uint rowSize, uint n, vector<double>& A, uint iaMin, uint jaMin, vector<double>& B, uint ibMin, uint jbMin, vector<double>& C, uint icMin, uint jcMin, vector<double>& W, uint iWMin, uint jWMin){
    //if (checkDimensions(A,n,n) + checkDimensions(B,n,n) + checkDimensions(C,n,n) != 0) return 1;
    if (n == 1) {
        //StdSubMatMult(rowSize,n,A,iaMin,jaMin,B,ibMin,jbMin,C,icMin,jcMin);
        C[icMin*rowSize + jcMin] = A[iaMin*rowSize + jaMin] * B[ibMin*rowSize + jbMin];
        return 3;
    }

    uint h = n/2;
    // Go over all indices again... need to consider the offset when nesting
    // should end up looking like the first AddMat()-line for A11+A22
    // M1 = (A11 + A22) (B11 + B22)
    AddMat(rowSize,A,iaMin+0,jaMin+0,A,iaMin+h,jaMin+h,W,iWMin+0,jWMin+0,h,h);          // A11+A22 -> W11
    AddMat(rowSize,B,ibMin+0,jbMin+0,B,ibMin+h,jbMin+h,W,iWMin+0,jWMin+h,h,h);          // B11 + B22 -> W12
    Strassen(rowSize,h,W,iWMin+0,jWMin+0,W,iWMin+0,jWMin+h,C,icMin+0,jcMin+0,W,iWMin+h,jWMin+0);    // C11 = M1 = W11 * W12
    SetMat(rowSize,C,icMin+0,jcMin+0,C,icMin+h,jcMin+h,h,h);                            // C22 = M1 = C11
    // printMat(C,rowSize,rowSize);
    // return 0;

    // M2 = (A21 + A22) B11
    AddMat(rowSize,A,iaMin+h,jaMin+0,A,iaMin+h,jaMin+h,W,iWMin+0,jWMin+0,h,h);          // A21+A22 -> W11
    Strassen(rowSize,h,W,iWMin+0,jWMin+0,B,ibMin+0,jbMin+0,C,icMin+h,jcMin+0,W,iWMin+h,jWMin+0);    // C21 = M2 = W11 * B11
    SubMat(rowSize,C,icMin+h,jcMin+h,C,icMin+h,jcMin+0,C,icMin+h,jcMin+h,h,h);          // C22 -= M2 -= C21

    // M3 = A11 (B12 − B22)
    SubMat(rowSize,B,ibMin+0,jbMin+h,B,ibMin+h,jbMin+h,W,iWMin+0,jWMin+0,h,h);          // B12 - B22 -> W11
    Strassen(rowSize,h,A,iaMin+0,jaMin+0,W,iWMin+0,jWMin+0,C,icMin+0,jcMin+h,W,iWMin+h,jWMin+0);    // C12 = M3 = A11 * W11
    AddMat(rowSize,C,icMin+h,jcMin+h,C,icMin+0,jcMin+h,C,icMin+h,jcMin+h,h,h);          // C22 += M3 += C12

    // M4 = A22 (B21 − B11)
    SubMat(rowSize,B,ibMin+h,jbMin+0,B,ibMin+0,jbMin+0,W,iWMin+0,jWMin+0,h,h);          // B21 - B11 -> W11
    Strassen(rowSize,h,A,iaMin+h,jaMin+h,W,iWMin+0,jWMin+0,W,iWMin+h,jWMin+h,W,iWMin+h,jWMin+0);    // M4 = A22*W11 -> W22
    AddMat(rowSize,C,icMin+0,jcMin+0,W,iWMin+h,jWMin+h,C,icMin+0,jcMin+0,h,h);          // C11 += M4 += W22
    AddMat(rowSize,C,icMin+h,jcMin+0,W,iWMin+h,jWMin+h,C,icMin+h,jcMin+0,h,h);          // C21 += M4 += W22

    // M5 = (A11 + A12) B22
    AddMat(rowSize,A,iaMin+0,jaMin+0,A,iaMin+0,jaMin+h,W,iWMin+0,jWMin+0,h,h);          // A11 + A12 -> W11
    Strassen(rowSize,h,W,iWMin+0,jWMin+0,B,ibMin+h,jbMin+h,W,iWMin+h,jWMin+h,W,iWMin+h,jWMin+0);    // M5 = W11*B22 -> W22
    SubMat(rowSize,C,icMin+0,jcMin+0,W,iWMin+h,jWMin+h,C,icMin+0,jcMin+0,h,h);          // C11 -= M5 -= W22
    AddMat(rowSize,C,icMin+0,jcMin+h,W,iWMin+h,jWMin+h,C,icMin+0,jcMin+h,h,h);          // C12 += M5 += W22

    // M6 = (A21 − A11) (B11 + B12)
    SubMat(rowSize,A,iaMin+h,jaMin+0,A,iaMin+0,jaMin+0,W,iWMin+0,jWMin+0,h,h);          // A21 - A11 -> W11
    AddMat(rowSize,B,ibMin+0,jbMin+0,B,ibMin+0,jbMin+h,W,iWMin+0,jWMin+h,h,h);          // B11 + B12 -> W12
    Strassen(rowSize,h,W,iWMin+0,jWMin+0,W,iWMin+0,jWMin+h,W,iWMin+h,jWMin+h,W,iWMin+h,jWMin+0);    // M6 = W11*W12 -> W22
    AddMat(rowSize,C,icMin+h,jcMin+h,W,iWMin+h,jWMin+h,C,icMin+h,jcMin+h,h,h);          // C22 += M6 += W22

    // M7 = (A12 − A22) (B21 + B22)
    SubMat(rowSize,A,iaMin+0,jaMin+h,A,iaMin+h,jaMin+h,W,iWMin+0,jWMin+0,h,h);          // A12 - A22 -> W11
    AddMat(rowSize,B,ibMin+h,jbMin+0,B,ibMin+h,jbMin+h,W,iWMin+0,jWMin+h,h,h);          // B21 + B22 -> W12
    Strassen(rowSize,h,W,iWMin+0,jWMin+0,W,iWMin+0,jWMin+h,W,iWMin+h,jWMin+h,W,iWMin+h,jWMin+0);    // M6 = W11*W12 -> W22
    AddMat(rowSize,C,icMin+0,jcMin+0,W,iWMin+h,jWMin+h,C,icMin+0,jcMin+0,h,h);          // C11 += M7 += W22
    
    // C11 = M1 + M4 − M5 + M7
    // C12 = M3 + M5
    // C21 = M2 + M4
    // C22 = M1 − M2 + M3 + M6

    return 0;
}

int main(void)
{
    /* Testing the Strassen algorithm vs standard matrix multiplication*/
    cout << "Starting test of Strassen algorithm" << endl;
    uint n = 8;
    vector<double> A(n*n, 1);
    vector<double> B(n*n, 0);
    vector<double> C(n*n, 0);
    vector<double> CStr(n*n, 1);
    vector<double> CRef(n*n, 0);
    vector<double> W(n*n, 0);

    uint sMS = n/2;
    for (uint i = 0; i < sMS; i++){
        for (uint j = 0; j < sMS; j++) B[(sMS+j)*n + i] += i+1.1;
    }
    
    AddMat(n,A,0,0,B,sMS,0,C,sMS,sMS,sMS,sMS);
    cout << "-----A-----" << endl;;
    printMat(A,n,n);
    cout << "-----B-----" << endl;;
    printMat(B,n,n);
    cout << "-----C-----" << endl;;
    printMat(C,n,n);
    cout << "-----Multiplication-----" << endl;
    cout << "-----A-----" << endl;;
    printMat(A,n,n);
    cout << "-----B-----" << endl;;
    printMat(B,n,n);
    StdSubMatMult(n,sMS,A,0,0,B,sMS,0,C,sMS,sMS);
    cout << "-----C-----" << endl;;
    printMat(C,n,n);
    cout << endl;
    cout << "_______________TEST START_______________" << endl;
    cout << "Testing the Strassen algorithm vs standard (naive) Matrix-Matrix-Multiplication" << endl;
    StdMatMult(n,n,A,B,CRef);
    Strassen(n,n,A,0,0,B,0,0,CStr,0,0,W,0,0);
    cout << "-----CRef-----(standard (naive) Matrix-Matrix-Multiplication)" << endl;;
    printMat(CRef,n,n);
    cout << "-----CStr-----(Strassen algorithm)" << endl;;
    printMat(CStr,n,n);

    //StdError = StdMatMult(n,n,A,B,C);
    //StrassenError = Strassen(n,n,A,B,C,W);
    //cout << 'Strassen algorithm terminated due to error: ' << StrassenError << endl;
}