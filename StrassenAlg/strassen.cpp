#include "strassen.hpp"

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

double getElement(vector<double>& A, uint i, uint j){
    if (i >= A.size() or j >= A.size()) return -10000;
    return A[i*A.size() + j];
}

int MatEye(vector<double>& A, uint rowSize, double alpha, uint k){
    if (k >= 0){
        for (uint i = 0; i < (rowSize-k); i++){
            A[(i+k)*rowSize + i] = alpha;
        }
    }
    else {
        for (uint i = 0; i < (rowSize-k); i++){
            A[i*rowSize + i+k] = alpha;
        }
    }
    return 0;
}

int MatRand(vector<double>& A, double min, double max, uint seed){
    mt19937_64 RNG(seed);
    uniform_real_distribution<double> dice(0.0, 1.0);
    for (uint i=0; i < A.size(); i++){
        A[i] = min + (max-min)*dice(RNG);
    }
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
    for (uint i = 0; i < n; i++)
        {
            for (uint k = 0; k < n; k++)
            {
                double sum = 0;
                for (uint j = 0; j < n; j++)
                {
                    sum += A[(i+iaMin)*rowSize + (j+jaMin)] * B[(j+ibMin)*rowSize + (k+jbMin)];
                }
                C[(i+icMin)*rowSize + (k+jcMin)] = sum;
            }
        }
	return 0;
}

int StrassenR(uint rowSize, uint n, vector<double>& A, uint iaMin, uint jaMin, vector<double>& B, uint ibMin, uint jbMin, vector<double>& C, uint icMin, uint jcMin, vector<double>& W, uint iWMin, uint jWMin, uint min_size){
    // Matrices of size (min_size)x(min_size) will be processed using the standard naive matrix-matrix multiplication
    if (n <= min_size) {
        if (min_size==1){
            C[icMin*rowSize + jcMin] = A[iaMin*rowSize + jaMin] * B[ibMin*rowSize + jbMin];
            return 3;  
        }
        StdSubMatMult(rowSize,n,A,iaMin,jaMin,B,ibMin,jbMin,C,icMin,jcMin);
        //C[icMin*rowSize + jcMin] = A[iaMin*rowSize + jaMin] * B[ibMin*rowSize + jbMin];
        return 3;
    }

    uint h = n/2;
    // Go over all indices again... need to consider the offset when nesting
    // should end up looking like the first AddMat()-line for A11+A22
    // M1 = (A11 + A22) (B11 + B22)
    AddMat(rowSize,A,iaMin,jaMin,A,iaMin+h,jaMin+h,W,iWMin,jWMin,h,h);          // A11+A22 -> W11
    AddMat(rowSize,B,ibMin,jbMin,B,ibMin+h,jbMin+h,W,iWMin,jWMin+h,h,h);          // B11 + B22 -> W12
    StrassenR(rowSize,h,W,iWMin,jWMin,W,iWMin,jWMin+h,C,icMin,jcMin,W,iWMin+h,jWMin,min_size);    // C11 = M1 = W11 * W12
    SetMat(rowSize,C,icMin,jcMin,C,icMin+h,jcMin+h,h,h);                            // C22 = M1 = C11

    // M2 = (A21 + A22) B11
    AddMat(rowSize,A,iaMin+h,jaMin,A,iaMin+h,jaMin+h,W,iWMin,jWMin,h,h);          // A21+A22 -> W11
    StrassenR(rowSize,h,W,iWMin,jWMin,B,ibMin,jbMin,C,icMin+h,jcMin,W,iWMin+h,jWMin,min_size);    // C21 = M2 = W11 * B11
    SubMat(rowSize,C,icMin+h,jcMin+h,C,icMin+h,jcMin,C,icMin+h,jcMin+h,h,h);          // C22 -= M2 -= C21

    // M3 = A11 (B12 − B22)
    SubMat(rowSize,B,ibMin,jbMin+h,B,ibMin+h,jbMin+h,W,iWMin,jWMin,h,h);          // B12 - B22 -> W11
    StrassenR(rowSize,h,A,iaMin,jaMin,W,iWMin,jWMin,C,icMin,jcMin+h,W,iWMin+h,jWMin,min_size);    // C12 = M3 = A11 * W11
    AddMat(rowSize,C,icMin+h,jcMin+h,C,icMin,jcMin+h,C,icMin+h,jcMin+h,h,h);          // C22 += M3 += C12

    // M4 = A22 (B21 − B11)
    SubMat(rowSize,B,ibMin+h,jbMin,B,ibMin,jbMin,W,iWMin,jWMin,h,h);          // B21 - B11 -> W11
    StrassenR(rowSize,h,A,iaMin+h,jaMin+h,W,iWMin,jWMin,W,iWMin+h,jWMin+h,W,iWMin+h,jWMin,min_size);    // M4 = A22*W11 -> W22
    AddMat(rowSize,C,icMin,jcMin,W,iWMin+h,jWMin+h,C,icMin,jcMin,h,h);          // C11 += M4 += W22
    AddMat(rowSize,C,icMin+h,jcMin,W,iWMin+h,jWMin+h,C,icMin+h,jcMin,h,h);          // C21 += M4 += W22

    // M5 = (A11 + A12) B22
    AddMat(rowSize,A,iaMin,jaMin,A,iaMin,jaMin+h,W,iWMin,jWMin,h,h);          // A11 + A12 -> W11
    StrassenR(rowSize,h,W,iWMin,jWMin,B,ibMin+h,jbMin+h,W,iWMin+h,jWMin+h,W,iWMin+h,jWMin,min_size);    // M5 = W11*B22 -> W22
    SubMat(rowSize,C,icMin,jcMin,W,iWMin+h,jWMin+h,C,icMin,jcMin,h,h);          // C11 -= M5 -= W22
    AddMat(rowSize,C,icMin,jcMin+h,W,iWMin+h,jWMin+h,C,icMin,jcMin+h,h,h);          // C12 += M5 += W22

    // M6 = (A21 − A11) (B11 + B12)
    SubMat(rowSize,A,iaMin+h,jaMin,A,iaMin,jaMin,W,iWMin,jWMin,h,h);          // A21 - A11 -> W11
    AddMat(rowSize,B,ibMin,jbMin,B,ibMin,jbMin+h,W,iWMin,jWMin+h,h,h);          // B11 + B12 -> W12
    StrassenR(rowSize,h,W,iWMin,jWMin,W,iWMin,jWMin+h,W,iWMin+h,jWMin+h,W,iWMin+h,jWMin,min_size);    // M6 = W11*W12 -> W22
    AddMat(rowSize,C,icMin+h,jcMin+h,W,iWMin+h,jWMin+h,C,icMin+h,jcMin+h,h,h);          // C22 += M6 += W22

    // M7 = (A12 − A22) (B21 + B22)
    SubMat(rowSize,A,iaMin,jaMin+h,A,iaMin+h,jaMin+h,W,iWMin,jWMin,h,h);          // A12 - A22 -> W11
    AddMat(rowSize,B,ibMin+h,jbMin,B,ibMin+h,jbMin+h,W,iWMin,jWMin+h,h,h);          // B21 + B22 -> W12
    StrassenR(rowSize,h,W,iWMin,jWMin,W,iWMin,jWMin+h,W,iWMin+h,jWMin+h,W,iWMin+h,jWMin,min_size);    // M6 = W11*W12 -> W22
    AddMat(rowSize,C,icMin,jcMin,W,iWMin+h,jWMin+h,C,icMin,jcMin,h,h);          // C11 += M7 += W22
    
    // End state:
    // C11 = M1 + M4 − M5 + M7
    // C12 = M3 + M5
    // C21 = M2 + M4
    // C22 = M1 − M2 + M3 + M6

    return 0;
}

int Strassen(uint n, uint m, vector<double>& A, vector<double>& B, vector<double>& C, vector<double>& W, uint min_size)
{
    return StrassenR(n,n,A,0,0,B,0,0,C,0,0,W,0,0,min_size);
}
