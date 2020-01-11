#include "strassen.hpp"

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
    return 0;
}