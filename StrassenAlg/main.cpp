#include "strassen.hpp"
#include <ctime>

int main(void)//int argc, char *argv[])
{
    /* Testing the Strassen algorithm vs standard matrix multiplication*/
    cout << "Starting test of Strassen algorithm" << endl << endl;
    uint n = 4096;
    vector<double> A(n*n, 0);
    vector<double> B(n*n, 0);
    vector<double> C(n*n, 0);
    vector<double> CStr(n*n, 1);
    vector<double> CRef(n*n, 0);
    vector<double> W(n*n, 0);
    
    uint sMS = n/2;
/*
    for (uint i = 0; i < sMS; i++){
        for (uint j = 0; j < sMS; j++) B[(sMS+j)*n + i] += i+1.1;
    }*/
    uint seed = 67;
    cout << "Seed for generating the matrices: " << seed << endl << endl;
    MatRand(B, 1.0, 9.9, seed);
    MatRand(A, 1.0, 9.9, seed+1);

    cout << "-----A-----" << endl;;
    printMat(A,n,n);
    cout << "-----B-----" << endl;;
    printMat(B,n,n);
    cout << "-----C-----" << endl;;
    printMat(C,n,n);

    // AddMat(n,A,0,0,B,sMS,0,C,sMS,sMS,sMS,sMS);
    // cout << "-----A-----" << endl;;
    // printMat(A,n,n);
    // cout << "-----B-----" << endl;;
    // printMat(B,n,n);
    // cout << "-----C-----" << endl;;
    // printMat(C,n,n);
    // cout << "-----Multiplication-----" << endl;
    // cout << "-----A-----" << endl;;
    // printMat(A,n,n);
    // cout << "-----B-----" << endl;;
    // printMat(B,n,n);
    // StdSubMatMult(n,sMS,A,0,0,B,sMS,0,C,sMS,sMS);
    // cout << "-----C-----" << endl;;
    // printMat(C,n,n);
    cout << endl;
    cout << "_______________TEST START_______________" << endl;
    cout << "Testing the Strassen algorithm vs standard (naive) Matrix-Matrix-Multiplication" << endl << endl;
    clock_t startN = clock();
    StdMatMult(n,n,A,B,CRef);
    clock_t endN = clock();
    uint strassen_min_size = pow(2,4);
    clock_t startS = clock();
    Strassen(n,n,A,0,0,B,0,0,CStr,0,0,W,0,0,strassen_min_size);
    clock_t endS = clock();
    cout << "-----CRef-----(standard (naive) Matrix-Matrix-Multiplication)" << endl;
    printMat(CRef,n,n);
    cout << endl << "Standard took " << difftime(endN,startN)*1000.0/CLOCKS_PER_SEC << " millisec" << endl << endl << endl;
    cout << "-----CStr-----(Strassen algorithm)";
    printMat(CStr,n,n);
    cout << endl << "Strassen (min_size=" << strassen_min_size << ") took " << difftime(endS,startS)*1000.0/CLOCKS_PER_SEC << " millisec" << endl << endl;

    //StdError = StdMatMult(n,n,A,B,C);
    //StrassenError = Strassen(n,n,A,B,C,W);
    //cout << 'Strassen algorithm terminated due to error: ' << StrassenError << endl;
    return 0;
}
