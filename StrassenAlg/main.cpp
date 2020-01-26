#include "strassen.hpp"
#include <ctime>

int main(int argc, char *argv[]) // m is given by command line argument
{
    /* INPUT */
    assert(argc==3+1);
    string str = argv[1];
    uint m = stoi(str);
    str = argv[2];
    uint strassen_min_size = pow(2,stoi(str));
    str = argv[3];
    uint strassen_min_size2 = pow(2,stoi(str));

    /* Testing the Strassen algorithm vs standard matrix multiplication*/
    cout << "Starting test of Strassen algorithm with m= "<< m << endl << endl;
    uint n = pow(2,m);              // matrix size
    vector<double> A(n*n, 0);
    vector<double> B(n*n, 0);
    vector<double> C(n*n, 0);
    vector<double> CStr(n*n, 3);
    vector<double> CStr2(n*n, 3);
    vector<double> CRef(n*n, 0);
    vector<double> W(n*n, 0);
    
    //uint sMS = n/2;                 // submatrix size
    /*
    for (uint i = 0; i < sMS; i++){
        for (uint j = 0; j < sMS; j++) B[(sMS+j)*n + i] += i+1.1;
    }*/
    uint seed = 67;
    cout << "Seed for generating the matrices: " << seed << endl << endl;
    MatRand(B, 0.0, 1.0, seed);
    MatRand(A, 0.0, 1.0, seed+1);

    // cout << "-----A-----" << endl;;
    // printMat(A,n,n);
    // cout << "-----B-----" << endl;;
    // printMat(B,n,n);
    // cout << "-----C-----" << endl;;
    // printMat(C,n,n);

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
    /* MAIN EXECUTION CODE */
    clock_t startN = clock();
    StdMatMult(n,n,A,B,CRef);
    clock_t endN = clock();
    clock_t startS = clock();
    Strassen(n,n,A,B,CStr,W,1);
    clock_t endS = clock();
    clock_t startS2 = clock();
    Strassen(n,n,A,B,CStr2,W,strassen_min_size);
    clock_t endS2 = clock();
    clock_t startS3 = clock();
    Strassen(n,n,A,B,CStr2,W,strassen_min_size2);
    clock_t endS3 = clock();
    cout << "-----CRef-----(standard (naive) Matrix-Matrix-Multiplication)" << endl;
    /**/
    //printMat(CRef,n,n);
    cout << endl << "Standard took " << difftime(endN,startN)*1000.0/CLOCKS_PER_SEC << " millisec" << endl << endl << endl;
    cout << "-----CStr-----(Strassen algorithm)";
    //printMat(CStr,n,n);
    cout << endl << "Strassen (min_size=" << 1 << ") took " << difftime(endS,startS)*1000.0/CLOCKS_PER_SEC << " millisec" << endl << endl;
    cout << endl << "Strassen (min_size=" << strassen_min_size << ") took " << difftime(endS2,startS2)*1000.0/CLOCKS_PER_SEC << " millisec" << endl << endl;
    cout << endl << "Strassen (min_size=" << strassen_min_size2 << ") took " << difftime(endS3,startS3)*1000.0/CLOCKS_PER_SEC << " millisec" << endl << endl;
    return 0;
}
