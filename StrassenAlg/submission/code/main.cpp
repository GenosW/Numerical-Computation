#include "strassen.hpp"
#include <ctime>

int main(int argc, char *argv[]) // m is given by command line argument
{
    /* INPUT */
    /* called by: ./strassenTest m [min_size-for-hybrid1 [min_size-for-hybrid2]]
    Needed argument:
    m ... gives size of matrices A,B,C (each n x n matrices, where n=2^m)
    Optional arguments: []
    min_size-for-hybrid1 ... minimal matrix size for hybrid 1. If reached, recursion is resolved by standard matrix multiplication
    min_size-for-hybrid2 ... minimal matrix size for hybrid 2. If reached, recursion is resolved by standard matrix multiplication */
    assert(argc>1);
    string str = argv[1];
    uint m = stoi(str);
    uint strassen_min_size = 0;
    uint strassen_min_size2 = 0;
    if (argc>2)
    {
        str = argv[2];
        strassen_min_size = pow(2,stoi(str));
        if (argc>3)
        {
            str = argv[3];
            strassen_min_size2 = pow(2,stoi(str));
        }   
    }

    /* Testing the Strassen algorithm vs standard matrix multiplication*/
    cout << "Starting test of Strassen algorithm with m= "<< m << endl << endl;
    uint n = pow(2,m);              // matrix size
    vector<double> A(n*n, 0);
    vector<double> B(n*n, 0);
    vector<double> CStr(n*n, 3);
    vector<double> CRef(n*n, 0);
    vector<double> W(n*n, 0);
    
    /* Fill A,B with random numbers */
    uint seed = 67;
    cout << "Seed for generating the matrices: " << seed << endl << endl;
    
    MatRand(A, 0.0, 1.0, seed+1);
    MatRand(B, 0.0, 1.0, seed);

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
    /* Additional output to console if matrices are small */
    if (n<20){
        cout << "-----C-----" << endl;;
        printMat(CRef,n,n);
        cout << "-----CStr-----" << endl;;
        printMat(CStr,n,n);
    }

    /* Declaration of clocks for hybrid Strassen algorithms
        (not needed if only given 1 command line argument) */
    clock_t startS2;
    clock_t endS2;
    clock_t startS3;
    clock_t endS3;
    if (argc>2)
    {   
        vector<double> CStr2(n*n, 3);
        startS2 = clock();
        Strassen(n,n,A,B,CStr2,W,strassen_min_size);
        endS2 = clock();
        /* Additional output to console if matrices are small */
        if (n<20){
            cout << "-----CStr2-----" << endl;;
            printMat(CStr2,n,n);
        }
        
        if (argc>3)
        {  
            vector<double> CStr3(n*n, 3);
            startS3 = clock();
            Strassen(n,n,A,B,CStr3,W,strassen_min_size2);
            endS3 = clock();
            /* Additional output to console if matrices are small */
            if (n<20){
                cout << "-----CStr3-----" << endl;;
                printMat(CStr3,n,n);
            }
        }   
    }
    
    /* OUTPUT */
    cout << "-----CRef-----(standard (naive) Matrix-Matrix-Multiplication)" << endl;
    cout << endl << "Standard took " << difftime(endN,startN)*1000.0/CLOCKS_PER_SEC << " millisec" << endl << endl << endl;
    cout << "-----CStr-----(Strassen algorithm)";
    cout << endl << "Strassen (min_size=" << 1 << ") took " << difftime(endS,startS)*1000.0/CLOCKS_PER_SEC << " millisec" << endl << endl;
    if (argc>2)
    {  
        cout << endl << "Strassen (min_size=" << strassen_min_size << ") took " << difftime(endS2,startS2)*1000.0/CLOCKS_PER_SEC << " millisec" << endl << endl;
        if (argc>3)
        {  
            cout << endl << "Strassen (min_size=" << strassen_min_size2 << ") took " << difftime(endS3,startS3)*1000.0/CLOCKS_PER_SEC << " millisec" << endl << endl;
        }
    }
}
