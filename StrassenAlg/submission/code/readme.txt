Numerical Computation Project
Topic: Strassen Algorithm
Group: 2
Members:    Christian   GOLLMANN,   01435044
            Peter       HOLZNER,    01426733

1) Compilation:
On Linux (and Mac), the project can be compiled by navigating to the folder "/code" and executing the makefile (no special arguments needed).
We did not test on Windows, but it should still work since we did not use any special libraries outside of the std-libraries.

2) Running our test program
On Linux (and Mac), it can be executed with
    ./strassenTest m [min_size-for-hybrid1 [min_size-for-hybrid2]]
e.g. for m=6
    ./strassenTest 6
    ./strassenTest 6 3      # 1 hybrid implementation with min_size=2^3 enabled
    ./strassenTest 6 3 4    # 2nd hybrid implementation with min_size=2^5   
                            # enabled
The command line arguments are simply seperated by spaces and the last two arguments (inside the square brackets []) are both optional.
The optional arguments trigger up to two additional implementations of hybrd algorithms (Strassen algorithm that resolve the recursion once
a minimal submatrix size is reached by using the standard matrix-matrix multiplication algorithm. More information in the documentation.)

3) If questions arise:
We can be reached via our TU-Email addresses:
    e1435044@student.tuwien.ac.at
    peter.anton.holzner@student.tuwien.ac.at

