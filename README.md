# Maintainable Matrix Commitment Scheme with	Constant-Size Public Parameters and Incremental	Aggregation

### Dependency

Our code was built on the library of flint 2.8.0, which support large integer operation based on GMP 6.2.1, and use relic 0.5.0 to implement the bilinear pairing. Our code successfully excecuted on a computer with Ubuntu 20.04.5 LTS and the compiler is g++ with version 9.4.0.

### How to run the code

First compile the test file using 
`g++ -o test main.cpp proof.cpp Vc.cpp util.cpp -lgmp -lflint`
then run the executable file.
