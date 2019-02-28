## Compile project

`mkdir build`
`cd build`
`cmake ..`
`make`

## Clean generated files

`rm -r build`

## Launch program

Usage : `./bench [-i <iter-nb> | -p <precision-exp>]
                          -e <eigvec-nb>
                          [-n <matrix-size> |  -f <matrix-filename>]
                          -m <Krylov-size>`

Options : 

* -i `<iter-nb>` :                  Number of iteration of the method
* -p `<precision-exp>` :    	  Base 10 exposant for minimum precision of the eigenvectors (100 iterations max)
* -e `<eigvec-nb>` :                Number of desired eigeinvectors
* -n `<matrix-size>` :              Generates a n x n random matrix as input
* -f `<matrix-filename>` :  	  Use a `.mtx` matrix as input
* -m `<Krylov-size>` :              Size of the desired Krylov space

You may not use -i and -p options together, but at least one of both is required.

You may not use -n and -f options together, but at least one of both is required.

-e and -m options are mandatory.


## Output

The program generate an `output.csv` file which contains execution time, eigenvectors and eigenvalues for each iteration.
All data is in a single table with the following column :

* iteration : Iteration number when the data was collected
* eigindex : Index of the eigenvector or eigenvalue among the `k` asked for
* value : Value of the data
* type : type of data. `vr` means real part of eigenvalue, `vi` being the imaginary part. `t` stand for time in seconds.


