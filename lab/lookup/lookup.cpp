#include <iostream>
#include <mpi.h>
#include <cmath>
using std :: cout;
using std :: endl;

// Function to interpolate values
double lookupVal(int n, double *x, double *y, double xval) {
    for (int i = 0; i < n - 1; ++i) {
        if (xval >= x[i] && xval <= x[i + 1]) {
            return y[i] + (xval - x[i]) * (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
        }
    }
    return 0.;
}

int main(int argc, char *argv[]) {

    int numPE, myPE;

    MPI_Init(&argc, &argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &numPE);
    MPI_Comm_rank(MPI_COMM_WORLD, &myPE);

    int m = 20;
    int n = 100;

    double x[n], y[n];  
    double xVal[m] = {0};
    double yVal[m] = {0}; 

    for (int i = 0; i < n; ++i) {
        x[i] = i;
        y[i] = i * i;
    }

    for (int i = 0; i < m; ++i) {
        xVal[i] = 2. * i;
    }

    int num_per_process = m / numPE; 
    int start_index = num_per_process * myPE;
    int end_index = start_index + num_per_process;

    for (int i = start_index; i < end_index; ++i) {
        yVal[i] = lookupVal(n, x, y, xVal[i]);
    }

    double test[m] = {0};  // Changed to double

    // Debug: Print out yVal on each process
    /*
    cout << "Process " << myPE << " yVal: ";
    for (int i = start_index; i < end_index; ++i) {
        cout << yVal[i] << " ";
    }
    cout << endl;
    */

    MPI_Gather(yVal + start_index, num_per_process, MPI_DOUBLE, test, num_per_process, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Finalize();

    if (myPE == 0) {
        cout << "======================================\n";
        for (int i = 0; i < m; ++i) {
            cout << test[i] << endl;
        }
    }
    
    return 0;
}
