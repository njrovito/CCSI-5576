// #include "main.h"
#include <mpi.h>

void printArray(string arrayName, int n, int *array, in myPE){
    for(int=0;, i<n; i++){
        cout << "myPE: " << myPE << " " << arrayName << " " << 
    }
}

int main(int argc, char *argv[] ){
    printf("Hello, world!");

    MPI_Init(&argc, &argv); // & = pass the locations in memory, not the value 
    int numPE = 4;  // This needs to match the number of PEs on which the code is running
    int myPE;

    MPI_Comm_size(MPI_COMM_WORLD, &numPE);
    MPI_Comm_rank(MPI_COMM_WORLD, &myPE);

    cout << "myPE:" << myPE << "Hello, World" << endl;

    MPI_Barrier(MPI_COMM_WORLD); // want until all processes are here

    // Broadcast 
    in n = 10;
    in a[n];

    for(int i=0; i<10; i++) a[i] = 0;

    if(myPE == 0){
        for(int i=0; i<10; i++) a[i] = 333;
    }

    MPI_Bcast(a, n, MPI_INT, 0, MPI_COMM_WORLD);  //bcast a from 0 to all 


    MPI_Finalize();

}

/*
| takes output and sends to another command (grep)
| grep "search_str"

'rookiehpc.org
*/
