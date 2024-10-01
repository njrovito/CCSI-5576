void Exit()
{
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  exit(0);
}


//  ==
//  ||
//  ||   C L A S S:    m p i I n f o
//  ||
//  ==

class mpiInfo
{
 public:

  int myPE;
  int numPE;
  int nRealx,  nRealy;
  int nPEx, nPEy;
  int iPE , jPE;
  int iMin, iMax, jMin, jMax ; // The global i-j numbers on this processor
  int nei_n, nei_s, nei_e, nei_w;
  int countx, county;

  double *phiL, *phiR;
  double *phiT, *phiB;
    
  double *phiSend_n,  *phiSend_s;
  double *phiSend_e,  *phiSend_w;
  double *phiRecv_n,  *phiRecv_s, *phiRecv_e,  *phiRecv_w;
  
  MPI_Status  status;
  int         err;
  int         tag;
  MPI_Request request;

  //  -
  //  |
  //  |   GridDecomposition: Set up PE numbering system in figure below and
  //  |                      establish communication arrays.
  //  |
  //  |                      nPEx -- number of PEs in the x-direction
  //  |                      nPEy -- number of PEs in the y-direction
  //  |                      numPE = total number of PEs
  //  |
  //  |                       +-------+-------+ . . .   +-------+
  //  |                       |       |       |         |       |
  //  |                       |       |       |         | numPE |
  //  |                       |       |       |         |       |
  //  |                       +-------+-------+ . . .   +-------+
  //  |                       .       .       .         .       .
  //  |                       .       .       .         .       .
  //  |                       .       .       .         .       .
  //  |                       +-------+-------+ . . .   +-------+
  //  |                       |       |       |         |       |
  //  |                       | nPEx  | nPEx+1|         |       |
  //  |                       |       |       |         |       |
  //  |                       +-------+-------+ . . .   +-------+
  //  |                       |       |       |         |       |
  //  |                       |   0   |   1   |         | nPEx-1|
  //  |                       |       |       |         |       |
  //  |                       +-------+-------+ . . .   +-------+
  //  |
  //  |
  //  -
  

  void GridDecomposition(int _nPEx, int _nPEy, int nCellx , int nCelly)
  {
    nRealx = nCellx;
    nRealy = nCelly;
    // Store and check incoming processor counts
    nPEx = _nPEx;
    nPEy = _nPEy;
    if (nPEx*nPEy != numPE){
        if ( myPE == 0 ) cout << "Fatal Error:  Number of PEs in x-y directions do not add up to numPE" << endl;
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        exit(0);
      }
    
    // Get the i-j location of this processor, given its number.  See figure above:
    jPE = int(myPE/nPEx);
    iPE = myPE - jPE*nPEx;

    // Set neighbor values
    nei_n = nei_s = nei_e = nei_w = -1;  //init to -1 (if no neighbor)
    
    if(iPE > 0){
	      nei_w = myPE - 1;  // West neighbor
      }

    if(iPE < nPEx - 1){
        nei_e = myPE + 1;  // East neighbor
    }

    if(jPE > 0){
        nei_s = myPE - nPEx;  // South neighbor
    }

    if(jPE < nPEy - 1){
        nei_n = myPE + nPEx;  // North neighbor
    }
    
    // For rank debugging
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::string filename = "neighbors_rank_" + std::to_string(rank - 1) + ".dat";
    std::ofstream file;
    file.open(filename, std::ios::out | std::ios::trunc);
    if (file.is_open()) {
            file << "Rank = " << rank << std::endl;
            file << "North neighbor: " << nei_n << std::endl;
            file << "South neighbor: " << nei_s << std::endl;
            file << "East  neighbor: " << nei_e << std::endl;
            file << "West  neighbor: " << nei_w << std::endl;
            file.close();
    }
    
    countx = nRealx + 2;
    county = nRealy + 2;
    
    phiL = new double [ county ];	
    phiR = new double [ county ];	
    phiT = new double [ countx ];	
    phiB = new double [ countx ];
    
    phiSend_n = new double [ countx ];
    phiSend_s = new double [ countx ];
    phiSend_e = new double [ county ];
    phiSend_w = new double [ county ];
    phiRecv_n = new double [ countx ];
    phiRecv_s = new double [ countx ];
    phiRecv_e = new double [ county ];
    phiRecv_w = new double [ county ];

    tag = 0;
  }

  void ExchangeBoundaryInfo(VD &Solution, VD &b){
    sLOOP phiSend_n[s] = 0.;
    sLOOP phiSend_s[s] = 0.;
    tLOOP phiSend_e[t] = 0.;
    tLOOP phiSend_w[t] = 0.;
    
    // ----------------------------------------------
    // Parallel communication on PE Boundaries       ** See fd.h for tLOOP and sLOOP macros **
    // ----------------------------------------------

    // (1.1) Put values into communication arrays
    sLOOP phiSend_n[s] = Solution[pid(s, nRealy - 1)]; 
    sLOOP phiSend_s[s] = Solution[pid(s, 2)]; 
    tLOOP phiSend_w[t] = Solution[pid(2, t)];  
    tLOOP phiSend_e[t] = Solution[pid(nRealx - 1, t)];  
    // (1.2) Send them to neighboring PEs

    // int MPI_Isend(const void* buffer, int count, MPI_Datatype datatype, int recipient, int tag, MPI_Comm communicator, MPI_Request* request); 

    if ( nei_n >= 0 )  err = MPI_Isend(phiSend_n, countx, MPI_DOUBLE, nei_n, 0, MPI_COMM_WORLD, &request);
    if ( nei_s >= 0 )  err = MPI_Isend(phiSend_s, countx, MPI_DOUBLE, nei_s, 0, MPI_COMM_WORLD, &request);
    if ( nei_e >= 0 )  err = MPI_Isend(phiSend_e, county, MPI_DOUBLE, nei_e, 0, MPI_COMM_WORLD, &request);
    if ( nei_w >= 0 )  err = MPI_Isend(phiSend_w, county, MPI_DOUBLE, nei_w, 0, MPI_COMM_WORLD, &request);

    // (1.3) Receive values from neighobring PEs' physical boundaries.
    if ( nei_n >= 0 ) { err = MPI_Irecv(phiRecv_n, countx, MPI_DOUBLE, nei_n, 0, MPI_COMM_WORLD, &request);   MPI_Wait(&request,&status); } 
    if ( nei_s >= 0 ) { err = MPI_Irecv(phiRecv_s, countx, MPI_DOUBLE, nei_s, 0, MPI_COMM_WORLD, &request);   MPI_Wait(&request,&status); }
    if ( nei_e >= 0 ) { err = MPI_Irecv(phiRecv_e, county, MPI_DOUBLE, nei_e, 0, MPI_COMM_WORLD, &request);   MPI_Wait(&request,&status); }
    if ( nei_w >= 0 ) { err = MPI_Irecv(phiRecv_w, county, MPI_DOUBLE, nei_w, 0, MPI_COMM_WORLD, &request);   MPI_Wait(&request,&status); }
    
    // (1.4) Update BCs using the exchanged information
    // Set b vector from adjoining processes
    if (nei_n >= 0) sLOOP b[pid(s, nRealy + 1)] = phiRecv_n[s] ;  // recieving from n
    if (nei_s >= 0) sLOOP b[pid(s, 0)]          = phiRecv_s[s] ;  // recieving from s
    if (nei_w >= 0) tLOOP b[pid(0, t)]          = phiRecv_w[t] ;  // recieving from w
    if (nei_e >= 0) tLOOP b[pid(nRealx + 1, t)] = phiRecv_e[t] ;  // recieving from e
    
  }
  
  int pid(int i,int j) { return (i+1) + (j)*(nRealx+2); }  

};
