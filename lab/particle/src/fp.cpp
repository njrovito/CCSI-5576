//  ================================================================================
//  || ||
//  ||              fp ||
//  ||              ------------------------------------------- ||
//  ||              F R E E   P A R T I C L E ||
//  || ||
//  ||              D E M O N S T R A T I O N   C O D E ||
//  ||              ------------------------------------------- ||
//  || ||
//  ||       Developed by: Scott R. Runnels, Ph.D. ||
//  ||                     University of Colorado Boulder ||
//  || ||
//  ||                For: CU Boulder CSCI 4576/5576 and associated labs ||
//  || ||
//  ||           Copyright 2020 Scott Runnels ||
//  || ||
//  ||                     Not for distribution or use outside of the ||
//  ||                     this course. ||
//  || ||
//  ================================================================================

#include "fp.h"

#include "mpi.h"
#include "mpiInfo.h"
// #include "particles.h"

// ==
// ||
// ||      C L A S S:   M E S H
// ||
// ||
// ||
// ||  (21)----(22)----(23)----(24)-----(25)  <--- nRealx+1 , nRealy+1
// ||   |       |       |       |        |
// ||   |       |       |       |        |
// ||   |       |       |       |        |
// ||  (16)----[17]----[18]----[19]-----(20)
// ||   |    1,3|    2,3|    3,3|        |
// ||   |       |       |       |        |
// ||   |       |       |       |        |
// ||  (11)----[12]----[13]----[14]-----(15)
// ||   |    1,2|    2,2|    3,2|        |
// ||   |       |       |       |        |
// ||   |       |       |       |        |
// ||  (6)-----[7]-----[8]-----[9]------(10)
// ||   |    1,1|    2,1|    3,1|        |
// ||   |       |       |       |        |
// ||   |       |       |       |        |
// ||  (1)-----(2)-----(3)-----(4)------(5)
// || 0,0                              4,0
// ||
// ==

class Mesh {
 public:
  double x0, x1, y0, y1;
  VD x, y;
  int nRealx, nRealy, nField, nReal;
  double lengthx, lengthy, dx, dy;
  int myPE;
  VD Qval;

  //  ==
  //  ||
  //  ||  Constructor:  Initialize values
  //  ||
  //  ==

  Mesh(double _x0, double _x1, double _y0, double _y1, int ncell_x, int ncell_y,
       mpiInfo &myMPI) {
    // Copy incoming values

    x0 = _x0;
    x1 = _x1;
    y0 = _y0;
    y1 = _y1;
    myPE = myMPI.myPE;

    // Compute number of real (physical) nodes, and their spacing

    nRealx = ncell_x + 1;
    nRealy = ncell_y + 1;
    nReal = nRealx * nRealy;
    dx = (x1 - x0) / ncell_x;
    dy = (y1 - y0) / ncell_y;

    // Compute the size of the field variables, which also lie on ghost nodes

    nField = (nRealx + 2) * (nRealy + 2);

    // Allocate memory -- Note that the node numbers and field variable numbers
    // must match.  So even though we only will be caring about real nodes,
    // their node numbers are naturally ordered, so must be of size nField.

    x.resize(nField + 1);
    y.resize(nField + 1);
    Qval.resize(nField + 1);

    for (int i = 1; i <= nRealx; ++i)
      for (int j = 1; j <= nRealy; ++j) {
        int p = pid(i, j);
        x[p] = x0 + (i - 1) * dx;
        y[p] = y0 + (j - 1) * dy;
      }
  }

  //  ==
  //  ||
  //  ||  ParticlesOnMesh
  //  ||
  //  ||  Quantify the particles' interaction with the mesh
  //  ||
  //  ||  References:  [1] https://www.particleincell.com/2010/es-pic-method/
  //  ||
  //  ==

  void ParticlesOnMesh(particles &PTCL, mpiInfo &myMPI) {
    double hx, hy;
    double w[5];
    int p[5];
    int iL, iR, jB, jT;
    int iPEnew, jPEnew;  // These store the i-j indicies of the processor
                         // receiving a particle,
                         //    if that particle is leaving the mesh.
    VI ptcl_send_list,
        ptcl_send_PE;  // These collect information about particles that
                       //    have left this processor and are heading onto
                       //    another processor.

    // -
    // |
    // | Determine which particles are still on this mesh and which have left
    // |
    // -

    for (int k = 1; k <= PTCL.n; ++k) {
      // First, check to be sure the particle is still in the mesh.  If it is
      // not, set its "active" flag to zero, and note the processor to which it
      // is going for MPI exchange.

      if (PTCL.active[k] == 1) {
        iPEnew = myMPI.iPE;
        jPEnew = myMPI.jPE;

        if (PTCL.x[k] < x0) {  // leaving the left boundary
          PTCL.active[k] = -1;
          iPEnew = myMPI.iPE - 1;
        }
        if (PTCL.x[k] > x1) {  // leaving the right boundary
          PTCL.active[k] = -1;
          iPEnew = myMPI.iPE + 1;
        }
        if (PTCL.y[k] < y0) {  // leaving the bottom boundary
          PTCL.active[k] = -1;
          jPEnew = myMPI.jPE - 1;
        }
        if (PTCL.y[k] > y1) {  // leaving the top boundary
          PTCL.active[k] = -1;
          jPEnew = myMPI.jPE + 1;
        }
      }

      // The particle is not in the mesh.  Collect this particle into a holding
      // array that will be sent to the neighboring processor.

      if (PTCL.active[k] == -1) {
        if (iPEnew >= 0 && iPEnew < myMPI.nPEx)  // if our new PE is in the
                                                 // boundaries of our total mesh
          if (jPEnew >= 0 && jPEnew < myMPI.nPEy) {
            ptcl_send_list.push_back(k);
            ptcl_send_PE.push_back(
                iPEnew + jPEnew * myMPI.nPEx);  // calculate which PE it needs
                                                // to go to based on iPE and jPE
          }

        PTCL.active[k] = 0;  // Remove it from the list of active particles
      }
    }

    // -
    // |
    // | Give and receive particles to/with other processors
    // |
    // -

    myMPI.ParticleExchange(ptcl_send_list, ptcl_send_PE, PTCL);

    // -
    // |
    // | Accumulate particles to the nodes (to be completed in the esPIC code,
    // next week)
    // |
    // -

    for (int k = 1; k <= nField; ++k) Qval[k] = 0.;

    // -
    // |
    // | Compute forces on particles
    // |
    // -

    for (int k = 1; k <= PTCL.n; ++k) PTCL.xf[k] = PTCL.yf[k] = 0.;

    for (int k = 1; k <= PTCL.n; ++k) {
      if (PTCL.active[k] == 1) {
        PTCL.xf[k] = 0.;
        PTCL.yf[k] = -.4;
      }
    }
  }

#include "mesh_plotter.h"

  int pid(int i, int j) { return (i + 1) + (j) * (nRealx + 2); }
};

//  ==
//  ||
//  ||
//  ||  Main Program
//  ||
//  ||
//  ==
//

int main(int argc, char *argv[]) {
  mpiInfo myMPI;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &myMPI.numPE);
  MPI_Comm_rank(MPI_COMM_WORLD, &myMPI.myPE);

  int nPEx, nPEy, nCellx, nCelly;
  double tEnd, dt;
  double flux, npHat, vx_bdy;

  // -
  // |
  // | Banner and Input
  // |
  // -

  if (myMPI.myPE == 0) {
    cout << "\n";
    cout << "---------------------------------------------\n";
    cout << "\n";
    cout << " F R E E   P A R T I C L E                   \n";
    cout << " D E M O   C O D E                           \n";
    cout << "\n";
    cout << " Running on " << myMPI.numPE << " processors \n";
    cout << "\n";
    cout << "---------------------------------------------\n";
    cout << "\n";
  }

  for (int count = 0; count < argc; ++count) {
    if (!strcmp(argv[count], "-nPEx")) nPEx = atoi(argv[count + 1]);
    if (!strcmp(argv[count], "-nPEy")) nPEy = atoi(argv[count + 1]);
    if (!strcmp(argv[count], "-nCellx")) nCellx = atoi(argv[count + 1]);
    if (!strcmp(argv[count], "-nCelly")) nCelly = atoi(argv[count + 1]);

    if (!strcmp(argv[count], "-flux")) flux = atof(argv[count + 1]);
    if (!strcmp(argv[count], "-vx_bdy")) vx_bdy = atof(argv[count + 1]);
    if (!strcmp(argv[count], "-npHat")) npHat = atoi(argv[count + 1]);

    if (!strcmp(argv[count], "-tEnd")) tEnd = atof(argv[count + 1]);
    if (!strcmp(argv[count], "-dt")) dt = atof(argv[count + 1]);
  }

  if (myMPI.myPE == 0) {
    cout << endl;
    cout << "Input Summary: " << endl;
    cout << "--------------------------- " << endl;
    cout << "No. PE   in  x-direction: " << nPEx << endl;
    cout << "             y-direction: " << nPEy << endl;
    cout << "No. Cells in x-direction: " << nCellx << endl;
    cout << "             y-direction: " << nCelly << endl;
    cout << "Flux density            : " << flux << endl;
    cout << "End Time                : " << tEnd << endl;
    cout << "Time Step               : " << dt << endl;
    cout << endl;
  }

  // -
  // |
  // | MPI / Processor ID
  // |
  // -

  myMPI.GridDecomposition(nPEx, nPEy, nCellx, nCelly);

  // -
  // |
  // | Parallel Grid Generation
  // |
  // -

  double totalLength = 1.;
  double eachPElength_x = totalLength / nPEx;
  double eachPElength_y = totalLength / nPEy;

  double x0 = eachPElength_x * myMPI.iPE;
  double x1 = x0 + eachPElength_x;
  double y0 = eachPElength_y * myMPI.jPE;
  double y1 = y0 + eachPElength_y;

  Mesh MESH(x0, x1, y0, y1, nCellx, nCelly, myMPI);

  // -
  // |
  // | Set up Particles
  // |
  // -

  particles PTCL(vx_bdy, flux, npHat);
  int count = 0;

  // -
  // |
  // | Time Marching Loop
  // |
  // -

  for (double t = 0.; t <= tEnd; t += dt) {
    if (myMPI.myPE == 0) cout << "Time = " << t << endl;

    // Inject particles

    if (myMPI.iPE == 0) PTCL.addFlux(t, dt, MESH.y0, MESH.y1, myMPI.myPE);

    // Move particles

    PTCL.move(dt);

    // Map between particles and the mesh

    MESH.ParticlesOnMesh(PTCL, myMPI);

    // Plot

    PTCL.plot("ptcl", count, myMPI.myPE);
    MESH.plot("mesh", count, myMPI);

    ++count;
  }

  // -
  // |
  // | Wrap-Up
  // |
  // -

  if (myMPI.myPE == 0) cout << "\n\n ** Successful Completion ** \n\n";

  MPI_Finalize();

  return 0;
}
