//  ================================================================================
//  ||                                                                            ||
//  ||              fd                                                            ||
//  ||              -------------------------------------------                   ||
//  ||              F I N I T E   D I F F E R E N C E                             ||
//  ||                                                                            ||
//  ||              D E M O N S T R A T I O N   C O D E                           ||
//  ||              -------------------------------------------                   ||
//  ||                                                                            ||
//  ||       Developed by: Scott R. Runnels, Ph.D.                                ||
//  ||                     University of Colorado Boulder                         ||
//  ||                                                                            ||
//  ||                For: CU Boulder CSCI 4576/5576 and associated labs          ||
//  ||                                                                            ||
//  ||           Copyright 2020 Scott Runnels                                     ||
//  ||                                                                            ||
//  ||                     Not for distribution or use outside of the             ||
//  ||                     this course.                                           ||
//  ||                                                                            ||
//  ================================================================================

  //  ==
  //  ||
  //  ||  Perform Jacobiiterations on nField X nField linear
  //  ||  system:  A*Solution = RHS
  //  ||
  //  ==

void SolveLinearSystem(int max_iter , VD RHS, VD &Solution )
  {
    int converged, it_converged;
    int iter = 0;
    double newval;
    double cur_delta = 0.;
    double max_delta = 0.;
    double tol;

    int        err;
    int        global_converged;
    int        zero = 0;
    int        one  = 1;
    
    VD SolutionNew; SolutionNew.resize(Solution.size());

    rLOOP Solution[r] = 0.;
    tol       = 1.e-04;
    
    // ====================================================
    // Begin Iterations
    // ====================================================
    
    converged = 0;

    
    while ( converged == 0 )
      {

	// ----------------------------------------------
	// (1) Convergence test
	// ----------------------------------------------
	
	if ( ++iter > max_iter )
	  {
	    cout << "Jacobi: max_iter of " << max_iter << " reached.\n";
	    return;
	  }

	max_delta    = 0.;  it_converged = 1;

	// ----------------------------------------------
	// (2) One Jacobi Iteration
	// ----------------------------------------------
	
	rLOOP
	  {

	    // (2.1) Compute new guess for row r

	    newval = b[r];
	    cLOOP  if ( r != c ) newval -=  A[r][c] * Solution[c];
	    newval /= A[r][r];

	    // (2.2) Convergence check

	    cur_delta  = fabs(Solution[r] - newval);

	    if ( cur_delta > tol ) it_converged = 0;

	    // (2.3) Record new value in solution

	    SolutionNew[r]       = newval;
	      
	  }

	rLOOP Solution[r] = SolutionNew[r];

	// ----------------------------------------------
	// (3) Make note of the convergence state
	// ----------------------------------------------
	
	converged = it_converged;

	// ----------------------------------------------
	// (4) Gather convergence information from PEs
	// ----------------------------------------------

	cout << "Jacobi iteration " << iter << endl;
      }
	
    cout << "Jacobi iteration converged " << endl;
  }


