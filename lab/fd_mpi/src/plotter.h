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

void plot(string descriptor, VD &field, mpiInfo &myMPI)
  {
    // Convert my processor number to a string:
    
    std::stringstream myPE_str;  myPE_str << myMPI.myPE;
    
    // Create plot filename:
    //    Filename is the following concatenated information:
    //    caller-provided descriptor + myPE + ".plt"

    string filename = descriptor + "_" + myPE_str.str() + ".plt";

    // Open plot file
    
    std::fstream f;
    f.open(filename.c_str(),std::ios::out);

    // Write to plot file

    int p[6];
    
    for ( int i = 1 ; i <= nRealx-1 ; ++i )
      {
	for ( int j = 1 ; j <= nRealy-1; ++j )
	  {
	    p[1] = pid( i   , j   );
	    p[2] = pid( i+1 , j   );
	    p[3] = pid( i+1 , j+1 );
	    p[4] = pid( i   , j+1 );
	    p[5] = p[1]            ;
	    
	    for ( int k = 1 ; k <= 5 ; ++k )  f << x[p[k]] << " " << y[p[k]] << " " << field[ p[k] ] << endl;
	  }
	
	f << endl;
	f << x[ p[1] ] << " " << y[ p[1] ] << " " << field[ p[1] ] << endl;
	f << endl;

      }

    f.close();


    // ----------------------------------------------
    // Write the gnuplot plot command file 'pc'
    // ----------------------------------------------

    if ( myMPI.myPE == 0 )
      {
	cout << "(o) Writing gnuplot command file 'pc'" << endl;
	
	f.open("pc",std::ios::out);
	f << "#" << endl;
	f << "# This file is written by fd_mpi.  See plotter.h" << endl;
	f << "# To run this file in gnuplot: " << endl;
	f << "#   (1) Start gnuplot" << endl;
	f << "#   (2) At the gnuplot prompt, enter the command:" << endl;
	f << "#          load 'pc'" << endl;
	f << "#" << endl;
	f << "# The current version of this plot command file will display" << endl;
	f << "# the plot to the screen.  To write the plot to a png file," << endl;
	f << "# uncomment the two 'set' commands in the following lines." << endl;
	f << "#" << endl;
	f << "#set term png" << endl;
	f << "#set output 'plot.png'" << endl;

	f << endl;
	f << endl;

	string plotCmd = "splot ";
	for ( int i = 0 ; i < myMPI.numPE ; ++i )
	  {
	    std::stringstream myPE_str2;  myPE_str2 << i;
	    plotCmd += "'F_phi_" + myPE_str2.str() + ".plt' w l, ";
	  }
	f << plotCmd << endl;
	f << endl;
	f.close();
      }
    
    
  }

