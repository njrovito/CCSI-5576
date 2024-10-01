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

void plot(string descriptor, VD &field)
  {
    // Convert my processor number to a string:
    
    std::stringstream myPE_str;  myPE_str << 0;
    
    // Create plot filename:
    //    Filename is the following concatenated information:
    //    caller-provided descriptor + myPE + ".plt"

    string filename = descriptor + "_" + myPE_str.str() + ".vtk";

    // Open plot file
    
    std::fstream f;
    f.open(filename.c_str(),std::ios::out);

    // Write to plot file

    
    f << "# vtk DataFile Version 2.0	" << endl;
    f << "Temperature Data"               << endl;
    f << "ASCII           "               << endl;
    f << "DATASET STRUCTURED_POINTS"      << endl;
    f << "SPACING " << dx << " " << dy << " " << dz << endl;
    f << "DIMENSIONS " << ncell_x << " " << ncell_y << " " << ncell_z << endl;
    f << "ORIGIN 0 0 0 " << endl;
    f << "POINT_DATA " << nField << endl;
    f << "SCALARS Temperature float " << endl;
    f << "LOOKUP_TABLE Default " << endl;

    rLOOP f << field[r] << endl;
    
    f.close();
  }

