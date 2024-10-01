//  ================================================================================
//  ||                                                                            ||
//  ||              esPIC                                                         ||
//  ||              ------------------------------------------------------        ||
//  ||              E L E C T R O S T I C   P A R T I C L E - I N - C E L L       ||
//  ||                                                                            ||
//  ||              D E M O N S T R A T I O N   C O D E                           ||
//  ||              ------------------------------------------------------        ||
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


class particles
{

 public:
  int n;
  VD x;
  VD y;
  VD vx;
  VD vy;
  VI active;
  VD xf, yf;
  double fluxTimer;

  // Particle physics
  
  double Qp;                     // Extensive physical quantity represented by each particle (C/ptcl)  
  double rho;                    // Charge Density (C/m^3)
  double f;                      // Flux Density (C/m^2)
  double vx_bdy;                 // x-velocity of injected particles (m/sec)
  double vy_bdy;                 // y-velocity of injected particles (m/sec)
  int    np_dot;                 // Particle flux rate (ptcls/m^2/sec)
  double nHat;                   // Number of particles/area to be injected, when injection occurs.
                                 //   - This user input ensures a smooth spread of particles when
                                 //     and injection occurs.
  double injectionInterval;      // Injection interval (sec)
  
  double previousInjectionTime;

  // ==
  // ||
  // ||  Constructor
  // ||
  // ==

  particles(double _vx_bdy, double _f , double _nHat )
    {
      n  = 5000;
      x.resize(n+1);
      y.resize(n+1);
      vx.resize(n+1);
      vy.resize(n+1);
      xf.resize(n+1);
      yf.resize(n+1);
      active.resize(n+1);

      for ( int i = 1 ; i < active.size() ; ++i ) active[i] = 0;
      fluxTimer = 0.;


      // Capture particle physics
      
      Qp     = 1.;           
      vx_bdy = _vx_bdy;
      vy_bdy = 0.;
      f      = _f;            
      nHat   = _nHat;

      // Compute remaining physical variables

      np_dot            = f/Qp;
      injectionInterval = nHat / np_dot;

      // Initialize particle injection timer
      
      previousInjectionTime = 0.;
    }

  // ==
  // ||
  // ||  Adds new particles to the list of particles
  // ||
  // ==

  void add(VD &xAdd , VD &yAdd, VD &vxAdd , VD &vyAdd)
  {

    // Add new particles in spaces left inactive (active[j] != 1)
    // to the existing arrays.

    int count = 0;
    for ( int j = 1 ; j <= x.size()-1 ; ++j )
      {
	if ( active[j] == 0 )
	  {
	    ++count;
	    
	    if ( count > xAdd.size()-1 )
	      {
		return;
	      }
	    
	    x [j] = xAdd [count];
	    y [j] = yAdd [count];
	    vx[j] = vxAdd[count];
	    vy[j] = vyAdd[count];
	    active[j] = 1;

	  }
      }
    
    // If we are out of space in the existing arrays, resize
    // them and try again.

    FatalError("Out of room");

  }

  
  // ==
  // ||
  // ||  Moves particles the (provided) dt
  // ||
  // ==

  void move(double dt)
  {
    for ( int i = 1 ; i < x.size() ; ++i )
      {
      if ( active[i] == 1 )
	{

	  // Apply force to update velocity with acceleration

	  vx[i] += dt * xf[i];
	  vy[i] += dt * yf[i];

	  // Apply velocity to update position
	  
	  x [i] += dt * vx[i];
	  y [i] += dt * vy[i];
	}
      }
  }

  // ==
  // ||
  // ||  ------------------------
  // ||  Function addFlux
  // ||  ------------------------
  // ||
  // ||  Regarding particles:
  // ||
  // ||
  // ||
  // ==

  void addFlux(double time, double dt , double y0 , double y1 , int myPE)
  {

    if ( injectionInterval < dt )
      {
	cout << "injectionInterval = " << injectionInterval << endl;
	FatalError("injectionInterval is too small.  Increase the particle injection density or decrease dt.");
      }

    double howLateWeAre = (time - previousInjectionTime) - injectionInterval;

    if ( howLateWeAre < 0. ) return;

    // So we are late!  Calculate how far into the mesh the particles have traveled
    // since entering.

    double distanceIn = howLateWeAre * vx_bdy;

    // Calculate the number of particles to be injected

    int Np = int ( nHat * (y1-y0) ) ;

    // Calculcate the vertical spacing of the particles

    double dy = (y1-y0) / Np;  
                               
    // Allocate memory to store the new particles
    
    VD xAdd ;  xAdd.resize(Np+1);    VD yAdd ;  yAdd.resize(Np+1);
    VD vxAdd; vxAdd.resize(Np+1);    VD vyAdd; vyAdd.resize(Np+1);

    // Compute their initial position

    for ( int i = 1; i <=  Np ; ++i )
      {
	xAdd[i] = distanceIn; 
	yAdd[i] = y0 + dy/2. + dy*(i-1);
      }

    // Store their initial velocity

    for ( int i = 1 ; i <= Np ; ++i )
      {
	vxAdd[i] = vx_bdy;
	vyAdd[i] = vy_bdy;
      }
    
    // Record the time that this set of particles entered the mesh.  Note that
    // we are recording the actual time, which is the current time minus how
    // late we are.

    previousInjectionTime = time - howLateWeAre;

    // Add the particles

    add(xAdd , yAdd, vxAdd , vyAdd);
    
  }

  


#include "particle_plotter.h"


};


