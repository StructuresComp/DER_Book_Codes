//
// Khalid Jawed
// khalidjm@seas.ucla.edu
//

#include "ProblemBase.hh"
#include "stdio.h"

#include "BASim/src/externalForces/dampingForce.hh"
#include "BASim/src/externalForces/terminalLoads.hh"

#include <iostream>
#include <fstream>
using namespace std;


class ch1_helixTerminalMoment : public Problem
{
public:
    ch1_helixTerminalMoment();
    virtual ~ch1_helixTerminalMoment();
    
protected:
    
    void Setup();
    void AtEachTimestep();
    void AfterEachTimestep();
    void Render();
    void OpenFile(); //A utility subroutine to write the rod configuration
    
    ElasticRod* rod;
    RodTimeStepper* stepper;
    std::vector<Vec3d> vertices, undeformed;
    
    Scalar Phi; // Total twist
    Scalar dTheta; 
    Scalar dt; // time step
    int nv, ne;
    Scalar l, a, Y, G;
    Scalar v1_i, v1_f; // curvature 1
    Scalar v2_i, v2_f; // curvature 2
    Scalar v3_i, v3_f; // twist at two ends (per length)
    
    std::ofstream m_outputfile;
    
    // Forces
    //ViscousForce* m_viscous;
    dampingForce* m_dampingforce;
    terminalLoads* tLoads;
    Scalar endKappa, endTwist;
    Scalar EI, GJ;

    Scalar m_viscosity;
    
    Scalar object_density; //Density of rod
    std::ostringstream name;
    std::ofstream out_file;
    
    Scalar ttime; // time
    bool saveData; // should a file be opened to save the data
    Scalar endTime; // simulation exits at this time
};

