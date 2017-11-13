//
// Khalid Jawed
// khalidjm@seas.ucla.edu
//

#include "ProblemBase.hh"
#include "stdio.h"

#include <iostream>
#include <fstream>
using namespace std;

class ch5_helixUncoiling : public Problem
{
public:
    ch5_helixUncoiling();
    virtual ~ch5_helixUncoiling();
    
protected:
    
    void Setup();
    void AtEachTimestep();
    void AfterEachTimestep();
    void Render();
    void OpenFile(); //A utility subroutine to write the rod configuration
    
    ElasticRod* rod;
    RodTimeStepper* stepper;
    std::vector<Vec3d> vertices, undeformed;
    
    int nv, ne;
    Scalar l, a, Y, G, dt;
    
    Scalar naturalRad, alpha; // helix parameters
    
    Scalar object_density; //Density of rod
    std::ostringstream name;
    std::ofstream out_file;
    
    Scalar ttime; // time
    bool saveData; // should a file be opened to save the data
    
    Scalar endTime; // simulation exits at this time
};

