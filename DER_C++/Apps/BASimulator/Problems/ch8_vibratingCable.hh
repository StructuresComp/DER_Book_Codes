#include "ProblemBase.hh"
#include "stdio.h"
#include "BASim/src/externalForces/dampingForce.hh"

#include <iostream>
#include <fstream>
using namespace std;

class ch8_vibratingCable : public Problem
{
public:
    ch8_vibratingCable();
    virtual ~ch8_vibratingCable();
    
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
    Scalar lastnodeX, lastnode_1X;
    
    // Forces
    //ViscousForce* m_viscous;
    dampingForce* m_dampingforce;
    Scalar m_viscosity;
    Vec3d m_gravity;
    
    Scalar object_density; //Density of rod
    std::ostringstream name;
    std::ofstream out_file;
    std::ostringstream nameEnergy;
    std::ofstream out_fileEnergy;
    
    Scalar ttime; // time
    bool saveData; // should a file be opened to save the data
    
    Scalar endTime; // simulation exits at this time
    
    Scalar amplitude, frequencyStart, frequencyEnd, tStep; // vibration
};

