#include "ch1_cantilever.hh"

ch1_cantilever::ch1_cantilever() : Problem("Cantilever under gravity", "Cantilever under gravity")
, rod(NULL)
, stepper(NULL)
{
    addDynamicsProps();
    addRodOptions();
    addRodTimeStepperOptions();
    
    AddOption("length", "length of rod", 0.1);
    AddOption("save-data", "save data", false);
    AddOption("end-time", "Sim. exits at this time", -1.0);
}

ch1_cantilever::~ch1_cantilever()
{
    if (rod != NULL) delete rod;
    if (stepper != NULL) delete stepper;
    if (saveData) out_file.close();
}

void ch1_cantilever::Setup()
{
    loadDynamicsProps();
    
    RodOptions opts;
    getRodOptions(opts);
    
    object_density = GetScalarOpt("density");
    dt = GetScalarOpt("dt");
    
    l = GetScalarOpt("length");
    a = GetScalarOpt("major-radius");
    Y = GetScalarOpt("youngs-modulus");
    G = GetScalarOpt("shear-modulus");
    endTime = GetScalarOpt("end-time");
    cout << "Y = " << GetScalarOpt("youngs-modulus") << "  G = " << GetScalarOpt("shear-modulus") << std::endl;
    
    saveData = GetBoolOpt("save-data");
    cout << "save data " << saveData << endl;
    
    Scalar edgeLen = l / ((double)(opts.numVertices - 2));
	for (int i=0; i < opts.numVertices-1; i++)
	{
	  Scalar x;
	  if (i==0)
	  {
	    x = - edgeLen;
        vertices.push_back(Vec3d(x, 0, 0));
        undeformed.push_back(Vec3d(x, 0, 0));
	  }
      x = ((double)i) * edgeLen;
	  vertices.push_back(Vec3d(x, 0, 0));
      undeformed.push_back(Vec3d(x, 0, 0));
	}

    rod = setupRod(opts, vertices, undeformed);
    nv = rod->nv();
    ne = rod->ne();
    
    stepper = getRodTimeStepper(*rod);
    m_world->addObject(rod);
    m_world->addController(stepper);

    // Set BC
    stepper->getBoundaryCondition()->setDesiredVertexPosition(0, rod->getVertex(0));
    stepper->getBoundaryCondition()->setDesiredVertexPosition(1, rod->getVertex(1));
	stepper->getBoundaryCondition()->setDesiredEdgeAngle(0, rod->getTheta(0));
    
    if (saveData) OpenFile();
}


void ch1_cantilever::AtEachTimestep()
{
	static int counter = 0;
	ttime = getTime();

    if (saveData && counter%100 == 0)
    {
      for (int i = 0; i < nv; i++)
      {
	    Vec3d pos = rod->getVertex(i);
        out_file << ttime << " " << pos.x() << " " << pos.y() << " " <<
          pos.z() << " " << endl;
	  }
	  if (endTime>0 && ttime > endTime) exit(1); // exit if condition met
    }
    counter++;
}

void ch1_cantilever::AfterEachTimestep()
{
    ;
}


void ch1_cantilever::Render()
{
	;
}

void ch1_cantilever::OpenFile()
{
    name<<"ch1_cantileverData.txt";
    
    out_file.open(name.str().c_str());
    out_file.precision(10);        
    out_file <<"#t(s) x y z\n";
}
