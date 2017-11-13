#include "ch5_helixUncoiling.hh"

ch5_helixUncoiling::ch5_helixUncoiling() : Problem("Helical rod under gravity", "Helical rod uncoils due to gravity")
, rod(NULL)
, stepper(NULL)
{
    addDynamicsProps();
    addRodOptions();
    addRodTimeStepperOptions();
    
    AddOption("length", "length of rod", 0.1);
    AddOption("save-data", "save data", false);
    AddOption("end-time", "Sim. exits at this time", -1.0);
    AddOption("natural-radius", "helix parameter", naturalRad);
    AddOption("alpha", "helix parameter", alpha);
}

ch5_helixUncoiling::~ch5_helixUncoiling()
{
    if (rod != NULL) delete rod;
    if (stepper != NULL) delete stepper;
    if (saveData) out_file.close();
}

void ch5_helixUncoiling::Setup()
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
    naturalRad = GetScalarOpt("natural-radius");
    alpha = GetScalarOpt("alpha");    
    endTime = GetScalarOpt("end-time");
    cout << "Y = " << GetScalarOpt("youngs-modulus") << "  G = " << GetScalarOpt("shear-modulus") << endl;
    cout << "Natural Radius = " << naturalRad << "  alpha = " << alpha << endl;
    
    saveData = GetBoolOpt("save-data");
    cout << "save data " << saveData << endl;
    
	for (int i=0; i < opts.numVertices; i++)
	{
	  Scalar thTmp = ((double)i) / ((double) opts.numVertices - 1 ) * l / naturalRad;
	  Scalar x = naturalRad * (1.0 - cos(thTmp));
	  Scalar y = naturalRad * sin(thTmp);
	  Scalar z = - alpha * naturalRad * thTmp;
	  
	  vertices.push_back(Vec3d(x, y, z));
      undeformed.push_back(Vec3d(x, y, z));
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


void ch5_helixUncoiling::AtEachTimestep()
{
	ttime = getTime();

    if (saveData)
    {
      for (int i = 0; i < nv; i++)
      {
		Scalar theta = 0.0;
		if (i < nv-1) theta = rod->getTheta(i);
		Scalar refTwist = 0.0;
		if (i>=1 && i < nv-1) refTwist = rod->getReferenceTwist(i);
		
	    Vec3d pos = rod->getVertex(i);
        out_file << ttime << " " << pos.x() << " " << pos.y() << " " <<
          pos.z() << " " << theta << " " << refTwist << endl;
	  }
	  if (endTime>0 && ttime > endTime) exit(1); // exit if condition met
    }
}

void ch5_helixUncoiling::AfterEachTimestep()
{
    ;
}


void ch5_helixUncoiling::Render()
{
	;
}

void ch5_helixUncoiling::OpenFile()
{
    name<<"ch5_helixUncoilingData.txt";
    
    out_file.open(name.str().c_str());
    out_file.precision(10);        
    out_file <<"#t(s) x y z theta refTwist\n";
}
