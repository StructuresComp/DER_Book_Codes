#include "ch1_helixTerminalMoment.hh"

ch1_helixTerminalMoment::ch1_helixTerminalMoment() : Problem("Rod under terminal moments", "A straight turns into a helical one")
, rod(NULL)
, stepper(NULL)
{
    addDynamicsProps();
    addRodOptions();
    addRodTimeStepperOptions();
    
    AddOption("length", "length of rod", 0.1);
    AddOption("end-curvature1", "end curvature 1", 0.0);
    AddOption("end-curvature2", "end curvature 2", 0.0);
    AddOption("end-twist", "end twist", 0.0);
    AddOption("end-time", "Sim. exits at this time", -1.0);
    AddOption("medium-viscosity", "medium viscosity", 0.0);
    AddOption("save-data", "save data", false);
}

ch1_helixTerminalMoment::~ch1_helixTerminalMoment()
{
    if (rod != NULL) delete rod;
    if (stepper != NULL) delete stepper;
    if (saveData) out_file.close();
}

void ch1_helixTerminalMoment::Setup()
{
    loadDynamicsProps();
    
    RodOptions opts;
    getRodOptions(opts);
    
    m_viscosity = GetScalarOpt("medium-viscosity");
    object_density = GetScalarOpt("density");
    dt = GetScalarOpt("dt");
    v1_i = GetScalarOpt("end-curvature1"); // unit: 1/length
    v2_i = GetScalarOpt("end-curvature2"); // unit: 1/length
    v3_i = GetScalarOpt("end-twist"); // unit: 1/length
    
    l = GetScalarOpt("length");
    a = GetScalarOpt("major-radius");
    Y = GetScalarOpt("youngs-modulus");
    G = GetScalarOpt("shear-modulus");
    endTime = GetScalarOpt("end-time");
    cout << "Y = " << GetScalarOpt("youngs-modulus") << "  G = " << GetScalarOpt("shear-modulus") << std::endl;
    cout << "End curvature = " << v1_i << "," << v2_i << " end-twist = " << v3_i << endl;
    
    saveData = GetBoolOpt("save-data");
    cout << "save data " << saveData << endl;
    
	for (int i=0; i < opts.numVertices; i++)
	{
	  Scalar x = ((double)i) * l / ((double)(opts.numVertices - 1));
	  vertices.push_back(Vec3d(x, 0, 0));
      undeformed.push_back(Vec3d(x, 0, 0));
	}

    rod = setupRod(opts, vertices, undeformed);
    nv = rod->nv();
    ne = rod->ne();
    
    stepper = getRodTimeStepper(*rod);
    
    //Add damping
    if(m_viscosity > 0)
    {
        m_dampingforce = new dampingForce(m_viscosity, ttime);
        stepper->addExternalForce(m_dampingforce);
    }
    
    EI = Y * M_PI * pow(a,4) / 4.0;
    GJ = G * M_PI * pow(a,4) / 2.0;
    // Compute twist and curvature at the other end
    v3_f = v3_i;    
    Scalar beta = 2.0 * G / Y;
    Scalar deltaEps = (1 - beta) * v3_i * l;
    v1_f = cos(deltaEps) * v1_i + sin(deltaEps) * v2_i;
    v2_f = - sin(deltaEps) * v1_i + cos(deltaEps) * v2_i;  
    
    tLoads = new terminalLoads(v1_i, v2_i, v3_i, v1_f, v2_f, v3_f, EI, GJ, ttime);
    stepper->addExternalForce(tLoads);
    
    m_world->addObject(rod);
    m_world->addController(stepper);

    // Set BC
    stepper->getBoundaryCondition()->setDesiredVertexPosition(0, rod->getVertex(0));
    
    if (saveData) OpenFile();
}


void ch1_helixTerminalMoment::AtEachTimestep()
{
	static int counter = 0;
	ttime = getTime();
	
	//calculate twist.
    Phi = 0.0;
    Scalar twist;
    for (int i = 0; i < nv; i++)
    {
	  Vec3d pos = rod->getVertex(i);
	  if (i==0 || i==nv-1)
	    twist = 0.0;
	  else
	    twist = rod->getReferenceTwist(i) + rod->getTheta(i) - rod->getTheta(i - 1);

      Phi += twist;
      if (saveData && counter%10 == 0)
      {
        out_file << ttime << " " << pos.x() << " " << pos.y() << " " <<
          pos.z() << " " << twist << endl;
	  }
    }

    if (endTime>0 && ttime > endTime) exit(1); // exit if condition met
    counter++;
}

void ch1_helixTerminalMoment::AfterEachTimestep()
{
    ;
}


void ch1_helixTerminalMoment::Render()
{
	;
}

void ch1_helixTerminalMoment::OpenFile()
{
    name << "ch1_helixTerminalMomentData.txt";

    out_file.open(name.str().c_str());
    out_file.precision(10);        
    out_file <<"#t(s) x y z twist\n";
}
