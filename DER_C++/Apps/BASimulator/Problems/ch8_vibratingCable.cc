#include "ch8_vibratingCable.hh"

ch8_vibratingCable::ch8_vibratingCable() : Problem("End vibration", "A rod with imposed vibration at the end")
, rod(NULL)
, stepper(NULL)
{
    addDynamicsProps();
    addRodOptions();
    addRodTimeStepperOptions();
    
    AddOption("length", "length of rod", 0.1);
    AddOption("medium-viscosity", "medium viscosity", 0.0);
    AddOption("save-data", "save data", false);
    AddOption("end-time", "Sim. exits at this time", -1.0);
    AddOption("frequency-start", "Frequency of vibration", 1.0);
    AddOption("frequency-end", "Frequency of vibration", -1.0);
    AddOption("t-step", "time over which the frequency varies from start to end", -1.0);
    AddOption("amplitude", "Amplitude of vibration", -1.0);
}

ch8_vibratingCable::~ch8_vibratingCable()
{
    if (rod != NULL) delete rod;
    if (stepper != NULL) delete stepper;
    if (saveData) 
    {
      out_file.close();
      out_fileEnergy.close();
    }
}

void ch8_vibratingCable::Setup()
{
    loadDynamicsProps();
    
    RodOptions opts;
    getRodOptions(opts);
    
    m_viscosity = GetScalarOpt("medium-viscosity");
    object_density = GetScalarOpt("density");
    dt = GetScalarOpt("dt");
    
    l = GetScalarOpt("length");
    a = GetScalarOpt("major-radius");
    Y = GetScalarOpt("youngs-modulus");
    G = GetScalarOpt("shear-modulus");
    m_gravity = GetVecOpt("gravity");
    amplitude = GetScalarOpt("amplitude");
    frequencyStart = GetScalarOpt("frequency-start");
    frequencyEnd = GetScalarOpt("frequency-end");
    if (frequencyEnd < 0)
    {
	  tStep = 1.0; // doesn't matter
      frequencyEnd = frequencyStart;
    }

    tStep = GetScalarOpt("t-step");
    
    endTime = GetScalarOpt("end-time");
    cout << "l = " << l << endl;
    cout << "Y = " << Y << "  G = " << G << std::endl;
    cout << "amplitude = " << amplitude << "  frequency = " << frequencyStart 
      << " frequency-end = " << frequencyEnd << " tstep = " << tStep << endl;
    
    // Output area moment of inertia
    Scalar AreaMoment = M_PI * pow(a, 4.0) / 4.0;
    cout << "I = " << AreaMoment << endl;
    // Output force per length
    Vec3d gVector = GetVecOpt("gravity");
    Scalar rho = GetScalarOpt("density");
    Scalar g = gVector.norm();
    Scalar w = rho * M_PI * square(a) * g;
    cout << "load per length: " << w << endl;
    
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
    
    //Add damping
    if(m_viscosity > 0)
    {
        m_dampingforce = new dampingForce(m_viscosity, ttime);
        stepper->addExternalForce(m_dampingforce);
    }
        
    m_world->addObject(rod);
    m_world->addController(stepper);

    // Set BC
    stepper->getBoundaryCondition()->setDesiredVertexPosition(0, rod->getVertex(0));
    stepper->getBoundaryCondition()->setDesiredVertexPosition(1, rod->getVertex(1));
	stepper->getBoundaryCondition()->setDesiredEdgeAngle(0, rod->getTheta(0));
    
    if (saveData) OpenFile();
    
    lastnodeX = (rod->getVertex(nv-1)).x();
    lastnode_1X = (rod->getVertex(nv-2)).x();
}


void ch8_vibratingCable::AtEachTimestep()
{
	static int counter = 0;
	ttime = getTime();

    //
    // Apply vibration
    //
    Scalar translationEndTime = 1.0;
    Scalar fEffective = frequencyStart;
    if (ttime <= translationEndTime)
      fEffective = 0.0;
    else if (ttime <= tStep + translationEndTime)
      fEffective = frequencyStart + ( (ttime - translationEndTime)/tStep) * (frequencyEnd - frequencyStart);
    else
      fEffective = frequencyEnd;
    Scalar verticalPos = amplitude * sin(2.0 * M_PI * fEffective * ttime);
    
    fEffective = round(fEffective);

    // last node 
    Vec3d endNode = rod->getVertex(nv-1);
    Scalar desiredX;
    if (ttime < translationEndTime)
      desiredX = lastnodeX - 0.1 * ttime * l;
    else
      desiredX = lastnodeX - 0.1 * translationEndTime * l;
    Vec3d desiredPos = Vec3d(desiredX, verticalPos, endNode.z());
    stepper->getBoundaryCondition()->setDesiredVertexPosition(nv-1, desiredPos);

    // second-last node
    Vec3d endNode2 = rod->getVertex(nv-2);
    if (ttime < translationEndTime)
      desiredX = lastnode_1X - 0.1 * ttime * l;
    else
      desiredX = lastnode_1X - 0.1 * translationEndTime * l;
    Vec3d desiredPos2 = Vec3d(desiredX, verticalPos, endNode2.z());
    stepper->getBoundaryCondition()->setDesiredVertexPosition(nv-2, desiredPos2);

    if (saveData)
    {
      // Energy computation
      Scalar bendingE, stretchingE;
      std::vector<RodForce*>& forces = rod->getForces();
      std::vector<RodForce*>::iterator fIt;
      for (fIt = forces.begin(); fIt != forces.end(); ++fIt) 
      {
        if (dynamic_cast<RodBendingForceSym*> (*fIt) != NULL ) 
        {
          RodBendingForceSym *f = dynamic_cast<RodBendingForceSym*> (*fIt);
          if (!f->viscous()) bendingE = f->globalEnergy();
        }
        if (dynamic_cast<RodStretchingForce*> (*fIt) != NULL ) 
        {
          RodStretchingForce *f = dynamic_cast<RodStretchingForce*> (*fIt);
          if (!f->viscous()) stretchingE = f->globalEnergy();
        }
      }
      // Potential energy
      Scalar gravityE = 0;
      for (int i = 0; i < nv; i++)
      {
	    Vec3d pos = rod->getVertex(i);
        gravityE -= pos.dot(m_gravity) * rod->getVertexMass(i);
      }
      // Kinetic energy
      Scalar kineticE = 0;
      for (int i = 0; i < nv; i++)
      {
	    Vec3d v = rod->getVelocity(i);
        kineticE += rod->getVertexMass(i) * square(v.norm());
      }
      
      for (int i = 0; i < nv; i++)
      {
	    Vec3d pos = rod->getVertex(i);
        out_file << ttime << " " << fEffective << " " << pos.x() << 
          " " << pos.y() << " " << pos.z() << endl;
	  }
	  
     out_fileEnergy << ttime << " " << fEffective <<
        " " << bendingE << " " << stretchingE << 
        " " << gravityE << " " << kineticE << endl;
          
	  if (endTime>0 && ttime > endTime) exit(1); // exit if condition met
    }
    counter++;
}

void ch8_vibratingCable::AfterEachTimestep()
{
    ;
}


void ch8_vibratingCable::Render()
{
	;
}

void ch8_vibratingCable::OpenFile()
{
    name << "ch8_vibratingCableData.txt";
    
    nameEnergy << "ch8_vibratingCableEnergyData.txt";

    out_file.open(name.str().c_str());
    out_file.precision(10);        
    out_file <<"#t(s) f(Hz) x y z\n";
    
    out_fileEnergy.open(nameEnergy.str().c_str());
    out_fileEnergy.precision(10);        
    out_fileEnergy <<"#t(s) f(Hz) Eb Es Eg Ek (unit:Joules)\n";
}
