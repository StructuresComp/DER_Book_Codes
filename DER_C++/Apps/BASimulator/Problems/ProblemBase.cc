/**
 * \file ProblemBase.cc
 *
 * \author miklos@cs.columbia.edu (based on problem-setup.cc from sar2120@columbia.edu)
 * \date 09/09/2009
 */

#include "ProblemBase.hh"
#include <BASim/BASim>
#include <fstream>

using namespace std;
using namespace BASim;

Problem::Problem(const string& name, const string& desc)
  : m_problemName(name)
  , m_problemDesc(desc)
  , m_dynamicsProps(false)
{
msg_val = 0.0;
  m_world = new World();
  
  update_display = false;

}

Problem::~Problem()
{
  assert( m_world != NULL );
  if( m_world != NULL ) delete m_world;
  m_world = NULL;
}

Option* Problem::GetOption(const string& name)
{
  if (m_options.find(name) == m_options.end()) {
    cerr << "Option " << name << " does not exist" << endl;
    //exit(-1);
  }

  return &(m_options.find(name)->second);
}

bool& Problem::GetBoolOpt(const string& name)
{
  return GetOption(name)->b;
}

int& Problem::GetIntOpt(const string& name)
{
  return GetOption(name)->i;
}

Scalar& Problem::GetScalarOpt(const string& name)
{
  return GetOption(name)->r;
}

Vec3d& Problem::GetVecOpt(const string& name)
{
  return GetOption(name)->v;
}

string& Problem::GetStringOpt(const string& name)
{
  return GetOption(name)->s;
}

int Problem::LoadOptions(const char* filename)
{
  ifstream input(filename);
  if (!input.is_open()) {
    cerr << "ERROR: File " << filename << " not found" << endl;
    return -1;
  }

  string line, option;
  istringstream sIn;
  string tmp;
  for (getline(input, line); !input.eof(); getline(input, line)) {
    sIn.clear();
    option.clear();
    sIn.str(line);
    sIn >> option;
    if (option.size() == 0 || option.c_str()[0] == '#') continue;
    OptionMap::iterator itr;
    itr = m_options.find(option);
    if (itr == m_options.end()) {
      cerr << "Invalid option: " << option << endl;
      continue;
    }
    if (itr->second.type == Option::BOOL) {
      sIn >> tmp;
      if (tmp == "true" || tmp == "1") itr->second.b = true;
      else if (tmp == "false" || tmp == "0") itr->second.b = false;
    } else if (itr->second.type == Option::INT) {
      sIn >> itr->second.i;
    } else if (itr->second.type == Option::SCALAR) {
      sIn >> itr->second.r;
    } else if (itr->second.type == Option::VEC) {
      Vec3d& v = itr->second.v;
      sIn >> v[0];
      sIn >> v[1];
      sIn >> v[2];
    } else if (itr->second.type == Option::STRING) {
      sIn >> itr->second.s;
    } else {
      cerr << "Invalid option type" << endl;
    }
  }
  input.close();

  return 0;
}

int Problem::LoadOptions(int argc, char** argv)
{
  string option, tmp;
  int start = 0;
  while (start < argc && string(argv[start]) != "--") ++start;
  for (int i = start + 1; i < argc; ++i) {
    option = argv[i];
    OptionMap::iterator itr;
    itr = m_options.find(option);
    if (itr == m_options.end()) {
      cerr << "Invalid option on command line: " << option << endl;
      continue;
    }
    if (i == argc - 1) {
      cerr << "Too few arguments on command line" << endl;
      break;
    }
    if (itr->second.type == Option::BOOL) {
      tmp = argv[i+1]; ++i;
      if (tmp == "true" || tmp == "1") itr->second.b = true;
      if (tmp == "false" || tmp == "0") itr->second.b = false;
    } else if (itr->second.type == Option::INT) {
      itr->second.i = atoi(argv[i+1]); ++i;
    } else if (itr->second.type == Option::SCALAR) {
      itr->second.r = atof(argv[i+1]); ++i;
    } else if (itr->second.type == Option::VEC) {
      if (i >= argc - 3) {
        cerr << "Too few arguments on command line" << endl;
        break;
      }
      Vec3d& v = itr->second.v;
      v[0] = atof(argv[i+1]); ++i;
      v[1] = atof(argv[i+1]); ++i;
      v[2] = atof(argv[i+1]); ++i;
    } else if (itr->second.type == Option::STRING) {
      itr->second.s = argv[i+1]; ++i;
    } else {
      //cerr << "Invalid option type" << endl;
    }
  }
  return 0;
}

void Problem::BaseSetup(int argc, char** argv)
{
#ifdef HAVE_PETSC
  PetscUtils::initializePetsc(&argc, &argv);
#endif // HAVE_PETSC
  Setup();
}

void Problem::BaseFinalize()
{
#ifdef HAVE_PETSC
  PetscUtils::finalizePetsc();
#endif // HAVE_PETSC
}

void Problem::BaseAtEachTimestep()
{
//  cout << "Problem::BaseAtEachTimestep()\n";
  
  AtEachTimestep();
  m_world->execute();

  if (m_dynamicsProps) {
    setTime(getTime() + getDt());
  }
  AfterStep();

//  cout << "Problem::BaseAtEachTimestep() comp \n";

}

void Problem::BaseRender()
{
  Render();
}

void Problem::BaseHUDRender(int w, int h)
{
  HUDRender(w, h);
}

void Problem::PrintOptions(ostream& os)
{
  OptionMap::const_iterator itr;
  for (itr = m_options.begin(); itr != m_options.end(); ++itr) {
    const Option* opt = &itr->second;
    os << opt->name << " ";
    switch(opt->type) {
    case Option::BOOL:
      if (opt->b) {
        os << "true";
      } else {
        os << "false";
      }
      break;
    case Option::INT:
      os << opt->i;
      break;
    case Option::SCALAR:
      os << opt->r;
      break;
    case Option::VEC:
      os << opt->v.format(EIGEN_SPACES_ONLY_IO);
      break;
    case Option::STRING:
      os << opt->s;
      break;
    }
    os << "\t# " << opt->label << endl;
  }
}

void Problem::addDynamicsProps()
{
  if (m_dynamicsProps) {
    cerr << "Dynamics properties already exist" << endl;
    return;
  }
  m_dynamicsProps = true;

  m_world->add_property(m_time, "time", 0.0);

  m_world->add_property(m_dt, "time-step", 0.01);
  AddOption("dt", "time-step size", getDt());

  m_world->add_property(m_gravity, "gravity", Vec3d(0, -981, 0));
  AddOption("gravity", "gravity", getGravity());
}

void Problem::loadDynamicsProps()
{
  assert(dynamicsPropsLoaded());

  setDt(GetScalarOpt("dt"));
  setGravity(GetVecOpt("gravity"));
}

const Scalar& Problem::getTime() const
{
  assert(dynamicsPropsLoaded());

  return m_world->property(m_time);
}

void Problem::setTime(const Scalar& time)
{
  assert(dynamicsPropsLoaded());

  m_world->property(m_time) = time;
}

const Scalar& Problem::getDt() const
{
  assert(dynamicsPropsLoaded());

  return m_world->property(m_dt);
}

void Problem::setDt(const Scalar& dt)
{
  assert(dynamicsPropsLoaded());

  m_world->property(m_dt) = dt;
}

const Vec3d& Problem::getGravity() const
{
  assert(dynamicsPropsLoaded());

  return m_world->property(m_gravity);
}

void Problem::setGravity(const Vec3d& gravity)
{
  assert(dynamicsPropsLoaded());

  m_world->property(m_gravity) = gravity;
}

void Problem::addRodOptions()
{
  RodOptions opts;
  AddOption("nv", "number of vertices in the rod", opts.numVertices);
  AddOption("density", "volumetric density of the rod", opts.density);
  AddOption("youngs-modulus", "Young's modulus for the rod", opts.YoungsModulus);
  AddOption("shear-modulus", "Shear modulus for the rod", opts.ShearModulus);
  AddOption("viscosity", "Dynamic viscosity for the rod", opts.viscosity);
  AddOption("major-radius", "radius of major axis of cross-section", opts.radiusA);
  AddOption("minor-radius", "radius of minor axis of cross-section", opts.radiusB);
  AddOption("radius-scale", "scaling for rendering and collisions", opts.radiusScale);
  AddOption("reference-frame", "type of reference frame to use (time/space)",
            opts.refFrame == ElasticRod::TimeParallel ? "time" : "space");
  AddOption("quasistatic", "update material frame quasistatically", opts.quasistatic);
}

void Problem::getRodOptions(RodOptions& opts)
{
  opts.numVertices = GetIntOpt("nv");
  opts.density = GetScalarOpt("density");
  opts.YoungsModulus = GetScalarOpt("youngs-modulus");
  opts.ShearModulus = GetScalarOpt("shear-modulus");
  opts.viscosity = GetScalarOpt("viscosity");
  opts.radiusA = GetScalarOpt("major-radius");
  opts.radiusB = GetScalarOpt("minor-radius");
  opts.radiusScale = GetScalarOpt("radius-scale");
  string& frame = GetStringOpt("reference-frame");
  if (frame == "time") opts.refFrame = ElasticRod::TimeParallel;
  else if (frame == "space") opts.refFrame = ElasticRod::SpaceParallel;
  else {
    cerr << "Unknown reference frame type " << frame << ". Using default"
         << endl;
  }
  opts.quasistatic = GetBoolOpt("quasistatic");
}

void Problem::addRodTimeStepperOptions()
{
  AddOption("mass-damping", "mass damping for the rod", 0.0);
  AddOption("integrator", "type of integrator to use for the rod", "implicit");
  AddOption("iterations", "maximum number of iterations for implicit method",
            (int) 50);
  AddOption("atol", "absolute convergence tolerance", 1e-8);
  AddOption("rtol", "relative convergence tolerance", 1e-8);
  AddOption("stol", "convergence tolerance in terms of the norm of the change in the solution between steps", 1e-8);
  AddOption("inftol", "infinity norm convergence tolerance", 1e-8);
  AddOption("velocity-solve", "solve for velocity increments in implicit Euler (if false, the solve is for position increments)", false);
}

RodTimeStepper* Problem::getRodTimeStepper(ElasticRod& rod)
{
  RodTimeStepper* stepper = new RodTimeStepper(rod);

  string& integrator = GetStringOpt("integrator");
  if (integrator == "symplectic") {
    stepper->setDiffEqSolver(RodTimeStepper::SYMPL_EULER);

  } else if (integrator == "implicit") {
    stepper->setDiffEqSolver(RodTimeStepper::IMPL_EULER);

  } else if (integrator == "statics") {
    stepper->setDiffEqSolver(RodTimeStepper::STATICS);

  } else {
    cerr << "Unknown integrator " << integrator
         << ". Using default instead." << endl;
  }

  stepper->setTimeStep(getDt());

  Scalar massDamping = GetScalarOpt("mass-damping");
  if (massDamping != 0) {
    stepper->addExternalForce(new RodMassDamping(massDamping));
  }

  if (getGravity().norm() > 0) {
    stepper->addExternalForce(new RodGravity(getGravity()));
  }

  stepper->setMaxIterations(GetIntOpt("iterations"));
  stepper->set_stol(GetScalarOpt("stol"));
  stepper->set_atol(GetScalarOpt("atol"));
  stepper->set_rtol(GetScalarOpt("rtol"));
  stepper->set_inftol(GetScalarOpt("inftol"));

//   cout << "stol " << stepper->get_stol() << endl
//        << "atol " << stepper->get_atol() << endl
//        << "rtol " << stepper->get_rtol() << endl
//        << "inftol " << stepper->get_inftol() << endl
//        << "maxit " << stepper->getMaxIterations() << endl;

  return stepper;
}

MultipleRodTimeStepper* Problem::getMultipleRodTimeStepper()
{
  MultipleRodTimeStepper* stepper = new MultipleRodTimeStepper();

  string& integrator = GetStringOpt("integrator");

  if (integrator == "symplectic") {
    stepper->setDiffEqSolver(MultipleRodTimeStepper::SYMPL_EULER);
  } 
  if (integrator == "implicit")
  {
    stepper->setDiffEqSolver(MultipleRodTimeStepper::IMPL_EULER);
  } 
  else 
  {
    cerr << "Unknown integrator " << integrator
    << ". Using default instead." << endl;
  }
  //else if (integrator == "statics") {
  //  stepper->setDiffEqSolver(RodTimeStepper::STATICS);
  //} 
  
  stepper->setTimeStep(getDt());
  
  Scalar massDamping = GetScalarOpt("mass-damping");
  if (massDamping != 0) {
    stepper->addExternalForce(new RodMassDamping(massDamping));
  }
  
  if (getGravity().norm() > 0) {
    stepper->addExternalForce(new RodGravity(getGravity()));
  }
  
  stepper->setMaxIterations(GetIntOpt("iterations"));
  stepper->set_stol(GetScalarOpt("stol"));
  stepper->set_atol(GetScalarOpt("atol"));
  stepper->set_rtol(GetScalarOpt("rtol"));
  stepper->set_inftol(GetScalarOpt("inftol"));
  
  //   cout << "stol " << stepper->get_stol() << endl
  //        << "atol " << stepper->get_atol() << endl
  //        << "rtol " << stepper->get_rtol() << endl
  //        << "inftol " << stepper->get_inftol() << endl
  //        << "maxit " << stepper->getMaxIterations() << endl;
  
  return stepper;  
}



/*Problem::Problem()
  : m_time(0.0)
  , m_dt(0.1)
  {}
  {

  Problem::Problem()
  : _integrator(NULL)
  , _projection(NULL)
  , m_time(0.0)
  {
  AddOption("integrator", "the type of time integration", "dynamic");
  AddOption("dt", "the timestep", 1.0e-2);
  AddOption("max-adaptive-level", "maximum level for apative timestepping", 3);

  RodOpts opts;
  AddOption("nv", "number of vertices", opts.vertices);
  AddOption("kb", "bending stiffness", opts.bendStiffness);
  AddOption("kt", "twisting stiffness", opts.twistStiffness);
  AddOption("ks", "stretching stiffness", opts.stretchStiffness);
  AddOption("surface-tension", "surface tension", opts.surfaceTension);
  AddOption("generalized-stiffness", "assume inputs are physical parameters",
  opts.generalizedStiffness);
  AddOption("radius", "radius of rod", opts.radius);
  AddOption("rod-density", "density of rod", opts.density);
  AddOption("linear-density", "use linear or volumetric density",
  opts.linearDensity);
  AddOption("inextensible", "The rod is inextensible", opts.inextensible);
  AddOption("aniso", "The rod is anisotropic", opts.anisotropic);
  AddOption("closed", "The rod is a loop", opts.isClosed);
  AddOption("quasistatic", "material-frame updated quasistatically",
  opts.quasistatic);
  AddOption("viscous", "Viscous or elastic rod", opts.viscous);
  AddOption("damping", "mass damping", opts.damping);
  AddOption("inertia-scale-factor", "number to divide rotational inertia by",
  opts.inertiaScaleFactor);
  AddOption("use-time-parallel", "use time or space parallel reference frame",
  opts.useTimeParallel);

  AddOption("twist", "initial twist of the rod", 25.0);
  AddOption("gravity", "gravity", -9.81);
  AddOption("gravity-rod", "apply gravity to rod", true);
  AddOption("clampEnds", "don't know", false);

  AddOption("detect-rod-collisions", "perform rod collision detection", true);
  AddOption("detect-body-collisions", "perform body collision detection", true);

  AddOption("video-logging", "log data required to create a video", false);
  AddOption("video-log-file", "log file for video data", "video.log");
  AddOption("video-ignore-pausing", "Don't produce frames when the simulation "
  "is paused", false);
  AddOption("save-images", "turn a video log into a sequence of images",
  false);
  AddOption("save-pngs", "saves pngs in simulation time", false);

  AddOption("friction", "the surface friction of rigid bodies", 10.0);
  AddOption("body-density", "The density of the rigid bodies", 0.01);

  AddOption("render", "render the simulation", true);
  AddOption("window-width", "The width in pixels of the render window", 640);
  AddOption("window-height", "The height in pixels of the render window", 480);
  AddOption("render-rod", "render the rod nicely", true);
  AddOption("texture-width", "width of texture on rod", 1.0);
  AddOption("render-slices", "resolution of the segments", 12);
  AddOption("color1", "render color 1", Vec3d(255, 125, 75));
  AddOption("color2", "render color 2", Vec3d(75, 125, 255));
  AddOption("twist-texture", "twist texture or vertices", false);
  AddOption("use-striped-texture", "use candy cane stripes for the texture",
  true);
  AddOption("draw-arrows", "draw arrows for the frames", false);

  AddOption("read-camera", "read camera parameters from config file", false);
  AddOption("eye", "eye", Vec3d(0,0,0));
  AddOption("up", "up", Vec3d(0,1,0));
  AddOption("center", "view center", Vec3d(0,0,-1));
  AddOption("fps", "frames to capture per second", 30.0);
  AddOption("max-frames", "The maximum number of frames to output", -1);
  AddOption("start-paused", "Start simulation paused", true);

  AddOption("save-rods-as-objs", "Output the rod mesh to obj files in "
  "simulation time", false);
  AddOption("save-rods-as-text", "Outputs the state information of the rod.  "
  "The rod object can then be reconstructed later and interpolated "
  "to any vertex resolution", false);
  AddOption("save-bodies", "Outputs the bodies to files.  Trimesh bodies are "
  "output as objs and other rigid bodies are output in a textual "
  "format", false);
  AddOption("save-all-in-one-obj", "Saves all rods and bodies for a timestep "
  "in one obj", false);
  AddOption("data-directory", "The directory in which to store data",
  "render_data");
  AddOption("print-time", "Prints time after each step", false);

  AddOption("youngs-modulus", "For ropes only", 20000000);
  AddOption("quaternion-unit-constraint", "for ropes only", 0.0);
  AddOption("quaternion-parallel-constraint", "for ropes only", 100000.0);
  AddOption("stretch-damping", "for ropes only", 0.0);
  AddOption("twist-damping", "for ropes only", 0.0);
  AddOption("inertia-tensor-factor", "for ropes only", 0.01);
  AddOption("director-damping", "for ropes only", 0.0);
  AddOption("point-damping", "for ropes only", 0.0);
  AddOption("stretch-stiffness", "for ropes", 1000000);
  AddOption("quaternion-iterations", "for ropes", 10);
  }



  Problem::~Problem()
  {
  for (uint i = 0; i < _forces.size(); ++i) delete _forces[i];
  if (_integrator != NULL) delete _integrator;
  }

*/

/*
  void Problem::BaseAfterLoad() {
  m_dt = GetRealOpt("dt");
  AfterLoad();
  // Set defaults from the options file
  RigidBody::setFriction(GetRealOpt("friction"));
  RigidBody::setDensity(GetRealOpt("body-density"));
  RigidBody::setGravity(GetRealOpt("gravity"));
  RigidBody::setCollide(GetBoolOpt("detect-body-collisions"));
  }

  void Problem::BaseSetup() {
  Real g = GetRealOpt("gravity");
  if (g != 0.0 && GetBoolOpt("gravity-rod")) {
  _forces.push_back(new Gravity(0, g, 0));
  }

  Setup();

  string type = GetStringOpt("integrator");

  if (type == "static") {
  _integrator = new Static(_rods, _bodies, _forces);
  } else if (type == "viscous") {
  _integrator = new ViscousDynamics(_rods, _bodies, _forces);
  } else if (type == "dynamic") {
  _integrator = new Dynamic(_rods, _bodies, _forces);
  } else if (type == "smart") {
  _integrator = new SmartIntegrator(_rods, _bodies, _forces);
  } else if (type == "implicit") {
  _integrator = new ImplicitIntegrator(_rods, _bodies, _forces);
  } else {
  cerr << "Unknown integrator type \'" << type << "\'" << endl;
  assert(0);
  exit(1);
  }

  if (_constraints.size() > 0) {
  _projection = new FastProjection(_objs, _constraints);
  _projection->project();
  }
  }

  bool Problem::AddOption(const OptionKey& name, const string& desc,
  const bool& def)
  {
  if(m_options.find(name) != m_options.end()) return false;

  Option opt(name, desc, def);
  m_options.insert(make_pair(name, opt));
  return true;
  }


  bool Problem::AddOption(const OptionKey& name, const string& desc,
  const int& def)
  {
  if(m_options.find(name) != m_options.end()) return false;

  Option opt(name, desc, def);
  m_options.insert(make_pair(name, opt));
  return true;
  }

  bool Problem::AddOption(const OptionKey& name, const string& desc,
  const uint& def)
  {
  int iDef = def;
  return AddOption(name, desc, iDef);
  }

  bool Problem::AddOption(const OptionKey& name, const string& desc,
  const Real& def)
  {
  if(m_options.find(name) != m_options.end()) return false;

  Option opt(name, desc, def);
  m_options.insert(make_pair(name, opt));
  return true;
  }

  bool Problem::AddOption(const OptionKey& name, const string& desc,
  const Vec3d& def)
  {
  if(m_options.find(name) != m_options.end()) return false;

  Option opt(name, desc, def);
  m_options.insert(make_pair(name, opt));
  return true;
  }

  bool Problem::AddOption(const OptionKey& name, const string& desc,
  const string& def)
  {
  if(m_options.find(name) != m_options.end()) return false;

  Option opt(name, desc, def);
  m_options.insert(make_pair(name, opt));
  return true;
  }

  // added by basile
  // for options defined by C strings, we explicitly cast to char, as implicit cast yields wrong answer (i.e. converts def to bool)
  bool Problem::AddOption(const OptionKey& name, const string& desc,
  const char* def)
  {
  // in case of ambiguity, issue a warning
  if(  !strcmp(def, "true") || !strcmp(def, "1")  || !strcmp(def, "false") || !strcmp(def, "0") ) {
  cout << " ************ WARNING !! ************************************" << endl;
  cout << " ************ ambiguous option initialization found in AddOption(..., ..., \"" << def << "\");" << endl;
  cout << " ************ option assumed to be of type string" << endl;
  cout << " ************ if bool is intended instead, make last argument of AddOption a bool, not a string " << endl;
  cout << " ************************************************************" << endl;
  }
  string strdef(def);
  return AddOption(name, desc, strdef);
  }

  void Problem::AddRod(Rod* rod)
  {
  assert(rod != NULL);

  rod->setDamping(GetRealOpt("damping"));

  rod->setup();

  if (rod->ks() == 0.0 && rod->inextensible()) {
  Constraints* con = new Inextensibility(*rod);
  AddConstraint(con);
  }

  RodController* controller = new RodController(*rod);
  string type = GetStringOpt("integrator");
  if (type == "dynamic" || type == "smart" || type == "explicit") {
  controller->setIntegrator(RodController::EXPLICIT);
  } else if (type == "implicit") {
  controller->setIntegrator(RodController::IMPLICIT);
  }
  for (uint i = 0; i < _forces.size(); ++i) {
  controller->addExternalForce(_forces[i]);
  }
  rod->setController(controller);

  _rods.push_back(rod);
  _objs.push_back(rod);

  if (GetBoolOpt("detect-rod-collisions")) {
  simulation::ArrayDistribution<Real>* masses
  = new simulation::ArrayDistribution<Real>(*rod);
  _collisions.AddRod(rod, *masses, rod->nv(), _rods.size()-1);
  }
  }



  void Problem::AddBody(RigidBody* body)
  {
  assert(body != NULL);
  _bodies.push_back(body);
  _objs.push_back(body);
  }



  void Problem::AddConstraint(Constraints* con)
  {
  static uint cid = 0;
  IdxArray cidx;
  for (uint i = 0; i < con->nConstraints(); i++) cidx.push_back(cid++);
  con->setCidx(cidx);
  _constraints.push_back(con);
  }

  Option* Problem::GetOption(OptionKey s)
  {
  map<string, Option>::iterator itr = m_options.find(s);
  if (itr != m_options.end()) {
  return &(itr->second);
  } else {
  cerr << "No such option: " << s << endl;
  exit(1);
  }
  }

  int& Problem::GetIntOpt(OptionKey s) {
  return GetOption(s)->i;
  }

  Real& Problem::GetRealOpt(OptionKey s) {
  return GetOption(s)->r;
  }

  bool& Problem::GetBoolOpt(OptionKey s) {
  return GetOption(s)->b;
  }

  Vec3d& Problem::GetVecOpt(OptionKey s) {
  return GetOption(s)->v;
  }

  string& Problem::GetStringOpt(OptionKey s) {
  return GetOption(s)->s;
  }

  void Problem::BaseAtEachTimestep()
  {
  START_TIMER("TOTAL");

  static bool collisions = GetBoolOpt("detect-rod-collisions");
  static int max_adaptive_level = GetIntOpt("max-adaptive-level");

  bool accept = false;
  bool cutTimeStep = false;
  static int level = 0; // for adaptive time stepping;

  for (uint i = 0; i < _rods.size(); ++i) {
  _rods[i]->getController()->prepareStep();
  }

  Real beginTime = m_time;
  Real divisor = 1.;
  while (!accept) {
  m_time = beginTime;
  accept = true;
  divisor = pow(2.0, level);
  Real dt = m_dt / divisor;
  AtEachTimestep();

  START_TIMER("CONSTRAINTS");
  // transfers torque from twisting of rod to rigid bodies (if any)
  for (uint i = 0; i < _constraints.size(); i++)
  _constraints[i]->beforeStep();
  STOP_TIMER("CONSTRAINTS");

  START_TIMER("STEP");
  // steps the simulation forward in time
  for (uint i = 0; i < _rods.size(); ++i) {
  if (_rods[i]->getController()->step(dt) < 0) {
  accept = false;
  for (uint j = 0; j <= i; ++j)
  _rods[j]->getController()->rollBackStep(dt);
  break;
  }
  }
  STOP_TIMER("STEP");

  START_TIMER("CONSTRAINTS");
  // fixes positions so constraints are satisfied
  if (_projection != NULL && accept) {
  //_projection->project();
  }
  STOP_TIMER("CONSTRAINTS");*/
/*
  START_TIMER("STEP");
  // fixes velocities to match new positions
  if (_integrator != NULL) _integrator->cleanup(dt);
  START_TIMER("STEP");
*/
/*START_TIMER("COLLISIONS");
// fixes rod self-collisions
// TODO: rod/rigid body collision detection
if (collisions && accept) {
for (uint i = 0; i < _rods.size(); ++i) {
_collisions.ChangeRodArray(i, _rods[i]->nv());
}
_collisions.DetectAndRespond(dt);
}

STOP_TIMER("COLLISIONS");

START_TIMER("CONSTRAINTS");
if (accept) {
// matches material frame on end of rod to orientation of rigid body
for (uint i = 0; i < _constraints.size(); i++) {
//_constraints[i]->cleanup();
}
STOP_TIMER("CONSTRAINTS");

m_time += dt;
AfterStep();

START_TIMER("QUASISTATIC");
// update rod quantities
for (uint i = 0; i < _rods.size(); ++i) {
//_rods[i]->updateQuantities();
_rods[i]->updateQuantities();
}
STOP_TIMER("QUASISTATIC");
}
*/
/*
  if (level < max_adaptive_level) {
  for (uint i = 0; i < _rods.size(); ++i) {
  //_rods[i]->computeInternalForce(afterForce[i], afterTorque[i]);
  Real max = 0.0; int idx = -1;
  for (uint j = 0; j < afterForce[i].size(); ++j) {
  if (afterForce[i][j].norm() > max) {
  max = afterForce[i][j].norm();
  idx = (int) j;
  }
  }
  if (idx >= 0 && idx < (int)beforeForce.size()
  && beforeForce[i][idx].norm() > 0.0)
  {
  Vec3d diff = afterForce[i][idx] - beforeForce[i][idx];
  Real ratio = diff.norm() / beforeForce[i][idx].norm();
  if (ratio >= 1.1) {
  accept = false;
  break;
  }
  }
  }
  } else accept = true;
*/
/*
  for (uint i = 0; i < _rods.size(); ++i) {
  if (_rods[i]->getController()->checkStableStep(dt) < 0) {
  accept = false;
  for (uint j = 0; j < _rods.size(); ++j) {
  _rods[j]->getController()->rollBackStep(dt);
  }
  break;
  }
  }

  if (!accept) {
  cutTimeStep = true;
  ++level;
  //for (uint i = 0; i < _rods.size(); ++i) {
  //_rods[i]->restoreSavedState();
  //}
  } else {
  for (uint i = 0; i < _rods.size(); ++i) {
  _rods[i]->getController()->finalizeStep(dt);
  }
  }
  }

  START_TIMER("TIMER");
  STOP_TIMER("TIMER");

  if (cutTimeStep) {
  cout << "time step divisor increased to " << divisor << " (2^"
  << level << ")" << endl;
  } else if (level > 0) {
  --level;
  divisor = pow(2.0, level);
  cout << "time step divisor decreased to " << divisor << " (2^"
  << level << ")" << endl;
  }

  STOP_TIMER("TOTAL");
  }



  void Problem::GetRodOpts(RodOpts& opts)
  {
  opts.vertices = GetIntOpt("nv");
  opts.bendStiffness = GetRealOpt("kb");
  opts.twistStiffness = GetRealOpt("kt");
  opts.stretchStiffness = GetRealOpt("ks");
  opts.surfaceTension = GetRealOpt("surface-tension");
  opts.generalizedStiffness = GetBoolOpt("generalized-stiffness");
  opts.radius = GetRealOpt("radius");
  opts.density = GetRealOpt("rod-density");
  opts.linearDensity = GetBoolOpt("linear-density");
  opts.inextensible = GetBoolOpt("inextensible");
  opts.anisotropic = GetBoolOpt("aniso");
  opts.isClosed = GetRealOpt("closed");
  opts.quasistatic = GetBoolOpt("quasistatic");
  opts.viscous = GetBoolOpt("viscous");
  opts.damping = GetRealOpt("damping");
  opts.inertiaScaleFactor = GetRealOpt("inertia-scale-factor");
  opts.useTimeParallel = GetBoolOpt("use-time-parallel");
  }



  Rod* Problem::CreateRod()
  {
  RodOpts opts;
  GetRodOpts(opts);
  return new Rod(opts);
  }



  Rod* Problem::CreateHelicalRodWithNiceEnds(const Vec3d& start_pos,
  const Vec3d& direction,
  Real radius,
  Real height,
  Real windings,
  bool wrap_ccw)
  {
  Rod* rod = CreateRod();
  int ne = rod->ne();

  Vec3d dir = direction;
  dir.normalize();

  // Get a vector orthogonal to the direction vector
  Vec3d radius_vector;
  Utils::findOrthogonal(dir, radius_vector);
  radius_vector *= radius;

  // Some portion of the height goes to making the nice ends. let's
  // say it's half a winding in height on each end plus an additional
  // quarter winding of straight part
  //int ne_straight = ceil(0.15*ne);
  int ne_straight = 1;
  Real edges_per_winding = (ne - 2*ne_straight) / windings;
  int ne_curve = (int) (0.5*(edges_per_winding+1));
  int ne_twist = ne - 2*ne_curve - 2*ne_straight;

  Vec3d center_line = start_pos;
  Vec3d dir_increment = (height / ne) * dir;

  rod->vertex(0) = rod->vertex0(0) = center_line;
  int i = 0; // current edge

  // Draw the straight part at the start
  for (i = 0; i < ne_straight; ++i) {
  center_line += dir_increment;
  rod->vertex(i+1) = rod->vertex0(i+1) = center_line;
  }

  // The theta for the small semicircles on the ends
  Real theta1 = M_PI / ne_curve;
  if (!wrap_ccw) {
  theta1 = -theta1;
  }

  // Get the rotation matrix for the small semicircles
  Vec3d center_offset = radius_vector / 2;
  Vec3d mini_radius = -1.0 * center_offset;

  Real b = 0.7;
  Real a = 0.5 * ne_curve;
  a /= (1.0 / b) * exp(ne_curve * b) - ne_curve - 1.0 / b;

  // Draw the end semicircle
  int i_curr = i;
  for (; i < ne_curve + i_curr; ++i) {
  Utils::rotateAxisAngle(mini_radius, dir, theta1);
  Real len = 0.5 - a + a * exp((ne_curve - i + i_curr) * b);
  center_line += len * dir_increment;
  //center_line += dir_increment;
  rod->vertex(i+1) = rod->vertex0(i+1)
  = center_line + center_offset + mini_radius;*/
/*
  positions[3*i] = undeformed[3*i] = v[0];
  positions[3*i+1] = undeformed[3*i+1] = v[1];
  positions[3*i+2] = undeformed[3*i+2] = v[2];
*/
/*}

// The angle around the direction axis between sequential vertices
Real theta2 = windings * 2 * M_PI / ne_twist;
if (!wrap_ccw) {
theta2 = -theta2;
}

i_curr = i;
for (; i < ne_twist + i_curr; i++) {
center_line += dir_increment;
Utils::rotateAxisAngle(radius_vector, dir, theta2);
rod->vertex(i+1) = rod->vertex0(i+1) = center_line + radius_vector;
}

center_offset = radius_vector / 2;
mini_radius = center_offset;

i_curr = i;
for (; i < ne_curve + i_curr; ++i) {
Real len = 0.5 - a + a * exp((i - i_curr) * b);
center_line += len * dir_increment;
//center_line += dir_increment;
Utils::rotateAxisAngle(mini_radius, dir, theta1);
rod->vertex(i+1) = rod->vertex0(i+1) = center_line + center_offset + mini_radius;
}

i_curr = i;
for (; i < ne_straight + i_curr; ++i) {
center_line += dir_increment;
rod->vertex(i+1) = rod->vertex0(i+1) = center_line;
}

return rod;

}

Rod* Problem::CreateHelicalRod(const Vec3d& start_pos,
const Vec3d& direction,
Real radius,
Real height,
Real windings,
bool wrap_ccw) {
Rod* rod = CreateRod();
int nv = rod->nv();

Vec3d dir = direction;
dir.normalize();

// Get a vector orthogonal to the direction vector
Vec3d radius_vector;
Utils::findOrthogonal(dir, radius_vector);
radius_vector *= radius;

// The angle around the direction axis between sequential vertices
Real theta = windings * 2 * M_PI / nv;
if (!wrap_ccw) {
theta = -theta;
}
Mat3d m;
Utils::matrixAxisAngle(m, dir, theta);

Vec3d center_line = start_pos;
Vec3d dir_increment = (height / nv) * dir;
for (int i = 0; i < nv; i++) {
Vec3d v = radius_vector + center_line;
rod->vertex(i) = rod->vertex0(i) = v;

radius_vector = m*radius_vector;
center_line += dir_increment;
}

return rod;
}

Rod* Problem::CreateStraightRod(const Vec3d& start_pos,
const Vec3d& end_pos) {
Rod* rod = CreateRod();
int nv = rod->nv();

Vec3d dir = end_pos - start_pos;
dir = dir / (nv - 1);
Vec3d pos = start_pos;
for (int i = 0; i < nv; ++i) {
rod->vertex(i) = rod->vertex0(i) = pos;
pos = pos + dir;
}
return rod;
}



Rod* Problem::CreateCircularRod(const Vec3d& center,
const Vec3d& normal,
Real radius)
{
Rod* rod = CreateRod();
int nv = rod->nv();

Vec3d axis1; Utils::findOrthogonal(normal, axis1);
Vec3d axis2 = axis1.cross(normal);
axis2.normalize();

for (int i = 0; i < nv; ++i) {
Real theta = 2.0*M_PI * i / rod->ne();
Vec3d v = radius * cos(theta) * axis1 + radius * sin(theta) * axis2;
rod->vertex(i) = rod->vertex0(i) = v;
}

return rod;
}



RigidBody* Problem::AddBox(const Vec3d& len,
const Vec3d& pos,
const Quaternion& orn,
const Vec3d& lvel,
const Vec3d& avel) {
RigidBody::Descriptor desc;
desc.type = RigidBody::BOX;
desc.params[0] = len[0];
desc.params[1] = len[1];
desc.params[2] = len[2];

RigidBody* body = new RigidBody(desc);
body->setPosition(pos);
body->setOrientation(orn);
body->setVelocity(lvel);
body->setAngularVelocity(avel);
AddBody(body);
return body;
}

RigidBody* Problem::AddBox(const Vec3d& len,
const Vec3d& pos,
const Quaternion& orn,
bool fixed) {
RigidBody::Descriptor desc;
desc.type = RigidBody::BOX;
desc.params[0] = len[0];
desc.params[1] = len[1];
desc.params[2] = len[2];
desc.fixed = fixed;

RigidBody* body = new RigidBody(desc);
body->setPosition(pos);
body->setOrientation(orn);
AddBody(body);
return body;
}

RigidBody* Problem::AddPlane(Real a, Real b, Real c, Real d) {
RigidBody::Descriptor desc;
desc.type = RigidBody::PLANE;
desc.params[0] = a;
desc.params[1] = b;
desc.params[2] = c;
desc.params[3] = d;
desc.fixed = true;

RigidBody* body = new RigidBody(desc);
AddBody(body);
return body;
}

RigidBody* Problem::AddPlane(const Real* params) {
RigidBody::Descriptor desc;
desc.type = RigidBody::PLANE;
desc.params[0] = params[0];
desc.params[1] = params[1];
desc.params[2] = params[2];
desc.params[3] = params[3];
desc.fixed = true;

RigidBody* body = new RigidBody(desc);
AddBody(body);
return body;
}

RigidBody* Problem::AddSphere(Real radius,
const Vec3d& pos,
const Quaternion& orn,
const Vec3d& lvel,
const Vec3d& avel) {
RigidBody::Descriptor desc;
desc.type = RigidBody::SPHERE;
desc.params[0] = radius;

RigidBody* body = new RigidBody(desc);
body->setPosition(pos);
body->setOrientation(orn);
body->setVelocity(lvel);
body->setAngularVelocity(avel);
AddBody(body);
return body;
}

RigidBody* Problem::AddSphere(Real radius,
const Vec3d& pos,
const Quaternion& orn,
bool fixed) {
RigidBody::Descriptor desc;
desc.type = RigidBody::SPHERE;
desc.params[0] = radius;
desc.fixed = fixed;

RigidBody* body = new RigidBody(desc);
body->setPosition(pos);
body->setOrientation(orn);
AddBody(body);
return body;
}

RigidBody* Problem::AddTrimesh(const string& file,
const Vec3d& pos,
const Quaternion& orn,
const Vec3d& lvel,
const Vec3d& avel) {
RigidBody::Descriptor desc;
desc.type = RigidBody::TRIMESH;
desc.mesh_file = file;

RigidBody* body = new RigidBody(desc);
body->setPosition(pos);
body->setOrientation(orn);
body->setVelocity(lvel);
body->setAngularVelocity(avel);
AddBody(body);
return body;
}

RigidBody* Problem::AddTrimesh(const string& file,
const Vec3d& pos,
const Quaternion& orn,
bool fixed) {
RigidBody::Descriptor desc;
desc.type = RigidBody::TRIMESH;
desc.mesh_file = file;
desc.fixed = fixed;

RigidBody* body = new RigidBody(desc);
body->setPosition(pos);
body->setOrientation(orn);
AddBody(body);
return body;
}

RigidBody* Problem::AddCylinder(Real radius,
Real height,
const Vec3d& pos,
const Quaternion& orn,
const Vec3d& lvel,
const Vec3d& avel) {
RigidBody::Descriptor desc;
desc.type = RigidBody::CYLINDER;
desc.params[0] = radius;
desc.params[1] = height;

RigidBody* body = new RigidBody(desc);
body->setPosition(pos);
body->setOrientation(orn);
body->setVelocity(lvel);
body->setAngularVelocity(avel);
AddBody(body);
return body;
}

RigidBody* Problem::AddCylinder(Real radius,
Real height,
const Vec3d& pos,
const Quaternion& orn,
bool fixed) {
RigidBody::Descriptor desc;
desc.type = RigidBody::CYLINDER;
desc.params[0] = radius;
desc.params[1] = height;
desc.fixed = fixed;

RigidBody* body = new RigidBody(desc);
body->setPosition(pos);
body->setOrientation(orn);
AddBody(body);
return body;
}
*/
