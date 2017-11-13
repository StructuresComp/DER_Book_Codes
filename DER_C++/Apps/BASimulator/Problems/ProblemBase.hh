/**
 * \file ProblemBase.hh
 *
 * \author miklos@cs.columbia.edu (based on problem-setup.h from sar2120@columbia.edu)
 * \date 09/09/2009
 */

#ifndef PROBLEMBASE_HH
#define PROBLEMBASE_HH

#include "Option.hh"

#include <map>
#include <string>
#include <BASim/BASim>

#include "BASim/src/Physics/ElasticRods/MultipleRodTimeStepper.hh"

class Problem
{
public:

  typedef std::map<std::string, Option> OptionMap;

  Problem(const std::string& name = "", const std::string& desc = "");
  virtual ~Problem();

  void BaseSetup(int argc, char** argv);
  void BaseFinalize();
  void BaseAtEachTimestep();
  void BaseAfterLoad();
  void BaseRender();
  void BaseHUDRender(int w, int h);
  
  const std::string& ProblemName() const { return m_problemName; }
  const std::string& ProblemDescription() const { return m_problemDesc; }

  template <typename T>
  int AddOption(const std::string& name, const std::string& desc, const T& def);

  Option* GetOption(const std::string& name);
  bool& GetBoolOpt(const std::string& name);
  int& GetIntOpt(const std::string& name);
  Scalar& GetScalarOpt(const std::string& name);
  Vec3d& GetVecOpt(const std::string& name);
  std::string& GetStringOpt(const std::string& name);

  int LoadOptions(const char* filename);
  int LoadOptions(const std::string& filename)
  {
    return LoadOptions(filename.c_str());
  }
  int LoadOptions(int argc, char** argv);

  void PrintOptions(std::ostream& os);

  const World& getWorld() const { return *m_world; }
  World& getWorld() { return *m_world; }

  /**
   * Adds properties to the World that are generally used in dynamics
   * problems.
   */
  void addDynamicsProps();

  /**
   * Loads in dynamics properties from the options. setupDynamics must
   * be called before this can be used.
   */
  void loadDynamicsProps();

  bool dynamicsPropsLoaded() const { return m_dynamicsProps; }

  const Scalar& getTime() const;
  void setTime(const Scalar& time);

  const Scalar& getDt() const;
  void setDt(const Scalar& dt);

  const Vec3d& getGravity() const;
  void setGravity(const Vec3d& gravity);

  void addRodOptions();
  void getRodOptions(RodOptions& opts);
  void addRodTimeStepperOptions();
  RodTimeStepper* getRodTimeStepper(ElasticRod& rod);
  MultipleRodTimeStepper* getMultipleRodTimeStepper();

	double msg_val;
	std::string screenshot_filename;
	
	bool update_display;
	
protected:

  virtual void Setup() {}
  virtual void AtEachTimestep() {}
  virtual void AfterLoad() {}
  virtual void AfterStep() {}
  virtual void Render() {}
  virtual void HUDRender(int w, int h) {}
  
  std::string m_problemName;
  std::string m_problemDesc;
  OptionMap m_options;

  World* m_world;

  bool m_dynamicsProps;
  ObjPropHandle<Scalar> m_time;
  ObjPropHandle<Scalar> m_dt;
  ObjPropHandle<Vec3d> m_gravity;
};

#include "ProblemBase.tcc"

#endif // PROBLEMBASE_HH
