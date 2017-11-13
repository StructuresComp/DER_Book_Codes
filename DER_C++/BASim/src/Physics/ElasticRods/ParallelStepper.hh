/**
 * \file ParallelStepper.hh
 *
 * \author smith@cs.columbia.edu
 * \date 05/23/2010
 */


#ifndef PARALLELSTEPPER_HH
#define PARALLELSTEPPER_HH

#include "BASim/src/Core/ObjectControllerBase.hh"
#include "BASim/src/Physics/ElasticRods/ElasticRod.hh"
#include "BASim/src/Physics/ElasticRods/RodTimeStepper.hh"
#include "BASim/src/Physics/ElasticRods/RodMassDamping.hh"
#include "BASim/src/Physics/ElasticRods/RodGravity.hh"
#include "BASim/src/Math/Math.hh"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include <limits>

namespace BASim {

  /**
   * Class to execute a number of controllers in parallel. This assumes
   * that the controllers can safely execute in parallel. 
   */
  class ParallelStepper : public ObjectControllerBase
  {

  public:
    /**
     * Default constructor.
     */
    ParallelStepper();

    /**
     * Destructor.
     */
    virtual ~ParallelStepper();

    /**
     * Adds a rod that will be evolved in time using this BridsonStepper.
     */
    void addController( ObjectControllerBase* controller );

    /**
     * Executes all inserted controllers in parallel. 
     */
    void execute();
    
  private:
    std::vector<ObjectControllerBase*> m_controllers;

  };
  
} // namespace BASim

#endif // PARALLELSTEPPER_HH

