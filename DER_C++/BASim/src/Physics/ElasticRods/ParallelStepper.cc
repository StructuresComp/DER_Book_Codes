/**
 * \file ParallelStepper.cc
 *
 * \author smith@cs.columbia.edu
 * \date 05/23/2010
 */

#include "ParallelStepper.hh"


namespace BASim 
{
  
ParallelStepper::ParallelStepper()
: m_controllers()
{}

ParallelStepper::~ParallelStepper()
{}
  
  
void ParallelStepper::addController( ObjectControllerBase* controller )
{
  assert( controller != NULL );
  m_controllers.push_back( controller );
}

void ParallelStepper::execute()
{
  #ifdef HAVE_OPENMP
  #pragma omp parallel for
  #endif
  for( int i = 0; i < (int) m_controllers.size(); ++i ) 
  {
    assert( m_controllers[i] != NULL );
    m_controllers[i]->execute();
  }
}


}
