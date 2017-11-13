/**
 * \file ViscousRod.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 08/29/2009
 */

#ifndef VISCOUSROD_HH
#define VISCOUSROD_HH

#include "BASim/src/Physics/ElasticRods/ElasticRod.hh"

namespace BASim {

/** Class for simulating viscous rods. */
class ViscousRod : public ElasticRod
{
public:

  ViscousRod(int numVertices = 3);

  void setup();

  void setViscosity(Scalar viscosity);

  void updateReferenceProperties();

protected:

  void computeEdgeVolumes();

  EPropHandle<Scalar> m_edgeVolume;
};

} // namespace BASim

#endif // VISCOUSROD_HH
