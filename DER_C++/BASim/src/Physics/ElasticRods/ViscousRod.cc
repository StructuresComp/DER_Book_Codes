/**
 * \file ViscousRod.cc
 *
 * \author miklos@cs.columbia.edu
 * \date 08/29/2009
 */

#include "BASim/src/Physics/ElasticRods/ViscousRod.hh"
#include "BASim/src/Physics/ElasticRods/RodStretchingForce.hh"
#include "BASim/src/Physics/ElasticRods/RodBendingForce.hh"
#include "BASim/src/Physics/ElasticRods/RodTwistingForce.hh"
#include "BASim/src/Physics/ElasticRods/RodAnisoForce.hh"

namespace BASim {

ViscousRod::ViscousRod(int numVertices)
  : ElasticRod(numVertices)
{
  setRefFrameType(TimeParallel);
}

void ViscousRod::setup()
{
  ElasticRod::setup();
  computeEdgeVolumes();

  addForce(new RodStretchingForce(*this));
  addForce(new RodTwistingForce(*this));
  if (refFrameType() == TimeParallel) addForce(new RodBendingForce(*this));
  else addForce(new RodAnisoForce(*this));
}

void ViscousRod::computeEdgeVolumes()
{
  edge_iter eit, end = edges_end();
  for (eit = edges_begin(); eit != end; ++eit) {
    edge_handle& eh = *eit;
    Scalar a = radiusA(eh);
    Scalar b = radiusB(eh);
    Scalar h = getEdgeLength(eh);
    PhysObject::property(m_edgeVolume)[eh] = computeMass(1.0, a, b, h);
  }
}

void ViscousRod::setViscosity(Scalar viscosity)
{
  setYoungsModulus(3.0 * viscosity);
  setShearModulus(2.0 * viscosity);
}

void ViscousRod::updateReferenceProperties()
{
  // preservation of volume for edges means radius must change
  // according to change in length of edge
  edge_iter eit, end = edges_end();
  for (eit = edges_begin(); eit != end; ++eit) {
    edge_handle& eh = *eit;
    const Scalar& v = property(m_edgeVolume)[eh];
    const Scalar& h = getEdgeLength(eh);
    setRadius(eh, sqrt(v / (M_PI * h)));
  }

  ElasticRod::updateReferenceProperties();
}

} // namespace BASim
