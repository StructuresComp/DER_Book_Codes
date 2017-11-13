/**
 * \file RodUtils.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/03/2009
 */

#ifndef RODUTILS_HH
#define RODUTILS_HH

#include "BASim/src/Physics/ElasticRods/AnisotropicRod.hh"

namespace BASim {

struct RodOptions {
  int numVertices;
  Scalar density;
  Scalar radiusA;
  Scalar radiusB;
  Scalar radiusScale;
  Scalar YoungsModulus;
  Scalar ShearModulus;
  Scalar viscosity;
  Scalar radiusInner;  // for hollow tube
  bool anisotropic;
  bool elastic;
  bool quasistatic;
  bool inextensible;
  bool hollow;
  ElasticRod::RefFrameType refFrame;

  RodOptions()
    : numVertices(50)
    , density(1.3)
    , radiusA(5e-3)
    , radiusB(5e-3)
    , radiusScale(1.0)
    , YoungsModulus(1e10)
    , ShearModulus(3.4e9)
    , viscosity(50)
    , radiusInner(0.0)
    , anisotropic(true)
    , elastic(true)
    , quasistatic(true)
    , inextensible(false)
    , hollow(false)
    , refFrame(ElasticRod::TimeParallel)
  {}
};

inline ElasticRod* setupRod(const RodOptions& opts,
                            const std::vector<Vec3d>& initialPosition,
                            const std::vector<Vec3d>& undeformedPosition)
{
  assert(opts.numVertices == (int) initialPosition.size());
  assert(opts.numVertices == (int) undeformedPosition.size());

  ElasticRod* rod = NULL;
  if (opts.anisotropic) rod = new AnisotropicRod(opts.numVertices);
  else rod = new AnisotropicRod(opts.numVertices);

  rod->setRadius(opts.radiusA, opts.radiusB);
  rod->setRadiusScale(opts.radiusScale);
  rod->setDensity(opts.density);
  rod->setYoungsModulus(opts.YoungsModulus);
  rod->setShearModulus(opts.ShearModulus);
  rod->setViscosity(opts.viscosity);
  rod->setQuasistatic(opts.quasistatic);
  rod->setRefFrameType(opts.refFrame);
  
  rod->setHollow(opts.hollow, opts.radiusInner);

  // set up using undeformed positions
  for (int i = 0; i < rod->nv(); ++i)
    rod->setVertex(i, undeformedPosition[i]);
  rod->setup();

  // update to initial positions
  for (int i = 0; i < rod->nv(); ++i)
    rod->setVertex(i, initialPosition[i]);
  rod->updateProperties();

  return rod;
}

inline Vec3d calculateObjectCenter(const ElasticRod& rod)
{
  Vec3d center = Vec3d::Zero();

  ElasticRod::vertex_iter vit, end = rod.vertices_end();
  for (vit = rod.vertices_begin(); vit != end; ++vit) {
    center += rod.getVertex(*vit);
  }

  center /= rod.nv();

  return center;
}

inline Scalar calculateObjectBoundingRadius(const ElasticRod& rod,
                                            const Vec3d& center)
{
  Scalar radius = 0.0;

  ElasticRod::vertex_iter vit, end = rod.vertices_end();
  for (vit = rod.vertices_begin(); vit != end; ++vit) {
    radius = std::max(radius, (rod.getVertex(*vit) - center).norm());
  }

  return radius;
}

} // namespace BASim

#endif // RODUTILS_HH
