/**
 * \file AnisotropicRod.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 08/29/2009
 */

#ifndef ANISOTROPICROD_HH
#define ANISOTROPICROD_HH

#include "BASim/src/Physics/ElasticRods/ElasticRod.hh"

namespace BASim {

/** Class for simulating anisotropic rods. */
class AnisotropicRod : public ElasticRod
{
public:

  AnisotropicRod(int numVertices = 3);

  void setup();
};

} // namespace BASim

#endif // ANISOTROPICROD_HH
