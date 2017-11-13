#ifndef TRANSLATOR_HH
#define TRANSLATOR_HH

#include "BASim/src/Render/Camera.hh"

namespace BASim {

/** Converts mouse motions into world translations. Modifies the
    camera's eye and view center directly. */
class Translator {
public:

  /// Default constructor
  Translator(Camera* c, const Scalar scale = 1.0);

  /// Set a particular camera to use.
  void setCamera(Camera* c);

  /// Set the scaling for mouse motions.
  void setScale(const Scalar s);

  /// Start a translation mouse motion.
  /// Position in [-1,1] x [-1,1].
  void start(const Vec2d& p);

  /// Update a translation mouse motion with a new mouse point.
  /// Position in [-1,1] x [-1,1].
  void update(const Vec2d& p);

  /// Stop a translation mouse motion.
  void stop();

  /// Get the current translation.
  const Vec3d& getTranslation() const;

private:
  Camera* m_camera;             ///< The current camera to obey.
  bool m_translating;           ///< Whether we are translating.
  Vec2d m_startPos;             ///< Start of a translation drag.
  Vec3d m_translation;          ///< Current user translation.
  Scalar m_scale;               ///< Scaling of mouse motions to translations.
};

} // namespace BASim

#endif // TRANSLATOR_HH
