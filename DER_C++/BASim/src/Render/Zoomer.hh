#ifndef ZOOMER_H
#define ZOOMER_H

namespace BASim {

/** Converts mouse motions into world translations. Modifies the
    camera's eye and view center directly. */
class Zoomer
{
public:

  /// Default constructor
  Zoomer(Camera* c, const Scalar scale = 1.0);

  /// Set a particular camera to use.
  void setCamera(Camera* c);

  /// Set the scaling for mouse motions.
  void setScale(const Scalar s);

  /// Start a mouse motion.
  /// Position in [-1,1] x [-1,1].
  void start(const Vec2d& p);

  /// Update a mouse motion with a new mouse point.
  /// Position in [-1,1] x [-1,1].
  void update(const Vec2d& p);

  /// Stop a mouse motion.
  void stop();

protected:

  Camera* m_camera;                       ///< The current camera to obey.
  bool m_translating;                     ///< Whether we are translating.
  Vec2d m_startPos;                       ///< Start of a translation drag.
  Scalar m_scale;                         ///< Scaling of mouse motions.
};

} // namespace BASim

#endif // ZOOMER_H
