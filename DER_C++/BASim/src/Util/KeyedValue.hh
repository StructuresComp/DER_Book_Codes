/**
 * \file KeyedValue.hh
 *
 * \author Nicolas Clauvelin, Basile Audoly, and miklos@cs.columbia.edu
 * \date 09/16/2009
 */

#ifndef KEYEDVALUE_HH
#define KEYEDVALUE_HH

namespace BASim {

/** A simple keyframe class for specifying a certain value at a
    certain time. */
template <class T>
class KeyFrame
{
public:

  KeyFrame()
    : m_time(0.0)
    , m_value(T())
  {}

  KeyFrame(const Scalar& time, const T& value)
    : m_time(time)
    , m_value(value)
  {}

  ~KeyFrame() {}

  const Scalar& getTime() const
  {
    return m_time;
  }

  const T& getValue() const
  {
    return m_value;
  }

  friend bool operator<(const KeyFrame<T>& frame1, const KeyFrame<T>& frame2)
  {
    return (frame1.m_time < frame2.m_time);
  }

  friend bool operator>(const KeyFrame& frame1, const KeyFrame& frame2)
  {
    return (frame1.m_time > frame2.m_time);
  }

  friend std::ostream& operator<< (std::ostream& os, const KeyFrame& frame)
  {
    os << "{ " << frame.m_time << ", " << frame.m_value << " }";
    return os;
  }

private:
  Scalar m_time;
  T m_value;
};

/** A class for specifying a value that changes over time using
    keyframes. */
template <class T>
class KeyedValue
{
public:

  typedef std::vector< KeyFrame<T> > KeyFrameArray;

  KeyedValue() {}

  KeyedValue(const std::string& specs)
  {
    setKeyFramesFromString(specs);
  }

  ~KeyedValue() {}

  size_t getNumFrames() const
  {
    return m_keyFrames.size();
  }

  const Scalar& getFrameTime(int i) const
  {
    assert(i >= 0 && i < m_keyFrames.size());
    return m_keyFrames[i].getTime();
  }

  const T& getFrameValue(int i) const
  {
    assert(i >= 0 && i < m_keyFrames.size());
    return m_keyFrames[i].getValue();
  }

  void addKeyFrame(const KeyFrame<T>& key)
  {
    m_keyFrames.push_back(key);
    std::sort(m_keyFrames.begin(), m_keyFrames.end());
  }

  void addKeyFrame(const Scalar& time, const T& value)
  {
    addKeyFrame(KeyFrame<T>(time, value));
  }

  const T& getValueAtTime(const Scalar& time)
  {
    assert(m_keyFrames.size() > 0);

    if (m_keyFrames.size() == 1) {
      m_currentValue = m_keyFrames[0].getValue();
      return m_currentValue;
    }

    typename KeyFrameArray::iterator pos;
    pos = std::upper_bound(m_keyFrames.begin(),
                           m_keyFrames.end(),
                           KeyFrame<T>(time, T()));

    if (pos == m_keyFrames.begin()) {
      m_currentValue = (*pos).getValue();
      return m_currentValue;
    }

    if (pos == m_keyFrames.end()) {
      m_currentValue = (*(pos - 1)).getValue();
      return m_currentValue;
    }

    const T& valueA = (*(pos - 1)).getValue();
    const T& valueB = (*pos).getValue();

    Scalar timeA = (*(pos - 1)).getTime();
    Scalar timeB = (*pos).getTime();

    Scalar alpha = (time - timeA) / (timeB - timeA);
    assert(alpha >= 0. && alpha <= 1.);

    m_currentValue = (1. - alpha) * valueA + alpha * valueB;
    return m_currentValue;
  }

  void setTime(const Scalar& time)
  {
    m_currentValue = getValueAtTime(time);
  }

  const T& getValue() const
  {
    return m_currentValue;
  }

  void clear()
  {
    m_keyFrames.clear();
  }

  /** Sets the keyframes up from a string. */
  void setKeyFramesFromString(const std::string& specs)
  {
    clear();
    m_specs = specs;

    // groups of keyframes
    std::vector<std::string> groups;
    stringTokenizer(groups, specs, "_");

    for (size_t i = 0; i < groups.size(); ++i) {
      // specify format of keyframe
      std::vector<std::string> format;
      stringTokenizer(format, groups[i], ":");
      if (format.size() == 0) continue;
      if (format.size() == 1) {
        addKeyFrame(0, fromString<T>(format[0]));
        continue;
      }

      if (format.size() != 2) {
        std::cerr << "Keyframe group syntax error: "
                  << groups[i] << std::endl;
        continue;
      }

      std::vector<std::string> values;
      stringTokenizer(values, format[1], ",");

      // time,value
      if (format[0] == "t") {
        if (values.size() % 2 != 0) {
          std::cerr << "Keyframe format syntax error: "
                    << format[1] << std::endl;
          continue;
        }

        for (size_t j = 0; j < values.size(); j += 2) {
          Scalar time = fromString<Scalar>(values[j]);
          T value = fromString<T>(values[j + 1]);
          addKeyFrame(time, value);
        }

        // value,speed,value
      } else if (format[0] == "s") {
        if (values.size() < 2 || values.size() % 2 != 0) {
          std::cerr << "Keyframe value syntax error: "
                    << format[1] << std::endl;
          continue;
        }

        Scalar time = fromString<Scalar>(values[0]);
        T value = fromString<T>(values[1]);
        addKeyFrame(time, value);
        for (size_t j = 2; j < values.size(); j += 2) {
          Scalar speed = fromString<Scalar>(values[j]);
          T nextValue = fromString<T>(values[j + 1]);
          time += (nextValue - value) / speed;
          value = nextValue;
          addKeyFrame(time, value);
        }

      } else {
        std::cerr << "Unknown keyframe format: " << groups[i] << std::endl;
        continue;
      }
    }
  }

  const std::string& getSpecString() const
  {
    return m_specs;
  }

private:
  KeyFrameArray m_keyFrames;
  T m_currentValue;
  std::string m_specs;
};

} // namespace BASim

#endif // KEYEDVALUE_HH
