/**
 * \file Tokenizer.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 11/17/2009
 */

#ifndef TOKENIZER_HH
#define TOKENIZER_HH

namespace BASim {

/** Class to tokenize an input stream */
class Tokenizer
{
public:

  Tokenizer(std::istream& stream) : m_stream(stream) {}

  int read(std::string& str)
  {
    while (m_stream.good()) {
      std::string tok;
      m_stream >> tok;

      if (tok.size() == 0) continue;

      if (tok.compare(0, 2, "//") == 0) {
        m_stream.ignore(256,'\n');
        continue;
      }

      str = tok;
      return 0;
    }

    return -1;
  }



  int read(bool& b)
  {
    std::string str;
    if (read(str) != 0) return -1;

    b = (str == "1") || (str == "true");
    return 0;
  }



  template<class T> int read(T& t)
  {
    std::string str;
    if (read(str) != 0) return -1;

    std::istringstream iss(str);
    return (iss >> t).fail() ? -1 : 0;
  }



  template <class Vec> int readVec(Vec& v)
  {
    for (int i = 0; i < (int) v.size(); ++i) {
      if (read(v(i)) != 0) return -1;
    }

    return 0;
  }

private:

  std::istream& m_stream;
};

} // namespace BASim

#endif // TOKENIZER_HH
