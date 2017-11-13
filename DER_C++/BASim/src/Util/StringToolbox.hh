/**
 * \file StringToolbox.hh
 *
 * \author Nicolas Clauvelin and Basile Audoly
 * \date June 2009
 */

#ifndef STRINGTOOLBOX_HH
#define STRINGTOOLBOX_HH

namespace BASim {

/// string converter
template <class T>
inline T fromString(const std::string& ToConvert)
{
  std::stringstream streamIn(ToConvert);
  T ReturnValue = T();
  streamIn >> ReturnValue;
  return ReturnValue;
}

/// string tokenizer
inline void stringTokenizer(std::vector<std::string>& tokens,
                            const std::string& str,
                            const std::string& delimiters)
{
  // skip delimiters at beginning
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);

  // find first "non-delimiter"
  std::string::size_type pos = str.find_first_of(delimiters, lastPos);

  while (std::string::npos != pos || std::string::npos != lastPos) {
    // found a token, add it to the vector
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    // skip delimiters
    lastPos = str.find_first_not_of(delimiters, pos);
    // find next non-delimiter
    pos = str.find_first_of(delimiters, lastPos);
  }
}

} // namespace BASim

#endif // STRINGTOOLBOX_HH
