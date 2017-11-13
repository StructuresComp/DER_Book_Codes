/**
 * \file ProblemBase.tcc
 *
 * \author miklos@cs.columbia.edu
 * \date 09/09/2009
 */

#include <iostream>

template <typename T> int
Problem::AddOption(const std::string& name, const std::string& desc,
                   const T& def)
{
  if (m_options.find(name) != m_options.end()) {
    std::cerr << "Option " << name << " already exists" << std::endl;
    return -1;
  }
  m_options.insert(std::make_pair(name, Option(name, desc, def)));
  return 0;
}
