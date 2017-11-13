/**
 * \file XMLSceneParser.hh
 *
 * \author smith@cs.columbia.edu
 * \date 04/17/2010
 */

#ifndef XMLSCENEPARSER_HH
#define XMLSCENEPARSER_HH

#include "rapidxml.hpp"
#include "RodTextFileParser.hh"

namespace BASim
{

class XMLSceneParser
{

public:

  XMLSceneParser( const std::string& filename )
  {    
    
    // Attempt to read the text from the user-specified xml file
    std::string filecontents;
    if( !loadTextFileIntoString(filename,filecontents) )
    {
      std::cout << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m XML scene file " << filename << ". Failed to locate <sim> node." << std::endl;
      exit(1);
    }
    
    // Copy string into an array of characters for the xml parser
    for( int i = 0; i < (int) filecontents.size(); ++i ) m_xml_chars.push_back(filecontents[i]);
    m_xml_chars.push_back('\0');
    
    // Initialize the xml parser with the character vector
    m_doc.parse<0>(&m_xml_chars[0]);
  }
  
  ~XMLSceneParser()
  {
  }
  
  void loadRodFileNames( std::vector<std::string>& rodfilenames )
  {
    rapidxml::xml_node<>* node = m_doc.first_node("sim");
    if( node == NULL ) 
    {
      std::cout << "\033[31;1mERROR IN XMLSCENEPARSER:\033[m Failed to parse xml scene file. Failed to locate <sim> node." << std::endl;
      exit(1);
    }

    for( rapidxml::xml_node<>* extrnlobjct = node->first_node("externalobject"); extrnlobjct; extrnlobjct = extrnlobjct->next_sibling() )
    {
      if( extrnlobjct->first_attribute("type") ) 
      {
        if( std::string(extrnlobjct->first_attribute("type")->value()) == "rod" )
        {
          if( extrnlobjct->first_attribute("file") )
          {
            std::string rodfilename(extrnlobjct->first_attribute("file")->value());
            if( rodfilename.size() != 0 ) rodfilenames.push_back(rodfilename);
          }
        }
      }
    }
  }
  
private:
  
  bool loadTextFileIntoString( const std::string& filename, std::string& filecontents )
  {
    // Attempt to open the text file for reading
    std::ifstream textfile(filename.c_str(),std::ifstream::in);
    if(!textfile) return false;
    
    // Read the entire file into a single string
    std::string line;
    while(getline(textfile,line)) filecontents.append(line);

    textfile.close();

    return true;
  }

  void loadRodFileNames()
  {
  }
  
  std::vector<char> m_xml_chars;
  rapidxml::xml_document<> m_doc;
};
  
}

#endif


