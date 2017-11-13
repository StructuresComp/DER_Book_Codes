/**
 * \file RodTextFileParser.hh
 *
 * \author smith@cs.columbia.edu
 * \date 04/18/2010
 */

#ifndef RODTEXFILEPARSER_HH
#define RODTEXFILEPARSER_HH

#include <string>

namespace BASim
{

class RodTextFileParser
{

public:
  
  RodTextFileParser( const std::string& filename )
  :m_numverts_specified(false)
  ,m_undeformed_config_specified(false)
  ,m_deformed_config_specified(false)
  ,m_num_verts(-1)
  ,m_deformed()
  {
    parseRodFile( filename );
  }

  int getNumVerts()
  {
    return m_num_verts;
  }

  void getDeformedVertices( std::vector<Vec3d>& deformed )
  {
    for( int i = 0; i < (int) m_deformed.size(); ++i ) deformed.push_back(m_deformed[i]);
  }
  
private:
  
  void parseNumVerts( std::istream& command_stream )
  {    
    if( command_stream >> m_num_verts  ) 
    {
      // Skip trailing whitespace
      command_stream >> std::ws;
      
      // Extra commands are not expected
      if( !command_stream.eof() )
      {
        std::cout << "\033[31;1mERROR IN RODTEXTFILEPARSER:\033[m Invalid value for numverts." << std::endl;
        exit(1);
      }
      
      // Ensure number of verts is positive
      if( m_num_verts <= 0 )
      {
        std::cout << "\033[31;1mERROR IN RODTEXTFILEPARSER:\033[m Invalid value for numverts." << std::endl;
        exit(1);
      }
    }
    else
    {
      std::cout << "\033[31;1mERROR IN RODTEXTFILEPARSER:\033[m Invalid value for numverts." << std::endl;
      exit(1);
    }	
  }
  
  void parseDeformed( std::istream& command_stream, std::istream& input_file_stream  )
  {    
    command_stream >> std::ws;

    if( !command_stream.eof() )
    {
      std::cout << "\033[31;1mERROR IN RODTEXTFILEPARSER:\033[m Error, deformed command does not take immediate arguements." << std::endl;
      exit(1);
    }

    while( !input_file_stream.eof() )
    {
      std::string position_command;
      getline( input_file_stream, position_command );

      std::istringstream position_stream(position_command);

      std::string command_part;
      if( position_stream >> command_part )
      {
        double x;
        double y;
        double z;
        
        if( command_part == "enddeformed" )
        {
          return;
        }
        
        std::istringstream vertex_stream(position_command);
        
        if( (vertex_stream >> x >> y >> z >> std::ws) && vertex_stream.eof() )
        {
          m_deformed.push_back(Vec3d(x,y,z));
        }
        else
        {
          std::cout << "\033[31;1mERROR IN RODTEXTFILEPARSER:\033[m Error, invalid vertex position specified." << std::endl;
          exit(1);
        }
      }
    }
    
    std::cout << "\033[31;1mERROR IN RODTEXTFILEPARSER:\033[m Error, deformed command never terminated with enddeformed." << std::endl;
    exit(1);
  }
  
  void parseRodFile( const std::string& filename )
  {
    std::ifstream inputrod(filename.c_str());
    
    if( inputrod.is_open() )
    {
      std::string rod_command;
      while( !inputrod.eof() )
      {
        // Read a single line at a time
        getline( inputrod, rod_command );

        // Use a string stream for easy tokenizing
        std::istringstream commandstream(rod_command);

        // First element of a command is the command's name
        std::string command;

        commandstream >> command;

        // Convert to lower case
        for( unsigned int i = 0; i < command.length(); i++ )
        {
          command[i] = (char) tolower( command[i] );
        }

        if( command == "numverts" )
        {
          m_numverts_specified = true;
          this->parseNumVerts( commandstream );
        }
        else if( command == "deformed" )
        {
          m_deformed_config_specified = true;
          // deformed is a multiline command. Thus, pass it string stream with the command name 
          // and the remaining text.
          this->parseDeformed( commandstream, inputrod );
        }
        else
        {
          // If the line wasn't a known command,
          // it could have been a comment (started
          // with #) or whitespace, which are acceptable.
          // Otherwise, throw an exception and generate an
          // error message.
          
          std::istringstream test_stream(command);
          
          // Strip whitespace, see what remains
          test_stream >> std::ws;
          
          // If line isn't empty or isn't a comment, we have
          // an invalid command
          if( !test_stream.eof() && command.c_str()[0] != '#' )
          {
            std::cout << "\033[31;1mERROR IN RODTEXTFILEPARSER:\033[m Invalid command: " << rod_command << std::endl;
            exit(1);
          }
        }        
      }
      inputrod.close();
    }
    else
    {
      std::cout << "\033[31;1mERROR IN RODTEXTFILEPARSER:\033[m Failed to load simple rod file." << std::endl;
      exit(1);
    }
    
  }
  
  void parseDeformedConfig( std::istream& command_stream, std::istream& input_file_stream );
  
  // Whether or not certain properties were specified in rod file
  bool m_numverts_specified;
  bool m_undeformed_config_specified;
  bool m_deformed_config_specified;
  
  int m_num_verts;
  std::vector<Vec3d> m_deformed;
};

}

#endif


