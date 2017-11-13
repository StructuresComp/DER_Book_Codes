/**
 * \file ObjParser.hh
 *
 * \author smith@cs.columbia.edu
 * \date 08/29/2009
 */

#ifndef OBJPARSER_HH
#define OBJPARSER_HH

#include <string>
#include "BASim/src/Core/TriangleMesh.hh"

namespace BASim 
{
  /** 
   * Loades meshes from OBJ format files. 
   */
  class ObjParser
  {
  public:
    /**
     * Loads a triangle mesh into tri_mesh from OBJ format file obj_file_name.
     *
     * \param[in] obj_file_name String filename of OBJ format file.
     * \param[in] tri_mesh Empty tirangle mesh to load OBJ file into.
     */
    bool loadTriangularMesh( const std::string& obj_file_name, TriangleMesh& tri_mesh );
  private:
    void parsevCommand( std::istringstream& commandstream, TriangleMesh& tri_mesh );
    void parsetrianglefCommand( std::istringstream& commandstream, TriangleMesh& tri_mesh );
    void tokenize( const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters = " " );
    
    bool m_vn_vt_unsupported_printed;
    bool m_vt_unsupported_printed;
    bool m_vn_unsupported_printed;
    bool m_o_unsupported_printed;
    bool m_g_unsupported_printed;
    bool m_mtllib_unsupported;
    bool m_usemtl_unsupported;
    bool m_s_unsupported;
    
    bool m_successful_parse;
    
    int m_line_num;
  };
}

#endif
