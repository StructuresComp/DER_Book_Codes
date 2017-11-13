/**
 * \file XMLSceneOutputter.hh
 *
 * \author smith@cs.columbia.edu
 * \date 04/18/2010
 */

#ifndef XMLSCENEOUTPUTTER_HH
#define XMLSCENEOUTPUTTER_HH

#include "rapidxml.hpp"
#include "rapidxml_print.hpp"

#include "Apps/BASimulator/Problems/ProblemBase.hh"
#include <iomanip>
#include <fstream>
#include <typeinfo>

// In Windows, these probably won't work :'(
#include <sys/stat.h>
#include <sys/types.h>

namespace BASim
{

class XMLSceneOutputter
{
  
public:
  
  // Saves the current scene (xml description, rods, obj files) to outputdirectory. xmlfilename is the name 
  // of the output xml scene file. 
  void outputScene( const std::string& outputdirectory, const std::string& xmlfilename, const Problem& problem )
  {
    // TODO: should proably add some sanity checks to append / and stuff to outputdirectory.

    // Try to create the user-requested directory
    //if( mkdir(outputdirectory.c_str(), 0755) )
    //{
    //  std::cout << "\033[31;1mWARNING IN XMLSceneOutputter:\033[m Failed to create NEW output folder. Output folder might already exist?" << std::endl;
    //}
    
    const World& world = problem.getWorld();
    const World::Objects& objects = world.getObjects();

    rapidxml::xml_document<> doc;

    // Create an xml header
    rapidxml::xml_node<>* decl = doc.allocate_node(rapidxml::node_declaration);
    decl->append_attribute(doc.allocate_attribute("version", "1.0"));
    decl->append_attribute(doc.allocate_attribute("encoding", "utf-8"));
    doc.append_node(decl);

    // Create the root simulation node
    rapidxml::xml_node<>* root = doc.allocate_node(rapidxml::node_element, "sim");
    doc.append_node(root);
    
    // Create a child node for each rod (save as external files, for now)
    int num_rods = 0;
    int num_objs = 0;
    for( World::Objects::const_iterator it = objects.begin(); it != objects.end(); ++it ) 
    {
      assert( (*it) != NULL );

      ObjectBase* objptr = *it;

      // Elastic rods
      if( typeid(*objptr)==typeid(AnisotropicRod) )
      {
        ElasticRod* rod = dynamic_cast<ElasticRod*>(objptr);
        if (rod == NULL) continue;

        // Create a node in the xml hieararchy
        std::stringstream name;
        name << std::setfill('0');
        name << "rod" << std::setw(10) << num_rods << ".txt";
        
        rapidxml::xml_node<>* child = doc.allocate_node(rapidxml::node_element, "externalobject");
        child->append_attribute(doc.allocate_attribute("type", "rod"));
        child->append_attribute(doc.allocate_attribute("file",doc.allocate_string(name.str().c_str())));
        root->append_node(child);
        
        // Write the file containing the rod
        std::stringstream filename;
        filename << std::setfill('0');
        filename << outputdirectory << "/rod" << std::setw(10) << num_rods << ".txt";
        outputRod( *rod, filename.str() );
        
        ++num_rods;
      }
      // Obj meshes
      else if( typeid(*objptr)==typeid(TriangleMesh) )
      {
        TriangleMesh* triobj = dynamic_cast<TriangleMesh*>(objptr);
        if( triobj == NULL ) continue;
        
        // Create a node in the xml hieararchy
        std::stringstream name;
        name << std::setfill('0');
        name << "obj" << std::setw(10) << num_objs << ".obj";        
        
        rapidxml::xml_node<>* child = doc.allocate_node(rapidxml::node_element, "externalobject");
        child->append_attribute(doc.allocate_attribute("type", "obj"));
        child->append_attribute(doc.allocate_attribute("file",doc.allocate_string(name.str().c_str())));
        root->append_node(child);        
        
        // Write the file containing the obj
        std::stringstream filename;
        filename << std::setfill('0');
        filename << outputdirectory << "/obj" << std::setw(10) << num_objs << ".obj";    
        outputObj( *triobj, filename.str() );
        
        ++num_objs;
      }
      else
      {
        std::cout << "\033[31;1mWARNING IN XMLSceneOutputter:\033[m Encountered nonexistant object type!" << std::endl;
      }
    }
    
    // Write the scene representation to a string
    std::string xml_as_string;
    rapidxml::print(std::back_inserter(xml_as_string), doc);

    // Output the xml to the user-requested file
    std::string xmlfullfilename(outputdirectory);
    xmlfullfilename.append("/");
    xmlfullfilename.append(xmlfilename);
    std::ofstream xmlfileout(xmlfullfilename.c_str());
    if( xmlfileout.is_open() )
    {
      xmlfileout << xml_as_string;
      
      xmlfileout.close();
    }
    else 
    {
      std::cout << "\033[31;1mWARNING IN XMLSceneOutputter:\033[m Failed to open xml output file for writting." << std::endl;
    }
  }
  
  void outputRod( const ElasticRod& rod, const std::string& filename )
  {
    std::ofstream rodoutput(filename.c_str());
    
    if( rodoutput.is_open() )
    {
      ////////////////////////////////
      // Output the number of vertices
      rodoutput << "numverts " << rod.nv() << std::endl << std::endl;

      /////////////////////////////////////////////////////////
      // Output the radii. Assumes constant across rod, for now
      rodoutput << "radiusScale " << rod.getRadiusScale() << std::endl;
      rodoutput << "radiusA " << rod.radiusA(0) << std::endl;
      rodoutput << "radiusB " << rod.radiusB(0) << std::endl << std::endl;
      
      
      ////////////////////////////////////
      // Output the deformed configuration
      rodoutput << "vertices" << std::endl;
      for( int i = 0; i < rod.nv(); ++i )
      {
        rodoutput << "   " << rod.getVertex(i).x() << " " << rod.getVertex(i).y() << " " << rod.getVertex(i).z() << std::endl;
      }
      rodoutput << "endvertices" << std::endl << std::endl;

      ////////////////////////////////////
      // Output the vertex  configuration
      rodoutput << "vertexvelocities" << std::endl;
      for( int i = 0; i < rod.nv(); ++i )
      {
        rodoutput << "   " << rod.getVelocity(i).x() << " " << rod.getVelocity(i).y() << " " << rod.getVelocity(i).z() << std::endl;
      }
      rodoutput << "endvertexvelocities" << std::endl << std::endl;

      ////////////////////////////
      // Output the material frame
      rodoutput << "material1" << std::endl;
      for( int i = 0; i < rod.ne(); ++i )
      {
        rodoutput << "   " << rod.getMaterial1(i).x() << " " << rod.getMaterial1(i).y() << " " << rod.getMaterial1(i).z() << std::endl;
      }
      rodoutput << "endmaterial1" << std::endl << std::endl;
      
      rodoutput.close();
    }
    else
    {
      std::cout << "\033[31;1mWARNING IN XMLSceneOutputter:\033[m Failed to open text rod file for output." << std::endl;
    }
  }
  
  void outputObj( const TriangleMesh& trimesh, const std::string& filename )
  {
    std::ofstream objoutput(filename.c_str());
    
    if( objoutput.is_open() )
    {
      objoutput << "# vertices (1-indexed)" << std::endl;
      for( TriangleMesh::vertex_iter itr = trimesh.vertices_begin(); itr != trimesh.vertices_end(); ++itr )
      {
        objoutput << "v " << trimesh.getVertex(*itr).x() << " " << trimesh.getVertex(*itr).y() << " " << trimesh.getVertex(*itr).z() << std::endl;
      }

      objoutput << "# faces" << std::endl;
      for( TriangleMesh::face_iter itr = trimesh.faces_begin(); itr != trimesh.faces_end(); ++itr )
      {
        objoutput << "f ";
        for( TriangleMesh::FaceVertexIter fvit = trimesh.fv_iter(*itr); fvit; ++fvit )
        {
          objoutput << (1+(*fvit).idx()) << " ";
        }
        objoutput << std::endl;
      }

      objoutput.close();
    }
    else
    {
      std::cout << "\033[31;1mWARNING IN XMLSceneOutputter:\033[m Failed to open obj file for output." << std::endl;
    }
  }
  
};
  
}

#endif


