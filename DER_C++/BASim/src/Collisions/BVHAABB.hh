/**
 * \file BVHAABB.hh
 *
 * \author smith@cs.columbia.edu
 * \date 03/21/2010
 */

#ifndef BVHAABB_HH
#define BVHAABB_HH

#include <set>

namespace BASim 
{
  
class IntPair
{
public:

  IntPair( const int& i , const int& j )
  {
    m_i = std::min(i,j);
    m_j = std::max(i,j);
  }

	bool operator==( const IntPair& rhs ) const
  {
    assert( m_i <= m_j );
    assert( rhs.m_i <= rhs.m_j );
    
    return m_i == rhs.m_i && m_j == rhs.m_j;
  }
  
  bool operator!=( const IntPair& rhs ) const
  {
    assert( m_i <= m_j );
    assert( rhs.m_i <= rhs.m_j );
    
    return !(*this==rhs);
  }
  
  bool operator<( const IntPair& rhs ) const
  {
    assert( m_i <= m_j );
    assert( rhs.m_i <= rhs.m_j );
    
    if( m_i != rhs.m_i ) return m_i < rhs.m_i;
    else return m_j < rhs.m_j;
  }
  
  int first() const { return m_i; }
  int second() const { return m_j; }

private:
  int m_i;
  int m_j;
};
  
/**
 * Struct to store face indices.
 */
struct TriangularFace
{
  int idx[3];
};

/**
 * Struct to store information needed to resolve a "proximity" collision between two edges.
 */
struct EdgeEdgeProximityCollision
{
  // Vertices involved in collision
  int e0_v0;
  int e0_v1;
  int e1_v0;
  int e1_v1;
  // Radii of edge 0 and edge 1
  double r0;
  double r1;
  // Barycentric coordinates of closest points between edges
  double s;
  double t;
  // Penetration depth (sum of radii - minimum distance)
  double pen;
  // Collision normal
  Vec3d n;
};

/**
 * Struct to store information needed to resolve a "proximity" collision between a vertex and a face.
 */
struct VertexFaceProximityCollision
{
  // Index of vertex
  int v0;
  // Index of face vertices
  int f0;
  int f1;
  int f2;
  // Radii of vertex and face
  double r0;
  double r1;
  // Barycentric coordinates of closest point on triangle
  double u;
  double v;
  double w;
  // Penetration depth (sum of radii - minimum distance)
  double pen;
  // Collision normal, points TOWARDS the vertex
  Vec3d n;
};  

/**
 * Struct to store information needed to resolve a continuous collision between two edges.
 */
struct EdgeEdgeContinuousTimeCollision
{
  // Vertices involved in collision
  int e0_v0;
  int e0_v1;
  int e1_v0;
  int e1_v1;
  // Barycentric coordinates of closest points between edges
  double s;
  double t;
  // Penetration depth (sum of radii - minimum distance)
  //double pen;
  // Collision normal
  Vec3d n;
};

/**
 * Struct to store information needed to resolve a continuous time collision between a vertex
 * and a face. 
 */
struct VertexFaceContinuousTimeCollision
{
  // Index of vertex
  int v0;
  // Index of face vertices
  int f0;
  int f1;
  int f2;
  // Barycentric coordinates of closest point on triangle
  double u;
  double v;
  double w;
  // Collision normal, points TOWARDS the vertex
  Vec3d n;
};

// TODO: Rework this so I only store face, edge indcs for leaf nodes
struct BVHAABBNode
{
  BVHAABBNode* children[2];
  double lwr_bnds[3];
  double upr_bnds[3];
  
  std::vector<int> edge_indcs;
  std::vector<int> face_indcs;
};

class BVHAABB
{
public:
  BVHAABB( const VecXd& x, const VecXd& v, const std::vector<std::pair<int,int> >& edges, const std::vector<TriangularFace>& faces, const std::vector<double>& vert_radii, const double& timestep )
  :m_root_node(NULL)
  ,m_x(x)
  ,m_v(v)
  ,m_edges(edges)
  ,m_faces(faces)
  ,m_vert_radii(vert_radii)
  ,m_h(timestep)
  {
    assert( m_x.size() == (int) (3*m_vert_radii.size()) );
    assert( m_x.size() == m_v.size() );
    regenerateBVH();
  }
  
  ~BVHAABB()
  {
    if( m_root_node != NULL )
    {
      deleteTree( m_root_node );
      m_root_node = NULL;
    }
  }
  
  //void regenerateBVH( const VecXd& x, const std::vector<std::pair<int,int> >& edges, const std::vector<TriangularFace>& faces, const std::vector<double>& vert_radii )
  void regenerateBVH()
  {
    //m_x = x;
    //m_edges = edges;
    //m_faces = faces;
    //m_vert_radii = vert_radii;
    
    // TODO: Add sanity checks here. 
    
    //std::cout << "Regenerating BVH" << std::endl;
    
    if( m_x.size() == 0 ) return;
    
    // Create the root node
    if( m_root_node != NULL ) deleteTree( m_root_node );
    m_root_node = new BVHAABBNode;
    m_root_node->children[0] = m_root_node->children[1] = NULL;
    
    // Compute the bounds on the root node
    std::vector<int> edge_indcs; for( int i = 0; i < (int) m_edges.size(); ++i ) edge_indcs.push_back(i);
    std::vector<int> face_indcs; for( int i = 0; i < (int) m_faces.size(); ++i ) face_indcs.push_back(i);
    computeProximityBounds( edge_indcs, face_indcs, m_root_node );
    
    // Create the remainder of the tree
    divideNode( edge_indcs, face_indcs, m_root_node );
    //recomputeProximityBoundsBottomUp( m_root_node );
    
    
    //std::cout << "m_x.size(): " << m_x.size() << std::endl;
    //std::cout << "    x: " << m_x.transpose() << std::endl;
    //std::cout << "m_vert_radii.size(): " << m_vert_radii.size() << std::endl;
    //std::cout << "m_edges.size(): " << m_edges.size() << std::endl;
    //std::cout << "m_faces.size(): " << m_faces.size() << std::endl;
    //printTreeInfo( m_root_node, 0, "Root node" );

    assert( treeBuiltCorrectly() );
  }
  
  void getProximityCollisions( std::vector<EdgeEdgeProximityCollision>& edg_edg, std::vector<VertexFaceProximityCollision>& vrt_fce )
  {
    if( m_x.size() == 0 ) return;
    
    assert( m_root_node != NULL );
    assert( treeBuiltCorrectly() );
    
    recomputeProximityBoundsBottomUp( m_root_node );

    //std::set<int> edges;
    //std::set<int> faces;
    //getContainedObjects( m_root_node, edges, faces );
    //std::cout << edges.size() << std::endl;

    //std::cout << std::endl;
    //std::cout << "Number of levels: " << computeDepth( m_root_node ) << std::endl;
    //printTreeInfo( m_root_node, 0, "Root Node" );
    //std::cout << std::endl;

    std::set<IntPair> edg_edg_set;
    std::set<std::pair<int,int> > vrt_fce_set;
    std::set<std::pair<int,int> > edg_fce_set;
    std::set<IntPair> fce_fce_set;
    computePossibleProximityCollisions( edg_edg_set, vrt_fce_set, edg_fce_set, fce_fce_set, m_root_node, m_root_node );

    // TODO: Run some tests to ensure no duplicates in the four sets above
    
    for( std::set<IntPair>::const_iterator it = edg_edg_set.begin(); it != edg_edg_set.end(); ++it )
    {
      EdgeEdgeProximityCollision edg_edg_prxm;
      edg_edg_prxm.e0_v0 = m_edges[it->first()].first;
      edg_edg_prxm.e0_v1 = m_edges[it->first()].second;
      edg_edg_prxm.e1_v0 = m_edges[it->second()].first;
      edg_edg_prxm.e1_v1 = m_edges[it->second()].second;
      
      edg_edg_prxm.r0 = m_vert_radii[edg_edg_prxm.e0_v0];
      edg_edg_prxm.r1 = m_vert_radii[edg_edg_prxm.e1_v0];

      edg_edg.push_back(edg_edg_prxm);
    }
    
    for( std::set<std::pair<int,int> >::const_iterator it = vrt_fce_set.begin(); it != vrt_fce_set.end(); ++it )
    {
      VertexFaceProximityCollision vrt_fce_prxm;

      assert( it->first >= 0 ); assert( it->first < (int) m_vert_radii.size() );
      vrt_fce_prxm.v0 = it->first;

      assert( it->second >= 0 ); assert( it->second < (int) m_faces.size() );
      
      assert( m_faces[it->second].idx[0] >= 0 ); assert( m_faces[it->second].idx[0] < (int) m_vert_radii.size() );
      vrt_fce_prxm.f0 = m_faces[it->second].idx[0];
      assert( m_faces[it->second].idx[1] >= 0 ); assert( m_faces[it->second].idx[1] < (int) m_vert_radii.size() );
      vrt_fce_prxm.f1 = m_faces[it->second].idx[1];
      assert( m_faces[it->second].idx[2] >= 0 ); assert( m_faces[it->second].idx[2] < (int) m_vert_radii.size() );
      vrt_fce_prxm.f2 = m_faces[it->second].idx[2];

      vrt_fce_prxm.r0 = m_vert_radii[vrt_fce_prxm.v0];
      vrt_fce_prxm.r1 = m_vert_radii[vrt_fce_prxm.f0];
      
      vrt_fce.push_back(vrt_fce_prxm);
    }
    
    for( std::set<std::pair<int,int> >::const_iterator it = edg_fce_set.begin(); it != edg_fce_set.end(); ++it )
    {
      EdgeEdgeProximityCollision edg_edg_prxm0;
      edg_edg_prxm0.e0_v0 = m_edges[it->first].first;
      edg_edg_prxm0.e0_v1 = m_edges[it->first].second;
      edg_edg_prxm0.e1_v0 = m_faces[it->second].idx[0];
      edg_edg_prxm0.e1_v1 = m_faces[it->second].idx[1];
      
      edg_edg_prxm0.r0 = m_vert_radii[edg_edg_prxm0.e0_v0];
      edg_edg_prxm0.r1 = m_vert_radii[edg_edg_prxm0.e1_v0];
      
      edg_edg.push_back(edg_edg_prxm0);
      
      
      EdgeEdgeProximityCollision edg_edg_prxm1;
      edg_edg_prxm1.e0_v0 = m_edges[it->first].first;
      edg_edg_prxm1.e0_v1 = m_edges[it->first].second;
      edg_edg_prxm1.e1_v0 = m_faces[it->second].idx[1];
      edg_edg_prxm1.e1_v1 = m_faces[it->second].idx[2];
      
      edg_edg_prxm1.r0 = m_vert_radii[edg_edg_prxm1.e0_v0];
      edg_edg_prxm1.r1 = m_vert_radii[edg_edg_prxm1.e1_v0];
      
      edg_edg.push_back(edg_edg_prxm1);
      
      
      EdgeEdgeProximityCollision edg_edg_prxm2;
      edg_edg_prxm2.e0_v0 = m_edges[it->first].first;
      edg_edg_prxm2.e0_v1 = m_edges[it->first].second;
      edg_edg_prxm2.e1_v0 = m_faces[it->second].idx[2];
      edg_edg_prxm2.e1_v1 = m_faces[it->second].idx[0];
      
      edg_edg_prxm2.r0 = m_vert_radii[edg_edg_prxm2.e0_v0];
      edg_edg_prxm2.r1 = m_vert_radii[edg_edg_prxm2.e1_v0];
      
      edg_edg.push_back(edg_edg_prxm2);
    }
  }

  void getContinuousTimeCollisions( std::vector<EdgeEdgeContinuousTimeCollision>& edg_edg, std::vector<VertexFaceContinuousTimeCollision>& vrt_fce )
  {
    if( m_x.size() == 0 ) return;
    
    assert( m_root_node != NULL );
    assert( treeBuiltCorrectly() );
    
    recomputeContinuousTimeBoundsBottomUp( m_root_node );
    
    std::set<IntPair> edg_edg_set;
    std::set<std::pair<int,int> > vrt_fce_set;
    std::set<std::pair<int,int> > edg_fce_set;
    std::set<IntPair> fce_fce_set;
    computePossibleProximityCollisions( edg_edg_set, vrt_fce_set, edg_fce_set, fce_fce_set, m_root_node, m_root_node );
    
    // TODO: Run some tests to ensure no duplicates in the four sets above
    
    for( std::set<IntPair>::const_iterator it = edg_edg_set.begin(); it != edg_edg_set.end(); ++it )
    {
      EdgeEdgeContinuousTimeCollision edg_edg_ct;
      edg_edg_ct.e0_v0 = m_edges[it->first()].first;
      edg_edg_ct.e0_v1 = m_edges[it->first()].second;
      edg_edg_ct.e1_v0 = m_edges[it->second()].first;
      edg_edg_ct.e1_v1 = m_edges[it->second()].second;
      
      edg_edg.push_back(edg_edg_ct);
    }

    //if( vrt_fce_set.size() != 0 ) std::cout << "vrt_fce_set size: " << vrt_fce_set.size() << std::endl;
    
    for( std::set<std::pair<int,int> >::const_iterator it = vrt_fce_set.begin(); it != vrt_fce_set.end(); ++it )
    {
      VertexFaceContinuousTimeCollision vrt_fce_ct;
      
      vrt_fce_ct.v0 = it->first;
      vrt_fce_ct.f0 = m_faces[it->second].idx[0];
      vrt_fce_ct.f1 = m_faces[it->second].idx[1];
      vrt_fce_ct.f2 = m_faces[it->second].idx[2];
      
      vrt_fce.push_back(vrt_fce_ct);
    }

    for( std::set<std::pair<int,int> >::const_iterator it = edg_fce_set.begin(); it != edg_fce_set.end(); ++it )
    {
      EdgeEdgeContinuousTimeCollision edg_edg_ct0;
      edg_edg_ct0.e0_v0 = m_edges[it->first].first;
      edg_edg_ct0.e0_v1 = m_edges[it->first].second;
      edg_edg_ct0.e1_v0 = m_faces[it->second].idx[0];
      edg_edg_ct0.e1_v1 = m_faces[it->second].idx[1];
      edg_edg.push_back(edg_edg_ct0);

      EdgeEdgeContinuousTimeCollision edg_edg_ct1;
      edg_edg_ct1.e0_v0 = m_edges[it->first].first;
      edg_edg_ct1.e0_v1 = m_edges[it->first].second;
      edg_edg_ct1.e1_v0 = m_faces[it->second].idx[1];
      edg_edg_ct1.e1_v1 = m_faces[it->second].idx[2];
      edg_edg.push_back(edg_edg_ct1);

      EdgeEdgeContinuousTimeCollision edg_edg_ct2;
      edg_edg_ct2.e0_v0 = m_edges[it->first].first;
      edg_edg_ct2.e0_v1 = m_edges[it->first].second;
      edg_edg_ct2.e1_v0 = m_faces[it->second].idx[2];
      edg_edg_ct2.e1_v1 = m_faces[it->second].idx[0];
      edg_edg.push_back(edg_edg_ct2);
    }
  }
  
private:
  
  /////////////////////////////////////////////////////////////////////////////
  // DEBUG METHODS
  
  void printTreeInfo( BVHAABBNode* node, int depth, const std::string left_or_right )
  {
    assert( node != NULL );
    
    if( isLeaf(node) ) std::cout << "LEAF  ";
    
    std::cout << left_or_right << "     Depth: " << depth << std::endl;
    if( node->children[0] == NULL && node->children[1] == NULL )
    {
      std::cout << "     x bounds: " << node->lwr_bnds[0] << " " << node->upr_bnds[0] << std::endl;
      std::cout << "     y bounds: " << node->lwr_bnds[1] << " " << node->upr_bnds[1] << std::endl;
      std::cout << "     z bounds: " << node->lwr_bnds[2] << " " << node->upr_bnds[2] << std::endl;
      int numedges = 0; int numfaces = 0; computeNumEnclosedObjects(node,numedges,numfaces);
      std::cout << "     Num objects: " << numedges << " " << numfaces << std::endl;
      std::cout << "     Edges: "; for( int i = 0; i < (int) node->edge_indcs.size(); ++i ) std::cout << node->edge_indcs[i] << " "; std::cout << std::endl;
      std::cout << "     Faces: "; for( int i = 0; i < (int) node->face_indcs.size(); ++i ) std::cout << node->face_indcs[i] << " "; std::cout << std::endl;
      std::cout << "     No children " << std::endl;
      return;
    }
    else if( node->children[0] != NULL && node->children[1] == NULL )
    {
      std::cout << "     Left child only " << std::endl;
      printTreeInfo( node->children[0], depth+1, "Left Child" );
      return;
    }
    else if( node->children[0] == NULL && node->children[1] != NULL )
    {
      std::cout << "     Right child only " << std::endl;
      printTreeInfo( node->children[1], depth+1, "Right Child" );
      return;
    }
    else 
    {
      std::cout << "     Two children " << std::endl;
      printTreeInfo( node->children[0], depth+1, "Left Child" );
      printTreeInfo( node->children[1], depth+1, "Right Child" );
      return;
    }
  }
  
  void getContainedObjects( BVHAABBNode* node, std::set<int>& edges, std::set<int>& faces )
  {
    assert( node != NULL );
    
    if( isLeaf(node) )
    {
      for( int i = 0; i < (int) node->edge_indcs.size(); ++i ) edges.insert(node->edge_indcs[i]);
      for( int i = 0; i < (int) node->face_indcs.size(); ++i ) faces.insert(node->face_indcs[i]);
    }
    if( hasLeftChild(node) ) 
    {
      getContainedObjects( node->children[0], edges, faces );
    }
    if( hasRightChild(node) ) 
    {
      getContainedObjects( node->children[1], edges, faces );
    }
  }
  
  int computeDepth( BVHAABBNode* node )
  {
    assert( node != NULL );

    if( isLeaf(node) ) return 1;
    else if(  hasLeftChild(node) && !hasRightChild(node) ) return computeDepth(node->children[0])+1;
    else if( !hasLeftChild(node) &&  hasRightChild(node) ) return computeDepth(node->children[1])+1;
    else return std::max(computeDepth(node->children[0]),computeDepth(node->children[1]))+1;
  }
  
  bool hasProperChildren( BVHAABBNode* node )
  {
    assert( node != NULL );
    
    // If we found a leaf, leaf is correct
    if( isLeaf(node) ) return true;
    // If we found a node that only has one child, incorrect
    if( (!hasLeftChild(node) && hasRightChild(node)) || (hasLeftChild(node) && !hasRightChild(node)) ) return false;
    // Otherwise, node is correct if both children are correct
    return hasProperChildren( node->children[0] ) && hasProperChildren( node->children[1] );
  }
  
  bool onlyLeavesHaveObjects( BVHAABBNode* node )
  {
    assert( node != NULL );

    // If we found a leaf, it should have at least one object
    if( isLeaf(node) ) return (node->edge_indcs.size()+node->face_indcs.size())>0;

    // Otherwise check that this nodes children have proper number of objects
    if( hasLeftChild(node) && hasRightChild(node) ) return onlyLeavesHaveObjects( node->children[0] ) && onlyLeavesHaveObjects( node->children[1] );
    // In a properly built tree, shouldn't his these cases. But throw them in so the recursion terminates nicely if called on
    // a badly built tree.
    if( hasLeftChild(node) )  return onlyLeavesHaveObjects( node->children[0] );
    if( hasRightChild(node) ) return onlyLeavesHaveObjects( node->children[1] );
    
    assert( false ); 
    return false;
  }
  
  void computeNumEnclosedObjects( BVHAABBNode* node, int& num_edges, int& num_faces )
  {
    assert( node != NULL );

    // If we found a leaf, add to the total number of objects
    if( isLeaf(node) )
    {
      num_edges += node->edge_indcs.size();
      num_faces += node->face_indcs.size();
    }
    // Otherwise recurse on the children
    if( hasLeftChild(node)  ) computeNumEnclosedObjects( node->children[0], num_edges, num_faces );
    if( hasRightChild(node) ) computeNumEnclosedObjects( node->children[1], num_edges, num_faces );
  }
  
  bool parentsContainChildren( BVHAABBNode* node )
  {
    assert( node != NULL );

    // Check the left child
    if( hasLeftChild(node) )
    {
      if( node->lwr_bnds[0] > node->children[0]->lwr_bnds[0] ) return false;
      if( node->lwr_bnds[1] > node->children[0]->lwr_bnds[1] ) return false;
      if( node->lwr_bnds[2] > node->children[0]->lwr_bnds[2] ) return false;
      
      if( node->upr_bnds[0] < node->children[0]->upr_bnds[0] ) return false;
      if( node->upr_bnds[1] < node->children[0]->upr_bnds[1] ) return false;
      if( node->upr_bnds[2] < node->children[0]->upr_bnds[2] ) return false;
      
      if( !parentsContainChildren( node->children[0] ) ) return false;
    }
    // Check the right child
    if( hasRightChild(node) )
    {
      if( node->lwr_bnds[0] > node->children[1]->lwr_bnds[0] ) return false;
      if( node->lwr_bnds[1] > node->children[1]->lwr_bnds[1] ) return false;
      if( node->lwr_bnds[2] > node->children[1]->lwr_bnds[2] ) return false;
      
      if( node->upr_bnds[0] < node->children[1]->upr_bnds[0] ) return false;
      if( node->upr_bnds[1] < node->children[1]->upr_bnds[1] ) return false;
      if( node->upr_bnds[2] < node->children[1]->upr_bnds[2] ) return false;
      
      if( !parentsContainChildren( node->children[1] ) ) return false;
    }

    return true;
  }

  bool childrenPointersSetupCorrectly( BVHAABBNode* node )
  {
    assert( node != NULL );
    
    // Leaves have no children
    if( isLeaf(node) ) return true;
    
    // If the children are equal, this is incorect
    if( node->children[0] == node->children[1] ) return false;
    
    if( hasLeftChild(node) ) if( !childrenPointersSetupCorrectly(node->children[0]) ) return false;
    
    if( hasRightChild(node) ) if( !childrenPointersSetupCorrectly(node->children[1]) ) return false;

    return true;
  }
  
  bool treeBuiltCorrectly()
  {
    // Ensure the children of each node are different
    bool children_different = childrenPointersSetupCorrectly(m_root_node);
    
    // Ensure that the each node either has two children or is a leaf
    bool proper_children = hasProperChildren( m_root_node );
    
    // Ensure that only leaves have objects
    bool leaves_have_objects = onlyLeavesHaveObjects( m_root_node );
    
    // Ensure that proper number of objects are in the BVH (can tighten this up a bit more by checking that all objects are accounted for)
    int num_edges = 0; int num_faces = 0;
    computeNumEnclosedObjects( m_root_node, num_edges, num_faces );
    bool proper_num_enclosed = (num_edges==(int)m_edges.size())&&(num_faces==(int)m_faces.size());
    
    // Ensure all objects are enclosed in the BVH
    std::set<int> edges;
    std::set<int> faces;
    getContainedObjects( m_root_node, edges, faces );
    bool all_objects_enclosed = (edges.size()==m_edges.size())&&(faces.size()==m_faces.size());
    // Ensure all objects enclosed in the BVH should be there :)
    bool all_object_valid = true;
    for( std::set<int>::const_iterator i = edges.begin(); i != edges.end(); ++i ) if( *i >= (int) m_edges.size() ) all_object_valid = false;
    for( std::set<int>::const_iterator i = faces.begin(); i != faces.end(); ++i ) if( *i >= (int) m_faces.size() ) all_object_valid = false;
    
    // Ensure the children are contained in the parents
    bool children_in_parents = parentsContainChildren( m_root_node );
    
    return children_different && proper_children && leaves_have_objects && proper_num_enclosed && all_objects_enclosed && all_object_valid && children_in_parents;
  }

  
  /////////////////////////////////////////////////////////////////////////////
  // Methods (M-E-T-H-O-D MAN, M-E-T-H-O-D MAN)

  bool isLeaf( BVHAABBNode* node )
  {
    assert( node != NULL );
    return node->children[0]==NULL && node->children[1]==NULL;
  }
  
  bool hasLeftChild( BVHAABBNode* node )
  {
    assert( node != NULL );
    return node->children[0]!=NULL;
  }
  
  bool hasRightChild( BVHAABBNode* node )
  {
    assert( node != NULL );
    return node->children[1]!=NULL;
  }
  
  int numContainedObjects( BVHAABBNode* node )
  {
    assert( node != NULL );
    return node->edge_indcs.size() + node->face_indcs.size();
  }
  
  double computeNodeVolume( BVHAABBNode* node )
  {
    assert( node != NULL );
    return (node->upr_bnds[0]-node->lwr_bnds[0])*(node->upr_bnds[1]-node->lwr_bnds[1])*(node->upr_bnds[2]-node->lwr_bnds[2]);
  }

  
  void computePossibleProximityCollisions( std::set<IntPair>& edg_edg, std::set<std::pair<int,int> >& vrt_fce, std::set<std::pair<int,int> >& edg_fce, std::set<IntPair>& fce_fce, BVHAABBNode* node_a, BVHAABBNode* node_b )
  {
    assert( node_a != NULL ); assert( node_b != NULL );
    // If the bounding volumes do not overlap, there are no possible collisions between their objects
    if( !AABBsOverlap( node_a, node_b ) ) return;
    // If both bounding volumes are leaves, add their contents to list potential collisions
    if( isLeaf(node_a) && isLeaf(node_b) ) compareContentsProximity(edg_edg,vrt_fce,edg_fce,fce_fce,node_a,node_b);
    // If one bounding volume is a leaf, we must recurse on the other volume
    else if( isLeaf(node_a) )
    {
      computePossibleProximityCollisions( edg_edg, vrt_fce, edg_fce, fce_fce, node_a, node_b->children[0] );
      computePossibleProximityCollisions( edg_edg, vrt_fce, edg_fce, fce_fce, node_a, node_b->children[1] );
    }
    else if( isLeaf(node_b) )
    {
      computePossibleProximityCollisions( edg_edg, vrt_fce, edg_fce, fce_fce, node_a->children[0], node_b );
      computePossibleProximityCollisions( edg_edg, vrt_fce, edg_fce, fce_fce, node_a->children[1], node_b );
    }
    else 
    {
      // Recurse on the larger volume node
      if( computeNodeVolume( node_a ) >= computeNodeVolume( node_b ) ) 
      {
        computePossibleProximityCollisions( edg_edg, vrt_fce, edg_fce, fce_fce, node_a->children[0], node_b );
        computePossibleProximityCollisions( edg_edg, vrt_fce, edg_fce, fce_fce, node_a->children[1], node_b );
      }
      else
      {
        computePossibleProximityCollisions( edg_edg, vrt_fce, edg_fce, fce_fce, node_a, node_b->children[0] );
        computePossibleProximityCollisions( edg_edg, vrt_fce, edg_fce, fce_fce, node_a, node_b->children[1] );
      }
    }
  }
  
  bool AABBsOverlap( BVHAABBNode* node_a, BVHAABBNode* node_b )
  {
    assert( node_a != NULL ); assert( node_b != NULL );
    
    // Attempt to find a separating hyperplane
    if( node_a->upr_bnds[0] < node_b->lwr_bnds[0] ) return false; 
    if( node_b->upr_bnds[0] < node_a->lwr_bnds[0] ) return false; 
    if( node_a->upr_bnds[1] < node_b->lwr_bnds[1] ) return false; 
    if( node_b->upr_bnds[1] < node_a->lwr_bnds[1] ) return false; 
    if( node_a->upr_bnds[2] < node_b->lwr_bnds[2] ) return false; 
    if( node_b->upr_bnds[2] < node_a->lwr_bnds[2] ) return false; 
    
    // Otherwise the bounding volumes overlap!
    return true;
  }
  
  void compareContentsProximity( std::set<IntPair>& edg_edg, std::set<std::pair<int,int> >& vrt_fce, std::set<std::pair<int,int> >& edg_fce, std::set<IntPair>& fce_fce, BVHAABBNode* node_a, BVHAABBNode* node_b )
  {
    assert( node_a != NULL ); assert( node_b != NULL );
    assert( isLeaf(node_a) ); assert( isLeaf(node_b) );
    assert( numContainedObjects(node_a) > 0 ); assert( numContainedObjects(node_b) > 0 );
    
    // Compare the edges in node_a to the edges in node_b
    compareEdgeLists( node_a->edge_indcs, node_b->edge_indcs, edg_edg );

    // Compare the edges in node_a to the faces in node_b
    compareEdgeFaceLists( node_a->edge_indcs, node_b->face_indcs, edg_edg, vrt_fce, edg_fce );
    
    // Compare the edges in node_b to the faces in node_a
    compareEdgeFaceLists( node_b->edge_indcs, node_a->face_indcs, edg_edg, vrt_fce, edg_fce );
    
    // Compare the faces in node_a to the faces in node_b
    compareFaceLists( node_a->face_indcs, node_b->face_indcs, edg_edg, vrt_fce, fce_fce );
  }
  
  void compareEdgeLists( std::vector<int>& edges_a, std::vector<int>& edges_b, std::set<IntPair>& edg_edg )
  {
    for( int i = 0; i < (int) edges_a.size(); ++i ) for( int j = 0; j < (int) edges_b.size(); ++j )
    {
      assert( edges_a[i] >= 0 ); assert( edges_a[i] < (int) m_edges.size() );
      assert( m_edges[edges_a[i]].first  < (int) m_vert_radii.size() );
      assert( m_edges[edges_a[i]].second < (int) m_vert_radii.size() );
      assert( edges_b[j] >= 0 ); assert( edges_b[j] < (int) m_edges.size() );
      assert( m_edges[edges_b[i]].first  < (int) m_vert_radii.size() );
      assert( m_edges[edges_b[i]].second < (int) m_vert_radii.size() );

      if( edges_a[i] == edges_b[j] ) continue;
      
      edg_edg.insert(IntPair(edges_a[i],edges_b[j]));
    }
  }
  
  void compareEdgeFaceLists( std::vector<int>& edges, std::vector<int>& faces, std::set<IntPair>& edg_edg, std::set<std::pair<int,int> >& vrt_fce, std::set<std::pair<int,int> >& edg_fce )
  {
    for( int i = 0; i < (int) edges.size(); ++i ) for( int j = 0; j < (int) faces.size(); ++j )
    {
      assert( edges[i] >= 0 ); assert( edges[i] < (int) m_edges.size() );
      assert( m_edges[edges[i]].first  >= 0 ); assert( m_edges[edges[i]].first  < (int) m_vert_radii.size() );
      assert( m_edges[edges[i]].second >= 0 ); assert( m_edges[edges[i]].second < (int) m_vert_radii.size() );
      assert( faces[j] >= 0 ); assert( faces[j] < (int) m_faces.size() );
      assert( m_faces[faces[j]].idx[0] >= 0 ); assert( m_faces[faces[j]].idx[0] < (int) m_vert_radii.size() );
      assert( m_faces[faces[j]].idx[1] >= 0 ); assert( m_faces[faces[j]].idx[1] < (int) m_vert_radii.size() );
      assert( m_faces[faces[j]].idx[2] >= 0 ); assert( m_faces[faces[j]].idx[2] < (int) m_vert_radii.size() );

      // Each edge has two vertices that possibly collide with the face
      vrt_fce.insert(std::pair<int,int>(m_edges[edges[i]].first,faces[j]));
      vrt_fce.insert(std::pair<int,int>(m_edges[edges[i]].second,faces[j]));
      
      // Each edge could possibly collide with each edge of the face
      edg_fce.insert(std::pair<int,int>(edges[i],faces[j]));
    }
  }
  
  void compareFaceLists( std::vector<int>& faces_a, std::vector<int>& faces_b, std::set<IntPair>& edg_edg, std::set<std::pair<int,int> >& vrt_fce, std::set<IntPair>& fce_fce )
  {
    // TODO: Implement this method!
    //for( int i = 0; i < (int) faces_a.size(); ++i ) for( int j = 0; j < (int) faces_b.size(); ++j )
    //{
    //  assert( faces_a[i] >= 0 ); assert( faces_a[i] < (int) m_faces.size() );
    //  assert( m_faces[faces_a[i]].idx[0] >= 0 ); assert( m_faces[faces_a[i]].idx[0] < (int) m_vert_radii.size() );
    //  assert( m_faces[faces_a[i]].idx[1] >= 0 ); assert( m_faces[faces_a[i]].idx[1] < (int) m_vert_radii.size() );
    //  assert( m_faces[faces_a[i]].idx[2] >= 0 ); assert( m_faces[faces_a[i]].idx[2] < (int) m_vert_radii.size() );
    //  assert( faces_b[j] >= 0 ); assert( faces_b[j] < (int) m_faces.size() );
    //  assert( m_faces[faces_b[j]].idx[0] >= 0 ); assert( m_faces[faces_b[j]].idx[0] < (int) m_vert_radii.size() );
    //  assert( m_faces[faces_b[j]].idx[1] >= 0 ); assert( m_faces[faces_b[j]].idx[1] < (int) m_vert_radii.size() );
    //  assert( m_faces[faces_b[j]].idx[2] >= 0 ); assert( m_faces[faces_b[j]].idx[2] < (int) m_vert_radii.size() );
    //  
    //  if( faces_a[i] == faces_b[j] ) continue;
    //  
    //  // Compare all of face a's verts to face b
    //  //std::pair<int,int>(m_faces[faces_a[i]].idx[0],faces_b[j]);
    //  //std::pair<int,int>(m_faces[faces_a[i]].idx[1],faces_b[j]);
    //  //std::pair<int,int>(m_faces[faces_a[i]].idx[2],faces_b[j]);
    //  
    //  // Compare all of face b's verts to face a
    //  //std::pair<int,int>(m_faces[faces_b[j]].idx[0],faces_a[i]);
    //  //std::pair<int,int>(m_faces[faces_b[j]].idx[1],faces_a[i]);
    //  //std::pair<int,int>(m_faces[faces_b[j]].idx[2],faces_a[i]);
    //}
  }
  
  
  
  
  void deleteTree( BVHAABBNode* node )
  {
    assert( node != NULL );
    
    // Delete the left sub-tree
    if( hasLeftChild(node) )
    {
      deleteTree( node->children[0] );
      node->children[0] = NULL;
    }
    
    // Delete the right sub-tree
    if( hasRightChild(node) )
    {
      deleteTree( node->children[1] );
      node->children[1] = NULL;
    }
    
    // Delete the node itself
    delete node;
    node = NULL;
  }
  
  void recomputeProximityBoundsBottomUp( BVHAABBNode* node )
  {
    assert( node != NULL );

    // If node is a leaf, update the bounding box for that node
    if( isLeaf(node) ) computeProximityBounds( node->edge_indcs, node->face_indcs, node );
    // Otherwise update the node's children bounding box's and compute based off them
    else 
    {
      assert( hasLeftChild(node) ); 
      assert( hasRightChild(node) );
      
      recomputeProximityBoundsBottomUp( node->children[0] );
      recomputeProximityBoundsBottomUp( node->children[1] );
      
      node->lwr_bnds[0] = std::min(node->children[0]->lwr_bnds[0],node->children[1]->lwr_bnds[0]);
      node->lwr_bnds[1] = std::min(node->children[0]->lwr_bnds[1],node->children[1]->lwr_bnds[1]);
      node->lwr_bnds[2] = std::min(node->children[0]->lwr_bnds[2],node->children[1]->lwr_bnds[2]);
      
      node->upr_bnds[0] = std::max(node->children[0]->upr_bnds[0],node->children[1]->upr_bnds[0]);
      node->upr_bnds[1] = std::max(node->children[0]->upr_bnds[1],node->children[1]->upr_bnds[1]);
      node->upr_bnds[2] = std::max(node->children[0]->upr_bnds[2],node->children[1]->upr_bnds[2]);
      
      // Inflate the bounds to account for FPA errors
      node->lwr_bnds[0] -= 1.0e-6;
      node->lwr_bnds[1] -= 1.0e-6;
      node->lwr_bnds[2] -= 1.0e-6;
      
      node->upr_bnds[0] += 1.0e-6;
      node->upr_bnds[1] += 1.0e-6;
      node->upr_bnds[2] += 1.0e-6;
    }

  }

  void computeProximityBounds( const std::vector<int>& edge_indcs, const std::vector<int>& face_indcs, BVHAABBNode* box )
  {
    assert( box != NULL );

    box->lwr_bnds[0] = box->lwr_bnds[1] = box->lwr_bnds[2] = std::numeric_limits<double>::infinity();  
    box->upr_bnds[0] = box->upr_bnds[1] = box->upr_bnds[2] = -std::numeric_limits<double>::infinity();    
    
    // Compute bounds including all edges in the box
    for( int i = 0; i < (int) edge_indcs.size(); ++i )
    {
      assert( edge_indcs[i] >= 0 ); assert( edge_indcs[i] < (int) m_edges.size() );
      assert( m_edges[edge_indcs[i]].first  >= 0 ); assert( m_edges[edge_indcs[i]].first  < (int) m_vert_radii.size() );
      assert( m_edges[edge_indcs[i]].second >= 0 ); assert( m_edges[edge_indcs[i]].second < (int) m_vert_radii.size() );

      const double r0 = m_vert_radii[m_edges[edge_indcs[i]].first]; assert( r0 >= 0.0 );
      const double r1 = m_vert_radii[m_edges[edge_indcs[i]].second]; assert( r1 >= 0.0 );
      const Vec3d& x0 = m_x.segment<3>(3*m_edges[edge_indcs[i]].first);
      const Vec3d& x1 = m_x.segment<3>(3*m_edges[edge_indcs[i]].second);
      
      // Update lower bounds
      box->lwr_bnds[0] = std::min(box->lwr_bnds[0],x0.x()-r0);
      box->lwr_bnds[1] = std::min(box->lwr_bnds[1],x0.y()-r0);
      box->lwr_bnds[2] = std::min(box->lwr_bnds[2],x0.z()-r0);
      
      box->lwr_bnds[0] = std::min(box->lwr_bnds[0],x1.x()-r1);
      box->lwr_bnds[1] = std::min(box->lwr_bnds[1],x1.y()-r1);
      box->lwr_bnds[2] = std::min(box->lwr_bnds[2],x1.z()-r1);

      // Update upper bounds
      box->upr_bnds[0] = std::max(box->upr_bnds[0],x0.x()+r0);
      box->upr_bnds[1] = std::max(box->upr_bnds[1],x0.y()+r0);
      box->upr_bnds[2] = std::max(box->upr_bnds[2],x0.z()+r0);
      
      box->upr_bnds[0] = std::max(box->upr_bnds[0],x1.x()+r1);
      box->upr_bnds[1] = std::max(box->upr_bnds[1],x1.y()+r1);
      box->upr_bnds[2] = std::max(box->upr_bnds[2],x1.z()+r1);
    }
    
    // Compute bounds including all faces in the box
    for( int i = 0; i < (int) face_indcs.size(); ++i )
    {
      assert( face_indcs[i] >= 0 ); assert( face_indcs[i] < (int) m_faces.size() );
      assert( m_faces[face_indcs[i]].idx[0] >= 0 ); assert( m_faces[face_indcs[i]].idx[0] < (int) m_vert_radii.size() );
      assert( m_faces[face_indcs[i]].idx[1] >= 0 ); assert( m_faces[face_indcs[i]].idx[1] < (int) m_vert_radii.size() );
      assert( m_faces[face_indcs[i]].idx[2] >= 0 ); assert( m_faces[face_indcs[i]].idx[2] < (int) m_vert_radii.size() );
      
      const double r0 = m_vert_radii[m_faces[face_indcs[i]].idx[0]]; assert( r0 >= 0.0 );
      const double r1 = m_vert_radii[m_faces[face_indcs[i]].idx[1]]; assert( r1 >= 0.0 );
      const double r2 = m_vert_radii[m_faces[face_indcs[i]].idx[2]]; assert( r2 >= 0.0 );
      const Vec3d& x0 = m_x.segment<3>(3*m_faces[face_indcs[i]].idx[0]);
      const Vec3d& x1 = m_x.segment<3>(3*m_faces[face_indcs[i]].idx[1]);
      const Vec3d& x2 = m_x.segment<3>(3*m_faces[face_indcs[i]].idx[2]);

      // Update lower bounds
      box->lwr_bnds[0] = std::min(box->lwr_bnds[0],x0.x()-r0);
      box->lwr_bnds[1] = std::min(box->lwr_bnds[1],x0.y()-r0);
      box->lwr_bnds[2] = std::min(box->lwr_bnds[2],x0.z()-r0);
      
      box->lwr_bnds[0] = std::min(box->lwr_bnds[0],x1.x()-r1);
      box->lwr_bnds[1] = std::min(box->lwr_bnds[1],x1.y()-r1);
      box->lwr_bnds[2] = std::min(box->lwr_bnds[2],x1.z()-r1);
      
      box->lwr_bnds[0] = std::min(box->lwr_bnds[0],x2.x()-r2);
      box->lwr_bnds[1] = std::min(box->lwr_bnds[1],x2.y()-r2);
      box->lwr_bnds[2] = std::min(box->lwr_bnds[2],x2.z()-r2);

      // Update lower bounds
      box->upr_bnds[0] = std::max(box->upr_bnds[0],x0.x()+r0);
      box->upr_bnds[1] = std::max(box->upr_bnds[1],x0.y()+r0);
      box->upr_bnds[2] = std::max(box->upr_bnds[2],x0.z()+r0);
      
      box->upr_bnds[0] = std::max(box->upr_bnds[0],x1.x()+r1);
      box->upr_bnds[1] = std::max(box->upr_bnds[1],x1.y()+r1);
      box->upr_bnds[2] = std::max(box->upr_bnds[2],x1.z()+r1);
      
      box->upr_bnds[0] = std::max(box->upr_bnds[0],x2.x()+r2);
      box->upr_bnds[1] = std::max(box->upr_bnds[1],x2.y()+r2);
      box->upr_bnds[2] = std::max(box->upr_bnds[2],x2.z()+r2);
    }

    // Enlarge the bounds to account for FPA errors
    box->lwr_bnds[0] -= 1.0e-6;
    box->lwr_bnds[1] -= 1.0e-6;
    box->lwr_bnds[2] -= 1.0e-6;

    box->upr_bnds[0] += 1.0e-6;
    box->upr_bnds[1] += 1.0e-6;
    box->upr_bnds[2] += 1.0e-6;
  }

  double computeEdgeCentroidX( const int& edg_idx )
  {
    assert( edg_idx >= 0 ); assert( edg_idx < (int) m_edges.size() );
    assert( m_edges[edg_idx].first  >= 0 ); assert( 3*m_edges[edg_idx].first < m_x.size() );
    assert( m_edges[edg_idx].second >= 0 ); assert( 3*m_edges[edg_idx].second < m_x.size() );

    return 0.5*(m_x(3*m_edges[edg_idx].first)+m_x(3*m_edges[edg_idx].second));
  }
  
  double computeFaceCentroidX( const int& fce_idx )
  {
    assert( fce_idx >= 0 ); assert( fce_idx < (int) m_faces.size() );
    assert( m_faces[fce_idx].idx[0] >= 0 ); assert( 3*m_faces[fce_idx].idx[0] < m_x.size() );
    assert( m_faces[fce_idx].idx[1] >= 0 ); assert( 3*m_faces[fce_idx].idx[1] < m_x.size() );
    assert( m_faces[fce_idx].idx[2] >= 0 ); assert( 3*m_faces[fce_idx].idx[2] < m_x.size() );
    
    return (m_x(3*m_faces[fce_idx].idx[0])+m_x(3*m_faces[fce_idx].idx[1])+m_x(3*m_faces[fce_idx].idx[2]))/3.0;
  }
  
  // TODO: In all of these methods (in all of response?) we should precompute the medians before getting proximity stuff
  bool splitOnXMedian( const std::vector<int>& edge_indcs, const std::vector<int>& face_indcs, BVHAABBNode* to_divide )
  {
    assert( to_divide != NULL );

    // Compute the median
    std::vector<double> xmedians;
    for( int i = 0; i < (int) edge_indcs.size(); ++i ) xmedians.push_back(computeEdgeCentroidX(edge_indcs[i]));
    for( int i = 0; i < (int) face_indcs.size(); ++i ) xmedians.push_back(computeFaceCentroidX(face_indcs[i]));
    std::sort(xmedians.begin(),xmedians.end());
    double xmedian;
    if( xmedians.size()%2 == 1 ) xmedian = xmedians[(xmedians.size()-1)/2];
    else xmedian = 0.5*(xmedians[xmedians.size()/2]+xmedians[(xmedians.size()/2)-1]);
    //std::cout << "xmedian: " << xmedian << std::endl;
    
    // Divide the objects to the "left" and to the "right" of the median
    std::vector<int> xlftedgs;
    std::vector<int> xrhtedgs;
    for( int i = 0; i < (int) edge_indcs.size(); ++i ) 
    {
      if( computeEdgeCentroidX(edge_indcs[i]) <= xmedian ) 
      {
        xlftedgs.push_back(edge_indcs[i]);
      }
      else 
      {
        xrhtedgs.push_back(edge_indcs[i]);
      }
    }
    std::vector<int> xlftfaces;
    std::vector<int> xrhtfaces;
    for( int i = 0; i < (int) face_indcs.size(); ++i ) 
    {
      if( computeFaceCentroidX(face_indcs[i]) <= xmedian ) 
      {
        xlftfaces.push_back(face_indcs[i]);
      }
      else 
      {
        xrhtfaces.push_back(face_indcs[i]);
      }
    }
    
    assert( edge_indcs.size() == xlftedgs.size()+xrhtedgs.size() );
    assert( face_indcs.size() == xlftfaces.size()+xrhtfaces.size() );
    
    //std::cout << "xmedian: " << xmedian << std::endl;
    //for( int i = 0; i < (int) edge_indcs.size(); ++i ) std::cout << computeEdgeCentroidX(i) << std::endl;
    //std::cout << "xleftsize: " << xlftedgs.size() << "     xrhtsize: " << xrhtedgs.size() << std::endl;
    
    //#ifdef DEBUG
      std::set<int> testedges; 
      for( int i = 0; i < (int) xlftedgs.size(); ++i ) testedges.insert(xlftedgs[i]);
      for( int i = 0; i < (int) xrhtedgs.size(); ++i ) testedges.insert(xrhtedgs[i]);
      assert( testedges.size() == edge_indcs.size() );

      std::set<int> testfaces;
      for( int i = 0; i < (int) xlftfaces.size(); ++i ) testfaces.insert(xlftfaces[i]);
      for( int i = 0; i < (int) xrhtfaces.size(); ++i ) testfaces.insert(xrhtfaces[i]);
      assert( testfaces.size() == face_indcs.size() );
    //#endif
    
    // If we were not able to split about the median
    if( (xlftedgs.size()+xlftfaces.size()) == 0 || (xrhtedgs.size()+xrhtfaces.size()) == 0 ) return false;

    //std::cout << "Splitting tree into (along x): " << std::endl;
    //std::cout << "    "; for( int i = 0; i < (int) xlftedgs.size(); ++i ) { std::cout << xlftedgs[i] << " "; } std::cout << std::endl;
    //std::cout << "    "; for( int i = 0; i < (int) xrhtedgs.size(); ++i ) { std::cout << xrhtedgs[i] << " "; } std::cout << std::endl;
    
    // Otherwise build children of the node, recurse on children
    to_divide->children[0] = new BVHAABBNode;
    to_divide->children[0]->children[0] = to_divide->children[0]->children[1] = NULL;
    computeProximityBounds( xlftedgs, xlftfaces, to_divide->children[0] );
    divideNode( xlftedgs, xlftfaces, to_divide->children[0] );

    to_divide->children[1] = new BVHAABBNode;
    to_divide->children[1]->children[0] = to_divide->children[1]->children[1] = NULL;
    computeProximityBounds( xrhtedgs, xrhtfaces, to_divide->children[1] );
    divideNode( xrhtedgs, xrhtfaces, to_divide->children[1] );
    
    return true;
  }
  
  double computeEdgeCentroidY( const int& edg_idx )
  {
    assert( edg_idx >= 0 ); assert( edg_idx < (int) m_edges.size() );
    assert( m_edges[edg_idx].first  >= 0 ); assert( 3*m_edges[edg_idx].first+1 < m_x.size() );
    assert( m_edges[edg_idx].second >= 0 ); assert( 3*m_edges[edg_idx].second+1 < m_x.size() );
  
    return 0.5*(m_x(3*m_edges[edg_idx].first+1)+m_x(3*m_edges[edg_idx].second+1));
  }
  
  double computeFaceCentroidY( const int& fce_idx )
  {    
    assert( fce_idx >= 0 ); assert( fce_idx < (int) m_faces.size() );
    assert( m_faces[fce_idx].idx[0] >= 0 ); assert( 3*m_faces[fce_idx].idx[0]+1 < m_x.size() );
    assert( m_faces[fce_idx].idx[1] >= 0 ); assert( 3*m_faces[fce_idx].idx[1]+1 < m_x.size() );
    assert( m_faces[fce_idx].idx[2] >= 0 ); assert( 3*m_faces[fce_idx].idx[2]+1 < m_x.size() );
    
    return (m_x(3*m_faces[fce_idx].idx[0]+1)+m_x(3*m_faces[fce_idx].idx[1]+1)+m_x(3*m_faces[fce_idx].idx[2]+1))/3.0;
  }
  
  bool splitOnYMedian( const std::vector<int>& edge_indcs, const std::vector<int>& face_indcs, BVHAABBNode* to_divide )
  {
    assert( to_divide != NULL );

    // Compute the median
    std::vector<double> ymedians;
    for( int i = 0; i < (int) edge_indcs.size(); ++i ) ymedians.push_back(computeEdgeCentroidY(edge_indcs[i]));
    for( int i = 0; i < (int) face_indcs.size(); ++i ) ymedians.push_back(computeFaceCentroidY(face_indcs[i]));
    std::sort(ymedians.begin(),ymedians.end());
    double ymedian;
    if( ymedians.size()%2 == 1 ) ymedian = ymedians[(ymedians.size()-1)/2];
    else ymedian = 0.5*(ymedians[ymedians.size()/2]+ymedians[(ymedians.size()/2)-1]);
    
    // Divide the objects to the "left" and to the "right" of the median
    std::vector<int> ylftedgs;
    std::vector<int> yrhtedgs;
    for( int i = 0; i < (int) edge_indcs.size(); ++i ) 
    {
      if( computeEdgeCentroidY(edge_indcs[i]) <= ymedian ) 
      {
        ylftedgs.push_back(edge_indcs[i]); 
      }
      else 
      {
        yrhtedgs.push_back(edge_indcs[i]);
      }
    }
    std::vector<int> ylftfaces;
    std::vector<int> yrhtfaces;
    for( int i = 0; i < (int) face_indcs.size(); ++i ) 
    {
      if( computeFaceCentroidY(face_indcs[i]) <= ymedian ) 
      {
        ylftfaces.push_back(face_indcs[i]);
      }
      else 
      {
        yrhtfaces.push_back(face_indcs[i]);
      }
    }
    
    assert( edge_indcs.size() == ylftedgs.size()+yrhtedgs.size() );
    assert( face_indcs.size() == ylftfaces.size()+yrhtfaces.size() );
    
    //#ifdef DEBUG
      std::set<int> testedges; 
      for( int i = 0; i < (int) ylftedgs.size(); ++i ) testedges.insert(ylftedgs[i]);
      for( int i = 0; i < (int) yrhtedgs.size(); ++i ) testedges.insert(yrhtedgs[i]);
      assert( testedges.size() == edge_indcs.size() );
    
      std::set<int> testfaces;
      for( int i = 0; i < (int) ylftfaces.size(); ++i ) testfaces.insert(ylftfaces[i]);
      for( int i = 0; i < (int) yrhtfaces.size(); ++i ) testfaces.insert(yrhtfaces[i]);
      assert( testfaces.size() == face_indcs.size() );
    //#endif

    // If we were not able to split about the median
    if( (ylftedgs.size()+ylftfaces.size()) == 0 || (yrhtedgs.size()+yrhtfaces.size()) == 0 ) return false;

    //std::cout << "Splitting tree into (along y): "  << std::endl;
    //std::cout << "    "; for( int i = 0; i < (int) ylftedgs.size(); ++i ) { std::cout << ylftedgs[i] << " "; } std::cout << std::endl;
    //std::cout << "    "; for( int i = 0; i < (int) yrhtedgs.size(); ++i ) { std::cout << yrhtedgs[i] << " "; } std::cout << std::endl;
    
    // Otherwise build children of the node, recurse on children
    to_divide->children[0] = new BVHAABBNode;
    to_divide->children[0]->children[0] = to_divide->children[0]->children[1] = NULL;
    computeProximityBounds( ylftedgs, ylftfaces, to_divide->children[0] );
    divideNode( ylftedgs, ylftfaces, to_divide->children[0] );

    to_divide->children[1] = new BVHAABBNode;
    to_divide->children[1]->children[0] = to_divide->children[1]->children[1] = NULL;
    computeProximityBounds( yrhtedgs, yrhtfaces, to_divide->children[1] );
    divideNode( yrhtedgs, yrhtfaces, to_divide->children[1] );
    
    return true;
  }
  
  double computeEdgeCentroidZ( const int& edg_idx )
  {
    assert( edg_idx >= 0 ); assert( edg_idx < (int) m_edges.size() );
    assert( m_edges[edg_idx].first  >= 0 ); assert( 3*m_edges[edg_idx].first+2 < m_x.size() );
    assert( m_edges[edg_idx].second >= 0 ); assert( 3*m_edges[edg_idx].second+2 < m_x.size() );

    return 0.5*(m_x(3*m_edges[edg_idx].first+2)+m_x(3*m_edges[edg_idx].second+2));
  }
  
  double computeFaceCentroidZ( const int& fce_idx )
  {
    assert( fce_idx >= 0 ); assert( fce_idx < (int) m_faces.size() );
    assert( m_faces[fce_idx].idx[0] >= 0 ); assert( 3*m_faces[fce_idx].idx[0]+2 < m_x.size() );
    assert( m_faces[fce_idx].idx[1] >= 0 ); assert( 3*m_faces[fce_idx].idx[1]+2 < m_x.size() );
    assert( m_faces[fce_idx].idx[2] >= 0 ); assert( 3*m_faces[fce_idx].idx[2]+2 < m_x.size() );
    
    return (m_x(3*m_faces[fce_idx].idx[0]+2)+m_x(3*m_faces[fce_idx].idx[1]+2)+m_x(3*m_faces[fce_idx].idx[2]+2))/3.0;
  }
  
  bool splitOnZMedian( const std::vector<int>& edge_indcs, const std::vector<int>& face_indcs, BVHAABBNode* to_divide )
  {
    assert( to_divide != NULL );
    
    // Compute the median
    std::vector<double> zmedians;
    for( int i = 0; i < (int) edge_indcs.size(); ++i ) zmedians.push_back(computeEdgeCentroidZ(edge_indcs[i]));
    for( int i = 0; i < (int) face_indcs.size(); ++i ) zmedians.push_back(computeFaceCentroidZ(face_indcs[i]));
    std::sort(zmedians.begin(),zmedians.end());
    double zmedian;
    if( zmedians.size()%2 == 1 ) zmedian = zmedians[(zmedians.size()-1)/2];
    else zmedian = 0.5*(zmedians[zmedians.size()/2]+zmedians[(zmedians.size()/2)-1]);
    
    // Divide the objects to the "left" and to the "right" of the median
    std::vector<int> zlftedgs;
    std::vector<int> zrhtedgs;
    for( int i = 0; i < (int) edge_indcs.size(); ++i ) 
    {
      if( computeEdgeCentroidZ(edge_indcs[i]) <= zmedian ) 
      {
        zlftedgs.push_back(edge_indcs[i]); 
      }
      else 
      {
        zrhtedgs.push_back(edge_indcs[i]);
      }
    }
    std::vector<int> zlftfaces;
    std::vector<int> zrhtfaces;
    for( int i = 0; i < (int) face_indcs.size(); ++i ) 
    {
      if( computeFaceCentroidZ(face_indcs[i]) <= zmedian ) 
      {
        zlftfaces.push_back(face_indcs[i]);
      }
      else 
      {
        zrhtfaces.push_back(face_indcs[i]);
      }
    }
    
    assert( edge_indcs.size() == zlftedgs.size()+zrhtedgs.size() );
    assert( face_indcs.size() == zlftfaces.size()+zrhtfaces.size() );
    
    //#ifdef DEBUG
      std::set<int> testedges; 
      for( int i = 0; i < (int) zlftedgs.size(); ++i ) testedges.insert(zlftedgs[i]);
      for( int i = 0; i < (int) zrhtedgs.size(); ++i ) testedges.insert(zrhtedgs[i]);
      assert( testedges.size() == edge_indcs.size() );
    
      std::set<int> testfaces;
      for( int i = 0; i < (int) zlftfaces.size(); ++i ) testfaces.insert(zlftfaces[i]);
      for( int i = 0; i < (int) zrhtfaces.size(); ++i ) testfaces.insert(zrhtfaces[i]);
      assert( testfaces.size() == face_indcs.size() );
    //#endif

    // If we were not able to split about the median
    if( (zlftedgs.size()+zlftfaces.size()) == 0 || (zrhtedgs.size()+zrhtfaces.size()) == 0 ) return false;
    
    //std::cout << "Splitting tree into (along z): " << std::endl;
    //std::cout << "    "; for( int i = 0; i < (int) zlftedgs.size(); ++i ) { std::cout << zlftedgs[i] << " "; } std::cout << std::endl;
    //std::cout << "    "; for( int i = 0; i < (int) zrhtedgs.size(); ++i ) { std::cout << zrhtedgs[i] << " "; } std::cout << std::endl;
    
    // Otherwise build children of the node, recurse on children
    to_divide->children[0] = new BVHAABBNode;
    to_divide->children[0]->children[0] = to_divide->children[0]->children[1] = NULL;
    computeProximityBounds( zlftedgs, zlftfaces, to_divide->children[0] );
    divideNode( zlftedgs, zlftfaces, to_divide->children[0] );

    to_divide->children[1] = new BVHAABBNode;
    to_divide->children[1]->children[0] = to_divide->children[1]->children[1] = NULL;
    computeProximityBounds( zrhtedgs, zrhtfaces, to_divide->children[1] );
    divideNode( zrhtedgs, zrhtfaces, to_divide->children[1] );
    
    return true;
  }
  
  // Lots of possible cleanups and optimizations here. Precompute the medians, and use O(n) not O(nlgn) median method. 
  void divideNode( const std::vector<int>& edge_indcs, const std::vector<int>& face_indcs, BVHAABBNode* to_divide )
  {
    //std::cout << "Entering divideNode" << std::endl;
    
    assert( to_divide != NULL );
    assert( to_divide->children[0] == NULL && to_divide->children[1] == NULL );

    // If this node only has 1 object, its a leaf. 
    if( edge_indcs.size()+face_indcs.size() == 1 )
    {
      //std::cout << "Base case, 1 object" << std::endl;

      // Record the objects in the node
      //to_divide->edge_indcs = edge_indcs;
      //to_divide->face_indcs = face_indcs;
      to_divide->edge_indcs.insert(to_divide->edge_indcs.begin(),edge_indcs.begin(),edge_indcs.end());
      to_divide->face_indcs.insert(to_divide->face_indcs.begin(),face_indcs.begin(),face_indcs.end());
      //computeProximityBounds(to_divide->edge_indcs,to_divide->face_indcs,to_divide);
      return;
    }

    // Try to split on a cartesian axis, starting with the longest
    std::pair<double,int> xpair = std::pair<double,int>(to_divide->upr_bnds[0]-to_divide->lwr_bnds[0],0);
    std::pair<double,int> ypair = std::pair<double,int>(to_divide->upr_bnds[1]-to_divide->lwr_bnds[1],1);
    std::pair<double,int> zpair = std::pair<double,int>(to_divide->upr_bnds[2]-to_divide->lwr_bnds[2],2);
    std::vector<std::pair<double,int> > axis_extents;
    axis_extents.push_back(xpair);
    axis_extents.push_back(ypair);
    axis_extents.push_back(zpair);
    std::sort(axis_extents.begin(),axis_extents.end());
    for( int i = axis_extents.size()-1; i >= 0; --i )
    {
      if     ( axis_extents[i].second == 0 ) { if( splitOnXMedian( edge_indcs, face_indcs, to_divide ) ) {return;} }
      else if( axis_extents[i].second == 1 ) { if( splitOnYMedian( edge_indcs, face_indcs, to_divide ) ) {return;} }
      else if( axis_extents[i].second == 2 ) { if( splitOnZMedian( edge_indcs, face_indcs, to_divide ) ) {return;} }
      else assert( false );
    }
    
    // Couldn't split on one of the axis - put all of the objects into this node and make it a leaf
    std::cout << "\033[31;1mBVHAABB WARNING:\033[m Failed to split a node, storing all objects into one node." << std::endl;
    // Record the objects in the node
    to_divide->edge_indcs.insert(to_divide->edge_indcs.begin(),edge_indcs.begin(),edge_indcs.end());
    to_divide->face_indcs.insert(to_divide->face_indcs.begin(),face_indcs.begin(),face_indcs.end());
    return;
  }
  
  void recomputeContinuousTimeBoundsBottomUp( BVHAABBNode* node )
  {
    assert( node != NULL );
    
    // If node is a leaf, update the bounding box for that node
    if( isLeaf(node) ) computeContinuousBounds( node->edge_indcs, node->face_indcs, node );
    // Otherwise update the node's children bounding box's and compute based off them
    else 
    {
      assert( hasLeftChild(node) ); 
      assert( hasRightChild(node) );
      
      recomputeContinuousTimeBoundsBottomUp( node->children[0] );
      recomputeContinuousTimeBoundsBottomUp( node->children[1] );
      
      node->lwr_bnds[0] = std::min(node->children[0]->lwr_bnds[0],node->children[1]->lwr_bnds[0]);
      node->lwr_bnds[1] = std::min(node->children[0]->lwr_bnds[1],node->children[1]->lwr_bnds[1]);
      node->lwr_bnds[2] = std::min(node->children[0]->lwr_bnds[2],node->children[1]->lwr_bnds[2]);
      
      node->upr_bnds[0] = std::max(node->children[0]->upr_bnds[0],node->children[1]->upr_bnds[0]);
      node->upr_bnds[1] = std::max(node->children[0]->upr_bnds[1],node->children[1]->upr_bnds[1]);
      node->upr_bnds[2] = std::max(node->children[0]->upr_bnds[2],node->children[1]->upr_bnds[2]);
      
      // Inflate the bounds to account for FPA errors
      node->lwr_bnds[0] -= 1.0e-6;
      node->lwr_bnds[1] -= 1.0e-6;
      node->lwr_bnds[2] -= 1.0e-6;
      
      node->upr_bnds[0] += 1.0e-6;
      node->upr_bnds[1] += 1.0e-6;
      node->upr_bnds[2] += 1.0e-6;
    }
  }
  
  void computeContinuousBounds( const std::vector<int>& edge_indcs, const std::vector<int>& face_indcs, BVHAABBNode* box )
  {
    assert( box != NULL );
  
    // Compute bounds on current timestep position
    computeProximityBounds(edge_indcs,face_indcs,box);

    // Also include next timestep positions when computing bounds
    // Compute bounds including all edges in the box
    for( int i = 0; i < (int) edge_indcs.size(); ++i )
    {
      assert( edge_indcs[i] >= 0 ); assert( edge_indcs[i] < (int) m_edges.size() );
      assert( m_edges[edge_indcs[i]].first  >= 0 ); assert( m_edges[edge_indcs[i]].first  < (int) m_vert_radii.size() );
      assert( m_edges[edge_indcs[i]].second >= 0 ); assert( m_edges[edge_indcs[i]].second < (int) m_vert_radii.size() );

      const double r0 = m_vert_radii[m_edges[edge_indcs[i]].first]; assert( r0 >= 0.0 );
      const double r1 = m_vert_radii[m_edges[edge_indcs[i]].second]; assert( r1 >= 0.0 );
      const Vec3d& x0 = m_x.segment<3>(3*m_edges[edge_indcs[i]].first)  + m_h*m_v.segment<3>(3*m_edges[edge_indcs[i]].first);
      const Vec3d& x1 = m_x.segment<3>(3*m_edges[edge_indcs[i]].second) + m_h*m_v.segment<3>(3*m_edges[edge_indcs[i]].second);

      // Update lower bounds
      box->lwr_bnds[0] = std::min(box->lwr_bnds[0],x0.x()-r0);
      box->lwr_bnds[1] = std::min(box->lwr_bnds[1],x0.y()-r0);
      box->lwr_bnds[2] = std::min(box->lwr_bnds[2],x0.z()-r0);

      box->lwr_bnds[0] = std::min(box->lwr_bnds[0],x1.x()-r1);
      box->lwr_bnds[1] = std::min(box->lwr_bnds[1],x1.y()-r1);
      box->lwr_bnds[2] = std::min(box->lwr_bnds[2],x1.z()-r1);

      // Update upper bounds
      box->upr_bnds[0] = std::max(box->upr_bnds[0],x0.x()+r0);
      box->upr_bnds[1] = std::max(box->upr_bnds[1],x0.y()+r0);
      box->upr_bnds[2] = std::max(box->upr_bnds[2],x0.z()+r0);

      box->upr_bnds[0] = std::max(box->upr_bnds[0],x1.x()+r1);
      box->upr_bnds[1] = std::max(box->upr_bnds[1],x1.y()+r1);
      box->upr_bnds[2] = std::max(box->upr_bnds[2],x1.z()+r1);
    }
    
    // Compute bounds including all faces in the box
    for( int i = 0; i < (int) face_indcs.size(); ++i )
    {
      assert( face_indcs[i] >= 0 ); assert( face_indcs[i] < (int) m_faces.size() );
      assert( m_faces[face_indcs[i]].idx[0] >= 0 ); assert( m_faces[face_indcs[i]].idx[0] < (int) m_vert_radii.size() );
      assert( m_faces[face_indcs[i]].idx[1] >= 0 ); assert( m_faces[face_indcs[i]].idx[1] < (int) m_vert_radii.size() );
      assert( m_faces[face_indcs[i]].idx[2] >= 0 ); assert( m_faces[face_indcs[i]].idx[2] < (int) m_vert_radii.size() );
      
      const double r0 = m_vert_radii[m_faces[face_indcs[i]].idx[0]]; assert( r0 >= 0.0 );
      const double r1 = m_vert_radii[m_faces[face_indcs[i]].idx[1]]; assert( r1 >= 0.0 );
      const double r2 = m_vert_radii[m_faces[face_indcs[i]].idx[2]]; assert( r2 >= 0.0 );
      const Vec3d& x0 = m_x.segment<3>(3*m_faces[face_indcs[i]].idx[0]) + m_h*m_v.segment<3>(3*m_faces[face_indcs[i]].idx[0]);
      const Vec3d& x1 = m_x.segment<3>(3*m_faces[face_indcs[i]].idx[1]) + m_h*m_v.segment<3>(3*m_faces[face_indcs[i]].idx[1]);
      const Vec3d& x2 = m_x.segment<3>(3*m_faces[face_indcs[i]].idx[2]) + m_h*m_v.segment<3>(3*m_faces[face_indcs[i]].idx[2]);
      
      // Update lower bounds
      box->lwr_bnds[0] = std::min(box->lwr_bnds[0],x0.x()-r0);
      box->lwr_bnds[1] = std::min(box->lwr_bnds[1],x0.y()-r0);
      box->lwr_bnds[2] = std::min(box->lwr_bnds[2],x0.z()-r0);
      
      box->lwr_bnds[0] = std::min(box->lwr_bnds[0],x1.x()-r1);
      box->lwr_bnds[1] = std::min(box->lwr_bnds[1],x1.y()-r1);
      box->lwr_bnds[2] = std::min(box->lwr_bnds[2],x1.z()-r1);
      
      box->lwr_bnds[0] = std::min(box->lwr_bnds[0],x2.x()-r2);
      box->lwr_bnds[1] = std::min(box->lwr_bnds[1],x2.y()-r2);
      box->lwr_bnds[2] = std::min(box->lwr_bnds[2],x2.z()-r2);
      
      // Update lower bounds
      box->upr_bnds[0] = std::max(box->upr_bnds[0],x0.x()+r0);
      box->upr_bnds[1] = std::max(box->upr_bnds[1],x0.y()+r0);
      box->upr_bnds[2] = std::max(box->upr_bnds[2],x0.z()+r0);
      
      box->upr_bnds[0] = std::max(box->upr_bnds[0],x1.x()+r1);
      box->upr_bnds[1] = std::max(box->upr_bnds[1],x1.y()+r1);
      box->upr_bnds[2] = std::max(box->upr_bnds[2],x1.z()+r1);
      
      box->upr_bnds[0] = std::max(box->upr_bnds[0],x2.x()+r2);
      box->upr_bnds[1] = std::max(box->upr_bnds[1],x2.y()+r2);
      box->upr_bnds[2] = std::max(box->upr_bnds[2],x2.z()+r2);
    }
    
    // Enlarge the bounds to account for FPA errors
    box->lwr_bnds[0] -= 1.0e-6;
    box->lwr_bnds[1] -= 1.0e-6;
    box->lwr_bnds[2] -= 1.0e-6;
    
    box->upr_bnds[0] += 1.0e-6;
    box->upr_bnds[1] += 1.0e-6;
    box->upr_bnds[2] += 1.0e-6;
  }

  BVHAABBNode* m_root_node;
  const VecXd& m_x;
  const VecXd& m_v;
  const std::vector<std::pair<int,int> >& m_edges;
  const std::vector<TriangularFace>& m_faces;
  const std::vector<double>& m_vert_radii;
  double m_h;
};

}

#endif

// TODO:
//   o Pull code out of header file :)
//   o A number of speedups possible in constructing stuff by precomputing medians, etc
//   o Implement face-face in case we decide to do surface-surface collisions

