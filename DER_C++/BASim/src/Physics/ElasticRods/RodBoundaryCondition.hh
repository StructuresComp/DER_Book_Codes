/**
 * \file RodBoundaryCondition.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 12/09/2009
 */

#ifndef RODBOUNDARYCONDITION_HH
#define RODBOUNDARYCONDITION_HH

#include "BASim/src/Core/ObjectHandles.hh"
#include "BASim/src/Core/TopologicalObject/TopObjHandles.hh"

namespace BASim {
    
    class ElasticRod;
    
    /** Class for managing fixed/scripted vertices and edges of a rod */
    // TODO: For consistency, methods in RodBoundaryCondition should acccept
    //       edge and/or vertex iterators to be consistent with the rest of
    //       BASim.
    class RodBoundaryCondition
    {
    public:
        
        typedef std::vector<int> BCList;
        
        RodBoundaryCondition(ElasticRod& rod);
        
        const BCList& scriptedVertices() const;
        bool isVertexScripted(int vertIdx) const;
        void setDesiredVertexPosition(int vertIdx, const Vec3d& position);
        const Vec3d& getDesiredVertexPosition(int vertIdx);
        void releaseVertex(int vertIdx);
        
        void adjustListAfterDeleteVertex(int vertIdx);
        
        const BCList& scriptedEdges() const;
        bool isEdgeScripted(int edgeIdx) const;
        void setDesiredEdgeAngle(int edgeIdx, const Scalar& theta);
        const Scalar& getDesiredEdgeAngle(int edgeIdx);
        void releaseEdge(int edgeIdx);
        
        void adjustListAfterDeleteEdge(int edgeIdx);
        
        void setVerticalConstraint(int vertex_id, Scalar desired_pos, Scalar desired_vel, Scalar uncons_pos, bool is_static)
        {
            scriptedVerticalVertices.push_back(vertex_id);
            desiredVerticalPositions.push_back(desired_pos);
            // desiredVerticalVelocities.push_back(desired_vel);
            // unconstrainedVerticalPositions.push_back(uncons_pos);
            // isStaticMode.push_back(is_static);
        }
        
        // For CNT project (Khalid)
        // Vertical (Y) constraint (The following 5 functions)
        void clearVerticalConstraints()
        {
            scriptedVerticalVertices.clear();
            desiredVerticalPositions.clear();
            // desiredVerticalVelocities.clear();
            // unconstrainedVerticalPositions.clear();
            // isStaticMode.clear();
        }
        int getNumberOfVerticallyScriptedVertices() const {
            return (int)scriptedVerticalVertices.size();
        }
        void setDesiredVerticalPosition(int vertex_id, Scalar desired_pos)
        {
            scriptedVerticalVertices.push_back(vertex_id);
            desiredVerticalPositions.push_back(desired_pos);
        }
        Scalar getDesiredVerticalPosition(int id) const
        {
            return desiredVerticalPositions[id];
        }
        int getVerticalScriptedVertId(int id) const
        {
            // 'id' means the order in the vector of scripted vertices
            return scriptedVerticalVertices[id];
        }
        // Horizontal (X) constraint (The following 5 functions)
        void clearHorizontalXConstraints()
        {
            scriptedHorizontalXVertices.clear();
            desiredHorizontalXPositions.clear();
        }
        int getNumberOfHorizontallyXScriptedVertices() const {
            return (int)scriptedHorizontalXVertices.size();
        }
        void setDesiredHorizontalXPosition(int vertex_id, Scalar desired_pos)
        {
            scriptedHorizontalXVertices.push_back(vertex_id);
            desiredHorizontalXPositions.push_back(desired_pos);
        }
        Scalar getDesiredHorizontalXPosition(int id) const
        {
            return desiredHorizontalXPositions[id];
        }
        int getHorizontalXScriptedVertId(int id) const
        {
            return scriptedHorizontalXVertices[id];
        }
        // Horizontal (Z) constraint (The following 5 functions)
        void clearHorizontalZConstraints()
        {
            scriptedHorizontalZVertices.clear();
            desiredHorizontalZPositions.clear();
        }
        int getNumberOfHorizontallyZScriptedVertices() const {
            return (int)scriptedHorizontalZVertices.size();
        }
        void setDesiredHorizontalZPosition(int vertex_id, Scalar desired_pos)
        {
            scriptedHorizontalZVertices.push_back(vertex_id);
            desiredHorizontalZPositions.push_back(desired_pos);
        }
        Scalar getDesiredHorizontalZPosition(int id) const
        {
            return desiredHorizontalZPositions[id];
        }
        int getHorizontalZScriptedVertId(int id) const
        {
            return scriptedHorizontalZVertices[id];
        }
        
        Scalar getDesiredVerticalVelocity(int id) const
        {
            return desiredVerticalVelocities[id];
        }
        
        Scalar getUnconstrainedVerticalPosition(int id) const
        {
            return unconstrainedVerticalPositions[id];
        }
        
        bool isStatic(int id) const
        {
            return isStaticMode[id];
        }
        
        void setStatic(int id, bool is_static)
        {
            isStaticMode[id] = is_static;
        }
        
        
    protected:
        
        ElasticRod& m_rod;
        
        ObjPropHandle<BCList> m_scriptedVerts;
        VPropHandle<Vec3d> m_desiredPositions;
        VPropHandle<bool> m_isVertexScripted;
        
        ObjPropHandle<BCList> m_scriptedEdges;
        EPropHandle<Scalar> m_desiredTheta;
        EPropHandle<bool> m_isMaterialScripted;
        
        std::vector<Scalar> desiredVerticalVelocities;
        std::vector<Scalar> unconstrainedVerticalPositions;
        std::vector<bool> isStaticMode;
        
        // Constrain node DOF only x,y or z direction
        std::vector<int> scriptedVerticalVertices;
        std::vector<Scalar> desiredVerticalPositions;
        std::vector<int> scriptedHorizontalXVertices;
        std::vector<Scalar> desiredHorizontalXPositions;      
        std::vector<int> scriptedHorizontalZVertices;
        std::vector<Scalar> desiredHorizontalZPositions;
    };
    
} // namespace BASim

#endif // RODBOUNDARYCONDITION_HH
