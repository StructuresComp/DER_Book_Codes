#include "RodBoundaryCondition.hh"
#include "BASim/src/Physics/ElasticRods/ElasticRod.hh"

namespace BASim {
    
    RodBoundaryCondition::RodBoundaryCondition(ElasticRod& rod)
    : m_rod(rod)
    {
        m_rod.add_property(m_scriptedVerts, "list of scripted vertices");
        m_rod.add_property(m_desiredPositions, "desired vertex positions", Vec3d::Zero().eval());
        m_rod.add_property(m_isVertexScripted, "is vertex scripted", false);
        
        m_rod.add_property(m_scriptedEdges, "list of scripted edges");
        m_rod.add_property(m_desiredTheta, "desired theta values", 0.0);
        m_rod.add_property(m_isMaterialScripted, "is material scripted", false);
    }
    
    const RodBoundaryCondition::BCList&
    RodBoundaryCondition::scriptedVertices() const
    {
        return m_rod.property(m_scriptedVerts);
    }
    
    bool RodBoundaryCondition::isVertexScripted(int vertIdx) const
    {
        return m_rod.property(m_isVertexScripted)[vertIdx];
    }
    
    void RodBoundaryCondition::setDesiredVertexPosition(int vertIdx,
                                                        const Vec3d& position)
    {
        assert(vertIdx >= 0);
        assert(vertIdx < m_rod.nv());
        
        m_rod.property(m_isVertexScripted)[vertIdx] = true;
        
        BCList& verts = m_rod.property(m_scriptedVerts);
        BCList::iterator result = std::find(verts.begin(), verts.end(), vertIdx);
        if (result == verts.end()) verts.push_back(vertIdx);
        m_rod.property(m_desiredPositions)[vertIdx] = position;
    }
    
    const Vec3d& RodBoundaryCondition::getDesiredVertexPosition(int vertIdx)
    {
        return m_rod.property(m_desiredPositions)[vertIdx];
    }
    
    void RodBoundaryCondition::releaseVertex(int vertIdx)
    {
        assert(vertIdx >= 0);
        assert(vertIdx < m_rod.nv());
        
        m_rod.property(m_isVertexScripted)[vertIdx] = false;
        
        BCList& verts = m_rod.property(m_scriptedVerts);
        
        //  std::remove(verts.begin(), verts.end(), vertIdx);   // this doesn't change the size of vector.
        for(int i=0; i<(int)verts.size(); i++) {
            if (verts[i] == vertIdx) {
                verts.erase(verts.begin() + i);
                break;
            }
        }
        
    }
    
    const RodBoundaryCondition::BCList& RodBoundaryCondition::scriptedEdges() const
    {
        return m_rod.property(m_scriptedEdges);
    }
    
    bool RodBoundaryCondition::isEdgeScripted(int edgeIdx) const
    {
        return m_rod.property(m_isMaterialScripted)[edgeIdx];
    }
    
    void RodBoundaryCondition::setDesiredEdgeAngle(int edgeIdx, const Scalar& theta)
    {
        assert(edgeIdx >= 0);
        assert(edgeIdx < m_rod.ne());
        
        m_rod.property(m_isMaterialScripted)[edgeIdx] = true;
        
        BCList& edges = m_rod.property(m_scriptedEdges);
        BCList::iterator result
        = std::find(edges.begin(), edges.end(), edgeIdx);
        if (result == edges.end()) edges.push_back(edgeIdx);
        
        m_rod.property(m_desiredTheta)[edgeIdx] = theta;
    }
    
    const Scalar& RodBoundaryCondition::getDesiredEdgeAngle(int edgeIdx)
    {
        return m_rod.property(m_desiredTheta)[edgeIdx];
    }
    
    void RodBoundaryCondition::releaseEdge(int edgeIdx)
    {
        assert(edgeIdx >= 0);
        assert(edgeIdx < m_rod.ne());
        
        m_rod.property(m_isMaterialScripted)[edgeIdx] = false;
        
        BCList& edges = m_rod.property(m_scriptedEdges);
        //  std::remove(edges.begin(), edges.end(), edgeIdx);	// this doesn't change the size of vector.
        
        for(int i=0; i<(int)edges.size(); i++) {
            if (edges[i] == edgeIdx) {
                edges.erase(edges.begin() + i);
                break;
            }
        }
    }
    
    void RodBoundaryCondition::adjustListAfterDeleteVertex(int vertIdx){
        BCList& verts = m_rod.property(m_scriptedVerts);
        
        for(int i=0; i<(int)verts.size(); i++) {
            if (verts[i] == vertIdx) {
                verts.erase(verts.begin() + i);
                i--;
            } else if (verts[i] > vertIdx) {
                verts[i]--;
            }
        }
        
    }
    
    void RodBoundaryCondition::adjustListAfterDeleteEdge(int edgeIdx) {
        BCList& edges = m_rod.property(m_scriptedEdges);
        
        for(int i=0; i<(int)edges.size(); i++) {
            if (edges[i] == edgeIdx) {
                edges.erase(edges.begin() + i);
                i--;
            } else if (edges[i] > edgeIdx) {
                edges[i]--;
            }
        }
        
    }
    
    
    
} // namespace BASim
