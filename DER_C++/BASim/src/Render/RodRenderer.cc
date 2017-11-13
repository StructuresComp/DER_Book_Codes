/**
 * \file RodRenderer.cc
 *
 * \author miklos@cs.columbia.edu
 * \date 08/30/2009
 */

#include "BASim/src/Render/RodRenderer.hh"
#include "BASim/src/Render/RenderUtils.hh"
#include "TriangleMeshRenderer.hh"

namespace BASim {

RodRenderer * RodRenderer::s_last_instance = NULL;
    
RodRenderer::RodRenderer(ElasticRod& rod)
  : m_rod(rod)
  , m_tube(m_rod)
  , m_mode(SMOOTH)
  , m_drawMaterial(false)
  , m_drawReference(false)
  , m_scaleToRadius(false)
  , m_drawArrows(false)
{
  // material colors
  m_palette.push_back(Color(107,144,212));
  m_palette.push_back(Color(16,71,169));
  //m_palette.push_back(Color(255, 125, 75));
  //m_palette.push_back(Color(55, 125, 255));

  // reference frame colors
  m_palette.push_back(Color(200, 200, 200));
  m_palette.push_back(Color(50, 50, 50));

  // curvature binormal color
  m_palette.push_back(Color(102, 204, 51));
    
  s_last_instance = this;
}

void RodRenderer::render()
{
  if (m_mode == SMOOTH) drawSmoothRod();
  else if (m_mode == SIMPLE) drawSimpleRod();
  else if (m_mode == PHOTOREALISTIC) drawPRRod();
    
  if (m_mode != PHOTOREALISTIC)
  {
      if (m_drawMaterial) drawMaterialFrame();
      if (m_drawReference) drawReferenceFrame();
  }
}

void RodRenderer::drawSimpleRod()
{
  glDisable(GL_LIGHTING);
  
  glLineWidth(2);
  glBegin(GL_LINES);

  const Color& edgeColor = m_palette[1];
  const Color& inverseEdgeColor = edgeColor.inverse();

  int i;
  ElasticRod::edge_iter eit;
  for (eit = m_rod.edges_begin(), i = 0; eit != m_rod.edges_end(); ++eit, ++i ) 
  {
    if( m_rod.getBoundaryCondition()->isEdgeScripted(i) ) OpenGL::color(inverseEdgeColor);
    else OpenGL::color(edgeColor);

    ElasticRod::EdgeVertexIter evit = m_rod.ev_iter(*eit);
    for (evit = m_rod.ev_iter(*eit); evit; ++evit) {
      Vec3d x = m_rod.getVertex(*evit);
      OpenGL::vertex(x);
    }
  }

  glEnd();

  glPointSize(5);
  glBegin(GL_POINTS);

  ElasticRod::vertex_iter vit;
  for (vit = m_rod.vertices_begin(), i = 0; vit != m_rod.vertices_end(); ++vit, ++i ) 
  {
    if( m_rod.getBoundaryCondition()->isVertexScripted(i) ) OpenGL::color(inverseEdgeColor);
    else OpenGL::color(edgeColor);

    const Vec3d& x = m_rod.getVertex(*vit);
    OpenGL::vertex(x);
  }

  glEnd();
  
  
  glEnable(GL_LIGHTING);
}

void RodRenderer::drawSmoothRod()
{
  m_tube.buildTube();

  glEnable(GL_LIGHTING);
  glEnable(GL_COLOR_MATERIAL);
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

  const Color& color1 = m_palette[0];
  const Color& color2 = m_palette[1];

  int last_vertex_fixed = -1;
  for(int i=m_rod.nv() - 5; i>= 0; i--) {
    if (m_rod.getBoundaryCondition()->isVertexScripted(i)) {          
      last_vertex_fixed = i; break;
    }
  }
  
  

  int slices = m_tube.getSlices();
  for (int j = 0; j < m_rod.ne(); ++j) {
    const Vec3d& v0 = m_rod.getVertex(j);
    const Vec3d& v1 = m_rod.getVertex((j + 1) % m_rod.nv());
    glBegin(GL_QUAD_STRIP);
    
    double twist = (m_rod.m_twisting_energy(j) + m_rod.m_twisting_energy(j+1)) / 2.0;
    twist *= 5.0;
    
    Color color1t, color2t;
    
    if (twist > 0) {
      if (twist > 1) twist = 1.0;
      color1t = Color(0, 0, (int)(255 * twist));
      color2t = Color(0, 0, (int)(255 * twist));
    } else {
      if (twist < -1.0) twist = -1.0;
      color1t = Color((int)(255 * -twist),0,0);
      color2t = Color((int)(255 * -twist),0,0);
    }
    
    /*
    double twist_level = (m_rod.m_twisting_energy(j) + m_rod.m_twisting_energy(j+1)) / 5.0;
    if (twist_level > 1.0) twist_level=1.0;
    
    Color color1t = Color(0, 0, (int)(255 * twist_level));
    Color color2t = Color(0, 0, (int)(200 * twist_level));
    */
    
    
    for (int k = 0; k <= slices; ++k) {
      if ((0 <= k && k < (slices + 3) / 4)
          || (slices / 2 <= k && k < (3 * slices + 3) / 4)) {
        OpenGL::color(color1t);
      } else {
        OpenGL::color(color2t);
      }
      
      if (j + 1 <= last_vertex_fixed)
      {
        glColor3f(1, 0, 0);
      }

      const Vec3d& x0 = m_tube.getPointAtVert(j, k);
      Vec3d n = (x0 - v0).normalized();
      OpenGL::normal(n);
      OpenGL::vertex(x0);

      const Vec3d& x1 = m_tube.getPointAtVert(j + 1, k);
      n = (x1 - v1).normalized();
      OpenGL::normal(n);
      OpenGL::vertex(x1);
    }

    glEnd();
  }

  glDisable(GL_COLOR_MATERIAL);
  glDisable(GL_LIGHTING);
}

void RodRenderer::drawMaterialFrame()
{
  glDisable(GL_LIGHTING);
  
  glLineWidth(2);
  glBegin(GL_LINES);

  const Color& color1 = m_palette[0];
  const Color& color2 = m_palette[1];
  Scalar r = 1;

  Color colorr = Color(0, 100, 200);
      
  {
    std::vector<Vec3d> pts;
    
    ElasticRod::edge_iter eit, end = m_rod.edges_end();
    for (eit = m_rod.edges_begin(); eit != end; ++eit) {
      ElasticRod::edge_handle& eh = *eit;
      ElasticRod::vertex_handle vh0 = m_rod.fromVertex(eh);
      ElasticRod::vertex_handle vh1 = m_rod.toVertex(eh);
  
      Vec3d x = (m_rod.getVertex(vh0) + m_rod.getVertex(vh1)) / 2.0;
      if (m_scaleToRadius) r = m_rod.radiusA(eh) * 1.5;
  
      OpenGL::color(color1);
      Vec3d y = x + r * m_rod.getMaterial1(eh);
      OpenGL::vertex(x);
      OpenGL::vertex(y);
      
      pts.push_back(y);
    }
    
//    glLineWidth(3);
    OpenGL::color(colorr);
    
    for(int i=0; i<pts.size() - 1; i++) {
            
      OpenGL::vertex(pts[i]);
      OpenGL::vertex(pts[i+1]);
  
    }
    
    glEnd();
  
    glEnable(GL_LIGHTING);
    return;
  } 
  
  ElasticRod::edge_iter eit, end = m_rod.edges_end();
  for (eit = m_rod.edges_begin(); eit != end; ++eit) {
    ElasticRod::edge_handle& eh = *eit;
    ElasticRod::vertex_handle vh0 = m_rod.fromVertex(eh);
    ElasticRod::vertex_handle vh1 = m_rod.toVertex(eh);

    Vec3d x = (m_rod.getVertex(vh0) + m_rod.getVertex(vh1)) / 2.0;
    if (m_scaleToRadius) r = m_rod.radiusA(eh);

    OpenGL::color(color1);
    if (m_drawArrows) drawArrow(x, m_rod.getMaterial1(eh), r);
    else {
      Vec3d y = x + r * m_rod.getMaterial1(eh);
      OpenGL::vertex(x);
      OpenGL::vertex(y);
    }

    if (m_scaleToRadius) r = m_rod.radiusB(eh);

    OpenGL::color(color2);
    if (m_drawArrows) drawArrow(x, m_rod.getMaterial2(eh), r);
    else {
      Vec3d y = x + r * m_rod.getMaterial2(eh);
      OpenGL::vertex(x);
      OpenGL::vertex(y);
    }
  }

  glEnd();
  
  glEnable(GL_LIGHTING);
}

void RodRenderer::drawReferenceFrame()
{
  glDisable(GL_LIGHTING);

  glLineWidth(2);
  glBegin(GL_LINES);

  const Color& color1 = m_palette[2];
  const Color& color2 = m_palette[3];

  ElasticRod::edge_iter eit, end = m_rod.edges_end();
  for (eit = m_rod.edges_begin(); eit != end; ++eit) {
    ElasticRod::edge_handle& eh = *eit;
    ElasticRod::vertex_handle vh0 = m_rod.fromVertex(eh);
    ElasticRod::vertex_handle vh1 = m_rod.toVertex(eh);

    Vec3d x = (m_rod.getVertex(vh0) + m_rod.getVertex(vh1)) / 2.0;

    OpenGL::color(color1);
    if (m_drawArrows) drawArrow(x, m_rod.getReferenceDirector1(eh));
    else {
      Vec3d y = x + m_rod.getReferenceDirector1(eh);
      OpenGL::vertex(x);
      OpenGL::vertex(y);
    }

    OpenGL::color(color2);
    if (m_drawArrows) drawArrow(x, m_rod.getReferenceDirector2(eh));
    else {
      Vec3d y = x + m_rod.getReferenceDirector2(eh);
      OpenGL::vertex(x);
      OpenGL::vertex(y);
    }
  }

  glEnd();

  glEnable(GL_LIGHTING);
}

Vec3d RodRenderer::calculateObjectCenter()
{
  Vec3d center = Vec3d::Zero();
  
  ElasticRod::vertex_iter vit, end = m_rod.vertices_end();
  for (vit = m_rod.vertices_begin(); vit != end; ++vit) {
    center += m_rod.getVertex(*vit);
  }
  
  center /= m_rod.nv();
  
  return center;
}

Scalar RodRenderer::calculateObjectBoundingRadius(const Vec3d& center)
{
  Scalar radius = 0.0;
  
  ElasticRod::vertex_iter vit, end = m_rod.vertices_end();
  for (vit = m_rod.vertices_begin(); vit != end; ++vit) {
    radius = std::max(radius, (m_rod.getVertex(*vit) - center).norm());
  }
  
  return radius;
}
    
    class PerlinNoiseGenerator
    {
        // reference: http://mrl.nyu.edu/~perlin/doc/oscar.html#noise
    public:
        static const int B = 0x100;
        static const int BM = 0xff;
        
        static const int N = 0x1000;
        static const int NP = 12;   // 2^N
        static const int NM = 0xfff;
        
        static int p[B + B + 2];
        static float g2[B + B + 2][2];
        static int start;
        
#define s_curve(t) ( t * t * (3. - 2. * t) )
        
#define lerp(t, a, b) ( a + t * (b - a) )
        
#define setup(i,b0,b1,r0,r1)\
t = vec[i] + N;\
b0 = ((int)t) & BM;\
b1 = (b0+1) & BM;\
r0 = t - (int)t;\
r1 = r0 - 1.;
        
        static float noise2(float vec[2])
        {
            int bx0, bx1, by0, by1, b00, b10, b01, b11;
            float rx0, rx1, ry0, ry1, *q, sx, sy, a, b, t, u, v;
            int i, j;
            
            if (start) {
                start = 0;
                init();
            }
            
            setup(0, bx0,bx1, rx0,rx1);
            setup(1, by0,by1, ry0,ry1);
            
            i = p[ bx0 ];
            j = p[ bx1 ];
            
            b00 = p[ i + by0 ];
            b10 = p[ j + by0 ];
            b01 = p[ i + by1 ];
            b11 = p[ j + by1 ];
            
            sx = s_curve(rx0);
            sy = s_curve(ry0);
            
#define at2(rx,ry) ( rx * q[0] + ry * q[1] )
            
            q = g2[ b00 ] ; u = at2(rx0,ry0);
            q = g2[ b10 ] ; v = at2(rx1,ry0);
            a = lerp(sx, u, v);
            
            q = g2[ b01 ] ; u = at2(rx0,ry1);
            q = g2[ b11 ] ; v = at2(rx1,ry1);
            b = lerp(sx, u, v);
            
            return lerp(sy, a, b);
        }
        
        static void normalize2(float v[2])
        {
            float s;
            
            s = sqrt(v[0] * v[0] + v[1] * v[1]);
            v[0] = v[0] / s;
            v[1] = v[1] / s;
        }
        
        static void init(void)
        {
            int i, j, k;
            
            for (i = 0 ; i < B ; i++) {
                p[i] = i;
                
                for (j = 0 ; j < 2 ; j++)
                    g2[i][j] = (float)((random() % (B + B)) - B) / B;
                normalize2(g2[i]);
            }
            
            while (--i) {
                k = p[i];
                p[i] = p[j = random() % B];
                p[j] = k;
            }
            
            for (i = 0 ; i < B + 2 ; i++) {
                p[B + i] = p[i];
                for (j = 0 ; j < 2 ; j++)
                    g2[B + i][j] = g2[i][j];
            }
        }
    };
    
    int PerlinNoiseGenerator::p[PerlinNoiseGenerator::B + PerlinNoiseGenerator::B + 2];
    float PerlinNoiseGenerator::g2[PerlinNoiseGenerator::B + PerlinNoiseGenerator::B + 2][2];
    int PerlinNoiseGenerator::start = 1;
    
void RodRenderer::drawPRRod()
{/*
    m_tube.buildTube();

    static GLuint fboBelt = 0;
    static GLuint texBelt = 0;
    static GLuint fboReflection = 0;
    static GLuint texReflection = 0;
    static GLuint texWall = 0;
    static int window_width = 0;
    static int window_height = 0;
    int ww = glutGet(GLUT_WINDOW_WIDTH);
    int wh = glutGet(GLUT_WINDOW_HEIGHT);
    if (ww != window_width || wh != window_height)
    {
        window_width = ww;
        window_height = wh;
        
        // create textures to render to
        glEnable(GL_TEXTURE_2D);
        
        glGenTextures(1, &texBelt);
        glBindTexture(GL_TEXTURE_2D, texBelt);
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, window_width, window_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);

        glGenTextures(1, &texReflection);
        glBindTexture(GL_TEXTURE_2D, texReflection);
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, window_width, window_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);
        
        glBindTexture(GL_TEXTURE_2D, 0);
        glDisable(GL_TEXTURE_2D);

        // render buffer for the depth attachment
        GLuint rboBeltDepth;
        glGenRenderbuffers(1, &rboBeltDepth);
        glBindRenderbuffer(GL_RENDERBUFFER, rboBeltDepth);
        glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, window_width, window_height);

        GLuint rboReflectionDepth;
        glGenRenderbuffers(1, &rboReflectionDepth);
        glBindRenderbuffer(GL_RENDERBUFFER, rboReflectionDepth);
        glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, window_width, window_height);
        
        glBindRenderbuffer(GL_RENDERBUFFER, 0);
        
        // create the FBO and assign the attachments
        glGenFramebuffers(1, &fboBelt);
        glBindFramebuffer(GL_FRAMEBUFFER, fboBelt);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texBelt, 0);    
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rboBeltDepth);
        
        glGenFramebuffers(1, &fboReflection);
        glBindFramebuffer(GL_FRAMEBUFFER, fboReflection);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texReflection, 0);    
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rboReflectionDepth);
        
        // check FBO status
        assert(glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE);
        assert(glGetError() == GL_NO_ERROR);
        
        // create a noise texture
        glGenTextures(1, &texWall);
        glBindTexture(GL_TEXTURE_2D, texWall);
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

        int N = 1024;
        char buffer[N * N * 4];
        float st1[2];
        float st2[2];
        float st3[2];
        float st4[2];
        float mul = 15;
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
            {
                st1[0] = (float)i / (N - 1) * 2 * mul;
                st1[1] = (float)j / (N - 1) * 2 * mul;
                st2[0] = (float)i / (N - 1) * 4 * mul;
                st2[1] = (float)j / (N - 1) * 4 * mul;
                st3[0] = (float)i / (N - 1) * 8 * mul;
                st3[1] = (float)j / (N - 1) * 8 * mul;
                st4[0] = (float)i / (N - 1) * 16 * mul;
                st4[1] = (float)j / (N - 1) * 16 * mul;
                float v = PerlinNoiseGenerator::noise2(st1) + PerlinNoiseGenerator::noise2(st2) + PerlinNoiseGenerator::noise2(st3) + PerlinNoiseGenerator::noise2(st4);
                char * p = buffer + (i * 512 + j) * 4;
                *(p + 0) = (int)(255 * (0.35 + v * 0.02));
                *(p + 1) = (int)(255 * (0.4 + v * 0.02));
                *(p + 2) = (int)(255 * (0.4 + v * 0.02));
                *(p + 3) = (int)(255);
            }
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, 512, 512, 0, GL_RGBA, GL_UNSIGNED_BYTE, buffer);
    }
    
    glShadeModel(GL_SMOOTH);
    glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
    
    glEnable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    
    // render the rod reflection into the belt texture (with shading)
    glBindFramebuffer(GL_FRAMEBUFFER, fboReflection);
    glBindTexture(GL_TEXTURE_2D, 0);

    glClearColor(0.3, 0.3, 0.3, 0);
    glClearDepth(1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    for (int j = 0; j < m_rod.ne(); j++) 
    {
        const Vec3d & v0 = m_rod.getVertex(j);
        const Vec3d & v1 = m_rod.getVertex((j + 1) % m_rod.nv());
        
        Color color1t, color2t;        

        color1t = Color(0.8, 1.0, 0.5, 1.0);
        color2t = Color(0.8, 1.0, 0.5, 1.0);
        
        glBegin(GL_QUAD_STRIP);
        for (int k = 0; k <= m_tube.getSlices(); k++) 
        {
            OpenGL::color(color1t);
            
            Vec3d x0 = m_tube.getPointAtVert(j, k);
            Vec3d n = (x0 - v0).normalized();
            x0.y() = -x0.y();
            n.y() = -n.y();
            OpenGL::normal(n);
            OpenGL::vertex(x0);
            
            Vec3d x1 = m_tube.getPointAtVert(j + 1, k);
            n = (x1 - v1).normalized();
            x1.y() = -x1.y();
            n.y() = -n.y();
            OpenGL::normal(n);
            OpenGL::vertex(x1);
        }
        glEnd();
    }
    
    // render the environment map into the belt texture
    ////
    
    // render the belt markers
    Scalar r = -12;
    Scalar l = 12;
    Scalar tmin = TriangleMeshRenderer::lastInstance()->mesh().getVertex(*(TriangleMeshRenderer::lastInstance()->mesh().vertices_begin())).x();
    Scalar tmax = tmin + 2000;

    glDisable(GL_LIGHTING);
    glDisable(GL_COLOR_MATERIAL);

//    glLineWidth(1);
//    glColor4f(0, 0, 0, 1);
//    glBegin(GL_LINES);
//    for (Scalar t = tmin; t < tmax; t += 5.0)
//    {
//        glVertex3f(t, 0.0, r);
//        glVertex3f(t, 0.0, r + 5.0);
//        
//        glVertex3f(t, 0.0, l - 5.0);
//        glVertex3f(t, 0.0, l);
//    }
//    glEnd();
        
    // render the rod shadow
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    glDisable(GL_DEPTH_TEST);
    glDepthMask(GL_FALSE);
    
    Vec3d light_direction = Vec3d(1, 2, 1).normalized();
    Vec3d ground_normal = Vec3d(0, 1, 0);
    Mat3d light_ground_projection = Mat3d::Identity() - light_direction * ground_normal.transpose() / light_direction.dot(ground_normal);
    
    Scalar light_source_area_factor = 2.0;
    Scalar alpha_scalar = 0.6;
    for (int j = 0; j < m_rod.ne(); j++) 
    {
        const Vec3d & v0 = m_rod.getVertex(j);
        const Vec3d & v1 = m_rod.getVertex((j + 1) % m_rod.nv());
        
        glBegin(GL_QUAD_STRIP);
        for (int k = 0; k <= m_tube.getSlices(); k++) 
        {
            Vec3d x0 = m_tube.getPointAtVert(j, k);
            Vec3d x0_protruded = v0 + (x0 - v0) * (x0.y() + light_source_area_factor) / light_source_area_factor;
            Vec3d x0_projected = light_ground_projection * x0_protruded;
            Vec3d n0 = (x0 - v0).normalized();
            double alpha0 = pow(fabs(n0.dot(light_direction)), 2.0) / ((x0.y() + light_source_area_factor) / light_source_area_factor);
            
            Vec3d x1 = m_tube.getPointAtVert(j + 1, k);
            Vec3d x1_protruded = v1 + (x1 - v1) * (x1.y() + light_source_area_factor) / light_source_area_factor;
            Vec3d x1_projected = light_ground_projection * x1_protruded;
            Vec3d n1 = (x1 - v1).normalized();
            double alpha1 = pow(fabs(n1.dot(light_direction)), 2.0) / ((x1.y() + light_source_area_factor) / light_source_area_factor);
            
            glColor4f(0, 0, 0, alpha0 * alpha_scalar);
            OpenGL::vertex(x0_projected);
            glColor4f(0, 0, 0, alpha1 * alpha_scalar);
            OpenGL::vertex(x1_projected);
            
        }
        glEnd();
    }
    
    glDepthMask(GL_TRUE);
    glEnable(GL_DEPTH_TEST);


    // render the scene
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glDisable(GL_BLEND);

    // render the belt (for depth only; color will be overwritten later)
    glColor4f(0, 0, 0, 1);
    glBegin(GL_QUADS);
    glVertex3f(-100, 0, r);
    glVertex3f(-100, 0, l);
    glVertex3f(18, 0, l);
    glVertex3f(18, 0, r);
    for (int i = 0; i < 12; i++)
    {
        Scalar t1 = (i + 0) * 2 * M_PI / 12;
        Scalar t2 = (i + 1) * 2 * M_PI / 12;
        glVertex3f(18 + 1 * sin(t1), -1 + 1 * cos(t1), r);
        glVertex3f(18 + 1 * sin(t1), -1 + 1 * cos(t1), l);
        glVertex3f(18 + 1 * sin(t2), -1 + 1 * cos(t2), l);
        glVertex3f(18 + 1 * sin(t2), -1 + 1 * cos(t2), r);
    }
    glEnd();    
    
    // render the reflection texture onto the belt
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, texReflection);

    glColor4f(0.2, 0.225, 0.25, 1.0);
//    glColor4f(1, 1, 1, 1);
    glDepthMask(GL_FALSE);
    double z = 0;
    glBegin(GL_QUADS);
    glTexCoord2f(0.0, 0.0); glVertex3f(-1.0, -1.0, z);
    glTexCoord2f(0.0, 1.0); glVertex3f(-1.0, 1.0, z);
    glTexCoord2f(1.0, 1.0); glVertex3f(1.0, 1.0, z);
    glTexCoord2f(1.0, 0.0); glVertex3f(1.0, -1.0, z);
    glEnd();
    glDepthMask(GL_TRUE);

//    // the red lines indicating the belt region
//    glDisable(GL_TEXTURE_2D);
//    glDisable(GL_DEPTH_TEST);
//    glColor4f(1, 0, 0, 1);
//    glBegin(GL_LINE_STRIP);
//    float w = 915;
//    float h = 481;
//    glVertex2f(-1, 1 - 2 * 427 / h);
//    glVertex2f(-1 + 2 * 416 / w, 1 - 2 * 157 / h);
//    glVertex2f(-1 + 2 * 881 / w, 1 - 2 * 254 / h);
//    glVertex2f(-1 + 2 * 781 / w, 1 - 2 * 479 / h);
//    glEnd();
//    glEnable(GL_DEPTH_TEST);

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    
    // the white bar besides the belt
    glDisable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, 0);
    
    glColor4f(1, 1, 1, 1);
    glBegin(GL_QUADS);
    glVertex3f(-100, 0, r);
    glVertex3f(16, 0, r);
    glVertex3f(16, 0, r - 2.0);
    glVertex3f(-100, 0, r - 2.0);
    glVertex3f(-100, 0, l);
    glVertex3f(16, 0, l);
    glVertex3f(16, 0, l + 2.0);
    glVertex3f(-100, 0, l + 2.0);
    glEnd();
    
    // render the real rod
    glEnable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);
    
    for (int j = 0; j < m_rod.ne(); j++) 
    {
        const Vec3d & v0 = m_rod.getVertex(j);
        const Vec3d & v1 = m_rod.getVertex((j + 1) % m_rod.nv());
        
        Color color1t, color2t;        
        color1t = Color(1.6, 2.0, 1.0, 1.0);
        color2t = Color(1.6, 2.0, 1.0, 1.0);
        
        glBegin(GL_QUAD_STRIP);
        for (int k = 0; k <= m_tube.getSlices(); k++) 
        {
            OpenGL::color(color1t);
            
            Vec3d x0 = m_tube.getPointAtVert(j, k);
            Vec3d n = (x0 - v0).normalized();
            OpenGL::normal(n);
            OpenGL::vertex(x0);
            
            Vec3d x1 = m_tube.getPointAtVert(j + 1, k);
            n = (x1 - v1).normalized();
            OpenGL::normal(n);
            OpenGL::vertex(x1);
        }
        glEnd();
    }
    
    // render the walls around the setup
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, texWall);
    
    glColor4f(1, 1, 1, 1);
    glBegin(GL_QUADS);
    glNormal3f(0, 0, 1);
    glTexCoord2f(0.01, 0.01); glVertex3f(-40, -5, -15);
    glTexCoord2f(0.99, 0.01); glVertex3f(19, -5, -15);
    glTexCoord2f(0.99, 0.50); glVertex3f(19, 30, -15);
    glTexCoord2f(0.01, 0.50); glVertex3f(-40, 30, -15);
    glNormal3f(-1, 0, 0);
    glTexCoord2f(0.01, 0.50); glVertex3f(19, -5, -15);
    glTexCoord2f(0.99, 0.50); glVertex3f(19, -5, 20);
    glTexCoord2f(0.99, 0.99); glVertex3f(19, 30, 20);
    glTexCoord2f(0.01, 0.99); glVertex3f(19, 30, -15);
    glNormal3f(0, 1, 0);
    glTexCoord2f(0.01, 0.50); glVertex3f(19, -5, -15);
    glTexCoord2f(0.99, 0.50); glVertex3f(19, -5, 20);
    glTexCoord2f(0.99, 0.99); glVertex3f(-40, -5, 20);
    glTexCoord2f(0.01, 0.99); glVertex3f(-40, -5, -15);
    glEnd();
    
    glDisable(GL_BLEND);
    glDisable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, 0);
    glShadeModel(GL_FLAT);
    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);
    glClearColor(161.0/255.0,161.0/255.0,161.0/255.0,0.0); */
}
  
} // namespace BASim
