/**
 * \file BASimulator.cc
 *
 * \author miklos@cs.columbia.edu (based on problem-setup.cc from sar2120@columbia.edu)
 * \date 09/09/2009
 */

#include <BASim/BASim>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <sys/time.h>
#include <tclap/CmdLine.h>
#include <typeinfo>

#include "BASimulator.hh"

#include "BASim/src/Render/TriangleMeshRenderer.hh"

#include "BASim/src/IO/XMLSceneOutputter.hh"

// Might not work in windows? :'(
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>


std::vector<Problem*> problems;
#include "CreateProblemVector.inl"

Problem* current_problem = NULL;
ViewController controller;

// Renderable objects include RodRenderers and TriangleMeshRenderers.
std::vector<RenderBase*> renderable_objects;
std::vector<RodRenderer*> rod_renderers;
std::vector<TriangleMeshRenderer*> triangle_mesh_renderers;

bool render = true;
bool paused = true;
bool continuous = true;
int window_width = 512;
int window_height = 512;
int max_frames = -1;
Scalar max_time = -1;
bool progress_indicator = false;
bool generate_movie = false;

// Display the simulation time or not
bool g_dsp_sim_tm = true;

// For outputting time-stamped simulation captures,
// for matching current simulation run to output directory.
time_t g_rawtime;
struct tm* g_timeinfo;

// For dumping movies, etc
Scalar fps = 24;
double frame_period;
int current_frame;
std::string outputdirectory;

void cleanup()
{
  Timer::report();
  if (current_problem != NULL) current_problem->BaseFinalize();
  for( size_t i = 1; i < problems.size(); ++i )
  {
    assert( problems[i] != NULL );
    delete problems[i];
  }
  
  for( int i = 0; i < (int) renderable_objects.size(); ++i )
  {
    // RenderBases should only be deleted here
    assert( renderable_objects[i] != NULL );
    if( renderable_objects[i] != NULL ) delete renderable_objects[i];
    renderable_objects[i] = NULL;
  }
}

void SetLighting()
{
  // Create a directional white light with a small ambient component
  glEnable(GL_LIGHT0);
  GLfloat white_ambient[] = {0.1,0.1,0.1,1.0};
  glLightfv(GL_LIGHT0,GL_AMBIENT,white_ambient);
  GLfloat white_diffuse[] = {0.55,0.55,0.55,1.0};
  glLightfv(GL_LIGHT0,GL_DIFFUSE,white_diffuse);
  GLfloat upper_corner[] = {1.0,1.0,1.0,0.0};
  glLightfv(GL_LIGHT0,GL_POSITION,upper_corner);
  
  // Create a much weaker direction light
  glEnable(GL_LIGHT1);
  GLfloat weak_white_diffuse[] = {0.3,0.3,0.3,1.0};
  glLightfv(GL_LIGHT1,GL_DIFFUSE,weak_white_diffuse);
  GLfloat negative_z[] = {0.0,0.0,1.0,0.0};
  glLightfv(GL_LIGHT1,GL_POSITION,negative_z);
  
  glShadeModel(GL_FLAT);
}

//void SetMaterial()
//{
//  GLfloat mat_specular[] = {1.0, 1.0, 1.0, 1.0};
//  GLfloat mat_shininess[] = {50.0};
//  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
//  glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
//  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
//  glEnable(GL_COLOR_MATERIAL);
//}

void InitMenu()
{
  const int centerMenu = glutCreateMenu(menu);
  glutAddMenuEntry("Object (c)", 'c');
  glutAddMenuEntry("World origin (C)", 'C');

  const int showMenu = glutCreateMenu(menu);
  glutAddMenuEntry("Material (1)", '1');
  glutAddMenuEntry("Reference (2)", '2');
  //glutAddMenuEntry("Bishop (3)", '3');
  //glutAddMenuEntry("Show u axis (4)", '4');
  //glutAddMenuEntry("Show v axis (5)", '5');
  //glutAddMenuEntry("Curvature (6)", '6');
  glutAddMenuEntry("Mode (m)", 'm');
  glutAddMenuEntry("Draw arrows (a)", 'a');
  //glutAddMenuEntry("Force (f)", 'f');
  glutAddMenuEntry("Camera (p)", 'p');
  glutAddMenuEntry("Scale to radius (r)", 'r');

  glutCreateMenu(menu);
  glutAddMenuEntry("Quit (q)", 'q');
  glutAddMenuEntry("Save screenshot (s)", 's');

  glutAddSubMenu("View center", centerMenu);
  glutAddSubMenu("Display", showMenu);
  //glutAddMenuEntry("Write to file (w)", 'w');
  //glutAddMenuEntry("Print vertex velocities (v)", 'v');

  glutAttachMenu(GLUT_RIGHT_BUTTON);
}

Vec3d calcSimCenter()
{
  Vec3d center = Vec3d::Zero();
  int n = 0;

  for( int i = 0; i < (int) renderable_objects.size(); ++i ) 
  {
    center += renderable_objects[i]->calculateObjectCenter();
    ++n;
  }
  
  if( n != 0 ) center /= ((double)n);

  return center;
}

Scalar calcViewRadius(Vec3d& simCenter)
{
  Scalar radius = 0.0;

  for( int i = 0; i < (int) renderable_objects.size(); ++i ) 
  {
    Vec3d center = renderable_objects[i]->calculateObjectCenter();
    Scalar r = renderable_objects[i]->calculateObjectBoundingRadius(center);
    radius = std::max(radius, r + (center - simCenter).norm());
  }
  
  if( radius == 0.0 ) radius += 0.1;
  
  return radius;
}

void centerObject()
{
  Vec3d simCenter = calcSimCenter();
  controller.setCenterMode(ViewController::CENTER_OBJECT);
  controller.setViewCenter(simCenter);

  Scalar radius = calcViewRadius(simCenter);
  controller.setBoundingRadius(radius / 1.5);
}

void InitCamera()
{
  controller.setViewDirection(Vec3d(0, 0, -2));
  centerObject();
  /*
    if (current_problem->GetBoolOpt("read-camera")) {
    Vec3d& eye = current_problem->GetVecOpt("eye");
    Vec3d& up = current_problem->GetVecOpt("up");
    Vec3d& center = current_problem->GetVecOpt("center");

    SMV::Vec3 e(eye(0), eye(1), eye(2));
    SMV::Vec3 u(up(0), up(1), up(2));
    SMV::Vec3 c(center(0), center(1), center(2));
    SMV::Camera& cam = controller.getCamera();
    cam.setEye(e);
    cam.setUp(u);
    cam.setViewCenter(c);

    } else {
    SMV::Vec3 view_dir(0.0, 0.0, -2.0);
    controller.setViewDirection(view_dir);
    centerObject();
    }*/
	

}

void renderBitmapString( float x, float y, float z, void *font, std::string s ) 
{
  glDisable(GL_LIGHTING);
  glColor3f(0.0,0.0,0.0);
	glRasterPos3f(x, y, z);
	for (std::string::iterator i = s.begin(); i != s.end(); ++i)
	{
		char c = *i;
		glutBitmapCharacter(font, c);
	}
  glEnable(GL_LIGHTING);
}

void drawSimpleAxis()
{
  glDisable(GL_LIGHTING);
  
  // Label x axis
  GLvoid *font_style = GLUT_BITMAP_HELVETICA_18;
  glColor3f(1.0,0.0,0.0);
  glRasterPos3f(1.0,0.0,0.0);
  glutBitmapCharacter(font_style,'x');
  
  // Label y axis
  glColor3f(0.0,1.0,0.0);
  glRasterPos3f(0.0,1.0,0.0);
  glutBitmapCharacter(font_style,'y');
  
  // Label z axis
  glColor3f(0.0,0.0,1.0);
  glRasterPos3f(0.0,0.0,1.0);
  glutBitmapCharacter(font_style,'z');
  
  glLineWidth(1.0);
  
  // Red x axis
  glColor3f(1.0,0.0,0.0);
  glBegin(GL_LINES);
  glVertex3f(0.0, 0.0, 0.0);
  glVertex3f(1.0, 0.0, 0.0);
  glEnd();
  
  // Green y axis
  glColor3f(0.0,1.0,0.0);
  glBegin(GL_LINES);
  glVertex3f(0.0, 0.0, 0.0);
  glVertex3f(0.0, 1.0, 0.0);
  glEnd();
  
  // Blue z axis
  glColor3f(0.0,0.0,1.0);
  glBegin(GL_LINES);
  glVertex3f(0.0, 0.0, 0.0);
  glVertex3f(0.0, 0.0, 1.0);
  glEnd();
  
  glEnable(GL_LIGHTING);
}

void display()
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
  glPushMatrix();

  controller.ApplyCamera();
  
  //controller.setDefault2D();
  //drawSimpleAxis();
  
  if (current_problem->update_display && render) {
    current_problem->update_display = false;
    
    for( int i = 0; i < (int) renderable_objects.size(); ++i )
    {
      // RenderBases should only be deleted here
      assert( renderable_objects[i] != NULL );
      if( renderable_objects[i] != NULL ) delete renderable_objects[i];
      renderable_objects[i] = NULL;
    }  
    renderable_objects.clear();
  
    World& world = current_problem->getWorld();
    World::Objects& objects = world.getObjects();

    World::Objects::iterator it;
    for (it = objects.begin(); it != objects.end(); ++it) 
    {
      assert( (*it) != NULL );
      
      ObjectBase* objptr = *it;

      if( typeid(*objptr)==typeid(AnisotropicRod) )
      {
        ElasticRod* rod = dynamic_cast<ElasticRod*>(objptr);
        if (rod == NULL) continue;
        rod_renderers.push_back(new RodRenderer(*rod));
        renderable_objects.push_back(rod_renderers.back());
        
		    bool& b = rod_renderers.back()->drawMaterial();
		    b = true;
        
      }
      else if( typeid(*objptr)==typeid(TriangleMesh) )
      {
        TriangleMesh* tripobj = dynamic_cast<TriangleMesh*>(objptr);
        if( tripobj == NULL ) continue;
        triangle_mesh_renderers.push_back(new TriangleMeshRenderer(*tripobj));
        renderable_objects.push_back(triangle_mesh_renderers.back());
      }
      else
      {
        std::cout << "Unknown object encountered" << std::endl;
      }
    }  
  }
  
  for( int i = 0; i < (int) renderable_objects.size(); ++i ) renderable_objects[i]->render();
  current_problem->BaseRender();
  
  // set HUD transformation
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	
	gluOrtho2D(0, window_width, 0, window_height);
	
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

  glColor3f(0.0,0.0,0.0);
  //if( g_dsp_sim_tm ) renderBitmapString( 5, window_height-20, 0.0, GLUT_BITMAP_HELVETICA_18, current_problem->ProblemName() );
  if( g_dsp_sim_tm ) renderBitmapString( 5, window_height-40, 0.0, GLUT_BITMAP_HELVETICA_18, toString(current_problem->getTime()) );
  
  glDisable(GL_LIGHTING);
  current_problem->BaseHUDRender(window_width, window_height);
  glEnable(GL_LIGHTING);
  
  static double msg_max = 0;
  
  if (!RodRenderer::lastInstance() || RodRenderer::lastInstance()->getMode() != RodRenderer::PHOTOREALISTIC)
    renderBitmapString( 5, window_height-80, 0.0, GLUT_BITMAP_HELVETICA_18, toString(current_problem->msg_val) );
  
  if (msg_max < current_problem->msg_val) msg_max = current_problem->msg_val;
  
  if (!RodRenderer::lastInstance() || RodRenderer::lastInstance()->getMode() != RodRenderer::PHOTOREALISTIC)
    renderBitmapString( 5, window_height-120, 0.0, GLUT_BITMAP_HELVETICA_18, toString(msg_max) );
  
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();

  glutSwapBuffers();
}

void reshape(int w, int h)
{
  window_width = w;
  window_height = h;
  
  Camera& c = controller.getCamera();
  if (!RodRenderer::lastInstance() || RodRenderer::lastInstance()->getMode() != RodRenderer::PHOTOREALISTIC)
    c.setPerspective(60, 1);
  else
    c.setPerspective(20, 1);
	
//	c.setOrthographic(-50, 50, -50, 50);

  Scalar radius = controller.getBoundingRadius();
  c.setZClipping(0.1 * radius, 3 * radius);
  c.setViewport(w, h);

//		Camera& c = controller.getCamera();

	if (0) {
		c.setEye(Vec3d(62.4, 46.94, 26.45));
		c.setUp(Vec3d(-0.17, 0.964, -0.2));
		c.setViewCenter(Vec3d(32.873, 4.5159, -7.44));
	}

}

void idle()
{
  if (max_frames != -1) {
    if (current_problem->dynamicsPropsLoaded()) {
      if (current_problem->getTime() * fps >= max_frames) {
        std::cout << "Computed " << max_frames << " frames. Exiting"
                  << std::endl;
        exit(0);
      }
    }
  }

  if (max_time != -1) {
    if (current_problem->dynamicsPropsLoaded()) {
      if (current_problem->getTime() >= max_time) {
        std::cout << "Computed up to time " << max_time << ". Exiting"
                  << std::endl;
        exit(0);
      }
    }
  }

  if (!paused) {
    if (!continuous) paused = true;
    current_problem->BaseAtEachTimestep();
    if (render) glutPostRedisplay();
    //if(render&&generate_movie) 
    //{
    //  int file_width = 8;
    //  std::stringstream name;
    //  name << std::setfill('0');
    //  name << "screencap/frame_" << std::setw(file_width) << int(current_problem->getTime()/current_problem->getDt()) << ".png";
    //  saveScreen(name.str());
    //  std::cout << "Saving screencap as: " << name.str() << std::endl;
    //}    
  }
  
  if( render && generate_movie )
  {
    if( floor(current_problem->getTime()/frame_period) >= current_frame ) 
    {
      //std::cout << outputdirectory << std::endl;
      mkdir(outputdirectory.c_str(), 0755);
  
      int file_width = 20;

      //if( generate_movie == 2 || generate_movie == 3 )
      //{
        //std::stringstream framedirname;
        //framedirname << std::setfill('0');
        //framedirname << outputdirectory << "/frame" << std::setw(file_width) << current_frame;
        //
        ////std::cout << framedirname.str() << std::endl;
        //mkdir(framedirname.str().c_str(), 0755);
        //
        //XMLSceneOutputter scenesaver;
        //scenesaver.outputScene( framedirname.str(), "scene.xml", *current_problem );
        //
        //std::cout << "Frame: " << current_frame << "   Time: " << current_problem->getTime() << "   Scenedir: " << framedirname.str() << std::endl;
      //}
      
      //if( generate_movie == 1 || generate_movie == 3 )
      //{
        glutPostRedisplay();
        
        std::stringstream name;
        name << std::setfill('0');
        name << outputdirectory << "/frame" << std::setw(file_width) << current_frame << ".png";

        saveScreen(name.str());
        std::cout << "Frame: " << current_frame << "   Time: " << current_problem->getTime() << "   Screencapture: " << name.str() << std::endl;
      //}
      
      ++current_frame;
    }
  }

  if (progress_indicator && current_problem->dynamicsPropsLoaded()) {
    std::cout << "\rtime = " << current_problem->getTime() << std::flush;
  }
}

void scaleMousePos(const int x, const int y, Scalar& xx, Scalar& yy)
{
  int w, h;
  controller.getCamera().getViewport(&w, &h);

  xx = 2 *           x / (Scalar) (w - 1) - 1.0;
  yy = 2 * (h - y - 1) / (Scalar) (h - 1) - 1.0;
}

void mouse(int button, int state, int x, int y)
{
  const bool zooming     = (button == GLUT_MIDDLE_BUTTON) ||
    ((button == GLUT_LEFT_BUTTON) && (glutGetModifiers() & GLUT_ACTIVE_CTRL));
  const bool translating = (button == GLUT_LEFT_BUTTON)
    && (glutGetModifiers() & GLUT_ACTIVE_SHIFT);
  const bool rotating    = (button == GLUT_LEFT_BUTTON)
    && (glutGetModifiers() == 0);

  Scalar xx, yy;
  scaleMousePos(x, y, xx, yy);
  if (state == GLUT_DOWN) {
    if (translating)
      controller.beginTranslationDrag(xx, yy);

    if (zooming)
      controller.beginZoomDrag(xx, yy);

    if (rotating)
      controller.beginRotationDrag(xx, yy);

  } else {
    controller.endTranslationDrag(xx, yy);
    controller.endRotationDrag(xx, yy);
    controller.endZoomDrag(xx, yy);
  }

  glutPostRedisplay();
}

void motion(int x, int y)
{
  Scalar xx, yy;
  scaleMousePos(x, y, xx, yy);
  controller.updateDrag(xx, yy);
  glutPostRedisplay();
}

void keyboard(unsigned char key, int, int)
{
  menu(key);
}

void menu(int id)
{
  switch(id) {
  case 'q':
    exit(0);
    break;

  case ' ':
    paused = !paused;
    break;

  case 'i':
    continuous = !continuous;
      break;
      
  case '1': {
    for( int i = 0; i < (int) rod_renderers.size(); ++i )
    {
      bool& b = rod_renderers[i]->drawMaterial();
      b = !b;
    }
    glutPostRedisplay();
    break;
  }

  case '2': {
    for( int i = 0; i < (int) rod_renderers.size(); ++i )
    {
      bool& b = rod_renderers[i]->drawReference();
      b = !b;
    }
    glutPostRedisplay();
    break;
  }

  case 'm': {
    for( int i = 0; i < (int) rod_renderers.size(); ++i )
    {
      RodRenderer::DrawMode mode = rod_renderers[i]->getMode();
      mode = (RodRenderer::DrawMode) ((mode + 1) % RodRenderer::NONE);
      rod_renderers[i]->setMode(mode);
    }
    for( int i = 0; i < (int) triangle_mesh_renderers.size(); ++i )
    {
      TriangleMeshRenderer::DrawMode mode = triangle_mesh_renderers[i]->getMode();
      mode = (TriangleMeshRenderer::DrawMode) ((mode + 1) % TriangleMeshRenderer::NONE);
      triangle_mesh_renderers[i]->setMode(mode);
    }
    if (RodRenderer::lastInstance() && RodRenderer::lastInstance()->getMode() == RodRenderer::PHOTOREALISTIC)
    {
      int w = glutGet(GLUT_WINDOW_WIDTH);
      int h = glutGet(GLUT_WINDOW_HEIGHT);
      if (w > h * 1.909)
        glutReshapeWindow(w, w / 1.909);
      else
        glutReshapeWindow(h * 1.909, h);
      controller.getCamera().setExperimentCamera();
      controller.getCamera().setPerspective(20, 1);
    } else
    {
      controller.getCamera().setPerspective(60, 1);
    }
    glutPostRedisplay();
    break;
  }

  case 'r': {
    for( int i = 0; i < (int) rod_renderers.size(); ++i )
    {
      bool& b = rod_renderers[i]->scaleToRadius();
      b = !b;
    }
    glutPostRedisplay();
    break;
  }

  case 'a': {
    for( int i = 0; i < (int) rod_renderers.size(); ++i )
    {
      bool& b = rod_renderers[i]->drawArrows();
      b = !b;
    }
    glutPostRedisplay();
    break;
  }

  case 's':
//    saveScreen("screen.png");

    saveScreen(current_problem->screenshot_filename);
    break;

  case 'd':
      generate_movie = !generate_movie;
      std::cout << "Dump movie: " << generate_movie << std::endl;
      
      if (generate_movie) {
        current_frame = floor(current_problem->getTime()/frame_period);
      }
    break;
      
  case 'c':
    centerObject();
    controller.setCenterMode(ViewController::CENTER_OBJECT);
    glutPostRedisplay();
    break;

  case 'C':
    controller.setCenterMode(ViewController::CENTER_WORLD_ORIGIN);
    glutPostRedisplay();
    break;

  case 'f':
    current_problem->BaseAtEachTimestep();
    glutPostRedisplay();
    break;

  case 'g': {
    Camera& c = controller.getCamera();
    std::cout << "    c.setEye(Vec3d( " << c.getEye()(0) << ", " << c.getEye()(1) << ", " << c.getEye()(2) << " ));\n";
    std::cout << "    c.setUp(Vec3d( " << c.getUp()(0) << ", " << c.getUp()(1) << ", " << c.getUp()(2) << " ));\n";
    std::cout << "    c.setViewCenter(Vec3d( " << c.getViewCenter()(0) << ", " << c.getViewCenter()(1) << ", " << c.getViewCenter()(2) << " ));\n";

    c.setEye(Vec3d( 389.406, 87.4058, 54.0277 ));
    c.setUp(Vec3d( -0.227028, 0.942404, -0.245628 ));
    c.setViewCenter(Vec3d( 328.1, -11.0899, -12.3007 ));


    glutPostRedisplay();
    
    break;
  }              
  case 'p':
  Camera& c = controller.getCamera();
  std::cout << "eye: " << c.getEye() << std::endl;
  std::cout << "up: " << c.getUp() << std::endl;
  std::cout << "center: " << c.getViewCenter() << std::endl;

//  c.setEye(Vec3d(-17.82, 79.06, 57.16));
//  c.setUp(Vec3d(0.121, 0.9819, -0.1449));
//  c.setViewCenter(Vec3d(25.936,11.850,6.12));

//  c.setEye(Vec3d(69.5, 46.887, -25.85));
//  c.setUp(Vec3d(-0.156, 0.9766, 0.147));
//  c.setViewCenter(Vec3d(38.42, 2.382, 2.754));

//  c.setEye(Vec3d(269, 109, 36));
//  c.setUp(Vec3d(-0.35, 0.86, -0.36));
//  c.setViewCenter(Vec3d(232, 23, -1.75));

//  c.setEye(Vec3d(62.4, 46.94, 26.45));
//  c.setUp(Vec3d(-0.17, 0.964, -0.2));
//  c.setViewCenter(Vec3d(32.873, 4.5159, -7.44));

  glutPostRedisplay();
  
//eye: {{-16.8287}, {79.062983}, {57.164423}}
//up: {{0.12142154}, {0.98196629}, {-0.14491036}}
//center: {{25.936425}, {11.850348}, {6.126444}}

  break;
            
  }
  /*
  switch(id) {
      case '2':
      {
      bool& showBishop = rod_renderer->GetOptions().show_bishop;
      showBishop = !showBishop;
      glutPostRedisplay();
      break;
      }
      case '4':
      {
      bool& showU = rod_renderer->GetOptions().show_u;
      showU = !showU;
      glutPostRedisplay();
      break;
      }
      case '5':
      {
      bool& showV = rod_renderer->GetOptions().show_v;
      showV = !showV;
      glutPostRedisplay();
      break;
      }
      case '6':
      {
      bool& showKappa = rod_renderer->GetOptions().show_kappa;
      showKappa = !showKappa;
      glutPostRedisplay();
      break;
      }
      case 'f':
      {
      bool& showForce = rod_renderer->GetOptions().show_force;
      showForce = !showForce;
      glutPostRedisplay();
      break;
      }
      case 'p':
      {
      SMV::Camera& c = controller.getCamera();
      std::cout << "eye: " << c.getEye() << std::endl;
      std::cout << "up: " << c.getUp() << std::endl;
      std::cout << "center: " << c.getViewCenter() << std::endl;
      break;
      }
      }
      case 'v':
      {
      std::vector<Rod*>& rods = current_problem->GetRods();
      for (uint i = 0; i < rods.size(); ++i) {
      Rod& rod = *rods[i];
      for (uint j = 0; j < rod.nv(); ++j) {
      std::cout << rod.velocity(j).norm() << std::endl;
      }
      }
      break;
      }

      case 'w':
      std::vector<Rod*>& rods = current_problem->GetRods();
      for (uint i = 0; i < rods.size(); ++i) {
      std::ostringstream filename;
      filename << "rod_" << i << ".txt";

      std::ofstream out_file;
      out_file.open(filename.str().c_str());
      if (!out_file.is_open()) {
      std::cout << "Failed to write rod to file " << filename << std::endl;

      } else {
      out_file << *rods[i];
      out_file.close();
      }
      }
      break;
      }
  */
}

/// Copy the current OpenGL render buffer and save it to a PNG-format file.
bool saveScreen(const std::string& filename)
{
#ifdef HAVE_PNG
  int w, h;
  controller.getCamera().getViewport(&w, &h);

  YImage image;
  image.resize(w, h);

  glReadPixels(0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, image.data());
  image.flip();
  image.setAllAlpha(255);

  std::cerr << "Screenshot was stored to " << filename 
            << std::endl;

  return image.save(filename.c_str());
#else
  std::cerr << "Not compiled with PNG support, can't save images."
            << std::endl;
  return false;
#endif
}

void initializeOpenGL(int argc, char** argv)
{
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_ALPHA);
  glutInitWindowPosition(100, 100);
  glutInitWindowSize(window_width, window_height);
  glutCreateWindow(argv[0]);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);
  glClearColor(161.0/255.0,161.0/255.0,161.0/255.0,0.0);
  //glClearColor(0.0/255.0,0.0/255.0,0.0/255.0,0.0);

  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);
  glutKeyboardFunc(keyboard);
  glutIdleFunc(idle);

  InitMenu();
  InitCamera();

  SetLighting();
  //SetMaterial();
}

void setOptions()
{
  current_problem->AddOption("render", "display output in OpenGL", render);
  current_problem->AddOption("paused", "start simulation paused", paused);
  current_problem->AddOption("continuous", "pause after every step", continuous);
  current_problem->AddOption("window-width", "width of the window", 512);
  current_problem->AddOption("window-height", "height of the window", 512);
  current_problem->AddOption("fps", "frames per second for output movie", fps);
  current_problem->AddOption("max-frames", "number of frames to compute", max_frames);
  current_problem->AddOption("max-time", "maximum (simulation) time", max_time);
  current_problem->AddOption("progress-indicator", "prints out time", progress_indicator);
  current_problem->AddOption("generate-movie", "generate_movie", generate_movie);
  
}

void getOptions()
{
  render = current_problem->GetBoolOpt("render");
  paused = render && current_problem->GetBoolOpt("paused");
  continuous = render && current_problem->GetBoolOpt("continuous");

  window_width = current_problem->GetIntOpt("window-width");
  window_height = current_problem->GetIntOpt("window-height");

  fps = current_problem->GetScalarOpt("fps");
  max_frames = current_problem->GetIntOpt("max-frames");
  max_time = current_problem->GetScalarOpt("max-time");
  progress_indicator = current_problem->GetBoolOpt("progress-indicator");
  
  if( fps > 0 ) 
  {
    frame_period = 1.0/((double)fps);
    current_frame = 0;
  }
    
  generate_movie = current_problem->GetBoolOpt("generate-movie");
}

void RunProblem(int argc, char** argv)
{
  current_problem->BaseSetup(argc, argv);
  
  // TODO: This rendering is a mess! dynamic_cast of pointers, tyepid, the horror!
  if (render) {
    initializeOpenGL(argc, argv);
    World& world = current_problem->getWorld();
    World::Objects& objects = world.getObjects();

    World::Objects::iterator it;
    for (it = objects.begin(); it != objects.end(); ++it) 
    {
      assert( (*it) != NULL );
      
      ObjectBase* objptr = *it;

      if( typeid(*objptr)==typeid(AnisotropicRod) )
      {
        ElasticRod* rod = dynamic_cast<ElasticRod*>(objptr);
        if (rod == NULL) continue;
        rod_renderers.push_back(new RodRenderer(*rod));
        renderable_objects.push_back(rod_renderers.back());
        
		    bool& b = rod_renderers.back()->drawMaterial();
		    b = true;
        
      }
      else if( typeid(*objptr)==typeid(TriangleMesh) )
      {
        TriangleMesh* tripobj = dynamic_cast<TriangleMesh*>(objptr);
        if( tripobj == NULL ) continue;
        triangle_mesh_renderers.push_back(new TriangleMeshRenderer(*tripobj));
        renderable_objects.push_back(triangle_mesh_renderers.back());
      }
      else
      {
        std::cout << "Unknown object encountered" << std::endl;
      }
    }
    
    centerObject();

    glutMainLoop();

  } else {
    paused = false;
    continuous = true;
    while(1) idle();
  }
}

void PrintProblemTypes()
{
  for (size_t i = 1; i < problems.size(); ++i) {
    Problem& prob = *problems[i];
    std::cout << i << ": " << prob.ProblemName() << std::endl
              << "Description: " << prob.ProblemDescription() << std::endl
              << std::endl;
  }
}

void CreateOptionsFile(int idx, const std::string& filename)
{
  std::ofstream file;

  file.open(filename.c_str());
  if (!file.is_open()) {
    std::cerr << "Failed to open file " << filename << std::endl;
    return;
  }

  problems[idx]->PrintOptions(file);

  file.close();

  std::cout << "Generated options file for "
            << problems[idx]->ProblemName() << ": "
            << filename << std::endl;
}

template <typename T>
class ProblemConstraint : public TCLAP::Constraint<T>
{
public:

  ProblemConstraint(const T& lower, const T& upper)
    : m_lower(lower)
    , m_upper(upper)
  {
    m_typeDesc = toString(lower) + "-" + toString(upper);
  }

  virtual std::string description() const { return m_typeDesc; }
  virtual std::string shortID() const { return m_typeDesc; }

  virtual bool check(const T& value) const
  {
    return ((value >= m_lower) && (value <= m_upper));
  }

protected:

  T m_lower;
  T m_upper;
  std::string m_typeDesc;
};

int parseCommandLine(int argc, char** argv)
{
  ProblemConstraint<int> allowedProblems(1, problems.size() - 1);

  try {
    std::string version(BASIM_VERSION);
    TCLAP::CmdLine cmd("Send inquiries to basim.support@gmail.com", ' ',
                       version);

    TCLAP::SwitchArg
      print("p", "print", "Print available problems", cmd);

    TCLAP::ValueArg<int>
      run("r", "run", "Run a problem", false, 1, &allowedProblems, cmd);

    TCLAP::ValueArg<int>
      opts("o", "options", "Print a problem's options", false, 1,
           &allowedProblems, cmd);

    TCLAP::ValueArg<std::string>
      file("f", "file", "Options file for a problem", false, "",
           "string", cmd);

    TCLAP::ValueArg<int>
      gen("g", "generate","Generate a default options file for a problem",
          false, 1, &allowedProblems, cmd);

    TCLAP::ValueArg<std::string>
      solver("s", "solver", "File describing options for solver", false, "",
             "string", cmd);

    cmd.parse(argc, argv);

    if (solver.isSet()) {
      if (readSolverFile(solver.getValue()) != 0) return -1;
    }

    if (print.getValue()) {
      PrintProblemTypes();
      return -1;
    }

    if (opts.isSet()) {
      int idx = opts.getValue();
      std::cout << "# Options for " << problems[idx]->ProblemName() << std::endl;
      problems[idx]->PrintOptions(std::cout);
      return -1;
    }

    if (gen.isSet()) {
      if (!file.isSet()) {
        TCLAP::ArgException e("Must also specify -" + file.getFlag(),
                              "-" + gen.getFlag());
        throw(e);
      }

      CreateOptionsFile(gen.getValue(), file.getValue());
      return -1;
    }

    if (run.isSet()) {
      int idx = run.getValue();
      current_problem = problems[idx];
      setOptions();

      if (file.isSet()) {
        if (current_problem->LoadOptions(file.getValue()) == -1) return -1;
      } else {
        std::cout << "No options file specified. Using default options."
                  << std::endl;
      }
      if (current_problem->LoadOptions(argc, argv) == -1) return -1;
      getOptions();
      return idx;
    }

    std::cerr << cmd.getProgramName() << ": missing operand" << std::endl
              << "Try `" << cmd.getProgramName() << " --help'"
              << " for more information" << std::endl;

  } catch (TCLAP::ArgException& e) {
    std::cerr << "ERROR: " << e.argId() << std::endl
              << "       " << e.error() << std::endl;
  }

  return -1;
}


void printCommandLineSplashScreen()
{
  std::cout << "\033[32;1m";
  std::cout << "----------------------------------" << std::endl;
  std::cout << "  ____    _    ____  _" << std::endl;
  std::cout << " | __ )  / \\  / ___|(_)_ __ ___" << std::endl;
  std::cout << " |  _ \\ / _ \\ \\___ \\| | '_ ` _ \\" << std::endl;
  std::cout << " | |_) / ___ \\ ___) | | | | | | |" << std::endl;
  std::cout << " |____/_/   \\_\\____/|_|_| |_| |_|" << std::endl;
  std::cout << "----------------------------------" << std::endl;
  std::cout << "\033[m";
  std::cout << std::endl;
  
  #ifdef DEBUG
  std::cout << " Build mode: DEBUG" << std::endl;
  #else
  std::cout << " Build mode: RELEASE" << std::endl;
  #endif
  std::cout << std::setfill('0');
  std::cout << " Timestamp:  " << (1900+g_timeinfo->tm_year) << "/" << std::setw(2) << (1+g_timeinfo->tm_mon) << "/";
  std::cout << (g_timeinfo->tm_mday) << "  " << std::setw(2) << (g_timeinfo->tm_hour) << ":" << std::setw(2) << (g_timeinfo->tm_min);
  std::cout << ":" << std::setw(2) << (g_timeinfo->tm_sec) << std::endl;
  
  std::cout << " Simulation: " << current_problem->ProblemName() << std::endl;
  
  std::cout << " Linear Solver:   " << SolverUtils::instance()->getSolverName() << std::endl;

  #ifdef HAVE_OPENMP
  std::cout << " OpenMP: Enabled" << std::endl;
  std::cout << " OpenMP Max Threads: " << omp_get_max_threads() << std::endl;
  #else
  std::cout << " OpenMP: Disabled" << std::endl;
  #endif
  
  std::cout << std::endl;
}

std::string generateOutputDirName()
{
  assert( g_timeinfo != NULL );

  std::stringstream datestream;
  datestream.fill('0');
  datestream << std::setw(4) << (1900+g_timeinfo->tm_year) << "_" << std::setw(2) << (1+g_timeinfo->tm_mon) << "_";
  datestream << std::setw(2) << (g_timeinfo->tm_mday) << "_" << std::setw(2) << (g_timeinfo->tm_hour) << "_" << std::setw(2) << (g_timeinfo->tm_min) << "_";
  datestream << std::setw(2) << (g_timeinfo->tm_sec) << "_" << "simulation_capture";

  return datestream.str();
}

int main(int argc, char** argv)
{
  // Generate a directory name with the date and time
  time(&g_rawtime);
  g_timeinfo = localtime(&g_rawtime);
  outputdirectory = generateOutputDirName();
  
  // Some stuff for dumping movies
  if( fps > 0 ) 
  {
    frame_period = 1.0/((double)fps);
    current_frame = 0;
  }

  CreateProblemVector();
  atexit(cleanup);
  if (parseCommandLine(argc, argv) < 0) return -1;
  
  printCommandLineSplashScreen();
  
  RunProblem(argc, argv);
  return 0;
}
