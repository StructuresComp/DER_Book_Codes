/**
 * \file BASimulator.hh
 *
 * \author miklos@cs.columbia.edu (based on problem-execution.h from sar2120@columbia.edu)
 * \date 09/09/2009
 */

#ifndef BASIMULATOR_HH
#define BASIMULATOR_HH

#include <BASim/BASim>
#include "Problems/ProblemBase.hh"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

using namespace BASim;

void CreateProblemVector();
void PrintProblemTypes();

void RunProblem(int argc, char** argv);

void SetLighting();
void InitMenu();
void InitCamera();

void menu(int id);

void scaleMousePos(const int x, const int y, Scalar* xx, Scalar* yy);
bool saveScreen(const std::string& filename);

void LoadNextFrame();

void drawBox(Scalar lx, Scalar ly, Scalar lz);
//void drawMesh(SMV::Mesh* mesh);
void drawPlane(Scalar a, Scalar b, Scalar c, Scalar d);

void display(void);
void reshape(int w, int h);
void idle();
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void keyboard(unsigned char key, int, int);

#endif // BASIMULATOR_HH
