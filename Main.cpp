#include <string>
#include <iostream>
#include <stdlib.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "MyWorld.h"
#include "Timer.h"
#include <Math.h>

using namespace std;

// opengl setup related variables
unsigned int window_width = 640, window_height = 640;

// ui related variables
bool mouse_down = false;
bool leftClick = false;
bool rightClick = false;
int mouseX;
int mouseY;
bool viewVelocity = false;

// simulation related variables
MyWorld mySimulator;
bool simulating = false;
int frame_number = 0;
Timer timer;
int numCells = 64; // Number os cells in a row/column

bool screenSaverOn = 0;

// simulation functions
int getNumCells() {
  return numCells;
}

// opengl functions
void myGlutResize(int w, int h);

void myGlutIdle(void);

void myGlutDisplay(void);

void myGlutKeyboard(unsigned char key, int x, int y);

void myGlutMouse(int button, int state, int x, int y);

void myGlutMotion(int x, int y);

void drawVelocity();

void drawDensity();

void initializeFields();

// main function
int main(int argc, char *argv[])
{
    mySimulator.initialize(numCells, 0.1, 0.0, 0.0);

    if (argc == 2 && strcmp(argv[1], "-p") == 0)
      screenSaverOn = 1;
    
    if (screenSaverOn)    
      initializeFields();
    
    
    glutInit(&argc, argv);
    glutInitWindowSize(window_width, window_height);
    glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
    glutCreateWindow("Fluid Sim");
    glutIdleFunc(myGlutIdle);
    glutDisplayFunc(myGlutDisplay);
    glutReshapeFunc(myGlutResize);
    glutKeyboardFunc(myGlutKeyboard);
    glutMouseFunc(myGlutMouse);
    glutMotionFunc(myGlutMotion);
    
    // anti aliasing
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    glutMainLoop();
    return 0;
}

void myGlutResize(int w, int h)
{
    window_width = w;
    window_height = h;
    glViewport(0,0,w,h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glutPostRedisplay();
}

void myGlutIdle(void) {
    if (simulating) {
        timer.stop();
        double time_diff_in_sec = timer.getLastElapsedTime();
        if (time_diff_in_sec > 0.01) {
            while (time_diff_in_sec > 0.01) {
                mySimulator.simulate();
                frame_number++;
                time_diff_in_sec -= 0.01;
            }
            timer.start();
        }
    }
  
    glutPostRedisplay();
}

void myGlutDisplay(void) {
    glClearColor(1.f , 1.f, 1.f ,1.0f);
    ::glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_COLOR_MATERIAL);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
    
    glViewport ( 0, 0, window_width, window_height);
    glMatrixMode(GL_PROJECTION);    // opengl matrix for camera
    glLoadIdentity();
    gluOrtho2D ( 0.0, 1.0, 0.0, 1.0 );

    // lighting
    glEnable(GL_LIGHTING);
    float ambient[4] = {0.5, 0.5, 0.5, 1};
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);
    float diffuse[4] = {0.5, 0.5, 0.5, 1};
    float position[4] = {10, 10, 10, 0};
    glEnable(GL_LIGHT0);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
    glLightfv(GL_LIGHT0, GL_POSITION, position);
    

    if (viewVelocity)
      drawVelocity();
    else
      drawDensity();

    glutSwapBuffers();
}

void myGlutKeyboard(unsigned char key, int x, int y) {
    switch (key) {
        case 27:    // esc
            exit(0);
            break;
        case ' ':   // toggle simulation
            simulating =  !simulating;
            if (simulating) timer.start();
            break;
        case 'v': // toggle between velocity view and density view
            viewVelocity = !viewVelocity;
            break;
                case 'c': // toggle the color being drawn
                    mySimulator.toggleColor();
        default:
            break;
    }
    
    glutPostRedisplay();
}

void myGlutMouse(int button, int state, int x, int y) {

  if(screenSaverOn)
    return;

  mouse_down = (state == GLUT_DOWN);
        
    if(mouse_down){
        mouseX = x;
        mouseY = y;
        int i = (int)((mouseX / (double)window_width) * mySimulator.getNumCells() + 1);
        int j = (int)(((window_height - mouseY) / (double)window_height) * mySimulator.getNumCells() + 1);
	
        if (i < 1 || i > mySimulator.getNumCells() || j < 1 || j > mySimulator.getNumCells())
            return;
        
        if (button == GLUT_LEFT_BUTTON) {
            leftClick = true;
            mySimulator.setDensity(i, j, 100.0);
     
        } else if (button == GLUT_RIGHT_BUTTON || button == GLUT_MIDDLE_BUTTON) {
            rightClick = true;
            mySimulator.setU(i, j, 5.0);
            mySimulator.setV(i, j, 5.0);
        }
    } else {
        leftClick = false;
        rightClick = false;
    }
    glutPostRedisplay();
}

void myGlutMotion(int x, int y) {    

  if(screenSaverOn)
      return;

    int i = (int)((x / (double)window_width) * mySimulator.getNumCells() + 1);
    int j = (int)(((window_height - y) / (double)window_height) * mySimulator.getNumCells() + 1);
	
    if (i < 1 || i > mySimulator.getNumCells() || j < 1 || j > mySimulator.getNumCells())
        return;
        
    if (leftClick) {
        mySimulator.setDensity(i, j, 100.0);

    } else if (rightClick) {        
        mySimulator.setU(i, j, x - mouseX);
        mySimulator.setV(i, j, mouseY - y);
    }

    mouseX = x;
    mouseY = y;
    glutPostRedisplay();
}

void RenderBitmapString(float x, float y, void *font,char *string)
{
    char *c;
    ::glRasterPos2f(x, y);
    for (c=string; *c != '\0'; c++) {
        ::glutBitmapCharacter(font, *c);
    }
    ::glRasterPos2f(x+1, y);
    for (c=string; *c != '\0'; c++) {
        ::glutBitmapCharacter(font, *c);
    }
}


void drawVelocity() {
    double h = 1.0 / mySimulator.getNumCells();

    glColor3f ( 1.0f, 1.0f, 1.0f );
    glLineWidth ( 1.0f );

    glBegin ( GL_LINES );
    for (int i=1 ; i <= mySimulator.getNumCells(); i++) {
        double x = (i - 0.5) * h;
        for (int j = 1; j <= mySimulator.getNumCells(); j++) {
            double y = (j - 0.5) * h;

            glVertex2f(x, y );
            glVertex2f (x + mySimulator.getVelocityU(IX(i,j)), y + mySimulator.getVelocityV(IX(i,j)));
        }
    }
    glEnd ();
}

void drawDensity() {
      double h = 1.0 / mySimulator.getNumCells();
    glBegin(GL_QUADS);
    for (int i = 0; i <= mySimulator.getNumCells(); i++) {
        double x = (i - 0.5) * h;
        for (int j = 0; j <= mySimulator.getNumCells(); j++) {
            double y = (j - 0.5) * h;

            double d00r = mySimulator.getDensity(0, IX(i, j));
            double d01r = mySimulator.getDensity(0, IX(i, j+1));
            double d10r = mySimulator.getDensity(0, IX(i+1, j));
            double d11r = mySimulator.getDensity(0, IX(i+1, j+1));   
            double d00 = mySimulator.getDensity(IX(i, j));
            double d01 = mySimulator.getDensity(IX(i, j+1));
            double d10 = mySimulator.getDensity(IX(i+1, j));
            double d11 = mySimulator.getDensity(IX(i+1, j+1));
            double d00g = mySimulator.getDensity(1, IX(i, j));
			double d01g = mySimulator.getDensity(1, IX(i, j + 1));
			double d10g = mySimulator.getDensity(1, IX(i + 1, j));
			double d11g = mySimulator.getDensity(1, IX(i + 1, j + 1));
			double d00b = mySimulator.getDensity(2, IX(i, j));
			double d01b = mySimulator.getDensity(2, IX(i, j + 1));
			double d10b = mySimulator.getDensity(2, IX(i + 1, j));
			double d11b = mySimulator.getDensity(2, IX(i + 1, j + 1));

            glColor3d(d00r, d00g, d00b); 
            glVertex3f(x, y, 0);
            glColor3d(d10r, d10g, d10b); 
            glVertex3f(x + h, y, 0);            
            glColor3d(d11r, d11g, d11b); 
            glVertex3f(x + h, y + h, 0);
            glColor3d(d01r, d01g, d01b); 
            glVertex3f(x, y + h, 0);
        }
    }
    glEnd();

}

void initializeFields() {
  //timer.start();
  mySimulator.setDensity(32, 32, 100.0);
  mySimulator.setU(32, 32, 10);
  mySimulator.setV(32, 32, 200);
  for (int i = 0; i < 64; i++){
      for(int j = 0;j< 64;j++) {
          if( (i*j) % 12 == 0) {
               mySimulator.setColor(j);
               mySimulator.setDensity(i, j, i*j %100);
               mySimulator.setU(i, j, 3*i);
               mySimulator.setV(i, j, 2*j);
          }
      }
  }
}
