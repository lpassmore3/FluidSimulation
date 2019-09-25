#include <string>
#include <iostream>
#include <stdlib.h>
#include <atlimage.h>
#include <atlstr.h>

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
unsigned int window_width = 320, window_height = 320;

// ui related variables
bool mouse_down = false;
bool leftClick = false;
bool rightClick = false;
int mouseX;
int mouseY;
bool viewVelocity = false;
bool viewDivergence = false;

char color = 'r'; // Sets the color being painted, default is red

CImage myImage; // Stores the image to be loaded
double R = 0, G = 0, B = 0; //These are the channel values for loading the image
COLORREF PixColor = 0; //This is a color data of the image

// simulation related variables
MyWorld mySimulator;
bool simulating = false;
int frame_number = 0;
Timer timer;
int numCells = 4; // Number os cells in a row/column
int numCells_I = 4;
int numCells_J = 4;

bool screenSaverOn = 0;


// simulation functions
int getNumCells() {
  return numCells;
}

int getNumCells_I() {
	return numCells_I;
}

int getNumCells_J() {
	return numCells_J;
}

// opengl functions
void myGlutResize(int w, int h);

void myGlutIdle(void);

void myGlutDisplay(void);

void myGlutKeyboard(unsigned char key, int x, int y);

void myGlutMouse(int button, int state, int x, int y);

void myGlutMotion(int x, int y);

void drawVelocity();

void drawDivergence();

void drawDensity();

void initializeFields();

// main function
int main(int argc, char *argv[])
{

	// Edit window size
	float cellRatio = numCells_I / numCells_J;
	/*if (numCells_I > numCells_J) {
		window_width = window_width * cellRatio;
	}
	else if (numCells_I < numCells_J) {
		window_height = window_height * cellRatio;
	}*/
	window_width = numCells_I * 10;
	window_height = numCells_J * 10;

	//std:cout << color << endl;

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
	else if (viewDivergence)
		drawDivergence();
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
            simulating = !simulating;
            if (simulating) timer.start();
            break;
        case 'v': // toggle between velocity view and density view
			viewDivergence = false;
			viewVelocity = !viewVelocity;
            break;
		case 'd': // toggle between divergence and density view
			viewVelocity = false;
			viewDivergence = !viewDivergence;
			break;
		case 'p': // print the average divergence
			cout << "Average Divergence: " << mySimulator.getAverageDivergence() << "\n";
			cout << "Average Absolute Divergence: " << mySimulator.getAverageAbsDivergence() << "\n\n";
		case 'r': // toggle to editing Red velocity field
			//std::cout << "r was pressed" << endl;
			color = 'r';
			break;
		case 'g': // toggle to editing Green velocity field
			//std::cout << "g was pressed" << endl;
			color = 'g';
			break;
		case 'b': // toggle to editing Blue velocity field
			//std::cout << "b was pressed" << endl;
			color = 'b';
			break;
		//case 'e': // toggle erasing (filling with Black)
		//	color = 'e';
		//	break;
		case 'i': // loads the saved image
			if (!screenSaverOn) {
				// draws the image file on screen
				if (numCells_I <= 64 && numCells_J <= 64) {
					myImage.Load(_T("moana_64x64.bmp"));
					for (int i = 0; i < numCells_I; i++) {
						for (int j = 0; j < numCells_J; j++) {
							PixColor = myImage.GetPixel(i, numCells - 1 - j);
							R = GetRValue(PixColor);//This macro extracts red channel value
							G = GetGValue(PixColor);//This macro extracts green channel value
							B = GetBValue(PixColor);//This extracts blue pixel color

							mySimulator.setDensity_R(i + 1, j + 1, (R / 255.0) * 10.0);
							mySimulator.setDensity_G(i + 1, j + 1, (G / 255.0) * 10.0);
							mySimulator.setDensity_B(i + 1, j + 1, (B / 255.0) * 10.0);
						}
					}
				}
				else {
					cout << "Number of cells is too large to load image.";
				}
			}
			break;
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
        int i = (int)((mouseX / (double)window_width) * mySimulator.getNumCells_I() + 1);
        int j = (int)(((window_height - mouseY) / (double)window_height) * mySimulator.getNumCells_J() + 1);
	
        if (i < 1 || i > mySimulator.getNumCells_I() || j < 1 || j > mySimulator.getNumCells_J())
            return;
        
        if (button == GLUT_LEFT_BUTTON) {
            leftClick = true;
			//mySimulator.setDensity(i, j, 100.0);
			if (color == 'r') {
				mySimulator.setDensity_R(i, j, 100.0);
			}
			else if (color == 'g') {
				mySimulator.setDensity_G(i, j, 100.0);
			}
			else if (color == 'b') {
				mySimulator.setDensity_B(i, j, 100.0);
			}
			/*else if (color == 'e') {
				mySimulator.setDensity_R(i, j, 0.0);
				mySimulator.setDensity_G(i, j, 0.0);
				mySimulator.setDensity_B(i, j, 0.0);
			}*/
     
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

    int i = (int)((x / (double)window_width) * mySimulator.getNumCells_I() + 1);
    int j = (int)(((window_height - y) / (double)window_height) * mySimulator.getNumCells_J() + 1);
	
    if (i < 1 || i > mySimulator.getNumCells_I() || j < 1 || j > mySimulator.getNumCells_J())
        return;
        
    if (leftClick) {
        //mySimulator.setDensity(i, j, 100.0);
		if (color == 'r') {
			mySimulator.setDensity_R(i, j, 100.0);
		}
		else if (color == 'g') {
			mySimulator.setDensity_G(i, j, 100.0);
		}
		else if (color == 'b') {
			mySimulator.setDensity_B(i, j, 100.0);
		}

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
    double hx = 1.0 / mySimulator.getNumCells_I();
	double hy = 1.0 / mySimulator.getNumCells_J();

    glColor3f ( 1.0f, 1.0f, 1.0f );
    glLineWidth ( 1.0f );

    glBegin ( GL_LINES );
    for (int i = 0; i <= mySimulator.getNumCells_I(); i++) {
        double x = (i - 0.5) * hx;
        for (int j = 0; j <= mySimulator.getNumCells_J(); j++) {
            double y = (j - 0.5) * hy;

            glVertex2f(x, y);
            //glVertex2f (x + mySimulator.getVelocityU(IX(i,j)), y + mySimulator.getVelocityV(IX(i,j)));
			glVertex2f(x + mySimulator.getVelocityU(i, j), y + mySimulator.getVelocityV(i, j));
        }
    }
    glEnd ();
}

// TODO: Fix
void drawDivergence() {
	double hx = 1.0 / mySimulator.getNumCells_I();
	double hy = 1.0 / mySimulator.getNumCells_J();

	glBegin(GL_QUADS);
	for (int i = 0; i <= mySimulator.getNumCells_I(); i++) {
		double x = (i - 0.5) * hx;
		for (int j = 0; j <= mySimulator.getNumCells_J(); j++) {
			double y = (j - 0.5) * hy;

			double divergence = mySimulator.getDivergence(i, j);

			if (abs(divergence) > 0.01) {
				//glColor3d(1.0, 1.0, 1.0);
				if (divergence < 0) {
					glColor3d(1.0, 0.0, 0.0);
				}
				else if (divergence > 0) {
					glColor3d(0.0, 0.0, 1.0);
				}
			}
			else if (abs(divergence) > 0.005) {
				if (divergence < 0) {
					glColor3d(0.5, 0.0, 0.0);
				}
				else if (divergence > 0) {
					glColor3d(0.0, 0.0, 0.5);
				}
			}
			else {
				glColor3d(0.0, 0.0, 0.0);
			}
			glVertex3f(x, y, 0);
			glVertex3f(x + hx, y, 0);
			glVertex3f(x + hx, y + hy, 0);
			glVertex3f(x, y + hy, 0);

		}
	}
	glEnd();

}

void drawDensity() {
	double hx = 1.0 / mySimulator.getNumCells_I();
	double hy = 1.0 / mySimulator.getNumCells_J();

    glBegin(GL_QUADS);
    for (int i = 0; i <= mySimulator.getNumCells_I(); i++) {
        double x = (i - 0.5) * hx;
        for (int j = 0; j <= mySimulator.getNumCells_J(); j++) {
            double y = (j - 0.5) * hy;

			/*double d00_R = mySimulator.getDensity_R(IX(i, j));
			double d00_G = mySimulator.getDensity_G(IX(i, j));
			double d00_B = mySimulator.getDensity_B(IX(i, j));

			double d01_R = mySimulator.getDensity_R(IX(i, j + 1));
			double d01_G = mySimulator.getDensity_G(IX(i, j + 1));
			double d01_B = mySimulator.getDensity_B(IX(i, j + 1));

			double d10_R = mySimulator.getDensity_R(IX(i + 1, j));
			double d10_G = mySimulator.getDensity_G(IX(i + 1, j));
			double d10_B = mySimulator.getDensity_B(IX(i + 1, j));

			double d11_R = mySimulator.getDensity_R(IX(i + 1, j + 1));
			double d11_G = mySimulator.getDensity_G(IX(i + 1, j + 1));
			double d11_B = mySimulator.getDensity_B(IX(i + 1, j + 1));*/

			double d00_R = mySimulator.getDensity_R(i, j);
			double d00_G = mySimulator.getDensity_G(i, j);
			double d00_B = mySimulator.getDensity_B(i, j);

			double d01_R = mySimulator.getDensity_R(i, j + 1);
			double d01_G = mySimulator.getDensity_G(i, j + 1);
			double d01_B = mySimulator.getDensity_B(i, j + 1);

			double d10_R = mySimulator.getDensity_R(i + 1, j);
			double d10_G = mySimulator.getDensity_G(i + 1, j);
			double d10_B = mySimulator.getDensity_B(i + 1, j);

			double d11_R = mySimulator.getDensity_R(i + 1, j + 1);
			double d11_G = mySimulator.getDensity_G(i + 1, j + 1);
			double d11_B = mySimulator.getDensity_B(i + 1, j + 1);

			glColor3d(d00_R, d00_G, d00_B);
			glVertex3f(x, y, 0);
			glColor3d(d10_R, d10_G, d10_B);
			glVertex3f(x + hx, y, 0);
			glColor3d(d11_R, d11_G, d11_B);
			glVertex3f(x + hx, y + hy, 0);
			glColor3d(d01_R, d01_G, d01_B);
			glVertex3f(x, y + hy, 0);

        }
    }
    glEnd();

}

void initializeFields() {

  for (int i = 1; i <= numCells_I; i++) {
	  for (int j = 1; j <= numCells_J; j++) {

		  if (j <= numCells_J / 3 || j > numCells_J * 7 / 9 && j < numCells_J * 8 / 9) {
			  mySimulator.setDensity_R(i, j, 100.0);
		  }
		  /*if (j >= numCells / 3 && j <= numCells * 2 / 3) {
			  mySimulator.setDensity_G(i, j, 100.0);
		  }*/
		  if (j > numCells_J * 4 / 9 && j < numCells_J * 5 / 9) {
			  mySimulator.setDensity_G(i, j, 100.0);
		  }
		  if (j >= numCells_J * 2 / 3 || j > numCells_J * 1 / 9 && j < numCells_J * 2 / 9) {
			  mySimulator.setDensity_B(i, j, 100.0);
		  }

		  if (i <= j && i + j >= numCells_J) {
			  mySimulator.setV(i, j, 2.0);
		  }
		  if (i >= j && i <= numCells_J - j) {
			  mySimulator.setV(i, j, -2.0);
		  }
		  if (i <= j && i + j <= numCells_I) {
			  mySimulator.setU(i, j, 2.0);
		  }
		  if (i >= j && i + j >= numCells_I) {
			  mySimulator.setU(i, j, -2.0);
		  }
	  }
  }

}

//// loads the red, green, and blue values of the image to the three density fields
//void drawImage(CImage myImage) {
//	for (int i = 0; i < numCells; i++) {
//		for (int j = 0; j < numCells; j++) {
//			myImage.GetPixel(i, j);
//		}
//	}
//}
