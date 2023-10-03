#include <GL/glut.h>

#include <iostream>

#include "error.hpp"
#include "parameters.hpp"

#include "fluid.hpp"

int mouseXlast = 0;
int mouseYlast = 0;

FluidParticle fluid = FluidParticle(parameters.N, parameters.dt, parameters.diff, parameters.visc);

void display()
{
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f); CHECK_GL_ERROR();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); CHECK_GL_ERROR();
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); CHECK_GL_ERROR();
    glEnable(GL_BLEND); CHECK_GL_ERROR();

    fluid.step();

    for (int i = 0; i < fluid.getSize(); ++i)
    {
        for (int j = 0; j < fluid.getSize(); ++j)
        {
            int x = i * parameters.scale;
            int y = j * parameters.scale;
            float d = fluid.getDensity()[INDEX(i, j)];

            glColor4f(d, d, d, 1.0f);
            
            glBegin(GL_QUADS);
            glVertex2f(x, y);
            glVertex2f(x + parameters.scale, y);
            glVertex2f(x + parameters.scale, y + parameters.scale);
            glVertex2f(x, y + parameters.scale);
            glEnd();
        }
    }

    fluid.fadeDensity();

    glutSwapBuffers(); CHECK_GL_ERROR();
}

void reshape(int width, int height)
{
    height = std::max(1, height);
	width = std::max(1, width);

	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0, width, 0, height);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void mouse_move(int x, int y)
{
    int x_diff = x - mouseXlast;
	int y_diff = mouseYlast - y;

	int x_vel = std::min(100 * x_diff, 1000);
	int y_vel = std::min(100 * y_diff, 1000);

	float amountX = std::max(x_vel, -1000);
	float amountY = std::max(y_vel, -1000);

	int x_actual = x / parameters.scale;
	int y_actual = parameters.N - (y / parameters.scale);

	fluid.addDensity(x_actual, y_actual, 40);
	fluid.addVelocity(x_actual, y_actual, amountX, amountY);

    mouseXlast = x;
	mouseYlast = y;
	glutPostRedisplay();
}

void mouse_click(int button, int state, int x, int y)
{
	mouseXlast = x;
	mouseYlast = y;
}

void timer(int value) {
	glutPostRedisplay();
	glutTimerFunc(value, timer, 0);
}

int main(int argc, char** argv)
{
    glutInit(&argc, argv); CHECK_GL_ERROR();
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH); CHECK_GL_ERROR();
    glutInitWindowSize(parameters.window_width, parameters.window_height); CHECK_GL_ERROR();
    glutInitWindowPosition(10, 10); CHECK_GL_ERROR();
    glutCreateWindow("ANI3D - Fluid Simulation"); CHECK_GL_ERROR();
    glutDisplayFunc(display); CHECK_GL_ERROR();
    glutReshapeFunc(reshape); CHECK_GL_ERROR();
    glutTimerFunc(10, timer, 0); CHECK_GL_ERROR();
    glutMotionFunc(mouse_move); CHECK_GL_ERROR();
    glutMouseFunc(mouse_click); CHECK_GL_ERROR();

    glutMainLoop();

    return EXIT_SUCCESS;
}