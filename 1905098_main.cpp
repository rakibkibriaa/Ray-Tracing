#ifdef __linux__
#include <GL/glut.h>
#elif WIN32
#include <windows.h>
#include <GL/glut.h>
#endif

#include <cmath>
#include <iostream>
#include <bits/stdc++.h>
#define pi (2 * acos(0.0))
#include "1905098_classes.h"
#include "1905098_bitmap_image.hpp"
using namespace std;
/* Global variables */
char title[] = "3D Shapes";

Point center(0, 0, 0);

Point camera(120, 150, 50);

Point forward_vector;

Point right_vector;

Point world_up_vector(0, 0, 1);

Point camera_up_vector;

void axes()
{
    glBegin(GL_LINES);
    glVertex3f(-5, 0, 0);
    glVertex3f(5, 0, 0);
    glVertex3f(0, -5, 0);
    glVertex3f(0, 5, 0);
    glVertex3f(0, 0, -5);
    glVertex3f(0, 0, 5);
    glEnd();
}
int recursion_level;
int img_dimension_pixel;
vector<Object *> objects;
vector<PointLight *> pointLights;
vector<SpotLight *> spotLights;

double viewAngle = 80.0;


void loadData()
{
    ifstream input_file("scene.txt");
    input_file >> recursion_level >> img_dimension_pixel;
    int totalObjects;
    input_file >> totalObjects;

    for (int i = 0; i < totalObjects; i++)
    {
        string obj_type;
        input_file >> obj_type;
        if (obj_type == "sphere")
        {
            double x, y, z, r;
            input_file >> x >> y >> z;
            input_file >> r;
            Sphere *sphere = new Sphere(Point(x, y, z), r);
            input_file >> sphere->color[0] >> sphere->color[1] >> sphere->color[2];
            input_file >> sphere->coEfficients[0] >> sphere->coEfficients[1] >> sphere->coEfficients[2] >> sphere->coEfficients[3];
            input_file >> sphere->shine;

            objects.push_back(sphere);
        }
        else if (obj_type == "triangle")
        {
            double x1, y1, z1, x2, y2, z2, x3, y3, z3;
            input_file >> x1 >> y1 >> z1;
            Point p1(x1, y1, z1);
            input_file >> x2 >> y2 >> z2;
            Point p2(x2, y2, z2);
            input_file >> x3 >> y3 >> z3;
            Point p3(x3, y3, z3);
            Triangle *triangle = new Triangle(p1, p2, p3);

            input_file >> triangle->color[0] >> triangle->color[1] >> triangle->color[2];
            input_file >> triangle->coEfficients[0] >> triangle->coEfficients[1] >> triangle->coEfficients[2] >> triangle->coEfficients[3];
            input_file >> triangle->shine;

            objects.push_back(triangle);
        }
        else if (obj_type == "general")
        {
            double A, B, C, D, E, F, G, H, I, J;
            input_file >> A >> B >> C >> D >> E >> F >> G >> H >> I >> J;
            General *general = new General(A, B, C, D, E, F, G, H, I, J);
            input_file >> general->reference_point.x >> general->reference_point.y >> general->reference_point.z;
            input_file >> general->length >> general->width >> general->height;
            input_file >> general->color[0] >> general->color[1] >> general->color[2];
            input_file >> general->coEfficients[0] >> general->coEfficients[1] >> general->coEfficients[2] >> general->coEfficients[3];
            input_file >> general->shine;

            objects.push_back(general);
        }
    }
    int totalPointLightSources;
    input_file >> totalPointLightSources;
    for (int i = 0; i < totalPointLightSources; i++)
    {
        double x, y, z;
        input_file >> x >> y >> z;
        Point p(x, y, z);
        PointLight *pointLight = new PointLight;
        pointLight->light_pos = p;
        input_file >> pointLight->color[0] >> pointLight->color[1] >> pointLight->color[2];

        pointLights.push_back(pointLight);
    }
    int totalSpotLightSources;
    input_file >> totalSpotLightSources;
    for (int i = 0; i < totalSpotLightSources; i++)
    {
        double x, y, z;

        input_file >> x >> y >> z;
        Point p(x, y, z);

        SpotLight *spotLight = new SpotLight;
        PointLight point_light;

        point_light.light_pos = p;

        spotLight->point_light = point_light;

        input_file >> point_light.color[0] >> point_light.color[1] >> point_light.color[2];

        input_file >> spotLight->light_direction.x >> spotLight->light_direction.y >> spotLight->light_direction.z;

        input_file >> spotLight->cutoff_angle;

        spotLights.push_back(spotLight);


    }
    Object *floor = new Floor(1000, 20);

    floor->coEfficients[0] = 0.4;
    floor->coEfficients[1] = 0.2;
    floor->coEfficients[2] = 0.2;
    floor->coEfficients[3] = 0.2;

    floor->shine = 0.6;

    objects.push_back(floor);
}
int window_height = 480;
int  window_width = 600;

int imageCount = 1;
bitmap_image image;

void capture()
{
	cout<<"Started to capture"<<endl;

    image.set_all_channels(0, 0, 0);

	double planeDistance = (window_height / 2.0) / tan((acos(-1) * viewAngle) / (2 * 180.0));

	Point topLeft = camera + (forward_vector * planeDistance) - (right_vector * (window_width / 2.0)) + (camera_up_vector * (window_height / 2.0));

	double du = 1.0 * window_width / img_dimension_pixel;
	double dv = 1.0 * window_height / img_dimension_pixel;

	// Choose middle of the grid cell
	topLeft = topLeft + (right_vector * du * 0.5) - (camera_up_vector * dv * 0.5);

	int nearestObjectIndex = -1;


    int cnt = 0;
	for(int i=0;i<img_dimension_pixel;i++)
	{
		for(int j=0;j<img_dimension_pixel;j++)
		{


			// calculate current pixel
			Point pixel = topLeft + (right_vector * du * i) - (camera_up_vector * dv * j);

			// cast ray from EYE to (curPixel-eye) direction ; eye is the position of the camera
			Ray ray;
			ray.start = camera;
			Point tt = pixel - camera;

			ray.dir = normalized_point(tt);

			double tmpcolor[3];
			double mini_t = INT_MAX;
			int index = -1;


			for(int k=0;k<objects.size();k++)
			{

				double t = objects[k]->intersect(ray,tmpcolor, 0);

                if(t > 0 && t < mini_t)
                {
                    mini_t = t;
                    index = k;
                }

			}

			// if nearest object is found, then shade the pixel
			if(index != -1)
			{

				double color[3];

				double t = objects[index]->intersect(ray,color, 1);

				if(color[0] > 1) color[0] = 1;
				if(color[1] > 1) color[1] = 1;
				if(color[2] > 1) color[2] = 1;

				if(color[0] < 0) color[0] = 0;
				if(color[1] < 0) color[1] = 0;
				if(color[2] < 0) color[2] = 0;

				image.set_pixel(i, j, 255*color[0],255*color[1],255*color[2]);

			}

		}
	}

	image.save_image("Output_1"+to_string(imageCount)+".bmp");
	cout<<"Capture Finished"<<endl;

	imageCount++;

}

void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW); // To operate on model-view matrix

    // Render a color-cube consisting of 6 quads with different colors
    glLoadIdentity();
    // Reset the model-view matrix
    // glTranslatef(1.5f, 0.0f, -7.0f);  // Move right and into the screen

    gluLookAt(camera.x, camera.y, camera.z,
              camera.x + forward_vector.x, camera.y + forward_vector.y, camera.z + forward_vector.z,
              camera_up_vector.x, camera_up_vector.y, camera_up_vector.z);

    glMatrixMode(GL_MODELVIEW);

    // Begin drawing the color cube with 6 quads
    // Top face (y = 1.0f)
    // Define vertices in counter-clockwise (CCW) order with normal pointing out

    for (int i = 0; i < objects.size(); i++)
    {
        objects[i]->draw();
    }
    for (int i = 0; i < pointLights.size(); i++)
    {
        pointLights[i]->draw();
    }
    for (int i = 0; i < spotLights.size(); i++)
    {
        spotLights[i]->draw();
    }


    glutSwapBuffers(); // Swap the front and back frame buffers (double buffering)
}

void calculateNecessaryVectors()
{
    forward_vector.x = center.x - camera.x;
    forward_vector.y = center.y - camera.y;
    forward_vector.z = center.z - camera.z;

    forward_vector = normalized_point(forward_vector);

    right_vector = cross_product(forward_vector, world_up_vector);
    right_vector = normalized_point(right_vector);

    camera_up_vector = cross_product(right_vector, forward_vector);
    camera_up_vector = normalized_point(camera_up_vector);
}

void reshape(GLsizei width, GLsizei height)
{ // GLsizei for non-negative integer
    // Compute aspect ratio of the new window
    if (height == 0)
        height = 1; // To prevent divide by 0
    GLfloat aspect = (GLfloat)width / (GLfloat)height;

    // Set the viewport to cover the new window
    glViewport(0, 0, width, height);

    // Set the aspect ratio of the clipping volume to match the viewport
    glMatrixMode(GL_PROJECTION); // To operate on the Projection matrix
    glLoadIdentity();            // Reset
    // Enable perspective projection with fovy, aspect, zNear and zFar
    gluPerspective(viewAngle, aspect, 1, 1000.0f);
}

void keyboardListener(unsigned char key, int x, int y)
{
    double v = 0.2;
    double rate = 0.2;
    // Point oldEye = eye;
    float theta = 1 * acos(-1) / 180.0;
    double s;
    // float v = 0.1;
    switch (key)
    {

    case '0':
        capture();
        break;

    case '1':


            right_vector = right_vector * cos(theta) + forward_vector * sin(theta);

            forward_vector = cross_product(camera_up_vector,right_vector);

            forward_vector = normalized_point(forward_vector);

            right_vector = normalized_point(right_vector);

            center = forward_vector * getLength(center - camera);



            break;

        case '2':

            forward_vector = forward_vector * cos(theta) + right_vector * sin(theta);

            forward_vector = normalized_point(forward_vector);

            right_vector = cross_product(forward_vector , camera_up_vector);

            right_vector = normalized_point(right_vector);

            center = forward_vector * getLength(center - camera);

            break;

        case '3':

            forward_vector = forward_vector * cos(theta) + camera_up_vector * sin(theta);

            forward_vector = normalized_point(forward_vector);

            camera_up_vector = cross_product(right_vector, forward_vector);

            camera_up_vector = normalized_point(camera_up_vector);

            center = forward_vector * getLength(center - camera);


            break;

        case '4':

            camera_up_vector = camera_up_vector * cos(theta) + forward_vector * sin(theta);

            camera_up_vector = normalized_point(camera_up_vector);

            forward_vector = cross_product(camera_up_vector, right_vector);

            forward_vector = normalized_point(forward_vector);

            center = forward_vector * getLength(center - camera);

            break;

        case '5':


            camera_up_vector = camera_up_vector * cos(theta) + right_vector * sin(theta);

            camera_up_vector = normalized_point(camera_up_vector);

            right_vector = cross_product(forward_vector,camera_up_vector);

            right_vector = normalized_point(right_vector);

            //center = forward_vector * getLength(center - camera);


            break;

        case '6':

            right_vector = right_vector * cos(theta) + camera_up_vector * sin(theta);

            right_vector = normalized_point(right_vector);

            camera_up_vector = cross_product(right_vector,forward_vector);

            camera_up_vector = normalized_point(camera_up_vector);

            break;

        // Control center (location where the eye is looking at)
        // control centerx
        // Control what is shown
        // Control exit
    case 27:     // ESC key
        exit(0); // Exit window
        break;
    }

    // look = look - eye + oldEye;

    glutPostRedisplay(); // Post a paint request to activate display()
}

void initGL()
{
    calculateNecessaryVectors();
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);              // Set background color to black and opaque
    glClearDepth(1.0f);                                // Set background depth to farthest
    glEnable(GL_DEPTH_TEST);                           // Enable depth testing for z-culling
    glDepthFunc(GL_LEQUAL);                            // Set the type of depth-test
    glShadeModel(GL_SMOOTH);                           // Enable smooth shading
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST); // Nice perspective corrections

    loadData();

    image = bitmap_image(img_dimension_pixel,img_dimension_pixel);
}

void specialKeyListener(int key, int x, int y)
{
    double r = 1;
    switch (key)
    {
    case GLUT_KEY_UP: // down arrow key
        // pos = pos + l;
        camera.x = camera.x + r * forward_vector.x;
        camera.y = camera.y + r * forward_vector.y;
        camera.z = camera.z + r * forward_vector.z;

        center.x = center.x + r * forward_vector.x;
        center.y = center.y + r * forward_vector.y;
        center.z = center.z + r * forward_vector.z;

        break;
    case GLUT_KEY_DOWN: // up arrow key
                        // pos = pos - l;
        camera.x = camera.x - r * forward_vector.x;
        camera.y = camera.y - r * forward_vector.y;
        camera.z = camera.z - r * forward_vector.z;

        center.x = center.x - r * forward_vector.x;
        center.y = center.y - r * forward_vector.y;
        center.z = center.z - r * forward_vector.z;

        break;

    case GLUT_KEY_RIGHT:
        // pos = pos + r;
        camera.x = camera.x + r * right_vector.x;
        camera.y = camera.y + r * right_vector.y;
        camera.z = camera.z + r * right_vector.z;

        center.x = center.x + r * right_vector.x;
        center.y = center.y + r * right_vector.y;
        center.z = center.z + r * right_vector.z;

        break;
    case GLUT_KEY_LEFT:
        // pos = pos - r;
        camera.x = camera.x - r * right_vector.x;
        camera.y = camera.y - r * right_vector.y;
        camera.z = camera.z - r * right_vector.z;

        center.x = center.x - r * right_vector.x;
        center.y = center.y - r * right_vector.y;
        center.z = center.z - r * right_vector.z;
        break;

    case GLUT_KEY_PAGE_UP:
        // pos = pos + u;
        camera.x = camera.x + r * camera_up_vector.x;
        camera.y = camera.y + r * camera_up_vector.y;
        camera.z = camera.z + r * camera_up_vector.z;

        center.x = center.x + r * camera_up_vector.x;
        center.y = center.y + r * camera_up_vector.y;
        center.z = center.z + r * camera_up_vector.z;
        break;
    case GLUT_KEY_PAGE_DOWN:
        camera.x = camera.x - r * camera_up_vector.x;
        camera.y = camera.y - r * camera_up_vector.y;
        camera.z = camera.z - r * camera_up_vector.z;

        center.x = center.x - r * camera_up_vector.x;
        center.y = center.y - r * camera_up_vector.y;
        center.z = center.z - r * camera_up_vector.z;
        break;

    default:
        break;
    }
    glutPostRedisplay();
}

/* Main function: GLUT runs as a console application starting at main() */
int main(int argc, char **argv)
{
    glutInit(&argc, argv); // Initialize GLUT
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(window_width, window_height); // Set the window's initial width & height
    glutInitWindowPosition(0, 0); // Position the window's initial top-left corner
    glutCreateWindow(title);      // Create window with the given title
    glutDisplayFunc(display);     // Register callback handler for window re-paint event
    glutReshapeFunc(reshape);     // Register callback handler for window re-size event
    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    initGL();       // Our own OpenGL initialization
    glutMainLoop(); // Enter the infinite event-processing loop

    objects.clear();
    objects.shrink_to_fit();

    pointLights.clear();
    pointLights.shrink_to_fit();

    spotLights.clear();
    spotLights.shrink_to_fit();
    return 0;
}
