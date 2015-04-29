#include <stdio.h>
#include <stdlib.h>
#include <gl\glut.h>
#include "udray.h"
#include "glm.h"

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

extern Camera *ray_cam;       // camera info
extern int image_i, image_j;  // current pixel being shaded
extern bool wrote_image;      // has the last pixel been shaded?

// reflection/refraction recursion control

extern int maxlevel;          // maximum depth of ray recursion 
extern double minweight;      // minimum fractional contribution to color

// these describe the scene

extern vector < GLMmodel * > model_list;
extern vector < Sphere * > sphere_list;
extern vector < Light * > light_list;

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

// intersect a ray with the entire scene (.obj models + spheres)

// x, y are in pixel coordinates with (0, 0) the upper-left hand corner of the image.
// color variable is result of this function--it carries back info on how to draw the pixel

void trace_ray(int level, double weight, Ray *ray, Vect color)
{
	Intersection *nearest_inter = NULL;
	Intersection *inter = NULL;
	int i;

	// test for intersection with all .obj models

	for (i = 0; i < model_list.size(); i++) {
		inter = intersect_ray_glm_object(ray, model_list[i]);
		update_nearest_intersection(&inter, &nearest_inter);
	}

	// test for intersection with all spheres

	for (i = 0; i < sphere_list.size(); i++) {
		inter = intersect_ray_sphere(ray, sphere_list[i]);
		update_nearest_intersection(&inter, &nearest_inter);
	}

	// "color" the ray according to intersecting surface properties

	// choose one of the simpler options below to debug or preview your scene more quickly.
	// another way to render faster is to decrease the image size.

	if (nearest_inter) {
		//shade_ray_false_color_normal(nearest_inter, color);
		//    shade_ray_intersection_mask(color);  
		shade_ray_diffuse(ray, nearest_inter, color);
		//   shade_ray_recursive(level, weight, ray, nearest_inter, color);
	}

	// color the ray using a default

	else
		shade_ray_background(ray, color); 
}

//----------------------------------------------------------------------------

// test for ray-sphere intersection; return details of intersection if true

#define SMALL_NUM  0.00000001 // anything that avoids division overflow
#define rayeps 1e-7
Intersection *intersect_ray_sphere(Ray *ray, Sphere *S)
{
	// FILL IN CODE (line below says "no" for all spheres, so replace it)
	Vect    w0, w;          // ray vectors
	float   a, b, c;           // params to calc ray-sphere intersect
	float t;
	Vect I;
	Intersection *inter;
	//float rayDeltaP = *(ray)->dir;
	Vect deltaP;
	//VectSub(S->P, ray->orig, deltaP);
	//VectSub(ray->orig, S->P, deltaP);
	VectSub(S->P, ray->orig, deltaP);
	//float mag = VectUnit(deltaP);
	//a = VectDotProd(deltaP, deltaP);
	//a = 1;
	//b = 2 * *(ray)->dir  
	b = VectDotProd(ray->dir, deltaP);
	//c = VectDotProd(deltaP, deltaP) - pow(S->radius, 2);
	c = VectDotProd(deltaP, deltaP) - pow(S->radius, 2);
	//c = VectDotProd(deltaP, deltaP) - (pow(mag, 2) - pow(S->radius, 2));
	//c = VectUnit(deltaP) - pow(S->radius, 2);
	float b2c = sqrt(pow(b, 2) - c); 
	float t01a = -b + b2c;
	float t02a = -b - b2c;
	float t01 = b + b2c;
	float t02 = b - b2c;
	float b2_c = pow(b, 2) - c;
	if(t01 < 0 && t02 < 0)
		return NULL;
	if(b2_c < 0)
		return NULL;
	else if(b2_c > 0){
		if(t01 < t02)
			t = t01;
		else
			t = t02;
	}
	else // b^2 - c = 0
		if(t01 = t02)
			t = t01;
		else
			t = (t01 + t02) /2;
	/**
	//float mag = VectUnit(b);
	float rayDeltaP = VectDotProd(ray->dir, deltaP);  //ray->dir * deltaP
	float tout = rayDeltaP; //ray->dir * deltaP;
	float tin = sqrt(pow((rayDeltaP), 2) - (pow(fabs(VectMag(deltaP)), 2) - pow(S->radius, 2)));
	// get triangle edge vectors and plane normal

	float t0 = tout + tin;
	float t1 = tout - tin;
	
	if((pow(b,2) - c) < 0)
		return NULL;
	else if((pow(b,2) - c) > 0){
		if(t01 < t02)
			t = t01;
		else
			t = t02;
	}
	else // b^2 - c = 0
		if(t01 = t02)
			t = t01;
		else
			t = (t01 + t02) /2;
	

	**/
	// get intersect point of ray with triangle plane

	//if (t < rayeps)                   // triangle is behind/too close to ray => no intersect
	//	return NULL;                 // for a segment, also test if (t > 1.0) => no intersect

	// intersect point of ray and plane

	VectAddS(t, ray->dir, ray->orig, I);        

	inter = make_intersection();
	inter->t = t;
	inter->N;
	inter->surf = S->surf;
	//inter->N = deltaP; // mag;
	VectCopy(inter->P, I);
	VectUnit(I);
	return inter;                      // I is in T
	//return NULL;
}

//----------------------------------------------------------------------------

// only local, ambient + diffuse lighting (no specular, shadows, reflections, or refractions)

void shade_ray_diffuse(Ray *ray, Intersection *inter, Vect color)
{
	Vect L;
	double diff_factor;

	// iterate over lights

	for (int i = 0; i < light_list.size(); i++) {

		// AMBIENT

		color[R] += inter->surf->amb[R] * light_list[i]->amb[R];
		color[G] += inter->surf->amb[G] * light_list[i]->amb[G];
		color[B] += inter->surf->amb[B] * light_list[i]->amb[B];

		// DIFFUSE 

		// FILL IN CODE
		Vect v;
		Vect N;
		VectCopy(N, inter->N);
		VectUnit(N);
		Vect P;
		VectCopy(P, light_list[i]->P);
		VectUnit(P);
		double normal = VectDotProd(N, P);
		VectAddS(1, inter->N, light_list[i]->P, v);
		//double normal = VectUnit(v);
		color[R] += inter->surf->diff[R] * light_list[i]->diff[R] * normal;//cos(
		color[G] += inter->surf->diff[G] * light_list[i]->diff[G] * normal;
		color[B] += inter->surf->diff[B] * light_list[i]->diff[B] * normal;
	}

	// clamp color to [0, 1]

	VectClamp(color, 0, 1);
}

//----------------------------------------------------------------------------

// same as shade_ray_diffuse(), but add specular lighting + shadow rays (i.e., full Phong illumination model)

void shade_ray_local(Ray *ray, Intersection *inter, Vect color)
{
	Vect L;
	double diff_factor;

	// iterate over lights

	for (int i = 0; i < light_list.size(); i++) {
	// FILL IN CODE
		Vect v;
		VectAddS(1, inter->N, light_list[i]->P, v);
		double normal = VectUnit(v);
		double shiny = inter->surf->spec_exp;
		color[R] += inter->surf->spec[R] * light_list[i]->spec[R] * pow(normal, shiny) ;//cos(
		color[G] += inter->surf->spec[G] * light_list[i]->spec[G] * pow(normal, shiny);
		color[B] += inter->surf->spec[B] * light_list[i]->spec[B] * pow(normal, shiny);
	}
}

//----------------------------------------------------------------------------

// full shading model: ambient/diffuse/specular lighting, shadow rays, recursion for reflection, refraction

// level = recursion level (only used for reflection/refraction)

void shade_ray_recursive(int level, double weight, Ray *ray, Intersection *inter, Vect color)
{
	Surface *surf;
	int i;

	// initialize color to Phong reflectance model

	shade_ray_local(ray, inter, color);

	// if not too deep, recurse

	if (level + 1 < maxlevel) {

		// add reflection component to color

		if (surf->reflectivity * weight > minweight) {

			// FILL IN CODE

		}

		// add refraction component to color

		if (surf->transparency * weight > minweight) {

			// GRAD STUDENTS -- FILL IN CODE

		}
	}
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

// ray trace another pixel if the image isn't finished yet

void idle()
{
	if (image_j < ray_cam->im->h) {

		raytrace_one_pixel(image_i, image_j);

		image_i++;

		if (image_i == ray_cam->im->w) {
			image_i = 0;
			image_j++;
		}    
	}

	// write rendered image to file when done

	else if (!wrote_image) {

		write_PPM("output.ppm", ray_cam->im);

		wrote_image = true;
	}

	glutPostRedisplay();
}

//----------------------------------------------------------------------------

// show the image so far

void display(void)
{
	// draw it!

	glPixelZoom(1, -1);
	glRasterPos2i(0, ray_cam->im->h);

	glDrawPixels(ray_cam->im->w, ray_cam->im->h, GL_RGBA, GL_FLOAT, ray_cam->im->data);

	glFlush ();
}

//----------------------------------------------------------------------------

void init()
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, ray_cam->im->w, 0.0, ray_cam->im->h);
}

//----------------------------------------------------------------------------

int main(int argc, char** argv)
{
	glutInit(&argc, argv);

	// initialize scene (must be done before scene file is parsed)

	init_raytracing();

	if (argc == 2)
		parse_scene_file(argv[1], ray_cam);
	else {
		printf("missing .scene file\n");
		exit(1);
	}

	// opengl business

	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(ray_cam->im->w, ray_cam->im->h);
	glutInitWindowPosition(500, 300);
	glutCreateWindow("hw3");
	init();

	glutDisplayFunc(display); 
	glutIdleFunc(idle);

	glutMainLoop();

	return 0; 
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
