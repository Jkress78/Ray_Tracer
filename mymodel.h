/*
* James Kress - 200429558
* CS 405 Term Project - RayTracing with reflections, refraction, and shadows
* April 20, 2023
*
* mymodel.h: contains all of the hard coded data for "rayTrace.cpp" as well as all of the function definitions.
*
*
*/

#pragma once
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>


using namespace std;

/*
* Super Class for objects to be rendered in a scene
*/
class Object {

	float kd; //diffuse coefficient 
	//percentage of shading for current point, reflected point, and refracted point
	float w1, w2, w3; // w1 + w2 + w3 = 1 
	int type;

public:


	Object(float nKd, float nW1, float nW2, float nW3, int nType) {
		this->kd = nKd;
		this->w1 = nW1;
		this->w2 = nW2;
		this->w3 = nW3;
		this->type = nType;
	}

	float getKd() { return this->kd; }

	float getW1() { return this->w1; }

	float getW2() { return this->w2; }

	float getW3() { return this->w3; }

	int getType() { return this->type; }

};

class Poly4 : public Object {
public:
	float v[4][3]; //list of vertices
	float n[3]; //normal of polygon

	Poly4(float nKd, float nW1, float nW2, float nW3, float nV[4][3], float nN[3], int nType) :Object(nKd, nW1, nW2, nW3, nType) {

		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 3; j++)
				this->v[i][j] = nV[i][j];
		}

		for (int i = 0; i < 3; i++)
			this->n[i] = nN[i];
	}


};

class Sphere : public Object {
private:
	float x, y, z;
	float rad;

public:
	Sphere(float nKd, float nW1, float nW2, float nW3, float nX, float nY, float nZ, float nRad, int nType)
		:Object(nKd, nW1, nW2, nW3, nType) {
		this->x = nX;
		this->y = nY;
		this->z = nZ;
		this->rad = nRad;
	}

	float getX() { return this->x; }
	float getY() { return this->y; }
	float getZ() { return this->z; }
	float getRad() { return this->rad; }

};

/* definition of the image buffer */
#define ROWS 512
#define COLS 512

unsigned char img[ROWS][COLS];

unsigned char img2[ROWS][COLS];

/* definition of window on the image plane in the camera coordinates */
/* They are used in mapping (j, i) in the screen coordinates into */
/* (x, y) on the image plane in the camera coordinates */
/* The window size used here simulates the 35 mm film. */
float xmin = 0.0175;
float ymin = -0.0175;
float xmax = -0.0175;
float ymax = 0.0175;


/* definition of the camera parameters */
float VRP[3] = { 1.0, 3.5, 6.5 };
float VPN[3] = { 0.0, -1, -2.5 };
float VUP[3] = { 0.0, 1.0, 0.0 };

float focal = 0.025;	/* focal length simulating 50 mm lens */


/* definition of light source */
float LPR[3] = { -10.0, 10.0, -2.0 };	/* light position */
float Ip = 300.0;	/* intensity of the point light source */

/* === transformation matrices (to be constructed) === */

/* Transformation from the world to the camera coordinates */
float Mwc[4][4];

/* Transformation from the camera to the world coordinates */
float Mcw[4][4];


//Vectors to store all objects used in the scene
vector<Poly4> poly_list; 
vector<Sphere> sphere_list;

/*-------- Functions --------*/
void rayConstruction(int i, int j, float(&p0)[3], float(&v0)[3]);
int rayTracing(float p0[3], float v[3], int level);
int rayObjIntrsct(float p0[3], float v0[3], float(&p)[3], float(&n)[3], float& kd, float& w1, float& w2, float& w3);
void reflectionRay(float v[3], float n[3], float(&r1)[3]);
void refractionRay(float v[3], float n[3], float(&r1)[3]);
bool isInShadow(float p[3]);
int shading(float p[3], float n[3], float kd);
float magnitude(float v[3]);
float dotProd(float v1[3], float v2[3]);
float raySphereIntersect(float p0[3], float v0[3], Sphere s, float(&n)[3]);
float rayPolyIntersect(float p0[3], float v0[3], Poly4 pl);
bool isInside(float p[], float v[][3], float n[], int nVert);
vector<vector<float>> copyVl(float vl[][3], int nVert, int index);
void copyP(float p[], float(&p2)[2], int index);
vector<float> findnuv(float v1[], float v2[], float p[]);
void findmWCmCW(float vpn[], float vup[], float vrp[]);

