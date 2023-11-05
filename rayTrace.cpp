/*
* James Kress - 200429558
* CS 405 Term Project - RayTracing with reflections, refraction, and shadows
* April 20, 2023
*
* rayTrace.cpp: takes data from hard coded objects in "mymodel.h" and performs ray tracing to generate a
* raw image and store it as a file which can be converted to a proper image using online resources.
*
*
*/

#include "mymodel.h"
//using namespace std;

int main() {

	findmWCmCW(VPN, VUP, VRP);
	int numShapes;
	std::cout << "How many shapes would you like to make?" << endl;
	std::cout << "# of Shapes: ";
	std::cin >> numShapes;

	char choice;
	while (numShapes > 0) {
		std::cout << "What kind of shape do you want to add?" << endl;
		std::cout << "p: polygon plane" << endl;
		std::cout << "s: sphere" << endl;
		std::cout << "Choice: ";
		std::cin >> choice;

		if (choice == 'p') {
			std::cout << "Enter the choords of the points of the plane." << endl;
			float tempV[4][3];
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 3; j++) {
					std::cout << "Point: " << i << endl;
					if (j == 0)
						std::cout << "X: ";
					else if (j == 1)
						std::cout << "Y: ";
					else
						std::cout << "Z: ";

					std::cin >> tempV[i][j];
				}
			}

			float tempN[3];
			std::cout << "Enter the values for the normal vector.";
			std::cout << "X: ";
			std::cin >> tempN[0];
			std::cout << "Y: ";
			std::cin >> tempN[1];
			std::cout << "Z: ";
			std::cin >> tempN[2];

			float tempKd;
			std::cout << "Enter the value of the defuse coefficient (between 0-1): ";
			std::cin >> tempKd;

			float tempW1, tempW2, tempW3;
			std::cout << "Now enter the values for the percentage of shading calculated at\n";
			std::cout << "the current point(w1), the reflected point(w2), the refracted point(w3).\n";
			std::cout << "Ensure that w1 + w2 + w3 = 1.\n";

			std::cout << "w1: ";
			std::cin >> tempW1;

			std::cout << "w2: ";
			std::cin >> tempW2;

			std::cout << "w3: ";
			std::cin >> tempW3;

			poly_list.push_back(Poly4(tempKd, tempW1, tempW2, tempW3, tempV, tempN, 0));
		}

		else if (choice == 's') {
			float t1, t2, t3;
			std::cout << "Enter the coordinates of the center of the Sphere.\n";
			std::cout << "X: ";
			std::cin >> t1;
			std::cout << "Y: ";
			std::cin >> t2;
			std::cout << "Z: ";
			std::cin >> t3;

			float tempR;
			std::cout << "Enter the radius size.\n";
			std::cout << "Radius: ";
			std::cin >> tempR;

			float tempKd;
			std::cout << "Enter the value of the defuse coefficient (between 0-1): ";
			std::cin >> tempKd;

			float tempW1, tempW2, tempW3;
			std::cout << "Now enter the values for the percentage of shading calculated at\n";
			std::cout << "the current point(w1), the reflected point(w2), the refracted point(w3).\n";
			std::cout << "Ensure that w1 + w2 + w3 = 1.\n";

			std::cout << "w1: ";
			std::cin >> tempW1;

			std::cout << "w2: ";
			std::cin >> tempW2;

			std::cout << "w3: ";
			std::cin >> tempW3;

			sphere_list.push_back(Sphere(tempKd, tempW1, tempW2, tempW3, t1, t2, t3, tempR, 1));
		}

		else {
			cout << "Invalid decision try again.\n\n";
			numShapes++;
		}

		numShapes--;
	}

	

	float P0[3] = { VRP[0], VRP[1], VRP[2] }; //camera VRP (origin of the rays)
	float V[3]; //var for generated rays
	int C;	//shading value
	int level = 2;
	for (int i = 0; i < ROWS; i++) {
		for (int j = 0; j < COLS; j++) {
			//construct Ray
			rayConstruction(i, j, P0, V);

			if ((C = rayTracing(P0, V, level)) != NULL) {
				img[i][j] = C;
			}
		}
	}
	//output the generated image to a raw image file
	//img buffer is being output from bottom up so the result is not upside down as windows 
	//generates imgs bottom up
	ofstream binaryFile("new2.raw", ios::out | ios::binary);
	for (int i = ROWS-1; i >= 0; i--) {
		for (int j = 0; j < COLS; j++) {
			binaryFile.write((char*)&img[i][j], sizeof(unsigned char));
		}
	}

	
	return 0;
}

void rayConstruction(int i, int j, float(&p0)[3], float(&v0)[3]) {
	//Convert point in image buffer to camera coords
	float x = (xmax - xmin) * j / (COLS - 1) + xmin;
	float y = (ymax - ymin) * i / (ROWS - 1) + ymin;


	//add focal length to make the converted image buffer coords a 3D point on the image plane
	float P1[4] = { x, y, focal, 1 };

	//Convert image plane coord into world coords storing in a temp variable
	//Temp var is in homogenious coords
	float temp[4] = { 0, 0, 0, 0 };
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			temp[i] += Mcw[i][j] * P1[j];
			
		}
		
	}

	//take converted point from temp and place into P1
	//also convert the coords in temp back into regular 3D coords from homogeneous

	for (int i = 0; i < 3; i++) {
		P1[i] = temp[i] / temp[3];
	}


	//Calculate vector from P0 (VPN) to P1 (point on image plane)
	for (int i = 0; i < 3; i++) {
		v0[i] = P1[i] - p0[i];
	}

	//get magnitude of new vector
	float vMag = magnitude(v0);

	//Normalize the new vector

	for (int i = 0; i < 3; i++) {
		v0[i] = v0[i] / vMag;
	}
	

}

int rayTracing(float p0[3], float v[3], int level) {
	float p[3]; //collision point
	float n[3];	//surface normal for collision point
	float kd;   //diffuse coefficient for collision point
	float w1, w2, w3;

	int found = rayObjIntrsct(p0, v, p, n, kd, w1, w2, w3);

	if (found == NULL)
		return 61;
	else {
		if (isInShadow(p)) //check for shadows, if there is skip further calcs and return a shading value of 0
			return 0;
		else {
			int c1 = 0, c2 = 0, c3 = 0;//holds shading value to be returned
			
			//Basic diffuse shading
			if (w1 > 0.0)
				c1 = shading(p, n, kd);
			//Reflection calculation
			if (w2 > 0.0 && level > 1) {
				float r1[3];
				reflectionRay(v, n, r1);
				c2 = rayTracing(p, r1, level - 1);
			}
			//refraction calculation
			if (w3 > 0 && level > 1) {
				float r2[3];
				refractionRay(v, n, r2);
			}

			return (c1 * w1 + c2 * w2 + c3 * w3);
		}
	}
}

bool isInShadow(float p[3]) {
	//Create a shadow ray from the intersection point to the light source
	//and normalize it
	float shad[3];
	for (int i = 0; i < 3; i++) {
		shad[i] = LPR[i] - p[i];
	}
	
	float vecMag = magnitude(shad);

	for (int i = 0; i < 3; i++)
		shad[i] = shad[i] / vecMag;
	//empty values used to fill in the intersection functions
	float n[3];

	//check for object intersection along the shadow ray
	float t1 = NULL;
	for (int i = 0; i < poly_list.size(); i++) {
		t1 = rayPolyIntersect(p, shad, poly_list.at(i));
		if (t1 != NULL)
			break;
	}
	float t2 = NULL;
	for (int i = 0; i < sphere_list.size(); i++) {
		t2 = raySphereIntersect(p, shad, sphere_list.at(i), n);
		if (t2 != NULL)
			break;
	}
	float epsilon = 0.1; //used so the ray does not intersect with the object where the original intersection occured

	/* if no intersections or the intersections found are less than epsilon then no shadow */
	if (t1 == NULL && t2 == NULL)
		return false;
	else if (t1 < epsilon && t2 < epsilon)
		return false;
	else 
		return true;
}

void reflectionRay(float v[3], float n[3], float(&r1)[3]) {
	/* Calculate the refelction ray using the incident ray V, surface normal N, and the 
		angle between the two */
	r1[0] = v[0] - 2 * dotProd(v, n) * n[0];
	r1[1] = v[1] - 2 * dotProd(v, n) * n[1];
	r1[2] = v[2] - 2 * dotProd(v, n) * n[2];

}

void refractionRay(float v[3], float n[3], float(&r2)[3]) {
	/* solve for the refraction ray equation:
	* R = (eta_r*(N dot V) - sqrt(1 - eta_r*eta_r * (1 - (N dot V)*(N dot V))) - N - eta_r*V
	* where eta_r is equal to the refraction index of the current medium / the refraction index of the next medium
	* V-> the current ray
	* N-> surface normal
	* sqrt-> square root
	*/
	float eta_t1 = 1.333; //water
	float eta_t2 = 1.52; //glass
	float eta_i = 1.0; //air

	float eta_r = eta_i / eta_t2;

	float cos_theta_i = dotProd(n, v);

	if (cos_theta_i < 0)
		cos_theta_i = -cos_theta_i;

	float root = sqrt(1 - eta_r * eta_r * (1 - cos_theta_i * cos_theta_i));

	r2[0] = (eta_r * cos_theta_i - root) * n[0] - eta_r * v[0];
	r2[1] = (eta_r * cos_theta_i - root) * n[0] - eta_r * v[1];
	r2[2] = (eta_r * cos_theta_i - root) * n[0] - eta_r * v[2];
	
}



int rayObjIntrsct(float p0[3], float v0[3], float(&p)[3], float(&n)[3], float& kd, float& w1, float& w2, float& w3) {
	//Check each obj list to see if empty
	bool is_poly_empty = poly_list.empty();
	bool is_sphere_empty = sphere_list.empty();

	int index_p;
	int index_s;
	float t_old_p = NULL;
	float t_new_p;
	float t_old_s = NULL;
	float t_new_s;
	float n_s[3];

	for (int i = 0; i < poly_list.size(); i++) {
		t_new_p = rayPolyIntersect(p0, v0, poly_list.at(i));
		if (t_new_p == NULL) continue;
		if (t_old_p == NULL || t_old_p > t_new_p) { t_old_p = t_new_p; index_p = i; }
	}

	for (int i = 0; i < sphere_list.size(); i++) {
		t_new_s = raySphereIntersect(p0, v0, sphere_list.at(i), n_s);
		if (t_new_s == NULL) continue;
		if (t_old_s == NULL || t_old_s > t_new_s) { t_old_s = t_new_s; index_s = i; }
	}

	if (t_old_s == NULL && t_old_p == NULL)
		return NULL;

	else if (t_old_s == NULL) {
		for (int i = 0; i < 3; i++) {
			p[i] = p0[i] + t_old_p * v0[i];

			n[i] = poly_list[index_p].n[i];
		}

		kd = poly_list[index_p].getKd();
		w1 = poly_list[index_p].getW1();
		w2 = poly_list[index_p].getW2();
		w3 = poly_list[index_p].getW3();
		return 1;
	}

	else if (t_old_p == NULL) {
		for (int i = 0; i < 3; i++) {
			p[i] = p0[i] + t_old_s * v0[i];

			n[i] = n_s[i];
		}

		kd = sphere_list[index_s].getKd();
		w1 = sphere_list[index_s].getW1();
		w2 = sphere_list[index_s].getW2();
		w3 = sphere_list[index_s].getW3();
		return 1;
	}

	else if (t_old_p < t_old_s) {
		for (int i = 0; i < 3; i++) {
			p[i] = p0[i] + t_old_p * v0[i];

			n[i] = poly_list[index_p].n[i];
		}

		kd = poly_list[index_p].getKd();
		w1 = poly_list[index_p].getW1();
		w2 = poly_list[index_p].getW2();
		w3 = poly_list[index_p].getW3();
		return 1;
	}

	else {
		for (int i = 0; i < 3; i++) {
			p[i] = p0[i] + t_old_s * v0[i];

			n[i] = n_s[i];
		}

		kd = sphere_list[index_s].getKd();
		w1 = sphere_list[index_s].getW1();
		w2 = sphere_list[index_s].getW2();
		w3 = sphere_list[index_s].getW3();
		return 1;
	}
}

float raySphereIntersect(float p0[3], float v0[3], Sphere s, float(&n)[3]) {

	float c[3] = { s.getX(), s.getY(), s.getZ()}; //center of sphere

	float L[3];      //vector from camera origin to center of sphere c
	for (int i = 0; i < 3; i++) { //get vector from O to C
		L[i] = c[i] - p0[i];
	}
	float tca = dotProd(L, v0);

	float lMag = magnitude(L);
	float d = sqrt((lMag * lMag) - (tca * tca));
	//check if distance from the center to tca is greater than radius
	if (d > s.getRad())
		return NULL;

	else {
		float thc = sqrt((s.getRad() * s.getRad()) - (d * d));
		float t = tca + thc; // distance from camera origin to intersection point

		//get the normal of the sphere
		for (int i = 0; i < 3; i++) {
			n[i] = c[i] - (p0[i] + t * v0[i]);
		}
		//Normalize the normal
		float nMag = magnitude(n);
		for (int i = 0; i < 3; i++)
			n[i] = n[i] / nMag;
		

		return t;
	}



}

float rayPolyIntersect(float p0[3], float v0[3], Poly4 pl) {

	//get the normal of the polygon from 2 sets of vertices defined by l1 and l2
	float l1[3];
	float l2[3];

	for (int i = 0; i < 3; i++) {
		l1[i] = pl.v[1][i] - pl.v[0][i]; //l1 = v1-v0
		l2[i] = pl.v[2][i] - pl.v[0][i]; //l2 = v2-v0
	}

	//cross product between vectors l1 and l2 gives the normal of the polygon
	pl.n[0] = l1[1] * l2[2] - l1[2] * l2[1];
	pl.n[1] = l1[2] * l2[0] - l1[0] * l2[2];
	pl.n[2] = l1[1] * l2[0] - l1[0] * l2[1];

	//normalize the normal
	float nMag = magnitude(pl.n);
	for (int i = 0; i < 3; i++)
		pl.n[i] = pl.n[i] / nMag;

	//Now we find the value D with the values of ABC given by the normal
	// A = Nx, B = Ny, C = Nz
	//D = -(A*v0x + B*v0y + C*v0z) 
	//use vertex 1 in pl.v to calculate D
	float D = 0 - ((pl.n[0] * pl.v[1][0]) + (pl.n[1] * pl.v[1][1]) + (pl.n[2] * pl.v[1][2]));

	//next check for intersection with the plane
	float dot1 = dotProd(pl.n, v0);
	if (dot1 == 0)
		return NULL;

	else {
		//get the distance from the point O to the plane intersection
		float t =  -(dotProd(pl.n, p0) + D) / dot1;

		//get the coords of the intersection point on the plane
		float p[3];
		for (int i = 0; i < 3; i++)
			p[i] = p0[i] + t * v0[i];
		//get number of vertices in the polygon to be tested
		int numVerts = sizeof(pl.v) / sizeof(pl.v[0]);

		//check if the intersection point is in the polygon or not
		bool inside = isInside(p, pl.v, pl.n, numVerts);

		//if the point is inside the polygon set the value of kd
		//then return t
		if (inside) { return t; }
		else { return NULL; }
	}

}

bool isInside(float p[], float v[][3], float n[], int nVert) {
	int index;
	float p2[2]; //holds 2D coords of p
	int counter = 0; //counts how many times an edge of a polygon crosses the +x-axis

	//check which coord to drop 
	if (n[0] > n[1] && n[0] > n[2])
		index = 0;
	else if (n[1] > n[0] && n[1] > n[2])
		index = 1;
	else
		index = 2;

	//copy the vertices of the polygon to new variable 
	//this also drops the x, y, or z coord specified by the index
	//converting it from 3D to 2D
	vector<vector<float>> vl2 = copyVl(v, nVert, index);

	//convert the point P from 3D to 2D and store in p2
	copyP(p, p2, index);

	//translate the 2D polygon so p2 is at the origin
	for (int i = 0; i < nVert; i++)
		for (int j = 0; j < 2; j++)
			vl2[i][j] -= p2[j];

	int j; //used to roll index back to 0 when needed
	for (int i = 0; i < nVert; i++) {
		j = i + 1;
		if (i == nVert - 1) //used to get v3-v0
			j = 0;

		//if both y coords of the points are above x-axis can't cross it
		if (vl2[i][1] > 0 && vl2[j][1] > 0)
			counter += 0;

		//if both y coords of the points are below x-axis can't cross it
		else if (vl2[i][1] < 0 && vl2[j][1] < 0)
			counter += 0;

		//if both x coords are negative can't cross +x-axis
		else if (vl2[i][0] < 0 && vl2[j][0] < 0)
			counter += 0;
		//if both x coords are positive a cross occurs
		else if (vl2[i][0] > 0 && vl2[j][0] > 0)
			counter++;

		//if one x is positve and the other isn't check where the intersection is
		else {
			float m = (vl2[j][1] - vl2[i][1]) / (vl2[j][0] - vl2[i][0]); //slope of the line from vi -> vi+1
			float b = vl2[i][1] - (m * vl2[i][0]); //y-intercept of the line from vi -> vi+1

			float xInt = -b / m; //x-intercept

			//if the intercept is negative no +x-axis cross
			if (xInt <= 0)
				counter += 0;
			else
				counter++;
		}
	}
	//if there are an odd number of crosses the point is inside the polygon if even it is not
	if (counter % 2 == 1)
		return true;
	else
		return false;
}

/*
* copyVl: converts 3D coords of a set of vertices to 2D and saves them in a vector to be returned
*/
vector<vector<float>> copyVl(float vl[][3], int nVert, int index) {
	vector<vector<float>> result(nVert, vector<float>(2, 0));

	if (index == 0) {
		for (int i = 0; i < nVert; i++) {
			result[i][0] = vl[i][1];
			result[i][1] = vl[i][2];
		}
		return result;
	}
	else if (index == 1) {
		for (int i = 0; i < nVert; i++) {
			result[i][0] = vl[i][0];
			result[i][1] = vl[i][2];
		}
		return result;
	}
	else {
		for (int i = 0; i < nVert; i++) {
			result[i][0] = vl[i][0];
			result[i][1] = vl[i][1];
		}
		return result;
	}



}

/*
* copyP: converts 3D coords of a point to 2D uses pass by reference so no return
* is necesary.
*/
void copyP(float p[], float(&p2)[2], int index) {

	if (index == 0) {
		p2[0] = p[1];
		p2[1] = p[2];
	}

	else if (index == 1) {
		p2[0] = p[0];
		p2[1] = p[2];
	}

	else {
		p2[0] = p[0];
		p2[1] = p[1];
	}
}

int shading(float p[3], float n[3], float kd) {
	float l[3]; //Vector L
	int c;		//shading value

	//calculate the vector L
	for (int i = 0; i < 3; i++)
		l[i] = p[i] - LPR[i];
	//Get the magnitude of the vector L
	float mag = magnitude(l);

	//Normalize the vector L
	for (int i = 0; i < 3; i++) {
		l[i] = l[i] / mag;
		
	}
	
	//Var to hold the dot product 
	float dot = dotProd(n, l);

	//compute the shading value for the current pixel
	c = Ip * kd * dot;
	
	if (c < 0)
		c = 0;
	//return the shading value
	
	return c;
}

float magnitude(float v[3]) {
	float res; //holds result of magnitude computation

	//compute the magnitude of the Vector
	res = sqrt(pow(v[0], 2) + pow(v[1], 2) + pow(v[2], 2));

	//retrun the computed magnitude
	return res;
}

float dotProd(float v1[3], float v2[3]) {
	float res = 0; //holds the dot product result

	//compute the dot product of the two vectors
	for (int i = 0; i < 3; i++) {
		res += v1[i] * v2[i];
	}

	//return the computed dot product
	return res;
}

vector<float> findnuv(float v1[], float v2[], float p[]) {
	float n[3]; //unit vector n
	float u[3]; //unit vector u
	float v[3]; //unit vector v
	vector<float> result;

	//------------- Calculate n ---------------


	float tM = magnitude(v1);

	for (int i = 0; i < 3; i++)
		n[i] = v1[i] / tM;
	//-------------------------------------------

	//------------- calculate u -----------------


	//cross product between vectors (P3-P1) and (P2-P1)
	float temp[3];
	temp[0] = v2[1] * v1[2] - v2[2] * v1[1];
	temp[1] = v2[2] * v1[0] - v2[0] * v1[2];
	temp[2] = v2[1] * v1[0] - v2[0] * v1[1];

	tM = magnitude(temp);

	for (int i = 0; i < 3; i++)
		u[i] = temp[i] / tM;
	//-----------------------------------------------

	//----------- Calculate v------------------------
	v[0] = n[1] * u[2] - n[2] * u[1];
	v[1] = n[2] * u[0] - n[0] * u[2];
	v[2] = n[0] * u[1] - n[1] * u[0];
	//-----------------------------------------------

	for (int i = 0; i < 3; i++)
		result.push_back(u[i]);

	for (int i = 0; i < 3; i++)
		result.push_back(v[i]);

	for (int i = 0; i < 3; i++)
		result.push_back(n[i]);


	return result;
}

void findmWCmCW(float vpn[], float vup[], float vrp[]) {
	vector<float> nuv = findnuv(vpn, vup, vrp);// holds the u, v, n, values for the given vectors

	//----------------------- World to Camera -------------------------------
	float t[4][4] = { {1, 0, 0, -vrp[0]},
					  {0, 1, 0, -vrp[1]},
					  {0, 0, 1, -vrp[2]},
					  {0, 0, 0, 1}
	};

	float r[4][4] = { {nuv[0], nuv[1], nuv[2], 0},
					  {nuv[3], nuv[4], nuv[5], 0},
					  {nuv[6], nuv[7], nuv[8], 0},
					  {     0,      0,      0, 1}
	};


	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {

				Mwc[i][j] += r[i][k] * t[k][j];
			}
		}
	}
	//----------------------------------------------------------------------------
	//-------------------- Camera to World -----------------------------------
	float invT[4][4] = {
						{1, 0, 0, vrp[0]},
						{0, 1, 0, vrp[1]},
						{0, 0, 1, vrp[2]},
						{0, 0, 0, 1}
	};

	float invR[4][4] = {
					  {nuv[0], nuv[3], nuv[6], 0},
					  {nuv[1], nuv[4], nuv[7], 0},
					  {nuv[2], nuv[5], nuv[8], 0},
					  {     0,      0,      0, 1}
	};



	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {

				Mcw[i][j] += invT[i][k] * invR[k][j];
			}
		}
	}

}
