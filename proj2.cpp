// Name: 			Ryan Watkins
// Quarter, Year:   Fall 2013
// Project:			2
//
// This file is to be modified by the student.
// main.cpp
////////////////////////////////////////////////////////////
#include <GL/glut.h>
#include <cmath>
#include <math.h>
#include <vector>
#include <algorithm>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <string>

#define _USE_MATH_DEFINES

using namespace std;

const int WINDOW_WIDTH = 800;
const int WINDOW_HEIGHT = 800;

const int VSIZE_T = 3;

#define PI 3.14159265
#define RAD PI/180

struct Point3D {
		double x;
		double y;
		double z;
		
		bool even;
		
		vector<int> adj;

		Point3D() : x(0), y(0), z(0), even(false)
		{}

		Point3D(double x, double y, double z, bool e = false)
			: x(x), y(y), z(z), even(e)
		{}

		Point3D(const Point3D & p)
			: x(p.x), y(p.y), z(p.z), adj(p.adj), even(p.even)
		{}
		
		Point3D operator+(const Point3D & rhs) const
		{ return Point3D(x + rhs.x, y + rhs.y, z + rhs.z, even); }
		Point3D operator-(const Point3D & rhs) const
		{ return Point3D(x - rhs.x, y - rhs.y, z - rhs.z, even); }
		Point3D operator*(double rhs) const
		{ return Point3D(x * rhs, y * rhs, z * rhs, even); }
		Point3D operator/(double rhs) const
		{ return Point3D(x / rhs, y / rhs, z / rhs, even); }
		Point3D operator+=(const Point3D & rhs)
		{ x += rhs.x; y += rhs.y; z += rhs.z; return *this; }
		Point3D operator-=(const Point3D & rhs)
		{ x -= rhs.x; y -= rhs.y; z -= rhs.z; return *this; }
		Point3D operator*=(double rhs)
		{ x *= rhs; y *= rhs; z *= rhs; return *this; }
		Point3D operator/=(double rhs)
		{ x /= rhs; y /= rhs; z /= rhs; return *this; }
		bool operator==(Point3D rhs) { 
			if (x == rhs.x && y == rhs.y && z == rhs.z) 
				return true; 
			return false;
		}

		double magnitude() const
		{ return sqrt(x * x + y * y + z * z); }
		void normalize()
		{ *this /= magnitude(); }
		double dot(const Point3D & rhs) const
		{
			return x * rhs.x + y * rhs.y + z * rhs.z;
		}
		Point3D cross(const Point3D & rhs) const
		{
			return Point3D(y * rhs.z - z * rhs.y,
					z * rhs.x - x * rhs.z,
					x * rhs.y - y * rhs.x, even);
		}
		
		void addAdj(int p) {
			bool exists = false;
			for (int i = 0; i < adj.size(); i++) {
				if (adj[i] == p)
					exists = true;
			}
			if (exists == false) {
				adj.push_back(p);
			}
		}
};

struct TriPoly {
		int a;
		int b;
		int c;

		TriPoly() : a(0), b(0), c(0)
		{}

		TriPoly(int a, int b, int c)
			: a(a), b(b), c(c)
		{}

		TriPoly(const TriPoly & t)
			: a(t.a), b(t.b), c(t.c)
		{}
};

struct HPair {
	int p;
	int f;
	
	HPair() : p(0), f(0)
	{}
	
	HPair(int p, int f)
		: p(p), f(f)
	{}
};

struct Vector3
{
	double x;
	double y;
	double z;

	Vector3() : x(0.0), y(0.0), z(0.0)
	{}

	Vector3(double x, double y, double z)
		: x(x), y(y), z(z)
	{}

	Vector3(const Vector3 & v)
		: x(v.x), y(v.y), z(v.z)
	{}

	Vector3 operator+(const Vector3 & rhs) const
	{ return Vector3(x + rhs.x, y + rhs.y, z + rhs.z); }
	Vector3 operator-(const Vector3 & rhs) const
	{ return Vector3(x - rhs.x, y - rhs.y, z - rhs.z); }
	Vector3 operator*(double rhs) const
	{ return Vector3(x * rhs, y * rhs, z * rhs); }
	Vector3 operator/(double rhs) const
	{ return Vector3(x / rhs, y / rhs, z / rhs); }
	Vector3 operator+=(const Vector3 & rhs)
	{ x += rhs.x; y += rhs.y; z += rhs.z; return *this; }
	Vector3 operator-=(const Vector3 & rhs)
	{ x -= rhs.x; y -= rhs.y; z -= rhs.z; return *this; }
	Vector3 operator*=(double rhs)
	{ x *= rhs; y *= rhs; z *= rhs; return *this; }
	Vector3 operator/=(double rhs)
	{ x /= rhs; y /= rhs; z /= rhs; return *this; }

	double magnitude() const
	{ return sqrt(x * x + y * y + z * z); }
	void normalize()
	{ *this /= magnitude(); }
	double dot(const Vector3 & rhs) const
	{
		return x * rhs.x + y * rhs.y + z * rhs.z;
	}
	Vector3 cross(const Vector3 & rhs) const
	{
		return Vector3(y * rhs.z - z * rhs.y,
					z * rhs.x - x * rhs.z,
					x * rhs.y - y * rhs.x);
	}
};

struct Ray {
	Vector3 position;
   	Vector3 direction;
   	
   	Ray() : position(Vector3(0,0,0)), direction(Vector3(0,0,0))
	{}

	Ray(Vector3 p, Vector3 d)
		: position(p), direction(d)
	{}

	Ray(const Ray & r)
		: position(r.position), direction(r.direction)
	{}
};

/////////////////////
// Global Storage
////////////////////
vector<Point3D> points;
vector<TriPoly> faces;
vector<Point3D> pointsD1;
vector<TriPoly> facesD1;
vector<Point3D> pointsOrig;
vector<TriPoly> facesOrig;
vector<Point3D> planePoints;
vector<TriPoly> planeFaces;
HPair header;
int view = 1;
int LSDCout = 0;
bool moveMode = false;
bool rotMode = false;
bool shaderMode = false;
double tREF = INFINITY;
Point3D sceneCamera(50, 20, -100);
Point3D sceneViewer(400, 400, -100);
int lightHorz = 800;
int lightVert = 800;
int lightDepth = -500;

////////////////Used for division
//find the point among existing to prevent redundency
int findRef1 (Point3D f) {
	int ref = -1;
	for (int i = 0; i < pointsD1.size(); i++) {
		if (pointsD1[i] == f)
			ref = i;
	}
	return ref;
}

// Renders a quad at cell (x, y) with dimensions CELL_LENGTH
void renderPixel(double x, double y, float r = 1.0, float g = 1.0, float b = 1.0)
{
	glBegin(GL_POINTS);
		glColor3f(r, g, b);
		glVertex2f(x, y);
	glEnd();
}

void DDA(double xStart, double yStart, double xEnd, double yEnd) { 
	int steps;
	float dx, dy, xInc, yInc, slope, x, y;

	dx = xEnd - xStart; 
	dy = yEnd - yStart;

	slope = dy/dx; 

	if (abs(dx) > abs(dy)) { 
		steps = (int)abs(dx); 
	} else { 
		steps = (int)abs(dy); 
	} 

	xInc = dx/steps; 
	yInc = dy/steps; 

	x = xStart; 
	y = yStart; 

	for (unsigned int i = 0; i < steps; i++) { 	
		x += xInc; 
		y += yInc; 
		renderPixel((double)x, (double)y);
	} 
}

//////////////////////
// File Processing
////////////////////
void parseFile(string filepath) {
	// open file
	ifstream fin;
	fin.open(filepath.c_str());
	if (!fin.good()) {
		cout << "FILE READ ERROR" << endl;
		return;
	}

	// locally store all data from file
	int i = 0;
	vector<int> data;
	string line;
	while (fin.good()) {
		getline(fin, line);
		istringstream os(line);

		while(os >> i) {
    		data.push_back(i);
    	}
    }
	
	//grab the header
	int pointCount = data[0];
	int faceCount = data[1];
	header.p = pointCount;
	header.f = faceCount;
	data.erase(data.begin());
	data.erase(data.begin());

	// parse out the points into a vector 
    Point3D p;
    for (int i = 0; i < pointCount*VSIZE_T; i+=VSIZE_T) {
 		p.x = data[i];
    	p.y = data[i+1];
    	p.z = data[i+2];
    	points.push_back(p);
    }
    
    // deque the points section from the local storage to access the faces
    for (int i = 0; i < pointCount*VSIZE_T; i++) {
    	data.erase(data.begin());
    }
    
    // parse out the faces into another vector 		
    TriPoly t;
    for (int i = 0; i < faceCount*VSIZE_T; i+=VSIZE_T) {
    	t.a = data[i];
    	t.b = data[i+1];
    	t.c = data[i+2];
    	
    	points[t.a].addAdj(t.b);
    	points[t.a].addAdj(t.c);
    	
    	points[t.b].addAdj(t.a);
    	points[t.b].addAdj(t.c);
    	
    	points[t.c].addAdj(t.a);
    	points[t.c].addAdj(t.b);
    	
    	faces.push_back(t);
    }
    facesOrig = faces;
    pointsOrig = points;
}

bool rayTraceTri(Ray r, TriPoly face) {
	Point3D pa = points[face.a];
	Point3D pb = points[face.b];
	Point3D pc = points[face.c];
	
	// Triangle edge Vectors
	Point3D ab = pb - pa;
	Point3D ac = pc - pa;
	
	// variables to build the matrices:
	// | a d g || beta  |   | j | 
	// | b e h || gamma | = | k | 
	// | c f i ||   t   |   | l |
	
	// In use with Crammer's Rule
	double a = -(ab.x);
	double b = -(ab.y);
	double c = -(ab.z);
	
	double d = -(ac.x);
	double e = -(ac.y);
	double f = -(ac.z);
	
	double g = r.direction.x;
	double h = r.direction.y;
	double i = r.direction.z;

	double j = pa.x - r.position.x;
	double k = pa.y - r.position.y;
	double l = pa.z - r.position.z;
	
	// Selecting certain reference rows we can reuse some calculations
	// so these condense those calculations
	double eihf = (e*i) - (h*f);
	double difg = (g*f) - (d*i);
	double dheg = (d*h) - (e*g);
	
	double akjb = (a*k) - (j*b);
	double alcj = (c*j) - (a*l);
	double blkc = (b*l) - (k*c);
	
	// A is the determinant of the matrix that divides each matrix
	double A = (a*((e*i) - (f*h))) - (b*((d*i) - (f*g))) + (c*((d*h) - (e*g)));
	// determinant that will be t
	double t = ((a*((e*l) - (f*k))) - (b*((d*l) - (f*j))) + (c*((d*k) - (e*j)))) / A;
	// determinant that will be gamma
	double gamma = ((a*((k*i) - (l*h))) - (b*((j*i) - (l*g))) + (c*((j*h) - (k*g)))) / A;
	// determinant that will be beta
	double beta = ((j*((e*i) - (f*h))) - (k*((d*i) - (f*g))) + (l*((d*h) - (e*g)))) / A;

	// Constraints of hitting an intersection
	if (t < 0 || t > INFINITY)
		return false;
	if (gamma < 0 || gamma > 1)
		return false;
	if (beta < 0 || beta > 1-gamma)
		return false;
	if (t > tREF)
		return false;
	
	tREF = t;
	return true;
}


////////////////Used for division
//places new verticies
Point3D createOdd1(Point3D a, Point3D b) {
	//mimics a winged edge
	int wing1 = -1;
	int wing2 = -1;
	Point3D odd;
	for (int i = 0; i < a.adj.size(); i++) {
		for (int j = 0; j < b.adj.size(); j++) {
			if (a.adj[i] == b.adj[j])
				wing1 = a.adj[i];
		}
	}
	for (int i = 0; i < a.adj.size(); i++) {
		for (int j = 0; j < b.adj.size(); j++) {
			if (a.adj[i] == b.adj[j] && a.adj[i] != wing1)
				wing2 = a.adj[i];
		}
	}
	if (wing2 < 0) 
		odd = (a + b) / 2;
	else {
		odd = a * 3/8 + b * 3/8 + points[wing1] * 1/8 + points[wing2] * 1/8;
	}
	return odd;
}

////////////////Used for division 
//repositions old verticies
void refineEven1() {
	//Update adj list by removing even to even connections
	for (int i = 0; i < pointsD1.size(); i++) {
		vector<int> newAdj;
		for (int j = 0; j < pointsD1[i].adj.size(); j++) {
			if (pointsD1[i].even == true && pointsD1[pointsD1[i].adj[j]].even == false) {
				newAdj.push_back(pointsD1[i].adj[j]);
			}

		}
		if (pointsD1[i].even == true) {
		//redundent checking for redundencies
			bool exists = false;
			vector<int> newnewAdj;
			for (int j = 0; j < newAdj.size(); j++) {
				for (int k = 0; k < newnewAdj.size(); k++) {
					if (pointsD1[newAdj[j]].x == pointsD1[newnewAdj[k]].x && pointsD1[newAdj[j]].y == pointsD1[newnewAdj[k]].y && pointsD1[newAdj[j]].z == pointsD1[newnewAdj[k]].z)
						exists = true;
				}
				if (!exists)
					newnewAdj.push_back(newAdj[j]);
				exists = false;
			}
			
			pointsD1[i].adj = newnewAdj;
			//redefine position of even verticies
			double n = pointsD1[i].adj.size();

			double b = 0;
			if (n > 3)
				b = 3/(8*n);
			else {
				b = 4;
				b /=16;
			}

			Point3D pSum;
			for (int j = 0; j < pointsD1[i].adj.size(); j++) {
				pSum += pointsD1[pointsD1[i].adj[j]] * b;
			}
			pointsD1[i] *= 1 - n*b;
			pointsD1[i] += pSum;
		}
		newAdj.clear();
	}
	//Update all verticies to odd
	for (int i = 0; i < pointsD1.size(); i++) {
		pointsD1[i].even = false;
	}
}

////////////////Used for division
//creates new verticies, calls createOdd to place them and calls 
//refine even to repos. old verts
void subdiv1() {
	for (int i = 0; i < faces.size(); i ++) {
		//Pull face reference
		Point3D a = points[faces[i].a];
		Point3D b = points[faces[i].b];
		Point3D c = points[faces[i].c];
		//Create odd points and indecies
		Point3D p1 = createOdd1(a,b);
		int	p1Ref = findRef1(p1);
			if (p1Ref < 0) {
				pointsD1.push_back(p1);
				p1Ref = pointsD1.size() - 1;
			}
		Point3D p2 = createOdd1(b,c);
		int	p2Ref = findRef1(p2);
			if (p2Ref < 0) {
				pointsD1.push_back(p2);
				p2Ref = pointsD1.size() - 1;
			}
		Point3D p3 = createOdd1(a,c);
		int	p3Ref = findRef1(p3);
			if (p3Ref < 0) {
				pointsD1.push_back(p3);
				p3Ref = pointsD1.size() - 1;
			}
		//Set new points adj lists
		pointsD1[p1Ref].addAdj(faces[i].a);
		pointsD1[p1Ref].addAdj(faces[i].b);
		pointsD1[p1Ref].addAdj(p2Ref);
		pointsD1[p1Ref].addAdj(p3Ref);
		
		pointsD1[p2Ref].addAdj(faces[i].b);
		pointsD1[p2Ref].addAdj(faces[i].c);
		pointsD1[p2Ref].addAdj(p1Ref);
		pointsD1[p2Ref].addAdj(p3Ref);
		
		pointsD1[p3Ref].addAdj(faces[i].a);
		pointsD1[p3Ref].addAdj(faces[i].c);
		pointsD1[p3Ref].addAdj(p1Ref);
		pointsD1[p3Ref].addAdj(p2Ref);
		
		//Update old adj lists
		pointsD1[faces[i].a].addAdj(p1Ref);
		pointsD1[faces[i].a].addAdj(p3Ref);
		
		pointsD1[faces[i].b].addAdj(p1Ref);
		pointsD1[faces[i].b].addAdj(p2Ref);
		
		pointsD1[faces[i].c].addAdj(p2Ref);
		pointsD1[faces[i].c].addAdj(p3Ref);
		
		//Create sub faces
		TriPoly oddFace(p1Ref, p2Ref, p3Ref);
		TriPoly evenFace1(faces[i].a, p1Ref, p3Ref);
		TriPoly evenFace2(faces[i].b, p1Ref, p2Ref);
		TriPoly evenFace3(faces[i].c, p2Ref, p3Ref);
		
		//Add new faces to face table
		facesD1.push_back(oddFace);
		facesD1.push_back(evenFace1);
		facesD1.push_back(evenFace2);
		facesD1.push_back(evenFace3);
		
		//Set even face ids
		pointsD1[faces[i].a].even = true;
		pointsD1[faces[i].b].even = true;
		pointsD1[faces[i].c].even = true;
	}
	//smoothe the mesh
	refineEven1();
}

////////////////Used for division
//driver for the subdivision algorithm for quick calls
void loopSubDivide1() {
	// Restore to default division after subdivide iteration 3
	if (LSDCout > 2) {
		points = pointsOrig;
		faces = facesOrig;
		pointsD1.clear();
		facesD1.clear();
		LSDCout = 0;
		return;
	}
	
	pointsD1 = points;
	subdiv1();
	points = pointsD1;
	faces = facesD1;
	LSDCout++;
}

// Translation handlers
// translate in given direction
void translate (unsigned char direction) {
	int dist = 10; 
	switch (direction) {
		case 'u':
			for (int i = 0; i < points.size(); i++) {
				points[i].y += dist;
				if (view == 3)
					points[i].z += dist;
			}
			pointsOrig = points;
			break;
		case 'd':
			dist *= -1;
			for (int i = 0; i < points.size(); i++) {
				points[i].y += dist;
				if (view == 3)
					points[i].z += dist;
			}
			pointsOrig = points;
			break;
		case 'l':
			dist *= -1;
			for (int i = 0; i < points.size(); i++) {
				points[i].x += dist;
				if (view == 2)
					points[i].z += dist;
			}
			pointsOrig = points;
			break;
		case 'r':
			for (int i = 0; i < points.size(); i++) {
				points[i].x += dist;
				if (view == 2)
					points[i].z += dist;
			}
			pointsOrig = points;
			break;
	}

}
// choose the direction with the ARROW keys
void moveObj(int key, int x, int y) {
	if (!moveMode)
		return;
	
	switch (key) {
		case GLUT_KEY_LEFT:
			translate ('l');
			break;
		case GLUT_KEY_UP:
			translate ('u');
			break;
		case GLUT_KEY_DOWN:
			translate ('d');
			break;
		case GLUT_KEY_RIGHT:
			translate ('r');
			break;
	}
	
	glutPostRedisplay();
}

// Rotation handlers
// helper to get the origin translation vector
double getCenterInc (char coord) {
	double min = INFINITY;
	double max = 0;
	switch (coord) {
		case 'x':
			for (int i = 0; i < points.size(); i++) {
				if (points[i].x < min)
					min = points[i].x;
				if (points[i].x > max)
					max = points[i].x;	
			}
			break;
		case 'y':
			for (int i = 0; i < points.size(); i++) {
				if (points[i].y < min)
					min = points[i].y;
				if (points[i].y > max)
					max = points[i].y;	
			}
			break;
		case 'z':
			for (int i = 0; i < points.size(); i++) {
				if (points[i].z < min)
					min = points[i].z;
				if (points[i].z > max)
					max = points[i].z;	
			}
			break;
	}
	
	return (min + max) / 2;
}
// Rotation handlers
// rotate in given direction
void rotate (unsigned char direction) {
	int deg = 10; 
	// Get the origin average incrementers to translate to origin
	double midx = getCenterInc('x');
	double midy = getCenterInc('y');
	double midz = getCenterInc('z');
	switch (direction) {
		case 'u':
			// Translate to origin
			for (int i = 0; i < points.size(); i++) {
				points[i].x -= midx;
				points[i].y -= midy;
				points[i].z -= midz;
			}
			// Rotate about axis
			for (int i = 0; i < points.size(); i++) {
				double newY = points[i].y * cos(deg*RAD) - points[i].z * sin(deg*RAD);
				double newZ = points[i].y * sin(deg*RAD) + points[i].z * cos(deg*RAD);
				points[i].y = newY;
				points[i].z = newZ;
			}
			// Translate back
			for (int i = 0; i < points.size(); i++) {
				points[i].x += midx;
				points[i].y += midy;
				points[i].z += midz;
			}
			//pointsOrig = points;
			break;
		case 'd':
			deg *= -1;
			// Translate to origin
			for (int i = 0; i < points.size(); i++) {
				points[i].x -= midx;
				points[i].y -= midy;
				points[i].z -= midz;
			}
			// Rotate about axis
			for (int i = 0; i < points.size(); i++) {
				double newY = points[i].y * cos(deg*RAD) - points[i].z * sin(deg*RAD);
				double newZ = points[i].y * sin(deg*RAD) + points[i].z * cos(deg*RAD);
				points[i].y = newY;
				points[i].z = newZ;
			}
			// Translate back
			for (int i = 0; i < points.size(); i++) {
				points[i].x += midx;
				points[i].y += midy;
				points[i].z += midz;
			}
			//pointsOrig = points;
			break;
		case 'l':
			deg *= -1;
			// Translate to origin
			for (int i = 0; i < points.size(); i++) {
				points[i].x -= midx;
				points[i].y -= midy;
				points[i].z -= midz;
			}
			// Rotate about axis
			for (int i = 0; i < points.size(); i++) {
				double newX = points[i].x * cos(deg*RAD) + points[i].z * sin(deg*RAD);
				double newZ = points[i].x * -(sin(deg*RAD)) + points[i].z * cos(deg*RAD);
				points[i].x = newX;
				points[i].z = newZ;
			}
			// Translate back
			for (int i = 0; i < points.size(); i++) {
				points[i].x += midx;
				points[i].y += midy;
				points[i].z += midz;
			}
			//pointsOrig = points;
			break;
		case 'r':
			// Translate to origin
			for (int i = 0; i < points.size(); i++) {
				points[i].x -= midx;
				points[i].y -= midy;
				points[i].z -= midz;
			}
			// Rotate about axis
			for (int i = 0; i < points.size(); i++) {
				double newX = points[i].x * cos(deg*RAD) + points[i].z * sin(deg*RAD);
				double newZ = points[i].x * -(sin(deg*RAD)) + points[i].z * cos(deg*RAD);
				points[i].x = newX;
				points[i].z = newZ;
			}
			// Translate back
			for (int i = 0; i < points.size(); i++) {
				points[i].x += midx;
				points[i].y += midy;
				points[i].z += midz;
			}
			//pointsOrig = points;
			break;
	}

}
// choose the direction with the ARROW keys
void rotObj(int key, int x, int y) {
	if (!rotMode)
		return;
	
	switch (key) {
		case GLUT_KEY_LEFT:
			rotate ('l');
			break;
		case GLUT_KEY_UP:
			rotate ('u');
			break;
		case GLUT_KEY_DOWN:
			rotate ('d');
			break;
		case GLUT_KEY_RIGHT:
			rotate ('r');
			break;
	}
	
	glutPostRedisplay();
}

Point3D phongShade (float r, float g, float b, Point3D p, int f) {
	// Ambient is set as base color
	Point3D ambient(r, g, b);
	// Specular reflection
	double reflect = 50;
	// Light source
	Point3D light(lightHorz,lightVert, lightDepth);
	renderPixel(lightHorz, lightVert, 1.0, 0.0, 0.0);
	// Direction to viewer
	Point3D viewer(200, 200, 0);
	// Get face normal
	Point3D pa = points[faces[f].a];
	Point3D pb = points[faces[f].b];
	Point3D pc = points[faces[f].c];
	Point3D v = pa - pc;
	Point3D w = pa - pb;
	Point3D fNorm = v.cross(w);
	
	double centerX = getCenterInc('x');
	double centerY = getCenterInc('y');
	double centerZ = getCenterInc('z');
	
	Point3D centerP(centerX, centerY, centerZ); 
	
	if (fNorm.dot(centerP) > 0) {
		fNorm = fNorm * -1;
	}

	fNorm.normalize();
	
	// Get light direction
	light -= p;
	light.normalize();

	// Calculate diffuse factor
	double difScale = fNorm.dot(light);
	// Diffuse
	Point3D diffuse = ambient;
	diffuse = diffuse * max(0.0, difScale);
	diffuse = diffuse * 2;

	// Calculate reflected ray
	Point3D rbounce = fNorm;
	rbounce = rbounce * difScale;
	rbounce = rbounce * 2;
	rbounce -= light; 
	viewer -= p;
	rbounce.normalize();
	viewer.normalize();
	
	// Calculate specular factor
	double specScale = viewer.dot(rbounce);
	Point3D specular(.5,.5,.5);
	
	// Specular
	specular = specular * pow(max(0.0, specScale), reflect);

	Point3D phong;
	phong += diffuse;
	phong += specular;
	return phong;
}


void iRender() {
	// For perspective projection
	Point3D camera(4*50, 4*20, -600);
	Point3D viewer(400, 400, -400);
	cout << "Raytracing..." << endl;
	for (int i = 0; i < WINDOW_WIDTH; i++) {
		for (int j = 0; j < WINDOW_HEIGHT; j++) {
			tREF = INFINITY;
			for (int k = 0; k < faces.size(); k++) {
				//ray construction
				Vector3 rp(i,j,0);
				Vector3 rd(0,0,1);
				Ray rz(rp,rd); 

				if (rayTraceTri(rz, faces[k])) {
					Point3D toRend(i, j, tREF);
					Point3D color = phongShade (0.0, 0.0, 1.0, toRend, k);
					renderPixel(i, j, color.x, color.y, color.z);
				}
			}
		}
	}
	cout << "Raytrace Complete" << endl;
}

bool rayTraceTriPlane(Ray r, TriPoly face) {
	Point3D pa = planePoints[face.a];
	Point3D pb = planePoints[face.b];
	Point3D pc = planePoints[face.c];

	// Triangle edge Vectors
	Point3D ab = pb - pa;
	Point3D ac = pc - pa;
	
	// In use with Crammer's Rule
	double a = -(ab.x);
	double b = -(ab.y);
	double c = -(ab.z);
	
	double d = -(ac.x);
	double e = -(ac.y);
	double f = -(ac.z);
	
	double g = r.direction.x;
	double h = r.direction.y;
	double i = r.direction.z;

	double j = pa.x - r.position.x;
	double k = pa.y - r.position.y;
	double l = pa.z - r.position.z;
	
	// Selecting certain reference rows we can reuse some calculations
	// so these condense those calculations
	double eihf = (e*i) - (h*f);
	double difg = (g*f) - (d*i);
	double dheg = (d*h) - (e*g);
	
	double akjb = (a*k) - (j*b);
	double alcj = (c*j) - (a*l);
	double blkc = (b*l) - (k*c);
	
	// A is the determinant of the matrix that divides each matrix
	double A = (a*((e*i) - (f*h))) - (b*((d*i) - (f*g))) + (c*((d*h) - (e*g)));
	// determinant that will be t
	double t = ((a*((e*l) - (f*k))) - (b*((d*l) - (f*j))) + (c*((d*k) - (e*j)))) / A;
	// determinant that will be gamma
	double gamma = ((a*((k*i) - (l*h))) - (b*((j*i) - (l*g))) + (c*((j*h) - (k*g)))) / A;
	// determinant that will be beta
	double beta = ((j*((e*i) - (f*h))) - (k*((d*i) - (f*g))) + (l*((d*h) - (e*g)))) / A;

	// Constraints of hitting an intersection
	if (t < 0 || t > INFINITY)
		return false;
	if (gamma < 0 || gamma > 1)
		return false;
	if (beta < 0 || beta > 1-gamma)
		return false;
	if (t > tREF)
		return false;
	
	tREF = t;
	return true;
}

void rendPlane() {
	for (int i = 0; i < WINDOW_WIDTH; i++) {
		for (int j = 0; j < WINDOW_HEIGHT; j++) {
			tREF = INFINITY;
			for (int k = 0; k < planeFaces.size(); k++) {
				//ray construction
				Vector3 rp(i,j,0);
				Vector3 rd(0,0,1);
				Ray rz(rp,rd); 

				if (rayTraceTriPlane(rz, planeFaces[k])) {
					Point3D color(0.5, 0.5, 0.5);
					rp.x = i;
					rp.y = j;
					rp.z = tREF;	
					rd.x = lightHorz;
					rd.y = lightVert;
					rd.z = lightDepth;
					rd -= rp;
					rz.position = rp;
					rz.direction = rd;
					if (rayTraceTri(rz, faces[k])) {
						color = color / 2;
					}
					renderPixel(i, j, color.x, color.y, color.z);
				}
			}
		}
	}
}

void wRender(char viewPort) {
	if (shaderMode)
		return;
		
	switch (viewPort) {
		case 'f':
			for (int i = 0; i < faces.size(); i++) {
				int a = faces[i].a;
				int b = faces[i].b;
				int c = faces[i].c;
				Point3D pa = points[a];
				Point3D pb = points[b];
				Point3D pc = points[c];

				DDA(pa.x, pa.y, pb.x, pb.y);
				DDA(pb.x, pb.y, pc.x, pc.y);
				DDA(pc.x, pc.y, pa.x, pa.y);
			}
			break;
		case 's':
			for (int i = 0; i < faces.size(); i++) {
				int a = faces[i].a;
				int b = faces[i].b;
				int c = faces[i].c;
				Point3D pa = points[a];
				Point3D pb = points[b];
				Point3D pc = points[c];

				DDA(pa.z, pa.y, pb.z, pb.y);
				DDA(pb.z, pb.y, pc.z, pc.y);
				DDA(pc.z, pc.y, pa.z, pa.y);
			}
			break;
		case 't':
			for (int i = 0; i < faces.size(); i++) {
				int a = faces[i].a;
				int b = faces[i].b;
				int c = faces[i].c;
				Point3D pa = points[a];
				Point3D pb = points[b];
				Point3D pc = points[c];

				DDA(pa.x, pa.z, pb.x, pb.z);
				DDA(pb.x, pb.z, pc.x, pc.z);
				DDA(pc.x, pc.z, pa.x, pa.z);
		
			}
			break;
	}
}

//////////////////
// FRONT VIEW
void renderFront() {
	wRender('f');
}

//////////////////
// SIDE VIEW
void renderSide() {
	wRender('s');
}

//////////////////
// TOP VIEW
void renderTop() {
	wRender('t');
}

void scaleObj(double amount) {
	for (int i = 0; i < points.size(); i++) {
		
		double midx = getCenterInc('x');
		double midy = getCenterInc('y');
		double midz = getCenterInc('z');

		// Translate to origin
//		points[i].x -= midx;
//		points[i].y -= midy;
//		points[i].z -= midz;
		// Scale
		points[i] *= amount;
		// Translate back
//		points[i].x += midx;
//		points[i].y += midy;
//		points[i].z += midz;
	}
	pointsOrig = points;
}

// Perspective View
void renderPerspective() {
	if (shaderMode)
		return;

	Point3D camera(4*50, 4*20, -600);
	Point3D viewer(400, 400, -400);
	for (int i = 0; i < faces.size(); i++) {
		int a = faces[i].a;
		int b = faces[i].b;
		int c = faces[i].c;
		Point3D pa = points[a];
		Point3D pb = points[b];
		Point3D pc = points[c];	

		// for point a
		Point3D da = pa - camera;
		Point3D db = pb - camera;
		Point3D dc = pc - camera;
		
		double baX = viewer.z / da.z * da.x - viewer.x;
		double baY = viewer.z / da.z * da.y - viewer.y;
		double bbX = viewer.z / db.z * db.x - viewer.x;
		double bbY = viewer.z / db.z * db.y - viewer.y;
		double bcX = viewer.z / dc.z * dc.x - viewer.x;
		double bcY = viewer.z / dc.z * dc.y - viewer.y;
		
		baX *= -1;
		baY *= -1;
		bbX *= -1;
		bbY *= -1;
		bcX *= -1;
		bcY *= -1;
		
		DDA(baX, baY, bbX, bbY);
		DDA(bbX, bbY, bcX, bcY);
		DDA(bcX, bcY, baX, baY);
	}
}


// VIEW SELECTION CONTROLS
// TRANSFORM CONTROLS
void viewSelect(unsigned char key, int x, int y) {
	switch (key) {
		case '1':
			view = 1;
			break;
		case '2':
			view = 2;
			break;
		case '3':
			view = 3;
			break;
		case 'l':
			loopSubDivide1();
			break;
		case 't':
			if (rotMode) {
				cout << "Rotate Mode [OFF]" << endl;
				rotMode = false;
			}
			if (!moveMode) {
				cout << "Translate Mode [ON]" << endl;
				moveMode = true;
				glutSpecialFunc(moveObj);
			} else {
				cout << "Translate Mode [OFF]" << endl;
				moveMode = false;
			}
			break;
		case 'r':
			if (moveMode) {
				cout << "Translate Mode [OFF]" << endl;
				moveMode = false;
			}
			if (!rotMode) {
				cout << "Rotate Mode [ON]" << endl;
				rotMode = true;
				glutSpecialFunc(rotObj);
			} else {
				cout << "Rotate Mode [OFF]" << endl;
				rotMode = false;
			}
			break;
		case 'i':
			if (!shaderMode) {
				cout << "Raytrace Mode [ON]" << endl;
				shaderMode = true;
			} else {
				cout << "Raytrace Mode [OFF]" << endl;
				shaderMode = false;
			}
			break;
	}
	glutPostRedisplay();
}



//Output function to OpenGL Buffer
void GL_render()
{
	glClear(GL_COLOR_BUFFER_BIT);
	
	switch (view) {
		case 1:
			renderFront();
			if (shaderMode) {
				rendPlane();
				iRender();
			}
			break;
		case 2:
			renderSide();
			if (shaderMode) {
				rendPlane();
				iRender();
			}
			break;
		case 3:
			renderTop();
			if (shaderMode) {
				rendPlane();
				iRender();
			}
			break;
		default:
			renderFront();
			break;
	}

	glutSwapBuffers();
}

//Initializes OpenGL attributes
void GLInit(int* argc, char** argv)
{
	glutInit(argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);

	// ...
	// Complete this function
	// ...
	glutCreateWindow("CS 130 - Ryan Watkins");

	// The default view coordinates is (-1.0, -1.0) bottom left & (1.0, 1.0) top right.
	// For the purposes of this lab, this is set to the number of pixels
	// in each dimension.
	glMatrixMode(GL_PROJECTION_MATRIX);
	glOrtho(0, WINDOW_WIDTH, 0, WINDOW_HEIGHT, -1, 1);

	glutDisplayFunc(GL_render);
	
}

int main(int argc, char** argv)
{	
	if (argc < 2) {
		cout << "Usage: render renderfile.txt" << endl;
		return -1;
	}

	parseFile(argv[1]);

	cout << "finished parse" << endl;
	for (int i = 0; i < points.size(); i++) {
		cout << i << ": " << points[i].x << " , " << points[i].y << " , " << points[i].z << endl;
	}
	cout << "---------------------" << endl;
	for (int i = 0; i < faces.size(); i++) {
		cout << faces[i].a << " , " << faces[i].b << " , " << faces[i].c << endl;
	}
	cout << "END" << endl;
	scaleObj(6);
	
	Point3D p1(50, 50, 500);
	Point3D p2;
	Point3D p3(WINDOW_WIDTH, 0, 0);
	Point3D p4(WINDOW_WIDTH - 50, 50, 500);
	p1 *= 3;
	p2 *= 3;
	p3 *= 3;
	p4 *= 3;
	planePoints.push_back(p1);
	planePoints.push_back(p2);
	planePoints.push_back(p3);
	planePoints.push_back(p4);
	TriPoly f1(0, 2, 1);
	TriPoly f2(3, 2, 0);
	planeFaces.push_back(f1);
	planeFaces.push_back(f2);
	
	GLInit(&argc, argv);
//	glutSpecialFunc(viewSelect);
	glutKeyboardFunc(viewSelect);
	glutMainLoop();

	return 0;
}


