/* Release code for program 1 CS 180 */

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <assert.h>

#include "tiny_obj_loader.h"
#include "Image.h"

// This allows you to skip the `std::` in front of C++ standard library
// functions. You can also say `using std::cout` to be more selective.
// You should never do this in a header file.
using namespace std;

struct Vertex {
   float x;
   float y;
   float z;
};

struct Color {
   float r;
   float g;
   float b;
};

struct Triangle {
   Vertex v1;
   Vertex v2;
   Vertex v3;
   Color c1;
   Color c2;
   Color c3;
};

struct BoundingBox {
   float xmin;
   float xmax;
   float ymin;
   float ymax;
};

//global variables for possible use with transforms, etc.
int g_width, g_height;

// Recommended values
float r_right;
float r_left;
float r_bottom;
float r_top;

// Bounding Box
struct BoundingBox bb;

float epsilon = 0.00001;

// Chosen colors
float global_r = 255;
float global_g = 0;
float global_b = 0;

// Store color mode
int colorMode;

// Track groupings of vertices
vector<vector<float>> parsedT;

// z buffer to store z values
vector<vector<float>> zBuf;

/*
   Helper function you will want all quarter
   Given a vector of shapes which has already been read from an obj file
   resize all vertices to the range [-1, 1]
   THIS IS NOT A WINDOW COORDINATE TRANSFORM
 */
void resize_obj(std::vector<tinyobj::shape_t> &shapes){
   float minX, minY, minZ;
   float maxX, maxY, maxZ;
   float scaleX, scaleY, scaleZ;
   float shiftX, shiftY, shiftZ;
   float epsilon = 0.001;

   minX = minY = minZ = 1.1754E+38F;
   maxX = maxY = maxZ = -1.1754E+38F;

   //Go through all vertices to determine min and max of each dimension
   for (size_t i = 0; i < shapes.size(); i++) {
      for (size_t v = 0; v < shapes[i].mesh.positions.size() / 3; v++) {
         if(shapes[i].mesh.positions[3*v+0] < minX) minX = shapes[i].mesh.positions[3*v+0];
         if(shapes[i].mesh.positions[3*v+0] > maxX) maxX = shapes[i].mesh.positions[3*v+0];

         if(shapes[i].mesh.positions[3*v+1] < minY) minY = shapes[i].mesh.positions[3*v+1];
         if(shapes[i].mesh.positions[3*v+1] > maxY) maxY = shapes[i].mesh.positions[3*v+1];

         if(shapes[i].mesh.positions[3*v+2] < minZ) minZ = shapes[i].mesh.positions[3*v+2];
         if(shapes[i].mesh.positions[3*v+2] > maxZ) maxZ = shapes[i].mesh.positions[3*v+2];
      }
   }

	//From min and max compute necessary scale and shift for each dimension
   float maxExtent, xExtent, yExtent, zExtent;
   xExtent = maxX-minX;
   yExtent = maxY-minY;
   zExtent = maxZ-minZ;
   if (xExtent >= yExtent && xExtent >= zExtent) {
      maxExtent = xExtent;
   }
   if (yExtent >= xExtent && yExtent >= zExtent) {
      maxExtent = yExtent;
   }
   if (zExtent >= xExtent && zExtent >= yExtent) {
      maxExtent = zExtent;
   }
   scaleX = 2.0 /maxExtent;
   shiftX = minX + (xExtent/ 2.0);
   scaleY = 2.0 / maxExtent;
   shiftY = minY + (yExtent / 2.0);
   scaleZ = 2.0/ maxExtent;
   shiftZ = minZ + (zExtent)/2.0;

   //Go through all verticies shift and scale them
   for (size_t i = 0; i < shapes.size(); i++) {
      for (size_t v = 0; v < shapes[i].mesh.positions.size() / 3; v++) {
         shapes[i].mesh.positions[3*v+0] = (shapes[i].mesh.positions[3*v+0] - shiftX) * scaleX;
         assert(shapes[i].mesh.positions[3*v+0] >= -1.0 - epsilon);
         assert(shapes[i].mesh.positions[3*v+0] <= 1.0 + epsilon);
         shapes[i].mesh.positions[3*v+1] = (shapes[i].mesh.positions[3*v+1] - shiftY) * scaleY;
         assert(shapes[i].mesh.positions[3*v+1] >= -1.0 - epsilon);
         assert(shapes[i].mesh.positions[3*v+1] <= 1.0 + epsilon);
         shapes[i].mesh.positions[3*v+2] = (shapes[i].mesh.positions[3*v+2] - shiftZ) * scaleZ;
         assert(shapes[i].mesh.positions[3*v+2] >= -1.0 - epsilon);
         assert(shapes[i].mesh.positions[3*v+2] <= 1.0 + epsilon);
      }
   }
}

void setBoundingBox(struct BoundingBox *bb, float xmin, float xmax, float ymin, float ymax) {
   bb->xmin = xmin;
   bb->xmax = xmax;
   bb->ymin = ymin;
   bb->ymax = ymax;
}

float findMin(float c1, float c2, float c3) {
    float minimum = min(c1, min(c2, c3));
    return minimum;
}

float findMax(float c1, float c2, float c3) {
    float maximum = max(c1, max(c2, c3));
    return maximum;
}

void calcBounds( struct Triangle t, struct BoundingBox *bb ) {
    int xmin = findMin(t.v1.x, t.v2.x, t.v3.x);
    int xmax = findMax(t.v1.x, t.v2.x, t.v3.x);
    int ymin = findMin(t.v1.y, t.v2.y, t.v3.y);
    int ymax = findMax(t.v1.y, t.v2.y, t.v3.y);
    setBoundingBox(bb, xmin, xmax, ymin, ymax);
}

bool inTriangle(float val) {
    return ((val >= -epsilon) && (val <= (1+epsilon)));
}

float scaleIntensity(float z) { // Returns scalar to multiply color by
    // 1 is 100%, 0.5 is 75%, 0 is 50%, -0.5 is 25%, -1 is 0%
    return (((50*z) + 50)/100);
}

void colorPixel(int x, int y, float alpha, float beta, float gamma, shared_ptr<Image> image, struct Triangle t) {
    if(inTriangle(alpha) && inTriangle(beta) && inTriangle(gamma)) {
        float currZ = (alpha*t.v1.z)+(beta*t.v2.z)+(gamma*t.v3.z);
        float r = (alpha*t.c1.r)+(beta*t.c2.r)+(gamma*t.c3.r);
        float g = (alpha*t.c1.g)+(beta*t.c2.g)+(gamma*t.c3.g);
        float b = (alpha*t.c1.b)+(beta*t.c2.b)+(gamma*t.c3.b);
        // float r = 255;
        // float g = 0;
        // float b = 0;
        // cout << "PIXEL COLORS" << endl;
        // cout << "X: " << x << endl;
        // cout << "Y: " << y << endl;
        if (zBuf[x][y] < currZ) {
            image->setPixel(x, y, r, g, b);
            zBuf[x][y] = currZ;
        }
    }
}

void calcBary(struct Triangle t, int oldX, int oldY, shared_ptr<Image> image, float tArea) {
    float xaxc = t.v1.x - t.v3.x;
    float yayc = t.v1.y - t.v3.y;
    float xbxa = t.v2.x - t.v1.x;
    float ybya = t.v2.y - t.v1.y;
    float bArea = (xaxc*(oldY-t.v3.y)) - ((oldX-t.v3.x)*yayc);
    float gArea = (xbxa*(oldY-t.v1.y)) - ((oldX-t.v1.x)*ybya);
    float beta = bArea/tArea;
    float gamma = gArea/tArea;
    float alpha = 1 - beta - gamma;
    colorPixel(oldX, oldY, alpha, beta, gamma, image, t);
}

// g_width = input width; g_height = input height
// r_right = w/h; r_left = -w/h
// r_bottom = -1; r_top = 1
// xw = current x; yw = current y (in world)
void transformCoords(float xw, float yw, float pixels[]) {
    // cout << "XW: " << xw << endl;
    // cout << "YW: " << yw << endl;
    float c = (g_width-1)/(r_right-r_left);
    float d = c * r_left * -1;
    float xp = (c * xw) + d;
    pixels[0] = xp;

    float e = (g_height-1)/(r_top-r_bottom);
    float f = e * r_bottom * -1;
    float yp = (e * yw) + f;
    pixels[1] = yp;
    // cout << "XP: " << xp << endl;
    // cout << "YP: " << yp << endl;
}

void setVertProps(struct Vertex *v, float x, float y, float z) {
    // cout << "Z: " << z << endl;
    v->x = x;
    v->y = y;
    v->z = z;
}

void setColorProps(struct Color *c, float r, float g, float b) {
    // cout << "RED: " << r << endl;
    // cout << "GREEN: " << g << endl;
    // cout << "BLUE: " << b << endl;
    c->r = r;
    c->g = g;
    c->b = b;
}

void handleColors(struct Triangle *t) {
    //cout << "HANDLING COLORS" << endl;
    struct Color c1;
    struct Color c2;
    struct Color c3;
    if (colorMode == 1) { // depth
        float scalar = scaleIntensity(t->v1.z);
        // cout << "Z1: " << t->v1.z << endl;
        // cout << "SCALAR 1: " << scalar << endl;
        setColorProps(&c1, (global_r*scalar), (global_g*scalar), (global_b*scalar));
        scalar = scaleIntensity(t->v2.z);
        // cout << "Z2: " << t->v2.z << endl;
        // cout << "SCALAR 2: " << scalar << endl;
        setColorProps(&c2, (global_r*scalar), (global_g*scalar), (global_b*scalar));
        scalar = scaleIntensity(t->v3.z);
        // cout << "Z3: " << t->v3.z << endl;
        // cout << "SCALAR 3: " << scalar << endl;
        setColorProps(&c3, (global_r*scalar), (global_g*scalar), (global_b*scalar));
    }
    t->c1 = c1;
    t->c2 = c2;
    t->c3 = c3;
}

Vertex handleSetup(unsigned int ind) {
    struct Vertex v;
    float pixels[2];
    vector<float> vCoords = parsedT[ind];
    transformCoords(vCoords[0], vCoords[1], pixels);
    setVertProps(&v, pixels[0], pixels[1], vCoords[2]);
    return v;
}

struct Triangle createTriangle(unsigned int ind1, unsigned int ind2, unsigned int ind3, vector<float> posBuf) {
    struct Triangle t;
    t.v1 = handleSetup(ind1);
    t.v2 = handleSetup(ind2);
    t.v3 = handleSetup(ind3);
    handleColors(&t);
    return t;
}

void drawTriangle(struct Triangle t, shared_ptr<Image> image) {
    calcBounds(t, &bb);

    // Calculate area of triangle
    float tArea = ((t.v2.x-t.v1.x)*(t.v3.y-t.v1.y)) - ((t.v3.x-t.v1.x)*(t.v2.y-t.v1.y));

	// Draw bounding box
	for(int y = bb.ymin; y <= bb.ymax; ++y) {
		for(int x = bb.xmin; x <= bb.xmax; ++x) {
            calcBary(t, x, y, image, tArea);
		}
	}
}

vector<vector<float>> parseTriangles(vector<float> posBuf) {
    vector<vector<float>> parsedT;
    for (int i=0; i < posBuf.size(); i += 3) {
        vector<float> tmp;
        tmp.push_back(posBuf[i]);
        tmp.push_back(posBuf[i+1]);
        tmp.push_back(posBuf[i+2]);
        parsedT.push_back(tmp);
    } 
    return parsedT;
}

/* basic main - modify */
int main(int argc, char **argv)
{
	if(argc < 6) {
      cerr << "Usage: ./raster meshfile imagefile image_width image_height coloring_mode" << endl;
      return 0;
    }
	// OBJ filename
	string meshName(argv[1]);
	string imgName(argv[2]);

	g_width = atoi(argv[3]);
    g_height = atoi(argv[4]);
    colorMode = atoi(argv[5]);

    if (g_width < g_height) {
        r_right = 1;
        r_left = -1;
        r_top = float(g_height)/float(g_width);
        r_bottom = r_top * -1;
    } else {
        r_right = float(g_width)/float(g_height);
        r_left = r_right * -1;
        r_top = 1;
        r_bottom = -1;
    }

   //create an image
	auto image = make_shared<Image>(g_width, g_height);

	// triangle buffer to store indices to which vertices together form a triangle
	vector<unsigned int> triBuf;
	// position buffer to store vertex positions {x, y, z}
	vector<float> posBuf;
    float inf = std::numeric_limits<float>::infinity();
    inf *= -1;
    for (int i = 0; i < g_width; i++) {
        vector<float> tmp;
        for (int j = 0; j < g_height; j++) {
            tmp.push_back(inf);
        }
        zBuf.push_back(tmp);
    }
	// Some obj files contain material information.
	// We'll ignore them for this assignment.
	vector<tinyobj::shape_t> shapes; // geometry
	vector<tinyobj::material_t> objMaterials; // material
	string errStr;
	
   //load in the obj file - shapes has the key information
   bool rc = tinyobj::LoadObj(shapes, objMaterials, errStr, meshName.c_str());

	/* error checking on read */
	if(!rc) {
		cerr << errStr << endl;
	} else {
 		//keep this code to resize your object to be within -1 -> 1. 
        //then you will still WANT a window transform
   	    resize_obj(shapes); 
        //slightly redundant but makes some consistency with later OGL code
        //access the mesh data via posBuf and triBuf 
        // position buffer to store vertex positions {x, y, z}
        // triangle buffer to store indices to which vertices together form a triangle
		posBuf = shapes[0].mesh.positions;
		triBuf = shapes[0].mesh.indices;
	}
	cout << "Number of vertices: " << posBuf.size()/3 << endl;
	cout << "Number of triangles: " << triBuf.size()/3 << endl;

	//TODO add code to iterate through each triangle and rasterize it 
    //BECAREFUL this step is the most 'bug' prone for students - data is PACKED into a single array

    parsedT = parseTriangles(posBuf); // Group vertices

    for (int i=0; i < triBuf.size(); i += 3) {
        Triangle tmp = createTriangle(triBuf[i], triBuf[i+1], triBuf[i+2], posBuf);
        drawTriangle(tmp, image);
    }
	
	//write out the image
    image->writeToFile(imgName);

	return 0;
}
