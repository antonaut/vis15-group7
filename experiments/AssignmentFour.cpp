<<<<<<< HEAD
//---------------------------------------------------------------------------
#include "stdafx.h"
//---------------------------------------------------------------------------
#include "AssignmentFour.h"
//---------------------------------------------------------------------------
#include "Properties.h"
#include "GLGeometryViewer3D.h"
#include "MarchingCubes3D.h"
//---------------------------------------------------------------------------
#include <fstream>

#define rnd01() (((double) rand()) / RAND_MAX)
#define SIDE_LENGTH 1.0f

IMPLEMENT_GEOX_CLASS( AssignmentFour ,0) {
    BEGIN_CLASS_INIT( AssignmentFour );
    
    ADD_STRING_PROP(Filename, 0)
    ADD_FLOAT32_PROP(IsoValue, 0)

    ADD_NOARGS_METHOD(AssignmentFour::DrawCubes);
    ADD_NOARGS_METHOD(AssignmentFour::showPerVertexLighting)
    ADD_NOARGS_METHOD(AssignmentFour::applyMarchingCubes)
}

QWidget *AssignmentFour::createViewer() {
   viewer = new GLGeometryViewer3D();
   return viewer;
}

AssignmentFour::AssignmentFour() {
    IsoValue = 0;
    Filename = "data/cubes.txt";
    viewer = NULL;
}

AssignmentFour::~AssignmentFour() {}

void AssignmentFour::loadCubes() {
    fstream in(Filename.c_str());

    size_t num_cubes;
    in >> num_cubes;

    cubes.clear();
    for (size_t i = 0; i < num_cubes; ++i) {
        cubes.push_back(Cube());
        in >> cubes[i];
    }
}

void AssignmentFour::setCorners(const Cube &c, const Vector3f &offset, Vector3f corners[8]) {
    float32 r = SIDE_LENGTH/2;
    float32 x1 = -r + offset[0], x2 = r + offset[0];
    float32 y1 = -r + offset[1], y2 = r + offset[1];
    float32 z1 = -r + offset[2], z2 = r + offset[2];

    corners[0] = makeVector3f(x1, y1, z2);
    corners[1] = makeVector3f(x2, y1, z2);
    corners[2] = makeVector3f(x2, y1, z1);
    corners[3] = makeVector3f(x1, y1, z1);
    corners[4] = makeVector3f(x1, y2, z2);
    corners[5] = makeVector3f(x2, y2, z2);
    corners[6] = makeVector3f(x2, y2, z1);
    corners[7] = makeVector3f(x1, y2, z1);
}

void AssignmentFour::drawCube(const Cube &c, const Vector3f &offset) {
    Vector3f corners[8];
    setCorners(c, offset, corners);

    const Vector3f &v1 = corners[0];
    const Vector3f &v2 = corners[1];
    const Vector3f &v3 = corners[2];
    const Vector3f &v4 = corners[3];
    const Vector3f &v5 = corners[4];
    const Vector3f &v6 = corners[5];
    const Vector3f &v7 = corners[6];
    const Vector3f &v8 = corners[7];

    viewer->addLine(v1, v2);
    viewer->addLine(v2, v3);
    viewer->addLine(v3, v4);
    viewer->addLine(v4, v1);

    viewer->addLine(v1, v5);
    viewer->addLine(v2, v6);
    viewer->addLine(v3, v7);
    viewer->addLine(v4, v8);

    viewer->addLine(v5, v6);
    viewer->addLine(v6, v7);
    viewer->addLine(v7, v8);
    viewer->addLine(v8, v5);
}

void AssignmentFour::drawContours(const Cube &c, const Vector3f &offset) {
    Vector3f corners[8];
    setCorners(c, offset, corners);

    vector<Vector3f> vertices;

    MarchingCubes3D::triangulate(corners, c.v, IsoValue, vertices);
    for (size_t i = 0; i < vertices.size() - 2; i += 3) {
        viewer->addTriangle(vertices[i], vertices[i+1], vertices[i+2]);
    }
}


void AssignmentFour::DrawCubes() {
    loadCubes();
    viewer->clear();

    size_t num_cubes = cubes.size();
    float32 start_oy, start_ox;
    if (num_cubes == 1) {
        start_oy = 0.0f;
        start_ox = 0.0f;
    }
    else {
        start_oy = ((int) (sqrt(num_cubes) + 0.5f)) * 0.5f;
        start_ox = -start_oy;
    }

    vector<Vector3f> offsets;
    float32 oy = start_oy;
    for (int i = 0; i < (int) (sqrt(num_cubes) + 0.5f); ++i) {
        float ox = start_ox;
        for (int j = 0; j < (int) (sqrt(num_cubes) + 0.5f); ++j) {
            offsets.push_back(makeVector3f(ox, oy, 0));
            ox += 1.5f;
        }
        oy -= 1.5f;
    }

    for (size_t i = 0; i < cubes.size(); ++i) {
        drawCube(cubes[i], offsets[i]);
    }

    for (size_t i = 0; i < cubes.size(); ++i) {
        drawContours(cubes[i], offsets[i]);
    }

    viewer->refresh();
}

void AssignmentFour::showPerVertexLighting() 
{
    viewer->clear();

    //llf
    Point3D p0;
    p0.position = makeVector3f(-2.0f,-2.0f,-2.0f);
    p0.color = makeVector4f(0.0f,0.0f,0.0f,0.0f);
    p0.normal = -normalize(makeVector3f(-1.0f,-1.0f,-1.0f));
    int p0Handle = viewer->addPoint(p0);
    
    //lrf
    Point3D p1;
    p1.position = makeVector3f(2.0f,-2.0f,-2.0f);
    p1.color = makeVector4f(1.0f,0.0f,0.0f,0.0f);
    p1.normal = -normalize(makeVector3f(1.0f,-1.0f,-1.0f));
    int p1Handle = viewer->addPoint(p1);
    
    //urf
    Point3D p2;
    p2.position = makeVector3f(2.0f,2.0f,-2.0f);
    p2.color = makeVector4f(1.0f,1.0f,0.0f,0.0f);
    p2.normal = -normalize(makeVector3f(1.0f,1.0f,-1.0f));
    int p2Handle = viewer->addPoint(p2);
    
    //ulf
    Point3D p3;
    p3.position = makeVector3f(-2.0f,2.0f,-2.0f);
    p3.color = makeVector4f(0.0f,1.0f,0.0f,0.0f);
    p3.normal = -normalize(makeVector3f(-1.0f,1.0f,-1.0f));
    int p3Handle = viewer->addPoint(p3);
    
    //llb
    Point3D p4;
    p4.position = makeVector3f(-2.0f,-2.0f,2.0f);
    p4.color = makeVector4f(0.0f,0.0f,1.0f,0.0f);
    p4.normal = -normalize(makeVector3f(-1.0f,-1.0f,1.0f));
    int p4Handle = viewer->addPoint(p4);

    //lrb
    Point3D p5;
    p5.position = makeVector3f(2.0f,-2.0f,2.0f);
    p5.color = makeVector4f(1.0f,0.0f,1.0f,0.0f);
    p5.normal = -normalize(makeVector3f(1.0f,-1.0f,1.0f));
    int p5Handle = viewer->addPoint(p5);

    //urb
    Point3D p6;
    p6.position = makeVector3f(2.0f,2.0f,2.0f);
    p6.color = makeVector4f(1.0f,1.0f,1.0f,0.0f);
    p6.normal = -normalize(makeVector3f(1.0f,1.0f,1.0f));
    int p6Handle = viewer->addPoint(p6);

    //ulb
    Point3D p7;
    p7.position = makeVector3f(-2.0f,2.0f,2.0f);
    p7.color = makeVector4f(0.0f,1.0f,1.0f,0.0f);
    p7.normal = -normalize(makeVector3f(-1.0f,1.0f,1.0f));
    int p7Handle = viewer->addPoint(p7);

    // front
    Triangle t0;
    t0.vertices[0] = p0Handle;
    t0.vertices[1] = p1Handle;
    t0.vertices[2] = p2Handle;
    int t0Handle = viewer->addTriangle(t0);
    
    Triangle t1;
    t1.vertices[0] = p2Handle;
    t1.vertices[1] = p3Handle;
    t1.vertices[2] = p0Handle;
    int t1Handle = viewer->addTriangle(t1);

    // right
    Triangle t2;
    t2.vertices[0] = p1Handle;
    t2.vertices[1] = p5Handle;
    t2.vertices[2] = p6Handle;
    int t2Handle = viewer->addTriangle(t2);

    Triangle t3;
    t3.vertices[0] = p6Handle;
    t3.vertices[1] = p2Handle;
    t3.vertices[2] = p1Handle;
    int t3Handle = viewer->addTriangle(t3);
    
    // left
    Triangle t4;
    t4.vertices[0] = p0Handle;
    t4.vertices[1] = p3Handle;
    t4.vertices[2] = p7Handle;
    int t4Handle = viewer->addTriangle(t4);

    Triangle t5;
    t5.vertices[0] = p7Handle;
    t5.vertices[1] = p4Handle;
    t5.vertices[2] = p0Handle;
    int t5Handle = viewer->addTriangle(t5);

    // back
    Triangle t6;
    t6.vertices[0] = p6Handle;
    t6.vertices[1] = p5Handle;
    t6.vertices[2] = p4Handle;
    int t6Handle = viewer->addTriangle(t6);

    Triangle t7;
    t7.vertices[0] = p4Handle;
    t7.vertices[1] = p7Handle;
    t7.vertices[2] = p6Handle;
    int t7Handle = viewer->addTriangle(t7);

    // top
    Triangle t8;
    t8.vertices[0] = p3Handle;
    t8.vertices[1] = p2Handle;
    t8.vertices[2] = p6Handle;
    int t8Handle = viewer->addTriangle(t8);

    Triangle t9;
    t9.vertices[0] = p6Handle;
    t9.vertices[1] = p7Handle;
    t9.vertices[2] = p3Handle;
    int t9Handle = viewer->addTriangle(t9);

    // bottom
    Triangle t10;
    t10.vertices[0] = p5Handle;
    t10.vertices[1] = p1Handle;
    t10.vertices[2] = p0Handle;
    int t10Handle = viewer->addTriangle(t10);

    Triangle t11;
    t11.vertices[0] = p0Handle;
    t11.vertices[1] = p4Handle;
    t11.vertices[2] = p5Handle;
    int t11Handle = viewer->addTriangle(t11);

    // display changes
    viewer->refresh();
}

void AssignmentFour::applyMarchingCubes()
{
    viewer->clear();

    Vector3f corners[8];
    corners[0] = makeVector3f(-2.0f,-2.0f,-2.0f);
    corners[1] = makeVector3f( 2.0f,-2.0f,-2.0f);
    corners[2] = makeVector3f( 2.0f, 2.0f,-2.0f);
    corners[3] = makeVector3f(-2.0f, 2.0f,-2.0f);
    corners[4] = makeVector3f(-2.0f,-2.0f, 2.0f);
    corners[5] = makeVector3f( 2.0f,-2.0f, 2.0f);
    corners[6] = makeVector3f( 2.0f, 2.0f, 2.0f);
    corners[7] = makeVector3f(-2.0f, 2.0f, 2.0f);

    float32 values[] = {
        -1, 3, 2, 3,
        1, 2, -2, 2, 
    };

    vector<Vector3f> vertices;

    MarchingCubes3D::triangulate(corners, values, 0, vertices);
    for (unsigned i = 0; i < vertices.size(); i++ )
        viewer->addPoint( Point3D(vertices[i]) );

    for (int j = 0; j < (int) (viewer->getNumberOfPoints() / 3); j++)
    {
        Triangle t;
        for(int k = 0; k < 3;k++)
            t.vertices[k] = j*3+k;
        viewer->addTriangle(t);
    }

    // display changes
    viewer->refresh();

}
=======
//---------------------------------------------------------------------------
#include "stdafx.h"
//---------------------------------------------------------------------------
#include "AssignmentFour.h"
//---------------------------------------------------------------------------
#include "Properties.h"
#include "GeoXOutput.h"
//---------------------------------------------------------------------------

#include <limits>
#include <vector>


IMPLEMENT_GEOX_CLASS( AssignmentFour, 0)
{
    BEGIN_CLASS_INIT( AssignmentFour );

	ADD_SEPARATOR("Start Position")
	ADD_FLOAT32_PROP(XStart, 0)
	ADD_FLOAT32_PROP(YStart, 0)
	ADD_FLOAT32_PROP(MaxDistance, 0)

	ADD_SEPARATOR("Options")
	ADD_BOOLEAN_PROP(DirectionFieldOnly, 0)

	ADD_SEPARATOR("Euler")
	ADD_FLOAT32_PROP(EulerStepSize, 0)
	ADD_INT32_PROP(EulerStep, 0)

	ADD_SEPARATOR("Runge Kutta")
	ADD_FLOAT32_PROP(RKStepSize, 0)
	ADD_INT32_PROP(RKStep, 0)

    ADD_SEPARATOR("Vectorfield")
    ADD_STRING_PROP(VectorfieldFilename, 0)
	ADD_FLOAT32_PROP(ArrowScale, 0)
	ADD_BOOLEAN_PROP(NormalizeVectorField, 0);

    ADD_SEPARATOR("Seeding of Stream Lines")
	ADD_INT32_PROP(NumStreamLines, 0)
	ADD_BOOLEAN_PROP(DrawField, 0)
	ADD_BOOLEAN_PROP(GridSeed, 0)
	ADD_INT32_PROP(GridPointsX, 0)
	ADD_INT32_PROP(GridPointsY, 0)
	
	ADD_NOARGS_METHOD(AssignmentFour::DrawVectorField)
	
	ADD_NOARGS_METHOD(AssignmentFour::UseEllipseField)
	ADD_NOARGS_METHOD(AssignmentFour::LoadandRefreshVectorField)
	ADD_NOARGS_METHOD(AssignmentFour::EulerStreamline)
	ADD_NOARGS_METHOD(AssignmentFour::RungeKuttaStreamline)
	ADD_NOARGS_METHOD(AssignmentFour::GoodStepSize)
	ADD_NOARGS_METHOD(AssignmentFour::SeedingStreamLines)
	ADD_NOARGS_METHOD(AssignmentFour::DistributionSeed)


}

QWidget* AssignmentFour::createViewer()
{
    viewer = new GLGeometryViewer();
    return viewer;
}

AssignmentFour::AssignmentFour()
{
    viewer = NULL;
   
	//VectorfieldFilename = "C:\\Users\\Eyob\\Desktop\\Sink.am";
	VectorfieldFilename = "/home/simon/Git/vis15-group7/data/assignment05/ANoise2CT4.am";
	XStart = 1;
	YStart = 0;
	MaxDistance = 5.3;
	
	EulerStepSize = 0.1;
	EulerStep = 100;
	
	RKStepSize = 0.3;
	RKStep = 30;

	ArrowScale = 0.1;
	MaxArchLength = -1.0f;

	DirectionFieldOnly = false;
	
	UseVectorField = false;

	VectorFieldAccessor = &AssignmentFour::ExampleFieldValue;
	NormalizeVectorField = false;

	NumStreamLines = 1600;
	GridPointsX = 40;
	GridPointsY = 40;
	DrawField = true;
	GridSeed = false;
}

AssignmentFour::~AssignmentFour() {}


void AssignmentFour::UseEllipseField() {
	VectorFieldAccessor = &AssignmentFour::ExampleFieldValue;
	UseVectorField = false;
}

void AssignmentFour::LoadVectorField() {

    if (!Field.load(VectorfieldFilename))
    {
        error("Error loading field file " + VectorfieldFilename + "\n");
    }
	VectorFieldAccessor = &AssignmentFour::FieldValue;
	UseVectorField = true;
}

void AssignmentFour::DrawVectorFieldHelper() {
	//Draw vector directions (constant length)
	for (float32 x = Field.boundMin()[0]; x <= Field.boundMax()[0]; x += 0.2)
	{
		for (float32 y = Field.boundMin()[1]; y <= Field.boundMax()[1]; y += 0.2)
		{
			Vector2f vec = Field.sample(x, y);

			if (NormalizeVectorField) {
				vec.normalize();
			}

			float32 x2 = x + ArrowScale*vec[0];
			float32 y2 = y + ArrowScale*vec[1];

			viewer->addLine(x, y, x2, y2);
			viewer->addPoint(Point2D(x2, y2));
		}
	}
}

void AssignmentFour::LoadandRefreshVectorField() {
    viewer->clear();

    LoadVectorField();

    viewer->refresh();
}

bool AssignmentFour::IsTooSlow(Vector2f vec) {
	float threshold = 1E-4;
	return vec.getSqrNorm() < threshold;
}

/**
	Takes two parameters. The number of steps and a method of integration (i.e. Euler or RK4).
	Returns a vector with a stream line path.
*/
vector<Vector2f> 
AssignmentFour::Integrator(
	int numberOfSteps, 
	Vector2f(AssignmentFour::*Method)(Vector2f),
	float32 xstart, float32 ystart
)
{
	Vector2f xi;
	vector<Vector2f> path;
	float32 arcLength = 0.0f;

	xi[0] = xstart;
	xi[1] = ystart;
	//xi.normalize();
	path.push_back(xi);

	for (int i = 0; i < numberOfSteps; i++)
	{
		Vector2f xp = (this->*Method)(xi);
		if (IsTooSlow(xp)) {
			output << "Stopped early after " << i << " steps. (Going too slow)\n";
			return path;
		}

		if (UseVectorField && !Field.insideBounds(xp)) {
			output << "Stopped early after " << i << " steps. (Outside bounds)\n";
			return path;
		}

		if (DirectionFieldOnly) {
			xp.normalize();
		}

		path.push_back(xp);
		
		arcLength += (xp - xi).getSqrNorm();

		if (arcLength > MaxDistance) {
			output << "Stopped early after " << i << " steps. (Maximum distance)\n";
			return path;
		}
		
		xi = xp;
	}

	return path;
}

Vector2f AssignmentFour::FieldValue(Vector2f xi) {
	StaticVector<float, 2U> vec = Field.sample(xi[0], xi[1]);
	return makeVector2f(vec[0], vec[1]);
}

Vector2f AssignmentFour::ExampleFieldValue(Vector2f vec) {
	return makeVector2f(-vec[1], vec[0] / 2.0f);
}

Vector2f AssignmentFour::Euler(Vector2f xi)
{
	Vector2f xp = xi + (this->*VectorFieldAccessor)(xi)*EulerStepSize;
	return xp;
}

Vector2f AssignmentFour::RK4(Vector2f xi)
{
	Vector2f v1, v2, v3, v4;
	v1 = (this->*VectorFieldAccessor)(xi);
	v2 = (this->*VectorFieldAccessor)(xi + (v1 * (RKStepSize / 2.0f)));
	v3 = (this->*VectorFieldAccessor)(xi + (v2 * (RKStepSize / 2.0f)));
	v4 = (this->*VectorFieldAccessor)(xi + (v3 * RKStepSize));
	return xi + (v1 + v2 * 2.0f + v3 * 2.0f + v4) * RKStepSize / 6.0f;
}


void AssignmentFour::EulerStreamline()
{
	Vector4f color = makeVector4f(1, 1, 0, 1);

	output << "Calculating integration points...";
	vector<Vector2f> path = Integrator(EulerStep, &AssignmentFour::Euler, XStart, YStart);
	output << "done\n";

	output << "Drawing stream line...";
	DrawStreamline(path, color);
	output << "done\n";
}

void AssignmentFour::RungeKuttaStreamline()
{
	Vector4f color = makeVector4f(1, 0, 1, 1);
	
	vector<Vector2f> path = Integrator(RKStep, &AssignmentFour::RK4, XStart, YStart);
	DrawStreamline(path, color);
}

float32 AssignmentFour::randomFloat(float32 a, float32 b) {
    float32 random = ((float32) rand()) / (float32) RAND_MAX;
    float32 diff = b - a;
    float32 r = random * diff;
    return a + r;
}

float32 AssignmentFour::getMagnitude(Vector2f xi) {
	return (this->*VectorFieldAccessor)(xi).getSqrNorm();
}

void AssignmentFour::magnitudeDistributionHelper(int n, float32 minX, float32 maxX, float32 minY, float32 maxY, vector<Vector2f> &points) {
	float32 spread = (float32) (sqrt(n) * 1.2);

	float32 dx = (maxX - minX) / spread;
	float32 dy = (maxY - minY) / spread;

	if (n <= 0) {
		return;
	}

	vector< pair<float32, Vector2f> > magnitudes;

	for (int xi = 0; xi < n; ++xi) {
		for (int yi = 0; yi < n; ++yi) {
			float32 x = minX + xi*dx;
			float32 y = minY + yi*dy;
			Vector2f v = makeVector2f(x, y);
			float32 magnitude = getMagnitude(v);
			magnitudes.push_back(make_pair(magnitude, v));
		}
	}

	sort(magnitudes.begin(), magnitudes.end(), [](const pair<float32, Vector2f> &m1, const pair<float32, Vector2f> &m2) -> bool {
		return m1.first > m2.first;
	});

	for (size_t i = 0; i < (size_t) n; ++i) {
		points.push_back(magnitudes[i].second);
	}
}

vector<Vector2f> AssignmentFour::magnitudeDistribution(int n) {
	float32 minX = Field.boundMin()[0]; 
	float32 maxX = Field.boundMax()[0];
	float32 minY = Field.boundMin()[1]; 
	float32 maxY = Field.boundMax()[1];

	vector<Vector2f> points;

	magnitudeDistributionHelper(n, minX, maxX, minY, maxY, points);

	return points;
}


void AssignmentFour::GoodStepSize()
{
	// TODO
}

void AssignmentFour::SeedingStreamLines()
{
	viewer->clear();

	if (DrawField) {
		DrawVectorFieldHelper();
	}

	Vector4f color = makeVector4f(0.8f, 0.2f, 0.0f, 0.40f);
	float32 minX = Field.boundMin()[0]; 
	float32 maxX = Field.boundMax()[0];
	float32 minY = Field.boundMin()[1]; 
	float32 maxY = Field.boundMax()[1];
	
	if (GridSeed) {
		float32 dx = (maxX - minX) / GridPointsX;
		float32 dy = (maxY - minY) / GridPointsY;
		float32 x = minX; 

		for (int ix = 0; ix < GridPointsX; ++ix) {
			float32 y = minY;
			for (int iy = 0; iy < GridPointsY; ++iy) {
				vector<Vector2f> path = Integrator(RKStep, &AssignmentFour::RK4, x, y);
				DrawStreamline(path, color);				

				y += dy;
			}
			x += dx;
		}
	}
	else {
		int n = NumStreamLines;

		for (int i = 0; i < n; ++i) {
			float32 x = randomFloat(minX, maxX);
			float32 y = randomFloat(minY, maxY);

			vector<Vector2f> path = Integrator(RKStep, &AssignmentFour::RK4, x, y);
			DrawStreamline(path, color);
		}
	}

	viewer->refresh();
}

void AssignmentFour::DistributionSeed() {
	viewer->clear();

	if (DrawField) {
		DrawVectorFieldHelper();
	}

	Vector4f color = makeVector4f(0.8f, 0.2f, 0.0f, 0.40f);

	vector<Vector2f> points = magnitudeDistribution(NumStreamLines);

	for (size_t i = 0; i < points.size(); ++i) {
		float32 x = points[i][0];
		float32 y = points[i][1];
		vector<Vector2f> path = Integrator(RKStep, &AssignmentFour::RK4, x, y);
		DrawStreamline(path, color);

		output << "x: " << x << "\ty:" << y << "\n";
	}

	viewer->refresh();
}

void AssignmentFour::DrawStreamline(vector<Vector2f> path, const Vector4f &color)
{
	int arraySize = path.size();
	if (arraySize < 2) {
		return;
	}
	Vector2f p1 = path[0];
	Vector2f p2;
	for(int i=1; i<arraySize; i++)
	{
		p2 = path[i];
		viewer->addLine(p1[0], p1[1], p2[0], p2[1], color);
		viewer->addPoint(p1);
		p1 = p2;
	}
	viewer->addPoint(path[arraySize-1]);
}

void AssignmentFour::DrawVectorField()
{
    viewer->clear();

    //Load the vector field
    VectorField2 field;
    if (!field.load(VectorfieldFilename))
    {
        output << "Error loading field file " << VectorfieldFilename << "\n";
        return;
    }

	DrawVectorFieldHelper();

    viewer->refresh();
}


namespace
{
    ///Returns the next power-of-two
    int32 NextPOT(int32 n)
    {
        n--;
        n |= n >> 1;   // Divide by 2^k for consecutive doublings of k up to 32,
        n |= n >> 2;   // and then or the results.
        n |= n >> 4;
        n |= n >> 8;
        n |= n >> 16;
        n++;           // The result is a number of 1 bits equal to the number
                       // of bits in the original number, plus 1. That's the
                       // next highest power of 2.
        return n;
    }
}


>>>>>>> 649b2f2aa6a1009e3dbfe8184b8951cfa393642b
