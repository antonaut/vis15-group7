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
	ADD_BOOLEAN_PROP(IntegrateBackwards, 0)
	ADD_BOOLEAN_PROP(ShowPoints, 0)

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
	VectorfieldFilename = "./vis15-group7/data/assignment05/ANoise2CT4.am";
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
	IntegrateBackwards = false;
	ShowPoints = false;

	UseVectorField = false;

	VectorFieldAccessor = &AssignmentFour::ExampleFieldValue;
	NormalizeVectorField = false;

	NumStreamLines = 1600;
	GridPointsX = 40;
	GridPointsY = 40;
	DrawField = true;
	GridSeed = false;

	arcLength = 0.0f;
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
			if (ShowPoints) {
				viewer->addPoint(Point2D(x2, y2));
			}
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
	arcLength = 0.0f;

	xi[0] = xstart;
	xi[1] = ystart;

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
	Vector2f v = makeVector2f(vec[0], vec[1]);
	if (DirectionFieldOnly) v.normalize();
	return IntegrateBackwards ? -v : v;
}

Vector2f AssignmentFour::ExampleFieldValue(Vector2f vec) {
	Vector2f v = makeVector2f(-vec[1], vec[0] / 2.0f);
	return IntegrateBackwards ? -v : v;
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

	//viewer->refresh();
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
		if (ShowPoints) {
			viewer->addPoint(p1);
		}
		p1 = p2;
		viewer->refresh();
	}
	if (ShowPoints) {
		viewer->addPoint(path[arraySize - 1]);
	}
	viewer->refresh();
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