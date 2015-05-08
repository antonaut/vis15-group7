//---------------------------------------------------------------------------
#include "stdafx.h"
//---------------------------------------------------------------------------
#include "AssignmentSix.h"
//---------------------------------------------------------------------------
#include "Properties.h"
#include "GeoXOutput.h"
//---------------------------------------------------------------------------

#include <limits>
#include <vector>


IMPLEMENT_GEOX_CLASS(AssignmentSix, 0)
{
	BEGIN_CLASS_INIT(AssignmentSix);

	ADD_SEPARATOR("Vectorfield")
	ADD_STRING_PROP(VectorfieldFilename, 0)

	ADD_SEPARATOR("Options")
	ADD_INT32_PROP(SampleX, 0)
	ADD_INT32_PROP(SampleY, 0)
	ADD_INT32_PROP(KernelSize, 0)
	ADD_INT32_PROP(Seed, 0)
	ADD_NOARGS_METHOD(AssignmentSix::LoadVectorFieldAndRefresh)
}

QWidget* AssignmentSix::createViewer()
{
	viewer = new GLGeometryViewer();
	return viewer;
}

AssignmentSix::AssignmentSix()
{
	viewer = NULL;

	//VectorfieldFilename = "C:\\Users\\Eyob\\Desktop\\Sink.am";
	VectorfieldFilename = "./vis15-group7/data/assignment05/ANoise2CT4.am";

	RKStepSize = 0.3;
	RKStep = 30;

	DirectionFieldOnly = false;
	IntegrateBackwards = false;
	
	VectorFieldAccessor = &AssignmentSix::FieldValue;

	arcLength = 0.0f;
}

AssignmentSix::~AssignmentSix() {}


void AssignmentSix::LoadVectorField() {

	if (!Field.load(VectorfieldFilename))
	{
		error("Error loading field file " + VectorfieldFilename + "\n");
	}
	VectorFieldAccessor = &AssignmentSix::FieldValue;
}


void AssignmentSix::LoadVectorFieldAndRefresh() {
	viewer->clear();

	LoadVectorField();

	viewer->refresh();
}

bool AssignmentSix::IsTooSlow(Vector2f vec) {
	float threshold = 1E-4;
	return vec.getSqrNorm() < threshold;
}

/**
Takes two parameters. The number of steps and a method of integration (i.e. Euler or RK4).
Returns a vector with a stream line path.
*/
vector<Vector2f>
AssignmentSix::Integrator(
int numberOfSteps,
Vector2f(AssignmentSix::*Method)(Vector2f),
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

		if (!Field.insideBounds(xp)) {
			output << "Stopped early after " << i << " steps. (Outside bounds)\n";
			return path;
		}

		path.push_back(xp);
		arcLength += (xp - xi).getSqrNorm();

		/*
		if (arcLength > MaxDistance) {
			output << "Stopped early after " << i << " steps. (Maximum distance)\n";
			return path;
		}*/

		xi = xp;
	}

	return path;
}

Vector2f AssignmentSix::FieldValue(Vector2f xi) {
	StaticVector<float, 2U> vec = Field.sample(xi[0], xi[1]);
	Vector2f v = makeVector2f(vec[0], vec[1]);
	return IntegrateBackwards ? -v : v;
}

Vector2f AssignmentSix::RK4(Vector2f xi)
{
	Vector2f v1, v2, v3, v4;
	v1 = (this->*VectorFieldAccessor)(xi);
	if (DirectionFieldOnly) v1.normalize();

	v2 = (this->*VectorFieldAccessor)(xi + (v1 * (RKStepSize / 2.0f)));
	if (DirectionFieldOnly) v2.normalize();

	v3 = (this->*VectorFieldAccessor)(xi + (v2 * (RKStepSize / 2.0f)));
	if (DirectionFieldOnly) v3.normalize();

	v4 = (this->*VectorFieldAccessor)(xi + (v3 * RKStepSize));
	if (DirectionFieldOnly) v4.normalize();

	return xi + (v1 + v2 * 2.0f + v3 * 2.0f + v4) * RKStepSize / 6.0f;
}


void AssignmentSix::RungeKuttaStreamline()
{
	Vector4f color = makeVector4f(1, 0, 1, 1);
	float xstart, ystart;
	vector<Vector2f> path = Integrator(RKStep, &AssignmentSix::RK4, xstart, ystart);
	/*DrawStreamline(path, color);*/
}

float32 AssignmentSix::randomFloat(float32 a, float32 b) {
	float32 random = ((float32)rand()) / (float32)RAND_MAX;
	float32 diff = b - a;
	float32 r = random * diff;
	return a + r;
}
