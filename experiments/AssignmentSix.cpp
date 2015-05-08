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

	ADD_SEPARATOR("Texture");
	ADD_STRING_PROP(TextureFilename, 0)

	ADD_SEPARATOR("Options")
	ADD_INT32_PROP(SampleX, 0)
	ADD_INT32_PROP(SampleY, 0)
	ADD_INT32_PROP(KernelSize, 0)
	ADD_INT32_PROP(Seed, 0)
	ADD_NOARGS_METHOD(AssignmentSix::LoadVectorFieldAndRefresh)

	ADD_NOARGS_METHOD(AssignmentSix::DrawTexture)
	ADD_NOARGS_METHOD(AssignmentSix::LIC)


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
	VectorfieldFilename = "/home/simon/Git/vis15-group7/data/assignment05/ANoise2CT4.am";
	TextureFilename = "/home/simon/Git/vis15-group7/data/assignment06/";

	RKStepSize = 0.3;
	RKStep = 30;

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

void AssignmentSix::DrawTexture() {
	viewer->clear();

	srand((unsigned) Seed);

	ScalarField2 field = getRandomField(makeVector2f(-5, -5), makeVector2f(5, 5), makeVector2ui(512, 512), false);
	viewer->setTextureGray(field.getData());

	viewer->refresh();
}

void AssignmentSix::LIC() {
	viewer->clear();

	srand((unsigned) Seed);
	LoadVectorField();

	vector< vector<Vector2f> > streamLines = getStreamLines(Field);
	ScalarField2 randomField = getRandomField(makeVector2f(-5, -5), makeVector2f(5, 5), makeVector2ui(128, 128), false);

	viewer->setTextureGray(Field.getData());

	viewer->refresh();
}

ScalarField2 AssignmentSix::getRandomField(const Vector2f &lowerBounds, const Vector2f &upperBounds,
										   const Vector2ui &dims, bool grayscale) {
	ScalarField2 field = ScalarField2();
	field.init(lowerBounds, upperBounds, dims);

	for (card32 i = 0; i < dims[0]; ++i) {
		for (card32 j = 0; j < dims[1]; ++j) {
			float32 value = randomFloat(0, 1);
			if (grayscale) {
				value = round(value);
			}

			field.setNodeScalar(i, j, value);
		}
	}

	return field;
}


vector<vector<Vector2f> > AssignmentSix::getStreamLines(const VectorField2 &field) {

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

Vector2f AssignmentSix::FieldValue(Vector2f xi, bool integrateBackwards) {
	StaticVector<float, 2U> vec = Field.sample(xi[0], xi[1]);
	Vector2f v = makeVector2f(vec[0], vec[1]);

	return integrateBackwards ? -v : v;
}

Vector2f AssignmentSix::RK4(Vector2f xi, bool integrateBackwards)
{
	Vector2f v1, v2, v3, v4;
	v1 = (this->*VectorFieldAccessor)(xi, integrateBackwards);
	v2 = (this->*VectorFieldAccessor)(xi + (v1 * (RKStepSize / 2.0f)), integrateBackwards);
	v3 = (this->*VectorFieldAccessor)(xi + (v2 * (RKStepSize / 2.0f)), integrateBackwards);
	v4 = (this->*VectorFieldAccessor)(xi + (v3 * RKStepSize), integrateBackwards);

	return xi + (v1 + v2 * 2.0f + v3 * 2.0f + v4) * RKStepSize / 6.0f;
}


void AssignmentSix::RungeKuttaStreamline(float32 xstart, float32 ystart)
{
	Vector4f color = makeVector4f(1, 0, 1, 1);
	vector<Vector2f> path = Integrator(32, AssignmentSix::RK4, xstart, ystart);
}

float32 AssignmentSix::randomFloat(float32 a, float32 b) {
	float32 random = ((float32)rand()) / (float32)RAND_MAX;
	float32 diff = b - a;
	float32 r = random * diff;
	return a + r;
}
