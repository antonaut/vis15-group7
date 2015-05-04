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
	//ADD_NOARGS_METHOD(AssignmentThree::DrawVectorField)
	
	ADD_NOARGS_METHOD(AssignmentFour::UseEllipseField)
	ADD_NOARGS_METHOD(AssignmentFour::LoadandRefreshVectorField)
	ADD_NOARGS_METHOD(AssignmentFour::EulerStreamlines)
	ADD_NOARGS_METHOD(AssignmentFour::RungeKuttaStreamlines)
	ADD_NOARGS_METHOD(AssignmentFour::GoodStepSize)


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
	VectorfieldFilename = "./data/assignment05/Sink.am";
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
	// The method used to get vector field data.
	VectorFieldAccessor = &AssignmentFour::ExampleFieldValue;
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
			vec.normalize();

			viewer->addLine(x, y, x + ArrowScale*vec[0], y + ArrowScale*vec[1]);
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
	Vector2f(AssignmentFour::*Method)(Vector2f)
)
{
	Vector2f xi;
	vector<Vector2f> path;
	float32 arcLength = 0.0f;

	xi[0] = XStart;
	xi[1] = YStart;
	xi.normalize();
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


void AssignmentFour::EulerStreamlines()
{
	vector<Vector2f> path = Integrator(EulerStep, &AssignmentFour::Euler);

	DrawStreamlines(path);
}
void AssignmentFour::RungeKuttaStreamlines()
{
	vector<Vector2f> path = Integrator(RKStep, &AssignmentFour::RK4);
	DrawStreamlines(path);
}

void AssignmentFour::GoodStepSize()
{
	/* TODO */
}

void AssignmentFour::DrawStreamlines(vector<Vector2f> path)
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
		viewer->addLine(p1[0],p1[1], p2[0], p2[1]);
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


