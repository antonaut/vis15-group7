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


	ADD_SEPARATOR("Euler")
	ADD_FLOAT32_PROP(EulerStepSize, 0)
	ADD_INT32_PROP(EulerStep, 0)

	ADD_SEPARATOR("Runge Kutta")
	ADD_FLOAT32_PROP(RKStepSize, 0)
	ADD_INT32_PROP(RKStep, 0)

    ADD_SEPARATOR("Vectorfield")
    ADD_STRING_PROP(VectorfieldFilename, 0)
   
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
   
	VectorfieldFilename = "C:\\Users\\Eyob\\Desktop\\vis15-group7\\data\\assignment05\\Sink.am";

	XStart = 1;
	YStart = 0;
	MaxDistance = 5.3;

	EulerStepSize = 0.1;
	EulerStep = 100;
	
	RKStepSize = 0.3;
	RKStep = 30;

}

AssignmentFour::~AssignmentFour() {}

void AssignmentFour::LoadVectorField() {
        //Load scalar field
    if (!Field.load(VectorfieldFilename))
    {
        error("Error loading field file " + VectorfieldFilename + "\n");
    }
}

void AssignmentFour::LoadandRefreshVectorField() {
    viewer->clear();

    LoadVectorField();

    viewer->refresh();
}

vector<Vector2f> AssignmentFour::Integrator(int step, Vector2f (AssignmentFour::*Method)(Vector2f))
{
	Vector2f xi;
	vector<Vector2f> path;

	xi[0] = XStart;
	xi[1] = YStart;
	xi.normalize();
	path.push_back(xi);

	for(int i=0; i<step; i++)
	{
		Vector2f xp = (this->*Method)(xi);
		path.push_back(xp);
		xi = xp;
	}

	return path;
}

Vector2f AssignmentFour::FieldValue(Vector2f xi) {
	StaticVector<float, 1U> vec = Field.sample(xi[0], xi[1]);
	float x = vec[0];
	float y = vec[1];
	return makeVector2f(x,y);
}

Vector2f AssignmentFour::Euler(Vector2f xi)
{
	Vector2f xp = xi + FieldValue(xi)*EulerStepSize;
	return xp;
}

Vector2f AssignmentFour::RK4(Vector2f xi)
{
	Vector2f v1, v2, v3, v4;
	v1 = FieldValue(xi)/6.0f;
	v2 = FieldValue(xi + (v1 * (RKStepSize/2.0f)))/3.0f;
	v3 = FieldValue(xi + (v2 * (RKStepSize/2.0f)))/3.0f;
	v4 = FieldValue(xi + (v3 * RKStepSize))/6.0f;
	return xi + (v1 + v2 + v3 + v4) * RKStepSize;
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

    //Draw vector directions (constant length)
    for(float32 x=field.boundMin()[0]; x<=field.boundMax()[0]; x+=0.2)
    {
        for(float32 y=field.boundMin()[1]; y<=field.boundMax()[1]; y+=0.2)
        {
            Vector2f vec = field.sample(x,y);
            vec.normalize();

            //viewer->addLine(x, y, x + ArrowScale*vec[0], y + ArrowScale*vec[1]);
        }
    }

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


