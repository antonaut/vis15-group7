//---------------------------------------------------------------------------
#ifndef AssignmentSixH
#define AssignmentSixH
//---------------------------------------------------------------------------
#include "Experiment.h"
#include "LinearAlgebra.h"
#include "Field2.hpp"
#include "GLGeometryViewer.h"
//---------------------------------------------------------------------------


/// This is an assignment experiment.
///
class AssignmentSix : public Experiment
{
	GEOX_CLASS(AssignmentSix)

private:
	//Methods
	void LoadVectorField();
	bool IsTooSlow(Vector2f);
	Vector2f RK4(Vector2f);

	float32 randomFloat(float32, float32);
	//Attrs

	///The method used to get vector field data.
	Vector2f(AssignmentSix::*VectorFieldAccessor)(Vector2f);

	float arcLength;

	//Constructor / Destructor
public:
	AssignmentSix();
	virtual ~AssignmentSix();

	//Methods
public:
	void LoadVectorFieldAndRefresh();
	void EulerStreamline();
	void RungeKuttaStreamline();

	Vector2f FieldValue(Vector2f);
	Vector2f Method(Vector2f);
	vector<Vector2f> Integrator(int, Vector2f(AssignmentSix::*)(Vector2f), float32 x, float32 y);
	virtual QWidget* createViewer();

	//Attributes
public:
	///The loaded field
	VectorField2 Field;

	///File name of the vector field
	string VectorfieldFilename;

	bool IntegrateBackwards;
	bool DirectionFieldOnly;

	//RK
	float RKStepSize;
	int RKStep;

	int32 SampleX;
	int32 SampleY;
	int32 KernelSize;
	int32 Seed;

protected:
	GLGeometryViewer* viewer;

};

#endif
