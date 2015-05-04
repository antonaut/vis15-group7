//---------------------------------------------------------------------------
#ifndef AssignmentFourH
#define AssignmentFourH
//---------------------------------------------------------------------------
#include "Experiment.h"
#include "LinearAlgebra.h"
#include "Field2.hpp"
#include "GLGeometryViewer.h"
//---------------------------------------------------------------------------


/// This is an assignment experiment.
///
class AssignmentFour : public Experiment
{
    GEOX_CLASS(AssignmentFour)

private:
    //Methods
    void LoadVectorField();
	Vector2f Euler(Vector2f);
	Vector2f RK4(Vector2f);
	Vector2f FieldValue(Vector2f);

    //Attrs

//Constructor / Destructor
public:
    AssignmentFour();
    virtual ~AssignmentFour();

//Methods
public:
    void LoadandRefreshVectorField();
	void EulerStreamlines();
	void RungeKuttaStreamlines();
	void GoodStepSize();

    void DrawVectorField();
	void DrawStreamlines(vector<Vector2f>);

	Vector2f Method(Vector2f);
	vector<Vector2f> Integrator(int, Vector2f (AssignmentFour::*)(Vector2f));
    virtual QWidget* createViewer();

//Attributes
public:
    ///The loaded field
    VectorField2 Field;

    ///File name of the vector field
    string VectorfieldFilename;

 
	float XStart, YStart, MaxDistance;
	//Euler
	float EulerStepSize;
	int EulerStep;
	//RK
	float RKStepSize;
	int RKStep;

	float ArrowScale;
	
protected:
    GLGeometryViewer* viewer;
};

#endif
