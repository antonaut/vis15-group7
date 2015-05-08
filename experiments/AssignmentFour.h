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
	Vector2f ExampleFieldValue(Vector2f);
	void DrawVectorFieldHelper();
	bool IsTooSlow(Vector2f);
	float32 randomFloat(float32, float32);
	vector<Vector2f> magnitudeDistribution(int n);
	void magnitudeDistributionHelper(int n, float32 minX, float32 maxX, float32 minY, float32 maxY, vector<Vector2f> &);
	float32 getMagnitude(Vector2f);
	bool magSort(const pair<float32, Vector2f> &m1, const pair<float32, Vector2f> &m2);

    //Attrs

	///The method used to get vector field data.
	Vector2f(AssignmentFour::*VectorFieldAccessor)(Vector2f);

	///Should we use the vector field or the example field
	bool UseVectorField;
	float arcLength;

//Constructor / Destructor
public:
    AssignmentFour();
    virtual ~AssignmentFour();

//Methods
public:
    void LoadandRefreshVectorField();
	void EulerStreamline();
	void RungeKuttaStreamline();
	void SeedingStreamLines();
	void DistributionSeed();

	void UseEllipseField();
    void DrawVectorField();
	void DrawStreamline(vector<Vector2f>, const Vector4f &color);

	Vector2f Method(Vector2f);
	vector<Vector2f> Integrator(int, Vector2f(AssignmentFour::*)(Vector2f), float32 x, float32 y);
    virtual QWidget* createViewer();

//Attributes
public:
    ///The loaded field
    VectorField2 Field;

    ///File name of the vector field
    string VectorfieldFilename;

	bool NormalizeVectorField;
	bool IntegrateBackwards;
	bool ShowPoints;

 
	float XStart, YStart, MaxDistance;
	//Euler
	float EulerStepSize;
	int EulerStep;
	//RK
	float RKStepSize;
	float MaxArchLength;
	int RKStep;

	bool DirectionFieldOnly;

	float ArrowScale;

	// Seeding Stream Lines
	int NumStreamLines;
	bool DrawField;
	bool GridSeed;
	int GridPointsX;
	int GridPointsY;
	
protected:
    GLGeometryViewer* viewer;

};

#endif
