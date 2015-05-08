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
	Vector2f RK4(Vector2f, bool integrateBackwards);
	ScalarField2 enhanceContrast(ScalarField2);
	float32 randomFloat(float32, float32);

	ScalarField2 getRandomField(const Vector2f &lowerBounds, const Vector2f &upperBounds,
								const Vector2ui &dims, bool grayscale);

	vector< vector<Vector2f> > getStreamLines(const VectorField2 &field);
	vector<Vector2ui> streamLineToPixels(const VectorField2 &field, const vector<Vector2f> &streamLine);
	float32 smear(const ScalarField2 &field, const vector<Vector2ui> &pixels);
	//Attrs
	ScalarField2 texture;
	///The method used to get vector field data.
	Vector2f(AssignmentSix::*VectorFieldAccessor)(Vector2f, bool);

	float arcLength;

	//Constructor / Destructor
public:
	AssignmentSix();
	virtual ~AssignmentSix();

	//Methods
public:
	void LoadVectorFieldAndRefresh();
	void RungeKuttaStreamline(float32 xstart, float32 ystart);
	void DrawTexture();
	void EnhanceTexture();
	void LIC();
	void Coloring();

	Vector2f FieldValue(Vector2f, bool integrateBackwards);
	vector<Vector2f> Integrator(int, Vector2f(AssignmentSix::*)(Vector2f, bool), float32 x, float32 y);
	virtual QWidget* createViewer();

	//Attributes
public:
	///The loaded field
	VectorField2 Field;

	///File name of the vector field
	string VectorfieldFilename;
	string TextureFilename;

	//RK
	float RKStepSize;
	int RKStep;

	int32 SampleX;
	int32 SampleY;
	int32 KernelSize;
	int32 Seed;

	///Whether to draw the texture in RGB or grayscale
	bool ColoredTexture;

protected:
	GLGeometryViewer* viewer;

};

#endif
