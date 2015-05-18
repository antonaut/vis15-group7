//---------------------------------------------------------------------------
#ifndef AssignmentSevenH
#define AssignmentSevenH
//---------------------------------------------------------------------------
#include "Experiment.h"
#include "LinearAlgebra.h"
#include "Field2.hpp"
#include "GLGeometryViewer.h"
#include <random>
//---------------------------------------------------------------------------


/// This is an assignment experiment.
///
class AssignmentSeven : public Experiment
{
	GEOX_CLASS(AssignmentSeven)

private:
	// Help methods
	float32 randomFloat(float32, float32);
	Vector2f rotate(const Vector2f &, float32 angle) const;

	// Drawing help methods
	void drawStreamline(vector<Vector2f> path, const Vector4f &color);

	// Field creators
	ScalarField2 getRandomField(const Vector2f &boundMin, const Vector2f &boundMax, const Vector2ui &dims, bool grayscale);
	VectorField2 getEllipseField(const Vector2f &boundMin, const Vector2f &boundMax, const Vector2ui &dims) const;
	void LoadVectorField();
	void LoadScalarField();

	// StreamLine helpers
	bool IsTooSlow(Vector2f);
	Vector2f RK4(Vector2f, bool integrateBackwards);
	ScalarField2 enhanceContrast(ScalarField2);
	Vector2f FieldValue(Vector2f, bool integrateBackwards);
	vector<Vector2f> Integrator(int, Vector2f(AssignmentSeven::*)(Vector2f, bool), float32 x, float32 y);

	// Fast LIC helpers
	vector<Vector2ui> streamLineToPixels(const ScalarField2 &field, const vector<Vector2f> &streamLine);
	vector<Vector2ui> lineToPixels(const ScalarField2 &field, const Vector2f &v0, const Vector2f &v1) const;
	vector<float32> smear(const ScalarField2 &field, const vector<Vector2ui> &pixels);

	//Attrs
	ScalarField2 texture;
	std::mt19937_64 mt;

	///The method used to get vector field data.
	Vector2f(AssignmentSeven::*VectorFieldAccessor)(Vector2f, bool);

	float arcLength;

	//Constructor / Destructor
public:
	AssignmentSeven();
	virtual ~AssignmentSeven();

	//Methods
public:
	void LoadVectorFieldAndRefresh();
	void LoadScalarFieldAndRefresh();
	void ShowExtremePoints();
	void DrawTexture();
	void EnhanceTexture();
	void FastLIC();
	void Coloring();

	virtual QWidget* createViewer();

	//Attributes
public:
	///The loaded fields
	VectorField2 Field;

	ScalarField2 SField;

	///File name of the vector- and scalar fields
	string VectorfieldFilename;
	string ScalarfieldFilename;
	string TextureFilename;

	// The number sample points of the texture in both directions
	Vector2ui TextureResolution;

	//RK
	float RKStepSize;
	int RKStep;

	card32 KernelSize;
	int32 Seed;

	///Whether to draw the texture in RGB or grayscale
	bool GrayScale;
	bool ColoredTexture;

	// Add StreamLines, for debugging
	bool AddStreamLines;

protected:
	GLGeometryViewer* viewer;

};

#endif
