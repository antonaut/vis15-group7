//---------------------------------------------------------------------------
#ifndef AssignmentThreeH
#define AssignmentThreeH
//---------------------------------------------------------------------------
#include "Experiment.h"
#include "LinearAlgebra.h"
#include "GLGeometryViewer.h"
//---------------------------------------------------------------------------


/// This is an example experiment.
///
/// The code is meant to demonstrate how
///  to use the GeoX framework
///
class AssignmentThree : public Experiment
{
    GEOX_CLASS(AssignmentThree)

private:
    void AddContours(Point2D&, float32, Point2D&, float32, Point2D&, float32, Point2D&, float32);
    void AddSingleContour(const Point2D&, float32, const Point2D&, float32, const Point2D&, float32, const Point2D&, float32);
    float32 Interpolate(float32, float32, float32, float32);
    void DrawLineFromPoints(const Point2D&, const Point2D&);
    int square_count;

//Constructor / Destructor
public:
    AssignmentThree();
    virtual ~AssignmentThree();

//Methods
public:
    void DrawMesh();
    void MarchingSquares();
    void DrawScalarField();
    void DrawVectorField();
    virtual QWidget* createViewer();

//Attributes
public:
    ///File name of the scalar field
    string ScalarfieldFilename;

    ///The iso value for marching squares
    float IsoValue;

    ///File name of the vector field
    string VectorfieldFilename;

    ///Length of the arrows
    float ArrowScale;

    bool UseMidPointDecider;
    bool ShowScalarPoints;
    bool ShowMesh;

protected:
    GLGeometryViewer* viewer;
};

bool xcomp(Point2D p1, Point2D p2);
#endif
