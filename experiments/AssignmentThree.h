//---------------------------------------------------------------------------
#ifndef AssignmentThreeH
#define AssignmentThreeH
//---------------------------------------------------------------------------
#include "Experiment.h"
#include "LinearAlgebra.h"
#include "Field2.hpp"
#include "GLGeometryViewer.h"
//---------------------------------------------------------------------------


/// This is an assignment experiment.
///
class AssignmentThree : public Experiment
{
    GEOX_CLASS(AssignmentThree)

private:
    //Methods
    void AddContours(Point2D&, float32, Point2D&, float32, Point2D&, float32, Point2D&, float32);
    void AddSingleContour(const Point2D&, float32, const Point2D&, float32, const Point2D&, float32, const Point2D&, float32);
    float32 Interpolate(float32, float32, float32, float32);
    void DrawLineFromPoints(const Point2D&, const Point2D&);
    void MarchingSquaresHelper();
    void DrawScalarFieldHelper();
    void DrawMeshHelper();
    void LoadScalarField();

    //Attrs
    int square_count;
    Vector4f isocolor;

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
    void IsoContours();
    virtual QWidget* createViewer();

//Attributes
public:
    ///File name of the scalar field
    string ScalarfieldFilename;
    ///The loaded field
    ScalarField2 Field;

    ///The iso value for marching squares
    float IsoValue;

    ///File name of the vector field
    string VectorfieldFilename;

    ///Length of the arrows
    float ArrowScale;

    bool UseMidPointDecider;
    bool ShowScalarPoints;
    bool ShowMesh;

    int32 NumberOfIsoContours;

protected:
    GLGeometryViewer* viewer;
};

bool xcomp(Point2D p1, Point2D p2);
#endif
