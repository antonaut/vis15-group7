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

    //Attrs
    int square_count;
    Vector4f isocolor;

//Constructor / Destructor
public:
    AssignmentFour();
    virtual ~AssignmentFour();

//Methods
public:
    void DrawMesh();
    void DrawVectorField();
    virtual QWidget* createViewer();

//Attributes
public:
    ///The loaded field
    ScalarField2 Field;

    ///File name of the vector field
    string VectorfieldFilename;

    ///Length of the arrows
    float ArrowScale;
	
protected:
    GLGeometryViewer* viewer;
};

bool xcomp(Point2D p1, Point2D p2);
#endif
