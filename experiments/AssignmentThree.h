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

//Constructor / Destructor
public:
    AssignmentThree();
    virtual ~AssignmentThree();

//Methods
public:
    void CreateData_Random();
    void DrawScatterPlot();
    void DrawLine();
    void DrawCircle();
    void BouncingBall();
    virtual QWidget* createViewer();

private:
    void MoveBall();

//Attributes
public:
    vector< Vector2f > Data;

    // Used for the circle
    float Radius;
    Vector2f Center;
    int NumSamples;

    //Parameters for Bouncy Ball
    Vector2f Velocity;

protected:
    GLGeometryViewer* viewer;
};


#endif
