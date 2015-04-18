//---------------------------------------------------------------------------
#include "stdafx.h"
//---------------------------------------------------------------------------
#include "AssignmentThree.h"
//---------------------------------------------------------------------------
#include "Properties.h"
#include "GLGeometryViewer.h"
#include "GeoXOutput.h"

//---------------------------------------------------------------------------

IMPLEMENT_GEOX_CLASS( AssignmentThree, 0)
{
    BEGIN_CLASS_INIT( AssignmentThree );
    ADD_NOARGS_METHOD(AssignmentThree::CreateData_Random);
    ADD_NOARGS_METHOD(AssignmentThree::DrawScatterPlot);
    ADD_NOARGS_METHOD(AssignmentThree::DrawLine);


    ADD_SEPARATOR("Circle");
    ADD_FLOAT32_PROP(Radius, 0);
    ADD_VECTOR2F_PROP(Center, 0);
    ADD_INT32_PROP(NumSamples, 0);
    ADD_NOARGS_METHOD(AssignmentThree::DrawCircle);

    ADD_SEPARATOR("BOUNCY BALL");
    ADD_VECTOR2F_PROP(Velocity, 0);
    ADD_NOARGS_METHOD(AssignmentThree::BouncingBall);
}

QWidget* AssignmentThree::createViewer()
{
    viewer = new GLGeometryViewer();
    return viewer;
}

AssignmentThree::AssignmentThree()
{
    viewer = NULL;
    Radius = 1.0f;
    Center = makeVector2f(0, 0);
    NumSamples = 100;
}

AssignmentThree::~AssignmentThree() {}


void AssignmentThree::CreateData_Random()
{
    Data.clear();
    Data.resize(1000);
    for(int i=0;i<(int)Data.size();i++)
    {
        Data[i][0] = (float)rand() / RAND_MAX;
        Data[i][1] = (float)rand() / RAND_MAX;
    }
}


void AssignmentThree::DrawScatterPlot()
{
    for(int i=0;i<(int)Data.size();i++)
    {
        Point2D NewPoint(Data[i][0], Data[i][1]);
        NewPoint.color = makeVector4f(Data[i][0], Data[i][1], 0.5, 1);
        viewer->addPoint(NewPoint);
    }

    //Axes
    Point2D Origin(0, 0);
    Origin.color = makeVector4f(1,1,1,1);
    Origin.size = 10;
    const int idOrigin = viewer->addPoint(Origin);

    Point2D XOff(1.1, 0);
    XOff.color = makeVector4f(1,0,0,1);
    XOff.size = 10;
    const int idXOff = viewer->addPoint(XOff);

    Point2D YOff(0, 1.1);
    YOff.color = makeVector4f(0,1,0,1);
    YOff.size = 10;
    const int idYOff = viewer->addPoint(YOff);

    //X-Axis
    Line Axis;
    Axis.vertices[0] = idOrigin;
    Axis.vertices[1] = idXOff;
    Axis.color = makeVector4f(1,1,1,1);
    Axis.thickness = 3;
    viewer->addLine(Axis);

    //Y-Axis
    Axis.vertices[1] = idYOff;
    viewer->addLine(Axis);

    // display changes
    viewer->refresh();
}

void AssignmentThree::DrawLine()
{
    viewer->addLine(0, 0, 10, 10);

    viewer->refresh();
}


void AssignmentThree::DrawCircle()
{
    for(int i=1; i <= NumSamples; ++i) {
        const float PI = 3.14159265359;
        const float x1 = Center[0] +
            Radius * cos(2 * PI * float(i)/float(NumSamples));
        const float y1 = Center[1] +
            Radius * sin(2 * PI * float(i)/float(NumSamples));
        const float x2 = Center[0] +
            Radius * cos(2 * PI * float(i+1)/float(NumSamples));
        const float y2 = Center[1] +
            Radius * sin(2 * PI * float(i+1)/float(NumSamples));

        viewer->addLine(x1, y1, x2, y2);
    }
    viewer->refresh();
}

void AssignmentThree::BouncingBall()
{
    for (int i=0; i<100; ++i)
    {
        MoveBall();
        DrawCircle();
    }
}

void AssignmentThree::MoveBall()
{
    Center += Velocity;
    if (Center.getSqrNorm() > 8.0) {
        Velocity = -Velocity;
    }
}
