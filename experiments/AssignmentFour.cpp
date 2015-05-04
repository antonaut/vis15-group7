//---------------------------------------------------------------------------
#include "stdafx.h"
//---------------------------------------------------------------------------
#include "AssignmentFour.h"
//---------------------------------------------------------------------------
#include "Properties.h"
#include "GeoXOutput.h"
//---------------------------------------------------------------------------

#include <limits>
#include <vector>


IMPLEMENT_GEOX_CLASS( AssignmentFour, 0)
{
    BEGIN_CLASS_INIT( AssignmentFour );

    ADD_NOARGS_METHOD(AssignmentFour::DrawMesh);

    ADD_SEPARATOR("Vectorfield")
    ADD_STRING_PROP(VectorfieldFilename, 0)
    ADD_FLOAT32_PROP(ArrowScale, 0)
    ADD_NOARGS_METHOD(AssignmentFour::DrawVectorField)
}

QWidget* AssignmentFour::createViewer()
{
    viewer = new GLGeometryViewer();
    return viewer;
}

AssignmentFour::AssignmentFour()
{
    viewer = NULL;
   
	VectorfieldFilename = "";
    ArrowScale = 0.1;
    square_count = 0;
}

AssignmentFour::~AssignmentFour() {}

void AssignmentFour::LoadVectorField() {
        //Load scalar field
    if (!Field.load(VectorfieldFilename))
    {
        error("Error loading field file " + VectorfieldFilename + "\n");
    }
}

void AssignmentFour::DrawMesh() {
    viewer->clear();

    LoadVectorField();

    viewer->refresh();
}

void AssignmentFour::DrawVectorField()
{
    viewer->clear();

    //Load the vector field
    VectorField2 field;
    if (!field.load(VectorfieldFilename))
    {
        output << "Error loading field file " << VectorfieldFilename << "\n";
        return;
    }

    //Draw vector directions (constant length)
    for(float32 x=field.boundMin()[0]; x<=field.boundMax()[0]; x+=0.2)
    {
        for(float32 y=field.boundMin()[1]; y<=field.boundMax()[1]; y+=0.2)
        {
            Vector2f vec = field.sample(x,y);
            vec.normalize();

            viewer->addLine(x, y, x + ArrowScale*vec[0], y + ArrowScale*vec[1]);
        }
    }

    viewer->refresh();
}


namespace
{
    ///Returns the next power-of-two
    int32 NextPOT(int32 n)
    {
        n--;
        n |= n >> 1;   // Divide by 2^k for consecutive doublings of k up to 32,
        n |= n >> 2;   // and then or the results.
        n |= n >> 4;
        n |= n >> 8;
        n |= n >> 16;
        n++;           // The result is a number of 1 bits equal to the number
                       // of bits in the original number, plus 1. That's the
                       // next highest power of 2.
        return n;
    }
}

///Used for sorting on x-value.
bool xcomp(Point2D p1, Point2D p2) {
    return p1.position[0] < p2.position[0];
}

