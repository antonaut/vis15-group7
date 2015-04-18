//---------------------------------------------------------------------------
#include "stdafx.h"
//---------------------------------------------------------------------------
#include "AssignmentThree.h"
//---------------------------------------------------------------------------
#include "Properties.h"
#include "GLGeometryViewer.h"
#include "GeoXOutput.h"
//---------------------------------------------------------------------------

#include <limits>
#include "Field2.hpp"

IMPLEMENT_GEOX_CLASS( AssignmentThree, 0)
{
    BEGIN_CLASS_INIT( AssignmentThree );

    ADD_NOARGS_METHOD(AssignmentThree::DrawMesh);

    ADD_SEPARATOR("Scalarfield")
    ADD_STRING_PROP(ScalarfieldFilename, 0)
    ADD_NOARGS_METHOD(AssignmentThree::DrawScalarField)

    ADD_SEPARATOR("Marching Squares")
    ADD_FLOAT32_PROP(IsoValue, 0)
    ADD_NOARGS_METHOD(AssignmentThree::MarchingSquares)

    ADD_SEPARATOR("Vectorfield")
    ADD_STRING_PROP(VectorfieldFilename, 0)
    ADD_FLOAT32_PROP(ArrowScale, 0)
    ADD_NOARGS_METHOD(AssignmentThree::DrawVectorField)
}

QWidget* AssignmentThree::createViewer()
{
    viewer = new GLGeometryViewer();
    return viewer;
}

AssignmentThree::AssignmentThree()
{
    viewer = NULL;
    ScalarfieldFilename = "/home/simon/Dropbox/KTH/5an/visualization/ass3/SimpleGrid.am";
    
    IsoValue = 3;

    VectorfieldFilename = "";
    ArrowScale = 0.1;
}

AssignmentThree::~AssignmentThree() {}

void AssignmentThree::DrawMesh() {
    viewer->clear();

    //Load scalar field
    ScalarField2 field;
    if (!field.load(ScalarfieldFilename))
    {
        output << "Error loading field file " << ScalarfieldFilename << "\n";
        return;
    }

    for(size_t i = 0; i < field.dims()[0] - 1; ++i) {
        for(size_t j = 0; j < field.dims()[1] - 1; ++j) {
            Point2D p1, p2, p3;
            p1.position = field.nodePosition(i, j);
            p2.position = field.nodePosition(i+1, j);
            p3.position = field.nodePosition(i, j+1);

            float32 x1 = p1.position[0];
            float32 y1 = p1.position[1];
            float32 x2 = p2.position[0];
            float32 y2 = p2.position[1];
            float32 x3 = p3.position[0];
            float32 y3 = p3.position[1];

            viewer->addLine(x1, y1, x2, y2);
            viewer->addLine(x1, y1, x3, y3);
        }
    }

    for(size_t i = 0, j = field.dims()[1] - 1; i < field.dims()[0] - 1; ++i) {
        Point2D p1, p2;
        p1.position = field.nodePosition(i, j);
        p2.position = field.nodePosition(i+1, j);
    

        float32 x1 = p1.position[0];
        float32 y1 = p1.position[1];
        float32 x2 = p2.position[0];
        float32 y2 = p2.position[1];
     
        viewer->addLine(x1, y1, x2, y2);
    }

    for(size_t j = 0, i = field.dims()[0] - 1; j < field.dims()[1] - 1; ++j) {
        Point2D p1, p2;
        p1.position = field.nodePosition(i, j);
        p2.position = field.nodePosition(i, j+1);
    

        float32 x1 = p1.position[0];
        float32 y1 = p1.position[1];
        float32 x2 = p2.position[0];
        float32 y2 = p2.position[1];
     
        viewer->addLine(x1, y1, x2, y2);
    }
    
    viewer->refresh();
}

void AssignmentThree::DrawMesh() {
    viewer->clear();
    DrawMesh();
    
    viewer->refresh();
}

void AssignmentThree::DrawScalarField()
{
    viewer->clear();

    //Load scalar field
    ScalarField2 field;
    if (!field.load(ScalarfieldFilename))
    {
        output << "Error loading field file " << ScalarfieldFilename << "\n";
        return;
    }

    //Get the minimum/maximum value in that field
    float32 min = std::numeric_limits<float32>::max();
    float32 max = -std::numeric_limits<float32>::max();
    for(size_t j=0; j<field.dims()[1]; j++)
    {
        for(size_t i=0; i< field.dims()[0]; i++)
        {
            const float32 val = field.nodeScalar(i,j);
            min = val < min ? val : min;
            max = val > max ? val : max;
        }
    }

    //Draw a point for each grid vertex.
    for(size_t j=0; j<field.dims()[1]; j++)
    {
        for(size_t i=0; i<field.dims()[0]; i++)
        {
            const float32 val = field.nodeScalar(i, j);
            const float32 c = (val - min) / (max - min);

            Point2D p;
            p.position  = field.nodePosition(i, j);
            p.size = 5;
            //Use a grayscale depending on the actual value
            p.color[0] = c; p.color[1] = c; p.color[2] = c;
            viewer->addPoint(p);
        }
    }

    viewer->refresh();
}

void AssignmentThree::DrawVectorField()
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
