//---------------------------------------------------------------------------
#include "stdafx.h"
//---------------------------------------------------------------------------
#include "AssignmentThree.h"
//---------------------------------------------------------------------------
#include "Properties.h"
#include "GeoXOutput.h"
//---------------------------------------------------------------------------

#include <limits>
#include <vector>


IMPLEMENT_GEOX_CLASS( AssignmentThree, 0)
{
    BEGIN_CLASS_INIT( AssignmentThree );

    ADD_NOARGS_METHOD(AssignmentThree::DrawMesh);

    ADD_SEPARATOR("Scalarfield")
    ADD_STRING_PROP(ScalarfieldFilename, 0)
    ADD_NOARGS_METHOD(AssignmentThree::DrawScalarField)

    ADD_SEPARATOR("Marching Squares")
    ADD_FLOAT32_PROP(IsoValue, 0)
    ADD_BOOLEAN_PROP(UseMidPointDecider, false)
    ADD_NOARGS_METHOD(AssignmentThree::MarchingSquares)

    ADD_SEPARATOR("IsoContours")
    ADD_INT32_PROP(NumberOfIsoContours, 0)
    ADD_NOARGS_METHOD(AssignmentThree::IsoContours)

    ADD_SEPARATOR("Drawing options")
    ADD_BOOLEAN_PROP(ShowScalarPoints, false)
    ADD_BOOLEAN_PROP(ShowMesh, false)

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
    ScalarfieldFilename = "./data/assignment03/SimpleGrid.am";

    IsoValue = 3.5;
    UseMidPointDecider = true;
    ShowScalarPoints = true;
    ShowMesh = true;

    NumberOfIsoContours = 5;
    isocolor = makeVector4f(1.0, 0.0, 0.0, 1); // Red

    VectorfieldFilename = "";
    ArrowScale = 0.1;
    square_count = 0;
}

AssignmentThree::~AssignmentThree() {}

void AssignmentThree::LoadScalarField() {
        //Load scalar field
    if (!Field.load(ScalarfieldFilename))
    {
        error("Error loading field file " + ScalarfieldFilename + "\n");
    }
}

void AssignmentThree::DrawMesh() {
    viewer->clear();

    LoadScalarField();

    viewer->refresh();
}

void AssignmentThree::DrawMeshHelper() {
    for(size_t i = 0; i < Field.dims()[0] - 1; ++i) {
        for(size_t j = 0; j < Field.dims()[1] - 1; ++j) {
            Point2D p1, p2, p3;
            p1.position = Field.nodePosition(i, j);
            p2.position = Field.nodePosition(i+1, j);
            p3.position = Field.nodePosition(i, j+1);

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

    for(size_t i = 0, j = Field.dims()[1] - 1; i < Field.dims()[0] - 1; ++i) {
        Point2D p1, p2;
        p1.position = Field.nodePosition(i, j);
        p2.position = Field.nodePosition(i+1, j);


        float32 x1 = p1.position[0];
        float32 y1 = p1.position[1];
        float32 x2 = p2.position[0];
        float32 y2 = p2.position[1];

        viewer->addLine(x1, y1, x2, y2);
    }

    for(size_t j = 0, i = Field.dims()[0] - 1; j < Field.dims()[1] - 1; ++j) {
        Point2D p1, p2;
        p1.position = Field.nodePosition(i, j);
        p2.position = Field.nodePosition(i, j+1);


        float32 x1 = p1.position[0];
        float32 y1 = p1.position[1];
        float32 x2 = p2.position[0];
        float32 y2 = p2.position[1];

        viewer->addLine(x1, y1, x2, y2);
    }
}

float32 AssignmentThree::Interpolate(float32 xy1, float32 v1, float32 xy2, float32 v2) {
    float32 t = (IsoValue - v1) / (v2 - v1);
    float32 res = (1-t)*xy1 + t*xy2;
    /*
    output << "Interpolation: (xy1, v1, xy2, v2) : " << "(" << xy1 << ", " << v1 << ", " << xy2 << ", " << v2 << ") : " << t << "\n";
    output << "Yields: " << res << "\n";
    */
    return res;
}


void AssignmentThree::AddSingleContour(const Point2D &p1, float32 v1, const Point2D &p2, float32 v2, const Point2D &p3, float32 v3, const Point2D &p4, float32 v4) {
    float32 x1 = p1.position[0];
    float32 y1 = p1.position[1];
    float32 x2 = p2.position[0];
    float32 y2 = p2.position[1];
    float32 x3 = p3.position[0];
    float32 y3 = p3.position[1];
    float32 x4 = p4.position[0];
    float32 y4 = p4.position[1];

    float32 line_x1 = Interpolate(x1, v1, x2, v2);
    float32 line_y1 = Interpolate(y1, v1, y2, v2);
    float32 line_x2 = Interpolate(x3, v3, x4, v4);
    float32 line_y2 = Interpolate(y3, v3, y4, v4);

    viewer->addLine(line_x1, line_y1, line_x2, line_y2, isocolor);
}

void AssignmentThree::AddContours(Point2D &p1, float32 v1, Point2D &p2, float32 v2, Point2D &p3, float32 v3, Point2D &p4, float32 v4) {
    /*
    output << "\n\nAdding contours (" << square_count << ") for " << p1 << " v:" << v1 << ", " << p2 << " v:" << v2 << ", " << p3 << " v:" << v3 << ", " << p4 << " v:" << v4 << "\n";

    */
    /*
     * No contours, return
     */

    if (v1 < IsoValue && v2 < IsoValue && v3 < IsoValue && v4 < IsoValue) {
        return;
    }

    if (v1 >= IsoValue && v2 >= IsoValue && v3 >= IsoValue && v4 >= IsoValue) {
        return;
    }

    /*
     * Single bottom left contour
     */

    if (v1 < IsoValue && v2 >= IsoValue && v3 >= IsoValue && v4 >= IsoValue) {
        AddSingleContour(p1, v1, p2, v2, p1, v1, p3, v3);
        return;
    }

    if (v1 > IsoValue && v2 < IsoValue && v3 < IsoValue && v4 < IsoValue) {
        AddSingleContour(p1, v1, p2, v2, p1, v1, p3, v3);
        return;
    }

    /*
     * Single bottom right contour
     */

    if (v2 < IsoValue && v1 >= IsoValue && v3 >= IsoValue && v4 >= IsoValue) {
        AddSingleContour(p2, v2, p1, v1, p2, v2, p4, v4);
        return;
    }

    if (v2 >= IsoValue && v1 < IsoValue && v3 < IsoValue && v4 < IsoValue) {
        AddSingleContour(p2, v2, p1, v1, p2, v2, p4, v4);
        return;
    }

    /*
     * Single top left contour
     */

    if (v3 < IsoValue && v1 >= IsoValue && v2 >= IsoValue && v4 >= IsoValue) {
        AddSingleContour(p3, v3, p1, v1, p3, v3, p4, v4);
        return;
    }

    if (v3 >= IsoValue && v1 < IsoValue && v2 < IsoValue && v4 < IsoValue) {
        AddSingleContour(p3, v3, p1, v1, p3, v3, p4, v4);
        return;
    }

    /*
     * Single top right contour
     */

    if (v4 < IsoValue && v1 >= IsoValue && v2 >= IsoValue && v3 >= IsoValue) {
        AddSingleContour(p4, v4, p2, v2, p4, v4, p3, v3);
        return;
    }

    if (v4 >= IsoValue && v1 < IsoValue && v2 < IsoValue && v3 < IsoValue) {
        AddSingleContour(p4, v4, p2, v2, p4, v4, p3, v3);
        return;
    }

    /*
     * Single vertical contour
     */

    if (v1 < IsoValue && v3 < IsoValue && v2 >= IsoValue && v4 >= IsoValue) {
        AddSingleContour(p1, v1, p2, v2, p3, v3, p4, v4);
        return;
    }

    if (v1 >= IsoValue && v3 >= IsoValue && v2 < IsoValue && v4 < IsoValue) {
        AddSingleContour(p1, v1, p2, v2, p3, v3, p4, v4);
        return;
    }

    /*
     * Single horizontal contour
     */

    if (v1 < IsoValue && v2 < IsoValue && v3 >= IsoValue && v4 >= IsoValue) {
        AddSingleContour(p1, v1, p3, v3, p2, v2, p4, v4);
        return;
    }

    if (v1 >= IsoValue && v2 >= IsoValue && v3 < IsoValue && v4 < IsoValue) {
        AddSingleContour(p1, v1, p3, v3, p2, v2, p4, v4);
        return;
    }

    /*
     * Double contours
     */

    float32 x1 = p1.position[0];
    float32 y1 = p1.position[1];

    float32 x2 = p2.position[0];
    float32 y2 = p2.position[1];
    
    float32 x3 = p3.position[0];
    float32 y3 = p3.position[1];
    
    float32 x4 = p4.position[0];
    float32 y4 = p4.position[1];


    float32 ip12 = Interpolate(x1, v1, x2, v2);
    float32 ip24 = Interpolate(y2, v2, y4, v4);
    float32 ip13 = Interpolate(y1, v1, y3, v3);
    float32 ip34 = Interpolate(x3, v3, x4, v4);

    Point2D pi12;
    pi12.position[0] = ip12;
    pi12.position[1] = y1;
    Point2D pi24;
    pi24.position[0] = x2;
    pi24.position[1] = ip24;
    Point2D pi13;
    pi13.position[0] = x1;
    pi13.position[1] = ip13;
    Point2D pi34;
    pi34.position[0] = ip34;
    pi34.position[1] = y3;
    /*
    output << "pi12: " << pi12 << " between " << p1 << " and " << p2 << ".\n";
    output << "pi24: " << pi24 << " between " << p2 << " and " << p4 << ".\n";
    output << "pi13: " << pi13 << " between " << p1 << " and " << p3 << ".\n";
    output << "pi34: " << pi34 << " between " << p3 << " and " << p4 << ".\n";
    */
   
    /*
     * Mid point decider.
     */

    if (UseMidPointDecider) {

        float32 mv = 0.25*(v1+v2+v3+v4);
        output << "mid value: " << mv << "\n";

        if (mv >= IsoValue) {
            DrawLineFromPoints(pi13, pi34);
            DrawLineFromPoints(pi12, pi24);
            return;
        }

        DrawLineFromPoints(pi13, pi12);
        DrawLineFromPoints(pi34, pi24);
        return;
    }

    /*
     * Asymptotic decider.
     */
    std::vector<Point2D> points;

    points.push_back(pi12);
    points.push_back(pi13);
    points.push_back(pi24);
    points.push_back(pi34);

    std::sort(points.begin(), points.end(), xcomp);

    p1 = points[0];
    p2 = points[1];
    p3 = points[2];
    p4 = points[3];

    DrawLineFromPoints(p1, p2);
    DrawLineFromPoints(p3, p4);
}

void AssignmentThree::DrawLineFromPoints(const Point2D& p1, const Point2D& p2) {
    /*
    output << "Drawing a line between " << p1 << " and " << p2 << ".\n";
    */
    viewer->addLine(p1.position[0], p1.position[1], p2.position[0], p2.position[1], isocolor);
}

void AssignmentThree::IsoContours() {
    viewer->clear();
    
    if (NumberOfIsoContours < 1) {
        error("Invalid number of IsoContours.");
        return;
    }

    LoadScalarField();

    float32 min_v = std::numeric_limits<float32>::max();
    float32 max_v = -std::numeric_limits<float32>::max();
    for(size_t i = 0; i < Field.dims()[0]; ++i) {
        for(size_t j = 0; j < Field.dims()[1]; ++j) {
            float32 v = Field.nodeScalar(i, j);

            min_v = min(min_v, v);
            max_v = max(max_v, v);
        }
    }

    output << "Min v: " << min_v << "\n";
    output << "Max v: " << max_v << "\n";
    
    // Red is highest iso value
    isocolor = makeVector4f(1.0, 0.0, 0.0, 1);

    float32 deltaColor = 1/((float32)(NumberOfIsoContours));
    float32 deltaIso = (max_v - min_v)*deltaColor;
    IsoValue = max_v;

    for (int i=NumberOfIsoContours; i>0;i--) {
        MarchingSquaresHelper();
        isocolor = makeVector4f(isocolor[0]-deltaColor, 0.0, isocolor[2]+deltaColor, 1.0);
        IsoValue -= deltaIso;
        viewer->refresh();
    }
}

void AssignmentThree::MarchingSquares() {
    viewer->clear();
    LoadScalarField();
    MarchingSquaresHelper();
    viewer->refresh();
}

void AssignmentThree::MarchingSquaresHelper() {
    if (ShowMesh) {
        DrawMeshHelper();
    }

    square_count = 0;

    for(size_t i = 0; i < Field.dims()[0] - 1; ++i) {
        for(size_t j = 0; j < Field.dims()[1] - 1; ++j) {
            Point2D p1, p2, p3, p4;
            p1.position = Field.nodePosition(i, j);
            p2.position = Field.nodePosition(i+1, j);
            p3.position = Field.nodePosition(i, j+1);
            p4.position = Field.nodePosition(i+1, j+1);

            float32 v1 = Field.nodeScalar(i, j);
            float32 v2 = Field.nodeScalar(i+1, j);
            float32 v3 = Field.nodeScalar(i, j+1);
            float32 v4 = Field.nodeScalar(i+1, j+1);

            ++square_count;
            AddContours(p1, v1, p2, v2, p3, v3, p4, v4);
        }
    }

    if (ShowScalarPoints) {
        DrawScalarFieldHelper();
    }
}

void AssignmentThree::DrawScalarField()
{
    viewer->clear();
    DrawScalarFieldHelper();
    viewer->refresh();
}

void AssignmentThree::DrawScalarFieldHelper() {
    
    //Draw a point for each grid vertex.
    for(size_t j=0; j<Field.dims()[1]; j++)
    {
        for(size_t i=0; i<Field.dims()[0]; i++)
        {
            const float32 val = Field.nodeScalar(i, j);
            const float32 c = val < IsoValue ? 0 : 1;

            Point2D p;
            p.position  = Field.nodePosition(i, j);
            p.size = 5;
            //Use a grayscale depending on the actual value
            p.color[0] = c; p.color[1] = c; p.color[2] = c;
            viewer->addPoint(p);
        }
    }
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

///Used for sorting on x-value.
bool xcomp(Point2D p1, Point2D p2) {
    return p1.position[0] < p2.position[0];
}

