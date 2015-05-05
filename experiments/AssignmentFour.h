//---------------------------------------------------------------------------
#ifndef AssignmentFourH
#define AssignmentFourH
//---------------------------------------------------------------------------
#include "Experiment.h"
#include "LinearAlgebra.h"
#include "GLGeometryViewer3D.h"
//---------------------------------------------------------------------------
#include <ostream>
#include <vector>

using namespace std;

struct Cube {
    float32 v[8];

    friend ostream& operator<< (ostream &out, const Cube &c) {
        return out << "(" << c.v[0] << ", " << c.v[1] << ", " << c.v[2] << ", " 
            << c.v[3] << ", " << c.v[4] << ", " << c.v[5] << ", " << c.v[6] << ", " << c.v[7] << ")";
    }

    friend istream& operator>> (istream &in, Cube &c) {
        return in >> c.v[0] >> c.v[1] >> c.v[2] >> c.v[3] >> c.v[4] >> c.v[5] >> c.v[6] >> c.v[7];
    }
};

///
/// This is an example experiment.
///
/// The code is meant to demonstrate how
///  to use the GeoX framework
///
class AssignmentFour : public Experiment
{
   GEOX_CLASS(AssignmentFour)
    private:
        vector<Cube> cubes;        

        GLGeometryViewer3D* viewer;

        void setCorners(const Cube &, const Vector3f &, Vector3f[8]);
        void loadCubes();
        void drawContours(const Cube &, const Vector3f &);
        void drawCube(const Cube &, const Vector3f &);

    public:
        string Filename;
        float32 IsoValue;

        AssignmentFour();
        ~AssignmentFour();

        void DrawCubes();
        
        void showPerVertexLighting();
        void applyMarchingCubes();
        virtual QWidget *createViewer();
};


#endif
