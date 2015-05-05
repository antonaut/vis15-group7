<<<<<<< HEAD
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
=======
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
	Vector2f Euler(Vector2f);
	Vector2f RK4(Vector2f);
	Vector2f FieldValue(Vector2f);
	Vector2f ExampleFieldValue(Vector2f);
	void DrawVectorFieldHelper();
	bool IsTooSlow(Vector2f);

    //Attrs

	///The method used to get vector field data.
	Vector2f(AssignmentFour::*VectorFieldAccessor)(Vector2f);

	///Should we use the vector field or the example field
	bool UseVectorField;

//Constructor / Destructor
public:
    AssignmentFour();
    virtual ~AssignmentFour();

//Methods
public:
    void LoadandRefreshVectorField();
	void EulerStreamline();
	void RungeKuttaStreamline();

	void GoodStepSize();

	void UseEllipseField();
    void DrawVectorField();
	void DrawStreamline(vector<Vector2f>);

	Vector2f Method(Vector2f);
	vector<Vector2f> Integrator(int, Vector2f(AssignmentFour::*)(Vector2f));
    virtual QWidget* createViewer();

//Attributes
public:
    ///The loaded field
    VectorField2 Field;

    ///File name of the vector field
    string VectorfieldFilename;

 
	float XStart, YStart, MaxDistance;
	//Euler
	float EulerStepSize;
	int EulerStep;
	//RK
	float RKStepSize;
	float MaxArchLength;
	int RKStep;

	bool DirectionFieldOnly;

	float ArrowScale;
	
protected:
    GLGeometryViewer* viewer;
};

#endif
>>>>>>> 649b2f2aa6a1009e3dbfe8184b8951cfa393642b
