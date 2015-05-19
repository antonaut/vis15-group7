//---------------------------------------------------------------------------
#include "stdafx.h"
//---------------------------------------------------------------------------
#include "AssignmentSeven.h"
//---------------------------------------------------------------------------
#include "Properties.h"
#include "GeoXOutput.h"
//---------------------------------------------------------------------------

#include <limits>
#include <vector>
#include <algorithm>
#include <random>

#ifndef IS_POW_2
#define IS_POW_2(x) (((x) != 0) && (((x) & ((x) - 1)) == 0))
#endif

#ifdef _WIN32
#define isnan(x) (_isnan(x) != 0)
#endif

IMPLEMENT_GEOX_CLASS(AssignmentSeven, 0)
{
	BEGIN_CLASS_INIT(AssignmentSeven);

	ADD_SEPARATOR("Vectorfield")
	ADD_STRING_PROP(VectorfieldFilename, 0)
	ADD_NOARGS_METHOD(AssignmentSeven::LoadVectorFieldAndRefresh)

	ADD_SEPARATOR("Scalarfield")
	ADD_STRING_PROP(ScalarfieldFilename, 0)
	ADD_NOARGS_METHOD(AssignmentSeven::LoadScalarFieldAndRefresh)

	ADD_SEPARATOR("Critical Points")
	ADD_NOARGS_METHOD(AssignmentSeven::ShowExtremePoints)

	ADD_SEPARATOR("Texture");
	ADD_STRING_PROP(TextureFilename, 0)
	ADD_VECTOR2I_PROP(TextureResolution, 0);
	ADD_BOOLEAN_PROP(GrayScale, 0);
	ADD_BOOLEAN_PROP(AddStreamLines, 0);
	ADD_NOARGS_METHOD(AssignmentSeven::DrawTexture)
	ADD_NOARGS_METHOD(AssignmentSeven::EnhanceTexture)
	
	ADD_SEPARATOR("LIC Options")
	ADD_CARD32_PROP(KernelSize, 0)
	ADD_INT32_PROP(Seed, 0)
	ADD_BOOLEAN_PROP(ColoredTexture, 0)
	ADD_NOARGS_METHOD(AssignmentSeven::FastLIC)

	ADD_NOARGS_METHOD(AssignmentSeven::ScalarSkeleton)

}

QWidget* AssignmentSeven::createViewer()
{
	viewer = new GLGeometryViewer();
	return viewer;
}

AssignmentSeven::AssignmentSeven()
{
	viewer = NULL;

	VectorfieldFilename = "../data/assignment07/Data/Vector/ANoise2CT4.am";

	VectorfieldFilename = "/home/simon/Git/vis15-group7/data/assignment07/Data/Vector/ANoise2CT4.am";
	ScalarfieldFilename = "/home/simon/Git/vis15-group7/data/assignment07/Data/Scalar/NoisyHill.am";

	TextureFilename = "../data/assignment06/";
	TextureResolution = makeVector2ui(64, 64);
	GrayScale = false;
	AddStreamLines = false;

	RKStepSize = 0.3;
	RKStep = 30;

	VectorFieldAccessor = &AssignmentSeven::FieldValue;

	arcLength = 0.0f;
	ArrowScale = 0.2f;
	KernelSize = 32;

	Seed = 1;

	ColoredTexture = false;

	texture = getRandomField(makeVector2f(-5, -5), makeVector2f(5, 5), TextureResolution, GrayScale);
}

AssignmentSeven::~AssignmentSeven() {}

void AssignmentSeven::DrawVectorFieldHelper() {
	//Draw vector directions (constant length)
	for (float32 x = Field.boundMin()[0]; x <= Field.boundMax()[0]; x += 0.2)
	{
		for (float32 y = Field.boundMin()[1]; y <= Field.boundMax()[1]; y += 0.2)
		{
			Vector2f vec = Field.sample(x, y);

			float32 x2 = x + ArrowScale*vec[0];
			float32 y2 = y + ArrowScale*vec[1];

			viewer->addLine(x, y, x2, y2);

			Point2D current = Point2D(x2,y2);
			current.color = makeVector4f(1, 1, 1, 1);
			viewer->addPoint(current);
			
		}
	}
}

void AssignmentSeven::ShowExtremePoints() {

	const float32 ZERO_THRESHOLD = 1E-8;

	vector< Vector2f > sources;
	vector< Vector2f > repfocs;
	vector< Vector2f > saddles;
	vector< Vector2f > centers;
	vector< Vector2f > sinks;
	vector< Vector2f > attfocs;
	//DrawVectorFieldHelper();
	float32 dx = abs(Field.boundMin()[0] - Field.boundMax()[0]) / 10.0;
	float32 dy = abs(Field.boundMin()[1] - Field.boundMax()[1]) / 10.0;

	for (float32 j = Field.boundMin()[1]; j <= Field.boundMax()[1] - dy; j += (dy)) {
		for (float32 i = Field.boundMin()[0]; i <= Field.boundMax()[0]-dx; i += (dx)) {
			//output << "Processing node(" << i << ", " << j << ").\n";
			Vector2f xy = Field.sample(i, j);
			Vector2f xxy = Field.sample(i + dx, j);
			Vector2f xyy = Field.sample(i, j + dy);
			Vector2f xxyy = Field.sample(i + dx, j + dy);
			bool extremePointInSquare = false;
			for (auto d = 0; d < 2; ++d) {
				float32 v1 = xy[d];
				float32 v2 = xxy[d];
				float32 v3 = xyy[d];
				float32 v4 = xxyy[d];
				bool b1 = signbit(v1);
				bool b2 = signbit(v2);
				bool b3 = signbit(v3);
				bool b4 = signbit(v4);

				//output << "[ " << (b1 ? "-" : "+") << " " << (b2 ? "-" : "+") << "\n  " << (b3 ? "-" : "+") << " " << (b4 ? "-" : "+") << " ]\n";
				int minusSigns = 0;
				if (b1) ++minusSigns;
				if (b2) ++minusSigns;
				if (b3) ++minusSigns;
				if (b4) ++minusSigns;

				output << "nsigns: " << minusSigns << "\n";
				if (minusSigns == 1 ||
					minusSigns == 3 /*||
					(b1 == b3 && b2 == b4 && b1 != b2) ||
					(b1 != b3 && b2 != b4 && b1 == b2)*/
					) {
					extremePointInSquare = true;
					break;
				}
			}

			if (extremePointInSquare) {
				// Let's find that extreme point.
				// Start from point xy and head towards middelkth.
				Vector2f point = makeVector2f(i, j);
				int32 iter = 0;
				Matrix2f jac;
				Matrix2f jacInv;
				//output << "Start: " << point << "\n";
				float32 norm = std::numeric_limits<float32>::max();

				Vector2f dp;
				while (norm > ZERO_THRESHOLD && iter < 10) {
					Point2D current = Point2D(point[0], point[1]);
					current.color = makeVector4f((10.0-iter)/10.0, iter/10.0, 0.3, 1);
					current.size = 10;
					//viewer->addPoint(current);
					//viewer->refresh();

					jac = Field.sampleJacobian(point);
					jacInv = invertMatrix(jac);
					dp = jacInv*Field.sample(point);
					point = point - dp;
					if (!Field.insideBounds(point)){
						iter = 20;
						break;
					}
					norm = point.getSqrNorm();
					iter++;
					//output << "point: " << point << ", dp: " << dp << ", det(J): " << jac.getDeterminant() << ", norm(dp): " << norm << "\n";
				}

				if (iter <= 10) {
					output << "Found extreme point: " << point << " after " << iter << " iterations.\n";

					Point2D pt(point[0], point[1]);

					// Time to classify!
					Matrix2f eigenVectors;
					Vector2f eigenRealValues;
					Vector2f eigenImagValues;
					jac.solveEigenProblem(eigenRealValues, eigenImagValues, eigenVectors);

					float32 r1 = eigenRealValues[0];
					float32 r2 = eigenRealValues[1];
					float32 i1 = eigenImagValues[0];
					float32 i2 = eigenImagValues[1];

					if ((r1 < 0 && r2 > 0) || (r1 > 0 && r2 < 0) && equals(i1, 0) && equals(i2, 0)) {
						// Saddle
						pt.color = makeVector4f(1, 1, 0, 1);
						saddles.push_back(point);
					}
					else if (r1 > 0 && r2 > 0 && equals(i1, 0) && equals(i2, 0)) {
						// Repelling node
						pt.color = makeVector4f(1, 0, 0, 1);
					}
					else if (equals(r1, r2) && r1 > 0 && equals(i1, -i2) && !equals(i1, 0)) {
						// Repelling focus
						pt.color = makeVector4f(1, 0.5f, 0, 1);
					}
					else if (equals(r1, r2) && equals(r1, 0) && equals(i1, -i2) && !equals(i1, 0)) {
						// Center
						pt.color = makeVector4f(0, 1, 0, 1);
					}
					else if (equals(r1, r2) && r1 < 0 && equals(i1, -i2) && !equals(i1, 0)) {
						// Attracting focus
						pt.color = makeVector4f(0.5f, 1, 0, 1);
					}
					else if (r1 < 0 && r2 < 0 && equals(i1, 0) && equals(i2, 0)) {
						// Attracting node
						pt.color = makeVector4f(0, 0, 1, 1);
					}
					else {
						output << "r1, r2, i1, i2: " << r1 << ", " << r2 << ", " << i1 << ", " << i2 << "\n";
						pt.color = makeVector4f(1, 1, 1, 1);
					}

					pt.size = 6;
					viewer->addPoint(pt);
				}
			}
		}
	}

	viewer->refresh();
}

bool AssignmentSeven::equals(float32 a, float32 b) const {
	return abs(a - b) < 0.001;
}


void AssignmentSeven::LoadVectorField() {

	if (!Field.load(VectorfieldFilename))
	{
		error("Error loading field file " + VectorfieldFilename + "\n");
	}

	VectorFieldAccessor = &AssignmentSeven::FieldValue;
}


void AssignmentSeven::LoadVectorFieldAndRefresh() {
	viewer->clear();

	LoadVectorField();

	viewer->refresh();
}

void AssignmentSeven::LoadScalarField() {
	if (!SField.load(ScalarfieldFilename)) {
		error("Error loading file" + ScalarfieldFilename + "\n");
	}
}

void AssignmentSeven::LoadScalarFieldAndRefresh() {
	viewer->clear();
	LoadScalarField();
	viewer->refresh();
}


void AssignmentSeven::DrawTexture() {
	viewer->clear();

	LoadVectorField();


	texture = getRandomField(Field.boundMin(), Field.boundMax(), TextureResolution, GrayScale);
	viewer->setTextureGray(texture.getData());

	viewer->refresh();
}

void AssignmentSeven::EnhanceTexture() {
	viewer->clear();
	texture = enhanceContrast(texture);
	viewer->setTextureGray(texture.getData());
	viewer->refresh();
}

void AssignmentSeven::FastLIC() {
	viewer->clear();

	mt.seed(Seed);
	LoadVectorField();

	const Vector2ui &textureResolution = TextureResolution;

	//const VectorField2 vectorField = getEllipseField(makeVector2f(-5, -5), makeVector2f(5, 5), makeVector2ui(16, 16));
	const VectorField2 &vectorField = Field;

	const Vector2ui &dims = vectorField.dims();
	const Vector2f &boundMin = vectorField.boundMin();
	const Vector2f &boundMax = vectorField.boundMax();

	ScalarField2 randomField = getRandomField(boundMin, boundMax, textureResolution, GrayScale);
	ScalarField2 smearedField(randomField);
	smearedField.setZero();

	vector< vector<card32> > timesRendered(textureResolution[0], vector<card32>(textureResolution[1], 0));
	vector< vector<float32> > valueSum(textureResolution[0], vector<float32>(textureResolution[1], 0));

	size_t numStreamLines = 0;

	const float32 dx = (boundMax[0] - boundMin[0]) / textureResolution[0];
	const float32 dy = (boundMax[1] - boundMin[1]) / textureResolution[1];
	for (card32 x = 0; x < textureResolution[0]; ++x) {
		for (card32 y = 0; y < textureResolution[1]; ++y) {
			if (timesRendered[x][y] > 0) {
				continue;
			}

			float32 xstart = boundMin[0] + dx * x + dx/2;
			float32 ystart = boundMin[1] + dy * y + dy/2;

			if (vectorField.sample(xstart, ystart).getSqrNorm() == 0) {
				timesRendered[x][y] += 1;
				valueSum[x][y] += 0;
				continue;
			}

			vector<Vector2f> streamLine = Integrator(128, &AssignmentSeven::RK4, xstart, ystart);

			numStreamLines += 1;

			if (AddStreamLines) {
				drawStreamline(streamLine, makeVector4f(1, 0, 0, 1));
			}

			if (streamLine.empty()) {
				timesRendered[x][y] += 1;
				valueSum[x][y] += 0;
				continue;
			}

			vector<Vector2ui> pixels = streamLineToPixels(randomField, streamLine);

			assert(std::find(pixels.begin(), pixels.end(), makeVector2ui(x, y)) != pixels.end());

			assert(!pixels.empty());

			vector<float32> smearValues = smear(randomField, pixels);

			assert(smearValues.size() == pixels.size());

			for (size_t i = 0; i < smearValues.size(); ++i) {
				const Vector2ui &p = pixels[i];
				card32 px = p[0];
				card32 py = p[1];

				assert(px >= 0 && py >= 0 && px < textureResolution[0] && py < textureResolution[1]);

				timesRendered[p[0]][p[1]] += 1;
				valueSum[p[0]][p[1]] += smearValues[i];
			}
		}
	}

	for (card32 x = 0; x < textureResolution[0]; ++x) {
		for (card32 y = 0; y < textureResolution[1]; ++y) {
			float32 v = valueSum[x][y] / timesRendered[x][y];
			if (isnan(v)) {
				v = randomField.nodeScalar(x, y);
				output << "nan (" << x << ", " << y << "): vs=" << valueSum[x][y] << ", " << "tr=" << timesRendered[x][y] << "\n";
			}
			smearedField.setNodeScalar(x, y, v);
		}
	}

	output << "Number of streamLines: " << numStreamLines << "\n";

	texture = smearedField;

	viewer->setTextureGray(smearedField.getData());

	viewer->refresh();
}

void AssignmentSeven::drawStreamline(vector<Vector2f> path, const Vector4f &color)
{
	size_t arraySize = path.size();
	if (arraySize < 2) {
		return;
	}
	Vector2f p1 = path[0];
	Vector2f p2;
	for(size_t i = 1; i < arraySize; ++i) {
		p2 = path[i];
		viewer->addLine(p1[0], p1[1], p2[0], p2[1], color);
		viewer->addPoint(p1);
		p1 = p2;
	}
	viewer->addPoint(path[arraySize - 1]);
}

ScalarField2 AssignmentSeven::getRandomField(const Vector2f &boundMin, const Vector2f &boundMax,
										   const Vector2ui &dims, bool grayscale) {
	ScalarField2 field = ScalarField2();
	field.init(boundMin, boundMax, dims);

	for (card32 i = 0; i < dims[0]; ++i) {
		for (card32 j = 0; j < dims[1]; ++j) {
			float32 value = randomFloat(0, 1);
			if (!grayscale) {
				value = round(value);
			}

			field.setNodeScalar(i, j, value);
		}
	}

	return field;
}

VectorField2 AssignmentSeven::getEllipseField(const Vector2f &boundMin, const Vector2f &boundMax, const Vector2ui &dims) const {
	VectorField2 field = VectorField2();
	field.init(boundMin, boundMax, dims);

	Vector2f c = (boundMax + boundMin) / 2;

	for (card32 x = 0; x < dims[0]; ++x) {
		for (card32 y = 0; y < dims[0]; ++y) {
			Vector2f p = field.nodePosition(x, y);
			Vector2f v = rotate(makeVector2f(0, 1), atan2(p[1] - c[1], p[0] - c[0]));
			field.setNode(x, y, v);
		}
	}

	return field;
}

vector<Vector2ui> AssignmentSeven::streamLineToPixels(const ScalarField2 &field, const vector<Vector2f> &streamLine) {
	vector<Vector2ui> pixels;

	for (size_t i = 0; i < streamLine.size() - 1; ++i) {
		const Vector2f &v0 = streamLine[i];
		const Vector2f &v1 = streamLine[i+1];

		vector<Vector2ui> line_pixels = lineToPixels(field, v0, v1);
		std::copy(line_pixels.begin(), line_pixels.end(), back_inserter(pixels));
	}

	return pixels;
}


vector<Vector2ui> AssignmentSeven::lineToPixels(const ScalarField2 &field, const Vector2f &v0, const Vector2f &v1) const {
	vector<Vector2ui> points;

	const Vector2f &bn = field.boundMin();
	const Vector2f &bx = field.boundMax();
	float32 dx = (bx[0] - bn[0]) / field.dims()[0];
	float32 dy = (bx[1] - bn[1]) / field.dims()[1];

	size_t steps = dx < dy ? abs(v1[0] - v0[0]) / dx : abs(v1[1] - v0[1]) / dy;
	float32 d = min(dx, dy);

	Vector2f direction = v1 - v0;
	direction.normalize();
	direction *= d;

	points.push_back(field.closestNode(v0));

	for (size_t i = 0; i <= steps; ++i) {
		Vector2ui node = field.closestNode(v0 + direction*i);
		if (node != points.back()) {
			points.push_back(node);
		}
	}

	Vector2ui last_node = field.closestNode(v1);
	if (points.back() != last_node) {
		points.push_back(last_node);
	}

	return points;
}

vector<float32> AssignmentSeven::smear(const ScalarField2 &field, const vector<Vector2ui> &pixels) {
	vector<float32> smears;

	vector<float32> vals;
	for (size_t i = 0; i < pixels.size(); ++i) {
		vals.push_back(field.nodeScalar(pixels[i][0], pixels[i][1]));
	}

	card32 kernelSize = std::min<card32>(KernelSize, (card32) pixels.size());
	card32 numVal = 0;
	float32 sumVal = 0.0f;
	for (size_t i = 0; i < kernelSize / 2; ++i) {
		sumVal += vals[i];
		numVal += 1;
	}

	smears.push_back(sumVal / numVal);

	for (size_t i = 1; i < pixels.size(); ++i) {

		if (i >= kernelSize / 2) {
			sumVal -= vals[i - kernelSize / 2];
			numVal -= 1;
		}

		if (i < pixels.size() - kernelSize / 2) {
			sumVal += vals[i + kernelSize / 2];
			numVal += 1;
		}

		smears.push_back(sumVal / numVal);
	}

	return smears;
}

Vector2f AssignmentSeven::rotate(const Vector2f &f, float32 angle) const {
	float32 x = f[0]*cos(angle) - f[1]*sin(angle);
	float32 y = f[0]*sin(angle) + f[1]*cos(angle);

	return makeVector2f(x, y);
}

ScalarField2 AssignmentSeven::enhanceContrast(ScalarField2 image) {

	Vector2ui dims = makeVector2ui(image.dims()[0], image.dims()[1]);
	float32 sum = 0.0f;
	float32 P = 0.0f;
	int nonBlack = 0;

	for (card32 i = 0; i < dims[0]; ++i){
		for (card32 j = 0; j < dims[1]; ++j) {
			float32 val = image.nodeScalar(i, j);
			if (val != 0.0f) {
				nonBlack++;
				sum += val;
				P += pow(val, 2);
			}
		}
	}
	float32 mean = sum / (dims[0] * dims[1]);
	float32 dev = sqrt((P-nonBlack*pow(mean, 2))/(nonBlack-1));
	float32 f = 0.1 / dev;
	const float32 MAX_STRETCH_FACTOR = 20;
	f = f > MAX_STRETCH_FACTOR ? f : MAX_STRETCH_FACTOR;
	
	ScalarField2 enhanced = ScalarField2();
	enhanced.init(makeVector2f(image.boundMin()[0], image.boundMin()[1]), makeVector2f(image.boundMax()[0], image.boundMax()[1]), dims);

	for (card32 i = 0; i < dims[0]; ++i) {
		for (card32 j = 0; j < dims[1]; ++j) {
			float32 val = image.nodeScalar(i, j);
			enhanced.setNodeScalar(i, j, 0.5 + f*(val - mean));
		}
	}

	return enhanced;
}

bool AssignmentSeven::IsTooSlow(Vector2f vec) {
	float threshold = 1E-4;
	return vec.getSqrNorm() < threshold;
}

void AssignmentSeven::DrawScalarField(const ScalarField2 &field) {
	float32 minValue = field.nodeScalar(0, 0);
	float32 maxValue = field.nodeScalar(0, 0);

	ScalarField2 renderField = ScalarField2();
	renderField.init(field.boundMin(), field.boundMax(), field.dims());

	for (card32 row = 0; row < field.dims()[0]; ++row) {
		for (card32 col = 0; col < field.dims()[1]; ++col) {
			minValue = min(minValue, field.nodeScalar(row, col));
			maxValue = max(maxValue, field.nodeScalar(row, col));
		}
	}

	for (card32 row = 0; row < field.dims()[0]; ++row) {
		for (card32 col = 0; col < field.dims()[1]; ++col) {
			float32 data = (field.nodeScalar(row, col) - minValue) / (maxValue - minValue);
			renderField.setNodeScalar(row, col, data);
		}
	}

	viewer->setTextureGray(renderField.getData());
}

void AssignmentSeven::ScalarSkeleton() {
	viewer->clear();

	LoadScalarField();

	DrawScalarField(SField);

	Field.init(SField.boundMin(), SField.boundMax(), SField.dims());
	const card32 nx = Field.dims()[0];
	const card32 ny = Field.dims()[1];
	const float dx = (Field.boundMax()[0] - Field.boundMin()[0]) / nx;
	const float dy = (Field.boundMax()[1] - Field.boundMin()[1]) / ny;
	for (card32 xi = 0; xi < nx; ++xi) {
		float32 x = SField.boundMin()[0] + xi * dx;
		for (card32 yi = 0; yi < ny; ++yi) {
			float32 y = SField.boundMin()[1] + yi * dy;
			Vector2f grad = SField.sampleGradient(x, y);
			//output << "grad: " << grad[0] << ", " << grad[1] << "\n";
			Field.setNode(xi, yi, grad);
		}
	}

	viewer->refresh();

	ShowExtremePoints();
}

/**
Takes two parameters. The number of steps and a method of integration (i.e. Euler or RK4).
Returns a vector with a stream line path.
*/
vector<Vector2f>
AssignmentSeven::Integrator(
		int numberOfSteps,
		Vector2f(AssignmentSeven::*Method)(Vector2f, bool),
		float32 xstart, float32 ystart
)
{
	Vector2f xi;
	vector<Vector2f> path;
	arcLength = 0.0f;

	xi[0] = xstart;
	xi[1] = ystart;

	vector<Vector2f> bw;

	for (int i = 0; i < numberOfSteps; i++)
	{
		Vector2f xp = (this->*Method)(xi, true);
		if (IsTooSlow(xp)) {
			output << "Stopped early after " << i << " steps. (Going too slow)\n";
			break;
		}

		if (!Field.insideBounds(xp)) {
			bw.push_back(xp);
			break;
		}

		bw.push_back(xp);
		arcLength += (xp - xi).getSqrNorm();

		/*
		if (arcLength > MaxDistance) {
			output << "Stopped early after " << i << " steps. (Maximum distance)\n";
			return path;
		}*/

		xi = xp;
	}

	bw.push_back(xi);
	std::copy(bw.rbegin(), bw.rend(), back_inserter(path));

	xi[0] = xstart;
	xi[1] = ystart;

	path.push_back(xi);

	for (int i = 0; i < numberOfSteps; i++)
	{
		Vector2f xp = (this->*Method)(xi, false);
		if (IsTooSlow(xp)) {
			output << "Stopped early after " << i << " steps. (Going too slow)\n";
			return path;
		}

		if (!Field.insideBounds(xp)) {
			path.push_back(xp);
			return path;
		}

		path.push_back(xp);
		arcLength += (xp - xi).getSqrNorm();

		/*
		if (arcLength > MaxDistance) {
			output << "Stopped early after " << i << " steps. (Maximum distance)\n";
			return path;
		}*/

		xi = xp;
	}

	path.push_back(xi);

	return path;
}

Vector2f AssignmentSeven::FieldValue(Vector2f xi, bool integrateBackwards) {
	StaticVector<float, 2U> vec = Field.sample(xi[0], xi[1]);
	Vector2f v = makeVector2f(vec[0], vec[1]);

	return integrateBackwards ? -v : v;
}

Vector2f AssignmentSeven::RK4(Vector2f xi, bool integrateBackwards)
{
	Vector2f v1, v2, v3, v4;
	v1 = (this->*VectorFieldAccessor)(xi, integrateBackwards);
	v2 = (this->*VectorFieldAccessor)(xi + (v1 * (RKStepSize / 2.0f)), integrateBackwards);
	v3 = (this->*VectorFieldAccessor)(xi + (v2 * (RKStepSize / 2.0f)), integrateBackwards);
	v4 = (this->*VectorFieldAccessor)(xi + (v3 * RKStepSize), integrateBackwards);

	return xi + (v1 + v2 * 2.0f + v3 * 2.0f + v4) * RKStepSize / 6.0f;
}

float32 AssignmentSeven::randomFloat(float32 a, float32 b) {
	std::uniform_real_distribution<float32> dist(a, b);
	return dist(mt);
}
