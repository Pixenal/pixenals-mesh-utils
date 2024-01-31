#include "Types.h"
#include "math.h"

Vec3 vec3MultiplyScalar(Vec3 a, float b) {
	Vec3 c = {a.x * b, a.y * b, a.z * b};
	return c;
}

void vec3DivideEqualScalar(Vec3 *pA, float b) {
	pA->x /= b;
	pA->y /= b;
	pA->z /= b;
}

Vec3 vec3SubtractScalar(Vec3 a, float b) {
	Vec3 c = {a.x - b, a.y - b, a.z - b};
	return c;
}

Vec3 vec3Subtract(Vec3 a, Vec3 b) {
	Vec3 c = {a.x - b.x, a.y - b.y, a.z - b.z};
	return c;
}

Vec3 vec3AddScalar(Vec3 a, float b) {
	Vec3 c = {a.x + b, a.y + b, a.z + b};
	return c;
}

Vec3 vec3Add(Vec3 a, Vec3 b) {
	Vec3 c = {a.x + b.x, a.y + b.y, a.z + b.z};
	return c;
}

void vec3AddEqual(Vec3 *pA, Vec3 b) {
	pA->x += b.x;
	pA->y += b.y;
	pA->z += b.z;
}

int32_t vec3GreaterThan(Vec3 a, Vec3 b) {
	return (a.x > b.x) && (a.y > b.y) && (a.z > b.z);
}

int32_t vec3LessThan(Vec3 a, Vec3 b) {
	return (a.x < b.x) && (a.y < b.y) && (a.z < b.z);
}

int32_t vec3ApproxEqual(Vec3 a, Vec3 b) {
	return 0;
}

Vec2 vec2Multiply(Vec2 a, Vec2 b) {
	Vec2 c;
	c.x = a.x * b.x;
	c.y = a.y * b.y;
	return c;
}

void vec2MultiplyEqual(Vec2 *a, Vec2 b) {
	a->x *= b.x;
	a->y *= b.y;
}

Vec2 vec2DivideScalar(Vec2 a, float b) {
	Vec2 c;
	c.x = a.x / b;
	c.y = a.y / b;
	return c;
}

void vec2DivideEqualScalar(Vec2 *a, float b) {
	a->x /= b;
	a->y /= b;
}

Vec2 vec2Subtract(Vec2 a, Vec2 b) {
	Vec2 c = {a.x - b.x, a.y - b.y};
	return c;
}

void vec2SubtractEqual(Vec2 *a, Vec2 b) {
	a->x -= b.x;
	a->y -= b.y;
}

Vec2 vec2SubtractScalar(Vec2 a, float b) {
	Vec2 c = {a.x - b, a.y - b};
	return c;
}

Vec2 vec2Add(Vec2 a, Vec2 b) {
	Vec2 c = {a.x + b.x, a.y + b.y};
	return c;
}

Vec2 vec2AddScalar(Vec2 a, float b) {
	Vec2 c = {a.x + b, a.y + b};
	return c;
}

void vec2AddEqual(Vec2 *a, Vec2 b) {
	a->x += b.x;
	a->y += b.y;
}

void vec2AddEqualScalar(Vec2 *pA, float b) {
	pA->x += b;
	pA->y += b;
}

void vec2MultiplyEqualScalar(Vec2 *a, float b) {
	a->x *= b;
	a->y *= b;
}

Vec2 vec2MultiplyScalar(Vec2 a, float b) {
	Vec2 c = {a.x * b, a.y * b};
	return c;
}

float vec2Dot(Vec2 a, Vec2 b) {
	return a.x * b.x + a.y * b.y;
}

Vec2 vec2Cross(Vec2 a) {
	Vec2 b = {a.y, -a.x};
	return b;
}

Vec2 vec2ModScalar(Vec2 a, float b) {
	Vec2 c = {fmod(a.x, b), fmod(a.y, b)};
	return c;
}

void vec2ModEqualScalar(Vec2 *a, float b) {
	a->x = fmod(a->x, b);
	a->y = fmod(a->y, b);
}

int32_t vec2GreaterThan(Vec2 a, Vec2 b) {
	return (a.x > b.x) && (a.y > b.y);
}

int32_t vec2LessThan(Vec2 a, Vec2 b) {
	return (a.x < b.x) && (a.y < b.y);
}

int32_t vec2LessThanEqualTo(Vec2 a, Vec2 b) {
	return (a.x <= b.x) && (a.y <= b.y);
}

float customFloor(float a) {
	int32_t aTrunc = a;
	aTrunc -= ((float)aTrunc != a) && (a < .0f);
	return aTrunc;
}

iVec2 vec2FloorAssign(Vec2 *pA) {
	iVec2 c;
	c.x = pA->x = customFloor(pA->x);
	c.y = pA->y = customFloor(pA->y);
	return c;
}


Vec3 barycentricToCartesian(Vec3 *pTri, Vec3 *pPoint) {
	Vec3 pointCartesian;
	pointCartesian.x = (pPoint->x * pTri[0].x) + (pPoint->y * pTri[1].x) + (pPoint->z * pTri[2].x);
	pointCartesian.y = (pPoint->x * pTri[0].y) + (pPoint->y * pTri[1].y) + (pPoint->z * pTri[2].y);
	pointCartesian.z = (pPoint->x * pTri[0].z) + (pPoint->y * pTri[1].z) + (pPoint->z * pTri[2].z);
	return pointCartesian;
}

Vec3 cartesianToBarycentric(Vec2 *pTri, Vec3 *pPoint) {
	Vec3 pointBc;
	double derta = .0;
	double dertau = .0;
	double dertav = .0;

	/*Perform cramers rule*/

	derta = (pTri[0].x * pTri[1].y) - (pTri[0].x * pTri[2].y) -
			(pTri[1].x * pTri[0].y) + (pTri[1].x * pTri[2].y) +
			(pTri[2].x * pTri[0].y) - (pTri[2].x * pTri[1].y);
	/*Get determinate of Au*/
	dertau = (pPoint->x * pTri[1].y) - (pPoint->x * pTri[2].y) -
			 (pTri[1].x * pPoint->y) + (pTri[1].x * pTri[2].y) +
			 (pTri[2].x * pPoint->y) - (pTri[2].x * pTri[1].y);
	/*Get determinate of Av*/
	dertav = (pTri[0].x * pPoint->y) - (pTri[0].x * pTri[2].y) -
			 (pPoint->x * pTri[0].y) + (pPoint->x * pTri[2].y) +
			 (pTri[2].x * pTri[0].y) - (pTri[2].x * pPoint->y);

	/*u = dert(Au) / dert(A)*/
	pointBc.x = dertau / derta;
	/*u = dert(Av) / dert(A)*/
	pointBc.y = dertav / derta;
	/*w can be derived from u and v*/
	pointBc.z = 1.0 - pointBc.x - pointBc.y;

	return pointBc;
}
