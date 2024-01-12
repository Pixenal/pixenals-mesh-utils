#include "Types.h"

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
int32_t vec2GreaterThan(Vec2 a, Vec2 b) {
	return (a.x > b.x) && (a.y > b.y);
}
int32_t vec2LessThan(Vec2 a, Vec2 b) {
	return (a.x < b.x) && (a.y < b.y);
}

Vec3 barycentricToCartesian(Vec3 *triVert0, Vec3 *triVert1, Vec3 *triVert2,
                            Vec3 *point) {
	Vec3 pointCartesian;
	pointCartesian.x = (point->x * triVert0->x) + (point->y * triVert1->x) + (point->z * triVert2->x);
	pointCartesian.y = (point->x * triVert0->y) + (point->y * triVert1->y) + (point->z * triVert2->y);
	pointCartesian.z = (point->x * triVert0->z) + (point->y * triVert1->z) + (point->z * triVert2->z);
	return pointCartesian;
}

Vec3 cartesianToBarycentric(Vec3 *triVert0, Vec3 *triVert1, Vec3 *triVert2,
                            Vec3 *point) {
	Vec3 pointBc;
	double derta = .0;
	double dertau = .0;
	double dertav = .0;

	/*Perform cramers rule*/

	derta = (triVert0->x * triVert1->y) - (triVert0->x * triVert2->y) -
			(triVert1->x * triVert0->y) + (triVert1->x * triVert2->y) +
			(triVert2->x * triVert0->y) - (triVert2->x * triVert1->y);
	/*Get determinate of Au*/
	dertau = (point->x * triVert1->y) - (point->x * triVert2->y) -
			 (triVert1->x * point->y) + (triVert1->x * triVert2->y) +
			 (triVert2->x * point->y) - (triVert2->x * triVert1->y);
	/*Get determinate of Av*/
	dertav = (triVert0->x * point->y) - (triVert0->x * triVert2->y) -
			 (point->x * triVert0->y) + (point->x * triVert2->y) +
			 (triVert2->x * triVert0->y) - (triVert2->x * point->y);

	/*u = dert(Au) / dert(A)*/
	pointBc.x = dertau / derta;
	/*u = dert(Av) / dert(A)*/
	pointBc.y = dertav / derta;
	/*w can be derived from u and v*/
	pointBc.z = 1.0 - pointBc.x - pointBc.y;

	return pointBc;
}
