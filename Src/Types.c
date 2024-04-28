#include <math.h>
#include <float.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include <Types.h>

#define INDEX_ATTRIB(t, pD, i, v, c) ((t (*)[v])pD->pData)[i][c]

#define BLEND_REPLACE(t, pD, iD, pA, iA, pB, iB, v, c)\
	INDEX_ATTRIB(t, pD, iD, v, c) = INDEX_ATTRIB(t, pB, iB, v, c)

#define BLEND_MULTIPLY(t, pDest, iDest, pA, iA, pB, iB, v, c)\
	INDEX_ATTRIB(t,pD,iD,v,c) = INDEX_ATTRIB(t,pA,iA,v,c) * INDEX_ATTRIB(t,pB,iB,v,c)

#define BLEND_DIVIDE(t, pDest, iDest, pA, iA, pB, iB, v, c)\
	INDEX_ATTRIB(t,pD,iD,v,c) = INDEX_ATTRIB(t,pB,iB,v,c) != (t)0 ?\
		INDEX_ATTRIB(t,pA,iA,v,c) / INDEX_ATTRIB(t,pB,iB,v,c) : (t)0

#define BLEND_ADD(t, pDest, iDest, pA, iA, pB, iB, v, c)\
	INDEX_ATTRIB(t,pD,iD,v,c) = INDEX_ATTRIB(t,pA,iA,v,c) + INDEX_ATTRIB(t,pB,iB,v,c)

#define BLEND_SUBTRACT(t, pDest, iDest, pA, iA, pB, iB, v, c)\
	INDEX_ATTRIB(t,pD,iD,v,c) = INDEX_ATTRIB(t,pA,iA,v,c) - INDEX_ATTRIB(t,pB,iB,v,c)

#define BLEND_ADD_SUB(t, pDest, iDest, pA, iA, pB, iB, v, c)\
	INDEX_ATTRIB(t,pD,iD,v,c) = INDEX_ATTRIB(t,pA,iA,v,c) +\
		INDEX_ATTRIB(t,pB,iB,v,c) - ((t)1 - INDEX_ATTRIB(t,pB,iB,v,c))

#define BLEND_LIGHTEN(t, pDest, iDest, pA, iA, pB, iB, v, c)\
	INDEX_ATTRIB(t,pD,iD,v,c) = INDEX_ATTRIB(t,pA,iA,v,c) > INDEX_ATTRIB(t,pB,iB,v,c) ?\
		INDEX_ATTRIB(t,pA,iA,v,c) : INDEX_ATTRIB(t,pB,iB,v,c)

#define BLEND_DARKEN(t, pDest, iDest, pA, iA, pB, iB, v, c)\
	INDEX_ATTRIB(t,pD,iD,v,c) = INDEX_ATTRIB(t,pA,iA,v,c) < INDEX_ATTRIB(t,pB,iB,v,c) ?\
		INDEX_ATTRIB(t,pA,iA,v,c) : INDEX_ATTRIB(t,pB,iB,v,c)

#define BLEND_OVERLAY(t, pDest, iDest, pA, iA, pB, iB, v, c)\
	INDEX_ATTRIB(t,pD,iD,v,c) = INDEX_ATTRIB(t,pA,iA,v,c) > .5 ?\
		2.0 * INDEX_ATTRIB(t,pA,iA,v,c) * INDEX_ATTRIB(t,pB,iB,v,c) :\
		1.0 - 2.0 * (1.0 - INDEX_ATTRIB(t,pA,iA,v,c)) * (1.0 - INDEX_ATTRIB(t,pB,iB,v,c))

#define BLEND_SOFT_LIGHT(t, pDest, iDest, pA, iA, pB, iB, v, c)\
	INDEX_ATTRIB(t,pD,iD,v,c) = INDEX_ATTRIB(t,pB,iB,v,c) < .5 ?\
\
		2.0 * INDEX_ATTRIB(t,pA,iA,v,c) * INDEX_ATTRIB(t,pB,iB,v,c) +\
		INDEX_ATTRIB(t,pA,iA,v,c) * INDEX_ATTRIB(t,pA,iA,v,c) *\
		(1.0 - 2.0 * INDEX_ATTRIB(t,pB,iB,v,c)) :\
\
		2.0 * INDEX_ATTRIB(t,pA,iA,v,c) * (1.0 - INDEX_ATTRIB(t,pB,iB,v,c)) +\
		sqrt(INDEX_ATTRIB(t,pA,iA,v,c)) * (2.0 * INDEX_ATTRIB(t,pB,iB,v,c) - 1.0)

#define BLEND_COLOR_DODGE(t, pDest, iDest, pA, iA, pB, iB, v, c) {\
	t d =INDEX_ATTRIB(t,pA,iA,v,c) / (1.0 - INDEX_ATTRIB(t,pB,iB,v,c));\
	t e = d < 1.0 ? d : 1.0;\
	INDEX_ATTRIB(t,pD,iD,v,c) = INDEX_ATTRIB(t,pB,iB,v,c) == 1.0 ?\
		INDEX_ATTRIB(t,pB,iB,v,c) : e;\
}

#define INTERPOLATE_SCALAR(t, pD, iD, pS, iA, iB, iC, bc)\
	INDEX_ATTRIB(t, pD, iD, 1, 0) = INDEX_ATTRIB(t, pS, iA, 1, 0) * bc.x;\
 	INDEX_ATTRIB(t, pD, iD, 1, 0) += INDEX_ATTRIB(t, pS, iB, 1, 0) * bc.y;\
	INDEX_ATTRIB(t, pD, iD, 1, 0) += INDEX_ATTRIB(t, pS, iC, 1, 0) * bc.z;\
	INDEX_ATTRIB(t, pD, iD, 1, 0) /= bc.x + bc.y + bc.z;

#define INTERPOLATE_V2(t, pD, iD, pS, iA, iB, iC, bc) {\
	INDEX_ATTRIB(t,pD,iD,2,0) = INDEX_ATTRIB(t,pS,iA,2,0) * bc.x;\
	INDEX_ATTRIB(t,pD,iD,2,1) = INDEX_ATTRIB(t,pS,iA,2,1) * bc.x;\
	INDEX_ATTRIB(t,pD,iD,2,0) += INDEX_ATTRIB(t,pS,iB,2,0) * bc.y;\
	INDEX_ATTRIB(t,pD,iD,2,1) += INDEX_ATTRIB(t,pS,iB,2,1) * bc.y;\
	INDEX_ATTRIB(t,pD,iD,2,0) += INDEX_ATTRIB(t,pS,iC,2,0) * bc.z;\
	INDEX_ATTRIB(t,pD,iD,2,1) += INDEX_ATTRIB(t,pS,iC,2,1) * bc.z;\
	float sum = bc.x + bc.y + bc.z;\
	INDEX_ATTRIB(t,pD,iD,2,0) /= sum;\
	INDEX_ATTRIB(t,pD,iD,2,1) /= sum;\
}

#define INTERPOLATE_V3(t, pD, iD, pS, iA, iB, iC, bc) {\
	INDEX_ATTRIB(t,pD,iD,3,0) = INDEX_ATTRIB(t,pS,iA,3,0) * bc.x;\
	INDEX_ATTRIB(t,pD,iD,3,1) = INDEX_ATTRIB(t,pS,iA,3,1) * bc.x;\
	INDEX_ATTRIB(t,pD,iD,3,2) = INDEX_ATTRIB(t,pS,iA,3,2) * bc.x;\
	INDEX_ATTRIB(t,pD,iD,3,0) += INDEX_ATTRIB(t,pS,iB,3,0) * bc.y;\
	INDEX_ATTRIB(t,pD,iD,3,1) += INDEX_ATTRIB(t,pS,iB,3,1) * bc.y;\
	INDEX_ATTRIB(t,pD,iD,3,2) += INDEX_ATTRIB(t,pS,iB,3,2) * bc.y;\
	INDEX_ATTRIB(t,pD,iD,3,0) += INDEX_ATTRIB(t,pS,iC,3,0) * bc.z;\
	INDEX_ATTRIB(t,pD,iD,3,1) += INDEX_ATTRIB(t,pS,iC,3,1) * bc.z;\
	INDEX_ATTRIB(t,pD,iD,3,2) += INDEX_ATTRIB(t,pS,iC,3,2) * bc.z;\
	float sum = bc.x + bc.y + bc.z;\
	INDEX_ATTRIB(t,pD,iD,3,0) /= sum;\
	INDEX_ATTRIB(t,pD,iD,3,1) /= sum;\
	INDEX_ATTRIB(t,pD,iD,3,2) /= sum;\
}

#define INTERPOLATE_V4(t, pD, iD, pS, iA, iB, iC, bc) {\
	INDEX_ATTRIB(t,pD,iD,4,0) = INDEX_ATTRIB(t,pS,iA,4,0) * bc.x;\
	INDEX_ATTRIB(t,pD,iD,4,1) = INDEX_ATTRIB(t,pS,iA,4,1) * bc.x;\
	INDEX_ATTRIB(t,pD,iD,4,2) = INDEX_ATTRIB(t,pS,iA,4,2) * bc.x;\
	INDEX_ATTRIB(t,pD,iD,4,3) = INDEX_ATTRIB(t,pS,iA,4,3) * bc.x;\
	INDEX_ATTRIB(t,pD,iD,4,0) += INDEX_ATTRIB(t,pS,iB,4,0) * bc.y;\
	INDEX_ATTRIB(t,pD,iD,4,1) += INDEX_ATTRIB(t,pS,iB,4,1) * bc.y;\
	INDEX_ATTRIB(t,pD,iD,4,2) += INDEX_ATTRIB(t,pS,iB,4,2) * bc.y;\
	INDEX_ATTRIB(t,pD,iD,4,3) += INDEX_ATTRIB(t,pS,iB,4,3) * bc.y;\
	INDEX_ATTRIB(t,pD,iD,4,0) += INDEX_ATTRIB(t,pS,iC,4,0) * bc.z;\
	INDEX_ATTRIB(t,pD,iD,4,1) += INDEX_ATTRIB(t,pS,iC,4,1) * bc.z;\
	INDEX_ATTRIB(t,pD,iD,4,2) += INDEX_ATTRIB(t,pS,iC,4,2) * bc.z;\
	INDEX_ATTRIB(t,pD,iD,4,3) += INDEX_ATTRIB(t,pS,iC,4,3) * bc.z;\
	float sum = bc.x + bc.y + bc.z;\
	INDEX_ATTRIB(t,pD,iD,4,0) /= sum;\
	INDEX_ATTRIB(t,pD,iD,4,1) /= sum;\
	INDEX_ATTRIB(t,pD,iD,4,2) /= sum;\
	INDEX_ATTRIB(t,pD,iD,4,3) /= sum;\
}

Vec364 vec364MultiplyScalar(Vec364 a, double b) {
	Vec364 c = {a.x * b, a.y * b, a.z * b};
	return c;
}
void vec364DivideEqualScalar(Vec364 *pA, double b) {
	pA->x /= b;
	pA->y /= b;
	pA->z /= b;
}
Vec364 vec364DivideScalar(Vec364 a, double b) {
	Vec364 c = {a.x / b, a.y / b, a.z / b};
	return c;
}
void vec364AddEqual(Vec364 *pA, Vec364 b) {
	pA->x += b.x;
	pA->y += b.y;
	pA->z += b.z;
}

Vec264 vec264MultiplyScalar(Vec264 a, double b) {
	Vec264 c = {a.x * b, a.y * b};
	return c;
}
void vec264DivideEqualScalar(Vec264 *pA, double b) {
	pA->x /= b;
	pA->y /= b;
}
Vec264 vec264DivideScalar(Vec264 a, double b) {
	Vec264 c = {a.x / b, a.y / b};
	return c;
}
void vec264AddEqual(Vec264 *pA, Vec264 b) {
	pA->x += b.x;
	pA->y += b.y;
}

iVec38 veci38MultiplyScalar(iVec38 a, float b) {
	iVec38 c = {a.x * b, a.y * b, a.z * b};
	return c;
}
void vec38DivideEqualScalar(iVec38 *pA, float b) {
	pA->x /= b;
	pA->y /= b;
	pA->z /= b;
}
iVec38 veci38DivideScalar(iVec38 a, float b) {
	iVec38 c = {a.x / b, a.y / b, a.z / b};
	return c;
}
void veci38AddEqual(iVec38 *pA, iVec38 b) {
	pA->x += b.x;
	pA->y += b.y;
	pA->z += b.z;
}

iVec316 veci316MultiplyScalar(iVec316 a, float b) {
	iVec316 c = {a.x * b, a.y * b, a.z * b};
	return c;
}
void veci316DivideEqualScalar(iVec316 *pA, float b) {
	pA->x /= b;
	pA->y /= b;
	pA->z /= b;
}
iVec316 veci316DivideScalar(iVec316 a, float b) {
	iVec316 c = {a.x / b, a.y / b, a.z / b};
	return c;
}
void veci316AddEqual(iVec316 *pA, iVec316 b) {
	pA->x += b.x;
	pA->y += b.y;
	pA->z += b.z;
}

iVec3 veci3MultiplyScalar(iVec3 a, float b) {
	iVec3 c = {a.x * b, a.y * b, a.z * b};
	return c;
}
void veci3DivideEqualScalar(iVec3 *pA, float b) {
	pA->x /= b;
	pA->y /= b;
	pA->z /= b;
}
iVec3 veci3DivideScalar(iVec3 a, float b) {
	iVec3 c = {a.x / b, a.y / b, a.z / b};
	return c;
}
void veci3AddEqual(iVec3 *pA, iVec3 b) {
	pA->x += b.x;
	pA->y += b.y;
	pA->z += b.z;
}

iVec364 veci364MultiplyScalar(iVec364 a, float b) {
	iVec364 c = {a.x * b, a.y * b, a.z * b};
	return c;
}
void veci364DivideEqualScalar(iVec364 *pA, float b) {
	pA->x /= b;
	pA->y /= b;
	pA->z /= b;
}
iVec364 veci364DivideScalar(iVec364 a, float b) {
	iVec364 c = {a.x / b, a.y / b, a.z / b};
	return c;
}
void veci364AddEqual(iVec364 *pA, iVec364 b) {
	pA->x += b.x;
	pA->y += b.y;
	pA->z += b.z;
}

iVec28 veci28MultiplyScalar(iVec28 a, float b) {
	iVec28 c = {a.x * b, a.y * b};
	return c;
}
void veci28DivideEqualScalar(iVec28 *pA, float b) {
	pA->x /= b;
	pA->y /= b;
}
iVec28 veci28DivideScalar(iVec28 a, float b) {
	iVec28 c = {a.x / b, a.y / b};
	return c;
}
void veci28AddEqual(iVec28 *pA, iVec28 b) {
	pA->x += b.x;
	pA->y += b.y;
}

iVec216 veci216MultiplyScalar(iVec216 a, float b) {
	iVec216 c = {a.x * b, a.y * b};
	return c;
}
void veci216DivideEqualScalar(iVec216 *pA, float b) {
	pA->x /= b;
	pA->y /= b;
}
iVec216 veci216DivideScalar(iVec216 a, float b) {
	iVec216 c = {a.x / b, a.y / b};
	return c;
}
void veci216AddEqual(iVec216 *pA, iVec216 b) {
	pA->x += b.x;
	pA->y += b.y;
}

iVec2 veci2MultiplyScalar(iVec2 a, float b) {
	iVec2 c = {a.x * b, a.y * b};
	return c;
}
void veci2DivideEqualScalar(iVec2 *pA, float b) {
	pA->x /= b;
	pA->y /= b;
}
iVec2 veci2DivideScalar(iVec2 a, float b) {
	iVec2 c = {a.x / b, a.y / b};
	return c;
}
void veci2AddEqual(iVec2 *pA, iVec2 b) {
	pA->x += b.x;
	pA->y += b.y;
}

iVec264 veci264MultiplyScalar(iVec264 a, float b) {
	iVec264 c = {a.x * b, a.y * b};
	return c;
}
void veci264DivideEqualScalar(iVec264 *pA, float b) {
	pA->x /= b;
	pA->y /= b;
}
iVec264 veci264DivideScalar(iVec264 a, float b) {
	iVec264 c = {a.x / b, a.y / b};
	return c;
}
void veci264AddEqual(iVec264 *pA, iVec264 b) {
	pA->x += b.x;
	pA->y += b.y;
}

Vec3 vec3MultiplyScalar(Vec3 a, float b) {
	Vec3 c = {a.x * b, a.y * b, a.z * b};
	return c;
}

void vec3DivideEqualScalar(Vec3 *pA, float b) {
	pA->x /= b;
	pA->y /= b;
	pA->z /= b;
}

Vec3 vec3DivideScalar(Vec3 a, float b) {
	Vec3 c = {a.x / b, a.y / b, a.z / b};
	return c;
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

int32_t vec3AproxEqual(Vec3 a, Vec3 b) {
	Vec3 bLow = _(b V3SUBS FLOAT_EQUAL_MARGIN);
	Vec3 bHigh = _(b V3ADDS FLOAT_EQUAL_MARGIN);
	return _(a V3LESS bHigh) && _(a V3GREAT bLow);
}

Vec3 vec3Lerp(Vec3 a, Vec3 b, float alpha) {
	float alphaInverse = 1.0f - alpha;
	Vec3 c;
	c.x = a.x * alphaInverse + b.x * alpha;
	c.y = a.y * alphaInverse + b.y * alpha;
	c.z = a.z * alphaInverse + b.z * alpha;
	return c;
}

Vec3 vec3UnitFromPoints(Vec3 a, Vec3 b) {
	Vec3 dir = _(b V3SUB a); //direction
	float magnitude = sqrt(dir.x * dir.x + dir.y * dir.y);
	return _(dir V3DIVS magnitude);
}

Vec3 vec3Cross(Vec3 a, Vec3 b) {
	Vec3 c = {
		.x = a.y * b.z - a.z * b.y,
		.y = a.z * b.x - a.x * b.z,
		.z = a.x * b.y - a.y * b.x
	};
	return c;
}

Vec3 vec3MultiplyMat3x3(Vec3 a, Mat3x3 *pB) {
	Vec3 c;
	c.x = a.x * pB->d[0][0] + a.y * pB->d[1][0] + a.z * pB->d[2][0];
	c.y = a.x * pB->d[0][1] + a.y * pB->d[1][1] + a.z * pB->d[2][1];
	c.z = a.x * pB->d[0][2] + a.y * pB->d[1][2] + a.z * pB->d[2][2];
	return c;
}

Vec3 vec3Normalize(Vec3 a) {
	float magnitude = sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
	return _(a V3DIVS magnitude);
}

float vec3Dot(Vec3 a, Vec3 b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec2 vec2Abs(Vec2 a) {
	if (a.x < .0f) {
		a.x *= -1.0f;
	}
	if (a.y < .0f) {
		a.y *= -1.0f;
	}
	return a;
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

int32_t vec2GreaterThanEqualTo(Vec2 a, Vec2 b) {
	return (a.x >= b.x) && (a.y >= b.y);
}

int32_t vec2GreaterThan(Vec2 a, Vec2 b) {
	return (a.x > b.x) && (a.y > b.y);
}

int32_t vec2GreaterThanScalar(Vec2 a, float b) {
	return (a.x > b) && (a.y > b);
}

int32_t vec2LessThan(Vec2 a, Vec2 b) {
	return (a.x < b.x) && (a.y < b.y);
}

int32_t vec2LessThanScalar(Vec2 a, float b) {
	return (a.x < b) && (a.y < b);
}

int32_t vec2LessThanEqualTo(Vec2 a, Vec2 b) {
	return (a.x <= b.x) && (a.y <= b.y);
}

int32_t vec2NotEqual(Vec2 a, Vec2 b) {
	return a.x != b.x || a.y != b.y;
}

int32_t vec2Equal(Vec2 a, Vec2 b) {
	return a.x == b.x && a.y == b.y;
}

int32_t vec2AproxEqual(Vec2 a, Vec2 b) {
	Vec2 bLow = _(b V2SUBS FLOAT_EQUAL_MARGIN);
	Vec2 bHigh = _(b V2ADDS FLOAT_EQUAL_MARGIN);
	return _(a V2LESS bHigh) && _(a V2GREAT bLow);
}

int32_t vec2WindingCompare(Vec2 a, Vec2 b, Vec2 centre, int32_t fallBack) {
	Vec2 aDiff = _(a V2SUB centre);
	Vec2 bDiff = _(b V2SUB centre);
	float cross = aDiff.x * bDiff.y - bDiff.x * aDiff.y;
	if (cross != .0f) {
		return cross > .0f;
	}
	if (fallBack) {
	float aDist = aDiff.x * aDiff.x + aDiff.y * aDiff.y;
	float bDist = bDiff.x * bDiff.x + bDiff.y * bDiff.y;
	return bDist > aDist;
	}
	else {
		return 2;
	}
}

Mat2x2 mat2x2Adjugate(Mat2x2 a) {
	Mat2x2 c;
	c.d[0][0] = a.d[1][1];
	c.d[0][1] = -a.d[0][1];
	c.d[1][0] = -a.d[1][0];
	c.d[1][1] = a.d[0][0];
	return c;
}

float mat2x2Determinate(Mat2x2 a) {
	return a.d[0][0] * a.d[1][1] - a.d[0][1] * a.d[1][0];
}

void mat2x2MultiplyEqualScalar(Mat2x2 *pA, float b) {
	pA->d[0][0] *= b;
	pA->d[0][1] *= b;
	pA->d[1][0] *= b;
	pA->d[1][1] *= b;
}

Mat2x2 mat2x2Invert(Mat2x2 a) {
	float determinate = mat2x2Determinate(a);
	Mat2x2 inverse = mat2x2Adjugate(a);
	mat2x2MultiplyEqualScalar(&inverse, 1.0f / determinate);
	return inverse;
}

float mat3x3Determinate(Mat3x3 *pA) {
	float aDet = pA->d[1][1] * pA->d[2][2] - pA->d[2][1] * pA->d[1][2];
	float bDet = pA->d[0][1] * pA->d[2][2] - pA->d[2][1] * pA->d[0][2];
	float cDet = pA->d[0][1] * pA->d[1][2] - pA->d[1][1] * pA->d[0][2];
	return pA->d[0][0] * aDet - pA->d[1][0] * bDet + pA->d[2][0] * cDet;
}

Mat3x3 mat3x3Adjugate(Mat3x3 *pA) {
	Mat3x3 c;
	c.d[0][0] = pA->d[1][1] * pA->d[2][2] - pA->d[2][1] * pA->d[1][2];
	c.d[0][1] = pA->d[0][1] * pA->d[2][2] - pA->d[2][1] * pA->d[0][2];
	c.d[0][2] = pA->d[0][1] * pA->d[1][2] - pA->d[1][1] * pA->d[0][2];
	c.d[1][0] = pA->d[1][0] * pA->d[2][2] - pA->d[2][0] * pA->d[1][2];
	c.d[1][1] = pA->d[0][0] * pA->d[2][2] - pA->d[2][0] * pA->d[0][2];
	c.d[1][2] = pA->d[0][0] * pA->d[1][2] - pA->d[1][0] * pA->d[0][2];
	c.d[2][0] = pA->d[1][0] * pA->d[2][1] - pA->d[2][0] * pA->d[1][1];
	c.d[2][1] = pA->d[0][0] * pA->d[2][1] - pA->d[2][0] * pA->d[0][1];
	c.d[2][2] = pA->d[0][0] * pA->d[1][1] - pA->d[1][0] * pA->d[0][1];
	c.d[1][0] *= -1.0f;
	c.d[0][1] *= -1.0f;
	c.d[2][1] *= -1.0f;
	c.d[1][2] *= -1.0f;
	return c;
}

void mat3x3MultiplyEqualScalar(Mat3x3 *pA, float b) {
	for (int32_t i = 0; i < 3; ++i) {
		for (int32_t j = 0; j < 3; ++j) {
			pA->d[i][j] *= b;
		}
	}
}

Mat3x3 mat3x3FromVec3(Vec3 a, Vec3 b, Vec3 c) {
	Mat3x3 mat;
	mat.d[0][0] = a.x;
	mat.d[0][1] = a.y;
	mat.d[0][2] = a.z;
	mat.d[1][0] = b.x;
	mat.d[1][1] = b.y;
	mat.d[1][2] = b.z;
	mat.d[2][0] = c.x;
	mat.d[2][1] = c.y;
	mat.d[2][2] = c.z;
	return mat;
}

Mat3x3 mat3x3Invert(Mat3x3 *pA) {
	float determinate = mat3x3Determinate(pA);
	Mat3x3 inverse = mat3x3Adjugate(pA);
	mat3x3MultiplyEqualScalar(&inverse, 1.0f / determinate);
	return inverse;
}

Mat2x3 mat2x2MultiplyMat2x3(Mat2x2 a, Mat2x3 b) {
	Mat2x3 c;
	c.d[0][0] = a.d[0][0] * b.d[0][0] + a.d[0][1] * b.d[1][0];
	c.d[0][1] = a.d[0][0] * b.d[0][1] + a.d[0][1] * b.d[1][1];
	c.d[0][2] = a.d[0][0] * b.d[0][2] + a.d[0][1] * b.d[1][2];
	c.d[1][0] = a.d[1][0] * b.d[0][0] + a.d[1][1] * b.d[1][0];
	c.d[1][1] = a.d[1][0] * b.d[0][1] + a.d[1][1] * b.d[1][1];
	c.d[1][2] = a.d[1][0] * b.d[0][2] + a.d[1][1] * b.d[1][2];
	return c;
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

int32_t checkFaceIsInBounds(Vec2 min, Vec2 max, FaceInfo face, Mesh *pMesh) {
		/*
		int32_t isInside = 0;
		for (int32_t j = 0; j < face.size; ++j) {
			int32_t vertIndex = pMesh->pLoops[face.start + j];
			Vec2 *vert = (Vec2 *)(pMesh->pVerts + vertIndex);
			if (_(*vert V2GREATEQL min) && _(*vert V2LESSEQL max)) {
				isInside = 1;
				break;
			}
		}
		return isInside;
		*/
	Vec2 faceMin, faceMax;
	faceMin.x = faceMin.y = FLT_MAX;
	faceMax.x = faceMax.y = 0;
	for (int32_t i = 0; i < face.size; ++i) {
		int32_t vertIndex = pMesh->mesh.pLoops[face.start + i];
		Vec3 *pVert = pMesh->pVerts + vertIndex;
		if (pVert->x < faceMin.x) {
			faceMin.x = pVert->x;
		}
		if (pVert->y < faceMin.y) {
			faceMin.y = pVert->y;
		}
		if (pVert->x > faceMax.x) {
			faceMax.x = pVert->x;
		}
		if (pVert->y > faceMax.y) {
			faceMax.y = pVert->y;
		}
	}
	iVec2 inside;
	inside.x = (faceMin.x > min.x && faceMin.x < max.x) ||
	           (faceMax.x > min.x && faceMax.x < max.x) ||
			   (faceMin.x < min.x && faceMax.x > max.x);
	inside.y = (faceMin.y > min.y && faceMin.y < max.y) ||
	           (faceMax.y > min.y && faceMax.y < max.y) ||
			   (faceMin.y < min.y && faceMax.y > max.y);
	return inside.x && inside.y;
}

uint32_t ruvmFnvHash(uint8_t *value, int32_t valueSize, uint32_t size) {
	uint32_t hash = 2166136261;
	for (int32_t i = 0; i < valueSize; ++i) {
		hash ^= value[i];
		hash *= 16777619;
	}
	hash %= size;
	return hash;
}

void getFaceBounds(FaceBounds *pBounds, RuvmVec2 *pUvs, FaceInfo faceInfo) {
	pBounds->fMin.x = pBounds->fMin.y = FLT_MAX;
	pBounds->fMax.x = pBounds->fMax.y = .0f;
	for (int32_t i = 0; i < faceInfo.size; ++i) {
		Vec2 *uv = pUvs + faceInfo.start + i;
		pBounds->fMin.x = uv->x < pBounds->fMin.x ? uv->x : pBounds->fMin.x;
		pBounds->fMin.y = uv->y < pBounds->fMin.y ? uv->y : pBounds->fMin.y;
		pBounds->fMax.x = uv->x > pBounds->fMax.x ? uv->x : pBounds->fMax.x;
		pBounds->fMax.y = uv->y > pBounds->fMax.y ? uv->y : pBounds->fMax.y;
	}
}

int32_t checkIfEdgeIsSeam(int32_t edgeIndex, FaceInfo face, int32_t loop,
                          Mesh *pMesh, EdgeVerts *pEdgeVerts) {
	int32_t *pVerts = pEdgeVerts[edgeIndex].verts;
	if (pVerts[1] < 0) {
		return 2;
	}
	else {
		int32_t whichLoop = pVerts[0] == face.start + loop;
		int32_t otherLoop = pVerts[whichLoop];
		int32_t iNext = (loop + 1) % face.size;
		int32_t nextBaseLoop = face.start + iNext;
		Vec2 uv = pMesh->pUvs[nextBaseLoop];
		Vec2 uvOther = pMesh->pUvs[otherLoop];
		int32_t isSeam = _(uv V2NOTEQL uvOther);
		if (isSeam) {
			return 1;
		}
	}
	return 0;
}

int32_t checkIfEdgeIsPreserve(Mesh* pMesh, int32_t edge) {
	return pMesh->pEdgePreserve ? pMesh->pEdgePreserve[edge] : 0;
}

static int32_t getOtherVert(int32_t i, int32_t faceSize, int8_t *pVertsRemoved) {
	int32_t ib = (i + 1) % faceSize;
	//search from i + 1 to facesize, and if non found,
	//then run again from 0 to facesize. If non found then,
	//return error
	int32_t attempts = 0;
	do {
		attempts++;
		for (; ib < faceSize; ++ib) {
			if (!pVertsRemoved[ib]) {
				return ib;
			}
		}
		ib = 0;
	} while (attempts == 1);
	return -1;
}

FaceTriangulated triangulateFace(RuvmAllocator alloc, FaceInfo baseFace, Mesh *pMesh) {
	FaceTriangulated outMesh = {0};
	int32_t triCount = baseFace.size - 2;
	outMesh.pTris = alloc.pMalloc(sizeof(int32_t) * triCount);
	int32_t loopCount = triCount * 3;
	outMesh.pLoops = alloc.pMalloc(sizeof(int32_t) * loopCount);
	
	int8_t *pVertsRemoved = alloc.pCalloc(baseFace.size, 1);
	int32_t loopsLeft = baseFace.size;
	do {
		//loop through ears, and find one with shortest edge
		float minDist = FLT_MAX; //min distance
		int32_t nextEar[3];
		for (int32_t i = 0; i < baseFace.size; ++i) {
			if (pVertsRemoved[i]) {
				continue;
			}
			int32_t ib = getOtherVert(i, baseFace.size, pVertsRemoved);
			int32_t ic = getOtherVert(ib, baseFace.size, pVertsRemoved);
			Vec2 uva = pMesh->pUvs[baseFace.start + i];
			Vec2 uvb = pMesh->pUvs[baseFace.start + ib];
			Vec2 uvc = pMesh->pUvs[baseFace.start + ic];
			int32_t windingDir = vec2WindingCompare(uva, uvb, uvc, 0);
			if (!windingDir) {
				continue;
			}
			Vec2 vDist = vec2Abs(_(uvc V2SUB uva));
			float dist = sqrt(vDist.x * vDist.x + vDist.y * vDist.y); //distance
			if (dist < minDist) {
				minDist = dist;
				nextEar[0] = i;
				nextEar[1] = ib;
				nextEar[2] = ic;
			}
		}
		outMesh.pTris[outMesh.triCount] = outMesh.loopCount;
		for (int32_t i = 0; i < 3; ++i) {
			//set to equal loop index, rather than vert index
			outMesh.pLoops[outMesh.loopCount] = nextEar[i];
			outMesh.loopCount++;
		}
		outMesh.triCount++;
		pVertsRemoved[nextEar[1]] = 1;
		loopsLeft--;
	} while (loopsLeft >= 3);
	alloc.pFree(pVertsRemoved);
	return outMesh;
}

int32_t getAttribSize(RuvmAttribType type) {
	switch (type) {
		case RUVM_ATTRIB_I8:
			return 1;
		case RUVM_ATTRIB_I16:
			return 2;
		case RUVM_ATTRIB_I32:
			return 4;
		case RUVM_ATTRIB_I64:
			return 8;
		case RUVM_ATTRIB_F32:
			return 4;
		case RUVM_ATTRIB_F64:
			return 8;
		case RUVM_ATTRIB_V2_I8:
			return 2;
		case RUVM_ATTRIB_V2_I16:
			return 4;
		case RUVM_ATTRIB_V2_I32:
			return 8;
		case RUVM_ATTRIB_V2_I64:
			return 16;
		case RUVM_ATTRIB_V2_F32:
			return 8;
		case RUVM_ATTRIB_V2_F64:
			return 16;
		case RUVM_ATTRIB_V3_I8:
			return 3;
		case RUVM_ATTRIB_V3_I16:
			return 6;
		case RUVM_ATTRIB_V3_I32:
			return 12;
		case RUVM_ATTRIB_V3_I64:
			return 24;
		case RUVM_ATTRIB_V3_F32:
			return 12;
		case RUVM_ATTRIB_V3_F64:
			return 24;
		case RUVM_ATTRIB_V4_I8:
			return 4;
		case RUVM_ATTRIB_V4_I16:
			return 8;
		case RUVM_ATTRIB_V4_I32:
			return 16;
		case RUVM_ATTRIB_V4_I64:
			return 32;
		case RUVM_ATTRIB_V4_F32:
			return 16;
		case RUVM_ATTRIB_V4_F64:
			return 32;
		case RUVM_ATTRIB_STRING:
			return RUVM_ATTRIB_STRING_MAX_LEN;
	}
}

RuvmAttrib *getAttrib(char *pName, RuvmAttrib *pAttribs, int32_t attribCount) {
	for (int32_t i = 0; i < attribCount; ++i) {
		if (0 == strncmp(pName, pAttribs[i].name, RUVM_ATTRIB_NAME_MAX_LEN)) {
			return pAttribs + i;
		}
	}
	return NULL;
}

Vec3 *attribAsV3(RuvmAttrib *pAttrib, int32_t index) {
	return (Vec3 *)pAttrib->pData + index;
}

Vec2 *attribAsV2(RuvmAttrib *pAttrib, int32_t index) {
	return (Vec2 *)pAttrib->pData + index;
}

int32_t *attribAsI32(RuvmAttrib *pAttrib, int32_t index) {
	return (int32_t *)pAttrib->pData + index;
}

void *attribAsVoid(RuvmAttrib *pAttrib, int32_t index) {
	switch (pAttrib->type) {
		case RUVM_ATTRIB_I8:
			return ((int8_t *)pAttrib->pData) + index;
		case RUVM_ATTRIB_I16:
			return ((int16_t *)pAttrib->pData) + index;
		case RUVM_ATTRIB_I32:
			return ((int32_t *)pAttrib->pData) + index;
		case RUVM_ATTRIB_I64:
			return ((int64_t *)pAttrib->pData) + index;
		case RUVM_ATTRIB_F32:
			return ((float *)pAttrib->pData) + index;
		case RUVM_ATTRIB_F64:
			return ((double *)pAttrib->pData) + index;
		case RUVM_ATTRIB_V2_I8:
			return ((int8_t (*)[2])pAttrib->pData) + index;
		case RUVM_ATTRIB_V2_I16:
			return ((int16_t (*)[2])pAttrib->pData) + index;
		case RUVM_ATTRIB_V2_I32:
			return ((int32_t (*)[2])pAttrib->pData) + index;
		case RUVM_ATTRIB_V2_I64:
			return ((int64_t (*)[2])pAttrib->pData) + index;
		case RUVM_ATTRIB_V2_F32:
			return ((float (*)[2])pAttrib->pData) + index;
		case RUVM_ATTRIB_V2_F64:
			return ((double (*)[2])pAttrib->pData) + index;
		case RUVM_ATTRIB_V3_I8:
			return ((int8_t (*)[3])pAttrib->pData) + index;
		case RUVM_ATTRIB_V3_I16:
			return ((int16_t (*)[3])pAttrib->pData) + index;
		case RUVM_ATTRIB_V3_I32:
			return ((int32_t (*)[3])pAttrib->pData) + index;
		case RUVM_ATTRIB_V3_I64:
			return ((int64_t (*)[3])pAttrib->pData) + index;
		case RUVM_ATTRIB_V3_F32:
			return ((float (*)[3])pAttrib->pData) + index;
		case RUVM_ATTRIB_V3_F64:
			return ((double (*)[3])pAttrib->pData) + index;
		case RUVM_ATTRIB_V4_I8:
			return ((int8_t (*)[4])pAttrib->pData) + index;
		case RUVM_ATTRIB_V4_I16:
			return ((int16_t (*)[4])pAttrib->pData) + index;
		case RUVM_ATTRIB_V4_I32:
			return ((int32_t (*)[4])pAttrib->pData) + index;
		case RUVM_ATTRIB_V4_I64:
			return ((int64_t (*)[4])pAttrib->pData) + index;
		case RUVM_ATTRIB_V4_F32:
			return ((float (*)[4])pAttrib->pData) + index;
		case RUVM_ATTRIB_V4_F64:
			return ((double (*)[4])pAttrib->pData) + index;
		case RUVM_ATTRIB_STRING:
			return ((char (*)[RUVM_ATTRIB_STRING_MAX_LEN])pAttrib->pData) + index;
	}
}

int32_t copyAttrib(RuvmAttrib *pDest, int32_t iDest,
                   RuvmAttrib *pSrc, int32_t iSrc) {
	if (pSrc->type != pDest->type) {
		return 1;
	}
	switch (pSrc->type) {
		case RUVM_ATTRIB_I8:
			((int8_t *)pSrc->pData)[iSrc] = ((int8_t *)pDest->pData)[iDest];
			break;
		case RUVM_ATTRIB_I16:
			((int16_t *)pSrc->pData)[iSrc] = ((int16_t *)pDest->pData)[iDest];
			break;
		case RUVM_ATTRIB_I32:
			((int32_t *)pSrc->pData)[iSrc] = ((int32_t *)pDest->pData)[iDest];
			break;
		case RUVM_ATTRIB_I64:
			((int64_t *)pSrc->pData)[iSrc] = ((int64_t *)pDest->pData)[iDest];
			break;
		case RUVM_ATTRIB_F32:
			((float *)pSrc->pData)[iSrc] = ((float *)pDest->pData)[iDest];
			break;
		case RUVM_ATTRIB_F64:
			((double *)pSrc->pData)[iSrc] = ((double *)pDest->pData)[iDest];
			break;
		case RUVM_ATTRIB_V2_I8:
			memcpy(((int8_t (*)[2])pDest->pData)[iDest],
			       ((int8_t (*)[2])pSrc->pData)[iSrc], sizeof(int8_t[2]));
			break;
		case RUVM_ATTRIB_V2_I16:
			memcpy(((int16_t (*)[2])pDest->pData)[iDest],
			       ((int16_t (*)[2])pSrc->pData)[iSrc], sizeof(int16_t[2]));
			break;
		case RUVM_ATTRIB_V2_I32:
			memcpy(((int32_t (*)[2])pDest->pData)[iDest],
			       ((int32_t (*)[2])pSrc->pData)[iSrc], sizeof(int32_t[2]));
			break;
		case RUVM_ATTRIB_V2_I64:
			memcpy(((int64_t (*)[2])pDest->pData)[iDest],
			       ((int64_t (*)[2])pSrc->pData)[iSrc], sizeof(int64_t[2]));
			break;
		case RUVM_ATTRIB_V2_F32:
			memcpy(((float (*)[2])pDest->pData)[iDest],
			       ((float (*)[2])pSrc->pData)[iSrc], sizeof(float[2]));
			break;
		case RUVM_ATTRIB_V2_F64:
			memcpy(((double (*)[2])pDest->pData)[iDest],
			       ((double (*)[2])pSrc->pData)[iSrc], sizeof(double[2]));
			break;
		case RUVM_ATTRIB_V3_I8:
			memcpy(((int8_t (*)[3])pDest->pData)[iDest],
			       ((int8_t (*)[3])pSrc->pData)[iSrc], sizeof(int8_t[3]));
			break;
		case RUVM_ATTRIB_V3_I16:
			memcpy(((int16_t (*)[3])pDest->pData)[iDest],
			       ((int16_t (*)[3])pSrc->pData)[iSrc], sizeof(int16_t[3]));
			break;
		case RUVM_ATTRIB_V3_I32:
			memcpy(((int32_t (*)[3])pDest->pData)[iDest],
			       ((int32_t (*)[3])pSrc->pData)[iSrc], sizeof(int32_t[3]));
			break;
		case RUVM_ATTRIB_V3_I64:
			memcpy(((int64_t (*)[3])pDest->pData)[iDest],
			       ((int64_t (*)[3])pSrc->pData)[iSrc], sizeof(int64_t[3]));
			break;
		case RUVM_ATTRIB_V3_F32:
			memcpy(((float (*)[3])pDest->pData)[iDest],
			       ((float (*)[3])pSrc->pData)[iSrc], sizeof(float[3]));
			break;
		case RUVM_ATTRIB_V3_F64:
			memcpy(((double (*)[3])pDest->pData)[iDest],
			       ((double (*)[3])pSrc->pData)[iSrc], sizeof(double[3]));
			break;
		case RUVM_ATTRIB_V4_I8:
			memcpy(((int8_t (*)[4])pDest->pData)[iDest],
			       ((int8_t (*)[4])pSrc->pData)[iSrc], sizeof(int8_t[4]));
			break;
		case RUVM_ATTRIB_V4_I16:
			memcpy(((int16_t (*)[4])pDest->pData)[iDest],
			       ((int16_t (*)[4])pSrc->pData)[iSrc], sizeof(int16_t[4]));
			break;
		case RUVM_ATTRIB_V4_I32:
			memcpy(((int32_t (*)[4])pDest->pData)[iDest],
			       ((int32_t (*)[4])pSrc->pData)[iSrc], sizeof(int32_t[4]));
			break;
		case RUVM_ATTRIB_V4_I64:
			memcpy(((int64_t (*)[4])pDest->pData)[iDest],
			       ((int64_t (*)[4])pSrc->pData)[iSrc], sizeof(int64_t[4]));
			break;
		case RUVM_ATTRIB_V4_F32:
			memcpy(((float (*)[4])pDest->pData)[iDest],
			       ((float (*)[4])pSrc->pData)[iSrc], sizeof(float[4]));
			break;
		case RUVM_ATTRIB_V4_F64:
			memcpy(((double (*)[4])pDest->pData)[iDest],
			       ((double (*)[4])pSrc->pData)[iSrc], sizeof(double[4]));
			break;
		case RUVM_ATTRIB_STRING:
			memcpy(((char (*)[RUVM_ATTRIB_STRING_MAX_LEN])pDest->pData)[iDest],
			       ((char (*)[RUVM_ATTRIB_STRING_MAX_LEN])pSrc->pData)[iSrc],
				   sizeof(double[RUVM_ATTRIB_STRING_MAX_LEN]));
			break;
	}
	return 0;
}

void copyAllAttribs(RuvmAttrib *pDest, int32_t iDest,
                    RuvmAttrib *pSrc, int32_t iSrc,
                    int32_t attribCount) {
	for (int32_t i = 0; i < attribCount; ++i) {
		copyAttrib(pDest + i, iDest, pSrc + i, iSrc);
	}
}

RuvmTypeDefault *getTypeDefaultConfig(RuvmTypeDefaultConfig *pConfig,
                                      RuvmAttribType type) {
	switch (type) {
		case RUVM_ATTRIB_I8:
			return &pConfig->i8;
		case RUVM_ATTRIB_I16:
			return &pConfig->i16;
		case RUVM_ATTRIB_I32:
			return &pConfig->i32;
		case RUVM_ATTRIB_I64:
			return &pConfig->i64;
		case RUVM_ATTRIB_F32:
			return &pConfig->f32;
		case RUVM_ATTRIB_F64:
			return &pConfig->f64;
		case RUVM_ATTRIB_V2_I8:
			return &pConfig->v2_i8;
		case RUVM_ATTRIB_V2_I16:
			return &pConfig->v2_i16;
		case RUVM_ATTRIB_V2_I32:
			return &pConfig->v2_i32;
		case RUVM_ATTRIB_V2_I64:
			return &pConfig->v2_i64;
		case RUVM_ATTRIB_V2_F32:
			return &pConfig->v2_f32;
		case RUVM_ATTRIB_V2_F64:
			return &pConfig->v2_f64;
		case RUVM_ATTRIB_V3_I8:
			return &pConfig->v3_i8;
		case RUVM_ATTRIB_V3_I16:
			return &pConfig->v3_i16;
		case RUVM_ATTRIB_V3_I32:
			return &pConfig->v3_i32;
		case RUVM_ATTRIB_V3_I64:
			return &pConfig->v3_i64;
		case RUVM_ATTRIB_V3_F32:
			return &pConfig->v3_f32;
		case RUVM_ATTRIB_V3_F64:
			return &pConfig->v3_f64;
		case RUVM_ATTRIB_V4_I8:
			return &pConfig->v4_i8;
		case RUVM_ATTRIB_V4_I16:
			return &pConfig->v4_i16;
		case RUVM_ATTRIB_V4_I32:
			return &pConfig->v4_i32;
		case RUVM_ATTRIB_V4_I64:
			return &pConfig->v4_i64;
		case RUVM_ATTRIB_V4_F32:
			return &pConfig->v4_f32;
		case RUVM_ATTRIB_V4_F64:
			return &pConfig->v4_f64;
		case RUVM_ATTRIB_STRING:
			return &pConfig->string;
	}
}

RuvmCommonAttrib *getCommonAttrib(RuvmCommonAttrib *pAttribs, int32_t attribCount,
                                  char *pName) {
	//this is it's own function (rather than a general getAttrib function),
	//because of the below todo:
	//TODO replace linear search with hash table.
	for (int32_t i = 0; i < attribCount; ++i) {
		if (0 == strncmp(pName, pAttribs[i].name, RUVM_ATTRIB_NAME_MAX_LEN)) {
			return pAttribs + i;
		}
	}
	return NULL;
}

void interpolateAttrib(RuvmAttrib *pDest, int32_t iDest, RuvmAttrib *pSrc,
                       int32_t iSrcA, int32_t iSrcB, int32_t iSrcC, Vec3 bc) {
	if (pDest->type != pSrc->type) {
		printf("Type mismatch in interpolateAttrib\n");
		//TODO remove all uses of abort(), and add proper exception handling
		abort();
	}
	RuvmAttribType type = pDest->type;
	switch (type) {
		case RUVM_ATTRIB_I8:
			INTERPOLATE_SCALAR(int8_t, pDest, iDest, pSrc, iSrcA, iSrcB, iSrcC, bc);
			break;
		case RUVM_ATTRIB_I16:
			INTERPOLATE_SCALAR(int16_t, pDest, iDest, pSrc, iSrcA, iSrcB, iSrcC, bc);
			break;
		case RUVM_ATTRIB_I32:
			INTERPOLATE_SCALAR(int32_t, pDest, iDest, pSrc, iSrcA, iSrcB, iSrcC, bc);
			break;
		case RUVM_ATTRIB_I64:
			INTERPOLATE_SCALAR(int64_t, pDest, iDest, pSrc, iSrcA, iSrcB, iSrcC, bc);
			break;
		case RUVM_ATTRIB_F32:
			INTERPOLATE_SCALAR(float, pDest, iDest, pSrc, iSrcA, iSrcB, iSrcC, bc);
			break;
		case RUVM_ATTRIB_F64:
			INTERPOLATE_SCALAR(double, pDest, iDest, pSrc, iSrcA, iSrcB, iSrcC, bc);
			break;
		case RUVM_ATTRIB_V2_I8:
			INTERPOLATE_V2(int8_t, pDest, iDest, pSrc, iSrcA, iSrcB, iSrcC, bc);
			break;
		case RUVM_ATTRIB_V2_I16:
			INTERPOLATE_V2(int16_t, pDest, iDest, pSrc, iSrcA, iSrcB, iSrcC, bc);
			break;
		case RUVM_ATTRIB_V2_I32:
			INTERPOLATE_V2(int32_t, pDest, iDest, pSrc, iSrcA, iSrcB, iSrcC, bc);
			break;
		case RUVM_ATTRIB_V2_I64:
			INTERPOLATE_V2(int64_t, pDest, iDest, pSrc, iSrcA, iSrcB, iSrcC, bc);
			break;
		case RUVM_ATTRIB_V2_F32:
			INTERPOLATE_V2(float, pDest, iDest, pSrc, iSrcA, iSrcB, iSrcC, bc);
			break;
		case RUVM_ATTRIB_V2_F64:
			INTERPOLATE_V2(double, pDest, iDest, pSrc, iSrcA, iSrcB, iSrcC, bc);
			break;
		case RUVM_ATTRIB_V3_I8:
			INTERPOLATE_V3(int8_t, pDest, iDest, pSrc, iSrcA, iSrcB, iSrcC, bc);
			break;
		case RUVM_ATTRIB_V3_I16:
			INTERPOLATE_V3(int16_t, pDest, iDest, pSrc, iSrcA, iSrcB, iSrcC, bc);
			break;
		case RUVM_ATTRIB_V3_I32:
			INTERPOLATE_V3(int32_t, pDest, iDest, pSrc, iSrcA, iSrcB, iSrcC, bc);
			break;
		case RUVM_ATTRIB_V3_I64:
			INTERPOLATE_V3(int64_t, pDest, iDest, pSrc, iSrcA, iSrcB, iSrcC, bc);
			break;
		case RUVM_ATTRIB_V3_F32:
			INTERPOLATE_V3(float, pDest, iDest, pSrc, iSrcA, iSrcB, iSrcC, bc);
			break;
		case RUVM_ATTRIB_V3_F64:
			INTERPOLATE_V3(double, pDest, iDest, pSrc, iSrcA, iSrcB, iSrcC, bc);
			break;
		case RUVM_ATTRIB_V4_I8:
			//TODO using unsigned here temporarily for vert color
			//make a proper set of attrib types for unsigned types
			INTERPOLATE_V4(uint8_t, pDest, iDest, pSrc, iSrcA, iSrcB, iSrcC, bc);
			break;
		case RUVM_ATTRIB_V4_I16:
			INTERPOLATE_V4(int16_t, pDest, iDest, pSrc, iSrcA, iSrcB, iSrcC, bc);
			break;
		case RUVM_ATTRIB_V4_I32:
			INTERPOLATE_V4(int32_t, pDest, iDest, pSrc, iSrcA, iSrcB, iSrcC, bc);
			break;
		case RUVM_ATTRIB_V4_I64:
			INTERPOLATE_V4(int64_t, pDest, iDest, pSrc, iSrcA, iSrcB, iSrcC, bc);
			break;
		case RUVM_ATTRIB_V4_F32:
			INTERPOLATE_V4(float, pDest, iDest, pSrc, iSrcA, iSrcB, iSrcC, bc);
			break;
		case RUVM_ATTRIB_V4_F64:
			INTERPOLATE_V4(double, pDest, iDest, pSrc, iSrcA, iSrcB, iSrcC, bc);
			break;
		case RUVM_ATTRIB_STRING:
			break;
	}
}

static void appendOnNonString() {
	printf("Blend mode append reached on non string attrib!\n");
	abort();
}

//TODO replace with a callback function in blendConfig. This many switch cases is absurd
void blendAttribs(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                  RuvmAttrib *pB, int32_t iB, RuvmBlendConfig blendConfig) {
	RuvmAttribType type = pD->type;
	switch (type) {
		case RUVM_ATTRIB_I8:
			switch (blendConfig.blend) {
				case RUVM_BLEND_REPLACE:
					BLEND_REPLACE(int8_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_MULTIPLY:
					BLEND_MULTIPLY(int8_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_DIVIDE:
					BLEND_DIVIDE(int8_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_ADD:
					BLEND_ADD(int8_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_SUBTRACT:
					BLEND_SUBTRACT(int8_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_ADD_SUB:
					BLEND_ADD_SUB(int8_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_LIGHTEN:
					BLEND_LIGHTEN(int8_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_DARKEN:
					BLEND_DARKEN(int8_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_OVERLAY:
					BLEND_OVERLAY(int8_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_SOFT_LIGHT:
					BLEND_SOFT_LIGHT(int8_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_COLOR_DODGE:
					BLEND_COLOR_DODGE(int8_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_APPEND:
					appendOnNonString();
					break;
			}
			break;
		case RUVM_ATTRIB_I16:
			switch (blendConfig.blend) {
				case RUVM_BLEND_REPLACE:
					BLEND_REPLACE(int16_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_MULTIPLY:
					BLEND_MULTIPLY(int16_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_DIVIDE:
					BLEND_DIVIDE(int16_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_ADD:
					BLEND_ADD(int16_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_SUBTRACT:
					BLEND_SUBTRACT(int16_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_ADD_SUB:
					BLEND_ADD_SUB(int16_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_LIGHTEN:
					BLEND_LIGHTEN(int16_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_DARKEN:
					BLEND_DARKEN(int16_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_OVERLAY:
					BLEND_OVERLAY(int16_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_SOFT_LIGHT:
					BLEND_SOFT_LIGHT(int16_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_COLOR_DODGE:
					BLEND_COLOR_DODGE(int16_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_APPEND:
					appendOnNonString();
					break;
			}
			break;
		case RUVM_ATTRIB_I32:
			switch (blendConfig.blend) {
				case RUVM_BLEND_REPLACE:
					BLEND_REPLACE(int32_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_MULTIPLY:
					BLEND_MULTIPLY(int32_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_DIVIDE:
					BLEND_DIVIDE(int32_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_ADD:
					BLEND_ADD(int32_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_SUBTRACT:
					BLEND_SUBTRACT(int32_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_ADD_SUB:
					BLEND_ADD_SUB(int32_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_LIGHTEN:
					BLEND_LIGHTEN(int32_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_DARKEN:
					BLEND_DARKEN(int32_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_OVERLAY:
					BLEND_OVERLAY(int32_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_SOFT_LIGHT:
					BLEND_SOFT_LIGHT(int32_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_COLOR_DODGE:
					BLEND_COLOR_DODGE(int32_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_APPEND:
					appendOnNonString();
					break;
			}
			break;
		case RUVM_ATTRIB_I64:
			switch (blendConfig.blend) {
				case RUVM_BLEND_REPLACE:
					BLEND_REPLACE(int64_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_MULTIPLY:
					BLEND_MULTIPLY(int64_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_DIVIDE:
					BLEND_DIVIDE(int64_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_ADD:
					BLEND_ADD(int64_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_SUBTRACT:
					BLEND_SUBTRACT(int64_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_ADD_SUB:
					BLEND_ADD_SUB(int64_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_LIGHTEN:
					BLEND_LIGHTEN(int64_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_DARKEN:
					BLEND_DARKEN(int64_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_OVERLAY:
					BLEND_OVERLAY(int64_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_SOFT_LIGHT:
					BLEND_SOFT_LIGHT(int64_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_COLOR_DODGE:
					BLEND_COLOR_DODGE(int64_t, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_APPEND:
					appendOnNonString();
					break;
			}
			break;
		case RUVM_ATTRIB_F32:
			switch (blendConfig.blend) {
				case RUVM_BLEND_REPLACE:
					BLEND_REPLACE(float, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_MULTIPLY:
					BLEND_MULTIPLY(float, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_DIVIDE:
					BLEND_DIVIDE(float, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_ADD:
					BLEND_ADD(float, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_SUBTRACT:
					BLEND_SUBTRACT(float, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_ADD_SUB:
					BLEND_ADD_SUB(float, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_LIGHTEN:
					BLEND_LIGHTEN(float, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_DARKEN:
					BLEND_DARKEN(float, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_OVERLAY:
					BLEND_OVERLAY(float, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_SOFT_LIGHT:
					BLEND_SOFT_LIGHT(float, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_COLOR_DODGE:
					BLEND_COLOR_DODGE(float, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_APPEND:
					appendOnNonString();
					break;
			}
			break;
		case RUVM_ATTRIB_F64:
			switch (blendConfig.blend) {
				case RUVM_BLEND_REPLACE:
					BLEND_REPLACE(double, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_MULTIPLY:
					BLEND_MULTIPLY(double, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_DIVIDE:
					BLEND_DIVIDE(double, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_ADD:
					BLEND_ADD(double, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_SUBTRACT:
					BLEND_SUBTRACT(double, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_ADD_SUB:
					BLEND_ADD_SUB(double, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_LIGHTEN:
					BLEND_LIGHTEN(double, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_DARKEN:
					BLEND_DARKEN(double, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_OVERLAY:
					BLEND_OVERLAY(double, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_SOFT_LIGHT:
					BLEND_SOFT_LIGHT(double, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_COLOR_DODGE:
					BLEND_COLOR_DODGE(double, pD, iD, pA, iA, pB, iB, 1, 0);
					break;
				case RUVM_BLEND_APPEND:
					appendOnNonString();
					break;
			}
			break;
		case RUVM_ATTRIB_V2_I8:
			switch (blendConfig.blend) {
				case RUVM_BLEND_REPLACE:
					BLEND_REPLACE(int8_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_REPLACE(int8_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_MULTIPLY:
					BLEND_MULTIPLY(int8_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_MULTIPLY(int8_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_DIVIDE:
					BLEND_DIVIDE(int8_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_DIVIDE(int8_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_ADD:
					BLEND_ADD(int8_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_ADD(int8_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_SUBTRACT:
					BLEND_SUBTRACT(int8_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_SUBTRACT(int8_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_ADD_SUB:
					BLEND_ADD_SUB(int8_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_ADD_SUB(int8_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_LIGHTEN:
					BLEND_LIGHTEN(int8_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_LIGHTEN(int8_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_DARKEN:
					BLEND_DARKEN(int8_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_DARKEN(int8_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_OVERLAY:
					BLEND_OVERLAY(int8_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_OVERLAY(int8_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_SOFT_LIGHT:
					BLEND_SOFT_LIGHT(int8_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_SOFT_LIGHT(int8_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_COLOR_DODGE:
					BLEND_COLOR_DODGE(int8_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_COLOR_DODGE(int8_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_APPEND:
					appendOnNonString();
					break;
			}
			break;
		case RUVM_ATTRIB_V2_I16:
			switch (blendConfig.blend) {
				case RUVM_BLEND_REPLACE:
					BLEND_REPLACE(int16_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_REPLACE(int16_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_MULTIPLY:
					BLEND_MULTIPLY(int16_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_MULTIPLY(int16_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_DIVIDE:
					BLEND_DIVIDE(int16_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_DIVIDE(int16_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_ADD:
					BLEND_ADD(int16_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_ADD(int16_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_SUBTRACT:
					BLEND_SUBTRACT(int16_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_SUBTRACT(int16_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_ADD_SUB:
					BLEND_ADD_SUB(int16_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_ADD_SUB(int16_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_LIGHTEN:
					BLEND_LIGHTEN(int16_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_LIGHTEN(int16_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_DARKEN:
					BLEND_DARKEN(int16_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_DARKEN(int16_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_OVERLAY:
					BLEND_OVERLAY(int16_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_OVERLAY(int16_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_SOFT_LIGHT:
					BLEND_SOFT_LIGHT(int16_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_SOFT_LIGHT(int16_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_COLOR_DODGE:
					BLEND_COLOR_DODGE(int16_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_COLOR_DODGE(int16_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_APPEND:
					appendOnNonString();
					break;
			}
			break;
		case RUVM_ATTRIB_V2_I32:
			switch (blendConfig.blend) {
				case RUVM_BLEND_REPLACE:
					BLEND_REPLACE(int32_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_REPLACE(int32_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_MULTIPLY:
					BLEND_MULTIPLY(int32_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_MULTIPLY(int32_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_DIVIDE:
					BLEND_DIVIDE(int32_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_DIVIDE(int32_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_ADD:
					BLEND_ADD(int32_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_ADD(int32_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_SUBTRACT:
					BLEND_SUBTRACT(int32_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_SUBTRACT(int32_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_ADD_SUB:
					BLEND_ADD_SUB(int32_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_ADD_SUB(int32_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_LIGHTEN:
					BLEND_LIGHTEN(int32_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_LIGHTEN(int32_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_DARKEN:
					BLEND_DARKEN(int32_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_DARKEN(int32_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_OVERLAY:
					BLEND_OVERLAY(int32_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_OVERLAY(int32_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_SOFT_LIGHT:
					BLEND_SOFT_LIGHT(int32_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_SOFT_LIGHT(int32_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_COLOR_DODGE:
					BLEND_COLOR_DODGE(int32_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_COLOR_DODGE(int32_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_APPEND:
					appendOnNonString();
					break;
			}
			break;
		case RUVM_ATTRIB_V2_I64:
			switch (blendConfig.blend) {
				case RUVM_BLEND_REPLACE:
					BLEND_REPLACE(int64_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_REPLACE(int64_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_MULTIPLY:
					BLEND_MULTIPLY(int64_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_MULTIPLY(int64_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_DIVIDE:
					BLEND_DIVIDE(int64_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_DIVIDE(int64_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_ADD:
					BLEND_ADD(int64_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_ADD(int64_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_SUBTRACT:
					BLEND_SUBTRACT(int64_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_SUBTRACT(int64_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_ADD_SUB:
					BLEND_ADD_SUB(int64_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_ADD_SUB(int64_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_LIGHTEN:
					BLEND_LIGHTEN(int64_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_LIGHTEN(int64_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_DARKEN:
					BLEND_DARKEN(int64_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_DARKEN(int64_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_OVERLAY:
					BLEND_OVERLAY(int64_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_OVERLAY(int64_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_SOFT_LIGHT:
					BLEND_SOFT_LIGHT(int64_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_SOFT_LIGHT(int64_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_COLOR_DODGE:
					BLEND_COLOR_DODGE(int64_t, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_COLOR_DODGE(int64_t, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_APPEND:
					appendOnNonString();
					break;
			}
			break;
		case RUVM_ATTRIB_V2_F32:
			switch (blendConfig.blend) {
				case RUVM_BLEND_REPLACE:
					BLEND_REPLACE(float, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_REPLACE(float, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_MULTIPLY:
					BLEND_MULTIPLY(float, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_MULTIPLY(float, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_DIVIDE:
					BLEND_DIVIDE(float, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_DIVIDE(float, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_ADD:
					BLEND_ADD(float, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_ADD(float, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_SUBTRACT:
					BLEND_SUBTRACT(float, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_SUBTRACT(float, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_ADD_SUB:
					BLEND_ADD_SUB(float, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_ADD_SUB(float, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_LIGHTEN:
					BLEND_LIGHTEN(float, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_LIGHTEN(float, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_DARKEN:
					BLEND_DARKEN(float, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_DARKEN(float, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_OVERLAY:
					BLEND_OVERLAY(float, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_OVERLAY(float, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_SOFT_LIGHT:
					BLEND_SOFT_LIGHT(float, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_SOFT_LIGHT(float, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_COLOR_DODGE:
					BLEND_COLOR_DODGE(float, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_COLOR_DODGE(float, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_APPEND:
					appendOnNonString();
					break;
			}
			break;
		case RUVM_ATTRIB_V2_F64:
			switch (blendConfig.blend) {
				case RUVM_BLEND_REPLACE:
					BLEND_REPLACE(double, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_REPLACE(double, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_MULTIPLY:
					BLEND_MULTIPLY(double, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_MULTIPLY(double, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_DIVIDE:
					BLEND_DIVIDE(double, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_DIVIDE(double, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_ADD:
					BLEND_ADD(double, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_ADD(double, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_SUBTRACT:
					BLEND_SUBTRACT(double, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_SUBTRACT(double, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_ADD_SUB:
					BLEND_ADD_SUB(double, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_ADD_SUB(double, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_LIGHTEN:
					BLEND_LIGHTEN(double, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_LIGHTEN(double, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_DARKEN:
					BLEND_DARKEN(double, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_DARKEN(double, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_OVERLAY:
					BLEND_OVERLAY(double, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_OVERLAY(double, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_SOFT_LIGHT:
					BLEND_SOFT_LIGHT(double, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_SOFT_LIGHT(double, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_COLOR_DODGE:
					BLEND_COLOR_DODGE(double, pD, iD, pA, iA, pB, iB, 2, 0);
					BLEND_COLOR_DODGE(double, pD, iD, pA, iA, pB, iB, 2, 1);
					break;
				case RUVM_BLEND_APPEND:
					appendOnNonString();
					break;
			}
			break;
		case RUVM_ATTRIB_V3_I8:
			switch (blendConfig.blend) {
				case RUVM_BLEND_REPLACE:
					BLEND_REPLACE(int8_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_REPLACE(int8_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_REPLACE(int8_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_MULTIPLY:
					BLEND_MULTIPLY(int8_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_MULTIPLY(int8_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_MULTIPLY(int8_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_DIVIDE:
					BLEND_DIVIDE(int8_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_DIVIDE(int8_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_DIVIDE(int8_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_ADD:
					BLEND_ADD(int8_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_ADD(int8_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_ADD(int8_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_SUBTRACT:
					BLEND_SUBTRACT(int8_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_SUBTRACT(int8_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_SUBTRACT(int8_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_ADD_SUB:
					BLEND_ADD_SUB(int8_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_ADD_SUB(int8_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_ADD_SUB(int8_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_LIGHTEN:
					BLEND_LIGHTEN(int8_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_LIGHTEN(int8_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_LIGHTEN(int8_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_DARKEN:
					BLEND_DARKEN(int8_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_DARKEN(int8_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_DARKEN(int8_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_OVERLAY:
					BLEND_OVERLAY(int8_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_OVERLAY(int8_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_OVERLAY(int8_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_SOFT_LIGHT:
					BLEND_SOFT_LIGHT(int8_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_SOFT_LIGHT(int8_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_SOFT_LIGHT(int8_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_COLOR_DODGE:
					BLEND_COLOR_DODGE(int8_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_COLOR_DODGE(int8_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_COLOR_DODGE(int8_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_APPEND:
					appendOnNonString();
					break;
			}
			break;
		case RUVM_ATTRIB_V3_I16:
			switch (blendConfig.blend) {
				case RUVM_BLEND_REPLACE:
					BLEND_REPLACE(int16_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_REPLACE(int16_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_REPLACE(int16_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_MULTIPLY:
					BLEND_MULTIPLY(int16_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_MULTIPLY(int16_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_MULTIPLY(int16_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_DIVIDE:
					BLEND_DIVIDE(int16_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_DIVIDE(int16_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_DIVIDE(int16_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_ADD:
					BLEND_ADD(int16_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_ADD(int16_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_ADD(int16_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_SUBTRACT:
					BLEND_SUBTRACT(int16_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_SUBTRACT(int16_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_SUBTRACT(int16_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_ADD_SUB:
					BLEND_ADD_SUB(int16_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_ADD_SUB(int16_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_ADD_SUB(int16_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_LIGHTEN:
					BLEND_LIGHTEN(int16_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_LIGHTEN(int16_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_LIGHTEN(int16_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_DARKEN:
					BLEND_DARKEN(int16_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_DARKEN(int16_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_DARKEN(int16_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_OVERLAY:
					BLEND_OVERLAY(int16_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_OVERLAY(int16_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_OVERLAY(int16_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_SOFT_LIGHT:
					BLEND_SOFT_LIGHT(int16_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_SOFT_LIGHT(int16_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_SOFT_LIGHT(int16_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_COLOR_DODGE:
					BLEND_COLOR_DODGE(int16_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_COLOR_DODGE(int16_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_COLOR_DODGE(int16_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_APPEND:
					appendOnNonString();
					break;
			}
			break;
		case RUVM_ATTRIB_V3_I32:
			switch (blendConfig.blend) {
				case RUVM_BLEND_REPLACE:
					BLEND_REPLACE(int32_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_REPLACE(int32_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_REPLACE(int32_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_MULTIPLY:
					BLEND_MULTIPLY(int32_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_MULTIPLY(int32_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_MULTIPLY(int32_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_DIVIDE:
					BLEND_DIVIDE(int32_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_DIVIDE(int32_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_DIVIDE(int32_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_ADD:
					BLEND_ADD(int32_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_ADD(int32_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_ADD(int32_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_SUBTRACT:
					BLEND_SUBTRACT(int32_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_SUBTRACT(int32_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_SUBTRACT(int32_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_ADD_SUB:
					BLEND_ADD_SUB(int32_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_ADD_SUB(int32_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_ADD_SUB(int32_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_LIGHTEN:
					BLEND_LIGHTEN(int32_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_LIGHTEN(int32_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_LIGHTEN(int32_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_DARKEN:
					BLEND_DARKEN(int32_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_DARKEN(int32_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_DARKEN(int32_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_OVERLAY:
					BLEND_OVERLAY(int32_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_OVERLAY(int32_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_OVERLAY(int32_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_SOFT_LIGHT:
					BLEND_SOFT_LIGHT(int32_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_SOFT_LIGHT(int32_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_SOFT_LIGHT(int32_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_COLOR_DODGE:
					BLEND_COLOR_DODGE(int32_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_COLOR_DODGE(int32_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_COLOR_DODGE(int32_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_APPEND:
					appendOnNonString();
					break;
			}
			break;
		case RUVM_ATTRIB_V3_I64:
			switch (blendConfig.blend) {
				case RUVM_BLEND_REPLACE:
					BLEND_REPLACE(int64_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_REPLACE(int64_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_REPLACE(int64_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_MULTIPLY:
					BLEND_MULTIPLY(int64_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_MULTIPLY(int64_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_MULTIPLY(int64_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_DIVIDE:
					BLEND_DIVIDE(int64_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_DIVIDE(int64_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_DIVIDE(int64_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_ADD:
					BLEND_ADD(int64_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_ADD(int64_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_ADD(int64_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_SUBTRACT:
					BLEND_SUBTRACT(int64_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_SUBTRACT(int64_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_SUBTRACT(int64_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_ADD_SUB:
					BLEND_ADD_SUB(int64_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_ADD_SUB(int64_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_ADD_SUB(int64_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_LIGHTEN:
					BLEND_LIGHTEN(int64_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_LIGHTEN(int64_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_LIGHTEN(int64_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_DARKEN:
					BLEND_DARKEN(int64_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_DARKEN(int64_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_DARKEN(int64_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_OVERLAY:
					BLEND_OVERLAY(int64_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_OVERLAY(int64_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_OVERLAY(int64_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_SOFT_LIGHT:
					BLEND_SOFT_LIGHT(int64_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_SOFT_LIGHT(int64_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_SOFT_LIGHT(int64_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_COLOR_DODGE:
					BLEND_COLOR_DODGE(int64_t, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_COLOR_DODGE(int64_t, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_COLOR_DODGE(int64_t, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_APPEND:
					appendOnNonString();
					break;
			}
			break;
		case RUVM_ATTRIB_V3_F32:
			switch (blendConfig.blend) {
				case RUVM_BLEND_REPLACE:
					BLEND_REPLACE(float, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_REPLACE(float, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_REPLACE(float, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_MULTIPLY:
					BLEND_MULTIPLY(float, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_MULTIPLY(float, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_MULTIPLY(float, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_DIVIDE:
					BLEND_DIVIDE(float, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_DIVIDE(float, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_DIVIDE(float, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_ADD:
					BLEND_ADD(float, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_ADD(float, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_ADD(float, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_SUBTRACT:
					BLEND_SUBTRACT(float, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_SUBTRACT(float, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_SUBTRACT(float, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_ADD_SUB:
					BLEND_ADD_SUB(float, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_ADD_SUB(float, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_ADD_SUB(float, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_LIGHTEN:
					BLEND_LIGHTEN(float, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_LIGHTEN(float, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_LIGHTEN(float, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_DARKEN:
					BLEND_DARKEN(float, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_DARKEN(float, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_DARKEN(float, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_OVERLAY:
					BLEND_OVERLAY(float, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_OVERLAY(float, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_OVERLAY(float, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_SOFT_LIGHT:
					BLEND_SOFT_LIGHT(float, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_SOFT_LIGHT(float, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_SOFT_LIGHT(float, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_COLOR_DODGE:
					BLEND_COLOR_DODGE(float, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_COLOR_DODGE(float, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_COLOR_DODGE(float, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_APPEND:
					appendOnNonString();
					break;
			}
			break;
		case RUVM_ATTRIB_V3_F64:
			switch (blendConfig.blend) {
				case RUVM_BLEND_REPLACE:
					BLEND_REPLACE(double, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_REPLACE(double, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_REPLACE(double, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_MULTIPLY:
					BLEND_MULTIPLY(double, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_MULTIPLY(double, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_MULTIPLY(double, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_DIVIDE:
					BLEND_DIVIDE(double, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_DIVIDE(double, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_DIVIDE(double, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_ADD:
					BLEND_ADD(double, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_ADD(double, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_ADD(double, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_SUBTRACT:
					BLEND_SUBTRACT(double, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_SUBTRACT(double, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_SUBTRACT(double, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_ADD_SUB:
					BLEND_ADD_SUB(double, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_ADD_SUB(double, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_ADD_SUB(double, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_LIGHTEN:
					BLEND_LIGHTEN(double, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_LIGHTEN(double, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_LIGHTEN(double, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_DARKEN:
					BLEND_DARKEN(double, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_DARKEN(double, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_DARKEN(double, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_OVERLAY:
					BLEND_OVERLAY(double, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_OVERLAY(double, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_OVERLAY(double, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_SOFT_LIGHT:
					BLEND_SOFT_LIGHT(double, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_SOFT_LIGHT(double, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_SOFT_LIGHT(double, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_COLOR_DODGE:
					BLEND_COLOR_DODGE(double, pD, iD, pA, iA, pB, iB, 3, 0);
					BLEND_COLOR_DODGE(double, pD, iD, pA, iA, pB, iB, 3, 1);
					BLEND_COLOR_DODGE(double, pD, iD, pA, iA, pB, iB, 3, 2);
					break;
				case RUVM_BLEND_APPEND:
					appendOnNonString();
					break;
			}
			break;
		case RUVM_ATTRIB_V4_I8:
			switch (blendConfig.blend) {
				case RUVM_BLEND_REPLACE:
					BLEND_REPLACE(int8_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_REPLACE(int8_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_REPLACE(int8_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_REPLACE(int8_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_MULTIPLY:
					BLEND_MULTIPLY(int8_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_MULTIPLY(int8_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_MULTIPLY(int8_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_MULTIPLY(int8_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_DIVIDE:
					BLEND_DIVIDE(int8_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_DIVIDE(int8_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_DIVIDE(int8_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_DIVIDE(int8_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_ADD:
					BLEND_ADD(int8_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_ADD(int8_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_ADD(int8_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_ADD(int8_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_SUBTRACT:
					BLEND_SUBTRACT(int8_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_SUBTRACT(int8_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_SUBTRACT(int8_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_SUBTRACT(int8_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_ADD_SUB:
					BLEND_ADD_SUB(int8_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_ADD_SUB(int8_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_ADD_SUB(int8_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_ADD_SUB(int8_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_LIGHTEN:
					BLEND_LIGHTEN(int8_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_LIGHTEN(int8_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_LIGHTEN(int8_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_LIGHTEN(int8_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_DARKEN:
					BLEND_DARKEN(int8_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_DARKEN(int8_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_DARKEN(int8_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_DARKEN(int8_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_OVERLAY:
					BLEND_OVERLAY(int8_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_OVERLAY(int8_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_OVERLAY(int8_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_OVERLAY(int8_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_SOFT_LIGHT:
					BLEND_SOFT_LIGHT(int8_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_SOFT_LIGHT(int8_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_SOFT_LIGHT(int8_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_SOFT_LIGHT(int8_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_COLOR_DODGE:
					BLEND_COLOR_DODGE(int8_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_COLOR_DODGE(int8_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_COLOR_DODGE(int8_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_COLOR_DODGE(int8_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_APPEND:
					appendOnNonString();
					break;
			}
			break;
		case RUVM_ATTRIB_V4_I16:
			switch (blendConfig.blend) {
				case RUVM_BLEND_REPLACE:
					BLEND_REPLACE(int16_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_REPLACE(int16_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_REPLACE(int16_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_REPLACE(int16_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_MULTIPLY:
					BLEND_MULTIPLY(int16_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_MULTIPLY(int16_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_MULTIPLY(int16_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_MULTIPLY(int16_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_DIVIDE:
					BLEND_DIVIDE(int16_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_DIVIDE(int16_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_DIVIDE(int16_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_DIVIDE(int16_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_ADD:
					BLEND_ADD(int16_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_ADD(int16_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_ADD(int16_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_ADD(int16_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_SUBTRACT:
					BLEND_SUBTRACT(int16_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_SUBTRACT(int16_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_SUBTRACT(int16_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_SUBTRACT(int16_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_ADD_SUB:
					BLEND_ADD_SUB(int16_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_ADD_SUB(int16_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_ADD_SUB(int16_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_ADD_SUB(int16_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_LIGHTEN:
					BLEND_LIGHTEN(int16_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_LIGHTEN(int16_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_LIGHTEN(int16_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_LIGHTEN(int16_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_DARKEN:
					BLEND_DARKEN(int16_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_DARKEN(int16_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_DARKEN(int16_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_DARKEN(int16_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_OVERLAY:
					BLEND_OVERLAY(int16_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_OVERLAY(int16_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_OVERLAY(int16_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_OVERLAY(int16_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_SOFT_LIGHT:
					BLEND_SOFT_LIGHT(int16_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_SOFT_LIGHT(int16_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_SOFT_LIGHT(int16_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_SOFT_LIGHT(int16_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_COLOR_DODGE:
					BLEND_COLOR_DODGE(int16_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_COLOR_DODGE(int16_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_COLOR_DODGE(int16_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_COLOR_DODGE(int16_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_APPEND:
					appendOnNonString();
					break;
			}
			break;
		case RUVM_ATTRIB_V4_I32:
			switch (blendConfig.blend) {
				case RUVM_BLEND_REPLACE:
					BLEND_REPLACE(int32_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_REPLACE(int32_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_REPLACE(int32_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_REPLACE(int32_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_MULTIPLY:
					BLEND_MULTIPLY(int32_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_MULTIPLY(int32_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_MULTIPLY(int32_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_MULTIPLY(int32_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_DIVIDE:
					BLEND_DIVIDE(int32_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_DIVIDE(int32_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_DIVIDE(int32_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_DIVIDE(int32_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_ADD:
					BLEND_ADD(int32_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_ADD(int32_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_ADD(int32_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_ADD(int32_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_SUBTRACT:
					BLEND_SUBTRACT(int32_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_SUBTRACT(int32_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_SUBTRACT(int32_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_SUBTRACT(int32_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_ADD_SUB:
					BLEND_ADD_SUB(int32_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_ADD_SUB(int32_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_ADD_SUB(int32_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_ADD_SUB(int32_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_LIGHTEN:
					BLEND_LIGHTEN(int32_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_LIGHTEN(int32_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_LIGHTEN(int32_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_LIGHTEN(int32_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_DARKEN:
					BLEND_DARKEN(int32_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_DARKEN(int32_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_DARKEN(int32_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_DARKEN(int32_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_OVERLAY:
					BLEND_OVERLAY(int32_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_OVERLAY(int32_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_OVERLAY(int32_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_OVERLAY(int32_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_SOFT_LIGHT:
					BLEND_SOFT_LIGHT(int32_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_SOFT_LIGHT(int32_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_SOFT_LIGHT(int32_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_SOFT_LIGHT(int32_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_COLOR_DODGE:
					BLEND_COLOR_DODGE(int32_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_COLOR_DODGE(int32_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_COLOR_DODGE(int32_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_COLOR_DODGE(int32_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_APPEND:
					appendOnNonString();
					break;
			}
			break;
		case RUVM_ATTRIB_V4_I64:
			switch (blendConfig.blend) {
				case RUVM_BLEND_REPLACE:
					BLEND_REPLACE(int64_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_REPLACE(int64_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_REPLACE(int64_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_REPLACE(int64_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_MULTIPLY:
					BLEND_MULTIPLY(int64_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_MULTIPLY(int64_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_MULTIPLY(int64_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_MULTIPLY(int64_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_DIVIDE:
					BLEND_DIVIDE(int64_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_DIVIDE(int64_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_DIVIDE(int64_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_DIVIDE(int64_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_ADD:
					BLEND_ADD(int64_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_ADD(int64_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_ADD(int64_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_ADD(int64_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_SUBTRACT:
					BLEND_SUBTRACT(int64_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_SUBTRACT(int64_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_SUBTRACT(int64_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_SUBTRACT(int64_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_ADD_SUB:
					BLEND_ADD_SUB(int64_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_ADD_SUB(int64_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_ADD_SUB(int64_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_ADD_SUB(int64_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_LIGHTEN:
					BLEND_LIGHTEN(int64_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_LIGHTEN(int64_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_LIGHTEN(int64_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_LIGHTEN(int64_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_DARKEN:
					BLEND_DARKEN(int64_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_DARKEN(int64_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_DARKEN(int64_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_DARKEN(int64_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_OVERLAY:
					BLEND_OVERLAY(int64_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_OVERLAY(int64_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_OVERLAY(int64_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_OVERLAY(int64_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_SOFT_LIGHT:
					BLEND_SOFT_LIGHT(int64_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_SOFT_LIGHT(int64_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_SOFT_LIGHT(int64_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_SOFT_LIGHT(int64_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_COLOR_DODGE:
					BLEND_COLOR_DODGE(int64_t, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_COLOR_DODGE(int64_t, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_COLOR_DODGE(int64_t, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_COLOR_DODGE(int64_t, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_APPEND:
					appendOnNonString();
					break;
			}
			break;
		case RUVM_ATTRIB_V4_F32:
			switch (blendConfig.blend) {
				case RUVM_BLEND_REPLACE:
					BLEND_REPLACE(float, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_REPLACE(float, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_REPLACE(float, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_REPLACE(float, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_MULTIPLY:
					BLEND_MULTIPLY(float, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_MULTIPLY(float, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_MULTIPLY(float, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_MULTIPLY(float, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_DIVIDE:
					BLEND_DIVIDE(float, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_DIVIDE(float, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_DIVIDE(float, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_DIVIDE(float, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_ADD:
					BLEND_ADD(float, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_ADD(float, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_ADD(float, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_ADD(float, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_SUBTRACT:
					BLEND_SUBTRACT(float, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_SUBTRACT(float, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_SUBTRACT(float, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_SUBTRACT(float, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_ADD_SUB:
					BLEND_ADD_SUB(float, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_ADD_SUB(float, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_ADD_SUB(float, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_ADD_SUB(float, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_LIGHTEN:
					BLEND_LIGHTEN(float, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_LIGHTEN(float, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_LIGHTEN(float, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_LIGHTEN(float, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_DARKEN:
					BLEND_DARKEN(float, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_DARKEN(float, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_DARKEN(float, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_DARKEN(float, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_OVERLAY:
					BLEND_OVERLAY(float, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_OVERLAY(float, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_OVERLAY(float, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_OVERLAY(float, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_SOFT_LIGHT:
					BLEND_SOFT_LIGHT(float, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_SOFT_LIGHT(float, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_SOFT_LIGHT(float, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_SOFT_LIGHT(float, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_COLOR_DODGE:
					BLEND_COLOR_DODGE(float, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_COLOR_DODGE(float, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_COLOR_DODGE(float, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_COLOR_DODGE(float, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_APPEND:
					appendOnNonString();
					break;
			}
			break;
		case RUVM_ATTRIB_V4_F64:
			switch (blendConfig.blend) {
				case RUVM_BLEND_REPLACE:
					BLEND_REPLACE(double, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_REPLACE(double, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_REPLACE(double, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_REPLACE(double, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_MULTIPLY:
					BLEND_MULTIPLY(double, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_MULTIPLY(double, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_MULTIPLY(double, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_MULTIPLY(double, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_DIVIDE:
					BLEND_DIVIDE(double, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_DIVIDE(double, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_DIVIDE(double, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_DIVIDE(double, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_ADD:
					BLEND_ADD(double, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_ADD(double, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_ADD(double, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_ADD(double, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_SUBTRACT:
					BLEND_SUBTRACT(double, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_SUBTRACT(double, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_SUBTRACT(double, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_SUBTRACT(double, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_ADD_SUB:
					BLEND_ADD_SUB(double, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_ADD_SUB(double, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_ADD_SUB(double, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_ADD_SUB(double, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_LIGHTEN:
					BLEND_LIGHTEN(double, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_LIGHTEN(double, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_LIGHTEN(double, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_LIGHTEN(double, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_DARKEN:
					BLEND_DARKEN(double, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_DARKEN(double, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_DARKEN(double, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_DARKEN(double, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_OVERLAY:
					BLEND_OVERLAY(double, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_OVERLAY(double, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_OVERLAY(double, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_OVERLAY(double, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_SOFT_LIGHT:
					BLEND_SOFT_LIGHT(double, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_SOFT_LIGHT(double, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_SOFT_LIGHT(double, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_SOFT_LIGHT(double, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_COLOR_DODGE:
					BLEND_COLOR_DODGE(double, pD, iD, pA, iA, pB, iB, 4, 0);
					BLEND_COLOR_DODGE(double, pD, iD, pA, iA, pB, iB, 4, 1);
					BLEND_COLOR_DODGE(double, pD, iD, pA, iA, pB, iB, 4, 2);
					BLEND_COLOR_DODGE(double, pD, iD, pA, iA, pB, iB, 4, 3);
					break;
				case RUVM_BLEND_APPEND:
					appendOnNonString();
					break;
			}
			break;
		case RUVM_ATTRIB_STRING:
			//TODO add string append
			break;
	}
}
