#include <math.h>
#include <float.h>

#include <Types.h>

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

Vec3 vec3Lerp(Vec3 a, Vec3 b, float alpha) {
	float alphaInverse = 1.0f - alpha;
	Vec3 c;
	c.x = a.x * alphaInverse + b.x * alpha;
	c.y = a.y * alphaInverse + b.y * alpha;
	c.z = a.z * alphaInverse + b.z * alpha;
	return c;
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
		int32_t vertIndex = pMesh->pLoops[face.start + i];
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

void getFaceBounds(FaceBounds *pBounds, Vec2 *pUvs, FaceInfo faceInfo) {
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
                          RuvmMesh *pMesh, EdgeVerts *pEdgeVerts) {
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
