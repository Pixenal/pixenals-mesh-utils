#pragma once
#include <stdint.h>

typedef struct {
	float x;
	float y;
} Vec2;

typedef struct {
	float x;
	float y;
	float z;
} Vec3;

typedef struct {
	Vec3 a;
	Vec3 b;
	Vec3 c;
} TriXyz;

typedef struct {
	Vec2 a;
	Vec2 b;
	Vec2 c;
} TriUv;

typedef struct {
	int32_t x;
	int32_t y;
	int32_t z;
	int32_t w;
} iVec4;

typedef struct {
	int32_t x;
	int32_t y;
} iVec2;

typedef struct BoundaryVert{
	struct BoundaryVert *pNext;
	int32_t faceIndex;
	int32_t tile;
	int32_t firstVert;
	int32_t face;
	int32_t valid;
} BoundaryVert;

typedef struct EdgeTable {
	struct EdgeTable *pNext;
	int32_t ruvmVert;
	int32_t vert;
	int32_t tile;
	int32_t loops;
} EdgeTable;

typedef struct {
	int32_t vertSize;
	int32_t boundaryVertSize;
	Vec3 *pVerts;
	int32_t loopSize;
	int32_t boundaryLoopSize;
	int32_t *pLoops;
	Vec3 *pNormals;
	int32_t faceSize;
	int32_t boundaryFaceSize;
	int32_t *pFaces;
	Vec2 *pUvs;
} MeshData;

typedef struct Cell {
	uint32_t localIndex;
	uint32_t initialized;
	struct Cell *pChildren;
	int32_t faceSize;
	int32_t *pFaces;
	int32_t edgeFaceSize;
	int32_t *pEdgeFaces;
	int32_t cellIndex;
	Vec2 boundsMin;
	Vec2 boundsMax;
} Cell;

typedef struct {
	Cell **pCells;
	int32_t *pCellType;
	int32_t cellSize;
	int32_t faceSize;
} FaceCellsInfo;

typedef struct VertAdj{
	struct VertAdj *pNext;
	int32_t vert;
	int32_t ruvmVert;
	int32_t loopSize;
} VertAdj;

typedef struct {
	VertAdj *pVertTable;
	BoundaryVert *pBoundaryTable;
	EdgeTable *pEdgeTable;
} MergeTables;

typedef struct {
	Vec3 loop;
	int32_t index;
	int32_t sort;
	int32_t baseLoop;
	Vec2 uv;
} LoopBuffer;

typedef struct {
	int32_t id;
	int32_t *pJobsCompleted;
	int32_t averageVertAdjDepth;
	int32_t *pBoundaryFaceStart;
	BoundaryVert *pBoundaryBuffer;
	MeshData localMesh;
	int32_t bufferSize;
	int32_t totalVerts;
	int32_t totalLoops;
	int32_t totalFaces;
	BoundaryVert **pFinalBoundary;
	int32_t totalBoundaryFaces;
	int32_t loopBufferSize;
	int32_t *pBoundaryVerts;
	int32_t boundaryVertSize;
	MeshData mesh;
} ThreadArg;

Vec3 vec3SubtractScalar(Vec3 a, float b);
Vec3 vec3AddScalar(Vec3 a, float b);
int32_t vec3GreaterThan(Vec3 a, Vec3 b);
int32_t vec3LessThan(Vec3 a, Vec3 b);
int32_t vec3ApproxEqual(Vec3 a, Vec3 b);

Vec2 vec2Multiply(Vec2 a, Vec2 b);
void vec2MultiplyEqual(Vec2 *pA, Vec2 b);
Vec2 vec2DivideScalar(Vec2 a, float b);
void vec2DivideEqualScalar(Vec2 *pA, float b);
Vec2 vec2Subtract(Vec2 a, Vec2 b);
void vec2SubtractEqual(Vec2 *pA, Vec2 b);
Vec2 vec2SubtractScalar(Vec2 a, float b);
Vec2 vec2Add(Vec2 a, Vec2 b);
Vec2 vec2AddScalar(Vec2 a, float b);
void vec2AddEqual(Vec2 *pA, Vec2 b);
void vec2MultiplyEqualScalar(Vec2 *pA, float b);
Vec2 vec2MultiplyScalar(Vec2 a, float b);
float vec2Dot(Vec2 a, Vec2 b);
Vec2 vec2Cross(Vec2 a);
Vec2 vec2ModScalar(Vec2 a, float b);
void vec2ModEqualScalar(Vec2 *pA, float b);
int32_t vec2GreaterThan(Vec2 a, Vec2 b);
int32_t vec2LessThan(Vec2 a, Vec2 b);
int32_t vec2LessThanEqualTo(Vec2 a, Vec2 b);
Vec3 cartesianToBarycentric(Vec2 *pTri, Vec3 *pPoint);
Vec3 barycentricToCartesian(Vec3 *pTri, Vec3 *pPoint);


#define V3SUBS ,3SubtractScalar,
#define V3ADDS ,3AddScalar,
#define V3GREAT ,3GreaterThan,
#define V3LESS ,3LessThan,
#define V3APXEQL ,3ApproxEqual,

#define V2MUL ,2Multiply,
#define V2MULEQL ,2MultiplyEqual,
#define V2DIVS ,2DivideScalar,
#define V2DIVSEQL ,2DivideEqualScalar,
#define V2SUB ,2Subtract,
#define V2SUBEQL ,2SubtractEqual,
#define V2SUBS ,2SubtractScalar,
#define V2ADD ,2Add,
#define V2ADDS ,2AddScalar,
#define V2ADDEQL ,2AddEqual,
#define V2MULSEQL ,2MultiplyEqualScalar,
#define V2MULS ,2MultiplyScalar,
#define V2DOT ,2Dot,
#define V2MODS ,2ModScalar,
#define V2MODEQLS ,2ModEqualScalar,
#define V2GREAT ,2GreaterThan,
#define V2LESS ,2LessThan,
#define V2LESSEQL ,2LessThanEqualTo,
#define INFIX(a,o,b) vec##o((a),(b))
#define _(a) INFIX(a)
