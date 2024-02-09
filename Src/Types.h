#pragma once
#include <stdint.h>
#include <Mesh.h>
#include <RUVM.h>

typedef RuvmVec2 Vec2;
typedef RuvmVec3 Vec3;

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

typedef struct {
	int32_t start;
	int32_t end;
	int32_t size;
	int32_t index;
} FaceInfo;

typedef struct BoundaryVert{
	struct BoundaryVert *pNext;
	int32_t faceIndex;
	int32_t tile;
	int8_t firstVert;
	int8_t lastVert;
	int32_t face;
	int8_t type;
	int32_t job;
	int8_t hasPreservedEdge;
	int32_t baseLoop;
	int8_t baseLoops[8];
} BoundaryVert;

typedef struct EdgeTable {
	struct EdgeTable *pNext;
	int32_t ruvmVert;
	int32_t ruvmVertNext;
	int32_t vert;
	int32_t tile;
	int32_t loops;
	int32_t baseEdge;
	int32_t baseVert;
} EdgeTable;

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
	int32_t linkEdgeSize;
	int32_t *pLinkEdges;
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
	int8_t isBaseLoop;
} LoopBuffer;

typedef struct {
	LoopBuffer buf[12];
	int32_t size;
} LoopBufferWrap;

typedef struct BoundaryDir {
	struct BoundaryDir *pNext;
	BoundaryVert *pEntry;
} BoundaryDir;

typedef struct {
	int32_t maxLoopSize;
	RuvmAllocator alloc;
	RuvmMap pMap;
	int32_t id;
	int32_t averageVertAdjDepth;
	int32_t *pBoundaryFaceStart;
	BoundaryDir *pBoundaryBuffer;
	int32_t boundaryBufferSize;
	int32_t totalBoundaryFaces;
	WorkMesh localMesh;
	int32_t bufferSize;
	int32_t totalVerts;
	int32_t totalLoops;
	int32_t totalFaces;
	int32_t loopBufferSize;
	int32_t *pBoundaryVerts;
	int32_t boundaryVertSize;
	Mesh mesh;
	int32_t vertBase;
	int64_t averageRuvmFacesPerFace;
	int8_t *pInVertTable;
} ThreadArg;

typedef struct {
	RuvmContext pContext;
	RuvmMap pMap;
	void *pMutex;
	int32_t id;
	int32_t *pJobsCompleted;
	BoundaryDir *pBoundaryBuffer;
	int32_t boundaryBufferSize;
	int32_t averageVertAdjDepth;
	Mesh mesh;
	int64_t averageRuvmFacesPerFace;
	WorkMesh localMesh;
	int32_t vertBase;
	int32_t totalBoundaryFaces;
	int32_t totalVerts;
	int32_t totalLoops;
	int32_t totalFaces;
	int8_t *pInVertTable;
} SendOffArgs;

typedef struct {
	iVec2 min, max;
	Vec2 fMin, fMax;
	Vec2 fMinSmall, fMaxSmall;
} FaceBounds;

typedef struct {
	int32_t index;
	int32_t indexNext;
	int8_t localIndex;
	int8_t localIndexNext;
	Vec2 vert;
	Vec2 vertNext;
	Vec2 dir;
	Vec2 dirBack;
} LoopInfo;

typedef struct {
	int32_t loopStart;
	int32_t boundaryLoopStart;
	int32_t firstRuvmVert, lastRuvmVert;
	int32_t ruvmLoops;
	int32_t vertIndex;
	int32_t loopIndex;
} AddClippedFaceVars;

typedef struct {
	uint64_t timeSpent[3];
	int32_t maxDepth;
} DebugAndPerfVars;

typedef struct {
	Vec2 uv[3];
	Vec3 xyz[3];
	Vec3 *pNormals;
} BaseTriVerts;

Vec3 vec3MultiplyScalar(Vec3 a, float b);
void vec3DivideEqualScalar(Vec3 *pA, float b);
Vec3 vec3SubtractScalar(Vec3 a, float b);
Vec3 vec3AddScalar(Vec3 a, float b);
Vec3 vec3Add(Vec3 a, Vec3 b);
Vec3 vec3Subtract(Vec3 a, Vec3 b);
void vec3AddEqual(Vec3 *pA, Vec3 b);
int32_t vec3GreaterThan(Vec3 a, Vec3 b);
int32_t vec3LessThan(Vec3 a, Vec3 b);

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
void vec2AddEqualScalar(Vec2 *pA, float b);
void vec2MultiplyEqualScalar(Vec2 *pA, float b);
Vec2 vec2MultiplyScalar(Vec2 a, float b);
float vec2Dot(Vec2 a, Vec2 b);
Vec2 vec2Cross(Vec2 a);
Vec2 vec2ModScalar(Vec2 a, float b);
void vec2ModEqualScalar(Vec2 *pA, float b);
int32_t vec2GreaterThan(Vec2 a, Vec2 b);
int32_t vec2GreaterThanScalar(Vec2 a, float b);
int32_t vec2GreaterThanEqualTo(Vec2 a, Vec2 b);
int32_t vec2LessThan(Vec2 a, Vec2 b);
int32_t vec2LessThanScalar(Vec2 a, float b);
int32_t vec2LessThanEqualTo(Vec2 a, Vec2 b);
float customFloor(float a);
iVec2 vec2FloorAssign(Vec2 *pA);
Vec3 cartesianToBarycentric(Vec2 *pTri, Vec3 *pPoint);
Vec3 barycentricToCartesian(Vec3 *pTri, Vec3 *pPoint);

int32_t checkFaceIsInBounds(Vec2 min, Vec2 max, FaceInfo face, Mesh *pMesh);
void getFaceBounds(FaceBounds *pBounds, Vec2 *pUvs, FaceInfo faceInfo);

uint32_t ruvmFnvHash(uint8_t *value, int32_t valueSize, uint32_t size);

#define V3MULS ,3MultiplyScalar,
#define V3DIVEQLS ,3DivideEqualScalar,
#define V3SUB ,3Subtract,
#define V3SUBS ,3SubtractScalar,
#define V3ADDS ,3AddScalar,
#define V3ADD ,3Add,
#define V3ADDEQL ,3AddEqual,
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
#define V2ADDEQLS ,2AddEqualScalar,
#define V2MULSEQL ,2MultiplyEqualScalar,
#define V2MULS ,2MultiplyScalar,
#define V2DOT ,2Dot,
#define V2MODS ,2ModScalar,
#define V2MODEQLS ,2ModEqualScalar,
#define V2GREAT ,2GreaterThan,
#define V2GREATS ,2GreaterThanScalar,
#define V2GREATEQL ,2GreaterThanEqualTo,
#define V2LESS ,2LessThan,
#define V2LESSS ,2LessThanScalar,
#define V2LESSEQL ,2LessThanEqualTo,
#define INFIX(a,o,b) vec##o((a),(b))
#define _(a) INFIX(a)
