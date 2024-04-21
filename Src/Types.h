#pragma once
#include <stdint.h>
#include <Mesh.h>
#include <RUVM.h>

#define FLOAT_EQUAL_MARGIN .0001f

typedef RuvmVec2 Vec2;
typedef RuvmVec3 Vec3;

typedef struct {
	double x;
	double y;
	double z;
} Vec364;

typedef struct {
	double x;
	double y;
} Vec264;

typedef struct {
	int8_t x;
	int8_t y;
	int8_t z;
} iVec38;

typedef struct {
	int16_t x;
	int16_t y;
	int16_t z;
} iVec316;

typedef struct {
	int32_t x;
	int32_t y;
	int32_t z;
} iVec3;

typedef struct {
	int64_t x;
	int64_t y;
	int64_t z;
} iVec364;

typedef struct {
	int8_t x;
	int8_t y;
} iVec28;

typedef struct {
	int16_t x;
	int16_t y;
} iVec216;

typedef struct {
	int64_t x;
	int64_t y;
} iVec264;

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
	float d[2][3];
} Mat2x3;

typedef struct {
	float d[2][2];
} Mat2x2;

typedef struct {
	float d[3][3];
} Mat3x3;

typedef struct {
	int32_t start;
	int32_t end;
	int32_t size;
	int32_t index;
} FaceInfo;

typedef struct {
	int32_t triCount;
	int32_t *pTris;
	int32_t loopCount;
	int32_t *pLoops;
} FaceTriangulated;

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
	int32_t baseFace;
	int8_t baseLoops[8];
	int8_t onLine;
	int8_t isRuvm;
	int8_t temp;
	int8_t seam;
	int8_t seams;
	int16_t fSorts[8];
} BoundaryVert;

typedef struct EdgeTable {
	struct EdgeTable *pNext;
	int32_t ruvmFace;
	int32_t ruvmEdge;
	int32_t vert;
	int32_t tile;
	int32_t loops;
	int32_t baseEdge;
	int32_t baseVert;
	int8_t keepBaseLoop;
	int32_t job;
	int32_t loopIndex;
} EdgeTable;

typedef struct RealEdgeTable {
	struct RealEdgeTable *pNext;
	int32_t edge;
	int32_t inEdge;
	int32_t mapFace;
	int32_t valid;
} SeamEdgeTable;

typedef struct VertSeamTable {
	struct VertSeamTable *pNext;
	int32_t seams;
} VertSeamTable;

typedef struct {
	int32_t verts[2];
} EdgeVerts;

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
	iVec2 min, max;
	Vec2 fMin, fMax;
	Vec2 fMinSmall, fMaxSmall;
} FaceBounds;

typedef struct {
	Cell **pCells;
	int8_t *pCellType;
	int32_t cellSize;
	int32_t faceSize;
	FaceBounds faceBounds;
} FaceCellsInfo;

typedef struct VertAdj{
	struct VertAdj *pNext;
	int32_t vert;
	int32_t loopSize;
	int32_t baseFace;
	int32_t mapVert;
} VertAdj;

typedef struct MeshBufEdgeTable  {
	struct MeshBufEdgeTable *pNext;
	int32_t edge;
	int32_t refFace;
	int32_t refEdge;
	int32_t loopCount;
} MeshBufEdgeTable;

typedef struct {
	VertAdj *pVertTable;
	BoundaryVert *pBoundaryTable;
	EdgeTable *pEdgeTable;
} MergeTables;

typedef struct {
	Vec3 loop;
	int32_t index;
	int32_t sort;
	int32_t fSort;
	int32_t baseLoop;
	int32_t refEdge;
	Vec2 uv;
	int8_t onLine;
	int32_t onLineBase;
	int8_t isBaseLoop;
	int8_t seam;
	int8_t preserve;
	Vec3 normal;
	Vec3 bc; //bearycentric coords
	int32_t triLoops[3];
} LoopBuffer;

typedef struct OnLineTable {
	struct OnLineTable *pNext;
	int32_t baseEdgeOrLoop;
	int32_t ruvmVert;
	int32_t outVert;
	int32_t type;
} OnLineTable;

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
	int32_t totalBoundaryEdges;
	BufMesh bufMesh;
	int32_t bufferSize;
	int32_t totalVerts;
	int32_t totalLoops;
	int32_t totalEdges;
	int32_t totalFaces;
	int32_t loopBufferSize;
	int32_t *pBoundaryVerts;
	int32_t boundaryVertSize;
	Mesh mesh;
	int32_t vertBase;
	int64_t averageRuvmFacesPerFace;
	EdgeVerts *pEdgeVerts;
	RuvmCommonAttribList *pCommonAttribList;
} ThreadArg;

typedef struct {
	RuvmContext pContext;
	RuvmMap pMap;
	int32_t bufferSize;
	int32_t loopBufferSize;
	void *pMutex;
	int32_t id;
	int32_t *pJobsCompleted;
	BoundaryDir *pBoundaryBuffer;
	int32_t boundaryBufferSize;
	int32_t averageVertAdjDepth;
	Mesh mesh;
	int64_t averageRuvmFacesPerFace;
	BufMesh bufMesh;
	int32_t vertBase;
	int32_t edgeBase;
	int32_t totalBoundaryFaces;
	int32_t totalBoundaryEdges;
	int32_t totalVerts;
	int32_t totalLoops;
	int32_t totalEdges;
	int32_t totalFaces;
	int8_t *pInVertTable;
	int8_t *pVertSeamTable;
	EdgeVerts *pEdgeVerts;
	RuvmCommonAttribList *pCommonAttribList;
} SendOffArgs;

typedef struct {
	int32_t index;
	int32_t edgeIndex;
	int32_t edgeIndexNext;
	int8_t edgeIsSeam;
	int8_t edgeNextIsSeam;
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
	int32_t edgeIndex;
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

Vec364 vec364MultiplyScalar(Vec364 a, double b);
void vec364DivideEqualScalar(Vec364 *pA, double b);
Vec364 vec364DivideScalar(Vec364 a, double b);
void vec364AddEqual(Vec364 *pA, Vec364 b);

Vec264 vec264MultiplyScalar(Vec264 a, double b);
void vec264DivideEqualScalar(Vec264 *pA, double b);
Vec264 vec264DivideScalar(Vec264 a, double b);
void vec264AddEqual(Vec264 *pA, Vec264 b);

iVec38 veci38MultiplyScalar(iVec38 a, float b);
void vec38DivideEqualScalar(iVec38 *pA, float b);
iVec38 veci38DivideScalar(iVec38 a, float b);
void veci38AddEqual(iVec38 *pA, iVec38 b);

iVec316 veci316MultiplyScalar(iVec316 a, float b);
void veci316DivideEqualScalar(iVec316 *pA, float b);
iVec316 veci316DivideScalar(iVec316 a, float b);
void veci316AddEqual(iVec316 *pA, iVec316 b);

iVec3 veci3MultiplyScalar(iVec3 a, float b);
void veci3DivideEqualScalar(iVec3 *pA, float b);
iVec3 veci3DivideScalar(iVec3 a, float b);
void veci3AddEqual(iVec3 *pA, iVec3 b);

iVec364 veci364MultiplyScalar(iVec364 a, float b);
void veci364DivideEqualScalar(iVec364 *pA, float b);
iVec364 veci364DivideScalar(iVec364 a, float b);
void veci364AddEqual(iVec364 *pA, iVec364 b);

iVec28 veci28MultiplyScalar(iVec28 a, float b);
void veci28DivideEqualScalar(iVec28 *pA, float b);
iVec28 veci28DivideScalar(iVec28 a, float b);
void veci28AddEqual(iVec28 *pA, iVec28 b);

iVec216 veci216MultiplyScalar(iVec216 a, float b);
void veci216DivideEqualScalar(iVec216 *pA, float b);
iVec216 veci216DivideScalar(iVec216 a, float b);
void veci216AddEqual(iVec216 *pA, iVec216 b);

iVec2 veci2MultiplyScalar(iVec2 a, float b);
void veci2DivideEqualScalar(iVec2 *pA, float b);
iVec2 veci2DivideScalar(iVec2 a, float b);
void veci2AddEqual(iVec2 *pA, iVec2 b);

iVec264 veci264MultiplyScalar(iVec264 a, float b);
void veci264DivideEqualScalar(iVec264 *pA, float b);
iVec264 veci264DivideScalar(iVec264 a, float b);
void veci264AddEqual(iVec264 *pA, iVec264 b);

Vec3 vec3MultiplyScalar(Vec3 a, float b);
void vec3DivideEqualScalar(Vec3 *pA, float b);
Vec3 vec3DivideScalar(Vec3 a, float b);
Vec3 vec3SubtractScalar(Vec3 a, float b);
Vec3 vec3AddScalar(Vec3 a, float b);
Vec3 vec3Add(Vec3 a, Vec3 b);
Vec3 vec3Subtract(Vec3 a, Vec3 b);
void vec3AddEqual(Vec3 *pA, Vec3 b);
int32_t vec3GreaterThan(Vec3 a, Vec3 b);
int32_t vec3LessThan(Vec3 a, Vec3 b);
int32_t vec3AproxEqual(Vec3 a, Vec3 b);
Vec3 vec3Lerp(Vec3 a, Vec3 b, float alpha);
Vec3 vec3Cross(Vec3 a, Vec3 b);
Vec3 vec3UnitFromPoints(Vec3 a, Vec3 b);
Vec3 vec3MultiplyMat3x3(Vec3 a, Mat3x3 *pB);
Vec3 vec3Normalize(Vec3 a);

Vec2 vec2Abs(Vec2 a);
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
int32_t vec2NotEqual(Vec2 a, Vec2 b);
int32_t vec2Equal(Vec2 a, Vec2 b);
int32_t vec2AproxEqual(Vec2 a, Vec2 b);
int32_t vec2WindingCompare(Vec2 a, Vec2 b, Vec2 centre, int32_t fallBack);

float customFloor(float a);
iVec2 vec2FloorAssign(Vec2 *pA);
Vec3 cartesianToBarycentric(Vec2 *pTri, Vec3 *pPoint);
Vec3 barycentricToCartesian(Vec3 *pTri, Vec3 *pPoint);

int32_t checkFaceIsInBounds(Vec2 min, Vec2 max, FaceInfo face, Mesh *pMesh);
void getFaceBounds(FaceBounds *pBounds, RuvmAttrib *pUvs, FaceInfo faceInfo);
int32_t checkIfEdgeIsSeam(int32_t edgeIndex, FaceInfo face, int32_t loop,
                          Mesh *pMesh, EdgeVerts *pEdgeVerts);

Mat2x2 mat2x2Adjugate(Mat2x2 a);
float mat2x2Determinate(Mat2x2 a);
void mat2x2MultiplyEqualScalar(Mat2x2 *pA, float b);
Mat2x2 mat2x2Invert(Mat2x2 a);

Mat2x3 mat2x2MultiplyMat2x3(Mat2x2 a, Mat2x3 b);

uint32_t ruvmFnvHash(uint8_t *value, int32_t valueSize, uint32_t size);


int32_t checkIfEdgeIsPreserve(Mesh* pMesh, int32_t edge);
FaceTriangulated triangulateFace(RuvmAllocator alloc, FaceInfo baseFace, Mesh *pMesh);

void blendAttribReplace_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                           RuvmAttrib *pB, int32_t iB);
void blendAttribReplace_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribReplace_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribReplace_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribReplace_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribReplace_F64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribReplace_V2_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribReplace_V2_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribReplace_V2_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribReplace_V2_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribReplace_V2_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribReplace_V2_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribReplace_V3_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribReplace_V3_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribReplace_V3_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribReplace_V3_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribReplace_V3_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribReplace_V3_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribReplace_V4_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribReplace_V4_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribReplace_V4_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribReplace_V4_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribReplace_V4_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribReplace_V4_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribMultiply_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                           RuvmAttrib *pB, int32_t iB);
void blendAttribMultiply_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribMultiply_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribMultiply_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribMultiply_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribMultiply_F64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribMultiply_V2_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribMultiply_V2_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribMultiply_V2_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribMultiply_V2_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribMultiply_V2_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribMultiply_V2_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribMultiply_V3_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribMultiply_V3_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribMultiply_V3_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribMultiply_V3_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribMultiply_V3_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribMultiply_V3_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribMultiply_V4_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribMultiply_V4_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribMultiply_V4_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribMultiply_V4_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribMultiply_V4_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribMultiply_V4_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribDivide_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                           RuvmAttrib *pB, int32_t iB);
void blendAttribDivide_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribDivide_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribDivide_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribDivide_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribDivide_F64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribDivide_V2_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribDivide_V2_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribDivide_V2_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribDivide_V2_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribDivide_V2_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribDivide_V2_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribDivide_V3_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribDivide_V3_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribDivide_V3_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribDivide_V3_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribDivide_V3_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribDivide_V3_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribDivide_V4_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribDivide_V4_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribDivide_V4_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribDivide_V4_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribDivide_V4_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribDivide_V4_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                           RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_F64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V2_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V2_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V2_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V2_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V2_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V2_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V3_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V3_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V3_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V3_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V3_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V3_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V4_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V4_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V4_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V4_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V4_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V4_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                           RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_F64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V2_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V2_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V2_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V2_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V2_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V2_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V3_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V3_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V3_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V3_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V3_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V3_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V4_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V4_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V4_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V4_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V4_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAdd_V4_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribSubtract_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                           RuvmAttrib *pB, int32_t iB);
void blendAttribSubtract_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribSubtract_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribSubtract_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribSubtract_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribSubtract_F64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribSubtract_V2_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribSubtract_V2_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribSubtract_V2_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribSubtract_V2_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribSubtract_V2_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribSubtract_V2_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribSubtract_V3_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribSubtract_V3_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribSubtract_V3_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribSubtract_V3_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribSubtract_V3_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribSubtract_V3_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribSubtract_V4_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribSubtract_V4_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribSubtract_V4_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribSubtract_V4_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribSubtract_V4_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribSubtract_V4_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAddSub_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                           RuvmAttrib *pB, int32_t iB);
void blendAttribAddSub_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribAddSub_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribAddSub_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribAddSub_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribAddSub_F64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribAddSub_V2_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribAddSub_V2_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAddSub_V2_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAddSub_V2_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAddSub_V2_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAddSub_V2_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAddSub_V3_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribAddSub_V3_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAddSub_V3_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAddSub_V3_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAddSub_V3_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAddSub_V3_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAddSub_V4_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribAddSub_V4_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAddSub_V4_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAddSub_V4_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAddSub_V4_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAddSub_V4_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribLighten_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                           RuvmAttrib *pB, int32_t iB);
void blendAttribLighten_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribLighten_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribLighten_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribLighten_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribLighten_F64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribLighten_V2_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribLighten_V2_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribLighten_V2_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribLighten_V2_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribLighten_V2_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribLighten_V2_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribLighten_V3_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribLighten_V3_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribLighten_V3_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribLighten_V3_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribLighten_V3_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribLighten_V3_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribLighten_V4_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribLighten_V4_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribLighten_V4_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribLighten_V4_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribLighten_V4_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribLighten_V4_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribDarken_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                           RuvmAttrib *pB, int32_t iB);
void blendAttribDarken_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribDarken_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribDarken_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribDarken_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribDarken_F64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribDarken_V2_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribDarken_V2_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribDarken_V2_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribDarken_V2_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribDarken_V2_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribDarken_V2_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribDarken_V3_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribDarken_V3_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribDarken_V3_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribDarken_V3_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribDarken_V3_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribDarken_V3_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribDarken_V4_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribDarken_V4_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribDarken_V4_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribDarken_V4_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribDarken_V4_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribDarken_V4_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribOverlay_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                           RuvmAttrib *pB, int32_t iB);
void blendAttribOverlay_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribOverlay_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribOverlay_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribOverlay_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribOverlay_F64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribOverlay_V2_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribOverlay_V2_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribOverlay_V2_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribOverlay_V2_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribOverlay_V2_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribOverlay_V2_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribOverlay_V3_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribOverlay_V3_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribOverlay_V3_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribOverlay_V3_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribOverlay_V3_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribOverlay_V3_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribOverlay_V4_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribOverlay_V4_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribOverlay_V4_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribOverlay_V4_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribOverlay_V4_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribOverlay_V4_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribSoftLight_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                           RuvmAttrib *pB, int32_t iB);
void blendAttribSoftLight_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribSoftLight_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribSoftLight_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribSoftLight_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribSoftLight_F64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribSoftLight_V2_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribSoftLight_V2_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribSoftLight_V2_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribSoftLight_V2_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribSoftLight_V2_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribSoftLight_V2_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribSoftLight_V3_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribSoftLight_V3_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribSoftLight_V3_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribSoftLight_V3_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribSoftLight_V3_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribSoftLight_V3_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribSoftLight_V4_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribSoftLight_V4_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribSoftLight_V4_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribSoftLight_V4_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribSoftLight_V4_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribSoftLight_V4_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribColorDodge_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                           RuvmAttrib *pB, int32_t iB);
void blendAttribColorDodge_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribColorDodge_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribColorDodge_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribColorDodge_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribColorDodge_F64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                            RuvmAttrib *pB, int32_t iB);
void blendAttribColorDodge_V2_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribColorDodge_V2_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribColorDodge_V2_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribColorDodge_V2_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribColorDodge_V2_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribColorDodge_V2_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribColorDodge_V3_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribColorDodge_V3_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribColorDodge_V3_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribColorDodge_V3_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribColorDodge_V3_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribColorDodge_V3_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribColorDodge_V4_I8(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                              RuvmAttrib *pB, int32_t iB);
void blendAttribColorDodge_V4_I16(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribColorDodge_V4_I32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribColorDodge_V4_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribColorDodge_V4_F32(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribColorDodge_V4_I64(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                               RuvmAttrib *pB, int32_t iB);
void blendAttribAppend_Str(RuvmAttrib *pD, int32_t iD, RuvmAttrib *pA, int32_t iA,
                           RuvmAttrib *pB, int32_t iB);

int32_t getAttribSize(RuvmAttribType type);
RuvmAttrib *getAttrib(char *pName, RuvmAttrib *pAttribs, int32_t attribCount);
Vec3 *attribAsV3(RuvmAttrib *pAttrib, int32_t index);
Vec2 *attribAsV2(RuvmAttrib *pAttrib, int32_t index);
int32_t *attribAsI32(RuvmAttrib *pAttrib, int32_t index);
void *attribAsVoid(RuvmAttrib *pAttrib, int32_t index);
int32_t copyAttrib(RuvmAttrib *pDest, int32_t iDest,
                   RuvmAttrib *pSrc, int32_t iSrc);
void copyAllAttribs(RuvmAttrib *pDest, int32_t iDest,
                    RuvmAttrib *pSrc, int32_t iSrc,
                    int32_t attribCount);
RuvmTypeDefault *getTypeDefaultConfig(RuvmTypeDefaultConfig *pConfig,
                                      RuvmAttribType type);
RuvmCommonAttrib *getCommonAttrib(RuvmCommonAttrib *pAttribs, int32_t attribCount,
                                  char *pName);
void interpolateAttrib(RuvmAttrib *pDest, int32_t iDest, RuvmAttrib *pSrc,
                       int32_t iSrcA, int32_t iSrcB, int32_t iSrcC, Vec3 bc);
void blendAttribs(RuvmAttrib *pDest, int32_t iDest, RuvmAttrib *pA, int32_t iA,
                  RuvmAttrib *pB, int32_t iB, RuvmBlendConfig blendConfig);

#define V3_F64_MULS 	,364MultiplyScalar,
#define V3_F64_DIVEQLS 	,364DivideEqualScalar,
#define V3_F64_DIVS 	,364DivideScalar,
#define V3_F64_ADDEQL 	,364AddEqual,

#define V2_F64_MUL  	,264MultiplyScalar,
#define V2_F64_DIVEQLS	,264DivideEqualScalar,
#define V2_F64_DIVS 	,264DivideScalar,
#define V2_F64_ADDEQL 	,264AddEqual,

#define V3_I8_MULS		,i38MultiplyScalar,
#define V3_I8_DIVEQLS	,i38DivideEqualScalar,
#define V3_I8_DIVS		,i38DivideScalar,
#define V3_I8_ADDEQL	,i38AddEqual,

#define V3_I16_MULS		,i316MultiplyScalar,
#define V3_I16_DIVEQLS	,i316DivideEqualScalar,
#define V3_I16_DIVS		,i316DivideScalar,
#define V3_I16_ADDEQL	,i316AddEqual,

#define V3_I32_MULS		,i3MultiplyScalar,
#define V3_I32_DIVEQLS	,i3DivideEqualScalar,
#define V3_I32_DIVS		,i3DivideScalar,
#define V3_I32_ADDEQL	,i3AddEqual,

#define V3_I64_MULS		,i364MultiplyScalar,
#define V3_I64_DIVEQLS	,i364DivideEqualScalar,
#define V3_I64_DIVS		,i364DivideScalar,
#define V3_I64_ADDEQL	,i364AddEqual,

#define V3_I64_MULS		,i364MultiplyScalar,
#define V3_I64_DIVEQLS	,i364DivideEqualScalar,
#define V3_I64_DIVS		,i364DivideScalar,
#define V3_I64_ADDEQL	,i364AddEqual,

#define V2_I8_MULS		,i28MultiplyScalar,
#define V2_I8_DIVEQLS	,i28DivideEqualScalar,
#define V2_I8_DIVS		,i28DivideScalar,
#define V2_I8_ADDEQL	,i28AddEqual,

#define V2_I16_MULS		,i216MultiplyScalar,
#define V2_I16_DIVEQLS	,i216DivideEqualScalar,
#define V2_I16_DIVS		,i216DivideScalar,
#define V2_I16_ADDEQL	,i216AddEqual,

#define V2_I32_MULS		,i2MultiplyScalar,
#define V2_I32_DIVEQLS	,i2DivideEqualScalar,
#define V2_I32_DIVS		,i2DivideScalar,
#define V2_I32_ADDEQL	,i2AddEqual,

#define V2_I64_MULS		,i264MultiplyScalar,
#define V2_I64_DIVEQLS	,i264DivideEqualScalar,
#define V2_I64_DIVS		,i264DivideScalar,
#define V2_I64_ADDEQL	,i264AddEqual,

#define V3MULS ,3MultiplyScalar,
#define V3DIVEQLS ,3DivideEqualScalar,
#define V3DIVS ,3DivideScalar,
#define V3SUB ,3Subtract,
#define V3SUBS ,3SubtractScalar,
#define V3ADDS ,3AddScalar,
#define V3ADD ,3Add,
#define V3ADDEQL ,3AddEqual,
#define V3GREAT ,3GreaterThan,
#define V3LESS ,3LessThan,
#define V3CROSS ,3Cross,
#define V3APROXEQL ,3AproxEqual,
#define V3MULM3X3 ,3MultiplyMat3x3,

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
#define V2NOTEQL ,2NotEqual,
#define V2EQL ,2Equal,
#define V2APROXEQL ,2AproxEqual,

#define INFIX(a,o,b) vec##o((a),(b))
#define _(a) INFIX(a)
