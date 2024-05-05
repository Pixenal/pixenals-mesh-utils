#pragma once
#include <stdint.h>
#include <Mesh.h>
#include <MathUtils.h>

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
	uint64_t baseLoop : 22;
	uint64_t ruvmLoop : 33;
	uint64_t job : 8;
	uint64_t seam : 1;

	uint64_t tile : 24;
	uint64_t isRuvm : 11;
	uint64_t onLine : 11;
	uint64_t hasPreservedEdge : 1;

	int32_t faceIndex;
	int32_t face;
	int32_t baseFace;
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
	int32_t loop;
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
	V2_F32 boundsMin;
	V2_F32 boundsMax;
	int32_t linkEdgeSize;
	int32_t *pLinkEdges;
} Cell;

typedef struct {
	V2_I32 min, max;
	V2_F32 fMin, fMax;
	V2_F32 fMinSmall, fMaxSmall;
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
	V3_F32 loop;
	V3_F32 normal;
	V3_F32 bc; //barycentric coords
	V2_F32 uv;
	int8_t baseLoop;
	uint8_t ruvmLoop : 4;
	uint8_t onLine : 1;
	uint8_t isBaseLoop : 1;
	uint8_t preserve : 1;
	uint8_t isRuvm : 1;
	int8_t triLoops[3];
} LoopBuffer;

typedef struct OnLineTable {
	struct OnLineTable *pNext;
	int32_t baseEdgeOrLoop;
	int32_t ruvmVert;
	int32_t outVert;
	int32_t type;
} OnLineTable;

typedef struct {
	LoopBuffer buf[11];
	int32_t size;
} LoopBufferWrap;

typedef struct BoundaryDir {
	struct BoundaryDir *pNext;
	BoundaryVert *pEntry;
} BoundaryDir;

typedef struct {
	BoundaryDir *pTable;
	int32_t size;
} BoundaryTable;

typedef struct {
	VertAdj *pRuvmVertAdj;
	MeshBufEdgeTable *pEdgeTable;
	uint32_t vertAdjSize;
	int32_t edgeTableSize;
} MappingTables;

typedef struct {
	RuvmMap pMap;
	BufMesh bufMesh;
	Mesh mesh;
	RuvmAllocator alloc;
	MappingTables mTables;
	Mat3x3 tbn;
	RuvmCommonAttribList *pCommonAttribList;
	BoundaryTable boundsTable;
	EdgeVerts *pEdgeVerts;
	int8_t *pInVertTable;
	int32_t id;
	int32_t bufferSize;
	int32_t loopBufferSize;
} MappingJobVars;

typedef struct {
	int32_t vertBase;
	int32_t edgeBase;
} JobBases;

typedef struct {
	RuvmMap pMap;
	Mesh mesh;
	BufMesh bufMesh;
	RuvmContext pContext;
	BoundaryTable boundsTable;
	EdgeVerts *pEdgeVerts;
	RuvmCommonAttribList *pCommonAttribList;
	int32_t *pJobsCompleted;
	void *pMutex;
	int8_t *pInVertTable;
	int32_t id;
	int32_t totalBoundaryFaces;
	int32_t totalBoundaryEdges;
	int32_t totalVerts;
	int32_t totalLoops;
	int32_t totalEdges;
	int32_t totalFaces;
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
	V2_F32 vert;
	V2_F32 vertNext;
	V2_F32 dir;
	V2_F32 dirBack;
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
	V2_F32 uv[4];
	V3_F32 xyz[4];
} BaseTriVerts;
