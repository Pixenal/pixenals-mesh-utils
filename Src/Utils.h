#pragma once

#include <RuvmInternal.h>
#include <MathUtils.h>
#include <Mesh.h>

typedef struct {
	int32_t triCount;
	int32_t *pTris;
	int32_t loopCount;
	int32_t *pLoops;
} FaceTriangulated;

typedef struct {
	int32_t start;
	int32_t end;
	int32_t size;
	int32_t index;
} FaceRange;

typedef struct {
	V2_I32 min, max;
	V2_F32 fMin, fMax;
	V2_F32 fMinSmall, fMaxSmall;
} FaceBounds;

int32_t checkFaceIsInBounds(V2_F32 min, V2_F32 max, FaceRange face, Mesh *pMesh);
void getFaceBounds(FaceBounds *pBounds, V2_F32 *pUvs, FaceRange face);
int32_t checkIfEdgeIsSeam(int32_t edgeIndex, FaceRange face, int32_t loop,
                          Mesh *pMesh, EdgeVerts *pEdgeVerts);

uint32_t ruvmFnvHash(uint8_t *value, int32_t valueSize, uint32_t size);

int32_t checkIfEdgeIsPreserve(Mesh* pMesh, int32_t edge);
int32_t checkIfEdgeIsReceive(Mesh* pMesh, int32_t edge);
FaceTriangulated triangulateFace(RuvmAlloc alloc, FaceRange baseFace, void *pVerts,
                                 int32_t *pLoops, int32_t useUvs);
V3_F32 getBarycentricInFace(V2_F32 *pTriUvs, int8_t *pTriLoops,
                            int32_t loopCount, V2_F32 vert);
void waitForJobs(RuvmContext pContext, int32_t *pJobsCompleted, void *pMutex);
FaceRange getFaceRange(const RuvmMesh *pMesh, int32_t index, _Bool border);
void buildEdgeList(RuvmContext pContext, Mesh* pMesh);
_Bool isMeshInvalid(Mesh* pMesh);