#pragma once

#include <RuvmInternal.h>
#include <MathUtils.h>

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
} FaceInfo;

typedef struct {
	V2_I32 min, max;
	V2_F32 fMin, fMax;
	V2_F32 fMinSmall, fMaxSmall;
} FaceBounds;

int32_t checkFaceIsInBounds(V2_F32 min, V2_F32 max, FaceInfo face, Mesh *pMesh);
void getFaceBounds(FaceBounds *pBounds, V2_F32 *pUvs, FaceInfo faceInfo);
int32_t checkIfEdgeIsSeam(int32_t edgeIndex, FaceInfo face, int32_t loop,
                          Mesh *pMesh, EdgeVerts *pEdgeVerts);

uint32_t ruvmFnvHash(uint8_t *value, int32_t valueSize, uint32_t size);

int32_t checkIfEdgeIsPreserve(Mesh* pMesh, int32_t edge);
FaceTriangulated triangulateFace(RuvmAlloc alloc, FaceInfo baseFace, Mesh *pMesh);
V3_F32 getBarycentricInFace(V2_F32 *pTriUvs, int8_t *pTriLoops,
                            int32_t loopCount, V2_F32 vert);
void waitForJobs(RuvmContext pContext, int32_t *pJobsCompleted, void *pMutex);
