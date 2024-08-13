#pragma once

#include <RuvmInternal.h>
#include <MathUtils.h>
#include <Mesh.h>

typedef struct {
	int32_t* pLoops;
	int32_t triCount;
	int32_t loopCount;
} FaceTriangulated;

typedef struct {
	V2_I32 min, max;
	V2_F32 fMin, fMax;
	V2_F32 fMinSmall, fMaxSmall;
} FaceBounds;

typedef struct {
	V2_F32 uv[4];
	V3_F32 xyz[4];
	float scale[4];
} BaseTriVerts;

int32_t checkFaceIsInBounds(V2_F32 min, V2_F32 max, FaceRange face, Mesh *pMesh);
void getFaceBounds(FaceBounds *pBounds, V2_F32 *pUvs, FaceRange face);
int32_t checkIfEdgeIsSeam(int32_t edgeIndex, FaceRange face, int32_t loop,
                          Mesh *pMesh, EdgeVerts *pEdgeVerts);

uint32_t ruvmFnvHash(uint8_t *value, int32_t valueSize, uint32_t size);

bool checkIfEdgeIsPreserve(Mesh* pMesh, int32_t edge);
bool checkIfVertIsPreserve(Mesh* pMesh, int32_t vert);
int32_t checkIfEdgeIsReceive(Mesh* pMesh, int32_t edge);
FaceTriangulated triangulateFace(RuvmAlloc alloc, FaceRange baseFace, void *pVerts,
                                 int32_t *pLoops, int32_t useUvs);
V3_F32 getBarycentricInFace(V2_F32 *pTriUvs, int8_t *pTriLoops,
                            int32_t loopCount, V2_F32 vert);
void waitForJobs(RuvmContext pContext, int32_t *pJobsCompleted, void *pMutex);
void buildEdgeList(RuvmContext pContext, Mesh* pMesh);
bool isMeshInvalid(Mesh* pMesh);
void progressBarClear();
void progressBarPrint(RuvmContext pContext, int32_t progress);
void stageBegin(void *pContext, RuvmStageReport* pReport, const char * pName);
void stageProgress(void *pContext, RuvmStageReport* pReport, int32_t progress);
void stageEnd(void *pContext, RuvmStageReport* pReport);
void stageBeginWrap(RuvmContext pContext, const char* pName, int32_t max);
void stageProgressWrap(RuvmContext pContext, int32_t progress);
void stageEndWrap(RuvmContext pContext);
void setStageName(RuvmContext pContext, const char* pName);
Mat3x3 buildFaceTbn(FaceRange face, Mesh *pMesh, int32_t *pLoopOveride);
void getTriScale(int32_t size, BaseTriVerts *pTri);
bool calcIntersection(V3_F32 a, V3_F32 b, V2_F32 c, V2_F32 cd,
                      V3_F32 *pPoint, float *pt, float *pt2);
bool indexBitArray(UBitField8 *pArr, int32_t index);
void setBitArr(UBitField8 *pArr, int32_t index, bool value);
void fInsertionSort(int32_t *pIndexTable, int32_t count, float *pSort);
Mat3x3 getInterpolatedTbn(Mesh *pMesh, FaceRange *pFace,
                          int8_t *pTriLoops, V3_F32 bc);
int32_t calcFaceOrientation(Mesh *pMesh, FaceRange *pFace, bool useUvs);