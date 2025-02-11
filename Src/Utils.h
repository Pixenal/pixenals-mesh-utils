#pragma once

#include <UvStuccoIntern.h>
#include <MathUtils.h>
#include <Mesh.h>

typedef struct {
	int32_t* pCorners;
	int32_t triCount;
	int32_t cornerCount;
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

typedef struct {
	UBitField8 *pBaseCorner;
	UBitField8 *pStucCorner;
	UBitField8 *pSegment;
	UBitField8 *pIsStuc;
	UBitField8 *pOnLine;
	UBitField8 *pOnInVert;
} BorderFaceBitArrs;

int32_t checkFaceIsInBounds(V2_F32 min, V2_F32 max, FaceRange face, Mesh *pMesh);
void getFaceBounds(FaceBounds *pBounds, V2_F32 *pStuc, FaceRange face);
int32_t checkIfEdgeIsSeam(int32_t edgeIdx, FaceRange face, int32_t corner,
                          Mesh *pMesh, EdgeVerts *pEdgeVerts);

uint32_t stucFnvHash(uint8_t *value, int32_t valueSize, uint32_t size);

bool checkIfEdgeIsPreserve(Mesh* pMesh, int32_t edge);
bool checkIfVertIsPreserve(Mesh* pMesh, int32_t vert);
int32_t checkIfEdgeIsReceive(Mesh* pMesh, int32_t edge);
FaceTriangulated triangulateFace(StucAlloc alloc, FaceRange baseFace, void *pVerts,
                                 int32_t *pCorners, int32_t useStuc);
V3_F32 getBarycentricInFace(V2_F32 *pTriStuc, int8_t *pTriCorners,
                            int32_t cornerCount, V2_F32 vert);
void buildEdgeList(StucContext pContext, Mesh* pMesh);
bool isMeshInvalid(Mesh* pMesh);
void progressBarClear();
void progressBarPrint(StucContext pContext, int32_t progress);
void stageBegin(void *pContext, StucStageReport* pReport, const char * pName);
void stageProgress(void *pContext, StucStageReport* pReport, int32_t progress);
void stageEnd(void *pContext, StucStageReport* pReport);
void stageBeginWrap(StucContext pContext, const char* pName, int32_t max);
void stageProgressWrap(StucContext pContext, int32_t progress);
void stageEndWrap(StucContext pContext);
void setStageName(StucContext pContext, const char* pName);
Mat3x3 buildFaceTbn(FaceRange face, Mesh *pMesh, int32_t *pCornerOveride);
void getTriScale(int32_t size, BaseTriVerts *pTri);
bool calcIntersection(V3_F32 a, V3_F32 b, V2_F32 c, V2_F32 cd,
                      V3_F32 *pPoint, float *pt, float *pt2);
int32_t idxBitArray(UBitField8 *pArr, int32_t idx, int32_t len);
void setBitArr(UBitField8 *pArr, int32_t idx, int32_t value, int32_t len);
void insertionSort(int32_t *pIdxTable, int32_t count, int32_t *pSort);
void fInsertionSort(int32_t *pIdxTable, int32_t count, float *pSort);
Mat3x3 getInterpolatedTbn(Mesh *pMesh, FaceRange *pFace,
                          int8_t *pTriCorners, V3_F32 bc);
int32_t calcFaceOrientation(Mesh *pMesh, FaceRange *pFace, bool useStuc);
int32_t getBorderFaceMemType(int32_t mapFaceSize, int32_t bufFaceSize);
int32_t getBorderFaceSize(int32_t memType);
void getBorderFaceBitArrs(BorderFace *pEntry, BorderFaceBitArrs *pArrs);
Result stucValidateAndDestroyJobs(StucContext pContext, int32_t jobCount,
                                  void ***pppJobHandles);