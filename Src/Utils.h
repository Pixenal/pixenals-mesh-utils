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

int32_t stucCheckFaceIsInBounds(V2_F32 min, V2_F32 max, FaceRange face, Mesh *pMesh);
void stucGetFaceBounds(FaceBounds *pBounds, V2_F32 *pStuc, FaceRange face);
int32_t stucCheckIfEdgeIsSeam(int32_t edgeIdx, FaceRange face, int32_t corner,
                          Mesh *pMesh, EdgeVerts *pEdgeVerts);

uint32_t stucFnvHash(uint8_t *value, int32_t valueSize, uint32_t size);

bool stucCheckIfEdgeIsPreserve(Mesh* pMesh, int32_t edge);
bool stucCheckIfVertIsPreserve(Mesh* pMesh, int32_t vert);
int32_t stucCheckIfEdgeIsReceive(Mesh* pMesh, int32_t edge);
FaceTriangulated stucTriangulateFace(StucAlloc alloc, FaceRange *pInFace, void *pVerts,
                                 int32_t *pCorners, int32_t useStuc);
V3_F32 stucGetBarycentricInFace(V2_F32 *pTriStuc, int8_t *pTriCorners,
                            int32_t cornerCount, V2_F32 vert);
void stucBuildEdgeList(StucContext pContext, Mesh* pMesh);
bool stucIsMeshInvalid(Mesh* pMesh);
void stucProgressBarClear();
void stucProgressBarPrint(StucContext pContext, int32_t progress);
void stucStageBegin(void *pContext, StucStageReport* pReport, const char * pName);
void stucStageProgress(void *pContext, StucStageReport* pReport, int32_t progress);
void stucStageEnd(void *pContext, StucStageReport* pReport);
void stucStageBeginWrap(StucContext pContext, const char* pName, int32_t max);
void stucStageProgressWrap(StucContext pContext, int32_t progress);
void stucStageEndWrap(StucContext pContext);
void stucSetStageName(StucContext pContext, const char* pName);
Mat3x3 stucBuildFaceTbn(FaceRange face, Mesh *pMesh, int32_t *pCornerOveride);
void stucGetTriScale(int32_t size, BaseTriVerts *pTri);
bool stucCalcIntersection(V3_F32 a, V3_F32 b, V2_F32 c, V2_F32 cd,
                      V3_F32 *pPoint, float *pt, float *pt2);
int32_t stucIdxBitArray(UBitField8 *pArr, int32_t idx, int32_t len);
void stucSetBitArr(UBitField8 *pArr, int32_t idx, int32_t value, int32_t len);
void stucInsertionSort(int32_t *pIdxTable, int32_t count, int32_t *pSort);
void stucFInsertionSort(int32_t *pIdxTable, int32_t count, float *pSort);
Mat3x3 stucGetInterpolatedTbn(Mesh *pMesh, FaceRange *pFace,
                          int8_t *pTriCorners, V3_F32 bc);
int32_t stucCalcFaceOrientation(Mesh *pMesh, FaceRange *pFace, bool useStuc);
int32_t stucGetBorderFaceMemType(int32_t mapFaceSize, int32_t bufFaceSize);
int32_t stucGetBorderFaceSize(int32_t memType);
void stucGetBorderFaceBitArrs(BorderFace *pEntry, BorderFaceBitArrs *pArrs);