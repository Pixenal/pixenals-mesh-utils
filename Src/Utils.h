#pragma once

#include <UvStuccoIntern.h>
#include <MathUtils.h>
#include <Mesh.h>
#include <Types.h>

typedef struct {
	I32* pCorners;
	I32 triCount;
	I32 cornerCount;
} FaceTriangulated;

typedef struct {
	V2_I32 min, max;
	V2_F32 fMin, fMax;
	V2_F32 fMinSmall, fMaxSmall;
} FaceBounds;

typedef struct {
	V2_F32 uv[4];
	V3_F32 xyz[4];
	F32 scale[4];
} BaseTriVerts;

typedef struct {
	UBitField8 *pBaseCorner;
	UBitField8 *pStucCorner;
	UBitField8 *pSegment;
	UBitField8 *pIsStuc;
	UBitField8 *pOnLine;
	UBitField8 *pOnInVert;
} BorderFaceBitArrs;

I32 stucCheckFaceIsInBounds(V2_F32 min, V2_F32 max, FaceRange face, Mesh *pMesh);
void stucGetFaceBounds(FaceBounds *pBounds, V2_F32 *pStuc, FaceRange face);
I32 stucCheckIfEdgeIsSeam(I32 edgeIdx, FaceRange face, I32 corner,
                          Mesh *pMesh, EdgeVerts *pEdgeVerts);

U32 stucFnvHash(U8 *value, I32 valueSize, U32 size);

bool stucCheckIfEdgeIsPreserve(Mesh* pMesh, I32 edge);
bool stucCheckIfVertIsPreserve(Mesh* pMesh, I32 vert);
I32 stucCheckIfEdgeIsReceive(Mesh* pMesh, I32 edge);
FaceTriangulated stucTriangulateFace(StucAlloc alloc, FaceRange *pInFace, void *pVerts,
                                 I32 *pCorners, I32 useStuc);
V3_F32 stucGetBarycentricInFace(V2_F32 *pTriStuc, I8 *pTriCorners,
                            I32 cornerCount, V2_F32 vert);
void stucBuildEdgeList(StucContext pContext, Mesh* pMesh);
bool stucIsMeshInvalid(Mesh* pMesh);
void stucProgressBarClear();
void stucProgressBarPrint(StucContext pContext, I32 progress);
void stucStageBegin(void *pContext, StucStageReport* pReport, const char * pName);
void stucStageProgress(void *pContext, StucStageReport* pReport, I32 progress);
void stucStageEnd(void *pContext, StucStageReport* pReport);
void stucStageBeginWrap(StucContext pContext, const char* pName, I32 max);
void stucStageProgressWrap(StucContext pContext, I32 progress);
void stucStageEndWrap(StucContext pContext);
void stucSetStageName(StucContext pContext, const char* pName);
Mat3x3 stucBuildFaceTbn(FaceRange face, Mesh *pMesh, I32 *pCornerOveride);
void stucGetTriScale(I32 size, BaseTriVerts *pTri);
bool stucCalcIntersection(V3_F32 a, V3_F32 b, V2_F32 c, V2_F32 cd,
                      V3_F32 *pPoint, F32 *pt, F32 *pt2);
I32 stucIdxBitArray(UBitField8 *pArr, I32 idx, I32 len);
void stucSetBitArr(UBitField8 *pArr, I32 idx, I32 value, I32 len);
void stucInsertionSort(I32 *pIdxTable, I32 count, I32 *pSort);
void stucFInsertionSort(I32 *pIdxTable, I32 count, F32 *pSort);
Mat3x3 stucGetInterpolatedTbn(Mesh *pMesh, FaceRange *pFace,
                          I8 *pTriCorners, V3_F32 bc);
I32 stucCalcFaceOrientation(Mesh *pMesh, FaceRange *pFace, bool useStuc);
I32 stucGetBorderFaceMemType(I32 mapFaceSize, I32 bufFaceSize);
I32 stucGetBorderFaceSize(I32 memType);
void stucGetBorderFaceBitArrs(BorderFace *pEntry, BorderFaceBitArrs *pArrs);