#pragma once

#include <uv_stucco_intern.h>
#include <math_utils.h>
#include <mesh.h>
#include <types.h>

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
	V3_F32 xyz[4];
	V2_F32 uv[4];
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

I32 stucCheckFaceIsInBounds(V2_F32 min, V2_F32 max, FaceRange face, const Mesh *pMesh);
void stucGetFaceBounds(FaceBounds *pBounds, const V2_F32 *pUvs, FaceRange face);
I32 stucCheckIfEdgeIsSeam(
	I32 edgeIdx,
	FaceRange face,
	I32 corner,
	const Mesh *pMesh
);
U32 stucFnvHash(const U8 *value, I32 valueSize, U32 size);
bool stucCheckIfEdgeIsPreserve(const Mesh *pMesh, I32 edge);
bool stucCheckIfVertIsPreserve(const Mesh *pMesh, I32 vert);
bool stucCheckIfEdgeIsReceive(const Mesh *pMesh, I32 edge, F32 receiveLen);
FaceTriangulated stucTriangulateFace(
	const StucAlloc alloc,
	const FaceRange *pInFace,
	const void *pVerts,
	const I32 *pCorners,
	I32 useStuc
);
V3_F32 stucGetBarycentricInFace(
	const V2_F32 *pTriStuc,
	I8 *pTriCorners,
	I32 cornerCount,
	V2_F32 vert
);
StucResult stucBuildEdgeList(StucContext pCtx, Mesh *pMesh);
void stucProgressBarClear();
void stucProgressBarPrint(StucContext pCtx, I32 progress);
void stucStageBegin(void *pCtx, StucStageReport *pReport, const char *pName);
void stucStageProgress(void *pCtx, StucStageReport *pReport, I32 progress);
void stucStageEnd(void *pCtx, StucStageReport *pReport);
void stucStageBeginWrap(StucContext pCtx, const char *pName, I32 max);
void stucStageProgressWrap(StucContext pCtx, I32 progress);
void stucStageEndWrap(StucContext pCtx);
void stucSetStageName(StucContext pCtx, const char *pName);
Mat3x3 stucBuildFaceTbn(FaceRange face, const Mesh *pMesh, const I32 *pCornerOveride);
void stucGetTriScale(I32 size, BaseTriVerts *pTri);
bool stucCalcIntersection(
	V3_F32 a,
	V3_F32 b,
	V2_F32 c,
	V2_F32 cd,
	V3_F32 *pPoint,
	F32 *pt,
	F32 *pt2
);
I32 stucIdxBitArray(UBitField8 *pArr, I32 idx, I32 len);
void stucSetBitArr(UBitField8 *pArr, I32 idx, I32 value, I32 len);
void stucInsertionSort(I32 *pIdxTable, I32 count, I32 *pSort);
void stucFInsertionSort(I32 *pIdxTable, I32 count, F32 *pSort);
Mat3x3 stucGetInterpolatedTbn(
	const Mesh *pMesh,
	const FaceRange *pFace,
	const I8 *pTriCorners,
	V3_F32 bc
);
I32 stucCalcFaceOrientation(const Mesh *pMesh, const FaceRange *pFace, bool useStuc);
I32 stucGetBorderFaceMemType(I32 mapFaceSize, I32 bufFaceSize);
I32 stucGetBorderFaceSize(I32 memType);
Result stucAllocBorderFace(I32 memType, BorderTableAlloc *pHandles, void **ppOut);
void stucGetBorderFaceBitArrs(BorderFace *pEntry, BorderFaceBitArrs *pArrs);
void stucBorderTableDestroyAlloc(BorderTableAlloc *pTableAlloc);