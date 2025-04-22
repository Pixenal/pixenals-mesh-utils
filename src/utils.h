/* 
SPDX-FileCopyrightText: 2025 Caleb Dawson
SPDX-License-Identifier: Apache-2.0
*/

#pragma once

#include <float.h>

#include <uv_stucco_intern.h>
#include <math_utils.h>
#include <mesh.h>
#include <types.h>
#include <thread_pool.h>

typedef Stuc_String String;

typedef struct BBox {
	V2_F32 min;
	V2_F32 max;
} BBox;

typedef struct FaceBounds {
	V2_I32 min, max;
	BBox fBBox;
	BBox fBBoxSmall;
} FaceBounds;

typedef struct BaseTriVerts {
	V3_F32 xyz[4];
	V2_F32 uv[4];
	F32 scale[4];
} BaseTriVerts;

typedef enum InsideStatus {
	STUC_INSIDE_STATUS_NONE,
	STUC_INSIDE_STATUS_OUTSIDE,
	STUC_INSIDE_STATUS_INSIDE,
	STUC_INSIDE_STATUS_ON_LINE
} InsideStatus;

typedef enum StucCompare {
	STUC_COMPARE_LESS,
	STUC_COMPARE_EQUAL,
	STUC_COMPARE_GREAT
} StucCompare;

//TODO remove unused fields
typedef struct HalfPlane {
	V2_F32 uv;
	V2_F32 uvNext;
	V2_F32 dir;
	V2_F32 dirUnit;
	V2_F32 halfPlane;
	F32 len;
	I32 edge;
	I8 idx;
	I8 idxPrev;
	I8 idxNext;
	bool flipEdgeDir;
} HalfPlane;

static inline
void initHalfPlaneLookup(
	const Mesh *pMesh,
	const FaceRange *pInFace,
	HalfPlane *pCache
) {
	for (I32 i = 0; i < pInFace->size; ++i) {
		pCache[i].idx = (I8)i;
		pCache[i].idxNext = (I8)stucGetCornerNext(i, pInFace);
		pCache[i].idxPrev = (I8)stucGetCornerPrev(i, pInFace);

		pCache[i].edge = stucGetMeshEdge(
			&pMesh->core,
			(FaceCorner) {.face = pInFace->idx, .corner = i}
		);

		pCache[i].uv = pMesh->pUvs[pInFace->start + i];
		pCache[i].uvNext = pMesh->pUvs[pInFace->start + pCache[i].idxNext];
		
		pCache[i].dir = _(pCache[i].uvNext V2SUB pCache[i].uv);
		pCache[i].len = pixmV2F32Len(pCache[i].dir);
		pCache[i].dirUnit = _(pCache[i].dir V2DIVS pCache[i].len);
		pCache[i].halfPlane = pixmV2F32LineNormal(pCache[i].dir);
	}
}

BBox stucBBoxGet(const Mesh *pMesh, FaceRange *pFace);
static inline
bool stucIsBBoxInBBox(BBox bboxA, BBox bboxB) {
	V2_I32 inside = {0};
	inside.d[0] =
		(bboxB.min.d[0] >= bboxA.min.d[0] && bboxB.min.d[0] < bboxA.max.d[0]) ||
		(bboxB.max.d[0] >= bboxA.min.d[0] && bboxB.max.d[0] < bboxA.max.d[0]) ||
		(bboxB.min.d[0] < bboxA.min.d[0] && bboxB.max.d[0] >= bboxA.max.d[0]);
	inside.d[1] =
		(bboxB.min.d[1] >= bboxA.min.d[1] && bboxB.min.d[1] < bboxA.max.d[1]) ||
		(bboxB.max.d[1] >= bboxA.min.d[1] && bboxB.max.d[1] < bboxA.max.d[1]) ||
		(bboxB.min.d[1] < bboxA.min.d[1] && bboxB.max.d[1] >= bboxA.max.d[1]);
	return inside.d[0] && inside.d[1];
}
void stucGetInFaceBounds(FaceBounds *pBounds, const V2_F32 *pUvs, FaceRange face);
I32 stucIsEdgeSeam(const Mesh *pMesh, I32 edge);
U32 stucFnvHash(const U8 *value, I32 valueSize, U32 size);
bool stucGetIfPreserveEdge(const Mesh *pMesh, I32 edge);
bool stucCheckIfVertIsPreserve(const Mesh *pMesh, I32 vert);
bool stucCheckIfEdgeIsReceive(const Mesh *pMesh, I32 edge, F32 receiveLen);

typedef struct Ear {
	struct Ear *pNext;
	struct Ear *pPrev;
	I32 corner;
	F32 len;
} Ear;

typedef struct TriangulateState {
	PixalcLinAlloc earAlloc;
	Ear *pEarList;
	const Mesh *pMesh;
	const FaceRange *pFace;
	V2_F32 (* fpGetPoint)(const Mesh *, const FaceRange *, I32);
	FaceTriangulated *pTris;
	bool *pRemoved;
	I32 triCount;
	bool wind;
} TriangulateState;

InsideStatus stucIsPointInHalfPlane(
	V2_F32 point,
	V2_F32 lineA,
	V2_F32 halfPlane,
	bool wind
);

STUC_FORCE_INLINE
bool doesEarIntersectFace(
	const TriangulateState *pState,
	I32 a, I32 b, I32 c,
	V2_F32 aPos, V2_F32 bPos, V2_F32 cPos
) {
	bool wind = pState->wind;
	V2_F32 abHalfPlane = pixmV2F32LineNormal(_(bPos V2SUB aPos));
	V2_F32 bcHalfPlane = pixmV2F32LineNormal(_(cPos V2SUB bPos));
	V2_F32 caHalfPlane = pixmV2F32LineNormal(_(aPos V2SUB cPos));
	for (I32 i = 0; i < pState->pFace->size; ++i) {
		if (i == a || i == b || i == c) {
			continue;
		}
		V2_F32 point = pState->fpGetPoint(pState->pMesh, pState->pFace, i);
		if (stucIsPointInHalfPlane(point, aPos, abHalfPlane, wind) == STUC_INSIDE_STATUS_OUTSIDE
		) {
			continue;
		}
		if (stucIsPointInHalfPlane(point, bPos, bcHalfPlane, wind) == STUC_INSIDE_STATUS_OUTSIDE
		) {
			continue;
		}
		if (stucIsPointInHalfPlane(point, cPos, caHalfPlane, wind) != STUC_INSIDE_STATUS_OUTSIDE
		) {
			return true;
		}
	}
	return false;
}

static inline
I32 stucGetNextRemaining(
	const TriangulateState *pState,
	I32 corner,
	const FaceRange *pFace
) {
	PIX_ERR_ASSERT("", corner >= 0 && corner < pFace->size);
	I32 start = corner;
	while (corner = stucGetCornerNext(corner, pFace), corner != start) {
		if (!pState->pRemoved[corner]) {
			return corner;
		}
	}
	PIX_ERR_ASSERT("no corners remain", false);
	return -1;
}

static inline
I32 stucGetPrevRemaining(
	const TriangulateState *pState,
	I32 corner,
	const FaceRange *pFace
) {
	PIX_ERR_ASSERT("", corner >= 0 && corner < pFace->size);
	I32 start = corner;
	while (corner = stucGetCornerPrev(corner, pFace), corner != start) {
		if (!pState->pRemoved[corner]) {
			return corner;
		}
	}
	PIX_ERR_ASSERT("no corners remain", false);
	return -1;
}

STUC_FORCE_INLINE
Ear *addEarCandidate(TriangulateState *pState, I32 corner) {
	I32 cornerPrev = stucGetPrevRemaining(pState, corner, pState->pFace);
	I32 cornerNext = stucGetNextRemaining(pState, corner, pState->pFace);
	V2_F32 a = pState->fpGetPoint(pState->pMesh, pState->pFace, cornerPrev);
	V2_F32 b = pState->fpGetPoint(pState->pMesh, pState->pFace, corner);
	V2_F32 c = pState->fpGetPoint(pState->pMesh, pState->pFace, cornerNext);
	V2_F32 ac = _(c V2SUB a);
	F32 dot = _(ac V2DOT pixmV2F32LineNormal(_(b V2SUB a)));
	if (dot >= 0 || // ear is concave or degenerate
		doesEarIntersectFace(pState, cornerPrev, corner, cornerNext, a, b, c)
	) {
		return NULL;
	}
	F32 len = pixmV2F32Len(ac);
	Ear *pNewEar = NULL;
	if (!pState->pEarList) {
		pixalcLinAlloc(&pState->earAlloc, (void **)&pState->pEarList, 1);
		pNewEar = pState->pEarList;
	}
	else {
		Ear *pEar = pState->pEarList;
		while(pEar->pNext && len > pEar->pNext->len) {
			pEar = pEar->pNext;
		}
		pixalcLinAlloc(&pState->earAlloc, (void **)&pNewEar, 1);
		if (len < pEar->len) {
			pEar->pPrev = pNewEar;
			pNewEar->pNext = pEar;
			pState->pEarList = pNewEar;
		}
		else {
			if (pEar->pNext) {
				pEar->pNext->pPrev = pNewEar;
				pNewEar->pNext = pEar->pNext;
			}
			pNewEar->pPrev = pEar;
			pEar->pNext = pNewEar;
		}
	}
	pNewEar->corner = corner;
	pNewEar->len = len;
	return pNewEar;
}

STUC_FORCE_INLINE
void addAdjEarCandidates(TriangulateState *pState, Ear *pEar) {
	I32 cornerNext = stucGetNextRemaining(pState, pEar->corner, pState->pFace);
	I32 cornerPrev = stucGetPrevRemaining(pState, pEar->corner, pState->pFace);
	addEarCandidate(pState, cornerNext);
	addEarCandidate(pState, cornerPrev);
}

static inline
Ear *addEar(TriangulateState *pState) {
	Ear *pEar = pState->pEarList;
	U8 *pTri = pState->pTris->pTris + pState->triCount * 3;
	I32 cornerPrev = stucGetPrevRemaining(pState, pEar->corner, pState->pFace);
	I32 cornerNext = stucGetNextRemaining(pState, pEar->corner, pState->pFace);
	pTri[0] = cornerPrev;
	pTri[1] = pEar->corner;
	pTri[2] = cornerNext;
	pState->triCount++;
	
	pState->pRemoved[pEar->corner] = true;
	return pEar;
}

/*
STUC_FORCE_INLINE
void cacheHalfPlanes(TriangulateState *pState) {
	PIX_ERR_ASSERT("", pState->pHalfPlanes);
	for (I32 i = 0; i < pState->pFace->size; ++i) {
		I32 iNext = stucGetCornerNext(i, pState->pFace);
		V2_F32 a = pState->fpGetPoint(pState->pMesh, pState->pFace, i);
		V2_F32 b = pState->fpGetPoint(pState->pMesh, pState->pFace, iNext);
		V2_F32 ab = _(b V2SUB a);
		pState->pHalfPlanes[i] = v2F32LineNormal(ab);
	}
}
*/

static inline
void stucRemoveEar(TriangulateState *pState) {
	Ear *pEar = pState->pEarList;
	if (pEar->pNext) {
		pEar->pNext->pPrev = NULL;
	}
	pState->pEarList = pEar->pNext;
}

static inline
bool isMarkedSkip(I32 *pSkip, I32 skipCount, I32 idx) {
	for (I32 i = 0; i < skipCount; ++i) {
		if (idx == pSkip[i]) {
			return true;
		}
	}
	return false;
}
//finds corner on convex hull of face, & determines wind direction from that
//returns 0 for clockwise, 1 for counterclockwise, & 2 if degenerate
STUC_FORCE_INLINE
I32 stucCalcFaceWind(
	const FaceRange *pFace,
	const Mesh *pMesh,
	V2_F32 (* fpGetPoint) (const Mesh *, const FaceRange *, I32)
) {
	PIX_ERR_ASSERT("", pFace->start >= 0 && pFace->size >= 3);
	I32 skip[32] = {0};
	I32 skipCount = 0;
	do {
		I32 lowestCorner = 0;
		V2_F32 lowestCoord = { FLT_MAX, FLT_MAX };
		for (I32 i = 0; i < pFace->size; ++i) {
			if (isMarkedSkip(skip, skipCount, i)) {
				continue;
			}
			V2_F32 pos = fpGetPoint(pMesh, pFace, i);
			if (pos.d[0] > lowestCoord.d[0]) {
				continue;
			}
			else if (pos.d[0] == lowestCoord.d[0] &&
			         pos.d[1] >= lowestCoord.d[1]) {
				continue;
			}
			lowestCorner = i;
			lowestCoord = pos;
		}
		I32 prev = lowestCorner == 0 ? pFace->size - 1 : lowestCorner - 1;
		I32 next = (lowestCorner + 1) % pFace->size;
		V2_F32 a = fpGetPoint(pMesh, pFace, prev);
		V2_F32 b = fpGetPoint(pMesh, pFace, lowestCorner);
		V2_F32 c = fpGetPoint(pMesh, pFace, next);
		//alt formula for determinate,
		//shorter and less likely to cause numerical error
		F32 det =
			(b.d[0] - a.d[0]) * (c.d[1] - a.d[1]) -
			(c.d[0] - a.d[0]) * (b.d[1] - a.d[1]);
		if (det) {
			return det > .0f;
		}
		//abc is degenerate, find another corner
		skip[skipCount] = lowestCorner;
		skipCount++;
	} while(skipCount < pFace->size);
	return 2;
}
static inline
I32 stucCalcFaceWindFromVerts(const FaceRange *pFace, const Mesh *pMesh) {
	return stucCalcFaceWind(pFace, pMesh, stucGetVertPosAsV2);
}
static inline
I32 stucCalcFaceWindFromUvs(const FaceRange *pFace, const Mesh *pMesh) {
	return stucCalcFaceWind(pFace, pMesh, stucGetUvPos);
}

STUC_FORCE_INLINE
void stucTriangulateFace(
	const StucAlloc *pAlloc,
	const FaceRange *pFace,
	const Mesh *pMesh,
	V2_F32 (* fpGetPoint)(const Mesh *, const FaceRange *, I32),
	FaceTriangulated *pTris
) {
	PIX_ERR_ASSERT("", pTris->pTris);
	TriangulateState state = {
		.pMesh = pMesh,
		.fpGetPoint = fpGetPoint,
		.pFace = pFace,
		.wind = stucCalcFaceWind(pFace, pMesh, fpGetPoint),
		.pTris = pTris,
		.pRemoved = pAlloc->fpCalloc(pFace->size, sizeof(bool))
	};
	PIX_ERR_ASSERT(
		"degenerate faces shouldn't be passed to this func",
		state.wind % 2 == state.wind
	);

	pixalcLinAllocInit(pAlloc, &state.earAlloc, sizeof(Ear), pFace->size, true);

	//add initial ears
	for (I32 i = 0; i < pFace->size; ++i) {
		addEarCandidate(&state, i);
	}
	PIX_ERR_ASSERT("", state.pEarList);
	do {
		Ear *pAddedEar = NULL;
		if (!state.pRemoved[state.pEarList->corner]) {
			pAddedEar = addEar(&state);
		}
		stucRemoveEar(&state);
		if (pAddedEar) {
			addAdjEarCandidates(&state, pAddedEar);
		}
	} while(state.triCount < pFace->size - 2);
	pixalcLinAllocDestroy(&state.earAlloc);
	pAlloc->fpFree(state.pRemoved);
	PIX_ERR_ASSERT("", state.triCount == pFace->size - 2);
	return;
}

static inline
void stucTriangulateFaceFromVerts(
	const StucAlloc *pAlloc,
	const FaceRange *pFace,
	const Mesh *pMesh,
	FaceTriangulated *pTris
) {
	stucTriangulateFace(pAlloc, pFace, pMesh, stucGetVertPosAsV2, pTris);
}

static inline
void stucTriangulateFaceFromUvs(
	const StucAlloc *pAlloc,
	const FaceRange *pFace,
	const Mesh *pMesh,
	FaceTriangulated *pTris
) {
	stucTriangulateFace(pAlloc, pFace, pMesh, stucGetUvPos, pTris);
}

STUC_FORCE_INLINE
V3_F32 stucGetBarycentricInTri(
	const Mesh *pMesh,
	const FaceRange *pFace,
	V2_F32 (* fpGetPoint)(const Mesh *, const FaceRange *, I32),
	const U8 *pTriCorners,
	V2_F32 vert
) {
	V2_F32 tri[3] = {0};
	for (I32 i = 0; i < 3; ++i) {
		tri[i] = fpGetPoint(pMesh, pFace, pTriCorners[i]);
	}
	return pixmCartesianToBarycentric(tri, &vert);

}

static inline
V3_F32 stucGetBarycentricInTriFromVerts(
	const Mesh *pMesh,
	const FaceRange *pFace,
	const U8 *pTriCorners,
	V2_F32 vert
) {
	return stucGetBarycentricInTri(pMesh, pFace, stucGetVertPosAsV2, pTriCorners, vert);
}

static inline
V3_F32 stucGetBarycentricInTriFromUvs(
	const Mesh *pMesh,
	const FaceRange *pFace,
	const U8 *pTriCorners,
	V2_F32 vert
) {
	return stucGetBarycentricInTri(pMesh, pFace, stucGetUvPos, pTriCorners, vert);
}

//Caller must check for nan in return value
STUC_FORCE_INLINE
V3_F32 stucGetBarycentricInFace(
	const Mesh *pMesh,
	const FaceRange *pFace,
	V2_F32 (* fpGetPoint)(const Mesh *, const FaceRange *, I32),
	I8 *pTriCorners,
	V2_F32 vert
) {
	PIX_ERR_ASSERT("", pixmV2F32IsFinite(vert));
	PIX_ERR_ASSERT("", (pFace->size == 3 || pFace->size == 4) && pTriCorners);
	V2_F32 triA[3] = {0};
	for (I32 i = 0; i < 3; ++i) {
		triA[i] = fpGetPoint(pMesh, pFace, i);
	}
	V3_F32 vertBc = pixmCartesianToBarycentric(triA, &vert);
	if (pFace->size == 4 && pixmV3F32IsFinite(vertBc) && vertBc.d[1] < 0) {
		//base face is a quad, and vert is outside first tri,
		//so use the second tri
		
		V2_F32 triB[3] = {triA[2], fpGetPoint(pMesh, pFace, 3), triA[0]};
		vertBc = pixmCartesianToBarycentric(triB, &vert);
		pTriCorners[0] = 2;
		pTriCorners[1] = 3;
	}
	else {
		for (I32 k = 0; k < 3; ++k) {
			pTriCorners[k] = (I8)k;
		}
	}
	return vertBc;
}

static inline
V3_F32 stucGetBarycentricInFaceFromVerts(
	const Mesh *pMesh,
	const FaceRange *pFace,
	I8 *pTriCorners,
	V2_F32 vert
) {
	return stucGetBarycentricInFace(
		pMesh,
		pFace,
		stucGetVertPosAsV2,
		pTriCorners,
		vert
	);
}

static inline
V3_F32 stucGetBarycentricInFaceFromUvs(
	const Mesh *pMesh,
	const FaceRange *pFace,
	I8 *pTriCorners,
	V2_F32 vert
) {
	return stucGetBarycentricInFace(
		pMesh,
		pFace,
		stucGetUvPos,
		pTriCorners,
		vert
	);
}


StucErr stucBuildEdgeList(StucContext pCtx, Mesh *pMesh);
void stucProgressBarClear();
void stucProgressBarPrint(StucContext pCtx, I32 progress);
void stucStageBegin(void *pCtx, StucStageReport *pReport, const char *pName);
void stucStageProgress(void *pCtx, StucStageReport *pReport, I32 progress);
void stucStageEnd(void *pCtx, StucStageReport *pReport);
void stucStageBeginWrap(StucContext pCtx, const char *pName, I32 max);
void stucStageProgressWrap(StucContext pCtx, I32 progress);
void stucStageEndWrap(StucContext pCtx);
void stucSetStageName(StucContext pCtx, const char *pName);
M3x3 stucBuildFaceTbn(FaceRange face, const Mesh *pMesh, const I32 *pCornerOveride);
void stucGetTriScale(I32 size, BaseTriVerts *pTri);
F32 stucGetT(V2_F32 point, V2_F32 lineA, V2_F32 lineUnit, F32 lineLen);
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
STUC_FORCE_INLINE
void stucInsertionSort(
	I32 *pIdxArr,
	I32 count,
	const void *pData,
	StucCompare (*fpCompare)(const void *, I32, I32)
) {
	bool order = fpCompare(pData, 0, 1) == STUC_COMPARE_LESS;
	pIdxArr[0] = !order;
	pIdxArr[1] = order;
	I32 bufSize = 2;
	for (I32 i = bufSize; i < count; ++i) {
		bool insert = false;
		I32 j;
		for (j = bufSize - 1; j >= 0; --j) {
			insert =
				fpCompare(pData, i, pIdxArr[j]) == STUC_COMPARE_LESS &&
				fpCompare(pData, i, pIdxArr[j - 1]) == STUC_COMPARE_GREAT;
			if (insert) {
				break;
			}
		}
		if (!insert) {
			pIdxArr[bufSize] = i;
		}
		else {
			for (I32 m = bufSize; m > j; --m) {
				pIdxArr[m] = pIdxArr[m - 1];
				PIX_ERR_ASSERT("", m <= bufSize && m > j);
			}
			pIdxArr[j] = i;
		}
		bufSize++;
	}
}
//void stucFInsertionSort(I32 *pIdxTable, I32 count, F32 *pSort);
M3x3 stucGetInterpolatedTbn(
	const Mesh *pMesh,
	const FaceRange *pFace,
	const I8 *pTriCorners,
	V3_F32 bc
);

StucErr stucDoJobInParallel(
	const MapToMeshBasic *pBasic,
	I32 jobCount, void *pJobArgs, I32 argStructSize,
	StucErr (* func)(void *)
);
//U32 stucGetEncasedFaceHash(I32 mapFace, V2_I16 tile, I32 tableSize);

typedef struct HTableEntryCore {
	struct HTableEntryCore *pNext;
} HTableEntryCore;

typedef struct HTableBucket {
	HTableEntryCore *pList;
} HTableBucket;

#define STUC_HTABLE_ALLOC_HANDLES_MAX 2

typedef struct HTable {
	const StucAlloc *pAlloc;
	PixalcLinAlloc allocHandles[STUC_HTABLE_ALLOC_HANDLES_MAX];
	HTableBucket *pTable;
	void *pUserData;
	I32 size;
} HTable;

typedef enum SearchResult {
	STUC_SEARCH_FOUND,
	STUC_SEARCH_NOT_FOUND,
	STUC_SEARCH_ADDED
} SearchResult;

void stucHTableInit(
	const StucAlloc *pAlloc,
	HTable *pHandle,
	I32 targetSize,
	I32Arr allocTypeSizes,
	void *pUserData
);
void stucHTableDestroy(HTable *pHandle);
PixalcLinAlloc *stucHTableAllocGet(HTable *pHandle, I32 idx);
const PixalcLinAlloc *stucHTableAllocGetConst(const HTable *pHandle, I32 idx);
static inline
HTableBucket *stucHTableBucketGet(HTable *pHandle, U64 key) {
	U64 hash = stucFnvHash((U8 *)&key, sizeof(key), pHandle->size);
	return pHandle->pTable + hash;
}
STUC_FORCE_INLINE
SearchResult stucHTableGet(
	HTable *pHandle,
	I32 alloc,
	const void *pKeyData,
	void **ppEntry,
	bool addEntry,
	void *pInitInfo,
	U64 (* fpMakeKey)(const void *),
	bool (* fpAddPredicate)(const void *, const void *, const void *),
	void (* fpInitEntry)(void *, HTableEntryCore *, const void *, void *, I32),
	bool (* fpCompareEntry)(const HTableEntryCore *, const void *, const void *)
) {
	PIX_ERR_ASSERT("", pHandle->pTable && pHandle->size);
	PIX_ERR_ASSERT(
		"",
		alloc < STUC_HTABLE_ALLOC_HANDLES_MAX && pHandle->allocHandles[alloc].valid
	);
	PIX_ERR_ASSERT("", (!addEntry || fpInitEntry) && fpCompareEntry);
	HTableBucket *pBucket = stucHTableBucketGet(pHandle, fpMakeKey(pKeyData));
	if (!pBucket->pList) {
		if (!addEntry ||
			fpAddPredicate && !fpAddPredicate(pHandle->pUserData, pKeyData, pInitInfo)
		) {
			return STUC_SEARCH_NOT_FOUND;
		}
		I32 linIdx =
			pixalcLinAlloc(pHandle->allocHandles + alloc, (void **)&pBucket->pList, 1);
		fpInitEntry(pHandle->pUserData, pBucket->pList, pKeyData, pInitInfo, linIdx);
		if (ppEntry) {
			*ppEntry = pBucket->pList;
		}
		return STUC_SEARCH_ADDED;
	}
	HTableEntryCore *pEntry = pBucket->pList;
	do {
		if (fpCompareEntry(pEntry, pKeyData, pInitInfo)) {
			if (ppEntry) {
				*ppEntry = pEntry;
			}
			return STUC_SEARCH_FOUND;
		}
		if (!pEntry->pNext) {
			if (!addEntry ||
				fpAddPredicate && !fpAddPredicate(pHandle->pUserData, pKeyData, pInitInfo)
			) {
				return STUC_SEARCH_NOT_FOUND;
			}
			I32 linIdx =
				pixalcLinAlloc(pHandle->allocHandles + alloc, (void **)&pEntry->pNext, 1);
			fpInitEntry(pHandle->pUserData, pEntry->pNext, pKeyData, pInitInfo, linIdx);
			if (ppEntry) {
				*ppEntry = pEntry->pNext;
			}
			return STUC_SEARCH_ADDED;
		}
		pEntry = pEntry->pNext;
	} while(true);
}

STUC_FORCE_INLINE
SearchResult stucHTableGetConst(
	HTable *pHandle,
	I32 alloc,
	const void *pKeyData,
	void **ppEntry,
	bool addEntry,
	const void *pInitInfo,
	U64 (* fpMakeKey)(const void *),
	bool (* fpAddPredicate)(const void *, const void *, const void *),
	void (* fpInitEntry)(void *, HTableEntryCore *, const void *, void *, I32),
	bool (* fpCompareEntry)(const HTableEntryCore *, const void *, const void *)
) {
	return stucHTableGet(
		pHandle,
		alloc,
		pKeyData,
		ppEntry,
		addEntry,
		(void *)pInitInfo,
		fpMakeKey,
		fpAddPredicate,
		fpInitEntry,
		fpCompareEntry
	);
}

static inline
bool stucHTableCmpFalse(
	const HTableEntryCore *pEntry,
	const void *pKeyData,
	const void *pInitInfo
) {
	return false;
}

U64 stucKeyFromI32(const void *pKeyData);

typedef struct EncasingInFace {
	U32 idx : 31;
	U32 wind : 1;
} EncasingInFace;

typedef struct EncasingInFaceArr {
	EncasingInFace *pArr;
	I32 size;
	I32 count;
} EncasingInFaceArr;

/*
typedef struct EncasedMapCorner {
	U32 inFace : 30; //index withing encasing in-face arr
	InsideStatus status : 2;
} EncasedMapCorner;

typedef struct EncasedMapCornerArr {
	EncasedMapCorner *pArr;
	I32 faceSize;
	I32 count;
} EncasedMapCornerArr;
*/

typedef struct EncasedMapFace {
	HTableEntryCore core;
	EncasingInFaceArr inFaces;
	I32 mapFace;
	V2_I16 tile;
} EncasedMapFace;

typedef struct EncasedMapFaceTableState {
	const MapToMeshBasic *pBasic;
} EncasedMapFaceTableState;

/*
typedef struct {
	EncasedMapFace *pList;
} EncasedMapFaceBucket;

typedef struct {
	EncasedMapFaceBucket *pTable;
	I32 size;
	I32 entryCount;
} EncasedFacesTable;
*/

typedef struct EncasedEntryIdx {
	struct EncasedEntryIdx *pNext;
	I32 mapFace;
	V2_I16 tile;
	I32 entryIdx;
} EncasedEntryIdx;

/*
typedef struct EncasedEntryIdxBucket {
	EncasedEntryIdx *pList;
} EncasedEntryIdxBucket;

typedef struct {
	EncasedEntryIdxBucket *pTable;
	I32 size;
} EncasedEntryIdxTable;
*/

typedef struct Border {
	FaceCorner start;
	I32 len;
} Border;

typedef struct BorderArr {
	Border *pArr;
	I32 size;
	I32 count;
} BorderArr;

typedef struct BufFace {
	I32 start;
	I32 size;
	I32 inPiece;
} BufFace;

typedef struct InPiece {
	EncasedMapFace *pList;
	BorderArr borderArr;
	I32 faceCount;
	//BufFace bufFace;
} InPiece;

typedef enum BufVertType {
	STUC_BUF_VERT_IN_OR_MAP,
	STUC_BUF_VERT_ON_EDGE,
	STUC_BUF_VERT_OVERLAP,
	STUC_BUF_VERT_INTERSECT,
	STUC_BUF_VERT_SUB_TYPE_IN,
	STUC_BUF_VERT_SUB_TYPE_MAP,
	STUC_BUF_VERT_SUB_TYPE_EDGE_IN,
	STUC_BUF_VERT_SUB_TYPE_EDGE_MAP
} BufVertType;


typedef struct MapVert {
	I8 type;// BufVertType
	I8 mapCorner;
	I32 inFace;
} MapVert;

typedef struct InVert {
	I8 type;
	I8 inCorner;
	I8 tri;
	I32 inFace;
} InVert;

typedef union InOrMapVert {
	InVert in;
	MapVert map;
} InOrMapVert;

typedef struct EdgeMapVert {
	I8 type;// BufVertType
	I8 inCorner;
	I8 mapCorner;
	U32 inFace;
	F32 tInEdge;
} EdgeMapVert;

typedef struct EdgeInVert {
	I8 type;
	I8 mapCorner;
	I8 inCorner;
	U32 inFace;
	F32 tMapEdge;
} EdgeInVert;

typedef union BufVertOnEdge {
	EdgeInVert in;//currently unused
	EdgeMapVert map;
} BufVertOnEdge;

typedef struct OverlapVert {
	I32 inFace;
	I8 inCorner;
	I8 mapCorner;
} OverlapVert;

typedef struct BufVertOverlapArr {
	OverlapVert *pArr;
	I32 size;
	I32 count;
} BufVertOverlapArr;

typedef struct IntersectVert {
	V2_F32 pos;
	F32 tInEdge;
	F32 tMapEdge;
	I32 inFace;
	I32 mapCorner;
	I8 inCorner;
	I8 tri;
} IntersectVert;

typedef struct BufVertOnEdgeArr {
	BufVertOnEdge *pArr;
	I32 size;
	I32 count;
} BufVertOnEdgeArr;

typedef struct BufVertInOrMapArr {
	InOrMapVert *pArr;
	I32 size;
	I32 count;
} BufVertInOrMapArr;

typedef struct BufVertIntersectArr {
	IntersectVert *pArr;
	I32 size;
	I32 count;
} BufVertIntersectArr;

typedef struct BufCorner {
	U32 vert : 30;
	U32 type : 2;// BufVertType 
} BufCorner;

typedef struct BufCornerArr {
	BufCorner *pArr;
	I32 size;
	I32 count;
} BufCornerArr;

typedef struct BufFaceArr {
	BufFace *pArr;
	I32 size;
	I32 count;
} BufFaceArr;

typedef struct BufMesh {
	BufFaceArr faces;
	BufCornerArr corners;
	BufVertInOrMapArr inOrMapVerts;
	BufVertOnEdgeArr onEdgeVerts;
	BufVertOverlapArr overlapVerts;
	BufVertIntersectArr intersectVerts;
} BufMesh;

typedef struct BufMeshArr {
	BufMesh arr[PIX_THREAD_MAX_SUB_MAPPING_JOBS];
	I32 count;
} BufMeshArr;

typedef struct InPieceArr {
	BufMeshArr *pBufMeshes;
	InPiece *pArr;
	I32 size;
	I32 count;
} InPieceArr;

typedef struct FindEncasedFacesJobArgs {
	JobArgs core;
	HTable encasedFaces;
	InPieceArr inPiecesMono;
} FindEncasedFacesJobArgs;

typedef struct SplitInPiecesAlloc {
	PixalcLinAlloc encased;
	PixalcLinAlloc inFace;
	PixalcLinAlloc border;
} SplitInPiecesAlloc;

typedef struct SplitInPiecesJobArgs {
	JobArgs core;
	const InPieceArr *pInPieceArr;
	InPieceArr newInPieces;
	InPieceArr newInPiecesClip;
	SplitInPiecesAlloc alloc;
} SplitInPiecesJobArgs;

typedef struct InVertKey {
	I32 inVert;
	I32 mapFace;
} InVertKey;

typedef struct MapVertKey {
	I32 mapVert;
	I32 inFace;
} MapVertKey;

typedef struct EdgeInVertKey {
	I32 inVert;
	I32 mapEdge;
} EdgeInVertKey;

typedef struct EdgeMapVertKey {
	I32 mapVert;
	I32 inEdge;
} EdgeMapVertKey;

typedef struct OverlapVertKey {
	I32 inVert;
	I32 mapVert;
} OverlapVertKey;

typedef struct IntersectVertKey {
	I32 inEdge;
	I32 mapEdge;
} IntersectVertKey;

typedef union MergeTableInOrMapKey {
	InVertKey in;
	MapVertKey map;
} MergeTableInOrMapKey;

typedef union MergeTableOnEdgeKey {
	EdgeInVertKey in;
	EdgeMapVertKey map;
} MergeTableOnEdgeKey;

typedef union MergeTableVertKey {
	MergeTableInOrMapKey inOrMap;
	MergeTableOnEdgeKey onEdge;
	OverlapVertKey overlap;
	IntersectVertKey intersect;
} MergeTableVertKey;

typedef struct MergeTableKey {
	BufVertType type;
	V2_I16 tile;
	MergeTableVertKey key;
} MergeTableKey;

typedef union MergeTableInitInfoVert {
	InOrMapVert inOrMap;
	BufVertOnEdge onEdge;
	OverlapVert overlap;
	IntersectVert intersect;
} MergeTableInitInfoVert;

typedef struct VertMergeCorner {
	FaceCorner corner;
	I8 bufMesh;
	bool clipped;
} VertMergeCorner;

typedef struct VertMergeTransform {
	M3x3 tbn;
} VertMergeTransform;

typedef struct VertMerge {
	HTableEntryCore core;
	MergeTableKey key;
	VertMergeCorner bufCorner;
	VertMergeTransform transform;
	U32 cornerCount : 31;
	U32 preserve : 1;
	U32 removed : 1;
	I32 outVert;
	I32 linIdx;
} VertMerge;

typedef struct VertMergeIntersect {
	VertMerge core;
	VertMerge *pSnapTo;
} VertMergeIntersect;

/*
typedef struct VertMergeBucket {
	VertMergeCore *pList;
} VertMergeBucket;

typedef struct VertMergeTable {
	PixalcLinAlloc alloc;
	PixalcLinAlloc intersectAlloc;
	VertMergeBucket *pTable;
	I32 size;
} VertMergeTable;
*/

bool couldInEdgeIntersectMapFace(const Mesh *pInMesh, I32 edge);

typedef struct InPieceKey {
	I32 mapFace;
	V2_I16 tile;
} InPieceKey;

static inline
U64 stucInPieceMakeKey(const void *pKeyData) {
	const InPieceKey *pKey = pKeyData;
	PIX_ERR_ASSERT("tile isn't 16 bit anymore?", sizeof(pKey->tile.d[0]) == 2);
	return
		(U64)pKey->mapFace << 32 |
		(U64)pKey->tile.d[0] << 16 |
		(U64)pKey->tile.d[1];
}

static inline
const U8 *stucGetTri(const FaceTriangulated *pFaceTris, I32 face, I32 idx) {
	return pFaceTris ? pFaceTris[face].pTris + idx * 3 : NULL;
}

typedef struct CachedBc {
	V3_F32 bc;
	bool valid;
} CachedBc;

typedef struct SrcFaces {
	I32 in;
	I32 map;
} SrcFaces;

SrcFaces stucGetSrcFacesForBufCorner(
	const MapToMeshBasic *pBasic,
	const InPiece *pInPiece,
	const BufMesh *pBufMesh,
	FaceCorner corner
);

void stucGetBufMeshForVertMergeEntry(
	const InPieceArr *pInPieces,
	const InPieceArr *pInPiecesClip,
	const VertMerge *pVert,
	const InPiece **ppInPiece,
	const BufMesh **ppBufMesh
);

static inline
I32 stucRangeGetSize(Range range) {
	I32 size = range.end - range.start;
	PIX_ERR_ASSERT("'Range' type doesn't support empty range", size > 0);
	return size;
}

typedef struct TPieceBuf {
	U32 mergedWith : 31;
	U32 merged : 1;
	U32 idx : 31;
	U32 added : 1;
} TPieceBuf;

typedef struct TPieceBufArr {
	TPieceBuf *pArr;
	I32 size;
	I32 count;
} TPieceBufArr;

typedef struct TPieceInFace {
	U32 idx : 29;
	U32 size : 3; //used for in-faces, so max face size of 4
} TPieceInFace;

typedef struct TPieceInFaceArr {
	TPieceInFace *pArr;
	I32 size;
	I32 count;
} TPieceInFaceArr;

typedef struct TPiece {
	 TPieceInFaceArr inFaces;
} TPiece;

typedef struct TPieceArr {
	TPiece *pArr;
	I32 size;
	I32 count;
	I32 *pInFaces;
	I32 faceCount;
} TPieceArr;

typedef struct TPieceVert {
	HTableEntryCore core;
	I32 vert;
	I32 tPiece;
} TPieceVert;

typedef struct TangentJobArgs {
	JobArgs core;
	const TPieceArr *pTPieces;
	I32Arr faces;
	I32 cornerCount;
	V3_F32 *pTangents;
	F32 *pTSigns;
} TangentJobArgs;

StucErr stucBuildTangents(void *pArgsVoid);

void stucThreadPoolSetDefault(StucContext context);
void stucAllocSetCustom(PixalcFPtrs *pAlloc, PixalcFPtrs *pCustomAlloc);
void stucAllocSetDefault(PixalcFPtrs *pAlloc);
