/* 
SPDX-FileCopyrightText: 2025 Caleb Dawson
SPDX-License-Identifier: Apache-2.0
*/

#pragma once

#include <float.h>

#include <pixenals_thread_utils.h>
#include <pixenals_math_utils.h>

#include <uv_stucco_intern.h>

typedef struct BaseTriVerts {
	V3_F32 xyz[4];
	V2_F32 uv[4];
	F32 scale[4];
} BaseTriVerts;

typedef enum InsideStatus {
	STUC_INSIDE_STATUS_NONE,
	STUC_INSIDE_STATUS_OUTSIDE,
	STUC_INSIDE_STATUS_INSIDE,
	STUC_INSIDE_STATUS_ON_LINE,
	STUC_INSIDE_STATUS_ON_VERT
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
	V2_I16 tile,
	HalfPlane *pCache
) {
	V2_F32 fTile = {.d = {(F32)tile.d[0], (F32)tile.d[1]}};
	for (I32 i = 0; i < pInFace->size; ++i) {
		pCache[i] = (HalfPlane){
			.idx = (I8)i,
			.idxNext = (I8)stucGetCornerNext(i, pInFace),
			.idxPrev = (I8)stucGetCornerPrev(i, pInFace),
			.edge = stucGetMeshEdge(
				&pMesh->core,
				(FaceCorner) {.face = pInFace->idx, .corner = i}
			),
			.uv = _(pMesh->pUvs[pInFace->start + i] V2SUB fTile)
		};
		pCache[i].uvNext = _(pMesh->pUvs[pInFace->start + pCache[i].idxNext] V2SUB fTile);
		
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
	V3_F32 (* fpGetPoint)(const Mesh *, const FaceRange *, I32);
	FaceTriangulated *pTris;
	bool *pRemoved;
	I32 triCount;
	V3_F32 normal;
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
	const I32 *pTriIdx,
	const V3_F32 *pTri,
	const V3_F32 *pTriNormal
) {
	V3_F32 normal = pixmV3F32Normalize(*pTriNormal);
	F32 triDistFromOrigin = _(normal V3DOT pTri[0]) * -1.0f;
	for (I32 i = 0; i < pState->pFace->size; ++i) {
		if (i == pTriIdx[0] || i == pTriIdx[1] || i == pTriIdx[2]) {
			continue;
		}
		V3_F32 point = pState->fpGetPoint(pState->pMesh, pState->pFace, i);
		F32 distFromPlane = _(normal V3DOT point) + triDistFromOrigin;
		V3_F32 projPoint = _(point V3ADD _(normal V3MULS (distFromPlane * -1.0f)));
		V3_F32 bc = pixmCartesianToBarycentric(pTri, &projPoint, &normal);
		if (_(bc V3GREAT (V3_F32){0})) {
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
	V3_F32 a = pState->fpGetPoint(pState->pMesh, pState->pFace, cornerPrev);
	V3_F32 b = pState->fpGetPoint(pState->pMesh, pState->pFace, corner);
	V3_F32 c = pState->fpGetPoint(pState->pMesh, pState->pFace, cornerNext);
	V3_F32 ac = _(c V3SUB a);
	V3_F32 cross = _(_(b V3SUB a) V3CROSS ac);
	if (_(cross V3DOT pState->normal) <= 0 || // ear is concave or degenerate
		doesEarIntersectFace(
			pState,
			(I32[]){cornerPrev, corner, cornerNext},
			(V3_F32[]){a, b, c},
			&cross
		)
	) {
		return NULL;
	}
	F32 len = pixmV3F32Len(ac);
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

static inline
void stucRemoveEar(TriangulateState *pState) {
	Ear *pEar = pState->pEarList;
	if (pEar->pNext) {
		pEar->pNext->pPrev = NULL;
	}
	pState->pEarList = pEar->pNext;
}

static inline
bool isMarkedSkip(I32Arr *pSkip, I32 idx) {
	for (I32 i = 0; i < pSkip->count; ++i) {
		if (idx == pSkip->pArr[i]) {
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
	I32 skipArr[32] = {0};
	I32Arr skip = {.pArr = skipArr};
	do {
		I32 lowCorner = 0;
		V2_F32 lowCoord = { FLT_MAX, FLT_MAX };
		for (I32 i = 0; i < pFace->size; ++i) {
			if (isMarkedSkip(&skip, i)) {
				continue;
			}
			V2_F32 pos = fpGetPoint(pMesh, pFace, i);
			if (pos.d[0] > lowCoord.d[0] ||

				pos.d[0] == lowCoord.d[0] &&
				pos.d[1] >= lowCoord.d[1]
			) {
				continue;
			}
			lowCorner = i;
			lowCoord = pos;
		}
		I32 prev = lowCorner == 0 ? pFace->size - 1 : lowCorner - 1;
		I32 next = (lowCorner + 1) % pFace->size;
		V2_F32 a = fpGetPoint(pMesh, pFace, prev);
		V2_F32 b = fpGetPoint(pMesh, pFace, lowCorner);
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
		skip.pArr[skip.count] = lowCorner;
		++skip.count;
	} while(skip.count < pFace->size);
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

static
V3_F32 getTriNormal(
	const Mesh *pMesh,
	const FaceRange *pFace,
	I32Arr *pSkip,
	I32 corner,
	V3_F32 (* fpGetPoint) (const Mesh *, const FaceRange *, I32)
) {
	I32 prev = corner == 0 ? pFace->size - 1 : corner - 1;
	I32 next = (corner + 1) % pFace->size;
	V3_F32 a = fpGetPoint(pMesh, pFace, prev);
	V3_F32 b = fpGetPoint(pMesh, pFace, corner);
	V3_F32 c = fpGetPoint(pMesh, pFace, next);
	V3_F32 cross = _(_(b V3SUB a) V3CROSS _(c V3SUB a));
	if (_(cross V3EQL (V3_F32){0})) {
		//tri is degenerate, find another corner
		pSkip->pArr[pSkip->count] = corner;
		++pSkip->count;
	}
	return cross;
}

static inline
V3_F32 stucCalcFaceNormal(
	const FaceRange *pFace,
	const Mesh *pMesh,
	V3_F32 (* fpGetPoint) (const Mesh *, const FaceRange *, I32)
) {
	PIX_ERR_ASSERT("", pFace->start >= 0 && pFace->size >= 3);
	I32 skipArr[32] = {0};
	I32Arr skip = {.pArr = skipArr};
	do {
		I32 highCorner = 0;
		I32 lowCorner = 0;
		V3_F32 highCoord = { -FLT_MAX, -FLT_MAX, -FLT_MAX };
		V3_F32 lowCoord = { FLT_MAX, FLT_MAX, FLT_MAX };
		for (I32 i = 0; i < pFace->size; ++i) {
			if (isMarkedSkip(&skip, i)) {
				continue;
			}
			V3_F32 pos = fpGetPoint(pMesh, pFace, i);
			if (!(
				pos.d[0] > lowCoord.d[0] ||

				pos.d[0] == lowCoord.d[0] &&
				pos.d[1] > lowCoord.d[1] ||

				pos.d[0] == lowCoord.d[0] &&
				pos.d[1] == lowCoord.d[1] &&
				pos.d[2] >= lowCoord.d[2]
			)) {
				lowCorner = i;
				lowCoord = pos;
			}
			if (!(
				pos.d[0] < highCoord.d[0] ||

				pos.d[0] == highCoord.d[0] &&
				pos.d[1] < highCoord.d[1] ||

				pos.d[0] == highCoord.d[0] &&
				pos.d[1] == highCoord.d[1] &&
				pos.d[2] <= highCoord.d[2]
			)) {
				highCorner = i;
				highCoord = pos;
			}
		}
		V3_F32 lowNormal = getTriNormal(pMesh, pFace, &skip, lowCorner, fpGetPoint);
		if (_(lowNormal V3EQL (V3_F32){0})) {
			continue;
		}
		V3_F32 highNormal = getTriNormal(pMesh, pFace, &skip, highCorner, fpGetPoint);
		if (_(highNormal V3EQL (V3_F32){0})) {
			continue;
		}
		if (_(lowNormal V3EQL highNormal)) {
			return lowNormal;
		}
		return _(
			_(pixmV3F32Normalize(highNormal) V3ADD pixmV3F32Normalize(lowNormal)) V3DIVS
			2.0f
		);
	} while(skip.count < pFace->size);
	PIX_ERR_ASSERT("", false);
	return (V3_F32){0};
}

STUC_FORCE_INLINE
void stucTriangulateFace(
	const StucAlloc *pAlloc,
	const FaceRange *pFace,
	const Mesh *pMesh,
	V3_F32 (* fpGetPoint)(const Mesh *, const FaceRange *, I32),
	FaceTriangulated *pTris
) {
	PIX_ERR_ASSERT("", pTris->pTris);
	TriangulateState state = {
		.pMesh = pMesh,
		.fpGetPoint = fpGetPoint,
		.pFace = pFace,
		.pTris = pTris,
		.pRemoved = pAlloc->fpCalloc(pFace->size, sizeof(bool)),
		.normal = stucCalcFaceNormal(pFace, pMesh, fpGetPoint)
	};

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
	stucTriangulateFace(pAlloc, pFace, pMesh, stucGetVertPos, pTris);
}

STUC_FORCE_INLINE
V3_F32 stucGetBarycentricInTri(
	const Mesh *pMesh,
	const FaceRange *pFace,
	V3_F32 (* fpGetPoint)(const Mesh *, const FaceRange *, I32),
	const U8 *pTriCorners,
	V2_F32 vert
) {
	V3_F32 tri[3] = {0};
	for (I32 i = 0; i < 3; ++i) {
		tri[i] = fpGetPoint(pMesh, pFace, pTriCorners[i]);
	}
	return pixmCartesianToBarycentric(
		tri,
		&(V3_F32){.d = {vert.d[0], vert.d[1]}},
		&(V3_F32){.d = {.0f, .0f, 1.0f}}
	);

}

static inline
V3_F32 stucGetBarycentricInTriFromVerts(
	const Mesh *pMesh,
	const FaceRange *pFace,
	const U8 *pTriCorners,
	V2_F32 vert
) {
	return stucGetBarycentricInTri(pMesh, pFace, stucGetVertPos, pTriCorners, vert);
}

//Caller must check for nan in return value
STUC_FORCE_INLINE
V3_F32 stucGetBarycentricInFace(
	const Mesh *pMesh,
	const FaceRange *pFace,
	V2_I16 tile,
	V3_F32 (* fpGetPoint)(const Mesh *, const FaceRange *, I32),
	I8 *pTriCorners,
	V2_F32 vertV2
) {
	PIX_ERR_ASSERT("", pixmV2F32IsFinite(vertV2));
	PIX_ERR_ASSERT("", (pFace->size == 3 || pFace->size == 4) && pTriCorners);
	V3_F32 vert = {.d = {vertV2.d[0], vertV2.d[1]}};
	V3_F32 fTile = {.d = {(F32)tile.d[0], (F32)tile.d[1]}};
	V3_F32 triA[3] = {0};
	for (I32 i = 0; i < 3; ++i) {
		triA[i] = _(fpGetPoint(pMesh, pFace, i) V3SUB fTile);
	}
	V3_F32 up = {.d = {.0f, .0f, 1.0f}};
	V3_F32 vertBc = pixmCartesianToBarycentric(triA, &vert, &up);
	if (pFace->size == 4 && pixmV3F32IsFinite(vertBc) && vertBc.d[1] < 0) {
		//base face is a quad, and vert is outside first tri,
		//so use the second tri
		
		V3_F32 triB[3] = {
			triA[2],
			_(fpGetPoint(pMesh, pFace, 3) V3SUB fTile),
			triA[0]
		};
		vertBc = pixmCartesianToBarycentric(triB, &vert, &up);
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
		(V2_I16) {0},
		stucGetVertPos,
		pTriCorners,
		vert
	);
}

static inline
V3_F32 stucGetBarycentricInFaceFromUvs(
	const Mesh *pMesh,
	const FaceRange *pFace,
	V2_I16 tile,
	I8 *pTriCorners,
	V2_F32 vert
) {
	return stucGetBarycentricInFace(
		pMesh,
		pFace,
		tile,
		stucGetUvPosAsV3,
		pTriCorners,
		vert
	);
}

StucErr stucBuildEdgeList(StucContext pCtx, StucMesh *pMesh);
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
	StucContext pCtx,
	I32 jobCount, void *pJobArgs, I32 argStructSize,
	StucErr (* func)(void *)
);
//U32 stucGetEncasedFaceHash(I32 mapFace, V2_I16 tile, I32 tableSize);

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

typedef struct InPieceKey {
	I32 mapFace;
	V2_I16 tile;
} InPieceKey;

static inline
PixuctKey stucInPieceMakeKey(const void *pKeyData) {
	return (PixuctKey){.pKey = pKeyData, .size = sizeof(InPieceKey)};
}

static inline
const U8 *stucGetTri(const FaceTriangulated *pFaceTris, I32 face, I32 idx) {
	return pFaceTris ? pFaceTris[face].pTris + idx * 3 : NULL;
}

typedef struct CachedBc {
	V3_F32 bc;
	bool valid;
} CachedBc;

static inline
I32 stucRangeGetSize(Range range) {
	I32 size = range.end - range.start;
	PIX_ERR_ASSERT("'Range' type doesn't support empty range", size > 0);
	return size;
}

void stucThreadPoolSetDefault(StucContext context);
void stucAllocSetCustom(PixalcFPtrs *pAlloc, PixalcFPtrs *pCustomAlloc);
void stucAllocSetDefault(PixalcFPtrs *pAlloc);
