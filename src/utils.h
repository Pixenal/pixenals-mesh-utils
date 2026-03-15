/* 
SPDX-FileCopyrightText: 2025 Caleb Dawson
SPDX-License-Identifier: Apache-2.0
*/

#pragma once

#include <float.h>

#include <pixenals_thread_utils.h>
#include <pixenals_math_utils.h>
#include <pixenals_structs.h>

#include <uv_stucco_intern.h>

//TODO replace Clust prefix with Stuc
typedef struct ClustFaceCorner {
	I32 face;
	I32 corner;
} ClustFaceCorner;

typedef struct ClustBorderBuf {
	BorderArr arr;
} ClustBorderBuf;

typedef struct ClustSplitFaceIdx {
	PixuctHTableEntryCore core;
	I32 face;
	bool removed;
	bool pendingRemove;
	bool preserve[4];
} ClustSplitFaceIdx;

typedef struct ClustSplitFaceBuf {
	ClustSplitFaceIdx **ppArr;
	I32 size;
	I32 count;
} ClustSplitFaceBuf;

typedef struct ClustBorderEdgeTableEntry {
	PixuctHTableEntryCore core;
	ClustFaceCorner corner;
	bool checked;
} ClustBorderEdgeTableEntry;

typedef struct SplitMesh {
	void *pUserData;
	PixtyRange (*fpFaceRange)(const void *, int32_t);
	int32_t (*fpVert)(const void *, int32_t);
	PixtyV2_F32 (*fpPos)(const void *, int32_t);
	ClustFaceCorner (*fpAdjCorner)(const void *, ClustFaceCorner);
	int32_t (*fpEdge)(const void *, ClustFaceCorner);
	int32_t faceCount;
} SplitMesh;

typedef struct ClustSplitCallbacks {
	PixErr (*fpIdxTableBuild)(
		void *,
		void *,
		void (*)(void *, const PixtyI32Arr *)
	);
	bool (*fpSplitPredicate)(const void *, PixtyRange, ClustFaceCorner);
	PixErr (*fpBorderAdd)(void *, const ClustBorderBuf *);
	PixErr (*fpFacesAdd)(void *, const ClustSplitFaceBuf *);
	bool (* fpIsEdgeIntern)(const void *, const ClustSplitFaceIdx *, int32_t);
	bool (* fpBorderPredicate)(const void *, int32_t);
} ClustSplitCallbacks;

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

typedef struct HalfPlane {
	V2_F32 uv;
	V2_F32 dirUnit;
	F32 len;
	I32 edge;
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
			.edge = stucGetMeshEdge(
				&pMesh->core,
				(FaceCorner) {.face = pInFace->idx, .corner = i}
			),
			.uv = _(pMesh->pUvs[pInFace->start + i] V2SUB fTile)
		};
	}
	for (I32 i = 0; i < pInFace->size; ++i) {
		I32 iNext = stucGetCornerNext(i, pInFace);
		V2_F32 dir = _(pCache[iNext].uv V2SUB pCache[i].uv);
		pCache[i].len = pixmV2F32Len(dir);
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
	I32 cornerPrev;
	I32 corner;
	I32 cornerNext;
	F32 len;
} Ear;

typedef struct TriangulateState {
	PixalcLinAlloc earAlloc;
	Ear *pEarList;
	const Mesh *pMesh;
	const FaceRange *pFace;
	V3_F32 (* fpGetPoint)(const Mesh *, const FaceRange *, I32);
	I8 *pRemoved;
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
	pNewEar->cornerPrev = cornerPrev;
	pNewEar->corner = corner;
	pNewEar->cornerNext = cornerNext;
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
Ear *addEar(TriangulateState *pState, I32 *pCount, U8 *pTris) {
	Ear *pEar = pState->pEarList;
	U8 *pTri = pTris + *pCount * 3;
	I32 cornerPrev = stucGetPrevRemaining(pState, pEar->corner, pState->pFace);
	I32 cornerNext = stucGetNextRemaining(pState, pEar->corner, pState->pFace);
	if (cornerPrev != pEar->cornerPrev || cornerNext != pEar->cornerNext) {
		return NULL;//ear entry is outdated
	}
	pTri[0] = cornerPrev;
	pTri[1] = pEar->corner;
	pTri[2] = cornerNext;
	++*pCount;
	
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

STUC_FORCE_INLINE
I32 stucGetNonDegenBoundCorner(
	const FaceRange *pFace,
	const Mesh *pMesh,
	V2_F32 (* fpGetPoint) (const Mesh *, const FaceRange *, I32),
	bool useMin,
	I32Arr *pExternSkip,
	F32 *pDet
) {
	PIX_ERR_ASSERT("", pFace->start >= 0 && pFace->size >= 3);
	I32 skipArr[STUC_NGON_MAX_SIZE] = {0};
	I32Arr skip = {.pArr = skipArr};
	do {
		I32 corner = 0;
		V2_F32 boundPos = {FLT_MAX, FLT_MAX};
		boundPos = useMin ? boundPos : _(boundPos V2MULS -1.0f);
		for (I32 i = 0; i < pFace->size; ++i) {
			if (isMarkedSkip(&skip, i) || pExternSkip && isMarkedSkip(pExternSkip, i)) {
				continue;
			}
			V2_F32 pos = fpGetPoint(pMesh, pFace, i);
			if (useMin) {
				if (pos.d[0] > boundPos.d[0] ||

					pos.d[0] == boundPos.d[0] &&
					pos.d[1] >= boundPos.d[1]
				) {
					continue;
				}
			}
			else {
				if (pos.d[0] < boundPos.d[0] ||

					pos.d[0] == boundPos.d[0] &&
					pos.d[1] <= boundPos.d[1]
				) {
					continue;
				}
			}
			corner = i;
			boundPos = pos;
		}
		I32 prev = corner == 0 ? pFace->size - 1 : corner - 1;
		I32 next = (corner + 1) % pFace->size;
		V2_F32 a = fpGetPoint(pMesh, pFace, prev);
		V2_F32 b = fpGetPoint(pMesh, pFace, corner);
		V2_F32 c = fpGetPoint(pMesh, pFace, next);
		//alt formula for determinate,
		//shorter and less likely to cause numerical error
		F32 det =
			(b.d[0] - a.d[0]) * (c.d[1] - a.d[1]) -
			(c.d[0] - a.d[0]) * (b.d[1] - a.d[1]);
		if (det) {
			if (pDet) {
				*pDet = det;
			}
			return corner;
		}
		//abc is degenerate, find another corner
		skip.pArr[skip.count] = corner;
		++skip.count;
	} while(skip.count < pFace->size);
	return -1;
}

//finds corner on convex hull of face, & determines wind direction from that
//returns 0 for clockwise, 1 for counterclockwise, & 2 if degenerate
STUC_FORCE_INLINE
I32 stucCalcFaceWind(
	const FaceRange *pFace,
	const Mesh *pMesh,
	V2_F32 (* fpGetPoint) (const Mesh *, const FaceRange *, I32)
) {
	F32 det = .0f;
	I32 corner = stucGetNonDegenBoundCorner(pFace, pMesh, fpGetPoint, true, NULL, &det);
	return corner != -1 ? det > .0f : 2;
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
	I32 corner,
	V3_F32 (* fpGetPoint) (const Mesh *, const FaceRange *, I32)
) {
	I32 prev = corner == 0 ? pFace->size - 1 : corner - 1;
	I32 next = (corner + 1) % pFace->size;
	V3_F32 a = fpGetPoint(pMesh, pFace, prev);
	V3_F32 b = fpGetPoint(pMesh, pFace, corner);
	V3_F32 c = fpGetPoint(pMesh, pFace, next);
	return _(_(b V3SUB a) V3CROSS _(c V3SUB a));
}

typedef struct AxisBounds{
	I32 minIdx;
	I32 maxIdx;
	F32 min;
	F32 max;
	F32 len;
} AxisBounds;

static inline
void axisBoundsCalcLen(AxisBounds *pBounds) {
	pBounds->len = pBounds->max - pBounds->min;
	PIX_ERR_ASSERT("", pBounds->len >= .0f);
}

static inline
void axisBoundsCmp(AxisBounds *pBounds, float pos, I32 idx) {
	if (pos < pBounds->min) {
		pBounds->min = pos;
		pBounds->minIdx = idx;
	}
	if (pos > pBounds->max) {
		pBounds->max = pos;
		pBounds->maxIdx = idx;
	}
}

//returns the longest axis
static inline
I32 axisBoundsMake(
	const FaceRange *pFace,
	const Mesh *pMesh,
	V3_F32 (* fpGetPoint) (const Mesh *, const FaceRange *, I32),
	I32Arr *pSkip,
	AxisBounds *pBounds
) {
	for (I32 i = 0; i < 3; ++i) {
		pBounds[i].max = -FLT_MAX;
		pBounds[i].min = FLT_MAX;
	}
	for (I32 i = 0; i < pFace->size; ++i) {
		if (pSkip && isMarkedSkip(pSkip, i)) {
			continue;
		}
		V3_F32 pos = fpGetPoint(pMesh, pFace, i);
		axisBoundsCmp(pBounds + 0, pos.d[0], i);
		axisBoundsCmp(pBounds + 1, pos.d[1], i);
		axisBoundsCmp(pBounds + 2, pos.d[2], i);
	}
	axisBoundsCalcLen(pBounds + 0);
	axisBoundsCalcLen(pBounds + 1);
	axisBoundsCalcLen(pBounds + 2);
	I32 lowAxis = pBounds[0].len > pBounds[1].len;
	if (pBounds[2].len < pBounds[1].len && pBounds[2].len < pBounds[0].len) {
		lowAxis = 2;
	}
	return lowAxis;
}

static inline
void markSkip(I32Arr *pSkip, I32 idx) {
	PIX_ERR_ASSERT("", pSkip->count < STUC_NGON_MAX_SIZE);
	pSkip->pArr[pSkip->count] = idx;
	++pSkip->count;
}

static inline
V3_F32 stucCalcFaceNormal(
	const FaceRange *pFace,
	const Mesh *pMesh,
	V3_F32 (* fpGetPoint) (const Mesh *, const FaceRange *, I32)
) {
	PIX_ERR_ASSERT(
		"invalid face size",
		pFace->start >= 0 && pFace->size >= 3 && pFace->size <= STUC_NGON_MAX_SIZE
	);
	I32 skipArr[STUC_NGON_MAX_SIZE] = {0};
	I32Arr skip = {.pArr = skipArr};
	AxisBounds bounds[3] = {0};
	I32 axis = axisBoundsMake(pFace, pMesh, fpGetPoint, NULL, bounds);
	V2_F32 (*fpPos)(const Mesh *, const FaceRange *, I32) =
		axis == 2 ? stucVertPosXy : axis ? stucVertPosXz : stucVertPosYz;
	do {
		I32 minIdx = stucGetNonDegenBoundCorner(pFace, pMesh, fpPos, true, &skip, NULL);
		I32 maxIdx = stucGetNonDegenBoundCorner(pFace, pMesh, fpPos, false, &skip, NULL);
		if (minIdx == -1 || maxIdx == -1) {
			return (V3_F32){0};
		}
		V3_F32 minNormal = getTriNormal(pMesh, pFace, minIdx, fpGetPoint);
		V3_F32 maxNormal = getTriNormal(pMesh, pFace, maxIdx, fpGetPoint);
		if (_(minNormal V3DOT maxNormal) <= .0f) {
			markSkip(&skip, minIdx);
			markSkip(&skip, maxIdx);
			continue;
		}
		if (_(minNormal V3EQL maxNormal)) {
			return minNormal;
		}
		return _(
			_(pixmV3F32Normalize(maxNormal) V3ADD pixmV3F32Normalize(minNormal)) V3DIVS
			2.0f
		);
	} while(skip.count < pFace->size);
	return (V3_F32){0};
}

//returns tri count (may be less than size - 2 if face is degen)
STUC_FORCE_INLINE
I32 stucTriangulateFace(
	const StucAlloc *pAlloc,
	const FaceRange *pFace,
	const Mesh *pMesh,
	V3_F32 (* fpGetPoint)(const Mesh *, const FaceRange *, I32),
	U8 *pTris
) {
	PIX_ERR_ASSERT("", pTris);
	TriangulateState state = {
		.pMesh = pMesh,
		.fpGetPoint = fpGetPoint,
		.pFace = pFace,
		.pRemoved = pAlloc->fpCalloc(pFace->size, 1),
		.normal = stucCalcFaceNormal(pFace, pMesh, fpGetPoint)
	};
	if (_(state.normal V3EQL (V3_F32){0})) {
		return 0;
	}

	pixalcLinAllocInit(pAlloc, &state.earAlloc, sizeof(Ear), pFace->size, true);

	//add initial ears
	for (I32 i = 0; i < pFace->size; ++i) {
		addEarCandidate(&state, i);
	}
	I32 triCount = 0;
	while (state.pEarList) {
		Ear *pAddedEar = NULL;
		if (!state.pRemoved[state.pEarList->cornerPrev] &&
			!state.pRemoved[state.pEarList->corner] &&
			!state.pRemoved[state.pEarList->cornerNext]
		) {
			pAddedEar = addEar(&state, &triCount, pTris);
		}
		stucRemoveEar(&state);
		if (pAddedEar) {
			addAdjEarCandidates(&state, pAddedEar);
		}
	}
	pixalcLinAllocDestroy(&state.earAlloc);
	pAlloc->fpFree(state.pRemoved);
	PIX_ERR_ASSERT("", triCount <= pFace->size - 2);
	return triCount;
}

static inline
I32 stucTriangulateFaceFromVerts(
	const StucAlloc *pAlloc,
	const FaceRange *pFace,
	const Mesh *pMesh,
	U8 *pTris
) {
	return stucTriangulateFace(pAlloc, pFace, pMesh, stucGetVertPos, pTris);
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
M3x3 stucGetInterpolatedTbn(
	const Mesh *pMesh,
	const FaceRange *pFace,
	const I8 *pTriCorners,
	V3_F32 bc
);

StucErr stucDoJobInParallel(
	StucContext pCtx,
	I32 threadId,
	I32 jobCount, void *pJobArgs, I32 argStructSize,
	StucErr (* func)(void *)
);

typedef struct InPieceKey {
	I32 cluster;
	V2_I16 tile;
	bool clip;
} InPieceKey;

static inline
PixuctKey stucInPieceMakeKey(const void *pKeyData) {
	return (PixuctKey){.pKey = pKeyData, .size = sizeof(InPieceKey)};
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

static inline
I32 stucCouldInEdgeIntersectMapFace(const Mesh *pMesh, I32 edge) {
	bool preserve = stucGetIfPreserveEdge(pMesh, edge);
	bool ret = stucGetIfSeamEdge(pMesh, edge) || stucGetIfMatBorderEdge(pMesh, edge);
	return preserve && !ret ? 2 : preserve || ret;
}

static inline
PixtyV2_F32 stucClustPos(const void *pMeshRaw, I32 vert) {
	const Mesh *pMesh = pMeshRaw;
	PIX_ERR_ASSERT("", pMesh->pPos && vert >= 0 && vert < pMesh->core.vertCount);
	return *(PixtyV2_F32 *)&pMesh->pPos[vert];
}

static inline
I32 stucClustVert(const void *pMeshRaw, I32 corner) {
	const Mesh *pMesh = pMeshRaw;
	PIX_ERR_ASSERT("", corner >= 0 && corner < pMesh->core.cornerCount);
	return pMesh->core.pCorners[corner];
}

static inline
PixtyRange stucClustFaceRange(const void *pMeshRaw, I32 face) {
	const Mesh *pMesh = pMeshRaw;
	PIX_ERR_ASSERT("", face >= 0 && face < pMesh->core.faceCount);
	return (PixtyRange) {
		.start = pMesh->core.pFaces[face],
		.end = pMesh->core.pFaces[face + 1]
	};
}

#if false
void clustBuildFaceIdxTable(void *pTableRaw, const PixtyI32Arr *pFaces);
SearchResult clustFaceIdxTableGet(PixuctHTable *pTable, int32_t face, void **ppEntry);
ClustBorderEdgeTableEntry *clustBorderEdgeAddOrGet(
	PixuctHTable *pBorderTable,
	ClustFaceCorner corner,
	bool add
);

CLUST_FORCE_INLINE
ClustSplitFaceIdx *clustGetAdjFace(
	const SplitMesh *pMesh,
	bool (* fpIsEdgeIntern)(const void *, const ClustSplitFaceIdx *, int32_t),
	PixuctHTable *pIdxTable,
	ClustFaceCorner corner,
	int32_t *pAdjCorner
) {
	ClustFaceCorner adj = pMesh->fpAdjCorner(pMesh->pUserData, corner);
	if (adj.corner == -1) {
		return NULL;
	}
	PIX_ERR_ASSERT("", adj.corner >= 0);
	ClustSplitFaceIdx *pAdjIdxEntry = NULL;
	clustFaceIdxTableGet(pIdxTable, adj.face, (void **)&pAdjIdxEntry);
	if (!pAdjIdxEntry) {
		if (pAdjCorner) {
			*pAdjCorner = -1;
		}
		return NULL;
	}
	int32_t edge = pMesh->fpEdge(pMesh->pUserData, corner);
	if (pAdjIdxEntry->removed || fpIsEdgeIntern(pMesh->pUserData, pAdjIdxEntry, edge)) {
		if (pAdjCorner) {
			*pAdjCorner = -1;
		}
		return NULL;
	}
	if (pAdjCorner) {
		PixtyRange adjFace = pMesh->fpFaceRange(pMesh->pUserData, pAdjIdxEntry->face);
		adj.corner = (adj.corner + 1) % (adjFace.end - adjFace.start);
		if (pAdjIdxEntry->preserve[adj.corner]) {
			if (pAdjCorner) {
				*pAdjCorner = -1;
			}
			return NULL;
		}
		*pAdjCorner = adj.corner;
	}
	return pAdjIdxEntry;
}

CLUST_FORCE_INLINE
void clustAddAdjFaces(
	const SplitMesh *pMesh,
	void *pUserData,
	const ClustSplitCallbacks *pCallbacks,
	ClustSplitFaceBuf *pInFaceBuf,
	PixuctHTable *pIdxTable,
	PixuctHTable *pBorderEdges,
	ClustSplitFaceIdx *pFace
) {
	PixtyRange inFace = pMesh->fpFaceRange(pMesh->pUserData, pFace->face);
	int32_t faceSize = inFace.end - inFace.start;
	for (int32_t i = 0; i < faceSize; ++i) {
		ClustFaceCorner corner = {.face = pFace->face, .corner = i};
		int32_t adjCorner = -1;
		ClustSplitFaceIdx *pAdjFace = clustGetAdjFace(
			pMesh,
			pCallbacks->fpIsEdgeIntern,
			pIdxTable,
			corner,
			&adjCorner
		);
		if (!pAdjFace) {
			clustBorderEdgeAddOrGet(pBorderEdges, corner, true);
			continue;
		}
		else if (pAdjFace->pendingRemove) {
			//already added to this piece
			continue; 
		}
		else if (pCallbacks->fpSplitPredicate(pUserData, inFace, corner)) {
			pFace->preserve[i] = true;
			PIX_ERR_ASSERT("", adjCorner != -1);
			pAdjFace->preserve[adjCorner] = true;
			clustBorderEdgeAddOrGet(pBorderEdges, corner, true);
			continue;
		}

		pAdjFace->pendingRemove = true;
		PIX_ERR_ASSERT("", pInFaceBuf->count < pInFaceBuf->size);
		pInFaceBuf->ppArr[pInFaceBuf->count] = pAdjFace;
		pInFaceBuf->count++;
	}
}

static inline
ClustSplitFaceIdx *clustGetFirstRemainingFace(PixuctHTable *pIdxTable) {
	PixalcLinAlloc *pTableAlloc = pixuctHTableAllocGet(pIdxTable, 0);
	PixalcLinAllocIter iter = {0};
	pixalcLinAllocIterInit(pTableAlloc, (PixtyRange) { 0, INT32_MAX }, &iter);
	for (; !pixalcLinAllocIterAtEnd(&iter); pixalcLinAllocIterInc(&iter)) {
		ClustSplitFaceIdx *pEntry = pixalcLinAllocGetItem(&iter);
		PIX_ERR_ASSERT("", pEntry);
		if (!pEntry->removed) {
			return pEntry;
		}
	}
	PIX_ERR_ASSERT("this func shouldn't have been called if no faces remained", false);
	return NULL;
}

static inline
void clustAddBorderToArr(const PixalcFPtrs *pAlloc, BorderArr *pArr, Border border) {
	PIX_ERR_ASSERT("", pArr->count <= pArr->size);
	if (!pArr->size) {
		pArr->size = 2;
		pArr->pArr = pAlloc->fpMalloc(pArr->size * sizeof(Border));
	}
	else if (pArr->count == pArr->size) {
		pArr->size *= 2;
		pArr->pArr = pAlloc->fpRealloc(pArr->pArr, pArr->size * sizeof(Border));
	}
	pArr->pArr[pArr->count] = border;
	pArr->count++;
}

CLUST_FORCE_INLINE
bool clustFindAndAddBorder(
	const PixalcFPtrs *pAlloc,
	const SplitMesh *pMesh,
	bool (* fpIsEdgeIntern)(const void *, const ClustSplitFaceIdx *, int32_t),
	ClustBorderBuf *pBorderBuf,
	PixuctHTable *pIdxTable,
	PixuctHTable *pBorderEdgeTable,
	int32_t edgesMax,
	ClustBorderEdgeTableEntry *pStart
) {
	Border border = {.start = pStart->corner};
	ClustFaceCorner corner = pStart->corner;
	ClustBorderEdgeTableEntry *pEntry  = pStart;
	do {
		if (border.len != 0) {//dont run this on first edge
			if (
				corner.face == pStart->corner.face &&
				corner.corner == pStart->corner.corner
			) {
				break;//full loop
			}
			pEntry = clustBorderEdgeAddOrGet(pBorderEdgeTable, corner, false);
		}
		if (pEntry) {
			PIX_ERR_ASSERT("", pEntry->checked == false);
			pEntry->checked = true;
			border.len++;
		}
		PIX_ERR_ASSERT("", border.len <= edgesMax);
		int32_t adjCorner = 0;
		//this is using the table for the pre-split piece.
		//this is fine, as faces arn't marked removed until the end of this func
		const ClustSplitFaceIdx *pAdjFace = clustGetAdjFace(
			pMesh,
			fpIsEdgeIntern,
			pIdxTable,
			corner,
			&adjCorner
		);
		PIX_ERR_ASSERT(
			"if edge isn't in border arr, there should be an adj face",
			!pEntry ^ !pAdjFace
		);
		if (pAdjFace) {
			//edge is internal, move to next adjacent
			corner.face = pAdjFace->face;
			corner.corner = adjCorner;
			//wind = pAdjFace->pInFace->wind;
		}
		else {
			PixtyRange face = pMesh->fpFaceRange(pMesh->pUserData, corner.face);
			corner.corner = (corner.corner + 1) % (face.end - face.start);
		}
	} while(true);
	clustAddBorderToArr(pAlloc, &pBorderBuf->arr, border);
	return true;
}

CLUST_FORCE_INLINE
void clustFillBorderBuf(
	const PixalcFPtrs *pAlloc,
	const SplitMesh *pMesh,
	const ClustSplitCallbacks *pCallbacks,
	ClustBorderBuf *pBorderBuf,
	PixuctHTable *pIdxTable,
	PixuctHTable *pBorderEdges
) {
	PixalcLinAlloc *pTableAlloc = pixuctHTableAllocGet(pBorderEdges, 0);
	pBorderBuf->arr.count = 0;
	int32_t edgeCount = pixalcLinAllocGetCount(pTableAlloc);
	PixalcLinAllocIter iter = {0};
	pixalcLinAllocIterInit(pTableAlloc, (PixtyRange) { 0, INT32_MAX }, &iter);
	for (; !pixalcLinAllocIterAtEnd(&iter); pixalcLinAllocIterInc(&iter)) {
		ClustBorderEdgeTableEntry *pEntry = pixalcLinAllocGetItem(&iter);
		if (pEntry->checked) {
			continue;
		}
		int32_t adjCorner = 0;
		//this is using the table for the pre-split piece.
		//this is fine, as faces arn't marked removed until the end of this func
		if (clustGetAdjFace(
			pMesh,
			pCallbacks->fpIsEdgeIntern,
			pIdxTable,
			pEntry->corner,
			&adjCorner
		)) {
			continue;//we want an exterior edge to start, so skip this one
		}
		int32_t edge = pMesh->fpEdge(pMesh->pUserData, pEntry->corner);
		if (pCallbacks->fpBorderPredicate(pMesh->pUserData, edge)) {
			clustFindAndAddBorder(
				pAlloc,
				pMesh,
				pCallbacks->fpIsEdgeIntern,
				pBorderBuf,
				pIdxTable,
				pBorderEdges,
				edgeCount,
				pEntry
			);
		}
	}
}

CLUST_FORCE_INLINE
PixErr clustSplitAdjFaces(
	const PixalcFPtrs *pAlloc,
	const SplitMesh *pMesh,
	void *pUserData,
	const ClustSplitCallbacks *pCallbacks,
	ClustSplitFaceBuf *pInFaceBuf,
	ClustBorderBuf *pBorderBuf,
	PixuctHTable *pIdxTable,
	int32_t *pFacesRemaining
) {
	PixErr err = PIX_ERR_SUCCESS;
	PixuctHTable borderEdges = {0};
	pixuctHTableInit(
		pAlloc,
		&borderEdges,
		*pFacesRemaining / 2 + 1,
		(PixtyI32Arr) {
			.pArr = (int32_t[]) {sizeof(ClustBorderEdgeTableEntry)},
			.count = 1
		},
		NULL,
		NULL,
		true
	);
	pInFaceBuf->count = 0;
	{
		ClustSplitFaceIdx *pStartFace = clustGetFirstRemainingFace(pIdxTable);
		pInFaceBuf->ppArr[0] = pStartFace;
		pStartFace->pendingRemove = true;
		pInFaceBuf->count++;
		int32_t i = 0;
		do {
			ClustSplitFaceIdx *pIdxEntry = NULL;
			clustFaceIdxTableGet(
				pIdxTable,
				pInFaceBuf->ppArr[i]->face,
				(void **)&pIdxEntry
			);
			clustAddAdjFaces(
				pMesh,
				pUserData,
				pCallbacks,
				pInFaceBuf,
				pIdxTable,
				&borderEdges,
				pIdxEntry
			);
		} while (i++, i < pInFaceBuf->count);
	}
	clustFillBorderBuf(pAlloc, pMesh, pCallbacks, pBorderBuf, pIdxTable, &borderEdges);
	if (pBorderBuf->arr.count) {
		err = pCallbacks->fpBorderAdd(pUserData, pBorderBuf);
		PIX_ERR_RETURN_IFNOT(err, "");
	}
	pixuctHTableDestroy(&borderEdges);
	
	// copy buf into new in-piece, & mark in-faces as removed in idx-table
	err = pCallbacks->fpFacesAdd(pUserData, pInFaceBuf);
	PIX_ERR_RETURN_IFNOT(err, "");
	for (int32_t i = 0; i < pInFaceBuf->count; ++i) {
		PIX_ERR_ASSERT("", pInFaceBuf->ppArr[i]->pendingRemove);
		pInFaceBuf->ppArr[i]->removed = true;
		pInFaceBuf->ppArr[i]->pendingRemove = false;
	}
	*pFacesRemaining -= pInFaceBuf->count;
	return err;
}

STUC_FORCE_INLINE
PixErr clustSplitIslands(
	const PixalcFPtrs *pAlloc,
	const SplitMesh *pMesh,
	void *pUserData,
	const ClustSplitCallbacks *pCallbacks,
	int32_t faceCount,
	ClustSplitFaceBuf *pInFaceBuf,
	ClustBorderBuf *pBorderBuf
) {
	PixErr err = PIX_ERR_SUCCESS;
	PixuctHTable idxTable = {0};
	pixuctHTableInit(
		pAlloc,
		&idxTable,
		faceCount / 4 + 1,
		(PixtyI32Arr) {.pArr = (int32_t[]) {sizeof(ClustSplitFaceIdx)}, .count = 1},
		NULL,
		NULL,
		true
	);
	err = pCallbacks->fpIdxTableBuild(pUserData, &idxTable, clustBuildFaceIdxTable);
	PIX_ERR_THROW_IFNOT(err, "", 0);
	int32_t facesRemaining = faceCount;
	do {
		err = clustSplitAdjFaces(
			pAlloc,
			pMesh,
			pUserData,
			pCallbacks,
			pInFaceBuf,
			pBorderBuf,
			&idxTable,
			&facesRemaining
		);
		PIX_ERR_THROW_IFNOT(err, "", 0);
		PIX_ERR_ASSERT("", facesRemaining >= 0 && facesRemaining < faceCount);
	} while(facesRemaining);
	PIX_ERR_CATCH(0, err, ;);
	pixuctHTableDestroy(&idxTable);
	return err;
}
#endif

typedef struct StucIdxTable {
	U32 idx : 31;
	U32 valid : 1;
} StucIdxTable;

typedef struct StucIdxTableArr {
	StucIdxTable *pArr;
	I32 size;
	I32 count;
} StucIdxTableArr;

struct StucBorderNode;

typedef struct StucBorderLink {
	StucBorderNode *pNode;
} StucBorderLink;

typedef struct StucBorderNode {
	FaceCorner corners[2];
	I32 idx;
	StucIdxTable seen[2];
	bool intern;
} StucBorderNode;

typedef struct StucBorderNodeArr {
	StucBorderNode *pArr;
	I32 size;
	I32 count;
} StucBorderNodeArr;

typedef struct StucBufIsland {
	struct StucBufIsland *pNext;
	//PixtyI32Arr faces;
	I32 idx;
} StucBufIsland;

typedef struct EdgeCorners {
	FaceCorner corners[2];
} EdgeCorners;

typedef struct FaceBuf {
	PixtyI32Arr faces;
	I32 island;
} FaceBuf;

typedef struct FaceBufArr {
	FaceBuf *pArr;
	I32 size;
	I32 count;
} FaceBufArr;

typedef struct StucIdxRedir {
	U32 idx : 31;
	U32 redir : 1;
} StucIdxRedir;

typedef struct StucIdxRedirArr {
	StucIdxRedir *pArr;
	I32 size;
	I32 count;
} StucIdxRedirArr;

typedef struct StucBorderBb {
	PixtyV2_F32 min;
	PixtyV2_F32 max;
	I32 border;
} StucBorderBb;

typedef struct StucBorderBbArr {
	StucBorderBb *pArr;
	I32 size;
	I32 count;
} StucBorderBbArr;

typedef struct StucSplitMem {
	FaceBufArr faceBuf;
	StucIdxRedirArr redirArr;
	StucIdxTableArr faceTable;
	StucIdxTableArr edgeTable;
	StucBorderNodeArr edges;
	StucBorderBbArr bb;
} StucSplitMem;

typedef struct StucSplitMesh {
	const void *pUserData;
	PixtyRange (*fpFaceRange)(const void *, I32);
	I32 (*fpEdge)(const void *, FaceCorner);
	I32 (*fpVert)(const void *, I32);
	PixtyV2_F32 (*fpPos)(const void *, I32);
	EdgeCorners (*fpEdgeCorners)(const void *, I32);
	FaceCorner (*fpAdjCorner)(const void *pMeshRaw, FaceCorner corner);
	I32 faceCount;
} StucSplitMesh;

typedef struct StucIslands {
	void *pUserData;
	StucErr (*fpIslandAdd)(const PixalcFPtrs *, void *, I32, I32 *);
	StucErr (*fpRangeSet)(void *, I32, PixtyRange);
	StucErr (*fpFacesInit)(const PixalcFPtrs *, void *, I32, I32 **);
	StucErr (*fpBorderInit)(const PixalcFPtrs *, void *, I32, I32 *);
	StucErr (*fpBorderAddEdge)(const PixalcFPtrs *, void *, I32, I32, I32, FaceCorner, I32);
	StucErr (*fpBorderMarkAsOuter)(void *, I32, I32, const ClutreBb *);
} StucIslands;

StucErr stucSplitToIslands(
	const PixalcFPtrs *pAlloc,
	StucSplitMem *pMem,
	const StucSplitMesh *pMesh,
	StucIslands *pIslands,
	bool (*fpSplitPredicate)(const void *, I32)
);

//TODO move this out
static
FaceCorner callGetAdjCorner(const void *pMeshRaw, FaceCorner corner) {
	FaceCorner adj = {0};
	stucGetAdjCorner(pMeshRaw, corner, &adj);
	return adj;
}

void stucSplitMemDestroy(const PixalcFPtrs *pAlloc, StucSplitMem *pMem);