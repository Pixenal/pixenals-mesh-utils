/* 
SPDX-FileCopyrightText: 2025 Caleb Dawson
SPDX-License-Identifier: Apache-2.0
*/

#pragma once

#include <float.h>

#include <pixenals_math_utils.h>
#include <pixenals_structs.h>

#define PIXMSH_NGON_MAX_SIZE 256

#ifndef PIX_FORCE_INLINE
#ifdef NDEBUG
#ifdef WIN32
#define PIX_FORCE_INLINE __forceinline
#else
#define PIX_FORCE_INLINE __attribute__((always_inline)) static inline
#endif
#else
#define PIX_FORCE_INLINE static inline
#endif
#endif


typedef struct PixmshFaceRange {
	I32 start;
	I32 size;
} PixmshFaceRange;

typedef struct PixmshFaceCorner {
	int32_t face;
	int32_t corner;
} PixmshFaceCorner;

typedef struct PixmshBaseTriVerts {
	PixtyV3_F32 xyz[4];
	PixtyV2_F32 uv[4];
	float scale[4];
} PixmshBaseTriVerts;

typedef struct PixmshV2Bb {
	PixtyV2_F32 min;
	PixtyV2_F32 max;
} PixmshV2Bb;

PixmshV2Bb pixmshV2BbGet(
	const void *pMesh,
	PixtyV2_F32 (*fpPos)(const void *, PixmshFaceRange, int32_t),
	PixmshFaceRange face
);
static inline
bool pixmshV2BbOverlap(const PixmshV2Bb *pB, const PixmshV2Bb *pA) {
	return
		pA->min.d[0] <= pB->max.d[0] && pA->max.d[0] >= pB->min.d[0] &&
		pA->min.d[1] <= pB->max.d[1] && pA->max.d[1] >= pB->min.d[1];
}

typedef struct PixmshEar {
	struct PixmshEar *pNext;
	struct PixmshEar *pPrev;
	int32_t cornerPrev;
	int32_t corner;
	int32_t cornerNext;
	float len;
} PixmshEar;

typedef struct PixmshTriangulateState {
	PixalcLinAlloc earAlloc;
	PixmshEar *pEarList;
	const void *pMesh;
	PixtyV3_F32 (* fpPos)(const void *, PixmshFaceRange, int32_t);
	int8_t *pRemoved;
	PixmshFaceRange face;
	PixtyV3_F32 normal;
} PixmshTriangulateState;

PIX_FORCE_INLINE
bool pixmshDoesEarIntersectFace(
	const PixmshTriangulateState *pState,
	const int32_t *pTriIdx,
	const PixtyV3_F32 *pTri,
	const PixtyV3_F32 *pTriNormal
) {
	PixtyV3_F32 normal = pixmV3F32Normalize(*pTriNormal);
	float triDistFromOrigin = _(normal V3DOT pTri[0]) * -1.0f;
	for (int32_t i = 0; i < pState->face.size; ++i) {
		if (i == pTriIdx[0] || i == pTriIdx[1] || i == pTriIdx[2]) {
			continue;
		}
		PixtyV3_F32 point = pState->fpPos(pState->pMesh, pState->face, i);
		float distFromPlane = _(normal V3DOT point) + triDistFromOrigin;
		PixtyV3_F32 projPoint = _(point V3ADD _(normal V3MULS (distFromPlane * -1.0f)));
		PixtyV3_F32 bc = pixmCartesianToBarycentric(pTri, &projPoint, &normal);
		if (_(bc V3GREAT (PixtyV3_F32){0})) {
			return true;
		}
	}
	return false;
}

static inline
I32 pixmshGetCornerPrev(I32 corner, PixmshFaceRange face) {
	PIX_ERR_ASSERT("", corner >= 0 && corner < face.size);
	I32 prev = corner ? corner - 1 : face.size - 1;
	return prev;
}

static inline
I32 pixmshGetCornerNext(I32 corner, PixmshFaceRange face) {
	PIX_ERR_ASSERT("", corner >= 0 && corner < face.size);
	I32 next = (corner + 1) % face.size;
	return next;
}

static inline
int32_t pixmshGetNextRemaining(
	const PixmshTriangulateState *pState,
	int32_t corner,
	PixmshFaceRange face
) {
	PIX_ERR_ASSERT("", corner >= 0 && corner < face.size);
	int32_t start = corner;
	while (corner = pixmshGetCornerNext(corner, face), corner != start) {
		if (!pState->pRemoved[corner]) {
			return corner;
		}
	}
	PIX_ERR_ASSERT("no corners remain", false);
	return -1;
}

static inline
int32_t pixmshGetPrevRemaining(
	const PixmshTriangulateState *pState,
	int32_t corner,
	PixmshFaceRange face
) {
	PIX_ERR_ASSERT("", corner >= 0 && corner < face.size);
	int32_t start = corner;
	while (corner = pixmshGetCornerPrev(corner, face), corner != start) {
		if (!pState->pRemoved[corner]) {
			return corner;
		}
	}
	PIX_ERR_ASSERT("no corners remain", false);
	return -1;
}

PIX_FORCE_INLINE
PixmshEar *pixmshAddEarCandidate(PixmshTriangulateState *pState, int32_t corner) {
	int32_t cornerPrev = pixmshGetPrevRemaining(pState, corner, pState->face);
	int32_t cornerNext = pixmshGetNextRemaining(pState, corner, pState->face);
	PixtyV3_F32 a = pState->fpPos(pState->pMesh, pState->face, cornerPrev);
	PixtyV3_F32 b = pState->fpPos(pState->pMesh, pState->face, corner);
	PixtyV3_F32 c = pState->fpPos(pState->pMesh, pState->face, cornerNext);
	PixtyV3_F32 ac = _(c V3SUB a);
	PixtyV3_F32 cross = _(_(b V3SUB a) V3CROSS ac);
	if (_(cross V3DOT pState->normal) <= 0 || // ear is concave or degenerate
		pixmshDoesEarIntersectFace(
			pState,
			(int32_t[]){cornerPrev, corner, cornerNext},
			(PixtyV3_F32[]){a, b, c},
			&cross
		)
	) {
		return NULL;
	}
	float len = pixmV3F32Len(ac);
	PixmshEar *pNewEar = NULL;
	if (!pState->pEarList) {
		pixalcLinAlloc(&pState->earAlloc, (void **)&pState->pEarList, 1);
		pNewEar = pState->pEarList;
	}
	else {
		PixmshEar *pEar = pState->pEarList;
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

PIX_FORCE_INLINE
void pixmshAddAdjEarCandidates(PixmshTriangulateState *pState, PixmshEar *pEar) {
	int32_t cornerNext = pixmshGetNextRemaining(pState, pEar->corner, pState->face);
	int32_t cornerPrev = pixmshGetPrevRemaining(pState, pEar->corner, pState->face);
	pixmshAddEarCandidate(pState, cornerNext);
	pixmshAddEarCandidate(pState, cornerPrev);
}

static inline
PixmshEar *pixmshAddEar(PixmshTriangulateState *pState, int32_t *pCount, uint8_t *pTris) {
	PixmshEar *pEar = pState->pEarList;
	uint8_t *pTri = pTris + *pCount * 3;
	int32_t cornerPrev = pixmshGetPrevRemaining(pState, pEar->corner, pState->face);
	int32_t cornerNext = pixmshGetNextRemaining(pState, pEar->corner, pState->face);
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
void pixmshRemoveEar(PixmshTriangulateState *pState) {
	PixmshEar *pEar = pState->pEarList;
	if (pEar->pNext) {
		pEar->pNext->pPrev = NULL;
	}
	pState->pEarList = pEar->pNext;
}

static inline
bool pixmshIsMarkedSkip(PixtyI32Arr *pSkip, int32_t idx) {
	for (int32_t i = 0; i < pSkip->count; ++i) {
		if (idx == pSkip->pArr[i]) {
			return true;
		}
	}
	return false;
}

PIX_FORCE_INLINE
int32_t pixmshGetNonDegenBoundCorner(
	const PixmshFaceRange face,
	const void *pMesh,
	PixtyV2_F32 (* fpPos) (const void *, const PixmshFaceRange, int32_t),
	bool useMin,
	PixtyI32Arr *pExternSkip,
	float *pDet
) {
	PIX_ERR_ASSERT("", face.start >= 0 && face.size >= 3);
	int32_t skipArr[PIXMSH_NGON_MAX_SIZE] = {0};
	PixtyI32Arr skip = {.pArr = skipArr};
	do {
		int32_t corner = 0;
		PixtyV2_F32 boundPos = {FLT_MAX, FLT_MAX};
		boundPos = useMin ? boundPos : _(boundPos V2MULS -1.0f);
		for (int32_t i = 0; i < face.size; ++i) {
			if (pixmshIsMarkedSkip(&skip, i) ||
				pExternSkip && pixmshIsMarkedSkip(pExternSkip, i)
			) {
				continue;
			}
			PixtyV2_F32 pos = fpPos(pMesh, face, i);
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
		int32_t prev = corner == 0 ? face.size - 1 : corner - 1;
		int32_t next = (corner + 1) % face.size;
		PixtyV2_F32 a = fpPos(pMesh, face, prev);
		PixtyV2_F32 b = fpPos(pMesh, face, corner);
		PixtyV2_F32 c = fpPos(pMesh, face, next);
		//alt formula for determinate,
		//shorter and less likely to cause numerical error
		float det =
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
	} while(skip.count < face.size);
	return -1;
}

//finds corner on convex hull of face, & determines wind direction from that
//returns 0 for clockwise, 1 for counterclockwise, & 2 if degenerate
PIX_FORCE_INLINE
int32_t pixmshCalcFaceWind(
	PixmshFaceRange face,
	const void *pMesh,
	PixtyV2_F32 (* fpPos) (const void *, const PixmshFaceRange, int32_t)
) {
	float det = .0f;
	int32_t corner = pixmshGetNonDegenBoundCorner(face, pMesh, fpPos, true, NULL, &det);
	return corner != -1 ? det > .0f : 2;
}

static
PixtyV3_F32 pixmshGetTriNormal(
	const void *pMesh,
	PixmshFaceRange face,
	int32_t corner,
	PixtyV3_F32 (* fpPos) (const void *, PixmshFaceRange, int32_t)
) {
	int32_t prev = corner == 0 ? face.size - 1 : corner - 1;
	int32_t next = (corner + 1) % face.size;
	PixtyV3_F32 a = fpPos(pMesh, face, prev);
	PixtyV3_F32 b = fpPos(pMesh, face, corner);
	PixtyV3_F32 c = fpPos(pMesh, face, next);
	return _(_(b V3SUB a) V3CROSS _(c V3SUB a));
}

typedef struct AxisBounds{
	int32_t minIdx;
	int32_t maxIdx;
	float min;
	float max;
	float len;
} AxisBounds;

static inline
void pixmshAxisBoundsCalcLen(AxisBounds *pBounds) {
	pBounds->len = pBounds->max - pBounds->min;
	PIX_ERR_ASSERT("", pBounds->len >= .0f);
}

static inline
void pixmshAxisBoundsCmp(AxisBounds *pBounds, float pos, int32_t idx) {
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
int32_t pixmshAxisBoundsMake(
	PixmshFaceRange face,
	const void *pMesh,
	PixtyV3_F32 (* fpPos) (const void *, PixmshFaceRange, int32_t),
	PixtyI32Arr *pSkip,
	AxisBounds *pBounds
) {
	for (int32_t i = 0; i < 3; ++i) {
		pBounds[i].max = -FLT_MAX;
		pBounds[i].min = FLT_MAX;
	}
	for (int32_t i = 0; i < face.size; ++i) {
		if (pSkip && pixmshIsMarkedSkip(pSkip, i)) {
			continue;
		}
		PixtyV3_F32 pos = fpPos(pMesh, face, i);
		pixmshAxisBoundsCmp(pBounds + 0, pos.d[0], i);
		pixmshAxisBoundsCmp(pBounds + 1, pos.d[1], i);
		pixmshAxisBoundsCmp(pBounds + 2, pos.d[2], i);
	}
	pixmshAxisBoundsCalcLen(pBounds + 0);
	pixmshAxisBoundsCalcLen(pBounds + 1);
	pixmshAxisBoundsCalcLen(pBounds + 2);
	int32_t lowAxis = pBounds[0].len > pBounds[1].len;
	if (pBounds[2].len < pBounds[1].len && pBounds[2].len < pBounds[0].len) {
		lowAxis = 2;
	}
	return lowAxis;
}

static inline
void pixmshMarkSkip(PixtyI32Arr *pSkip, int32_t idx) {
	PIX_ERR_ASSERT("", pSkip->count < PIXMSH_NGON_MAX_SIZE);
	pSkip->pArr[pSkip->count] = idx;
	++pSkip->count;
}

typedef struct PixmshTriIntf {
	const void *pMesh;
	PixtyV3_F32 (* fpPos) (const void *, PixmshFaceRange, int32_t);
} PixmshTriIntf;

static inline
PixtyV2_F32 pixmshVertPosXy(const PixmshTriIntf *pInterf, PixmshFaceRange face, I32 corner) {
	PIX_ERR_ASSERT("", corner >= 0 && corner < face.size);
	PixtyV3_F32 pos = pInterf->fpPos(pInterf->pMesh, face, corner);
	return (PixtyV2_F32){pos.d[0], pos.d[1]};
}
static inline
PixtyV2_F32 pixmshVertPosXz(const PixmshTriIntf *pInterf, PixmshFaceRange face, I32 corner) {
	PIX_ERR_ASSERT("", corner >= 0 && corner < face.size);
	PixtyV3_F32 pos = pInterf->fpPos(pInterf->pMesh, face, corner);
	return (PixtyV2_F32){pos.d[0], pos.d[2]};
}
static inline
PixtyV2_F32 pixmshVertPosYz(const PixmshTriIntf *pInterf, PixmshFaceRange face, I32 corner) {
	PIX_ERR_ASSERT("", corner >= 0 && corner < face.size);
	PixtyV3_F32 pos = pInterf->fpPos(pInterf->pMesh, face, corner);
	return (PixtyV2_F32){pos.d[1], pos.d[2]};
}

static inline
PixtyV3_F32 pixmshCalcFaceNormal(
	PixmshFaceRange face,
	const void *pMesh,
	PixtyV3_F32 (* fpPos) (const void *, PixmshFaceRange, int32_t)
) {
	PIX_ERR_ASSERT(
		"invalid face size",
		face.start >= 0 && face.size >= 3 && face.size <= PIXMSH_NGON_MAX_SIZE
	);
	int32_t skipArr[PIXMSH_NGON_MAX_SIZE] = {0};
	PixtyI32Arr skip = {.pArr = skipArr};
	AxisBounds bounds[3] = {0};
	int32_t axis = pixmshAxisBoundsMake(face, pMesh, fpPos, NULL, bounds);
	PixtyV2_F32 (*fpAxisPos)(const void *, PixmshFaceRange, int32_t) =
		axis == 2 ? pixmshVertPosXy : axis ? pixmshVertPosXz : pixmshVertPosYz;
	PixmshTriIntf wrap = {.pMesh = pMesh, .fpPos = fpPos};
	do {
		int32_t minIdx =
			pixmshGetNonDegenBoundCorner(face, &wrap, fpAxisPos, true, &skip, NULL);
		int32_t maxIdx =
			pixmshGetNonDegenBoundCorner(face, &wrap, fpAxisPos, false, &skip, NULL);
		if (minIdx == -1 || maxIdx == -1) {
			return (PixtyV3_F32){0};
		}
		PixtyV3_F32 minNormal = pixmshGetTriNormal(pMesh, face, minIdx, fpPos);
		PixtyV3_F32 maxNormal = pixmshGetTriNormal(pMesh, face, maxIdx, fpPos);
		if (_(minNormal V3DOT maxNormal) <= .0f) {
			pixmshMarkSkip(&skip, minIdx);
			pixmshMarkSkip(&skip, maxIdx);
			continue;
		}
		if (_(minNormal V3EQL maxNormal)) {
			return minNormal;
		}
		return _(
			_(pixmV3F32Normalize(maxNormal) V3ADD pixmV3F32Normalize(minNormal)) V3DIVS
			2.0f
		);
	} while(skip.count < face.size);
	return (PixtyV3_F32){0};
}

//returns tri count (may be less than size - 2 if face is degen)
PIX_FORCE_INLINE
int32_t pixmshTriangulateFace(
	const PixalcFPtrs *pAlloc,
	const PixmshFaceRange face,
	const void *pMesh,
	PixtyV3_F32 (* fpPos)(const void *, PixmshFaceRange, int32_t),
	uint8_t *pTris
) {
	PIX_ERR_ASSERT("", pTris);
	PixmshTriangulateState state = {
		.pMesh = pMesh,
		.fpPos = fpPos,
		.face = face,
		.pRemoved = pAlloc->fpCalloc(face.size, 1),
		.normal = pixmshCalcFaceNormal(face, pMesh, fpPos)
	};
	if (_(state.normal V3EQL (PixtyV3_F32){0})) {
		return 0;
	}

	pixalcLinAllocInit(pAlloc, &state.earAlloc, sizeof(PixmshEar), face.size, true);

	//add initial ears
	for (int32_t i = 0; i < face.size; ++i) {
		pixmshAddEarCandidate(&state, i);
	}
	int32_t triCount = 0;
	while (state.pEarList) {
		PixmshEar *pAddedEar = NULL;
		if (!state.pRemoved[state.pEarList->cornerPrev] &&
			!state.pRemoved[state.pEarList->corner] &&
			!state.pRemoved[state.pEarList->cornerNext]
		) {
			pAddedEar = pixmshAddEar(&state, &triCount, pTris);
		}
		pixmshRemoveEar(&state);
		if (pAddedEar) {
			pixmshAddAdjEarCandidates(&state, pAddedEar);
		}
	}
	pixalcLinAllocDestroy(&state.earAlloc);
	pAlloc->fpFree(state.pRemoved);
	PIX_ERR_ASSERT("", triCount <= face.size - 2);
	return triCount;
}

PIX_FORCE_INLINE
PixtyV3_F32 pixmshGetBarycentricInTri(
	const void *pMesh,
	PixmshFaceRange face,
	PixtyV3_F32 (* fpPos)(const void *, PixmshFaceRange, int32_t),
	const uint8_t *pTriCorners,
	PixtyV2_F32 vert
) {
	PixtyV3_F32 tri[3] = {0};
	for (int32_t i = 0; i < 3; ++i) {
		tri[i] = fpPos(pMesh, face, pTriCorners[i]);
	}
	return pixmCartesianToBarycentric(
		tri,
		&(PixtyV3_F32){.d = {vert.d[0], vert.d[1]}},
		&(PixtyV3_F32){.d = {.0f, .0f, 1.0f}}
	);
}

//Caller must check for nan in return value
PIX_FORCE_INLINE
PixtyV3_F32 pixmshGetBarycentricInFace(
	const void *pMesh,
	PixmshFaceRange face,
	PixtyV2_I16 tile,
	PixtyV3_F32 (* fpPos)(const void *, PixmshFaceRange, int32_t),
	int8_t *pTriCorners,
	PixtyV2_F32 vertV2
) {
	PIX_ERR_ASSERT("", pixmV2F32IsFinite(vertV2));
	PIX_ERR_ASSERT("", (face.size == 3 || face.size == 4) && pTriCorners);
	PixtyV3_F32 vert = {.d = {vertV2.d[0], vertV2.d[1]}};
	PixtyV3_F32 fTile = {.d = {(F32)tile.d[0], (F32)tile.d[1]}};
	PixtyV3_F32 triA[3] = {0};
	for (int32_t i = 0; i < 3; ++i) {
		triA[i] = _(fpPos(pMesh, face, i) V3SUB fTile);
	}
	PixtyV3_F32 up = {.d = {.0f, .0f, 1.0f}};
	PixtyV3_F32 vertBc = pixmCartesianToBarycentric(triA, &vert, &up);
	if (face.size == 4 && pixmV3F32IsFinite(vertBc) && vertBc.d[1] < 0) {
		//base face is a quad, and vert is outside first tri,
		//so use the second tri
		
		PixtyV3_F32 triB[3] = {
			triA[2],
			_(fpPos(pMesh, face, 3) V3SUB fTile),
			triA[0]
		};
		vertBc = pixmCartesianToBarycentric(triB, &vert, &up);
		pTriCorners[0] = 2;
		pTriCorners[1] = 3;
	}
	else {
		for (int32_t k = 0; k < 3; ++k) {
			pTriCorners[k] = (int8_t)k;
		}
	}
	return vertBc;
}

PixtyM3x3 pixmshBuildFaceTbn(
	PixmshFaceRange face,
	const void *pMesh,
	PixtyV3_F32 (*fpPos)(const void *, PixmshFaceRange, int32_t),
	PixtyV2_F32 (*fpUv)(const void *, PixmshFaceRange, int32_t)
);
void pixmshGetTriScale(int32_t size, PixmshBaseTriVerts *pTri);

typedef struct PixmshSplitIdxTable {
	uint32_t idx : 31;
	uint32_t valid : 1;
} PixmshSplitIdxTable;

typedef struct PixmshSplitIdxTableArr {
	PixmshSplitIdxTable *pArr;
	int32_t size;
	int32_t count;
} PixmshSplitIdxTableArr;

struct PixmshBorderNode;

typedef struct PixmshBorderNode {
	PixmshFaceCorner corners[2];
	int32_t idx;
	PixmshSplitIdxTable seen[2];
	bool intern;
} PixmshBorderNode;

typedef struct PixmshBorderLink {
	PixmshBorderNode *pNode;
} PixmshBorderLink;

typedef struct PixmshBorderNodeArr {
	PixmshBorderNode *pArr;
	int32_t size;
	int32_t count;
} PixmshBorderNodeArr;

typedef struct PixmshBufIsland {
	struct PixmshBufIsland *pNext;
	//PixtyI32Arr faces;
	int32_t idx;
} PixmshBufIsland;

typedef struct PixmshEdgeCorners {
	PixmshFaceCorner corners[2];
} PixmshEdgeCorners;

typedef struct PixmshFaceBuf {
	PixtyI32Arr faces;
	int32_t island;
} PixmshFaceBuf;

typedef struct PixmshFaceBufArr {
	PixmshFaceBuf *pArr;
	int32_t size;
	int32_t count;
} PixmshFaceBufArr;

typedef struct PixmshIdxRedir {
	uint32_t idx : 31;
	uint32_t redir : 1;
} PixmshIdxRedir;

typedef struct PixmshIdxRedirArr {
	PixmshIdxRedir *pArr;
	int32_t size;
	int32_t count;
} PixmshIdxRedirArr;

typedef struct PixmshBorderBb {
	PixtyV2_F32 min;
	PixtyV2_F32 max;
	int32_t border;
} PixmshBorderBb;

typedef struct PixmshBorderBbArr {
	PixmshBorderBb *pArr;
	int32_t size;
	int32_t count;
} PixmshBorderBbArr;

typedef struct PixmshSplitMem {
	PixmshFaceBufArr faceBuf;
	PixmshIdxRedirArr redirArr;
	PixmshSplitIdxTableArr faceTable;
	PixmshSplitIdxTableArr edgeTable;
	PixmshBorderNodeArr edges;
	PixmshBorderBbArr bb;
} PixmshSplitMem;

//TODO replace all func ptrs in param lists with typedefs? maybe?
typedef struct PixmshSplitIntfIn {
	const void *pUserData;
	PixmshFaceRange (*fpFaceRange)(const void *, int32_t);
	int32_t (*fpEdge)(const void *, PixmshFaceCorner);
	PixtyV2_F32 (*fpPos)(const void *, int32_t);
	PixmshEdgeCorners (*fpEdgeCorners)(const void *, int32_t);
	PixmshFaceCorner (*fpAdjCorner)(const void *pMeshRaw, PixmshFaceCorner corner);
	int32_t faceCount;
} PixmshSplitIntfIn;

typedef struct PixmshSplitIntfOut {
	void *pUserData;
	PixErr (*fpIslandAdd)(const PixalcFPtrs *, void *, int32_t, int32_t *);
	PixErr (*fpRangeSet)(void *, int32_t, PixtyRange);
	PixErr (*fpFacesInit)(const PixalcFPtrs *, void *, int32_t, int32_t **);
	PixErr (*fpBorderInit)(const PixalcFPtrs *, void *, int32_t, int32_t *);
	PixErr (*fpBorderAddEdge)(const PixalcFPtrs *, void *, int32_t, int32_t, int32_t, PixmshFaceCorner, int32_t);
	PixErr (*fpBorderMarkAsOuter)(void *, int32_t, int32_t, const PixmshV2Bb *);
} PixmshSplitIntfOut;

PixErr pixmshSplitToIslands(
	const PixalcFPtrs *pAlloc,
	PixmshSplitMem *pMem,
	const PixmshSplitIntfIn *pMesh,
	PixmshSplitIntfOut *pIslands,
	bool (*fpSplitPredicate)(const void *, int32_t)
);

void pixmshSplitMemDestroy(const PixalcFPtrs *pAlloc, PixmshSplitMem *pMem);
