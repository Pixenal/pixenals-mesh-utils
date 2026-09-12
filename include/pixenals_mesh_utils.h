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

static inline
bool pixmshCornerCmp(PixmshFaceCorner a, PixmshFaceCorner b) {
	return a.face == b.face && a.corner == b.corner;
}

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

typedef struct PixmshTriMem {
	PixtyI8Arr removedArr;
	PixalcLinAlloc earAlloc;
} PixmshTriMem;

typedef struct PixmshTriangulateState {
	PixmshTriMem *pMem;
	PixmshEar *pEarList;
	const void *pMesh;
	PixtyV3_F32 (* fpPos)(const void *, PixmshFaceRange, int32_t);
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
		if (!pState->pMem->removedArr.pArr[corner]) {
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
		if (!pState->pMem->removedArr.pArr[corner]) {
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
		pixalcLinAlloc(&pState->pMem->earAlloc, (void **)&pState->pEarList, 1);
		pNewEar = pState->pEarList;
	}
	else {
		PixmshEar *pEar = pState->pEarList;
		while(pEar->pNext && len > pEar->pNext->len) {
			pEar = pEar->pNext;
		}
		pixalcLinAlloc(&pState->pMem->earAlloc, (void **)&pNewEar, 1);
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
	
	pState->pMem->removedArr.pArr[pEar->corner] = true;
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

static inline
PixtyV3_F32 pixmshCalcFaceNormal(
	PixmshFaceRange face,
	const void *pMesh,
	PixtyV3_F32 (* fpPos) (const void *, PixmshFaceRange, int32_t),
	bool normalize
) {
	PixtyV3_F32 normal = {0};
	if (face.size == 3) {
		 normal = pixmshGetTriNormal(pMesh, face, 0, fpPos);
	}
	else if (face.size == 4) {
		PixtyV3_F32 a = pixmshGetTriNormal(pMesh, face, 0, fpPos);
		PixtyV3_F32 b = pixmshGetTriNormal(pMesh, face, 2, fpPos);
		normal = _(a V3ADD b);
	}
	else {
		PIX_ERR_ASSERT("invalid face", face.start >= 0 && face.size > 4);
		//"Newell's method for computing the plane equation of a polygon" Filippo Tampieri 1992:
		//https://dl.acm.org/doi/10.5555/130745.130783
		//see also: https://wikis.khronos.org/opengl/Calculating_a_Surface_Normal
		for (I32 i = 0; i < face.size; ++i) {
			PixtyV3_F32 pos = fpPos(pMesh, face, i);
			PixtyV3_F32 posNext = fpPos(pMesh, face, (i + 1) % face.size);
			normal.d[0] += (pos.d[1] - posNext.d[1]) * (pos.d[2] + posNext.d[2]);
			normal.d[1] += (pos.d[2] - posNext.d[2]) * (pos.d[0] + posNext.d[0]);
			normal.d[2] += (pos.d[0] - posNext.d[0]) * (pos.d[1] + posNext.d[1]);
		}
	}
	return normalize ? pixmV3F32Normalize(normal) : normal;
}

static inline
void pixmshTriMemClear(PixmshTriMem *pMem) {
	if (pMem->removedArr.count) {
		memset(
			pMem->removedArr.pArr,
			0,
			pMem->removedArr.count * PIXALC_ITEMSIZE(pMem->removedArr.pArr)
		);
		pMem->removedArr.count = 0;
	}
	if (pMem->earAlloc.valid) {
		pixalcLinAllocClear(&pMem->earAlloc);
	}
}

//returns tri count (may be less than size - 2 if face is degen)
//TODO take in PixtyU8Arr and pass final tri-count as .count,
//this way func can return error
PIX_FORCE_INLINE
int32_t pixmshTriangulateFace(
	const PixalcFPtrs *pAlloc,
	PixmshTriMem *pMem,
	const PixmshFaceRange face,
	const void *pMesh,
	PixtyV3_F32 (* fpPos)(const void *, PixmshFaceRange, int32_t),
	uint8_t *pTris
) {
	PIX_ERR_ASSERT("", pTris);
	pixmshTriMemClear(pMem);
	PixmshTriangulateState state = {
		.pMem = pMem,
		.pMesh = pMesh,
		.fpPos = fpPos,
		.face = face,
		.normal = pixmshCalcFaceNormal(face, pMesh, fpPos, false)
	};
	if (_(state.normal V3EQL (PixtyV3_F32){0})) {
		return 0;
	}
	pMem->removedArr.count = face.size;
	PIXALC_DYN_ARR_RESIZE_ZERO(pAlloc, &pMem->removedArr, pMem->removedArr.count);
	if (!pMem->earAlloc.valid) {
		pixalcLinAllocInit(pAlloc, &pMem->earAlloc, sizeof(PixmshEar), face.size, true);
	}

	//add initial ears
	for (int32_t i = 0; i < face.size; ++i) {
		pixmshAddEarCandidate(&state, i);
	}
	int32_t triCount = 0;
	while (state.pEarList) {
		PixmshEar *pAddedEar = NULL;
		if (!pMem->removedArr.pArr[state.pEarList->cornerPrev] &&
			!pMem->removedArr.pArr[state.pEarList->corner] &&
			!pMem->removedArr.pArr[state.pEarList->cornerNext]
		) {
			pAddedEar = pixmshAddEar(&state, &triCount, pTris);
		}
		pixmshRemoveEar(&state);
		if (pAddedEar) {
			pixmshAddAdjEarCandidates(&state, pAddedEar);
		}
	}
	PIX_ERR_ASSERT("", triCount <= face.size - 2);
	return triCount;
}

static inline
void pixmshTriMemDestroy(const PixalcFPtrs *pAlloc, PixmshTriMem *pMem) {
	PIXALC_DYN_ARR_DESTROY(pAlloc, &pMem->removedArr);
	if (pMem->earAlloc.valid) {
		pixalcLinAllocDestroy(&pMem->earAlloc);
	}
	*pMem = (PixmshTriMem){0};
}

PIX_FORCE_INLINE
PixtyV3_F32 pixmshGetBarycentricInTri(
	const void *pMesh,
	PixmshFaceRange face,
	PixtyV3_F32 (* fpPos)(const void *, PixmshFaceRange, int32_t),
	const int8_t *pTriCorners,
	PixtyV2_F32 vert
) {
	PixtyV3_F32 tri[3] = {0};
	for (int32_t i = 0; i < 3; ++i) {
		tri[i] = fpPos(pMesh, face, (int32_t)pTriCorners[i]);
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
	if (face.size == 4 && (!pixmV3F32IsFinite(vertBc) || vertBc.d[1] < 0)) {
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

typedef struct PixmshSplitEdgeInfo {
	PixmshBorderNode *pNode;
	PixmshFaceCorner corner;
	int32_t edge;
} PixmshSplitEdgeInfo;

typedef struct PixmshSplitEdgeInfoArr {
	PixmshSplitEdgeInfo *pArr;
	int32_t size;
	int32_t count;
} PixmshSplitEdgeInfoArr;

typedef struct PixmshBorderNodePtrArr {
	PixmshBorderNode **pArr;
	int32_t size;
} PixmshBorderNodePtrArr;

typedef struct PixmshSplitMem {
	PixmshFaceBufArr faceBuf;
	PixmshIdxRedirArr redirArr;
	PixmshSplitIdxTableArr faceTable;
	PixmshSplitIdxTableArr edgeTable;
	PixmshBorderNodeArr edges;
	PixmshBorderBbArr bb;
	PixmshSplitEdgeInfoArr edgeBuf;
	PixmshBorderNodePtrArr fallbacks;
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

void pixmshSplitMemInit(const PixalcFPtrs *pAlloc, PixmshSplitMem *pMem, int32_t faceCount);
PixErr pixmshSplitConstructIslandsFromAdj(
	const PixalcFPtrs *pAlloc,
	PixmshSplitMem *pMem,
	const PixmshSplitIntfIn *pMesh,
	PixmshSplitIntfOut *pIslands,
	int32_t splitTotal,
	int32_t *pIslandCount
);
PixErr pixmshSplitBordersMake(
	const PixalcFPtrs *pAlloc,
	PixmshSplitMem *pMem,
	const PixmshSplitIntfIn *pMesh,
	PixmshSplitIntfOut *pIslands,
	int32_t islandCount
);
PixmshBorderNode *pixmshBorderNodeInit(
	const PixalcFPtrs *pAlloc,
	PixmshBorderNodeArr *pEdges,
	const PixmshEdgeCorners *pCorners
);
void pixmshSplitIslandIdxInit(
	const PixalcFPtrs *pAlloc,
	PixmshSplitIdxTableArr *pFaceTable,
	PixmshIdxRedirArr *pArr,
	int32_t face
);
void pixmshSplitEdgeTableAdd(
	const PixalcFPtrs *pAlloc,
	PixmshSplitIdxTableArr *pTable,
	int32_t edge,
	int32_t idx
);
void pixmshSplitHandleAdj(
	const PixalcFPtrs *pAlloc,
	PixmshSplitMem *pMem,
	int32_t *pFaces
);

PIX_FORCE_INLINE
PixErr pixmshSplitFindAdjForCorner(
	const PixalcFPtrs *pAlloc,
	PixmshSplitMem *pMem,
	const PixmshSplitIntfIn *pMesh,
	bool (*fpSplitPredicate)(const void *, int32_t),
	PixmshFaceCorner corner,
	int32_t *pSplitTotal
) {
	PixErr err = PIX_ERR_SUCCESS;
	int32_t edge = pMesh->fpEdge(pMesh->pUserData, corner);
	if (edge < pMem->edgeTable.count && pMem->edgeTable.pArr[edge].valid) {
		return err;
	}
	PixmshEdgeCorners corners = pMesh->fpEdgeCorners(pMesh->pUserData, edge);
	int32_t faces[2] = {corners.corners[0].face, corners.corners[1].face};
	bool borderEdge = faces[0] == -1 || faces[1] == -1;
	if (borderEdge || fpSplitPredicate && fpSplitPredicate(pMesh->pUserData, edge)) {
		++*pSplitTotal;
		if (faces[0] != -1 && !pMem->faceTable.pArr[faces[0]].valid) {
			pixmshSplitIslandIdxInit(pAlloc, &pMem->faceTable, &pMem->redirArr, faces[0]);
		}
		if (faces[1] != -1 && !pMem->faceTable.pArr[faces[1]].valid) {
			pixmshSplitIslandIdxInit(pAlloc, &pMem->faceTable, &pMem->redirArr, faces[1]);
		}
		if (edge >= pMem->edgeTable.size || !pMem->edgeTable.pArr[edge].valid) {
			PixmshBorderNode *pNode = pixmshBorderNodeInit(pAlloc, &pMem->edges, &corners);
			pixmshSplitEdgeTableAdd(pAlloc, &pMem->edgeTable, edge, pNode->idx);
		}
		return err;
	}
	pixmshSplitHandleAdj(pAlloc, pMem, faces);
	return err;
}

PIX_FORCE_INLINE
PixErr pixmshSplitToIslands(
	const PixalcFPtrs *pAlloc,
	PixmshSplitMem *pMem,
	const PixmshSplitIntfIn *pMesh,
	PixmshSplitIntfOut *pIslands,
	bool makeBorders,
	bool (*fpSplitPredicate)(const void *, int32_t)
) {
	PixErr err = PIX_ERR_SUCCESS;
	pixmshSplitMemInit(pAlloc, pMem, pMesh->faceCount);
	int32_t splitTotal = 0;
	//determine adjacency
	for (int32_t i = 0; i < pMesh->faceCount; ++i) {
		PixmshFaceRange face = pMesh->fpFaceRange(pMesh->pUserData, i);
		for (int32_t j = 0; j < face.size; ++j) {
			PixmshFaceCorner corner = {.face = i, .corner = j};
			err = pixmshSplitFindAdjForCorner(
				pAlloc,
				pMem,
				pMesh,
				fpSplitPredicate,
				corner,
				&splitTotal
			);
			PIX_ERR_THROW_IFNOT(err, "", 0);
		}
	}
	PIX_ERR_THROW_IFNOT_COND(err, pMem->redirArr.count, "failed to split mesh", 0);

	int32_t islandCount = 0;
	err = pixmshSplitConstructIslandsFromAdj(
		pAlloc,
		pMem,
		pMesh,
		pIslands,
		splitTotal,
		&islandCount
	);
	PIX_ERR_THROW_IFNOT(err, "", 0);
	
	if (makeBorders) {
		err = pixmshSplitBordersMake(pAlloc, pMem, pMesh, pIslands, islandCount);
		PIX_ERR_THROW_IFNOT(err, "failed to make border(s)", 0);
	}
	PIX_ERR_CATCH(0, err, ;);
	return err;
}

void pixmshSplitMemDestroy(const PixalcFPtrs *pAlloc, PixmshSplitMem *pMem);
