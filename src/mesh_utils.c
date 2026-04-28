/* 
SPDX-FileCopyrightText: 2025 Caleb Dawson
SPDX-License-Identifier: Apache-2.0
*/

#include <float.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include <pixenals_mesh_utils.h>
#include <pixenals_math_utils.h>
#include <pixenals_error_utils.h>

typedef float F32;
typedef double F64;
typedef uint8_t U8;
typedef int8_t I8;
typedef int32_t I32;

PixmshV2Bb pixmshV2BbGet(
	const void *pMesh,
	PixtyV2_F32 (*fpPos)(const void *, PixmshFaceRange, int32_t),
	PixmshFaceRange face
) {
	PixmshV2Bb bbox = {
		.min = {.d = {FLT_MAX, FLT_MAX}},
		.max = {.d = {-FLT_MAX, -FLT_MAX}}
	};
	for (I32 i = 0; i < face.size; ++i) {
		PixtyV2_F32 pos = fpPos(pMesh, face, i);
		if (pos.d[0] < bbox.min.d[0]) {
			bbox.min.d[0] = pos.d[0];
		}
		if (pos.d[1] < bbox.min.d[1]) {
			bbox.min.d[1] = pos.d[1];
		}
		if (pos.d[0] > bbox.max.d[0]) {
			bbox.max.d[0] = pos.d[0];
		}
		if (pos.d[1] > bbox.max.d[1]) {
			bbox.max.d[1] = pos.d[1];
		}
	}
	//>= not >,
	//because faces can be flat (they may be facing sideways in a map for instance)
	PIX_ERR_ASSERT("", _(bbox.max V2GREATEQL bbox.min));
	return bbox;
}

PixtyM3x3 pixmshBuildFaceTbn(
	PixmshFaceRange face,
	const void *pMesh,
	PixtyV3_F32 (*fpPos)(const void *, PixmshFaceRange, I32),
	PixtyV2_F32 (*fpUv)(const void *, PixmshFaceRange, I32)
) {
	I32 corner = face.start;
	PixtyV2_F32 uv = fpUv(pMesh, face, corner);
	PixtyV3_F32 vert = fpPos(pMesh, face, corner);
	I32 next = face.start + 1;
	PixtyV2_F32 uvNext = fpUv(pMesh, face, next);
	PixtyV3_F32 vertNext = fpPos(pMesh, face, next);
	I32 prev = face.start + face.size - 1;
	PixtyV2_F32 uvPrev = fpUv(pMesh, face, prev);
	PixtyV3_F32 vertPrev = fpPos(pMesh, face, prev);
	//uv space direction vectors,
	//forming the coefficient matrix
	PixtyM2x2 coeffMat = {0};
	*(PixtyV2_F32 *)&coeffMat.d[0] = _(uvNext V2SUB uv);
	*(PixtyV2_F32 *)&coeffMat.d[1] = _(uvPrev V2SUB uv);
	//object space direction vectors,
	//forming the variable matrix
	PixtyM2x3 varMat = {0};
	PixtyV3_F32 osDirA = _(vertNext V3SUB vert);
	PixtyV3_F32 osDirB = _(vertPrev V3SUB vert);
	*(PixtyV3_F32 *)&varMat.d[0] = osDirA;
	*(PixtyV3_F32 *)&varMat.d[1] = osDirB;
	PixtyM2x2 coeffMatInv = pixmM2x2Invert(coeffMat);
	PixtyM2x3 tb = pixmM2x2MultiplyM2x3(coeffMatInv, varMat);
	PixtyM3x3 tbn = {0};
	*(PixtyV3_F32 *)&tbn.d[0] = pixmV3F32Normalize(*(PixtyV3_F32 *)&tb.d[0]);
	*(PixtyV3_F32 *)&tbn.d[1] = pixmV3F32Normalize(*(PixtyV3_F32 *)&tb.d[1]);
	PixtyV3_F32 normal = _(osDirA V3CROSS osDirB);
	*(PixtyV3_F32 *)&tbn.d[2] = pixmV3F32Normalize(normal);
	return tbn;
}

void pixmshGetTriScale(I32 size, PixmshBaseTriVerts *pTri) {
	for (I32 i = 0; i < size; ++i) {
		I32 iLast = i == 0 ? size - 1 : i - 1;
		I32 iNext = (i + 1) % size;
		F32 uvArea = pixmV2F32TriArea(pTri->uv[iLast], pTri->uv[i], pTri->uv[iNext]);
		F32 xyzArea = pixmV3F32TriArea(pTri->xyz[iLast], pTri->xyz[i], pTri->xyz[iNext]);
		pTri->scale[i] = xyzArea / uvArea;
	}
}

static
void islandIdxInit(
	const PixalcFPtrs *pAlloc,
	PixmshSplitIdxTableArr *pFaceTable,
	PixmshIdxRedirArr *pArr,
	I32 face
) {
	I32 newIdx = 0;
	PIXALC_DYN_ARR_ADD(PixmshIdxRedir, pAlloc, pArr, newIdx);
	pArr->pArr[newIdx] = (PixmshIdxRedir){.idx = newIdx};
	pFaceTable->pArr[face].idx = newIdx;
	pFaceTable->pArr[face].valid = true;
}

static
PixmshBorderNode *pixmshBorderNodeGet(
	const PixmshSplitIntfIn *pMesh,
	const PixmshSplitIdxTable *pTable,
	PixalcLinAlloc *pEdgeAlloc,
	I32 vert
) {
	return pTable[vert].valid ? pixalcLinAllocIdx(pEdgeAlloc, pTable[vert].idx) : NULL;
}

static
PixmshBorderNode *pixmshBorderNodeInit(
	const PixalcFPtrs *pAlloc,
	PixmshBorderNodeArr *pEdges,
	const PixmshEdgeCorners *pCorners
) {
	I32 idx = 0;
	PIXALC_DYN_ARR_ADD(PixmshBorderNode, pAlloc, pEdges, idx);
	pEdges->pArr[idx] = (PixmshBorderNode) {
		.idx = idx,
		.corners[0] = pCorners->corners[0],
		.corners[1] = pCorners->corners[1]
	};
	return pEdges->pArr + idx;
}

#define STUC_IDX_REDIR_THRES 8

static
PixmshIdxRedir *getBuf(const PixmshSplitIdxTableArr *pFaceTable, PixmshIdxRedirArr *pArr, I32 face) {
	PixmshIdxRedir *pId = NULL;
	I32 idx = pFaceTable->pArr[face].idx;
	I32 i = 0;
	do {
		PIX_ERR_ASSERT("", idx < pArr->count && i < pArr->count);
		pId = pArr->pArr + idx;
		idx = pId->idx;
	} while(++i, pId->redir);
	if (i > STUC_IDX_REDIR_THRES) {
		I32 finalIdx = idx;
		idx = pFaceTable->pArr[face].idx;
		do {
			pId = pArr->pArr + idx;
			idx = pId->idx;
			pId->idx = finalIdx;
		} while(pId->redir);
	}
	return pId;
}

static
I32 getEdgeIslandSingle(
	PixmshSplitMem *pMem,
	const PixmshBorderNode *pEdge,
	I32 side
) {
	PixmshFaceCorner corner = pEdge->corners[side];
	if (corner.face == -1) {
		return -1;
	}
	else {
		I32 face = getBuf(&pMem->faceTable, &pMem->redirArr, corner.face)->idx;
		return pMem->faceBuf.pArr[face].island;
	}
}

static
void getEdgeIslands(
	PixmshSplitMem *pMem,
	const PixmshBorderNode *pEdge,
	I32 *pIslands
) {
	pIslands[0] = getEdgeIslandSingle(pMem, pEdge, 0);
	pIslands[1] = getEdgeIslandSingle(pMem, pEdge, 1);
	PIX_ERR_ASSERT("floating edge", pIslands[0] != -1 || pIslands[1] != -1);
}

static
bool isEdgeIntern(
	PixmshSplitMem *pMem,
	PixmshBorderNode *pEdge,
	I32 *pIslands
) {
	if (pEdge->intern) {
		return true;
	}
	I32 islands[2] = {0};
	getEdgeIslands(pMem, pEdge, islands);
	pEdge->intern = islands[0] == islands[1];
	if (pIslands) {
		pIslands[0] = islands[0];
		pIslands[1] = islands[1];
	}
	return pEdge->intern;
}

static
bool seenThisEdge(const PixmshBorderNode *pEdge, I32 island) {
	return
		pEdge->seen[0].valid && pEdge->seen[0].idx == island ||
		pEdge->seen[1].valid && pEdge->seen[1].idx == island;
}

static
void markEdgeSeen(PixmshBorderNode *pEdge, PixmshFaceCorner corner, I32 island) {
	PIX_ERR_ASSERT("", !pEdge->seen[0].valid || !pEdge->seen[1].valid);
	bool side = corner.face == pEdge->corners[1].face;
	PIX_ERR_ASSERT(
		"",
		pEdge->corners[side].face == corner.face &&
		pEdge->corners[0].face != pEdge->corners[1].face
	);
	pEdge->seen[side] = (PixmshSplitIdxTable){.idx = island, .valid = true};
}

static
void pixmshBorderBbCmp(
	const PixmshSplitIntfIn *pMesh,
	PixmshBorderBb *pBb,
	I32 corner,
	I32 border
) {
	PixtyV2_F32 pos = pMesh->fpPos(pMesh->pUserData, corner);
	pBb->min.d[0] = pos.d[0] < pBb->min.d[0] ? pos.d[0] : pBb->min.d[0];
	pBb->min.d[1] = pos.d[1] < pBb->min.d[1] ? pos.d[1] : pBb->min.d[1];
	pBb->max.d[0] = pos.d[0] > pBb->max.d[0] ? pos.d[0] : pBb->max.d[0];
	pBb->max.d[1] = pos.d[1] > pBb->max.d[1] ? pos.d[1] : pBb->max.d[1];
	if (pBb->border == -1) {
		pBb->border = border;
	}
}

static
PixErr walkAndAddBorder(
	const PixalcFPtrs *pAlloc,
	PixmshSplitMem *pMem,
	const PixmshSplitIntfIn *pMesh,
	PixmshSplitIntfOut *pIslands,
	PixmshBorderNode *pStart,
	I32 *pIslandIdx,
	I32 idx
) {
	PixErr err = PIX_ERR_SUCCESS;
	I32 islandIdx = pIslandIdx[idx];
	PIX_ERR_ASSERT("", islandIdx >= 0);
	if (seenThisEdge(pStart, islandIdx)) {
		return err;
	}
	PixmshBorderNode *pNode = pStart;
	I32 borderIdx = 0;
	err = pIslands->fpBorderInit(pAlloc, pIslands->pUserData, islandIdx, &borderIdx);
	PIX_ERR_RETURN_IFNOT(err, "");
	PixmshFaceCorner corner = pStart->corners[idx];
	PixmshFaceRange face = pMesh->fpFaceRange(pMesh->pUserData, corner.face);
	I32 edge = pMesh->fpEdge(pMesh->pUserData, corner);
	do {
		PIX_ERR_ASSERT("", !seenThisEdge(pNode, islandIdx));
		markEdgeSeen(pNode, corner, islandIdx);
		pIslands->fpBorderAddEdge(
			pAlloc,
			pIslands->pUserData,
			islandIdx,
			pIslandIdx[!idx],
			borderIdx,
			corner,
			edge
		);
		PIX_ERR_RETURN_IFNOT(err, "");
		pixmshBorderBbCmp(
			pMesh,
			pMem->bb.pArr + islandIdx,
			face.start + corner.corner,
			borderIdx
		);

		corner.corner = (corner.corner + 1) % (face.size);
		edge = pMesh->fpEdge(pMesh->pUserData, corner);
		if (edge < pMem->edgeTable.size && pMem->edgeTable.pArr[edge].valid ||
			pMesh->fpAdjCorner(pMesh->pUserData, corner).face == -1
		) {
			pNode = pMem->edges.pArr + pMem->edgeTable.pArr[edge].idx;
			if (!isEdgeIntern(pMem, pNode, NULL)) {
				continue;
			}
		}
		do {
			corner = pMesh->fpAdjCorner(pMesh->pUserData, corner);
			face = pMesh->fpFaceRange(pMesh->pUserData, corner.face);
			corner.corner = (corner.corner + 1) % (face.size);
			edge = pMesh->fpEdge(pMesh->pUserData, corner);
			pNode = pMem->edgeTable.pArr[edge].valid ?
				pMem->edges.pArr + pMem->edgeTable.pArr[edge].idx : NULL;
		} while(!pNode || isEdgeIntern(pMem, pNode, NULL));
	} while(pNode != pStart);
	return err;
}

static
void edgeTableAdd(const PixalcFPtrs *pAlloc, PixmshSplitIdxTableArr *pTable, I32 edge, I32 idx) {
	I32 oldSize = pTable->size;
	PIXALC_DYN_ARR_RESIZE(PixmshSplitIdxTable, pAlloc, pTable, edge + 1);
	if (oldSize <= edge) {
		memset(
			pTable->pArr + oldSize,
			0,
			sizeof(PixmshSplitIdxTable) * (pTable->size - oldSize)
		);
	}
	pTable->pArr[edge] = (PixmshSplitIdxTable){.idx = idx, .valid = true};
}

static
PixErr findAdjForCorner(
	const PixalcFPtrs *pAlloc,
	PixmshSplitMem *pMem,
	const PixmshSplitIntfIn *pMesh,
	bool (*fpSplitPredicate)(const void *, I32),
	PixmshFaceCorner corner,
	I32 *pSplitTotal
) {
	PixErr err = PIX_ERR_SUCCESS;
	I32 edge = pMesh->fpEdge(pMesh->pUserData, corner);
	if (edge < pMem->edgeTable.count && pMem->edgeTable.pArr[edge].valid) {
		return err;
	}
	PixmshEdgeCorners corners = pMesh->fpEdgeCorners(pMesh->pUserData, edge);
	I32 faces[2] = {corners.corners[0].face, corners.corners[1].face};
	bool borderEdge = faces[0] == -1 || faces[1] == -1;
	if (borderEdge || fpSplitPredicate && fpSplitPredicate(pMesh->pUserData, edge)) {
		++*pSplitTotal;
		if (faces[0] != -1 && !pMem->faceTable.pArr[faces[0]].valid) {
			islandIdxInit(pAlloc, &pMem->faceTable, &pMem->redirArr, faces[0]);
		}
		if (faces[1] != -1 && !pMem->faceTable.pArr[faces[1]].valid) {
			islandIdxInit(pAlloc, &pMem->faceTable, &pMem->redirArr, faces[1]);
		}
		if (edge >= pMem->edgeTable.size || !pMem->edgeTable.pArr[edge].valid) {
			PixmshBorderNode *pNode = pixmshBorderNodeInit(pAlloc, &pMem->edges, &corners);
			edgeTableAdd(pAlloc, &pMem->edgeTable, edge, pNode->idx);
		}
		return err;
	}
	bool joinTo;
	if (!pMem->faceTable.pArr[faces[0]].valid && !pMem->faceTable.pArr[faces[1]].valid) {
		joinTo = 0;
		islandIdxInit(pAlloc, &pMem->faceTable, &pMem->redirArr, faces[joinTo]);
	}
	else {
		joinTo = pMem->faceTable.pArr[faces[1]].valid;
	}
	if (!pMem->faceTable.pArr[faces[!joinTo]].valid) {
		islandIdxInit(pAlloc, &pMem->faceTable, &pMem->redirArr, faces[!joinTo]);
	}
	PixmshIdxRedir *pIdTo = getBuf(&pMem->faceTable, &pMem->redirArr, faces[joinTo]);
	PixmshIdxRedir *pIdFrom = getBuf(&pMem->faceTable, &pMem->redirArr, faces[!joinTo]);
	if (pIdTo != pIdFrom) {
		pIdFrom->idx = pIdTo->idx;
		pIdFrom->redir = true;
	}
	return err;
}

static
void pixmshSplitMemInit(const PixalcFPtrs *pAlloc, PixmshSplitMem *pMem, I32 faceCount) {
	PIXALC_DYN_ARR_RESIZE(PixmshSplitIdxTable, pAlloc, &pMem->faceTable, faceCount);
	if (pMem->faceTable.pArr) {
		memset(pMem->faceTable.pArr, 0, sizeof(PixmshSplitIdxTable) * faceCount);
	}
	if (pMem->edgeTable.pArr) {
		memset(pMem->edgeTable.pArr, 0, sizeof(PixmshSplitIdxTable) * pMem->edgeTable.size);
	}
	pMem->faceBuf.count = 0;
	pMem->redirArr.count = 0;
	pMem->faceTable.count = 0;
	pMem->edgeTable.count = 0;
	pMem->edges.count = 0;
	pMem->bb.count = 0;
}

//TODO reuse memory across multiple calls for tables, buffers
PixErr pixmshSplitToIslands(
	const PixalcFPtrs *pAlloc,
	PixmshSplitMem *pMem,
	const PixmshSplitIntfIn *pMesh,
	PixmshSplitIntfOut *pIslands,
	bool (*fpSplitPredicate)(const void *, I32)
) {
	PixErr err = PIX_ERR_SUCCESS;
	pixmshSplitMemInit(pAlloc, pMem, pMesh->faceCount);
	I32 splitTotal = 0;
	for (I32 i = 0; i < pMesh->faceCount; ++i) {
		PixmshFaceRange face = pMesh->fpFaceRange(pMesh->pUserData, i);
		for (I32 j = 0; j < face.size; ++j) {
			PixmshFaceCorner corner = {.face = i, .corner = j};
			err = findAdjForCorner(
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
	{
		I32 oldSize = pMem->faceBuf.size;
		PIXALC_DYN_ARR_RESIZE(PixmshFaceBuf, pAlloc, &pMem->faceBuf, pMem->redirArr.count);
		if (pMem->faceBuf.size > oldSize) {
			memset(
				pMem->faceBuf.pArr + oldSize,
				0,
				sizeof(PixmshFaceBuf) * (pMem->faceBuf.size - oldSize)
			);
		}
	}
	for (I32 i = 0; i < pMem->redirArr.count; ++i) {
		pMem->faceBuf.pArr[i].faces.count = 0;
	}
	for (I32 i = 0; i < pMesh->faceCount; ++i) {
		PIX_ERR_THROW_IFNOT_COND(err, pMem->faceTable.pArr[i].valid, "", 0);
		PixmshIdxRedir *pId = getBuf(&pMem->faceTable, &pMem->redirArr, i);
		PixmshFaceBuf *pBuf = pMem->faceBuf.pArr + pId->idx;
		I32 newIdx = 0;
		PIXALC_DYN_ARR_ADD(I32, pAlloc, &pBuf->faces, newIdx);
		pBuf->faces.pArr[newIdx] = i;
	}
	I32 *pFaces = NULL;
	err = pIslands->fpFacesInit(pAlloc, pIslands->pUserData, pMesh->faceCount, &pFaces);
	PIX_ERR_THROW_IFNOT_COND(err, pFaces, "", 0);
	I32 offset = 0;
	I32 islandCount = 0;
	for (I32 i = 0; i < pMem->redirArr.count; ++i) {
		if (!pMem->faceBuf.pArr[i].faces.count) {
			continue;
		}
		I32 newIdx = 0;
		err = pIslands->fpIslandAdd(pAlloc, pIslands->pUserData, splitTotal, &newIdx);
		PIX_ERR_THROW_IFNOT(err, "", 0);
		pMem->faceBuf.pArr[i].island = newIdx;
		PixtyRange range = {.start = offset};
		memcpy(
			pFaces + offset,
			pMem->faceBuf.pArr[i].faces.pArr,
			sizeof(I32) * pMem->faceBuf.pArr[i].faces.count
		);
		offset += pMem->faceBuf.pArr[i].faces.count;
		range.end = offset;
		err = pIslands->fpRangeSet(pIslands->pUserData, newIdx, range);
		PIX_ERR_THROW_IFNOT(err, "", 0);
		++islandCount;
	}
	PIX_ERR_ASSERT("", islandCount >= 0 && offset == pMesh->faceCount);
	PIXALC_DYN_ARR_RESIZE(PixmshBorderBb, pAlloc, &pMem->bb, islandCount);
	for (I32 i = 0; i < islandCount; ++i) {
		pMem->bb.pArr[i] = (PixmshBorderBb){
			.min = {FLT_MAX, FLT_MAX},
			.max = {-FLT_MAX, -FLT_MAX},
			.border = -1
		};
	}
	PIX_ERR_ASSERT("", pMem->edges.count > 0);
	for (I32 i = 0; i < pMem->edges.count; ++i) {
		PixmshBorderNode *pStart = pMem->edges.pArr + i;
		I32 islands[2] = {0};
		if (isEdgeIntern(pMem, pStart, islands)) {
			continue;
		}
		for (I32 j = 0; j < 2; ++j) {
			if (islands[j] == -1) {
				continue;
			}
			err = walkAndAddBorder(pAlloc, pMem, pMesh, pIslands, pStart, islands, j);
			PIX_ERR_THROW_IFNOT(err, "", 0);
		}
	}
	for (I32 i = 0; i < islandCount; ++i) {
		PIX_ERR_ASSERT("", pMem->bb.pArr[i].border != -1);
		if (pIslands->fpBorderMarkAsOuter) {
			err = pIslands->fpBorderMarkAsOuter(
				pIslands->pUserData,
				i,
				pMem->bb.pArr[i].border,
				&(PixmshV2Bb){.min = pMem->bb.pArr[i].min, .max = pMem->bb.pArr[i].max}
			);
		}
	}
	PIX_ERR_CATCH(0, err, ;);
	return err;
}

void pixmshSplitMemDestroy(const PixalcFPtrs *pAlloc, PixmshSplitMem *pMem) {
	if (pMem->faceBuf.pArr) {
		for (I32 i = 0; i < pMem->faceBuf.size; ++i) {
			if (pMem->faceBuf.pArr[i].faces.pArr) {
				pAlloc->fpFree(pMem->faceBuf.pArr[i].faces.pArr);
			}
		}
		pAlloc->fpFree(pMem->faceBuf.pArr);
	}
	if (pMem->redirArr.pArr) {
		pAlloc->fpFree(pMem->redirArr.pArr);
	}
	if (pMem->faceTable.pArr) {
		pAlloc->fpFree(pMem->faceTable.pArr);
	}
	if (pMem->edgeTable.pArr) {
		pAlloc->fpFree(pMem->edgeTable.pArr);
	}
	if (pMem->edges.pArr) {
		pAlloc->fpFree(pMem->edges.pArr);
	}
	if (pMem->bb.pArr) {
		pAlloc->fpFree(pMem->bb.pArr);
	}
	*pMem = (PixmshSplitMem){0};
}
