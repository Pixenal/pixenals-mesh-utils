/* 
SPDX-FileCopyrightText: 2025 Caleb Dawson
SPDX-License-Identifier: Apache-2.0
*/

#include <float.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include <mikktspace.h>

#include <utils.h>
#include <math_utils.h>
#include <context.h>
#include <attrib_utils.h>
#include <error.h>
#include <thread_pool.h>

BBox stucBBoxGet(const Mesh *pMesh, FaceRange *pFace) {
	BBox bbox = {.min = {.d = {FLT_MAX, FLT_MAX}}, .max = {.d = {-FLT_MAX, -FLT_MAX}}};
	for (I32 i = 0; i < pFace->size; ++i) {
		I32 vertIdx = pMesh->core.pCorners[pFace->start + i];
		V3_F32 pos = pMesh->pPos[vertIdx];
		//STUC_ASSERT("", pVert && v3F32IsFinite(*pVert)); //TODO validate this earlier
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
	STUC_ASSERT("", _(bbox.max V2GREATEQL bbox.min));
	return bbox;
}

U32 stucFnvHash(const U8 *value, I32 valueSize, U32 size) {
	STUC_ASSERT("", value && valueSize > 0 && size > 0);
	U32 hash = 2166136261;
	for (I32 i = 0; i < valueSize; ++i) {
		hash ^= value[i];
		hash *= 16777619;
	}
	hash %= size;
	STUC_ASSERT("", hash >= 0);
	return hash;
}

void stucGetInFaceBounds(FaceBounds *pBounds, const V2_F32 *pUvs, FaceRange face) {
	STUC_ASSERT("", pBounds && pUvs);
	STUC_ASSERT("", face.size >= 3 && face.start >= 0);
	STUC_ASSERT("", face.end >= 0 && face.idx >= 0);
	pBounds->fBBox.min.d[0] = pBounds->fBBox.min.d[1] = FLT_MAX;
	pBounds->fBBox.max.d[0] = pBounds->fBBox.max.d[1] = -FLT_MAX;
	for (I32 i = 0; i < face.size; ++i) {
		const V2_F32 uv = pUvs[face.start + i];
		STUC_ASSERT("", v2F32IsFinite(uv));
		if (uv.d[0] < pBounds->fBBox.min.d[0]) {
			pBounds->fBBox.min.d[0] = uv.d[0];
		}
		if (uv.d[1] < pBounds->fBBox.min.d[1]) {
			pBounds->fBBox.min.d[1] = uv.d[1];
		}
		if (uv.d[0] > pBounds->fBBox.max.d[0]) {
			pBounds->fBBox.max.d[0] = uv.d[0];
		}
		if (uv.d[1] > pBounds->fBBox.max.d[1]) {
			pBounds->fBBox.max.d[1] = uv.d[1];
		}
	}
	STUC_ASSERT("", _(pBounds->fBBox.max V2GREATEQL pBounds->fBBox.min));
}

I32 stucIsEdgeSeam(const Mesh *pMesh, I32 edge) {
	STUC_ASSERT("", pMesh && pMesh->pEdgeFaces && pMesh->pEdgeCorners);
	V2_I32 faces = pMesh->pEdgeFaces[edge];
	if (faces.d[1] == -1) {
		return true;
	}
	FaceRange faceA = stucGetFaceRange(&pMesh->core, faces.d[0]);
	bool windA = stucCalcFaceWindFromUvs(&faceA, pMesh);
	FaceRange faceB = stucGetFaceRange(&pMesh->core, faces.d[1]);
	bool windB = stucCalcFaceWindFromUvs(&faceB, pMesh);
	if (windA != windB) {
		return true; //marking wind borders as seam for now
	}
	V2_I8 corners = pMesh->pEdgeCorners[edge];
	I32 aA = corners.d[0];
	I32 bA = corners.d[1];
	I32 bB = stucGetCornerNext(bA, &faceB);
	V2_F32 uvAA = pMesh->pUvs[faceA.start + aA];
	V2_F32 uvBB = pMesh->pUvs[faceB.start + bB];
	if (!_(uvAA V2EQL uvBB)) {
		return true;
	}
	I32 aB = stucGetCornerNext(aA, &faceA);
	V2_F32 uvAB = pMesh->pUvs[faceA.start + aB];
	V2_F32 uvBA = pMesh->pUvs[faceB.start + bA];
	return !_(uvAB V2EQL uvBA);
	V2_F32 halfPlane = v2F32LineNormal(_(uvAB V2SUB uvAA));
	V2_F32 uvAC = pMesh->pUvs[faceA.start + stucGetCornerNext(aB, &faceA)];
	V2_F32 uvBC = pMesh->pUvs[faceB.start + stucGetCornerNext(bB, &faceB)];
	bool refSign = _(_(uvAC V2SUB uvAA) V2DOT halfPlane) > 0;
	if (_(_(uvBC V2SUB uvAA) V2DOT halfPlane) > 0 != refSign) {
		return true;
	}
	if (faceA.size == 4) {
		V2_F32 uvAD = pMesh->pUvs[faceA.start + stucGetCornerPrev(aA, &faceA)];
		if (_(_(uvAD V2SUB uvAA) V2DOT halfPlane) > 0 != refSign) {
			return true;
		}
	}
	if (faceB.size == 4) {
		V2_F32 uvBD = pMesh->pUvs[faceB.start + stucGetCornerPrev(bA, &faceB)];
		if (_(_(uvBD V2SUB uvAA) V2DOT halfPlane) > 0 != refSign) {
			return true;
		}
	}
	return false;
}

bool stucGetIfPreserveEdge(const Mesh *pMesh, I32 edge) {
	STUC_ASSERT("", pMesh && edge >= 0);
	if (pMesh->pEdgePreserve) {
		STUC_ASSERT("", pMesh->pEdgePreserve[edge] % 2 == pMesh->pEdgePreserve[edge]);
	}
	return pMesh->pEdgePreserve ? pMesh->pEdgePreserve[edge] : false;
}

bool stucCheckIfVertIsPreserve(const Mesh *pMesh, I32 vert) {
	STUC_ASSERT("", pMesh && vert >= 0);
	bool preserveVert = pMesh->pVertPreserve ? pMesh->pVertPreserve[vert] : false;
	STUC_ASSERT("", pMesh->pNumAdjPreserve);
	I32 numAdjSeam = pMesh->pNumAdjPreserve[vert] & 0xf;
	I32 numAdjPreserve = pMesh->pNumAdjPreserve[vert] >> 4 & 0xf;
	STUC_ASSERT("", numAdjSeam <= 3 && numAdjPreserve <= 3);
	return
		preserveVert ||
		//if a vert is adj to both a seam & a preserve edge, we keep it.
		// this avoids split edges in the final mesh.
		//note that an edge is only preserve if not a seam.
		numAdjSeam && numAdjPreserve ||
		//verts with 1, 3 or more (not 2) adj preserve edges are also kept,
		//(they're junction points, removing them would cause holes)
		numAdjPreserve == 1 || numAdjPreserve == 3;
}

bool stucCheckIfEdgeIsReceive(const Mesh *pMesh, I32 edge, F32 receiveLen) {
	STUC_ASSERT("", pMesh && edge >= 0);
	if (pMesh->pEdgeReceive) {
		STUC_ASSERT("", pMesh->pEdgeReceive[edge] % 2 == pMesh->pEdgeReceive[edge]);
	}
	if (receiveLen >= .0f) {
		STUC_ASSERT("", pMesh->pEdgeLen);
		return pMesh->pEdgeLen[edge] >= receiveLen;
	}
	else if (pMesh->pEdgeReceive) {
		return pMesh->pEdgeReceive[edge];
	}
	return true;
}

/*
static
I32 getOtherVert(I32 i, I32 faceSize, I8 *pVertsRemoved) {
	STUC_ASSERT("", i >= 0 && faceSize >= 3 && i < faceSize && pVertsRemoved);
	I32 ib = (i + 1) % faceSize;
	//search from i + 1 to facesize, and if non found,
	//then run again from 0 to facesize. If non found then,
	//return error
	I32 attempts = 0;
	do {
		attempts++;
		for (; ib < faceSize; ++ib) {
			STUC_ASSERT("", pVertsRemoved[ib] >= 0);
			if (!pVertsRemoved[ib]) {
				return ib;
			}
		}
		ib = 0;
	} while (attempts == 1);
	return -1;
}
*/

/*
typedef struct TriEdge {
	struct TriEdge* pNext;
	I32 tris[2];
	I32 verts[2];
	bool valid;
} TriEdge;

static
void initTriEdgeEntry(TriEdge* pEntry, I32 verta, I32 vertb, I32 tri) {
	pEntry->tris[0] = tri;
	pEntry->verts[0] = verta;
	pEntry->verts[1] = vertb;
	pEntry->valid = true;
}
*/

/*
static
void addTriEdgeToTable(
	const StucAlloc *pAlloc,
	I32 tableSize,
	TriEdge *pEdgeTable,
	I32 verta, I32 vertb,
	I32 tri
) {
	U32 sum = verta + vertb;
	I32 hash = stucFnvHash((U8 *)&sum, sizeof(sum), tableSize);
	TriEdge *pEntry = pEdgeTable + hash;
	if (!pEntry->valid) {
		initTriEdgeEntry(pEntry, verta, vertb, tri);
	}
	else {
		do {
			if ((pEntry->verts[0] == verta || pEntry->verts[0] == vertb) &&
				(pEntry->verts[1] == verta || pEntry->verts[1] == vertb)) {

				pEntry->tris[1] = tri;
				break;
			}
			if (!pEntry->pNext) {
				pEntry = pEntry->pNext = pAlloc->fpCalloc(1, sizeof(TriEdge));
				initTriEdgeEntry(pEntry, verta, vertb, tri);
				break;
			}
			pEntry = pEntry->pNext;
		} while(pEntry);
	}
}
*/

typedef struct AdjEntry {
	I32 face;
	I32 corner;
} AdjEntry;

typedef struct AdjBucket {
	AdjEntry* pArr;
	I32 count;
	I32 size;
} AdjBucket;

static
void adjTableDestroyBuckets(const StucAlloc *pAlloc, I32 count, AdjBucket *pAdjTable) {
	for (I32 i = 0; i < count; ++i) {
		if (pAdjTable[i].pArr) {
			pAlloc->fpFree(pAdjTable[i].pArr);
		}
	}
}

static
Result buildCornerAdjTable(
	const StucAlloc *pAlloc,
	const Mesh* pMesh,
	AdjBucket *pAdjTable
) {
	Result err = STUC_SUCCESS;
	STUC_ASSERT("", pAdjTable);
	for (I32 i = 0; i < pMesh->core.faceCount; ++i) {
		FaceRange face = stucGetFaceRange(&pMesh->core, i);
		for (I32 j = 0; j < face.size; ++j) {
			AdjBucket* pBucket = pAdjTable + pMesh->core.pCorners[face.start + j];
			STUC_ASSERT("", pBucket->count <= pBucket->size);
			if (!pBucket->pArr) {
				pBucket->size = 2;
				pBucket->pArr =
					pAlloc->fpMalloc(sizeof(AdjEntry) * pBucket->size);
			}
			else if (pBucket->count == pBucket->size) {
				STUC_ASSERT("tried to realloc null arr", pBucket->pArr);
				pBucket->size *= 2;
				pBucket->pArr = pAlloc->fpRealloc(
					pBucket->pArr,
					sizeof(AdjEntry) * pBucket->size
				);
			}
			pBucket->pArr[pBucket->count].face = i;
			pBucket->pArr[pBucket->count].corner = j;
			pBucket->count++;
		}
	}
	STUC_CATCH(0, err,
		adjTableDestroyBuckets(pAlloc, pMesh->core.vertCount, pAdjTable);
	;);
	return err;
}

static
Result findEdgesForFace(Mesh* pMesh, AdjBucket* pAdjTable, I32 idx) {
	Result err = STUC_SUCCESS;
	FaceRange face = stucGetFaceRange(&pMesh->core, idx);
	for (I32 j = 0; j < face.size; ++j) {
		if (pMesh->core.pEdges[face.start + j] >= 0) {
			continue; //Already set
		}
		I32 edge = pMesh->core.edgeCount;
		pMesh->core.edgeCount++;
		AdjBucket* pBucket = pAdjTable + pMesh->core.pCorners[face.start + j];
		STUC_ASSERT("", pBucket->count > 0 && pBucket->size >= pBucket->count);
		for (I32 k = 0; k < pBucket->count; ++k) {
			AdjEntry* pEntry = pBucket->pArr + k;
			if (pEntry->face == idx) {
				STUC_RETURN_ERR_IFNOT_COND(
					err,
					pEntry->corner == j,
					"Invalid mesh, 2 corners in this face share 1 vert"
				);
				continue;
			}
			FaceRange otherFace = stucGetFaceRange(&pMesh->core, pEntry->face);
			I32 nextCorner = (j + 1) % face.size;
			I32 otherPrevCorner = pEntry->corner ?
				pEntry->corner - 1 : otherFace.size - 1;
			if (pMesh->core.pEdges[otherFace.start + otherPrevCorner] >= 0) {
				continue; //Already set
			}
			if (pMesh->core.pCorners[face.start + nextCorner] !=
				pMesh->core.pCorners[otherFace.start + otherPrevCorner]) {
				continue; //Not connected
			}
			pMesh->core.pEdges[otherFace.start + otherPrevCorner] = edge;
			break;
		}
		pMesh->core.pEdges[face.start + j] = edge;
	}
	return err;
}

static
Result findEdges(Mesh* pMesh, AdjBucket* pAdjTable) {
	Result err = STUC_SUCCESS;
	for (I32 i = 0; i < pMesh->core.faceCount; ++i) {
		err = findEdgesForFace(pMesh, pAdjTable, i);
		STUC_RETURN_ERR_IFNOT(err, "");
	}
	return err;
}

Result stucBuildEdgeList(StucContext pCtx, Mesh* pMesh) {
	Result err = STUC_SUCCESS;
	STUC_RETURN_ERR_IFNOT_COND(err, !pMesh->core.pEdges, "");
	const StucAlloc *pAlloc = &pCtx->alloc;
	STUC_ASSERT("", pMesh->core.vertCount);
	AdjBucket* pAdjTable =
		pAlloc->fpCalloc(pMesh->core.vertCount, sizeof(AdjBucket));
	err = buildCornerAdjTable(pAlloc, pMesh, pAdjTable);
	STUC_THROW_IFNOT(err, "", 0);
	{
		STUC_ASSERT("", pMesh->core.cornerCount);
		I32 dataSize = sizeof(I32) * pMesh->core.cornerCount;
		pMesh->core.pEdges = pAlloc->fpMalloc(dataSize);
		memset(pMesh->core.pEdges, -1, dataSize);
		err = findEdges(pMesh, pAdjTable);
		STUC_THROW_IFNOT(err, "'findEdges' returned error", 1);
		STUC_CATCH(1, err,
			pAlloc->fpFree(pMesh->core.pEdges);
		);
	}
	STUC_CATCH(0, err, ;);
	adjTableDestroyBuckets(pAlloc, pMesh->core.vertCount, pAdjTable);
	pAlloc->fpFree(pAdjTable);
	return err;
}

void stucProgressBarClear() {
	printf("\n");
	printf("\x1b[1F");
	printf("\x1b[2K");
}

void stucProgressBarPrint(StucContext pCtx, I32 progress) {
	printf("	");
	for (I32 i = 0u; i < pCtx->stageReport.outOf; ++i) {
		char character;
		if (i < progress) {
			character = '#';
		}
		else {
			character = '-';
		}
		printf("%c", character);
	}
}

void stucStageBegin(void *pCtx, StucStageReport *pReport, const char* pName) {
	return;
	stucSetStageName(pCtx, pName);
}
void stucStageProgress(void *pCtx, StucStageReport *pReport, I32 progress) {
	return;
	if (progress) {
		stucProgressBarClear();
	}
	printf("%s", pReport->stage);
	stucProgressBarPrint(pCtx, progress);
}
void stucStageEnd(void *pCtx, StucStageReport *pReport) {
	return;
	memset(pReport->stage, 0, STUC_STAGE_NAME_LEN);
	stucProgressBarClear();
}

void stucStageBeginWrap(StucContext pCtx, const char* pName, I32 max) {
	pCtx->stageReport.fpBegin(pCtx, &pCtx->stageReport, pName);
	//Only needed if using default stage report functions,
	//it's just used for the progress bar
	pCtx->stageInterval = max <= pCtx->stageReport.outOf ?
		1 : max / pCtx->stageReport.outOf;
}

void stucStageProgressWrap(StucContext pCtx, I32 progress) {
	if (pCtx->stageInterval != 1 && progress % pCtx->stageInterval) {
		return;
	}
	//Normalize progress within stageReport.outOf
	I32 normProgress = progress / pCtx->stageInterval;
	pCtx->stageReport.fpProgress(pCtx, &pCtx->stageReport, normProgress);
}

void stucStageEndWrap(StucContext pCtx) {
	pCtx->stageReport.fpEnd(pCtx, &pCtx->stageReport);
}

void stucSetStageName(StucContext pCtx, const char* pName) {
	strncpy(pCtx->stageReport.stage, pName, STUC_STAGE_NAME_LEN);
}

Mat3x3 stucBuildFaceTbn(FaceRange face, const Mesh *pMesh, const I32 *pCornerOveride) {
	I32 corner = pCornerOveride ? face.start + pCornerOveride[1] : face.start;
	I32 vertIdx = pMesh->core.pCorners[corner];
	V2_F32 uv = pMesh->pUvs[corner];
	V3_F32 vert = pMesh->pPos[vertIdx];
	I32 next = pCornerOveride ? face.start + pCornerOveride[2] : face.start + 1;
	I32 vertIdxNext = pMesh->core.pCorners[next];
	V2_F32 uvNext = pMesh->pUvs[next];
	V3_F32 vertNext = pMesh->pPos[vertIdxNext];
	I32 prev = pCornerOveride ? face.start + pCornerOveride[0] : face.end - 1;
	I32 vertIdxPrev = pMesh->core.pCorners[prev];
	V2_F32 uvPrev = pMesh->pUvs[prev];
	V3_F32 vertPrev = pMesh->pPos[vertIdxPrev];
	//uv space direction vectors,
	//forming the coefficient matrix
	Mat2x2 coeffMat = {0};
	*(V2_F32 *)&coeffMat.d[0] = _(uvNext V2SUB uv);
	*(V2_F32 *)&coeffMat.d[1] = _(uvPrev V2SUB uv);
	//object space direction vectors,
	//forming the variable matrix
	Mat2x3 varMat = {0};
	V3_F32 osDirA = _(vertNext V3SUB vert);
	V3_F32 osDirB = _(vertPrev V3SUB vert);
	*(V3_F32 *)&varMat.d[0] = osDirA;
	*(V3_F32 *)&varMat.d[1] = osDirB;
	Mat2x2 coeffMatInv = mat2x2Invert(coeffMat);
	Mat2x3 tb = mat2x2MultiplyMat2x3(coeffMatInv, varMat);
	Mat3x3 tbn = {0};
	*(V3_F32 *)&tbn.d[0] = v3F32Normalize(*(V3_F32 *)&tb.d[0]);
	*(V3_F32 *)&tbn.d[1] = v3F32Normalize(*(V3_F32 *)&tb.d[1]);
	V3_F32 normal = _(osDirA V3CROSS osDirB);
	*(V3_F32 *)&tbn.d[2] = v3F32Normalize(normal);
	return tbn;
}

void stucGetTriScale(I32 size, BaseTriVerts *pTri) {
	for (I32 i = 0; i < size; ++i) {
		I32 iLast = i == 0 ? size - 1 : i - 1;
		I32 iNext = (i + 1) % size;
		F32 uvArea = v2F32TriArea(pTri->uv[iLast], pTri->uv[i], pTri->uv[iNext]);
		F32 xyzArea = v3F32TriArea(pTri->xyz[iLast], pTri->xyz[i], pTri->xyz[iNext]);
		pTri->scale[i] = xyzArea / uvArea;
	}
}

F32 stucGetT(V2_F32 point, V2_F32 lineA, V2_F32 lineUnit, F32 lineLen) {
	V2_F32 dir = _(point V2SUB lineA);
	return _(dir V2DOT lineUnit) / lineLen;
}

//kind of specific to this lib,
//a and b are v3, while and d are v2.
//cd is also taken as a param, while ab is calced
bool stucCalcIntersection(
	V3_F32 a,
	V3_F32 b,
	V2_F32 c,
	V2_F32 cd,
	V3_F32 *pPoint,
	F32 *pt,
	F32 *pt2
) {
	V3_F32 ab = _(b V3SUB a);
	V2_F32 ac = _(c V2SUB *(V2_F32 *)&a);
	F32 det2 = _(*(V2_F32 *)&ab V2DET cd);
	if (det2 == .0f) {
		return false;
	}
	F32 t = _(ac V2DET cd) / det2;
	if (pPoint) {
		*pPoint = _(a V3ADD _(ab V3MULS t));
	}
	if (pt) {
		*pt = t;
	}
	if (pt2) {
		det2 = _(cd V2DET *(V2_F32 *)&ab);
		if (det2 == .0f) {
			return false;
		}
		V2_F32 ca = _(*(V2_F32 *)&a V2SUB c);
		*pt2 = _(ca V2DET *(V2_F32 *)&ab) / det2;
	}
	return true;
}

//does not bounds check
I32 stucIdxBitArray(UBitField8 *pArr, I32 idx, I32 len) {
	idx *= len;
	I32 byte = idx / 8;
	I32 bit = idx % 8;
	I32 mask = (0x1 << len) - 1;
	if (bit + len > 8) {
		//bit spans byte boundary
		return *(UBitField16 *)&pArr[byte] >> bit & mask;
	}
	else {
		return pArr[byte] >> bit & mask;
	}
}

//does not bounds check.
//Also, if value is 0, only 1 bit will be set, len is ignored
void stucSetBitArr(UBitField8 *pArr, I32 idx, I32 value, I32 len) {
	STUC_ASSERT("", (value & (0x1 << len) - 1) == value);
	idx *= len;
	I32 byte = idx / 8;
	I32 bit = idx % 8;
	if (value) {
		if (bit + len > 8) {
			//cast to 16 bit as value spans across byte boundary
			*(UBitField16 *)&pArr[byte] |= value << bit;
		}
		else {
			pArr[byte] |= value << bit;
		}
	}
	else {
		UBitField8 mask = -0x1 ^ (0x1 << bit);
		pArr[byte] &= mask;
	}
}
#ifndef TEMP_DISABLE
void stucSetBorderFaceMapAttrib(
	BorderFace *pEntry,
	UBitField8 *pArr,
	I32 corner,
	I32 value
) {
	I32 len = 3 + pEntry->memType;
	stucSetBitArr(pArr, corner, value, len);
}
#endif

/*
void stucFInsertionSort(I32 *pIdxTable, I32 count, F32 *pSort) {
	//insertion sort
	F32 a = pSort[0];
	F32 b = pSort[1];
	I32 order = a < b;
	pIdxTable[0] = !order;
	pIdxTable[1] = order;
	I32 bufSize = 2;
	for (I32 i = bufSize; i < count; ++i) {
		bool insert = false;
		I32 j;
		for (j = bufSize - 1; j >= 0; --j) {
			insert = pSort[i] < pSort[pIdxTable[j]] &&
			         pSort[i] > pSort[pIdxTable[j - 1]];
			if (insert) {
				break;
			}
		}
		if (!insert) {
			pIdxTable[bufSize] = i;
		}
		else {
			for (I32 m = bufSize; m > j; --m) {
				pIdxTable[m] = pIdxTable[m - 1];
				STUC_ASSERT("", m <= bufSize && m > j);
			}
			pIdxTable[j] = i;
		}
		bufSize++;
	}
}
*/

Mat3x3 stucGetInterpolatedTbn(
	const Mesh *pMesh,
	const FaceRange *pFace,
	const I8 *pTriCorners,
	V3_F32 bc
) {
	//TODO replace interpolation in this func with the attrib
	//     interpolation funcions or macros
	V3_F32 *pNormals = pMesh->pNormals;
	V3_F32 normal = _(pNormals[pFace->start + pTriCorners[0]] V3MULS bc.d[0]);
	_(&normal V3ADDEQL _(pNormals[pFace->start + pTriCorners[1]] V3MULS bc.d[1]));
	_(&normal V3ADDEQL _(pNormals[pFace->start + pTriCorners[2]] V3MULS bc.d[2]));
	_(&normal V3DIVEQLS bc.d[0] + bc.d[1] + bc.d[2]);
	V3_F32 *pTangents = pMesh->pTangents;
	V3_F32 tangent = _(pTangents[pFace->start + pTriCorners[0]] V3MULS bc.d[0]);
	_(&tangent V3ADDEQL _(pTangents[pFace->start + pTriCorners[1]] V3MULS bc.d[1]));
	_(&tangent V3ADDEQL _(pTangents[pFace->start + pTriCorners[2]] V3MULS bc.d[2]));
	_(&tangent V3DIVEQLS bc.d[0] + bc.d[1] + bc.d[2]);
	//TODO should this be interpolated? Or are such edge cases invalid?
	F32 tSign = pMesh->pTSigns[pFace->start + pTriCorners[0]];
	V3_F32 bitangent = _(_(normal V3CROSS tangent) V3MULS tSign);
	Mat3x3 tbn = {0};
	*(V3_F32 *)&tbn.d[0] = tangent;
	*(V3_F32 *)&tbn.d[1] = bitangent;
	*(V3_F32 *)&tbn.d[2] = normal;
	return tbn;
}

I32 stucGetBorderFaceMemType(I32 mapFaceSize, I32 bufFaceSize) {
	STUC_ASSERT("", bufFaceSize >= 0);
	if (bufFaceSize <= 14 && mapFaceSize <= 8) {
		return 0;
	}
	else if (bufFaceSize <= 26 && mapFaceSize <= 16) {
		return 1;
	}
	else if (bufFaceSize <= 50 && mapFaceSize <= 32) {
		return 2;
	}
	STUC_ASSERT("Border face size > 64", false);
	return 0;
}

#ifndef TEMP_DISABLE
I32 stucGetBorderFaceSize(I32 memType) {
	STUC_ASSERT("", memType >= 0 && memType <= 3);
	switch (memType) {
	case 0:
		return sizeof(BorderFaceSmall);
	case 1:
		return sizeof(BorderFaceMid);
	case 2:
		return sizeof(BorderFaceLarge);
	}
	STUC_ASSERT("This shouldn't be hit", false);
	return 0;
}

Result stucAllocBorderFace(I32 memType, BorderTableAlloc *pHandles, void **ppOut) {
	Result err = STUC_SUCCESS;
	void *pHandle = NULL;
	switch (memType) {
		case 0:
			pHandle = pHandles->pSmall;
			break;
		case 1:
			pHandle = pHandles->pMid;
			break;
		case 2:
			pHandle = pHandles->pLarge;
			break;
		default:
			STUC_THROW(err, "invalid memtype", 0);
			break;
	}
	err = stucLinAlloc(pHandle, ppOut, 1);
	STUC_THROW_IFNOT(err, "error allocating border face entry", 0);
	STUC_CATCH(0, err, ;);
	return err;
}

void stucGetBorderFaceBitArrs(BorderFace *pEntry, BorderFaceBitArrs *pArrs) {
	switch (pEntry->memType) {
		case 0: {
			BorderFaceSmall *pCast = (BorderFaceSmall *)pEntry;
			pArrs->pBaseCorner = pCast->baseCorner;
			pArrs->pStucCorner = pCast->stucCorner;
			pArrs->pSegment = pCast->segment;
			pArrs->pIsStuc = pCast->isStuc;
			pArrs->pOnLine = pCast->onLine;
			pArrs->pOnInVert = pCast->onInVert;
			return;
		}
		case 1: {
			BorderFaceMid *pCast = (BorderFaceMid *)pEntry;
			pArrs->pBaseCorner = pCast->baseCorner;
			pArrs->pStucCorner = pCast->stucCorner;
			pArrs->pSegment = pCast->segment;
			pArrs->pIsStuc = pCast->isStuc;
			pArrs->pOnLine = pCast->onLine;
			pArrs->pOnInVert = pCast->onInVert;
			return;
		}
		case 2: {
			BorderFaceLarge *pCast = (BorderFaceLarge *)pEntry;
			pArrs->pBaseCorner = pCast->baseCorner;
			pArrs->pStucCorner = pCast->stucCorner;
			pArrs->pSegment = pCast->segment;
			pArrs->pIsStuc = pCast->isStuc;
			pArrs->pOnLine = pCast->onLine;
			pArrs->pOnInVert = pCast->onInVert;
			return;
		}
	}
}

void stucBorderTableDestroyAlloc(BorderTableAlloc *pTableAlloc) {
	if (pTableAlloc->pSmall) {
		stucLinAllocDestroy(pTableAlloc->pSmall);
	}
	if (pTableAlloc->pMid) {
		stucLinAllocDestroy(pTableAlloc->pMid);
	}
	if (pTableAlloc->pLarge) {
		stucLinAllocDestroy(pTableAlloc->pLarge);
	}
}

V3_F32 stucGetBufCornerUvw(
	const MapToMeshBasic *pBasic,
	const BufMesh *pMesh,
	I32 corner
) {
	V3_F32 bufUvw = { .d = {
		0, 0,
		pMesh->pW[corner] * pBasic->wScale}
	};
	*(V2_F32 *)&bufUvw = pMesh->mesh.pUvs[corner];
	return bufUvw;
}
#endif

static
Result sendOffJobs(
	const MapToMeshBasic *pBasic,
	I32 jobCount,
	void *pJobArgs, I32 argStructSize,
	Result (* func)(void *),
	void ***pppJobHandles
) {
	Result err = STUC_SUCCESS;
	StucContext pCtx = pBasic->pCtx;
	void *jobArgPtrs[MAX_THREADS] = {0};
	for (I32 i = 0; i < jobCount; ++i) {
		jobArgPtrs[i] = (U8 *)pJobArgs + i * argStructSize;
	}
	*pppJobHandles = pCtx->alloc.fpCalloc(jobCount, sizeof(void *));
	err = pCtx->threadPool.pJobStackPushJobs(
		pCtx->pThreadPoolHandle,
		jobCount,
		*pppJobHandles,
		func,
		jobArgPtrs
	);
	STUC_RETURN_ERR_IFNOT(err, "");
	return err;
}

Result stucDoJobInParallel(
	const MapToMeshBasic *pBasic,
	I32 jobCount, void *pJobArgs, I32 argStructSize,
	Result (* func)(void *)
) {
	Result err = STUC_SUCCESS;
	STUC_ASSERT("", jobCount >= 0);
	if (!jobCount) {
		return err;
	}
	void **ppJobHandles = NULL;
	err = sendOffJobs(pBasic, jobCount, pJobArgs, argStructSize, func, &ppJobHandles);
	STUC_THROW_IFNOT(err, "", 0);
	err = pBasic->pCtx->threadPool.fpWaitForJobs(
		pBasic->pCtx->pThreadPoolHandle,
		jobCount,
		ppJobHandles,
		true,
		NULL
	);
	STUC_THROW_IFNOT(err, "", 0);
	err = stucJobGetErrs(pBasic->pCtx, jobCount, &ppJobHandles);
	STUC_THROW_IFNOT(err, "", 0);
	STUC_CATCH(0, err, ;);
	if (ppJobHandles) {
		stucJobDestroyHandles(pBasic->pCtx, jobCount, ppJobHandles);
		pBasic->pCtx->alloc.fpFree(ppJobHandles);
	}
	return err;
}

/*
U32 stucGetEncasedFaceHash(I32 mapFace, V2_I16 tile, I32 tableSize) {
	U32 key = (U32)mapFace + stucFnvHash((U8 *)tile.d, sizeof(tile.d), UINT32_MAX);
	return stucFnvHash((U8 *)&key, sizeof(key), tableSize);
}
*/

InsideStatus stucIsPointInHalfPlane(
	V2_F32 point,
	V2_F32 lineA,
	V2_F32 halfPlane,
	bool wind
) {
	V2_F32 dir = _(point V2SUB lineA);
	F32 dot = _(halfPlane V2DOT dir);
	if (dot == .0f) {
		return STUC_INSIDE_STATUS_ON_LINE;
	}
	else {
		return (dot > .0f) ^ wind ? STUC_INSIDE_STATUS_INSIDE : STUC_INSIDE_STATUS_OUTSIDE;
	}
}

static
I32 getNearbyPrime(I32 num) {
	I32 primes[] = {
		1,
		3,
		5,
		11,
		17,
		37,
		67,
		131,
		257,
		521,
		1031,
		2053,
		4099,
		8209,
		16411,
		32771,
		65537,
		131101,
		262147,
		524309,
		1048583,
		2097169,
		4194319,
		8388617,
		16777259,
		33554467,
		67108879,
		134217757,
		268435459
	};
	F32 exp = log2f((F32)num);
	I32 expRound = roundf(exp);
	STUC_ASSERT("a value this high shouldn't've been passed", expRound <= 28);
	return primes[expRound];
}

void stucHTableInit(
	const StucAlloc *pAlloc,
	HTable *pHandle,
	I32 targetSize,
	I32Arr allocTypeSizes,
	void *pUserData
) {
	STUC_ASSERT("", targetSize > 0);
	I32 size = getNearbyPrime(targetSize);
	*pHandle = (HTable){
		.pAlloc = pAlloc,
		.pUserData = pUserData,
		.size = size,
		.pTable = pAlloc->fpCalloc(size, sizeof(HTableBucket))
	};
	STUC_ASSERT(
		"",
		allocTypeSizes.count > 0 && allocTypeSizes.count <= STUC_HTABLE_ALLOC_HANDLES_MAX
	);
	I32 allocInitSize = size / allocTypeSizes.count / 2 + 1;
	for (I32 i = 0; i < allocTypeSizes.count; ++i) {
		stucLinAllocInit(
			pAlloc,
			pHandle->allocHandles + i,
			allocTypeSizes.pArr[i],
			allocInitSize,
			true
		);
	}
}

void stucHTableDestroy(HTable *pHandle) {
	if (pHandle->pTable) {
		STUC_ASSERT("", pHandle->size);
		pHandle->pAlloc->fpFree(pHandle->pTable);
	}
	STUC_ASSERT(
		"at least 1 lin alloc handle should have been initialized",
		pHandle->allocHandles[0].valid
	);
	for (I32 i = 0; i < STUC_HTABLE_ALLOC_HANDLES_MAX; ++i) {
		if (pHandle->allocHandles[i].valid) {
			stucLinAllocDestroy(pHandle->allocHandles + i);
		}
	}
	*pHandle = (HTable) {0};
}

LinAlloc *stucHTableAllocGet(HTable *pHandle, I32 idx) {
	STUC_ASSERT("", idx >= 0 && idx < STUC_HTABLE_ALLOC_HANDLES_MAX);
	return pHandle->allocHandles + idx;
}

const LinAlloc *stucHTableAllocGetConst(const HTable *pHandle, I32 idx) {
	STUC_ASSERT("", idx >= 0 && idx < STUC_HTABLE_ALLOC_HANDLES_MAX);
	return pHandle->allocHandles + idx;
}

U64 stucKeyFromI32(const void *pKeyData) {
	return *(I32 *)pKeyData;
}

static
int mikktGetNumFaces(const SMikkTSpaceContext *pCtx) {
	TangentJobArgs *pArgs = pCtx->m_pUserData;
	return pArgs->faces.count;
}

static
I32 mikktGetInFace(const TangentJobArgs *pArgs, I32 iFace) {
	return pArgs->pTPieces->pInFaces[pArgs->core.range.start + iFace];
}

static
int mikktGetNumVertsOfFace(const SMikkTSpaceContext *pCtx, const int iFace) {
	TangentJobArgs *pArgs = pCtx->m_pUserData;
	I32 inFaceIdx = mikktGetInFace(pArgs, iFace);
	return stucGetFaceRange(&pArgs->core.pBasic->pInMesh->core, inFaceIdx).size;
}

static
I32 mikktGetInFaceStart(const TangentJobArgs *pArgs, I32 inFaceIdx) {
	return pArgs->core.pBasic->pInMesh->core.pFaces[inFaceIdx];
}

static
void mikktGetPos(
	const SMikkTSpaceContext *pCtx,
	F32 *pFvPosOut,
	const int iFace,
	const int iVert
) {
	TangentJobArgs *pArgs = pCtx->m_pUserData;
	I32 inFaceIdx = mikktGetInFace(pArgs, iFace);
	I32 inFaceStart = mikktGetInFaceStart(pArgs, inFaceIdx);
	I32 vertIdx = pArgs->core.pBasic->pInMesh->core.pCorners[inFaceStart + iVert];
	*(V3_F32 *)pFvPosOut = pArgs->core.pBasic->pInMesh->pPos[vertIdx];
}

static
void mikktGetNormal(
	const SMikkTSpaceContext *pCtx,
	F32 *pFvNormOut,
	const int iFace,
	const int iVert
) {
	TangentJobArgs *pArgs = pCtx->m_pUserData;
	I32 inFaceIdx = mikktGetInFace(pArgs, iFace);
	I32 inFaceStart = mikktGetInFaceStart(pArgs, inFaceIdx);
	*(V3_F32 *)pFvNormOut = pArgs->core.pBasic->pInMesh->pNormals[inFaceStart + iVert];
}

static
void mikktGetTexCoord(
	const SMikkTSpaceContext *pCtx,
	F32 *pFvTexcOut,
	const int iFace,
	const int iVert
) {
	TangentJobArgs *pArgs = pCtx->m_pUserData;
	I32 inFaceIdx = mikktGetInFace(pArgs, iFace);
	I32 inFaceStart = mikktGetInFaceStart(pArgs, inFaceIdx);
	*(V2_F32 *)pFvTexcOut = pArgs->core.pBasic->pInMesh->pUvs[inFaceStart + iVert];
}

static
void mikktSetTSpaceBasic(
	const SMikkTSpaceContext *pCtx,
	const F32 *pFvTangent,
	const F32 fSign,
	const int iFace,
	const int iVert
) {
	TangentJobArgs *pArgs = pCtx->m_pUserData;
	I32 corner = pArgs->faces.pArr[iFace] + iVert;
	pArgs->pTangents[corner] = *(V3_F32 *)pFvTangent;
	pArgs->pTSigns[corner] = fSign;
}

Result stucBuildTangents(void *pArgsVoid) {
	Result err = STUC_SUCCESS;
	TangentJobArgs *pArgs = pArgsVoid;
	const StucAlloc *pAlloc = &pArgs->core.pBasic->pCtx->alloc;

	pArgs->pTangents = pAlloc->fpCalloc(pArgs->cornerCount, sizeof(V3_F32));
	pArgs->pTSigns = pAlloc->fpCalloc(pArgs->cornerCount, sizeof(F32));

	SMikkTSpaceInterface mikktInterface = {
		.m_getNumFaces = mikktGetNumFaces,
		.m_getNumVerticesOfFace = mikktGetNumVertsOfFace,
		.m_getPosition = mikktGetPos,
		.m_getNormal = mikktGetNormal,
		.m_getTexCoord = mikktGetTexCoord,
		.m_setTSpaceBasic = mikktSetTSpaceBasic
	};
	SMikkTSpaceContext mikktContext = {
		.m_pInterface = &mikktInterface,
		.m_pUserData = pArgs
	};
	if (!genTangSpaceDefault(&mikktContext)) {
		STUC_RETURN_ERR(err, "mikktspace func 'genTangSpaceDefault' returned error");
	}
	return err;
}
