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
#include <pixenals_math_utils.h>
#include <context.h>
#include <attrib_utils.h>
#include <pixenals_error_utils.h>
#include <pixenals_thread_utils.h>

BBox stucBBoxGet(const Mesh *pMesh, FaceRange *pFace) {
	BBox bbox = {.min = {.d = {FLT_MAX, FLT_MAX}}, .max = {.d = {-FLT_MAX, -FLT_MAX}}};
	for (I32 i = 0; i < pFace->size; ++i) {
		I32 vertIdx = pMesh->core.pCorners[pFace->start + i];
		V3_F32 pos = pMesh->pPos[vertIdx];
		//PIX_ERR_ASSERT("", pVert && pixmV3F32IsFinite(*pVert)); //TODO validate this earlier
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

void stucGetInFaceBounds(FaceBounds *pBounds, const V2_F32 *pUvs, FaceRange face) {
	PIX_ERR_ASSERT("", pBounds && pUvs);
	PIX_ERR_ASSERT("", face.size >= 3 && face.start >= 0);
	PIX_ERR_ASSERT("", face.end >= 0 && face.idx >= 0);
	pBounds->fBBox.min.d[0] = pBounds->fBBox.min.d[1] = FLT_MAX;
	pBounds->fBBox.max.d[0] = pBounds->fBBox.max.d[1] = -FLT_MAX;
	for (I32 i = 0; i < face.size; ++i) {
		const V2_F32 uv = pUvs[face.start + i];
		PIX_ERR_ASSERT("", pixmV2F32IsFinite(uv));
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
	PIX_ERR_ASSERT("", _(pBounds->fBBox.max V2GREATEQL pBounds->fBBox.min));
}

I32 stucIsEdgeSeam(const Mesh *pMesh, I32 edge) {
	PIX_ERR_ASSERT("", pMesh && pMesh->pEdgeFaces && pMesh->pEdgeCorners);
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
	V2_F32 halfPlane = pixmV2F32LineNormal(_(uvAB V2SUB uvAA));
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
	PIX_ERR_ASSERT("", pMesh && edge >= 0);
	if (pMesh->pEdgePreserve) {
		PIX_ERR_ASSERT("", pMesh->pEdgePreserve[edge] % 2 == pMesh->pEdgePreserve[edge]);
	}
	return pMesh->pEdgePreserve ? pMesh->pEdgePreserve[edge] : false;
}

bool stucCheckIfVertIsPreserve(const Mesh *pMesh, I32 vert) {
	PIX_ERR_ASSERT("", pMesh && vert >= 0);
	bool preserveVert = pMesh->pVertPreserve ? pMesh->pVertPreserve[vert] : false;
	PIX_ERR_ASSERT("", pMesh->pNumAdjPreserve);
	I32 numAdjSeam = pMesh->pNumAdjPreserve[vert] & 0xf;
	I32 numAdjPreserve = pMesh->pNumAdjPreserve[vert] >> 4 & 0xf;
	PIX_ERR_ASSERT("", numAdjSeam <= 3 && numAdjPreserve <= 3);
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
	PIX_ERR_ASSERT("", pMesh && edge >= 0);
	if (pMesh->pEdgeReceive) {
		PIX_ERR_ASSERT("", pMesh->pEdgeReceive[edge] % 2 == pMesh->pEdgeReceive[edge]);
	}
	if (receiveLen >= .0f) {
		PIX_ERR_ASSERT("", pMesh->pEdgeLen);
		return pMesh->pEdgeLen[edge] >= receiveLen;
	}
	else if (pMesh->pEdgeReceive) {
		return pMesh->pEdgeReceive[edge];
	}
	return true;
}

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
StucErr buildCornerAdjTable(
	const StucAlloc *pAlloc,
	const StucMesh *pMesh,
	AdjBucket *pAdjTable
) {
	StucErr err = PIX_ERR_SUCCESS;
	PIX_ERR_ASSERT("", pAdjTable);
	for (I32 i = 0; i < pMesh->faceCount; ++i) {
		FaceRange face = stucGetFaceRange(pMesh, i);
		for (I32 j = 0; j < face.size; ++j) {
			AdjBucket* pBucket = pAdjTable + pMesh->pCorners[face.start + j];
			PIX_ERR_ASSERT("", pBucket->count <= pBucket->size);
			if (!pBucket->pArr) {
				pBucket->size = 2;
				pBucket->pArr =
					pAlloc->fpMalloc(sizeof(AdjEntry) * pBucket->size);
			}
			else if (pBucket->count == pBucket->size) {
				PIX_ERR_ASSERT("tried to realloc null arr", pBucket->pArr);
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
	PIX_ERR_CATCH(0, err,
		adjTableDestroyBuckets(pAlloc, pMesh->vertCount, pAdjTable);
	;);
	return err;
}

static
StucErr findEdgesForFace(StucMesh *pMesh, AdjBucket* pAdjTable, I32 idx) {
	StucErr err = PIX_ERR_SUCCESS;
	FaceRange face = stucGetFaceRange(pMesh, idx);
	for (I32 j = 0; j < face.size; ++j) {
		if (pMesh->pEdges[face.start + j] >= 0) {
			continue; //Already set
		}
		I32 edge = pMesh->edgeCount;
		pMesh->edgeCount++;
		AdjBucket* pBucket = pAdjTable + pMesh->pCorners[face.start + j];
		PIX_ERR_ASSERT("", pBucket->count > 0 && pBucket->size >= pBucket->count);
		for (I32 k = 0; k < pBucket->count; ++k) {
			AdjEntry* pEntry = pBucket->pArr + k;
			if (pEntry->face == idx) {
				PIX_ERR_RETURN_IFNOT_COND(
					err,
					pEntry->corner == j,
					"Invalid mesh, 2 corners in this face share 1 vert"
				);
				continue;
			}
			FaceRange otherFace = stucGetFaceRange(pMesh, pEntry->face);
			I32 nextCorner = (j + 1) % face.size;
			I32 otherPrevCorner = pEntry->corner ?
				pEntry->corner - 1 : otherFace.size - 1;
			if (pMesh->pEdges[otherFace.start + otherPrevCorner] >= 0) {
				continue; //Already set
			}
			if (pMesh->pCorners[face.start + nextCorner] !=
				pMesh->pCorners[otherFace.start + otherPrevCorner]) {
				continue; //Not connected
			}
			pMesh->pEdges[otherFace.start + otherPrevCorner] = edge;
			break;
		}
		pMesh->pEdges[face.start + j] = edge;
	}
	return err;
}

static
StucErr findEdges(StucMesh* pMesh, AdjBucket* pAdjTable) {
	StucErr err = PIX_ERR_SUCCESS;
	for (I32 i = 0; i < pMesh->faceCount; ++i) {
		err = findEdgesForFace(pMesh, pAdjTable, i);
		PIX_ERR_RETURN_IFNOT(err, "");
	}
	return err;
}

StucErr stucBuildEdgeList(StucContext pCtx, StucMesh *pMesh) {
	StucErr err = PIX_ERR_SUCCESS;
	PIX_ERR_RETURN_IFNOT_COND(err, !pMesh->pEdges, "");
	const StucAlloc *pAlloc = &pCtx->alloc;
	PIX_ERR_ASSERT("", pMesh->vertCount);
	AdjBucket* pAdjTable =
		pAlloc->fpCalloc(pMesh->vertCount, sizeof(AdjBucket));
	err = buildCornerAdjTable(pAlloc, pMesh, pAdjTable);
	PIX_ERR_THROW_IFNOT(err, "", 0);
	{
		PIX_ERR_ASSERT("", pMesh->cornerCount);
		I32 dataSize = sizeof(I32) * pMesh->cornerCount;
		pMesh->pEdges = pAlloc->fpMalloc(dataSize);
		memset(pMesh->pEdges, -1, dataSize);
		err = findEdges(pMesh, pAdjTable);
		PIX_ERR_THROW_IFNOT(err, "'findEdges' returned error", 1);
		PIX_ERR_CATCH(1, err,
			pAlloc->fpFree(pMesh->pEdges);
		);
	}
	PIX_ERR_CATCH(0, err, ;);
	adjTableDestroyBuckets(pAlloc, pMesh->vertCount, pAdjTable);
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

M3x3 stucBuildFaceTbn(FaceRange face, const Mesh *pMesh, const I32 *pCornerOveride) {
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
	M2x2 coeffMat = {0};
	*(V2_F32 *)&coeffMat.d[0] = _(uvNext V2SUB uv);
	*(V2_F32 *)&coeffMat.d[1] = _(uvPrev V2SUB uv);
	//object space direction vectors,
	//forming the variable matrix
	M2x3 varMat = {0};
	V3_F32 osDirA = _(vertNext V3SUB vert);
	V3_F32 osDirB = _(vertPrev V3SUB vert);
	*(V3_F32 *)&varMat.d[0] = osDirA;
	*(V3_F32 *)&varMat.d[1] = osDirB;
	M2x2 coeffMatInv = pixmM2x2Invert(coeffMat);
	M2x3 tb = pixmM2x2MultiplyM2x3(coeffMatInv, varMat);
	M3x3 tbn = {0};
	*(V3_F32 *)&tbn.d[0] = pixmV3F32Normalize(*(V3_F32 *)&tb.d[0]);
	*(V3_F32 *)&tbn.d[1] = pixmV3F32Normalize(*(V3_F32 *)&tb.d[1]);
	V3_F32 normal = _(osDirA V3CROSS osDirB);
	*(V3_F32 *)&tbn.d[2] = pixmV3F32Normalize(normal);
	return tbn;
}

void stucGetTriScale(I32 size, BaseTriVerts *pTri) {
	for (I32 i = 0; i < size; ++i) {
		I32 iLast = i == 0 ? size - 1 : i - 1;
		I32 iNext = (i + 1) % size;
		F32 uvArea = pixmV2F32TriArea(pTri->uv[iLast], pTri->uv[i], pTri->uv[iNext]);
		F32 xyzArea = pixmV3F32TriArea(pTri->xyz[iLast], pTri->xyz[i], pTri->xyz[iNext]);
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
	PIX_ERR_ASSERT("", (value & (0x1 << len) - 1) == value);
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

M3x3 stucGetInterpolatedTbn(
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
	M3x3 tbn = {0};
	*(V3_F32 *)&tbn.d[0] = tangent;
	*(V3_F32 *)&tbn.d[1] = bitangent;
	*(V3_F32 *)&tbn.d[2] = normal;
	return tbn;
}

I32 stucGetBorderFaceMemType(I32 mapFaceSize, I32 bufFaceSize) {
	PIX_ERR_ASSERT("", bufFaceSize >= 0);
	if (bufFaceSize <= 14 && mapFaceSize <= 8) {
		return 0;
	}
	else if (bufFaceSize <= 26 && mapFaceSize <= 16) {
		return 1;
	}
	else if (bufFaceSize <= 50 && mapFaceSize <= 32) {
		return 2;
	}
	PIX_ERR_ASSERT("Border face size > 64", false);
	return 0;
}

static
StucErr sendOffJobs(
	StucContext pCtx,
	I32 jobCount,
	void *pJobArgs, I32 argStructSize,
	StucErr (* func)(void *),
	void ***pppJobHandles
) {
	StucErr err = PIX_ERR_SUCCESS;
	void *jobArgPtrs[PIX_THREAD_MAX_THREADS] = {0};
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
	PIX_ERR_RETURN_IFNOT(err, "");
	return err;
}

StucErr stucDoJobInParallel(
	StucContext pCtx,
	I32 jobCount, void *pJobArgs, I32 argStructSize,
	StucErr (* func)(void *)
) {
	StucErr err = PIX_ERR_SUCCESS;
	PIX_ERR_ASSERT("", jobCount >= 0);
	if (!jobCount) {
		return err;
	}
	void **ppJobHandles = NULL;
	err = sendOffJobs(pCtx, jobCount, pJobArgs, argStructSize, func, &ppJobHandles);
	PIX_ERR_THROW_IFNOT(err, "", 0);
	err = pCtx->threadPool.fpWaitForJobs(
		pCtx->pThreadPoolHandle,
		jobCount,
		ppJobHandles,
		true,
		NULL
	);
	PIX_ERR_THROW_IFNOT(err, "", 0);
	err = stucJobGetErrs(pCtx, jobCount, &ppJobHandles);
	PIX_ERR_THROW_IFNOT(err, "", 0);
	PIX_ERR_CATCH(0, err, ;);
	if (ppJobHandles) {
		stucJobDestroyHandles(pCtx, jobCount, ppJobHandles);
		pCtx->alloc.fpFree(ppJobHandles);
	}
	return err;
}

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

StucErr stucThreadPoolSetCustom(
	StucContext pCtx,
	const StucThreadPool *pThreadPool
) {
	if (!pThreadPool->fpInit || !pThreadPool->fpDestroy || !pThreadPool->fpMutexGet ||
	    !pThreadPool->fpMutexLock || !pThreadPool->fpMutexUnlock || !pThreadPool->fpMutexDestroy ||
	    !pThreadPool->fpJobStackGetJob || !pThreadPool->pJobStackPushJobs) {
		printf("Failed to set custom thread pool. One or more functions were NULL");
		return PIX_ERR_ERROR;
	}
	pCtx->threadPool.fpDestroy(pCtx);
	pCtx->threadPool = *pThreadPool;
	return PIX_ERR_SUCCESS;
}

void stucThreadPoolSetDefault(StucContext pCtx) {
	pCtx->threadPool.fpInit = pixthThreadPoolInit;
	pCtx->threadPool.fpWaitForJobs = pixthWaitForJobsIntern;
	pCtx->threadPool.fpGetJobErr = pixthGetJobErr;
	pCtx->threadPool.fpJobHandleDestroy = pixthJobHandleDestroy;
	pCtx->threadPool.fpDestroy = pixthThreadPoolDestroy;
	pCtx->threadPool.fpMutexGet = pixthMutexGet;
	pCtx->threadPool.fpMutexLock = pixthMutexLock;
	pCtx->threadPool.fpMutexUnlock = pixthMutexUnlock;
	pCtx->threadPool.fpMutexDestroy = pixthMutexDestroy;
	/*
	pCtx->threadPool.fpBarrierGet = stucBarrierGet;
	pCtx->threadPool.fpBarrierWait = stucBarrierWait;
	pCtx->threadPool.fpBarrierDestroy = stucBarrierDestroy;
	*/
	pCtx->threadPool.fpJobStackGetJob = pixthJobStackGetJob;
	pCtx->threadPool.pJobStackPushJobs = pixthJobStackPushJobs;
	pCtx->threadPool.fpGetAndDoJob = pixthGetAndDoJob;
}

void stucAllocSetCustom(PixalcFPtrs *pAlloc, PixalcFPtrs *pCustomAlloc) {
	PIX_ERR_ASSERT("", pAlloc && pCustomAlloc);
	if (!pCustomAlloc->fpMalloc || !pCustomAlloc->fpCalloc || !pCustomAlloc->fpFree) {
		printf("Failed to set custom alloc. One or more functions were NULL");
		return;
	}
	*pAlloc = *pCustomAlloc;
}

void stucAllocSetDefault(PixalcFPtrs *pAlloc) {
	PIX_ERR_ASSERT("", pAlloc);
	pAlloc->fpMalloc = malloc;
	pAlloc->fpCalloc = calloc;
	pAlloc->fpFree = free;
	pAlloc->fpRealloc = realloc;
}
