#include <float.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include <utils.h>
#include <math_utils.h>
#include <context.h>
#include <attrib_utils.h>
#include <error.h>
#include <thread_pool.h>

I32 stucCheckFaceIsInBounds(V2_F32 min, V2_F32 max, FaceRange face, const Mesh *pMesh) {
	STUC_ASSERT("", pMesh && pMesh->pVerts && pMesh->core.pCorners);
	STUC_ASSERT("", face.size >= 3 && face.start >= 0 && face.end >= 0 && face.idx >= 0);
	STUC_ASSERT("", v2F32IsFinite(min) && v2F32IsFinite(max));
	V2_F32 faceMin = {0};
	V2_F32 faceMax = {0};
	faceMin.d[0] = faceMin.d[1] = FLT_MAX;
	faceMax.d[0] = faceMax.d[1] = 0;
	for (I32 i = 0; i < face.size; ++i) {
		I32 vertIdx = pMesh->core.pCorners[face.start + i];
		V3_F32 *pVert = pMesh->pVerts + vertIdx;
		STUC_ASSERT("", pVert && v3F32IsFinite(*pVert));
		if (pVert->d[0] < faceMin.d[0]) {
			faceMin.d[0] = pVert->d[0];
		}
		if (pVert->d[1] < faceMin.d[1]) {
			faceMin.d[1] = pVert->d[1];
		}
		if (pVert->d[0] > faceMax.d[0]) {
			faceMax.d[0] = pVert->d[0];
		}
		if (pVert->d[1] > faceMax.d[1]) {
			faceMax.d[1] = pVert->d[1];
		}
	}
	//Faces can be flat (they may be facing sideways in a map for instance)
	STUC_ASSERT("", _(faceMax V2GREATEQL faceMin));
	V2_I32 inside = {0};
	inside.d[0] =
		(faceMin.d[0] >= min.d[0] && faceMin.d[0] < max.d[0]) ||
		(faceMax.d[0] >= min.d[0] && faceMax.d[0] < max.d[0]) ||
		(faceMin.d[0] < min.d[0] && faceMax.d[0] >= max.d[0]);
	inside.d[1] =
		(faceMin.d[1] >= min.d[1] && faceMin.d[1] < max.d[1]) ||
		(faceMax.d[1] >= min.d[1] && faceMax.d[1] < max.d[1]) ||
		(faceMin.d[1] < min.d[1] && faceMax.d[1] >= max.d[1]);
	return inside.d[0] && inside.d[1];
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

void stucGetFaceBounds(FaceBounds *pBounds, const V2_F32 *pUvs, FaceRange face) {
	STUC_ASSERT("", pBounds && pUvs);
	STUC_ASSERT("", face.size >= 3 && face.start >= 0);
	STUC_ASSERT("", face.end >= 0 && face.idx >= 0);
	pBounds->fMin.d[0] = pBounds->fMin.d[1] = FLT_MAX;
	pBounds->fMax.d[0] = pBounds->fMax.d[1] = -FLT_MAX;
	for (I32 i = 0; i < face.size; ++i) {
		const V2_F32 *uv = pUvs + face.start + i;
		STUC_ASSERT("", uv && v2F32IsFinite(*uv));
		pBounds->fMin.d[0] = uv->d[0] < pBounds->fMin.d[0] ?
			uv->d[0] : pBounds->fMin.d[0];
		pBounds->fMin.d[1] = uv->d[1] < pBounds->fMin.d[1] ?
			uv->d[1] : pBounds->fMin.d[1];
		pBounds->fMax.d[0] = uv->d[0] > pBounds->fMax.d[0] ?
			uv->d[0] : pBounds->fMax.d[0];
		pBounds->fMax.d[1] = uv->d[1] > pBounds->fMax.d[1] ?
			uv->d[1] : pBounds->fMax.d[1];
	}
	STUC_ASSERT("", _(pBounds->fMax V2GREATEQL pBounds->fMin));
}

I32 stucCheckIfEdgeIsSeam(
	I32 edgeIdx,
	FaceRange face,
	I32 corner,
	const Mesh *pMesh
) {
	STUC_ASSERT("", pMesh);
	STUC_ASSERT("", face.size >= 3 && face.start >= 0 && face.end >= 0 && face.size >= 0);
	STUC_ASSERT("", edgeIdx >= 0 && corner >= 0 && corner < face.size);
	V2_I32 corners = pMesh->pEdgeCorners[edgeIdx];
	if (corners.d[1] < 0) {
		return 2;
	}
	else {
		STUC_ASSERT(
			"",
			corners.d[0] == face.start + corner ||
			corners.d[1] == face.start + corner
		);
		I32 whichCorner = corners.d[0] == face.start + corner;
		I32 otherCorner = corners.d[whichCorner];
		I32 iNext = (corner + 1) % face.size;
		I32 nextBaseCorner = face.start + iNext;
		V2_F32 uv = pMesh->pUvs[nextBaseCorner];
		V2_F32 uvOther = pMesh->pUvs[otherCorner];
		STUC_ASSERT("", v2F32IsFinite(uv) && v2F32IsFinite(uvOther));
		I32 isSeam = !_(uv V2APROXEQL uvOther);
		if (isSeam) {
			return 1;
		}
	}
	return 0;
}

bool stucCheckIfEdgeIsPreserve(const Mesh *pMesh, I32 edge) {
	STUC_ASSERT("", pMesh && edge >= 0);
	if (pMesh->pEdgePreserve) {
		STUC_ASSERT("", pMesh->pEdgePreserve[edge] % 2 == pMesh->pEdgePreserve[edge]);
	}
	return pMesh->pEdgePreserve ? pMesh->pEdgePreserve[edge] : false;
}

bool stucCheckIfVertIsPreserve(const Mesh *pMesh, I32 vert) {
	STUC_ASSERT("", pMesh && vert >= 0);
	if (pMesh->pVertPreserve) {
		STUC_ASSERT("", pMesh->pVertPreserve[vert] % 2 == pMesh->pVertPreserve[vert]);
	}
	return pMesh->pVertPreserve ? pMesh->pVertPreserve[vert] : false;
}

bool stucCheckIfEdgeIsReceive(const Mesh *pMesh, I32 edge, F32 receiveLen) {
	STUC_ASSERT("", pMesh && edge >= 0);
	if (pMesh->pEdgeReceive) {
		STUC_ASSERT("", pMesh->pEdgeReceive[edge] % 2 == pMesh->pEdgeReceive[edge]);
	}
	if (receiveLen >= .0f) {
		STUC_ASSERT("", pMesh->pEdgeLen);
		return pMesh->pEdgeLen[edge] <= receiveLen;
	}
	else if (pMesh->pEdgeReceive) {
		return pMesh->pEdgeReceive[edge];
	}
	return false;
}

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

static
void addTriEdgeToTable(
	const StucAlloc *pAlloc,
	I32 tableSize,
	TriEdge *pEdgeTable,
	I32 verta, I32 vertb,
	I32 tri
) {
	U32 sum = verta + vertb;
	I32 hash = stucFnvHash((U8 *)&sum, 4, tableSize);
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
				pEntry = pEntry->pNext = pAlloc->pCalloc(1, sizeof(TriEdge));
				initTriEdgeEntry(pEntry, verta, vertb, tri);
				break;
			}
			pEntry = pEntry->pNext;
		} while(pEntry);
	}
}

//This gives really long tris, where short tris are possible.
//Re-add search to find short tris, and prefer those.
FaceTriangulated stucTriangulateFace(
	const StucAlloc alloc,
	const FaceRange *pInFace,
	const void *pVerts,
	const I32 *pCorners,
	I32 useStuc
) {
	FaceTriangulated outMesh = {0};
	outMesh.triCount = pInFace->size - 2;
	I32 cornerCount = outMesh.triCount * 3;
	outMesh.pCorners = alloc.pMalloc(sizeof(I32) * cornerCount);
	TriEdge *pEdgeTable = alloc.pCalloc(pInFace->size, sizeof(TriEdge));
	I8 *pVertsRemoved = alloc.pCalloc(pInFace->size, 1);
	I32 cornersLeft = pInFace->size;
	I32 start = 0;
	do {
		I32 ear[3] = {0};
		bool earIsValid = false;
		I32 earFallback[3] = {0};
		bool fallback = false;
		F32 shortestLen = FLT_MAX;
		for (I32 i = 0; i < pInFace->size; ++i) {
			//i %= pInFace->size;
			if (pVertsRemoved[i]) {
				continue;
			}
			I32 ib = getOtherVert(i, pInFace->size, pVertsRemoved);
			I32 ic = getOtherVert(ib, pInFace->size, pVertsRemoved);
			F32 height;
			F32 len;
			if (useStuc) {
				const V2_F32 *pUvs = pVerts;
				V2_F32 verta = pUvs[pInFace->start + i];
				V2_F32 vertb = pUvs[pInFace->start + ib];
				V2_F32 vertc = pUvs[pInFace->start + ic];
				height = v2F32TriHeight(verta, vertb, vertc);
				V2_F32 ac = _(vertc V2SUB verta);
				len = v2F32Len(ac);
			}
			else {
				const V3_F32 *pVertsCast = pVerts;
				V3_F32 verta = pVertsCast[pCorners[pInFace->start + i]];
				V3_F32 vertb = pVertsCast[pCorners[pInFace->start + ib]];
				V3_F32 vertc = pVertsCast[pCorners[pInFace->start + ic]];
				height = v3F32TriHeight(verta, vertb, vertc);
				V3_F32 ac = _(vertc V3SUB verta);
				len = v3F32Len(ac);
			}
			//If ear is not degenerate, then add.
			//Or, if skipped == i, the corner has wrapped back around,
			//without finding a non degenerate ear.
			//In this case, add the ear to avoid an infinite corner
			if (height > .000001f && len < shortestLen) {
				ear[0] = i;
				ear[1] = ib;
				ear[2] = ic;
				earIsValid = true;
				shortestLen = len;
			}
			else if (!fallback) {
				earFallback[0] = i;
				earFallback[1] = ib;
				earFallback[2] = ic;
				fallback = true;
			}
		}
		if (!earIsValid) {
			STUC_ASSERT("", fallback);
			ear[0] = earFallback[0];
			ear[1] = earFallback[1];
			ear[2] = earFallback[2];
		}
		addTriEdgeToTable(
			&alloc,
			pInFace->size,
			pEdgeTable,
			ear[0], ear[1],
			outMesh.cornerCount
		);
		addTriEdgeToTable(
			&alloc,
			pInFace->size,
			pEdgeTable,
			ear[1], ear[2],
			outMesh.cornerCount
		);
		addTriEdgeToTable(
			&alloc,
			pInFace->size,
			pEdgeTable,
			ear[2], ear[0],
			outMesh.cornerCount
		);
		for (I32 i = 0; i < 3; ++i) {
			outMesh.pCorners[outMesh.cornerCount] = ear[i];
			outMesh.cornerCount++;
		}
		start = ear[2];
		pVertsRemoved[ear[1]] = 1;
		cornersLeft--;
	} while (cornersLeft >= 3);
	alloc.pFree(pVertsRemoved);

	//spin

	//free tri edge table
	for (I32 i = 0; i < pInFace->size; ++i) {
		TriEdge* pEntry = pEdgeTable[i].pNext;
		while (pEntry) {
			TriEdge *pNext = pEntry->pNext;
			alloc.pFree(pEntry);
			pEntry = pNext;
		}
	}
	alloc.pFree(pEdgeTable);
	return outMesh;
}

//Caller must check for nan in return value
V3_F32 stucGetBarycentricInFace(
	const V2_F32 *pTriStuc,
	I8 *pTriCorners,
	I32 cornerCount,
	V2_F32 vert
) {
	STUC_ASSERT("", pTriStuc && v2F32IsFinite(*pTriStuc) && v2F32IsFinite(vert));
	STUC_ASSERT("", cornerCount >= 3 && pTriCorners);
	V3_F32 vertBc = cartesianToBarycentric(pTriStuc, &vert);
	if (cornerCount == 4 && v3F32IsFinite(vertBc) && vertBc.d[1] < 0) {
		//base face is a quad, and vert is outside first tri,
		//so use the second tri
		
		//regarding the above condition,
		//because triangulation uses ear clipping,
		//and ngons never hit this block of code,
		//we only need to compare y. As it will always
		//be the point opposite the dividing edge in the quad.
		//This avoids us needing to worry about cases where verts
		//are slightly outside of the quad, by a margin of error.
		//A vert will always end up in one or the other tri.
		V2_F32 triBuf[3] =
			{pTriStuc[2], pTriStuc[3], pTriStuc[0]};
		vertBc = cartesianToBarycentric(triBuf, &vert);
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

typedef struct {
	I32 face;
	I32 corner;
} AdjEntry;

typedef struct {
	AdjEntry* pArr;
	I32 count;
	I32 size;
} AdjBucket;

static
void adjTableDestroyBuckets(const StucAlloc *pAlloc, I32 count, AdjBucket *pAdjTable) {
	for (I32 i = 0; i < count; ++i) {
		if (pAdjTable[i].pArr) {
			pAlloc->pFree(pAdjTable[i].pArr);
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
		FaceRange face = stucGetFaceRange(&pMesh->core, i, false);
		for (I32 j = 0; j < face.size; ++j) {
			AdjBucket* pBucket = pAdjTable + pMesh->core.pCorners[face.start + j];
			STUC_ASSERT("", pBucket->count <= pBucket->size);
			if (!pBucket->pArr) {
				pBucket->size = 2;
				pBucket->pArr =
					pAlloc->pMalloc(sizeof(AdjEntry) * pBucket->size);
			}
			else if (pBucket->count == pBucket->size) {
				STUC_ASSERT("tried to realloc null arr", pBucket->pArr);
				pBucket->size *= 2;
				pBucket->pArr = pAlloc->pRealloc(
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
	FaceRange face = stucGetFaceRange(&pMesh->core, idx, false);
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
			FaceRange otherFace = stucGetFaceRange(&pMesh->core, pEntry->face, false);
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
		pAlloc->pCalloc(pMesh->core.vertCount, sizeof(AdjBucket));
	err = buildCornerAdjTable(pAlloc, pMesh, pAdjTable);
	STUC_THROW_IFNOT(err, "", 0);
	{
		STUC_ASSERT("", pMesh->core.cornerCount);
		I32 dataSize = sizeof(I32) * pMesh->core.cornerCount;
		pMesh->core.pEdges = pAlloc->pMalloc(dataSize);
		memset(pMesh->core.pEdges, -1, dataSize);
		err = findEdges(pMesh, pAdjTable);
		STUC_THROW_IFNOT(err, "'findEdges' returned error", 1);
		STUC_CATCH(1, err,
			pAlloc->pFree(pMesh->core.pEdges);
		);
	}
	STUC_CATCH(0, err, ;);
	adjTableDestroyBuckets(pAlloc, pMesh->core.vertCount, pAdjTable);
	pAlloc->pFree(pAdjTable);
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
	pCtx->stageReport.pBegin(pCtx, &pCtx->stageReport, pName);
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
	pCtx->stageReport.pProgress(pCtx, &pCtx->stageReport, normProgress);
}

void stucStageEndWrap(StucContext pCtx) {
	pCtx->stageReport.pEnd(pCtx, &pCtx->stageReport);
}

void stucSetStageName(StucContext pCtx, const char* pName) {
	strncpy(pCtx->stageReport.stage, pName, STUC_STAGE_NAME_LEN);
}

Mat3x3 stucBuildFaceTbn(FaceRange face, const Mesh *pMesh, const I32 *pCornerOveride) {
	I32 corner = pCornerOveride ? face.start + pCornerOveride[1] : face.start;
	I32 vertIdx = pMesh->core.pCorners[corner];
	V2_F32 uv = pMesh->pUvs[corner];
	V3_F32 vert = pMesh->pVerts[vertIdx];
	I32 next = pCornerOveride ? face.start + pCornerOveride[2] : face.start + 1;
	I32 vertIdxNext = pMesh->core.pCorners[next];
	V2_F32 uvNext = pMesh->pUvs[next];
	V3_F32 vertNext = pMesh->pVerts[vertIdxNext];
	I32 prev = pCornerOveride ? face.start + pCornerOveride[0] : face.end - 1;
	I32 vertIdxPrev = pMesh->core.pCorners[prev];
	V2_F32 uvPrev = pMesh->pUvs[prev];
	V3_F32 vertPrev = pMesh->pVerts[vertIdxPrev];
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
	STUC_ASSERT("", det2 != .0f);
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
		STUC_ASSERT("", det2 != .0f);
		*pt2 = _(ac V2DET *(V2_F32 *)&ab) / det2;
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

void stucSetBorderFaceMapAttrib(
	BorderFace *pEntry,
	UBitField8 *pArr,
	I32 corner,
	I32 value
) {
	I32 len = 3 + pEntry->memType;
	stucSetBitArr(pArr, corner, value, len);
}

void stucInsertionSort(I32 *pIdxTable, I32 count, I32 *pSort) {
	//insertion sort
	I32 a = pSort[0];
	I32 b = pSort[1];
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
			STUC_ASSERT("", j < bufSize && j >= 0);
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
		STUC_ASSERT("", i >= bufSize && i < count);
		bufSize++;
	}
}

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
			STUC_ASSERT("", j < bufSize && j >= 0);
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
		STUC_ASSERT("", i >= bufSize && i < count);
		bufSize++;
	}
}

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

static
bool isMarkedSkip(I32 *pSkip, I32 skipCount, I32 idx) {
	for (I32 i = 0; i < skipCount; ++i) {
		if (idx == pSkip[i]) {
			return true;
		}
	}
	return false;
}

//0 for clockwise, returns 1 for counterclockwise, & 2 if degenerate
I32 stucCalcFaceOrientation(const Mesh *pMesh, const FaceRange *pFace, bool useStuc) {
	STUC_ASSERT("", pFace->start >= 0 && pFace->size >= 3);
	I32 skip[16] = {0};
	I32 skipCount = 0;
	do {
		I32 lowestCorner = 0;
		V2_F32 lowestCoord = { FLT_MAX, FLT_MAX };
		for (I32 i = 0; i < pFace->size; ++i) {
			if (isMarkedSkip(skip, skipCount, i)) {
				continue;
			}
			I32 corner = pFace->start + i;
			V2_F32 pos;
			if (useStuc) {
				pos = pMesh->pUvs[corner];
			}
			else {
				I32 vert = pMesh->core.pCorners[corner];
				pos = *(V2_F32 *)&pMesh->pVerts[vert];
			}
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
		V2_F32 a;
		V2_F32 b;
		V2_F32 c;
		if (useStuc) {
			a = pMesh->pUvs[pFace->start + prev];
			b = pMesh->pUvs[pFace->start + lowestCorner];
			c = pMesh->pUvs[pFace->start + next];
		}
		else {
			I32 vertPrev = pMesh->core.pCorners[pFace->start + prev];
			I32 vert = pMesh->core.pCorners[pFace->start + lowestCorner];
			I32 vertNext = pMesh->core.pCorners[pFace->start + next];
			a = *(V2_F32 *)&pMesh->pVerts[vertPrev];
			b = *(V2_F32 *)&pMesh->pVerts[vert];
			c = *(V2_F32 *)&pMesh->pVerts[vertNext];
		}
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