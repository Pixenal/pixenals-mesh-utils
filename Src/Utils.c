#include <float.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include <Utils.h>
#include <MathUtils.h>
#include <Context.h>
#include <AttribUtils.h>
#include <Clock.h>
#include <Error.h>

int32_t checkFaceIsInBounds(V2_F32 min, V2_F32 max, FaceRange face, Mesh *pMesh) {
	RUVM_ASSERT("", pMesh && pMesh->pVerts && pMesh->mesh.pLoops);
	RUVM_ASSERT("", face.size >= 3 && face.start >= 0 && face.end >= 0 && face.index >= 0);
	RUVM_ASSERT("", v2IsFinite(min) && v2IsFinite(max));
	V2_F32 faceMin, faceMax;
	faceMin.d[0] = faceMin.d[1] = FLT_MAX;
	faceMax.d[0] = faceMax.d[1] = 0;
	for (int32_t i = 0; i < face.size; ++i) {
		int32_t vertIndex = pMesh->mesh.pLoops[face.start + i];
		V3_F32 *pVert = pMesh->pVerts + vertIndex;
		RUVM_ASSERT("", pVert && v3IsFinite(*pVert));
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
	RUVM_ASSERT("", _(faceMax V2GREATEQL faceMin));
	V2_I32 inside;
	inside.d[0] = (faceMin.d[0] >= min.d[0] && faceMin.d[0] < max.d[0]) ||
	           (faceMax.d[0] >= min.d[0] && faceMax.d[0] < max.d[0]) ||
			   (faceMin.d[0] < min.d[0] && faceMax.d[0] >= max.d[0]);
	inside.d[1] = (faceMin.d[1] >= min.d[1] && faceMin.d[1] < max.d[1]) ||
	           (faceMax.d[1] >= min.d[1] && faceMax.d[1] < max.d[1]) ||
			   (faceMin.d[1] < min.d[1] && faceMax.d[1] >= max.d[1]);
	return inside.d[0] && inside.d[1];
}

uint32_t ruvmFnvHash(uint8_t *value, int32_t valueSize, uint32_t size) {
	RUVM_ASSERT("", value && valueSize > 0 && size > 0);
	uint32_t hash = 2166136261;
	for (int32_t i = 0; i < valueSize; ++i) {
		hash ^= value[i];
		hash *= 16777619;
	}
	hash %= size;
	RUVM_ASSERT("", hash >= 0);
	return hash;
}

void getFaceBounds(FaceBounds *pBounds, V2_F32 *pUvs, FaceRange face) {
	RUVM_ASSERT("", pBounds && pUvs && v2IsFinite(*pUvs));
	RUVM_ASSERT("", face.size >= 3 && face.start >= 0);
	RUVM_ASSERT("", face.end >= 0 && face.index >= 0);
	pBounds->fMin.d[0] = pBounds->fMin.d[1] = FLT_MAX;
	pBounds->fMax.d[0] = pBounds->fMax.d[1] = -FLT_MAX;
	for (int32_t i = 0; i < face.size; ++i) {
		V2_F32 *uv = pUvs + face.start + i;
		RUVM_ASSERT("", uv && v2IsFinite(*uv));
		pBounds->fMin.d[0] = uv->d[0] < pBounds->fMin.d[0] ?
			uv->d[0] : pBounds->fMin.d[0];
		pBounds->fMin.d[1] = uv->d[1] < pBounds->fMin.d[1] ?
			uv->d[1] : pBounds->fMin.d[1];
		pBounds->fMax.d[0] = uv->d[0] > pBounds->fMax.d[0] ?
			uv->d[0] : pBounds->fMax.d[0];
		pBounds->fMax.d[1] = uv->d[1] > pBounds->fMax.d[1] ?
			uv->d[1] : pBounds->fMax.d[1];
	}
	//Faces can be flat (they may be facing sideways in a map for instance)
	RUVM_ASSERT("", _(pBounds->fMax V2GREATEQL pBounds->fMin));
}

void getFaceBoundsVert(FaceBounds* pBounds, V3_F32* pVerts, FaceRange face) {
	RUVM_ASSERT("", pBounds && v2IsFinite(pBounds->fMin) && v2IsFinite(pBounds->fMax));
	RUVM_ASSERT("", pVerts && v3IsFinite(*pVerts));
	RUVM_ASSERT("", face.size >= 3 && face.start >= 0);
	RUVM_ASSERT("", face.end >= 0 && face.index >= 0);
	pBounds->fMin.d[0] = pBounds->fMin.d[1] = FLT_MAX;
	pBounds->fMax.d[0] = pBounds->fMax.d[1] = .0f;
	for (int32_t i = 0; i < face.size; ++i) {
		V3_F32* vert = pVerts + face.start + i;
		RUVM_ASSERT("", vert && v3IsFinite(*vert));
		pBounds->fMin.d[0] = vert->d[0] < pBounds->fMin.d[0] ?
			vert->d[0] : pBounds->fMin.d[0];
		pBounds->fMin.d[1] = vert->d[1] < pBounds->fMin.d[1] ?
			vert->d[1] : pBounds->fMin.d[1];
		pBounds->fMax.d[0] = vert->d[0] > pBounds->fMax.d[0] ?
			vert->d[0] : pBounds->fMax.d[0];
		pBounds->fMax.d[1] = vert->d[1] > pBounds->fMax.d[1] ?
			vert->d[1] : pBounds->fMax.d[1];
	}
	//Faces can be flat (they may be facing sideways in a map for instance)
	RUVM_ASSERT("", _(pBounds->fMax V2GREATEQL pBounds->fMin));
}

int32_t checkIfEdgeIsSeam(int32_t edgeIndex, FaceRange face, int32_t loop,
                          Mesh *pMesh, EdgeVerts *pEdgeVerts) {
	RUVM_ASSERT("", pMesh && pEdgeVerts);
	RUVM_ASSERT("", face.size >= 3 && face.start >= 0 && face.end >= 0 && face.size >= 0);
	RUVM_ASSERT("", edgeIndex >= 0 && loop >= 0 && loop < face.size);
	int32_t *pVerts = pEdgeVerts[edgeIndex].verts;
	RUVM_ASSERT("", pVerts);
	if (pVerts[1] < 0) {
		return 2;
	}
	else {
		RUVM_ASSERT("", pVerts[0] == face.start + loop || pVerts[1] == face.start + loop);
		int32_t whichLoop = pVerts[0] == face.start + loop;
		int32_t otherLoop = pVerts[whichLoop];
		int32_t iNext = (loop + 1) % face.size;
		int32_t nextBaseLoop = face.start + iNext;
		V2_F32 uv = pMesh->pUvs[nextBaseLoop];
		V2_F32 uvOther = pMesh->pUvs[otherLoop];
		RUVM_ASSERT("", v2IsFinite(uv) && v2IsFinite(uvOther));
		int32_t isSeam = !_(uv V2APROXEQL uvOther);
		if (isSeam) {
			return 1;
		}
	}
	return 0;
}

bool checkIfEdgeIsPreserve(Mesh* pMesh, int32_t edge) {
	RUVM_ASSERT("", pMesh && edge >= 0);
	if (pMesh->pEdgePreserve) {
		RUVM_ASSERT("", pMesh->pEdgePreserve[edge] % 2 == pMesh->pEdgePreserve[edge]);
	}
	return pMesh->pEdgePreserve ? pMesh->pEdgePreserve[edge] : false;
}

bool checkIfVertIsPreserve(Mesh* pMesh, int32_t vert) {
	RUVM_ASSERT("", pMesh && vert >= 0);
	if (pMesh->pVertPreserve) {
		RUVM_ASSERT("", pMesh->pVertPreserve[vert] % 2 == pMesh->pVertPreserve[vert]);
	}
	return pMesh->pVertPreserve ? pMesh->pVertPreserve[vert] : false;
}

int32_t checkIfEdgeIsReceive(Mesh* pMesh, int32_t edge) {
	RUVM_ASSERT("", pMesh && edge >= 0);
	if (pMesh->pEdgeReceive) {
		RUVM_ASSERT("", pMesh->pEdgeReceive[edge] == 0 || pMesh->pEdgeReceive[edge] == 1);
	}
	return pMesh->pEdgeReceive ? pMesh->pEdgeReceive[edge] : 0;
}

static
int32_t getOtherVert(int32_t i, int32_t faceSize, int8_t *pVertsRemoved) {
	RUVM_ASSERT("", i >= 0 && faceSize >= 3 && i < faceSize && pVertsRemoved);
	int32_t ib = (i + 1) % faceSize;
	//search from i + 1 to facesize, and if non found,
	//then run again from 0 to facesize. If non found then,
	//return error
	int32_t attempts = 0;
	do {
		attempts++;
		for (; ib < faceSize; ++ib) {
			RUVM_ASSERT("", pVertsRemoved[ib] >= 0);
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
	int32_t tris[2];
	int32_t verts[2];
	bool valid;
} TriEdge;

static
void initTriEdgeEntry(TriEdge* pEntry, int32_t verta, int32_t vertb, int32_t tri) {
	pEntry->tris[0] = tri;
	pEntry->verts[0] = verta;
	pEntry->verts[1] = vertb;
	pEntry->valid = true;
}

static
void addTriEdgeToTable(RuvmAlloc *pAlloc, int32_t tableSize, TriEdge *pEdgeTable, int32_t verta, int32_t vertb, int32_t tri) {
	uint32_t sum = verta + vertb;
	int32_t hash = ruvmFnvHash((uint8_t *)&sum, 4, tableSize);
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
FaceTriangulated triangulateFace(RuvmAlloc alloc, FaceRange baseFace, void *pVerts,
                                 int32_t *pLoops, int32_t useUvs) {
	FaceTriangulated outMesh = {0};
	outMesh.triCount = baseFace.size - 2;
	int32_t loopCount = outMesh.triCount * 3;
	outMesh.pLoops = alloc.pMalloc(sizeof(int32_t) * loopCount);
	TriEdge *pEdgeTable = alloc.pCalloc(baseFace.size, sizeof(TriEdge));
	int8_t *pVertsRemoved = alloc.pCalloc(baseFace.size, 1);
	int32_t loopsLeft = baseFace.size;
	int32_t start = 0;
	int32_t end = baseFace.size;
	do {
		int32_t ear[3] = {0};
		bool earIsValid = false;
		int32_t earFallback[3] = {0};
		bool fallback = false;
		float shortestLen = FLT_MAX;
		for (int32_t i = 0; i < baseFace.size; ++i) {
			//i %= baseFace.size;
			if (pVertsRemoved[i]) {
				continue;
			}
			int32_t ib = getOtherVert(i, baseFace.size, pVertsRemoved);
			int32_t ic = getOtherVert(ib, baseFace.size, pVertsRemoved);
			float height;
			float len;
			if (useUvs) {
				V2_F32 *pUvs = pVerts;
				V2_F32 verta = pUvs[baseFace.start + i];
				V2_F32 vertb = pUvs[baseFace.start + ib];
				V2_F32 vertc = pUvs[baseFace.start + ic];
				height = v2TriHeight(verta, vertb, vertc);
				V2_F32 ac = _(vertc V2SUB verta);
				len = v2Len(ac);
			}
			else {
				V3_F32 *pVertsCast = pVerts;
				V3_F32 verta = pVertsCast[pLoops[baseFace.start + i]];
				V3_F32 vertb = pVertsCast[pLoops[baseFace.start + ib]];
				V3_F32 vertc = pVertsCast[pLoops[baseFace.start + ic]];
				height = v3TriHeight(verta, vertb, vertc);
				V3_F32 ac = _(vertc V3SUB verta);
				len = v3Len(ac);
			}
			//If ear is not degenerate, then add.
			//Or, if skipped == i, the loop has wrapped back around,
			//without finding a non degenerate ear.
			//In this case, add the ear to avoid an infinite loop
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
			RUVM_ASSERT("", fallback);
			ear[0] = earFallback[0];
			ear[1] = earFallback[1];
			ear[2] = earFallback[2];
		}
		addTriEdgeToTable(&alloc, baseFace.size, pEdgeTable, ear[0], ear[1],
		                  outMesh.loopCount);
		addTriEdgeToTable(&alloc, baseFace.size, pEdgeTable, ear[1], ear[2],
		                  outMesh.loopCount);
		addTriEdgeToTable(&alloc, baseFace.size, pEdgeTable, ear[2], ear[0],
		                  outMesh.loopCount);
		for (int32_t i = 0; i < 3; ++i) {
			outMesh.pLoops[outMesh.loopCount] = ear[i];
			outMesh.loopCount++;
		}
		start = ear[2];
		pVertsRemoved[ear[1]] = 1;
		loopsLeft--;
	} while (loopsLeft >= 3);
	alloc.pFree(pVertsRemoved);

	//spin

	//free tri edge table
	for (int32_t i = 0; i < baseFace.size; ++i) {
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
V3_F32 getBarycentricInFace(V2_F32 *pTriUvs, int8_t *pTriLoops,
                          int32_t loopCount, V2_F32 vert) {
	RUVM_ASSERT("", pTriUvs && v2IsFinite(*pTriUvs) && v2IsFinite(vert));
	RUVM_ASSERT("", loopCount >= 3 && pTriLoops);
	V3_F32 vertBc = cartesianToBarycentric(pTriUvs, &vert);
	if (loopCount == 4 && v3IsFinite(vertBc) && vertBc.d[1] < 0) {
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
			{pTriUvs[2], pTriUvs[3], pTriUvs[0]};
		vertBc = cartesianToBarycentric(triBuf, &vert);
		pTriLoops[0] = 2;
		pTriLoops[1] = 3;
	}
	else {
		for (int32_t k = 0; k < 3; ++k) {
			pTriLoops[k] = k;
		}
	}
	return vertBc;
}

//TODO replace custom barrier with system barrier?
void waitForJobs(RuvmContext pContext, int32_t *pJobsCompleted, void *pMutex) {
	RUVM_ASSERT("", pContext && pJobsCompleted && *pJobsCompleted >= 0);
	RUVM_ASSERT("", pContext->pThreadPoolHandle && pMutex);
	int32_t waiting;
	do  {
		void (*pJob)(void *) = NULL;
		void *pArgs = NULL;
		pContext->threadPool.pJobStackGetJob(pContext->pThreadPoolHandle,
		                                     &pJob, &pArgs);
		if (pJob) {
			RUVM_ASSERT("", pArgs);
			pJob(pArgs);
		}
		pContext->threadPool.pMutexLock(pContext->pThreadPoolHandle, pMutex);
		RUVM_ASSERT("", pContext->threadCount >= 0);
		waiting = *pJobsCompleted < pContext->threadCount;
		pContext->threadPool.pMutexUnlock(pContext->pThreadPoolHandle, pMutex);
	} while(waiting);
}

typedef struct {
	int32_t face;
	int32_t loop;
} AdjEntry;

typedef struct {
	AdjEntry* pArr;
	int32_t count;
	int32_t size;
} AdjBucket;

static
void buildLoopAdjTable(RuvmAlloc *pAlloc,
                       Mesh* pMesh, AdjBucket *pAdjTable) {
	for (int32_t i = 0; i < pMesh->mesh.faceCount; ++i) {
		FaceRange face = getFaceRange(&pMesh->mesh, i, false);
		for (int32_t j = 0; j < face.size; ++j) {
			AdjBucket* pBucket = pAdjTable + pMesh->mesh.pLoops[face.start + j];
			RUVM_ASSERT("", pBucket->count <= pBucket->size);
			if (!pBucket->pArr) {
				pBucket->size = 2;
				pBucket->pArr =
					pAlloc->pMalloc(sizeof(AdjEntry) * pBucket->size);
			}
			else if (pBucket->count == pBucket->size) {
				pBucket->size *= 2;
				pBucket->pArr =
					pAlloc->pRealloc(pBucket->pArr, sizeof(AdjEntry) *
						pBucket->size);
			}
			pBucket->pArr[pBucket->count].face = i;
			pBucket->pArr[pBucket->count].loop = j;
			pBucket->count++;
		}
	}
}

static
void findEdges(Mesh* pMesh, AdjBucket* pAdjTable) {
	for (int32_t i = 0; i < pMesh->mesh.faceCount; ++i) {
		FaceRange face = getFaceRange(&pMesh->mesh, i, false);
		for (int32_t j = 0; j < face.size; ++j) {
			if (pMesh->mesh.pEdges[face.start + j] >= 0) {
				continue; //Already set
			}
			int32_t edge = pMesh->mesh.edgeCount;
			pMesh->mesh.edgeCount++;
			AdjBucket* pBucket = pAdjTable + pMesh->mesh.pLoops[face.start + j];
			RUVM_ASSERT("", pBucket->count > 0 &&
				pBucket->size >= pBucket->count);
			for (int32_t k = 0; k < pBucket->count; ++k) {
				AdjEntry* pEntry = pBucket->pArr + k;
				if (pEntry->face == i) {
					RUVM_ASSERT("Invalid mesh, 2 loops in this face share 1 vert",
						pEntry->loop == j);
					continue;
				}
				FaceRange otherFace = getFaceRange(&pMesh->mesh, pEntry->face, false);
				int32_t nextLoop = (j + 1) % face.size;
				int32_t otherPrevLoop = pEntry->loop ?
					pEntry->loop - 1 : otherFace.size - 1;
				if (pMesh->mesh.pEdges[otherFace.start + otherPrevLoop] >= 0) {
					continue; //Already set
				}
				if (pMesh->mesh.pLoops[face.start + nextLoop] !=
					pMesh->mesh.pLoops[otherFace.start + otherPrevLoop]) {
					continue; //Not connected
				}
				pMesh->mesh.pEdges[otherFace.start + otherPrevLoop] = edge;
				break;
			}
			pMesh->mesh.pEdges[face.start + j] = edge;
		}
	}
}

void buildEdgeList(RuvmContext pContext, Mesh* pMesh) {
	RUVM_ASSERT("", !pMesh->mesh.pEdges);
	RuvmAlloc* pAlloc = &pContext->alloc;
	RUVM_ASSERT("", pMesh->mesh.vertCount);
	AdjBucket* pAdjTable =
		pAlloc->pCalloc(pMesh->mesh.vertCount, sizeof(AdjBucket));
	buildLoopAdjTable(pAlloc, pMesh, pAdjTable);

	RUVM_ASSERT("", pMesh->mesh.loopCount);
	int32_t dataSize = sizeof(int32_t) * pMesh->mesh.loopCount;
	pMesh->mesh.pEdges = pAlloc->pMalloc(dataSize);
	memset(pMesh->mesh.pEdges, -1, dataSize);
	findEdges(pMesh, pAdjTable);

	for (int32_t i = 0; i < pMesh->mesh.vertCount; ++i) {
		pAlloc->pFree(pAdjTable[i].pArr);
	}
	pAlloc->pFree(pAdjTable);
}

bool isMeshInvalid(Mesh* pMesh) {
	for (int32_t i = 0; i < pMesh->mesh.faceCount; ++i) {
		FaceRange face = getFaceRange(&pMesh->mesh, i, false);
		if (face.size < 3) {
			return true;
		}
	}
	return false;
}

void progressBarClear() {
	printf("\n");
	printf("\x1b[1F");
	printf("\x1b[2K");
}

void progressBarPrint(RuvmContext pContext, int32_t progress) {
	printf("	");
	for (int32_t i = 0u; i < pContext->stageReport.outOf; ++i) {
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

void stageBegin(void *pContext, RuvmStageReport *pReport, const char* pName) {
	return;
	setStageName(pContext, pName);
}
void stageProgress(void *pContext, RuvmStageReport *pReport, int32_t progress) {
	return;
	if (progress) {
		progressBarClear();
	}
	printf("%s", pReport->stage);
	progressBarPrint(pContext, progress);
}
void stageEnd(void *pContext, RuvmStageReport *pReport) {
	return;
	memset(pReport->stage, 0, RUVM_STAGE_NAME_LEN);
	progressBarClear();
}

void stageBeginWrap(RuvmContext pContext, const char* pName, int32_t max) {
	pContext->stageReport.pBegin(pContext, &pContext->stageReport, pName);
	//Only needed if using default stage report functions,
	//it's just used for the progress bar
	pContext->stageInterval = max <= pContext->stageReport.outOf ?
		1 : max / pContext->stageReport.outOf;
}

void stageProgressWrap(RuvmContext pContext, int32_t progress) {
	if (pContext->stageInterval != 1 && progress % pContext->stageInterval) {
		return;
	}
	//Normalize progress within stageReport.outOf
	int32_t normProgress = progress / pContext->stageInterval;
	pContext->stageReport.pProgress(pContext, &pContext->stageReport, normProgress);
}

void stageEndWrap(RuvmContext pContext) {
	pContext->stageReport.pEnd(pContext, &pContext->stageReport);
}

void setStageName(RuvmContext pContext, const char* pName) {
	strncpy(pContext->stageReport.stage, pName, RUVM_STAGE_NAME_LEN);
}

Mat3x3 buildFaceTbn(FaceRange face, Mesh *pMesh, int32_t *pLoopOveride) {
	int32_t loop = pLoopOveride ? face.start + pLoopOveride[1] : face.start;
	int32_t vertIndex = pMesh->mesh.pLoops[loop];
	V2_F32 uv = pMesh->pUvs[loop];
	V3_F32 vert = pMesh->pVerts[vertIndex];
	int32_t next = pLoopOveride ? face.start + pLoopOveride[2] : face.start + 1;
	int32_t vertIndexNext = pMesh->mesh.pLoops[next];
	V2_F32 uvNext = pMesh->pUvs[next];
	V3_F32 vertNext = pMesh->pVerts[vertIndexNext];
	int32_t prev = pLoopOveride ? face.start + pLoopOveride[0] : face.end - 1;
	int32_t vertIndexPrev = pMesh->mesh.pLoops[prev];
	V2_F32 uvPrev = pMesh->pUvs[prev];
	V3_F32 vertPrev = pMesh->pVerts[vertIndexPrev];
	//uv space direction vectors,
	//forming the coefficient matrix
	Mat2x2 coeffMat;
	*(V2_F32 *)&coeffMat.d[0] = _(uvNext V2SUB uv);
	*(V2_F32 *)&coeffMat.d[1] = _(uvPrev V2SUB uv);
	//object space direction vectors,
	//forming the variable matrix
	Mat2x3 varMat;
	V3_F32 osDirA = _(vertNext V3SUB vert);
	V3_F32 osDirB = _(vertPrev V3SUB vert);
	*(V3_F32 *)&varMat.d[0] = osDirA;
	*(V3_F32 *)&varMat.d[1] = osDirB;
	Mat2x2 coeffMatInv = mat2x2Invert(coeffMat);
	Mat2x3 tb = mat2x2MultiplyMat2x3(coeffMatInv, varMat);
	Mat3x3 tbn;
	*(V3_F32 *)&tbn.d[0] = v3Normalize(*(V3_F32 *)&tb.d[0]);
	*(V3_F32 *)&tbn.d[1] = v3Normalize(*(V3_F32 *)&tb.d[1]);
	V3_F32 normal = _(osDirA V3CROSS osDirB);
	*(V3_F32 *)&tbn.d[2] = v3Normalize(normal);
	return tbn;
}

void getTriScale(int32_t size, BaseTriVerts *pTri) {
	for (int32_t i = 0; i < size; ++i) {
		int32_t iLast = i == 0 ? size - 1 : i - 1;
		int32_t iNext = (i + 1) % size;
		float uvArea = v2TriArea(pTri->uv[iLast], pTri->uv[i], pTri->uv[iNext]);
		float xyzArea = v3TriArea(pTri->xyz[iLast], pTri->xyz[i], pTri->xyz[iNext]);
		pTri->scale[i] = xyzArea / uvArea;
	}
}

//kind of specific to this lib,
//a and b are v3, while and d are v2.
//cd is also taken as a param, while ab is calced
bool calcIntersection(V3_F32 a, V3_F32 b, V2_F32 c, V2_F32 cd,
                      V3_F32 *pPoint, float *pt, float *pt2) {
	V3_F32 ab = _(b V3SUB a);
	V2_F32 ac = _(c V2SUB *(V2_F32 *)&a);
	float det2 = _(*(V2_F32 *)&ab V2DET cd);
	if (det2 == .0f) {
		return false;
	}
	RUVM_ASSERT("", det2 != .0f);
	float t = _(ac V2DET cd) / det2;
	if (pPoint) {
		*pPoint = _(a V3ADD _(ab V3MULS t));
	}
	if (pt) {
		*pt = fabs(t);
	}
	if (pt2) {
		det2 = _(cd V2DET *(V2_F32 *)&ab);
		if (det2 == .0f) {
			return false;
		}
		RUVM_ASSERT("", det2 != .0f);
		*pt2 = fabs(_(ac V2DET *(V2_F32 *)&ab) / det2);
	}
	return true;
}

bool indexBitArray(UBitField8 *pArr, int32_t index) {
	int32_t byte = index / 8;
	int32_t bit = index % 8;
	return pArr[byte] >> bit & 0x1;
}

void setBitArr(UBitField8 *pArr, int32_t index, bool value) {
	int32_t byte = index / 8;
	int32_t bit = index % 8;
	if (value) {
		pArr[byte] |= 0x1 << bit;
	}
	else {
		UBitField8 mask = -0x1 ^ (0x1 << bit);
		pArr[byte] &= mask;
	}
}

void fInsertionSort(int32_t *pIndexTable, int32_t count, float *pSort) {
	//insertion sort
	float a = pSort[0];
	float b = pSort[1];
	int32_t order = a < b;
	pIndexTable[0] = !order;
	pIndexTable[1] = order;
	int32_t bufSize = 2;
	for (int32_t i = bufSize; i < count; ++i) {
		bool insert;
		int32_t j;
		for (j = bufSize - 1; j >= 0; --j) {
			insert = pSort[i] < pSort[pIndexTable[j]] &&
			         pSort[i] > pSort[pIndexTable[j - 1]];
			if (insert) {
				break;
			}
			RUVM_ASSERT("", j < bufSize && j >= 0);
		}
		if (!insert) {
			pIndexTable[bufSize] = i;
		}
		else {
			for (int32_t m = bufSize; m > j; --m) {
				pIndexTable[m] = pIndexTable[m - 1];
				RUVM_ASSERT("", m <= bufSize && m > j);
			}
			pIndexTable[j] = i;
		}
		RUVM_ASSERT("", i >= bufSize && i < count);
		bufSize++;
	}
}

Mat3x3 getInterpolatedTbn(Mesh *pMesh, FaceRange *pFace,
                          int8_t *pTriLoops, V3_F32 bc) {
	//TODO replace interpolation in this func with the attrib
	//     interpolation funcions or macros
	V3_F32 *pNormals = pMesh->pNormals;
	V3_F32 normal = _(pNormals[pFace->start + pTriLoops[0]] V3MULS bc.d[0]);
	_(&normal V3ADDEQL _(pNormals[pFace->start + pTriLoops[1]] V3MULS bc.d[1]));
	_(&normal V3ADDEQL _(pNormals[pFace->start + pTriLoops[2]] V3MULS bc.d[2]));
	_(&normal V3DIVEQLS bc.d[0] + bc.d[1] + bc.d[2]);
	V3_F32 *pTangents = pMesh->pTangents;
	V3_F32 tangent = _(pTangents[pFace->start + pTriLoops[0]] V3MULS bc.d[0]);
	_(&tangent V3ADDEQL _(pTangents[pFace->start + pTriLoops[1]] V3MULS bc.d[1]));
	_(&tangent V3ADDEQL _(pTangents[pFace->start + pTriLoops[2]] V3MULS bc.d[2]));
	_(&tangent V3DIVEQLS bc.d[0] + bc.d[1] + bc.d[2]);
	//TODO should this be interpolated? Or are such edge cases invalid?
	float tSign = pMesh->pTSigns[pFace->start + pTriLoops[0]];
	V3_F32 bitangent = _(_(normal V3CROSS tangent) V3MULS tSign);
	Mat3x3 tbn;
	*(V3_F32 *)&tbn.d[0] = tangent;
	*(V3_F32 *)&tbn.d[1] = bitangent;
	*(V3_F32 *)&tbn.d[2] = normal;
	return tbn;
}

static
bool isMarkedSkip(int32_t *pSkip, int32_t skipCount, int32_t index) {
	for (int32_t i = 0; i < skipCount; ++i) {
		if (index == pSkip[i]) {
			return true;
		}
	}
	return false;
}

//0 for clockwise, returns 1 for counterclockwise, & 2 if degenerate
int32_t calcFaceOrientation(Mesh *pMesh, FaceRange *pFace, bool useUvs) {
	RUVM_ASSERT("", pFace->start >= 0 && pFace->size >= 3);
	int32_t skip[16] = {0};
	int32_t skipCount = 0;
	do {
		int32_t lowestLoop = 0;
		V2_F32 lowestCoord = { FLT_MAX, FLT_MAX };
		for (int32_t i = 0; i < pFace->size; ++i) {
			if (isMarkedSkip(skip, skipCount, i)) {
				continue;
			}
			int32_t loop = pFace->start + i;
			V2_F32 pos;
			if (useUvs) {
				pos = pMesh->pUvs[loop];
			}
			else {
				int32_t vert = pMesh->mesh.pLoops[loop];
				pos = *(V2_F32 *)&pMesh->pVerts[vert];
			}
			if (pos.d[0] < lowestCoord.d[0] ||
				pos.d[1] < lowestCoord.d[1]) {

				lowestLoop = i;
				lowestCoord = pos;
			}
		}
		int32_t prev = lowestLoop == 0 ? pFace->size - 1 : lowestLoop - 1;
		int32_t next = (lowestLoop + 1) % pFace->size;
		V2_F32 a;
		V2_F32 b;
		V2_F32 c;
		if (useUvs) {
			a = pMesh->pUvs[prev];
			b = pMesh->pUvs[lowestLoop];
			c = pMesh->pUvs[next];
		}
		else {
			int32_t vertPrev = pMesh->mesh.pLoops[pFace->start + prev];
			int32_t vert = pMesh->mesh.pLoops[pFace->start + lowestLoop];
			int32_t vertNext = pMesh->mesh.pLoops[pFace->start + next];
			a = *(V2_F32 *)&pMesh->pVerts[vertPrev];
			b = *(V2_F32 *)&pMesh->pVerts[vert];
			c = *(V2_F32 *)&pMesh->pVerts[vertNext];
		}
		//alt formula for determinate,
		//shorter and less likely to cause numerical error
		float det = (b.d[0] - a.d[0]) * (c.d[1] - a.d[1]) -
		            (c.d[0] - a.d[0]) * (b.d[1] - a.d[1]);
		if (det) {
			return det > .0f;
		}
		//abc is degenerate, find another corner
		skip[skipCount] = lowestLoop;
		skipCount++;
	} while(skipCount < pFace->size);
	RUVM_ASSERT("face is degenerate", skipCount == pFace->size);
	return 2;
}