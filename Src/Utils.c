#include <float.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include <Utils.h>
#include <MathUtils.h>
#include <Context.h>
#include <AttribUtils.h>
#include <Clock.h>

int32_t checkFaceIsInBounds(V2_F32 min, V2_F32 max, FaceRange face, Mesh *pMesh) {
	assert(pMesh && pMesh->pVerts && pMesh->mesh.pLoops);
	assert(face.size >= 3 && face.start >= 0 && face.end >= 0 && face.index >= 0);
	assert(v2IsFinite(min) && v2IsFinite(max));
	V2_F32 faceMin, faceMax;
	faceMin.d[0] = faceMin.d[1] = FLT_MAX;
	faceMax.d[0] = faceMax.d[1] = 0;
	for (int32_t i = 0; i < face.size; ++i) {
		int32_t vertIndex = pMesh->mesh.pLoops[face.start + i];
		V3_F32 *pVert = pMesh->pVerts + vertIndex;
		assert(pVert && v3IsFinite(*pVert));
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
	assert(_(faceMax V2GREATEQL faceMin));
	V2_I32 inside;
	inside.d[0] = (faceMin.d[0] > min.d[0] && faceMin.d[0] < max.d[0]) ||
	           (faceMax.d[0] > min.d[0] && faceMax.d[0] < max.d[0]) ||
			   (faceMin.d[0] < min.d[0] && faceMax.d[0] > max.d[0]);
	inside.d[1] = (faceMin.d[1] > min.d[1] && faceMin.d[1] < max.d[1]) ||
	           (faceMax.d[1] > min.d[1] && faceMax.d[1] < max.d[1]) ||
			   (faceMin.d[1] < min.d[1] && faceMax.d[1] > max.d[1]);
	return inside.d[0] && inside.d[1];
}

uint32_t ruvmFnvHash(uint8_t *value, int32_t valueSize, uint32_t size) {
	assert(value && valueSize > 0 && size > 0);
	uint32_t hash = 2166136261;
	for (int32_t i = 0; i < valueSize; ++i) {
		hash ^= value[i];
		hash *= 16777619;
	}
	hash %= size;
	assert(hash >= 0);
	return hash;
}

void getFaceBounds(FaceBounds *pBounds, V2_F32 *pUvs, FaceRange face) {
	assert(pBounds && pUvs && v2IsFinite(*pUvs));
	assert(face.size >= 3 && face.start >= 0);
	assert(face.end >= 0 && face.index >= 0);
	pBounds->fMin.d[0] = pBounds->fMin.d[1] = FLT_MAX;
	pBounds->fMax.d[0] = pBounds->fMax.d[1] = .0f;
	for (int32_t i = 0; i < face.size; ++i) {
		V2_F32 *uv = pUvs + face.start + i;
		assert(uv && v2IsFinite(*uv));
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
	assert(_(pBounds->fMax V2GREATEQL pBounds->fMin));
}

void getFaceBoundsVert(FaceBounds* pBounds, V3_F32* pVerts, FaceRange face) {
	assert(pBounds && v2IsFinite(pBounds->fMin) && v2IsFinite(pBounds->fMax));
	assert(pVerts && v3IsFinite(*pVerts));
	assert(face.size >= 3 && face.start >= 0);
	assert(face.end >= 0 && face.index >= 0);
	pBounds->fMin.d[0] = pBounds->fMin.d[1] = FLT_MAX;
	pBounds->fMax.d[0] = pBounds->fMax.d[1] = .0f;
	for (int32_t i = 0; i < face.size; ++i) {
		V3_F32* vert = pVerts + face.start + i;
		assert(vert && v3IsFinite(*vert));
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
	assert(_(pBounds->fMax V2GREATEQL pBounds->fMin));
}

int32_t checkIfEdgeIsSeam(int32_t edgeIndex, FaceRange face, int32_t loop,
                          Mesh *pMesh, EdgeVerts *pEdgeVerts) {
	assert(pMesh && pEdgeVerts);
	assert(face.size >= 3 && face.start >= 0 && face.end >= 0 && face.size >= 0);
	assert(edgeIndex >= 0 && loop >= 0 && loop < face.size);
	int32_t *pVerts = pEdgeVerts[edgeIndex].verts;
	assert(pVerts);
	if (pVerts[1] < 0) {
		return 2;
	}
	else {
		assert(pVerts[0] == face.start + loop || pVerts[1] == face.start + loop);
		int32_t whichLoop = pVerts[0] == face.start + loop;
		int32_t otherLoop = pVerts[whichLoop];
		int32_t iNext = (loop + 1) % face.size;
		int32_t nextBaseLoop = face.start + iNext;
		V2_F32 uv = pMesh->pUvs[nextBaseLoop];
		V2_F32 uvOther = pMesh->pUvs[otherLoop];
		assert(v2IsFinite(uv) && v2IsFinite(uvOther));
		int32_t isSeam = _(uv V2NOTEQL uvOther);
		if (isSeam) {
			return 1;
		}
	}
	return 0;
}

int32_t checkIfEdgeIsPreserve(Mesh* pMesh, int32_t edge) {
	assert(pMesh && edge >= 0);
	if (pMesh->pEdgePreserve) {
		assert(pMesh->pEdgePreserve[edge] % 2 == pMesh->pEdgePreserve[edge]);
	}
	return pMesh->pEdgePreserve ? pMesh->pEdgePreserve[edge] : 0;
}

int32_t checkIfEdgeIsReceive(Mesh* pMesh, int32_t edge) {
	assert(pMesh && edge >= 0);
	if (pMesh->pEdgeReceive) {
		assert(pMesh->pEdgeReceive[edge] == 0 || pMesh->pEdgeReceive[edge] == 1);
	}
	return pMesh->pEdgeReceive ? pMesh->pEdgeReceive[edge] : 0;
}

static
int32_t getOtherVert(int32_t i, int32_t faceSize, int8_t *pVertsRemoved) {
	assert(i >= 0 && faceSize >= 3 && i < faceSize && pVertsRemoved);
	int32_t ib = (i + 1) % faceSize;
	//search from i + 1 to facesize, and if non found,
	//then run again from 0 to facesize. If non found then,
	//return error
	int32_t attempts = 0;
	do {
		attempts++;
		for (; ib < faceSize; ++ib) {
			assert(pVertsRemoved[ib] >= 0);
			if (!pVertsRemoved[ib]) {
				return ib;
			}
		}
		ib = 0;
	} while (attempts == 1);
	return -1;
}

FaceTriangulated triangulateFace(RuvmAlloc alloc, FaceRange baseFace, void *pVerts,
                                 int32_t *pLoops, int32_t useUvs) {
	FaceTriangulated outMesh = {0};
	int32_t triCount = baseFace.size - 2;
	outMesh.pTris = alloc.pMalloc(sizeof(int32_t) * triCount);
	int32_t loopCount = triCount * 3;
	outMesh.pLoops = alloc.pMalloc(sizeof(int32_t) * loopCount);
	
	int8_t *pVertsRemoved = alloc.pCalloc(baseFace.size, 1);
	int32_t loopsLeft = baseFace.size;
	do {
		//loop through ears, and find one with shortest edge
		float minDist = FLT_MAX; //min distance
		int32_t nextEar[3] = {0};
		for (int32_t i = 0; i < baseFace.size; ++i) {
			if (pVertsRemoved[i]) {
				continue;
			}
			int32_t ib = getOtherVert(i, baseFace.size, pVertsRemoved);
			int32_t ic = getOtherVert(ib, baseFace.size, pVertsRemoved);
			V2_F32 verta;
			V2_F32 vertb;
			V2_F32 vertc;
			if (useUvs) {
				V2_F32 *pUvs = pVerts;
				verta = pUvs[baseFace.start + i];
				vertb = pUvs[baseFace.start + ib];
				vertc = pUvs[baseFace.start + ic];
			}
			else {
				V3_F32 *pVertsCast = pVerts;
				verta = *(V2_F32 *)(pVertsCast + pLoops[baseFace.start + i]);
				vertb = *(V2_F32 *)(pVertsCast + pLoops[baseFace.start + ib]);
				vertc = *(V2_F32 *)(pVertsCast + pLoops[baseFace.start + ic]);
			}
			if (v2DegenerateTri(verta, vertb, vertc, .00001f)) {
				continue;
			}
			V2_F32 vDist = v2Abs(_(vertc V2SUB verta));
			float dist =
				sqrt(vDist.d[0] * vDist.d[0] + vDist.d[1] * vDist.d[1]);
			if (dist < minDist) {
				minDist = dist;
				nextEar[0] = i;
				nextEar[1] = ib;
				nextEar[2] = ic;
			}
		}
		outMesh.pTris[outMesh.triCount] = outMesh.loopCount;
		for (int32_t i = 0; i < 3; ++i) {
			//set to equal loop index, rather than vert index
			outMesh.pLoops[outMesh.loopCount] = nextEar[i];
			outMesh.loopCount++;
		}
		outMesh.triCount++;
		pVertsRemoved[nextEar[1]] = 1;
		loopsLeft--;
	} while (loopsLeft >= 3);
	alloc.pFree(pVertsRemoved);
	return outMesh;
}

//Caller must check for nan in return value
V3_F32 getBarycentricInFace(V2_F32 *pTriUvs, int8_t *pTriLoops,
                          int32_t loopCount, V2_F32 vert) {
	assert(pTriUvs && v2IsFinite(*pTriUvs) && v2IsFinite(vert));
	assert(loopCount >= 3 && pTriLoops);
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

void waitForJobs(RuvmContext pContext, int32_t *pJobsCompleted, void *pMutex) {
	assert(pContext && pJobsCompleted && *pJobsCompleted >= 0);
	assert(pContext->pThreadPoolHandle && pMutex);
	int32_t waiting;
	do  {
		void (*pJob)(void *) = NULL;
		void *pArgs = NULL;
		pContext->threadPool.pJobStackGetJob(pContext->pThreadPoolHandle,
		                                     &pJob, &pArgs);
		if (pJob) {
			assert(pArgs);
			pJob(pArgs);
		}
		pContext->threadPool.pMutexLock(pContext->pThreadPoolHandle, pMutex);
		assert(pContext->threadCount >= 0);
		waiting = *pJobsCompleted < pContext->threadCount;
		pContext->threadPool.pMutexUnlock(pContext->pThreadPoolHandle, pMutex);
	} while(waiting);
	pContext->threadPool.pMutexDestroy(pContext->pThreadPoolHandle, pMutex);
}

FaceRange getFaceRange(const RuvmMesh *pMesh,
                      const int32_t index, const _Bool border) {
	assert(border % 2 == border);
	int32_t realIndex;
	int32_t direction;
	if (!border) {
		realIndex = index;
		direction = 1;
		assert(index >= 0 && index < pMesh->faceCount);
	}
	else {
		BufMeshIndex bufMeshIndex
			= convertBorderFaceIndex(((BufMesh*)pMesh), index);
		realIndex = bufMeshIndex.realIndex;
		direction = -1;
	}
	FaceRange face = {
		.index = realIndex,
		.start = pMesh->pFaces[realIndex],
		.end = pMesh->pFaces[realIndex + direction]
	};
	if (!border) {
		assert(face.start >= 0 && face.end <= pMesh->loopCount);
		assert(face.start < face.end);
		face.size = face.end - face.start;
	}
	else {
		BufMeshIndex start =
			convertBorderLoopIndex(((BufMesh *)pMesh), face.start);
		BufMeshIndex end =
			convertBorderLoopIndex(((BufMesh *)pMesh), face.end);
		face.start = start.realIndex;
		face.end = end.realIndex;
		assert(face.end >=
		       ((BufMesh *)pMesh)->loopBufSize - 1 -
			   ((BufMesh *)pMesh)->borderLoopCount);
		assert(face.end < face.start);
		assert(face.start < ((BufMesh *)pMesh)->loopBufSize);
		face.size = face.start - face.end;
	}
	assert(face.size >= 3);
	return face;
}

typedef struct {
	int32_t *pBufSize;
	int32_t **ppList[2];
	int32_t *pCount;
	int32_t *pBorderCount;
	AttribArray *pAttribArr;
} BufMeshDomain;

static
void reallocBufMesh(const RuvmAlloc *pAlloc,
                    BufMesh *pMesh, BufMeshDomain *pDomain) {
	int32_t realBorderEnd = *pDomain->pBufSize - 1 - *pDomain->pBorderCount;
	int32_t oldSize = *pDomain->pBufSize;
	*pDomain->pBufSize *= 2;
	int32_t diff = *pDomain->pBufSize - oldSize;
	assert(*pDomain->pBufSize > oldSize);
	for (int32_t i = 0; i < 2; ++i) {
		if (!pDomain->ppList[i]) {
			continue;
		}
		int32_t oldFirstElement = (*pDomain->ppList[i])[realBorderEnd + 1];
		int32_t oldLastElement = (*pDomain->ppList[i])[oldSize - 1];
		*pDomain->ppList[i] =
			pAlloc->pRealloc(*pDomain->ppList[i],
									 sizeof(int32_t) * *pDomain->pBufSize);
		int32_t *pStart = *pDomain->ppList[i] + realBorderEnd + 1;
		memmove(pStart + diff, pStart, sizeof(int32_t) * *pDomain->pBorderCount);
		int32_t newFirstElement = (*pDomain->ppList[i])[realBorderEnd + 1 + diff];
		int32_t newLastElement = (*pDomain->ppList[i])[*pDomain->pBufSize - 1];
		assert(newFirstElement == oldFirstElement);
		assert(newLastElement == oldLastElement);
	}
	reallocAndMoveAttribs(pAlloc, pMesh, pDomain->pAttribArr, realBorderEnd + 1,
	                      diff, *pDomain->pBorderCount, *pDomain->pBufSize);
}

static
BufMeshIndex getNewBufMeshIndex(const RuvmAlloc *pAlloc, BufMesh *pMesh,
                                BufMeshDomain *pDomain, const _Bool border,
								DebugAndPerfVars *pDbVars) {
	//TODO assertions like these need to be converted to release exceptions
	int32_t realBorderEnd = *pDomain->pBufSize - 1 - *pDomain->pBorderCount;
	assert(*pDomain->pCount <= realBorderEnd);
	if (*pDomain->pCount == realBorderEnd) {
		CLOCK_INIT;
		CLOCK_START;
		reallocBufMesh(pAlloc, pMesh, pDomain);
		CLOCK_STOP_NO_PRINT;
		pDbVars->reallocTime += CLOCK_TIME_DIFF(start, stop);
	}
	BufMeshIndex index = {0};
	if (border){
		index.index = *pDomain->pBorderCount;
		index.realIndex = *pDomain->pBufSize - 1 - index.index;
		++*pDomain->pBorderCount;
	}
	else {
		index.index = *pDomain->pCount;
		index.realIndex = index.index;
		++*pDomain->pCount;
	}
	return index;
}

BufMeshIndex bufMeshAddFace(const RuvmAlloc *pAlloc, BufMesh *pMesh,
                            _Bool border, DebugAndPerfVars *pDpVars) {
	BufMeshDomain domain = {
		.pBufSize = &pMesh->faceBufSize,
		.ppList = {&pMesh->mesh.pFaces, NULL},
		.pCount = &pMesh->mesh.faceCount,
		.pBorderCount = &pMesh->borderFaceCount,
		.pAttribArr = &pMesh->mesh.faceAttribs
	};
	BufMeshIndex index = getNewBufMeshIndex(pAlloc, pMesh, &domain, border,
	                                        pDpVars);
	return index;
}

BufMeshIndex bufMeshAddLoop(const RuvmAlloc *pAlloc, BufMesh *pMesh,
                            _Bool border, DebugAndPerfVars *pDpVars) {
	BufMeshDomain domain = {
		.pBufSize = &pMesh->loopBufSize,
		.ppList = {&pMesh->mesh.pLoops, &pMesh->mesh.pEdges},
		.pCount = &pMesh->mesh.loopCount,
		.pBorderCount = &pMesh->borderLoopCount,
		.pAttribArr = &pMesh->mesh.loopAttribs
	};
	BufMeshIndex index = getNewBufMeshIndex(pAlloc, pMesh, &domain, border,
	                                        pDpVars);
	return index;
}

BufMeshIndex bufMeshAddEdge(const RuvmAlloc *pAlloc, BufMesh *pMesh,
                            _Bool border, DebugAndPerfVars *pDpVars) {
	BufMeshDomain domain = {
		.pBufSize = &pMesh->edgeBufSize,
		.ppList = {NULL, NULL},
		.pCount = &pMesh->mesh.edgeCount,
		.pBorderCount = &pMesh->borderEdgeCount,
		.pAttribArr = &pMesh->mesh.edgeAttribs
	};
	BufMeshIndex index = getNewBufMeshIndex(pAlloc, pMesh, &domain, border,
	                                        pDpVars);
	return index;
}

BufMeshIndex bufMeshAddVert(const RuvmAlloc *pAlloc, BufMesh *pMesh,
                            _Bool border, DebugAndPerfVars *pDpVars) {
	BufMeshDomain domain = {
		.pBufSize = &pMesh->vertBufSize,
		.ppList = {NULL, NULL},
		.pCount = &pMesh->mesh.vertCount,
		.pBorderCount = &pMesh->borderVertCount,
		.pAttribArr = &pMesh->mesh.vertAttribs
	};
	BufMeshIndex index = getNewBufMeshIndex(pAlloc, pMesh, &domain, border,
	                                        pDpVars);
	return index;
}

BufMeshIndex convertBorderFaceIndex(const BufMesh *pMesh, int32_t face) {
	assert(face >= 0 && face <= pMesh->borderFaceCount);
	BufMeshIndex index = {
		.index = face,
		.realIndex = pMesh->faceBufSize - 1 - face
	};
	return index;
}
BufMeshIndex convertBorderLoopIndex(const BufMesh *pMesh, int32_t loop) {
	assert(loop >= 0 && loop <= pMesh->borderLoopCount);
	BufMeshIndex index = {
		.index = loop,
		.realIndex = pMesh->loopBufSize - 1 - loop
	};
	return index;
}
BufMeshIndex convertBorderEdgeIndex(const BufMesh *pMesh, int32_t edge) {
	assert(edge >= 0 && edge <= pMesh->borderEdgeCount);
	BufMeshIndex index = {
		.index = edge,
		.realIndex = pMesh->edgeBufSize - 1 - edge
	};
	return index;
}
BufMeshIndex convertBorderVertIndex(const BufMesh *pMesh, int32_t vert) {
	assert(vert >= 0 && vert <= pMesh->borderVertCount);
	BufMeshIndex index = {
		.index = vert,
		.realIndex = pMesh->vertBufSize - 1 - vert
	};
	return index;
}
