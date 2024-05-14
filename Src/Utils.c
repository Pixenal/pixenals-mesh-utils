#include <float.h>
#include <math.h>

#include <Utils.h>
#include <MathUtils.h>
#include <Context.h>

int32_t checkFaceIsInBounds(V2_F32 min, V2_F32 max, FaceInfo face, Mesh *pMesh) {
		/*
		int32_t isInside = 0;
		for (int32_t j = 0; j < face.size; ++j) {
			int32_t vertIndex = pMesh->pLoops[face.start + j];
			V2_F32 *vert = (V2_F32 *)(pMesh->pVerts + vertIndex);
			if (_(*vert V2GREATEQL min) && _(*vert V2LESSEQL max)) {
				isInside = 1;
				break;
			}
		}
		return isInside;
		*/
	V2_F32 faceMin, faceMax;
	faceMin.d[0] = faceMin.d[1] = FLT_MAX;
	faceMax.d[0] = faceMax.d[1] = 0;
	for (int32_t i = 0; i < face.size; ++i) {
		int32_t vertIndex = pMesh->mesh.pLoops[face.start + i];
		V3_F32 *pVert = pMesh->pVerts + vertIndex;
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
	uint32_t hash = 2166136261;
	for (int32_t i = 0; i < valueSize; ++i) {
		hash ^= value[i];
		hash *= 16777619;
	}
	hash %= size;
	return hash;
}

void getFaceBounds(FaceBounds *pBounds, V2_F32 *pUvs, FaceInfo faceInfo) {
	pBounds->fMin.d[0] = pBounds->fMin.d[1] = FLT_MAX;
	pBounds->fMax.d[0] = pBounds->fMax.d[1] = .0f;
	for (int32_t i = 0; i < faceInfo.size; ++i) {
		V2_F32 *uv = pUvs + faceInfo.start + i;
		pBounds->fMin.d[0] = uv->d[0] < pBounds->fMin.d[0] ?
			uv->d[0] : pBounds->fMin.d[0];
		pBounds->fMin.d[1] = uv->d[1] < pBounds->fMin.d[1] ?
			uv->d[1] : pBounds->fMin.d[1];
		pBounds->fMax.d[0] = uv->d[0] > pBounds->fMax.d[0] ?
			uv->d[0] : pBounds->fMax.d[0];
		pBounds->fMax.d[1] = uv->d[1] > pBounds->fMax.d[1] ?
			uv->d[1] : pBounds->fMax.d[1];
	}
}

void getFaceBoundsVert(FaceBounds* pBounds, V3_F32* pVerts, FaceInfo faceInfo) {
	pBounds->fMin.d[0] = pBounds->fMin.d[1] = FLT_MAX;
	pBounds->fMax.d[0] = pBounds->fMax.d[1] = .0f;
	for (int32_t i = 0; i < faceInfo.size; ++i) {
		V3_F32* vert = pVerts + faceInfo.start + i;
		pBounds->fMin.d[0] = vert->d[0] < pBounds->fMin.d[0] ?
			vert->d[0] : pBounds->fMin.d[0];
		pBounds->fMin.d[1] = vert->d[1] < pBounds->fMin.d[1] ?
			vert->d[1] : pBounds->fMin.d[1];
		pBounds->fMax.d[0] = vert->d[0] > pBounds->fMax.d[0] ?
			vert->d[0] : pBounds->fMax.d[0];
		pBounds->fMax.d[1] = vert->d[1] > pBounds->fMax.d[1] ?
			vert->d[1] : pBounds->fMax.d[1];
	}
}

int32_t checkIfEdgeIsSeam(int32_t edgeIndex, FaceInfo face, int32_t loop,
                          Mesh *pMesh, EdgeVerts *pEdgeVerts) {
	int32_t *pVerts = pEdgeVerts[edgeIndex].verts;
	if (pVerts[1] < 0) {
		return 2;
	}
	else {
		int32_t whichLoop = pVerts[0] == face.start + loop;
		int32_t otherLoop = pVerts[whichLoop];
		int32_t iNext = (loop + 1) % face.size;
		int32_t nextBaseLoop = face.start + iNext;
		V2_F32 uv = pMesh->pUvs[nextBaseLoop];
		V2_F32 uvOther = pMesh->pUvs[otherLoop];
		int32_t isSeam = _(uv V2NOTEQL uvOther);
		if (isSeam) {
			return 1;
		}
	}
	return 0;
}

int32_t checkIfEdgeIsPreserve(Mesh* pMesh, int32_t edge) {
	return pMesh->pEdgePreserve ? pMesh->pEdgePreserve[edge] : 0;
}

int32_t checkIfEdgeIsReceive(Mesh* pMesh, int32_t edge) {
	return pMesh->pEdgeReceive ? pMesh->pEdgeReceive[edge] : 0;
}

static
int32_t getOtherVert(int32_t i, int32_t faceSize, int8_t *pVertsRemoved) {
	int32_t ib = (i + 1) % faceSize;
	//search from i + 1 to facesize, and if non found,
	//then run again from 0 to facesize. If non found then,
	//return error
	int32_t attempts = 0;
	do {
		attempts++;
		for (; ib < faceSize; ++ib) {
			if (!pVertsRemoved[ib]) {
				return ib;
			}
		}
		ib = 0;
	} while (attempts == 1);
	return -1;
}

FaceTriangulated triangulateFace(RuvmAlloc alloc, FaceInfo baseFace, void *pVerts,
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

V3_F32 getBarycentricInFace(V2_F32 *pTriUvs, int8_t *pTriLoops,
                          int32_t loopCount, V2_F32 vert) {
	V3_F32 vertBc = cartesianToBarycentric(pTriUvs, &vert);
	if (loopCount == 4 && vertBc.d[1] < 0) {
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
	int32_t waiting;
	do  {
		void (*pJob)(void *) = NULL;
		void *pArgs = NULL;
		pContext->threadPool.pJobStackGetJob(pContext->pThreadPoolHandle,
		                                     &pJob, &pArgs);
		if (pJob) {
			pJob(pArgs);
		}
		pContext->threadPool.pMutexLock(pContext->pThreadPoolHandle, pMutex);
		waiting = *pJobsCompleted < pContext->threadCount;
		pContext->threadPool.pMutexUnlock(pContext->pThreadPoolHandle, pMutex);
	} while(waiting);
	pContext->threadPool.pMutexDestroy(pContext->pThreadPoolHandle, pMutex);
}
