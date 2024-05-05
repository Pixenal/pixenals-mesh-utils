#include <float.h>
#include <math.h>

#include <Utils.h>
#include <MathUtils.h>

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

FaceTriangulated triangulateFace(RuvmAlloc alloc, FaceInfo baseFace, Mesh *pMesh) {
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
		int32_t nextEar[3];
		for (int32_t i = 0; i < baseFace.size; ++i) {
			if (pVertsRemoved[i]) {
				continue;
			}
			int32_t ib = getOtherVert(i, baseFace.size, pVertsRemoved);
			int32_t ic = getOtherVert(ib, baseFace.size, pVertsRemoved);
			V2_F32 uva = pMesh->pUvs[baseFace.start + i];
			V2_F32 uvb = pMesh->pUvs[baseFace.start + ib];
			V2_F32 uvc = pMesh->pUvs[baseFace.start + ic];
			int32_t windingDir = v2WindingCompare(uva, uvb, uvc, 0);
			if (!windingDir) {
				continue;
			}
			V2_F32 vDist = v2Abs(_(uvc V2SUB uva));
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
