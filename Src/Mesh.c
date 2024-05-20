#include <stdint.h>
#include <string.h>
#include <assert.h>

#include <Clock.h>
#include <AttribUtils.h>
#include <Mesh.h>

typedef struct {
	int32_t *pBufSize;
	int32_t **ppList[2]; //2 so both pLoops and pEdges can be realloced at once
	int32_t *pCount;
	AttribArray *pAttribArr;
} MeshDomain;

typedef struct {
	MeshDomain domain;
	int32_t *pBorderCount;
} BufMeshDomain;

static
void reallocBufMesh(const RuvmAlloc *pAlloc,
                    BufMesh *pMesh, BufMeshDomain *pBufDomain) {
	MeshDomain *pDomain = (MeshDomain *)pBufDomain;
	int32_t realBorderEnd = *pDomain->pBufSize - 1 - *pBufDomain->pBorderCount;
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
		memmove(pStart + diff, pStart, sizeof(int32_t) * *pBufDomain->pBorderCount);
		int32_t newFirstElement = (*pDomain->ppList[i])[realBorderEnd + 1 + diff];
		int32_t newLastElement = (*pDomain->ppList[i])[*pDomain->pBufSize - 1];
		assert(newFirstElement == oldFirstElement);
		assert(newLastElement == oldLastElement);
	}
	reallocAndMoveAttribs(pAlloc, pMesh, pDomain->pAttribArr, realBorderEnd + 1,
	                      diff, *pBufDomain->pBorderCount, *pDomain->pBufSize);
}

static
BufMeshIndex getNewBufMeshIndex(const RuvmAlloc *pAlloc, BufMesh *pMesh,
                                BufMeshDomain *pBufDomain, const _Bool border,
								DebugAndPerfVars *pDbVars) {
	MeshDomain *pDomain = (MeshDomain *)pBufDomain;
	int32_t realBorderEnd = *pDomain->pBufSize - 1 - *pBufDomain->pBorderCount;
	//TODO assertions like these need to be converted to release exceptions
	assert(*pDomain->pCount <= realBorderEnd);
	if (*pDomain->pCount == realBorderEnd) {
		CLOCK_INIT;
		CLOCK_START;
		reallocBufMesh(pAlloc, pMesh, pBufDomain);
		CLOCK_STOP_NO_PRINT;
		pDbVars->reallocTime += CLOCK_TIME_DIFF(start, stop);
	}
	BufMeshIndex index = {0};
	if (border){
		index.index = *pBufDomain->pBorderCount;
		index.realIndex = *pDomain->pBufSize - 1 - index.index;
		++*pBufDomain->pBorderCount;
	}
	else {
		index.index = *pDomain->pCount;
		index.realIndex = index.index;
		++*pDomain->pCount;
	}
	return index;
}

static
void getFaceDomain(Mesh *pMesh, MeshDomain *pDomain) {
	pDomain->pBufSize = &pMesh->faceBufSize;
	pDomain->ppList[0] = &pMesh->mesh.pFaces;
	pDomain->pCount = &pMesh->mesh.faceCount;
	pDomain->pAttribArr = &pMesh->mesh.faceAttribs;
}

static
void getLoopDomain(Mesh *pMesh, MeshDomain *pDomain) {
	pDomain->pBufSize = &pMesh->loopBufSize;
	pDomain->ppList[0] = &pMesh->mesh.pLoops;
	pDomain->ppList[1] = &pMesh->mesh.pEdges;
	pDomain->pCount = &pMesh->mesh.loopCount;
	pDomain->pAttribArr = &pMesh->mesh.loopAttribs;
}

static
void getEdgeDomain(Mesh *pMesh, MeshDomain *pDomain) {
	pDomain->pBufSize = &pMesh->edgeBufSize;
	pDomain->ppList[0] = NULL;
	pDomain->pCount = &pMesh->mesh.edgeCount;
	pDomain->pAttribArr = &pMesh->mesh.edgeAttribs;
}

static
void getVertDomain(Mesh *pMesh, MeshDomain *pDomain) {
	pDomain->pBufSize = &pMesh->vertBufSize;
	pDomain->ppList[0] = NULL;
	pDomain->pCount = &pMesh->mesh.vertCount;
	pDomain->pAttribArr = &pMesh->mesh.vertAttribs;
}

static
void reallocMesh(const RuvmAlloc *pAlloc, Mesh *pMesh, MeshDomain *pDomain) {
	*pDomain->pBufSize *= 2;
	for (int32_t i = 0; i < 2; ++i) {
		if (!pDomain->ppList[i]) {
			continue;
		}
		*pDomain->ppList[i] =
			pAlloc->pRealloc(*pDomain->ppList[i],
		                     sizeof(int32_t) * *pDomain->pBufSize);
	}
	reallocAttribs(pAlloc, pMesh, pDomain->pAttribArr, *pDomain->pBufSize);
}

static
int32_t getNewMeshIndex(const RuvmAlloc *pAlloc,
                        Mesh *pMesh, MeshDomain *pDomain) {
	assert(*pDomain->pCount >= 0 && *pDomain->pBufSize > 0);
	assert(*pDomain->pCount <= *pDomain->pBufSize);
	if (*pDomain->pCount == *pDomain->pBufSize) {
		reallocMesh(pAlloc, pMesh, pDomain);
	}
	int32_t index = *pDomain->pCount;
	++*pDomain->pCount;
	return index;
}

BufMeshIndex bufMeshAddFace(const RuvmAlloc *pAlloc, BufMesh *pMesh,
                            _Bool border, DebugAndPerfVars *pDpVars) {
	BufMeshDomain domain = {0};
	getFaceDomain((Mesh *)pMesh, (MeshDomain *)&domain);
	domain.pBorderCount = &pMesh->borderFaceCount;
	BufMeshIndex index = getNewBufMeshIndex(pAlloc, pMesh, &domain, border,
	                                        pDpVars);
	return index;
}

BufMeshIndex bufMeshAddLoop(const RuvmAlloc *pAlloc, BufMesh *pMesh,
                            _Bool border, DebugAndPerfVars *pDpVars) {
	BufMeshDomain domain = {0};
	getLoopDomain((Mesh *)pMesh, (MeshDomain *)&domain);
	domain.pBorderCount = &pMesh->borderLoopCount;
	BufMeshIndex index = getNewBufMeshIndex(pAlloc, pMesh, &domain, border,
	                                        pDpVars);
	return index;
}

BufMeshIndex bufMeshAddEdge(const RuvmAlloc *pAlloc, BufMesh *pMesh,
                            _Bool border, DebugAndPerfVars *pDpVars) {
	BufMeshDomain domain = {0};
	getEdgeDomain((Mesh *)pMesh, (MeshDomain *)&domain);
	domain.pBorderCount = &pMesh->borderEdgeCount;
	BufMeshIndex index = getNewBufMeshIndex(pAlloc, pMesh, &domain, border,
	                                        pDpVars);
	return index;
}

BufMeshIndex bufMeshAddVert(const RuvmAlloc *pAlloc, BufMesh *pMesh,
                            _Bool border, DebugAndPerfVars *pDpVars) {
	BufMeshDomain domain = {0};
	getVertDomain((Mesh *)pMesh, (MeshDomain *)&domain);
	domain.pBorderCount = &pMesh->borderVertCount;
	BufMeshIndex index = getNewBufMeshIndex(pAlloc, pMesh, &domain, border,
	                                        pDpVars);
	return index;
}

BufMeshIndex convertBorderFaceIndex(const BufMesh *pMesh, int32_t face) {
	assert(face >= 0 && face <= pMesh->borderFaceCount);
	BufMeshIndex index = {
		.index = face,
		.realIndex = ((Mesh *)pMesh)->faceBufSize - 1 - face
	};
	return index;
}
BufMeshIndex convertBorderLoopIndex(const BufMesh *pMesh, int32_t loop) {
	assert(loop >= 0 && loop <= pMesh->borderLoopCount);
	BufMeshIndex index = {
		.index = loop,
		.realIndex = ((Mesh *)pMesh)->loopBufSize - 1 - loop
	};
	return index;
}
BufMeshIndex convertBorderEdgeIndex(const BufMesh *pMesh, int32_t edge) {
	assert(edge >= 0 && edge <= pMesh->borderEdgeCount);
	BufMeshIndex index = {
		.index = edge,
		.realIndex = ((Mesh *)pMesh)->edgeBufSize - 1 - edge
	};
	return index;
}
BufMeshIndex convertBorderVertIndex(const BufMesh *pMesh, int32_t vert) {
	assert(vert >= 0 && vert <= pMesh->borderVertCount);
	BufMeshIndex index = {
		.index = vert,
		.realIndex = ((Mesh *)pMesh)->vertBufSize - 1 - vert
	};
	return index;
}

int32_t meshAddFace(const RuvmAlloc *pAlloc, Mesh *pMesh) {
	MeshDomain domain = {0};
	getFaceDomain(pMesh, &domain);
	return getNewMeshIndex(pAlloc, pMesh, &domain);
}

int32_t meshAddLoop(const RuvmAlloc *pAlloc, Mesh *pMesh) {
	MeshDomain domain = {0};
	getLoopDomain(pMesh, &domain);
	return getNewMeshIndex(pAlloc, pMesh, &domain);
}

int32_t meshAddEdge(const RuvmAlloc *pAlloc, Mesh *pMesh) {
	MeshDomain domain = {0};
	getEdgeDomain(pMesh, &domain);
	return getNewMeshIndex(pAlloc, pMesh, &domain);
}

int32_t meshAddVert(const RuvmAlloc *pAlloc, Mesh *pMesh) {
	MeshDomain domain = {0};
	getVertDomain(pMesh, &domain);
	return getNewMeshIndex(pAlloc, pMesh, &domain);
}

void reallocMeshToFit(const RuvmAlloc *pAlloc, Mesh *pMesh) {
	RuvmMesh *pCore = &pMesh->mesh;
	pMesh->faceBufSize = pCore->faceCount + 1;
	int32_t newLen = sizeof(int32_t) * pMesh->faceBufSize;
	pCore->pFaces = pAlloc->pRealloc(pCore->pFaces, newLen);
	pMesh->loopBufSize = pCore->loopCount;
	newLen = sizeof(int32_t) * pMesh->loopBufSize;
	pCore->pLoops = pAlloc->pRealloc(pCore->pLoops, newLen);
	pCore->pEdges = pAlloc->pRealloc(pCore->pEdges, newLen);
	reallocAttribs(pAlloc, pMesh, &pCore->faceAttribs, pMesh->faceBufSize);
	reallocAttribs(pAlloc, pMesh, &pCore->loopAttribs, pMesh->loopBufSize);
	pMesh->edgeBufSize = pCore->edgeCount;
	reallocAttribs(pAlloc, pMesh, &pCore->edgeAttribs, pMesh->edgeBufSize);
	pMesh->vertBufSize = pCore->vertCount;
	reallocAttribs(pAlloc, pMesh, &pCore->vertAttribs, pMesh->vertBufSize);
}

void meshSetLastFace(const RuvmAlloc *pAlloc, Mesh *pMesh) {
	int32_t lastFace = meshAddFace(pAlloc, pMesh);
	pMesh->mesh.pFaces[lastFace] = pMesh->mesh.loopCount;
	//meshAddFace() increments this, so we need to undo that
	pMesh->mesh.faceCount--; 
}

void bufMeshSetLastFaces(const RuvmAlloc *pAlloc, BufMesh *pBufMesh,
                         DebugAndPerfVars *pDpVars) {
	Mesh *pMesh = asMesh(pBufMesh);
	BufMeshIndex lastFace = bufMeshAddFace(pAlloc, pBufMesh, false, pDpVars);
	pMesh->mesh.pFaces[lastFace.realIndex] = pMesh->mesh.loopCount;
	//bufMeshAddFace() increments this, so we need to undo that
	pMesh->mesh.faceCount--; 

	lastFace = bufMeshAddFace(pAlloc, pBufMesh, true, pDpVars);
	pMesh->mesh.pFaces[lastFace.realIndex] = pBufMesh->borderLoopCount;
	pBufMesh->borderFaceCount--;
}
