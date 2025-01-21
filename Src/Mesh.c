#include <stdint.h>
#include <string.h>

#include <mikktspace.h>

#include <Error.h>
#include <Clock.h>
#include <AttribUtils.h>
#include <Mesh.h>
#include <Context.h>

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

void createMesh(RuvmContext pContext, RuvmObject *pObj, RuvmObjectType type) {
	int32_t size = 0;
	switch (type) {
		case RUVM_OBJECT_DATA_MESH:
			size = sizeof(RuvmMesh);
			break;
		case RUVM_OBJECT_DATA_MESH_INTERN:
			size = sizeof(Mesh);
			break;
		case RUVM_OBJECT_DATA_MESH_BUF:
			size = sizeof(BufMesh);
			break;
		default:
			RUVM_ASSERT("Invalid object data type", false);
			return;
	}
	pObj->pData = pContext->alloc.pCalloc(1, size);
	((RuvmMesh *)pObj->pData)->type.type = type;
}

static
void reallocBufMesh(const RuvmAlloc *pAlloc,
                    BufMesh *pMesh, BufMeshDomain *pBufDomain) {
	MeshDomain *pDomain = (MeshDomain *)pBufDomain;
	int32_t realBorderEnd = *pDomain->pBufSize - 1 - *pBufDomain->pBorderCount;
	int32_t oldSize = *pDomain->pBufSize;
	*pDomain->pBufSize *= 2;
	int32_t diff = *pDomain->pBufSize - oldSize;
	RUVM_ASSERT("", *pDomain->pBufSize > oldSize);
	for (int32_t i = 0; i < 2; ++i) {
		if (!pDomain->ppList[i]) {
			continue;
		}
		int32_t oldFirstElement = (*pDomain->ppList[i])[realBorderEnd + 1];
		int32_t oldLastElement = (*pDomain->ppList[i])[oldSize - 1];
		*pDomain->ppList[i] =
			pAlloc->pRealloc(*pDomain->ppList[i],
									 sizeof(int32_t) * *pDomain->pBufSize);
		if (*pBufDomain->pBorderCount) {
			int32_t *pStart = *pDomain->ppList[i] + realBorderEnd + 1;
			memmove(pStart + diff, pStart, sizeof(int32_t) * *pBufDomain->pBorderCount);
			int32_t newFirstElement = (*pDomain->ppList[i])[realBorderEnd + 1 + diff];
			int32_t newLastElement = (*pDomain->ppList[i])[*pDomain->pBufSize - 1];
			RUVM_ASSERT("", newFirstElement == oldFirstElement);
			RUVM_ASSERT("", newLastElement == oldLastElement);
		}
	}
	reallocAndMoveAttribs(pAlloc, pMesh, pDomain->pAttribArr, realBorderEnd + 1,
	                      diff, *pBufDomain->pBorderCount, *pDomain->pBufSize);
}

static
BufMeshIndex getNewBufMeshIndex(const RuvmAlloc *pAlloc, BufMesh *pMesh,
                                BufMeshDomain *pBufDomain, const bool border,
								DebugAndPerfVars *pDbVars, bool *pRealloced) {
	MeshDomain *pDomain = (MeshDomain *)pBufDomain;
	int32_t realBorderEnd = *pDomain->pBufSize - 1 - *pBufDomain->pBorderCount;
	//TODO assertions like these need to be converted to release exceptions
	RUVM_ASSERT("", *pDomain->pCount <= realBorderEnd);
	if (*pDomain->pCount == realBorderEnd) {
		CLOCK_INIT;
		CLOCK_START;
		reallocBufMesh(pAlloc, pMesh, pBufDomain);
		if (pBufDomain->domain.pAttribArr == &pMesh->mesh.mesh.loopAttribs) {
			pMesh->pW = pMesh->pWAttrib->pData;
			pMesh->pInNormal = pMesh->pInNormalAttrib->pData;
			pMesh->pInTangent = pMesh->pInTangentAttrib->pData;
			pMesh->pAlpha = pMesh->pAlphaAttrib->pData;
			pMesh->pInTSign = pMesh->pInTSignAttrib->pData;
		}
		CLOCK_STOP_NO_PRINT;
		pDbVars->reallocTime += CLOCK_TIME_DIFF(start, stop);
		*pRealloced = true;
	}
	else {
		*pRealloced = false;
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
                        Mesh *pMesh, MeshDomain *pDomain, bool *pRealloced) {
	RUVM_ASSERT("", *pDomain->pCount >= 0 && *pDomain->pBufSize > 0);
	RUVM_ASSERT("", *pDomain->pCount <= *pDomain->pBufSize);
	if (*pDomain->pCount == *pDomain->pBufSize) {
		reallocMesh(pAlloc, pMesh, pDomain);
		*pRealloced = true;
	}
	else {
		*pRealloced = false;
	}
	int32_t index = *pDomain->pCount;
	++*pDomain->pCount;
	return index;
}

BufMeshIndex bufMeshAddFace(const RuvmAlloc *pAlloc, BufMesh *pMesh,
                            bool border, DebugAndPerfVars *pDpVars,
                            bool *pRealloced) {
	BufMeshDomain domain = {0};
	getFaceDomain((Mesh *)pMesh, (MeshDomain *)&domain);
	domain.pBorderCount = &pMesh->borderFaceCount;
	BufMeshIndex index = getNewBufMeshIndex(pAlloc, pMesh, &domain, border,
	                                        pDpVars, pRealloced);
	return index;
}

BufMeshIndex bufMeshAddLoop(const RuvmAlloc *pAlloc, BufMesh *pMesh,
                            bool border, DebugAndPerfVars *pDpVars,
	                        bool *pRealloced) {
	BufMeshDomain domain = {0};
	getLoopDomain((Mesh *)pMesh, (MeshDomain *)&domain);
	domain.pBorderCount = &pMesh->borderLoopCount;
	BufMeshIndex index = getNewBufMeshIndex(pAlloc, pMesh, &domain, border,
	                                        pDpVars, pRealloced);
	return index;
}

BufMeshIndex bufMeshAddEdge(const RuvmAlloc *pAlloc, BufMesh *pMesh,
                            bool border, DebugAndPerfVars *pDpVars,
	                        bool *pRealloced) {
	BufMeshDomain domain = {0};
	getEdgeDomain((Mesh *)pMesh, (MeshDomain *)&domain);
	domain.pBorderCount = &pMesh->borderEdgeCount;
	BufMeshIndex index = getNewBufMeshIndex(pAlloc, pMesh, &domain, border,
	                                        pDpVars, pRealloced);
	return index;
}

BufMeshIndex bufMeshAddVert(const RuvmAlloc *pAlloc, BufMesh *pMesh,
                            bool border, DebugAndPerfVars *pDpVars,
	                        bool *pRealloced) {
	BufMeshDomain domain = {0};
	getVertDomain((Mesh *)pMesh, (MeshDomain *)&domain);
	domain.pBorderCount = &pMesh->borderVertCount;
	BufMeshIndex index = getNewBufMeshIndex(pAlloc, pMesh, &domain, border,
	                                        pDpVars, pRealloced);
	return index;
}

BufMeshIndex convertBorderFaceIndex(const BufMesh *pMesh, int32_t face) {
	RUVM_ASSERT("", face >= 0 && face <= pMesh->borderFaceCount);
	BufMeshIndex index = {
		.index = face,
		.realIndex = ((Mesh *)pMesh)->faceBufSize - 1 - face
	};
	return index;
}
BufMeshIndex convertBorderLoopIndex(const BufMesh *pMesh, int32_t loop) {
	RUVM_ASSERT("", loop >= 0 && loop <= pMesh->borderLoopCount);
	BufMeshIndex index = {
		.index = loop,
		.realIndex = ((Mesh *)pMesh)->loopBufSize - 1 - loop
	};
	return index;
}
BufMeshIndex convertBorderEdgeIndex(const BufMesh *pMesh, int32_t edge) {
	RUVM_ASSERT("", edge >= 0 && edge <= pMesh->borderEdgeCount);
	BufMeshIndex index = {
		.index = edge,
		.realIndex = ((Mesh *)pMesh)->edgeBufSize - 1 - edge
	};
	return index;
}
BufMeshIndex convertBorderVertIndex(const BufMesh *pMesh, int32_t vert) {
	RUVM_ASSERT("", vert >= 0 && vert <= pMesh->borderVertCount);
	BufMeshIndex index = {
		.index = vert,
		.realIndex = ((Mesh *)pMesh)->vertBufSize - 1 - vert
	};
	return index;
}

int32_t meshAddFace(const RuvmAlloc *pAlloc, Mesh *pMesh, bool *pRealloced) {
	MeshDomain domain = {0};
	getFaceDomain(pMesh, &domain);
	return getNewMeshIndex(pAlloc, pMesh, &domain, pRealloced);
}

int32_t meshAddLoop(const RuvmAlloc *pAlloc, Mesh *pMesh, bool *pRealloced) {
	MeshDomain domain = {0};
	getLoopDomain(pMesh, &domain);
	return getNewMeshIndex(pAlloc, pMesh, &domain, pRealloced);
}

int32_t meshAddEdge(const RuvmAlloc *pAlloc, Mesh *pMesh, bool *pRealloced) {
	MeshDomain domain = {0};
	getEdgeDomain(pMesh, &domain);
	return getNewMeshIndex(pAlloc, pMesh, &domain, pRealloced);
}

int32_t meshAddVert(const RuvmAlloc *pAlloc, Mesh *pMesh, bool *pRealloced) {
	MeshDomain domain = {0};
	getVertDomain(pMesh, &domain);
	return getNewMeshIndex(pAlloc, pMesh, &domain, pRealloced);
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
	bool realloced = false;
	int32_t lastFace = meshAddFace(pAlloc, pMesh, &realloced);
	pMesh->mesh.pFaces[lastFace] = pMesh->mesh.loopCount;
	//meshAddFace() increments this, so we need to undo that
	pMesh->mesh.faceCount--; 
}

void bufMeshSetLastFaces(const RuvmAlloc *pAlloc, BufMesh *pBufMesh,
                         DebugAndPerfVars *pDpVars) {
	Mesh *pMesh = asMesh(pBufMesh);
	bool realloced = false;
	BufMeshIndex lastFace = bufMeshAddFace(pAlloc, pBufMesh, false, pDpVars, &realloced);
	pMesh->mesh.pFaces[lastFace.realIndex] = pMesh->mesh.loopCount;
	//bufMeshAddFace() increments this, so we need to undo that
	pMesh->mesh.faceCount--; 

	lastFace = bufMeshAddFace(pAlloc, pBufMesh, true, pDpVars, &realloced);
	pMesh->mesh.pFaces[lastFace.realIndex] = pBufMesh->borderLoopCount;
	pBufMesh->borderFaceCount--;
}

bool checkIfMesh(RuvmMesh *pMesh) {
	switch (pMesh->type.type) {
		case RUVM_OBJECT_DATA_MESH:
			return true;
		case RUVM_OBJECT_DATA_MESH_INTERN:
			return true;
		case RUVM_OBJECT_DATA_MESH_BUF:
			return true;
		default:
			return false;
	}
}

static
void bulkCopyAttribs(AttribArray *pSrc, int32_t SrcOffset,
                     AttribArray *pDest, int32_t dataLen) {
	for (int32_t i = 0; i < pDest->count; ++i) {
		Attrib *pSrcAttrib = getAttrib(pDest->pArr[i].name, pSrc);
		if (!pSrcAttrib) {
			continue;
		}
		void *attribDestStart = attribAsVoid(pDest->pArr + i, SrcOffset);
		int32_t attribTypeSize = getAttribSize(pDest->pArr[i].type);
		memcpy(attribDestStart, pSrcAttrib->pData, attribTypeSize * dataLen);
	}
}

void addToMeshCounts(RuvmContext pContext, MeshCounts *pCounts,
                     MeshCounts *pBoundsCounts, Mesh *pMeshSrc) {
	//TODO maybe replace *Counts vars in Mesh to use MeshCounts,
	//     so we can just do:
	//     meshCountsAdd(totalCount, pBufMesh->mesh.meshCounts);
	//     or something.
	pCounts->faces += pMeshSrc->mesh.faceCount;
	pCounts->loops += pMeshSrc->mesh.loopCount;
	pCounts->edges += pMeshSrc->mesh.edgeCount;
	pCounts->verts += pMeshSrc->mesh.vertCount;
	if (((RuvmObjectType *)pMeshSrc) == RUVM_OBJECT_DATA_MESH_BUF) {
		BufMesh *pBufMesh = (BufMesh *)pMeshSrc;
		pBoundsCounts->faces += pBufMesh->borderFaceCount;
		pBoundsCounts->loops += pBufMesh->borderLoopCount;
		pBoundsCounts->edges += pBufMesh->borderEdgeCount;
		pBoundsCounts->verts += pBufMesh->borderVertCount;
	}
}

void copyMesh(RuvmMesh *pDestMesh, RuvmMesh *pSrcMesh) {
	RUVM_ASSERT("", checkIfMesh(pDestMesh));
	RUVM_ASSERT("", checkIfMesh(pSrcMesh));
	int32_t faceBase = pDestMesh->faceCount;
	int32_t loopBase = pDestMesh->loopCount;
	int32_t edgeBase = pDestMesh->edgeCount;
	int32_t vertBase = pDestMesh->vertCount;
	int32_t *facesStart = pDestMesh->pFaces + faceBase;
	int32_t *loopsStart = pDestMesh->pLoops + loopBase;
	int32_t *edgesStart = pDestMesh->pEdges + loopBase;
	memcpy(facesStart, pSrcMesh->pFaces,
		sizeof(int32_t) * pSrcMesh->faceCount);
	bulkCopyAttribs(&pSrcMesh->faceAttribs, pDestMesh->faceCount,
		&pDestMesh->faceAttribs, pSrcMesh->faceCount);
	pDestMesh->faceCount += pSrcMesh->faceCount;
	pDestMesh->pFaces[pDestMesh->faceCount] = pSrcMesh->pFaces[pSrcMesh->faceCount];
	memcpy(loopsStart, pSrcMesh->pLoops,
		sizeof(int32_t) * pSrcMesh->loopCount);
	bulkCopyAttribs(&pSrcMesh->loopAttribs, pDestMesh->loopCount,
		&pDestMesh->loopAttribs, pSrcMesh->loopCount);
	pDestMesh->loopCount += pSrcMesh->loopCount;
	memcpy(edgesStart, pSrcMesh->pEdges,
		sizeof(int32_t) * pSrcMesh->loopCount);
	bulkCopyAttribs(&pSrcMesh->edgeAttribs, pDestMesh->edgeCount,
		&pDestMesh->edgeAttribs, pSrcMesh->edgeCount);
	pDestMesh->edgeCount += pSrcMesh->edgeCount;
	bulkCopyAttribs(&pSrcMesh->vertAttribs, pDestMesh->vertCount,
		&pDestMesh->vertAttribs, pSrcMesh->vertCount);
	pDestMesh->vertCount += pSrcMesh->vertCount;
	for (int32_t i = faceBase; i < pDestMesh->faceCount; ++i) {
		pDestMesh->pFaces[i] += loopBase;
	}
	pDestMesh->pFaces[pDestMesh->faceCount] += loopBase;
	for (int32_t i = loopBase; i < pDestMesh->loopCount; ++i) {
		pDestMesh->pLoops[i] += vertBase;
		pDestMesh->pEdges[i] += edgeBase;
	}
}

//move these to a separate Obj.c if more functions are made?

void applyObjTransform(RuvmObject *pObj) {
	Mesh *pMesh = pObj->pData;
	for (int32_t i = 0; i < pMesh->mesh.vertCount; ++i) {
		V3_F32 *pV3 = pMesh->pVerts + i;
		V4_F32 v4 = {pV3->d[0], pV3->d[1], pV3->d[2], 1.0f};
		_(&v4 V4MULEQLM4X4 &pObj->transform);
		*pV3 = *(V3_F32 *)&v4;
	}
	pObj->transform = identM4x4;
}

void mergeObjArr(RuvmContext pContext, Mesh *pMesh,
                 int32_t objCount, RuvmObject *pObjArr, bool setCommon) {
	//TODO allocate map mesh based on all meshes in obj arr
	Mesh **ppSrcs = pContext->alloc.pCalloc(objCount, sizeof(void *));
	MeshCounts totalCount = {0};
	for (int32_t i = 0; i < objCount; ++i) {
		ppSrcs[i] = pObjArr[i].pData;
		addToMeshCounts(pContext, &totalCount, NULL, (Mesh *)pObjArr[i].pData);
	}
	pMesh->faceBufSize = totalCount.faces + 1; //+1 for last face index
	pMesh->loopBufSize = totalCount.loops;
	pMesh->edgeBufSize = totalCount.edges;
	pMesh->vertBufSize = totalCount.verts;
	pMesh->mesh.pFaces =
		pContext->alloc.pMalloc(sizeof(int32_t) * pMesh->faceBufSize);
	pMesh->mesh.pLoops =
		pContext->alloc.pMalloc(sizeof(int32_t) * pMesh->loopBufSize);
	pMesh->mesh.pEdges =
		pContext->alloc.pMalloc(sizeof(int32_t) * pMesh->loopBufSize);
	allocAttribsFromMeshArr(&pContext->alloc, pMesh, objCount, ppSrcs, setCommon);
	for (int32_t i = 0; i < objCount; ++i) {
		copyMesh(pMesh, pObjArr[i].pData);
	}
}

void destroyObjArr(RuvmContext pContext, int32_t objCount, RuvmObject *pObjArr) {
	for (int32_t i = 0; i < objCount; ++i) {
		ruvmMeshDestroy(pContext, pObjArr[i].pData);
		pContext->alloc.pFree(pObjArr[i].pData);
	}
	pContext->alloc.pFree(pObjArr);
}

FaceRange getFaceRange(const RuvmMesh *pMesh,
                      const int32_t index, const bool border) {
	RUVM_ASSERT("", border % 2 == border);
	int32_t realIndex;
	int32_t direction;
	if (!border) {
		realIndex = index;
		direction = 1;
		RUVM_ASSERT("", index >= 0 && index < pMesh->faceCount);
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
		RUVM_ASSERT("", face.start >= 0 && face.end <= pMesh->loopCount);
		RUVM_ASSERT("", face.start < face.end);
		face.size = face.end - face.start;
	}
	else {
		BufMeshIndex start =
			convertBorderLoopIndex(((BufMesh *)pMesh), face.start);
		BufMeshIndex end =
			convertBorderLoopIndex(((BufMesh *)pMesh), face.end);
		face.start = start.realIndex;
		face.end = end.realIndex;
		RUVM_ASSERT("", face.end >=
		       ((Mesh *)pMesh)->loopBufSize - 1 -
			   ((BufMesh *)pMesh)->borderLoopCount);
		RUVM_ASSERT("", face.end < face.start);
		RUVM_ASSERT("", face.start < ((Mesh *)pMesh)->loopBufSize);
		face.size = face.start - face.end;
	}
	RUVM_ASSERT("", face.size >= 3);
	return face;
}

static
int mikktGetNumFaces(const SMikkTSpaceContext *pContext) {
	Mesh *pMesh = pContext->m_pUserData;
	return pMesh->mesh.faceCount;
}

static
int mikktGetNumVertsOfFace(const SMikkTSpaceContext *pContext, const int iFace) {
	Mesh *pMesh = pContext->m_pUserData;
	FaceRange face = getFaceRange(&pMesh->mesh, iFace, false);
	return face.size;
}

static
void mikktGetPos(const SMikkTSpaceContext *pContext, float *pFvPosOut,
                 const int iFace, const int iVert) {
	Mesh *pMesh = pContext->m_pUserData;
	FaceRange face = getFaceRange(&pMesh->mesh, iFace, false);
	int32_t vertIndex = pMesh->mesh.pLoops[face.start + iVert];
	*(V3_F32 *)pFvPosOut = pMesh->pVerts[vertIndex];
}

static
void mikktGetNormal(const SMikkTSpaceContext *pContext, float *pFvNormOut,
                    const int iFace, const int iVert) {
	Mesh *pMesh = pContext->m_pUserData;
	FaceRange face = getFaceRange(&pMesh->mesh, iFace, false);
	*(V3_F32 *)pFvNormOut = pMesh->pNormals[face.start + iVert];;
}

static
void mikktGetTexCoord(const SMikkTSpaceContext *pContext, float *pFvTexcOut,
                      const int iFace, const int iVert) {
	Mesh *pMesh = pContext->m_pUserData;
	FaceRange face = getFaceRange(&pMesh->mesh, iFace, false);
	*(V2_F32 *)pFvTexcOut = pMesh->pUvs[face.start + iVert];
}

static
void mikktSetTSpaceBasic(const SMikkTSpaceContext *pContext, const float *pFvTangent,
                         const float fSign, const int iFace, const int iVert) {
	Mesh *pMesh = pContext->m_pUserData;
	FaceRange face = getFaceRange(&pMesh->mesh, iFace, false);
	int32_t loop = face.start + iVert;
	pMesh->pTangents[loop] = *(V3_F32 *)pFvTangent;
	pMesh->pTSigns[loop] = fSign;
}

void buildTangents(Mesh *pMesh) {
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
		.m_pUserData = pMesh
	};
	genTangSpaceDefault(&mikktContext);
}