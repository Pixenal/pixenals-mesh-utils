#include <stdint.h>
#include <string.h>

#include <mikktspace.h>

#include <Error.h>
#include <Clock.h>
#include <AttribUtils.h>
#include <Mesh.h>
#include <Context.h>

typedef struct {
	I32 *pBufSize;
	I32 **ppList[2]; //2 so both pCorners and pEdges can be realloced at once
	I32 *pCount;
	AttribArray *pAttribArr;
} MeshDomain;

typedef struct {
	MeshDomain domain;
	I32 *pBorderCount;
} BufMeshDomain;

void stucCreateMesh(StucContext pContext, StucObject *pObj, StucObjectType type) {
	I32 size = 0;
	switch (type) {
		case STUC_OBJECT_DATA_MESH:
			size = sizeof(StucMesh);
			break;
		case STUC_OBJECT_DATA_MESH_INTERN:
			size = sizeof(Mesh);
			break;
		case STUC_OBJECT_DATA_MESH_BUF:
			size = sizeof(BufMesh);
			break;
		default:
			STUC_ASSERT("Invalid object data type", false);
			return;
	}
	pObj->pData = pContext->alloc.pCalloc(1, size);
	((StucMesh *)pObj->pData)->type.type = type;
}

static
void reallocBufMesh(const StucAlloc *pAlloc,
                    BufMesh *pMesh, BufMeshDomain *pBufDomain) {
	MeshDomain *pDomain = (MeshDomain *)pBufDomain;
	I32 realBorderEnd = *pDomain->pBufSize - 1 - *pBufDomain->pBorderCount;
	I32 oldSize = *pDomain->pBufSize;
	*pDomain->pBufSize *= 2;
	I32 diff = *pDomain->pBufSize - oldSize;
	STUC_ASSERT("", *pDomain->pBufSize > oldSize);
	for (I32 i = 0; i < 2; ++i) {
		if (!pDomain->ppList[i]) {
			continue;
		}
		I32 oldFirstElement = (*pDomain->ppList[i])[realBorderEnd + 1];
		I32 oldLastElement = (*pDomain->ppList[i])[oldSize - 1];
		*pDomain->ppList[i] =
			pAlloc->pRealloc(*pDomain->ppList[i], sizeof(I32) * *pDomain->pBufSize);
		if (*pBufDomain->pBorderCount) {
			I32 *pStart = *pDomain->ppList[i] + realBorderEnd + 1;
			memmove(pStart + diff, pStart, sizeof(I32) * *pBufDomain->pBorderCount);
			I32 newFirstElement = (*pDomain->ppList[i])[realBorderEnd + 1 + diff];
			I32 newLastElement = (*pDomain->ppList[i])[*pDomain->pBufSize - 1];
			STUC_ASSERT("", newFirstElement == oldFirstElement);
			STUC_ASSERT("", newLastElement == oldLastElement);
		}
	}
	stucReallocAndMoveAttribs(pAlloc, pMesh, pDomain->pAttribArr, realBorderEnd + 1,
	                          diff, *pBufDomain->pBorderCount, *pDomain->pBufSize);
}

static
BufMeshIdx getNewBufMeshIdx(const StucAlloc *pAlloc, BufMesh *pMesh,
                            BufMeshDomain *pBufDomain, const bool border,
                            DebugAndPerfVars *pDbVars, bool *pRealloced) {
	MeshDomain *pDomain = (MeshDomain *)pBufDomain;
	I32 realBorderEnd = *pDomain->pBufSize - 1 - *pBufDomain->pBorderCount;
	//TODO assertions like these need to be converted to release exceptions
	STUC_ASSERT("", *pDomain->pCount <= realBorderEnd);
	if (*pDomain->pCount == realBorderEnd) {
		CLOCK_INIT;
		CLOCK_START;
		reallocBufMesh(pAlloc, pMesh, pBufDomain);
		CLOCK_STOP_NO_PRINT;
		pDbVars->reallocTime += CLOCK_TIME_DIFF(start, stop);
		*pRealloced = true;
	}
	else {
		*pRealloced = false;
	}
	BufMeshIdx idx = {0};
	if (border){
		idx.idx = *pBufDomain->pBorderCount;
		idx.realIdx = *pDomain->pBufSize - 1 - idx.idx;
		++*pBufDomain->pBorderCount;
	}
	else {
		idx.idx = *pDomain->pCount;
		idx.realIdx = idx.idx;
		++*pDomain->pCount;
	}
	return idx;
}

static
void getFaceDomain(Mesh *pMesh, MeshDomain *pDomain) {
	pDomain->pBufSize = &pMesh->faceBufSize;
	pDomain->ppList[0] = &pMesh->core.pFaces;
	pDomain->pCount = &pMesh->core.faceCount;
	pDomain->pAttribArr = &pMesh->core.faceAttribs;
}

static
void getCornerDomain(Mesh *pMesh, MeshDomain *pDomain) {
	pDomain->pBufSize = &pMesh->cornerBufSize;
	pDomain->ppList[0] = &pMesh->core.pCorners;
	pDomain->ppList[1] = &pMesh->core.pEdges;
	pDomain->pCount = &pMesh->core.cornerCount;
	pDomain->pAttribArr = &pMesh->core.cornerAttribs;
}

static
void getEdgeDomain(Mesh *pMesh, MeshDomain *pDomain) {
	pDomain->pBufSize = &pMesh->edgeBufSize;
	pDomain->ppList[0] = NULL;
	pDomain->pCount = &pMesh->core.edgeCount;
	pDomain->pAttribArr = &pMesh->core.edgeAttribs;
}

static
void getVertDomain(Mesh *pMesh, MeshDomain *pDomain) {
	pDomain->pBufSize = &pMesh->vertBufSize;
	pDomain->ppList[0] = NULL;
	pDomain->pCount = &pMesh->core.vertCount;
	pDomain->pAttribArr = &pMesh->core.vertAttribs;
}

static
void reallocMesh(const StucAlloc *pAlloc, Mesh *pMesh, MeshDomain *pDomain) {
	*pDomain->pBufSize *= 2;
	for (I32 i = 0; i < 2; ++i) {
		if (!pDomain->ppList[i]) {
			continue;
		}
		*pDomain->ppList[i] =
			pAlloc->pRealloc(*pDomain->ppList[i], sizeof(I32) * *pDomain->pBufSize);
	}
	stucReallocAttribArr(pAlloc, pMesh, pDomain->pAttribArr, *pDomain->pBufSize);
}

static
I32 getNewMeshIdx(const StucAlloc *pAlloc,
                      Mesh *pMesh, MeshDomain *pDomain, bool *pRealloced) {
	STUC_ASSERT("", *pDomain->pCount >= 0 && *pDomain->pBufSize > 0);
	STUC_ASSERT("", *pDomain->pCount <= *pDomain->pBufSize);
	if (*pDomain->pCount == *pDomain->pBufSize) {
		reallocMesh(pAlloc, pMesh, pDomain);
		*pRealloced = true;
	}
	else {
		*pRealloced = false;
	}
	I32 idx = *pDomain->pCount;
	++*pDomain->pCount;
	return idx;
}

BufMeshIdx stucBufMeshAddFace(const StucAlloc *pAlloc, BufMesh *pMesh,
                              bool border, DebugAndPerfVars *pDpVars,
                              bool *pRealloced) {
	BufMeshDomain domain = {0};
	getFaceDomain((Mesh *)pMesh, (MeshDomain *)&domain);
	domain.pBorderCount = &pMesh->borderFaceCount;
	BufMeshIdx idx = getNewBufMeshIdx(pAlloc, pMesh, &domain, border,
	                                  pDpVars, pRealloced);
	return idx;
}

BufMeshIdx stucBufMeshAddCorner(const StucAlloc *pAlloc, BufMesh *pMesh,
                                bool border, DebugAndPerfVars *pDpVars,
                                bool *pRealloced) {
	BufMeshDomain domain = {0};
	getCornerDomain((Mesh *)pMesh, (MeshDomain *)&domain);
	domain.pBorderCount = &pMesh->borderCornerCount;
	BufMeshIdx idx = getNewBufMeshIdx(pAlloc, pMesh, &domain, border,
	                                  pDpVars, pRealloced);
	return idx;
}

BufMeshIdx stucBufMeshAddEdge(const StucAlloc *pAlloc, BufMesh *pMesh,
                              bool border, DebugAndPerfVars *pDpVars,
                              bool *pRealloced) {
	BufMeshDomain domain = {0};
	getEdgeDomain((Mesh *)pMesh, (MeshDomain *)&domain);
	domain.pBorderCount = &pMesh->borderEdgeCount;
	BufMeshIdx idx = getNewBufMeshIdx(pAlloc, pMesh, &domain, border,
	                                  pDpVars, pRealloced);
	return idx;
}

BufMeshIdx stucBufMeshAddVert(const StucAlloc *pAlloc, BufMesh *pMesh,
                              bool border, DebugAndPerfVars *pDpVars,
                              bool *pRealloced) {
	BufMeshDomain domain = {0};
	getVertDomain((Mesh *)pMesh, (MeshDomain *)&domain);
	domain.pBorderCount = &pMesh->borderVertCount;
	BufMeshIdx idx = getNewBufMeshIdx(pAlloc, pMesh, &domain, border,
	                                  pDpVars, pRealloced);
	return idx;
}

BufMeshIdx stucConvertBorderFaceIdx(const BufMesh *pMesh, I32 face) {
	STUC_ASSERT("", face >= 0 && face <= pMesh->borderFaceCount);
	BufMeshIdx idx = {
		.idx = face,
		.realIdx = ((Mesh *)pMesh)->faceBufSize - 1 - face
	};
	return idx;
}
BufMeshIdx stucConvertBorderCornerIdx(const BufMesh *pMesh, I32 corner) {
	STUC_ASSERT("", corner >= 0 && corner <= pMesh->borderCornerCount);
	BufMeshIdx idx = {
		.idx = corner,
		.realIdx = ((Mesh *)pMesh)->cornerBufSize - 1 - corner
	};
	return idx;
}
BufMeshIdx stucConvertBorderEdgeIdx(const BufMesh *pMesh, I32 edge) {
	STUC_ASSERT("", edge >= 0 && edge <= pMesh->borderEdgeCount);
	BufMeshIdx idx = {
		.idx = edge,
		.realIdx = ((Mesh *)pMesh)->edgeBufSize - 1 - edge
	};
	return idx;
}
BufMeshIdx stucConvertBorderVertIdx(const BufMesh *pMesh, I32 vert) {
	STUC_ASSERT("", vert >= 0 && vert <= pMesh->borderVertCount);
	BufMeshIdx idx = {
		.idx = vert,
		.realIdx = ((Mesh *)pMesh)->vertBufSize - 1 - vert
	};
	return idx;
}

I32 stucMeshAddFace(const StucAlloc *pAlloc, Mesh *pMesh, bool *pRealloced) {
	MeshDomain domain = {0};
	getFaceDomain(pMesh, &domain);
	return getNewMeshIdx(pAlloc, pMesh, &domain, pRealloced);
}

I32 stucMeshAddCorner(const StucAlloc *pAlloc, Mesh *pMesh, bool *pRealloced) {
	MeshDomain domain = {0};
	getCornerDomain(pMesh, &domain);
	return getNewMeshIdx(pAlloc, pMesh, &domain, pRealloced);
}

I32 stucMeshAddEdge(const StucAlloc *pAlloc, Mesh *pMesh, bool *pRealloced) {
	MeshDomain domain = {0};
	getEdgeDomain(pMesh, &domain);
	return getNewMeshIdx(pAlloc, pMesh, &domain, pRealloced);
}

I32 stucMeshAddVert(const StucAlloc *pAlloc, Mesh *pMesh, bool *pRealloced) {
	MeshDomain domain = {0};
	getVertDomain(pMesh, &domain);
	return getNewMeshIdx(pAlloc, pMesh, &domain, pRealloced);
}

void stucReallocMeshToFit(const StucAlloc *pAlloc, Mesh *pMesh) {
	StucMesh *pCore = &pMesh->core;
	pMesh->faceBufSize = pCore->faceCount + 1;
	I32 newLen = sizeof(I32) * pMesh->faceBufSize;
	pCore->pFaces = pAlloc->pRealloc(pCore->pFaces, newLen);
	pMesh->cornerBufSize = pCore->cornerCount;
	newLen = sizeof(I32) * pMesh->cornerBufSize;
	pCore->pCorners = pAlloc->pRealloc(pCore->pCorners, newLen);
	pCore->pEdges = pAlloc->pRealloc(pCore->pEdges, newLen);
	stucReallocAttribArr(pAlloc, pMesh, &pCore->faceAttribs, pMesh->faceBufSize);
	stucReallocAttribArr(pAlloc, pMesh, &pCore->cornerAttribs, pMesh->cornerBufSize);
	pMesh->edgeBufSize = pCore->edgeCount;
	stucReallocAttribArr(pAlloc, pMesh, &pCore->edgeAttribs, pMesh->edgeBufSize);
	pMesh->vertBufSize = pCore->vertCount;
	stucReallocAttribArr(pAlloc, pMesh, &pCore->vertAttribs, pMesh->vertBufSize);
}

void stucMeshSetLastFace(const StucAlloc *pAlloc, Mesh *pMesh) {
	bool realloced = false;
	I32 lastFace = stucMeshAddFace(pAlloc, pMesh, &realloced);
	pMesh->core.pFaces[lastFace] = pMesh->core.cornerCount;
	//meshAddFace() increments this, so we need to undo that
	pMesh->core.faceCount--; 
}

void stucBufMeshSetLastFaces(const StucAlloc *pAlloc, BufMesh *pBufMesh,
                             DebugAndPerfVars *pDpVars) {
	Mesh *pMesh = &pBufMesh->mesh;
	bool realloced = false;
	BufMeshIdx lastFace = stucBufMeshAddFace(pAlloc, pBufMesh, false, pDpVars, &realloced);
	pMesh->core.pFaces[lastFace.realIdx] = pMesh->core.cornerCount;
	//bufMeshAddFace() increments this, so we need to undo that
	pMesh->core.faceCount--; 

	lastFace = stucBufMeshAddFace(pAlloc, pBufMesh, true, pDpVars, &realloced);
	pMesh->core.pFaces[lastFace.realIdx] = pBufMesh->borderCornerCount;
	pBufMesh->borderFaceCount--;
}

bool stucCheckIfMesh(StucObjectData type) {
	switch (type.type) {
		case STUC_OBJECT_DATA_MESH:
			return true;
		case STUC_OBJECT_DATA_MESH_INTERN:
			return true;
		case STUC_OBJECT_DATA_MESH_BUF:
			return true;
		default:
			return false;
	}
}

static
void bulkCopyAttribs(AttribArray *pSrc, I32 SrcOffset,
                     AttribArray *pDest, I32 dataLen) {
	for (I32 i = 0; i < pDest->count; ++i) {
		Attrib *pSrcAttrib = stucGetAttribIntern(pDest->pArr[i].core.name, pSrc);
		if (!pSrcAttrib) {
			continue;
		}
		void *attribDestStart = stucAttribAsVoid(&pDest->pArr[i].core, SrcOffset);
		I32 attribTypeSize = stucGetAttribSizeIntern(pDest->pArr[i].core.type);
		memcpy(attribDestStart, pSrcAttrib->core.pData, attribTypeSize * dataLen);
	}
}

void stucAddToMeshCounts(StucContext pContext, MeshCounts *pCounts,
                         MeshCounts *pBoundsCounts, Mesh *pMeshSrc) {
	//TODO maybe replace *Counts vars in Mesh to use MeshCounts,
	//     so we can just do:
	//     meshCountsAdd(totalCount, pBufMesh->mesh.meshCounts);
	//     or something.
	pCounts->faces += pMeshSrc->core.faceCount;
	pCounts->corners += pMeshSrc->core.cornerCount;
	pCounts->edges += pMeshSrc->core.edgeCount;
	pCounts->verts += pMeshSrc->core.vertCount;
	if (pMeshSrc->core.type.type == STUC_OBJECT_DATA_MESH_BUF) {
		BufMesh *pBufMesh = (BufMesh *)pMeshSrc;
		pBoundsCounts->faces += pBufMesh->borderFaceCount;
		pBoundsCounts->corners += pBufMesh->borderCornerCount;
		pBoundsCounts->edges += pBufMesh->borderEdgeCount;
		pBoundsCounts->verts += pBufMesh->borderVertCount;
	}
}

void stucCopyMesh(StucMesh *pDestMesh, StucMesh *pSrcMesh) {
	printf("pSrcMesh->type.type			%d\n", pSrcMesh->type.type);
	printf("pSrcMesh->faceCount			%d\n", pSrcMesh->faceCount);
	printf("pSrcMesh->cornerCount		%d\n", pSrcMesh->cornerCount);
	printf("pSrcMesh->vertCount			%d\n", pSrcMesh->vertCount);
	printf("pSrcMesh->pFaces			%p\n", pSrcMesh->pFaces);
	if (pSrcMesh->type.type == STUC_OBJECT_DATA_NULL) {
		return;
	}
	STUC_ASSERT("", stucCheckIfMesh(pDestMesh->type));
	STUC_ASSERT("", stucCheckIfMesh(pSrcMesh->type));
	I32 faceBase = pDestMesh->faceCount;
	I32 cornerBase = pDestMesh->cornerCount;
	I32 edgeBase = pDestMesh->edgeCount;
	I32 vertBase = pDestMesh->vertCount;
	I32 *facesStart = pDestMesh->pFaces + faceBase;
	I32 *cornersStart = pDestMesh->pCorners + cornerBase;
	I32 *edgesStart = pDestMesh->pEdges + cornerBase;
	memcpy(facesStart, pSrcMesh->pFaces, sizeof(I32) * pSrcMesh->faceCount);
	bulkCopyAttribs(&pSrcMesh->faceAttribs, pDestMesh->faceCount,
	                &pDestMesh->faceAttribs, pSrcMesh->faceCount);
	pDestMesh->faceCount += pSrcMesh->faceCount;
	pDestMesh->pFaces[pDestMesh->faceCount] = pSrcMesh->pFaces[pSrcMesh->faceCount];
	memcpy(cornersStart, pSrcMesh->pCorners, sizeof(I32) * pSrcMesh->cornerCount);
	bulkCopyAttribs(&pSrcMesh->cornerAttribs, pDestMesh->cornerCount,
	                &pDestMesh->cornerAttribs, pSrcMesh->cornerCount);
	pDestMesh->cornerCount += pSrcMesh->cornerCount;
	memcpy(edgesStart, pSrcMesh->pEdges, sizeof(I32) * pSrcMesh->cornerCount);
	bulkCopyAttribs(&pSrcMesh->edgeAttribs, pDestMesh->edgeCount,
	                &pDestMesh->edgeAttribs, pSrcMesh->edgeCount);
	pDestMesh->edgeCount += pSrcMesh->edgeCount;
	bulkCopyAttribs(&pSrcMesh->vertAttribs, pDestMesh->vertCount,
	                &pDestMesh->vertAttribs, pSrcMesh->vertCount);
	pDestMesh->vertCount += pSrcMesh->vertCount;
	for (I32 i = faceBase; i < pDestMesh->faceCount; ++i) {
		pDestMesh->pFaces[i] += cornerBase;
	}
	pDestMesh->pFaces[pDestMesh->faceCount] += cornerBase;
	for (I32 i = cornerBase; i < pDestMesh->cornerCount; ++i) {
		pDestMesh->pCorners[i] += vertBase;
		pDestMesh->pEdges[i] += edgeBase;
	}
}

//move these to a separate Obj.c if more functions are made?

void stucApplyObjTransform(StucObject *pObj) {
	Mesh *pMesh = (Mesh *)pObj->pData;
	for (I32 i = 0; i < pMesh->core.vertCount; ++i) {
		V3_F32 *pV3 = pMesh->pVerts + i;
		V4_F32 v4 = {pV3->d[0], pV3->d[1], pV3->d[2], 1.0f};
		_(&v4 V4MULEQLM4X4 &pObj->transform);
		*pV3 = *(V3_F32 *)&v4;
	}
	for (I32 i = 0; i < pMesh->core.cornerCount; ++i) {
		Mat3x3 mat3x3 = Mat3x3FromMat4x4(&pObj->transform);
		V3_F32 *pNormal = pMesh->pNormals + i;
		_(pNormal V3MULEQLM3X3 &mat3x3);
		*pNormal = v3Normalize(*pNormal);
	}
	pObj->transform = identM4x4;
}

void stucMergeObjArr(StucContext pContext, Mesh *pMesh,
                     I32 objCount, StucObject *pObjArr, bool setCommon) {
	Mesh **ppSrcs = pContext->alloc.pCalloc(objCount, sizeof(void *));
	MeshCounts totalCount = {0};
	for (I32 i = 0; i < objCount; ++i) {
		ppSrcs[i] = (Mesh *)pObjArr[i].pData;
		stucAddToMeshCounts(pContext, &totalCount, NULL, (Mesh *)pObjArr[i].pData);
	}
	pMesh->faceBufSize = totalCount.faces + 1; //+1 for last face index
	pMesh->cornerBufSize = totalCount.corners;
	pMesh->edgeBufSize = totalCount.edges;
	pMesh->vertBufSize = totalCount.verts;
	pMesh->core.pFaces =
		pContext->alloc.pMalloc(sizeof(I32) * pMesh->faceBufSize);
	pMesh->core.pCorners =
		pContext->alloc.pMalloc(sizeof(I32) * pMesh->cornerBufSize);
	pMesh->core.pEdges =
		pContext->alloc.pMalloc(sizeof(I32) * pMesh->cornerBufSize);
	stucAllocAttribsFromMeshArr(&pContext->alloc, pMesh, objCount, ppSrcs, setCommon);
	for (I32 i = 0; i < objCount; ++i) {
		stucCopyMesh(&pMesh->core, (StucMesh *)pObjArr[i].pData);
	}
}

Result stucDestroyObjArr(StucContext pContext, I32 objCount, StucObject *pObjArr) {
	StucResult err = STUC_NOT_SET;
	for (I32 i = 0; i < objCount; ++i) {
		err = stucMeshDestroy(pContext, (StucMesh *)pObjArr[i].pData);
		STUC_THROW_IF(err, true, "", 0);
		pContext->alloc.pFree(pObjArr[i].pData);
	}
	pContext->alloc.pFree(pObjArr);
	STUC_CATCH(0, err, ;)
	return err;
}

FaceRange stucGetFaceRange(const StucMesh *pMesh,
                           const I32 idx, const bool border) {
	STUC_ASSERT("", border % 2 == border);
	I32 realIdx;
	I32 direction;
	if (!border) {
		realIdx = idx;
		direction = 1;
		STUC_ASSERT("", idx >= 0 && idx < pMesh->faceCount);
	}
	else {
		BufMeshIdx bufMeshIdx = stucConvertBorderFaceIdx(((BufMesh*)pMesh), idx);
		realIdx = bufMeshIdx.realIdx;
		direction = -1;
	}
	FaceRange face = {
		.idx = realIdx,
		.start = pMesh->pFaces[realIdx],
		.end = pMesh->pFaces[realIdx + direction]
	};
	if (!border) {
		STUC_ASSERT("", face.start >= 0 && face.end <= pMesh->cornerCount);
		STUC_ASSERT("", face.start < face.end);
		face.size = face.end - face.start;
	}
	else {
		BufMeshIdx start = stucConvertBorderCornerIdx(((BufMesh *)pMesh), face.start);
		BufMeshIdx end = stucConvertBorderCornerIdx(((BufMesh *)pMesh), face.end);
		face.start = start.realIdx;
		face.end = end.realIdx;
		STUC_ASSERT("", face.end >=
		                ((Mesh *)pMesh)->cornerBufSize - 1 -
		                ((BufMesh *)pMesh)->borderCornerCount);
		STUC_ASSERT("", face.end < face.start);
		STUC_ASSERT("", face.start < ((Mesh *)pMesh)->cornerBufSize);
		face.size = face.start - face.end;
	}
	STUC_ASSERT("", face.size >= 3);
	return face;
}

static
int mikktGetNumFaces(const SMikkTSpaceContext *pContext) {
	Mesh *pMesh = pContext->m_pUserData;
	return pMesh->core.faceCount;
}

static
int mikktGetNumVertsOfFace(const SMikkTSpaceContext *pContext, const int iFace) {
	Mesh *pMesh = pContext->m_pUserData;
	FaceRange face = stucGetFaceRange(&pMesh->core, iFace, false);
	return face.size;
}

static
void mikktGetPos(const SMikkTSpaceContext *pContext, F32 *pFvPosOut,
                 const int iFace, const int iVert) {
	Mesh *pMesh = pContext->m_pUserData;
	FaceRange face = stucGetFaceRange(&pMesh->core, iFace, false);
	I32 vertIdx = pMesh->core.pCorners[face.start + iVert];
	*(V3_F32 *)pFvPosOut = pMesh->pVerts[vertIdx];
}

static
void mikktGetNormal(const SMikkTSpaceContext *pContext, F32 *pFvNormOut,
                    const int iFace, const int iVert) {
	Mesh *pMesh = pContext->m_pUserData;
	FaceRange face = stucGetFaceRange(&pMesh->core, iFace, false);
	*(V3_F32 *)pFvNormOut = pMesh->pNormals[face.start + iVert];;
}

static
void mikktGetTexCoord(const SMikkTSpaceContext *pContext, F32 *pFvTexcOut,
                      const int iFace, const int iVert) {
	Mesh *pMesh = pContext->m_pUserData;
	FaceRange face = stucGetFaceRange(&pMesh->core, iFace, false);
	*(V2_F32 *)pFvTexcOut = pMesh->pUvs[face.start + iVert];
}

static
void mikktSetTSpaceBasic(const SMikkTSpaceContext *pContext, const F32 *pFvTangent,
                         const F32 fSign, const int iFace, const int iVert) {
	Mesh *pMesh = pContext->m_pUserData;
	FaceRange face = stucGetFaceRange(&pMesh->core, iFace, false);
	I32 corner = face.start + iVert;
	pMesh->pTangents[corner] = *(V3_F32 *)pFvTangent;
	pMesh->pTSigns[corner] = fSign;
}

void stucBuildTangents(Mesh *pMesh) {
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