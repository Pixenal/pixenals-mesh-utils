/* 
SPDX-FileCopyrightText: 2025 Caleb Dawson
SPDX-License-Identifier: Apache-2.0
*/

#include <stdint.h>
#include <string.h>

#include <error.h>
#include <attrib_utils.h>
#include <mesh.h>
#include <context.h>

typedef struct MeshDomain {
	StucDomain domain;
	I32 *pBufSize;
	I32 **ppList[2]; //2 so both pCorners and pEdges can be realloced at once
	I32 *pCount;
	AttribArray *pAttribArr;
} MeshDomain;

typedef struct BufMeshDomain {
	MeshDomain domain;
	I32 *pBorderCount;
} BufMeshDomain;

void stucCreateMesh(const StucContext pCtx, StucObject *pObj, StucObjectType type) {
	I32 size = 0;
	switch (type) {
		case STUC_OBJECT_DATA_MESH:
			size = sizeof(StucMesh);
			break;
		case STUC_OBJECT_DATA_MESH_INTERN:
			size = sizeof(Mesh);
			break;
		case STUC_OBJECT_DATA_MESH_BUF:
			size = sizeof(Mesh);
			break;
		default:
			STUC_ASSERT("Invalid object data type", false);
			return;
	}
	pObj->pData = pCtx->alloc.fpCalloc(1, size);
	((StucMesh *)pObj->pData)->type.type = type;
}

#ifndef TEMP_DISABLE
static
void reallocBufMesh(
	StucContext pCtx,
	BufMesh *pMesh,
	BufMeshDomain *pBufDomain
) {
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
			pCtx->alloc.fpRealloc(*pDomain->ppList[i], sizeof(I32) * *pDomain->pBufSize);
		if (*pBufDomain->pBorderCount) {
			I32 *pStart = *pDomain->ppList[i] + realBorderEnd + 1;
			memmove(pStart + diff, pStart, sizeof(I32) * *pBufDomain->pBorderCount);
			I32 newFirstElement = (*pDomain->ppList[i])[realBorderEnd + 1 + diff];
			I32 newLastElement = (*pDomain->ppList[i])[*pDomain->pBufSize - 1];
			STUC_ASSERT("", newFirstElement == oldFirstElement);
			STUC_ASSERT("", newLastElement == oldLastElement);
		}
	}
	stucReallocAndMoveAttribs(
		&pCtx->alloc,
		pMesh,
		pDomain->pAttribArr,
		realBorderEnd + 1,
		diff,
		*pBufDomain->pBorderCount,
		*pDomain->pBufSize
	);
	//0xffffffff for set all
	stucAssignActiveAliases(pCtx, &pMesh->mesh, 0xffffffff, pBufDomain->domain.domain);
}

static
BufMeshIdx getNewBufMeshIdx(
	StucContext pCtx,
	BufMesh *pMesh,
	BufMeshDomain *pBufDomain,
	const bool border,
	bool *pRealloced
) {
	MeshDomain *pDomain = (MeshDomain *)pBufDomain;
	I32 realBorderEnd = *pDomain->pBufSize - 1 - *pBufDomain->pBorderCount;
	//TODO assertions like these need to be converted to release exceptions
	STUC_ASSERT("", *pDomain->pCount <= realBorderEnd);
	if (*pDomain->pCount == realBorderEnd) {
		reallocBufMesh(pCtx, pMesh, pBufDomain);
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
#endif

static
void getFaceDomain(Mesh *pMesh, MeshDomain *pDomain) {
	pDomain->domain = STUC_DOMAIN_FACE;
	pDomain->pBufSize = &pMesh->faceBufSize;
	pDomain->ppList[0] = &pMesh->core.pFaces;
	pDomain->pCount = &pMesh->core.faceCount;
	pDomain->pAttribArr = &pMesh->core.faceAttribs;
}

static
void getCornerDomain(Mesh *pMesh, MeshDomain *pDomain) {
	pDomain->domain = STUC_DOMAIN_CORNER;
	pDomain->pBufSize = &pMesh->cornerBufSize;
	pDomain->ppList[0] = &pMesh->core.pCorners;
	pDomain->ppList[1] = &pMesh->core.pEdges;
	pDomain->pCount = &pMesh->core.cornerCount;
	pDomain->pAttribArr = &pMesh->core.cornerAttribs;
}

static
void getEdgeDomain(Mesh *pMesh, MeshDomain *pDomain) {
	pDomain->domain = STUC_DOMAIN_EDGE;
	pDomain->pBufSize = &pMesh->edgeBufSize;
	pDomain->ppList[0] = NULL;
	pDomain->pCount = &pMesh->core.edgeCount;
	pDomain->pAttribArr = &pMesh->core.edgeAttribs;
}

static
void getVertDomain(Mesh *pMesh, MeshDomain *pDomain) {
	pDomain->domain = STUC_DOMAIN_VERT;
	pDomain->pBufSize = &pMesh->vertBufSize;
	pDomain->ppList[0] = NULL;
	pDomain->pCount = &pMesh->core.vertCount;
	pDomain->pAttribArr = &pMesh->core.vertAttribs;
}

static
void reallocMesh(StucContext pCtx, Mesh *pMesh, MeshDomain *pDomain) {
	*pDomain->pBufSize *= 2;
	for (I32 i = 0; i < 2; ++i) {
		if (!pDomain->ppList[i]) {
			continue;
		}
		*pDomain->ppList[i] =
			pCtx->alloc.fpRealloc(*pDomain->ppList[i], *pDomain->pBufSize * sizeof(I32));
		/*
		if (pMesh->core.type.type == STUC_OBJECT_DATA_MESH_BUF &&
			pDomain->domain == STUC_DOMAIN_VERT
		) {
			BufVertInfoArr *pVertInfo = &((BufMesh *)pMesh)->vertInfo;
			pVertInfo->pArr = pCtx->alloc.fpRealloc(
				pVertInfo->pArr, pMesh->vertBufSize * sizeof(BufVertInfo)
			);
			pVertInfo->pTypeArr =
				pCtx->alloc.fpRealloc(pVertInfo->pTypeArr, pMesh->vertBufSize);
		}
		*/
	}
	stucReallocAttribArr(
		pCtx,
		pDomain->domain,
		pMesh,
		pDomain->pAttribArr,
		*pDomain->pBufSize
	);
}

static
I32 getNewMeshIdx(
	StucContext pCtx,
	Mesh *pMesh,
	MeshDomain *pDomain,
	bool *pRealloced
) {
	STUC_ASSERT("", *pDomain->pCount >= 0 && *pDomain->pBufSize > 0);
	STUC_ASSERT("", *pDomain->pCount <= *pDomain->pBufSize);
	if (*pDomain->pCount == *pDomain->pBufSize) {
		reallocMesh(pCtx, pMesh, pDomain);
		if (pRealloced) {
			*pRealloced = true;
		}
	}
	else {
		if (pRealloced) {
			*pRealloced = false;
		}
	}
	I32 idx = *pDomain->pCount;
	++*pDomain->pCount;
	return idx;
}

#ifndef TEMP_DISABLE
BufMeshIdx stucBufMeshAddFace(
	StucContext pCtx,
	BufMesh *pMesh,
	bool border,
	bool *pRealloced
) {
	BufMeshDomain domain = {0};
	getFaceDomain((Mesh *)pMesh, (MeshDomain *)&domain);
	domain.pBorderCount = &pMesh->borderFaceCount;
	BufMeshIdx idx =
		getNewBufMeshIdx(pCtx, pMesh, &domain, border, pRealloced);
	return idx;
}

BufMeshIdx stucBufMeshAddCorner(
	StucContext pCtx,
	BufMesh *pMesh,
	bool border,
	bool *pRealloced
) {
	BufMeshDomain domain = {0};
	getCornerDomain((Mesh *)pMesh, (MeshDomain *)&domain);
	domain.pBorderCount = &pMesh->borderCornerCount;
	BufMeshIdx idx =
		getNewBufMeshIdx(pCtx, pMesh, &domain, border, pRealloced);
	return idx;
}

BufMeshIdx stucBufMeshAddEdge(
	StucContext pCtx,
	BufMesh *pMesh,
	bool border,
	bool *pRealloced
) {
	BufMeshDomain domain = {0};
	getEdgeDomain((Mesh *)pMesh, (MeshDomain *)&domain);
	domain.pBorderCount = &pMesh->borderEdgeCount;
	BufMeshIdx idx =
		getNewBufMeshIdx(pCtx, pMesh, &domain, border, pRealloced);
	return idx;
}

BufMeshIdx stucBufMeshAddVert(
	StucContext pCtx,
	BufMesh *pMesh,
	bool border,
	bool *pRealloced
) {
	BufMeshDomain domain = {0};
	getVertDomain((Mesh *)pMesh, (MeshDomain *)&domain);
	domain.pBorderCount = &pMesh->borderVertCount;
	BufMeshIdx idx =
		getNewBufMeshIdx(pCtx, pMesh, &domain, border, pRealloced);
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
#endif

I32 stucMeshAddFace(const StucContext pCtx, Mesh *pMesh, bool *pRealloced) {
	MeshDomain domain = {0};
	getFaceDomain(pMesh, &domain);
	return getNewMeshIdx(pCtx, pMesh, &domain, pRealloced);
}

I32 stucMeshAddCorner(const StucContext pCtx, Mesh *pMesh, bool *pRealloced) {
	MeshDomain domain = {0};
	getCornerDomain(pMesh, &domain);
	return getNewMeshIdx(pCtx, pMesh, &domain, pRealloced);
}

I32 stucMeshAddEdge(const StucContext pCtx, Mesh *pMesh, bool *pRealloced) {
	MeshDomain domain = {0};
	getEdgeDomain(pMesh, &domain);
	return getNewMeshIdx(pCtx, pMesh, &domain, pRealloced);
}

I32 stucMeshAddVert(const StucContext pCtx, Mesh *pMesh, bool *pRealloced) {
	MeshDomain domain = {0};
	getVertDomain(pMesh, &domain);
	return getNewMeshIdx(pCtx, pMesh, &domain, pRealloced);
}

void stucReallocMeshToFit(const StucContext pCtx, Mesh *pMesh) {
	StucMesh *pCore = &pMesh->core;
	pMesh->faceBufSize = pCore->faceCount + 1;
	I32 newLen = sizeof(I32) * pMesh->faceBufSize;
	pCore->pFaces = pCtx->alloc.fpRealloc(pCore->pFaces, newLen);
	pMesh->cornerBufSize = pCore->cornerCount;
	newLen = sizeof(I32) * pMesh->cornerBufSize;
	pCore->pCorners = pCtx->alloc.fpRealloc(pCore->pCorners, newLen);
	pCore->pEdges = pCtx->alloc.fpRealloc(pCore->pEdges, newLen);
	stucReallocAttribArr(
		pCtx,
		STUC_DOMAIN_FACE,
		pMesh,
		&pCore->faceAttribs,
		pMesh->faceBufSize
	);
	stucReallocAttribArr(
		pCtx,
		STUC_DOMAIN_CORNER,
		pMesh,
		&pCore->cornerAttribs,
		pMesh->cornerBufSize
	);
	pMesh->edgeBufSize = pCore->edgeCount;
	stucReallocAttribArr(
		pCtx,
		STUC_DOMAIN_EDGE,
		pMesh,
		&pCore->edgeAttribs,
		pMesh->edgeBufSize
	);
	pMesh->vertBufSize = pCore->vertCount;
	stucReallocAttribArr(
		pCtx,
		STUC_DOMAIN_VERT,
		pMesh,
		&pCore->vertAttribs,
		pMesh->vertBufSize
	);
}

void stucMeshSetLastFace(const StucContext pCtx, Mesh *pMesh) {
	I32 lastFace = stucMeshAddFace(pCtx, pMesh, NULL);
	pMesh->core.pFaces[lastFace] = pMesh->core.cornerCount;
	//meshAddFace() increments this, so we need to undo that
	pMesh->core.faceCount--; 
}

#ifndef TEMP_DISABLE
void stucBufMeshSetLastFaces(StucContext pCtx, BufMesh *pBufMesh) {
	Mesh *pMesh = &pBufMesh->mesh;
	bool realloced = false;
	BufMeshIdx lastFace =
		stucBufMeshAddFace(pCtx, pBufMesh, false, &realloced);
	pMesh->core.pFaces[lastFace.realIdx] = pMesh->core.cornerCount;
	//bufMeshAddFace() increments this, so we need to undo that
	pMesh->core.faceCount--; 

	lastFace = stucBufMeshAddFace(pCtx, pBufMesh, true, &realloced);
	pMesh->core.pFaces[lastFace.realIdx] = pBufMesh->borderCornerCount;
	pBufMesh->borderFaceCount--;
}
#endif

bool stucCheckIfMesh(const StucObjectData type) {
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
void bulkCopyAttribs(
	const AttribArray *pSrc,
	I32 SrcOffset,
	AttribArray *pDest,
	I32 dataLen
) {
	for (I32 i = 0; i < pDest->count; ++i) {
		const Attrib *pSrcAttrib = 
			stucGetAttribInternConst(pDest->pArr[i].core.name, pSrc, false, NULL, NULL);
		if (!pSrcAttrib) {
			continue;
		}
		void *attribDestStart = stucAttribAsVoid(&pDest->pArr[i].core, SrcOffset);
		I32 attribTypeSize = stucGetAttribSizeIntern(pDest->pArr[i].core.type);
		memcpy(attribDestStart, pSrcAttrib->core.pData, attribTypeSize * dataLen);
	}
}

void stucAddToMeshCounts(
	MeshCounts *pCounts,
	MeshCounts *pBoundsCounts,
	const Mesh *pMeshSrc
) {
	//TODO maybe replace *Counts vars in Mesh to use MeshCounts,
	//     so we can just do:
	//     meshCountsAdd(totalCount, pBufMesh->mesh.meshCounts);
	//     or something.
	pCounts->faces += pMeshSrc->core.faceCount;
	pCounts->corners += pMeshSrc->core.cornerCount;
	pCounts->edges += pMeshSrc->core.edgeCount;
	pCounts->verts += pMeshSrc->core.vertCount;
}

Result stucCopyMesh(StucMesh *pDestMesh, const StucMesh *pSrcMesh) {
	Result err = STUC_SUCCESS;
	if (pSrcMesh->type.type == STUC_OBJECT_DATA_NULL) {
		//TODO why doesn't this return STUC_ERROR?
		return err;
	}
	STUC_RETURN_ERR_IFNOT_COND(err, stucCheckIfMesh(pDestMesh->type), "");
	STUC_RETURN_ERR_IFNOT_COND(err, stucCheckIfMesh(pSrcMesh->type), "");
	I32 faceBase = pDestMesh->faceCount;
	I32 cornerBase = pDestMesh->cornerCount;
	I32 edgeBase = pDestMesh->edgeCount;
	I32 vertBase = pDestMesh->vertCount;
	I32 *facesStart = pDestMesh->pFaces + faceBase;
	I32 *cornersStart = pDestMesh->pCorners + cornerBase;
	I32 *edgesStart = pDestMesh->pEdges + cornerBase;
	memcpy(facesStart, pSrcMesh->pFaces, sizeof(I32) * pSrcMesh->faceCount);
	bulkCopyAttribs(
		&pSrcMesh->faceAttribs,
		pDestMesh->faceCount,
		&pDestMesh->faceAttribs,
		pSrcMesh->faceCount
	);
	pDestMesh->faceCount += pSrcMesh->faceCount;
	pDestMesh->pFaces[pDestMesh->faceCount] = pSrcMesh->pFaces[pSrcMesh->faceCount];
	memcpy(cornersStart, pSrcMesh->pCorners, sizeof(I32) * pSrcMesh->cornerCount);
	bulkCopyAttribs(
		&pSrcMesh->cornerAttribs,
		pDestMesh->cornerCount,
		&pDestMesh->cornerAttribs,
		pSrcMesh->cornerCount
	);
	pDestMesh->cornerCount += pSrcMesh->cornerCount;
	memcpy(edgesStart, pSrcMesh->pEdges, sizeof(I32) * pSrcMesh->cornerCount);
	bulkCopyAttribs(
		&pSrcMesh->edgeAttribs,
		pDestMesh->edgeCount,
		&pDestMesh->edgeAttribs,
		pSrcMesh->edgeCount
	);
	pDestMesh->edgeCount += pSrcMesh->edgeCount;
	bulkCopyAttribs(
		&pSrcMesh->vertAttribs,
		pDestMesh->vertCount,
		&pDestMesh->vertAttribs,
		pSrcMesh->vertCount
	);
	pDestMesh->vertCount += pSrcMesh->vertCount;
	for (I32 i = faceBase; i < pDestMesh->faceCount; ++i) {
		pDestMesh->pFaces[i] += cornerBase;
	}
	pDestMesh->pFaces[pDestMesh->faceCount] += cornerBase;
	for (I32 i = cornerBase; i < pDestMesh->cornerCount; ++i) {
		pDestMesh->pCorners[i] += vertBase;
		pDestMesh->pEdges[i] += edgeBase;
	}
	return err;
}

//move these to a separate Obj.c if more functions are made?

void stucApplyObjTransform(StucObject *pObj) {
	Mesh *pMesh = (Mesh *)pObj->pData;
	for (I32 i = 0; i < pMesh->core.vertCount; ++i) {
		V3_F32 *pV3 = pMesh->pPos + i;
		V4_F32 v4 = {pV3->d[0], pV3->d[1], pV3->d[2], 1.0f};
		_(&v4 V4MULEQLM4X4 &pObj->transform);
		*pV3 = *(V3_F32 *)&v4;
	}
	if (pMesh->pNormals) {
		for (I32 i = 0; i < pMesh->core.cornerCount; ++i) {
			Mat3x3 mat3x3 = Mat3x3FromMat4x4(&pObj->transform);
			V3_F32 *pNormal = pMesh->pNormals + i;
			_(pNormal V3MULEQLM3X3 & mat3x3);
			*pNormal = v3F32Normalize(*pNormal);
		}
	}
	pObj->transform = STUC_IDENT_MAT4X4;
}

Result stucMergeObjArr(
	StucContext pCtx,
	Mesh *pMesh,
	I32 objCount,
	const StucObject *pObjArr,
	bool setCommon
) {
	Result err = STUC_SUCCESS;
	const Mesh **ppSrcs = pCtx->alloc.fpCalloc(objCount, sizeof(void *));
	MeshCounts totalCount = {0};
	for (I32 i = 0; i < objCount; ++i) {
		ppSrcs[i] = (Mesh *)pObjArr[i].pData;
		stucAddToMeshCounts(&totalCount, NULL, (Mesh *)pObjArr[i].pData);
	}
	pMesh->faceBufSize = totalCount.faces + 1; //+1 for last face index
	pMesh->cornerBufSize = totalCount.corners;
	pMesh->edgeBufSize = totalCount.edges;
	pMesh->vertBufSize = totalCount.verts;
	pMesh->core.pFaces =
		pCtx->alloc.fpMalloc(sizeof(I32) * pMesh->faceBufSize);
	pMesh->core.pCorners =
		pCtx->alloc.fpMalloc(sizeof(I32) * pMesh->cornerBufSize);
	pMesh->core.pEdges =
		pCtx->alloc.fpMalloc(sizeof(I32) * pMesh->cornerBufSize);
	err = stucAllocAttribsFromMeshArr(pCtx, pMesh, objCount, ppSrcs, -1, setCommon, true, false);
	STUC_THROW_IFNOT(err, "", 0);
	for (I32 i = 0; i < objCount; ++i) {
		stucCopyMesh(&pMesh->core, (StucMesh *)pObjArr[i].pData);
	}
	STUC_CATCH(0, err,
		stucMeshDestroy(pCtx, &pMesh->core);
	);
	return err;
}

Result stucDestroyObjArr(StucContext pCtx, I32 objCount, StucObject *pObjArr) {
	StucResult err = STUC_NOT_SET;
	for (I32 i = 0; i < objCount; ++i) {
		err = stucMeshDestroy(pCtx, (StucMesh *)pObjArr[i].pData);
		STUC_THROW_IFNOT(err, "", 0);
		pCtx->alloc.fpFree(pObjArr[i].pData);
	}
	pCtx->alloc.fpFree(pObjArr);
	STUC_CATCH(0, err, ;)
	return err;
}

Result stucValidateMesh(const StucMesh *pMesh, bool checkEdges) {
	Result err = STUC_SUCCESS;
	STUC_RETURN_ERR_IFNOT_COND(err, pMesh->faceCount && pMesh->pFaces, "");
	STUC_RETURN_ERR_IFNOT_COND(err, pMesh->cornerCount && pMesh->pCorners, "");
	STUC_RETURN_ERR_IFNOT_COND(err, !checkEdges || (pMesh->edgeCount && pMesh->pEdges), "");
	STUC_RETURN_ERR_IFNOT_COND(err, pMesh->vertCount, "");
	STUC_RETURN_ERR_IFNOT_COND(
		err,
		(pMesh->meshAttribs.pArr && pMesh->meshAttribs.count) ||
		!pMesh->meshAttribs.count,
		""
	);
	STUC_RETURN_ERR_IFNOT_COND(
		err,
		(pMesh->faceAttribs.pArr && pMesh->faceAttribs.count) ||
		!pMesh->faceAttribs.count,
		""
	);
	STUC_RETURN_ERR_IFNOT_COND(
		err,
		pMesh->cornerAttribs.pArr && pMesh->cornerAttribs.count,
		""
	);
	STUC_RETURN_ERR_IFNOT_COND(
		err,
		(pMesh->edgeAttribs.pArr && pMesh->edgeAttribs.count) ||
		!pMesh->edgeAttribs.count,
		""
	);
	STUC_RETURN_ERR_IFNOT_COND(
		err,
		pMesh->vertAttribs.pArr && pMesh->vertAttribs.count,
		""
	);
	for (I32 i = 0; i < pMesh->faceCount; ++i) {
		FaceRange face = stucGetFaceRange(pMesh, i);
		STUC_RETURN_ERR_IFNOT_COND(err, face.size >= 3, "");
		for (I32 j = 0; j < face.size; ++j) {
			I32 corner = face.start + j;
			STUC_RETURN_ERR_IFNOT_COND(err, corner < pMesh->cornerCount, "");
			STUC_RETURN_ERR_IFNOT_COND(err, pMesh->pCorners[corner] < pMesh->vertCount, "");
		}
	}
	for (I32 i = 1; i < STUC_ATTRIB_USE_ENUM_COUNT; ++i) {
		if (i == STUC_ATTRIB_USE_SP_ENUM_COUNT) {
			continue;
		}
		AttribActive idx = pMesh->activeAttribs[i];
		if (!idx.active) {
			STUC_RETURN_ERR_IFNOT_COND(err,
				!stucIsAttribUseRequired(i),
				"in-mesh must have active attribs for pos, uv, normal, and idx"
			);
			continue;
		}
		STUC_RETURN_ERR_IFNOT_COND(err, idx.idx >= 0, "invalid active attrib index");
		STUC_RETURN_ERR_IFNOT_COND(
			err,
			idx.domain >= 0 && idx.domain <= STUC_DOMAIN_VERT,
			"invalid active attrib domain"
		);
	}
	return err;
}

void stucAliasMeshCoreNoAttribs(StucMesh *pDest, StucMesh *pSrc) {
	pDest->type = pSrc->type;
	pDest->pFaces = pSrc->pFaces;
	pDest->pCorners = pSrc->pCorners;
	pDest->pEdges = pSrc->pEdges;
	pDest->faceCount = pSrc->faceCount;
	pDest->cornerCount = pSrc->cornerCount;
	pDest->edgeCount = pSrc->edgeCount;
	pDest->vertCount = pSrc->vertCount;
}

I32 stucGetDomainSize(const Mesh *pMesh, StucDomain domain) {
	switch (domain) {
		case STUC_DOMAIN_FACE:
			return pMesh->faceBufSize;
		case STUC_DOMAIN_CORNER:
			return pMesh->cornerBufSize;
		case STUC_DOMAIN_EDGE:
			return pMesh->edgeBufSize;
		case STUC_DOMAIN_VERT:
			return pMesh->vertBufSize;
		default:
			STUC_ASSERT("invalid domain", false);
	}
	return 0;
}

I32 stucDomainCountGetIntern(const StucMesh *pMesh, StucDomain domain) {
	switch (domain) {
		case STUC_DOMAIN_FACE:
			return pMesh->faceCount;
		case STUC_DOMAIN_CORNER:
			return pMesh->cornerCount;
		case STUC_DOMAIN_EDGE:
			return pMesh->edgeCount;
		case STUC_DOMAIN_VERT:
			return pMesh->vertCount;
		default:
			STUC_ASSERT("invalid domain", false);
	}
	STUC_ASSERT("invalid domain", false);
	return 0;
}

#ifndef TEMP_DISABLE
I32 stucGetVirtualBufIdx(BufMesh *pBufMesh, I32 corner) {
	I32 idxVirtual = pBufMesh->mesh.cornerBufSize - corner - 1;
	STUC_ASSERT(
		"",
		idxVirtual >= 0 && idxVirtual < pBufMesh->borderCornerCount
	);
	return idxVirtual;
}
#endif

I32 stucGetCornerPrev(I32 corner, const FaceRange *pFace) {
	I32 prev = corner ? corner - 1 : pFace->size - 1;
	STUC_ASSERT("", prev >= 0 && prev < pFace->size);
	return prev;
}

I32 stucGetCornerNext(I32 corner, const FaceRange *pFace) {
	STUC_ASSERT("", corner < pFace->size);
	I32 next = (corner + 1) % pFace->size;
	STUC_ASSERT("", next >= 0);
	return next;
}

bool stucGetIfSeamEdge(const Mesh *pMesh, I32 edge) {
	STUC_ASSERT("", pMesh->pSeamEdge);
	return pMesh->pSeamEdge[edge];
}

bool stucGetIfMatBorderEdge(const Mesh *pMesh, I32 edge) {
	STUC_ASSERT("", pMesh->pEdgeFaces && pMesh->pMatIdx);
	V2_I32 faces = pMesh->pEdgeFaces[edge];
	if (faces.d[1] == -1) {
		return false; //no adj face
	}
	return pMesh->pMatIdx[faces.d[0]] != pMesh->pMatIdx[faces.d[1]];
}

//does not increment the other corner, so it may not actually be adjacent.
//Do so after calling this func (if wind is equal between faces)
void stucGetAdjCorner(const Mesh *pMesh, FaceCorner corner, FaceCorner *pAdjCorner) {
	STUC_ASSERT("", pMesh->pEdgeFaces);
	FaceRange face = stucGetFaceRange(&pMesh->core, corner.face);
	I32 edge = pMesh->core.pEdges[face.start + corner.corner];
	V2_I8 corners = pMesh->pEdgeCorners[edge];
	V2_I32 faces = pMesh->pEdgeFaces[edge];
	bool which = corners.d[0] == corner.corner && faces.d[0] == corner.face;
	STUC_ASSERT(
		"",
		faces.d[!which] == -1 && corners.d[!which] == -1 ||
		faces.d[!which] == corner.face && corners.d[!which] == corner.corner
	);
	pAdjCorner->corner = corners.d[which];
	pAdjCorner->face = faces.d[which];
}

V2_F32 stucGetVertPosAsV2(const Mesh *pMesh, const FaceRange *pFace, I32 corner) {
	STUC_ASSERT("", corner >= 0 && corner < pFace->size);
	return *(V2_F32 *)&pMesh->pPos[pMesh->core.pCorners[pFace->start + corner]];
}

V2_F32 stucGetUvPos(const Mesh *pMesh, const FaceRange *pFace, I32 corner) {
	STUC_ASSERT("", corner >= 0 && corner < pFace->size);
	return pMesh->pUvs[pFace->start + corner];
}

I32 stucGetMeshVert(const StucMesh *pMesh, FaceCorner corner) {
	STUC_ASSERT("", pMesh && corner.face >= 0 && corner.face < pMesh->faceCount);
	FaceRange face = stucGetFaceRange(pMesh, corner.face);
	STUC_ASSERT("", corner.corner >= 0 && corner.corner < face.size);
	I32 vert = pMesh->pCorners[face.start + corner.corner];
	STUC_ASSERT("", vert >= 0 && vert < pMesh->vertCount);
	return vert;
}

I32 stucGetMeshEdge(const StucMesh *pMesh, FaceCorner corner) {
	STUC_ASSERT("", pMesh && corner.face >= 0 && corner.face < pMesh->faceCount);
	FaceRange face = stucGetFaceRange(pMesh, corner.face);
	STUC_ASSERT("", corner.corner >= 0 && corner.corner < face.size);
	I32 edge = pMesh->pEdges[face.start + corner.corner];
	STUC_ASSERT("", edge >= 0 && edge < pMesh->edgeCount);
	return edge;
}

bool checkForNgonsInMesh(const StucMesh *pMesh) {
	for (I32 i = 0; i < pMesh->faceCount; ++i) {
		FaceRange face = stucGetFaceRange(pMesh, i);
		if (face.size > 4) {
			return true;
		}
	}
	return false;
}
