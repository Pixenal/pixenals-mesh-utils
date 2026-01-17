/* 
SPDX-FileCopyrightText: 2025 Caleb Dawson
SPDX-License-Identifier: Apache-2.0
*/

#include <stdint.h>
#include <string.h>

#include <pixenals_error_utils.h>
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
			PIX_ERR_ASSERT("Invalid object data type", false);
			return;
	}
	pObj->pData = pCtx->alloc.fpCalloc(1, size);
	((StucMesh *)pObj->pData)->type.type = type;
}

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
	PIX_ERR_ASSERT("", *pDomain->pCount >= 0 && *pDomain->pBufSize > 0);
	PIX_ERR_ASSERT("", *pDomain->pCount <= *pDomain->pBufSize);
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
			stucGetAttribInternConst(pDest->pArr[i].core.name, pSrc, false, NULL, NULL, NULL);
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

StucErr stucCopyMesh(StucContext pCtx, StucMesh *pDest, const StucMesh *pSrc) {
	StucErr err = PIX_ERR_SUCCESS;
	if (pSrc->type.type == STUC_OBJECT_DATA_NULL) {
		//TODO why doesn't this return PIX_ERR_ERROR?
		return err;
	}
	PIX_ERR_RETURN_IFNOT_COND(err, stucCheckIfMesh(pDest->type), "");
	PIX_ERR_RETURN_IFNOT_COND(err, stucCheckIfMesh(pSrc->type), "");
	I32 faceBase = pDest->faceCount;
	I32 cornerBase = pDest->cornerCount;
	I32 edgeBase = pDest->edgeCount;
	I32 vertBase = pDest->vertCount;
	I32 *facesStart = pDest->pFaces + faceBase;
	I32 *cornersStart = pDest->pCorners + cornerBase;
	memcpy(facesStart, pSrc->pFaces, sizeof(I32) * pSrc->faceCount);
	bulkCopyAttribs(
		&pSrc->faceAttribs,
		pDest->faceCount,
		&pDest->faceAttribs,
		pSrc->faceCount
	);
	pDest->faceCount += pSrc->faceCount;
	pDest->pFaces[pDest->faceCount] = pSrc->pFaces[pSrc->faceCount];
	memcpy(cornersStart, pSrc->pCorners, sizeof(I32) * pSrc->cornerCount);
	bulkCopyAttribs(
		&pSrc->cornerAttribs,
		pDest->cornerCount,
		&pDest->cornerAttribs,
		pSrc->cornerCount
	);
	pDest->cornerCount += pSrc->cornerCount;
	if (pDest->pEdges) {
		I32 *edgesStart = pDest->pEdges + cornerBase;
		memcpy(edgesStart, pSrc->pEdges, sizeof(I32) * pSrc->cornerCount);
	}
	bulkCopyAttribs(
		&pSrc->edgeAttribs,
		pDest->edgeCount,
		&pDest->edgeAttribs,
		pSrc->edgeCount
	);
	pDest->edgeCount += pSrc->edgeCount;
	bulkCopyAttribs(
		&pSrc->vertAttribs,
		pDest->vertCount,
		&pDest->vertAttribs,
		pSrc->vertCount
	);
	pDest->vertCount += pSrc->vertCount;
	if (cornerBase) {
		for (I32 i = faceBase; i < pDest->faceCount; ++i) {
			pDest->pFaces[i] += cornerBase;
		}
		pDest->pFaces[pDest->faceCount] += cornerBase;
	}
	if (vertBase) {
		for (I32 i = cornerBase; i < pDest->cornerCount; ++i) {
			pDest->pCorners[i] += vertBase;
		}
	}
	if (edgeBase) {
		if (pDest->pEdges) {
			for (I32 i = cornerBase; i < pDest->cornerCount; ++i) {
				pDest->pEdges[i] += edgeBase;
			}
		}
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
			M3x3 m3x3 = pixmM3x3FromM4x4(&pObj->transform);
			V3_F32 *pNormal = pMesh->pNormals + i;
			_(pNormal V3MULEQLM3X3 & m3x3);
			*pNormal = pixmV3F32Normalize(*pNormal);
		}
	}
	pObj->transform = PIX_MATH_IDENT_MAT4X4;
}

StucErr stucMergeObjArr(
	StucContext pCtx,
	Mesh *pMesh,
	const StucObjArr *pObjArr,
	bool setCommon
) {
	StucErr err = PIX_ERR_SUCCESS;
	const Mesh **ppSrcs = pCtx->alloc.fpCalloc(pObjArr->count, sizeof(void *));
	MeshCounts totalCount = {0};
	for (I32 i = 0; i < pObjArr->count; ++i) {
		ppSrcs[i] = (Mesh *)pObjArr->pArr[i].pData;
		stucAddToMeshCounts(&totalCount, NULL, (Mesh *)pObjArr->pArr[i].pData);
	}
	pMesh->faceBufSize = totalCount.faces + 1; //+1 for last face index
	pMesh->cornerBufSize = totalCount.corners;
	pMesh->edgeBufSize = totalCount.edges;
	pMesh->vertBufSize = totalCount.verts;
	pMesh->core.pFaces = pCtx->alloc.fpMalloc(sizeof(I32) * pMesh->faceBufSize);
	pMesh->core.pCorners = pCtx->alloc.fpMalloc(sizeof(I32) * pMesh->cornerBufSize);
	pMesh->core.pEdges = pCtx->alloc.fpMalloc(sizeof(I32) * pMesh->cornerBufSize);
	err = stucAllocAttribsFromMeshArr(
		pCtx,
		pMesh,
		pObjArr->count,
		ppSrcs,
		-1,
		setCommon,
		true,
		false,
		false
	);
	PIX_ERR_THROW_IFNOT(err, "", 0);
	for (I32 i = 0; i < pObjArr->count; ++i) {
		err = stucCopyMesh(pCtx, &pMesh->core, (StucMesh *)pObjArr->pArr[i].pData);
		PIX_ERR_THROW_IFNOT(err, "", 0);
	}
	PIX_ERR_CATCH(0, err,
		stucMeshDestroy(pCtx, &pMesh->core);
	);
	return err;
}

StucErr stucMeshAllocCopy(
	StucContext pCtx,
	StucMesh *pDest,
	const StucMesh *pSrc,
	bool activeOnly
) {
	StucErr err = PIX_ERR_SUCCESS;
	Mesh destWrap = {
		.core.type.type = STUC_OBJECT_DATA_MESH,
		.faceBufSize = pSrc->faceCount + 1,
		.cornerBufSize = pSrc->cornerCount,
		.edgeBufSize = pSrc->edgeCount,
		.vertBufSize = pSrc->vertCount
	};
	if (pSrc->faceCount) {
		destWrap.core.pFaces = pCtx->alloc.fpMalloc(sizeof(I32) * destWrap.faceBufSize);
	}
	if (pSrc->cornerCount) {
		destWrap.core.pCorners =
			pCtx->alloc.fpMalloc(sizeof(I32) * destWrap.cornerBufSize);
		if (pSrc->pEdges) {
			destWrap.core.pEdges =
				pCtx->alloc.fpMalloc(sizeof(I32) * destWrap.cornerBufSize);
		}
	}
	Mesh srcWrap = {.core = *pSrc};
	err = stucAttemptToSetMissingActiveDomains(&srcWrap.core);
	PIX_ERR_RETURN_IFNOT(err, "");
	err = stucAssignActiveAliases(
		pCtx,
		&srcWrap,
		0xffffffff,
		STUC_DOMAIN_NONE
	);
	PIX_ERR_RETURN_IFNOT(err, "");
	const Mesh *pSrcWrap = &srcWrap;

	err = stucAllocAttribsFromMeshArr(
		pCtx,
		&destWrap,
		1,
		&pSrcWrap,
		-1,
		false, true, false, activeOnly
	);
	PIX_ERR_THROW_IFNOT(err, "", 0);
	*pDest = destWrap.core;
	err = stucCopyMesh(pCtx, pDest, pSrc);
	PIX_ERR_CATCH(0, err,
		stucMeshDestroy(pCtx, pDest);
	);
	return err;
}

StucErr stucObjArrDestroy(const StucContext pCtx, StucObjArr *pArr) {
	StucErr err = PIX_ERR_SUCCESS;
	PIX_ERR_RETURN_IFNOT_COND(err, pCtx && pArr, "");
	for (I32 i = 0; i < pArr->count; ++i) {
		if (pArr->pArr[i].pData) {
			StucMesh *pMesh = (StucMesh *)pArr->pArr[i].pData;
			stucMeshDestroy(pCtx, pMesh);
			pCtx->alloc.fpFree(pMesh);
		}
	}
	pCtx->alloc.fpFree(pArr->pArr);
	*pArr = (StucObjArr){0};
	return err;
}

StucErr stucValidateMesh(
	const StucAlloc *pAlloc,
	const StucMesh *pMesh,
	bool checkEdges,
	bool posOnly
) {
	StucErr err = PIX_ERR_SUCCESS;
	PIX_ERR_RETURN_IFNOT_COND(err, pMesh->faceCount && pMesh->pFaces, "");
	PIX_ERR_RETURN_IFNOT_COND(err, pMesh->cornerCount && pMesh->pCorners, "");
	PIX_ERR_RETURN_IFNOT_COND(err, !checkEdges || (pMesh->edgeCount && pMesh->pEdges), "");
	PIX_ERR_RETURN_IFNOT_COND(err, pMesh->vertCount, "");
	PIX_ERR_RETURN_IFNOT_COND(
		err,
		(pMesh->meshAttribs.pArr && pMesh->meshAttribs.count) ||
		!pMesh->meshAttribs.count,
		""
	);
	PIX_ERR_RETURN_IFNOT_COND(
		err,
		(pMesh->faceAttribs.pArr && pMesh->faceAttribs.count) ||
		!pMesh->faceAttribs.count,
		""
	);
	PIX_ERR_RETURN_IFNOT_COND(
		err,
		pMesh->cornerAttribs.pArr && pMesh->cornerAttribs.count,
		""
	);
	PIX_ERR_RETURN_IFNOT_COND(
		err,
		(pMesh->edgeAttribs.pArr && pMesh->edgeAttribs.count) ||
		!pMesh->edgeAttribs.count,
		""
	);
	PIX_ERR_RETURN_IFNOT_COND(
		err,
		pMesh->vertAttribs.pArr && pMesh->vertAttribs.count,
		""
	);
	I8 *pCornerUsed = pAlloc->fpCalloc(pMesh->cornerCount, 1);
	I8 *pVertUsed = pAlloc->fpCalloc(pMesh->vertCount, 1);
	for (I32 i = 0; i < pMesh->faceCount; ++i) {
		FaceRange face = stucGetFaceRange(pMesh, i);
		PIX_ERR_RETURN_IFNOT_COND(err, face.size >= 3, "");
		for (I32 j = 0; j < face.size; ++j) {
			I32 corner = face.start + j;
			PIX_ERR_RETURN_IFNOT_COND(err, corner < pMesh->cornerCount, "");
			PIX_ERR_RETURN_IFNOT_COND(err, pMesh->pCorners[corner] < pMesh->vertCount, "");
			pCornerUsed[corner] = true;
			pVertUsed[pMesh->pCorners[corner]] = true;
		}
	}
	for (I32 i = 0; i < pMesh->cornerCount; ++i) {
		PIX_ERR_THROW_IFNOT_COND(
			err,
			pCornerUsed[i],
			"in mesh contains floating corner(s)",
			0
		);
	}
	for (I32 i = 0; i < pMesh->vertCount; ++i) {
		PIX_ERR_THROW_IFNOT_COND(
			err,
			pVertUsed[i],
			"in mesh contains floating vert(s)",
			0
		);
	}
	PIX_ERR_CATCH(0, err, ;);
	pAlloc->fpFree(pCornerUsed);
	pCornerUsed = NULL;
	pAlloc->fpFree(pVertUsed);
	pVertUsed = NULL;
	PIX_ERR_RETURN_IFNOT(err, "");

	for (I32 i = 1; i < STUC_ATTRIB_USE_ENUM_COUNT; ++i) {
		if (i == STUC_ATTRIB_USE_SP_ENUM_COUNT) {
			continue;
		}
		AttribActive idx = pMesh->activeAttribs[i];
		if (!idx.active) {
			PIX_ERR_RETURN_IFNOT_COND(err,
				!stucIsAttribUseRequired(i) || posOnly && i != STUC_ATTRIB_USE_POS,
				"mesh must have active attribs for pos, uv, normal, and idx"
			);
			continue;
		}
		PIX_ERR_RETURN_IFNOT_COND(err, idx.idx >= 0, "invalid active attrib index");
		PIX_ERR_RETURN_IFNOT_COND(
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
			PIX_ERR_ASSERT("invalid domain", false);
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
			PIX_ERR_ASSERT("invalid domain", false);
	}
	PIX_ERR_ASSERT("invalid domain", false);
	return 0;
}

bool stucGetIfSeamEdge(const Mesh *pMesh, I32 edge) {
	PIX_ERR_ASSERT("", pMesh->pSeamEdge);
	return pMesh->pSeamEdge[edge];
}

bool stucGetIfMatBorderEdge(const Mesh *pMesh, I32 edge) {
	PIX_ERR_ASSERT("", pMesh->pEdgeFaces && pMesh->pMatIdx);
	V2_I32 faces = pMesh->pEdgeFaces[edge];
	if (faces.d[1] == -1) {
		return false; //no adj face
	}
	return pMesh->pMatIdx[faces.d[0]] != pMesh->pMatIdx[faces.d[1]];
}

//does not increment the other corner, so it may not actually be adjacent.
//Do so after calling this func (if wind is equal between faces)
void stucGetAdjCorner(const Mesh *pMesh, FaceCorner corner, FaceCorner *pAdjCorner) {
	PIX_ERR_ASSERT("", pMesh->pEdgeFaces);
	FaceRange face = stucGetFaceRange(&pMesh->core, corner.face);
	I32 edge = pMesh->core.pEdges[face.start + corner.corner];
	V2_I8 corners = pMesh->pEdgeCorners[edge];
	V2_I32 faces = pMesh->pEdgeFaces[edge];
	bool which = corners.d[0] == corner.corner && faces.d[0] == corner.face;
	PIX_ERR_ASSERT(
		"",
		faces.d[!which] == -1 && corners.d[!which] == -1 ||
		faces.d[!which] == corner.face && corners.d[!which] == corner.corner
	);
	pAdjCorner->corner = corners.d[which];
	pAdjCorner->face = faces.d[which];
}

I32 stucGetMeshVert(const StucMesh *pMesh, FaceCorner corner) {
	PIX_ERR_ASSERT("", pMesh && corner.face >= 0 && corner.face < pMesh->faceCount);
	FaceRange face = stucGetFaceRange(pMesh, corner.face);
	PIX_ERR_ASSERT("", corner.corner >= 0 && corner.corner < face.size);
	I32 vert = pMesh->pCorners[face.start + corner.corner];
	PIX_ERR_ASSERT("", vert >= 0 && vert < pMesh->vertCount);
	return vert;
}

I32 stucGetMeshEdge(const StucMesh *pMesh, FaceCorner corner) {
	PIX_ERR_ASSERT("", pMesh && corner.face >= 0 && corner.face < pMesh->faceCount);
	FaceRange face = stucGetFaceRange(pMesh, corner.face);
	PIX_ERR_ASSERT("", corner.corner >= 0 && corner.corner < face.size);
	I32 edge = pMesh->pEdges[face.start + corner.corner];
	PIX_ERR_ASSERT("", edge >= 0 && edge < pMesh->edgeCount);
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

bool stucQuickCmpMesh(StucContext pCtx, const StucMesh *pA, const StucMesh *pB) {
	if (pA->vertCount != pB->vertCount ||
		memcmp(
			pA->activeAttribs,
			pB->activeAttribs,
			sizeof(StucAttribActive) * STUC_ATTRIB_USE_ENUM_COUNT
	)) {
		return false;
	}
	const StucAttrib *pVerts = stucGetActiveAttribConst(pCtx, pA, STUC_ATTRIB_USE_POS);
	const StucAttrib *pRefVerts = stucGetActiveAttribConst(pCtx, pA, STUC_ATTRIB_USE_POS);
	return !(
		memcmp(
			pVerts->core.pData,
			pRefVerts->core.pData,
			stucGetAttribSizeIntern(pVerts->core.type) * pA->vertCount
		) ||
		memcmp(pA->pCorners, pB->pCorners, sizeof(I32) * pA->cornerCount) ||
		memcmp(pA->pFaces, pB->pFaces, sizeof(I32) * (pA->faceCount + 1))
	);
}

bool stucQuickCmpObj(StucContext pCtx, const StucObject *pA, const StucObject *pB) {
	if (memcmp(&pA->transform, &pB->transform, sizeof(Stuc_M4x4))) {
		return false;
	}
	const StucMesh *pAMesh = (StucMesh *)pA->pData;
	const StucMesh *pBMesh = (StucMesh *)pB->pData;
	return stucQuickCmpMesh(pCtx, pAMesh, pBMesh);
}
