/* 
SPDX-FileCopyrightText: 2025 Caleb Dawson
SPDX-License-Identifier: Apache-2.0
*/

#pragma once

#include <stdint.h>

#include <pixenals_alloc_utils.h>
#include <pixenals_math_utils.h>

#include <uv_stucco.h>
#include <types.h>

//is it worth just putting pData at the start of Attrib,
//so you can just have pNormalAttrib and remove pNormal.
//or would having to cast all the time be too annoying?
//probably
//TODO maybe combine I8 flags into a single 'StucEdgeFlags' or 'StucVertFlags' attrib?
//Like you could have edge preserve, edge receive and seam edge in a
//single attrib. Is it worth it? They don't take up much atm. If more flags were added
//it'd might be worth considering
typedef struct Mesh {
	StucMesh core;
	//aliases to special attribs
	Stuc_V3_F32 *pPos;
	Stuc_V3_F32 *pNormals;
	Stuc_V3_F32 *pTangents;
	F32 *pTSigns;
	F32 *pWScale;
	Stuc_V2_F32 *pUvs;
	I8 *pEdgePreserve;
	I8 *pVertPreserve;
	I8 *pEdgeReceive;
	I32 *pUsg;
	//it'd be ideal if this were unsigned, do dcc's usually store mats unsigned though?
	//probably not worth it tbh, i've never see an object with 128 mat slots let alone 256
	I8 *pMatIdx;
	F32 *pEdgeLen;
	I8 *pSeamEdge;
	I8 *pSeamVert;
	I8 *pNumAdjPreserve;
	Stuc_V2_I32 *pEdgeFaces;
	Stuc_V2_I8 *pEdgeCorners;
	Stuc_V3_F32 *pVertNormals;
	I32 faceBufSize;
	I32 cornerBufSize;
	I32 edgeBufSize;
	I32 vertBufSize;
} Mesh;

typedef struct MeshCounts {
	I32 faces;
	I32 corners;
	I32 edges;
	I32 verts;
} MeshCounts;

typedef struct FaceRange {
	I32 start;
	I32 end;
	I32 size;
	I32 idx;
} FaceRange;

typedef struct FaceCorner {
	I32 face;
	I32 corner;
} FaceCorner;

typedef struct FaceTriangulated {
	U8 *pTris;
} FaceTriangulated;

typedef struct TriCache {
	FaceTriangulated *pArr;
	PixalcLinAlloc alloc;
} TriCache;

void stucCreateMesh(const StucContext pCtx, StucObject *pObj, StucObjectType type);
I32 stucMeshAddFace(const StucContext pCtx, Mesh *pMesh, bool *pRealloced);
I32 stucMeshAddCorner(const StucContext pCtx, Mesh *pMesh, bool *pRealloced);
I32 stucMeshAddEdge(const StucContext pCtx, Mesh *pMesh, bool *pRealloced);
I32 stucMeshAddVert(const StucContext pCtx, Mesh *pMesh, bool *pRealloced);
void stucReallocMeshToFit(const StucContext pCtx, Mesh *pMesh);
void stucMeshSetLastFace(const StucContext pCtx, Mesh *pMesh);
bool stucCheckIfMesh(const StucObjectData type);
void stucAddToMeshCounts(
	MeshCounts *pCounts,
	MeshCounts *pBoundsCounts,
	const Mesh *pMeshSrc
);
StucErr stucCopyMesh(StucMesh *pDestMesh, const StucMesh *pSrcMesh);
void stucApplyObjTransform(StucObject *pObj);
StucErr stucMergeObjArr(
	StucContext pCtx,
	Mesh *pMesh,
	const StucObjArr *pObjArr,
	bool setCommon
);
static inline
FaceRange stucGetFaceRange(const StucMesh *pMesh, I32 idx) {
	PIX_ERR_ASSERT("", idx >= 0 && idx < pMesh->faceCount);
	FaceRange face = {
		.idx = idx,
		.start = pMesh->pFaces[idx],
		.end = pMesh->pFaces[idx + 1]
	};
	PIX_ERR_ASSERT("", face.start >= 0 && face.end <= pMesh->cornerCount);
	PIX_ERR_ASSERT("", face.start < face.end);
	face.size = face.end - face.start;
	PIX_ERR_ASSERT("", face.size >= 3);
	return face;
}
StucErr stucValidateMesh(const StucMesh *pMesh, bool checkEdges, bool posOnly);
void stucAliasMeshCoreNoAttribs(StucMesh *pDest, StucMesh *pSrc);
I32 stucGetDomainSize(const Mesh *pMesh, StucDomain domain);
I32 stucDomainCountGetIntern(const StucMesh *pMesh, StucDomain domain);
static inline
I32 stucGetCornerPrev(I32 corner, const FaceRange *pFace) {
	I32 prev = corner ? corner - 1 : pFace->size - 1;
	PIX_ERR_ASSERT("", prev >= 0 && prev < pFace->size);
	return prev;
}
static inline
I32 stucGetCornerNext(I32 corner, const FaceRange *pFace) {
	PIX_ERR_ASSERT("", corner < pFace->size);
	I32 next = (corner + 1) % pFace->size;
	PIX_ERR_ASSERT("", next >= 0);
	return next;
}
bool stucGetIfSeamEdge(const Mesh *pMesh, I32 edge);
bool stucGetIfMatBorderEdge(const Mesh *pMesh, I32 edge);
void stucGetAdjCorner(const Mesh *pMesh, FaceCorner corner, FaceCorner *pAdjCorner);
static inline
V3_F32 stucGetVertPos(const Mesh *pMesh, const FaceRange *pFace, I32 corner) {
	PIX_ERR_ASSERT("", corner >= 0 && corner < pFace->size);
	return pMesh->pPos[pMesh->core.pCorners[pFace->start + corner]];
}
static inline
V2_F32 stucGetVertPosAsV2(const Mesh *pMesh, const FaceRange *pFace, I32 corner) {
	PIX_ERR_ASSERT("", corner >= 0 && corner < pFace->size);
	return *(V2_F32 *)&pMesh->pPos[pMesh->core.pCorners[pFace->start + corner]];
}
static inline
V2_F32 stucGetUvPos(const Mesh *pMesh, const FaceRange *pFace, I32 corner) {
	PIX_ERR_ASSERT("", corner >= 0 && corner < pFace->size);
	return pMesh->pUvs[pFace->start + corner];
}
I32 stucGetMeshVert(const StucMesh *pMesh, FaceCorner corner);
I32 stucGetMeshEdge(const StucMesh *pMesh, FaceCorner corner);
bool checkForNgonsInMesh(const StucMesh *pMesh);
bool stucQuickCmpMesh(StucContext pCtx, const StucMesh *pA, const StucMesh *pB);
bool stucQuickCmpObj(StucContext pCtx, const StucObject *pA, const StucObject *pB);
