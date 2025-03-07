#pragma once

#include <stdint.h>

#include <uv_stucco.h>
#include <debug_and_perf.h>
#include <types.h>

typedef struct {
	I32 idx;
	I32 realIdx;
} BufMeshIdx;

//is it worth just putting pData at the start of Attrib,
//so you can just have pNormalAttrib and remove pNormal.
//or would having to cast all the time be too annoying?
//probably
//TODO maybe combine I8 flags into a single 'StucEdgeFlags' or 'StucVertFlags' attrib?
//Like you could have edge preserve, edge receive and seam edge in a
//single attrib. Is it worth it? They don't take up much atm. If more flags were added
//it'd might be worth considering
typedef struct {
	StucMesh core;
	//aliases to special attribs
	Stuc_V3_F32 *pVerts;
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
	Stuc_V2_I32 *pEdgeCorners;
	I32 faceBufSize;
	I32 cornerBufSize;
	I32 edgeBufSize;
	I32 vertBufSize;
} Mesh;

typedef struct {
	Mesh mesh;
	StucAttrib *pWAttrib;
	StucAttrib *pInNormalAttrib;
	StucAttrib *pInTangentAttrib;
	StucAttrib *pAlphaAttrib;
	StucAttrib *pInTSignAttrib;
	F32 *pW;
	Stuc_V3_F32 *pInNormal;
	Stuc_V3_F32 *pInTangent;
	F32 *pAlpha;
	F32 *pInTSign;
	I32 borderVertCount;
	I32 borderCornerCount;
	I32 borderEdgeCount;
	I32 borderFaceCount;
} BufMesh;

typedef struct {
	I32 faces;
	I32 corners;
	I32 edges;
	I32 verts;
} MeshCounts;

typedef struct {
	I32 start;
	I32 end;
	I32 size;
	I32 idx;
} FaceRange;

void stucCreateMesh(StucContext pCtx, StucObject *pObj, StucObjectType type);
BufMeshIdx stucBufMeshAddFace(
	StucContext pCtx,
	BufMesh *pMesh,
	bool border,
	bool *pRealloced
);
BufMeshIdx stucBufMeshAddCorner(
	StucContext pCtx,
	BufMesh *pMesh,
	bool border,
	bool *pRealloced
);
BufMeshIdx stucBufMeshAddEdge(
	StucContext pCtx,
	BufMesh *pMesh,
	bool border,
	bool *pRealloced
);
BufMeshIdx stucBufMeshAddVert(
	StucContext pCtx,
	BufMesh *pMesh,
	bool border,
	bool *pRealloced
);
BufMeshIdx stucConvertBorderFaceIdx(const BufMesh *pMesh, I32 face);
BufMeshIdx stucConvertBorderCornerIdx(const BufMesh *pMesh, I32 corner);
BufMeshIdx stucConvertBorderEdgeIdx(const BufMesh *pMesh, I32 edge);
BufMeshIdx stucConvertBorderVertIdx(const BufMesh *pMesh, I32 vert);
I32 stucMeshAddFace(StucContext pCtx, Mesh *pMesh, bool *pRealloced);
I32 stucMeshAddCorner(StucContext pCtx, Mesh *pMesh, bool *pRealloced);
I32 stucMeshAddEdge(StucContext pCtx, Mesh *pMesh, bool *pRealloced);
I32 stucMeshAddVert(StucContext pCtx, Mesh *pMesh, bool *pRealloced);
void stucReallocMeshToFit(StucContext pCtx, Mesh *pMesh);
void stucMeshSetLastFace(StucContext pCtx, Mesh *pMesh);
void stucBufMeshSetLastFaces(StucContext pCtx, BufMesh *pBufMesh);
bool stucCheckIfMesh(StucObjectData type);
void stucAddToMeshCounts(
	MeshCounts *pCounts,
	MeshCounts *pBoundsCounts,
	const Mesh *pMeshSrc
);
StucResult stucCopyMesh(StucMesh *pDestMesh, const StucMesh *pSrcMesh);
void stucApplyObjTransform(StucObject *pObj);
StucResult stucMergeObjArr(
	StucContext pCtx,
	Mesh *pMesh,
	I32 objCount,
	const StucObject *pObjArr,
	bool setCommon
);
StucResult stucDestroyObjArr(StucContext pCtx, I32 objCount, StucObject *pObjArr);
FaceRange stucGetFaceRange(const StucMesh *pMesh, I32 idx, bool border);
StucResult stucBuildTangents(Mesh *pMesh);
StucResult stucValidateMesh(const StucMesh *pMesh, bool checkEdges);
void stucAliasMeshCoreNoAttribs(StucMesh *pDest, StucMesh *pSrc);
I32 stucGetDomainSize(const Mesh *pMesh, StucDomain domain);
I32 stucDomainCountGetIntern(const StucMesh *pMesh, StucDomain domain);