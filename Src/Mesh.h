#pragma once

#include <stdint.h>

#include <UvStucco.h>
#include <DebugAndPerf.h>
#include <Types.h>

typedef struct {
	I32 idx;
	I32 realIdx;
} BufMeshIdx;

//is it worth just putting pData at the start of Attrib,
//so you can just have pNormalAttrib and remove pNormal.
//or would having to cast all the time be too annoying?
//probably
typedef struct {
	StucMesh core;
	StucAttrib *pVertAttrib;
	StucAttrib *pUvAttrib;
	StucAttrib *pNormalAttrib;
	StucAttrib *pTangentAttrib;
	StucAttrib *pTSignAttrib;
	StucAttrib *pEdgePreserveAttrib;
	StucAttrib *pVertPreserveAttrib;
	StucAttrib *pEdgeReceiveAttrib;
	StucAttrib *pUsgAttrib;
	StucAttrib *pWScaleAttrib;
	StucAttrib *pMatIdxAttrib;
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
	I32 faceBufSize;
	I32 cornerBufSize;
	I32 edgeBufSize;
	I32 vertBufSize;
} Mesh;

typedef struct {
	Mesh mesh;
	I32 borderVertCount;
	I32 borderCornerCount;
	I32 borderEdgeCount;
	I32 borderFaceCount;
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
	const StucAlloc *pAlloc,
	BufMesh *pMesh,
	bool border,
	DebugAndPerfVars *pDpVars,
	bool *pRealloced
);
BufMeshIdx stucBufMeshAddCorner(
	const StucAlloc *pAlloc,
	BufMesh *pMesh,
	bool border,
	DebugAndPerfVars *pDpVars,
	bool *pRealloced
);
BufMeshIdx stucBufMeshAddEdge(
	const StucAlloc *pAlloc,
	BufMesh *pMesh,
	bool border,
	DebugAndPerfVars *pDpVars,
	bool *pRealloced
);
BufMeshIdx stucBufMeshAddVert(
	const StucAlloc *pAlloc,
	BufMesh *pMesh,
	bool border,
	DebugAndPerfVars *pDpVars,
	bool *pRealloced
);
BufMeshIdx stucConvertBorderFaceIdx(const BufMesh *pMesh, I32 face);
BufMeshIdx stucConvertBorderCornerIdx(const BufMesh *pMesh, I32 corner);
BufMeshIdx stucConvertBorderEdgeIdx(const BufMesh *pMesh, I32 edge);
BufMeshIdx stucConvertBorderVertIdx(const BufMesh *pMesh, I32 vert);
I32 stucMeshAddFace(const StucAlloc *pAlloc, Mesh *pMesh, bool *pRealloced);
I32 stucMeshAddCorner(const StucAlloc *pAlloc, Mesh *pMesh, bool *pRealloced);
I32 stucMeshAddEdge(const StucAlloc *pAlloc, Mesh *pMesh, bool *pRealloced);
I32 stucMeshAddVert(const StucAlloc *pAlloc, Mesh *pMesh, bool *pRealloced);
void stucReallocMeshToFit(const StucAlloc *pAlloc, Mesh *pMesh);
void stucMeshSetLastFace(const StucAlloc *pAlloc, Mesh *pMesh);
void stucBufMeshSetLastFaces(
	const StucAlloc *pAlloc,
	BufMesh *pBufMesh,
	DebugAndPerfVars *pDpVars
);
bool stucCheckIfMesh(StucObjectData type);
void stucAddToMeshCounts(
	MeshCounts *pCounts,
	MeshCounts *pBoundsCounts,
	const Mesh *pMeshSrc
);
StucResult stucCopyMesh(StucMesh *pDestMesh, const StucMesh *pSrcMesh);
void stucApplyObjTransform(StucObject *pObj);
void stucMergeObjArr(
	StucContext pCtx,
	Mesh *pMesh,
	I32 objCount,
	const StucObject *pObjArr,
	bool setCommon
);
StucResult stucDestroyObjArr(StucContext pCtx, I32 objCount, StucObject *pObjArr);
FaceRange stucGetFaceRange(const StucMesh *pMesh, I32 idx, bool border);
StucResult stucBuildTangents(Mesh *pMesh);
StucResult stucValidateMesh(StucMesh *pMesh, bool checkEdges);