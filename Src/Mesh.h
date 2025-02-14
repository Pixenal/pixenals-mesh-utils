#pragma once

#include <stdint.h>

#include <UvStucco.h>
#include <DebugAndPerf.h>

typedef struct {
	int32_t idx;
	int32_t realIdx;
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
	float *pTSigns;
	float *pWScale;
	Stuc_V2_F32 *pUvs;
	int8_t *pEdgePreserve;
	int8_t *pVertPreserve;
	int8_t *pEdgeReceive;
	int32_t *pUsg;
	//it'd be ideal if this were unsigned, do dcc's usually store mats unsigned though?
	//probably not worth it tbh, i've never see an object with 128 mat slots let alone 256
	int8_t *pMatIdx;
	int32_t faceBufSize;
	int32_t cornerBufSize;
	int32_t edgeBufSize;
	int32_t vertBufSize;
} Mesh;

typedef struct {
	Mesh mesh;
	int32_t borderVertCount;
	int32_t borderCornerCount;
	int32_t borderEdgeCount;
	int32_t borderFaceCount;
	StucAttrib *pWAttrib;
	StucAttrib *pInNormalAttrib;
	StucAttrib *pInTangentAttrib;
	StucAttrib *pAlphaAttrib;
	StucAttrib *pInTSignAttrib;
	float *pW;
	Stuc_V3_F32 *pInNormal;
	Stuc_V3_F32 *pInTangent;
	float *pAlpha;
	float *pInTSign;
} BufMesh;

typedef struct {
	int32_t faces;
	int32_t corners;
	int32_t edges;
	int32_t verts;
} MeshCounts;

typedef struct {
	int32_t start;
	int32_t end;
	int32_t size;
	int32_t idx;
} FaceRange;

void stucCreateMesh(StucContext pContext, StucObject *pObj, StucObjectType type);
BufMeshIdx stucBufMeshAddFace(const StucAlloc *pAlloc, BufMesh *pMesh, bool border,
                              DebugAndPerfVars *pDpVars, bool *pRealloced);
BufMeshIdx stucBufMeshAddCorner(const StucAlloc *pAlloc, BufMesh *pMesh, bool border,
                                DebugAndPerfVars *pDpVars, bool *pRealloced);
BufMeshIdx stucBufMeshAddEdge(const StucAlloc *pAlloc, BufMesh *pMesh, bool border,
                              DebugAndPerfVars *pDpVars, bool *pRealloced);
BufMeshIdx stucBufMeshAddVert(const StucAlloc *pAlloc, BufMesh *pMesh, bool border,
                              DebugAndPerfVars *pDpVars, bool *pRealloced);
BufMeshIdx stucConvertBorderFaceIdx(const BufMesh *pMesh, int32_t face);
BufMeshIdx stucConvertBorderCornerIdx(const BufMesh *pMesh, int32_t corner);
BufMeshIdx stucConvertBorderEdgeIdx(const BufMesh *pMesh, int32_t edge);
BufMeshIdx stucConvertBorderVertIdx(const BufMesh *pMesh, int32_t vert);
int32_t stucMeshAddFace(const StucAlloc *pAlloc, Mesh *pMesh, bool *pRealloced);
int32_t stucMeshAddCorner(const StucAlloc *pAlloc, Mesh *pMesh, bool *pRealloced);
int32_t stucMeshAddEdge(const StucAlloc *pAlloc, Mesh *pMesh, bool *pRealloced);
int32_t stucMeshAddVert(const StucAlloc *pAlloc, Mesh *pMesh, bool *pRealloced);
void stucReallocMeshToFit(const StucAlloc *pAlloc, Mesh *pMesh);
void stucMeshSetLastFace(const StucAlloc *pAlloc, Mesh *pMesh);
void stucBufMeshSetLastFaces(const StucAlloc *pAlloc, BufMesh *pBufMesh,
                             DebugAndPerfVars *pDpVars);
bool stucCheckIfMesh(StucObjectData type);
void stucAddToMeshCounts(StucContext pContext, MeshCounts *pCounts,
                         MeshCounts *pBoundsCounts, Mesh *pMeshSrc);
void stucCopyMesh(StucMesh *pDestMesh, StucMesh *pSrcMesh);
void stucApplyObjTransform(StucObject *pObj);
void stucMergeObjArr(StucContext pContext, Mesh *pMesh,
                 int32_t objCount, StucObject *pObjArr, bool setCommon);
StucResult stucDestroyObjArr(StucContext pContext, int32_t objCount, StucObject *pObjArr);
FaceRange stucGetFaceRange(const StucMesh *pMesh, int32_t idx, bool border);
void stucBuildTangents(Mesh *pMesh);
