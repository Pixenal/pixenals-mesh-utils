#pragma once

#include <stdint.h>

#include <RUVM.h>
#include <DebugAndPerf.h>

typedef struct {
	int32_t idx;
	int32_t realIdx;
} BufMeshIdx;

typedef struct {
	RuvmMesh mesh;
	RuvmAttrib *pVertAttrib;
	RuvmAttrib *pUvAttrib;
	RuvmAttrib *pNormalAttrib;
	RuvmAttrib *pTangentAttrib;
	RuvmAttrib *pTSignAttrib;
	RuvmAttrib *pEdgePreserveAttrib;
	RuvmAttrib *pVertPreserveAttrib;
	RuvmAttrib *pEdgeReceiveAttrib;
	RuvmAttrib *pUsgAttrib;
	RuvmAttrib *pWScaleAttrib;
	Ruvm_V3_F32 *pVerts;
	Ruvm_V3_F32 *pNormals;
	Ruvm_V3_F32 *pTangents;
	float *pTSigns;
	float *pWScale;
	Ruvm_V2_F32 *pUvs;
	int8_t *pEdgePreserve;
	int8_t *pVertPreserve;
	int8_t *pEdgeReceive;
	int32_t *pUsg;
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
	RuvmAttrib *pWAttrib;
	RuvmAttrib *pInNormalAttrib;
	RuvmAttrib *pInTangentAttrib;
	RuvmAttrib *pAlphaAttrib;
	RuvmAttrib *pInTSignAttrib;
	float *pW;
	Ruvm_V3_F32 *pInNormal;
	Ruvm_V3_F32 *pInTangent;
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

void createMesh(RuvmContext pContext, RuvmObject *pObj, RuvmObjectType type);
BufMeshIdx bufMeshAddFace(const RuvmAlloc *pAlloc, BufMesh *pMesh, bool border,
                            DebugAndPerfVars *pDpVars, bool *pRealloced);
BufMeshIdx bufMeshAddCorner(const RuvmAlloc *pAlloc, BufMesh *pMesh, bool border,
                            DebugAndPerfVars *pDpVars, bool *pRealloced);
BufMeshIdx bufMeshAddEdge(const RuvmAlloc *pAlloc, BufMesh *pMesh, bool border,
                            DebugAndPerfVars *pDpVars, bool *pRealloced);
BufMeshIdx bufMeshAddVert(const RuvmAlloc *pAlloc, BufMesh *pMesh, bool border,
                            DebugAndPerfVars *pDpVars, bool *pRealloced);
BufMeshIdx convertBorderFaceIdx(const BufMesh *pMesh, int32_t face);
BufMeshIdx convertBorderCornerIdx(const BufMesh *pMesh, int32_t corner);
BufMeshIdx convertBorderEdgeIdx(const BufMesh *pMesh, int32_t edge);
BufMeshIdx convertBorderVertIdx(const BufMesh *pMesh, int32_t vert);
int32_t meshAddFace(const RuvmAlloc *pAlloc, Mesh *pMesh, bool *pRealloced);
int32_t meshAddCorner(const RuvmAlloc *pAlloc, Mesh *pMesh, bool *pRealloced);
int32_t meshAddEdge(const RuvmAlloc *pAlloc, Mesh *pMesh, bool *pRealloced);
int32_t meshAddVert(const RuvmAlloc *pAlloc, Mesh *pMesh, bool *pRealloced);
void reallocMeshToFit(const RuvmAlloc *pAlloc, Mesh *pMesh);
void meshSetLastFace(const RuvmAlloc *pAlloc, Mesh *pMesh);
void bufMeshSetLastFaces(const RuvmAlloc *pAlloc, BufMesh *pBufMesh,
                         DebugAndPerfVars *pDpVars);
bool checkIfMesh(RuvmMesh *pMesh);
void addToMeshCounts(RuvmContext pContext, MeshCounts *pCounts,
                     MeshCounts *pBoundsCounts, Mesh *pMeshSrc);
void copyMesh(RuvmMesh *pDestMesh, RuvmMesh *pSrcMesh);
void applyObjTransform(RuvmObject *pObj);
void mergeObjArr(RuvmContext pContext, Mesh *pMesh,
                 int32_t objCount, RuvmObject *pObjArr, bool setCommon);
void destroyObjArr(RuvmContext pContext, int32_t objCount, RuvmObject *pObjArr);
FaceRange getFaceRange(const RuvmMesh *pMesh, int32_t idx, bool border);
void buildTangents(Mesh *pMesh);
//TODO remove this, it's unecessary, just do &pBufMesh->mesh instead.
//     Obviously also rename mesh in Mesh to core, so it would be
//     pBufMesh->mesh.core.faceCount, for instance.
static inline
Mesh *asMesh(BufMesh *pMesh) {
	return (Mesh *)pMesh;
}
