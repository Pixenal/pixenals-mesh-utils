#pragma once

#include <stdint.h>

#include <RUVM.h>

typedef struct {
	RuvmMesh mesh;
	RuvmAttrib *pVertAttrib;
	RuvmAttrib *pUvAttrib;
	RuvmAttrib *pNormalAttrib;
	RuvmAttrib *pEdgePreserveAttrib;
	Ruvm_V3_F32 *pVerts;
	Ruvm_V3_F32 *pNormals;
	Ruvm_V2_F32 *pUvs;
	int32_t borderVertCount;
	int32_t borderLoopCount;
	int32_t borderEdgeCount;
	int32_t borderFaceCount;
	int32_t faceBufSize;
	int32_t loopBufSize;
	int32_t edgeBufSize;
	int32_t vertBufSize;
} BufMesh;

typedef struct {
	RuvmMesh mesh;
	RuvmAttrib *pVertAttrib;
	RuvmAttrib *pUvAttrib;
	RuvmAttrib *pNormalAttrib;
	RuvmAttrib *pEdgePreserveAttrib;
	RuvmAttrib *pEdgeReceiveAttrib;
	Ruvm_V3_F32 *pVerts;
	Ruvm_V3_F32 *pNormals;
	Ruvm_V2_F32 *pUvs;
	int8_t *pEdgePreserve;
	int8_t *pEdgeReceive;
} Mesh;
