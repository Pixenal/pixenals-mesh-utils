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
	int32_t boundaryVertSize;
	int32_t boundaryLoopSize;
	int32_t boundaryEdgeSize;
	int32_t boundaryFaceSize;
} BufMesh;

typedef struct {
	RuvmMesh mesh;
	RuvmAttrib *pVertAttrib;
	RuvmAttrib *pUvAttrib;
	RuvmAttrib *pNormalAttrib;
	RuvmAttrib *pEdgePreserveAttrib;
	Ruvm_V3_F32 *pVerts;
	Ruvm_V3_F32 *pNormals;
	Ruvm_V2_F32 *pUvs;
	int8_t *pEdgePreserve;
} Mesh;
