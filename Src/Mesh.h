#pragma once

#include <stdint.h>

#include <RUVM.h>

typedef struct {
	RuvmMesh mesh;
	RuvmAttrib *pVerts;
	RuvmAttrib *pUvs;
	RuvmAttrib *pNormals;
	int32_t boundaryVertSize;
	int32_t boundaryLoopSize;
	int32_t boundaryEdgeSize;
	int32_t boundaryFaceSize;
} BufMesh;

typedef struct {
	RuvmMesh mesh;
	RuvmAttrib *pVerts;
	RuvmAttrib *pUvs;
	RuvmAttrib *pNormals;
	RuvmAttrib *pEdgePreserve;
} Mesh;
