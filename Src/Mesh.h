#pragma once

#include <stdint.h>

#include <RUVM.h>

typedef RuvmMesh Mesh;

typedef struct {
	int32_t vertCount;
	int32_t boundaryVertSize;
	RuvmVec3 *pVerts;
	int32_t loopCount;
	int32_t boundaryLoopSize;
	int32_t *pLoops;
	RuvmVec3 *pNormals;
	int32_t faceCount;
	int32_t boundaryFaceSize;
	int32_t *pFaces;
	RuvmVec2 *pUvs;
} WorkMesh;
