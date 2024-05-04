#pragma once

#include <Types.h>

int32_t checkFaceIsInBounds(V2_F32 min, V2_F32 max, FaceInfo face, Mesh *pMesh);
void getFaceBounds(FaceBounds *pBounds, V2_F32 *pUvs, FaceInfo faceInfo);
int32_t checkIfEdgeIsSeam(int32_t edgeIndex, FaceInfo face, int32_t loop,
                          Mesh *pMesh, EdgeVerts *pEdgeVerts);

uint32_t ruvmFnvHash(uint8_t *value, int32_t valueSize, uint32_t size);

int32_t checkIfEdgeIsPreserve(Mesh* pMesh, int32_t edge);
FaceTriangulated triangulateFace(RuvmAllocator alloc, FaceInfo baseFace, Mesh *pMesh);
