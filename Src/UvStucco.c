#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include <Io.h>
#include <MapToJobMesh.h>
#include <CombineJobMeshes.h>
#include <MapFile.h>
#include <Context.h>
#include <Alloc.h>
#include <ThreadPool.h>
#include <UvStucco.h>
#include <Clock.h>
#include <AttribUtils.h>
#include <Utils.h>
#include <ImageUtils.h>
#include <Error.h>

// TODO
// - Reduce the bits written to the UVGP file for vert and corner indices, based on the total amount, in order to save space.
//   No point storing them as 32 bit if there's only like 4,000 verts
// - Split compressed data into chunks maybe?
//
// - Add blending options to interface, that control how MeshIn attributes blend with
//   those from the Map. Also add an option to disable or enable interpolation.
//   Add these to the StucAttrib struct.
//
//TODO a highly distorted meshIn can cause invalid geometry
//(enough to crash blender). When meshIn is quads atleast
//(I've not tested with tris). Find out why
//TODO add option to vary z projection depth with uv stretch (for wires and such)
//TODO Add an option for subdivision like smoothing (for instances where
//the map is higher res than the base mesh). So that the surface can be
//smoothed, without needing to induce the perf cost of actually subdividing
//the base mesh. Is this possible?
//TODO add user define void * args to custom callbacks
//TODO allow for layered mapping. eg, map-faces assigned layer 0 are only mapped
//to in-faces with a layer attribute of 0

static
void setDefaultStageReport(StucContext pContext) {
	pContext->stageReport.outOf = 50,
	pContext->stageReport.pBegin = stageBegin;
	pContext->stageReport.pProgress = stageProgress;
	pContext->stageReport.pEnd = stageEnd;
}

StucResult stucContextInit(StucContext *pContext, StucAlloc *pAlloc,
                           StucThreadPool *pThreadPool, StucIo *pIo,
					       StucTypeDefaultConfig *pTypeDefaultConfig,
                           StucStageReport *pStageReport) {
	StucAlloc alloc;
	if (pAlloc) {
		stucAllocSetCustom(&alloc, pAlloc);
	}
	else {
		stucAllocSetDefault(&alloc);
	}
	*pContext = alloc.pCalloc(1, sizeof(StucContextInternal));
	(*pContext)->alloc = alloc;
	if (pThreadPool) {
		stucThreadPoolSetCustom(*pContext, pThreadPool);
	}
	else {
		stucThreadPoolSetDefault(*pContext);
	}
	if (pIo) {
		stucIoSetCustom(*pContext, pIo);
	}
	else {
		stucIoSetDefault(*pContext);
	}
	(*pContext)->threadPool.pInit(&(*pContext)->pThreadPoolHandle,
	                              &(*pContext)->threadCount,
	                              &(*pContext)->alloc);
	if (pTypeDefaultConfig) {
		(*pContext)->typeDefaults = *pTypeDefaultConfig;
	}
	else {
		stucSetTypeDefaultConfig(*pContext);
	}
	if (pStageReport) {
		(*pContext)->stageReport = *pStageReport;
	}
	else {
		setDefaultStageReport(*pContext);
	}
	//TODO add ability to set custom specialAttrib names
	setDefaultSpecialAttribNames(*pContext);
	return STUC_SUCCESS;
}

StucResult stucContextDestroy(StucContext pContext) {
	pContext->threadPool.pDestroy(pContext->pThreadPoolHandle);
	pContext->alloc.pFree(pContext);
	return STUC_SUCCESS;
}

StucResult stucMapFileExport(StucContext pContext, char *pName,
                             int32_t objCount, StucObject* pObjArr,
                             int32_t usgCount, StucUsg* pUsgArr,
                             StucAttribIndexedArr *pIndexedAttribs) {
	return stucWriteStucFile(pContext, pName, objCount, pObjArr,
	                         usgCount, pUsgArr, pIndexedAttribs);
}

//TODO replace these with StucUsg and StucObj arr structs, that combine arr and count
StucResult stucMapFileLoadForEdit(StucContext pContext, char *filePath,
                                  int32_t *pObjCount, StucObject **ppObjArr,
                                  int32_t *pUsgCount, StucUsg **ppUsgArr,
                                  int32_t *pFlatCutoffCount, StucObject **ppFlatCutoffArr,
                                  StucAttribIndexedArr *pIndexedAttribs) {
	return stucLoadStucFile(pContext, filePath, pObjCount, ppObjArr, pUsgCount,
	                        ppUsgArr, pFlatCutoffCount, ppFlatCutoffArr, true, pIndexedAttribs);
}

StucResult stucMapFileLoad(StucContext pContext, StucMap *pMapHandle,
                           char *filePath) {
	StucResult err = STUC_NOT_SET;
	StucMap pMap = pContext->alloc.pCalloc(1, sizeof(MapFile));
	int32_t objCount = 0;
	StucObject *pObjArr = NULL;
	StucUsg *pUsgArr = NULL;
	int32_t flatCutoffCount = 0;
	StucObject *pFlatCutoffArr = NULL;
	err = stucLoadStucFile(pContext, filePath, &objCount, &pObjArr,
	                          &pMap->usgArr.count, &pUsgArr, &flatCutoffCount,
	                          &pFlatCutoffArr, false, &pMap->indexedAttribs);
	//TODO validate meshes, ensure pMatIdx is within mat range, faces are within max corner limit,
	//float values are valid, etc.
	STUC_THROW_IF(err, true, "failed to load file from disk", 0);

	for (int32_t i = 0; i < objCount; ++i) {
		setSpecialAttribs(pContext, (Mesh *)pObjArr[i].pData, 0x8ae); //10101110 - all except for preserve
		applyObjTransform(pObjArr + i);
	}
	pMap->mesh.core.type.type = STUC_OBJECT_DATA_MESH_INTERN;
	mergeObjArr(pContext, &pMap->mesh, objCount, pObjArr, false);

	setAttribOrigins(&pMap->mesh.core.meshAttribs, STUC_ATTRIB_ORIGIN_MAP);
	setAttribOrigins(&pMap->mesh.core.faceAttribs, STUC_ATTRIB_ORIGIN_MAP);
	setAttribOrigins(&pMap->mesh.core.cornerAttribs, STUC_ATTRIB_ORIGIN_MAP);
	setAttribOrigins(&pMap->mesh.core.edgeAttribs, STUC_ATTRIB_ORIGIN_MAP);
	setAttribOrigins(&pMap->mesh.core.vertAttribs, STUC_ATTRIB_ORIGIN_MAP);
	setAttribToDontCopy(pContext, &pMap->mesh, 0x7f0);

	setSpecialAttribs(pContext, &pMap->mesh, 0x8ae);

	//TODO some form of heap corruption when many objects
	//test with address sanitizer on CircuitPieces.stuc
	destroyObjArr(pContext, objCount, pObjArr);

	if (pMap->mesh.pUvAttrib) {
		//TODO as with all special attributes, allow user to define what should be considered
		//     the primary UV channel. This especially important for integration with other DCCs
		if (!strncmp(pMap->mesh.pUvAttrib->core.name, "UVMap", STUC_ATTRIB_NAME_MAX_LEN)) {
			char newName[STUC_ATTRIB_NAME_MAX_LEN] = "Map_UVMap";
			memcpy(pMap->mesh.pUvAttrib->core.name, newName, STUC_ATTRIB_NAME_MAX_LEN);
		}
	}

	//set corner attribs to interpolate by default
	//TODO make this an option in ui, even for non common attribs
	for (int32_t i = 0; i < pMap->mesh.core.cornerAttribs.count; ++i) {
		pMap->mesh.core.cornerAttribs.pArr[i].interpolate = true;
	}

	//the quadtree is created before USGs are assigned to verts,
	//as the tree's used to speed up the process
	printf("File loaded. Creating quad tree\n");
	err = stucCreateQuadTree(pContext, pMap);
	STUC_THROW_IF(err, true, "failed to create quadtree", 0);

	if (pMap->usgArr.count) {
		pMap->usgArr.pArr = pContext->alloc.pCalloc(pMap->usgArr.count, sizeof(Usg));
		for (int32_t i = 0; i < pMap->usgArr.count; ++i) {
			setSpecialAttribs(pContext, (Mesh *)pUsgArr[i].obj.pData, 0x02); //000010 - set only vert pos
			pMap->usgArr.pArr[i].origin = *(V2_F32 *)&pUsgArr[i].obj.transform.d[3];
			pMap->usgArr.pArr[i].pMesh = (Mesh *)pUsgArr[i].obj.pData;
			applyObjTransform(&pUsgArr[i].obj);
			if (pUsgArr[i].pFlatCutoff) {
				pMap->usgArr.pArr[i].pFlatCutoff = (Mesh *)pUsgArr[i].pFlatCutoff->pData;
				setSpecialAttribs(pContext, (Mesh *)pUsgArr[i].pFlatCutoff->pData, 0x02); //000010 - set only vert pos
				applyObjTransform(pUsgArr[i].pFlatCutoff);
			}
		}
		//TODO remove duplicate uses of alloc where pContext is present
		//like this
		allocUsgSquaresMesh(pContext, &pContext->alloc, pMap);
		fillUsgSquaresMesh(pMap, pUsgArr);
		assignUsgsToVerts(&pContext->alloc, pMap, pUsgArr);
		pMap->usgArr.pMemArr = pUsgArr;
	}

	*pMapHandle = pMap;
	//TODO add proper checks, and return STUC_ERROR if fails.
	//Do for all public functions (or internal ones as well)
	STUC_CATCH(0, err, stucMapFileUnload(pContext, pMap);)
	return err;
}

StucResult stucMapFileUnload(StucContext pContext, StucMap pMap) {
	stucDestroyQuadTree(pContext, &pMap->quadTree);
	stucMeshDestroy(pContext, &pMap->mesh.core);
	pContext->alloc.pFree(pMap);
	return STUC_SUCCESS;
}

static
void initCommonAttrib(StucContext pContext, StucCommonAttrib *pEntry,
                      StucAttrib *pAttrib) {
	memcpy(pEntry->name, pAttrib->core.name, STUC_ATTRIB_NAME_MAX_LEN);
	StucTypeDefault *pDefault = 
		getTypeDefaultConfig(&pContext->typeDefaults, pAttrib->core.type);
	pEntry->blendConfig = pDefault->blendConfig;
}

static
void getCommonAttribs(StucContext pContext, AttribArray *pMapAttribs,
					  AttribArray *pMeshAttribs,
					  int32_t *pCommonAttribCount,
					  StucCommonAttrib **ppCommonAttribs) {
	//TODO ignore special attribs like StucTangent or StucTSign
	if (!pMeshAttribs || !pMapAttribs) {
		return;
	}
	int32_t count = 0;
	for (int32_t i = 0; i < pMeshAttribs->count; ++i) {
		for (int32_t j = 0; j < pMapAttribs->count; ++j) {
			if (!strncmp(pMeshAttribs->pArr[i].core.name,
			             pMapAttribs->pArr[j].core.name,
			             STUC_ATTRIB_NAME_MAX_LEN)) {
				count++;
			}
		}
	}
	*ppCommonAttribs = count ?
		pContext->alloc.pMalloc(sizeof(StucCommonAttrib) * count) : NULL;
	count = 0;
	for (int32_t i = 0; i < pMeshAttribs->count; ++i) {
		for (int32_t j = 0; j < pMapAttribs->count; ++j) {
			if (!strncmp(pMeshAttribs->pArr[i].core.name,
			             pMapAttribs->pArr[j].core.name,
			             STUC_ATTRIB_NAME_MAX_LEN)) {
				initCommonAttrib(pContext, *ppCommonAttribs + count,
				                 pMeshAttribs->pArr + i);
				count++;
			}
		}
	}
	*pCommonAttribCount = count;
}

//TODO handle edge case, where attribute share the same name,
//but have incompatible types. Such as a float and a string.
StucResult stucQueryCommonAttribs(StucContext pContext, StucMap pMap, StucMesh *pMesh,
                            StucCommonAttribList *pCommonAttribs) {
	getCommonAttribs(pContext, &pMap->mesh.core.meshAttribs, &pMesh->meshAttribs,
					 &pCommonAttribs->meshCount, &pCommonAttribs->pMesh);
	getCommonAttribs(pContext, &pMap->mesh.core.faceAttribs, &pMesh->faceAttribs,
					 &pCommonAttribs->faceCount, &pCommonAttribs->pFace);
	getCommonAttribs(pContext, &pMap->mesh.core.cornerAttribs, &pMesh->cornerAttribs,
					 &pCommonAttribs->cornerCount, &pCommonAttribs->pCorner);
	getCommonAttribs(pContext, &pMap->mesh.core.edgeAttribs, &pMesh->edgeAttribs,
	                 &pCommonAttribs->edgeCount, &pCommonAttribs->pEdge);
	getCommonAttribs(pContext, &pMap->mesh.core.vertAttribs, &pMesh->vertAttribs,
					 &pCommonAttribs->vertCount, &pCommonAttribs->pVert);
	return STUC_SUCCESS;
}

StucResult stucDestroyCommonAttribs(StucContext pContext,
                              StucCommonAttribList *pCommonAttribs) {
	if (pCommonAttribs->pMesh) {
		pContext->alloc.pFree(pCommonAttribs->pMesh);
	}
	if (pCommonAttribs->pFace) {
		pContext->alloc.pFree(pCommonAttribs->pFace);
	}
	if (pCommonAttribs->pCorner) {
		pContext->alloc.pFree(pCommonAttribs->pCorner);
	}
	if (pCommonAttribs->pEdge) {
		pContext->alloc.pFree(pCommonAttribs->pEdge);
	}
	if (pCommonAttribs->pVert) {
		pContext->alloc.pFree(pCommonAttribs->pVert);
	}
	return STUC_SUCCESS;
}

static
void sendOffJobs(StucContext pContext, StucMap pMap, int32_t *pJobCount,
                 void ***pppJobHandles, SendOffArgs **ppJobArgs, Mesh *pMesh,
                 EdgeVerts *pEdgeVerts, int8_t *pInVertTable,
                 StucCommonAttribList *pCommonAttribList, bool getInFaces,
                 float wScale, int8_t maskIdx) {
	//struct timeval start, stop;
	//CLOCK_START;
	*pJobCount = MAX_SUB_MAPPING_JOBS;
	*pJobCount += *pJobCount == 0;
	int32_t facesPerThread = pMesh->core.faceCount / *pJobCount;
	bool singleThread = !facesPerThread;
	*pJobCount = singleThread ? 1 : *pJobCount;
	void *jobArgPtrs[MAX_THREADS];
	int32_t borderTableSize = pMap->mesh.core.faceCount / 5 + 2; //+ 2 incase is 0
	*ppJobArgs = pContext->alloc.pCalloc(*pJobCount, sizeof(SendOffArgs));
	printf("fromjobsendoff: BorderTableSize: %d\n", borderTableSize);
	for (int32_t i = 0; i < *pJobCount; ++i) {
		int32_t meshStart = facesPerThread * i;
		int32_t meshEnd = i == *pJobCount - 1 ?
			pMesh->core.faceCount : meshStart + facesPerThread;
		(*ppJobArgs)[i].inFaceOffset = meshStart;
		(*ppJobArgs)[i].pInVertTable = pInVertTable;
		(*ppJobArgs)[i].pEdgeVerts = pEdgeVerts;
		(*ppJobArgs)[i].pMap = pMap;
		(*ppJobArgs)[i].borderTable.size = borderTableSize;
		(*ppJobArgs)[i].mesh = *pMesh;
		(*ppJobArgs)[i].inFaceRange.start = meshStart;
		(*ppJobArgs)[i].inFaceRange.end = meshEnd;
		(*ppJobArgs)[i].pActiveJobs = pJobCount;
		(*ppJobArgs)[i].id = i;
		(*ppJobArgs)[i].pContext = pContext;
		(*ppJobArgs)[i].pCommonAttribList = pCommonAttribList;
		(*ppJobArgs)[i].getInFaces = getInFaces;
		(*ppJobArgs)[i].wScale = wScale;
		(*ppJobArgs)[i].maskIdx = maskIdx;
		jobArgPtrs[i] = *ppJobArgs + i;
	}
	*pppJobHandles = pContext->alloc.pCalloc(*pJobCount, sizeof(void *));
	pContext->threadPool.pJobStackPushJobs(pContext->pThreadPoolHandle, *pJobCount,
	                                       *pppJobHandles, stucMapToJobMesh, jobArgPtrs);
	//CLOCK_STOP("send off jobs");
}

static
void buildEdgeVertsTable(StucContext pContext, EdgeVerts **ppEdgeVerts,
                         StucMesh *pMesh) {
	*ppEdgeVerts = pContext->alloc.pMalloc(sizeof(EdgeVerts) * pMesh->edgeCount);
	memset(*ppEdgeVerts, -1, sizeof(EdgeVerts) * pMesh->edgeCount);
	for (int32_t i = 0; i < pMesh->cornerCount; ++i) {
		int32_t edge = pMesh->pEdges[i];
		int32_t whichVert = (*ppEdgeVerts)[edge].verts[0] >= 0;
		(*ppEdgeVerts)[edge].verts[whichVert] = i;
	}
}

typedef struct {
	int32_t d[2];
	int32_t idx;
} EdgeCache;

static
void addVertToTableEntry(Mesh *pMesh, FaceRange face, int32_t localCorner,
                         int32_t vert, int32_t edge, EdgeVerts *pEdgeVerts,
						 int8_t *pVertSeamTable, int8_t *pInVertTable,
						 EdgeCache *pEdgeCache, bool *pEdgeSeamTable) {
	//isSeam returns 2 if mesh border, and 1 if uv seam
	int32_t isSeam = checkIfEdgeIsSeam(edge, face, localCorner, pMesh,
									   pEdgeVerts);
	if (isSeam) {
		pVertSeamTable[vert] = isSeam;
		pEdgeSeamTable[edge] = true;
	}
	if (pInVertTable[vert] < 3 &&
		checkIfEdgeIsPreserve(pMesh, edge) &&
		pEdgeCache[vert].d[0] != edge + 1 &&
		pEdgeCache[vert].d[1] != edge + 1) {
		pInVertTable[vert]++;
		int32_t *pEdgeCacheIdx = &pEdgeCache[vert].idx;
		pEdgeCache[vert].d[*pEdgeCacheIdx] = edge + 1;
		++*pEdgeCacheIdx;
	}
}

static
void buildVertTables(StucContext pContext, Mesh *pMesh,
					 int8_t **ppInVertTable, int8_t **ppVertSeamTable,
					 EdgeVerts *pEdgeVerts, bool **ppEdgeSeamTable) {
	*ppInVertTable = pContext->alloc.pCalloc(pMesh->core.vertCount, 1);
	*ppVertSeamTable = pContext->alloc.pCalloc(pMesh->core.vertCount, 1);
	*ppEdgeSeamTable = pContext->alloc.pCalloc(pMesh->core.edgeCount, 1);
	EdgeCache *pEdgeCache =
		pContext->alloc.pCalloc(pMesh->core.vertCount, sizeof(EdgeCache));
	for (int32_t i = 0; i < pMesh->core.faceCount; ++i) {
		FaceRange face = {0};
		face.start = pMesh->core.pFaces[i];
		face.end = pMesh->core.pFaces[i + 1];
		face.size = face.end - face.start;
		face.idx = i;
		for (int32_t j = 0; j < face.size; ++j) {
			int32_t corner = face.start + j;
			int32_t vert = pMesh->core.pCorners[corner];
			int32_t edge = pMesh->core.pEdges[corner];
			addVertToTableEntry(pMesh, face, j, vert, edge, pEdgeVerts,
			                    *ppVertSeamTable, *ppInVertTable, pEdgeCache,
			                    *ppEdgeSeamTable);
			int32_t prevj = j == 0 ? face.size - 1 : j - 1;
			int32_t prevEdge = pMesh->core.pEdges[face.start + prevj];
			addVertToTableEntry(pMesh, face, prevj, vert, prevEdge, pEdgeVerts,
			                    *ppVertSeamTable, *ppInVertTable, pEdgeCache,
			                    *ppEdgeSeamTable);
		}
	}
	pContext->alloc.pFree(pEdgeCache);
}

static
bool checkIfNoFacesHaveMaskIdx(Mesh *pMesh, int8_t maskIdx) {
	if (!pMesh->pMatIdx) {
		return false;
	}
	for (int32_t i = 0; i < pMesh->core.faceCount; ++i) {
		if (pMesh->pMatIdx[i] == maskIdx) {
			return false;
		}
	}
	return true;
}

static
Result mapToMeshInternal(StucContext pContext, StucMap pMap, Mesh *pMeshIn,
                         StucMesh *pMeshOut, int8_t maskIdx,
                         StucCommonAttribList *pCommonAttribList,
                         InFaceArr **ppInFaceTable, float wScale) {
	StucResult err = STUC_SUCCESS;
	if (checkIfNoFacesHaveMaskIdx(pMeshIn, maskIdx)) {
		return STUC_SUCCESS;
	}
	CLOCK_INIT;
	CLOCK_START;
	EdgeVerts *pEdgeVerts = {0};
	printf("EdgeCount: %d\n", pMeshIn->core.edgeCount);
	buildEdgeVertsTable(pContext, &pEdgeVerts, &pMeshIn->core);
	int8_t *pInVertTable;
	int8_t *pVertSeamTable;
	bool *pEdgeSeamTable;
	buildVertTables(pContext, pMeshIn, &pInVertTable,
	                &pVertSeamTable, pEdgeVerts, &pEdgeSeamTable);
	CLOCK_STOP("Edge Table Time");

	CLOCK_START;
	int32_t jobCount = 0;
	SendOffArgs *pJobArgs = NULL;
	void **ppJobHandles = NULL;
	sendOffJobs(pContext, pMap, &jobCount, &ppJobHandles, &pJobArgs, pMeshIn, pEdgeVerts,
	            pInVertTable, pCommonAttribList, ppInFaceTable != NULL, wScale, maskIdx);
	if (!jobCount) {
		//no jobs sent
		//implement an STUC_CANCELLED status
		return STUC_SUCCESS;
	}
	CLOCK_STOP("Send Off Time");
	CLOCK_START;
	pContext->threadPool.pWaitForJobs(pContext->pThreadPoolHandle, jobCount,
	                                  ppJobHandles, true, NULL);
	CLOCK_STOP("Waiting Time");
	bool empty = true;
	for (int32_t i = 0; i < jobCount; ++i) {
		//STUC_ASSERT("", jobArgs[i].bufSize > 0);
		//you'll need to handle this properly when you re-enable multithreading
		if (pJobArgs[i].bufSize > 0) {
			empty = false;
		}
	}
	err = stucJobGetErrs(pContext, jobCount, &ppJobHandles);
	err = stucJobDestroyHandles(pContext, jobCount, &ppJobHandles);
	pContext->alloc.pFree(ppJobHandles);
	STUC_THROW_IF(err, true, "", 0);
	if (empty) {
		return err;
	}
	CLOCK_START;
	Mesh meshOutWrap = {0};
	err = stucCombineJobMeshes(pContext, pMap, &meshOutWrap, pJobArgs, pEdgeVerts,
	                          pVertSeamTable, pEdgeSeamTable, ppInFaceTable, wScale,
	                          pMeshIn, jobCount);
	CLOCK_STOP("Combine time");
	CLOCK_START;
	reallocMeshToFit(&pContext->alloc, &meshOutWrap);
	*pMeshOut = meshOutWrap.core;
	CLOCK_STOP("Realloc time");
	STUC_CATCH(0, err, ;);
	pContext->alloc.pFree(pEdgeVerts);
	pContext->alloc.pFree(pInVertTable);
	pContext->alloc.pFree(pVertSeamTable);
	pContext->alloc.pFree(pEdgeSeamTable);
	pContext->alloc.pFree(pJobArgs);
	return err;
}

static
void InFaceTableToHashTable(StucAlloc *pAlloc,
                            StucMap pMap, int32_t count, InFaceArr *pInFaceTable) {
	UsgInFace **ppHashTable = &pMap->usgArr.pInFaceTable;
	pMap->usgArr.tableSize = count * 2;
	*ppHashTable = pAlloc->pCalloc(pMap->usgArr.tableSize, sizeof(UsgInFace));
	for (int32_t i = 0; i < count; ++i) {
		for (int32_t j = 0; j < pInFaceTable[i].count; ++j) {
			uint32_t sum = pInFaceTable[i].usg + pInFaceTable[i].pArr[j];
			int32_t hash = stucFnvHash((uint8_t *)&sum, 4, pMap->usgArr.tableSize);
			UsgInFace *pEntry = *ppHashTable + hash;
			if (!pEntry->pEntry) {
				pEntry->pEntry = pInFaceTable + i;
				pEntry->face = pInFaceTable[i].pArr[j];
				continue;
			}
			do {
				if (!pEntry->pNext) {
					pEntry = pEntry->pNext = pAlloc->pCalloc(1, sizeof(UsgInFace));
					pEntry->pEntry = pInFaceTable + i;
					pEntry->face = pInFaceTable[i].pArr[j];
					break;
				}
				pEntry = pEntry->pNext;
			} while(true);
		}
	}
}

static
Result getMatsToAdd(StucContext pContext, StucMapArr *pMapArr, int32_t mapIdx,
                    CommonAttribList *pCommonAttribs,AttribIndexedArr *pInIndexedAttribs,
                    AttribIndexed **ppMatsToAdd) {
	AttribIndexedArr *pMapAttribArr = &pMapArr->ppArr[mapIdx]->indexedAttribs;
	CommonAttribList *pMapCommon = pCommonAttribs + mapIdx;
	char *pAttribName = pContext->spAttribs[STUC_ATTRIB_SP_MAT_IDX];
	CommonAttrib *pCommonAttrib = getCommonAttrib(pMapCommon->pFace, pMapCommon->faceCount,
	                                              pAttribName);
	BlendConfig *pConfig = {0};
	if (pCommonAttrib) {
		pConfig = &pCommonAttrib->blendConfig;
	}
	else {
		pConfig = &getTypeDefaultConfig(&pContext->typeDefaults,
		                                STUC_ATTRIB_STRING)->blendConfig;
	}
	char *pName = "StucMaterials";
	AttribIndexed *pInMats = getAttribIndexed(pInIndexedAttribs, pName);
	AttribIndexed *pMapMats = getAttribIndexed(pMapAttribArr, pName);
	if (!pInMats) {
		*ppMatsToAdd = pMapMats;
	}
	else if (!pMapMats) {
		*ppMatsToAdd = pInMats;
	}
	else {
		*ppMatsToAdd = pConfig->order ? pInMats : pMapMats;
	}
	return STUC_SUCCESS;
}

static
void appendToOutMatsBuf(StucAlloc *pAlloc, Mesh *pMesh,
                        AttribIndexed *pOutMats, char *pMatName) {
	STUC_ASSERT("", pOutMats->count <= pOutMats->size);
	if (pOutMats->count == pOutMats->size) {
		pOutMats->size *= 2;
		reallocAttrib(pAlloc, pMesh, &pOutMats->core, pOutMats->size);
	}
	memcpy(attribAsStr(&pOutMats->core, pOutMats->count), pMatName,
			STUC_ATTRIB_STRING_MAX_LEN);
	pOutMats->count++;
}

typedef struct {
	int8_t idx;
	bool hasRef;
} MatTableEntry;

static
void appendToOutMats(AttribIndexed *pOutMats, char *pOutMatBuf,
                     MatTableEntry *pEntry) {
	pEntry->hasRef = true;
	char *pDest = attribAsStr(&pOutMats->core, pOutMats->count);
	char *pSrc = pOutMatBuf + pEntry->idx * STUC_ATTRIB_STRING_MAX_LEN;
	memcpy(pDest, pSrc, STUC_ATTRIB_STRING_MAX_LEN);
	pEntry->idx = pOutMats->count;
	pOutMats->count++;
}

static
Result correctMatIndices(StucContext pContext, Mesh *pMeshArr, StucMapArr *pMapArr,
                         CommonAttribList *pCommonAttribs,
                         AttribIndexedArr *pInIndexedAttribs,
                         AttribIndexedArr *pOutIndexedAttribs) {
	Result err = STUC_SUCCESS;
	StucAlloc *pAlloc = &pContext->alloc;
	pOutIndexedAttribs->size = 1;
	pOutIndexedAttribs->count = 1;
	pOutIndexedAttribs->pArr = pAlloc->pCalloc(pOutIndexedAttribs->size, sizeof(AttribIndexed));
	AttribIndexed *pOutMats = pOutIndexedAttribs->pArr;
	pOutMats->size = pMapArr->count;
	pOutMats->size += (pOutMats->size % 2);
	initAttribCore(pAlloc, &pOutMats->core, "StucMaterials", pOutMats->size,
	               STUC_ATTRIB_STRING);
	MatTableEntry **ppMatTable = pAlloc->pCalloc(pMapArr->count, sizeof(void *));

	for (int32_t i = 0; i < pMapArr->count; ++i) {
		Mesh *pMesh = pMeshArr + i;
		setSpecialAttribs(pContext, pMesh, 0x800);//only set mat indices
		if (!pMesh->pMatIdx) {
			continue;
		}
		AttribIndexed *pMatsToAdd = NULL;
		err = getMatsToAdd(pContext, pMapArr, i, pCommonAttribs, pInIndexedAttribs,
		                   &pMatsToAdd);
		STUC_ASSERT("mesh faces have mat indices, but map has no materials?",
		            pMatsToAdd && pMatsToAdd->count);
		STUC_THROW_IF(err, true, "", 0);
		ppMatTable[i] = pAlloc->pCalloc(pMatsToAdd->count, sizeof(MatTableEntry));

		for (int32_t j = 0; j < pMesh->core.faceCount; ++j) {
			int32_t matIdx = pMesh->pMatIdx[j];
			STUC_THROW_IF(err, matIdx >= 0 && matIdx < pMatsToAdd->count, "", 0);
			MatTableEntry *pEntry = ppMatTable[i] + matIdx;
			if (!pEntry->hasRef) {
				pEntry->hasRef = true;
				//We're just through material slots, so linear search should be fine for now
				char *pMatName = attribAsStr(&pMatsToAdd->core, matIdx);
				int32_t idx = -1;
				for (int32_t k = 0; k < pOutMats->count; ++k) {
					if (!strncmp(attribAsStr(&pOutMats->core, k), pMatName, STUC_ATTRIB_STRING_MAX_LEN)) {
						idx = k;
						break;
					}
				}
				if (idx >= 0) {
					pEntry->idx = idx;
				}
				else {
					pEntry->idx = pOutMats->count;
					appendToOutMatsBuf(pAlloc, pMesh, pOutMats, pMatName);
				}
			}
			pMesh->pMatIdx[j] = pEntry->idx;
		}
		STUC_CATCH(0, err, ;)
		pAlloc->pFree(ppMatTable[i]);
		ppMatTable[i] = NULL;
		STUC_THROW_IF(err, true, "", 1);
	}
	STUC_CATCH(1, err, ;)
	pAlloc->pFree(ppMatTable);
	return err;
}

typedef struct {
	StucContext pContext;
	StucMapArr *pMapArr;
	StucMesh *pMeshIn;
	StucAttribIndexedArr *pInIndexedAttribs;
	StucMesh *pMeshOut;
	StucAttribIndexedArr *pOutIndexedAttribs;
	StucCommonAttribList *pCommonAttribList;
	float wScale;
} StucMapToMeshArgs;

static
Result mapToMeshFromJob(void *pArgsVoid) {
	StucMapToMeshArgs *pArgs = pArgsVoid;
	return stucMapToMesh(pArgs->pContext, pArgs->pMapArr, pArgs->pMeshIn,
	                     pArgs->pInIndexedAttribs, pArgs->pMeshOut,
	                     pArgs->pOutIndexedAttribs, pArgs->pCommonAttribList,
	                     pArgs->wScale);
}

Result stucQueueMapToMesh(StucContext pContext, void **ppJobHandle, StucMapArr *pMapArr,
                          StucMesh *pMeshIn, StucAttribIndexedArr *pInIndexedAttribs,
                          StucMesh *pMeshOut, StucAttribIndexedArr *pOutIndexedAttribs,
                          StucCommonAttribList *pCommonAttribList, float wScale) {
	StucMapToMeshArgs *pArgs = pContext->alloc.pCalloc(1, sizeof(StucMapToMeshArgs));
	pArgs->pContext = pContext;
	pArgs->pMapArr = pMapArr;
	pArgs->pMeshIn = pMeshIn;
	pArgs->pInIndexedAttribs = pInIndexedAttribs;
	pArgs->pMeshOut = pMeshOut;
	pArgs->pOutIndexedAttribs = pOutIndexedAttribs;
	pArgs->pCommonAttribList = pCommonAttribList;
	pArgs->wScale = wScale;
	pContext->threadPool.pJobStackPushJobs(pContext->pThreadPoolHandle, 1, ppJobHandle,
	                                       mapToMeshFromJob, &pArgs);
	return STUC_SUCCESS;
}

Result stucMapToMesh(StucContext pContext, StucMapArr *pMapArr,
                     StucMesh *pMeshIn, StucAttribIndexedArr *pInIndexedAttribs,
                     StucMesh *pMeshOut, StucAttribIndexedArr *pOutIndexedAttribs,
                     StucCommonAttribList *pCommonAttribList, float wScale) {
	//TODO replace vars called 'result' with 'err'
	StucResult err = STUC_NOT_SET;
	if (!pMeshIn) {
		printf("Stuc map to mesh failed, pMeshIn was null\n");
		return STUC_ERROR;
	}
	if (!pMapArr) {
		printf("Stuc map to mesh failed, pMap was null\n");
		return STUC_ERROR;
	}
	Mesh meshInWrap = {.core = *pMeshIn};

	setSpecialAttribs(pContext, &meshInWrap, 0xf0e); //don't set preserve yet

	setAttribOrigins(&meshInWrap.core.meshAttribs, STUC_ATTRIB_ORIGIN_MESH_IN);
	setAttribOrigins(&meshInWrap.core.faceAttribs, STUC_ATTRIB_ORIGIN_MESH_IN);
	setAttribOrigins(&meshInWrap.core.cornerAttribs, STUC_ATTRIB_ORIGIN_MESH_IN);
	setAttribOrigins(&meshInWrap.core.edgeAttribs, STUC_ATTRIB_ORIGIN_MESH_IN);
	setAttribOrigins(&meshInWrap.core.vertAttribs, STUC_ATTRIB_ORIGIN_MESH_IN);

	setAttribToDontCopy(pContext, &meshInWrap, 0x3f0);

	if (isMeshInvalid(&meshInWrap)) {
		return STUC_ERROR;
	}

	buildTangents(&meshInWrap);

	if (!meshInWrap.core.edgeCount) {
		STUC_ASSERT("", !meshInWrap.core.edgeAttribs.count);
		STUC_ASSERT("", !meshInWrap.core.edgeAttribs.pArr);
		buildEdgeList(pContext, &meshInWrap);
	}

	printf("meshInWrap.pMatIdx[5] = %d\n", meshInWrap.pMatIdx[5]);

	Mesh *pOutBufArr = pContext->alloc.pCalloc(pMapArr->count, sizeof(Mesh));
	StucObject *pOutObjWrapArr =
		pContext->alloc.pCalloc(pMapArr->count, sizeof(StucObject));
	for (int32_t i = 0; i < pMapArr->count; ++i) {
		pOutObjWrapArr[i].pData = (StucObjectData *)&pOutBufArr[i];
		StucMap pMap = pMapArr->ppArr[i];
		int8_t matIdx = pMapArr->pMatArr[i];
		InFaceArr *pInFaceTable = NULL;
		if (pMap->usgArr.count) {
			MapFile squares = { .mesh = pMap->usgArr.squares };
			err = stucCreateQuadTree(pContext, &squares);
			STUC_THROW_IF(err, true, "failed to create quadtree", 0);
			StucMesh squaresOut = {0};
			err = mapToMeshInternal(pContext, &squares, &meshInWrap, &squaresOut, matIdx,
			                        pCommonAttribList + i, &pInFaceTable, 1.0f);
			STUC_THROW_IF(err, true, "map to mesh usg failed", 0);
			sampleInAttribsAtUsgOrigins(pContext, pMap, &meshInWrap, &squaresOut, pInFaceTable);
			InFaceTableToHashTable(&pContext->alloc, pMap, squaresOut.faceCount, pInFaceTable);
			//*pMeshOut = squaresOut;
			//return STUC_SUCCESS;
			stucMeshDestroy(pContext, &squaresOut);
		}
		setSpecialAttribs(pContext, &meshInWrap, 0x50); //set perserve if present

		err = mapToMeshInternal(pContext, pMap, &meshInWrap, &pOutBufArr[i].core, matIdx,
		                        pCommonAttribList + i, NULL, wScale);
		STUC_THROW_IF(err, true, "map to mesh failed", 0);
		STUC_CATCH(0, err, stucMeshDestroy(pContext, pOutBufArr + i);)
		if (pMap->usgArr.count) {
			pContext->alloc.pFree(pMap->usgArr.pInFaceTable);
			pMap->usgArr.pInFaceTable = NULL;
			for (int32_t j = 0; j < pMap->usgArr.count; ++j) {
				//TODO uncomment this and fix any memory issue
				//pContext->alloc.pFree(pInFaceTable[i].pArr);
			}
			pContext->alloc.pFree(pInFaceTable);
		}
		STUC_THROW_IF(err, true, "", 1);
	}
	pMeshOut->type.type = STUC_OBJECT_DATA_MESH;
	Mesh meshOutWrap = {.core = *pMeshOut};
	printf("merging obj arr\n");
	err = correctMatIndices(pContext, pOutBufArr, pMapArr, pCommonAttribList,
	                        pInIndexedAttribs, pOutIndexedAttribs);
	STUC_THROW_IF(err, true, "", 1);
	mergeObjArr(pContext, &meshOutWrap, pMapArr->count, pOutObjWrapArr, false);
	printf("post-merging obj arr\n");
	*pMeshOut = meshOutWrap.core;
	printf("----------------------FINISHING IN MESH WITH FACE COUNT %d\n", pMeshIn->faceCount);
	STUC_CATCH(1, err, stucMeshDestroy(pContext, pMeshOut);)
	return err;
}

Result stucObjArrDestroy(StucContext pContext,
                         int32_t objCount, StucObject *pObjArr) {
	return destroyObjArr(pContext, objCount, pObjArr);
}

Result stucUsgArrDestroy(StucContext pContext,
                         int32_t count, StucUsg *pUsgArr) {
	Result err = STUC_NOT_SET;
	for (int32_t i = 0; i < count; ++i) {
		err = stucMeshDestroy(pContext, (StucMesh *)pUsgArr[i].obj.pData);
		STUC_THROW_IF(err, true, "", 0);
	}
	pContext->alloc.pFree(pUsgArr);
	STUC_CATCH(0, err, ;)
	return err;
}

StucResult stucMeshDestroy(StucContext pContext, StucMesh *pMesh) {
	for (int32_t i = 0; i < pMesh->meshAttribs.count; ++i) {
		if (pMesh->meshAttribs.pArr[i].core.pData) {
			pContext->alloc.pFree(pMesh->meshAttribs.pArr[i].core.pData);
		}
	}
	if (pMesh->meshAttribs.count && pMesh->meshAttribs.pArr) {
		pContext->alloc.pFree(pMesh->meshAttribs.pArr);
	}
	if(pMesh->pFaces) {
		pContext->alloc.pFree(pMesh->pFaces);
	}
	for (int32_t i = 0; i < pMesh->faceAttribs.count; ++i) {
		if (pMesh->faceAttribs.pArr[i].core.pData) {
			pContext->alloc.pFree(pMesh->faceAttribs.pArr[i].core.pData);
		}
	}
	if (pMesh->faceAttribs.count && pMesh->faceAttribs.pArr) {
		pContext->alloc.pFree(pMesh->faceAttribs.pArr);
	}
	if (pMesh->pCorners) {
		pContext->alloc.pFree(pMesh->pCorners);
	}
	for (int32_t i = 0; i < pMesh->cornerAttribs.count; ++i) {
		if (pMesh->cornerAttribs.pArr[i].core.pData) {
			pContext->alloc.pFree(pMesh->cornerAttribs.pArr[i].core.pData);
		}
	}
	if (pMesh->cornerAttribs.count && pMesh->cornerAttribs.pArr) {
		pContext->alloc.pFree(pMesh->cornerAttribs.pArr);
	}
	if (pMesh->pEdges) {
		pContext->alloc.pFree(pMesh->pEdges);
	}
	for (int32_t i = 0; i < pMesh->edgeAttribs.count; ++i) {
		if (pMesh->edgeAttribs.pArr[i].core.pData) {
			pContext->alloc.pFree(pMesh->edgeAttribs.pArr[i].core.pData);
		}
	}
	if (pMesh->edgeAttribs.count && pMesh->edgeAttribs.pArr) {
		pContext->alloc.pFree(pMesh->edgeAttribs.pArr);
	}
	for (int32_t i = 0; i < pMesh->vertAttribs.count; ++i) {
		if (pMesh->vertAttribs.pArr[i].core.pData) {
			pContext->alloc.pFree(pMesh->vertAttribs.pArr[i].core.pData);
		}
	}
	if (pMesh->vertAttribs.count && pMesh->vertAttribs.pArr) {
		pContext->alloc.pFree(pMesh->vertAttribs.pArr);
	}
	return STUC_SUCCESS;
}

StucResult stucGetAttribSize(StucAttrib *pAttrib, int32_t *pSize) {
	*pSize = getAttribSize(pAttrib->core.type);
	return STUC_SUCCESS;
}

StucResult stucGetAttrib(char *pName, StucAttribArray *pAttribs, StucAttrib **ppAttrib) {
	*ppAttrib = getAttrib(pName, pAttribs);
	return STUC_SUCCESS;
}

StucResult stucGetAttribIndexed(char *pName, StucAttribIndexedArr *pAttribs,
                                StucAttribIndexed **ppAttrib) {
	*ppAttrib = getAttribIndexed(pAttribs, pName);
	return STUC_SUCCESS;
}

typedef struct {
	StucImage imageBuf;
	StucMap pMap;
	StucContext pContext;
	int32_t *pActiveJobs;
	void *pMutex;
	int32_t bufOffset;
	int32_t pixelCount;
	int8_t id;
	V2_F32 zBounds;
} RenderArgs;

static void testPixelAgainstFace(RenderArgs *pVars, V2_F32 *pPos, FaceRange *pFace, Color *pColor) {
	Mesh* pMesh = &pVars->pMap->mesh;
	V2_F32 verts[4];
	for (int32_t i = 0; i < pFace->size; ++i) {
		verts[i] = *(V2_F32 *)(pMesh->pVerts + pMesh->core.pCorners[pFace->start + i]);
	}
	int8_t triCorners[4] = {0};
	V3_F32 bc = getBarycentricInFace(verts, triCorners, pFace->size, *pPos);
	if (bc.d[0] < -.0001f || bc.d[1] < -.0001f || bc.d[2] < -.0001f ||
		!v3IsFinite(bc)) {
		return;
	}
	V3_F32 vertsXyz[3];
	for (int32_t i = 0; i < 3; ++i) {
		int32_t vertIdx = pMesh->core.pCorners[pFace->start + triCorners[i]];
		vertsXyz[i] = pMesh->pVerts[vertIdx];
	}
	V3_F32 wsPos = barycentricToCartesian(vertsXyz, &bc);
	//the alpha channel is used as a depth buffer
	if (wsPos.d[2] <= pColor->d[3]) {
		return;
	}
	V3_F32 ab = _(vertsXyz[1] V3SUB vertsXyz[0]);
	V3_F32 ac = _(vertsXyz[2] V3SUB vertsXyz[0]);
	V3_F32 normal = v3Cross(ab, ac);
	float normalLen =
		sqrt(normal.d[0] * normal.d[0] +
		     normal.d[1] * normal.d[1] +
		     normal.d[2] * normal.d[2]);
	_(&normal V3DIVEQLS normalLen);
	V3_F32 up = { .0f, .0f, 1.0f };
	float dotUp = _(normal V3DOT up);
	if (dotUp < .0f) {
		return;
	}
	float depth = (wsPos.d[2] - pVars->zBounds.d[0]) / pVars->zBounds.d[1];
	float value = dotUp;
	value *= 1.0f - (1.0f - depth) * .5f;
	value *= .75;
	pColor->d[0] = value;
	pColor->d[1] = value;
	pColor->d[2] = value;
	pColor->d[3] = wsPos.d[2];
}

static void stucRenderJob(void *pArgs) {
	RenderArgs vars = *(RenderArgs *)pArgs;
	int32_t dataLen = vars.pixelCount * getPixelSize(vars.imageBuf.type);
	vars.imageBuf.pData = vars.pContext->alloc.pMalloc(dataLen);
	Mesh *pMesh = &vars.pMap->mesh;
	float pixelScale = 1.0 / (float)vars.imageBuf.res;
	float pixelHalfScale = pixelScale / 2.0f;
	FaceCells faceCells = {0};
	FaceCellsTable faceCellsTable = {.pFaceCells = &faceCells};
	QuadTreeSearch searchState = {0};
	stucInitQuadTreeSearch(&vars.pContext->alloc, vars.pMap, &searchState);
	for (int32_t i = 0; i < vars.pixelCount; ++i) {
		int32_t iOffset = vars.bufOffset + i;
		V2_F32 idx = {iOffset % vars.imageBuf.res,
		                iOffset / vars.imageBuf.res};
		V2_F32 pos = {pixelScale * idx.d[0] + pixelHalfScale,
		              pixelScale * idx.d[1] + pixelHalfScale};
		Color color = { 0 };
		color.d[3] = FLT_MAX * -1.0f;
		Range faceRange = {.start = 0, .end = 1};
		stucGetCellsForSingleFace(&searchState, 1, &pos, &faceCellsTable,
		                          NULL, 0, faceRange);
		int32_t leafIdx = faceCells.pCells[faceCells.cellSize - 1];
		Cell *pLeaf = vars.pMap->quadTree.cellTable.pArr + leafIdx;
		for (int32_t j = 0; j < faceCells.cellSize; ++j) {
			int32_t cellIdx = faceCells.pCells[j];
			Cell *pCell = vars.pMap->quadTree.cellTable.pArr + cellIdx;
			int32_t* pCellFaces;
			Range cellFaceRange = {0};
			if (faceCells.pCellType[j]) {
				pCellFaces = pCell->pEdgeFaces;
				cellFaceRange = pLeaf->pLinkEdgeRanges[j];
			}
			else if (faceCells.pCellType[j] != 1) {
				pCellFaces = pCell->pFaces;
				cellFaceRange.start = 0;
				cellFaceRange.end = pCell->faceSize;
			}
			else {
				continue;
			}
			for (int32_t k = cellFaceRange.start; k < cellFaceRange.end; ++k) {
				FaceRange face = {0};
				face.idx = pCellFaces[k];
				face.start = pMesh->core.pFaces[face.idx];
				face.end = pMesh->core.pFaces[face.idx + 1];
				face.size = face.end - face.start;
				FaceTriangulated faceTris = {0};
				if (face.size > 4) {
					faceTris = triangulateFace(vars.pContext->alloc, face,
					                           pMesh->pVerts,
					                           pMesh->core.pCorners, 0);
					for (int32_t l = 0; l < faceTris.triCount; ++l) {
						FaceRange tri = {0};
						tri.idx = face.idx;
						tri.start = l * 3;
						tri.end = tri.start + 3;
						tri.size = tri.end - tri.start;
						tri.start += face.start;
						tri.end += face.start;
						testPixelAgainstFace(&vars, &pos, &tri, &color);
					}
				}
				else {
					testPixelAgainstFace(&vars, &pos, &face, &color);
				}
			}
		}
		color.d[3] = color.d[3] != FLT_MAX * -1.0f;
		setPixelColor(&vars.imageBuf, i, &color);
	}
	stucDestroyQuadTreeSearch(&searchState);
	*(RenderArgs *)pArgs = vars;
	StucThreadPool *pThreadPool = &vars.pContext->threadPool;
	pThreadPool->pMutexLock(vars.pContext->pThreadPoolHandle, vars.pMutex);
	--*vars.pActiveJobs;
	pThreadPool->pMutexUnlock(vars.pContext->pThreadPoolHandle, vars.pMutex);
}

static
V2_F32 getZBounds(StucMap pMap) {
	Mesh* pMesh = &pMap->mesh;
	V2_F32 zBounds = {.d = {FLT_MAX, FLT_MIN}};
	for (int32_t i = 0; i < pMesh->core.vertCount; ++i) {
		if (pMesh->pVerts[i].d[2] < zBounds.d[0]) {
			zBounds.d[0] = pMesh->pVerts[i].d[2];
		}
		if (pMesh->pVerts[i].d[2] > zBounds.d[1]) {
			zBounds.d[1] = pMesh->pVerts[i].d[2];
		}
	}
	return zBounds;
}

StucResult stucMapFileGenPreviewImage(StucContext pContext, StucMap pMap, StucImage *pImage) {
	Result err = STUC_SUCCESS;
	V2_F32 zBounds = getZBounds(pMap);
	int32_t pixelCount = pImage->res * pImage->res;
	int32_t pixelsPerJob = pixelCount / pContext->threadCount;
	bool singleThread = !pixelsPerJob;
	void *pMutex = NULL;
	pContext->threadPool.pMutexGet(pContext->pThreadPoolHandle, &pMutex);
	int32_t activeJobs = 0;
	void *jobArgPtrs[MAX_THREADS];
	RenderArgs args[MAX_THREADS];
	activeJobs = singleThread ? 1 : pContext->threadCount;
	for (int32_t i = 0; i < activeJobs; ++i) {
		args[i].bufOffset = i * pixelsPerJob;
		args[i].imageBuf.res = pImage->res;
		args[i].imageBuf.type = pImage->type;
		args[i].pContext = pContext;
		args[i].pixelCount = i == activeJobs - 1 ?
			pixelCount - args[i].bufOffset : pixelsPerJob;
		args[i].pMap = pMap;
		args[i].zBounds = zBounds;
		args[i].pMutex = pMutex;
		args[i].pActiveJobs = &activeJobs;
		args[i].id = i;
		jobArgPtrs[i] = args + i;
	}
	void **ppJobHandles = pContext->alloc.pCalloc(activeJobs, sizeof(void *));
	pContext->threadPool.pJobStackPushJobs(pContext->pThreadPoolHandle, &ppJobHandles,
	                                       activeJobs, stucRenderJob, jobArgPtrs);
	stucWaitForJobsIntern(pContext->pThreadPoolHandle, activeJobs, ppJobHandles, true, NULL);
	err = stucJobGetErrs(pContext, activeJobs, &ppJobHandles);
	err = stucJobDestroyHandles(pContext, activeJobs, &ppJobHandles);
	pContext->alloc.pFree(ppJobHandles);
	STUC_THROW_IF(err, true, "", 0);
	int32_t pixelSize = getPixelSize(pImage->type);
	pImage->pData = pContext->alloc.pMalloc(pixelCount * pixelSize);
	for (int32_t i = 0; i < activeJobs; ++i) {
		void *pImageOffset = offsetImagePtr(pImage, i * pixelsPerJob);
		int32_t bytesToCopy = pixelSize * args[i].pixelCount;
		memcpy(pImageOffset, args[i].imageBuf.pData, bytesToCopy);
	}
	STUC_CATCH(0, err, ;);
	pContext->threadPool.pMutexDestroy(pContext->pThreadPoolHandle, pMutex);
	for (int32_t i = 0; i < activeJobs; ++i) {
		pContext->alloc.pFree(args[i].imageBuf.pData);
	}
	return err;
}

void stucMapIndexedAttribsGet(StucContext pContext, StucMap pMap,
                              StucAttribIndexedArr *pIndexedAttribs) {
	*pIndexedAttribs = pMap->indexedAttribs;
}

Result stucWaitForJobs(StucContext pContext, int32_t count, void **ppHandles,
                       bool wait, bool *pDone) {
	return pContext->threadPool.pWaitForJobs(pContext->pThreadPoolHandle, count,
	                                         ppHandles, wait, pDone);
}

Result stucJobGetErrs(StucContext pContext, int32_t jobCount,
                                  void ***pppJobHandles) {
	Result err = STUC_SUCCESS;
	STUC_THROW_IF(err, pContext && pppJobHandles, "", 0);
	STUC_THROW_IF(err, jobCount > 0, "", 0);
	for (int32_t i = 0; i < jobCount; ++i) {
		StucResult jobErr = STUC_NOT_SET;
		err = pContext->threadPool.pGetJobErr(pContext->pThreadPoolHandle,
		                                      (*pppJobHandles)[i], &jobErr);
		STUC_THROW_IF(err, jobErr == STUC_SUCCESS, "", 0);
	}
	STUC_CATCH(0, err, ;);
	return err;
}

Result stucJobDestroyHandles(StucContext pContext, int32_t jobCount,
                             void ***pppJobHandles) {
	Result err = STUC_SUCCESS;
	STUC_THROW_IF(err, pContext && pppJobHandles, "", 0);
	STUC_THROW_IF(err, jobCount > 0, "", 0);
	for (int32_t i = 0; i < jobCount; ++i) {
		pContext->threadPool.pJobHandleDestroy(pContext->pThreadPoolHandle,
		                                       *pppJobHandles + i);
	}
	STUC_CATCH(0, err, ;);
	return err;
}