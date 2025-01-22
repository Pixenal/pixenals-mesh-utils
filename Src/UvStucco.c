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
// - Split whole quadtree into chunks?
//
// - Add blending options to interface, that control how MeshIn attributes blend with
//   those from the Map. Also add an option to disable or enable interpolation.
//   Add these to the StucAttrib struct.
//
//TODO repalce localMesh with bufMesh.
//The old name is still present in some functions & vars
//TODO a highly distorted meshIn can cause invalid geometry
//(enough to crash blender). When meshIn is quads atleast
//(I've not tested with tris). Find out why
//TODO stucPreserve isn't working.
//TODO add option to vary z projection depth with uv stretch (for wires and such)
//Add option to mask stucpreserve by map edges. Where masking is defined
//per edge, not per face. An preserve meshin edge must pass through 2 map
//edges which are marked preserve, for the edge to cut that map face.
//This will cause there to be gaps in the corner, in cases of diagional meshin
//preserve edges. Just triangulate these faces, or leave as an ngon?
//TODO Add an option for subdivision like smoothing (for instances where
//the map is higher res than the base mesh). So that the surface can be
//smoothed, without needing to induce the perf cost of actually subdividing
//the base mesh. Is this possible?
//TODO add a cache (hash table?) for triangulated in faces, as your going to need to reference
//their edges when merging border faces
//TODO add user define void * args to custom callbacks
//TODO add the ability to open map files in dcc, to make edits to mesh, USGs etc

static
void stucSetTypeDefaultConfig(StucContext pContext) {
	StucTypeDefaultConfig config = {0};
	pContext->typeDefaults = config;
}

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
	return STUC_SUCCESS;
}

StucResult stucContextDestroy(StucContext pContext) {
	pContext->threadPool.pDestroy(pContext->pThreadPoolHandle);
	pContext->alloc.pFree(pContext);
	return STUC_SUCCESS;
}

StucResult stucMapFileExport(StucContext pContext, const char *pName,
                             int32_t objCount, StucObject* pObjArr,
                             int32_t usgCount, StucUsg* pUsgArr,
                             StucAttribIndexedArr indexedAttribs) {
	return stucWriteStucFile(pContext, pName, objCount, pObjArr,
	                         usgCount, pUsgArr, indexedAttribs);
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
	StucResult status = STUC_NOT_SET;
	StucMap pMap = pContext->alloc.pCalloc(1, sizeof(MapFile));
	int32_t objCount = 0;
	StucObject *pObjArr = NULL;
	StucUsg *pUsgArr = NULL;
	int32_t flatCutoffCount = 0;
	StucObject *pFlatCutoffArr = NULL;
	status = stucLoadStucFile(pContext, filePath, &objCount, &pObjArr,
	                          &pMap->usgArr.count, &pUsgArr, &flatCutoffCount,
	                          &pFlatCutoffArr, false, &pMap->indexedAttribs);
	if (status != STUC_SUCCESS) {
		return status;
	}

	for (int32_t i = 0; i < objCount; ++i) {
		setSpecialAttribs(pObjArr[i].pData, 0xae); //10101110 - all except for preserve
		applyObjTransform(pObjArr + i);
	}
	pMap->mesh.mesh.type.type = STUC_OBJECT_DATA_MESH_INTERN;
	mergeObjArr(pContext, &pMap->mesh, objCount, pObjArr, false);
	setSpecialAttribs(&pMap->mesh, 0xae);
	//TODO some form of heap corruption when many objects
	//test with address sanitizer on CircuitPieces.stuc
	destroyObjArr(pContext, objCount, pObjArr);

	if (pMap->mesh.pUvAttrib) {
		//TODO as with all special attributes, allow user to define what should be considered
		//     the primary UV channel. This especially important for integration with other DCCs
		if (!strncmp(pMap->mesh.pUvAttrib->name, "UVMap", STUC_ATTRIB_NAME_MAX_LEN)) {
			char newName[STUC_ATTRIB_NAME_MAX_LEN] = "Map_UVMap";
			memcpy(pMap->mesh.pUvAttrib->name, newName, STUC_ATTRIB_NAME_MAX_LEN);
		}
	}

	//set corner attribs to interpolate by default
	//TODO make this an option in ui, even for non common attribs
	for (int32_t i = 0; i < pMap->mesh.mesh.cornerAttribs.count; ++i) {
		pMap->mesh.mesh.cornerAttribs.pArr[i].interpolate = true;
	}

	//the quadtree is created before USGs are assigned to verts,
	//as the tree's used to speed up the process
	printf("File loaded. Creating quad tree\n");
	stucCreateQuadTree(pContext, pMap);

	if (pMap->usgArr.count) {
		pMap->usgArr.pArr = pContext->alloc.pCalloc(pMap->usgArr.count, sizeof(Usg));
		for (int32_t i = 0; i < pMap->usgArr.count; ++i) {
			setSpecialAttribs(pUsgArr[i].obj.pData, 0x02); //000010 - set only vert pos
			pMap->usgArr.pArr[i].origin = *(V2_F32 *)&pUsgArr[i].obj.transform.d[3];
			pMap->usgArr.pArr[i].pMesh = pUsgArr[i].obj.pData;
			applyObjTransform(&pUsgArr[i].obj);
			if (pUsgArr[i].pFlatCutoff) {
				pMap->usgArr.pArr[i].pFlatCutoff = pUsgArr[i].pFlatCutoff->pData;
				setSpecialAttribs(pUsgArr[i].pFlatCutoff->pData, 0x02); //000010 - set only vert pos
				applyObjTransform(pUsgArr[i].pFlatCutoff);
			}
		}
		allocUsgSquaresMesh(&pContext->alloc, pMap);
		fillUsgSquaresMesh(pMap, pUsgArr);
		assignUsgsToVerts(&pContext->alloc, pMap, pUsgArr);
		pMap->usgArr.pMemArr = pUsgArr;
	}

	*pMapHandle = pMap;
	//TODO add proper checks, and return STUC_ERROR if fails.
	//Do for all public functions (or internal ones as well)
	return STUC_SUCCESS;
}

StucResult stucMapFileUnload(StucContext pContext, StucMap pMap) {
	stucDestroyQuadTree(pContext, &pMap->quadTree);
	stucMeshDestroy(pContext, &pMap->mesh.mesh);
	pContext->alloc.pFree(pMap);
	return STUC_SUCCESS;
}

static
void initCommonAttrib(StucContext pContext, StucCommonAttrib *pEntry,
                      StucAttrib *pAttrib) {
	memcpy(pEntry->name, pAttrib->name, STUC_ATTRIB_NAME_MAX_LEN);
	StucTypeDefault *pDefault = 
		getTypeDefaultConfig(&pContext->typeDefaults, pAttrib->type);
	pEntry->blendConfig = pDefault->blendConfig;
}

static
void getCommonAttribs(StucContext pContext, AttribArray *pMapAttribs,
					  AttribArray *pMeshAttribs,
					  int32_t *pCommonAttribCount,
					  StucCommonAttrib **ppCommonAttribs) {
	if (!pMeshAttribs || !pMapAttribs) {
		return;
	}
	int32_t count = 0;
	for (int32_t i = 0; i < pMeshAttribs->count; ++i) {
		for (int32_t j = 0; j < pMapAttribs->count; ++j) {
			if (!strncmp(pMeshAttribs->pArr[i].name,
			             pMapAttribs->pArr[j].name,
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
			if (!strncmp(pMeshAttribs->pArr[i].name,
			             pMapAttribs->pArr[j].name,
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
	getCommonAttribs(pContext, &pMap->mesh.mesh.meshAttribs, &pMesh->meshAttribs,
					 &pCommonAttribs->meshCount, &pCommonAttribs->pMesh);
	getCommonAttribs(pContext, &pMap->mesh.mesh.faceAttribs, &pMesh->faceAttribs,
					 &pCommonAttribs->faceCount, &pCommonAttribs->pFace);
	getCommonAttribs(pContext, &pMap->mesh.mesh.cornerAttribs, &pMesh->cornerAttribs,
					 &pCommonAttribs->cornerCount, &pCommonAttribs->pCorner);
	getCommonAttribs(pContext, &pMap->mesh.mesh.edgeAttribs, &pMesh->edgeAttribs,
	                 &pCommonAttribs->edgeCount, &pCommonAttribs->pEdge);
	getCommonAttribs(pContext, &pMap->mesh.mesh.vertAttribs, &pMesh->vertAttribs,
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
void sendOffJobs(StucContext pContext, StucMap pMap, SendOffArgs *pJobArgs,
                 int32_t *pActiveJobs, int32_t *pMapJobsSent, Mesh *pMesh, void *pMutex,
                 EdgeVerts *pEdgeVerts, int8_t *pInVertTable,
				 StucCommonAttribList *pCommonAttribList,
	             bool getInFaces, float wScale) {
	//struct timeval start, stop;
	//CLOCK_START;
	int32_t facesPerThread = pMesh->mesh.faceCount / pContext->threadCount;
	bool singleThread = !facesPerThread;
	void *jobArgPtrs[MAX_THREADS];
	int32_t borderTableSize = pMap->mesh.mesh.faceCount / 5 + 2; //+ 2 incase is 0
	printf("fromjobsendoff: BorderTableSize: %d\n", borderTableSize);
	*pActiveJobs = singleThread ? 1 : pContext->threadCount;
	for (int32_t i = 0; i < *pActiveJobs; ++i) {
		int32_t meshStart = facesPerThread * i;
		int32_t meshEnd = i == *pActiveJobs - 1 ?
			pMesh->mesh.faceCount : meshStart + facesPerThread;
		Mesh meshPart = *pMesh;
		meshPart.mesh.pFaces += meshStart;
		meshPart.mesh.faceCount = meshEnd - meshStart;
		pJobArgs[i].inFaceOffset = meshStart;
		pJobArgs[i].pInVertTable = pInVertTable;
		pJobArgs[i].pEdgeVerts = pEdgeVerts;
		pJobArgs[i].pMap = pMap;
		pJobArgs[i].borderTable.size = borderTableSize;
		pJobArgs[i].mesh = meshPart;
		pJobArgs[i].pActiveJobs = pActiveJobs;
		pJobArgs[i].id = i;
		pJobArgs[i].pContext = pContext;
		pJobArgs[i].pMutex = pMutex;
		pJobArgs[i].pCommonAttribList = pCommonAttribList;
		pJobArgs[i].getInFaces = getInFaces;
		pJobArgs[i].wScale = wScale;
		jobArgPtrs[i] = pJobArgs + i;
	}
	*pMapJobsSent = *pActiveJobs;
	pContext->threadPool.pJobStackPushJobs(pContext->pThreadPoolHandle,
	                                       *pActiveJobs, stucMapToJobMesh, jobArgPtrs);
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
	*ppInVertTable = pContext->alloc.pCalloc(pMesh->mesh.vertCount, 1);
	*ppVertSeamTable = pContext->alloc.pCalloc(pMesh->mesh.vertCount, 1);
	*ppEdgeSeamTable = pContext->alloc.pCalloc(pMesh->mesh.edgeCount, 1);
	//TODO do we need to list number of unique preserve edges per vert?
	//I'm not doing so currently (hence why pEdgeCache is commented out),
	//and it seems to be working. (talking about split to pieces)
	EdgeCache *pEdgeCache =
		pContext->alloc.pCalloc(pMesh->mesh.vertCount, sizeof(EdgeCache));
	for (int32_t i = 0; i < pMesh->mesh.faceCount; ++i) {
		FaceRange face = {0};
		face.start = pMesh->mesh.pFaces[i];
		face.end = pMesh->mesh.pFaces[i + 1];
		face.size = face.end - face.start;
		face.idx = i;
		for (int32_t j = 0; j < face.size; ++j) {
			int32_t corner = face.start + j;
			int32_t vert = pMesh->mesh.pCorners[corner];
			int32_t edge = pMesh->mesh.pEdges[corner];
			addVertToTableEntry(pMesh, face, j, vert, edge, pEdgeVerts,
			                    *ppVertSeamTable, *ppInVertTable, pEdgeCache,
			                    *ppEdgeSeamTable);
			int32_t prevj = j == 0 ? face.size - 1 : j - 1;
			int32_t prevEdge = pMesh->mesh.pEdges[face.start + prevj];
			addVertToTableEntry(pMesh, face, prevj, vert, prevEdge, pEdgeVerts,
			                    *ppVertSeamTable, *ppInVertTable, pEdgeCache,
			                    *ppEdgeSeamTable);
		}
	}
	pContext->alloc.pFree(pEdgeCache);
}

static
void setAttribOrigins(AttribArray *pAttribs, StucAttribOrigin origin) {
	for (int32_t i = 0; i < pAttribs->count; ++i) {
		pAttribs->pArr[i].origin = origin;
	}
}

static
Result mapToMeshInternal(StucContext pContext, StucMap pMap, Mesh *pMeshIn,
                         StucMesh *pMeshOut, StucCommonAttribList *pCommonAttribList,
                         InFaceArr **ppInFaceTable, float wScale) {
	CLOCK_INIT;
	CLOCK_START;
	EdgeVerts *pEdgeVerts = {0};
	printf("EdgeCount: %d\n", pMeshIn->mesh.edgeCount);
	buildEdgeVertsTable(pContext, &pEdgeVerts, &pMeshIn->mesh);
	int8_t *pInVertTable;
	int8_t *pVertSeamTable;
	bool *pEdgeSeamTable;
	buildVertTables(pContext, pMeshIn, &pInVertTable,
	                &pVertSeamTable, pEdgeVerts, &pEdgeSeamTable);
	CLOCK_STOP("Edge Table Time");

	CLOCK_START;
	SendOffArgs jobArgs[MAX_THREADS] = {0};
	int32_t activeJobs = 0;
	int32_t mapJobsSent = 0;
	void *pMutex = NULL;
	pContext->threadPool.pMutexGet(pContext->pThreadPoolHandle, &pMutex);
	sendOffJobs(pContext, pMap, jobArgs, &activeJobs, &mapJobsSent, pMeshIn, pMutex,
	            pEdgeVerts, pInVertTable, pCommonAttribList,
	            ppInFaceTable != NULL, wScale);
	if (!mapJobsSent) {
		//no jobs sent
		//implement an STUC_CANCELLED status
		return STUC_SUCCESS;
	}
	CLOCK_STOP("Send Off Time");
	CLOCK_START;
	waitForJobs(pContext, &activeJobs, pMutex);
	pContext->threadPool.pMutexDestroy(pContext->pThreadPoolHandle, pMutex);
	CLOCK_STOP("Waiting Time");
	Result jobResult = STUC_SUCCESS;
	bool empty = true;
	for (int32_t i = 0; i < mapJobsSent; ++i) {
		//STUC_ASSERT("", jobArgs[i].bufSize > 0);
		//you'll need to handle this properly when you re-enable multithreading
		if (jobArgs[i].bufSize > 0) {
			empty = false;
		}
		if (jobArgs[i].result != STUC_SUCCESS) {
			jobResult = STUC_ERROR;
		}
	}
	if (empty || jobResult != STUC_SUCCESS) {
		return jobResult;
	}

	CLOCK_START;
	Mesh meshOutWrap = {0};
	stucCombineJobMeshes(pContext, pMap, &meshOutWrap, jobArgs, pEdgeVerts,
	                     pVertSeamTable, pEdgeSeamTable, ppInFaceTable, wScale,
	                     pMeshIn, mapJobsSent);
	CLOCK_STOP("Combine time");
	pContext->alloc.pFree(pEdgeVerts);
	pContext->alloc.pFree(pInVertTable);
	pContext->alloc.pFree(pVertSeamTable);
	pContext->alloc.pFree(pEdgeSeamTable);
	CLOCK_START;
	reallocMeshToFit(&pContext->alloc, &meshOutWrap);
	*pMeshOut = meshOutWrap.mesh;
	CLOCK_STOP("Realloc time");
	return STUC_SUCCESS;
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

Result stucMapToMesh(StucContext pContext, StucMap pMap, StucMesh *pMeshIn,
                     StucMesh *pMeshOut, StucCommonAttribList *pCommonAttribList,
                     float wScale) {
	if (!pMeshIn) {
		printf("Stuc map to mesh failed, pMeshIn was null\n");
		return STUC_ERROR;
	}
	if (!pMap) {
		printf("Stuc map to mesh failed, pMap was null\n");
		return STUC_ERROR;
	}
	Mesh meshInWrap = {.mesh = *pMeshIn};
	setSpecialAttribs(&meshInWrap, 0x70e); //don't set preserve yet
	if (isMeshInvalid(&meshInWrap)) {
		return STUC_ERROR;
	}
	buildTangents(&meshInWrap);
	//TODO remove this, I dont think it's necessary. Origin is only used in bufmesh
	//it doesn't matter what it's set to here
	setAttribOrigins(&meshInWrap.mesh.meshAttribs, STUC_ATTRIB_ORIGIN_MESH_IN);
	setAttribOrigins(&meshInWrap.mesh.faceAttribs, STUC_ATTRIB_ORIGIN_MESH_IN);
	setAttribOrigins(&meshInWrap.mesh.cornerAttribs, STUC_ATTRIB_ORIGIN_MESH_IN);
	setAttribOrigins(&meshInWrap.mesh.edgeAttribs, STUC_ATTRIB_ORIGIN_MESH_IN);
	setAttribOrigins(&meshInWrap.mesh.vertAttribs, STUC_ATTRIB_ORIGIN_MESH_IN);

	if (!meshInWrap.mesh.edgeCount) {
		STUC_ASSERT("", !meshInWrap.mesh.edgeAttribs.count);
		STUC_ASSERT("", !meshInWrap.mesh.edgeAttribs.pArr);
		buildEdgeList(pContext, &meshInWrap);
	}

	InFaceArr *pInFaceTable = NULL;
	if (pMap->usgArr.count) {
		MapFile squares = { .mesh = pMap->usgArr.squares };
		stucCreateQuadTree(pContext, &squares);
		StucMesh squaresOut = {0};
		mapToMeshInternal(pContext, &squares, &meshInWrap, &squaresOut, pCommonAttribList, &pInFaceTable, 1.0f);
		sampleInAttribsAtUsgOrigins(pMap, &meshInWrap, &squaresOut, pInFaceTable);
		InFaceTableToHashTable(&pContext->alloc, pMap, squaresOut.faceCount, pInFaceTable);
		//*pMeshOut = squaresOut;
		//return STUC_SUCCESS;
		stucMeshDestroy(pContext, &squaresOut);
	}
	setSpecialAttribs(&meshInWrap, 0x50); //set perserve if present

	mapToMeshInternal(pContext, pMap, &meshInWrap, pMeshOut, pCommonAttribList, NULL, wScale);
	if (pMap->usgArr.count) {
		pContext->alloc.pFree(pMap->usgArr.pInFaceTable);
		pMap->usgArr.pInFaceTable = NULL;
		for (int32_t i = 0; i < pMap->usgArr.count; ++i) {
			//pContext->alloc.pFree(pInFaceTable[i].pArr);
		}
		pContext->alloc.pFree(pInFaceTable);
	}
}

StucResult stucObjArrDestroy(StucContext pContext,
                             int32_t objCount, StucObject *pObjArr) {
	destroyObjArr(pContext, objCount, pObjArr);
}

StucResult stucUsgArrDestroy(StucContext pContext,
                                    int32_t count, StucUsg *pUsgArr) {
	for (int32_t i = 0; i < count; ++i) {
		stucMeshDestroy(pContext, pUsgArr[i].obj.pData);
	}
	pContext->alloc.pFree(pUsgArr);
	return STUC_SUCCESS;
}

StucResult stucMeshDestroy(StucContext pContext, StucMesh *pMesh) {
	for (int32_t i = 0; i < pMesh->meshAttribs.count; ++i) {
		if (pMesh->meshAttribs.pArr[i].pData) {
			pContext->alloc.pFree(pMesh->meshAttribs.pArr[i].pData);
		}
	}
	if (pMesh->meshAttribs.count && pMesh->meshAttribs.pArr) {
		pContext->alloc.pFree(pMesh->meshAttribs.pArr);
	}
	if(pMesh->pFaces) {
		pContext->alloc.pFree(pMesh->pFaces);
	}
	for (int32_t i = 0; i < pMesh->faceAttribs.count; ++i) {
		if (pMesh->faceAttribs.pArr[i].pData) {
			pContext->alloc.pFree(pMesh->faceAttribs.pArr[i].pData);
		}
	}
	if (pMesh->faceAttribs.count && pMesh->faceAttribs.pArr) {
		pContext->alloc.pFree(pMesh->faceAttribs.pArr);
	}
	if (pMesh->pCorners) {
		pContext->alloc.pFree(pMesh->pCorners);
	}
	for (int32_t i = 0; i < pMesh->cornerAttribs.count; ++i) {
		if (pMesh->cornerAttribs.pArr[i].pData) {
			pContext->alloc.pFree(pMesh->cornerAttribs.pArr[i].pData);
		}
	}
	if (pMesh->cornerAttribs.count && pMesh->cornerAttribs.pArr) {
		pContext->alloc.pFree(pMesh->cornerAttribs.pArr);
	}
	if (pMesh->pEdges) {
		pContext->alloc.pFree(pMesh->pEdges);
	}
	for (int32_t i = 0; i < pMesh->edgeAttribs.count; ++i) {
		if (pMesh->edgeAttribs.pArr[i].pData) {
			pContext->alloc.pFree(pMesh->edgeAttribs.pArr[i].pData);
		}
	}
	if (pMesh->edgeAttribs.count && pMesh->edgeAttribs.pArr) {
		pContext->alloc.pFree(pMesh->edgeAttribs.pArr);
	}
	for (int32_t i = 0; i < pMesh->vertAttribs.count; ++i) {
		if (pMesh->vertAttribs.pArr[i].pData) {
			pContext->alloc.pFree(pMesh->vertAttribs.pArr[i].pData);
		}
	}
	if (pMesh->vertAttribs.count && pMesh->vertAttribs.pArr) {
		pContext->alloc.pFree(pMesh->vertAttribs.pArr);
	}
	return STUC_SUCCESS;
}

StucResult stucGetAttribSize(StucAttrib *pAttrib, int32_t *pSize) {
	*pSize = getAttribSize(pAttrib->type);
	return STUC_SUCCESS;
}

StucResult stucGetAttrib(char *pName, AttribArray *pAttribs, Attrib **ppAttrib) {
	*ppAttrib = getAttrib(pName, pAttribs);
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
		verts[i] = *(V2_F32 *)(pMesh->pVerts + pMesh->mesh.pCorners[pFace->start + i]);
	}
	int8_t triCorners[4] = {0};
	V3_F32 bc = getBarycentricInFace(verts, triCorners, pFace->size, *pPos);
	if (bc.d[0] < -.0001f || bc.d[1] < -.0001f || bc.d[2] < -.0001f ||
		!v3IsFinite(bc)) {
		return;
	}
	V3_F32 vertsXyz[3];
	for (int32_t i = 0; i < 3; ++i) {
		int32_t vertIdx = pMesh->mesh.pCorners[pFace->start + triCorners[i]];
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
		stucGetCellsForSingleFace(&searchState, 1, &pos, &faceCellsTable, NULL, 0);
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
				face.start = pMesh->mesh.pFaces[face.idx];
				face.end = pMesh->mesh.pFaces[face.idx + 1];
				face.size = face.end - face.start;
				FaceTriangulated faceTris = {0};
				if (face.size > 4) {
					faceTris = triangulateFace(vars.pContext->alloc, face,
					                           pMesh->pVerts,
					                           pMesh->mesh.pCorners, 0);
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
	for (int32_t i = 0; i < pMesh->mesh.vertCount; ++i) {
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
	int32_t jobCount = activeJobs;
	pContext->threadPool.pJobStackPushJobs(pContext->pThreadPoolHandle,
	                                       activeJobs,
	                                       stucRenderJob, jobArgPtrs);
	CLOCK_INIT;
	CLOCK_START;
	waitForJobs(pContext, &activeJobs, pMutex);
	pContext->threadPool.pMutexDestroy(pContext->pThreadPoolHandle, pMutex);
	CLOCK_STOP("CREATE_IMAGE");
	int32_t pixelSize = getPixelSize(pImage->type);
	pImage->pData = pContext->alloc.pMalloc(pixelCount * pixelSize);
	for (int32_t i = 0; i < jobCount; ++i) {
		void *pImageOffset = offsetImagePtr(pImage, i * pixelsPerJob);
		int32_t bytesToCopy = pixelSize * args[i].pixelCount;
		memcpy(pImageOffset, args[i].imageBuf.pData, bytesToCopy);
		pContext->alloc.pFree(args[i].imageBuf.pData);
	}
	return STUC_SUCCESS;
}

void stucMapIndexedAttribsGet(StucContext pContext, StucMap pMap,
                              StucAttribIndexedArr *pIndexedAttribs) {
	*pIndexedAttribs = pMap->indexedAttribs;
}
