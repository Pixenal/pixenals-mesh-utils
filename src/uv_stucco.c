/* 
SPDX-FileCopyrightText: 2025 Caleb Dawson
SPDX-License-Identifier: Apache-2.0
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include <pixenals_alloc_utils.h>
#include <pixenals_thread_utils.h>
#include <pixenals_error_utils.h>

#include <io.h>
#include <attrib_utils.h>
#include <utils.h>
#include <interp_and_xform.h>
#include <merge_and_snap.h>

//TODO Reduce the bits written to the UVGP file for vert and corner indices, based on the total amount, in order to save space.
//	No point storing them as 32 bit if there's only like 4,000 verts
//TODO Split compressed data into chunks maybe?
//TODO add option to vary z projection depth with uv stretch (for wires and such)
//TODO Add an option for subdivision like smoothing (for instances where
//	the map is higher res than the base mesh). So that the surface can be
//	smoothed, without needing to induce the perf cost of actually subdividing
//	the base mesh. Is this possible?
//TODO add user define void * args to custom callbacks
//TODO allow for layered mapping. eg, map-faces assigned layer 0 are only mapped
//	to in-faces with a layer attribute of 0
//TODO make naming for MeshIn consistent

static
void setDefaultStageReport(StucContext pCtx) {
	pCtx->stageReport.outOf = 50,
	pCtx->stageReport.fpBegin = stucStageBegin;
	pCtx->stageReport.fpProgress = stucStageProgress;
	pCtx->stageReport.fpEnd = stucStageEnd;
}

StucErr stucContextInit(
	StucContext *pCtx,
	StucAlloc *pAlloc,
	StucThreadPool *pThreadPool,
	StucIo *pIo,
	StucTypeDefaultConfig *pTypeDefaultConfig,
	StucStageReport *pStageReport
) {
#ifndef NDEBUG
	stucIoDataTagValidate();
#endif

	StucAlloc alloc;
	if (pAlloc) {
		stucAllocSetCustom(&alloc, pAlloc);
	}
	else {
		stucAllocSetDefault(&alloc);
	}
	*pCtx = alloc.fpCalloc(1, sizeof(StucContextInternal));
	(*pCtx)->alloc = alloc;
	if (pThreadPool) {
		stucThreadPoolSetCustom(*pCtx, pThreadPool);
	}
	else {
		stucThreadPoolSetDefault(*pCtx);
	}
	if (pIo) {
		stucIoSetCustom(*pCtx, pIo);
	}
	else {
		stucIoSetDefault(*pCtx);
	}
	(*pCtx)->threadPool.fpInit(
		&(*pCtx)->pThreadPoolHandle,
		&(*pCtx)->threadCount,
		&(*pCtx)->alloc
	);
	if (pTypeDefaultConfig) {
		(*pCtx)->typeDefaults = *pTypeDefaultConfig;
	}
	else {
		stucSetTypeDefaultConfig(*pCtx);
	}
	if (pStageReport) {
		(*pCtx)->stageReport = *pStageReport;
	}
	else {
		setDefaultStageReport(*pCtx);
	}
	//TODO add ability to set custom specialAttrib names
	stucSetDefaultSpAttribNames(*pCtx);
	stucSetDefaultSpAttribDomains(*pCtx);
	stucSetDefaultSpAttribTypes(*pCtx);
	return PIX_ERR_SUCCESS;
}

StucErr stucContextDestroy(StucContext pCtx) {
	pCtx->threadPool.fpDestroy(pCtx->pThreadPoolHandle);
	pCtx->alloc.fpFree(pCtx);
	return PIX_ERR_SUCCESS;
}

//TODO replace these with StucUsg and StucObj arr structs, that combine arr and count
StucErr stucMapFileLoadForEdit(
	StucContext pCtx,
	const char *filePath,
	I32 *pObjCount,
	StucObject **ppObjArr,
	I32 *pUsgCount,
	StucUsg **ppUsgArr,
	I32 *pFlatCutoffCount,
	StucObject **ppFlatCutoffArr,
	StucAttribIndexedArr *pIndexedAttribs
) {
	//TODO reimplement
	return PIX_ERR_ERROR;
}

static
void buildEdgeLenList(StucContext pCtx, Mesh *pMesh) {
	PIX_ERR_ASSERT("", pMesh->pEdgeLen);
	V3_F32 *pPosCache = pCtx->alloc.fpMalloc(pMesh->core.edgeCount * sizeof(V3_F32));
	I8 *pSet = pCtx->alloc.fpCalloc(pMesh->core.edgeCount, 1);
	for (I32 i = 0; i < pMesh->core.cornerCount; ++i) {
		V3_F32 pos = pMesh->pPos[pMesh->core.pCorners[i]];
		I32 edge = pMesh->core.pEdges[i];
		if (!pSet[edge]) {
			pPosCache[edge] = pos;
			pSet[edge] = 1;
			continue;
		}
		//this occurs usually due to interior faces,
		// it shouldn't be an issue for for map-meshes, more so for in-meshes.
		//TODO remove this assert if no issues arise
		//PIX_ERR_ASSERT("more than 2 corners refernce 1 edge", pSet[edge] < 2);
		V3_F32 diff = _(pos V3SUB pPosCache[edge]);
		pMesh->pEdgeLen[edge] = pixmV3F32Len(diff);
		pSet[edge]++;
	}
	pCtx->alloc.fpFree(pSet);
	pCtx->alloc.fpFree(pPosCache);
}

static
StucErr attemptToSetMissingActiveDomains(StucMesh *pMesh) {
	StucErr err = PIX_ERR_SUCCESS;
	for (I32 i = 1; i < STUC_ATTRIB_USE_ENUM_COUNT; ++i) {
		if (i == STUC_ATTRIB_USE_SP_ENUM_COUNT) {
			continue;
		}
		AttribActive *pIdx = pMesh->activeAttribs + i;
		if (pIdx->domain != STUC_DOMAIN_NONE) {
			continue;
		}
		for (I32 j = STUC_DOMAIN_FACE; j <= STUC_DOMAIN_VERT; ++j) {
			const AttribArray *pAttribArr = stucGetAttribArrFromDomainConst(pMesh, j);
			if (pIdx->idx >= pAttribArr->count ||
				pAttribArr->pArr[pIdx->idx].core.use != i
			) {
				continue;
			}
			//the below is false, 2 domains have their own candidate.
			//the intended attrib is ambiguous, so return error
			PIX_ERR_RETURN_IFNOT_COND(
				err,
				pIdx->domain == STUC_DOMAIN_NONE,
				"Unable to determine active attrib domain"
			);
			pIdx->domain = j;
		}
	}
	return err;
}

static
void triCacheBuild(const StucAlloc *pAlloc, StucMap pMap) {
	bool ngons = checkForNgonsInMesh(&pMap->pMesh->core);
	if (ngons) {
		pMap->triCache.pArr =
			pAlloc->fpCalloc(pMap->pMesh->core.faceCount, sizeof(FaceTriangulated));
		pixalcLinAllocInit(pAlloc, &pMap->triCache.alloc, 3, 16, false);
		for (I32 i = 0; i < pMap->pMesh->core.faceCount; ++i) {
			FaceRange face = stucGetFaceRange(&pMap->pMesh->core, i);
			if (face.size > 4) {
				FaceTriangulated *pTris = pMap->triCache.pArr + i;
				pixalcLinAlloc(&pMap->triCache.alloc, (void **)&pTris->pTris, face.size - 2);
				stucTriangulateFaceFromVerts(pAlloc, &face, pMap->pMesh, pTris);
			}
		}
	}
}

static
void triCacheDestroy(const StucAlloc *pAlloc, StucMap pMap) {
	PIX_ERR_ASSERT("", !((pMap->triCache.pArr != NULL) ^ (pMap->triCache.alloc.valid)));
	if (pMap->triCache.pArr) {
		pAlloc->fpFree(pMap->triCache.pArr);
		pixalcLinAllocDestroy(&pMap->triCache.alloc);
		pMap->triCache = (TriCache) {0};
	}
}

static
void buildFaceBBoxes(const StucAlloc *pAlloc, StucMap pMap) {
	const Mesh *pMesh = pMap->pMesh;
	pMap->pFaceBBoxes = pAlloc->fpMalloc(pMesh->core.faceCount * sizeof(BBox));
	for (I32 i = 0; i < pMesh->core.faceCount; ++i) {
		FaceRange face = stucGetFaceRange(&pMesh->core, i);
		pMap->pFaceBBoxes[i] = stucBBoxGet(pMesh, &face);
	}
}

StucErr initFlatCutoff(
	StucContext pCtx,
	Usg *pUsg,
	StucObject *pCutoffObj
) {
	StucErr err = PIX_ERR_SUCCESS;
	pUsg->pFlatCutoff = pCtx->alloc.fpCalloc(1, sizeof(Mesh));
	pUsg->pFlatCutoff->core = *(StucMesh *)pCutoffObj->pData;

	err = attemptToSetMissingActiveDomains(&pUsg->pFlatCutoff->core);
	PIX_ERR_RETURN_IFNOT(err, "");
	err = stucAssignActiveAliases(
		pCtx,
		pUsg->pFlatCutoff,
		0x1 << STUC_ATTRIB_USE_POS,
		STUC_DOMAIN_NONE
	);
	PIX_ERR_RETURN_IFNOT(err, "");
	stucApplyObjTransform(
		&(StucObject){
			.pData = (StucObjectData *)pUsg->pFlatCutoff,
			.transform = pCutoffObj->transform
		}
	);
	return err;
}

struct MapDepEntry;

typedef struct MapDepPtrs {
	struct MapDepEntry **pArr;
	I32 size;
	I32 count;
} MapDepPtrs;

typedef struct MapDepEntry {
	HTableEntryCore core;
	MapDepPtrs deps;
	StucMap pMap;
	char *pName;
	char *pPath;
	bool onStack;
} MapDepEntry;

typedef struct MapDepStackEntry {
	MapDepEntry *pMap;
	I32 depIdx;
} MapDepStackEntry;

#define MAP_DEP_STACK_SIZE 64

typedef struct MapDepStack {
	MapDepStackEntry stack[MAP_DEP_STACK_SIZE];
	I32 ptr;
} MapDepStack;

typedef struct StrWithLen {
	const char *pStr;
	I32 len;
} StrWithLen;

static
void mapDepInit(
	void *pUserData,
	HTableEntryCore *pEntryCore,
	const void *pKeyData,
	void *pInitInfo,
	I32 linIdx
) {
	MapDepEntry *pEntry = (MapDepEntry *)pEntryCore;
	char **ppStrArr = pInitInfo;
	pEntry->pName = ppStrArr[0];
	pEntry->pPath = ppStrArr[1];
}

static
bool mapDepCmp(
	const HTableEntryCore *pEntryCore,
	const void *pKeyData,
	const void *pInitInfo
) {
	const MapDepEntry *pEntry = (MapDepEntry *)pEntryCore;
	const StrWithLen *pKey = pKeyData;
	return !strncmp(pEntry->pName, pKey->pStr, pixioPathMaxGet());
}

static
StucKey mapDepMakeKey(const void *pKeyData) {
	const StrWithLen *pStr = pKeyData;
	return (StucKey){.pKey = pStr->pStr, .size = pStr->len};
}

static
 SearchResult addMapDepEntry(
	 const StucAlloc *pAlloc,
	 HTable *pTable,
	 const char *pFilepath,
	 MapDepEntry **ppEntry
 ) {
	I32 nameLen = 0;
	I32 pathLen = 0;
	const char *pName = stucGetBasename(pFilepath, &nameLen, &pathLen);
	char *pNameCpy = pAlloc->fpMalloc(nameLen + 1);
	char *pPathCpy = pAlloc->fpMalloc(pathLen + 1);
	memcpy(pNameCpy, pName, nameLen + 1);
	memcpy(pPathCpy, pFilepath, pathLen + 1);
	return stucHTableGet(
		pTable,
		0,
		&(StrWithLen){.pStr = pNameCpy, .len = nameLen},
		ppEntry,
		true,
		(char *[]){pNameCpy, pPathCpy},
		mapDepMakeKey, NULL, mapDepInit, mapDepCmp
	);
}

static
StucErr addMapEntryDeps(
	const StucAlloc *pAlloc,
	MapDepStack *pStack,
	HTable *pTable,
	MapDepEntry *pEntry,
	StucMapDeps *pDeps
) {
	StucErr err = PIX_ERR_SUCCESS;
	for (I32 i = 0; i < pDeps->maps.count; ++i) {
		I32 nameLen = strnlen(pDeps->maps.pArr[i].pStr, pixioPathMaxGet());
		I32 idx = 0;
		PIXALC_DYN_ARR_ADD(void *, pAlloc, &pEntry->deps, idx);
		SearchResult result = addMapDepEntry(
			pAlloc,
			pTable,
			pDeps->maps.pArr[i].pStr,
			pEntry->deps.pArr + idx
		);
		PIX_ERR_RETURN_IFNOT_COND(
			err,
			pEntry->deps.pArr[idx]->onStack,
			"circular map dependency"
		);
	}
	return err;
}

static
StucErr mapDepStackPush(MapDepStack *pStack, MapDepEntry *pEntry) {
	StucErr err = PIX_ERR_SUCCESS;
	++pStack->ptr;
	PIX_ERR_RETURN_IFNOT_COND(
		err,
		pStack->ptr < MAP_DEP_STACK_SIZE,
		"too many map dependencies"
	);
	pStack->stack[pStack->ptr] = (MapDepStackEntry){.pMap = pEntry};
	return err;
}

static
StucErr mapDepStackPop(MapDepStack *pStack) {
	StucErr err = PIX_ERR_SUCCESS;
	PIX_ERR_RETURN_IFNOT_COND(err, pStack->ptr >= 0, "");
	--pStack->ptr;
	return err;
}

StucErr stucMapFileLoad(
	StucContext pCtx,
	const char *filePath,
	PixErr (* fpMapGet)(const char *, const char **, StucMap * const),
	PixErr (* fpMapStore)(const char *, StucMap)
) {
	StucErr err = PIX_ERR_SUCCESS;
	PIX_ERR_RETURN_IFNOT_COND(err, filePath, "provided filepath was NULL");

	HTable table;
	stucHTableInit(
		&pCtx->alloc,
		&table,
		64,
		(I32Arr){.pArr = (I32[]){sizeof(MapDepEntry)}, .count = 1},
		NULL
	);
	MapDepStack stack = {0};
	do {
		MapDepStackEntry *pStackEntry = stack.stack + stack.ptr;
		if (!filePath) {
			StucMap pMap = NULL;
			err = fpMapGet(pStackEntry->pMap->pName, &filePath, &pMap);
			PIX_ERR_THROW_IFNOT(err, "", 0);
			if (pMap) {
				mapDepStackPop(&stack);
				continue;
			}
			PIX_ERR_THROW_IFNOT_COND(err, filePath, "not path provided", 0);
		}
		StucMapDeps deps = {0};
		err = stucMapImportGetDep(pCtx, filePath, &deps);
		PIX_ERR_THROW_IFNOT(err, "", 0);
		if (addMapDepEntry(&pCtx->alloc, &table, filePath, &pStackEntry->pMap) ==
			STUC_SEARCH_ADDED &&
			deps.maps.count
		) {
			err = addMapEntryDeps(&pCtx->alloc, &stack, &table, pStackEntry->pMap, &deps);
			PIX_ERR_THROW_IFNOT(err, "", 0);
		}
		else if (pStackEntry->depIdx < pStackEntry->pMap->deps.count) {
			++pStackEntry->depIdx;
		}
		else {
			StucMap pMap = NULL;
			err = stucMapFileLoadIntern(
				pCtx,
				&pMap,
				pStackEntry->pMap
			);
			PIX_ERR_THROW_IFNOT_COND(err, pMap, "", 0);
			fpMapStore(pMap->pName, pMap);
			err = mapDepStackPop(&stack);
			PIX_ERR_THROW_IFNOT(err, "", 0);
			continue;
		}
		err = mapDepStackPush(&stack, pStackEntry->pMap->deps.pArr[pStackEntry->depIdx]);
		PIX_ERR_THROW_IFNOT(err, "", 0);
	} while(filePath = NULL, stack.ptr >= 0);

	PIX_ERR_CATCH(0, err, ;);
	stucHTableDestroy(&table);
	return err;
}

StucErr stucMapFileLoadIntern(
	StucContext pCtx,
	StucMap *pMapHandle,
	MapDepEntry *pEntry
) {
	StucErr err = PIX_ERR_SUCCESS;
	PIX_ERR_RETURN_IFNOT_COND(err, pEntry->pName, "file name is too long");
	StucMap pMap = pCtx->alloc.fpCalloc(1, sizeof(MapFile));
	pMap->pName = pEntry->pName;
	pMap->pPath = pEntry->pPath;
	StucObjArr objArr = {0};
	StucUsgArr usgArr = {0};
	StucObjArr cutoffArr = {0};
	ObjMapOptsArr mapOptsArr = {0};
	err = stucMapImport(
		pCtx, pEntry->pPath,
		&objArr,
		&mapOptsArr,
		&usgArr,
		&cutoffArr,
		NULL,
		&pMap->indexedAttribs,
		true
	);
	//TODO validate meshes, ensure pMatIdx is within mat range, faces are within max corner limit,
	//F32 values are valid, etc.
	PIX_ERR_THROW_IFNOT(err, "failed to load file from disk", 0);

	I32 targetIdx = 0;
	for (I32 i = 0; i < objArr.count; ++i) {
		Mesh *pMesh = (Mesh *)objArr.pArr[i].pData;
		
		err = attemptToSetMissingActiveDomains(&pMesh->core);
		PIX_ERR_THROW_IFNOT(err, "", 0);

		if (targetIdx < mapOptsArr.count &&
			mapOptsArr.pArr[targetIdx].obj == i
		) {
			StucMesh meshOut = {0};
			AttribIndexedArr outIdxAttribArr= {0};
			//TODO add an option to stucMapToMesh which retains all
			//original indexed attribs, rather than only keeping those
			//that are used
			stucMapToMesh(
				pCtx,
				&mapOptsArr.pArr[targetIdx].arr,
				&pMesh->core,
				&pMap->indexedAttribs,
				&meshOut,
				&outIdxAttribArr,
				1.0f,//TODO replace with actual wScale an receiveLen vars
				-1.0f,
				true
			);
			stucAttribIndexedArrDestroy(pCtx, &pMap->indexedAttribs);
			pMap->indexedAttribs = outIdxAttribArr;
			stucMeshDestroy(pCtx, &pMesh->core);
			pMesh->core = meshOut;
			++targetIdx;
		}
		err = stucAssignActiveAliases(
			pCtx,
			pMesh,
			STUC_ATTRIB_USE_FIELD(((StucAttribUse[]) {
				STUC_ATTRIB_USE_POS,
				STUC_ATTRIB_USE_UV,
				STUC_ATTRIB_USE_NORMAL,
				STUC_ATTRIB_USE_RECEIVE,
				STUC_ATTRIB_USE_IDX
			})),
			STUC_DOMAIN_NONE
		);
		PIX_ERR_THROW_IFNOT(err, "", 0);
		stucApplyObjTransform(objArr.pArr + i);
	}
	Mesh *pMapMesh = pCtx->alloc.fpCalloc(1, sizeof(Mesh));
	pMapMesh->core.type.type = STUC_OBJECT_DATA_MESH_INTERN;
	err = stucMergeObjArr(pCtx, pMapMesh, &objArr, false);
	PIX_ERR_THROW_IFNOT(err, "", 0);

	UBitField32 spToAppend = STUC_ATTRIB_USE_FIELD(((StucAttribUse[]) {
		STUC_ATTRIB_USE_EDGE_LEN
	}));
	stucAppendSpAttribsToMesh(
		pCtx,
		pMapMesh,
		spToAppend | (usgArr.count ? 0x1 << STUC_ATTRIB_USE_USG : 0x0),
		STUC_ATTRIB_ORIGIN_MAP
	);

	stucSetAttribOrigins(&pMapMesh->core.meshAttribs, STUC_ATTRIB_ORIGIN_MAP);
	stucSetAttribOrigins(&pMapMesh->core.faceAttribs, STUC_ATTRIB_ORIGIN_MAP);
	stucSetAttribOrigins(&pMapMesh->core.cornerAttribs, STUC_ATTRIB_ORIGIN_MAP);
	stucSetAttribOrigins(&pMapMesh->core.edgeAttribs, STUC_ATTRIB_ORIGIN_MAP);
	stucSetAttribOrigins(&pMapMesh->core.vertAttribs, STUC_ATTRIB_ORIGIN_MAP);

	stucSetAttribCopyOpt(
		pCtx,
		&pMapMesh->core,
		STUC_ATTRIB_DONT_COPY,
		~STUC_ATTRIB_USE_FIELD(((StucAttribUse[]) { //all except for
			STUC_ATTRIB_USE_POS,
			STUC_ATTRIB_USE_UV,
			STUC_ATTRIB_USE_NORMAL,
			STUC_ATTRIB_USE_IDX
		}))
	);
	err = stucAssignActiveAliases(
		pCtx,
		pMapMesh,
		STUC_ATTRIB_USE_FIELD(((StucAttribUse[]) {
			STUC_ATTRIB_USE_POS,
			STUC_ATTRIB_USE_UV,
			STUC_ATTRIB_USE_NORMAL,
			STUC_ATTRIB_USE_RECEIVE,
			STUC_ATTRIB_USE_USG,
			STUC_ATTRIB_USE_IDX,
			STUC_ATTRIB_USE_EDGE_LEN,
			STUC_ATTRIB_USE_NONE
		})),
		STUC_DOMAIN_NONE
	);
	PIX_ERR_THROW_IFNOT(err, "", 0);

	buildEdgeLenList(pCtx, pMapMesh);

	//TODO some form of heap corruption when many objects
	//test with address sanitizer on CircuitPieces.stuc
	stucObjArrDestroy(pCtx, &objArr);

	//set corner attribs to interpolate by default
	//TODO make this an option in ui, even for non common attribs
	for (I32 i = 0; i < pMapMesh->core.cornerAttribs.count; ++i) {
		pMapMesh->core.cornerAttribs.pArr[i].interpolate = true;
	}

	pMap->pMesh = pMapMesh;

	triCacheBuild(&pCtx->alloc, pMap);
	buildFaceBBoxes(&pCtx->alloc, pMap);

	//the quadtree is created before USGs are assigned to verts,
	//as the tree's used to speed up the process
	printf("File loaded. Creating quad tree\n");
	err = stucCreateQuadTree(pCtx, &pMap->quadTree, pMap->pMesh, pMap->pFaceBBoxes);
	PIX_ERR_THROW_IFNOT(err, "failed to create quadtree", 0);

	if (usgArr.count) {
		pMap->usgArr.count = usgArr.count;
		pMap->usgArr.pArr = pCtx->alloc.fpCalloc(pMap->usgArr.count, sizeof(Usg));
		for (I32 i = 0; i < pMap->usgArr.count; ++i) {
			Mesh *pUsgMesh = (Mesh *)usgArr.pArr[i].obj.pData;
			err = attemptToSetMissingActiveDomains(&pUsgMesh->core);
			PIX_ERR_THROW_IFNOT(err, "", 0);
			err = stucAssignActiveAliases(
				pCtx,
				pUsgMesh,
				0x1 << STUC_ATTRIB_USE_POS,
				STUC_DOMAIN_NONE
			);
			PIX_ERR_THROW_IFNOT(err, "", 0);
			Usg *pUsg = pMap->usgArr.pArr + i;
			pUsg->origin = *(V2_F32 *)&usgArr.pArr[i].obj.transform.d[3];
			pUsg->pMesh = pUsgMesh;
			stucApplyObjTransform(&usgArr.pArr[i].obj);
			if (usgArr.pArr[i].flatCutoff.enabled) {
				//TODO these shouldn't be duplicated for each usg,
				//store cutoffs in a separate arr
				I32 cutoffIdx = usgArr.pArr[i].flatCutoff.idx;
				initFlatCutoff(pCtx, pUsg, cutoffArr.pArr + cutoffIdx);
			}
		}
		Mesh *pSquares = pCtx->alloc.fpCalloc(1, sizeof(Mesh));
		stucAllocUsgSquaresMesh(pCtx, pMap, pSquares);
		stucFillUsgSquaresMesh(pMap, usgArr.pArr, pSquares);
		pMap->usgArr.pSquares = pSquares;
		stucAssignUsgsToVerts(&pCtx->alloc, pMap, usgArr.pArr);
		pMap->usgArr.pMemArr = usgArr.pArr;
	}

	*pMapHandle = pMap;
	PIX_ERR_CATCH(0, err, stucMapFileUnload(pCtx, pMap);)

	return err;
}

StucErr stucMapFileUnload(StucContext pCtx, StucMap pMap) {
	stucDestroyQuadTree(pCtx, &pMap->quadTree);
	if (pMap->pMesh) {
		stucMeshDestroy(pCtx, (StucMesh *)&pMap->pMesh->core);
		pCtx->alloc.fpFree((Mesh *)pMap->pMesh);
	}
	triCacheDestroy(&pCtx->alloc, pMap);
	if (pMap->pFaceBBoxes) {
		pCtx->alloc.fpFree(pMap->pFaceBBoxes);
	}
	if (pMap->usgArr.pSquares) {
		pCtx->alloc.fpFree((Mesh *)pMap->usgArr.pSquares);
	}
	pCtx->alloc.fpFree(pMap);
	return PIX_ERR_SUCCESS;
}

StucErr stucMapFileMeshGet(StucContext pCtx, StucMap pMap, const StucMesh **ppMesh) {
	StucErr err = PIX_ERR_SUCCESS;
	PIX_ERR_RETURN_IFNOT_COND(err, pCtx && pMap && ppMesh, "invalid args");
	*ppMesh = &pMap->pMesh->core;
	return err;
}

StucErr stucMapNameGet(StucContext pCtx, StucMap pMap, const char **ppName) {
	StucErr err = PIX_ERR_SUCCESS;
	PIX_ERR_RETURN_IFNOT_COND(err, pCtx && pMap && ppName, "");
	*ppName = pMap->pName;
	return err;
}

StucErr stucMapPathGet(StucContext pCtx, StucMap pMap, const char **ppPath) {
	StucErr err = PIX_ERR_SUCCESS;
	PIX_ERR_RETURN_IFNOT_COND(err, pCtx && pMap && ppPath, "");
	*ppPath = pMap->pPath;
	return err;
}

static
void initBlendOpt(
	StucContext pCtx,
	StucBlendOpt *pEntry,
	const StucAttrib *pAttrib,
	I32 attribIdx
) {
	StucTypeDefault *pDefault = 
		stucGetTypeDefaultConfig(&pCtx->typeDefaults, pAttrib->core.type);
	pEntry->blendConfig = pDefault->blendConfig;
	pEntry->attrib = attribIdx;
}

static
StucErr getCommonAttribs(
	StucContext pCtx,
	const StucMesh *pMapMesh,
	const AttribArray *pMapAttribs,
	const StucMesh *pMesh,
	const AttribArray *pMeshAttribs,
	StucBlendOptArr *pOptArr
) {
	StucErr err = PIX_ERR_SUCCESS;
	//TODO ignore special attribs like StucTangent or StucTSign
	for (I32 i = 0; i < pMeshAttribs->count; ++i) {
		Attrib *pAttrib = pMeshAttribs->pArr + i;
		if (pAttrib->core.use == STUC_ATTRIB_USE_POS) {
			continue;
		}
		const Attrib *pMapAttrib = NULL;
		err = stucGetMatchingAttribConst(
			pCtx,
			pMapMesh, pMapAttribs,
			pMesh, pAttrib,
			true,
			true,
			&pMapAttrib
		);
		PIX_ERR_THROW_IFNOT(err, "", 0);
		if (!pMapAttrib) {
			continue;
		}
		I32 newIdx = 0;
		PIXALC_DYN_ARR_ADD(StucBlendOpt, &pCtx->alloc, pOptArr, newIdx);
		initBlendOpt(pCtx, pOptArr->pArr + newIdx, pAttrib, i);
	}
	PIX_ERR_CATCH(0, err,
		pCtx->alloc.fpFree(pOptArr->pArr);
		*pOptArr = (StucBlendOptArr){0};
	);
	return err;
}

//TODO handle edge case, where attribute share the same name,
//but have incompatible types. Such as a F32 and a string.
StucErr stucQueryCommonAttribs(
	StucContext pCtx,
	const StucMap pMap,
	const StucMesh *pMesh,
	StucBlendOptArr *pOptArr
) {
	StucErr err = PIX_ERR_SUCCESS;
	PIX_ERR_RETURN_IFNOT_COND(err, pCtx && pMap && pMesh && pOptArr, "");
	const StucMesh *pMapMesh = &pMap->pMesh->core;
	StucMesh meshWrap = *pMesh;
	err = attemptToSetMissingActiveDomains(&meshWrap);
	PIX_ERR_RETURN_IFNOT(err, "");
	for (I32 domain = STUC_DOMAIN_FACE; domain <= STUC_DOMAIN_MESH; ++domain) {
		err = getCommonAttribs(
			pCtx,
			pMapMesh,
			stucGetAttribArrFromDomainConst(pMapMesh, domain),
			&meshWrap,
			stucGetAttribArrFromDomainConst(&meshWrap, domain),
			pOptArr + domain
		);
		PIX_ERR_THROW_IFNOT(err, "", 0);
	}
	PIX_ERR_CATCH(0, err,
		stucDestroyBlendOptArr(pCtx, pOptArr);
	);
	return err;
}

StucErr stucDestroyBlendOptArr(
	StucContext pCtx,
	StucBlendOptArr *pOptArr
) {
	StucErr err = PIX_ERR_SUCCESS;
	PIX_ERR_RETURN_IFNOT_COND(err, pCtx && pOptArr, "");
	for (I32 domain = STUC_DOMAIN_FACE; domain <= STUC_DOMAIN_MESH; ++domain) {
		StucBlendOptArr *pArr = pOptArr + domain;
		if (pArr) {
			pCtx->alloc.fpFree(pArr->pArr);
			pArr->pArr = NULL;
		}
		pArr->count = pArr->size = 0;
	}
	return err;
}

static
void buildEdgeAdj(Mesh *pMesh) {
	const StucMesh *pCore = &pMesh->core;
	memset(pMesh->pEdgeFaces, -1, sizeof(V2_I32) * pCore->edgeCount);
	memset(pMesh->pEdgeCorners, -1, sizeof(V2_I8) * pCore->edgeCount);
	for (I32 i = 0; i < pCore->faceCount; ++i) {
		FaceRange face = stucGetFaceRange(&pMesh->core, i);
		for (I32 j = 0; j < face.size; ++j) {
			I32 edge = pCore->pEdges[face.start + j];
			bool which = pMesh->pEdgeFaces[edge].d[0] >= 0;
			pMesh->pEdgeFaces[edge].d[which] = i;
			pMesh->pEdgeCorners[edge].d[which] = j;
		}
	}
}

static
void incNumAdjSeam(const Mesh *pMesh, I32 vert) {
	I32 numSeam = pMesh->pNumAdjPreserve[vert] & 0xf;
	if (numSeam < 3) {
		numSeam++;
		pMesh->pNumAdjPreserve[vert] &= 0xf0;
		pMesh->pNumAdjPreserve[vert] |= numSeam;
	}
}

static
void incNumAdjPreserve(const Mesh *pMesh, I32 vert) {
	I32 numPreserve = pMesh->pNumAdjPreserve[vert] >> 4 & 0xf;
	if (numPreserve < 3) { //only record up to 3
		numPreserve++;
		pMesh->pNumAdjPreserve[vert] &= 0xf;
		pMesh->pNumAdjPreserve[vert] |= numPreserve << 4;
	}
}

static
void buildSeamAndPreserveTables(Mesh *pMesh) {
	for (I32 i = 0; i < pMesh->core.edgeCount; ++i) {
		bool seam = stucIsEdgeSeam(pMesh, i);
		bool preserve = stucGetIfPreserveEdge(pMesh, i);
		if (seam || preserve) {
			V2_I32 faces = pMesh->pEdgeFaces[i];
			V2_I8 corners = pMesh->pEdgeCorners[i];
			I32 vert = stucGetMeshVert(
				&pMesh->core,
				(FaceCorner) {.face = faces.d[0], .corner = corners.d[0]}
			);
			if (seam) {
				pMesh->pSeamEdge[i] = seam;
				incNumAdjSeam(pMesh, vert);
			}
			else if (preserve) {
				incNumAdjPreserve(pMesh, vert);
			}
			if (faces.d[1] >= 0) {
				vert = stucGetMeshVert(
					&pMesh->core,
					(FaceCorner) {
					.face = faces.d[1], .corner = corners.d[1]
				}
				);
				if (seam) {
					incNumAdjSeam(pMesh, vert);
				}
				else if (preserve) {
					incNumAdjPreserve(pMesh, vert);
				}
			}
		}
	}
}

static
bool checkIfNoFacesHaveMaskIdx(const Mesh *pMesh, I8 maskIdx) {
	if (!pMesh->pMatIdx) {
		return false;
	}
	for (I32 i = 0; i < pMesh->core.faceCount; ++i) {
		if (pMesh->pMatIdx[i] == maskIdx) {
			return false;
		}
	}
	return true;
}

static
StucErr mapToMeshInternal(
	StucContext pCtx,
	const StucMap pMap,
	Mesh *pMeshIn,
	StucMesh *pOutMesh,
	I8 maskIdx,
	const StucBlendOptArr *pOptArr,
	InFaceTable *pInFaceTable,
	F32 wScale,
	F32 receiveLen
) {
	StucErr err = PIX_ERR_SUCCESS;
	if (checkIfNoFacesHaveMaskIdx(pMeshIn, maskIdx)) {
		return err;
	}
	MapToMeshBasic basic = {
		.pCtx = pCtx,
		.pMap = pMap,
		.pInMesh = pMeshIn,
		.pOptArr = pOptArr,
		.wScale = wScale,
		.receiveLen = receiveLen,
		.maskIdx = maskIdx,
		.pInFaceTable = pInFaceTable,
	};
	printf("A\n");
	if (pInFaceTable) {
		pixalcLinAllocInit(
			&pCtx->alloc,
			&pInFaceTable->alloc,
			sizeof(I32),
			pMeshIn->core.faceCount,
			true
		);
	}
	bool empty = false;
	InPieceArr inPieceArr = {0};
	I32 findEncasedJobCount = 0;
	FindEncasedFacesJobArgs findEncasedJobArgs[PIX_THREAD_MAX_SUB_MAPPING_JOBS] = {0};
	err = stucInPieceArrInit(
		&basic,
		&inPieceArr,
		&findEncasedJobCount, findEncasedJobArgs,
		&empty
	);
	PIX_ERR_RETURN_IFNOT(err, "");
	printf("B\n");
	if (!empty) {
		BufMeshArr bufMeshes = {0};
		BufMeshArr bufMeshesClip = {0};
		InPieceArr inPiecesSplit = {.pBufMeshes = &bufMeshes};
		InPieceArr inPiecesSplitClip = {.pBufMeshes = &bufMeshesClip};
		SplitInPiecesAllocArr splitAlloc = {
			.pArr = (SplitInPiecesAlloc[PIX_THREAD_MAX_SUB_MAPPING_JOBS]){0}
		};
		BufOutRangeTable bufOutTable = {0};
		OutBufIdxArr outBufIdxArr = {0};
		err = stucInPieceArrSplit(
			&basic,
			&inPieceArr,
			&inPiecesSplit, &inPiecesSplitClip,
			&splitAlloc
		);
		for (I32 i = 0; i < findEncasedJobCount; ++i) {
			PixalcLinAlloc *pEncasedAlloc =
				stucHTableAllocGet(&findEncasedJobArgs[i].encasedFaces, 0);
			PixalcLinAllocIter iter = {0};
			pixalcLinAllocIterInit(pEncasedAlloc, (Range) {0, INT32_MAX}, &iter);
			for (; !pixalcLinAllocIterAtEnd(&iter); pixalcLinAllocIterInc(&iter)) {
				EncasedMapFace *pEntry = pixalcLinAllocGetItem(&iter);
				if (pEntry->inFaces.pArr) {
					basic.pCtx->alloc.fpFree(pEntry->inFaces.pArr);
					pEntry->inFaces.pArr = NULL;
				}
			}
			stucHTableDestroy(&findEncasedJobArgs[i].encasedFaces);
		}
		PIX_ERR_RETURN_IFNOT(err, "");
		printf("C\n");
		
		err = stucInPieceArrInitBufMeshes(&basic, &inPiecesSplitClip, stucClipMapFace);
		PIX_ERR_RETURN_IFNOT(err, "");
		err = stucInPieceArrInitBufMeshes(&basic, &inPiecesSplit, stucAddMapFaceToBufMesh);
		PIX_ERR_RETURN_IFNOT(err, "");
		printf("D\n");

		HTable mergeTable = {0};
		stucVertMergeTableInit(&basic, &inPiecesSplit, &inPiecesSplitClip, &mergeTable);
		stucMergeVerts(&basic, &inPiecesSplit, false, &mergeTable);
		stucMergeVerts(&basic, &inPiecesSplitClip, true, &mergeTable);
		printf("E\n");

		I32 snappedVerts = 0;
		err = stucSnapIntersectVerts(
			&basic,
			&inPiecesSplit, &inPiecesSplitClip,
			&mergeTable,
			&snappedVerts
		);
		PIX_ERR_RETURN_IFNOT(err, "");
		printf("F\n");

		stucInitOutMesh(&basic, &mergeTable, snappedVerts);
		stucAddVertsToOutMesh(&basic, &mergeTable, 0);
		stucAddVertsToOutMesh(&basic, &mergeTable, 1);//intersect verts
		bufOutTable.size =
			inPiecesSplit.pBufMeshes->count + inPiecesSplitClip.pBufMeshes->count;
		bufOutTable.pArr = pCtx->alloc.fpCalloc(bufOutTable.size, sizeof(BufOutRange));
		stucAddFacesAndCornersToOutMesh(
			&basic,
			&inPiecesSplit,
			&mergeTable,
			&outBufIdxArr,
			&bufOutTable,
			false
		);
		stucAddFacesAndCornersToOutMesh(
			&basic,
			&inPiecesSplitClip,
			&mergeTable,
			&outBufIdxArr,
			&bufOutTable,
			true
		);
		if (!basic.outMesh.core.faceCount) {
			goto cleanUp;
		}
		stucMeshSetLastFace(pCtx, &basic.outMesh);
		printf("G\n");

		err = stucBuildTangentsForInPieces(
			&basic,
			pMeshIn,
			&inPiecesSplit, &inPiecesSplitClip,
			&mergeTable
		);
		PIX_ERR_RETURN_IFNOT(err, "");
		printf("H\n");
		
		err = stucXFormAndInterpVerts(&basic, &inPiecesSplit, &inPiecesSplitClip, &mergeTable, 0);
		PIX_ERR_RETURN_IFNOT(err, "");
		//intersect verts
		err = stucXFormAndInterpVerts(&basic, &inPiecesSplit, &inPiecesSplitClip, &mergeTable, 1);
		PIX_ERR_RETURN_IFNOT(err, "");
		err = stucInterpAttribs(
			&basic,
			&inPiecesSplit, &inPiecesSplitClip,
			&mergeTable,
			&bufOutTable,
			&outBufIdxArr,
			STUC_DOMAIN_FACE, stucInterpFaceAttribs
		);
		PIX_ERR_RETURN_IFNOT(err, "");
		//vert merge lin-idx is replaced with out-vert idx in corner-interp job,
		// so faces must be interpolated before corners
		err = stucInterpAttribs(
			&basic,
			&inPiecesSplit, &inPiecesSplitClip,
			&mergeTable,
			&bufOutTable,
			&outBufIdxArr,
			STUC_DOMAIN_CORNER, stucInterpCornerAttribs
		);
		PIX_ERR_RETURN_IFNOT(err, "");
		printf("I\n");

		stucReallocMeshToFit(pCtx, &basic.outMesh);
		*pOutMesh = basic.outMesh.core;
		printf("J\n");
	cleanUp:
		if (outBufIdxArr.pArr) {
			pCtx->alloc.fpFree(outBufIdxArr.pArr);
		}
		if (bufOutTable.pArr) {
			pCtx->alloc.fpFree(bufOutTable.pArr);
		}
		for (I32 i = 0; i < splitAlloc.count; ++i) {
			if (splitAlloc.pArr[i].encased.valid) {
				pixalcLinAllocDestroy(&splitAlloc.pArr[i].encased);
			}
			if (splitAlloc.pArr[i].inFace.valid) {
				pixalcLinAllocDestroy(&splitAlloc.pArr[i].inFace);
			}
			if (splitAlloc.pArr[i].border.valid) {
				pixalcLinAllocDestroy(&splitAlloc.pArr[i].border);
			}
		}
		stucHTableDestroy(&mergeTable);
		inPieceArrDestroy(pCtx, &inPiecesSplit);
		inPieceArrDestroy(pCtx, &inPiecesSplitClip);
		stucBufMeshArrDestroy(pCtx, &bufMeshes);
		stucBufMeshArrDestroy(pCtx, &bufMeshesClip);
		printf("K\n");
	}
	PIX_ERR_CATCH(0, err, ;);
	return err;
}

static
void addEntryToInFaceTable(
	const StucAlloc *pAlloc,
	UsgInFace **ppHashTable,
	StucMap pMap,
	InFaceArr *pInFaceTable,
	I32 squareIdx,
	I32 inFaceIdx
) {
	U32 sum = pInFaceTable[squareIdx].usg + pInFaceTable[squareIdx].pArr[inFaceIdx];
	I32 hash = stucFnvHash((U8 *)&sum, sizeof(sum), pMap->usgArr.tableSize);
	UsgInFace *pEntry = *ppHashTable + hash;
	if (!pEntry->pEntry) {
		pEntry->pEntry = pInFaceTable + squareIdx;
		pEntry->face = pInFaceTable[squareIdx].pArr[inFaceIdx];
		return;
	}
	do {
		if (!pEntry->pNext) {
			pEntry = pEntry->pNext = pAlloc->fpCalloc(1, sizeof(UsgInFace));
			pEntry->pEntry = pInFaceTable + squareIdx;
			pEntry->face = pInFaceTable[squareIdx].pArr[inFaceIdx];
			break;
		}
		pEntry = pEntry->pNext;
	} while(true);
}

static
void InFaceTableToHashTable(
	const StucAlloc *pAlloc,
	StucMap pMap,
	I32 count,
	InFaceArr *pInFaceTable
) {
	UsgInFace **ppHashTable = &pMap->usgArr.pInFaceTable;
	pMap->usgArr.tableSize = count * 2;
	*ppHashTable = pAlloc->fpCalloc(pMap->usgArr.tableSize, sizeof(UsgInFace));
	for (I32 i = 0; i < count; ++i) {
		for (I32 j = 0; j < pInFaceTable[i].count; ++j) {
			addEntryToInFaceTable(pAlloc, ppHashTable, pMap, pInFaceTable, i, j);
		}
	}
}

static
StucErr getOriginIndexedAttrib(
	StucContext pCtx,
	Attrib *pAttrib,
	I32 attribIdx,
	const StucMapArr *pMapArr,
	I32 mapIdx,
	const AttribIndexed *pMapIndexedAttrib,
	const AttribIndexed *pInIndexedAttrib,
	const AttribIndexed **ppMatsToAdd,
	StucDomain domain
) {
	StucErr err = PIX_ERR_SUCCESS;
	switch (pAttrib->origin) {
		case STUC_ATTRIB_ORIGIN_MAP:
			*ppMatsToAdd = pMapIndexedAttrib;
			break;
		case STUC_ATTRIB_ORIGIN_MESH_IN:
			*ppMatsToAdd = pInIndexedAttrib;
			break;
		case STUC_ATTRIB_ORIGIN_COMMON: {
			const BlendOptArr *pOptArr = pMapArr->pArr[mapIdx].blendOptArr;
			const BlendOpt *pOpts = stucGetBlendOpt(pOptArr, attribIdx, domain);
			BlendConfig config = {0};
			if (pOpts) {
				config = pOpts->blendConfig;
			}
			else {
				StucTypeDefault *pDefaultConfig =
					stucGetTypeDefaultConfig(&pCtx->typeDefaults, STUC_ATTRIB_STRING);
				config = pDefaultConfig->blendConfig;
			}
			*ppMatsToAdd = config.order ? pInIndexedAttrib : pMapIndexedAttrib;
			break;
		}
		default:
			PIX_ERR_ASSERT("invalid attrib origin for this function", false);
	}
	return err;
}

static
StucErr iterFacesAndCorrectIdxAttrib(
	StucContext pCtx,
	Attrib *pAttrib,
	Mesh *pMesh,
	AttribIndexed *pOutIndexedAttrib,
	const AttribIndexed *pOriginIndexedAttrib,
	StucDomain domain
) {
	StucErr err = PIX_ERR_SUCCESS;
	PIX_ERR_RETURN_IFNOT_COND(err, pAttrib->core.type == STUC_ATTRIB_I8, "");

	typedef struct TableEntry {
		I8 idx;
		bool hasRef;
	} TableEntry;
	TableEntry *pTable =
		pCtx->alloc.fpCalloc(pOriginIndexedAttrib->count, sizeof(TableEntry));

	I8 *pIdx = pAttrib->core.pData;
	I32 domainCount = stucDomainCountGetIntern(&pMesh->core, domain);
	for (I32 i = 0; i < domainCount; ++i) {
		I32 idx = pIdx[i];
		PIX_ERR_THROW_IFNOT_COND(err, idx >= 0 && idx < pOriginIndexedAttrib->count, "", 0);
		TableEntry *pEntry = pTable + idx;
		if (!pEntry->hasRef) {
			pEntry->hasRef = true;
			I32 newIdx = stucGetIdxInIndexedAttrib(
				pOutIndexedAttrib,
				pOriginIndexedAttrib,
				idx
			);
			if (newIdx >= 0) {
				pEntry->idx = (I8)newIdx;
			}
			else {
				pEntry->idx = (I8)pOutIndexedAttrib->count;
				stucAppendToIndexedAttrib(
					pCtx,
					pOutIndexedAttrib,
					&pOriginIndexedAttrib->core,
					idx
				);
			}
		}
		pIdx[i] = pEntry->idx;
	}
	PIX_ERR_CATCH(0, err, ;);
	pCtx->alloc.fpFree(pTable);
	return err;
}

static
StucErr getIndexedAttribInMaps(
	StucContext pCtx,
	const Mesh *pMesh,
	const StucMapArr *pMapArr,
	const Attrib *pAttrib,
	StucDomain domain,
	const AttribIndexed ***pppOut
) {
	StucErr err = PIX_ERR_SUCCESS;
	const AttribIndexed **ppAttribs = pCtx->alloc.fpCalloc(pMapArr->count, sizeof(void *));
	bool found = false;
	bool same = true;
	StucMap pMapCache = NULL;
	for (I32 i = 0; i < pMapArr->count; ++i) {
		const StucMap pMap = pMapArr->pArr[i].map.ptr;
		const char *pName = NULL;
		switch (pAttrib->origin) {
			case STUC_ATTRIB_ORIGIN_MAP:
				pName = pAttrib->core.name;
				break;
			case STUC_ATTRIB_ORIGIN_COMMON: {
				const AttribArray *pMapAttribArr =
					stucGetAttribArrFromDomainConst(&pMap->pMesh->core, domain);
				const Attrib *pMapAttrib = NULL;
				err = stucGetMatchingAttribConst(
					pCtx,
					&pMap->pMesh->core, pMapAttribArr,
					&pMesh->core, pAttrib,
					true,
					false,
					&pMapAttrib
				);
				PIX_ERR_THROW_IFNOT_COND(err, pMapAttrib, "", 0);
				pName = pMapAttrib->core.name;
				break;
			}
			default:
				PIX_ERR_ASSERT("invalid attrib origin", false);
		}
		const AttribIndexed *pIndexedAttrib =
			stucGetAttribIndexedInternConst(&pMap->indexedAttribs, pName);
		if (pIndexedAttrib) {
			found = true;
			ppAttribs[i] = pIndexedAttrib;
			if (!pMapCache) {
				pMapCache = pMapArr->pArr[i].map.ptr;
			}
			else if (same) {
				same = pMapCache == pMapArr->pArr[i].map.ptr;
			}
		}
	}
	if (found) {
		*pppOut = ppAttribs;
		return err;
	}
	PIX_ERR_CATCH(0, err, ;);
	pCtx->alloc.fpFree(ppAttribs);
	*pppOut = NULL;
	return err;
}

static
StucErr correctIdxIndices(
	StucContext pCtx,
	const char *pName,
	Mesh *pMeshArr,
	const StucMapArr *pMapArr,
	const AttribIndexedArr *pInIndexedAttribs,
	const AttribIndexed **ppMapAttribs,
	AttribIndexed *pOutIndexedAttrib,
	StucDomain domain
) {
	StucErr err = PIX_ERR_SUCCESS;
	const AttribIndexed *pInIndexedAttrib =
		stucGetAttribIndexedInternConst(pInIndexedAttribs, pName);
	for (I32 i = 0; i < pMapArr->count; ++i) {
		AttribArray *pAttribArr = stucGetAttribArrFromDomain(&pMeshArr[i].core, domain);
		I32 attribIdx = 0;
		Attrib *pAttrib =
			stucGetAttribIntern(pName, pAttribArr, false, NULL, NULL, &attribIdx);
		if (!ppMapAttribs[i] || !pAttrib) {
			continue;
		}
		Mesh *pMesh = pMeshArr + i;
		const AttribIndexed *pOriginIndexedAttrib = NULL;
		err = getOriginIndexedAttrib(
			pCtx,
			pAttrib,
			attribIdx,
			pMapArr,
			i,
			ppMapAttribs[i],
			pInIndexedAttrib,
			&pOriginIndexedAttrib,
			domain
		);
		PIX_ERR_RETURN_IFNOT(err, "");
		PIX_ERR_RETURN_IFNOT_COND(
			err,
			pOriginIndexedAttrib,
			"no indexed attrib found for idx attrib in mesh"
		);
		err = iterFacesAndCorrectIdxAttrib(
			pCtx,
			pAttrib,
			pMesh,
			pOutIndexedAttrib,
			pOriginIndexedAttrib,
			domain
		);
		PIX_ERR_RETURN_IFNOT(err, "");
	}
	return err;
}

static
StucErr mergeIndexedAttribs(
	StucContext pCtx,
	Mesh *pMeshArr,
	const StucMapArr *pMapArr,
	const AttribIndexedArr *pInIndexedAttribs,
	AttribIndexedArr *pOutIndexedAttribs,
	bool keepExisting
) {
	StucErr err = PIX_ERR_SUCCESS;
	const StucAlloc *pAlloc = &pCtx->alloc;
	pOutIndexedAttribs->size = pInIndexedAttribs->count;
	pOutIndexedAttribs->pArr =
		pAlloc->fpCalloc(pOutIndexedAttribs->size, sizeof(AttribIndexed));
	if (keepExisting) {
		for (I32 i = 0; i < pInIndexedAttribs->count; ++i) {
			err = stucAppendAndCopyIdxAttrib(
				pCtx,
				pInIndexedAttribs->pArr + i,
				pOutIndexedAttribs
			);
			PIX_ERR_THROW_IFNOT(err, "", 0);
		}
	}
	for (I32 i = 0; i < pMapArr->count; ++i) {
		Mesh *pMesh = pMeshArr + i;
		for (I32 j = STUC_DOMAIN_FACE; j <= STUC_DOMAIN_VERT; ++j) {
			AttribArray *pAttribArr = stucGetAttribArrFromDomain(&pMesh->core, j);
			for (I32 k = 0; k < pAttribArr->count; ++k) {
				Attrib *pAttrib = pAttribArr->pArr + k;
				if (pAttrib->core.use != STUC_ATTRIB_USE_IDX ||
					pAttrib->origin == STUC_ATTRIB_ORIGIN_MESH_OUT) {
					continue;
				}
				AttribIndexed *pIndexedAttrib = NULL;
				if (pAttrib->origin == STUC_ATTRIB_ORIGIN_MESH_IN) {
					if (!keepExisting) {
						err = stucAppendAndCopyIdxAttribFromName(
							pCtx,
							pAttrib->core.name,
							pInIndexedAttribs,
							pOutIndexedAttribs
						);
						PIX_ERR_THROW_IFNOT(err, "", 0);
					}
					continue;
				}
				const AttribIndexed **ppMapAttribs = NULL;
				stucGetAttribIndexed(pAttrib->core.name, pOutIndexedAttribs, &pIndexedAttrib);
				if (!pIndexedAttrib) {
					err = getIndexedAttribInMaps(
						pCtx,
						pMesh,
						pMapArr,
						pAttrib,
						j,
						&ppMapAttribs
					);
					PIX_ERR_THROW_IFNOT_COND(err, ppMapAttribs, "", 0);
					{
						const AttribIndexed *pRefAttrib = NULL;
						switch (pAttrib->origin) {
							case STUC_ATTRIB_ORIGIN_MAP:
								for (I32 l = 0; l < pMapArr->count; ++l) {
									if (ppMapAttribs[l]) {
										pRefAttrib = ppMapAttribs[l];
										break;
									}
								}
								PIX_ERR_ASSERT("", pRefAttrib);
								pIndexedAttrib = stucAppendIndexedAttrib(
									pCtx,
									pOutIndexedAttribs,
									pRefAttrib->core.name,
									0, //dont allocate pData
									pRefAttrib->core.type,
									pRefAttrib->core.use
								);
								break;
							case STUC_ATTRIB_ORIGIN_COMMON:
								pRefAttrib = stucGetAttribIndexedInternConst(
									pInIndexedAttribs,
									pAttrib->core.name
								);
								if (!keepExisting) {
									err = stucAppendAndCopyIdxAttrib(
										pCtx,
										pRefAttrib,
										pOutIndexedAttribs
									);
									PIX_ERR_THROW_IFNOT(err, "", 0);
								}
								break;
							default:
								PIX_ERR_ASSERT("invalid attrib origin", false);
						}
					}
				}
				err = correctIdxIndices(
					pCtx,
					pAttrib->core.name,
					pMeshArr,
					pMapArr,
					pInIndexedAttribs,
					ppMapAttribs,
					pIndexedAttrib,
					j
				);
				PIX_ERR_THROW_IFNOT(err, "", 0);
			}
		}
	}
	PIX_ERR_CATCH(0, err, ;)
	return err;
}

typedef struct StucMapToMeshArgs {
	StucContext pCtx;
	StucMapArr *pMapArr;
	StucMesh *pMeshIn;
	StucAttribIndexedArr *pInIndexedAttribs;
	StucMesh *pMeshOut;
	StucAttribIndexedArr *pOutIndexedAttribs;
	F32 wScale;
	F32 receiveLen;
} StucMapToMeshArgs;

static
StucErr mapToMeshFromJob(void *pArgsVoid) {
	StucMapToMeshArgs *pArgs = pArgsVoid;
	return stucMapToMesh(
		pArgs->pCtx,
		pArgs->pMapArr,
		pArgs->pMeshIn,
		pArgs->pInIndexedAttribs,
		pArgs->pMeshOut,
		pArgs->pOutIndexedAttribs,
		pArgs->wScale,
		pArgs->receiveLen,
		false
	);
}

StucErr stucQueueMapToMesh(
	StucContext pCtx,
	void **ppJobHandle,
	StucMapArr *pMapArr,
	StucMesh *pMeshIn,
	StucAttribIndexedArr *pInIndexedAttribs,
	StucMesh *pMeshOut,
	StucAttribIndexedArr *pOutIndexedAttribs,
	F32 wScale,
	F32 receiveLen
) {
	StucMapToMeshArgs *pArgs = pCtx->alloc.fpCalloc(1, sizeof(StucMapToMeshArgs));
	pArgs->pCtx = pCtx;
	pArgs->pMapArr = pMapArr;
	pArgs->pMeshIn = pMeshIn;
	pArgs->pInIndexedAttribs = pInIndexedAttribs;
	pArgs->pMeshOut = pMeshOut;
	pArgs->pOutIndexedAttribs = pOutIndexedAttribs;
	pArgs->wScale = wScale;
	pArgs->receiveLen = receiveLen;
	pCtx->threadPool.pJobStackPushJobs(
		pCtx->pThreadPoolHandle,
		1,
		ppJobHandle,
		mapToMeshFromJob,
		(void **)&pArgs
	);
	return PIX_ERR_SUCCESS;
}

static
StucErr mapMapArrToMesh(
	StucContext pCtx,
	const StucMapArr *pMapArr,
	Mesh *pMeshIn,
	const StucAttribIndexedArr *pInIndexedAttribs,
	StucMesh *pMeshOut,
	StucAttribIndexedArr *pOutIndexedAttribs,
	F32 wScale,
	F32 receiveLen,
	bool keepExistingIdxAttribs
) {
	StucErr err = PIX_ERR_SUCCESS;
	Mesh *pOutBufArr = pCtx->alloc.fpCalloc(pMapArr->count, sizeof(Mesh));
	StucObjArr outObjWrapArr = {0};
	outObjWrapArr.size = outObjWrapArr.count = pMapArr->count;
	outObjWrapArr.pArr = pCtx->alloc.fpCalloc(outObjWrapArr.size, sizeof(StucObject));
	for (I32 i = 0; i < pMapArr->count; ++i) {
		outObjWrapArr.pArr[i].pData = (StucObjectData *)&pOutBufArr[i];
		const StucMap pMap = pMapArr->pArr[i].map.ptr;
		I8 matIdx = pMapArr->pArr[i].matIdx;
		InFaceTable inFaceTable = {0};
		if (pMap->usgArr.count) {
			//set preserve to null to prevent usg squares from being split
			if (pMeshIn->pEdgePreserve || pMeshIn->pVertPreserve) {
				pMeshIn->pEdgePreserve = NULL;
				pMeshIn->pVertPreserve = NULL;
			}
			MapFile squares = { .pMesh = pMap->usgArr.pSquares };
			buildFaceBBoxes(&pCtx->alloc, &squares);
			err = stucCreateQuadTree(
				pCtx,
				&squares.quadTree,
				squares.pMesh,
				squares.pFaceBBoxes
			);
			PIX_ERR_THROW_IFNOT(err, "failed to create usg quadtree", 0);
			StucMesh squaresOut = { 0 };
			err = mapToMeshInternal(
				pCtx,
				&squares,
				pMeshIn,
				&squaresOut,
				matIdx,
				pMapArr->pArr[i].blendOptArr,
				&inFaceTable,
				1.0f,
				-1.0f
			);
			PIX_ERR_THROW_IFNOT(err, "map to mesh usg failed", 1);
			err = stucSampleInAttribsAtUsgOrigins(
				pCtx,
				pMap,
				pMeshIn,
				&squaresOut,
				inFaceTable.pArr
			);
			PIX_ERR_THROW_IFNOT(err, "", 1);
			InFaceTableToHashTable(&pCtx->alloc, pMap, squaresOut.faceCount, inFaceTable.pArr);
			//*pMeshOut = squaresOut;
			//return PIX_ERR_SUCCESS;
			stucMeshDestroy(pCtx, &squaresOut);
			stucAssignActiveAliases(
				pCtx,
				(Mesh *)pMeshIn,
				STUC_ATTRIB_USE_FIELD(((StucAttribUse[]) { //reassign preserve if present
					STUC_ATTRIB_USE_PRESERVE_EDGE,
					STUC_ATTRIB_USE_PRESERVE_VERT
				})),
				STUC_DOMAIN_NONE
			);
		}
		err = mapToMeshInternal(
			pCtx,
			pMap,
			pMeshIn,
			&pOutBufArr[i].core,
			matIdx,
			pMapArr->pArr[i].blendOptArr,
			NULL,
			wScale,
			receiveLen
		);
		PIX_ERR_THROW_IFNOT(err, "map to mesh failed", 1);
		PIX_ERR_CATCH(1, err, ;);
		if (pMap->usgArr.count) {
			pCtx->alloc.fpFree(pMap->usgArr.pInFaceTable);
			pMap->usgArr.pInFaceTable = NULL;
			pixalcLinAllocDestroy(&inFaceTable.alloc);
			inFaceTable = (InFaceTable) {0};
		}
		PIX_ERR_THROW_IFNOT(err, "", 0);
	}
	pMeshOut->type.type = STUC_OBJECT_DATA_MESH;
	Mesh meshOutWrap = {.core = *pMeshOut};
	err = mergeIndexedAttribs(
		pCtx,
		pOutBufArr,
		pMapArr,
		pInIndexedAttribs,
		pOutIndexedAttribs,
		keepExistingIdxAttribs
	);
	PIX_ERR_THROW_IFNOT(err, "", 0);
	err = stucMergeObjArr(pCtx, &meshOutWrap, &outObjWrapArr, false);
	PIX_ERR_THROW_IFNOT(err, "", 0);
	*pMeshOut = meshOutWrap.core;
	PIX_ERR_CATCH(0, err,
		stucMeshDestroy(pCtx, pMeshOut);
	);
	//meshes are stored on an arr buf, which we can't call stucObjArrDestroy
	for (I32 i = 0; i < pMapArr->count; ++i) {
		stucMeshDestroy(pCtx, &pOutBufArr[i].core);
	}
	pCtx->alloc.fpFree(pOutBufArr);
	stucObjArrDestroy(pCtx, &outObjWrapArr);
	return err;
}

static
StucErr appendSpAttribsToInMesh(
	const StucContext pCtx,
	Mesh *pWrap,
	const StucMesh *pMeshIn,
	UBitField32 flags
) {
	StucErr err = PIX_ERR_SUCCESS;
	UBitField32 has = 0;
	stucQuerySpAttribs(pCtx, pMeshIn, flags, &has);
	if (has) {
		PIX_ERR_RETURN(err, "in-mesh contains attribs it shouldn't");
	}
	Mesh meshInCpy = {.core = *pMeshIn};
	const Mesh *const pMeshInCpyPtr = &meshInCpy;
	stucAllocAttribsFromMeshArr(
		pCtx,
		pWrap,
		1, &pMeshInCpyPtr,
		0,
		false,
		false, //dont allocate data
		true //alias pMeshIn's data instead
	);
	stucAppendSpAttribsToMesh(
		pCtx,
		pWrap,
		flags, 
		STUC_ATTRIB_ORIGIN_MESH_IN
	);
	return err;
}

static
void destroyAppendedSpAttribs(StucContext pCtx, StucMesh *pMesh, UBitField32 flags) {
	for (I32 i = 1; i < STUC_ATTRIB_USE_SP_ENUM_COUNT; ++i) {
		if (!(flags >> i & 0x1)) {
			continue;
		}
		Attrib *pAttrib = stucGetActiveAttrib(pCtx, pMesh, i);
		if (pAttrib) {
			if (pAttrib->core.pData) {
				pCtx->alloc.fpFree(pAttrib->core.pData);
				pAttrib->core.pData = NULL;
			}
		}
	}
	if (pMesh->faceAttribs.pArr) {
		pCtx->alloc.fpFree(pMesh->faceAttribs.pArr);
		pMesh->faceAttribs.pArr = NULL;
	}
	if (pMesh->cornerAttribs.pArr) {
		pCtx->alloc.fpFree(pMesh->cornerAttribs.pArr);
		pMesh->cornerAttribs.pArr = NULL;
	}
	if (pMesh->edgeAttribs.pArr) {
		pCtx->alloc.fpFree(pMesh->edgeAttribs.pArr);
		pMesh->edgeAttribs.pArr = NULL;
	}
	if (pMesh->vertAttribs.pArr) {
		pCtx->alloc.fpFree(pMesh->vertAttribs.pArr);
		pMesh->vertAttribs.pArr = NULL;
	}
}

static
StucErr initMeshInWrap(
	StucContext pCtx,
	Mesh *pWrap,
	StucMesh meshIn, //passed by value so we can set active attrib domains if missing
	UBitField32 spAttribsToAppend,
	bool *pBuildEdges
) {
	StucErr err = PIX_ERR_SUCCESS;
	err = attemptToSetMissingActiveDomains(&meshIn);
	PIX_ERR_RETURN_IFNOT(err, "");
	stucAliasMeshCoreNoAttribs(&pWrap->core, &meshIn);
	*pBuildEdges = !meshIn.edgeCount;
	if (*pBuildEdges) {
		printf("no edge list found, building one\n");
		PIX_ERR_RETURN_IFNOT_COND(
			err,
			!meshIn.edgeAttribs.count,
			"in-mesh has edge attribs, yet no edge list"
		);
		err = stucBuildEdgeList(pCtx, pWrap);
		PIX_ERR_RETURN_IFNOT(err, "failed to build edge list");
		printf("finished building edge list\n");
	}
	err = appendSpAttribsToInMesh(pCtx, pWrap, &meshIn, spAttribsToAppend);
	PIX_ERR_RETURN_IFNOT(err, "");
	stucSetAttribOrigins(&pWrap->core.meshAttribs, STUC_ATTRIB_ORIGIN_MESH_IN);
	stucSetAttribOrigins(&pWrap->core.faceAttribs, STUC_ATTRIB_ORIGIN_MESH_IN);
	stucSetAttribOrigins(&pWrap->core.cornerAttribs, STUC_ATTRIB_ORIGIN_MESH_IN);
	stucSetAttribOrigins(&pWrap->core.edgeAttribs, STUC_ATTRIB_ORIGIN_MESH_IN);
	stucSetAttribOrigins(&pWrap->core.vertAttribs, STUC_ATTRIB_ORIGIN_MESH_IN);

	err = stucAssignActiveAliases(
		pCtx,
		pWrap,
		~STUC_ATTRIB_USE_FIELD(((StucAttribUse[]) { //all except for
			STUC_ATTRIB_USE_RECEIVE,
			STUC_ATTRIB_USE_USG,
			STUC_ATTRIB_USE_EDGE_LEN
		})),
		STUC_DOMAIN_NONE
	);
	PIX_ERR_RETURN_IFNOT(err, "");

	//err = stucBuildTangents(pWrap);
	PIX_ERR_RETURN_IFNOT(err, "failed to build tangents");
	buildEdgeAdj(pWrap);
	buildSeamAndPreserveTables(pWrap);

	//set sp
	stucSetAttribCopyOpt(
		pCtx,
		&pWrap->core, 
		STUC_ATTRIB_DONT_COPY,
		spAttribsToAppend | STUC_ATTRIB_USE_FIELD(((StucAttribUse[]) {
			STUC_ATTRIB_USE_WSCALE,
			STUC_ATTRIB_USE_NORMALS_VERT
		}))
	);
	//set required
	stucSetAttribCopyOpt(
		pCtx,
		&pWrap->core,
		STUC_ATTRIB_COPY,
		STUC_ATTRIB_USE_FIELD(((StucAttribUse[]) {
			STUC_ATTRIB_USE_POS,
			STUC_ATTRIB_USE_UV,
			STUC_ATTRIB_USE_NORMAL,
			STUC_ATTRIB_USE_IDX
		}))
	);

	return err;
}

StucErr stucMapToMesh(
	StucContext pCtx,
	const StucMapArr *pMapArr,
	const StucMesh *pMeshIn,
	const StucAttribIndexedArr *pInIndexedAttribs,
	StucMesh *pMeshOut,
	StucAttribIndexedArr *pOutIndexedAttribs,
	F32 wScale,
	F32 receiveLen,
	bool keepExistingIdxAttribs
) {
	StucErr err = PIX_ERR_SUCCESS;
	PIX_ERR_RETURN_IFNOT_COND(err, pMeshIn, "");
	err = stucValidateMesh(pMeshIn, false, false);
	PIX_ERR_RETURN_IFNOT(err, "invalid in-mesh");
	Mesh meshInWrap = {0};
	UBitField32 spAttribsToAppend = STUC_ATTRIB_USE_FIELD(((StucAttribUse[]) {
		STUC_ATTRIB_USE_TANGENT,
		STUC_ATTRIB_USE_TSIGN,
		STUC_ATTRIB_USE_SEAM_EDGE,
		STUC_ATTRIB_USE_SEAM_VERT,
		STUC_ATTRIB_USE_NUM_ADJ_PRESERVE,
		STUC_ATTRIB_USE_EDGE_FACES,
		STUC_ATTRIB_USE_EDGE_CORNERS
	}));
	bool builtEdges = false;
	err = initMeshInWrap(
		pCtx,
		&meshInWrap,
		*(StucMesh *)pMeshIn,
		spAttribsToAppend,
		&builtEdges
	);
	PIX_ERR_THROW_IFNOT(err, "", 0);

	PIX_ERR_THROW_IFNOT_COND(
		err,
		pMapArr && pMapArr->count && pMapArr->pArr,
		"", 0
	);

	err = mapMapArrToMesh(
		pCtx,
		pMapArr,
		&meshInWrap,
		pInIndexedAttribs,
		pMeshOut,
		pOutIndexedAttribs,
		wScale,
		receiveLen,
		keepExistingIdxAttribs
	);
	PIX_ERR_THROW_IFNOT(err, "mapMapArrToMesh returned error", 0);
	printf("----------------------FINISHING IN-MESH\n");
	PIX_ERR_CATCH(0, err, ;);
	if (builtEdges && meshInWrap.core.pEdges) {
		if (meshInWrap.core.pEdges) {
			pCtx->alloc.fpFree(meshInWrap.core.pEdges);
			meshInWrap.core.pEdges = NULL;
		}
	}
	destroyAppendedSpAttribs(pCtx, &meshInWrap.core, spAttribsToAppend);
	return err;
}

StucErr stucUsgArrDestroy(StucContext pCtx, I32 count, StucUsg *pUsgArr) {
	StucErr err = PIX_ERR_NOT_SET;
	for (I32 i = 0; i < count; ++i) {
		err = stucMeshDestroy(pCtx, (StucMesh *)pUsgArr[i].obj.pData);
		PIX_ERR_THROW_IFNOT(err, "", 0);
	}
	pCtx->alloc.fpFree(pUsgArr);
	PIX_ERR_CATCH(0, err, ;)
	return err;
}

StucErr stucMeshDestroy(StucContext pCtx, StucMesh *pMesh) {
	for (I32 i = 0; i < pMesh->meshAttribs.count; ++i) {
		if (pMesh->meshAttribs.pArr[i].core.pData) {
			pCtx->alloc.fpFree(pMesh->meshAttribs.pArr[i].core.pData);
		}
	}
	for (I32 i = 0; i < pMesh->faceAttribs.count; ++i) {
		if (pMesh->faceAttribs.pArr[i].core.pData) {
			pCtx->alloc.fpFree(pMesh->faceAttribs.pArr[i].core.pData);
		}
	}
	for (I32 i = 0; i < pMesh->cornerAttribs.count; ++i) {
		if (pMesh->cornerAttribs.pArr[i].core.pData) {
			pCtx->alloc.fpFree(pMesh->cornerAttribs.pArr[i].core.pData);
		}
	}
	for (I32 i = 0; i < pMesh->edgeAttribs.count; ++i) {
		if (pMesh->edgeAttribs.pArr[i].core.pData) {
			pCtx->alloc.fpFree(pMesh->edgeAttribs.pArr[i].core.pData);
		}
	}
	for (I32 i = 0; i < pMesh->vertAttribs.count; ++i) {
		if (pMesh->vertAttribs.pArr[i].core.pData) {
			pCtx->alloc.fpFree(pMesh->vertAttribs.pArr[i].core.pData);
		}
	}
	if (pMesh->meshAttribs.pArr) {
		pCtx->alloc.fpFree(pMesh->meshAttribs.pArr);
	}
	if (pMesh->faceAttribs.pArr) {
		pCtx->alloc.fpFree(pMesh->faceAttribs.pArr);
	}
	if (pMesh->edgeAttribs.pArr) {
		pCtx->alloc.fpFree(pMesh->edgeAttribs.pArr);
	}
	if (pMesh->cornerAttribs.pArr) {
		pCtx->alloc.fpFree(pMesh->cornerAttribs.pArr);
	}
	if (pMesh->vertAttribs.pArr) {
		pCtx->alloc.fpFree(pMesh->vertAttribs.pArr);
	}
	if(pMesh->pFaces) {
		pCtx->alloc.fpFree(pMesh->pFaces);
	}
	if (pMesh->pCorners) {
		pCtx->alloc.fpFree(pMesh->pCorners);
	}
	if (pMesh->pEdges) {
		pCtx->alloc.fpFree(pMesh->pEdges);
	}
	return PIX_ERR_SUCCESS;
}

StucErr stucGetAttribSize(StucAttribCore *pAttrib, I32 *pSize) {
	*pSize = stucGetAttribSizeIntern(pAttrib->type);
	return PIX_ERR_SUCCESS;
}

StucErr stucGetAttrib(const char *pName, StucAttribArray *pAttribs, StucAttrib **ppAttrib) {
	*ppAttrib = stucGetAttribIntern(pName, pAttribs, false, NULL, NULL, NULL);
	return PIX_ERR_SUCCESS;
}

StucErr stucAttribGetAsVoid(StucAttribCore *pAttrib, int32_t idx, void **ppOut) {
	StucErr err = PIX_ERR_SUCCESS;
	PIX_ERR_RETURN_IFNOT_COND(err, pAttrib && ppOut && idx >= 0, "");
	*ppOut = stucAttribAsVoid(pAttrib, idx);
	return PIX_ERR_SUCCESS;
}

StucErr stucAttribActiveGet(
	StucContext pCtx,
	StucMesh *pMesh,
	StucAttribUse use,
	Attrib **ppAttrib
) {
	*ppAttrib = stucGetActiveAttrib(pCtx, pMesh, use);
	return PIX_ERR_SUCCESS;
}

StucErr stucGetAttribIndexed(
	const char *pName,
	StucAttribIndexedArr *pAttribs,
	StucAttribIndexed **ppAttrib
) {
	*ppAttrib = stucGetAttribIndexedIntern(pAttribs, pName);
	return PIX_ERR_SUCCESS;
}

void stucMapIndexedAttribsGet(
	StucContext pCtx,
	StucMap pMap,
	StucAttribIndexedArr *pIndexedAttribs
) {
	*pIndexedAttribs = pMap->indexedAttribs;
}

StucErr stucWaitForJobs(
	StucContext pCtx,
	I32 count,
	void **ppHandles,
	bool wait,
	bool *pDone
) {
	return pCtx->threadPool.fpWaitForJobs(
		pCtx->pThreadPoolHandle,
		count,
		ppHandles,
		wait,
		pDone
	);
}

StucErr stucJobGetErrs(
	StucContext pCtx,
	I32 jobCount,
	void ***pppJobHandles
) {
	StucErr err = PIX_ERR_SUCCESS;
	PIX_ERR_ASSERT("", pCtx && pppJobHandles);
	PIX_ERR_ASSERT("", jobCount > 0);
	for (I32 i = 0; i < jobCount; ++i) {
		StucErr jobErr = PIX_ERR_NOT_SET;
		err = pCtx->threadPool.fpGetJobErr(
			pCtx->pThreadPoolHandle,
			(*pppJobHandles)[i],
			&jobErr
		);
		PIX_ERR_THROW_IFNOT_COND(err, jobErr == PIX_ERR_SUCCESS, "", 0);
	}
	PIX_ERR_CATCH(0, err, ;);
	return err;
}

void stucJobDestroyHandles(
	StucContext pCtx,
	I32 jobCount,
	void **ppJobHandles
) {
	PIX_ERR_ASSERT("", pCtx && ppJobHandles);
	PIX_ERR_ASSERT("", jobCount > 0);
	for (I32 i = 0; i < jobCount; ++i) {
		pCtx->threadPool.fpJobHandleDestroy(
			pCtx->pThreadPoolHandle,
			ppJobHandles + i
		);
	}
}

StucErr stucAttribSpTypesGet(StucContext pCtx, const AttribType **ppTypes) {
	StucErr err = PIX_ERR_SUCCESS;
	PIX_ERR_RETURN_IFNOT_COND(err, pCtx && ppTypes, "");
	*ppTypes = pCtx->spAttribTypes;
	return err;
}

StucErr stucAttribSpDomainsGet(StucContext pCtx, const StucDomain **ppDomains) {
	StucErr err = PIX_ERR_SUCCESS;
	PIX_ERR_RETURN_IFNOT_COND(err, pCtx && ppDomains, "");
	*ppDomains = pCtx->spAttribDomains;
	return err;
}

StucErr stucAttribSpIsValid(
	StucContext pCtx,
	const AttribCore *pCore,
	StucDomain domain
) {
	StucErr err = PIX_ERR_SUCCESS;
	PIX_ERR_RETURN_IFNOT_COND(err, pCtx && pCore, "");
	return 
		pCtx->spAttribTypes[pCore->use] == pCore->type &&
		pCtx->spAttribDomains[pCore->use] == domain;
}

StucErr stucAttribGetAllDomains(
	StucContext pCtx,
	StucMesh *pMesh,
	const char *pName,
	Attrib **ppAttrib,
	StucDomain *pDomain
) {
	for (I32 i = STUC_DOMAIN_FACE; i <= STUC_DOMAIN_VERT; ++i) {
		AttribArray *pArr = stucGetAttribArrFromDomain(pMesh, i);
		Attrib *pAttrib = stucGetAttribIntern(pName, pArr, false, NULL, NULL, NULL);
		if (pAttrib) {
			*ppAttrib = pAttrib;
			*pDomain = i;
			break;
		}
	}
	return PIX_ERR_SUCCESS;
}

StucErr stucAttribGetAllDomainsConst(
	StucContext pCtx,
	const StucMesh *pMesh,
	const char *pName,
	const StucAttrib **ppAttrib,
	StucDomain *pDomain
) {
	return stucAttribGetAllDomains(pCtx, (StucMesh *)pMesh, pName, (StucAttrib **)ppAttrib, pDomain);
}

StucErr stucAttribArrGet(
	StucContext pCtx,
	StucMesh *pMesh,
	StucDomain domain,
	StucAttribArray **ppArr
) {
	StucErr err = PIX_ERR_SUCCESS;
	PIX_ERR_RETURN_IFNOT_COND(err, pCtx && pMesh && ppArr, "");
	*ppArr = stucGetAttribArrFromDomain(pMesh, domain);
	return err;
}

StucErr stucAttribArrGetConst(
	StucContext pCtx,
	const StucMesh *pMesh,
	StucDomain domain,
	const StucAttribArray **ppArr
) {
	return stucAttribArrGet(
		pCtx,
		(StucMesh *)pMesh,
		domain,
		(AttribArray **)ppArr
	);
}

StucErr stucAttribGetCompType(
	StucContext pCtx,
	StucAttribType type,
	StucAttribType *pCompType
) {
	StucErr err = PIX_ERR_SUCCESS;
	PIX_ERR_RETURN_IFNOT_COND(err, pCtx && pCompType, "");
	*pCompType = stucAttribGetCompTypeIntern(type);
	return err;
}

StucErr stucAttribTypeGetVecSize(
	StucContext pCtx,
	StucAttribType type,
	I32 *pSize
) {
	StucErr err = PIX_ERR_SUCCESS;
	PIX_ERR_RETURN_IFNOT_COND(err, pCtx && pSize, "");
	*pSize = stucAttribTypeGetVecSizeIntern(type);
	return err;
}

StucErr stucDomainCountGet(
	StucContext pCtx,
	const StucMesh *pMesh,
	StucDomain domain,
	int32_t *pCount
) {
	StucErr err = PIX_ERR_SUCCESS;
	PIX_ERR_RETURN_IFNOT_COND(err, pCtx && pMesh && pCount, "");
	*pCount = stucDomainCountGetIntern(pMesh, domain);
	return err;
}

StucErr stucAttribIndexedArrDestroy(StucContext pCtx, StucAttribIndexedArr *pArr) {
	StucErr err = PIX_ERR_SUCCESS;
	PIX_ERR_RETURN_IFNOT_COND(err, pCtx && pArr, "");
	if (pArr->pArr) {
		for (I32 i = 0; i < pArr->count; ++i) {
			if (pArr->pArr[i].core.pData) {
				pCtx->alloc.fpFree(pArr->pArr[i].core.pData);
			}
		}
		pCtx->alloc.fpFree(pArr->pArr);
		pArr->pArr = NULL;
	}
	pArr->count = pArr->size = 0;
	return err;
}

StucErr stucMapArrDestroy(StucContext pCtx, StucMapArr *pMapArr) {
	StucErr err = PIX_ERR_SUCCESS;
	PIX_ERR_RETURN_IFNOT_COND(err, pCtx && pMapArr, "");
	if (pMapArr->pArr) {
		for (I32 i = 0; i < pMapArr->count; ++i) {
			stucDestroyBlendOptArr(pCtx, pMapArr->pArr[i].blendOptArr);
		}
		pCtx->alloc.fpFree(pMapArr->pArr);
	}
	*pMapArr = (StucMapArr){0};
	return err;
}

StucErr stucObjectInit(
	StucContext pCtx,
	StucObject *pObj,
	StucMesh *pMesh,
	const Stuc_M4x4 *pTransform
) {
	StucErr err = PIX_ERR_SUCCESS;
	PIX_ERR_RETURN_IFNOT_COND(err, pCtx && pObj, "");
	pObj->pData = (StucObjectData *)pMesh;
	if (pTransform) {
		pObj->transform = *pTransform;
	}
	else {
		pObj->transform = PIX_MATH_IDENT_MAT4X4;
	}
	return err;
}
