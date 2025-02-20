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
void setDefaultStageReport(StucContext pCtx) {
	pCtx->stageReport.outOf = 50,
	pCtx->stageReport.pBegin = stucStageBegin;
	pCtx->stageReport.pProgress = stucStageProgress;
	pCtx->stageReport.pEnd = stucStageEnd;
}

StucResult stucContextInit(
	StucContext *pCtx,
	StucAlloc *pAlloc,
	StucThreadPool *pThreadPool,
	StucIo *pIo,
	StucTypeDefaultConfig *pTypeDefaultConfig,
	StucStageReport *pStageReport
) {
	StucAlloc alloc;
	if (pAlloc) {
		stucAllocSetCustom(&alloc, pAlloc);
	}
	else {
		stucAllocSetDefault(&alloc);
	}
	*pCtx = alloc.pCalloc(1, sizeof(StucContextInternal));
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
	(*pCtx)->threadPool.pInit(
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
	return STUC_SUCCESS;
}

StucResult stucContextDestroy(StucContext pCtx) {
	pCtx->threadPool.pDestroy(pCtx->pThreadPoolHandle);
	pCtx->alloc.pFree(pCtx);
	return STUC_SUCCESS;
}

StucResult stucMapFileExport(
	StucContext pCtx,
	char *pName,
	I32 objCount,
	StucObject *pObjArr,
	I32 usgCount,
	StucUsg *pUsgArr,
	StucAttribIndexedArr *pIndexedAttribs
) {
	return stucWriteStucFile(
		pCtx,
		pName,
		objCount,
		pObjArr,
		usgCount,
		pUsgArr,
		pIndexedAttribs
	);
}

//TODO replace these with StucUsg and StucObj arr structs, that combine arr and count
StucResult stucMapFileLoadForEdit(
	StucContext pCtx,
	char *filePath,
	I32 *pObjCount,
	StucObject **ppObjArr,
	I32 *pUsgCount,
	StucUsg **ppUsgArr,
	I32 *pFlatCutoffCount,
	StucObject **ppFlatCutoffArr,
	StucAttribIndexedArr *pIndexedAttribs
) {
	return stucLoadStucFile(
		pCtx,
		filePath,
		pObjCount,
		ppObjArr,
		pUsgCount,
		ppUsgArr,
		pFlatCutoffCount,
		ppFlatCutoffArr, 
		true,
		pIndexedAttribs
	);
}

static
void buildEdgeLenList(StucContext pCtx, Mesh *pMesh) {
	STUC_ASSERT("", pMesh->pEdgeLen);
	V3_F32 *pPosCache = pCtx->alloc.pMalloc(pMesh->core.edgeCount * sizeof(V3_F32));
	I8 *pSet = pCtx->alloc.pCalloc(pMesh->core.edgeCount, 1);
	for (I32 i = 0; i < pMesh->core.cornerCount; ++i) {
		V3_F32 pos = pMesh->pVerts[pMesh->core.pCorners[i]];
		I32 edge = pMesh->core.pEdges[i];
		if (!pSet[edge]) {
			pPosCache[edge] = pos;
			pSet[edge] = 1;
			continue;
		}
		//STUC_ASSERT("more than 2 corners refernce 1 edge", pSet[edge] < 2);
		V3_F32 diff = _(pos V3SUB pPosCache[edge]);
		pMesh->pEdgeLen[edge] = v3Len(diff);
		pSet[edge]++;
	}
	pCtx->alloc.pFree(pSet);
	pCtx->alloc.pFree(pPosCache);
}

StucResult stucMapFileLoad(StucContext pCtx, StucMap *pMapHandle, char *filePath) {
	StucResult err = STUC_NOT_SET;
	StucMap pMap = pCtx->alloc.pCalloc(1, sizeof(MapFile));
	I32 objCount = 0;
	StucObject *pObjArr = NULL;
	StucUsg *pUsgArr = NULL;
	I32 flatCutoffCount = 0;
	StucObject *pFlatCutoffArr = NULL;
	err = stucLoadStucFile(
		pCtx, filePath,
		&objCount,
		&pObjArr,
		&pMap->usgArr.count,
		&pUsgArr,
		&flatCutoffCount,
		&pFlatCutoffArr,
		false,
		&pMap->indexedAttribs
	);
	//TODO validate meshes, ensure pMatIdx is within mat range, faces are within max corner limit,
	//F32 values are valid, etc.
	STUC_THROW_IFNOT(err, "failed to load file from disk", 0);

	for (I32 i = 0; i < objCount; ++i) {
		err = stucSetSpecialAttribs(pCtx, (Mesh *)pObjArr[i].pData, 0x8ae); //10101110 - all except for preserve
		STUC_THROW_IFNOT(err, "", 0);
		stucApplyObjTransform(pObjArr + i);
	}
	Mesh *pMapMesh = pCtx->alloc.pCalloc(1, sizeof(Mesh));
	pMapMesh->core.type.type = STUC_OBJECT_DATA_MESH_INTERN;
	stucMergeObjArr(pCtx, pMapMesh, objCount, pObjArr, false);

	//append edgeLen attrib
	stucAppendSpAttribsToMesh(pCtx, pMapMesh, 0x1000, STUC_ATTRIB_ORIGIN_MAP);

	stucSetAttribOrigins(&pMapMesh->core.meshAttribs, STUC_ATTRIB_ORIGIN_MAP);
	stucSetAttribOrigins(&pMapMesh->core.faceAttribs, STUC_ATTRIB_ORIGIN_MAP);
	stucSetAttribOrigins(&pMapMesh->core.cornerAttribs, STUC_ATTRIB_ORIGIN_MAP);
	stucSetAttribOrigins(&pMapMesh->core.edgeAttribs, STUC_ATTRIB_ORIGIN_MAP);
	stucSetAttribOrigins(&pMapMesh->core.vertAttribs, STUC_ATTRIB_ORIGIN_MAP);

	stucSetAttribToDontCopy(pCtx, pMapMesh, 0x17f0);
	err = stucSetSpecialAttribs(pCtx, pMapMesh, 0x18ae);
	STUC_THROW_IFNOT(err, "", 0);

	buildEdgeLenList(pCtx, pMapMesh);

	//TODO some form of heap corruption when many objects
	//test with address sanitizer on CircuitPieces.stuc
	stucDestroyObjArr(pCtx, objCount, pObjArr);

	if (pMapMesh->pUvAttrib) {
		//TODO as with all special attributes, allow user to define what should be considered
		//     the primary UV channel. This especially important for integration with other DCCs
		if (!strncmp(
			pMapMesh->pUvAttrib->core.name,
			pCtx->spAttribNames[STUC_ATTRIB_SP_UVS],
			STUC_ATTRIB_NAME_MAX_LEN
		)) {
			char newName[STUC_ATTRIB_NAME_MAX_LEN] = "Map_UVMap";
			memcpy(pMapMesh->pUvAttrib->core.name, newName, STUC_ATTRIB_NAME_MAX_LEN);
		}
	}

	//set corner attribs to interpolate by default
	//TODO make this an option in ui, even for non common attribs
	for (I32 i = 0; i < pMapMesh->core.cornerAttribs.count; ++i) {
		pMapMesh->core.cornerAttribs.pArr[i].interpolate = true;
	}

	pMap->pMesh = pMapMesh;

	//the quadtree is created before USGs are assigned to verts,
	//as the tree's used to speed up the process
	printf("File loaded. Creating quad tree\n");
	err = stucCreateQuadTree(pCtx, &pMap->quadTree, pMap->pMesh);
	STUC_THROW_IFNOT(err, "failed to create quadtree", 0);

	if (pMap->usgArr.count) {
		pMap->usgArr.pArr = pCtx->alloc.pCalloc(pMap->usgArr.count, sizeof(Usg));
		for (I32 i = 0; i < pMap->usgArr.count; ++i) {
			err = stucSetSpecialAttribs(pCtx, (Mesh *)pUsgArr[i].obj.pData, 0x02); //000010 - set only vert pos
			STUC_THROW_IFNOT(err, "", 0);
			pMap->usgArr.pArr[i].origin = *(V2_F32 *)&pUsgArr[i].obj.transform.d[3];
			pMap->usgArr.pArr[i].pMesh = (Mesh *)pUsgArr[i].obj.pData;
			stucApplyObjTransform(&pUsgArr[i].obj);
			if (pUsgArr[i].pFlatCutoff) {
				pMap->usgArr.pArr[i].pFlatCutoff = (Mesh *)pUsgArr[i].pFlatCutoff->pData;
				err = stucSetSpecialAttribs(pCtx, (Mesh *)pUsgArr[i].pFlatCutoff->pData, 0x02); //000010 - set only vert pos
				STUC_THROW_IFNOT(err, "", 0);
				stucApplyObjTransform(pUsgArr[i].pFlatCutoff);
			}
		}
		//TODO remove duplicate uses of alloc where pCtx is present
		//like this
		Mesh *pSquares = pCtx->alloc.pCalloc(1, sizeof(Mesh));
		stucAllocUsgSquaresMesh(pCtx, pMap, pSquares);
		stucFillUsgSquaresMesh(pMap, pUsgArr, pSquares);
		pMap->usgArr.pSquares = pSquares;
		stucAssignUsgsToVerts(&pCtx->alloc, pMap, pUsgArr);
		pMap->usgArr.pMemArr = pUsgArr;
	}

	*pMapHandle = pMap;
	//TODO add proper checks, and return STUC_ERROR if fails.
	//Do for all public functions (or internal ones as well)
	STUC_CATCH(0, err, stucMapFileUnload(pCtx, pMap);)
	return err;
}

StucResult stucMapFileUnload(StucContext pCtx, StucMap pMap) {
	stucDestroyQuadTree(pCtx, &pMap->quadTree);
	if (pMap->pMesh) {
		stucMeshDestroy(pCtx, (StucMesh *)&pMap->pMesh->core);
		pCtx->alloc.pFree((Mesh *)pMap->pMesh);
	}
	if (pMap->usgArr.pSquares) {
		pCtx->alloc.pFree((Mesh *)pMap->usgArr.pSquares);
	}
	pCtx->alloc.pFree(pMap);
	return STUC_SUCCESS;
}

static
void initCommonAttrib(StucContext pCtx, StucCommonAttrib *pEntry, const StucAttrib *pAttrib) {
	memcpy(pEntry->name, pAttrib->core.name, STUC_ATTRIB_NAME_MAX_LEN);
	StucTypeDefault *pDefault = 
		stucGetTypeDefaultConfig(&pCtx->typeDefaults, pAttrib->core.type);
	pEntry->blendConfig = pDefault->blendConfig;
}

static
void getCommonAttribs(
	StucContext pCtx,
	const AttribArray *pMapAttribs,
	const AttribArray *pMeshAttribs,
	I32 *pCommonAttribCount,
	StucCommonAttrib **ppCommonAttribs
) {
	//TODO ignore special attribs like StucTangent or StucTSign
	if (!pMeshAttribs || !pMapAttribs) {
		return;
	}
	I32 count = 0;
	for (I32 i = 0; i < pMeshAttribs->count; ++i) {
		for (I32 j = 0; j < pMapAttribs->count; ++j) {
			if (!strncmp(
				pMeshAttribs->pArr[i].core.name,
				pMapAttribs->pArr[j].core.name,
				STUC_ATTRIB_NAME_MAX_LEN)
				) {
				count++;
			}
		}
	}
	*ppCommonAttribs = count ?
		pCtx->alloc.pMalloc(sizeof(StucCommonAttrib) * count) : NULL;
	count = 0;
	for (I32 i = 0; i < pMeshAttribs->count; ++i) {
		for (I32 j = 0; j < pMapAttribs->count; ++j) {
			if (!strncmp(
				pMeshAttribs->pArr[i].core.name,
				pMapAttribs->pArr[j].core.name,
				STUC_ATTRIB_NAME_MAX_LEN
			)) {
				initCommonAttrib(
					pCtx,
					*ppCommonAttribs + count,
					pMeshAttribs->pArr + i
				);
				count++;
			}
		}
	}
	*pCommonAttribCount = count;
}

//TODO handle edge case, where attribute share the same name,
//but have incompatible types. Such as a F32 and a string.
StucResult stucQueryCommonAttribs(
	StucContext pCtx,
	StucMap pMap,
	StucMesh *pMesh,
	StucCommonAttribList *pCommonAttribs
) {
	getCommonAttribs(
		pCtx,
		&pMap->pMesh->core.meshAttribs,
		&pMesh->meshAttribs,
		&pCommonAttribs->meshCount,
		&pCommonAttribs->pMesh
	);
	getCommonAttribs(
		pCtx,
		&pMap->pMesh->core.faceAttribs,
		&pMesh->faceAttribs,
		&pCommonAttribs->faceCount,
		&pCommonAttribs->pFace
	);
	getCommonAttribs(
		pCtx,
		&pMap->pMesh->core.cornerAttribs,
		&pMesh->cornerAttribs,
		&pCommonAttribs->cornerCount,
		&pCommonAttribs->pCorner
	);
	getCommonAttribs(
		pCtx,
		&pMap->pMesh->core.edgeAttribs,
		&pMesh->edgeAttribs,
		&pCommonAttribs->edgeCount,
		&pCommonAttribs->pEdge
	);
	getCommonAttribs(
		pCtx,
		&pMap->pMesh->core.vertAttribs,
		&pMesh->vertAttribs,
		&pCommonAttribs->vertCount,
		&pCommonAttribs->pVert
	);
	return STUC_SUCCESS;
}

StucResult stucDestroyCommonAttribs(
	StucContext pCtx,
	StucCommonAttribList *pCommonAttribs
) {
	if (pCommonAttribs->pMesh) {
		pCtx->alloc.pFree(pCommonAttribs->pMesh);
	}
	if (pCommonAttribs->pFace) {
		pCtx->alloc.pFree(pCommonAttribs->pFace);
	}
	if (pCommonAttribs->pCorner) {
		pCtx->alloc.pFree(pCommonAttribs->pCorner);
	}
	if (pCommonAttribs->pEdge) {
		pCtx->alloc.pFree(pCommonAttribs->pEdge);
	}
	if (pCommonAttribs->pVert) {
		pCtx->alloc.pFree(pCommonAttribs->pVert);
	}
	return STUC_SUCCESS;
}

static
Result sendOffMappingJobs(
	MapToMeshBasic *pBasic,
	I32 *pJobCount,
	void ***pppJobHandles,
	MappingJobArgs **ppJobArgs
) {
	Result err = STUC_SUCCESS;
	*pJobCount = MAX_SUB_MAPPING_JOBS;
	*pJobCount += *pJobCount == 0;
	I32 facesPerThread = pBasic->pInMesh->core.faceCount / *pJobCount;
	bool singleThread = !facesPerThread;
	*pJobCount = singleThread ? 1 : *pJobCount;
	void *jobArgPtrs[MAX_THREADS] = {0};
	I32 borderTableSize = pBasic->pMap->pMesh->core.faceCount / 5 + 2; //+ 2 incase is 0
	*ppJobArgs = pBasic->pCtx->alloc.pCalloc(*pJobCount, sizeof(MappingJobArgs));
	printf("fromjobsendoff: BorderTableSize: %d\n", borderTableSize);
	for (I32 i = 0; i < *pJobCount; ++i) {
		I32 meshStart = facesPerThread * i;
		I32 meshEnd = i == *pJobCount - 1 ?
			pBasic->pInMesh->core.faceCount : meshStart + facesPerThread;
		(*ppJobArgs)[i].pBasic = pBasic;
		(*ppJobArgs)[i].borderTable.size = borderTableSize;
		(*ppJobArgs)[i].inFaceRange.start = meshStart;
		(*ppJobArgs)[i].inFaceRange.end = meshEnd;
		(*ppJobArgs)[i].pActiveJobs = pJobCount;
		(*ppJobArgs)[i].id = i;
		jobArgPtrs[i] = *ppJobArgs + i;
	}
	*pppJobHandles = pBasic->pCtx->alloc.pCalloc(*pJobCount, sizeof(void *));
	err = pBasic->pCtx->threadPool.pJobStackPushJobs(
		pBasic->pCtx->pThreadPoolHandle,
		*pJobCount,
		*pppJobHandles,
		stucMapToJobMesh,
		jobArgPtrs
	);
	STUC_RETURN_ERR_IFNOT(err, "");
	return err;
}

static
void buildEdgeCornersTable(Mesh *pMesh) {
	const StucMesh *pCore = &pMesh->core;
	memset(pMesh->pEdgeCorners, -1, sizeof(V2_I32) * pCore->edgeCount);
	for (I32 i = 0; i < pCore->cornerCount; ++i) {
		I32 edge = pCore->pEdges[i];
		bool which = (pMesh->pEdgeCorners)[edge].d[0] >= 0;
		(pMesh->pEdgeCorners)[edge].d[which] = i;
	}
}

typedef struct {
	I32 d[2];
	I32 idx;
} EdgeCache;

static
void addVertToTableEntry(
	const Mesh *pMesh,
	FaceRange face,
	I32 localCorner,
	I32 vert,
	I32 edge,
	EdgeCache *pEdgeCache
) {
	//isSeam returns 2 if mesh border, and 1 if uv seam
	I32 isSeam = stucCheckIfEdgeIsSeam(
		edge,
		face,
		localCorner,
		pMesh
	);
	if (isSeam) {
		pMesh->pSeamVert[vert] = (I8)isSeam;
		pMesh->pSeamEdge[edge] = true;
	}
	if (pMesh->pNumAdjPreserve[vert] < 3 &&
		stucCheckIfEdgeIsPreserve(pMesh, edge) &&
		pEdgeCache[vert].d[0] != edge + 1 &&
		pEdgeCache[vert].d[1] != edge + 1) {
		pMesh->pNumAdjPreserve[vert]++;
		I32 *pEdgeCacheIdx = &pEdgeCache[vert].idx;
		pEdgeCache[vert].d[*pEdgeCacheIdx] = edge + 1;
		++*pEdgeCacheIdx;
	}
}

static
void buildSeamAndPreserveTables(StucAlloc *pAlloc, Mesh *pMesh) {
	EdgeCache *pEdgeCache = pAlloc->pCalloc(pMesh->core.vertCount, sizeof(EdgeCache));
	for (I32 i = 0; i < pMesh->core.faceCount; ++i) {
		FaceRange face = {0};
		face.start = pMesh->core.pFaces[i];
		face.end = pMesh->core.pFaces[i + 1];
		face.size = face.end - face.start;
		face.idx = i;
		for (I32 j = 0; j < face.size; ++j) {
			I32 corner = face.start + j;
			I32 vert = pMesh->core.pCorners[corner];
			I32 edge = pMesh->core.pEdges[corner];
			addVertToTableEntry(pMesh, face, j, vert, edge, pEdgeCache);
			I32 prevj = j == 0 ? face.size - 1 : j - 1;
			I32 prevEdge = pMesh->core.pEdges[face.start + prevj];
			addVertToTableEntry(pMesh, face, prevj, vert, prevEdge, pEdgeCache);
		}
	}
	pAlloc->pFree(pEdgeCache);
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
Result mapToFaces( MapToMeshBasic *pBasic, I32 *pJobCount, MappingJobArgs **ppJobArgs, bool *pEmpty) {
	Result err = STUC_SUCCESS;
	void **ppJobHandles = NULL;
	err = sendOffMappingJobs(pBasic, pJobCount, &ppJobHandles, ppJobArgs);
	STUC_THROW_IFNOT(err, "", 0);
	if (!*pJobCount) {
		//no jobs sent
		//implement an STUC_CANCELLED status
		return err;
	}
	err = pBasic->pCtx->threadPool.pWaitForJobs(
		pBasic->pCtx->pThreadPoolHandle,
		*pJobCount,
		ppJobHandles,
		true,
		NULL
	);
	STUC_THROW_IFNOT(err, "", 0);
	*pEmpty = true;
	for (I32 i = 0; i < *pJobCount; ++i) {
		if ((*ppJobArgs)[i].bufSize > 0) {
			*pEmpty = false;
			break;
		}
	}
	err = stucJobGetErrs(pBasic->pCtx, *pJobCount, &ppJobHandles);
	STUC_THROW_IFNOT(err, "", 0);
	STUC_CATCH(0, err, ;);
	if (ppJobHandles) {
		stucJobDestroyHandles(pBasic->pCtx, *pJobCount, ppJobHandles);
		pBasic->pCtx->alloc.pFree(ppJobHandles);
	}
	return err;
}

static
void mappingJobArgsDestroy(MapToMeshBasic *pBasic, I32 jobCount, MappingJobArgs *pJobArgs) {
	for (I32 i = 0; i < jobCount; ++i) {
		BufMesh *pBufMesh = &pJobArgs[i].bufMesh;
		stucMeshDestroy(pBasic->pCtx, &pBufMesh->mesh.core);
		if (pJobArgs[i].borderTable.pTable) {
			pBasic->pCtx->alloc.pFree(pJobArgs[i].borderTable.pTable);
		}
		stucBorderTableDestroyAlloc(&pJobArgs[i].borderTableAlloc);
	}
}

static
Result mapToMeshInternal(
	StucContext pCtx,
	const StucMap pMap,
	const Mesh *pMeshIn,
	StucMesh *pOutMesh,
	I8 maskIdx,
	const StucCommonAttribList *pCommonAttribList,
	InFaceArr **ppInFaceTable,
	F32 wScale,
	F32 receiveLen
) {
	StucResult err = STUC_SUCCESS;
	if (checkIfNoFacesHaveMaskIdx(pMeshIn, maskIdx)) {
		return err;
	}
	MapToMeshBasic basic = {
		.pCtx = pCtx,
		.pMap = pMap,
		.pInMesh = pMeshIn,
		.pCommonAttribList = pCommonAttribList,
		.wScale = wScale,
		.receiveLen = receiveLen,
		.maskIdx = maskIdx,
		.ppInFaceTable = ppInFaceTable
	};
	I32 jobCount = 0;
	MappingJobArgs *pMappingJobArgs = NULL;

	bool empty = true;
	err = mapToFaces(&basic, &jobCount, &pMappingJobArgs, &empty);
	STUC_THROW_IFNOT(err, "", 0);
	if (!empty) {
		err = stucCombineJobMeshes(&basic, pMappingJobArgs, jobCount);
		STUC_THROW_IFNOT(err, "", 0);
		stucReallocMeshToFit(&pCtx->alloc, &basic.outMesh);
		*pOutMesh = basic.outMesh.core;
	}
	STUC_CATCH(0, err, ;);
	mappingJobArgsDestroy(&basic, jobCount, pMappingJobArgs);
	pCtx->alloc.pFree(pMappingJobArgs);
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
	I32 hash = stucFnvHash((U8 *)&sum, 4, pMap->usgArr.tableSize);
	UsgInFace *pEntry = *ppHashTable + hash;
	if (!pEntry->pEntry) {
		pEntry->pEntry = pInFaceTable + squareIdx;
		pEntry->face = pInFaceTable[squareIdx].pArr[inFaceIdx];
		return;
	}
	do {
		if (!pEntry->pNext) {
			pEntry = pEntry->pNext = pAlloc->pCalloc(1, sizeof(UsgInFace));
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
	*ppHashTable = pAlloc->pCalloc(pMap->usgArr.tableSize, sizeof(UsgInFace));
	for (I32 i = 0; i < count; ++i) {
		for (I32 j = 0; j < pInFaceTable[i].count; ++j) {
			addEntryToInFaceTable(pAlloc, ppHashTable, pMap, pInFaceTable, i, j);
		}
	}
}

static
Result getMatsToAdd(
	StucContext pCtx,
	const StucMapArr *pMapArr,
	I32 mapIdx,
	const CommonAttribList *pCommonAttribs,
	const AttribIndexedArr *pInIndexedAttribs,
	const AttribIndexed **ppMatsToAdd
) {
	const AttribIndexedArr *pMapAttribArr = &pMapArr->ppArr[mapIdx]->indexedAttribs;
	const CommonAttribList *pMapCommon = pCommonAttribs + mapIdx;
	char *pAttribName = pCtx->spAttribNames[STUC_ATTRIB_SP_MAT_IDX];
	const CommonAttrib *pCommonAttrib =
		stucGetCommonAttrib(pMapCommon->pFace, pMapCommon->faceCount, pAttribName);
	BlendConfig config = {0};
	if (pCommonAttrib) {
		config = pCommonAttrib->blendConfig;
	}
	else {
		config =
			stucGetTypeDefaultConfig(&pCtx->typeDefaults, STUC_ATTRIB_STRING)->blendConfig;
	}
	char *pName = "StucMaterials";
	const AttribIndexed *pInMats = stucGetAttribIndexedIntern(pInIndexedAttribs, pName);
	const AttribIndexed *pMapMats = stucGetAttribIndexedIntern(pMapAttribArr, pName);
	if (!pInMats) {
		*ppMatsToAdd = pMapMats;
	}
	else if (!pMapMats) {
		*ppMatsToAdd = pInMats;
	}
	else {
		*ppMatsToAdd = config.order ? pInMats : pMapMats;
	}
	return STUC_SUCCESS;
}

static
void appendToOutMatsBuf(
	const StucAlloc *pAlloc,
	Mesh *pMesh,
	AttribIndexed *pOutMats,
	const char *pMatName
) {
	STUC_ASSERT("", pOutMats->count <= pOutMats->size);
	if (pOutMats->count == pOutMats->size) {
		pOutMats->size *= 2;
		stucReallocAttrib(pAlloc, pMesh, &pOutMats->core, pOutMats->size);
	}
	memcpy(stucAttribAsStr(&pOutMats->core, pOutMats->count), pMatName,
	       STUC_ATTRIB_STRING_MAX_LEN);
	pOutMats->count++;
}

typedef struct {
	I8 idx;
	bool hasRef;
} MatTableEntry;

static
void appendToOutMats(
	AttribIndexed *pOutMats,
	char *pOutMatBuf,
	MatTableEntry *pEntry
) {
	pEntry->hasRef = true;
	char *pDest = stucAttribAsStr(&pOutMats->core, pOutMats->count);
	char *pSrc = pOutMatBuf + pEntry->idx * STUC_ATTRIB_STRING_MAX_LEN;
	memcpy(pDest, pSrc, STUC_ATTRIB_STRING_MAX_LEN);
	pEntry->idx = (I8)pOutMats->count;
	pOutMats->count++;
}

static
Result iterFacesAndCorrectMats(
	const StucAlloc *pAlloc,
	I32 i,
	Mesh *pMesh,
	AttribIndexed *pOutMats,
	MatTableEntry **ppMatTable,
	const AttribIndexed *pMatsToAdd
) {
	Result err = STUC_SUCCESS;
	for (I32 j = 0; j < pMesh->core.faceCount; ++j) {
		I32 matIdx = pMesh->pMatIdx[j];
		STUC_THROW_IFNOT_COND(err, matIdx >= 0 && matIdx < pMatsToAdd->count, "", 0);
		MatTableEntry *pEntry = ppMatTable[i] + matIdx;
		if (!pEntry->hasRef) {
			pEntry->hasRef = true;
			//We're just through material slots, so linear search should be fine for now
			const char *pMatName = stucAttribAsStrConst(&pMatsToAdd->core, matIdx);
			I32 idx = stucGetStrIdxInIndexedAttrib(pOutMats, pMatName);
			if (idx >= 0) {
				pEntry->idx = (I8)idx;
			}
			else {
				pEntry->idx = (I8)pOutMats->count;
				appendToOutMatsBuf(pAlloc, pMesh, pOutMats, pMatName);
			}
		}
		pMesh->pMatIdx[j] = pEntry->idx;
	}
	STUC_CATCH(0, err, ;);
	return err;
}

static
Result correctMatIndices(
	StucContext pCtx,
	Mesh *pMeshArr,
	const StucMapArr *pMapArr,
	const CommonAttribList *pCommonAttribs,
	const AttribIndexedArr *pInIndexedAttribs,
	AttribIndexedArr *pOutIndexedAttribs
) {
	Result err = STUC_SUCCESS;
	const StucAlloc *pAlloc = &pCtx->alloc;
	pOutIndexedAttribs->size = 1;
	pOutIndexedAttribs->count = 1;
	pOutIndexedAttribs->pArr =
		pAlloc->pCalloc(pOutIndexedAttribs->size, sizeof(AttribIndexed));
	AttribIndexed *pOutMats = pOutIndexedAttribs->pArr;
	pOutMats->size = pMapArr->count;
	pOutMats->size += (pOutMats->size % 2);
	stucInitAttribCore(
		pAlloc,
		&pOutMats->core,
		"StucMaterials",
		pOutMats->size,
		STUC_ATTRIB_STRING
	);
	MatTableEntry **ppMatTable = pAlloc->pCalloc(pMapArr->count, sizeof(void *));

	for (I32 i = 0; i < pMapArr->count; ++i) {
		Mesh *pMesh = pMeshArr + i;
		err = stucSetSpecialAttribs(pCtx, pMesh, 0x800);//only set mat indices
		STUC_THROW_IFNOT(err, "", 0);
		if (!pMesh->pMatIdx) {
			continue;
		}
		const AttribIndexed *pMatsToAdd = NULL;
		err = getMatsToAdd(
			pCtx,
			pMapArr,
			i,
			pCommonAttribs,
			pInIndexedAttribs,
			&pMatsToAdd
		);
		STUC_ASSERT("mesh faces have mat indices, but map has no materials?",
			pMatsToAdd && pMatsToAdd->count
		);
		STUC_THROW_IFNOT(err, "", 0);
		ppMatTable[i] = pAlloc->pCalloc(pMatsToAdd->count, sizeof(MatTableEntry));

		iterFacesAndCorrectMats(pAlloc, i, pMesh, pOutMats, ppMatTable, pMatsToAdd);
		STUC_CATCH(0, err, ;)
		pAlloc->pFree(ppMatTable[i]);
		ppMatTable[i] = NULL;
		STUC_THROW_IFNOT(err, "", 1);
	}
	STUC_CATCH(1, err, ;)
	pAlloc->pFree(ppMatTable);
	return err;
}

typedef struct {
	StucContext pCtx;
	StucMapArr *pMapArr;
	StucMesh *pMeshIn;
	StucAttribIndexedArr *pInIndexedAttribs;
	StucMesh *pMeshOut;
	StucAttribIndexedArr *pOutIndexedAttribs;
	StucCommonAttribList *pCommonAttribList;
	F32 wScale;
	F32 receiveLen;
} StucMapToMeshArgs;

static
Result mapToMeshFromJob(void *pArgsVoid) {
	StucMapToMeshArgs *pArgs = pArgsVoid;
	return stucMapToMesh(
		pArgs->pCtx,
		pArgs->pMapArr,
		pArgs->pMeshIn,
		pArgs->pInIndexedAttribs,
		pArgs->pMeshOut,
		pArgs->pOutIndexedAttribs,
		pArgs->pCommonAttribList,
		pArgs->wScale,
		pArgs->receiveLen
	);
}

Result stucQueueMapToMesh(
	StucContext pCtx,
	void **ppJobHandle,
	StucMapArr *pMapArr,
	StucMesh *pMeshIn,
	StucAttribIndexedArr *pInIndexedAttribs,
	StucMesh *pMeshOut,
	StucAttribIndexedArr *pOutIndexedAttribs,
	StucCommonAttribList *pCommonAttribList,
	F32 wScale,
	F32 receiveLen
) {
	StucMapToMeshArgs *pArgs = pCtx->alloc.pCalloc(1, sizeof(StucMapToMeshArgs));
	pArgs->pCtx = pCtx;
	pArgs->pMapArr = pMapArr;
	pArgs->pMeshIn = pMeshIn;
	pArgs->pInIndexedAttribs = pInIndexedAttribs;
	pArgs->pMeshOut = pMeshOut;
	pArgs->pOutIndexedAttribs = pOutIndexedAttribs;
	pArgs->pCommonAttribList = pCommonAttribList;
	pArgs->wScale = wScale;
	pArgs->receiveLen = receiveLen;
	pCtx->threadPool.pJobStackPushJobs(
		pCtx->pThreadPoolHandle,
		1,
		ppJobHandle,
		mapToMeshFromJob,
		&pArgs
	);
	return STUC_SUCCESS;
}

static
Result mapMapArrToMesh(
	StucContext pCtx,
	const StucMapArr *pMapArr,
	const Mesh *pMeshIn,
	const StucAttribIndexedArr *pInIndexedAttribs,
	StucMesh *pMeshOut,
	StucAttribIndexedArr *pOutIndexedAttribs,
	const StucCommonAttribList *pCommonAttribList,
	F32 wScale,
	F32 receiveLen
) {
	Result err = STUC_SUCCESS;
	Mesh *pOutBufArr = pCtx->alloc.pCalloc(pMapArr->count, sizeof(Mesh));
	StucObject *pOutObjWrapArr =
		pCtx->alloc.pCalloc(pMapArr->count, sizeof(StucObject));
	for (I32 i = 0; i < pMapArr->count; ++i) {
		pOutObjWrapArr[i].pData = (StucObjectData *)&pOutBufArr[i];
		const StucMap pMap = pMapArr->ppArr[i];
		I8 matIdx = pMapArr->pMatArr[i];
		InFaceArr *pInFaceTable = NULL;
		if (pMap->usgArr.count) {
			//set preserve to null to prevent usg squares from being split
			if (pMeshIn->pEdgePreserve || pMeshIn->pVertPreserve) {
				*(void **)&pMeshIn->pEdgePreserve = NULL;
				*(void **)&pMeshIn->pVertPreserve = NULL;
			}
			MapFile squares = {.pMesh = pMap->usgArr.pSquares};
			err = stucCreateQuadTree(pCtx, &squares.quadTree, squares.pMesh);
			STUC_THROW_IFNOT(err, "failed to create usg quadtree", 0);
			StucMesh squaresOut = {0};
			err = mapToMeshInternal(
				pCtx,
				&squares,
				pMeshIn,
				&squaresOut,
				matIdx,
				pCommonAttribList + i,
				&pInFaceTable,
				1.0f,
				-1.0f
			);
			STUC_THROW_IFNOT(err, "map to mesh usg failed", 1);
			err = stucSampleInAttribsAtUsgOrigins(
				pCtx,
				pMap,
				pMeshIn,
				&squaresOut,
				pInFaceTable
			);
			STUC_THROW_IFNOT(err, "", 1);
			InFaceTableToHashTable(&pCtx->alloc, pMap, squaresOut.faceCount, pInFaceTable);
			//*pMeshOut = squaresOut;
			//return STUC_SUCCESS;
			stucMeshDestroy(pCtx, &squaresOut);
			stucSetSpecialAttribs(pCtx, (Mesh *)pMeshIn, 0x50); //reassign preserve if present
		}
		err = mapToMeshInternal(
			pCtx,
			pMap,
			pMeshIn,
			&pOutBufArr[i].core,
			matIdx,
			pCommonAttribList + i,
			NULL,
			wScale,
			receiveLen
		);
		STUC_THROW_IFNOT(err, "map to mesh failed", 1);
		STUC_CATCH(1, err, ;);
		if (pMap->usgArr.count) {
			pCtx->alloc.pFree(pMap->usgArr.pInFaceTable);
			pMap->usgArr.pInFaceTable = NULL;
			for (I32 j = 0; j < pMap->usgArr.count; ++j) {
				//TODO uncomment this and fix any memory issue
				//pCtx->alloc.pFree(pInFaceTable[i].pArr);
			}
			pCtx->alloc.pFree(pInFaceTable);
		}
		STUC_THROW_IFNOT(err, "", 0);
	}
	pMeshOut->type.type = STUC_OBJECT_DATA_MESH;
	Mesh meshOutWrap = {.core = *pMeshOut};
	err = correctMatIndices(
		pCtx,
		pOutBufArr,
		pMapArr,
		pCommonAttribList,
		pInIndexedAttribs,
		pOutIndexedAttribs
	);
	STUC_THROW_IFNOT(err, "", 0);
	stucMergeObjArr(pCtx, &meshOutWrap, pMapArr->count, pOutObjWrapArr, false);
	*pMeshOut = meshOutWrap.core;
	STUC_CATCH(0, err,
		stucMeshDestroy(pCtx, pMeshOut);
	);
	//meshes are stored on an arr buf, which we can't call stucObjArrDestroy
	for (I32 i = 0; i < pMapArr->count; ++i) {
		stucMeshDestroy(pCtx, &pOutBufArr[i].core);
	}
	pCtx->alloc.pFree(pOutBufArr);
	pCtx->alloc.pFree(pOutObjWrapArr);
	return err;
}

static
Result appendSpAttribsToInMesh(
	const StucContext pCtx,
	Mesh *pWrap,
	const StucMesh *pMeshIn,
	UBitField32 flags
) {
	Result err = STUC_SUCCESS;
	UBitField32 has = 0;
	stucQuerySpAttribs(pCtx, pMeshIn, flags, &has);
	if (has) {
		STUC_RETURN_ERR(err, "in-mesh contains attribs it shouldn't");
	}
	Mesh meshInCpy = {.core = *pMeshIn};
	Mesh *pMeshInCpyPtr = &meshInCpy;
	stucAllocAttribsFromMeshArr(
		&pCtx->alloc,
		pWrap,
		1, &pMeshInCpyPtr,
		false,
		false, //dont allocation data
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
	for (I32 i = 1; i < STUC_ATTRIB_SP_COUNT; ++i) {
		if (!(flags >> i & 0x1)) {
			continue;
		}
		Attrib *pAttrib = stucGetSpAttrib(pCtx, pMesh, i);
		if (pAttrib) {
			if (pAttrib->core.pData) {
				pCtx->alloc.pFree(pAttrib->core.pData);
				pAttrib->core.pData = NULL;
			}
		}
	}
	if (pMesh->faceAttribs.pArr) {
		pCtx->alloc.pFree(pMesh->faceAttribs.pArr);
		pMesh->faceAttribs.pArr = NULL;
	}
	if (pMesh->cornerAttribs.pArr) {
		pCtx->alloc.pFree(pMesh->cornerAttribs.pArr);
		pMesh->cornerAttribs.pArr = NULL;
	}
	if (pMesh->edgeAttribs.pArr) {
		pCtx->alloc.pFree(pMesh->edgeAttribs.pArr);
		pMesh->edgeAttribs.pArr = NULL;
	}
	if (pMesh->vertAttribs.pArr) {
		pCtx->alloc.pFree(pMesh->vertAttribs.pArr);
		pMesh->vertAttribs.pArr = NULL;
	}
}

static
Result initMeshInWrap(
	StucContext pCtx,
	Mesh *pWrap,
	StucMesh *pMeshIn,
	UBitField32 spAttribsToAppend
) {
	Result err = STUC_SUCCESS;
	if (!pMeshIn->edgeCount) {
		STUC_RETURN_ERR_IFNOT_COND(
			err,
			!pMeshIn->edgeAttribs.count,
			"in-mesh has edge attribs, yet no edge list"
		);
		err = stucBuildEdgeList(pCtx, pWrap);
		STUC_RETURN_ERR_IFNOT(err, "failed to build edge list");
	}
	//const qualifiers removed here 
	stucAliasMeshCoreNoAttribs(&pWrap->core, pMeshIn);
	err = appendSpAttribsToInMesh(pCtx, pWrap, pMeshIn, spAttribsToAppend);
	STUC_RETURN_ERR_IFNOT(err, "");
	stucSetAttribOrigins(&pWrap->core.meshAttribs, STUC_ATTRIB_ORIGIN_MESH_IN);
	stucSetAttribOrigins(&pWrap->core.faceAttribs, STUC_ATTRIB_ORIGIN_MESH_IN);
	stucSetAttribOrigins(&pWrap->core.cornerAttribs, STUC_ATTRIB_ORIGIN_MESH_IN);
	stucSetAttribOrigins(&pWrap->core.edgeAttribs, STUC_ATTRIB_ORIGIN_MESH_IN);
	stucSetAttribOrigins(&pWrap->core.vertAttribs, STUC_ATTRIB_ORIGIN_MESH_IN);

	err = stucSetSpecialAttribs(pCtx, pWrap, 0x1ef5e);
	STUC_RETURN_ERR_IFNOT(err, "");

	err = stucBuildTangents(pWrap);
	STUC_RETURN_ERR_IFNOT(err, "failed to build tangents");
	buildEdgeCornersTable(pWrap);
	buildSeamAndPreserveTables(&pCtx->alloc, pWrap);

	//set all but pos, normal, uv, w-scale, and mat-idx
	stucSetAttribToDontCopy(pCtx, pWrap, 0x1f3f0);

	return err;
}

Result stucMapToMesh(
	StucContext pCtx,
	const StucMapArr *pMapArr,
	const StucMesh *pMeshIn,
	const StucAttribIndexedArr *pInIndexedAttribs,
	StucMesh *pMeshOut,
	StucAttribIndexedArr *pOutIndexedAttribs,
	const StucCommonAttribList *pCommonAttribList,
	F32 wScale,
	F32 receiveLen
) {
	StucResult err = STUC_SUCCESS;
	STUC_RETURN_ERR_IFNOT_COND(err, pMeshIn, "");
	err = stucValidateMesh(pMeshIn, false);
	STUC_RETURN_ERR_IFNOT(err, "invalid in-mesh");
	Mesh meshInWrap = {0};
	UBitField32 spAttribsToAppend = 0x1e300;
	err = initMeshInWrap(
		pCtx,
		&meshInWrap,
		(StucMesh *)pMeshIn,
		spAttribsToAppend
	);
	STUC_THROW_IFNOT(err, "", 0);

	STUC_THROW_IFNOT_COND(
		err,
		pMapArr && pMapArr->count && pMapArr->pMatArr && pMapArr->ppArr,
		"", 0
	);

	err = mapMapArrToMesh(
		pCtx,
		pMapArr,
		&meshInWrap,
		pInIndexedAttribs,
		pMeshOut,
		pOutIndexedAttribs,
		pCommonAttribList,
		wScale,
		receiveLen
	);
	STUC_THROW_IFNOT(err, "mapMapArrToMesh returned error", 0);
	printf("----------------------FINISHING IN-MESH\n");
	STUC_CATCH(0, err, ;);
	destroyAppendedSpAttribs(pCtx, &meshInWrap.core, spAttribsToAppend);
	return err;
}

Result stucObjArrDestroy(
	StucContext pCtx,
	I32 objCount,
	StucObject *pObjArr
) {
	return stucDestroyObjArr(pCtx, objCount, pObjArr);
}

Result stucUsgArrDestroy(StucContext pCtx, I32 count, StucUsg *pUsgArr) {
	Result err = STUC_NOT_SET;
	for (I32 i = 0; i < count; ++i) {
		err = stucMeshDestroy(pCtx, (StucMesh *)pUsgArr[i].obj.pData);
		STUC_THROW_IFNOT(err, "", 0);
	}
	pCtx->alloc.pFree(pUsgArr);
	STUC_CATCH(0, err, ;)
	return err;
}

StucResult stucMeshDestroy(StucContext pCtx, StucMesh *pMesh) {
	for (I32 i = 0; i < pMesh->meshAttribs.count; ++i) {
		if (pMesh->meshAttribs.pArr[i].core.pData) {
			pCtx->alloc.pFree(pMesh->meshAttribs.pArr[i].core.pData);
		}
	}
	if (pMesh->meshAttribs.count && pMesh->meshAttribs.pArr) {
		pCtx->alloc.pFree(pMesh->meshAttribs.pArr);
	}
	if(pMesh->pFaces) {
		pCtx->alloc.pFree(pMesh->pFaces);
	}
	for (I32 i = 0; i < pMesh->faceAttribs.count; ++i) {
		if (pMesh->faceAttribs.pArr[i].core.pData) {
			pCtx->alloc.pFree(pMesh->faceAttribs.pArr[i].core.pData);
		}
	}
	if (pMesh->faceAttribs.count && pMesh->faceAttribs.pArr) {
		pCtx->alloc.pFree(pMesh->faceAttribs.pArr);
	}
	if (pMesh->pCorners) {
		pCtx->alloc.pFree(pMesh->pCorners);
	}
	for (I32 i = 0; i < pMesh->cornerAttribs.count; ++i) {
		if (pMesh->cornerAttribs.pArr[i].core.pData) {
			pCtx->alloc.pFree(pMesh->cornerAttribs.pArr[i].core.pData);
		}
	}
	if (pMesh->cornerAttribs.count && pMesh->cornerAttribs.pArr) {
		pCtx->alloc.pFree(pMesh->cornerAttribs.pArr);
	}
	if (pMesh->pEdges) {
		pCtx->alloc.pFree(pMesh->pEdges);
	}
	for (I32 i = 0; i < pMesh->edgeAttribs.count; ++i) {
		if (pMesh->edgeAttribs.pArr[i].core.pData) {
			pCtx->alloc.pFree(pMesh->edgeAttribs.pArr[i].core.pData);
		}
	}
	if (pMesh->edgeAttribs.count && pMesh->edgeAttribs.pArr) {
		pCtx->alloc.pFree(pMesh->edgeAttribs.pArr);
	}
	for (I32 i = 0; i < pMesh->vertAttribs.count; ++i) {
		if (pMesh->vertAttribs.pArr[i].core.pData) {
			pCtx->alloc.pFree(pMesh->vertAttribs.pArr[i].core.pData);
		}
	}
	if (pMesh->vertAttribs.count && pMesh->vertAttribs.pArr) {
		pCtx->alloc.pFree(pMesh->vertAttribs.pArr);
	}
	return STUC_SUCCESS;
}

StucResult stucGetAttribSize(StucAttrib *pAttrib, I32 *pSize) {
	*pSize = stucGetAttribSizeIntern(pAttrib->core.type);
	return STUC_SUCCESS;
}

StucResult stucGetAttrib(char *pName, StucAttribArray *pAttribs, StucAttrib **ppAttrib) {
	*ppAttrib = stucGetAttribIntern(pName, pAttribs);
	return STUC_SUCCESS;
}

StucResult stucGetAttribIndexed(
	char *pName,
	StucAttribIndexedArr *pAttribs,
	StucAttribIndexed **ppAttrib
) {
	*ppAttrib = stucGetAttribIndexedIntern(pAttribs, pName);
	return STUC_SUCCESS;
}

typedef struct {
	StucImage imageBuf;
	StucMap pMap;
	StucContext pCtx;
	I32 *pActiveJobs;
	void *pMutex;
	I32 bufOffset;
	I32 pixelCount;
	I8 id;
	V2_F32 zBounds;
} RenderArgs;

static
void testPixelAgainstFace(
	RenderArgs *pArgs,
	V2_F32 *pPos,
	FaceRange *pFace,
	Color *pColor
) {
	const Mesh *pMesh = pArgs->pMap->pMesh;
	V2_F32 verts[4] = {0};
	for (I32 i = 0; i < pFace->size; ++i) {
		verts[i] = *(V2_F32 *)(pMesh->pVerts + pMesh->core.pCorners[pFace->start + i]);
	}
	I8 triCorners[4] = {0};
	V3_F32 bc = stucGetBarycentricInFace(verts, triCorners, pFace->size, *pPos);
	if (bc.d[0] < -.0001f || bc.d[1] < -.0001f || bc.d[2] < -.0001f ||
		!v3IsFinite(bc)) {
		return;
	}
	V3_F32 vertsXyz[3] = {0};
	for (I32 i = 0; i < 3; ++i) {
		I32 vertIdx = pMesh->core.pCorners[pFace->start + triCorners[i]];
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
	F32 normalLen = sqrtf(
		normal.d[0] * normal.d[0] +
		normal.d[1] * normal.d[1] +
		normal.d[2] * normal.d[2]
	);
	_(&normal V3DIVEQLS normalLen);
	V3_F32 up = { .0f, .0f, 1.0f };
	F32 dotUp = _(normal V3DOT up);
	if (dotUp < .0f) {
		return;
	}
	F32 depth = (wsPos.d[2] - pArgs->zBounds.d[0]) / pArgs->zBounds.d[1];
	F32 value = dotUp;
	value *= 1.0f - (1.0f - depth) * .5f;
	value *= .75;
	pColor->d[0] = value;
	pColor->d[1] = value;
	pColor->d[2] = value;
	pColor->d[3] = wsPos.d[2];
}

static
void testPixelAgainstCellFaces(
	RenderArgs *pArgs,
	const Mesh *pMesh,
	Cell *pLeaf,
	I32 j,
	FaceCells *pFaceCells,
	V2_F32 *pPos,
	Color *pColor
) {
	I32 cellIdx = pFaceCells->pCells[j];
	Cell *pCell = pArgs->pMap->quadTree.cellTable.pArr + cellIdx;
	I32* pCellFaces;
	Range cellFaceRange = {0};
	if (pFaceCells->pCellType[j]) {
		pCellFaces = pCell->pEdgeFaces;
		cellFaceRange = pLeaf->pLinkEdgeRanges[j];
	}
	else if (pFaceCells->pCellType[j] != 1) {
		pCellFaces = pCell->pFaces;
		cellFaceRange.start = 0;
		cellFaceRange.end = pCell->faceSize;
	}
	else {
		return;
	}
	for (I32 k = cellFaceRange.start; k < cellFaceRange.end; ++k) {
		FaceRange face = {0};
		face.idx = pCellFaces[k];
		face.start = pMesh->core.pFaces[face.idx];
		face.end = pMesh->core.pFaces[face.idx + 1];
		face.size = face.end - face.start;
		FaceTriangulated faceTris = {0};
		if (face.size > 4) {
			faceTris = stucTriangulateFace(
				pArgs->pCtx->alloc,
				&face,
				pMesh->pVerts,
				pMesh->core.pCorners,
				0
			);
			for (I32 l = 0; l < faceTris.triCount; ++l) {
				FaceRange tri = {0};
				tri.idx = face.idx;
				tri.start = l * 3;
				tri.end = tri.start + 3;
				tri.size = tri.end - tri.start;
				tri.start += face.start;
				tri.end += face.start;
				testPixelAgainstFace(pArgs, pPos, &tri, pColor);
			}
		}
		else {
			testPixelAgainstFace(pArgs, pPos, &face, pColor);
		}
	}
}

static
Result stucRenderJob(void *pArgs) {
	Result err = STUC_SUCCESS;
	RenderArgs args = *(RenderArgs *)pArgs;
	I32 dataLen = args.pixelCount * getPixelSize(args.imageBuf.type);
	args.imageBuf.pData = args.pCtx->alloc.pMalloc(dataLen);
	const Mesh *pMesh = args.pMap->pMesh;
	F32 pixelScale = 1.0f / (F32)args.imageBuf.res;
	F32 pixelHalfScale = pixelScale / 2.0f;
	FaceCells faceCells = {0};
	FaceCellsTable faceCellsTable = {.pFaceCells = &faceCells};
	QuadTreeSearch searchState = {.pAlloc = &args.pCtx->alloc, .pMap = args.pMap};
	stucInitQuadTreeSearch(&searchState);
	for (I32 i = 0; i < args.pixelCount; ++i) {
		I32 iOffset = args.bufOffset + i;
		V2_F32 idx = {
			(F32)(iOffset % args.imageBuf.res),
			(F32)(iOffset / args.imageBuf.res)
		};
		V2_F32 pos = {
			pixelScale * idx.d[0] + pixelHalfScale,
			pixelScale * idx.d[1] + pixelHalfScale
		};
		Color color = { 0 };
		color.d[3] = FLT_MAX * -1.0f;
		Range faceRange = {.start = 0, .end = 1};
		stucGetCellsForSingleFace(
			&searchState,
			1,
			&pos,
			&faceCellsTable,
			NULL,
			0,
			faceRange
		);
		I32 leafIdx = faceCells.pCells[faceCells.cellSize - 1];
		Cell *pLeaf = args.pMap->quadTree.cellTable.pArr + leafIdx;
		for (I32 j = 0; j < faceCells.cellSize; ++j) {
			testPixelAgainstCellFaces(&args, pMesh, pLeaf, j, &faceCells, &pos, &color);
		}
		color.d[3] = (F32)(color.d[3] != FLT_MAX) * -1.0f;
		setPixelColor(&args.imageBuf, i, &color);
	}
	stucDestroyQuadTreeSearch(&searchState);
	*(RenderArgs *)pArgs = args;
	StucThreadPool *pThreadPool = &args.pCtx->threadPool;
	pThreadPool->pMutexLock(args.pCtx->pThreadPoolHandle, args.pMutex);
	--*args.pActiveJobs;
	pThreadPool->pMutexUnlock(args.pCtx->pThreadPoolHandle, args.pMutex);
	return err;
}

static
V2_F32 getZBounds(StucMap pMap) {
	const Mesh *pMesh = pMap->pMesh;
	V2_F32 zBounds = {.d = {FLT_MAX, FLT_MIN}};
	for (I32 i = 0; i < pMesh->core.vertCount; ++i) {
		if (pMesh->pVerts[i].d[2] < zBounds.d[0]) {
			zBounds.d[0] = pMesh->pVerts[i].d[2];
		}
		if (pMesh->pVerts[i].d[2] > zBounds.d[1]) {
			zBounds.d[1] = pMesh->pVerts[i].d[2];
		}
	}
	return zBounds;
}

StucResult stucMapFileGenPreviewImage(
	StucContext pCtx,
	StucMap pMap,
	StucImage *pImage
) {
	Result err = STUC_SUCCESS;
	V2_F32 zBounds = getZBounds(pMap);
	I32 pixelCount = pImage->res * pImage->res;
	I32 pixelsPerJob = pixelCount / pCtx->threadCount;
	bool singleThread = !pixelsPerJob;
	void *pMutex = NULL;
	pCtx->threadPool.pMutexGet(pCtx->pThreadPoolHandle, &pMutex);
	I32 activeJobs = 0;
	void *jobArgPtrs[MAX_THREADS] = {0};
	RenderArgs args[MAX_THREADS] = {0};
	activeJobs = singleThread ? 1 : pCtx->threadCount;
	for (I32 i = 0; i < activeJobs; ++i) {
		args[i].bufOffset = i * pixelsPerJob;
		args[i].imageBuf.res = pImage->res;
		args[i].imageBuf.type = pImage->type;
		args[i].pCtx = pCtx;
		args[i].pixelCount = i == activeJobs - 1 ?
			pixelCount - args[i].bufOffset : pixelsPerJob;
		args[i].pMap = pMap;
		args[i].zBounds = zBounds;
		args[i].pMutex = pMutex;
		args[i].pActiveJobs = &activeJobs;
		args[i].id = (I8)i;
		jobArgPtrs[i] = args + i;
	}
	void **ppJobHandles = pCtx->alloc.pCalloc(activeJobs, sizeof(void *));
	pCtx->threadPool.pJobStackPushJobs(
		pCtx->pThreadPoolHandle,
		activeJobs,
		ppJobHandles,
		stucRenderJob,
		jobArgPtrs
	);
	stucWaitForJobsIntern(pCtx->pThreadPoolHandle, activeJobs, ppJobHandles, true, NULL);
	err = stucJobGetErrs(pCtx, activeJobs, &ppJobHandles);
	stucJobDestroyHandles(pCtx, activeJobs, ppJobHandles);
	pCtx->alloc.pFree(ppJobHandles);
	STUC_THROW_IFNOT(err, "", 0);
	I32 pixelSize = getPixelSize(pImage->type);
	pImage->pData = pCtx->alloc.pMalloc(pixelCount * pixelSize);
	for (I32 i = 0; i < activeJobs; ++i) {
		void *pImageOffset = offsetImagePtr(pImage, i * pixelsPerJob);
		I32 bytesToCopy = pixelSize * args[i].pixelCount;
		memcpy(pImageOffset, args[i].imageBuf.pData, bytesToCopy);
	}
	STUC_CATCH(0, err, ;);
	pCtx->threadPool.pMutexDestroy(pCtx->pThreadPoolHandle, pMutex);
	for (I32 i = 0; i < activeJobs; ++i) {
		pCtx->alloc.pFree(args[i].imageBuf.pData);
	}
	return err;
}

void stucMapIndexedAttribsGet(
	StucContext pCtx,
	StucMap pMap,
	StucAttribIndexedArr *pIndexedAttribs
) {
	*pIndexedAttribs = pMap->indexedAttribs;
}

Result stucWaitForJobs(
	StucContext pCtx,
	I32 count,
	void **ppHandles,
	bool wait,
	bool *pDone
) {
	return pCtx->threadPool.pWaitForJobs(
		pCtx->pThreadPoolHandle,
		count,
		ppHandles,
		wait,
		pDone
	);
}

Result stucJobGetErrs(
	StucContext pCtx,
	I32 jobCount,
	void ***pppJobHandles
) {
	Result err = STUC_SUCCESS;
	STUC_ASSERT("", pCtx && pppJobHandles);
	STUC_ASSERT("", jobCount > 0);
	for (I32 i = 0; i < jobCount; ++i) {
		StucResult jobErr = STUC_NOT_SET;
		err = pCtx->threadPool.pGetJobErr(
			pCtx->pThreadPoolHandle,
			(*pppJobHandles)[i],
			&jobErr
		);
		STUC_THROW_IFNOT_COND(err, jobErr == STUC_SUCCESS, "", 0);
	}
	STUC_CATCH(0, err, ;);
	return err;
}

void stucJobDestroyHandles(
	StucContext pCtx,
	I32 jobCount,
	void **ppJobHandles
) {
	STUC_ASSERT("", pCtx && ppJobHandles);
	STUC_ASSERT("", jobCount > 0);
	for (I32 i = 0; i < jobCount; ++i) {
		pCtx->threadPool.pJobHandleDestroy(
			pCtx->pThreadPoolHandle,
			ppJobHandles + i
		);
	}
}