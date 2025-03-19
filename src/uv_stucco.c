/* 
SPDX-FileCopyrightText: 2025 Caleb Dawson
SPDX-License-Identifier: Apache-2.0
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include <io.h>
#include <map_to_job_mesh.h>
#include <combine_job_meshes.h>
#include <map.h>
#include <context.h>
#include <alloc.h>
#include <thread_pool.h>
#include <uv_stucco.h>
#include <attrib_utils.h>
#include <utils.h>
#include <image_utils.h>
#include <error.h>

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
//TODO make naming for MeshIn consistent

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
	const char *pPath,
	I32 objCount,
	StucObject *pObjArr,
	I32 usgCount,
	StucUsg *pUsgArr,
	StucAttribIndexedArr *pIndexedAttribs
) {
	return stucWriteStucFile(
		pCtx,
		pPath,
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
	const char *filePath,
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
		//this occurs usually due to interior faces,
		// it shouldn't be an issue for for map-meshes, more so for in-meshes.
		//TODO remove this assert if no issues arise
		//STUC_ASSERT("more than 2 corners refernce 1 edge", pSet[edge] < 2);
		V3_F32 diff = _(pos V3SUB pPosCache[edge]);
		pMesh->pEdgeLen[edge] = v3F32Len(diff);
		pSet[edge]++;
	}
	pCtx->alloc.pFree(pSet);
	pCtx->alloc.pFree(pPosCache);
}

static
void TEMPsetSpFromAttribName(StucContext pCtx, StucMesh *pMesh, AttribArray *pArr) {
	for (I32 j = 0; j < pArr->count; ++j) {
		if (!strcmp(pArr->pArr[j].core.name, "StucMaterialIndices")) {
			strncpy(pArr->pArr[j].core.name, "StucMaterials", STUC_ATTRIB_NAME_MAX_LEN);
			pArr->pArr[j].core.use = STUC_ATTRIB_USE_IDX;
			pMesh->activeAttribs[STUC_ATTRIB_USE_IDX].active = true;
			pMesh->activeAttribs[STUC_ATTRIB_USE_IDX].idx = (I16)j;
		}
		else if (!strcmp(pArr->pArr[j].core.name, "Color")) {
			pArr->pArr[j].core.use = STUC_ATTRIB_USE_COLOR;
			pMesh->activeAttribs[STUC_ATTRIB_USE_COLOR].active = true;
			pMesh->activeAttribs[STUC_ATTRIB_USE_COLOR].idx = (I16)j;
		}
		else {
			for (I32 k = 1; k < STUC_ATTRIB_USE_SP_ENUM_COUNT; ++k) {
				if (!strncmp(pArr->pArr[j].core.name, pCtx->spAttribNames[k], STUC_ATTRIB_NAME_MAX_LEN)) {
					pArr->pArr[j].core.use = k;
					pMesh->activeAttribs[k].active = true;
					pMesh->activeAttribs[k].idx = (I16)j;
					break;
				}
			}
		}
	}
}

static
Result attemptToSetMissingActiveDomains(StucMesh *pMesh) {
	Result err = STUC_SUCCESS;
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
			STUC_RETURN_ERR_IFNOT_COND(
				err,
				pIdx->domain == STUC_DOMAIN_NONE,
				"Unable to determine active attrib domain"
			);
			pIdx->domain = j;
		}
	}
	return err;
}

StucResult stucMapFileLoad(StucContext pCtx, StucMap *pMapHandle, const char *filePath) {
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
		//TODO TEMP DELETE
		//TODO when you fix this, make sure map uvs arn't marked special,
		//or it'll conflit 
		Mesh *pMesh = (Mesh *)pObjArr[i].pData;
		TEMPsetSpFromAttribName(pCtx, &pMesh->core, &pMesh->core.faceAttribs);
		TEMPsetSpFromAttribName(pCtx, &pMesh->core, &pMesh->core.cornerAttribs);
		TEMPsetSpFromAttribName(pCtx, &pMesh->core, &pMesh->core.edgeAttribs);
		TEMPsetSpFromAttribName(pCtx, &pMesh->core, &pMesh->core.vertAttribs);

		//TODO this shouldn't be an issue.
		//check for this on export, so this doesn't throw and error on loaad
		attemptToSetMissingActiveDomains(&pMesh->core);

		err = stucAssignActiveAliases(
			pCtx,
			(Mesh *)pObjArr[i].pData,
			0x8ae, //10101110 - all except for preserve
			STUC_DOMAIN_NONE
		);
		STUC_THROW_IFNOT(err, "", 0);
		stucApplyObjTransform(pObjArr + i);
	}
	Mesh *pMapMesh = pCtx->alloc.pCalloc(1, sizeof(Mesh));
	pMapMesh->core.type.type = STUC_OBJECT_DATA_MESH_INTERN;
	err = stucMergeObjArr(pCtx, pMapMesh, objCount, pObjArr, false);
	STUC_THROW_IFNOT(err, "", 0);

	//append edgeLen attrib
	stucAppendSpAttribsToMesh(pCtx, pMapMesh, 0x1000, STUC_ATTRIB_ORIGIN_MAP);

	stucSetAttribOrigins(&pMapMesh->core.meshAttribs, STUC_ATTRIB_ORIGIN_MAP);
	stucSetAttribOrigins(&pMapMesh->core.faceAttribs, STUC_ATTRIB_ORIGIN_MAP);
	stucSetAttribOrigins(&pMapMesh->core.cornerAttribs, STUC_ATTRIB_ORIGIN_MAP);
	stucSetAttribOrigins(&pMapMesh->core.edgeAttribs, STUC_ATTRIB_ORIGIN_MAP);
	stucSetAttribOrigins(&pMapMesh->core.vertAttribs, STUC_ATTRIB_ORIGIN_MAP);

	stucSetAttribCopyOpt(pCtx, &pMapMesh->core, STUC_ATTRIB_DONT_COPY, 0x17f0);
	err = stucAssignActiveAliases(pCtx, pMapMesh, 0x18ae, STUC_DOMAIN_NONE);
	STUC_THROW_IFNOT(err, "", 0);

	buildEdgeLenList(pCtx, pMapMesh);

	//TODO some form of heap corruption when many objects
	//test with address sanitizer on CircuitPieces.stuc
	stucDestroyObjArr(pCtx, objCount, pObjArr);

	{
		Attrib *pUvAttrib = stucGetActiveAttrib(pCtx, &pMapMesh->core, STUC_ATTRIB_USE_UV);
		if (pUvAttrib) {
			//TODO as with all special attributes, allow user to define what should be considered
			//     the primary UV channel. This especially important for integration with other DCCs
			if (!strncmp(
				pUvAttrib->core.name,
				pCtx->spAttribNames[STUC_ATTRIB_USE_UV],
				STUC_ATTRIB_NAME_MAX_LEN
			)) {
				char newName[STUC_ATTRIB_NAME_MAX_LEN] = "Map_UVMap";
				memcpy(pUvAttrib->core.name, newName, STUC_ATTRIB_NAME_MAX_LEN);
			}
		}
		pMapMesh->core.activeAttribs[STUC_ATTRIB_USE_UV].active = false;
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
			Mesh *pUsgMesh = (Mesh *)pUsgArr[i].obj.pData;
			TEMPsetSpFromAttribName(pCtx, &pUsgMesh->core, &pUsgMesh->core.faceAttribs);
			TEMPsetSpFromAttribName(pCtx, &pUsgMesh->core, &pUsgMesh->core.cornerAttribs);
			TEMPsetSpFromAttribName(pCtx, &pUsgMesh->core, &pUsgMesh->core.edgeAttribs);
			TEMPsetSpFromAttribName(pCtx, &pUsgMesh->core, &pUsgMesh->core.vertAttribs);
			attemptToSetMissingActiveDomains(&pUsgMesh->core);
			err = stucAssignActiveAliases(
				pCtx,
				pUsgMesh,
				0x02, //000010 - set only vert pos
				STUC_DOMAIN_NONE
			);
			STUC_THROW_IFNOT(err, "", 0);
			pMap->usgArr.pArr[i].origin = *(V2_F32 *)&pUsgArr[i].obj.transform.d[3];
			pMap->usgArr.pArr[i].pMesh = pUsgMesh;
			stucApplyObjTransform(&pUsgArr[i].obj);
			if (pUsgArr[i].pFlatCutoff) {
				Mesh *pFlatCutoff = (Mesh *)pUsgArr[i].pFlatCutoff->pData;
				pMap->usgArr.pArr[i].pFlatCutoff = (Mesh *)pUsgArr[i].pFlatCutoff->pData;
				TEMPsetSpFromAttribName(pCtx, &pFlatCutoff->core, &pFlatCutoff->core.faceAttribs);
				TEMPsetSpFromAttribName(pCtx, &pFlatCutoff->core, &pFlatCutoff->core.cornerAttribs);
				TEMPsetSpFromAttribName(pCtx, &pFlatCutoff->core, &pFlatCutoff->core.edgeAttribs);
				TEMPsetSpFromAttribName(pCtx, &pFlatCutoff->core, &pFlatCutoff->core.vertAttribs);
				attemptToSetMissingActiveDomains(&pFlatCutoff->core);
				err = stucAssignActiveAliases(
					pCtx,
					pFlatCutoff,
					0x02, //000010 - set only vert pos
					STUC_DOMAIN_NONE
				);
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

Result stucMapFileMeshGet(StucContext pCtx, StucMap pMap, const StucMesh **ppMesh) {
	Result err = STUC_SUCCESS;
	STUC_RETURN_ERR_IFNOT_COND(err, pCtx && pMap && ppMesh, "invalid args");
	*ppMesh = &pMap->pMesh->core;
	return err;
}

static
void initCommonAttrib(
	StucContext pCtx,
	StucCommonAttrib *pEntry,
	const StucAttrib *pAttrib
) {
	memcpy(pEntry->name, pAttrib->core.name, STUC_ATTRIB_NAME_MAX_LEN);
	StucTypeDefault *pDefault = 
		stucGetTypeDefaultConfig(&pCtx->typeDefaults, pAttrib->core.type);
	pEntry->blendConfig = pDefault->blendConfig;
}

static
Result getCommonAttribs(
	StucContext pCtx,
	const StucMesh *pMapMesh,
	const AttribArray *pMapAttribs,
	const StucMesh *pMesh,
	const AttribArray *pMeshAttribs,
	StucCommonAttribArr *pCommonArr
) {
	Result err = STUC_SUCCESS;
	//TODO ignore special attribs like StucTangent or StucTSign
	pCommonArr->count = 0;
	pCommonArr->size = 2;
	pCommonArr->pArr = pCtx->alloc.pCalloc(pCommonArr->size, sizeof(StucCommonAttrib));
	for (I32 i = 0; i < pMeshAttribs->count; ++i) {
		Attrib *pAttrib = pMeshAttribs->pArr + i;
		Attrib *pMapAttrib = NULL;
		err = stucGetMatchingAttribConst(
			pCtx,
			pMapMesh, pMapAttribs,
			pMesh, pAttrib,
			true,
			true,
			&pMapAttrib
		);
		STUC_THROW_IFNOT(err, "", 0);
		if (!pMapAttrib) {
			continue;
		}
		STUC_ASSERT("", pCommonArr->count <= pCommonArr->size);
		if (pCommonArr->count == pCommonArr->size) {
			pCommonArr->size *= 2;
			pCommonArr->pArr = pCtx->alloc.pRealloc(pCommonArr->pArr, pCommonArr->size);
		}
		initCommonAttrib(pCtx, pCommonArr->pArr + pCommonArr->count, pAttrib);
		pCommonArr->count++;
	}
	STUC_CATCH(0, err,
		pCtx->alloc.pFree(pCommonArr->pArr);
		pCommonArr->count = pCommonArr->size = 0;
	);
	return err;
}

//TODO handle edge case, where attribute share the same name,
//but have incompatible types. Such as a F32 and a string.
StucResult stucQueryCommonAttribs(
	StucContext pCtx,
	const StucMap pMap,
	const StucMesh *pMesh,
	StucCommonAttribList *pCommonAttribs
) {
	Result err = STUC_SUCCESS;
	STUC_RETURN_ERR_IFNOT_COND(err, pCtx && pMap && pMesh && pCommonAttribs, "");
	const StucMesh *pMapMesh = &pMap->pMesh->core;
	StucMesh meshWrap = *pMesh;
	err = attemptToSetMissingActiveDomains(&meshWrap);
	STUC_RETURN_ERR_IFNOT(err, "");
	for (I32 i = STUC_DOMAIN_FACE; i <= STUC_DOMAIN_MESH; ++i) {
		StucCommonAttribArr *pCommonArr = NULL;
		stucCommonAttribArrGetFromDomain(pCtx, pCommonAttribs, i, &pCommonArr);
		err = getCommonAttribs(
			pCtx,
			pMapMesh,
			stucGetAttribArrFromDomainConst(pMapMesh, i),
			&meshWrap,
			stucGetAttribArrFromDomainConst(&meshWrap, i),
			pCommonArr
		);
		STUC_THROW_IFNOT(err, "", 0);
	}
	STUC_CATCH(0, err,
		stucDestroyCommonAttribs(pCtx, pCommonAttribs);
	);
	return err;
}

StucResult stucCommonAttribArrGetFromDomain(
	StucContext pCtx,
	StucCommonAttribList *pList,
	StucDomain domain,
	StucCommonAttribArr **ppArr
) {
	Result err = STUC_SUCCESS;
	switch (domain) {
	case STUC_DOMAIN_FACE:
		*ppArr = &pList->face;
		break;
	case STUC_DOMAIN_CORNER:
		*ppArr = &pList->corner;
		break;
	case STUC_DOMAIN_EDGE:
		*ppArr = &pList->edge;
		break;
	case STUC_DOMAIN_VERT:
		*ppArr = &pList->vert;
		break;
	case STUC_DOMAIN_MESH:
		*ppArr = &pList->mesh;
		break;
	default:
		STUC_RETURN_ERR(err, "invalid domain");
	}
	return err;
}

StucResult stucDestroyCommonAttribs(
	StucContext pCtx,
	StucCommonAttribList *pCommonAttribs
) {
	Result err = STUC_SUCCESS;
	STUC_RETURN_ERR_IFNOT_COND(err, pCtx && pCommonAttribs, "");
	for (I32 i = STUC_DOMAIN_FACE; i <= STUC_DOMAIN_MESH; ++i) {
		StucCommonAttribArr *pArr = NULL;
		stucCommonAttribArrGetFromDomain(pCtx, pCommonAttribs, i, &pArr);
		if (pArr) {
			pCtx->alloc.pFree(pArr->pArr);
			pArr->pArr = NULL;
		}
		pArr->count = pArr->size = 0;
	}
	return err;
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
Result mapToFaces(
	MapToMeshBasic *pBasic,
	I32 *pJobCount,
	MappingJobArgs **ppJobArgs,
	bool *pEmpty
) {
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
	printf("\n");
	for (I32 i = 0; i < *pJobCount; ++i) {
		printf("	job %d checked %d faces and used %d\n",
			i,
			(*ppJobArgs)[i].facesChecked,
			(*ppJobArgs)[i].facesUsed
		);
		if ((*ppJobArgs)[i].bufSize > 0) {
			*pEmpty = false;
			//break;
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
	InFaceTable *pInFaceTable,
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
		.pInFaceTable = pInFaceTable,
	};
	if (pInFaceTable) {
		stucLinAllocInit(
			&pCtx->alloc,
			&pInFaceTable->pAlloc,
			sizeof(I32),
			pMeshIn->core.faceCount
		);
	}
	I32 jobCount = 0;
	MappingJobArgs *pMappingJobArgs = NULL;

	bool empty = true;
	err = mapToFaces(&basic, &jobCount, &pMappingJobArgs, &empty);
	STUC_THROW_IFNOT(err, "", 0);
	if (!empty) {
		err = stucCombineJobMeshes(&basic, pMappingJobArgs, jobCount);
		STUC_THROW_IFNOT(err, "", 0);
		stucReallocMeshToFit(pCtx, &basic.outMesh);
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
Result getOriginIndexedAttrib(
	StucContext pCtx,
	Attrib *pAttrib,
	const StucMapArr *pMapArr,
	I32 mapIdx,
	const AttribIndexed *pMapIndexedAttrib,
	const AttribIndexed *pInIndexedAttrib,
	const AttribIndexed **ppMatsToAdd,
	StucDomain domain
) {
	Result err = STUC_SUCCESS;
	switch (pAttrib->origin) {
		case STUC_ATTRIB_ORIGIN_MAP:
			*ppMatsToAdd = pMapIndexedAttrib;
			break;
		case STUC_ATTRIB_ORIGIN_MESH_IN:
			*ppMatsToAdd = pInIndexedAttrib;
			break;
		case STUC_ATTRIB_ORIGIN_COMMON:
			const CommonAttribList *pMapCommon = pMapArr->pCommonAttribArr + mapIdx;
			const CommonAttrib *pCommonAttrib =
				stucGetCommonAttribFromDomain(pMapCommon, pAttrib->core.name, domain);
			BlendConfig config = {0};
			if (pCommonAttrib) {
				config = pCommonAttrib->blendConfig;
			}
			else {
				StucTypeDefault *pDefaultConfig =
					stucGetTypeDefaultConfig(&pCtx->typeDefaults, STUC_ATTRIB_STRING);
				config = pDefaultConfig->blendConfig;
			}
			*ppMatsToAdd = config.order ? pInIndexedAttrib : pMapIndexedAttrib;
			break;
		default:
			STUC_ASSERT("invalid attrib origin for this function", false);
	}
	return err;
}

static
Result iterFacesAndCorrectIdxAttrib(
	StucContext pCtx,
	Attrib *pAttrib,
	Mesh *pMesh,
	AttribIndexed *pOutIndexedAttrib,
	const AttribIndexed *pOriginIndexedAttrib,
	StucDomain domain
) {
	Result err = STUC_SUCCESS;
	STUC_RETURN_ERR_IFNOT_COND(err, pAttrib->core.type == STUC_ATTRIB_I8, "");

	typedef struct {
		I8 idx;
		bool hasRef;
	} TableEntry;
	TableEntry *pTable =
		pCtx->alloc.pCalloc(pOriginIndexedAttrib->count, sizeof(TableEntry));

	I8 *pIdx = pAttrib->core.pData;
	I32 domainCount = stucDomainCountGetIntern(&pMesh->core, domain);
	for (I32 i = 0; i < domainCount; ++i) {
		I32 idx = pIdx[i];
		STUC_THROW_IFNOT_COND(err, idx >= 0 && idx < pOriginIndexedAttrib->count, "", 0);
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
	STUC_CATCH(0, err, ;);
	pCtx->alloc.pFree(pTable);
	return err;
}

static
Result getIndexedAttribInMaps(
	StucContext pCtx,
	const Mesh *pMesh,
	const StucMapArr *pMapArr,
	const Attrib *pAttrib,
	bool *pSame,
	StucDomain domain,
	const AttribIndexed ***pppOut
) {
	Result err = STUC_SUCCESS;
	const AttribIndexed **ppAttribs = pCtx->alloc.pCalloc(pMapArr->count, sizeof(void *));
	bool found = false;
	*pSame = true;
	StucMap pMapCache = NULL;
	for (I32 i = 0; i < pMapArr->count; ++i) {
		const StucMap pMap = pMapArr->ppArr[i];
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
				STUC_THROW_IFNOT_COND(err, pMapAttrib, "", 0);
				pName = pMapAttrib->core.name;
				break;
			}
			default:
				STUC_ASSERT("invalid attrib origin", false);
		}
		const AttribIndexed *pIndexedAttrib =
			stucGetAttribIndexedInternConst(&pMap->indexedAttribs, pName);
		if (pIndexedAttrib) {
			found = true;
			ppAttribs[i] = pIndexedAttrib;
			if (!pMapCache) {
				pMapCache = pMapArr->ppArr[i];
			}
			else if (*pSame) {
				*pSame = pMapCache == pMapArr->ppArr[i];
			}
		}
	}
	if (found) {
		*pppOut = ppAttribs;
		return err;
	}
	STUC_CATCH(0, err, ;);
	pCtx->alloc.pFree(ppAttribs);
	*pppOut = NULL;
	return err;
}

static
Result correctIdxIndices(
	StucContext pCtx,
	const char *pName,
	Mesh *pMeshArr,
	const StucMapArr *pMapArr,
	const AttribIndexedArr *pInIndexedAttribs,
	const AttribIndexed **ppMapAttribs,
	AttribIndexed *pOutIndexedAttrib,
	StucDomain domain
) {
	Result err = STUC_SUCCESS;
	const AttribIndexed *pInIndexedAttrib =
		stucGetAttribIndexedInternConst(pInIndexedAttribs, pName);
	for (I32 i = 0; i < pMapArr->count; ++i) {
		AttribArray *pAttribArr = stucGetAttribArrFromDomain(&pMeshArr[i].core, domain);
		Attrib *pAttrib = stucGetAttribIntern(pName, pAttribArr, false, NULL, NULL);
		if (!ppMapAttribs[i] || !pAttrib) {
			continue;
		}
		Mesh *pMesh = pMeshArr + i;
		const AttribIndexed *pOriginIndexedAttrib = NULL;
		err = getOriginIndexedAttrib(
			pCtx,
			pAttrib,
			pMapArr,
			i,
			ppMapAttribs[i],
			pInIndexedAttrib,
			&pOriginIndexedAttrib,
			domain
		);
		STUC_RETURN_ERR_IFNOT(err, "");
		STUC_RETURN_ERR_IFNOT_COND(
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
		STUC_RETURN_ERR_IFNOT(err, "");
	}
	return err;
}

static
Result mergeIndexedAttribs(
	StucContext pCtx,
	Mesh *pMeshArr,
	const StucMapArr *pMapArr,
	const AttribIndexedArr *pInIndexedAttribs,
	AttribIndexedArr *pOutIndexedAttribs
) {
	Result err = STUC_SUCCESS;
	const StucAlloc *pAlloc = &pCtx->alloc;
	pOutIndexedAttribs->size = pInIndexedAttribs->count;
	pOutIndexedAttribs->pArr =
		pAlloc->pCalloc(pOutIndexedAttribs->size, sizeof(AttribIndexed));
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
				stucGetAttribIndexed(pAttrib->core.name, pOutIndexedAttribs, &pIndexedAttrib);
				if (pIndexedAttrib) {
					//already added
					continue;
				}
				if (pAttrib->origin == STUC_ATTRIB_ORIGIN_MESH_IN) {
					err = stucAppendAndCopyIndexedAttrib(
						pCtx,
						pAttrib->core.name,
						pOutIndexedAttribs,
						pInIndexedAttribs
					);
					STUC_THROW_IFNOT(err, "", 0);
					continue;
				}
				bool same = false;
				const AttribIndexed **ppMapAttribs = NULL;
				err = getIndexedAttribInMaps(
					pCtx,
					pMesh,
					pMapArr,
					pAttrib,
					&same,
					j,
					&ppMapAttribs
				);
				STUC_THROW_IFNOT_COND(err, ppMapAttribs, "", 0);
				//if (same) {
					//get map to add
					//append and copy into out
					//continue;
				//}
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
							break;
						case STUC_ATTRIB_ORIGIN_COMMON:
							STUC_ASSERT("", pAttrib->origin == STUC_ATTRIB_ORIGIN_COMMON);
							pRefAttrib = stucGetAttribIndexedInternConst(
								pInIndexedAttribs,
								pAttrib->core.name
							);
							break;
						default:
							STUC_ASSERT("invalid attrib origin", false);
					}
					STUC_ASSERT("", pRefAttrib);
					pIndexedAttrib = stucAppendIndexedAttrib(
						pCtx,
						pOutIndexedAttribs,
						pRefAttrib->core.name,
						0, //dont allocate pData
						pRefAttrib->core.type,
						pRefAttrib->core.use
					);
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
				STUC_THROW_IFNOT(err, "", 0);
			}
		}
	}
	STUC_CATCH(0, err, ;)
	return err;
}

typedef struct {
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
Result mapToMeshFromJob(void *pArgsVoid) {
	StucMapToMeshArgs *pArgs = pArgsVoid;
	return stucMapToMesh(
		pArgs->pCtx,
		pArgs->pMapArr,
		pArgs->pMeshIn,
		pArgs->pInIndexedAttribs,
		pArgs->pMeshOut,
		pArgs->pOutIndexedAttribs,
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
		InFaceTable inFaceTable = {0};
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
				pMapArr->pCommonAttribArr + i,
				&inFaceTable,
				1.0f,
				-1.0f
			);
			STUC_THROW_IFNOT(err, "map to mesh usg failed", 1);
			err = stucSampleInAttribsAtUsgOrigins(
				pCtx,
				pMap,
				pMeshIn,
				&squaresOut,
				inFaceTable.pArr
			);
			STUC_THROW_IFNOT(err, "", 1);
			InFaceTableToHashTable(&pCtx->alloc, pMap, squaresOut.faceCount, inFaceTable.pArr);
			//*pMeshOut = squaresOut;
			//return STUC_SUCCESS;
			stucMeshDestroy(pCtx, &squaresOut);
			stucAssignActiveAliases(
				pCtx,
				(Mesh *)pMeshIn,
				0x50, //reassign preserve if present
				STUC_DOMAIN_NONE
			);
		}
		err = mapToMeshInternal(
			pCtx,
			pMap,
			pMeshIn,
			&pOutBufArr[i].core,
			matIdx,
			pMapArr->pCommonAttribArr + i,
			NULL,
			wScale,
			receiveLen
		);
		STUC_THROW_IFNOT(err, "map to mesh failed", 1);
		STUC_CATCH(1, err, ;);
		if (pMap->usgArr.count) {
			pCtx->alloc.pFree(pMap->usgArr.pInFaceTable);
			pMap->usgArr.pInFaceTable = NULL;
			stucLinAllocDestroy(inFaceTable.pAlloc);
			inFaceTable.pAlloc = NULL;
			inFaceTable.pArr = NULL;
		}
		STUC_THROW_IFNOT(err, "", 0);
	}
	pMeshOut->type.type = STUC_OBJECT_DATA_MESH;
	Mesh meshOutWrap = {.core = *pMeshOut};
	err = mergeIndexedAttribs(
		pCtx,
		pOutBufArr,
		pMapArr,
		pInIndexedAttribs,
		pOutIndexedAttribs
	);
	STUC_THROW_IFNOT(err, "", 0);
	err = stucMergeObjArr(pCtx, &meshOutWrap, pMapArr->count, pOutObjWrapArr, false);
	STUC_THROW_IFNOT(err, "", 0);
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
	StucMesh meshIn, //passed by value so we can set active attrib domains if missing
	UBitField32 spAttribsToAppend,
	bool *pBuildEdges
) {
	Result err = STUC_SUCCESS;
	err = attemptToSetMissingActiveDomains(&meshIn);
	STUC_RETURN_ERR_IFNOT(err, "");
	stucAliasMeshCoreNoAttribs(&pWrap->core, &meshIn);
	*pBuildEdges = !meshIn.edgeCount;
	if (*pBuildEdges) {
		printf("no edge list found, building one\n");
		STUC_RETURN_ERR_IFNOT_COND(
			err,
			!meshIn.edgeAttribs.count,
			"in-mesh has edge attribs, yet no edge list"
		);
		err = stucBuildEdgeList(pCtx, pWrap);
		STUC_RETURN_ERR_IFNOT(err, "failed to build edge list");
		printf("finished building edge list\n");
	}
	err = appendSpAttribsToInMesh(pCtx, pWrap, &meshIn, spAttribsToAppend);
	STUC_RETURN_ERR_IFNOT(err, "");
	stucSetAttribOrigins(&pWrap->core.meshAttribs, STUC_ATTRIB_ORIGIN_MESH_IN);
	stucSetAttribOrigins(&pWrap->core.faceAttribs, STUC_ATTRIB_ORIGIN_MESH_IN);
	stucSetAttribOrigins(&pWrap->core.cornerAttribs, STUC_ATTRIB_ORIGIN_MESH_IN);
	stucSetAttribOrigins(&pWrap->core.edgeAttribs, STUC_ATTRIB_ORIGIN_MESH_IN);
	stucSetAttribOrigins(&pWrap->core.vertAttribs, STUC_ATTRIB_ORIGIN_MESH_IN);

	err = stucAssignActiveAliases(pCtx, pWrap, 0x1ef5e, STUC_DOMAIN_NONE);
	STUC_RETURN_ERR_IFNOT(err, "");

	err = stucBuildTangents(pWrap);
	STUC_RETURN_ERR_IFNOT(err, "failed to build tangents");
	buildEdgeCornersTable(pWrap);
	buildSeamAndPreserveTables(&pCtx->alloc, pWrap);

	//set sp
	stucSetAttribCopyOpt(pCtx, &pWrap->core, STUC_ATTRIB_DONT_COPY, spAttribsToAppend);
	//set required
	stucSetAttribCopyOpt(pCtx, &pWrap->core, STUC_ATTRIB_COPY, 0xc0e);

	return err;
}

Result stucMapToMesh(
	StucContext pCtx,
	const StucMapArr *pMapArr,
	const StucMesh *pMeshIn,
	const StucAttribIndexedArr *pInIndexedAttribs,
	StucMesh *pMeshOut,
	StucAttribIndexedArr *pOutIndexedAttribs,
	F32 wScale,
	F32 receiveLen
) {
	StucResult err = STUC_SUCCESS;
	STUC_RETURN_ERR_IFNOT_COND(err, pMeshIn, "");
	err = stucValidateMesh(pMeshIn, false);
	STUC_RETURN_ERR_IFNOT(err, "invalid in-mesh");
	Mesh meshInWrap = {0};
	UBitField32 spAttribsToAppend = 0x1e300;
	bool builtEdges = false;
	err = initMeshInWrap(
		pCtx,
		&meshInWrap,
		*(StucMesh *)pMeshIn,
		spAttribsToAppend,
		&builtEdges
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
		wScale,
		receiveLen
	);
	STUC_THROW_IFNOT(err, "mapMapArrToMesh returned error", 0);
	printf("----------------------FINISHING IN-MESH\n");
	STUC_CATCH(0, err, ;);
	if (builtEdges && meshInWrap.core.pEdges) {
		if (meshInWrap.core.pEdges) {
			pCtx->alloc.pFree(meshInWrap.core.pEdges);
			meshInWrap.core.pEdges = NULL;
		}
	}
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

StucResult stucGetAttribSize(StucAttribCore *pAttrib, I32 *pSize) {
	*pSize = stucGetAttribSizeIntern(pAttrib->type);
	return STUC_SUCCESS;
}

StucResult stucGetAttrib(const char *pName, StucAttribArray *pAttribs, StucAttrib **ppAttrib) {
	*ppAttrib = stucGetAttribIntern(pName, pAttribs, false, NULL, NULL);
	return STUC_SUCCESS;
}

StucResult stucAttribGetAsVoid(StucAttribCore *pAttrib, int32_t idx, void **ppOut) {
	StucResult err = STUC_SUCCESS;
	STUC_RETURN_ERR_IFNOT_COND(err, pAttrib && ppOut && idx >= 0, "");
	*ppOut = stucAttribAsVoid(pAttrib, idx);
	return STUC_SUCCESS;
}

StucResult stucGetAttribIndexed(
	const char *pName,
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
		!v3F32IsFinite(bc)) {
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
	V3_F32 normal = v3F32Cross(ab, ac);
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

Result stucAttribSpTypesGet(StucContext pCtx, const AttribType **ppTypes) {
	Result err = STUC_SUCCESS;
	STUC_RETURN_ERR_IFNOT_COND(err, pCtx && ppTypes, "");
	*ppTypes = pCtx->spAttribTypes;
	return err;
}

Result stucAttribSpDomainsGet(StucContext pCtx, const StucDomain **ppDomains) {
	Result err = STUC_SUCCESS;
	STUC_RETURN_ERR_IFNOT_COND(err, pCtx && ppDomains, "");
	*ppDomains = pCtx->spAttribDomains;
	return err;
}

Result stucAttribSpIsValid(
	StucContext pCtx,
	const AttribCore *pCore,
	StucDomain domain
) {
	Result err = STUC_SUCCESS;
	STUC_RETURN_ERR_IFNOT_COND(err, pCtx && pCore, "");
	return 
		pCtx->spAttribTypes[pCore->use] == pCore->type &&
		pCtx->spAttribDomains[pCore->use] == domain;
}

Result stucAttribGetAllDomains(
	StucContext pCtx,
	StucMesh *pMesh,
	const char *pName,
	Attrib **ppAttrib,
	StucDomain *pDomain
) {
	for (I32 i = STUC_DOMAIN_FACE; i <= STUC_DOMAIN_VERT; ++i) {
		AttribArray *pArr = stucGetAttribArrFromDomain(pMesh, i);
		Attrib *pAttrib = stucGetAttribIntern(pName, pArr, false, NULL, NULL);
		if (pAttrib) {
			*ppAttrib = pAttrib;
			*pDomain = i;
			break;
		}
	}
	return STUC_SUCCESS;
}

StucResult stucAttribGetAllDomainsConst(
	StucContext pCtx,
	const StucMesh *pMesh,
	const char *pName,
	const StucAttrib **ppAttrib,
	StucDomain *pDomain
) {
	return stucAttribGetAllDomains(pCtx, (StucMesh *)pMesh, pName, (StucAttrib **)ppAttrib, pDomain);
}

StucResult stucAttribArrGet(
	StucContext pCtx,
	StucMesh *pMesh,
	StucDomain domain,
	StucAttribArray **ppArr
) {
	Result err = STUC_SUCCESS;
	STUC_RETURN_ERR_IFNOT_COND(err, pCtx && pMesh && ppArr, "");
	*ppArr = stucGetAttribArrFromDomain(pMesh, domain);
	return err;
}

StucResult stucAttribArrGetConst(
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

StucResult stucAttribGetCompType(
	StucContext pCtx,
	StucAttribType type,
	StucAttribType *pCompType
) {
	Result err = STUC_SUCCESS;
	STUC_RETURN_ERR_IFNOT_COND(err, pCtx && pCompType, "");
	*pCompType = stucAttribGetCompTypeIntern(type);
	return err;
}

StucResult stucAttribTypeGetVecSize(
	StucContext pCtx,
	StucAttribType type,
	I32 *pSize
) {
	Result err = STUC_SUCCESS;
	STUC_RETURN_ERR_IFNOT_COND(err, pCtx && pSize, "");
	*pSize = stucAttribTypeGetVecSizeIntern(type);
	return err;
}

StucResult stucDomainCountGet(
	StucContext pCtx,
	const StucMesh *pMesh,
	StucDomain domain,
	int32_t *pCount
) {
	Result err = STUC_SUCCESS;
	STUC_RETURN_ERR_IFNOT_COND(err, pCtx && pMesh && pCount, "");
	*pCount = stucDomainCountGetIntern(pMesh, domain);
	return err;
}

StucResult stucAttribIndexedArrDestroy(StucContext pCtx, StucAttribIndexedArr *pArr) {
	Result err = STUC_SUCCESS;
	STUC_RETURN_ERR_IFNOT_COND(err, pCtx && pArr, "");
	if (pArr->pArr) {
		for (I32 i = 0; i < pArr->count; ++i) {
			if (pArr->pArr[i].core.pData) {
				pCtx->alloc.pFree(pArr->pArr[i].core.pData);
			}
		}
		pCtx->alloc.pFree(pArr->pArr);
		pArr->pArr = NULL;
	}
	pArr->count = pArr->size = 0;
	return err;
}

StucResult stucMapArrDestroy(StucContext pCtx, StucMapArr *pMapArr) {
	Result err = STUC_SUCCESS;
	STUC_RETURN_ERR_IFNOT_COND(err, pCtx && pMapArr, "");
	if (pMapArr->pCommonAttribArr) {
		for (I32 i = 0; i < pMapArr->count; ++i) {
			if (pMapArr->pCommonAttribArr[i].mesh.pArr) {
				pCtx->alloc.pFree(pMapArr->pCommonAttribArr[i].mesh.pArr);
			}
			if (pMapArr->pCommonAttribArr[i].face.pArr) {
				pCtx->alloc.pFree(pMapArr->pCommonAttribArr[i].face.pArr);
			}
			if (pMapArr->pCommonAttribArr[i].corner.pArr) {
				pCtx->alloc.pFree(pMapArr->pCommonAttribArr[i].corner.pArr);
			}
			if (pMapArr->pCommonAttribArr[i].edge.pArr) {
				pCtx->alloc.pFree(pMapArr->pCommonAttribArr[i].edge.pArr);
			}
			if (pMapArr->pCommonAttribArr[i].vert.pArr) {
				pCtx->alloc.pFree(pMapArr->pCommonAttribArr[i].vert.pArr);
			}
		}
		pCtx->alloc.pFree(pMapArr->pCommonAttribArr);
		pMapArr->pCommonAttribArr = NULL;
	}
	if (pMapArr->pMatArr) {
		pCtx->alloc.pFree(pMapArr->pMatArr);
		pMapArr->pMatArr = NULL;
	}
	if (pMapArr->ppArr) {
		pCtx->alloc.pFree(pMapArr->ppArr);
		pMapArr->ppArr = NULL;
	}
	return err;
}

StucResult stucObjectInit(
	StucContext pCtx,
	StucObject *pObj,
	StucMesh *pMesh,
	const Stuc_M4x4_F32 *pTransform
) {
	Result err = STUC_SUCCESS;
	STUC_RETURN_ERR_IFNOT_COND(err, pCtx && pObj, "");
	pObj->pData = (StucObjectData *)pMesh;
	if (pTransform) {
		pObj->transform = *pTransform;
	}
	else {
		pObj->transform = identM4x4;
	}
	return err;
}