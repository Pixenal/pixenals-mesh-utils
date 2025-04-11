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
	pCtx->stageReport.fpBegin = stucStageBegin;
	pCtx->stageReport.fpProgress = stucStageProgress;
	pCtx->stageReport.fpEnd = stucStageEnd;
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
	return STUC_SUCCESS;
}

StucResult stucContextDestroy(StucContext pCtx) {
	pCtx->threadPool.fpDestroy(pCtx->pThreadPoolHandle);
	pCtx->alloc.fpFree(pCtx);
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
		//STUC_ASSERT("more than 2 corners refernce 1 edge", pSet[edge] < 2);
		V3_F32 diff = _(pos V3SUB pPosCache[edge]);
		pMesh->pEdgeLen[edge] = v3F32Len(diff);
		pSet[edge]++;
	}
	pCtx->alloc.fpFree(pSet);
	pCtx->alloc.fpFree(pPosCache);
}

static
void TEMPsetSpFromAttribName(StucContext pCtx, StucMesh *pMesh, AttribArray *pArr) {
	for (I32 j = 0; j < pArr->count; ++j) {
		if (!strcmp(pArr->pArr[j].core.name, "StucMaterialIndices") ||
			!strcmp(pArr->pArr[j].core.name, "materials")
		) {
			strncpy(pArr->pArr[j].core.name, "materials", STUC_ATTRIB_NAME_MAX_LEN);
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

static
void triCacheBuild(const StucAlloc *pAlloc, StucMap pMap) {
	bool ngons = checkForNgonsInMesh(&pMap->pMesh->core);
	if (ngons) {
		pMap->triCache.pArr =
			pAlloc->fpCalloc(pMap->pMesh->core.faceCount, sizeof(FaceTriangulated));
		stucLinAllocInit(pAlloc, &pMap->triCache.alloc, 3, 16, false);
		for (I32 i = 0; i < pMap->pMesh->core.faceCount; ++i) {
			FaceRange face = stucGetFaceRange(&pMap->pMesh->core, i);
			FaceTriangulated *pTris = pMap->triCache.pArr + i;
			stucLinAlloc(&pMap->triCache.alloc, (void**)&pTris->pTris, face.size - 2);
			stucTriangulateFaceFromVerts(pAlloc, &face, pMap->pMesh, pTris);
		}
	}
}

static
void triCacheDestroy(const StucAlloc *pAlloc, StucMap pMap) {
	STUC_ASSERT("", !((pMap->triCache.pArr != NULL) ^ (pMap->triCache.alloc.valid)));
	if (pMap->triCache.pArr) {
		pAlloc->fpFree(pMap->triCache.pArr);
		stucLinAllocDestroy(&pMap->triCache.alloc);
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

StucResult stucMapFileLoad(StucContext pCtx, StucMap *pMapHandle, const char *filePath) {
	StucResult err = STUC_NOT_SET;
	StucMap pMap = pCtx->alloc.fpCalloc(1, sizeof(MapFile));
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

	strncpy(pMap->indexedAttribs.pArr[0].core.name, "materials", STUC_ATTRIB_NAME_MAX_LEN);

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
			STUC_ATTRIB_USE_FIELD(((StucAttribUse[]) {
				STUC_ATTRIB_USE_POS,
				STUC_ATTRIB_USE_UV,
				STUC_ATTRIB_USE_NORMAL,
				STUC_ATTRIB_USE_RECEIVE,
				STUC_ATTRIB_USE_USG,
				STUC_ATTRIB_USE_IDX
			})),
			STUC_DOMAIN_NONE
		);
		STUC_THROW_IFNOT(err, "", 0);
		stucApplyObjTransform(pObjArr + i);
	}
	Mesh *pMapMesh = pCtx->alloc.fpCalloc(1, sizeof(Mesh));
	pMapMesh->core.type.type = STUC_OBJECT_DATA_MESH_INTERN;
	err = stucMergeObjArr(pCtx, pMapMesh, objCount, pObjArr, false);
	STUC_THROW_IFNOT(err, "", 0);

	stucAppendSpAttribsToMesh(
		pCtx,
		pMapMesh,
		0x1 << STUC_ATTRIB_USE_EDGE_LEN,
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

	triCacheBuild(&pCtx->alloc, pMap);
	buildFaceBBoxes(&pCtx->alloc, pMap);

	//the quadtree is created before USGs are assigned to verts,
	//as the tree's used to speed up the process
	printf("File loaded. Creating quad tree\n");
	err = stucCreateQuadTree(pCtx, &pMap->quadTree, pMap->pMesh, pMap->pFaceBBoxes);
	STUC_THROW_IFNOT(err, "failed to create quadtree", 0);

	if (pMap->usgArr.count) {
		pMap->usgArr.pArr = pCtx->alloc.fpCalloc(pMap->usgArr.count, sizeof(Usg));
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
				0x1 << STUC_ATTRIB_USE_POS,
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
					0x1 << STUC_ATTRIB_USE_POS,
					STUC_DOMAIN_NONE
				);
				STUC_THROW_IFNOT(err, "", 0);
				stucApplyObjTransform(pUsgArr[i].pFlatCutoff);
			}
		}
		//TODO remove duplicate uses of alloc where pCtx is present
		//like this
		Mesh *pSquares = pCtx->alloc.fpCalloc(1, sizeof(Mesh));
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
	pCommonArr->pArr = pCtx->alloc.fpCalloc(pCommonArr->size, sizeof(StucCommonAttrib));
	for (I32 i = 0; i < pMeshAttribs->count; ++i) {
		Attrib *pAttrib = pMeshAttribs->pArr + i;
		const Attrib *pMapAttrib = NULL;
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
			pCommonArr->pArr = pCtx->alloc.fpRealloc(pCommonArr->pArr, pCommonArr->size);
		}
		initCommonAttrib(pCtx, pCommonArr->pArr + pCommonArr->count, pAttrib);
		pCommonArr->count++;
	}
	STUC_CATCH(0, err,
		pCtx->alloc.fpFree(pCommonArr->pArr);
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
			pCtx->alloc.fpFree(pArr->pArr);
			pArr->pArr = NULL;
		}
		pArr->count = pArr->size = 0;
	}
	return err;
}

#ifndef TEMP_DISABLE
static
Result sendOffMappingJobs(
	MapToMeshBasic *pBasic,
	I32 *pJobCount,
	void ***pppJobHandles,
	FindEncasedFacesJobArgs **ppJobArgs
) {
	Result err = STUC_SUCCESS;
	*pJobCount = MAX_SUB_MAPPING_JOBS;
	*pJobCount += *pJobCount == 0;
	I32 facesPerThread = pBasic->pInMesh->core.faceCount / *pJobCount;
	bool singleThread = !facesPerThread;
	*pJobCount = singleThread ? 1 : *pJobCount;
	void *jobArgPtrs[MAX_THREADS] = {0};
	I32 borderTableSize = pBasic->pMap->pMesh->core.faceCount / 5 + 2; //+ 2 incase is 0
	*ppJobArgs = pBasic->pCtx->alloc.fpCalloc(*pJobCount, sizeof(FindEncasedFacesJobArgs));
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
	*pppJobHandles = pBasic->pCtx->alloc.fpCalloc(*pJobCount, sizeof(void *));
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
#endif

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
void buildSeamAndPreserveTables(Mesh *pMesh) {
	for (I32 i = 0; i < pMesh->core.edgeCount; ++i) {
		pMesh->pSeamEdge[i] = stucIsEdgeSeam(pMesh, i);
		if (stucGetIfPreserveEdge(pMesh, i)) {
			V2_I32 faces = pMesh->pEdgeFaces[i];
			V2_I8 corners = pMesh->pEdgeCorners[i];
			I32 vert = stucGetMeshVert(
				&pMesh->core,
				(FaceCorner) {.face = faces.d[0], .corner = corners.d[0]}
			);
			if (pMesh->pNumAdjPreserve[vert] < 3) { //only record up to 3
				pMesh->pNumAdjPreserve[vert]++;
			}
			vert = stucGetMeshVert(
				&pMesh->core,
				(FaceCorner) {.face = faces.d[1], .corner = corners.d[1]}
			);
			if (pMesh->pNumAdjPreserve[vert] < 3) {
				pMesh->pNumAdjPreserve[vert]++;
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

#ifndef TEMP_DISABLE
static
void mappingJobArgsDestroy(MapToMeshBasic *pBasic, I32 jobCount, FindEncasedFacesJobArgs *pJobArgs) {
	for (I32 i = 0; i < jobCount; ++i) {
		BufMesh *pBufMesh = &pJobArgs[i].bufMesh;
		stucMeshDestroy(pBasic->pCtx, &pBufMesh->mesh.core);
		if (pJobArgs[i].borderTable.pTable) {
			pBasic->pCtx->alloc.fpFree(pJobArgs[i].borderTable.pTable);
		}
		stucBorderTableDestroyAlloc(&pJobArgs[i].borderTableAlloc);
	}
}
#endif

static
void setJobArgsCore(
	const MapToMeshBasic *pBasic,
	JobArgs *pCore,
	Range *pRanges,
	I32 jobIdx
) {
	pCore->pBasic = pBasic;
	pCore->range = pRanges[jobIdx];
	pCore->id = jobIdx;
}

static
void divideArrAmongstJobs(I32 arrSize, I32 *pJobCount, Range *pRanges) {
	if (!arrSize) {
		*pJobCount = 0;
		return;
	}
	STUC_ASSERT("", *pJobCount >= 0);
	I32 jobCount =  *pJobCount && *pJobCount < MAX_SUB_MAPPING_JOBS ?
		*pJobCount : MAX_SUB_MAPPING_JOBS;
	I32 piecesPerJob = arrSize / jobCount;
	jobCount = !piecesPerJob ? 1 : jobCount;
	for (I32 i = 0; i < jobCount; ++i) {
		pRanges[i].start = piecesPerJob * i;
		pRanges[i].end = i == jobCount - 1 ? arrSize : pRanges[i].start + piecesPerJob;
	}
	*pJobCount = jobCount;
}

static
void makeJobArgs(
	MapToMeshBasic *pBasic,
	I32 *pJobCount, void *pArgs, I32 argStructSize,
	void *pInitInfo,
	I32 (* fpGetArrCount)(const MapToMeshBasic *, void *),
	void (* fpInitArgEntry)(MapToMeshBasic *, void *, void *)
) {
	Range ranges[MAX_SUB_MAPPING_JOBS] = {0};
	divideArrAmongstJobs(fpGetArrCount(pBasic, pInitInfo), pJobCount, ranges);
	for (I32 i = 0; i < *pJobCount; ++i) {
		void *pArgEntry = (U8 *)pArgs + i * argStructSize;
		setJobArgsCore(pBasic, (JobArgs *)pArgEntry, ranges, i);
		if (fpInitArgEntry) {
			fpInitArgEntry(pBasic, pInitInfo, pArgEntry);
		}
	}
}

static
I32 encasedTableJobsGetRange(const MapToMeshBasic *pBasic, void *pInitInfo) {
	return pBasic->pInMesh->core.faceCount;
}

/*
static
EncasedEntryIdxBucket *getEntryIdxBucket(
	EncasedEntryIdxTable *pIdxTable,
	I32 mapFace,
	V2_I16 tile
) {
	U32 hash = stucGetEncasedFaceHash(mapFace, tile, pIdxTable->size);
	return pIdxTable->pTable + hash;
}
*/

typedef struct InPieceInitInfo {
	EncasedMapFace *pMapFace;
	InPieceArr *pInPieceArr;
} InPieceInitInfo;

static
void inPieceInit (
	void *pUserData,
	HTableEntryCore *pIdxEntryCore,
	const void *pKeyData,
	void *pInitInfoVoid,
	I32 linIdx
) {
	EncasedEntryIdx *pIdxEntry = (EncasedEntryIdx *)pIdxEntryCore;
	InPieceArr *pInPieceArr = ((InPieceInitInfo *)pInitInfoVoid)->pInPieceArr;
	const InPieceKey *pKey = pKeyData;
	pIdxEntry->mapFace = pKey->mapFace;
	pIdxEntry->tile = pKey->tile;
	pIdxEntry->entryIdx = pInPieceArr->count;

	EncasedMapFace *pMapFace = ((InPieceInitInfo *)pInitInfoVoid)->pMapFace;
	InPiece *pInPiece = pInPieceArr->pArr + pInPieceArr->count;
	pInPiece->pList = pMapFace;
	pInPiece->faceCount = pMapFace->inFaces.count;
	pInPieceArr->count++;
}

static
bool inPieceCmp(
	const HTableEntryCore *pIdxEntryCore,
	const void *pKeyData,
	const void *pInitInfo
) {
	const EncasedEntryIdx *pIdxEntry = (EncasedEntryIdx *)pIdxEntryCore;
	const InPieceKey *pMapFace = pKeyData;
	return
		pIdxEntry->mapFace == pMapFace->mapFace &&
		_(pIdxEntry->tile V2I16EQL pMapFace->tile);
}

static
void appendEncasedEntryToInPiece(
	EncasedMapFace *pEntry,
	EncasedEntryIdx *pIdxEntry,
	InPieceArr *pInPieceArr
) {
	pInPieceArr->pArr[pIdxEntry->entryIdx].faceCount += pEntry->inFaces.count;
	HTableEntryCore *pInPiece = (HTableEntryCore *)pInPieceArr->pArr[pIdxEntry->entryIdx].pList;
	while (pInPiece->pNext) {
		pInPiece = pInPiece->pNext;
	}
	pInPiece->pNext = (HTableEntryCore *)pEntry;
}

static
void addEncasedEntryToInPieceArr(
	const MapToMeshBasic *pBasic,
	HTable *pIdxTable,
	InPieceArr *pInPieceArr,
	EncasedMapFace *pMapFace
) {
	EncasedEntryIdx *pIdxEntry = NULL;
	SearchResult result = stucHTableGet(
		pIdxTable,
		0,
		&(InPieceKey) {.mapFace = pMapFace->mapFace, .tile = pMapFace->tile},
		(void **)&pIdxEntry,
		true, &(InPieceInitInfo) {.pMapFace = pMapFace, .pInPieceArr = pInPieceArr},
		stucInPieceMakeKey, NULL, inPieceInit, inPieceCmp
	);
	if (result == STUC_SEARCH_FOUND) {
		appendEncasedEntryToInPiece(pMapFace, pIdxEntry, pInPieceArr);
	}
}

static
void linkEncasedTableEntries(
	const MapToMeshBasic *pBasic,
	I32 jobCount,
	FindEncasedFacesJobArgs *pJobArgs,
	InPieceArr *pInPieceArr,
	bool *pEmpty
) {
	const StucAlloc *pAlloc = &pBasic->pCtx->alloc;
	pInPieceArr->size = pInPieceArr->count = 0;
	for (I32 i = 0; i < jobCount; ++i) {
		LinAlloc *pTableAlloc = stucHTableAllocGet(&pJobArgs[i].encasedFaces, 0);
		pInPieceArr->size += stucLinAllocGetCount(pTableAlloc);
	}
	if (pInPieceArr->size == 0) {
		*pEmpty = true;
		return;
	}
	pInPieceArr->pArr = pAlloc->fpCalloc(pInPieceArr->size, sizeof(InPiece));
	HTable idxTable = {0};
	stucHTableInit(
		pAlloc,
		&idxTable,
		pInPieceArr->size / 4 + 1,
		(I32Arr) {.pArr = (I32[]) {sizeof(EncasedEntryIdx)}, .count = 1},
		NULL
	);

	for (I32 i = 0; i < jobCount; ++i) {
		LinAlloc *pTableAlloc = stucHTableAllocGet(&pJobArgs[i].encasedFaces, 0);
		LinAllocIter iter = {0};
		stucLinAllocIterInit(pTableAlloc, (Range) {0, INT32_MAX}, &iter);
		for (; !stucLinAllocIterAtEnd(&iter); stucLinAllocIterInc(&iter)) {
			EncasedMapFace *pEntry = stucLinAllocGetItem(&iter);
			addEncasedEntryToInPieceArr(pBasic, &idxTable, pInPieceArr, pEntry);
			pEntry->core.pNext = NULL;
		}
	}
	stucHTableDestroy(&idxTable);
	*pEmpty = false;
}

I32 inPiecesJobsGetRange(const MapToMeshBasic *pBasic, void *pInitEntry) {
	return ((InPieceArr *)pInitEntry)->count;
}

void splitInPiecesJobInit(MapToMeshBasic *pBasic, void *pInitInfo, void *pEntryVoid) {
	SplitInPiecesJobArgs *pEntry = pEntryVoid;
	pEntry->pInPieceArr = pInitInfo;
}

typedef struct PieceFaceIdx {
	HTableEntryCore core;
	const EncasingInFace *pInFace;
	bool removed;
	bool pendingRemove;
	bool preserve[4];
} PieceFaceIdx;

/*
typedef struct {
	struct PieceFaceIdx *pList;
} PieceFaceIdxBucket;

typedef struct {
	PieceFaceIdxBucket *pTable;
	I32 size;
} PieceFaceIdxTable;
*/

static
void initPieceFaceIdxEntry(
	void *pUserData,
	HTableEntryCore *pEntry,
	const void *pKeyData,
	void *pInitInfo,
	I32 linIdx
) {
	EncasingInFace *pInFace = pInitInfo;
	((PieceFaceIdx *)pEntry)->pInFace = pInFace;
}

static
bool cmpPieceFaceIdxEntry(
	const HTableEntryCore *pEntry,
	const void *pKeyData,
	const void *pInitInfo
) {
	return ((PieceFaceIdx *)pEntry)->pInFace->idx == *(I32 *)pKeyData;
}

/*
static
U32 getPieceFaceIdxTableHash(PieceFaceIdxTable *pTable, I32 face) {
	U32 key = face;
	return stucFnvHash((U8 *)&key, sizeof(key), pTable->size);
}
*/

/*
static
PieceFaceIdx *getPieceFaceIdxEntry(HTable *pTable, I32 face) {
	PieceFaceIdx *pEntry = stucHTableBucketGet(pTable, face)->pList;
	while (pEntry) {
		if (pEntry->pInFace->idx == face) {
			return pEntry;
		}
		pEntry = pEntry->core.pNext;
	}
	return NULL;
}
*/

static
void buildPieceFaceIdxTable(
	HTable *pTable,
	const EncasingInFaceArr *pInFaces
) {
	for (I32 i = 0; i < pInFaces->count; ++i) {
		const EncasingInFace *pInFace = pInFaces->pArr + i;
		PieceFaceIdx *pEntry = NULL;
		I32 faceIdx = pInFace->idx;
		SearchResult result =
			stucHTableGetConst(
				pTable,
				0,
				&faceIdx,
				(void **)&pEntry,
				true, pInFace,
				stucKeyFromI32, NULL, initPieceFaceIdxEntry, cmpPieceFaceIdxEntry
			);
		STUC_ASSERT("", result == STUC_SEARCH_ADDED);
	}
}

static
SearchResult pieceFaceIdxTableGet(HTable *pTable, I32 face, void **ppEntry) {
	return stucHTableGet(
		pTable,
		0,
		&face,
		ppEntry,
		false, NULL,
		stucKeyFromI32, NULL, NULL, cmpPieceFaceIdxEntry
	);
}

static
PieceFaceIdx *getAdjFaceInPiece(
	const MapToMeshBasic *pBasic,
	HTable *pIdxTable,
	FaceCorner corner,
	I32 *pAdjCorner
) {
	FaceCorner adj = {0};
	stucGetAdjCorner(pBasic->pInMesh, corner, &adj);
	if (adj.corner == -1) {
		return NULL;
	}
	STUC_ASSERT("", adj.corner >= 0);
	PieceFaceIdx *pAdjIdxEntry = NULL;
	pieceFaceIdxTableGet(pIdxTable, adj.face, (void **)&pAdjIdxEntry);
	if (!pAdjIdxEntry) {
		if (pAdjCorner) {
			*pAdjCorner = -1;
		}
		return NULL;
	}
	I32 edge = stucGetMeshEdge(&pBasic->pInMesh->core, corner);
	if (pAdjIdxEntry->removed ||
		couldInEdgeIntersectMapFace(pBasic->pInMesh, edge)
	) {
		if (pAdjCorner) {
			*pAdjCorner = -1;
		}
		return NULL;
	}
	if (pAdjCorner) {
		FaceRange adjFaceRange =
			stucGetFaceRange(&pBasic->pInMesh->core, pAdjIdxEntry->pInFace->idx);
		adj.corner = stucGetCornerNext(adj.corner, &adjFaceRange);
		if (pAdjIdxEntry->preserve[adj.corner]) {
			if (pAdjCorner) {
				*pAdjCorner = -1;
			}
			return NULL;
		}
		*pAdjCorner = adj.corner;
	}
	return pAdjIdxEntry;
}

typedef struct BorderEdgeTableEntry {
	HTableEntryCore core;
	FaceCorner corner;
	bool checked;
} BorderEdgeTableEntry;

/*
typedef struct BorderEdgeBucket {
	BorderEdgeTableEntry *pList;
} BorderEdgeBucket;

typedef struct BorderEdgeTable {
	BorderEdgeBucket *pTable;
	LinAlloc alloc;
	I32 size;
} BorderEdgeTable;

typedef struct BorderEdgeArr {
	BorderEdgeTable table;
	FaceCorner *pCorners;
	I32 size;
	I32 count;
} BorderEdgeArr;
*/

static
void borderEdgeInit(
	void *pUserData,
	HTableEntryCore *pEntry,
	const void *pKeyData,
	void *pInitInfo,
	I32 linIdx
) {
	((BorderEdgeTableEntry *)pEntry)->corner = *(FaceCorner *)pKeyData;
}

static
bool borderEdgeCmp(
	const HTableEntryCore *pEntryCore,
	const void *pKeyData,
	const void *pInitInfo
) {
	BorderEdgeTableEntry *pEntry = (BorderEdgeTableEntry *)pEntryCore;
	return
		pEntry->corner.face == ((FaceCorner *)pKeyData)->face &&
		pEntry->corner.corner == ((FaceCorner *)pKeyData)->corner;
}

static
U64 borderEdgeMakeKey(const void *pKeyData) {
	return (U64)*(I32 *)pKeyData;
}

/*
static
BorderEdgeBucket *getBorderEdgeBucket(BorderEdgeTable *pTable, FaceCorner corner) {
	U32 key = corner.face + corner.corner;
	U32 hash = stucFnvHash((U8 *)&key, sizeof(key), pTable->size);
	return pTable->pTable + hash;
}
*/

static
BorderEdgeTableEntry *borderEdgeAddOrGet(
	HTable *pBorderTable,
	FaceCorner corner,
	bool add
) {
	BorderEdgeTableEntry *pEntry = NULL;
	SearchResult result = stucHTableGet(
		pBorderTable,
		0,
		&corner,
		(void **)&pEntry,
		add, &corner,
		borderEdgeMakeKey, NULL, borderEdgeInit, borderEdgeCmp
	);
	STUC_ASSERT(
		"there shouldn't be an existing entry if adding",
		!(add ^ (result == STUC_SEARCH_ADDED))
	);
	return pEntry;
}

static
void addBorderToArr(const MapToMeshBasic *pBasic, BorderArr *pArr, Border border) {
	StucAlloc *pAlloc = &pBasic->pCtx->alloc;
	STUC_ASSERT("", pArr->count <= pArr->size);
	if (!pArr->size) {
		pArr->size = 2;
		pArr->pArr = pBasic->pCtx->alloc.fpMalloc(pArr->size * sizeof(Border));
	}
	else if (pArr->count == pArr->size) {
		pArr->size *= 2;
		pArr->pArr = pAlloc->fpRealloc(pArr->pArr, pArr->size * sizeof(Border));
	}
	pArr->pArr[pArr->count] = border;
	pArr->count++;
}

bool couldInEdgeIntersectMapFace(const Mesh *pInMesh, I32 edge) {
	return
		stucGetIfSeamEdge(pInMesh, edge) ||
		stucGetIfPreserveEdge(pInMesh, edge) ||
		stucGetIfMatBorderEdge(pInMesh, edge);
}

typedef struct BorderBuf {
	BorderArr arr;
	BorderArr preserveRoots;
} BorderBuf;

static
bool findAndAddBorder(
	const MapToMeshBasic *pBasic,
	BorderBuf *pBorderBuf,
	HTable *pIdxTable,
	HTable *pBorderEdgeTable,
	I32 edgesMax,
	BorderEdgeTableEntry *pStart
) {
	const Mesh *pInMesh = pBasic->pInMesh;
	Border border = {.start = pStart->corner};
	FaceCorner corner = pStart->corner;
	BorderEdgeTableEntry *pEntry  = pStart;
	//bool wind = getPieceFaceIdxEntry(pIdxTable, corner.face)->pInFace->wind;
	do {
		if (border.len != 0) {//dont run this on first edge
			if (
				corner.face == pStart->corner.face &&
				corner.corner == pStart->corner.corner
			) {
				break;//full loop
			}
			pEntry = borderEdgeAddOrGet(pBorderEdgeTable, corner, false);
		}
		if (pEntry) {
			STUC_ASSERT("", pEntry->checked == false);
			pEntry->checked = true;
			border.len++;
		}
		STUC_ASSERT("", border.len <= edgesMax);
		I32 adjCorner = 0;
		//this is using the table for the pre-split piece.
		//this is fine, as faces arn't marked removed until the end of this func
		const PieceFaceIdx *pAdjFace = getAdjFaceInPiece(
			pBasic,
			pIdxTable,
			corner,
			&adjCorner
		);
		STUC_ASSERT(
			"if edge isn't in border arr, there should be an adj face",
			!pEntry ^ !pAdjFace
		);
		if (pAdjFace) {
			//edge is internal, move to next adjacent
			corner.face = pAdjFace->pInFace->idx;
			corner.corner = adjCorner;
			//wind = pAdjFace->pInFace->wind;
		}
		else {
			FaceRange faceRange = stucGetFaceRange(&pInMesh->core, corner.face);
			corner.corner = stucGetCornerNext(corner.corner, &faceRange);
		}
	} while(true);
	addBorderToArr(pBasic, &pBorderBuf->arr, border);
	return true;
}

typedef struct InFaceBuf {
	PieceFaceIdx **ppArr;
	I32 size;
	I32 count;
} InFaceBuf;

typedef enum ReceiveStatus {
	STUC_RECEIVE_NONE,
	STUC_RECEIVE_SOME,
	STUC_RECEIVE_ALL
} ReceiveStatus;

typedef struct MapCornerLookup {
	HalfPlane *pHalfPlanes;
	ReceiveStatus receive;
} MapCornerLookup;

typedef enum ReceiveIntersectResult {
	STUC_NO_INTERSECT,
	STUC_INTERSECTS_RECEIVE,
	STUC_INTERSECTS_NON_RECEIVE
} ReceiveIntersectResult;

static
ReceiveIntersectResult doesCornerIntersectReceive(
	const MapToMeshBasic *pBasic,
	const FaceRange *pMapFace, const MapCornerLookup *pMapCorners,
	const FaceRange *pInFace, FaceCorner inCorner
) {
	STUC_ASSERT(
		"check this before calling",
		pMapCorners->receive == STUC_RECEIVE_SOME ||
		pMapCorners->receive == STUC_RECEIVE_ALL
	);
	FaceCorner inCornerNext = {
		.face = pInFace->idx,
		.corner = stucGetCornerNext(inCorner.corner, pInFace)
	};
	V3_F32 inVert =
		pBasic->pInMesh->pPos[stucGetMeshVert(&pBasic->pInMesh->core, inCorner)];
	V3_F32 inVertNext =
		pBasic->pInMesh->pPos[stucGetMeshVert(&pBasic->pInMesh->core, inCornerNext)];
	for (I32 i = 0; i < pMapFace->size; ++i) {
		bool receive = true;
		if (pMapCorners->receive == STUC_RECEIVE_SOME &&
			!stucCheckIfEdgeIsReceive(
				pBasic->pMap->pMesh,
				pMapCorners->pHalfPlanes[i].edge,
				pBasic->receiveLen
		)) {
			receive = false;
		}
		F32 tMapEdge = 0;
		F32 tInEdge = 0;
		stucCalcIntersection(
			inVert, inVertNext,
			pMapCorners->pHalfPlanes[i].uv, pMapCorners->pHalfPlanes[i].dir,
			NULL,
			&tInEdge, &tMapEdge
		);
		if (tInEdge >= .0f && tInEdge <= 1.0f &&
			tMapEdge >= .0f && tMapEdge <= 1.0f
		) {
			return receive ? STUC_INTERSECTS_RECEIVE : STUC_INTERSECTS_NON_RECEIVE;
		}
	}
	return STUC_NO_INTERSECT;
}

static
bool isEdgeValidPreserve(
	const MapToMeshBasic *pBasic,
	const FaceRange *pMapFace,
	const MapCornerLookup *pMapCorners,
	const FaceRange *pInFace,
	FaceCorner inCorner
) {
	I32 edge = stucGetMeshEdge(&pBasic->pInMesh->core, inCorner);
	if (pMapCorners->receive != STUC_RECEIVE_NONE &&
		stucGetIfPreserveEdge(pBasic->pInMesh, edge)
	) {
		ReceiveIntersectResult result = doesCornerIntersectReceive(
			pBasic,
			pMapFace, pMapCorners,
			pInFace, inCorner
		);
		if (result == STUC_NO_INTERSECT || result == STUC_INTERSECTS_RECEIVE) {
			return true;
		}
		STUC_ASSERT("", result == STUC_INTERSECTS_NON_RECEIVE);
	}
	return false;
}

static
void addAdjFaces(
	const MapToMeshBasic *pBasic,
	const FaceRange *pMapFace,
	const MapCornerLookup *pMapCorners,
	InFaceBuf *pInFaceBuf,
	HTable *pIdxTable,
	HTable *pBorderEdges,
	PieceFaceIdx *pFace
) {
	FaceRange inFace = stucGetFaceRange(&pBasic->pInMesh->core, pFace->pInFace->idx);
	for (I32 i = 0; i < inFace.size; ++i) {
		FaceCorner corner = {.face = inFace.idx, .corner = i};
		I32 adjCorner = -1;
		PieceFaceIdx *pAdjFace = getAdjFaceInPiece(
			pBasic,
			pIdxTable,
			corner,
			&adjCorner
		);
		if (!pAdjFace) {
			borderEdgeAddOrGet(pBorderEdges, corner, true);
			continue;
		}
		else if (pAdjFace->pendingRemove) {
			//already added to this piece
			continue; 
		}
		else if (isEdgeValidPreserve(pBasic, pMapFace, pMapCorners, &inFace, corner)) {
			pFace->preserve[i] = true;
			STUC_ASSERT("", adjCorner != -1);
			pAdjFace->preserve[adjCorner] = true;
			borderEdgeAddOrGet(pBorderEdges, corner, true);
			continue;
		}
		STUC_ASSERT("", pAdjFace->pInFace);

		pAdjFace->pendingRemove = true;
		STUC_ASSERT("", pInFaceBuf->count < pInFaceBuf->size);
		pInFaceBuf->ppArr[pInFaceBuf->count] = pAdjFace;
		pInFaceBuf->count++;
	}
}

static
PieceFaceIdx *getFirstRemainingFace(HTable *pIdxTable) {
	LinAlloc *pTableAlloc = stucHTableAllocGet(pIdxTable, 0);
	LinAllocIter iter = {0};
	stucLinAllocIterInit(pTableAlloc, (Range) { 0, INT32_MAX }, &iter);
	for (; !stucLinAllocIterAtEnd(&iter); stucLinAllocIterInc(&iter)) {
		PieceFaceIdx *pEntry = stucLinAllocGetItem(&iter);
		STUC_ASSERT("", pEntry);
		if (!pEntry->removed) {
			return pEntry;
		}
	}
	STUC_ASSERT("this func shouldn't have been called if no faces remained", false);
	return NULL;
}

static
void fillBorderBuf(
	const MapToMeshBasic *pBasic,
	BorderBuf *pBorderBuf,
	HTable *pIdxTable,
	HTable *pBorderEdges
) {
	LinAlloc *pAlloc = stucHTableAllocGet(pBorderEdges, 0);
	pBorderBuf->arr.count = 0;
	pBorderBuf->preserveRoots.count = 0;
	I32 edgeCount = stucLinAllocGetCount(pAlloc);
	LinAllocIter iter = {0};
	stucLinAllocIterInit(pAlloc, (Range) { 0, INT32_MAX }, &iter);
	for (; !stucLinAllocIterAtEnd(&iter); stucLinAllocIterInc(&iter)) {
		BorderEdgeTableEntry *pEntry = stucLinAllocGetItem(&iter);
		if (pEntry->checked) {
			continue;
		}
		I32 edge = stucGetMeshEdge(&pBasic->pInMesh->core, pEntry->corner);
		if (couldInEdgeIntersectMapFace(pBasic->pInMesh, edge)) {
			findAndAddBorder(
				pBasic,
				pBorderBuf,
				pIdxTable,
				pBorderEdges,
				edgeCount,
				pEntry
			);
		}
	}
}

static
void splitAdjFacesIntoPiece(
	SplitInPiecesJobArgs *pArgs,
	const FaceRange *pMapFace,
	MapCornerLookup *pMapCorners,
	InFaceBuf *pInFaceBuf,
	BorderBuf *pBorderBuf,
	const InPiece *pInPiece,
	HTable *pIdxTable,
	InPiece *pNewInPiece,
	I32 *pFacesRemaining
) {
	const MapToMeshBasic *pBasic = pArgs->core.pBasic;
	HTable borderEdges = {0};
	stucHTableInit(
		&pBasic->pCtx->alloc,
		&borderEdges,
		*pFacesRemaining / 2 + 1,
		(I32Arr) {.pArr = (I32[]) {sizeof(BorderEdgeTableEntry)}, .count = 1},
		NULL
	);
	pInFaceBuf->count = 0;
	{
		PieceFaceIdx *pStartFace = getFirstRemainingFace(pIdxTable);
		pInFaceBuf->ppArr[0] = pStartFace;
		pStartFace->pendingRemove = true;
		pInFaceBuf->count++;
		I32 i = 0;
		do {
			PieceFaceIdx *pIdxEntry = NULL;
			pieceFaceIdxTableGet(
				pIdxTable,
				pInFaceBuf->ppArr[i]->pInFace->idx,
				(void **)&pIdxEntry
			);
			addAdjFaces(
				pBasic,
				pMapFace, pMapCorners,
				pInFaceBuf,
				pIdxTable,
				&borderEdges,
				pIdxEntry
			);
		} while (i++, i < pInFaceBuf->count);
	}
	pNewInPiece->faceCount = pInFaceBuf->count;
	fillBorderBuf(pBasic, pBorderBuf, pIdxTable, &borderEdges);
	if (pBorderBuf->arr.count) {
		pNewInPiece->borderArr.count = pBorderBuf->arr.count;
		stucLinAlloc(
			&pArgs->alloc.border,
			(void **)&pNewInPiece->borderArr.pArr,
			pNewInPiece->borderArr.count
		);
		memcpy(
			pNewInPiece->borderArr.pArr,
			pBorderBuf->arr.pArr,
			pBorderBuf->arr.count * sizeof(Border)
		);
	}
	stucHTableDestroy(&borderEdges);
	
	// copy buf into new in-piece, & mark in-faces as removed in idx-table
	pNewInPiece->pList->inFaces.count = pInFaceBuf->count;
	stucLinAlloc(
		&pArgs->alloc.inFace,
		(void **)&pNewInPiece->pList->inFaces.pArr,
		pNewInPiece->pList->inFaces.count
	);
	for (I32 i = 0; i < pInFaceBuf->count; ++i) {
		STUC_ASSERT("", pInFaceBuf->ppArr[i]->pendingRemove);
		pNewInPiece->pList->inFaces.pArr[i] = *pInFaceBuf->ppArr[i]->pInFace;
		pInFaceBuf->ppArr[i]->removed = true;
		pInFaceBuf->ppArr[i]->pendingRemove = false;
	}
	*pFacesRemaining -= pInFaceBuf->count;
}

static
ReceiveStatus getMapFaceReceiveStatus(
	const MapToMeshBasic *pBasic,
	const FaceRange *pFace
) {
	I32 count = 0;
	for (I32 i = 0; i < pFace->size; ++i) {
		I32 edge = stucGetMeshEdge(
			&pBasic->pMap->pMesh->core,
			(FaceCorner) {.face = pFace->idx, .corner = i}
		);
		if (stucCheckIfEdgeIsReceive(pBasic->pMap->pMesh, edge, pBasic->receiveLen)) {
			count++;
		}
	}
	if (!count) {
		return STUC_RECEIVE_NONE;
	}
	else if (count != pFace->size) {
		return STUC_RECEIVE_SOME;
	}
	return STUC_RECEIVE_ALL;
}

static
void splitInPieceEntry(
	SplitInPiecesJobArgs *pArgs,
	const InPiece *pInPiece,
	InFaceBuf *pInFaceBuf,
	BorderBuf *pBorderBuf
) {
	const MapToMeshBasic *pBasic = pArgs->core.pBasic;
	const StucAlloc *pAlloc = &pBasic->pCtx->alloc;

	HTable idxTable = {0};
	stucHTableInit(
		pAlloc,
		&idxTable,
		pInPiece->faceCount / 4 + 1,
		(I32Arr) {.pArr = (I32[]) {sizeof(PieceFaceIdx)}, .count = 1},
		NULL
	);
	{
		EncasedMapFace *pInFaces = pInPiece->pList;
		do {
			buildPieceFaceIdxTable(&idxTable, &pInFaces->inFaces);
			pInFaces = (EncasedMapFace *)pInFaces->core.pNext;
		} while (pInFaces);
	}

	STUC_ASSERT("", pInPiece->faceCount > 0);
	if (!pInFaceBuf->size) {
		pInFaceBuf->size = pInPiece->faceCount;
		pInFaceBuf->ppArr = pAlloc->fpMalloc(pInFaceBuf->size * sizeof(void *));
	}
	else if (pInFaceBuf->size < pInPiece->faceCount) {
		pInFaceBuf->size = pInPiece->faceCount;
		pInFaceBuf->ppArr =
			pAlloc->fpRealloc(pInFaceBuf->ppArr, pInFaceBuf->size * sizeof(void *));
	}

	FaceRange mapFace =
		stucGetFaceRange(&pBasic->pMap->pMesh->core, pInPiece->pList->mapFace);
	MapCornerLookup mapCorners = {
		.pHalfPlanes = pAlloc->fpCalloc(mapFace.size, sizeof(HalfPlane)),
		.receive = getMapFaceReceiveStatus(pBasic, &mapFace)
	};
	initHalfPlaneLookup(pBasic->pMap->pMesh, &mapFace, mapCorners.pHalfPlanes);

	I32 facesRemaining = pInPiece->faceCount;
	do {
		InPiece newInPiece = {0};
		stucLinAlloc(&pArgs->alloc.encased, (void **)&newInPiece.pList, 1);
		newInPiece.pList->mapFace = pInPiece->pList->mapFace;
		newInPiece.pList->tile = pInPiece->pList->tile;
		splitAdjFacesIntoPiece(
			pArgs,
			&mapFace,
			&mapCorners,
			pInFaceBuf,
			pBorderBuf,
			pInPiece,
			&idxTable,
			&newInPiece,
			&facesRemaining
		);
		InPieceArr *pNewInPieces = newInPiece.borderArr.count ?
			&pArgs->newInPiecesClip : &pArgs->newInPieces;
		STUC_ASSERT("", pNewInPieces->count <= pNewInPieces->size);
		if (pNewInPieces->count == pNewInPieces->size) {
			pNewInPieces->size *= 2;
			pNewInPieces->pArr = pAlloc->fpRealloc(
				pNewInPieces->pArr,
				pNewInPieces->size * sizeof(InPiece)
			);
		}
		pNewInPieces->pArr[pNewInPieces->count] = newInPiece;
		pNewInPieces->count++;
		STUC_ASSERT("", facesRemaining >= 0 && facesRemaining < pInPiece->faceCount);
	} while(facesRemaining);
	pAlloc->fpFree(mapCorners.pHalfPlanes);
	stucHTableDestroy(&idxTable);
}

static
Result splitInPieces(void *pArgsVoid) {
	Result err = STUC_SUCCESS;
	SplitInPiecesJobArgs *pArgs = pArgsVoid;
	const StucAlloc *pAlloc = &pArgs->core.pBasic->pCtx->alloc;
	I32 rangeSize = pArgs->core.range.end - pArgs->core.range.start;
	pArgs->newInPieces.size = rangeSize;
	pArgs->newInPiecesClip.size = rangeSize;
	pArgs->newInPieces.pArr = pAlloc->fpMalloc(pArgs->newInPieces.size * sizeof(InPiece));
	pArgs->newInPiecesClip.pArr =
		pAlloc->fpMalloc(pArgs->newInPiecesClip.size * sizeof(InPiece));
	stucLinAllocInit(pAlloc, &pArgs->alloc.encased, sizeof(EncasedMapFace), rangeSize, true);
	stucLinAllocInit(pAlloc, &pArgs->alloc.inFace, sizeof(EncasingInFace), rangeSize, true);
	stucLinAllocInit(pAlloc, &pArgs->alloc.border, sizeof(Border), rangeSize, true);
	InFaceBuf inFaceBuf = {0};
	BorderBuf borderBuf = {0};
	for (I32 i = pArgs->core.range.start; i < pArgs->core.range.end; ++i) {
		splitInPieceEntry(pArgs, pArgs->pInPieceArr->pArr + i, &inFaceBuf, &borderBuf);
	}
	if (inFaceBuf.ppArr) {
		pAlloc->fpFree(inFaceBuf.ppArr);
	}
	if (borderBuf.arr.pArr) {
		pAlloc->fpFree(borderBuf.arr.pArr);
	}
	if (borderBuf.preserveRoots.pArr) {
		pAlloc->fpFree(borderBuf.preserveRoots.pArr);
	}
	return err;
}

static
void appendNewPiecesToArr(
	const MapToMeshBasic *pBasic,
	InPieceArr *pInPiecesSplit,
	I32 jobCount,
	const SplitInPiecesJobArgs *pJobArgs,
	const InPieceArr *(* getNewInPieceArr) (const SplitInPiecesJobArgs *)
) {
	const StucAlloc *pAlloc = &pBasic->pCtx->alloc;
	pInPiecesSplit->count = pInPiecesSplit->size = 0;
	for (I32 i = 0; i < jobCount; ++i) {
		const InPieceArr *pNewInPieces = getNewInPieceArr(pJobArgs + i);
		pInPiecesSplit->size += pNewInPieces->count;
	}
	STUC_ASSERT("", pInPiecesSplit->size >= 0);
	pInPiecesSplit->pArr = pAlloc->fpCalloc(pInPiecesSplit->size, sizeof(InPiece));
	for (I32 i = 0; i < jobCount; ++i) {
		const InPieceArr *pNewInPieces = getNewInPieceArr(pJobArgs + i);
		memcpy(
			pInPiecesSplit->pArr + pInPiecesSplit->count,
			pNewInPieces->pArr,
			pNewInPieces->count * sizeof(InPiece)
		);
		pInPiecesSplit->count += pNewInPieces->count;
	}
}

static
void bufMeshArrDestroy(StucContext pCtx, BufMeshArr *pArr) {
	for (I32 i = 0; i < pArr->count; ++i) {
		if (pArr->arr[i].faces.pArr) {
			pCtx->alloc.fpFree(pArr->arr[i].faces.pArr);
		}
		if (pArr->arr[i].corners.pArr) {
			pCtx->alloc.fpFree(pArr->arr[i].corners.pArr);
		}
		if (pArr->arr[i].inOrMapVerts.pArr) {
			pCtx->alloc.fpFree(pArr->arr[i].inOrMapVerts.pArr);
		}
		if (pArr->arr[i].onEdgeVerts.pArr) {
			pCtx->alloc.fpFree(pArr->arr[i].onEdgeVerts.pArr);
		}
		if (pArr->arr[i].overlapVerts.pArr) {
			pCtx->alloc.fpFree(pArr->arr[i].overlapVerts.pArr);
		}
		if (pArr->arr[i].intersectVerts.pArr) {
			pCtx->alloc.fpFree(pArr->arr[i].intersectVerts.pArr);
		}
		pArr->arr[i] = (BufMesh) {0};
	}
}

void inPieceArrDestroy(const StucContext pCtx, InPieceArr *pArr) {
	if (pArr->pArr) {
		pCtx->alloc.fpFree(pArr->pArr);
	}
	*pArr = (InPieceArr) {0};
	//bufmeshes are stored on stack, so we don't free that
}

typedef struct BufMeshJobInitInfo {
	InPieceArr *pInPiecesSplit;
	Result (* fpAddPiece)(
		const MapToMeshBasic *,
		I32,
		const InPiece *,
		BufMesh *
	);
} BufMeshJobInitInfo;

static
I32 bufMeshInitJobsGetRange(const MapToMeshBasic *pBasic, void *pInitInfoVoid) {
	return ((BufMeshJobInitInfo *)pInitInfoVoid)->pInPiecesSplit->count;
}

static
void bufMeshInitJobInit(MapToMeshBasic *pBasic, void *pInitInfoVoid, void *pEntryVoid) {
	BufMeshInitJobArgs *pEntry = pEntryVoid;
	BufMeshJobInitInfo *pInitInfo = pInitInfoVoid;
	pEntry->pInPiecesSplit = pInitInfo->pInPiecesSplit;
	pEntry->fpAddPiece = pInitInfo->fpAddPiece;
}

static
I32 bufMeshArrGetVertCount(const BufMeshArr *pBufMeshes) {
	I32 total = 0;
	for (I32 i = 0; i < pBufMeshes->count; ++i) {
		total += pBufMeshes->arr[i].inOrMapVerts.count;
		total += pBufMeshes->arr[i].onEdgeVerts.count;
		total += pBufMeshes->arr[i].overlapVerts.count;
		total += pBufMeshes->arr[i].intersectVerts.count;
	}
	return total;
}

static
void vertMergeTableInit(
	const MapToMeshBasic *pBasic,
	const InPieceArr *pInPieces,
	const InPieceArr *pInPiecesClip,
	HTable *pTable
) {
	I32 vertTotal = bufMeshArrGetVertCount(pInPieces->pBufMeshes);
	vertTotal += bufMeshArrGetVertCount(pInPiecesClip->pBufMeshes);
	STUC_ASSERT(
		"mapToMesh should have returned before this func if empty",
		vertTotal > 0
	);
	stucHTableInit(
		&pBasic->pCtx->alloc,
		pTable,
		vertTotal / 4 + 1,
		(I32Arr){
			.pArr = (I32[]) {sizeof(VertMerge), sizeof(VertMergeIntersect)},
			.count = 2
		},
		NULL
	);
}

static
U64 mergeTableMakeKey(const void *pKeyData) {
	const MergeTableKey *pKey = pKeyData;
	U64 key = 0;
	switch (pKey->type) {
		case STUC_BUF_VERT_SUB_TYPE_IN: {
			const InVertKey *pInKey = &pKey->key.inOrMap.in;
			key = (U64)pInKey->inVert << 32 | (U64)pInKey->mapFace;
			break;
		}
		case STUC_BUF_VERT_SUB_TYPE_MAP: {
			const MapVertKey *pMapKey = &pKey->key.inOrMap.map;
			key = (U64)pMapKey->mapVert << 32 | (U64)pMapKey->inFace;
			break;
		}
		case STUC_BUF_VERT_SUB_TYPE_EDGE_IN: {
			const EdgeInVertKey *pInKey = &pKey->key.onEdge.in;
			key = (U64)pInKey->inVert << 32 | (U64)pInKey->mapEdge;
			break;
		}
		case STUC_BUF_VERT_SUB_TYPE_EDGE_MAP: {
			const EdgeMapVertKey *pMapKey = &pKey->key.onEdge.map;
			key = (U64)pMapKey->mapVert << 32 | (U64)pMapKey->inEdge;
			break;
		}
		case STUC_BUF_VERT_OVERLAP:
			key = (U64)pKey->key.overlap.inVert << 32 | (U64)pKey->key.overlap.mapVert;
			break;
		case STUC_BUF_VERT_INTERSECT:
			key = (U64)pKey->key.intersect.inEdge << 32 | (U64)pKey->key.intersect.mapEdge;
			break;
		default:
			STUC_ASSERT("invalid vert type", false);
	}
	return key + ((U64)pKey->tile.d[0] << 16 | (U64)pKey->tile.d[1]);
}

static
bool mergeTableEntryCmp(
	const HTableEntryCore *pEntryCore,
	const void *pKeyData,
	const void *pInitInfo
) {
	const VertMerge *pEntry = (VertMerge *)pEntryCore;
	const MergeTableKey *pKey = pKeyData;
	if (pKey->type != pEntry->key.type ||
		!_(pKey->tile V2I16EQL pEntry->key.tile)
	) {
		return false;
	}
	switch (pKey->type) {
		case STUC_BUF_VERT_SUB_TYPE_IN:
			return
				pKey->key.inOrMap.in.inVert == pEntry->key.key.inOrMap.in.inVert &&
				pKey->key.inOrMap.in.mapFace == pEntry->key.key.inOrMap.in.mapFace;
		case STUC_BUF_VERT_SUB_TYPE_MAP:
			return
				pKey->key.inOrMap.map.mapVert == pEntry->key.key.inOrMap.map.mapVert &&
				pKey->key.inOrMap.map.inFace == pEntry->key.key.inOrMap.map.inFace;
		case STUC_BUF_VERT_SUB_TYPE_EDGE_IN:
			return
				pKey->key.onEdge.in.inVert == pEntry->key.key.onEdge.in.inVert &&
				pKey->key.onEdge.in.mapEdge == pEntry->key.key.onEdge.in.mapEdge;
		case STUC_BUF_VERT_SUB_TYPE_EDGE_MAP:
			return
				pKey->key.onEdge.map.mapVert == pEntry->key.key.onEdge.map.mapVert &&
				pKey->key.onEdge.map.inEdge == pEntry->key.key.onEdge.map.inEdge;
		case STUC_BUF_VERT_OVERLAP:
			return
				pKey->key.overlap.inVert == pEntry->key.key.overlap.inVert &&
				pKey->key.overlap.mapVert == pEntry->key.key.overlap.mapVert;
		case STUC_BUF_VERT_INTERSECT:
			return
				pKey->key.intersect.inEdge == pEntry->key.key.intersect.inEdge &&
				pKey->key.intersect.mapEdge == pEntry->key.key.intersect.mapEdge;
		default:
			STUC_ASSERT("invalid vert type", false);
			return false;
	}
}

static
void mergeTableEntryInit(
	void *pUserData,
	HTableEntryCore *pEntryCore,
	const void *pKeyData,
	void *pInitInfo,
	I32 linIdx
) {
	VertMerge *pEntry = (VertMerge *)pEntryCore;
	const MergeTableKey *pKey = pKeyData;
	VertMergeCorner *pBufCorner = pInitInfo;
	pEntry->key = *pKey;
	pEntry->bufCorner = *pBufCorner;
	pEntry->linIdx = linIdx;
}

static
void initInOrMapKey(
	const MapToMeshBasic *pBasic,
	I32 mapFace,
	const MergeTableInitInfoVert *pInitInfoVert,
	MergeTableKey *pKey
) {
	switch (pInitInfoVert->inOrMap.in.type) {
		case STUC_BUF_VERT_SUB_TYPE_IN: {
			pKey->type = STUC_BUF_VERT_SUB_TYPE_IN;
			const InVert *pVert = &pInitInfoVert->inOrMap.in;
			FaceCorner inCorner = {.face = pVert->inFace, .corner = pVert->inCorner};
			pKey->key.inOrMap.in.inVert =
				stucGetMeshVert(&pBasic->pInMesh->core, inCorner);
			pKey->key.inOrMap.in.mapFace = mapFace;
			break;
		}
		case STUC_BUF_VERT_SUB_TYPE_MAP: {
			pKey->type = STUC_BUF_VERT_SUB_TYPE_MAP;
			const MapVert *pVert = &pInitInfoVert->inOrMap.map;
			FaceCorner mapCorner = {.face = mapFace, .corner = pVert->mapCorner};
			pKey->key.inOrMap.map.mapVert =
				stucGetMeshVert(&pBasic->pMap->pMesh->core, mapCorner);
			pKey->key.inOrMap.map.inFace = pVert->inFace;
			break;
		}
		default:
			STUC_ASSERT("invalid in-or-map buf vert type", false);
	}
}

static
void initOnEdgeKey(
	const MapToMeshBasic *pBasic,
	I32 mapFace,
	const MergeTableInitInfoVert *pInitInfoVert,
	MergeTableKey *pKey
) {
	switch (pInitInfoVert->onEdge.in.type) {
		case STUC_BUF_VERT_SUB_TYPE_EDGE_IN: {
			pKey->type = STUC_BUF_VERT_SUB_TYPE_EDGE_IN;
			const EdgeInVert *pVert = &pInitInfoVert->onEdge.in;
			FaceCorner inCorner = {.face = pVert->inFace, .corner = pVert->inCorner};
			pKey->key.onEdge.in.inVert =
				stucGetMeshVert(&pBasic->pInMesh->core, inCorner);
			FaceCorner mapCorner = {.face = mapFace, .corner = pVert->mapCorner};
			pKey->key.onEdge.in.mapEdge =
				(U64)stucGetMeshEdge(&pBasic->pMap->pMesh->core, mapCorner);
			break;
		}
		case STUC_BUF_VERT_SUB_TYPE_EDGE_MAP: {
			pKey->type = STUC_BUF_VERT_SUB_TYPE_EDGE_MAP;
			const EdgeMapVert *pVert = &pInitInfoVert->onEdge.map;
			FaceCorner mapCorner = {.face = mapFace, .corner = pVert->mapCorner};
			pKey->key.onEdge.map.mapVert =
				stucGetMeshVert(&pBasic->pMap->pMesh->core, mapCorner);
			FaceCorner inCorner = {.face = pVert->inFace, .corner = pVert->inCorner};
			pKey->key.onEdge.map.inEdge =
				(U64)stucGetMeshEdge(&pBasic->pInMesh->core, inCorner);
			break;
		}
		default:
			STUC_ASSERT("invalid edge buf vert type", false);
	}
}

static
void mergeTableInitKey(
	const MapToMeshBasic *pBasic,
	const InPiece *pInPiece,
	BufVertType type,
	const MergeTableInitInfoVert *pInitInfoVert,
	MergeTableKey *pKey
) {
	I32 mapFace = pInPiece->pList->mapFace;
	*pKey = (MergeTableKey) {
		.type = type,
		.tile = pInPiece->pList->tile
	};
	switch (type) {
		case STUC_BUF_VERT_IN_OR_MAP:
			initInOrMapKey(pBasic, mapFace, pInitInfoVert, pKey);
			break;
		case STUC_BUF_VERT_ON_EDGE:
			initOnEdgeKey(pBasic, mapFace, pInitInfoVert, pKey);
			break;
		case STUC_BUF_VERT_OVERLAP: {
			const OverlapVert *pVert = &pInitInfoVert->overlap;
			FaceCorner inCorner = {.face = pVert->inFace, .corner = pVert->inCorner};
			FaceCorner mapCorner = {.face = mapFace, .corner = pVert->mapCorner};
			pKey->key.overlap.inVert = stucGetMeshVert(&pBasic->pInMesh->core, inCorner);
			pKey->key.overlap.mapVert = stucGetMeshVert(&pBasic->pMap->pMesh->core, mapCorner);
			break;
		}
		case STUC_BUF_VERT_INTERSECT: {
			const IntersectVert *pVert = &pInitInfoVert->intersect;
			FaceCorner inCorner = {.face = pVert->inFace, .corner = pVert->inCorner};
			FaceCorner mapCorner = {.face = mapFace, .corner = pVert->mapCorner};
			pKey->key.intersect.inEdge = stucGetMeshEdge(&pBasic->pInMesh->core, inCorner);
			pKey->key.intersect.mapEdge = stucGetMeshEdge(&pBasic->pMap->pMesh->core, mapCorner);
			break;
		}
		default:
			STUC_ASSERT("invalid vert type", false);
	}
}

static
const InPiece *bufFaceGetInPiece(
	const BufMesh *pBufMesh,
	I32 face,
	const InPieceArr *pInPieces
) {
	I32 inPieceIdx = pBufMesh->faces.pArr[face].inPiece;
	return pInPieces->pArr + inPieceIdx;
}

static
BufVertType bufMeshGetType(const BufMesh *pBufMesh, FaceCorner corner) {
	BufFace bufFace = pBufMesh->faces.pArr[corner.face];
	BufCorner bufCorner = pBufMesh->corners.pArr[bufFace.start + corner.corner];
	return bufCorner.type;
}

static
MergeTableInitInfoVert mergeTableGetBufVert(const BufMesh *pBufMesh, FaceCorner corner) {
	BufFace bufFace = pBufMesh->faces.pArr[corner.face];
	BufCorner bufCorner = pBufMesh->corners.pArr[bufFace.start + corner.corner];
	switch (bufCorner.type) {
		case STUC_BUF_VERT_IN_OR_MAP:
			return (MergeTableInitInfoVert) {
				.inOrMap = pBufMesh->inOrMapVerts.pArr[bufCorner.vert]
			};
		case STUC_BUF_VERT_ON_EDGE:
			return (MergeTableInitInfoVert) {
				.onEdge = pBufMesh->onEdgeVerts.pArr[bufCorner.vert]
			};
		case STUC_BUF_VERT_OVERLAP:
			return (MergeTableInitInfoVert) {
				.overlap = pBufMesh->overlapVerts.pArr[bufCorner.vert]
			};
		case STUC_BUF_VERT_INTERSECT:
			return (MergeTableInitInfoVert) {
				.intersect = pBufMesh->intersectVerts.pArr[bufCorner.vert]
			};
		default:
			STUC_ASSERT("invalid buf vert type", false);
			return (MergeTableInitInfoVert) {0};
	}
}

static
void mergeTableAddVert(
	HTable *pTable,
	const MergeTableKey *pKey,
	VertMergeCorner *pInitInfo
) {
	void *pEntry = NULL;
	SearchResult result = stucHTableGet(
		pTable,
		pKey->type == STUC_BUF_VERT_INTERSECT,
		pKey,
		&pEntry,
		true, pInitInfo,
		mergeTableMakeKey, NULL, mergeTableEntryInit, mergeTableEntryCmp
	);
	STUC_ASSERT("", result == STUC_SEARCH_ADDED || result == STUC_SEARCH_FOUND);
}

static
void mergeTableGetVertKey(
	const MapToMeshBasic *pBasic,
	const InPiece *pInPiece,
	const BufMesh *pBufMesh,
	FaceCorner bufCorner,
	MergeTableKey *pKey
) {
	MergeTableInitInfoVert vertInfo = mergeTableGetBufVert(pBufMesh, bufCorner);
	mergeTableInitKey(
		pBasic,
		pInPiece,
		bufMeshGetType(pBufMesh, bufCorner),
		&vertInfo,
		pKey
	);
}

static
void mergeTableAddVerts(
	const MapToMeshBasic *pBasic,
	const InPieceArr *pInPieces,
	bool clipped,
	I32 bufMeshIdx,
	FaceCorner bufCorner,
	HTable *pTable
) {
	const BufMesh *pBufMesh = pInPieces->pBufMeshes->arr + bufMeshIdx;
	const InPiece *pInPiece = bufFaceGetInPiece(pBufMesh, bufCorner.face, pInPieces);
	VertMergeCorner initInfo = {
		.corner = bufCorner,
		.bufMesh = bufMeshIdx,
		.clipped = clipped
	};
	MergeTableKey key = {0};
	mergeTableGetVertKey(pBasic, pInPiece, pBufMesh, bufCorner, &key);
	mergeTableAddVert(pTable, &key, &initInfo);
}

static
void mergeVerts(
	const MapToMeshBasic *pBasic,
	const InPieceArr *pInPieces,
	bool clipped,
	HTable *pTable
) {
	for (I32 i = 0; i < pInPieces->pBufMeshes->count; ++i) {
		const BufMesh *pBufMesh = pInPieces->pBufMeshes->arr + i;
		for (I32 j = 0; j < pBufMesh->faces.count; ++j) {
			for (I32 k = 0; k < pBufMesh->faces.pArr[j].size; ++k) {
				mergeTableAddVerts(
					pBasic,
					pInPieces,
					clipped,
					i,
					(FaceCorner) {.face = j, .corner = k},
					pTable
				);
			}
		}
	}
}

static
const InPieceArr *getNewInPieces(const SplitInPiecesJobArgs *pJobArgs) {
	return &pJobArgs->newInPieces;
}

static
const InPieceArr *getNewInPiecesClip(const SplitInPiecesJobArgs *pJobArgs) {
	return &pJobArgs->newInPiecesClip;
}

static
void bufMeshArrMoveToInPieces(
	const InPieceArr *pInPieces,
	const BufMeshInitJobArgs *pJobArgs,
	I32 jobCount
) {
	BufMeshArr *pBufMeshes = pInPieces->pBufMeshes;
	pBufMeshes->count = jobCount;
	for (I32 i = 0; i < jobCount; ++i) {
		pBufMeshes->arr[i] = pJobArgs[i].bufMesh;
	}
}

static
Result inPieceArrInitBufMeshes(
	MapToMeshBasic *pBasic,
	InPieceArr *pInPieces,
	Result (* fpAddPiece)(
		const MapToMeshBasic *,
		I32,
		const InPiece *,
		BufMesh *
	)
) {
	Result err = STUC_SUCCESS;
	I32 jobCount = 0;
	BufMeshInitJobArgs jobArgs[MAX_SUB_MAPPING_JOBS] = {0};
	makeJobArgs(
		pBasic,
		&jobCount, jobArgs, sizeof(BufMeshInitJobArgs),
		&(BufMeshJobInitInfo) {.pInPiecesSplit = pInPieces, .fpAddPiece = fpAddPiece},
		bufMeshInitJobsGetRange, bufMeshInitJobInit);
	err = stucDoJobInParallel(
		pBasic,
		jobCount, jobArgs, sizeof(BufMeshInitJobArgs),
		stucBufMeshInit
	);
	STUC_RETURN_ERR_IFNOT(err, "");
	bufMeshArrMoveToInPieces(pInPieces, jobArgs, jobCount);
	return err;
}

typedef struct SplitInPiecesAllocArr {
	SplitInPiecesAlloc *pArr;
	I32 count;
} SplitInPiecesAllocArr;

//destroys in-piece arr after splitting
static
Result inPieceArrSplit(
	MapToMeshBasic *pBasic,
	InPieceArr *pInPieces,
	InPieceArr *pInPiecesSplit,
	InPieceArr *pInPiecesSplitClip,
	SplitInPiecesAllocArr *pSplitAlloc
) {
	Result err = STUC_SUCCESS;
	I32 jobCount = 0;
	SplitInPiecesJobArgs jobArgs[MAX_SUB_MAPPING_JOBS] = { 0 };
	makeJobArgs(
		pBasic,
		&jobCount, jobArgs, sizeof(SplitInPiecesJobArgs),
		pInPieces,
		inPiecesJobsGetRange, splitInPiecesJobInit
	);
	err = stucDoJobInParallel(
		pBasic,
		jobCount, jobArgs, sizeof(SplitInPiecesJobArgs),
		splitInPieces
	);
	STUC_RETURN_ERR_IFNOT(err, "");

	inPieceArrDestroy(pBasic->pCtx, pInPieces);
	*pInPieces = (InPieceArr) {0};

	appendNewPiecesToArr(pBasic, pInPiecesSplit, jobCount, jobArgs, getNewInPieces);
	appendNewPiecesToArr(pBasic, pInPiecesSplitClip, jobCount, jobArgs, getNewInPiecesClip);

	for (I32 i = 0; i < jobCount; ++i) {
		pSplitAlloc->pArr[i] = jobArgs[i].alloc;
		inPieceArrDestroy(pBasic->pCtx, &jobArgs[i].newInPieces);
		inPieceArrDestroy(pBasic->pCtx, &jobArgs[i].newInPiecesClip);
	}
	pSplitAlloc->count = jobCount;
	return err;
}

static
Result inPieceArrInit(
	MapToMeshBasic *pBasic,
	InPieceArr *pInPieces,
	I32 *pJobCount, FindEncasedFacesJobArgs *pJobArgs,
	bool *pEmpty
) {
	Result err = STUC_SUCCESS;
	makeJobArgs(
		pBasic,
		pJobCount, pJobArgs, sizeof(FindEncasedFacesJobArgs),
		NULL,
		encasedTableJobsGetRange, NULL
	);
	err = stucDoJobInParallel(
		pBasic,
		*pJobCount, pJobArgs, sizeof(FindEncasedFacesJobArgs),
		stucFindEncasedFaces
	);
	STUC_RETURN_ERR_IFNOT(err, "");

	linkEncasedTableEntries(
		pBasic,
		*pJobCount, pJobArgs,
		pInPieces,
		pEmpty
	);
	return err;
}

typedef struct SnapJobArgs {
	JobArgs core;
	LinAlloc *pIntersectAlloc;
	const InPieceArr *pInPieces;
	const InPieceArr *pInPiecesClip;
	I32 snappedCount;
} SnapJobArgs;

typedef struct SnapJobInitInfo {
	const InPieceArr *pInPieces;
	const InPieceArr *pInPiecesClip;
	HTable *pMergeTable;
} SnapJobInitInfo;

static
I32 snapJobsGetRange(const MapToMeshBasic *pBasic, void *pInitInfoVoid) {
	SnapJobInitInfo *pInitInfo = pInitInfoVoid;
	const LinAlloc *pIntersectAlloc = stucHTableAllocGet(pInitInfo->pMergeTable, 1);
	return stucLinAllocGetCount(pIntersectAlloc);
}

static
void snapJobInit(MapToMeshBasic *pBasic, void *pInitInfoVoid, void *pEntryVoid) {
	SnapJobArgs *pEntry = pEntryVoid;
	SnapJobInitInfo *pInitInfo = pInitInfoVoid;
	pEntry->pIntersectAlloc = stucHTableAllocGet(pInitInfo->pMergeTable, 1);
	pEntry->pInPieces = pInitInfo->pInPieces;
	pEntry->pInPiecesClip = pInitInfo->pInPiecesClip;
}

static
void snapIntersectVert(const SnapJobArgs *pArgs, VertMergeIntersect *pVert) {
	//snap
}

static
Result snapIntersectVertsInRange(void *pArgsVoid) {
	Result err = STUC_SUCCESS;
	SnapJobArgs *pArgs = pArgsVoid;
	LinAllocIter iter = {0};
	stucLinAllocIterInit(pArgs->pIntersectAlloc, pArgs->core.range, &iter);
	for (; !stucLinAllocIterAtEnd(&iter); stucLinAllocIterInc(&iter)) {
		VertMergeIntersect *pVert = stucLinAllocGetItem(&iter);
		snapIntersectVert(pArgs, pVert);
	}
	return err;
}

static
Result snapIntersectVerts(
	MapToMeshBasic *pBasic,
	const InPieceArr *pInPieces,
	const InPieceArr *pInPiecesClip,
	HTable *pMergeTable,
	I32 *pSnappedVerts
) {
	Result err = STUC_SUCCESS;
	I32 jobCount = 0;
	SnapJobArgs jobArgs[MAX_SUB_MAPPING_JOBS] = {0};
	makeJobArgs(
		pBasic,
		&jobCount, jobArgs, sizeof(SnapJobArgs),
		&(SnapJobInitInfo) {
			.pInPieces = pInPieces,
			.pInPiecesClip = pInPiecesClip,
			.pMergeTable = pMergeTable
		},
		snapJobsGetRange, snapJobInit);
	err = stucDoJobInParallel(
		pBasic,
		jobCount, jobArgs, sizeof(SnapJobArgs),
		snapIntersectVertsInRange
	);
	STUC_RETURN_ERR_IFNOT(err, "");
	*pSnappedVerts = 0;
	for (I32 i = 0; i < jobCount; ++i) {
		*pSnappedVerts += jobArgs[i].snappedCount;
	}
	return err;
}

typedef struct XformVertsJobInitInfo {
	const InPieceArr *pInPieces;
	const InPieceArr *pInPiecesClip;
	HTable *pMergeTable;
	I32 vertAllocIdx;
} XformVertsJobInitInfo;

static
I32 xformVertsJobsGetRange(const MapToMeshBasic *pBasic, void *pInitInfoVoid) {
	XformVertsJobInitInfo *pInitInfo = pInitInfoVoid;
	LinAlloc *pVertAlloc =
		stucHTableAllocGet(pInitInfo->pMergeTable, pInitInfo->vertAllocIdx);
	return stucLinAllocGetCount(pVertAlloc);
}

static
void xformVertsJobInit(MapToMeshBasic *pBasic, void *pInitInfoVoid, void *pEntryVoid) {
	xformAndInterpVertsJobArgs *pEntry = pEntryVoid;
	XformVertsJobInitInfo *pInitInfo = pInitInfoVoid;
	LinAlloc *pVertAlloc =
		stucHTableAllocGet(pInitInfo->pMergeTable, pInitInfo->vertAllocIdx);
	pEntry->pVertAlloc = pVertAlloc;
	pEntry->pOutMesh = &pBasic->outMesh;
	pEntry->pInPieces = pInitInfo->pInPieces;
	pEntry->pInPiecesClip = pInitInfo->pInPiecesClip;
	 //TODO again, make an enum or something for lin-alloc handles
	pEntry->intersect = pInitInfo->vertAllocIdx == 1;
}

static
Result xformAndInterpVerts(
	MapToMeshBasic *pBasic,
	const InPieceArr *pInPieces,
	const InPieceArr *pInPiecesClip,
	HTable *pMergeTable,
	I32 vertAllocIdx
) {
	Result err = STUC_SUCCESS;
	I32 jobCount = 0;
	xformAndInterpVertsJobArgs jobArgs[MAX_SUB_MAPPING_JOBS] = {0};
	makeJobArgs(
		pBasic,
		&jobCount, jobArgs, sizeof(xformAndInterpVertsJobArgs),
		&(XformVertsJobInitInfo) {
			.pInPieces = pInPieces,
			.pInPiecesClip = pInPiecesClip,
			.pMergeTable = pMergeTable,
			.vertAllocIdx = vertAllocIdx
		},
		xformVertsJobsGetRange, xformVertsJobInit
	);
	err = stucDoJobInParallel(
		pBasic,
		jobCount, jobArgs, sizeof(xformAndInterpVertsJobArgs),
		stucXformAndInterpVertsInRange
	);
	STUC_RETURN_ERR_IFNOT(err, "");
	return err;
}

typedef struct InterpAttribsJobInitInfo {
	const InPieceArr *pInPieces;
	const InPieceArr *pInPiecesClip;
	const HTable *pMergeTable;
	StucDomain domain;
} InterpAttribsJobInitInfo;

static
I32 interpAttribsJobsGetRange(const MapToMeshBasic *pBasic, void *pInitInfo) {
	return stucDomainCountGetIntern(
		&pBasic->outMesh.core,
		((InterpAttribsJobInitInfo *)pInitInfo)->domain
	);
}

static
void interpAttribsJobInit(MapToMeshBasic *pBasic, void *pInitInfoVoid, void *pEntryVoid) {
	InterpAttribsJobArgs *pEntry = pEntryVoid;
	InterpAttribsJobInitInfo *pInitInfo = pInitInfoVoid;
	pEntry->pOutMesh = &pBasic->outMesh;
	pEntry->pInPieces = pInitInfo->pInPieces;
	pEntry->pInPiecesClip = pInitInfo->pInPiecesClip;
	pEntry->pMergeTable = pInitInfo->pMergeTable;
}

static
Result interpAttribs(
	MapToMeshBasic *pBasic,
	const InPieceArr *pInPieces,
	const InPieceArr *pInPiecesClip,
	HTable *pMergeTable,
	StucDomain domain,
	Result (* job)(void *)
) {
	Result err = STUC_SUCCESS;
	I32 jobCount = 0;
	InterpAttribsJobArgs jobArgs[MAX_SUB_MAPPING_JOBS] = {0};
	makeJobArgs(
		pBasic,
		&jobCount, jobArgs, sizeof(InterpAttribsJobArgs),
		&(InterpAttribsJobInitInfo) {
			.pInPieces = pInPieces,
			.pInPiecesClip = pInPiecesClip,
			.pMergeTable = pMergeTable,
			.domain = domain
		},
		interpAttribsJobsGetRange, interpAttribsJobInit
	);
	err = stucDoJobInParallel(
		pBasic,
		jobCount, jobArgs, sizeof(InterpAttribsJobArgs),
		job
	);
	STUC_RETURN_ERR_IFNOT(err, "");
	return err;
}

static
void addVertsToOutMesh(
	MapToMeshBasic *pBasic,
	HTable *pMergeTable,
	I32 vertAllocIdx
) {
	LinAlloc *pVertAlloc = stucHTableAllocGet(pMergeTable, vertAllocIdx);
	LinAllocIter iter = {0};
	stucLinAllocIterInit(pVertAlloc, (Range) {0, INT32_MAX}, &iter);
	for (; !stucLinAllocIterAtEnd(&iter); stucLinAllocIterInc(&iter)) {
		VertMerge *pEntry = stucLinAllocGetItem(&iter);
		if (vertAllocIdx == 1) { //TODO make an enum for vert alloc types
			VertMergeIntersect *pIntersect = (VertMergeIntersect *)pEntry;
			if (pIntersect->pSnapTo) {
				continue; //vert has been snapped to another - skip
			}
		}
		I32 vert = stucMeshAddVert(pBasic->pCtx, &pBasic->outMesh, NULL);
		pEntry->outVert = vert;
	}
}

typedef struct OutCornerBuf {
	I32Arr buf;
	I32Arr final;
} OutCornerBuf;

static
void addBufFaceToOutMesh(
	MapToMeshBasic *pBasic,
	OutCornerBuf *pOutBuf,
	const InPiece *pInPiece,
	const BufMesh *pBufMesh,
	HTable *pMergeTable,
	I32 faceIdx
) {
	const StucAlloc *pAlloc = &pBasic->pCtx->alloc;
	pOutBuf->buf.count = 0;
	BufFace bufFace = pBufMesh->faces.pArr[faceIdx];
	for (I32 i = 0; i < bufFace.size; ++i) {
		STUC_ASSERT("", pOutBuf->buf.count <= pOutBuf->buf.size);
		if (pOutBuf->buf.count == pOutBuf->buf.size) {
			pOutBuf->buf.size *= 2;
			pOutBuf->buf.pArr =
				pAlloc->fpRealloc(pOutBuf->buf.pArr, pOutBuf->buf.size * sizeof(I32));
		}
		FaceCorner bufCorner = {.face = faceIdx, .corner = i};
		MergeTableKey key = { 0 };
		mergeTableGetVertKey(pBasic, pInPiece, pBufMesh, bufCorner, &key);
		VertMerge *pEntry = NULL;
		SearchResult result = stucHTableGet(
			pMergeTable,
			0,
			&key,
			(void **)&pEntry,
			false, NULL,
			mergeTableMakeKey, NULL, NULL, mergeTableEntryCmp
		);
		STUC_ASSERT("", pEntry && result == STUC_SEARCH_FOUND);
		BufVertType type = bufMeshGetType(pBufMesh, bufCorner);
		if (type == STUC_BUF_VERT_INTERSECT) {
			while (
				pEntry->key.type == STUC_BUF_VERT_INTERSECT &&
				((VertMergeIntersect *)pEntry)->pSnapTo
			) {
				pEntry = ((VertMergeIntersect *)pEntry)->pSnapTo;
			}
		}
		//using lin-idx (vert merge table) for now,
		//this is to allow easy access during fac eand corner attrib interp.
		//Will be set to out-vert later on.
		I32 linIdx = pEntry->linIdx;
		if (pEntry->key.type == STUC_BUF_VERT_INTERSECT) {
			linIdx |= -0x80000000;
		}
		pOutBuf->buf.pArr[pOutBuf->buf.count] = linIdx;
		pOutBuf->buf.count++;
	}
	pOutBuf->final.count = 0;
	for (I32 i = 0; i < pOutBuf->buf.count; ++i) {
		I32 vert = pOutBuf->buf.pArr[i];
		I32 iPrev = i ? i - 1 : pOutBuf->buf.count - 1;
		if (vert == pOutBuf->buf.pArr[iPrev]) {
			continue;
		}
		//not a dup, add
		STUC_ASSERT("", pOutBuf->final.count <= pOutBuf->final.size);
		if (pOutBuf->final.count == pOutBuf->final.size) {
			pOutBuf->final.size *= 2;
			pOutBuf->final.pArr = pAlloc->fpRealloc(
				pOutBuf->final.pArr,
				pOutBuf->final.size * sizeof(I32)
			);
		}
		pOutBuf->final.pArr[pOutBuf->final.count] = vert;
		pOutBuf->final.count++;
	}
	if (pOutBuf->final.count < 3) {
		return; //skip face
	}
	I32 outFace = stucMeshAddFace(pBasic->pCtx, &pBasic->outMesh, NULL);
	pBasic->outMesh.core.pFaces[outFace] = pBasic->outMesh.core.cornerCount;
	for (I32 i = 0; i < pOutBuf->final.count; ++i) {
		I32 outCorner = stucMeshAddCorner(pBasic->pCtx, &pBasic->outMesh, NULL);
		pBasic->outMesh.core.pCorners[outCorner] = pOutBuf->final.pArr[i];
	}
}

static
void addFacesAndCornersToOutMesh(
	MapToMeshBasic *pBasic,
	const InPieceArr *pInPieces,
	HTable *pMergeTable
) {
	OutCornerBuf outBuf = {.buf.size = 16, .final.size = 16};
	outBuf.buf.pArr = pBasic->pCtx->alloc.fpMalloc(outBuf.buf.size * sizeof(I32));
	outBuf.final.pArr = pBasic->pCtx->alloc.fpMalloc(outBuf.final.size * sizeof(I32));
	for (I32 i = 0; i < pInPieces->pBufMeshes->count; ++i) {
		const BufMesh *pBufMesh = pInPieces->pBufMeshes->arr + i;
		for (I32 j = 0; j < pBufMesh->faces.count; ++j) {
			addBufFaceToOutMesh(
				pBasic,
				&outBuf,
				pInPieces->pArr + pBufMesh->faces.pArr[j].inPiece,
				pBufMesh,
				pMergeTable,
				j
			);
		}
	}
	if (outBuf.buf.pArr) {
		pBasic->pCtx->alloc.fpFree(outBuf.buf.pArr);
	}
	if (outBuf.final.pArr) {
		pBasic->pCtx->alloc.fpFree(outBuf.final.pArr);
	}
}

static
Result initOutMesh(MapToMeshBasic *pBasic, HTable *pMergeTable, I32 snappedVerts) {
	Result err = STUC_SUCCESS;
	const StucAlloc *pAlloc = &pBasic->pCtx->alloc;
	Mesh *pMesh = &pBasic->outMesh;
	pMesh->core.type.type = STUC_OBJECT_DATA_MESH_BUF;
	LinAlloc *pVertAlloc = stucHTableAllocGet(pMergeTable, 0);
	LinAlloc *pIntersectVertAlloc = stucHTableAllocGet(pMergeTable, 1);
	I32 bufVertTotal = stucLinAllocGetCount(pVertAlloc);
	bufVertTotal += stucLinAllocGetCount(pIntersectVertAlloc);
	bufVertTotal -= snappedVerts;
	pMesh->faceBufSize = bufVertTotal;
	pMesh->cornerBufSize = pMesh->faceBufSize * 2;
	pMesh->edgeBufSize = pMesh->cornerBufSize;
	pMesh->vertBufSize = pMesh->faceBufSize;
	pMesh->core.pFaces = pAlloc->fpMalloc(sizeof(I32) * pMesh->faceBufSize);
	pMesh->core.pCorners = pAlloc->fpMalloc(sizeof(I32) * pMesh->cornerBufSize);
	//pMesh->core.pEdges = pAlloc->fpMalloc(sizeof(I32) * pMesh->edgeBufSize);

	//in-mesh is the active src,
	// unmatched active map attribs will not be marked active
	const Mesh *srcs[2] = {pBasic->pInMesh, pBasic->pMap->pMesh};
	err = stucAllocAttribsFromMeshArr(
		pBasic->pCtx,
		pMesh,
		2,
		srcs,
		0,
		true, true, false
	);
	STUC_THROW_IFNOT(err, "", 0);
	err = stucAssignActiveAliases(
		pBasic->pCtx,
		pMesh,
		STUC_ATTRIB_USE_FIELD(((StucAttribUse[]) {
			STUC_ATTRIB_USE_POS,
			STUC_ATTRIB_USE_UV,
			STUC_ATTRIB_USE_NORMAL
		})),
		STUC_DOMAIN_NONE
	);
	STUC_CATCH(0, err,
		stucMeshDestroy(pBasic->pCtx, &pMesh->core);
	;);
	return err;
}

static
void tPieceVertInit(
	void *pUserData,
	HTableEntryCore *pEntryCore,
	const void *pKey,
	void *pInitInfo,
	I32 linIdx
) {
	TPieceVert *pEntry = (TPieceVert *)pEntryCore;
	pEntry->vert = *(I32 *)pKey;
}

static
bool tPieceVertCmp(
	const HTableEntryCore *pEntryCore,
	const void *pKey,
	const void *pInitInfo
) {
	TPieceVert *pEntry = (TPieceVert *)pEntryCore;
	return pEntry->vert == *(I32 *)pKey;
}

static
U64 tPieceVertMakeKey(const void *pKey) {
	return *(I32 *)pKey;
}

typedef struct TPieceVertSearch {
	TPieceVert *pEntry;
	SearchResult result;
} TPieceVertSearch;


static
void getLowestTPiece(
	const TPieceBufArr *pTPieces,
	const TPieceVertSearch *pEntries,
	FaceRange *pFace,
	I32 *pLowestTPiece
) {
	for (I32 i = 0; i < pFace->size; ++i) {
		if (pEntries[i].result != STUC_SEARCH_FOUND) {
			continue;
		}
		const TPieceBuf *pTPiece = pTPieces->pArr + pEntries[i].pEntry->tPiece;
		I32 tPieceIdx = pEntries[i].pEntry->tPiece;
		while (pTPieces->pArr[tPieceIdx].merged) {
			tPieceIdx = pTPieces->pArr[tPieceIdx].mergedWith;
		}
		STUC_ASSERT(
			"if t-piece is merged, it should be with a piece of a lower idx",
			tPieceIdx <= pEntries[i].pEntry->tPiece
		);
		if (tPieceIdx < *pLowestTPiece) {
			*pLowestTPiece = tPieceIdx;
		}
	}
}

static
void setVertsAndMergeTPieces(
	TPieceBufArr *pTPieces,
	TPieceVertSearch *pEntries,
	FaceRange *pFace,
	I32 tPiece
) {
	for (I32 i = 0; i < pFace->size; ++i) {
		if (pEntries[i].result == STUC_SEARCH_FOUND &&
			pEntries[i].pEntry->tPiece != tPiece
		) {
			STUC_ASSERT(
				"we only merge to pieces with a lower idx",
				pEntries[i].pEntry->tPiece > tPiece
			);
			TPieceBuf *pTPiece = pTPieces->pArr + pEntries[i].pEntry->tPiece;
			pTPiece->merged = true;
			pTPiece->mergedWith = tPiece;
			//even if the vert's current tpiece was merged,
			// updating the vert entry reduces time spent walking merge chains
			pEntries[i].pEntry->tPiece = tPiece;
		}
		else if (pEntries[i].result == STUC_SEARCH_ADDED) {
			pEntries[i].pEntry->tPiece = tPiece;
		}
	}
}

static
void addOrMergeFaceTPieces(
	const MapToMeshBasic *pBasic,
	TPieceBufArr *pTPieces,
	HTable *pVertTable,
	I32 faceIdx,
	bool add
) {
	const StucAlloc *pAlloc = &pBasic->pCtx->alloc;
	const StucMesh *pMesh = &pBasic->pInMesh->core;
	FaceRange face = stucGetFaceRange(&pBasic->pInMesh->core, faceIdx);
	TPieceVertSearch vertEntries[4] = {0};
	for (I32 i = 0; i < face.size; ++i) {
		vertEntries[i].result = stucHTableGet(
			pVertTable,
			0,
			pMesh->pCorners + face.start + i,
			&vertEntries[i].pEntry,
			add, NULL,
			tPieceVertMakeKey, NULL, tPieceVertInit, tPieceVertCmp
		);
	}
	I32 lowestTPiece = INT32_MAX;
	getLowestTPiece(pTPieces, vertEntries, &face, &lowestTPiece);
	I32 tPiece = -1;
	if (lowestTPiece == INT32_MAX) {
		if (!add) {
			return; //no entries were found for this face
		}
		//all entries are new, so append new tPiece to arr
		STUC_DYN_ARR_ADD(TPieceBuf, pBasic, pTPieces, tPiece);
		pTPieces->pArr[tPiece] = (TPieceBuf) {0};
	}
	else {
		tPiece = lowestTPiece;
	}
	STUC_ASSERT("", tPiece >= 0 && tPiece < pTPieces->count);
	setVertsAndMergeTPieces(pTPieces, vertEntries, &face, tPiece);
}

static
void buildTPiecesForBufVerts(
	const MapToMeshBasic *pBasic,
	const InPieceArr *pInPieces, const InPieceArr *pInPiecesClip,
	LinAlloc *pMergeAlloc,
	TPieceBufArr *pTPieces,
	HTable *pVertTable,
	bool *pChecked
) {
	LinAllocIter iter = {0};
	stucLinAllocIterInit(pMergeAlloc, (Range) {0, INT32_MAX}, &iter);
	for (; !stucLinAllocIterAtEnd(&iter); stucLinAllocIterInc(&iter)) {
		VertMerge *pEntry = stucLinAllocGetItem(&iter);
		const InPiece *pInPiece = NULL;
		const BufMesh *pBufMesh = NULL;
		stucGetBufMeshForVertMergeEntry(
			pInPieces, pInPiecesClip,
			pEntry,
			&pInPiece,
			&pBufMesh
		);
		SrcFaces srcFaces = stucGetSrcFacesForBufCorner(
			pBasic,
			pInPiece,
			pBufMesh,
			pEntry->bufCorner.corner
		);
		STUC_ASSERT(
			"",
			srcFaces.in >= 0 && srcFaces.in < pBasic->pInMesh->core.faceCount
		);
		if (!pChecked[srcFaces.in]) {
			addOrMergeFaceTPieces(pBasic, pTPieces, pVertTable, srcFaces.in, true);
			pChecked[srcFaces.in] = true;
		}
	}
}

static
void buildTPieces(
	const MapToMeshBasic *pBasic,
	const InPieceArr *pInPieces, const InPieceArr *pInPiecesClip,
	HTable *pMergeTable,
	TPieceArr *pTPieces
) {
	LinAlloc *pMergeAlloc = stucHTableAllocGet(pMergeTable, 0);
	LinAlloc *pMergeAllocIntersect = stucHTableAllocGet(pMergeTable, 1);
	HTable vertTable = { 0 };
	stucHTableInit(
		&pBasic->pCtx->alloc,
		&vertTable,
		stucLinAllocGetCount(pMergeAlloc) + stucLinAllocGetCount(pMergeAllocIntersect),
		(I32Arr) {.pArr = (I32[]){sizeof(TPieceVert)}, .count = 1},
		NULL
	);
	bool *pChecked =
		pBasic->pCtx->alloc.fpCalloc(pBasic->pInMesh->core.faceCount, sizeof(bool));
	TPieceBufArr tPiecesBuf = {0};
	buildTPiecesForBufVerts(
		pBasic,
		pInPieces, pInPiecesClip,
		pMergeAlloc,
		&tPiecesBuf,
		&vertTable,
		pChecked
	);
	buildTPiecesForBufVerts(
		pBasic,
		pInPieces, pInPiecesClip,
		pMergeAllocIntersect,
		&tPiecesBuf,
		&vertTable,
		pChecked
	);
	STUC_ASSERT("map-to-mesh should have returned earlier if empty", tPiecesBuf.pArr);
	//not adding new t-pieces, only merging existing ones
	const StucMesh *pInMesh = &pBasic->pInMesh->core;
	for (I32 i = 0; i < pInMesh->faceCount; ++i) {
		if (!pChecked[i]) {
			addOrMergeFaceTPieces(pBasic, &tPiecesBuf, &vertTable, i, false);
			pChecked[i] = true;
		}
	}
	pBasic->pCtx->alloc.fpFree(pChecked);
	
	for (I32 i = 0; i < pInMesh->faceCount; ++i) {
		FaceRange face = stucGetFaceRange(pInMesh, i);
		for (I32 j = 0; j < face.size; ++j) {
			TPieceVert *pEntry = NULL;
			SearchResult result = stucHTableGet(
				&vertTable,
				0,
				&pInMesh->pCorners[face.start + j],
				&pEntry,
				false, NULL,
				tPieceVertMakeKey, NULL, NULL, tPieceVertCmp
			);
			if (result == STUC_SEARCH_NOT_FOUND) {
				continue;
			}
			I32 bufIdx = pEntry->tPiece;
			while (tPiecesBuf.pArr[bufIdx].merged) {
				bufIdx = tPiecesBuf.pArr[bufIdx].mergedWith;
			}
			//update so future searches don't need to walk merge chain
			pEntry->tPiece = bufIdx;
			TPieceBuf *pBuf = tPiecesBuf.pArr + bufIdx;
			if (!pBuf->added) {
				STUC_DYN_ARR_ADD(TPiece, pBasic, pTPieces, pBuf->idx);
				pTPieces->pArr[pBuf->idx] = (TPiece) {0};
				pBuf->added = true;
			}
			STUC_ASSERT("", pBuf->idx >= 0 && pBuf->idx < pTPieces->count);
			I32 faceArrIdx = -1;
			STUC_DYN_ARR_ADD(
				TPieceInFace,
				pBasic,
				(&pTPieces->pArr[pBuf->idx].inFaces),
				faceArrIdx
			);
			STUC_ASSERT("", faceArrIdx != -1);
			pTPieces->pArr[pBuf->idx].inFaces.pArr[faceArrIdx].idx = i;
			pTPieces->pArr[pBuf->idx].inFaces.pArr[faceArrIdx].size = face.size;
			pTPieces->faceCount++;
		}
	}
	stucHTableDestroy(&vertTable);
	pBasic->pCtx->alloc.fpFree(tPiecesBuf.pArr);
}

static
I32 tangentJobGetRange(const MapToMeshBasic *pBasic, void *pInitInfo) {
	return ((TPieceArr *)pInitInfo)->faceCount;
}

static
void tangentJobInit(MapToMeshBasic *pBasic, void *pInitInfo, void *pEntryVoid) {
	TangentJobArgs *pEntry = pEntryVoid;
	pEntry->pTPieces = (TPieceArr *)pInitInfo;
}

static
void copyTangentsFromJobFaces(
	Mesh *pInMesh,
	const TPieceArr *pTPieces,
	const TangentJobArgs *pArgs
) {
	for (I32 i = 0; i < pArgs->faces.count; ++i) {
		I32 jobFaceStart = pArgs->faces.pArr[i];
		I32 inFaceIdx = pTPieces->pInFaces[pArgs->core.range.start + i];
		FaceRange inFace = stucGetFaceRange(&pInMesh->core, inFaceIdx);
		for (I32 j = 0; j < inFace.size; ++j) {
			pInMesh->pTangents[inFace.start + j] = pArgs->pTangents[jobFaceStart + j];
			pInMesh->pTSigns[inFace.start + j] = pArgs->pTSigns[jobFaceStart + j];
		}
	}
}

static
Result buildTangentsForInPieces(
	MapToMeshBasic *pBasic,
	Mesh *pInMesh, //in-mesh is const in MapToMeshBasic
	const InPieceArr *pInPieces, const InPieceArr *pInPiecesClip,
	HTable *pMergeTable
) {
	Result err = STUC_SUCCESS;
	TPieceArr tPieces = {0};
	buildTPieces(pBasic, pInPieces, pInPiecesClip, pMergeTable, &tPieces);
	STUC_ASSERT("", tPieces.pArr);
	I32 jobCount = tPieces.count; //max jobs
	TangentJobArgs jobArgs[MAX_SUB_MAPPING_JOBS] = {0};
	makeJobArgs(
		pBasic,
		&jobCount,
		jobArgs, sizeof(TangentJobArgs),
		&tPieces,
		tangentJobGetRange, tangentJobInit
	);
	tPieces.pInFaces = pBasic->pCtx->alloc.fpCalloc(tPieces.faceCount, sizeof(I32));
	{
		tPieces.faceCount = 0;
		I32 job = 0;
		for (I32 i = 0; i < tPieces.count; ++i) {
			if (job < jobCount - 1 && tPieces.faceCount >= jobArgs[job].core.range.end) {
				jobArgs[job].core.range.end = tPieces.faceCount;
				job++;
				jobArgs[job].core.range.start = tPieces.faceCount;
			}
			STUC_ASSERT("", tPieces.pArr[i].inFaces.pArr);
			for (I32 j = 0; j < tPieces.pArr[i].inFaces.count; ++j) {
				I32 faceJobLocal = -1;
				STUC_DYN_ARR_ADD(I32, pBasic, (&jobArgs[job].faces), faceJobLocal);
				STUC_ASSERT("", faceJobLocal != -1);
				jobArgs[job].faces.pArr[faceJobLocal] = jobArgs[job].cornerCount;

				TPieceInFace face = tPieces.pArr[i].inFaces.pArr[j];
				tPieces.pInFaces[tPieces.faceCount] = face.idx;
				tPieces.faceCount++;
				jobArgs[job].cornerCount += face.size;
			}
			pBasic->pCtx->alloc.fpFree(tPieces.pArr[i].inFaces.pArr);
		}
		//last job may not match jobcount depending on num faces in each t-piece,
		//so update that here
		jobArgs[job].core.range.end = tPieces.faceCount;
		jobCount = job + 1;
	}
	pBasic->pCtx->alloc.fpFree(tPieces.pArr);
	tPieces = (TPieceArr) {.pInFaces = tPieces.pInFaces, .faceCount = tPieces.faceCount};
	err = stucDoJobInParallel(
		pBasic,
		jobCount, jobArgs, sizeof(TangentJobArgs),
		stucBuildTangents
	);
	STUC_THROW_IFNOT(err, "", 0);
	for (I32 i = 0; i < jobCount; ++i) {
		copyTangentsFromJobFaces(pInMesh, &tPieces, jobArgs + i);
	}
	STUC_CATCH(0, err, ;);
	for (I32 i = 0; i < jobCount; ++i) {
		pBasic->pCtx->alloc.fpFree(jobArgs[i].faces.pArr);
		pBasic->pCtx->alloc.fpFree(jobArgs[i].pTangents);
		pBasic->pCtx->alloc.fpFree(jobArgs[i].pTSigns);
	}
	pBasic->pCtx->alloc.fpFree(tPieces.pInFaces);
	return err;
}

static
Result mapToMeshInternal(
	StucContext pCtx,
	const StucMap pMap,
	Mesh *pMeshIn,
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
	printf("A\n");
	if (pInFaceTable) {
		stucLinAllocInit(
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
	FindEncasedFacesJobArgs findEncasedJobArgs[MAX_SUB_MAPPING_JOBS] = {0};
	err = inPieceArrInit(
		&basic,
		&inPieceArr,
		&findEncasedJobCount, findEncasedJobArgs,
		&empty
	);
	STUC_RETURN_ERR_IFNOT(err, "");
	printf("B\n");
	if (!empty) {
		BufMeshArr bufMeshes = {0};
		BufMeshArr bufMeshesClip = {0};
		InPieceArr inPiecesSplit = {.pBufMeshes = &bufMeshes};
		InPieceArr inPiecesSplitClip = {.pBufMeshes = &bufMeshesClip};
		SplitInPiecesAllocArr splitAlloc = {
			.pArr = (SplitInPiecesAlloc[MAX_SUB_MAPPING_JOBS]){0}
		};
		err = inPieceArrSplit(
			&basic,
			&inPieceArr,
			&inPiecesSplit, &inPiecesSplitClip,
			&splitAlloc
		);
		for (I32 i = 0; i < findEncasedJobCount; ++i) {
			LinAlloc *pEncasedAlloc =
				stucHTableAllocGet(&findEncasedJobArgs[i].encasedFaces, 0);
			LinAllocIter iter = {0};
			stucLinAllocIterInit(pEncasedAlloc, (Range) {0, INT32_MAX}, &iter);
			for (; !stucLinAllocIterAtEnd(&iter); stucLinAllocIterInc(&iter)) {
				EncasedMapFace *pEntry = stucLinAllocGetItem(&iter);
				if (pEntry->inFaces.pArr) {
					basic.pCtx->alloc.fpFree(pEntry->inFaces.pArr);
					pEntry->inFaces.pArr = NULL;
				}
			}
			stucHTableDestroy(&findEncasedJobArgs[i].encasedFaces);
		}
		STUC_RETURN_ERR_IFNOT(err, "");
		printf("C\n");
		
		err = inPieceArrInitBufMeshes(&basic, &inPiecesSplitClip, stucClipMapFace);
		STUC_RETURN_ERR_IFNOT(err, "");
		err = inPieceArrInitBufMeshes(&basic, &inPiecesSplit, stucAddMapFaceToBufMesh);
		STUC_RETURN_ERR_IFNOT(err, "");
		printf("D\n");

		HTable mergeTable = {0};
		vertMergeTableInit(&basic, &inPiecesSplit, &inPiecesSplitClip, &mergeTable);
		mergeVerts(&basic, &inPiecesSplit, false, &mergeTable);
		mergeVerts(&basic, &inPiecesSplitClip, true, &mergeTable);
		printf("E\n");

		I32 snappedVerts = 0;
		err = snapIntersectVerts(
			&basic,
			&inPiecesSplit, &inPiecesSplitClip,
			&mergeTable,
			&snappedVerts
		);
		STUC_RETURN_ERR_IFNOT(err, "");
		printf("F\n");

		initOutMesh(&basic, &mergeTable, snappedVerts);
		addVertsToOutMesh(&basic, &mergeTable, 0);
		addVertsToOutMesh(&basic, &mergeTable, 1);//intersect verts
		addFacesAndCornersToOutMesh(&basic, &inPiecesSplit, &mergeTable);
		addFacesAndCornersToOutMesh(&basic, &inPiecesSplitClip, &mergeTable);
		stucMeshSetLastFace(pCtx, &basic.outMesh);
		printf("G\n");

		err = buildTangentsForInPieces(
			&basic,
			pMeshIn,
			&inPiecesSplit, &inPiecesSplitClip,
			&mergeTable
		);
		STUC_RETURN_ERR_IFNOT(err, "");
		printf("H\n");
		
		err = xformAndInterpVerts(&basic, &inPiecesSplit, &inPiecesSplitClip, &mergeTable, 0);
		STUC_RETURN_ERR_IFNOT(err, "");
		//intersect verts
		err = xformAndInterpVerts(&basic, &inPiecesSplit, &inPiecesSplitClip, &mergeTable, 1);
		STUC_RETURN_ERR_IFNOT(err, "");
		err = interpAttribs(
			&basic,
			&inPiecesSplit, &inPiecesSplitClip,
			&mergeTable,
			STUC_DOMAIN_FACE, stucInterpFaceAttribs
		);
		STUC_RETURN_ERR_IFNOT(err, "");
		//vert merge lin-idx is replaced with out-vert idx in corner-interp job,
		// so faces must be interpolated before corners
		err = interpAttribs(
			&basic,
			&inPiecesSplit, &inPiecesSplitClip,
			&mergeTable,
			STUC_DOMAIN_CORNER, stucInterpCornerAttribs
		);
		STUC_RETURN_ERR_IFNOT(err, "");
		printf("I\n");

		for (I32 i = 0; i < splitAlloc.count; ++i) {
			if (splitAlloc.pArr[i].encased.valid) {
				stucLinAllocDestroy(&splitAlloc.pArr[i].encased);
			}
			if (splitAlloc.pArr[i].inFace.valid) {
				stucLinAllocDestroy(&splitAlloc.pArr[i].inFace);
			}
			if (splitAlloc.pArr[i].border.valid) {
				stucLinAllocDestroy(&splitAlloc.pArr[i].border);
			}
		}
		stucHTableDestroy(&mergeTable);
		inPieceArrDestroy(pCtx, &inPiecesSplit);
		inPieceArrDestroy(pCtx, &inPiecesSplitClip);
		bufMeshArrDestroy(pCtx, &bufMeshes);
		bufMeshArrDestroy(pCtx, &bufMeshesClip);
		printf("J\n");

		stucReallocMeshToFit(pCtx, &basic.outMesh);
		*pOutMesh = basic.outMesh.core;
		printf("K\n");
	}
	STUC_CATCH(0, err, ;);
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
		case STUC_ATTRIB_ORIGIN_COMMON: {
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
		}
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
	pCtx->alloc.fpFree(pTable);
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
	const AttribIndexed **ppAttribs = pCtx->alloc.fpCalloc(pMapArr->count, sizeof(void *));
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
	pCtx->alloc.fpFree(ppAttribs);
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
		pAlloc->fpCalloc(pOutIndexedAttribs->size, sizeof(AttribIndexed));
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
	return STUC_SUCCESS;
}

static
Result mapMapArrToMesh(
	StucContext pCtx,
	const StucMapArr *pMapArr,
	Mesh *pMeshIn,
	const StucAttribIndexedArr *pInIndexedAttribs,
	StucMesh *pMeshOut,
	StucAttribIndexedArr *pOutIndexedAttribs,
	F32 wScale,
	F32 receiveLen
) {
	Result err = STUC_SUCCESS;
	Mesh *pOutBufArr = pCtx->alloc.fpCalloc(pMapArr->count, sizeof(Mesh));
	StucObject *pOutObjWrapArr =
		pCtx->alloc.fpCalloc(pMapArr->count, sizeof(StucObject));
	for (I32 i = 0; i < pMapArr->count; ++i) {
		pOutObjWrapArr[i].pData = (StucObjectData *)&pOutBufArr[i];
		const StucMap pMap = pMapArr->ppArr[i];
		I8 matIdx = pMapArr->pMatArr[i];
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
			STUC_THROW_IFNOT(err, "failed to create usg quadtree", 0);
			StucMesh squaresOut = { 0 };
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
			pMapArr->pCommonAttribArr + i,
			NULL,
			wScale,
			receiveLen
		);
		STUC_THROW_IFNOT(err, "map to mesh failed", 1);
		STUC_CATCH(1, err, ;);
		if (pMap->usgArr.count) {
			pCtx->alloc.fpFree(pMap->usgArr.pInFaceTable);
			pMap->usgArr.pInFaceTable = NULL;
			stucLinAllocDestroy(&inFaceTable.alloc);
			inFaceTable = (InFaceTable) {0};
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
	pCtx->alloc.fpFree(pOutBufArr);
	pCtx->alloc.fpFree(pOutObjWrapArr);
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
	STUC_RETURN_ERR_IFNOT(err, "");

	//err = stucBuildTangents(pWrap);
	STUC_RETURN_ERR_IFNOT(err, "failed to build tangents");
	buildEdgeAdj(pWrap);
	buildSeamAndPreserveTables(pWrap);

	//set sp
	stucSetAttribCopyOpt(pCtx, &pWrap->core, STUC_ATTRIB_DONT_COPY, spAttribsToAppend);
	//set required
	stucSetAttribCopyOpt(
		pCtx,
		&pWrap->core,
		STUC_ATTRIB_COPY,
		STUC_ATTRIB_USE_FIELD(((StucAttribUse[]) {
			STUC_ATTRIB_USE_POS,
			STUC_ATTRIB_USE_UV,
			STUC_ATTRIB_USE_NORMAL,
			STUC_ATTRIB_USE_WSCALE,
			STUC_ATTRIB_USE_IDX
		}))
	);

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
			pCtx->alloc.fpFree(meshInWrap.core.pEdges);
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
	pCtx->alloc.fpFree(pUsgArr);
	STUC_CATCH(0, err, ;)
	return err;
}

StucResult stucMeshDestroy(StucContext pCtx, StucMesh *pMesh) {
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

typedef struct RenderArgs {
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
	I8 triCorners[4] = {0};
	V3_F32 bc = stucGetBarycentricInFaceFromVerts(pMesh, pFace, triCorners, *pPos);
	if (bc.d[0] < -.0001f || bc.d[1] < -.0001f || bc.d[2] < -.0001f ||
		!v3F32IsFinite(bc)) {
		return;
	}
	V3_F32 vertsXyz[3] = {0};
	for (I32 i = 0; i < 3; ++i) {
		I32 vertIdx = pMesh->core.pCorners[pFace->start + triCorners[i]];
		vertsXyz[i] = pMesh->pPos[vertIdx];
	}
	V3_F32 wsPos = stucBarycentricToCartesian(vertsXyz, bc);
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
	const StucAlloc *pAlloc = &pArgs->pCtx->alloc;
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
		if (face.size > 4) {
			I32 triCount = face.size - 2;
			FaceTriangulated tris = {.pTris = pAlloc->fpCalloc(triCount * 3, 1)};
			stucTriangulateFaceFromVerts(pAlloc, &face, pMesh, &tris);
			for (I32 l = 0; l < triCount; ++l) {
				FaceRange tri = {0};
				tri.idx = face.idx;
				tri.start = l * 3;
				tri.end = tri.start + 3;
				tri.size = tri.end - tri.start;
				tri.start += face.start;
				tri.end += face.start;
				testPixelAgainstFace(pArgs, pPos, &tri, pColor);
			}
			pAlloc->fpFree(tris.pTris);
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
	args.imageBuf.pData = args.pCtx->alloc.fpMalloc(dataLen);
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
	pThreadPool->fpMutexLock(args.pCtx->pThreadPoolHandle, args.pMutex);
	--*args.pActiveJobs;
	pThreadPool->fpMutexUnlock(args.pCtx->pThreadPoolHandle, args.pMutex);
	return err;
}

static
V2_F32 getZBounds(StucMap pMap) {
	const Mesh *pMesh = pMap->pMesh;
	V2_F32 zBounds = {.d = {FLT_MAX, FLT_MIN}};
	for (I32 i = 0; i < pMesh->core.vertCount; ++i) {
		if (pMesh->pPos[i].d[2] < zBounds.d[0]) {
			zBounds.d[0] = pMesh->pPos[i].d[2];
		}
		if (pMesh->pPos[i].d[2] > zBounds.d[1]) {
			zBounds.d[1] = pMesh->pPos[i].d[2];
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
	pCtx->threadPool.fpMutexGet(pCtx->pThreadPoolHandle, &pMutex);
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
	void **ppJobHandles = pCtx->alloc.fpCalloc(activeJobs, sizeof(void *));
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
	pCtx->alloc.fpFree(ppJobHandles);
	STUC_THROW_IFNOT(err, "", 0);
	I32 pixelSize = getPixelSize(pImage->type);
	pImage->pData = pCtx->alloc.fpMalloc(pixelCount * pixelSize);
	for (I32 i = 0; i < activeJobs; ++i) {
		void *pImageOffset = offsetImagePtr(pImage, i * pixelsPerJob);
		I32 bytesToCopy = pixelSize * args[i].pixelCount;
		memcpy(pImageOffset, args[i].imageBuf.pData, bytesToCopy);
	}
	STUC_CATCH(0, err, ;);
	pCtx->threadPool.fpMutexDestroy(pCtx->pThreadPoolHandle, pMutex);
	for (I32 i = 0; i < activeJobs; ++i) {
		pCtx->alloc.fpFree(args[i].imageBuf.pData);
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
	return pCtx->threadPool.fpWaitForJobs(
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
		err = pCtx->threadPool.fpGetJobErr(
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
		pCtx->threadPool.fpJobHandleDestroy(
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
				pCtx->alloc.fpFree(pArr->pArr[i].core.pData);
			}
		}
		pCtx->alloc.fpFree(pArr->pArr);
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
				pCtx->alloc.fpFree(pMapArr->pCommonAttribArr[i].mesh.pArr);
			}
			if (pMapArr->pCommonAttribArr[i].face.pArr) {
				pCtx->alloc.fpFree(pMapArr->pCommonAttribArr[i].face.pArr);
			}
			if (pMapArr->pCommonAttribArr[i].corner.pArr) {
				pCtx->alloc.fpFree(pMapArr->pCommonAttribArr[i].corner.pArr);
			}
			if (pMapArr->pCommonAttribArr[i].edge.pArr) {
				pCtx->alloc.fpFree(pMapArr->pCommonAttribArr[i].edge.pArr);
			}
			if (pMapArr->pCommonAttribArr[i].vert.pArr) {
				pCtx->alloc.fpFree(pMapArr->pCommonAttribArr[i].vert.pArr);
			}
		}
		pCtx->alloc.fpFree(pMapArr->pCommonAttribArr);
		pMapArr->pCommonAttribArr = NULL;
	}
	if (pMapArr->pMatArr) {
		pCtx->alloc.fpFree(pMapArr->pMatArr);
		pMapArr->pMatArr = NULL;
	}
	if (pMapArr->ppArr) {
		pCtx->alloc.fpFree(pMapArr->ppArr);
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
		pObj->transform = STUC_IDENT_MAT4X4;
	}
	return err;
}
