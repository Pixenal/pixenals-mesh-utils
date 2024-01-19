#pragma once
#include <stdint.h>

typedef struct {
	unsigned char *string;
	int32_t nextBitIndex;
	int32_t byteIndex;
} UvgpByteString;

typedef struct {
	float x;
	float y;
} Vec2;

typedef struct {
	float x;
	float y;
	float z;
} Vec3;

typedef struct {
	int32_t x;
	int32_t y;
} iVec2;

typedef struct BoundaryVert{
	struct BoundaryVert *next;
	int32_t vert;
	int32_t loopAmount;
} BoundaryVert;

typedef struct {
	int32_t vertAmount;
	Vec3 *vertBuffer;
	int32_t loopAmount;
	int32_t *loopBuffer;
	int32_t faceAmount;
	int32_t *faceBuffer;
	Vec2 *uvBuffer;
} BlenderMeshData;

typedef struct Cell {
	uint32_t localIndex;
	uint32_t initialized;
	struct Cell *children;
	int32_t faceAmount;
	int32_t *faces;
	int32_t cellIndex;
	Vec2 boundsMin;
	Vec2 boundsMax;
} Cell;

typedef struct {
	Cell **cells;
	int32_t cellAmount;
	int32_t faceAmount;
} FaceCellsInfo;

//TODO: reduce faces to 4 in order to lower this from 12
typedef struct VertAdj{
	struct VertAdj *next;
	int32_t vert;
	int32_t uvgpVert;
	int32_t loopAmount;
} VertAdj;

typedef struct {
	int32_t id;
	int32_t *jobsCompleted;
	int32_t averageVertAdjDepth;
	int32_t *boundaryFaceStart;
	BoundaryVert *boundaryBuffer;
	BlenderMeshData localMesh;
	int32_t *boundaryVerts;
	int32_t boundaryVertAmount;
	BlenderMeshData mesh;
} ThreadArg;

Vec3 vec3SubtractScalar(Vec3 a, float b);
Vec3 vec3AddScalar(Vec3 a, float b);
int32_t vec3GreaterThan(Vec3 a, Vec3 b);
int32_t vec3LessThan(Vec3 a, Vec3 b);

Vec2 vec2Multiply(Vec2 a, Vec2 b);
void vec2MultiplyEqual(Vec2 *a, Vec2 b);
Vec2 vec2DivideScalar(Vec2 a, float b);
void vec2DivideEqualScalar(Vec2 *a, float b);
Vec2 vec2Subtract(Vec2 a, Vec2 b);
void vec2SubtractEqual(Vec2 *a, Vec2 b);
Vec2 vec2SubtractScalar(Vec2 a, float b);
Vec2 vec2Add(Vec2 a, Vec2 b);
Vec2 vec2AddScalar(Vec2 a, float b);
void vec2AddEqual(Vec2 *a, Vec2 b);
void vec2MultiplyEqualScalar(Vec2 *a, float b);
Vec2 vec2MultiplyScalar(Vec2 a, float b);
float vec2Dot(Vec2 a, Vec2 b);
Vec2 vec2Cross(Vec2 a);
Vec2 vec2ModScalar(Vec2 a, float b);
void vec2ModEqualScalar(Vec2 *a, float b);
int32_t vec2GreaterThan(Vec2 a, Vec2 b);
int32_t vec2LessThan(Vec2 a, Vec2 b);
int32_t vec2LessThanEqualTo(Vec2 a, Vec2 b);
Vec3 cartesianToBarycentric(Vec3 *triVert0, Vec3 *triVert1, Vec3 *triVert2,
                            Vec3 *point);
Vec3 barycentricToCartesian(Vec3 *triVert0, Vec3 *triVert1, Vec3 *triVert2,
                            Vec3 *point);


#define V3SUBS ,3SubtractScalar,
#define V3ADDS ,3AddScalar,
#define V3GREAT ,3GreaterThan,
#define V3LESS ,3LessThan,

#define V2MUL ,2Multiply,
#define V2MULEQL ,2MultiplyEqual,
#define V2DIVS ,2DivideScalar,
#define V2DIVSEQL ,2DivideEqualScalar,
#define V2SUB ,2Subtract,
#define V2SUBEQL ,2SubtractEqual,
#define V2SUBS ,2SubtractScalar,
#define V2ADD ,2Add,
#define V2ADDS ,2AddScalar,
#define V2ADDEQL ,2AddEqual,
#define V2MULSEQL ,2MultiplyEqualScalar,
#define V2MULS ,2MultiplyScalar,
#define V2DOT ,2Dot,
#define V2MODS ,2ModScalar,
#define V2MODEQLS ,2ModEqualScalar,
#define V2GREAT ,2GreaterThan,
#define V2LESS ,2LessThan,
#define V2LESSEQL ,2LessThanEqualTo,
#define INFIX(a,o,b) vec##o((a),(b))
#define _(a) INFIX(a)
