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

Vec2 vec2Multiply(Vec2 a, Vec2 b);
void vec2MultiplyEqual(Vec2 *a, Vec2 b);
Vec2 vec2DivideScalar(Vec2 a, float b);
void vec2DivideEqualScalar(Vec2 *a, float b);
Vec2 vec2Subtract(Vec2 a, Vec2 b);
Vec2 vec2SubtractScalar(Vec2 a, float b);
Vec2 vec2Add(Vec2 a, Vec2 b);
Vec2 vec2AddScalar(Vec2 a, float b);
void vec2AddEqual(Vec2 *a, Vec2 b);
void vec2MultiplyEqualScalar(Vec2 *a, float b);
Vec2 vec2MultiplyScalar(Vec2 a, float b);
float vec2Dot(Vec2 a, Vec2 b);
Vec2 vec2Cross(Vec2 a);
int32_t vec2GreaterThan(Vec2 a, Vec2 b);
int32_t vec2LessThan(Vec2 a, Vec2 b);
Vec3 cartesianToBarycentric(Vec3 *triVert0, Vec3 *triVert1, Vec3 *triVert2,
                            Vec3 *point);
Vec3 barycentricToCartesian(Vec3 *triVert0, Vec3 *triVert1, Vec3 *triVert2,
                            Vec3 *point);

#define V2MUL ,Multiply,
#define V2MULEQL ,MultiplyEqual,
#define V2DIVS ,DivideScalar,
#define V2DIVSEQL ,DivideEqualScalar,
#define V2SUB ,Subtract,
#define V2SUBS ,SubtractScalar,
#define V2ADD ,Add,
#define V2ADDS ,AddScalar,
#define V2ADDEQL ,AddEqual,
#define V2MULSEQL ,MultiplyEqualScalar,
#define V2MULS ,MultiplyScalar,
#define V2DOT ,Dot,
#define V2GREAT ,GreaterThan,
#define V2LESS ,LessThan,
#define INFIX(a,o,b) vec2##o((a),(b))
#define _(a) INFIX(a)
