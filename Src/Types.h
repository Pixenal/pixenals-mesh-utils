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
	Vec3 pos;
} Vert;

typedef struct {
	int32_t vert;
	Vec3 normal;
} Loop;

typedef struct {
	int32_t loopAmount;
	Loop loops[4];
} Face;

Vec2 vec2Multiply(Vec2 a, Vec2 b);
void vec2MultiplyEqual(Vec2 *a, Vec2 b);
Vec2 vec2DivideScalar(Vec2 a, float b);
void vec2DivideEqualScalar(Vec2 *a, float b);
Vec2 vec2Subtract(Vec2 a, Vec2 b);
Vec2 vec2Add(Vec2 a, Vec2 b);
void vec2AddEqual(Vec2 *a, Vec2 b);
void vec2MultiplyEqualScalar(Vec2 *a, float b);
Vec2 vec2MultiplyScalar(Vec2 a, float b);

#define V2MUL ,Multiply,
#define V2MULEQL ,MultiplyEqual,
#define V2DIVS ,DivideScalar,
#define V2DIVSEQL ,DivideEqualScalar,
#define V2SUB ,Subtract,
#define V2ADD ,Add,
#define V2ADDEQL ,AddEqual,
#define V2MULSEQL ,MultiplyEqualScalar,
#define V2MULS ,MultiplyScalar,
#define INFIX(a,o,b) vec2##o((a),(b))
#define _(a) INFIX(a)
