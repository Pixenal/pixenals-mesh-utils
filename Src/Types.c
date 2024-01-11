#include "Types.h"

Vec2 vec2Multiply(Vec2 a, Vec2 b) {
	Vec2 c;
	c.x = a.x * b.x;
	c.y = a.y * b.y;
	return c;
}

void vec2MultiplyEqual(Vec2 *a, Vec2 b) {
	a->x *= b.x;
	a->y *= b.y;
}

Vec2 vec2DivideScalar(Vec2 a, float b) {
	Vec2 c;
	c.x = a.x / b;
	c.y = a.y / b;
	return c;
}

void vec2DivideEqualScalar(Vec2 *a, float b) {
	a->x /= b;
	a->y /= b;
}

Vec2 vec2Subtract(Vec2 a, Vec2 b) {
	Vec2 c = {.x = a.x - b.x, .y = a.y - b.y};
	return c;
}

Vec2 vec2Add(Vec2 a, Vec2 b) {
	Vec2 c = {.x = a.x + b.x, .y = a.y + b.y};
	return c;
}

void vec2AddEqual(Vec2 *a, Vec2 b) {
	a->x += b.x;
	a->y += b.y;
}

void vec2MultiplyEqualScalar(Vec2 *a, float b) {
	a->x *= b;
	a->y *= b;
}

Vec2 vec2MultiplyScalar(Vec2 a, float b) {
	Vec2 c = {.x = a.x * b, .y = a.y * b};
	return c;
}

float vec2Dot(Vec2 a, Vec2 b) {
	return a.x * b.x + a.y * b.y;
}
