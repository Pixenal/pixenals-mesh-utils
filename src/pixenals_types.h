/* 
SPDX-FileCopyrightText: 2025 Caleb Dawson
SPDX-License-Identifier: Apache-2.0
*/

#pragma once
#include <stdint.h>
#include <stdbool.h>

typedef struct PixtyRange {
	int32_t start;
	int32_t end;
} PixtyRange;

typedef struct PixtyI32Arr {
	int32_t *pArr;
	int32_t size;
	int32_t count;
} PixtyI32Arr;

typedef struct PixtyI8Arr {
	int8_t *pArr;
	int32_t size;
	int32_t count;
} PixtyI8Arr;


typedef struct PixtyV2_I8 {
	int8_t d[2];
} PixtyV2_I8;

typedef struct PixtyV2_I16 {
	int16_t d[2];
} PixtyV2_I16;

typedef struct PixtyV2_I32 {
	int32_t d[2];
} PixtyV2_I32;

typedef struct PixtyV2_I64 {
	int64_t d[2];
} PixtyV2_I64;

typedef struct PixtyV2_F32 {
	float d[2];
} PixtyV2_F32;

typedef struct PixtyV2_F64 {
	double d[2];
} PixtyV2_F64;

typedef struct PixtyV3_I8 {
	int8_t d[3];
} PixtyV3_I8;

typedef struct PixtyV3_I16 {
	int16_t d[3];
} PixtyV3_I16;

typedef struct PixtyV3_I32 {
	int32_t d[3];
} PixtyV3_I32;

typedef struct PixtyV3_I64 {
	int64_t d[3];
} PixtyV3_I64;

typedef struct PixtyV3_F32 {
	float d[3];
} PixtyV3_F32;

typedef struct PixtyV3_F64 {
	double d[3];
} PixtyV3_F64;

typedef struct PixtyV4_I8 {
	int8_t d[4];
} PixtyV4_I8;

typedef struct PixtyV4_I16 {
	int16_t d[4];
} PixtyV4_I16;

typedef struct PixtyV4_I32 {
	int32_t d[4];
} PixtyV4_I32;

typedef struct PixtyV4_I64 {
	int64_t d[4];
} PixtyV4_I64;

typedef struct PixtyV4_F32 {
	float d[4];
} PixtyV4_F32;

typedef struct PixtyV4_F64 {
	double d[4];
} PixtyV4_F64;

typedef struct PixtyM2x3 {
	float d[2][3];
} PixtyM2x3;

typedef struct PixtyM2x2 {
	float d[2][2];
} PixtyM2x2;

typedef struct PixtyM3x3 {
	float d[3][3];
} PixtyM3x3;

typedef struct PixtyM4x4 {
	float d[4][4];
} PixtyM4x4;
