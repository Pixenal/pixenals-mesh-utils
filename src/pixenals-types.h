/* 
SPDX-FileCopyrightText: 2025 Caleb Dawson
SPDX-License-Identifier: Apache-2.0
*/

#pragma once
#include <stdint.h>

typedef int8_t Pixtype_I8;
typedef int16_t Pixtype_I16;
typedef int32_t Pixtype_I32;
typedef int64_t Pixtype_I64;

typedef uint8_t Pixtype_U8;
typedef uint16_t Pixtype_U16;
typedef uint32_t Pixtype_U32;
typedef uint64_t Pixtype_U64;

typedef float Pixtype_F32;
typedef double Pixtype_F64;

typedef struct Pixtype_Range {
	Pixtype_I32 start;
	Pixtype_I32 end;
} Pixtype_Range;

typedef struct Pixtype_V2_I8 {
	int8_t d[2];
} Pixtype_V2_I8;

typedef struct Pixtype_V2_I16 {
	int16_t d[2];
} Pixtype_V2_I16;

typedef struct Pixtype_V2_I32 {
	int32_t d[2];
} Pixtype_V2_I32;

typedef struct Pixtype_V2_I64 {
	int64_t d[2];
} Pixtype_V2_I64;

typedef struct Pixtype_V2_F32 {
	float d[2];
} Pixtype_V2_F32;

typedef struct Pixtype_V2_F64 {
	double d[2];
} Pixtype_V2_F64;

typedef struct Pixtype_V3_I8 {
	int8_t d[3];
} Pixtype_V3_I8;

typedef struct Pixtype_V3_I16 {
	int16_t d[3];
} Pixtype_V3_I16;

typedef struct Pixtype_V3_I32 {
	int32_t d[3];
} Pixtype_V3_I32;

typedef struct Pixtype_V3_I64 {
	int64_t d[3];
} Pixtype_V3_I64;

typedef struct Pixtype_V3_F32 {
	float d[3];
} Pixtype_V3_F32;

typedef struct Pixtype_V3_F64 {
	double d[3];
} Pixtype_V3_F64;

typedef struct Pixtype_V4_I8 {
	int8_t d[4];
} Pixtype_V4_I8;

typedef struct Pixtype_V4_I16 {
	int16_t d[4];
} Pixtype_V4_I16;

typedef struct Pixtype_V4_I32 {
	int32_t d[4];
} Pixtype_V4_I32;

typedef struct Pixtype_V4_I64 {
	int64_t d[4];
} Pixtype_V4_I64;

typedef struct Pixtype_V4_F32 {
	float d[4];
} Pixtype_V4_F32;

typedef struct Pixtype_V4_F64 {
	double d[4];
} Pixtype_V4_F64;

typedef struct Pixtype_M4x4_F32 {
	float d[4][4];
} Pixtype_M4x4_F32;
