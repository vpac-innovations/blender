#ifndef __DNA_OBJECT_REFINER_H__
#define __DNA_OBJECT_REFINER_H__

#ifdef __cplusplus
extern "C" {
#endif

typedef enum RefineType {
	REFINE_NULL = 0,
	REFINE_POINT = 1,
	REFINE_POINTS = 2,
	REFINE_EDGES = 3,
	REFINE_FACES = 4,
} RefineType;

typedef enum SplitRatio {
	SPLIT2 = 2,
	SPLIT3 = 3,
	SPLIT9 = 9,
} SplitRatio;

typedef enum NSplits {
	REFINE_ONCE = 1,
	REFINE_TWICE = 2,
	REFINE_THREE_TIMES = 3,
	REFINE_FOUR_TIMES = 4,
	REFINE_FIVE_TIMES = 5
} NSplits;

typedef struct PartRefine {
	float radius;
	float max_mass;
	float min_mass;
	float falloff_grad;
	float falloff_offset;
	int split_ratio;
	int falloff_flag;
	short refine_type;
	short nsplits;
} PartRefine;

/* PartRefine->falloff_flag */
#define NO_FALLOFF	0
#define USE_FALLOFF	1

#ifdef __cplusplus
}
#endif

#endif

