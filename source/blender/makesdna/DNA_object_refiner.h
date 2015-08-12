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

typedef struct PartRefine {
	float radius;
	float max_mass;
	float min_mass;
	float falloff;
	short refine_type;
	short shape;
} PartRefine;

/* pr->shape */
#define REFINE_SHAPE_SPHERE		0
#define REFINE_SHAPE_BOX		1
#define REFINE_SHAPE_FALLOFF	2

#ifdef __cplusplus
}
#endif

#endif

