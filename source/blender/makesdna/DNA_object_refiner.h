#ifndef __DNA_OBJECT_REFINER_H__
#define __DNA_OBJECT_REFINER_H__

#ifdef __cplusplus
extern "C" {
#endif

typedef enum RefineType {
	REFINE_NULL = 0,
	REFINE_POINT  = 1,
	REFINE_SURFACE = 2,
} RefineType;

typedef struct PartRefine {
	short refine_type;
	short shape;
	float radius;
} PartRefine;

/* pr-> shape */
#define REFINE_SHAPE_SPHERE		0
#define REFINE_SHAPE_BOX		1
#define REFINE_SHAPE_FALLOFF	2

#ifdef __cplusplus
}
#endif

#endif

