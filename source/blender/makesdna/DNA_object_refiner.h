#ifndef __DNA_OBJECT_REFINER_H__
#define __DNA_OBJECT_REFINER_H__

#ifdef __cplusplus
extern "C" {
#endif

typedef enum RefineType {
	REFINE_NULL = 0,
	REFINE_OBJ  = 1,
} RefineType;

typedef struct PartRefine {
	short refine_type;
	short pad;
	float radius;
} PartRefine;

#ifdef __cplusplus
}
#endif

#endif

