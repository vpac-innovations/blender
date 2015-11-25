#ifndef __BKE_REFINE_H__
#define __BKE_REFINE_H__

/** \file BKE_refine.h
 *  \ingroup bke
 */

#include "DNA_modifier_types.h"

struct PartRefine *object_add_refiner(int type);
struct ListBase *prInitRefiners(struct Scene *scene, struct Object *ob_src);
void prEndRefiners(struct ListBase **refiners);
void add_refiner_custom_data_layers(struct Object *ob, int overwrite);
int closest_point_on_refiner(struct Object *ob, SurfaceModifierData *surmd, const float co[3], float surface_co[3], float surface_nor[3], int *mp_index);

typedef struct SPHRefiner {
	struct SPHRefiner *next, *prev;

	struct Scene *scene;
	struct Object *ob;
	struct ParticleSystem *psys;
	struct SurfaceModifierData *surmd;

	struct PartRefine *pr;

	float co[3];
	float nor[3];
	float vec_to_particle[3];
	float radius;\
} SPHRefiner;

typedef struct RefinerData {
	float v2p[3];
	float dist;
	int ratio;
} RefinerData;

#endif
