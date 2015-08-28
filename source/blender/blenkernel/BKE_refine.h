#ifndef __BKE_REFINE_H__
#define __BKE_REFINE_H__

/** \file BKE_refine.h
 *  \ingroup bke
 */

struct PartRefine *object_add_refiner(int type);
struct ListBase *prInitRefiners(struct Scene *scene, struct Object *ob_src);
void prEndRefiners(struct ListBase **refiners);

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
