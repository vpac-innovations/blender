#ifndef __BKE_REFINE_H__
#define __BKE_REFINE_H__

/** \file BKE_refine.h
 *  \ingroup bke
 */

struct PartRefine *object_add_refiner(int type);
struct ListBase *prInitRefiners(struct Scene *scene, struct Object *ob_src, struct ParticleSystem *psys_src);
void prEndRefiners(struct ListBase **refiners);

typedef struct SPHRefiner {
	struct SPHRefiner *next, *prev;
	struct PartRefine *pr;
	struct Object *ob;
	struct Scene *scene;
	float co[3];
	float radius;

	/* TODO: Add some way of flagging point refiner or surface refiner.
	 * Surface refiner will need to be associated with a parent object
	 * and will need to know which mesh surfaces apply refinement.
	 *
	 * Can point refiners be associated with a parent object?
	 *
	 * Expose all/some of this through Python API in order to develop
	 * an automatic feature detection script? */
} SPHRefiner;
#endif
