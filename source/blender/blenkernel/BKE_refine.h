#ifndef __BKE_REFINE_H__
#define __BKE_REFINE_H__

/** \file BKE_refine.h
 *  \ingroup bke
 */

struct PartRefine *object_add_refiner(int type);
struct ListBase *prInitRefiners(struct Scene *scene, struct Object *ob_src, struct ParticleSystem *psys_src);
void prEndRefiners(struct ListBase **refiners);
#endif
