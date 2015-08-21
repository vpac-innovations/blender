/*
 * ***** BEGIN GPL LICENSE BLOCK *****
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * The Original Code is: all of this file.
 *
 * Contributor(s): none yet.
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/blenkernel/intern/refine.c
 *  \ingroup bke
 */

#include "MEM_guardedalloc.h"

#include "DNA_listBase.h"
#include "DNA_object_refiner.h"
#include "DNA_object_types.h"
#include "DNA_particle_types.h"
#include "DNA_scene_types.h"

#include "BLI_math.h"
#include "BLI_listbase.h"

#include "BKE_modifier.h"
#include "BKE_refine.h"

#include "BPH_sph.h"

PartRefine *object_add_refiner(int type)
{
	PartRefine *pr;

	pr = MEM_callocN(sizeof(PartRefine), "PartRefine");

	pr->refine_type = type;
	pr->radius = 0.05f;
	pr->falloff = 1.f;
	pr->split_ratio = SPLIT2;
	pr->nsplits = REFINE_ONCE;
	pr->falloff_flag = NO_FALLOFF;
	//pr->max_mass = 11.f / (10.f * pow(pr->split_ratio, pr->nsplits));
	//pr->min_mass = 1.f / (pow(pr->split_ratio, pr->nsplits));

	return pr;
}

static SPHRefiner *new_sph_refiner(Scene *scene, Object *ob, PartRefine *pr)
{
	SPHRefiner *refiner = MEM_callocN(sizeof(SPHRefiner), "SPHRefiner");

	refiner->scene = scene;
	refiner->ob = ob;
	refiner->pr = pr;
	refiner->pr->min_mass = 1.f / (pow(refiner->pr->split_ratio, refiner->pr->nsplits));
	refiner->radius = pr->radius;
	copy_v3_v3(refiner->co, ob->loc);
	refiner->surmd = (SurfaceModifierData *)modifiers_findByType(refiner->ob, eModifierType_Surface);

	return refiner;
}

static void add_object_to_refiners(ListBase **refiners, Scene *scene, Object *ob, Object *ob_src)
{
	SPHRefiner *ref = NULL;

	if (ob == ob_src)
		return;

	/* derivedFinal may have been released in a separate function on
	 * another thread, hold here until it exists.
	 * TODO: Find a better solution. */
	while(ob->type == OB_MESH && !ob->derivedFinal){
	}

	if (*refiners == NULL)
		*refiners = MEM_callocN(sizeof(ListBase), "refiners list");

	ref = new_sph_refiner(scene, ob, ob->pr);

	/* make sure imat is up to date */
	invert_m4_m4(ob->imat, ob->obmat);

	BLI_addtail(*refiners, ref);
}

void prEndRefiners(ListBase **refiners)
{
	if (*refiners) {
		BLI_freelistN(*refiners);
		MEM_freeN(*refiners);
		*refiners = NULL;
	}
}

/* returns ListBase handle with objects taking part in refining */
ListBase *prInitRefiners(Scene *scene, Object *ob_src)
{
	Base *base;
	unsigned int layer= ob_src->lay;
	ListBase *refiners = NULL;

	for(base = scene->base.first; base; base = base->next) {
		if ( (base->lay & layer) ) {
			if (base->object->pr && base->object->pr->refine_type)
				add_object_to_refiners(&refiners, scene, base->object, ob_src);
		}
	}

	return refiners;
}
