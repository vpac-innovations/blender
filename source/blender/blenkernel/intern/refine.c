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
#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_object_refiner.h"
#include "DNA_object_types.h"
#include "DNA_particle_types.h"
#include "DNA_scene_types.h"

#include "BLI_math.h"
#include "BLI_listbase.h"

#include "BKE_bvhutils.h"
#include "BKE_cdderivedmesh.h"
#include "BKE_editmesh.h"
#include "BKE_modifier.h"
#include "BKE_refine.h"

#include "BPH_sph.h"

PartRefine *object_add_refiner(int type)
{
	PartRefine *pr;

	pr = MEM_callocN(sizeof(PartRefine), "PartRefine");

	pr->refine_type = type;
	pr->radius = 0.05f;
	pr->falloff_xn = 1.f;
	pr->falloff_xo = 0.f;
	pr->split_ratio = SPLIT2;
	pr->nsplits = REFINE_ONCE;
	pr->falloff_flag = NO_FALLOFF;

	return pr;
}

void add_refiner_custom_data_layers(Object *ob, int overwrite)
{
	PartRefine *pr = ob->pr;
	MFloatProperty *fp;
	CustomDataLayer *cdl;
	Mesh *me = ob->data;
	int index;
	int i;

	if (pr->refine_type == REFINE_FACES) {
		/* Retrieve existing custom data layer or create if not found */
		index = CustomData_get_named_layer_index(&me->pdata, CD_PROP_FLT, "radius");
		if (index < 0) {
			CustomData_add_layer_named(&me->pdata, CD_PROP_FLT, CD_DEFAULT, NULL, me->totpoly, "radius");
			index = CustomData_get_named_layer_index(&me->pdata, CD_PROP_FLT, "radius");
			overwrite = 1;
		}

		if (overwrite) {
			cdl = &me->pdata.layers[index];
			fp = ((MFloatProperty *)cdl->data);
			for (i = 0; i < me->totpoly; i++) {
				fp[i].f = ob->pr->radius;
			}
		}
	}
}

static SPHRefiner *new_sph_refiner(Scene *scene, Object *ob, PartRefine *pr)
{
	SPHRefiner *refiner = MEM_callocN(sizeof(SPHRefiner), "SPHRefiner");

	refiner->scene = scene;
	refiner->ob = ob;
	refiner->pr = pr;
	refiner->radius = pr->radius;
	copy_v3_v3(refiner->co, ob->loc);
	refiner->surmd = (SurfaceModifierData *)modifiers_findByType(refiner->ob, eModifierType_Surface);

	/* Need to account for non-uniform mass distribution of
	 * 9:1 splitting method */
	if (refiner->pr->split_ratio == SPLIT9)
		refiner->pr->min_mass = 2.f / (refiner->pr->split_ratio + 1.f);
	else
		refiner->pr->min_mass = 1.f / (pow(refiner->pr->split_ratio, refiner->pr->nsplits));

	return refiner;
}

static void add_object_to_refiners(ListBase **refiners, Scene *scene, Object *ob, Object *ob_src)
{
	SPHRefiner *ref = NULL;

	if (ob == ob_src)
		return;

	if (ob->type == OB_MESH) {
	/* Add refiner properties to a custom data layer for
	 * per face/edge/vert customisation. */
		add_refiner_custom_data_layers(ob, 0);
	}

	/* derivedFinal may have been released on another thread,
	 * hold here until it exists.
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

static int dm_tessface_to_poly_index(DerivedMesh *dm, int tessface_index)
{
	if (tessface_index != ORIGINDEX_NONE) {
		/* double lookup */
		const int *index_mf_to_mpoly;
		if ((index_mf_to_mpoly = dm->getTessFaceDataArray(dm, CD_ORIGINDEX))) {
			const int *index_mp_to_orig = dm->getPolyDataArray(dm, CD_ORIGINDEX);
			return DM_origindex_mface_mpoly(index_mf_to_mpoly, index_mp_to_orig, tessface_index);
		}
	}

	return ORIGINDEX_NONE;
}

int closest_point_on_refiner(Object *ob, SurfaceModifierData *surmd, const float co[3], float surface_co[3], float surface_nor[3], int *mp_index)
{
	BVHTreeNearest nearest;

	nearest.index = -1;
	nearest.dist_sq = FLT_MAX;

	BLI_bvhtree_find_nearest(surmd->bvhtree->tree, co, &nearest, surmd->bvhtree->nearest_callback, surmd->bvhtree);

	if (nearest.index != -1) {
	/* derivedFinal may have been released on another thread,
	 * hold here until it exists.
	 * TODO: Find a better solution. */
		while(ob->type == OB_MESH && !ob->derivedFinal){
		}

		*mp_index = dm_tessface_to_poly_index(ob->derivedFinal, nearest.index);
		copy_v3_v3(surface_co, nearest.co);

		if (surface_nor) {
			copy_v3_v3(surface_nor, nearest.no);
		}
		return 1;
	}

	return 0;
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
