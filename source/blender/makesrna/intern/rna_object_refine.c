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
 * Contributor(s): Blender Foundation (2015), Sean Loh
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/makesrna/intern/rna_object_refine.c
 *  \ingroup RNA
 */

#include <stdlib.h>

#include "DNA_object_types.h"
#include "DNA_object_refiner.h"
#include "DNA_particle_types.h"
#include "DNA_scene_types.h"

#include "RNA_define.h"

#include "rna_internal.h"

#include "WM_api.h"
#include "WM_types.h"

static EnumPropertyItem refiner_type_items[] = {
    {0, "NONE", 0, "None", ""},
    {REFINE_POINT, "POINT", 0, "From point", "Refiner for classical SPH adaptive resolution"},
    {REFINE_POINTS, "MESH_POINTS", 0, "From mesh points", "Refiner for classical SPH adaptive resolution"},
    {REFINE_EDGES, "MESH_EDGES", 0, "From mesh edges", "Refiner for classical SPH adaptive resolution"},
    {REFINE_FACES, "MESH_FACES", 0, "From mesh faces", "Refiner for classical SPH adaptive resolution"},
    {0, NULL, 0, NULL, NULL}
};

static EnumPropertyItem refiner_ratio_items[] = {
    {SPLIT2, "SPLIT2", 0, "2:1", "Split each parent particle into 2 daughter particles"},
    {SPLIT3, "SPLIT3", 0, "3:1", "Split each parent particle into 3 daughter particles"},
    {SPLIT9, "SPLIT9", 0, "9:1", "Split each parent particle into 9 daughter particles"},
    {0, NULL, 0, NULL, NULL}
};

static EnumPropertyItem refiner_split_items[] = {
    {REFINE_ONCE, "ONCE", 0, "1", "Allow particles to be split once."},
    {REFINE_TWICE, "TWICE", 0, "2", "Allow particles to be split twice"},
    {REFINE_THREE_TIMES, "THREE_TIMES", 0, "3", "Allow particles to be split three times"},
    {REFINE_FOUR_TIMES, "FOUR_TIMES", 0, "4", "Allow particles to be split four times"},
    {REFINE_FIVE_TIMES, "FIVE_TIMES", 0, "5", "Allow particles to be split five times"},
    {0, NULL, 0, NULL, NULL}
};
#ifdef RNA_RUNTIME

#include "BLI_math_base.h"

#include "MEM_guardedalloc.h"

#include "DNA_modifier_types.h"
#include "DNA_texture_types.h"

#include "BKE_context.h"
#include "BKE_DerivedMesh.h"
#include "BKE_modifier.h"
#include "BKE_pointcache.h"
#include "BKE_depsgraph.h"

#include "ED_object.h"

static int particle_id_check(PointerRNA *ptr)
{
	ID *id = ptr->id.data;

	return (GS(id->name) == ID_PA);
}

static void rna_RefinerSettings_update(Main *UNUSED(bmain), Scene *UNUSED(scene), PointerRNA *ptr)
{
	Object *ob = (Object *)ptr->id.data;

	DAG_id_tag_update(&ob->id, OB_RECALC_OB);
	WM_main_add_notifier(NC_OBJECT | ND_DRAW, ob);
}

static EnumPropertyItem point_type_items[] = {
	{0, "NONE", 0, "None", ""},
	{REFINE_POINT, "POINT", 0, "From Point", ""},
	{0, NULL, 0, NULL, NULL}
};

static EnumPropertyItem edge_type_items[] = {
	{0, "NONE", 0, "None", ""},
	{REFINE_POINTS, "MESH_POINTS", 0, "From mesh verts", ""},
    {REFINE_EDGES, "MESH_EDGES", 0, "From mesh edges", ""},
	{0, NULL, 0, NULL, NULL}
};
static EnumPropertyItem surface_type_items[] = {
	{0, "NONE", 0, "None", ""},
	{REFINE_POINTS, "MESH_POINTS", 0, "From mesh verts", ""},
    {REFINE_EDGES, "MESH_EDGES", 0, "From mesh edges", ""},
	{REFINE_FACES, "MESH_FACES", 0, "From mesh faces", ""},
	{0, NULL, 0, NULL, NULL}
};
static EnumPropertyItem volume_split_items[] = {
    {REFINE_ONCE, "ONCE", 0, "1", "Allow particles to be split once."},
    {REFINE_TWICE, "TWICE", 0, "2", "Allow particles to be split twice"},
    {REFINE_THREE_TIMES, "THREE_TIMES", 0, "3", "Allow particles to be split three times"},
    {0, NULL, 0, NULL, NULL}
};

static EnumPropertyItem *rna_Refiner_type_itemf(bContext *UNUSED(C), PointerRNA *ptr,
                                                 PropertyRNA *UNUSED(prop), bool *UNUSED(r_free))
{
	Object *ob = NULL;
	DerivedMesh *dm;

	ob = (Object *)ptr->id.data;

	if (ELEM(ob->type, OB_EMPTY))
		return point_type_items;

	if (ELEM(ob->type, OB_MESH)){
		dm = ob->derivedFinal;
		if (dm->getNumTessFaces(dm))
			return surface_type_items;
		else
			return edge_type_items;
	}

	return refiner_type_items;
}

static void rna_RefinerSettings_shape_update(Main *bmain, Scene *scene, PointerRNA *ptr)
{
	if (!particle_id_check(ptr)) {
		Object *ob = (Object *)ptr->id.data;
		ED_object_check_refiner_modifiers(bmain, scene, ob);
		WM_main_add_notifier(NC_OBJECT | ND_DRAW, ob);
		WM_main_add_notifier(NC_OBJECT | ND_MODIFIER, ob);
	}
}

static void rna_RefinerSettings_type_update(Main *bmain, Scene *scene, PointerRNA *ptr)
{
	if (!particle_id_check(ptr)) {
		Object *ob = (Object *)ptr->id.data;
		ED_object_check_refiner_modifiers(bmain, scene, ob);
		WM_main_add_notifier(NC_OBJECT | ND_DRAW, ob);
		WM_main_add_notifier(NC_OBJECT | ND_MODIFIER, ob);
	}
}

static void rna_RefinerSettings_ratio_update(Main *bmain, Scene *scene, PointerRNA *ptr)
{
	if (!particle_id_check(ptr)) {
		Object *ob = (Object *)ptr->id.data;
		ED_object_check_refiner_modifiers(bmain, scene, ob);
		WM_main_add_notifier(NC_OBJECT | ND_DRAW, ob);
		WM_main_add_notifier(NC_OBJECT | ND_MODIFIER, ob);
	}
}

static EnumPropertyItem *rna_Refiner_split_itemf(bContext *UNUSED(C), PointerRNA *ptr,
                                                 PropertyRNA *UNUSED(prop), bool *UNUSED(r_free))
{
	Object *ob = NULL;

	ob = (Object *)ptr->id.data;

	if (ob->pr->split_ratio != SPLIT9)
		return refiner_split_items;
	else
		return volume_split_items;
}

static void rna_RefinerSettings_split_update(Main *bmain, Scene *scene, PointerRNA *ptr)
{
	if (!particle_id_check(ptr)) {
		Object *ob = (Object *)ptr->id.data;
		ED_object_check_refiner_modifiers(bmain, scene, ob);
		WM_main_add_notifier(NC_OBJECT | ND_DRAW, ob);
		WM_main_add_notifier(NC_OBJECT | ND_MODIFIER, ob);
	}
}

static void rna_RefinerSettings_dependency_update(Main *bmain, Scene *scene, PointerRNA *ptr)
{
	if (particle_id_check(ptr)) {
		DAG_id_tag_update((ID *)ptr->id.data, OB_RECALC_OB | OB_RECALC_DATA | OB_RECALC_TIME | PSYS_RECALC_RESET);
	}
	else {
		Object *ob = (Object *)ptr->id.data;

		rna_RefinerSettings_shape_update(bmain, scene, ptr);

		DAG_relations_tag_update(bmain);

		DAG_id_tag_update(&ob->id, OB_RECALC_OB);

		WM_main_add_notifier(NC_OBJECT | ND_DRAW, ob);
	}
}

static char *rna_RefinerSettings_path(PointerRNA *ptr)
{
	PartRefine *pr = (PartRefine *)ptr->data;

	/* Check through all possible places the settings can be to find the right one */
	/* object refiner */
	Object *ob = (Object *)ptr->id.data;

	if (ob->pr == pr)
		return BLI_sprintfN("refiner");

	return NULL;
}
#else
static void rna_def_refiner(BlenderRNA *brna)
{
	StructRNA *srna;
	PropertyRNA *prop;

	srna = RNA_def_struct(brna, "RefinerSettings", NULL);
	RNA_def_struct_sdna(srna, "PartRefine");
	RNA_def_struct_path_func(srna, "rna_RefinerSettings_path");
	RNA_def_struct_ui_text(srna, "Refiner Settings", "Refiner settings for an object in physics simulation");
	RNA_def_struct_ui_icon(srna, ICON_PARTICLES);

	/* Enums */

	prop = RNA_def_property(srna, "type", PROP_ENUM, PROP_NONE);
	RNA_def_property_enum_sdna(prop, NULL, "refine_type");
	RNA_def_property_enum_items(prop, refiner_type_items);
	RNA_def_property_enum_funcs(prop, NULL, NULL, "rna_Refiner_type_itemf");
	RNA_def_property_ui_text(prop, "Refiner", "Measure distance from particle to refiner using objects' verts, edges or faces");
	RNA_def_property_update(prop, 0, "rna_RefinerSettings_dependency_update");
	//RNA_def_property_update(prop, 0, "rna_RefinerSettings_type_update");

	prop = RNA_def_property(srna, "ratio", PROP_ENUM, PROP_NONE);
	RNA_def_property_enum_sdna(prop, NULL, "split_ratio");
	RNA_def_property_enum_items(prop, refiner_ratio_items);
	RNA_def_property_ui_text(prop, "Splitting ratio", "Ratio of daughter to parent particles");
	RNA_def_property_update(prop, 0, "rna_RefinerSettings_ratio_update");

	prop = RNA_def_property(srna, "nsplit", PROP_ENUM, PROP_NONE);
	RNA_def_property_enum_sdna(prop, NULL, "nsplits");
	RNA_def_property_enum_items(prop, refiner_split_items);
	RNA_def_property_enum_funcs(prop, NULL, NULL, "rna_Refiner_split_itemf");
	RNA_def_property_ui_text(prop, "Number of splits allowed", "Max number of times a particle may be split");
	RNA_def_property_update(prop, 0, "rna_RefinerSettings_split_update");

	/* Floats */

	prop = RNA_def_property(srna, "radius", PROP_FLOAT, PROP_NONE);
	RNA_def_property_float_sdna(prop, NULL, "radius");
	RNA_def_property_range(prop, 0.0f, FLT_MAX);
	RNA_def_property_ui_range(prop, 0.0f, FLT_MAX, 1, 4);
	RNA_def_property_ui_text(prop, "Radius", "Radius within which adaptive resolution is active");
	RNA_def_property_update(prop, 0, "rna_RefinerSettings_update");

	prop = RNA_def_property(srna, "minimum_mass", PROP_FLOAT, PROP_NONE);
	RNA_def_property_float_sdna(prop, NULL, "min_mass");
	RNA_def_property_range(prop, 0.0f, 1.f);
	RNA_def_property_ui_text(prop, "Min mass fac", "Lower limit for particle mass in refiner region");
	RNA_def_property_update(prop, 0, "rna_RefinerSettings_update");

	prop = RNA_def_property(srna, "falloff_xo", PROP_FLOAT, PROP_NONE);
	RNA_def_property_float_sdna(prop, NULL, "falloff_xo");
	RNA_def_property_range(prop, 0.0f, FLT_MAX);
	RNA_def_property_ui_text(prop, "xmin", "Distance from refiner where target particle mass is minimum mass");
	RNA_def_property_update(prop, 0, "rna_RefinerSettings_update");

	prop = RNA_def_property(srna, "falloff_xn", PROP_FLOAT, PROP_NONE);
	RNA_def_property_float_sdna(prop, NULL, "falloff_xn");
	RNA_def_property_range(prop, 0.0f, FLT_MAX);
	RNA_def_property_ui_text(prop, "xmax", "Distance from refiner where target particle mass is the unrefined particle mass");
	RNA_def_property_update(prop, 0, "rna_RefinerSettings_update");

	/* Misc?*/
	prop = RNA_def_property(srna, "use_falloff", PROP_BOOLEAN, PROP_NONE);
	RNA_def_property_boolean_sdna(prop, NULL, "falloff_flag", NO_FALLOFF);
	RNA_def_property_ui_text(prop, "Use falloff function",
	                         "Use exponential falloff function to determine when splitting/merging happens");
	//RNA_def_property_update(prop, 0, "rna_Particle_reset");
}

void RNA_def_object_refine(BlenderRNA *brna)
{
	rna_def_refiner(brna);
}

#endif
