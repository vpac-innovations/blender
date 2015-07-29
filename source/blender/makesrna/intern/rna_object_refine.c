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
    {REFINE_POINT, "POINT", 0, "Point", "Refiner for classical SPH adaptive resolution"},
    {REFINE_SURFACE, "SURFACE", 0, "Surface", "Refiner for classical SPH adaptive resolution"},
    {0, NULL, 0, NULL, NULL}
};
#ifdef RNA_RUNTIME

#include "BLI_math_base.h"

#include "MEM_guardedalloc.h"

#include "DNA_modifier_types.h"
#include "DNA_texture_types.h"

#include "BKE_context.h"
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
	{REFINE_POINT, "POINT", 0, "Point", ""},
	{0, NULL, 0, NULL, NULL}
};

static EnumPropertyItem surface_type_items[] = {
	{0, "NONE", 0, "None", ""},
	{REFINE_POINT, "POINT", 0, "Point", ""},
	{REFINE_SURFACE, "SURFACE", 0, "Surface", ""},
	{0, NULL, 0, NULL, NULL}
};

static EnumPropertyItem *rna_Refiner_type_itemf(bContext *UNUSED(C), PointerRNA *ptr,
                                                 PropertyRNA *UNUSED(prop), bool *UNUSED(r_free))
{
	Object *ob = NULL;

	ob = (Object *)ptr->id.data;

	if (ELEM(ob->type, OB_EMPTY))
		return point_type_items;

	if (ELEM(ob->type, OB_MESH, OB_SURF))
			return surface_type_items;

	return refiner_type_items;
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
/*
static void rna_RefinerSettings_dependency_update(Main *bmain, Scene *scene, PointerRNA *ptr)
{
	if (particle_id_check(ptr)) {
		DAG_id_tag_update((ID *)ptr->id.data, OB_RECALC_OB | OB_RECALC_DATA | OB_RECALC_TIME | PSYS_RECALC_RESET);
	}
	else {
		Object *ob = (Object *)ptr->id.data;

		rna_FieldSettings_shape_update(bmain, scene, ptr);

		DAG_relations_tag_update(bmain);

		DAG_id_tag_update(&ob->id, OB_RECALC_OB);

		WM_main_add_notifier(NC_OBJECT | ND_DRAW, ob);
	}
}*/
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
	RNA_def_property_ui_text(prop, "Refiner", "Is object a refiner");
	RNA_def_property_update(prop, 0, "rna_RefinerSettings_type_update");

	/* Floats */

	prop = RNA_def_property(srna, "radius", PROP_FLOAT, PROP_NONE);
	RNA_def_property_float_sdna(prop, NULL, "radius");
	RNA_def_property_range(prop, 0.0f, FLT_MAX);
	RNA_def_property_ui_text(prop, "Radius", "Radius within which to apply splitting");
	RNA_def_property_update(prop, 0, "rna_RefinerSettings_update");

	prop = RNA_def_property(srna, "maximum_mass", PROP_FLOAT, PROP_NONE);
	RNA_def_property_float_sdna(prop, NULL, "max_mass");
	RNA_def_property_range(prop, 0.0f, FLT_MAX);
	RNA_def_property_ui_text(prop, "Maximum particle mass", "Upper limit for particle mass in refiner region");
	RNA_def_property_update(prop, 0, "rna_RefinerSettings_update");

	prop = RNA_def_property(srna, "minimum_mass", PROP_FLOAT, PROP_NONE);
	RNA_def_property_float_sdna(prop, NULL, "min_mass");
	RNA_def_property_range(prop, 0.0f, FLT_MAX);
	RNA_def_property_ui_text(prop, "Minimum particle mass", "Lower limit for particle mass in refiner region");
	RNA_def_property_update(prop, 0, "rna_RefinerSettings_update");
}

void RNA_def_object_refine(BlenderRNA *brna)
{
	rna_def_refiner(brna);
}

#endif
