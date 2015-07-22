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

static void rna_RefinerSettings_type_set(PointerRNA *ptr, int value)
{
	PartRefine *part_refine = (PartRefine *) ptr->data;

	part_refine->refine_type = value;
/*
	if (!particle_id_check(ptr)) {
		Object *ob = (Object *)ptr->id.data;
		ob->pd->forcefield = value;
		if (ELEM(value, PFIELD_WIND, PFIELD_VORTEX)) {
			ob->empty_drawtype = OB_SINGLE_ARROW;
		}
		else {
			ob->empty_drawtype = OB_PLAINAXES;
		}
	}*/
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

	static EnumPropertyItem refiner_type_items[] = {
	    {0, "NONE", 0, "None", ""},
		{REFINE_OBJ, "REFINER", ICON_PARTICLES, "Refiner", "Refiner for classical SPH adaptive resolution"},
	    {0, NULL, 0, NULL, NULL}
	};

	srna = RNA_def_struct(brna, "RefinerSettings", NULL);
	RNA_def_struct_sdna(srna, "PartRefine");
	RNA_def_struct_path_func(srna, "rna_RefinerSettings_path");
	RNA_def_struct_ui_text(srna, "Refiner Settings", "Refiner settings for an object in physics simulation");
	RNA_def_struct_ui_icon(srna, ICON_PARTICLES);

	/* Enums */

	prop = RNA_def_property(srna, "type", PROP_ENUM, PROP_NONE);
	RNA_def_property_enum_sdna(prop, NULL, "refine_type");
	RNA_def_property_enum_items(prop, refiner_type_items);
	RNA_def_property_enum_funcs(prop, NULL, "rna_RefinerSettings_type_set", NULL);
	RNA_def_property_ui_text(prop, "Refiner", "Is object a refiner");
	//RNA_def_property_update(prop, 0, "rna_RefinerSettings_dependency_update");
}

void RNA_def_object_refine(BlenderRNA *brna)
{
	rna_def_refiner(brna);
}

#endif
